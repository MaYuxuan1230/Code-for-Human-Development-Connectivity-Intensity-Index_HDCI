#!/usr/bin/env Rscript

required_pkgs <- c(
  "sf", "terra", "dplyr", "data.table",
  "ranger", "foreach", "doParallel", "yaml"
)

to_install <- required_pkgs[!required_pkgs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

invisible(lapply(required_pkgs, library, character.only = TRUE))

sf_use_s2(FALSE)
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript scripts/compute_reservoir_rfidw.R <config.yml>")
}

cfg <- yaml::read_yaml(args[1])

# ---------------------------
# Configuration
# ---------------------------
res_path   <- cfg$paths$reservoir
china_path <- cfg$paths$china_boundary
out_dir    <- cfg$paths$output_dir
tmp_dir    <- cfg$paths$temp_dir

area_field <- cfg$fields$area
stor_field <- cfg$fields$storage

res_m       <- cfg$raster$resolution_m
tile_size_m <- cfg$raster$tile_size_m

n_classes      <- cfg$model$n_classes
num_trees      <- cfg$model$num_trees
alpha_min      <- cfg$model$alpha_min
alpha_max      <- cfg$model$alpha_max
radius_min_km  <- cfg$model$radius_min_km
radius_max_km  <- cfg$model$radius_max_km
seed_value     <- cfg$model$seed

use_parallel    <- isTRUE(cfg$compute$use_parallel)
n_cores         <- cfg$compute$n_cores
terra_memfrac   <- cfg$compute$terra_memfrac
terra_memmax_gb <- cfg$compute$terra_memmax_gb

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

terraOptions(
  progress = 1,
  memfrac = terra_memfrac,
  memmax = terra_memmax_gb,
  tempdir = tmp_dir
)

# ---------------------------
# Helper functions
# ---------------------------
rescale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (isTRUE(all.equal(rng[1], rng[2]))) {
    rep(0.5, length(x))
  } else {
    (x - rng[1]) / (rng[2] - rng[1])
  }
}

winsorize_vec <- function(x, probs = c(0.01, 0.99)) {
  q <- quantile(x, probs = probs, na.rm = TRUE)
  x[x < q[1]] <- q[1]
  x[x > q[2]] <- q[2]
  x
}

make_tile_extents <- function(r, tile_size_m = 300000) {
  ex <- ext(r)
  xmin_r <- ex[1]; xmax_r <- ex[2]
  ymin_r <- ex[3]; ymax_r <- ex[4]

  xs <- unique(c(seq(xmin_r, xmax_r, by = tile_size_m), xmax_r))
  ys <- unique(c(seq(ymin_r, ymax_r, by = tile_size_m), ymax_r))

  tile_exts <- list()
  k <- 1
  for (i in 1:(length(xs) - 1)) {
    for (j in 1:(length(ys) - 1)) {
      tile_exts[[k]] <- c(xs[i], xs[i + 1], ys[j], ys[j + 1])
      k <- k + 1
    }
  }
  tile_exts
}

multiclass_logloss <- function(prob_mat, truth_factor) {
  idx <- cbind(seq_len(nrow(prob_mat)), match(truth_factor, colnames(prob_mat)))
  -mean(log(pmax(prob_mat[idx], 1e-15)))
}

fit_ranger_once <- function(dat, mtry, min.node.size, sample.fraction, num_trees = 1200, seed_value = 1L) {
  tab <- table(dat$size_class)
  class_weights <- as.numeric(median(tab) / tab)
  names(class_weights) <- names(tab)

  rf <- ranger(
    size_class ~ .,
    data = dat,
    probability = TRUE,
    num.trees = num_trees,
    mtry = mtry,
    min.node.size = min.node.size,
    sample.fraction = sample.fraction,
    importance = "permutation",
    class.weights = class_weights,
    seed = seed_value,
    num.threads = max(1, parallel::detectCores() - 1),
    respect.unordered.factors = "order"
  )

  p <- rf$predictions
  pred_class <- colnames(p)[max.col(p, ties.method = "first")]
  acc <- mean(pred_class == dat$size_class)
  ll  <- multiclass_logloss(p, dat$size_class)

  list(model = rf, oob_acc = acc, oob_logloss = ll)
}

run_cv <- function(dat, best_mtry, best_min_node, best_sample_fraction, v = 5, num_trees = 800, seed_value = 1L) {
  set.seed(seed_value)
  folds <- sample(rep(1:v, length.out = nrow(dat)))
  res_list <- vector("list", v)

  for (k in 1:v) {
    train_dat <- dat[folds != k, , drop = FALSE]
    test_dat  <- dat[folds == k, , drop = FALSE]

    tab <- table(train_dat$size_class)
    class_weights <- as.numeric(median(tab) / tab)
    names(class_weights) <- names(tab)

    rf <- ranger(
      size_class ~ .,
      data = train_dat,
      probability = TRUE,
      num.trees = num_trees,
      mtry = best_mtry,
      min.node.size = best_min_node,
      sample.fraction = best_sample_fraction,
      importance = "none",
      class.weights = class_weights,
      seed = seed_value + k,
      num.threads = max(1, parallel::detectCores() - 1),
      respect.unordered.factors = "order"
    )

    p <- predict(rf, data = test_dat)$predictions
    pred_class <- colnames(p)[max.col(p, ties.method = "first")]
    acc <- mean(pred_class == test_dat$size_class)
    ll  <- multiclass_logloss(p, test_dat$size_class)

    res_list[[k]] <- data.frame(fold = k, acc = acc, logloss = ll)
  }

  do.call(rbind, res_list)
}

process_one_tile <- function(tile_id, tile_ext_num, template_file, pts_file, max_radius_m, out_dir, tmp_dir,
                             terra_memfrac = 0.15, terra_memmax_gb = 4) {
  library(terra)

  terraOptions(
    progress = 0,
    memfrac = terra_memfrac,
    memmax = terra_memmax_gb,
    tempdir = tmp_dir
  )

  template <- rast(template_file)
  pts_v <- vect(pts_file)

  tile_ext <- ext(tile_ext_num[1], tile_ext_num[2], tile_ext_num[3], tile_ext_num[4])
  tile <- crop(template, tile_ext, snap = "out")
  if (is.null(tile)) return(NA_character_)

  tile_vals <- values(tile, mat = FALSE)
  if (length(tile_vals) == 0 || all(is.na(tile_vals))) return(NA_character_)

  tile_acc <- tile
  values(tile_acc) <- ifelse(!is.na(tile_vals), 0, NA_real_)

  exb <- ext(
    tile_ext_num[1] - max_radius_m,
    tile_ext_num[2] + max_radius_m,
    tile_ext_num[3] - max_radius_m,
    tile_ext_num[4] + max_radius_m
  )

  pts_tile <- crop(pts_v, exb)

  if (!is.null(pts_tile) && nrow(pts_tile) > 0) {
    for (i in 1:nrow(pts_tile)) {
      d <- distance(tile_acc, pts_tile[i], unit = "m")
      rad <- pts_tile$radius_m[i]
      alp <- pts_tile$alpha[i]
      strg <- pts_tile$strength[i]

      contrib <- ifel(d <= rad, strg / ((1 + d / 1000)^alp), 0)
      tile_acc <- tile_acc + contrib

      rm(d, contrib)
      gc()
    }
  }

  out_file <- file.path(out_dir, sprintf("tile_%04d.tif", tile_id))
  writeRaster(
    tile_acc, out_file,
    overwrite = TRUE,
    wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=LZW"))
  )

  rm(template, pts_v, tile, tile_vals, tile_acc, pts_tile)
  gc()

  out_file
}

# ---------------------------
# Read and prepare data
# ---------------------------
prov_sf <- st_read(china_path, quiet = TRUE)
res_sf  <- st_read(res_path, quiet = TRUE)

if (!all(c(area_field, stor_field) %in% names(res_sf))) {
  stop("Reservoir attribute fields not found. Please check config.yml.")
}

res_sf[[area_field]] <- as.numeric(res_sf[[area_field]])
res_sf[[stor_field]] <- as.numeric(res_sf[[stor_field]])

res_sf <- res_sf %>%
  filter(!is.na(.data[[area_field]]), !is.na(.data[[stor_field]]),
         .data[[area_field]] > 0, .data[[stor_field]] > 0)

target_crs <- "+proj=aea +lat_1=25 +lat_2=47 +lat_0=0 +lon_0=105 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

prov_sf <- st_make_valid(prov_sf)
res_sf  <- st_make_valid(res_sf)

prov_sf <- st_transform(prov_sf, target_crs)
res_sf  <- st_transform(res_sf, target_crs)

china_sf <- st_sf(geometry = st_union(st_geometry(prov_sf)))
china_sf <- st_make_valid(china_sf)

hit <- st_intersects(res_sf, china_sf, sparse = FALSE)[, 1]
res_sf <- res_sf[hit, ]

# ---------------------------
# Feature engineering
# ---------------------------
res_v <- vect(res_sf)

geom_area_m2 <- expanse(res_v, unit = "m")
perimeter_m  <- perim(res_v)

res_sf$geom_area_m2 <- as.numeric(geom_area_m2)
res_sf$perimeter_m  <- as.numeric(perimeter_m)

res_sf$compactness <- pmax(4 * pi * res_sf$geom_area_m2 / (res_sf$perimeter_m^2 + 1e-12), 0)
res_sf$shape_index <- res_sf$perimeter_m / (2 * sqrt(pi * res_sf$geom_area_m2 + 1e-12))
res_sf$stor_area_ratio <- res_sf[[stor_field]] / (res_sf[[area_field]] + 1e-12)

res_sf$Area_w            <- winsorize_vec(res_sf[[area_field]])
res_sf$STOR_w            <- winsorize_vec(res_sf[[stor_field]])
res_sf$geom_area_m2_w    <- winsorize_vec(res_sf$geom_area_m2)
res_sf$perimeter_m_w     <- winsorize_vec(res_sf$perimeter_m)
res_sf$stor_area_ratio_w <- winsorize_vec(res_sf$stor_area_ratio)

res_sf$log_area      <- log1p(res_sf$Area_w)
res_sf$log_stor      <- log1p(res_sf$STOR_w)
res_sf$log_geom_area <- log1p(res_sf$geom_area_m2_w)
res_sf$log_perimeter <- log1p(res_sf$perimeter_m_w)
res_sf$log_ratio     <- log1p(res_sf$stor_area_ratio_w)
res_sf$area_stor_inter <- res_sf$log_area * res_sf$log_stor

feature_names <- c(
  "log_area", "log_stor", "log_geom_area", "log_perimeter",
  "compactness", "shape_index", "log_ratio", "area_stor_inter"
)

feat_df <- st_drop_geometry(res_sf)[, feature_names]

# ---------------------------
# Initial classes by clustering
# ---------------------------
set.seed(seed_value)
km <- kmeans(scale(feat_df), centers = n_classes, nstart = 100, iter.max = 200)

res_sf$cluster_raw <- km$cluster

centers_df <- as.data.frame(km$centers)
centers_df$cluster_raw <- 1:n_classes
centers_df$order_score <- centers_df$log_area + centers_df$log_stor + centers_df$log_geom_area

ord <- centers_df$cluster_raw[order(centers_df$order_score)]
map_ord <- setNames(seq_len(n_classes), ord)

res_sf$size_id <- as.integer(map_ord[as.character(res_sf$cluster_raw)])
class_labels <- c("small", "medium", "large", "mega")[seq_len(n_classes)]
res_sf$size_class <- factor(res_sf$size_id, levels = 1:n_classes, labels = class_labels)

# ---------------------------
# RF tuning and training
# ---------------------------
model_df <- st_drop_geometry(res_sf)[, c("size_class", feature_names)]

grid <- expand.grid(
  mtry = unique(pmax(2, pmin(length(feature_names), c(2, 3, 4, floor(sqrt(length(feature_names))), ceiling(length(feature_names) / 2))))),
  min.node.size = c(1, 5, 10, 20),
  sample.fraction = c(0.60, 0.75, 0.90)
)

tune_list <- vector("list", nrow(grid))
for (i in 1:nrow(grid)) {
  rr <- fit_ranger_once(
    dat = model_df,
    mtry = grid$mtry[i],
    min.node.size = grid$min.node.size[i],
    sample.fraction = grid$sample.fraction[i],
    num_trees = num_trees,
    seed_value = seed_value
  )

  tune_list[[i]] <- data.frame(
    mtry = grid$mtry[i],
    min.node.size = grid$min.node.size[i],
    sample.fraction = grid$sample.fraction[i],
    oob_acc = rr$oob_acc,
    oob_logloss = rr$oob_logloss
  )
}

tune_res <- do.call(rbind, tune_list)
tune_res <- tune_res[order(tune_res$oob_logloss, -tune_res$oob_acc), ]
best_par <- tune_res[1, ]
write.csv(tune_res, file.path(out_dir, "RF_tuning_results.csv"), row.names = FALSE)

tab_final <- table(model_df$size_class)
class_weights_final <- as.numeric(median(tab_final) / tab_final)
names(class_weights_final) <- names(tab_final)

rf_fit <- ranger(
  size_class ~ .,
  data = model_df,
  probability = TRUE,
  num.trees = num_trees,
  mtry = best_par$mtry,
  min.node.size = best_par$min.node.size,
  sample.fraction = best_par$sample.fraction,
  importance = "permutation",
  class.weights = class_weights_final,
  seed = seed_value,
  num.threads = max(1, parallel::detectCores() - 1),
  respect.unordered.factors = "order"
)

oob_prob <- rf_fit$predictions
oob_pred <- colnames(oob_prob)[max.col(oob_prob, ties.method = "first")]
oob_acc  <- mean(oob_pred == model_df$size_class)
oob_ll   <- multiclass_logloss(oob_prob, model_df$size_class)

cv_res <- run_cv(
  dat = model_df,
  best_mtry = best_par$mtry,
  best_min_node = best_par$min.node.size,
  best_sample_fraction = best_par$sample.fraction,
  v = 5,
  num_trees = 800,
  seed_value = seed_value
)
write.csv(cv_res, file.path(out_dir, "RF_5foldCV_results.csv"), row.names = FALSE)

varimp_df <- data.frame(
  variable = names(rf_fit$variable.importance),
  importance = as.numeric(rf_fit$variable.importance)
) %>% arrange(desc(importance))
write.csv(varimp_df, file.path(out_dir, "RF_variable_importance.csv"), row.names = FALSE)

# ---------------------------
# Map RF output to decay parameters
# ---------------------------
pred_prob <- predict(rf_fit, data = model_df)$predictions
ord_num <- seq_len(n_classes)

expected_class <- as.vector(pred_prob %*% ord_num)
size_score_rf <- (expected_class - 1) / (n_classes - 1)

cont_size <- 0.45 * rescale01(res_sf$log_area) +
             0.45 * rescale01(res_sf$log_stor) +
             0.10 * rescale01(res_sf$log_geom_area)

res_sf$size_score <- 0.6 * size_score_rf + 0.4 * cont_size
res_sf$size_score <- pmin(pmax(res_sf$size_score, 0), 1)

res_sf$alpha <- alpha_max - res_sf$size_score * (alpha_max - alpha_min)
res_sf$radius_km <- radius_min_km + res_sf$size_score * (radius_max_km - radius_min_km)
res_sf$radius_m  <- res_sf$radius_km * 1000

res_sf$strength <- 0.3 * rescale01(res_sf$log_area) +
                   0.3 * rescale01(res_sf$log_stor) +
                   0.4 * res_sf$size_score
res_sf$strength <- pmin(pmax(res_sf$strength, 1e-6), 1)

res_out <- res_sf %>%
  select(
    all_of(c(area_field, stor_field)),
    geom_area_m2, perimeter_m,
    size_class, size_score, alpha, radius_km, strength
  )

st_write(
  res_out,
  file.path(out_dir, "reservoir_RF_parameters.gpkg"),
  delete_dsn = TRUE,
  quiet = TRUE
)

# ---------------------------
# Build template raster
# ---------------------------
china_v <- vect(china_sf)

template <- rast(ext(china_v), resolution = res_m, crs = crs(china_v))
values(template) <- 1
template <- mask(template, china_v)
names(template) <- "reservoir_hii_raw"

template_file <- file.path(out_dir, "template_1km.tif")
writeRaster(
  template,
  template_file,
  overwrite = TRUE,
  wopt = list(datatype = "INT1U", gdal = "COMPRESS=LZW")
)

# ---------------------------
# Reservoir representative points
# ---------------------------
pts_sf <- st_point_on_surface(res_sf)
pts_v <- vect(pts_sf)

pts_v$alpha    <- res_sf$alpha
pts_v$radius_m <- res_sf$radius_m
pts_v$strength <- res_sf$strength

pts_file <- file.path(out_dir, "reservoir_RF_points.gpkg")
writeVector(pts_v, pts_file, overwrite = TRUE)

# ---------------------------
# Tile-based interpolation
# ---------------------------
tile_exts <- make_tile_extents(template, tile_size_m = tile_size_m)
max_radius_m <- max(pts_v$radius_m, na.rm = TRUE)

tile_out_dir <- file.path(out_dir, "tiles")
dir.create(tile_out_dir, showWarnings = FALSE)

if (use_parallel) {
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  parallel::clusterEvalQ(cl, {
    library(terra)
    NULL
  })

  tile_files <- foreach(
    i = seq_along(tile_exts),
    .packages = c("terra")
  ) %dopar% {
    process_one_tile(
      tile_id = i,
      tile_ext_num = tile_exts[[i]],
      template_file = template_file,
      pts_file = pts_file,
      max_radius_m = max_radius_m,
      out_dir = tile_out_dir,
      tmp_dir = tmp_dir,
      terra_memfrac = terra_memfrac,
      terra_memmax_gb = terra_memmax_gb
    )
  }

  stopCluster(cl)
} else {
  tile_files <- vector("list", length(tile_exts))
  for (i in seq_along(tile_exts)) {
    tile_files[[i]] <- process_one_tile(
      tile_id = i,
      tile_ext_num = tile_exts[[i]],
      template_file = template_file,
      pts_file = pts_file,
      max_radius_m = max_radius_m,
      out_dir = tile_out_dir,
      tmp_dir = tmp_dir,
      terra_memfrac = terra_memfrac,
      terra_memmax_gb = terra_memmax_gb
    )
  }
}

tile_files <- unlist(tile_files)
tile_files <- tile_files[!is.na(tile_files)]
tile_files <- tile_files[file.exists(tile_files)]

if (length(tile_files) == 0) {
  stop("No valid tiles were produced.")
}

# ---------------------------
# Merge and normalize
# ---------------------------
template_r <- rast(template_file)
tile_rasters <- lapply(tile_files, rast)
s <- sprc(tile_rasters)

raw_r_vrt <- merge(
  s,
  algo = 3,
  first = TRUE,
  filename = file.path(out_dir, "merged_tiles_tmp.vrt"),
  overwrite = TRUE
)

raw_r <- resample(raw_r_vrt, template_r, method = "near")
raw_r <- mask(raw_r, china_v)
names(raw_r) <- "reservoir_hii_raw"

writeRaster(
  raw_r,
  file.path(out_dir, "reservoir_hii_raw_1km.tif"),
  overwrite = TRUE,
  wopt = list(datatype = "FLT4S", gdal = "COMPRESS=LZW")
)

rmin <- global(raw_r, "min", na.rm = TRUE)[1, 1]
rmax <- global(raw_r, "max", na.rm = TRUE)[1, 1]

if (isTRUE(all.equal(rmin, rmax))) {
  stop("Raw raster min and max are identical; normalization cannot be performed.")
}

hii_100 <- 100 * (raw_r - rmin) / (rmax - rmin)
hii_100 <- clamp(hii_100, lower = 0, upper = 100, values = TRUE)
names(hii_100) <- "reservoir_hii_0_100"

writeRaster(
  hii_100,
  file.path(out_dir, "reservoir_hii_0_100_1km.tif"),
  overwrite = TRUE,
  wopt = list(datatype = "FLT4S", gdal = "COMPRESS=LZW")
)

summary_df <- data.frame(
  metric = c(
    "n_reservoirs",
    "rf_oob_accuracy",
    "rf_oob_logloss",
    "cv_mean_accuracy",
    "cv_mean_logloss",
    "raw_min",
    "raw_max"
  ),
  value = c(
    nrow(res_sf),
    oob_acc,
    oob_ll,
    mean(cv_res$acc, na.rm = TRUE),
    mean(cv_res$logloss, na.rm = TRUE),
    rmin,
    rmax
  )
)

write.csv(summary_df, file.path(out_dir, "summary_metrics.csv"), row.names = FALSE)
