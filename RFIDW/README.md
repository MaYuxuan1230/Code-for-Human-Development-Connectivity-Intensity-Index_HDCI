# Reservoir RFIDW (Public Release)

Reproducible R code for building a 1-km reservoir influence surface using a random-forest-assisted inverse distance weighting workflow (RFIDW), based on reservoir polygons and reservoir attributes.

## Inputs

The script expects:

1. A reservoir polygon dataset with at least two numeric fields:
   - `Area`
   - `STOR`
2. A China boundary polygon dataset
3. User-defined output directory

## Outputs

The script writes:

- `reservoir_RF_parameters.gpkg`
- `reservoir_RF_points.gpkg`
- `reservoir_hii_raw_1km.tif`
- `reservoir_hii_0_100_1km.tif`
- `RF_tuning_results.csv`
- `RF_5foldCV_results.csv`
- `RF_variable_importance.csv`
- `summary_metrics.csv`

## Usage

Copy `config.example.yml` to `config.yml`, edit the paths, then run:

```bash
Rscript scripts/compute_reservoir_rfidw.R config.yml
```

## Notes

- The script uses a projected Albers equal-area CRS for distance-based computation.
- The default implementation keeps the original methodological logic but removes manuscript-specific formatting, local hard-coded paths, and non-essential console/report output.
- For reproducibility, key parameters are controlled through the YAML configuration file.

## Dependencies

The script will install missing packages if necessary:

- sf
- terra
- dplyr
- data.table
- ranger
- foreach
- doParallel
- yaml
