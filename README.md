# TMBstocks
Working examples of TMB models for salmon stock assessment

Repository to work on and run TMB models for stock-recruitment relationships. Examples include:
- Basic Ricker stock-recruit model
-...

## Folders:
- data_in: Raw data for models 
- output: Model results
- figures: Figures created by scripts
- TMB: contains .cpp TMB model files

## Files:
- `make.r`: Master file to source functions, load data, run TMB models, save output tables and figures
- `functions.r`: helper functions
- `load.r`: Load required packages and data
- `run_TMB.r`: run TMB model and save results
- `figures.r`: Functions to create figures
