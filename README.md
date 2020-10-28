# fishyTMB
Working examples of TMB models for salmon stock assessment

Repository to work on and run TMB models for stock-recruitment relationships. Examples include:
- Basic Ricker stock-recruit model
- (other examples to follow)

## Files:
### main directory: 
- `ricker.r`: Source functions, load data, run ricker TMB model, save output tables and figures
- `functions.r`: helper functions
- `load_data.r`: Download required data
- `figures.r`: Functions to create figures
### data_in: 
Raw data for models 
### output: 
Model results
### figures: 
Figures created by scripts
### TMB: 
Contains .cpp TMB model files
- `ricker.cpp`: simple Ricker model 
- `Aggregate_LRPs.cpp`: Combined Ricker and LRP model (Brooke Davis)
- `Aggregate_LRPs_Hier.cpp`: Combined Ricker and LRP model (Brooke Davis)


