./analysis_script.R 
  The file runs the entire analysis required to generate Figure 4 in the manuscript.
  Note however that the primary data is not included in this archive; hence, the 
  "Data processing" section in the script should be skipped. This section normally 
  generates the four data sets found in /data, which are read at a later point in the 
  analysis instead. It produces four fitted generalized additive models (included in
  ./analysis/fitted_models) which are used to interpolate between the data from the
  simulations (included in analysis/model_predictions). These predicted data are used
  to generate the .PDF and .PNG versions of Figure 4 (included in ./figure_files).

  NOTE: Fitting the GAMs takes an enormous amount of time, as does the generation of
  model predictions. For this reason, fitted GAMs and their predictons have been 
  saved and can be read in as is done in the file. This makes re-analysis MUCH less
  time-consuming.

# GITHUB NOTE: The /analysis folder is omitted due to the size of the files within,
but GAMs and model predictions can be regenerated with the analysis scripts as 
described above.