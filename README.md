# kif_golgi_fragmentation
matlab programs for analyzing Golgi fragmentation based on percentage of large Golgi objects

dv2tif.m converts dv files to individual tifs (requires bio-formats)
  dv_data_extract_metadata.m helper function for dv2tif.matlab

*.cpproj files provide baseline pipelines for cellprofiler analysis of noc or sirna effect

postcp_XXX_analysis.m function analyzes the .mat file produced by the cellprofiler pipeline
  postcp_XXX.m helper function for extracting data from the structure produced by cellprofiler
  thresholdstatistics.m helper function to discard cells not meeting the required thresholds