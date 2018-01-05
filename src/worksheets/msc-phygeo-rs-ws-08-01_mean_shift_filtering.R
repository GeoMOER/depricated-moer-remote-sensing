# rs-ws-08-1
# MOC - Data Analysis (T. Nauss, C. Reudenbach)
# Mean-shift filtering

# Set environment --------------------------------------------------------------
if(Sys.info()["sysname"] == "Windows"){
  source("D:/active/moc/msc-ui/scripts/msc-phygeo-ei/src/functions/set_environment.R")
} else {
  source("/media/permanent/active/moc/msc-ui/scripts/msc-phygeo-ei/src/functions/set_environment.R")
}


# Scale dataset ----------------------------------------------------------------
# rgb_ind_files <- list.files(path_muf_set1m_sample_segm,
#                             full.names = TRUE, pattern = glob2rx("*.tif"))
# muf <- stack(rgb_ind_files)
# 
# new_min = 0
# new_max = 255
# 
# muf_scaled <- lapply(seq(nlayers(muf)), function(i){
#   mfx <- getValues(muf[[i]])
#   old_min <- min(mfx)
#   old_max <- max(mfx)
#   
#   mfx_scaled <- new_min + (new_max - new_min) * (mfx - old_min) / (old_max - old_min)  
#   mfx_scaled <- setValues(muf[[i]], mfx_scaled)
#   return(mfx_scaled)
# })
# muf_scaled <- stack(muf_scaled)
# 
# writeRaster(muf_scaled, paste0(path_muf_set1m_sample_segm, "ortho_muf_rgb_idx_pca_scaled.tif"))


# Mean shift filtering ---------------------------------------------------------
spatialr_set = c(5)
range_set = c(15, 30)

muf <- paste0(path_muf_set1m_sample_segm, "ortho_muf_rgb_idx_pca_scaled.tif")

for(spatialr in spatialr_set){
  for(ranger in range_set){
    outfile_filter <- (paste0(path_muf_set1m_sample_segm, "ortho_muf_rgb_idx_pca_scaled_filter_r", 
                              spatialr, "_rng", ranger, ".tif"))
    
    outfile_spatial <- (paste0(path_muf_set1m_sample_segm, "ortho_muf_rgb_idx_pca_scaled_spatial_r", 
                               spatialr, "_rng", ranger, ".tif"))
    
    muf_msf <- otbcli_MeanShiftSmoothing(x = muf,
                                         outfile_filter = outfile_filter,
                                         outfile_spatial = outfile_spatial,
                                         return_raster = FALSE,
                                         spatialr = spatialr,
                                         ranger = ranger,
                                         thres = 0.1,
                                         maxiter = 100,
                                         rangeramp = 0,
                                         verbose=FALSE,
                                         ram="8192")
  }
}


# Segmentation -----------------------------------------------------------------
filter_files <- list.files(path_muf_set1m_sample_segm, 
                           pattern = glob2rx("ortho_muf_rgb_idx_pca_scaled_filter_r*.tif"),
                           full.names = TRUE)
spatial_files <- gsub("filter", "spatial", filter_files)


spatialr_set = c(1)
range_set = c(15, 30)

for(i in seq(length(filter_files))){
  print(i)
  for(spatialr in spatialr_set){
    for(ranger in range_set){
      outfile <- filter_files[[i]]
      outfile <- paste0(substr(outfile, 1, nchar(outfile)-4), "_segm_r", spatialr, "_rng", ranger, ".tif")
      muf_lsmss <- otbcli_ExactLargeScaleMeanShiftSegmentation(x = filter_files[[i]],
                                                               inpos = spatial_files[[i]],
                                                               out = outfile,
                                                               tmpdir = path_temp,
                                                               spatialr = spatialr,
                                                               ranger = ranger,
                                                               minsize = 0,
                                                               tilesizex = 1000,
                                                               tilesizey = 1000,
                                                               verbose = FALSE,
                                                               return_raster = FALSE)
    }
  }
}



# Aggregation of areas below minimum size --------------------------------------
segm_files <- list.files(path_muf_set1m_sample_segm, 
                         pattern = glob2rx("ortho_muf_rgb_idx_pca_scaled*segm*.tif"),
                         full.names = TRUE)
filter_files <- paste0(dirname(segm_files), "/",
                       substr(basename(segm_files), 1, regexpr("_segm", basename(segm_files))-1), ".tif")

minsize_set <- seq(40, 70, 10)

for(i in seq(length(segm_files))){
  for(minsize in minsize_set){
    outfile <- paste0(substr(segm_files[[i]], 1, nchar(segm_files[[i]])-4),
                      "_mins", minsize, ".tif")
    muf_lsmss <- otbcli_LSMSSmallRegionsMerging(x = filter_files[[i]],
                                                inseg = segm_files[[i]],
                                                out = outfile,
                                                minsize = minsize,
                                                tilesizex = 1000,
                                                tilesizey = 1000,
                                                verbose = FALSE,
                                                return_raster = FALSE,
                                                ram="8192")
  }
}


# Conversion of segments to polygon shape file ---------------------------------
segm_minsize_files <- list.files(path_muf_set1m_sample_segm, 
                                 pattern = glob2rx("ortho_muf_rgb_idx_pca_scaled*segm*_mins*.tif"),
                                 full.names = TRUE)
muf <- paste0(path_muf_set1m_sample_segm, "ortho_muf_rgb_idx_pca_scaled.tif")


for(i in seq(length(segm_minsize_files))){
  print(i)
  outfile <- paste0(substr(segm_minsize_files[[i]], 1, nchar(segm_minsize_files[[i]])-4),
                    "_vector", ".shp")
  if(!file.exists(outfile)){
    muf_lsmss <- otbcli_LSMSVectorization(x = muf,
                                          inseg = segm_minsize_files[[i]],
                                          out = outfile,
                                          tilesizex = 1000,
                                          tilesizey = 1000,
                                          verbose = FALSE,
                                          ram="8192")
  }
}
