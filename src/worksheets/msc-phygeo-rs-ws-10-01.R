# rs-ws-08-2
# MOC - Data Analysis (T. Nauss, C. Reudenbach)
# Predict LCC using gpm

# Set environment --------------------------------------------------------------
if(Sys.info()["sysname"] == "Windows"){
  source("D:/active/moc/msc-ui/scripts/msc-phygeo-ei/src/functions/set_environment.R")
} else {
  source("/media/permanent/active/moc/msc-ui/scripts/msc-phygeo-ei/src/functions/set_environment.R")
}


# Read land-cover prediction data and segments 
muf_lcc_prediction <- raster(paste0(path_muf_set1m_sample_rsws093, 
                                    "muf_lcc_prediction.tif"))
projection(muf_lcc_prediction) <- CRS("+init=epsg:25832")

segments <- readOGR(paste0(path_muf_set1m_sample_segm, 
                           "ortho_muf_rgb_idx_pca_scaled_segments.shp"),
                    "ortho_muf_rgb_idx_pca_scaled_segments")
projection(segments) <- CRS("+init=epsg:25832")

re <- setValues(muf_lcc_prediction, seq(ncell(muf_lcc_prediction)))
re_fields_nincl <- extract(re, segments)

# saveRDS(re_fields_nincl, paste0(path_muf_set1m_sample_rsws093, "re_fields_nincl.rds"))
# 
# re_fields_nincl <- readRDS(paste0(path_muf_set1m_sample_rsws093, "re_fields_nincl.rds"))


for(sgmt in re_fields_nincl){
  muf_lcc_prediction[sgmt] <- modal(muf_lcc_prediction[sgmt])
}

muf_lcc_prediction
writeRaster(muf_lcc_prediction, paste0(path_muf_set1m_sample_rsws093, 
                                       "muf_lcc_prediction_semgt_mod.tif"))
