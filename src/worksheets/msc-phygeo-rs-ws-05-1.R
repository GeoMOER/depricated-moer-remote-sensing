# rs-ws-05-1
# MOER - Remote Sensing
# Scaling

# Set environment --------------------------------------------------------------
if(Sys.info()["sysname"] == "Windows"){
  source("C:/Users/tnauss/permanent/edu/msc-phygeo-remote-sensing/moer-msc-phygeo-remote-sensing/src/functions/set_environment.R")
} else {
  source("/media/permanent/active/moc/msc-ui/scripts/msc-phygeo-ei/src/functions/set_environment.R")
}

# Start with ortho_muf_1m.tif
# Compute spectral indices -----------------------------------------------------
muf = stack(paste0(path_muf_set1m_sample_non_segm, "ortho_muf_1m.tif"))

idx = rgbIndices(muf, rgbi = c("GLI", "NGRDI", "TGI", "VVI"))
projection(idx) = CRS("+init=epsg:25832")

filepath = paste0(path_muf_set1m_sample_non_segm, 
                  "ortho_muf_", names(idx), ".tif")
writeRaster(idx, filename = filepath, bylayer = TRUE)


# Compute PCA ------------------------------------------------------------------
files_muf_rgb <- paste0(path_muf_set1m_sample_non_segm, "ortho_muf_1m.tif")
files_muf_idx <- list.files(path_muf_set1m_sample_non_segm, 
                            pattern = glob2rx("*i.tif"),full.names = TRUE)
files_muf_rgb_idx <- c(files_muf_rgb, files_muf_idx)

pca_data <- pca(stack(files_muf_rgb_idx))
saveRDS(pca_data, file = paste0(path_rdata, "ortho_muf_rgb_idx_pca.RDS"))
projection(pca_data$map) <- CRS("+init=epsg:25832")

filepath = paste0(path_muf_set1m_sample_non_segm, 
                  "ortho_muf_rgb_idx_", names(pca_data$map), ".tif")
writeRaster(pca_data$map, filename = filepath, bylayer = TRUE)

# Scale PCA
pca_scaled = scale(pca_data$map, center=TRUE, scale=TRUE)
names(pca_scaled) = paste0(names(pca_scaled), "_scaled")
filepath = paste0(path_muf_set1m_sample_non_segm, 
                  "ortho_muf_rgb_idx_", names(pca_scaled), ".tif")
writeRaster(pca_scaled, filename = filepath, bylayer = TRUE)

# Combine all scaled PCAs into one file
filepath = paste0(path_muf_set1m_sample_non_segm, 
                  "ortho_muf_rgb_idx_pca_scaled.tif")
writeRaster(pca_scaled, filename = filepath, bylayer = FALSE)


# Compute tecture filters on scaled PC1 using glcm -----------------------------
filename = paste0(path_muf_set1m_sample_non_segm, 
                  "ortho_muf_rgb_idx_PC1_scaled.tif")
pc1 <- raster(filename)

windows <- c(3, 5, 13, 21)
windows <- c(5, 13, 21)

for(win in windows){
  gt <- glcm(pc1, n_grey = 32, window = c(win, win),
             shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)),
             statistics = c("mean", "variance", "correlation"))
  
  filepath = paste0(path_muf_set1m_sample_non_segm, names(pc1), "_", names(gt), 
                    sprintf("_w%02d", win), ".tif")
  writeRaster(gt, filename = filepath, bylayer = TRUE)
}


# Compute haralick textures on scaled PC1 --------------------------------------
filename = paste0(path_muf_set1m_sample_non_segm, 
                  "ortho_muf_rgb_idx_PC1_scaled.tif")
pc1 <- raster(filename)

minv <- minValue(pc1)
maxv <- maxValue(pc1)

windows <- c(3, 5, 13, 21)
windows <- c(5, 13, 21)

for(win in windows){
  oth <- otbTexturesHaralick(x=pc1, 
                             path_output = path_muf_set1m_sample_non_segm, 
                             return_raster = TRUE, 
                             parameters.xyrad=list(c(win, win)),
                             parameters.xyoff=list(c(1,1)),
                             parameters.minmax=c(minv, maxv),
                             parameters.nbbin = 8,
                             texture="all",
                             channel = 1)
  filepath = paste0(path_muf_set1m_sample_non_segm, names(pc1), "_", 
                    substr(names(oth), 1, regexpr(".4.20637313219754.4780660526394", names(oth))),
                    "tif")
  writeRaster(oth, filename = filepath, bylayer = TRUE)
}
