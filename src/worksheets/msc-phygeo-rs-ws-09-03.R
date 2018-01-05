# rs-ws-08-2
# MOC - Data Analysis (T. Nauss, C. Reudenbach)
# Predict LCC using gpm

# Set environment --------------------------------------------------------------
if(Sys.info()["sysname"] == "Windows"){
  source("D:/active/moc/msc-ui/scripts/msc-phygeo-ei/src/functions/set_environment.R")
} else {
  source("/media/permanent/active/moc/msc-ui/scripts/msc-phygeo-ei/src/functions/set_environment.R")
}


# Prepare gpm data set used for remote sensing prediction study ----------------
# obsv <- readRDS(file = paste0(path_muf_set1m_lcc_ta, "muf_lc_ta_poly_large_tai_cmb_cc.rds"))
obsv <- readRDS(file = paste0(path_muf_set1m_lcc_ta, "muf_lc_ta_segm_traindata_clean.rds"))
# obsv <- obsv[which(obsv$ID %in% names(table(obsv$ID)[table(obsv$ID) > 1000])),]
p_val <- 0.05

obsv$LCC_ID <- as.factor(obsv$LCC_ID)

col_selector <- which(names(obsv) == "LCC_ID")
col_lc <- which(names(obsv) == "LCC_ID")
col_meta <- seq(which(names(obsv) == "LCC_NAME"), 
                which(names(obsv) == "PIXEL"))
col_precitors <- seq(which(names(obsv) == "band_1_ortho_muf_rgb_idx_pca_scaled_haralick__advanced_1515_11.1"),
                   which(names(obsv) == "ortho_muf_VVI"))
    
meta <- createGPMMeta(obsv, type = "input",
                      selector = col_selector, 
                      response = col_lc, 
                      predictor = col_precitors, 
                      meta = col_meta)

obsv <- gpm(obsv, meta, scale = FALSE)


# Clean predictor variables ----------------------------------------------------
if(compute){
  obsv <- cleanPredictors(x = obsv, nzv = TRUE, 
                          highcor = TRUE, cutoff = 0.80)
  saveRDS(obsv, file = paste0(path_muf_set1m_lcc, "gpm_muf_lc_ta_segm_clean_pred.rds"))
} else {
  obsv <- readRDS(file = paste0(path_muf_set1m_lcc, "gpm_muf_lc_ta_segm_clean_pred.rds"))
}


# Compile model training and evaluation dataset --------------------------------
obsv <- resamplingsByVariable(x = obsv,
                              use_selector = FALSE,
                              resample = 1)
    
# Split resamples into training and testing samples
obsv <- splitMultResp(x = obsv, 
                      p = p_val, 
                      use_selector = FALSE)


# Train and build model --------------------------------------------------------
cl <- makeCluster(detectCores())
registerDoParallel(cl)
obsv <- trainModel(x = obsv,
                   n_var = NULL, 
                   mthd = "rf",
                   mode = "rfe",
                   seed_nbr = 11, 
                   cv_nbr = 5,
                   var_selection = "indv", 
                   filepath_tmp = NULL)
saveRDS(obsv, file = paste0(path_muf_set1m_sample_rses093, "gpm_muf_lc_ta_segm_rf_model.rds"))
# obsv <- readRDS(file = paste0(path_muf_set1m_sample_rses093, "gpm_muf_lc_ta_segm_rf_model.rds"))

obsv@model$rf_rfe[[1]][[1]]$model
tstat <- compContTests(obsv@model$rf_rfe, mean = FALSE)

vi <- compVarImp(obsv@model$rf_rfe)
plotVarImp(vi)

# Predict lcc ------------------------------------------------------------------
library(randomForest)

muf_files <- list.files(path_muf_set1m_sample_non_segm, full.names = TRUE,
                        pattern = glob2rx("*.tif"))
muf_files
muf <- stack(muf_files)
muf_df <- getValues(muf)


gpm_muf_lc_ta_segm_rf_predict <- predict(obsv@model$rf_rfe[[1]][[1]]$model$fit, newdata = muf_df,
                              na.action = na.pass)

muf_lcc_prediction <- setValues(muf[[1]], as.numeric(as.character(gpm_muf_lc_ta_segm_rf_predict)))


writeRaster(muf_lcc_prediction, paste0(path_muf_set1m_sample_rses093, "muf_lcc_prediction.tif"), overwrite = TRUE)



