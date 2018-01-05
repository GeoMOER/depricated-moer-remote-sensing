# rs-ws-09-2
# MOC - Data Analysis (T. Nauss, C. Reudenbach)
# Extract training area information

# Set environment --------------------------------------------------------------
if(Sys.info()["sysname"] == "Windows"){
  source("D:/active/moc/msc-ui/scripts/msc-phygeo-ei/src/functions/set_environment.R")
} else {
  source("/media/permanent/active/moc/msc-ui/scripts/msc-phygeo-ei/src/functions/set_environment.R")
}


# Load/combine training areas --------------------------------------------------
if(compute){
  shp_names <- list.files(path_muf_set1m_lcc_ta_rsws091, 
                          pattern = glob2rx("*.shp"), full.names = TRUE)
  shift <- 0
  shps <- list()
  for(s in seq(length(shp_names))){
    act_shps <- readOGR(shp_names[s], ogrListLayers(shp_names[s]))
    projection(act_shps) <- CRS("+init=epsg:25832")
    shps[[s]] <- spChFIDs(act_shps, as.character(seq(nrow(act_shps)) + shift))
    shps[[s]] <- shps[[s]][, c(which("LCC_ID" == colnames(shps[[s]]@data)),
                               which("LCC_NAME" == colnames(shps[[s]]@data)))]
    shift <- shift + nrow(act_shps)
  }
  names(shps) <- basename(shp_names)
  
  # Compute areas
  areas <- lapply(shps, function(s){
    sum(area(s))
  }) 
  
  # Subsample fields, add information from original polygons and create a new
  # SpatialPointsDataFrame
  fields_subset <- spsample(shps$felder.shp, n = areas$Waldgebiete.shp, 
                            type = "regular", offset = c(0.5, 0.5))
  fields_information <- extract(shps$felder.shp, fields_subset)
  muf_lc_ta_segm_fields_incl <- SpatialPointsDataFrame(fields_subset, 
                                                       fields_information[, 3:4])
  writeOGR(muf_lc_ta_segm_fields_incl, 
           paste0(path_muf_set1m_lcc_ta_rsws091, "muf_lc_ta_segm_fields_incl.shp"), 
           "muf_lc_ta_segm_fields_incl", driver = "ESRI Shapefile")
  
  # Remove original fields polygons and combine the rest
  shps$felder.shp <- NULL
  muf_lc_ta_segm_fields_nincl <- do.call("rbind", shps)
  writeOGR(muf_lc_ta_segm_fields_nincl, 
           paste0(path_muf_set1m_lcc_ta_rsws091, "muf_lc_ta_segm_fields_nincl.shp"), 
           "muf_lc_ta_segm_fields_nincl", driver = "ESRI Shapefile")
} else {
  muf_lc_ta_segm_fields_incl <- 
    readOGR(paste0(path_muf_set1m_lcc_ta_rsws091, "muf_lc_ta_segm_fields_incl.shp"), 
            "muf_lc_ta_segm_fields_incl")
  muf_lc_ta_segm_fields_nincl <- 
    readOGR(paste0(path_muf_set1m_lcc_ta_rsws091, "muf_lc_ta_segm_fields_nincl.shp"), 
            "muf_lc_ta_segm_fields_nincl")
}


# Load raster datasets used as independent variables ---------------------------
muf_files <- list.files(path_muf_set1m_sample_non_segm, full.names = TRUE,
                        pattern = glob2rx("*.tif"))
muf <- stack(muf_files)


# Extract pixel IDs of training areas without fields ---------------------------
if(compute){
  re <- setValues(muf[[1]], seq(ncell(muf[[1]])))
  re_fields_nincl <- extract(re, muf_lc_ta_segm_fields_nincl)
  
  muf_lc_ta_segm_fields_nincl_df <- 
    lapply(seq(length(muf_lc_ta_segm_fields_nincl)), function(i){
      data.frame(LCC_ID = muf_lc_ta_segm_fields_nincl@data$LCC_ID[i],
                 LCC_NAME = muf_lc_ta_segm_fields_nincl@data$LCC_NAME[i],
                 PIXEL = re_fields_nincl[[i]])
    })
  muf_lc_ta_segm_fields_nincl_df <- do.call("rbind", muf_lc_ta_segm_fields_nincl_df)
  saveRDS(muf_lc_ta_segm_fields_nincl_df, 
          file = paste0(path_muf_set1m_lcc_ta_rsws091, "muf_lc_ta_segm_fields_nincl_df.rds"))
} else {
  muf_lc_ta_segm_fields_nincl_df <- readRDS(
    file = paste0(path_muf_set1m_lcc_ta_rsws091, "muf_lc_ta_segm_fields_nincl_df.rds"))
}


# Extract pixel IDs of training areas for fields -------------------------------
if(compute){
  re_fields_incl <- extract(re, muf_lc_ta_segm_fields_incl)
  
  muf_lc_ta_segm_fields_incl_df <- 
    data.frame(LCC_ID = muf_lc_ta_segm_fields_incl@data$LCC_ID,
               LCC_NAME = muf_lc_ta_segm_fields_incl@data$LCC_NAME,
               PIXEL = unlist(re_fields_incl))
  saveRDS(muf_lc_ta_segm_fields_incl_df, 
          file = paste0(path_muf_set1m_lcc_ta_rsws091, "muf_lc_ta_segm_fields_incl_df.rds"))
} else {
  muf_lc_ta_segm_fields_incl_df <- saveRDS(
    file = paste0(path_muf_set1m_lcc_ta_rsws091, "muf_lc_ta_segm_fields_incl_df.rds"))
}


# Combine extracted pixel information and extract raster values ----------------
if(compute){
  muf_lc_ta_segm_df <- rbind(muf_lc_ta_segm_fields_nincl_df,
                             muf_lc_ta_segm_fields_incl_df)
  
  # Extract raster information
  names_muf <- names(muf)
  muf_lc_ta_segm_traindata <- muf_lc_ta_segm_df
  for(i in seq(nlayers(muf))){
    muf_lc_ta_segm_traindata <- cbind(muf_lc_ta_segm_traindata, 
                                      muf[[i]][muf_lc_ta_segm_traindata$PIXEL])
    colnames(muf_lc_ta_segm_traindata)[i+3] <- names_muf[i]
  }
  saveRDS(muf_lc_ta_segm_traindata, file = paste0(path_muf_set1m_lcc_ta_rsws091, 
                                                  "muf_lc_ta_segm_traindata.rds"))
} else {
  muf_lc_ta_segm_traindata <-readRDS(
    file = paste0(path_muf_set1m_lcc_ta_rsws091, "muf_lc_ta_segm_traindata.rds"))
}

# Remove predictors with NA values ---------------------------------------------
na_sums <- sapply(muf_lc_ta_segm_traindata, function(x) sum(is.na(x)))
muf_lc_ta_segm_traindata_clean <- muf_lc_ta_segm_traindata[, names(na_sums[na_sums == 0])]
saveRDS(muf_lc_ta_segm_traindata_clean, file = paste0(path_muf_set1m_lcc_ta_rsws091, 
                                                "muf_lc_ta_segm_traindata_clean.rds"))


