# rs-ws-09-2
# MOC - Data Analysis (T. Nauss, C. Reudenbach)
# Extract training area information

# Set environment --------------------------------------------------------------
if(Sys.info()["sysname"] == "Windows"){
  source("D:/active/moc/msc-ui/scripts/msc-phygeo-ei/src/functions/set_environment.R")
} else {
  source("/media/permanent/active/moc/msc-ui/scripts/msc-phygeo-ei/src/functions/set_environment.R")
}


# Extract training area information --------------------------------------------
# Load/combine training areas
shp_names <- list.files(path_muf_set1m_lcc_ta, 
                        pattern = glob2rx("*.shp"), full.names = TRUE)
shift <- 0
shps <- list()
for(s in seq(length(shp_names))){
  act_shps <- readOGR(shp_names[s], ogrListLayers(shp_names[s]))
  shps[[s]] <- spChFIDs(act_shps, as.character(seq(nrow(act_shps)) + shift))
  shps[[s]] <- shps[[s]][, c(which("LCC_ID" == colnames(shps[[s]]@data)),
                             which("LCC_NAME" == colnames(shps[[s]]@data)))]
  shift <- shift + nrow(act_shps)
}
lc_ta <- do.call("rbind", shps)
projection(lc_ta) <- CRS("+init=epsg:25832")
writeOGR(lc_ta, paste0(path_muf_set1m_lcc_ta, "muf_lc_ta_segm.shp"), 
         "muf_lc_ta_segm", driver = "ESRI Shapefile")


# Load raster datasets used as independent variables
muf_files <- list.files(path_muf_set1m_sample_non_segm, full.names = TRUE,
                        pattern = glob2rx("*.tif"))
muf <- stack(muf_files)

# Extract pixel IDs of training areas
re <- setValues(muf[[1]], seq(ncell(muf[[1]])))
re_lc_ta <- extract(re, lc_ta)
saveRDS(re_lc_ta, file = paste0(path_muf_set1m_lcc_ta, "re_lc_ta.rds"))

re_lc_ta_df <- lapply(seq(length(re_lc_ta)), function(i){
  data.frame(LCC_ID = lc_ta@data$LCC_ID[i],
             LCC_NAME = lc_ta@data$LCC_NAME[i],
             PIXEL = re_lc_ta[[i]])
})
re_lc_ta_df <- do.call("rbind", re_lc_ta_df)
saveRDS(re_lc_ta_df, file = paste0(path_muf_set1m_lcc_ta, "re_lc_ta_df.rds"))

# Extract raster information
names_muf <- names(muf)
muf_lc_ta_segm_traindata <- readRDS(file = paste0(path_muf_set1m_lcc_ta, "re_lc_ta_df.rds"))
for(i in seq(nlayers(muf))){
  muf_lc_ta_segm_traindata <- cbind(muf_lc_ta_segm_traindata, 
                                 muf[[i]][muf_lc_ta_segm_traindata$PIXEL])
  colnames(muf_lc_ta_segm_traindata)[i+3] <- names_muf[i]
}
saveRDS(muf_lc_ta_segm_traindata, file = paste0(path_muf_set1m_lcc_ta, 
                                             "muf_lc_ta_segm_traindata.rds"))

# muf_lc_ta_segm_traindata <- readRDS(file = paste0(path_muf_set1m_lcc_ta, 
#                                                   "muf_lc_ta_segm_traindata.rds"))


# Extract sample from class "fields x" to reduce the data set ------------------
# Compute the number of observations (i.e. pixels) per land-cover group
lcc_numbers <- as.data.frame(table(muf_lc_ta_segm_traindata$LCC_ID))
colnames(lcc_numbers) <- c("LCC", "Freq")
lcc_numbers$LCC <- as.numeric(as.character(lcc_numbers$LCC))
lcc_numbers$Classes <- "Forest"
lcc_numbers$Classes[lcc_numbers$LCC >= 101] <- "Field"
lcc_numbers$Classes[lcc_numbers$LCC >= 201] <- "Infrastructure"
lcc_numbers$Classes[lcc_numbers$LCC >= 301] <- "Other"

lcc_groupsize <- aggregate(lcc_numbers$Freq, 
                           by = list(lcc_numbers$Classes), FUN = sum)

# Reduce observations from group "Field" to finally have the same number of 
# field observations than of forest observations
negative_sample_size <- lcc_groupsize$x[lcc_groupsize$Group.1 == "Field"] - 
  lcc_groupsize$x[lcc_groupsize$Group.1 == "Forest"]

rows_field <- which(muf_lc_ta_segm_traindata$LCC_ID %in% seq(101,200))
set.seed(2017)
negative_sample_rows <- sample(rows_field, negative_sample_size, replace = FALSE)

muf_lc_ta_segm_traindata_subset <- muf_lc_ta_segm_traindata[-negative_sample_rows,]

lcc_numbers_subset <- as.data.frame(table(muf_lc_ta_segm_traindata_subset$LCC_ID))
colnames(lcc_numbers_subset) <- c("LCC", "Freq")
lcc_numbers_subset$LCC <- as.numeric(as.character(lcc_numbers_subset$LCC))
lcc_numbers_subset$Classes <- "Forest"
lcc_numbers_subset$Classes[lcc_numbers_subset$LCC >= 101] <- "Field"
lcc_numbers_subset$Classes[lcc_numbers_subset$LCC >= 201] <- "Infrastructure"
lcc_numbers_subset$Classes[lcc_numbers_subset$LCC >= 301] <- "Other"

lcc_groupsize <- aggregate(lcc_numbers_subset$Freq, 
                           by = list(lcc_numbers_subset$Classes), FUN = sum)
lcc_groupsize

# Free memory
remove(muf_lc_ta_segm_traindata)
gc()

# Write output to file
muf_lc_ta_segm_traindata_subset <- readRDS(file = paste0(path_muf_set1m_lcc_ta, 
                                                         "muf_lc_ta_segm_traindata.rds"))
