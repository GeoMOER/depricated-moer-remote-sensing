# rs-ws-03-1
# MOC - Data Analysis (T. Nauss, C. Reudenbach)
# Merge "white" aerial images

# Set path ---------------------------------------------------------------------
if(Sys.info()["sysname"] == "Windows"){
  filepath_base <- "D:/active/moc/msc-phygeo-remote-sensing-2016/"
} else {
  filepath_base <- "/media/TOSHIBA\ EXT//GFO/BushEncroachment/"
}

path_data <- paste0(filepath_base, "data/")
path_aerial <- paste0(path_data, "aerial/")
path_aerial_merged <- paste0(path_data, "aerial_merged/")
path_temp <- paste0(filepath_base, "temp/")


# Libraries --------------------------------------------------------------------
library(raster)

rasterOptions(tmpdir = path_temp)


# Merge aerial files and write resulting raster to separate file ----------------
aerial_files <- list.files(path_aerial, full.names = TRUE, 
                           pattern = glob2rx("*_1.tif"))

for(name1 in aerial_files){
  name2 <- paste0(substr(name1, 1, nchar(name1)-6), ".tif")
  fn <- overlay(stack(name1), stack(name2), fun = min)
  dir.create(path_aerial_merged, showWarnings = FALSE)
  writeRaster(fn, filename = 
                paste0(paste0(path_aerial_merged, basename(name2))))
  
  # Rename original files to mark them
  file.rename(name1, paste0(name1, ".deprc"))
  file.rename(name2, paste0(name2, ".deprc"))
}



aerial_files <- list.files(path_aerial, full.names = TRUE, 
                           pattern = glob2rx("*.tif"))
for(i in aerial_files){
  s <- stack(i)
  crs(s) <- +proj=utm +zone=32 +ellps=GRS80 +units=m +no_defs 
  
  writeRaster(s, paste0(dirname(i), "/test/", basename(i)), format="GTiff")
}

