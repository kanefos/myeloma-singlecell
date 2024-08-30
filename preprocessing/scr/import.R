# Download from Zenodo

dir.create('preprocessing/data/import/',showWarnings=F)
options(timeout = max(1000, getOption("timeout")))
files=c('dataDescription.xlsx','metadata-donor.csv','Tcell-scTCR.csv','panImmune.h5ad','Tcell.h5ad')
for (file_name in files){
url <- paste0("https://zenodo.org/records/13199890/files/",file_name,"?download=1")
file_path <- "preprocessing/data/import/"
download.file(url, paste(file_path, file_name, sep = ""), mode = "wb", quiet=F)
}

