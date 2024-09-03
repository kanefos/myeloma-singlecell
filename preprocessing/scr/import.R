# Download from Zenodo

dir.create('data/import/',showWarnings=F)
options(timeout = max(1000, getOption("timeout")))
files=c('dataDescription.xlsx','metadata-donor.csv','Tcell-scTCR.csv','panImmune.h5ad','Tcell.h5ad')
for (file_name in files){
url <- paste0("https://zenodo.org/records/13646014/files/",file_name,"?download=1")
file_path <- "data/import/"
download.file(url, paste(file_path, file_name, sep = ""), mode = "wb", quiet=F)
}
