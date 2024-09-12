require(rmarkdown)
for (i in 1:5){rmarkdown::render(paste0('notebooks/',i,'.Rmd'),params=list(args = 'myarg'))}
for (i in c(1:3,5:9)){rmarkdown::render(paste0('notebooks/S',i,'.Rmd'),params=list(args = 'myarg'))}

