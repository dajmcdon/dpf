library(XML)

tmp = tempdir()
# put stuff here, fix the python script to write there

reticulate::source_python('getMeta.ph') # this "should" give R access to all the
# python objects, you may not need to save the csv

url = 'http://mazurka.org.uk/auto/earis/mazurka68-3/'
document = htmlParse(url)
all_links = as.vector(xpathSApply(document, "//a/@href"))
links = unique(all_links[substr(all_links, nchar(all_links)-8, nchar(all_links)) == "svaed.txt"])
for(i in 1:length(links)){
    download.file(paste(url, links[[i]], sep = ""), destfile = paste0(tmp, '/', links[[i]]))
}


source('read-data.R')

unlink(tmp) # deletes the tmp directory and its contents
