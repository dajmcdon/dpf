library(XML)
url = 'http://mazurka.org.uk/auto/earis/mazurka68-3/'
document = htmlParse(url)
all_links = as.vector(xpathSApply(document, "//a/@href"))
links = unique(all_links[substr(all_links, nchar(all_links)-8, nchar(all_links)) == "svaed.txt"])
for(i in 1:length(links)){
    download.file(paste(url, links[[i]], sep = ""), destfile = paste('extras/raw_data/', links[[i]], sep = ""))
}
