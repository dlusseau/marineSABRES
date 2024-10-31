##############################################################
#### R code for transforming WP4 graphs to df     ############
#### bmjv                                         ############
##############################################################
####   11 July 2024

install.packages("XML")
library(XML)

doc <- xmlParse("C:/Users/bmjv/OneDrive - Danmarks Tekniske Universitet/Postdoc/Marine SABRES/Models WP4/CLD_Arctic.xml")

xmlToDataFrame(nodes = getNodeSet(doc, "//attr"))[c("attrlabl", "attrdef", "attrtype", "attrdomv")]

xml_data <- xmlToList(doc)

install.packages("pysd2r")
library(pysd2r)

Arctic <- system.file(paste0("C:/Users/bmjv/OneDrive - Danmarks Tekniske Universitet/Postdoc/Marine SABRES/Models WP4/", "CLD_Arctic.xmile"), package = "pysd2r")
py <- pysd_connect()
Arctic <- read_vensim(py, Arctic)

Arctic <- read_xmile(py, Arctic)

get_doc(Arctic)

system.file()
