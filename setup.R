# setup ihc and qpcr data for the hd_endpoints R package

ihc = read.csv("from_jeff/ihc_endpoint_data.csv")
save(ihc , file="./hd_endpoints/data/ihc.rda")

qpcr = read.csv("from_jeff/qPCR_endpoint_data.csv")
save( qpcr , file="./hd_endpoints/data/qpcr.rda")



