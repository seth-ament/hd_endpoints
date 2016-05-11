make_endpoint_data_table = function(x=1) {

require(doBy)

# pre-process the qPCR data
data(qpcr)
targets = splitBy(~endpoint,qpcr)
nTargets = length(targets)
mice = unique( qpcr[,1] )
qpcr.mat = mice
for( i in 1:nTargets ) {
  tmp = targets[[i]][,c(1,6)]
  colnames(tmp)[2] = names(targets)[i]
  qpcr.mat = merge( qpcr.mat , tmp , by=1 , all.x = T  )
}
rownames(qpcr.mat) = qpcr.mat[,1]
qpcr.mat = qpcr.mat[,-1]


# pre-process the IHC data
# take the median value for each mouse
data(ihc)
ihc$genotype<-factor(ihc$genotype,levels=c("+/+","Q111/+"))
ihc.targets<-splitBy(~endpoint,ihc)


n = length(ihc.targets) - 1
mice = setdiff( unique( ihc[,2] ) , c("",NA) )
ihc.mat = mice
for( i in 1:n ) {
  tmp = ihc.targets[[i]][,c(2:5,6)]
  tmp2 = summaryBy( measure ~ mouse , tmp , FUN=median )
  colnames(tmp2)[2] = names(ihc.targets)[i]
  ihc.mat = merge( ihc.mat , tmp2 , by=1 , all.x = T )
}

# combine the qPCR and IHC data matrices
dat = merge( qpcr.mat , ihc.mat , by.x = 0 , by.y = 1 , all = T )
rownames(dat) = dat[,1]
endpoint_data = dat[,-1]

pdf( "qq.pdf" )
par( mfrow = c(5,5) , mar = c(2,2,1,1) , bty = "l" )
for( i in 1:ncol(endpoint_data) ) {
  qqnorm( endpoint_data[,i] , main = "" )
  mtext( side=3 , line=-1 , adj = 0.05 , cex = 1 , colnames(endpoint_data)[i] )
  if( i %in% c(1:14,20) ) {
    par( new = T )
    qqnorm( log( endpoint_data[,i] ) , col="red" , yaxt="n" , main = "" )
    axis( side = 4 )
    mtext( side=3 , line=-2.5 , adj=0.05 , "log" , col="red" )
  }
  if( i %in% c(18,21,23) ) {
    par( new = T )
    qqnorm( asin( endpoint_data[,i]/100 ) , col="red" , yaxt="n" , main="" )
    axis( side = 4 )
    mtext( side=3 , line=-2.5 , adj=0.05 , "asin" , col="red" )
  }
  if( i %in% c(19,22) ) {
    par( new = T )
    qqnorm( asin( endpoint_data[,i] ) , col="red" , yaxt="n" , main = "" )
    axis( side = 4 )
    mtext( side=3 , line=-2.5 , adj=0.05 , "asin" , col="red" )
  }
}
dev.off()

endpoint_data[ , c(1:14,20) ] = log( endpoint_data[ , c(1:14,20) ] )
endpoint_data[ , c(18,21,23) ] = asin( endpoint_data[ , c(18,21,23) ] / 100 )
endpoint_data[ , c(1:14,20) ] = asin( endpoint_data[ , c(19,22) ] )


# metadata for mice
tmp = unique( ihc[,2:5] )
matchRows = match( rownames(dat) , tmp[,1] )
row_metadata = tmp[ matchRows , ]

endpoint_data = list( measure = endpoint_data , row_metadata = row_metadata )

return( endpoint_data ) 

}











