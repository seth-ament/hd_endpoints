
require( hdendpoints )
require( doBy )
require( impute )
require( gplots )

data( endpoint_data )
data = endpoint_data


meta = endpoint_data$row_metadata

# find a subset of mice and phenotypes with mostly complete data
# impute the missing values for remaining features
drop.measure = which( colSums(is.na(data$measure)) > 70 )
drop.mouse = which( rowSums(is.na(data$measure)) > 3 )
x0 = as.matrix( data$measure[ -drop.mouse , -drop.measure ] )

# the resulting predictors and response variables are 
measure = impute.knn( x0 , k = 3 )$data
pheno = data$row_metadata[ -drop.mouse , ]
age = data$row_metadata$age[ -drop.mouse ]

measure.mean = colMeans( measure )
measure.sd = apply( measure , 2 , sd )
norm = t( ( t(measure) - measure.mean ) / measure.sd )

measure = norm

rownames( meta ) = meta$mouse
meta = meta[ rownames(measure) , ]

powers = matrix( NA , ncol=3 , nrow=20 )
effect.sizes = matrix( NA , ncol=3 , nrow=20 )
for( a in 1:3 ) {
  for( m in 1:20 ) {
     use = meta$age == c(3,9,12)[a]
     genotype = meta[ use , "genotype" ]
     y = measure[ use , m ]
     x = data.frame( y , genotype )
     x.sum<-summaryBy(y~genotype,x,FUN=c(mean,sd,length))
     delta<-x.sum[1,2]-x.sum[2,2]
     possible.effects = delta / 2
     sim = simulate_single_endpoint( 
        genotype = genotype ,
        measure = y , 
        n = 10 ,
        possible.effects = possible.effects )
     powers[ m , a ] = sim$powers[1]
     effect.sizes[ m , a ] = delta
  }
}


vars = c("Drd1a mRNA (striatum)","Drd2 mRNA (striatum)",
   "Darpp2 mRNA (striatum)","Scn4b mRNA (striatum)",
   "Cnr1 mRNA (striatum)","Homer1 mRNA (striatum)",
   "Abhd1 mRNA (liver)", "Insig2 mRNA (liver)",
   "Islr2 mRNA (liver)","H60b mRNA (liver)",
   "N4bp2 mRNA (liver)","Cell Ratio (striatum)","Cell Total (striatum)",
   "Fraction NeuN+ (striatum)","Darpp2 protein (striatum)",
   "mHTT Aggregates (striatum, MW8)" , "mHTT Area (striatum, MW8)",
   "mHTT Aggregates (striatum, p62)" , "mHTT Area (striatum, p62)",
   "synaptophysin (striatum)")


rownames(powers) = vars
colnames(powers) = c("3mo","9mo","12mo")

order = order( rowMeans( powers ) , decreasing = T )
powers = powers[ order , ]
tlab = 1:20*4-1.5
cex = 1.2

pdf("single_endpoints_power_50percentrescue.pdf" , width = 12 , height = 6 )
par( las = 3 , mar = c(15,8,2,2) )
barplot2( 
   height = t(powers) ,
   beside = T , axisnames=F , cex.axis = cex ,
   col = c("grey","darkgrey","black") )
mtext( side = 2 , line = 3 , "Power" , cex = cex )
abline( h = 0 )
axis( 1 , at=tlab , labels = F )
text( x = tlab , y = par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
labels = rownames(powers) , srt = 45 , adj = 1 , xpd = T , cex = cex )
legend( x = 65 , y = 0.9 ,
   fill = c("grey","darkgrey","black" ) ,
   bty = "n" , cex = cex ,
   legend = c("3 months","9 months","12 months") )
dev.off()




