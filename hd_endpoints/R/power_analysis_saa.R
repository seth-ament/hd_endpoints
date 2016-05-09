#####################################################################
#####################################################################
##########                                             
##########  Power analysis via simulated rescue of QRTPCR targets and
##inta########  histological measures for Q111/+ natural history paper
##########  
# Version 0.01: Jeff Carroll 
# Version 0.02: Seth Ament, created simulation functions for normals
#####################################################################
#####################################################################

require(googlesheets)
require(dplyr)
require(doBy)
require(ggplot2)
require(compute.es)
require(gridExtra)
require(doBy)

##################################################################################################
####### Code - Read in Data ######################################################################
##################################################################################################

# setwd("~/Google Drive/Carroll Lab Cloud Storage/Projects/Q111 Mouse Characterization/qRT-PCR/Data/qPCR data")
rt<-read.csv("master endpoint qPCR data.csv")
rt.targets<-splitBy(~endpoint,rt)

# setwd("~/Google Drive/Carroll Lab Cloud Storage/Projects/Q111 Mouse Characterization/Immunohistochemistry/Master Endpoint Data")
ihc<-read.csv("_Compiled.csv")
ihc$genotype<-factor(ihc$genotype,levels=c("+/+","Q111/+"))
ihc.targets<-splitBy(~endpoint,ihc)

##################################################################################################
####### Simulation
##################################################################################################

# simulate_power() is a function to calculate power of interaction effects for normally distributed variables
simulate_power = 
function( genotype , measure , alpha = 0.05 , sims = 500 , 
	n = NULL , possible.effects = NULL , floor = FALSE ) {

x = data.frame( genotype = genotype , y = measure )

if( length(unique(x$genotype)) != 2 ) {
  cat( "Error: x does not have exactly two genotypes\n")
  break 
}
if( is.numeric(x$y) == F ) {
  cat( "Error: y is not numeric\n" )
  break 
}

# get summary statistics for the distribution of x
x.sum<-summaryBy(y~genotype,x,FUN=c(mean,sd,length))
delta<-x.sum[1,2]-x.sum[2,2]
if( is.null(possible.effects) ) {
  possible.effects = seq( from=0 , to=delta , by = delta/100 )
}
if( is.null(n) ) {
  n = x.sum[1,4]
}

sim.pvals = matrix( 1 , ncol=sims , nrow = length(possible.effects) )
for( sim in 1:sims ) {
 for (i in 1:length(possible.effects)){
  
  # Increment the mean of the simulated treated HD mice by 1/100 of 
  # the true delta of their means to simulate increased rescue
  increment <- possible.effects[i]
  
  simulated.data <- data.frame(
	genotype = c(rep("+/+",n),rep("Q111/+",n)),
        y = c(
	  rnorm( n=n,mean=x.sum[1,2],sd=x.sum[1,3]),
          increment+rnorm(n=n,mean=x.sum[2,2],sd=x.sum[2,3]),
          rnorm( n=n,mean=x.sum[1,2],sd=x.sum[1,3]),
          rnorm( n=n,mean=x.sum[2,2],sd=x.sum[2,3])),
        treatment = c( 
	  rep("drug" , n*2 ) ,
    	  rep("control",n*2 ))
  )
  if( floor == T ) {
    negative.values = which( simulated.data$y < 0 )
    if( length(negative.values) > 0 ) {
      simulated.data$y[ negative.values ] = 0
    }
  }

  fit<-lm(y~genotype*treatment,simulated.data)
  aov = anova(fit)

  sim.pvals[i,sim] <- aov$'Pr(>F)'[3]
 }
 cat("done",sim,"of",sims,"simulations\n")
}

powers = rowSums( sim.pvals < alpha ) / sims

outp = data.frame( effect.size = possible.effects , powers )
return(outp)
}

##############
# calculate power for each gene at each time point (qPCR)

rt$endpoint = paste( rt$endpoint , rt$tissue , sep="_" )

endpoints = unique( rt$endpoint )

simulation.results = list()

for( i in 1:length(endpoints) ) {

  e = endpoints[i]

  cat( "**** Working on" , e , "****\n" )

  tmp = list()
 
  x1 = rt[ rt$endpoint == e & rt$age == 12 , ]
  tmp$m12 = simulate_power( 
	genotype = x1$genotype , 
	measure = log2( x1$measure ) ,
	sims = 100 , n = 10 )

  x2 = rt[ rt$endpoint == e & rt$age == 9 , ]
  tmp$m9 = simulate_power(
        genotype = x2$genotype ,   
        measure = log2( x2$measure ) ,
        sims = 100 , n = 10 )

  x3 = rt[ rt$endpoint == e & rt$age == 3 , ]
  tmp$m3 = simulate_power(
        genotype = x3$genotype ,   
        measure = log2( x3$measure ) ,
        sims = 100 , n = 10 )

  simulation.results[[endpoints[i]]] = tmp

}
names( simulation.results ) = endpoints

save( simulation.results , file="qpcr.power.n=10.effect=percentrescue.alpha=genoxtrt<0.05.RData" )

# plot the power to detect an interaction effect for each gene at each time point (qPCR)

pdf("power.qpcr_endpoints.pdf")
for( i in 1:length(simulation.results) ) {
  par( bty = "l" )
  # plot 12-month data
  y = simulation.results[[i]][[1]]$powers
  plot( x = c(0,length(y)) , 
	y = c(0,1) , 
	type = "n" ,
	xlab = "" , ylab = "" )
  points( y , type = "p" )
  lines( smooth.spline(  x = 1:length(y) , y = y ))

  # add the 9-month data
  y = simulation.results[[i]][[2]]$powers
  points( y , type = "p" , col = "blue" , pch = 22 )  
  lines( smooth.spline(  x = 1:length(y) , y = y ) , col = "blue" )

  # add the 3-month data
  y = simulation.results[[i]][[3]]$powers
  points( y , type = "p" , col = "red" , pch=25 )
  lines( smooth.spline(  x = 1:length(y) , y = y ) , col = "red" )

  abline( h = 0.8 , lty = 2 )
  mtext( side = 1 , line = 2.5 , "Percent Rescue" )
  mtext( side = 2 , line = 2.5 , "Power (alpha = 0.05)" )
  mtext( side = 3 , line = 0.5 , adj = 0 , font = 4 , endpoints[i] )

  legend( x = 1 , y = 1 ,
	legend = c("12m","9m","3m") ,
	col = c("black","blue","red") ,
	pch = c( 1 , 22 , 25 ) ,
	lty = 1 ,
	bty = "n" )
}
dev.off()






### Jeff's code from version 0.01




ggplot(simulation.results.scn,aes(x=i,y=interaction.pvalue))+
  geom_hline(yintercept=0.05,color='red',lwd=2)+
  ggtitle("Simulation results with Scn4b mRNA rescue")+
  xlab("Rescue, percent")+
  geom_point()

##################################################################################################
####### Simulation - MW8 at 12 months of age
##################################################################################################

# Pick one target (MW8 Aggregation) and run a simmulation study
# First, get summary at 12 months of age
agg<-ihc.targets[['mw8.aggregate.size']][,c(3,5,6,7)]
agg<-agg[agg$age=="12",c(1,3)]

# Pretend these are data from a control population
agg$treatment<-"control"
agg[,-1]->agg
agg.sum<-summaryBy(measure~treatment,data=agg,FUN=c(mean,sd,length))

simulation.results.mw8<-data.frame()
i=for (i in 1:100){
  
  # Increment the mean of the simulated treated HD mice by 1/100 of 
  # the true delta of their means to simulate increased rescue
  increment <- (agg.sum[1,2]/100)*i
  
  # Make a simulated data frame with aggregates in a "drug treated group"
  simulated.data <- data.frame(measure=c(rnorm(n=agg.sum[1,4],mean=agg.sum[1,2],sd=agg.sum[1,3])),
                               treatment=rep("drug",agg.sum[1,4]))
  
  # Decrement the aggregate size in the drug treatment group by the i % of the total aggregate size
  simulated.data$measure<-simulated.data$measure-increment
  joined.data<-rbind(agg,simulated.data)
  
  # Remove values below 0, which don't make sense for this endpoint
  joined.data[joined.data$measure<0,]$measure<-0
  
  # Run anova, in this case no genotype effect because wild-types have no aggregates
  simulated.result<-anova(lm(measure~treatment,joined.data))
  treatment.pvalue<-simulated.result$'Pr(>F)'[1]
  
  simulation.results.mw8<-rbind(simulation.results.mw8,data.frame(i=i,treatment.pvalue=treatment.pvalue))
  
}

ggplot(simulation.results.mw8,aes(x=i,y=treatment.pvalue))+
  geom_hline(yintercept=0.05,color='red',lwd=2)+
  ggtitle("Simulation results with MW8 Aggregate Size (12 Months)")+
  xlab("Rescue, percent")+
  geom_point()

##################################################################################################
####### Classifier Experiment, transcript data
##################################################################################################

# To do - bootstrap curves for each endpoint
# To do - MW8 size
# To do - confuse the classifier (RF)
