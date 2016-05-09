#####################################################################
#####################################################################
##########                                             
##########  Power analysis via simulated rescue of QRTPCR targets and
##inta########  histological measures for Q111/+ natural history paper
##########  
##########  Jeff Carroll 
##########
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
####### Simulation - Try once for Scn4b
##################################################################################################

# Pick one target (Scn4b) and run a simmulation study
# First, get summary
scn<-rt.targets[["Scn4b"]][,c(3,4,6)]
scn<-scn[scn$age=="12",c(1,3)]

# Pretend these are data from a control population
scn$treatment<-"control"
scn.sum<-summaryBy(measure~genotype,scn,FUN=c(mean,sd,length))
delta<-scn.sum[1,2]-scn.sum[2,2]

simulation.results.scn<-data.frame()
for (i in 1:100){
  
  # Increment the mean of the simulated treated HD mice by 1/100 of 
  # the true delta of their means to simulate increased rescue
  increment <- (delta/100) * i
  
  simulated.data <- data.frame(genotype=c(rep("+/+",scn.sum[1,4]),rep("Q111/+",scn.sum[2,4])),
                               measure=c(rnorm(n=scn.sum[1,4],mean=scn.sum[1,2],sd=scn.sum[1,3]),
                                             increment+rnorm(n=scn.sum[2,4],mean=scn.sum[2,2],sd=scn.sum[2,3])),
                               treatment=rep("drug",scn.sum[1,4]+scn.sum[2,4]))
    
  joined.data<-rbind(scn,simulated.data)
  
  simulated.result<-anova(lm(measure~genotype*treatment,joined.data))
  interaction.pvalue<-simulated.result$'Pr(>F)'[3]
  
  simulation.results.scn<-rbind(simulation.results.scn,data.frame(i=i,interaction.pvalue=interaction.pvalue))
  
}


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
