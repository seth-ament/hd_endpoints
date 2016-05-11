simulate_single_endpoint =
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


