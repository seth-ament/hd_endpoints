# construct a multivariate classifier and predict its power as an endpoint in HD therapeutics studies

# 1. assemble the endpoint data

require( hdendpoints )
require( glmnet )
require( impute )

data = endpoint_data

# find a subset of mice and phenotypes with mostly complete data
# impute the missing values for remaining features
drop.measure = which( colSums(is.na(data$measure)) > 70 )
drop.mouse = which( rowSums(is.na(data$measure)) > 3 )
x0 = as.matrix( data$measure[ -drop.mouse , -drop.measure ] )

# the resulting predictors and response variables are 
measure = impute.knn( x0 , k = 3 )$data
pheno = data$row_metadata[ -drop.mouse , ]
age = data$row_metadata$age[ -drop.mouse ]


# 2. now, construct multivariate classifiers by penalized logistic regression (elastic net)


target.ages = c(3,9,12)
classification = list()
beta.coefficients = matrix( NA , ncol = 3 , nrow = ncol(x) + 1 )
rownames(beta.coefficients) = c( "Intercept" , colnames(x) )
colnames(beta.coefficients) = c("3mo","9mo","12mo")
for( a in 1:3 ) { # begin classification loop

target.age = target.ages[a]

y = pheno[ age == target.age , "genotype" ]
x = measure[ age == target.age , ]

# training set performance
fit.full = cv.glmnet( x = x , y = y , family = "binomial" , alpha = 0.5 , nfolds = 5 ) 
pred.response.train = predict( fit.full , newx = x , type = "response" , s = fit.full$lambda.1se )[,1]
beta.coefficients[,a] = as.matrix(coef(fit.full))[,1]

# evaluate test set performance (leave-one-out cross-validation)
pred.class = rep( NA , length(y) )
pred.link = rep( NA , length(y) )
pred.response = rep( NA , length(y) )
for( i in 1:length(y) ) {
  xtrain = x[-i,]
  ytrain = y[-i]
  xtest = matrix( x[i,] , nrow = 1 )
  ytest = y[i]
  fit = cv.glmnet( 
	x = xtrain , 
	y = ytrain , 
	family = "binomial" , 
	alpha = 0.5 , 
	nfolds = 5 )
  pred.class[i] = predict( fit , newx = xtest , s = fit$lambda.1se , type = "class" )
  pred.link[i] = predict( fit , newx = xtest , s = fit$lambda.1se , type = "link" )
  pred.response[i] = predict( fit , newx = xtest  , s = fit$lambda.1se , type = "response" )
}

res = data.frame( pheno[ age == target.age , ] , pred.class , pred.link , pred.response , pred.response.train )
classification[[a]] = res[ order( res$pred.response ) , ]

} # end classification loop
# merge results across ages
classification = rbind( classification[[1]] , classification[[2]] , classification[[3]] )

save( classification, beta.coefficients , file="glmnet.classifiers.RData" )

# 3. use classifiers for power analysis


# create a function to simulate random correlated variables with known mean and standard deviation
simulate_x = function( x , sigma = cor(measure) , n , ntot = 1000 ) {
sim.x = mvrnorm( 
        n = ntot ,
        mu = colMeans(x) ,
        Sigma = sigma )
sd.simx = apply( sim.x , 2 , sd )
sd.x = apply( x , 2 , sd )
mu.simx = colMeans( sim.x )
resid.simx = t( t(sim.x) - mu.simx )
sim.x.scaled = t( ( mu.simx + t(resid.simx) / (sd.simx/sd.x) ) )
simx.final = sim.x.scaled[ sample(1:ntot,n) , ]
return( simx.final )
}

# create a function to calculate power for a given smaple size and age, across a range of effect sizes
simulate_power = 
function( measure , pheno , target.age , n = 10 , sims = 100 , alpha = 0.05 , verbose = F ) {

# select training data and construct classifier
possible.effects = 0:100 / 100 #seq( from=0 , to=delta , by = delta/100 )
y = pheno[ age == target.age , "genotype" ]
x = measure[ age == target.age , ]
fit = cv.glmnet( y = y , x = x , family = "binomial" , alpha = 0.5 , nfold = length(y) )

# start the simulation
# which variables have a floor (0)
floor = rep( FALSE , ncol(x) )
floor[ c(13:20) ] = TRUE

# summarize real data
wt.mean = colMeans( x[ y == "+/+", ] )
het.mean = colMeans( x[ y == "Q111/+", ] )
#wt.sd = apply( x[ y == "+/+", ] , 2 , sd )
#het.sd = apply( x[ y == "Q111/+", ] , 2 , sd )
delta=het.mean - wt.mean

sim.pvals = matrix( 1 , ncol=sims , nrow = length(possible.effects) )
for( sim in 1:sims ) {
 for (i in 1:length(possible.effects)){

  # Increment the mean of the simulated treated HD mice by 1/100 of 
  # the true delta of their means to simulate increased rescue
  increment <- possible.effects[i]

  treatment = c(
          rep("drug" , n*2 ) ,
          rep("control",n*2 ))
  geno.sim = rep( c(rep("+/+",n),rep("Q111/+",n)) , 2 )

# original version did not retain correlation structure among variables
#  x.sim = matrix( ncol = ncol(x) , nrow = 4*n )
#  for( j in 1:ncol(x) ) {
#        x.sim[,j] = c(
#          rnorm( n=n, mean=wt.mean[j] , sd=wt.sd[j] ),
#          -delta[j]*increment+rnorm(n=n , mean=het.mean[j] , sd=het.sd[j] ),
#          rnorm( n=n , mean = wt.mean[j] , sd = wt.sd[j] ),
#          rnorm( n=n , mean = het.mean[j] , sd = het.sd[j] ) )
#        if( floor[j] == T ) {
#    	  negative.values = which( x.sim[,j] < 0 )
#    	  if( length(negative.values) > 0 ) {
#            x.sim[ negative.values , j ] = 0
#	  }
#        }
#  }
 
  x.sim = rbind( simulate_x( x = x[ y == "+/+", ] , n = n ) ,
		t( -delta*increment + t(simulate_x( x = x[ y == "Q111/+", ] , n = n )) ) ,
                simulate_x( x = x[ y == "+/+", ] , n = n ) ,
                simulate_x( x = x[ y == "Q111/+", ] , n = n ) )
	for( j in 1:ncol(x) ) {
         if( floor[j] == T ) {
          negative.values = which( x.sim[,j] < 0 )
          if( length(negative.values) > 0 ) {
            x.sim[ negative.values , j ] = 0
          }
         }
	}

  pred.sim = predict( fit , newx = x.sim , type="response" , s = fit$lambda.1se ) 
  fit.scores = lm( pred.sim ~ geno.sim * treatment )
  sim.pvals[i,sim] = drop1( fit.scores , ~. , test="F" )[ 4 , 'Pr(>F)' ]
 
 }
 if( verbose == T ) cat("done",sim,"of",sims,"simulations\n")
}

# calculate power based on simulations
powers = rowSums( sim.pvals < alpha ) / sims
# output the results
param = list( n = n , sims = sims , alpha = alpha , age = target.age )
outp = list( simulated.pvals = sim.pvals , power = powers , param = param )
return( outp )

}


# 4. run power analysis for ages = 3,9,12 and n = 10,20,30

simulation_results2 = list()
for( i in 1:3 ) {
  res.n = list()
  for( j in 1:3 ) {
     cat("Working on: age =" , c(3,9,12)[j] , ", n =" , c(10,20,30)[i] , "\n" )
     res.age = simulate_power( 
	measure = measure , 
	pheno = pheno , 
	sims=100 , 
	target.age = c(3,9,12)[j] , 
	n = c(10,20,30)[i] )
     res.n[[j]] = res.age
  }
  names(res.n) = c("3m","9m","12m")
  simulation_results2[[i]] = res.n
}
names(simulation_results2) = c("n10","n20","n30")

save( simulation_results2 , file="simulation_results.elastic_net.RData")

# plot the power to detect an interaction effect for each gene at each time point (qPCR)

pdf("power.multivariate_endpoint.pdf")
for( n in c("n10","n20","n30") ) {
  par( bty = "l" )
  # plot 12-month data
  y = simulation_results2[[n]]$'12m'$power
  plot( x = c(0,length(y)) ,
        y = c(0,1) ,
        type = "n" ,
        xlab = "" , ylab = "" )
  points( y , type = "p" )
  lines( smooth.spline(  x = 1:length(y) , y = y , spar = 0.2 , nknots = 15 ))

  # add the 9-month data
  y = simulation_results2[[n]]$'9m'$power
  points( y , type = "p" , col = "blue" , pch = 22 )
  lines( smooth.spline(  x = 1:length(y) , y = y , spar = 0.2 , nknots = 15 ) , col = "blue" )

  # add the 3-month data
  y = simulation_results2[[n]]$'3m'$power
  points( y , type = "p" , col = "red" , pch=25 )
  lines( smooth.spline(  x = 1:length(y) , y = y , spar = 0.2 , nknots = 15 ) , col = "red" )

  abline( h = 0.8 , lty = 2 )
  mtext( side = 1 , line = 2.5 , "Percent Rescue" )
  mtext( side = 2 , line = 2.5 , "Power (alpha = 0.05)" )
  mtext( side = 3 , line = 0.5 , adj = 0 , font = 2 , 
	paste( "elastic net," , n ) )

  legend( x = 1 , y = 1 ,
        legend = c("12m","9m","3m") ,
        col = c("black","blue","red") ,
        pch = c( 1 , 22 , 25 ) ,
        lty = 1 ,
        bty = "n" )
}
dev.off()



