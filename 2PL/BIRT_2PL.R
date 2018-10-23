#############################################
# Author: Jinwen Luo
# Date: Oct 21, 2018
# Language: R
# Bayesian 2PL Model
# Reference : VanHoudnos(2012)
#############################################
# Set Working Directory
setwd("Desktop/learnC/B_IRT/")

#############################################
# Response Generation
source(file = "response_generate.r")
##*****************************************##

#############################################
# MCMC Algorithm Shell
source(file = "mcmc_alg.r")

# M-H Sampler for Parameters
source(file = "mh_sample.r")

# Generic Sampler for Normal Proposal Distribution
source(file = "generic_normal_prop.r")

## Generate proposal function
prop.th.abl <- generic.normal.proposal('th')
prop.a.disc <- generic.normal.proposal('a')
prop.b.diff <- generic.normal.proposal('b')

## Complete Conditional (CC) Distribution
source(file = "logprob.r")
### CC for the person ability parameters 
th.abl.cc <- function( U.data, state ) {
        return( ## Calculate the log-likelihood term
                apply( log.prob(U.data, state$th, state$a, state$b),1,sum)
                ## And add the prior log-density term 
                + dnorm(state$th, 0, sqrt(state$s2), log=TRUE))
}
### CC for the item discrimination parameters
a.disc.cc <- function( U.data, state ) {
        return( ## Calculate the log-likelihood term
                apply(log.prob(U.data, state$th, state$a, state$b),2,sum) 
                ## And add the prior log-density term
                + dlnorm(state$a, state$hyperpar$mu.a, sqrt(state$hyperpar$s2.a), log=TRUE))
}
### CCfor the item discrimination parameters
b.diff.cc <- function( U.data, state ) {
        return( ## Calculate the log-likelihood term
                apply(log.prob(U.data, state$th, state$a, state$b),2,sum)
                ## And add the prior log-density term
                + dnorm(state$b, 0, sqrt(state$hyperpar$s2.b), log=TRUE))
}

## Define the sampler
## Define the person ability sampler 
sample.th <- function( U.data, old ) {
        mh.sample( U.data, old,
                   cc.log.density    = th.abl.cc, 
                   proposal.function = prop.th.abl
        )
}


## Define the item discrimination sampler
sample.a <- function( U.data, old ) {
        mh.sample( U.data, old,
                   cc.log.density    = a.disc.cc, 
                   proposal.function = prop.a.disc
        )
}

## Define the item difficulty sampler
sample.b <- function( U.data, old ) {
        mh.sample( U.data, old,
                   cc.log.density    = b.diff.cc, 
                   proposal.function = prop.b.diff
        )
}
## Define the variance of person ability with Gibbs sampling
sample.s2 <- function(U.data, old) {
        ## Grab the appropriate values
        alpha.th  <- old$hyperpar$alpha.th
        beta.th   <- old$hyperpar$beta.th
        P.persons <- length(old$th)
        
        ## Update the chain 
        cur    <- old
        cur$s2 <- 1/rgamma(1, shape = alpha.th + P.persons/2,
                           rate  = beta.th + sum((old$th)^2)/2 )
        return(cur)
}
##*****************************************##

#############################################
# RUN The Chain

set.seed(314159)

# Specify the hyperparameters
hyperpars <- list( mu.a = 1.185,
                   s2.a = 1.185,
                   s2.b = 100,
                   alpha.th = 1,
                   beta.th = 1)

run.A <- run.chain.2pl(
        ## No burnin, keep 1000 iterations, do not thin
        M.burnin = 0, M.keep  = 1000, M.thin = 1,
        U.data  = U,
        hyperpars = hyperpars,
        ## Generate starting values for parameters
        th.init = runif( length(theta.abl), min=-5, max=5),
        a.init  = runif( length(a.disc), min=1/3, max=3),
        b.init  = runif( length(b.diff), min=-5, max=5),
        s2.init = runif(1, min = 0, max = 5),
#       MH.th=.85, MH.a=.15, MH.b=.15, verbose=TRUE)
# specify the MH tuning parameters
        MH.th=.80, MH.a=.12, MH.b=.12, verbose=TRUE)

#*****************************************##


#############################################
# Visualizing the chain
source(file = "get_2pl_params.r")
require(coda)
## Convert to mcmc type
run.A.mcmc <- coda::mcmc(t(run.A))
##  Trace Plot
plot( run.A.mcmc[ ,get.2pl.params(1,1,1,1)], density=FALSE, smooth=TRUE )
## burnin = 200
run.A.mcmc.conv <- window( run.A.mcmc, start=200)
plot( run.A.mcmc.conv[, get.2pl.params(1,1,1,1)], density=FALSE, smooth=TRUE )

## 9 person ability parameters
plot( run.A.mcmc.conv[, get.2pl.params(1:9,NULL,NULL,NULL)], 
      density=FALSE, smooth=TRUE, ylim=c(-5,5) )

## 9 item discrimination parameters
plot( run.A.mcmc.conv[, get.2pl.params(NULL,1:9,NULL,NULL)], 
      density=FALSE, smooth=TRUE, ylim=c(-2,2) )

## 9 item difficulty parameters
plot( run.A.mcmc.conv[, get.2pl.params(NULL,NULL,1:9,NULL)], 
      density=FALSE, smooth=TRUE, ylim=c(-4,-1) )

## sigma 
plot( run.A.mcmc.conv[, get.2pl.params(NULL,NULL,NULL,1)], 
      density=FALSE, smooth=TRUE, ylim=c(0,2) )

##############################################################
## Recovering parameters
## Visual Checkout
source(file = "check_eap.r")
all.eap <- apply( run.A.mcmc.conv, MARGIN=2, mean )
## visually compare
## Person Ability
check.sampler.graph(theta.abl, all.eap[ get.2pl.params(1:P.persons,NULL,NULL,NULL)],
        desc="Person Ability", ylab="EAP Estimates", col="blue" )
## Item discrimination 
check.sampler.graph(
        a.disc,
        all.eap[ get.2pl.params(NULL,1:I.items,NULL,NULL)],
        desc="Item discrimination", ylab="EAP Estimates", col="blue" )
## Item difficulty
check.sampler.graph(
        b.diff,
        all.eap[ get.2pl.params(NULL,NULL,1:I.items,NULL)],
        desc="Item difficulty", ylab="EAP Estimates", col="blue" )
## Sigma
check.sigma( run.A.mcmc.conv, c(1.2,2))






