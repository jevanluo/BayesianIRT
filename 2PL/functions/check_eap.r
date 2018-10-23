## We'll be doing these checks a lot, so let's make a function
check.sampler.graph <- function(true.values, estimated.values,
                                desc, ylab, col) {    
        par(mfrow=c(1,2))
        plot( true.values, estimated.values,
              xlab='True values', ylab=ylab,
              main=paste('Parameter recovery:',desc),
              col=col)
        abline(0,1)
        
        plot( ecdf( true.values ), xlab=deparse(substitute(true.values)),
              main=paste('Empirical CDFs:',desc) )
        plot( ecdf( estimated.values ), col=col, lty='dashed', add=TRUE )
        legend( 'topleft', c('True values', ylab), col=c('black', col),
                lty=c('solid','dashed'), bg='white')
        ## Add in the KS test
        legend( 
                'bottomright',
                paste("Kolmogorov-Smirnov p value:",
                      round( ks.test( true.values, estimated.values )$p.value, 2)),
                bg='white')
        par(mfrow=c(1,1))
}

check.sigma <- function( mcmc.conv , xylim) {
        par(mfrow=c(1,2))
        
        traceplot(mcmc.conv[, 'sigma.theta'],
                  smooth=TRUE, ylim=xylim,
                  main='Trace plot sig2.theta')
        abline(h=sig2.theta , col='blue')
        
        densplot( mcmc.conv[, 'sigma.theta'],
                  main='Density plot sig2.theta',
                  xlab='sig2.theta', xlim=xylim)
        abline(v=sig2.theta , col='blue')
        abline(v=var(theta.abl), col='purple')
        
        legend( "topright",
                c('true value','var(theta.abl)'),
                col=c('blue', 'purple'),
                lty='solid')
}
