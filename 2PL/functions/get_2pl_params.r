## get.2pl.params takes ranges of parameters and appends the 
## appropriate variable names. 
get.2pl.params <- function( theta.range, a.range, 
                            b.range, sigma.range) {
        ## Append the variable names to the ranges
        if( !is.null(theta.range) ) { 
                theta.range <- paste('theta.abl', theta.range) 
        }
        if( !is.null(a.range    ) ) { 
                a.range     <- paste('a.disc'   , a.range    ) 
        }
        if( !is.null(b.range    ) ) { 
                b.range     <- paste('b.diff'   , b.range    ) 
        }
        if( !is.null(sigma.range ) ) { 
                sigma.range  <-       'sigma.theta'              
        }
        
        return( c(theta.range, a.range, b.range, sigma.range ))
}