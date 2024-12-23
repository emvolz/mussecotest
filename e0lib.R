library(ape)
library( ggplot2 ) 
library( phydynR ) 
library( glue )
library( foreach ) 
library( doParallel )
library( mlesky ) 
library( diversitree )
library( musseco )

pancestral_mutsel_balance1 <- function( mu, gamma, alpha, omega, tol = 1e-3 ){
	proot <- function( p ){# 
		p*(gamma*(1-p)+mu*(1-p)-mu*alpha*p) - (omega*(1-p))*(gamma*p + mu*alpha*p - mu*(1-p))  
	}
	uniroot( proot, c(0,1), tol = tol )$root
}

dtR0 <- function( dtpars ) 
{
	pa = pancestral_mutsel_balance1( dtpars[5], dtpars[3] , dtpars[6]/dtpars[5], dtpars[1]/dtpars[2] )
	pa * dtpars[2] / dtpars[4] + (1-pa)*dtpars[1]/dtpars[3]
}

sapply( seq( 1, 4, length=12 ), function(x) 
{
	# λ0 , λ1 , µ0 , µ1 , q01 , q10 
	dtpars <- c( 0.00+x/100
		, 0.06-x/100/2
		, 0.03, 0.03
		, 0.003+6*x/1000
		, 0.01 #-x/1000
	) 
	print( '-------------' )
	print( x )
	print( dtpars )
	c( dtR0( dtpars ) , pancestral_mutsel_balance1( dtpars[5], dtpars[3] , dtpars[6]/dtpars[5], dtpars[1]/dtpars[2] ) ) |> print()  
	print( dtpars[1] / dtpars[2] )
	dtpars 
}
) -> parmmatrix 


simtree <- function( dtpars )
{
	set.seed( 1111 ) 
	dtphy <- tree.bisse(dtpars, max.taxa = 1e3, x0=1)
	while( is.null( dtphy )) 
		dtphy <- tree.bisse(dtpars, max.taxa = 1e3, x0=1)
	dtphy 
}


