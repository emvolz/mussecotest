
library(ape)
library( ggplot2 ) 
library( phydynR ) 
library( glue )
library( foreach ) 
library( doParallel )
library( mlesky ) 
library( diversitree )
library( musseco )



#' Theoretical frequency of the wild type assuming mutation & selection balance 
#' with a variant that has lower transmission fitness
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

# > pancestral_mutsel_balance1( .001, 52, 2000, .85 )
# [1] 0.7531224


# More likely to see accumulation of transmission-deleterious mutations in HIV than e.g. flu b/c of longer generation times (dtpars[3])
sapply( 1:9, function(x) 
{
	# zero is variant, one is ancestral 
	# λ0 , λ1 , µ0 , µ1 , q01 , q10 
	dtpars <- c( (1-x*.10)*2/10
		, 2/10
		, 1/10, 1/10  
		, 0.0015
		, 0.0015*(1*(1+x))
	) 
	print( '-------------' )
	print( x )
	print( dtpars )
	c( dtR0( dtpars ) , pancestral_mutsel_balance1( dtpars[5], dtpars[3] , dtpars[6]/dtpars[5], dtpars[1]/dtpars[2] ) ) |> print()  
	print( dtpars[1] / dtpars[2] )
	dtpars 
}
) -> parmmatrix 
rownames(parmmatrix) <- c( 'λ0' , 'λ1' , 'µ0' , 'µ1' , 'q01' , 'q10'  )

simtree <- function( dtpars )
{
	set.seed( 1111 ) 
	dtphy <- tree.bisse(dtpars, max.taxa = 1e3, x0=1)
	while( is.null( dtphy )) 
		dtphy <- tree.bisse(dtpars, max.taxa = 1e3, x0=1)
	dtphy 
}
#
# tr = simtree( parmmatrix[,9] )
# plot(tr)
# table( tr$tip.state )
