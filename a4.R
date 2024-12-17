invisible('
- uses Ne(t) and the musseco package 
- ?augmented likelihood, sample probability a / v 
- coalescentjl for tree simulation 
- optim model fits 
')
library(ape)
library( ggplot2 ) 
library( phydynR ) 
library( glue )
library( foreach ) 
library( doParallel )
library( mlesky ) 

# library( musseco )
source( '/home/erik/git/musseco/R/musseco.R' )


trids = 1:20 


read.csv( './musseco1.csv' , stringsAs=FALSE ) -> trdf 

mutselp <- function( s,mu ){
	mu / ( 1 - (1+s)*(1+mu)+mu )
}

trid <- 11 
# trid <- 18 

# proc_trid <- function(trid)
{

	st0 <- Sys.time() 
	tr = read.tree( glue( "trees1/{trid}.nwk" ) )
	tr$tip.label<- gsub( 'V', 'variant', tr$tip.label )
	tr$tip.label<- gsub( 'A', 'ancestral', tr$tip.label )
	demes = sapply( strsplit( tr$tip.label, split='\\.'), '[', 2 )

	sts <- node.depth.edgelength( tr )[1:Ntip(tr)] |> setNames( tr$tip.label )
	sts <- sts + (5e3 - max(sts))
	maxsts <- max(sts) 
	isv <-  grepl( names(sts) , patt = '.V' ) 
	vsts <- sts[isv]
	asts <- sts[ !isv  ]
	ssts <- matrix( 0, nrow = Ntip(tr), ncol = 2 )
	colnames(ssts) <- c( 'ancestral', 'variant' )
	rownames(ssts) <- tr$tip.label 
	ssts[ demes == 'ancestral' , 'ancestral'] <- 1.0 
	ssts[ demes == 'variant' , 'variant'] <- 1.0 

	f = mlskygrid( tr, sampleTimes = sts , tau = NULL, tau_lower = 0.1, tau_upper = 1e6, ncpu = 4 ) 
	Net <- cbind( f$time, f$ne )
	# plot(f)

	bdt <- DatedTree( tr, sts, sampleStates = ssts) 

	samplik <- function(s, mu)
	{# 
		pv = mutselp( s, mu ) 
		dbinom( isv, size = 1, prob = pv , log=TRUE ) |> sum()
	}

	lfun <- function( theta  )
	{
		muom = exp(theta) 
		s = muom[2]-1 
		mu = muom[1] 
		sl <- samplik( s,mu )

		l = loglikelihood_bisseco(  exp(theta) , bdt, Net  )

		print( c( l, sl, l + sl )); print( Sys.time() )
		round( l+sl, 2 ) # rounding to hasten convergence 

	}

	theta0 <- log( c( .04 , 0.75) )
	

	# st2 <- Sys.time() 
	 # l = lfun(theta0) 
	# st3 <- Sys.time() 
	# st3 - st2 
	 
	st1 = Sys.time() 

	o = optim( par = theta0, fn = lfun, method = 'Nelder-Mead', control = list(trace=-6, fnscale=-1))
	saveRDS( o, glue( 'a4-o-{trid}.rds') )

	print( trdf[trid, ] )

	muom <- exp( o$par ) 
	print( muom )
	print( muom[2]-1)

	o
}

