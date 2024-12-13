invisible('
- augmented likelihood, sample probability a / v 
- coalescentjl for tree simulation 
- optim model fits 
')
library(ape)
library( ggplot2 ) 
library( phydynR ) 
library( glue )
library( foreach ) 
library( doParallel )



trids = 1:20 


read.csv( './musseco1.csv' , stringsAs=FALSE ) -> trdf 

mutselp <- function( s,mu ){
	mu / ( 1 - (1+s)*(1+mu)+mu )
}

trid <- 11 
trid <- 18 

# proc_trid <- function(trid)
{

	st0 <- Sys.time() 
	
	m = 2 
	dnames = c( 'A', 'V' )

	bs <- matrix('0.0', nrow = m, ncol = m )
	rownames(bs) = colnames(bs) <- dnames 
	bs['A', 'A'] <- "parms$beta*A*(1-(A+V)/parms$K)"
	bs['V', 'V'] <-  "(1+parms$s)*parms$beta*V*(1-(A+V)/parms$K)"

	ms <- matrix('0.0', nrow = m, ncol = m )
	rownames(ms) = colnames(ms) <- dnames 
	ms[ 'A', 'V'] <- "parms$mu * A"

	ds <- rep("0.0", m ) |> setNames( dnames )
	ds["A"] <- "parms$gamma * A"
	ds["V"] <- "parms$gamma * V"

	trueparms <- list(
		beta = 1.5 
		, gamma = 1
		, K = 1e4 
		, s = -.20 
		, mu = .03 
	)

	M <- build.demographic.process( births = bs, migrations = ms, deaths = ds, parameterNames = names( trueparms ), rcpp=FALSE, sde=FALSE )

	tr = read.tree( glue( "trees1/{trid}.nwk" ) )
	demes = sapply( strsplit( tr$tip.label, split='\\.'), '[', 2 )

	sts <- node.depth.edgelength( tr )[1:Ntip(tr)] |> setNames( tr$tip.label )
	sts <- sts + (5e3 - max(sts))
	maxsts <- max(sts) 
	isv <-  grepl( names(sts) , patt = '.V' ) 
	vsts <- sts[isv]
	asts <- sts[ !isv  ]
	ssts <- matrix( 0, nrow = Ntip(tr), ncol = 2 )
	colnames(ssts) <- dnames 
	rownames(ssts) <- tr$tip.label 
	ssts[ demes == 'A' , 'A'] <- 1.0 
	ssts[ demes == 'V' , 'V'] <- 1.0 

	bdt <- DatedTree( tr, sts, sampleStates = ssts) 

	p = trueparms <- list(
		beta = 1.5 
		, gamma = 1
		, K = 1e4 
		, s = -.20 
		, mu = .03 
	)

	tfgy <- M( p, t0 = 0, t1 = maxsts, x0 = c(A = 1, V = 1 ))

	samplik <- function( bdt, M, p )
	{
		tfgy <- M( p, t0 = 0, t1 = maxsts
		, x0 = c( A = 1., V = 1 ) 
		, res = 1000 ) 
		ya <- sapply( tfgy$size, '[', 'A' )
		yv <- sapply( tfgy$size, '[', 'V' )
		at <- approxfun( rev(tfgy$time), rev(ya) ) 
		vt <- approxfun( rev(tfgy$time), rev(yv) ) 
		prv <- vt(sts) / (vt(sts) + at(sts)) 
		sum(log(prv[isv]) ) + sum( log(1-prv[!isv]))
	}

	lfun <- function( theta  )
	{
		s = theta[1] 
		logmu = theta[2] 
		mu = exp( logmu ) 


		p$s = s 
		p$mu = mu 
		l = phydynR::colik( tree = bdt, demographic.process.model=M , theta=p 
			, x0 = c( A = 1., V = 1 )
			, t0 = 0 
			, res = 100 
			, AgtY_penalty = 0 
			, likelihood='PL2'
			# , likelihood='QL'
		)
		sl <- samplik( bdt, M, p )
		print( c( l, sl, l + sl )); print( Sys.time() )
		round( l+sl, 2 ) # rounding to hasten convergence 
	}

	theta0 <- c( -.25, log( .03 ) )
	
	# st2 <- Sys.time() 
	l = lfun(theta0) 
	# st3 <- Sys.time() 
	
	st1 = Sys.time() 

	o = optim( par = theta0, fn = lfun, method = 'Nelder-Mead', control = list(trace=-6, fnscale=-1))
	saveRDS( o, glue( 'a3-o-{trid}.rds') )

	print( trdf[trid, ] )
	print( c( o$par[1], exp( o$par[2])) )

	o
}

