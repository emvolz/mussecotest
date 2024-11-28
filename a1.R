library(ape)
library( ggplot2 ) 
library( phydynR ) 
library( glue )
library( foreach ) 
library( doParallel )



trids = 1:20 

# show.demographic.process( M, x0 = c(A=1,V=1), theta = trueparms, t0 = 0, t1= 100 )


proc_trid <- function(trid)
{
	mus <- seq( 0.001, .15, length = 20 )
	ss <- seq( -.75, -.001, length = 20 )
	
	cl <- makeCluster( 24 )
	registerDoParallel( cl )
	ps <- expand.grid( ss, mus ) 
	ps$trid = trid 
	
	st0 <- Sys.time() 
	
	Z = foreach( i = 1:nrow(ps)
	    	    , .packages=c('phydynR', 'glue')
	    	    , .combine = c
	   	   ) %dopar% {
		trid = ps$trid[i] 
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
		ssts <- matrix( 0, nrow = Ntip(tr), ncol = 2 )
		colnames(ssts) <- dnames 
		rownames(ssts) <- tr$tip.label 
		ssts[ demes == 'A' , 'A'] <- 1.0 
		ssts[ demes == 'V' , 'V'] <- 1.0 

		bdt <- DatedTree( tr, sts) 

		lfun <- function( s, mu )
		{
			p = trueparms <- list(
				beta = 1.5 
				, gamma = 1
				, K = 1e4 
				, s = -.20 
				, mu = .03 
			)
			p$s = s 
			p$mu = mu 
			l = phydynR::colik( tree = bdt, demographic.process.model=M , theta=p 
				, x0 = c( A = 1., V = 1 )
				, t0 = 0 
				, res = 1000 
				, AgtY_penalty = 0 
				, likelihood='PL2'
			)
			l
		}

		lfun( ps[i,1], ps[i,2] ) 
	}
	
	st1 = Sys.time() 
	st1 - st0 
	stopCluster(cl) 

	pldf <- as.data.frame( ps )
	colnames(pldf) <- c( 's', 'mu' )
	pldf$loglik = Z 

	saveRDS( pldf, glue( 'a1-pldf-{trid}.rds') )
	pldf
}

pldfs <- lapply( trids, proc_trid )

saveRDS( pldfs, file = glue('a1-pldfs.rds'))
