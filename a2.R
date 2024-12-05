invisible('
a2: using phydynR for both tree simulation and likelihood 
	sample times based on Cjl trees 
')
library( ape )
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
		source('a2libs.R')

		s = ps[ i, 1] 
		mu = ps[i, 2 ]
		bdt <- simtree( trid, s, mu )

		lfun <- function( s, mu )
		{
			p <- list(
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

	saveRDS( pldf, glue( 'a2-pldf-{trid}.rds') )
	pldf
}

pldfs <- lapply( trids, proc_trid )

saveRDS( pldfs, file = glue('a2-pldfs.rds'))
