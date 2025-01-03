invisible('
Application of bisseco model to bisse birth/death trees simulated with diversitree 

')

library( ape )
library( ggplot2 ) 
library( phydynR ) 
library( glue )
library( foreach ) 
library( doParallel )
library( mlesky ) 
library( diversitree )
library( musseco )


source( 'e1lib.R' )

	# zero is variant, one is ancestral 
	# λ0 , λ1 , µ0 , µ1 , q01 , q10 

fits <- foreach( i = 1:ncol(parmmatrix), .combine = c ) %dopar% {
	source( 'e0lib.R' )
	dtpars <- parmmatrix[,i]
	dtphy <- simtree(dtpars)
	dtstates <- dtphy$tip.state
	dtisv <- dtstates==0
	dtfb = fitbisseco( dtphy, dtisv, Tg=1/dtpars[3], mu = dtpars[5], Net = NULL, theta0 = log(c(2,.75,1/2)), optim_parms = list(), mlesky_parms = list(tau = NULL, tau_lower = .1, tau_upper = 1e7, ncpu = 4, model = 1 ) )
	dtfb$theoralpha <- dtpars[6] / dtpars[5]
	dtfb$theoromega <- dtpars[1] / dtpars[2]
	dtfb$theorpa <-  pancestral_mutsel_balance1( dtpars[5],  dtpars[3] , coef(dtfb)['alpha'],  coef(dtfb)['omega'] )
	dtfb$empiricalpa <- 1-mean( dtisv )
	saveRDS( dtfb, file = glue::glue( 'e1-dtfb-{i}.rds' ) )
	dtfb 
}

saveRDS( fits, file = glue::glue( 'e1-fits.rds' ))

fitsalf <- foreach( i = 1:ncol(parmmatrix), .combine = c ) %dopar% {
	source( 'e0lib.R' )
	dtpars <- parmmatrix[,i]
	dtphy <- simtree(dtpars)
	dtstates <- dtphy$tip.state
	dtisv <- dtstates==0
	dtfb = fitbisseco( dtphy, dtisv, Tg=1/dtpars[3], mu = dtpars[5], Net = NULL, theta0 = log(c(2,.75,1/2)), augment_likelihood=FALSE, optim_parms = list(), mlesky_parms = list(tau = NULL, tau_lower = .1, tau_upper = 1e7, ncpu = 4, model = 1 ) )
	dtfb$theoralpha <- dtpars[6] / dtpars[5]
	dtfb$theoromega <- dtpars[1] / dtpars[2]
	dtfb$theorpa <-  pancestral_mutsel_balance1( dtpars[5],  dtpars[3] , coef(dtfb)['alpha'],  coef(dtfb)['omega'] )
	dtfb$empiricalpa <- 1-mean( dtisv )
	saveRDS( dtfb, file = glue::glue( 'e1-dtfbalf-{i}.rds' ) )
	dtfb 
}

saveRDS( fitsalf, file = glue::glue( 'e1-fitsalf.rds' ))

