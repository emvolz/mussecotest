library( phydynR ) 
library( ape ) 
library( glue )

read.csv( './musseco1.csv' , stringsAs=FALSE ) -> trdf 

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


simtree <- function(trid, s, mu )
{
	tr = read.tree( glue( "trees1/{trid}.nwk" ) )
	demes = sapply( strsplit( tr$tip.label, split='\\.'), '[', 2 )

	sts <- node.depth.edgelength( tr )[1:Ntip(tr)] |> setNames( tr$tip.label )
	sts <- sts + (5e3 - max(sts))
	ssts <- matrix( 0, nrow = Ntip(tr), ncol = 2 )
	colnames(ssts) <- dnames 
	rownames(ssts) <- tr$tip.label 
	ssts[ demes == 'A' , 'A'] <- 1.0 
	ssts[ demes == 'V' , 'V'] <- 1.0 

	p = trueparms 
	p$s = s 
	p$mu = mu 
	# a2:  resimulate tree here 
	tr = sim.co.tree( p , M
			, x0 = c( A = 1., V = 1 )
			, t0 = 0 
			, sampleTimes = sts 
			, sampleStates = ssts 
			, res = 50
			, finiteSizeCorrections=FALSE 
	)
	tr
}

# test
# tr = simtree( 11, trdf$s[11], trdf$mu[11] )
