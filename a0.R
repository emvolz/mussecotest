library(ape)
library( phydynR ) 
library( glue )


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


# show.demographic.process( M, x0 = c(A=1,V=1), theta = trueparms, t0 = 0, t1= 100 )


trid = 11 

tr = read.tree( glue( "trees/{trid}.nwk" ) )
demes = sapply( strsplit( tr$tip.label, split='\\.'), '[', 2 )

sts <- node.depth.edgelength( tr )[1:Ntip(tr)] |> setNames( tr$tip.label )
ssts <- matrix( 0, nrow = Ntip(tr), ncol = 2 )
colnames(ssts) <- dnames 
rownames(ssts) <- tr$tip.label 
ssts[ demes == 'A' , 'A'] <- 1.0 
ssts[ demes == 'V' , 'V'] <- 1.0 

bdt <- DatedTree( tr, sts) 

lfun <- function( s, mu )
{
	p = trueparms 
	p$s = s 
	p$mu = mu 
	l = phydynR::colik( tree = bdt, demographic.process.model=M , theta=p 
		      , x0 = c( A = 1., V = 1 )
			      , t0 = 0 
			      , res = 200 )
	l
}


ofun <- function( x ) 
{
	xx = exp(x) 
	s = xx[1] 
	mu = xx[2] 
	-lfun( s, mu )
}

# lfun( -.2, .04 ) |> system.time() 
# lfun( -.2, .2 ) 

mus <- seq( 0.001, .15, length = 20 )
ss <- seq( -.5, 0, length = 20 )

library( foreach ) 
library( doParallel )

ps <- expand.grid( ss, mus ) 

cl <- makeCluster( 24 )
registerDoParallel( cl )

# st0 = Sys.time() 
# Z = foreach( i = 1:4, .packages='phydynR', .combine = c) %dopar% lfun( ps[i,1], ps[i,2] ) 
# st1 = Sys.time() 
# st1 - st0 

st0 <- Sys.time() 
Z = foreach( i = 1:nrow(ps), .packages='phydynR', .combine = c) %dopar% lfun( ps[i,1], ps[i,2] ) 
st1 = Sys.time() 
st1 - st0 

# Zmat <- as.matrix( Z, nrow = 20 ) 

library( ggplot2 ) 
pldf <- as.data.frame( ps )
colnames(pldf) <- c( 's', 'mu' )
pldf$loglik = Z 
ggplot(  pldf , aes(s, mu, fill=loglik) ) + geom_tile() 

saveRDS( pldf, file = 'a0-pldf.rds' ) 
