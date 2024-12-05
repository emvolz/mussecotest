library( ggplot2 )

read.csv( './musseco1.csv' , stringsAs=FALSE ) -> trdf 

pldfs <- readRDS( './a2-pldfs.rds' ) 

pldf2aci <- function( d )
{
	imax <- which.max( d$loglik )
	smax <- d$s[ imax ]
	mumax <- d$mu[ imax ]
	d1 <- d[ d$loglik > (d$loglik[imax]-1.96), ]
	srange <- range( d1$s ) 
	murange <- range( d1$mu ) 
	list( c( smax, srange ) , c( mumax, murange ))
}

pldf2aci( pldfs[[1]] ) 

acis <- lapply( pldfs, pldf2aci )

sdf <- sapply( acis, function( aci ) aci[[1]] ) |> t() |> as.data.frame() 
colnames( sdf ) <- c( 'shat' , 'slb', 'sub' )
sdf$s <-  trdf$s 


mudf <- sapply( acis, function( aci ) aci[[2]] ) |> t() |> as.data.frame() 
colnames( mudf ) <- c( 'muhat' , 'mulb', 'muub' )
mudf$mu <-  trdf$mu 


pmu = ggplot( mudf, aes( x = 1:nrow(mudf), y = muhat, ymin = mulb, ymax=muub )) + geom_linerange()  + geom_point()  + geom_point( aes(x = 1:nrow(mudf), y = mu, colour = 'red') )
pmu 
ggsave( pmu, file='a2pl-pmu.png' )


ps = ggplot( sdf, aes( x = 1:nrow(sdf), y = shat, ymin = slb, ymax=sub )) + geom_linerange()  + geom_point()  + geom_point( aes(x = 1:nrow(sdf), y = s, colour = 'red') )
ps
ggsave( ps, file = 'a2pl-ps.png' )
