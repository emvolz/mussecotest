invisible('
- uses Ne(t) and the musseco package 
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
library( mlesky ) 

library( musseco )
# source( '/home/erik/git/musseco/R/musseco.R' )


trids = 1:20 


read.csv( './musseco1.csv' , stringsAs=FALSE ) -> trdf 

trid <- 11 

tr = read.tree( glue( "trees1/{trid}.nwk" ) )
isvariant <- grepl( tr$tip.label, patt = '.V' )
fb = fitbisseco( tr, isvariant, Net = NULL, theta0 = log(c(.05, .95, 1)), optim_parms = list(), mlesky_parms = list(tau = 1e5, tau_lower = 1e3, ncpu = 4) )

fb
#           coef(x)
# mu     0.03527616
# omega  0.77193676
# s     -0.22806324
# [1] "Likelihood: "

trdf[ trid, ]
#             s         mu     omega
# 11 -0.2421053 0.03768433 0.7578947


