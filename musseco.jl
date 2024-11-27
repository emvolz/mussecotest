using Coalescent 
using Plots
using Phylo 

sampsize = 1000
ntr= 20 

m = ModelFGY( "musseco.yaml" )
o = solveodes(m)
plot(o)

using Distributions
using StatsBase 
using Interpolations

staxis = range(95, maximum(o.t), length=20)  # o.t[ o.t .> 16 ] 
ouinterp1 = linear_interpolation( o.t , map(x->x[1],o.u) )
ouinterp2 = linear_interpolation( o.t , map(x->x[2],o.u) )
ouaxis = last( o.u, length(staxis))

Nwt =  sum( ouinterp1.(staxis) )
Ntfp = sum( ouinterp2.(staxis) )

nwt = sampsize*Nwt/(Ntfp+ Nwt ) |> floor |> Integer
ntfp = sampsize*Ntfp/(Ntfp + Nwt ) |> floor |> Integer


wtsamps  = sample( staxis, ouinterp1.(staxis)|>Weights , nwt; replace=true )
tfpsamps = sample( staxis, ouinterp2.(staxis)|>Weights , ntfp; replace=true )

s = SampleConfiguration( [ [ ("A", tt) for tt in wtsamps ]...
,[ ("V", tt) for tt in tfpsamps ]...
] )

tr = SimTree( m, s )

# tonewick(tr) |> parsenewick |> plot 

map(1:ntr) do i 
	tr = SimTree( m, s ) 
	otr = tonewick(tr )
	write( "trees/$(i).nwk" , otr ) 
end

