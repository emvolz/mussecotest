using Coalescent 
using Plots
using Phylo 
using Distributions
using StatsBase 
using Interpolations
using DataFrames
using CSV 

sampsize = 1000;
ntr = 20 ;
ss = range( -.5, -.01; length=ntr ) |> collect 
omegas = 1.0 .+ ss

function omegamutop( omega, mu )
	mu / ( 1 - omega*(1+mu)+mu )
end

function omegatomu( omega ; equilV = 0.15 )
	# p = (1 - omega*(1+mu)) / (1-omega*(1+mu) + mu ) )
	equilV*(1-omega)/ (1 + equilV*(omega-1))
end

mus = map( omegatomu, omegas )


m = ModelFGY( "musseco.yaml" )

function gensampconf(m::ModelFGY)

	o = solveodes(m)
	plot(o)

	staxis = range(4995, maximum(o.t), length=20)  # o.t[ o.t .> 16 ] 
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
	s

end


# tonewick(tr) |> parsenewick |> plot 

map(1:ntr) do i 
	s = ss[i] 
	mu = mus[i] 
	m.parameters["s"] = s 
	m.parameters["μ"] = mu 
	@show m.parameters 
	sconf = gensampconf( m )
	tr = SimTree( m, sconf ) 
	otr = tonewick(tr )
	write( "trees1/$(i).nwk" , otr ) 
end

odf = DataFrame( :s => ss, :mu => mus, :omega => omegas )
CSV.write( "musseco1.csv", odf )
