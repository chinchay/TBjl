module Integrations

# how to use it:
# include("Integrations.jl"); using .Integrations
# X = Integrations.get_X_for_GaussIntegration()

# you can use import, but I found more useful when combined with `Revise` (read )

# reload a module easily:
# https://timholy.github.io/Revise.jl/stable/
# using Revise
# include("Integrations.jl"); using .Integrations
# X = Integrations.get_X_for_GaussIntegration()
# ld = Integrations.ld

# in case you don't use Revise:
# include("Integrations.jl"); import .Integrations: get_X_for_GaussIntegration, Gauss5Quad
# instead of just `using .Integration`
# export get_listX_for_GaussIntegration, Gauss5Quad

# for get_listX_for_GaussIntegration()
const a  = 2.0 * sqrt( 10.0 / 7.0  )
const dx11 = -sqrt(5.0 + a) / 3.0
const dx22 = -sqrt(5.0 - a) / 3.0
const dx33 =  0.0
const dx44 = -dx22
const dx55 = -dx11
const xo = 0.0
const xf = 1.0
const nIntervals = 4 # enough or 5
const dx = (xf - xo) / nIntervals
const dx1 = dx * (0.5 + (dx11 / 2.0))
const dx2 = dx * (0.5 + (dx22 / 2.0))
const dx3 = dx * (0.5 + (dx33 / 2.0))
const dx4 = dx * (0.5 + (dx44 / 2.0))
const dx5 = dx * (0.5 + (dx55 / 2.0))
const ld = [dx1, dx2, dx3, dx4, dx5]

# for Gauss5Quad()
const c0 = 13.0 * sqrt(70.0)
const c1 = (322.0 - c0) / 900.0
const c2 = (322.0 + c0) / 900.0
const c3 = 128.0 / 225.0
const c4 = c2
const c5 = c1
const lc = [c1  c2  c3  c4  c5] # lc must be [a1 a2 a3 a4 a5] : 1×5 Matrix{Float32}
const fac = (xf - xo) / nIntervals / 2.0


function get_X_for_GaussIntegration()
    X = zeros(Float32, (nIntervals, 5))

    for i in 1:nIntervals
        D = xo + ((i - 1) * dx)
        X[i,:] .= D .+ ld
    end

    return X
end


function Gauss5Quad(Y, auxSum) # `Y` must be a nIntervals x 5 matrix
    
    # TODO: test if performance changes when using `.= @~ ` instead of `.=`
    auxSum  .= lc .* Y # it multplies c1 x a column in Y, c2 x second column of Y, ...
    
    integral = sum(auxSum) * fac
    return integral
end


end # end module

