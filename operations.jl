#*#############################################################################
"""
    smart_prod(A, B, C)

Compute the product of three matrices efficiently by determining the best 
factor order either `(A * B) * C` or `A * (B * C)`
"""
function smart_prod(A::Any, B::Any, C::Any)
    a1, a2 = size(A)
    b1, b2 = size(B)
    c1, c2 = size(C)
    @assert a2 == b1
    @assert b2 == c1
    
    # (A * B) * C
    complexity1 = (a1 * a2 * b2) + (a1 * b2 * c2)

    if a1 == a2 == b2 == c2
        return (A * B) * C, complexity1, -1, true
    end

    # A * (B * C)
    complexity2 = (b1 * b2 * c2) + (a1 * a2 * c2)

    if complexity1 < complexity2
        return (A * B) * C, complexity1, complexity2, true
    else
        return A * (B * C), complexity1, complexity2, false
    end
    
    # return (complexity1 <= complexity2) ? (A * B) * C : A * (B * C)
    
end

#*#############################################################################
# """
#     ×(A, B, C)

# Ternary operator definition for smart_prod(A, B, C)
# # Example
# ```julia-repl
# julia> A × B × C

# ```
# """
# ×(A::Any, B::Any, C::Any) = smart_prod(A, B, C)





#*#############################################################################

using LinearAlgebra #for  LinearAlgebra.BLAS.hemm(), and for diagind(out, 0)

"""
    ×(A, B)

Binary operator to substitute LinearAlgebra.BLAS.hemm('L', 'U', A, B)
for the product of hermitic matrices.

# References
- [BLAS functions](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#BLAS-functions)
- To know more aobut usage, see line 217 in [here](https://github.com/JuliaLang/julia/blob/ae8452a9e0b973991c30f27beb2201db1b0ea0d3/stdlib/LinearAlgebra/test/blas.jl)
  - Cnn, Cnm, Cmn = Matrix{ComplexF32}.(undef,((2,2), (2,1), (1,2)))
        

# Example
```julia-repl
julia> A × B

```
"""
×(A::Any, B::Any) = BLAS.hemm('L','U', A, B) 

#*#############################################################################
# func_up = 0.0
# func_dw = 0.0
# G_up    = 0.0
# G_dw    = 0.0

function get_Hubbard_energies(U, n)
    return U .* (n .- 0.5) # equivalent to U(site) * (n[site] - 0.5) for site in 1:nAtoms
end

function get_potential_energies(ϵ, V, hubbard)
    return ϵ .+ V .+ hubbard # ϵ_new   = ϵ(s) + V(s) + hubbard for s in 1:nAtoms
    
end

function get_diagIndices(G)
    return diagind(G, 0)
end

function get_init_Green!(E, G, diagIndices, n, U, ϵ, V, δ)
    G .= 0.0 # initialize
    hubbard = get_Hubbard_energies(U, n)
    ϵ_new   = get_potential_energies(ϵ, V, hubbard)

    # https://discourse.julialang.org/t/vectorized-variable-assignment-assign-matrix-off-diagonal-values/16363
    G[ diagIndices ] .= 1.0 ./ (E .- ϵ_new .+ δ*im)
    return G
end

# get_init_Green!(G, G_dw, U, ϵ, V, Δ) # modifies G_up and G_dw
# G = 0
function integrandComplexPlane(Efermi, y, G, n, U, ϵ, V, T00, T, TD, Ident) # G can be the initialized G_up or G_dw, and n can be n_dw, or n_up
    Δ = y / (1.0 - y)
    get_init_Green!(Efermi, G, diagIndices, n, U, ϵ, V, Δ)
    greenRenorm!(G, T00, T, TD, Ident) # renomalize G_up (now it is GRup)
    
    # defining the integrand
    d = (1.0 - y) ^ 2
    f = real( G[ diagIndices ] ) ./ d
    return f
end


function get_n(Efermi, G, n, n_, U, ϵ, V, T00, T, TD, Ident)
    using QuadGK
    yo = 0.0
    yf = 1.0
    integral, err = quadgk(
        y -> integrandComplexPlane(Efermi, y, G, n, U, ϵ, V, T00, T, TD, Ident),
        yo, yf, rtol=1e-4)
    # 
    n_ .= 0.5 .+ integral / π # integrate at each column
    return n_
end

nAtoms = 3
G_up = zeros(Float32, (nAtoms, nAtoms))
G_dw = zeros(Float32, (nAtoms, nAtoms))
n_up = zeros(Float32, nAtoms)
n_dw = zeros(Float32, nAtoms)
n_up_ = zeros(Float32, nAtoms)
n_dw_ = zeros(Float32, nAtoms)
U = zeros(Float32, nAtoms)
ϵ = zeros(Float32, nAtoms)
V = zeros(Float32, nAtoms)
T00 = zeros(Float32, (nAtoms, nAtoms))
T = zeros(Float32, (nAtoms, nAtoms))
TD = zeros(Float32, (nAtoms, nAtoms))
Ident = zeros(Float32, (nAtoms, nAtoms))


n_up_ = get_n(Efermi, G_up, n_dw, U, ϵ, V, T00, T, TD, Ident) #array of length=nAtoms
n_dw_ = get_n(Efermi, G_dw, n_up, U, ϵ, V, T00, T, TD, Ident)

function dQ_exceso(atomosT,Ndizimac,nAvgUPi,nAvgDOWNi,Efermi,hzeeman,Ue,Eoarray, Uarray,Vbordarr,T00,T,TD,Ident)
    IntegraPlanoCmplx(nAvgUPi,nAvgDOWNi,atomosT,hzeeman,Eoarray,Uarray,Vbordarr,Ue,Efermi,T00,T,TD,Ident,Ndizimac)
    dQ_exceso = Sum(  nAvgUPi + nAvgDOWNi   )  - atomosT ! dNe es el ***exceso del NUMERO DE ELECTRONES*** en el ribbon todo.
end

