using "operations.jl" # to use the × operator \times

function renorm(Ident, Zetas, Q)
    return inv(Ident - Zetas) × Q
    # TODO: do not use inv(), better to use some function in BLAS to solve a matrix equation: A * X = B
end

"""
    greenSuperf(G,T00,T,TD, Ident)

gL, gR  = greenSuperf(G, T00, T, TD, Ident)

"""
function greenSuperf(G,T00,T,TD, Ident) # Green renormalizado, de superficie
    Zetas = G × T00
    Q00   = renorm(Ident, Zetas, G)  # Q00 = inv(  Ident - G00*T00 )*G00
    # Qs0 = 0.0 used for Ndizimac not defined `if (j>1) continuaGR = continuaDizimando(Qs, Qs0, atomosT)`

    s₊ = Q00    # for Qs, Qs stands for superficie, for funcion de superficie RIGHT
    s₋ = Q00    # for QDs, for funcion de superficie LEFT
    τ₊ = T      # for TRs
    τ₋ = TD     # for TRDs

    b = Q00     # for Qb, Qb stands for bulk
    ν₊ = T      # for TR
    ν₋ = TD     # for TRD

    n = 10
    for i in 1:n
        # update central functions (not bulk)
        β₊ = b  × ν₊  # = Z
        β₋ = b  × ν₋  # = ZD
        σ₋ = s₋ × τ₋  # = ZDs
        σ₊ = s₊ × τ₊  # = Zs

        # update surface functions
        s₋ = σ₋ × β₊ × s₋  # "L" stands for LEFT
        s₊ = σ₊ × β₋ × s₊  # "R" stands for RIGHT

        # avoid computing at the last loop
        if i < n
            # update interactions τ and ν
            τ₋ = τ₋ × β₋
            τ₊ = τ₊ × β₊
            ν₊ = ν₊ × β₊
            ν₋ = ν₋ × β₋

            # update bulk function
            b  = ( (β₊ × β₋) + (β₋ × β₊) ) × b
        end
    end
    #
    return s₋, s₊ # LEFT, RIGHT surface functions
end


"""
    get_transmission(G, T00, T, TD, GRenorm)

transmission_up = get_transmission(G_up, T00, T, TD, GRenorm)

transmission_dw = get_transmission(G_dw, T00, T, TD, GRenorm)

"""
function get_transmission(G, T00, T, TD, GRenorm)
    # T, TD
    # gL means green Left (green correspondiente al sistemaaislado, es minuscula)
    # gR means green Right
    gL, gR  = greenSuperf(G, T00, T, TD, Ident)
    # TODO: define Ident as global variable perhaps. Ident is used to reduce memory in its creation

    Σ_L = TD × gL × T 
    # TODO: define a function that only multiplies only the imaginary with a real part (Sigma) to reduce computation.
    # TODO: T and TD are reals except in magnetic fields!
    # TODO: but greenSuperf() takes much more time, so better to solve there!
    
    Σ_R = T × gR × TD
    Γ_L = -2 * imag(Σ_L)
    Γ_R = -2 * imag(Σ_R)
    transmission = Γ_L × GRenorm × Γ_R × GRenorm' # ' means transpose conjugated: adjunta no es adjoint(otra definicion en Algebra2)
    return transmission

end
