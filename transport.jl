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

    b = Q00     # for Qb, Qb stands for bulk
    ν₊ = T      # for TR
    ν₋ = TD     # for TRD

    s₊ = Q00    # for Qs, Qs stands for superficie, for funcion de superficie RIGHT
    s₋ = Q00    # for QDs, for funcion de superficie LEFT
    τ₊ = T      # for TRs
    τ₋ = TD     # for TRDs


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
            Zetas = ( (β₊ × β₋) + (β₋ × β₊) ) × b
            b     = renorm(Ident, Zetas, b)  # Q00 = inv(  Ident - G00*T00 )*G00
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


function greenRenorm(G,T00,T,TD, Ident)
    Zetas = G × T00
    b  = renorm(Ident, Zetas, G)     # GR ## = inv(  Ident - G*T00 )*G
    ν₊ = T      # TR
    ν₋ = TD     # TRD

    n = 10
    for i in 1:n
        # update central functions (not bulk)
        β₊ = b × ν₊  # = Z ## Z(N-1)   = GR(N-1)*TR(N-1)
        β₋ = b × ν₋  # = ZD ## ZzD(N-1) = GR(N-1)*TRD(N-1)

        # avoid computing at the last loop, before updating `b`
        if i < n
            # update interactions τ and ν
            ν₊ = ν₊ × β₊   # TR(N) = TR(N-1)*GR(N-1)*TR(N-1)
            ν₋ = ν₋ × β₋   # TRD(N) = TRD(N-1)*GR(N-1)*TRD(N-1)
        end
        
        # update bulk function
        Zetas = (β₊ × β₋) + (β₋ × β₊)
        b = renorm(Ident, Zetas, b)  #GR(N)= inv(...)*GR(N-1)

    end
    #   
    return b
end


function greenRenorm_bulk_surface(G,T00,T,TD, Ident, is_bulk_or_surface)
    Zetas = G × T00
    # for Qb, Qb stands for bulk
    b  = renorm(Ident, Zetas, G)  # GR ## = inv(  Ident - G*T00 )*G
    ν₊ = T      # for TR
    ν₋ = TD     # for TRD

    if is_bulk_or_surface == "surface"
        s₊ = b    # for Qs, Qs stands for superficie, for funcion de superficie RIGHT
        s₋ = b    # for QDs, for funcion de superficie LEFT    
        τ₊ = T      # for TRs
        τ₋ = TD     # for TRDs
    end

    n = 10
    for i in 1:n
        # update central functions (not bulk)
        β₊ = b × ν₊  # = Z  ## Z(N-1)   = GR(N-1)*TR(N-1)
        β₋ = b × ν₋  # = ZD ## ZzD(N-1) = GR(N-1)*TRD(N-1)
        
        if is_bulk_or_surface == "surface"
            # update central functions (not bulk)
            σ₋ = s₋ × τ₋  # = ZDs
            σ₊ = s₊ × τ₊  # = Zs

            # update surface functions
            s₋ = σ₋ × β₊ × s₋  # "L" stands for LEFT
            s₊ = σ₊ × β₋ × s₊  # "R" stands for RIGHT

            # avoid computing at the last loop
            if i < n
                # update interactions τ
                τ₋ = τ₋ × β₋
                τ₊ = τ₊ × β₊
            end
        end

        # avoid computing at the last loop
        if i < n
            # update interactions ν
            ν₊ = ν₊ × β₊   # TR(N) = TR(N-1)*GR(N-1)*TR(N-1)
            ν₋ = ν₋ × β₋   # TRD(N) = TRD(N-1)*GR(N-1)*TRD(N-1)
        end

        # avoid computing at the last loop if surface green functions are asked:
        if (i < n) || (is_bulk_or_surface == "bulk")
            # update bulk function
            Zetas = (β₊ × β₋) + (β₋ × β₊)
            b = renorm(Ident, Zetas, b)  #GR(N)= inv(...)*GR(N-1)
        end

        if is_bulk_or_surface == "surface"
            return s₋, s₊ # LEFT, RIGHT surface functions
        else
            return b # bulk
        end
end