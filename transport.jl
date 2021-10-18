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
    Qs   = Q00  # Qs stands for superficie
    Qb   = Q00  # Qb stands for bulk
    TR   = T
    TRD  = TD
    TRs  = T
    TRDs = TD

    QDs   = Qs  # Qs stands for superficie

    n = 10
    for i in 1:n
        Z    = Qb  × TR
        ZD   = Qb  × TRD
        ZDs  = QDs × TRDs
        Zs   = Qs  × TRs

        if i < n # avoid computing at the last loop
            TRDs = TRDs × ZD
            TRs  = TRs  × Z
            TR   = TR   × Z
            TRD  = TRD  × ZD
        end

        QDs  = ZDs × Z  × QDs  # "L" stands for LEFT
        Qs   = Zs  × ZD × Qs  # "R" stands for RIGHT

        Qb   = ( (Z × ZD) + (ZD × Z) ) × Qb
    end
    #
    return QDs, Qs # LEFT, RIGHT

    # for _ in 1:10
    #     if is_L_or_R == "L" #"L" stands for LEFT
    #         Z    = Qb × TR
    #         ZD   = Qb × TRD
    #         ZDs  = Qs × TRDs
    #         # Zs   = Qs × TRs

    #         TRDs = TRDs × ZD
    #         # TRs  = TRs × Z
    #         TR   = TR   × Z
    #         TRD  = TRD  × ZD

    #         Qs   = ZDs  × Z  × Qs
    #         # Qs   = Zs × ZD × Qs
    #         Qb   = ( (Z × ZD) + (ZD × Z) ) × Qb

    #     elseif is_L_or_R == "R" # "R" stands for RIGHT
    #         Z    = Qb × TR
    #         ZD   = Qb × TRD
    #         # ZDs  = Qs × TRDs
    #         Zs   = Qs × TRs

    #         # TRDs = TRDs × ZD
    #         TRs  = TRs × Z
    #         TR   = TR  × Z
    #         TRD  = TRD × ZD

    #         # Qs   = ZDs  × Z  × Qs
    #         Qs   = Zs × ZD × Qs
    #         Qb   = ( (Z × ZD) + (ZD × Z) ) × Qb
    #     end
    # end
    # #
    # return Qs
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
