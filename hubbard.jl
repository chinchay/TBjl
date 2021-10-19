
function self_consistence()

end

# G_up = 0.0d0
# G_dw = 0.0d0
function get_init_Green(G_up, G_dw, U, ϵ, V, δ)
    for s in 1:nAtoms
        hubbard_up = U(s) * (n_dw - 0.5)
        hubbard_dw = U(s) * (n_up - 0.5)
        ϵ_new_up = ϵ(s) + V(s) + hubbard_up
        ϵ_new_dw = ϵ(s) + V(s) + hubbard_dw
        G_up[s, s] = 1.0 / (E - ϵ_new_up + δ*im)
        G_dw[s, s] = 1.0 / (E - ϵ_new_dw + δ*im)
    end
end
        
