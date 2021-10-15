



U = 1.0
eCC = 2.7
eCB = 2.0
eCN = 2.0
eNB = 2.0
atom_types = 3
pair_energies = [] #

atoms = getZGNR()

[n_up, n_dw] = getConvergedSpin(atoms, U)

plotBands(n_up, n_dw)

plotSpin(n_up, n_dw)

