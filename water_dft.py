from pyscf import gto, dft, scf, lib
import numpy as np

mol = gto.M(atom="""
O 0.000000 0.000000 0.118
H 0.000000 0.755000 -0.471
H 0.000000 -0.755000 -0.471
""", unit='Angstrom', basis='cc-pvtz', verbose=4)

mf = dft.RKS(mol)
mf.xc = 'B3LYP'
mf.diis_space = 8
mf.level_shift = 0.5
mf.conv_tol = 1e-8
mf.max_cycle = 200
mf.kernel()

# density and overlap checks
dm = mf.make_rdm1()
S = mol.intor('int1e_ovlp')
print("trace(P S) = ", np.trace(dm @ S))

# dipole
dip_ao = mol.intor('int1e_r', comp=3)
mu_elec_au = -np.array([np.einsum('ij,ij->', dip_ao[ax], dm) for ax in range(3)])

coords = mol.atom_coords()   # in Bohr
charges = np.array([mol.atom_charge(i) for i in range(mol.natm)])
# Fixed: proper numpy array multiplication
nuc_dip_au = np.sum(charges[:, np.newaxis] * coords, axis=0)

mu_total_au = nuc_dip_au + mu_elec_au
AU_TO_DEBYE = 2.541746286
print("Total dipole (Debye):", mu_total_au * AU_TO_DEBYE, "Magnitude:", np.linalg.norm(mu_total_au * AU_TO_DEBYE))
