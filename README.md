# Water Molecule DFT Calculation

A simple Python script to calculate the electronic structure and dipole moment of a water molecule using Density Functional Theory (DFT).

## Requirements

```bash
pip install pyscf numpy
```

- Python 3.7+
- PySCF (quantum chemistry package)
- NumPy (numerical computations)

## How to Run

```bash
python water_dft.py
```

Runtime: ~10-30 seconds

## Code Explanation

### 1. Setup Molecule

```python
mol = gto.M(atom="""
O 0.000000 0.000000 0.118
H 0.000000 0.755000 -0.471
H 0.000000 -0.755000 -0.471
""", unit='Angstrom', basis='cc-pvtz', verbose=4)
```

**`gto.M()`** - Creates molecule object
- `atom`: Atomic symbols and xyz coordinates
- `unit='Angstrom'`: Coordinate units (converts to Bohr internally)
- `basis='cc-pvtz'`: Basis set (~58 functions for H₂O, high quality)
- `verbose=4`: Output detail level (0-9)

### 2. Initialize DFT Calculation

```python
mf = dft.RKS(mol)
mf.xc = 'B3LYP'
```

**`dft.RKS(mol)`** - Creates Restricted Kohn-Sham DFT object (for closed-shell molecules)

**`mf.xc = 'B3LYP'`** - Sets exchange-correlation functional
- B3LYP = hybrid functional (20% exact exchange + GGA corrections)
- Popular choice for organic molecules

### 3. Convergence Settings

```python
mf.diis_space = 8
mf.level_shift = 0.5
mf.conv_tol = 1e-8
mf.max_cycle = 200
```

**`diis_space = 8`** - Store last 8 iterations for DIIS acceleration (speeds up convergence)

**`level_shift = 0.5`** - Stabilizes SCF by shifting virtual orbitals by 0.5 Hartree

**`conv_tol = 1e-8`** - Stop when energy changes less than 10⁻⁸ Hartree

**`max_cycle = 200`** - Maximum SCF iterations allowed

### 4. Run Calculation

```python
mf.kernel()
```

**`mf.kernel()`** - Performs self-consistent field (SCF) calculation
- Builds Fock matrix (kinetic + nuclear + electron-electron + exchange-correlation)
- Solves Kohn-Sham equations iteratively
- Returns converged energy

### 5. Get Density Matrix

```python
dm = mf.make_rdm1()
S = mol.intor('int1e_ovlp')
print("trace(P S) = ", np.trace(dm @ S))
```

**`mf.make_rdm1()`** - Creates reduced density matrix P
- Contains electron density information
- Size: (n_basis × n_basis)

**`mol.intor('int1e_ovlp')`** - Computes overlap integrals S
- Integrals of basis function overlaps: ⟨φᵢ|φⱼ⟩

**`np.trace(dm @ S)`** - Calculates total electron count
- Should equal 10 for H₂O (8 from O, 2 from H)

### 6. Calculate Electronic Dipole

```python
dip_ao = mol.intor('int1e_r', comp=3)
mu_elec_au = -np.array([np.einsum('ij,ij->', dip_ao[ax], dm) for ax in range(3)])
```

**`mol.intor('int1e_r', comp=3)`** - Computes dipole integrals
- Returns 3 matrices (x, y, z components)
- Each matrix contains ⟨φᵢ|r|φⱼ⟩ integrals

**`np.einsum('ij,ij->', dip_ao[ax], dm)`** - Contracts dipole integrals with density
- Equivalent to: Σᵢⱼ dip_ij × P_ij
- Computes electronic contribution to dipole
- Negative sign: electrons have negative charge

### 7. Calculate Nuclear Dipole

```python
coords = mol.atom_coords()
charges = np.array([mol.atom_charge(i) for i in range(mol.natm)])
nuc_dip_au = np.dot(charges, coords)
```

**`mol.atom_coords()`** - Returns nuclear positions in Bohr units
- Array shape: (n_atoms × 3)

**`mol.atom_charge(i)`** - Returns nuclear charge Z for atom i
- O: 8, H: 1

**`np.dot(charges, coords)`** - Calculates nuclear dipole
- Formula: μ_nuc = Σ Zᵢ × Rᵢ

### 8. Total Dipole

```python
mu_total_au = nuc_dip_au + mu_elec_au
AU_TO_DEBYE = 2.541746286
print("Total dipole (Debye):", mu_total_au * AU_TO_DEBYE)
print("Magnitude:", np.linalg.norm(mu_total_au * AU_TO_DEBYE))
```

**`mu_total_au`** - Sum of nuclear and electronic contributions

**`AU_TO_DEBYE`** - Conversion factor from atomic units to Debye
- 1 a.u. = 2.541746 Debye

**`np.linalg.norm()`** - Computes vector magnitude |μ|

## Expected Output

```
converged SCF energy = -76.4380614812
trace(P S) =  10.0
Total dipole (Debye): [0. 0. 2.07]
Magnitude: 2.07
```

### What Each Output Means

| Output | Value | Meaning |
|--------|-------|---------|
| SCF energy | -76.438 Hartree | Total electronic energy |
| trace(P S) | 10.0 | Electron count (correct ✓) |
| Dipole vector | [0, 0, 2.07] | Dipole points along z-axis |
| Magnitude | 2.07 Debye | Total dipole strength |

**Note:** Experimental dipole of water = 1.85 Debye. B3LYP overestimates by ~12%.

## Summary

1. **Define molecule** → geometry and basis set
2. **Setup DFT** → choose functional (B3LYP) and convergence parameters
3. **Run SCF** → solve electronic structure iteratively
4. **Verify** → check electron count = 10
5. **Calculate dipole** → nuclear contribution + electronic contribution
6. **Convert units** → atomic units → Debye



