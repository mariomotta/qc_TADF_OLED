#!/usr/bin/env python
# coding: utf-8

# In[1]:


#get_ipython().system(' python -V')
#get_ipython().system(' pip list | grep "qiskit"')


# In[2]:



# import common packages
import numpy as np
import pickle
import qiskit

from qiskit import Aer
# lib from Qiskit Aqua

from qiskit import QuantumRegister, QuantumCircuit, ClassicalRegister
from qiskit.providers.aer import (QasmSimulator, StatevectorSimulator, UnitarySimulator, noise)
from qiskit.aqua.utils.backend_utils import is_statevector_backend
from qiskit.aqua.utils import summarize_circuits
from qiskit.ignis.mitigation.measurement import complete_meas_cal, CompleteMeasFitter, MeasurementFilter
from qiskit.quantum_info import state_fidelity, process_fidelity

#from qiskit.aqua import Operator, QuantumInstance
from qiskit.aqua import QuantumInstance
from qiskit.aqua.operators import (TPBGroupedWeightedPauliOperator, WeightedPauliOperator, MatrixOperator, op_converter, Z2Symmetries)
from qiskit.aqua.algorithms import VQE, ExactEigensolver
from qiskit.aqua.operators import Z2Symmetries
from qiskit.aqua.components.optimizers import (Optimizer, COBYLA, SPSA, NELDER_MEAD, L_BFGS_B, TNC, CG)
from qiskit.aqua.components.variational_forms import (ry, ryrz, swaprz, variational_form)
from qiskit.aqua.components.initial_states import (zero, custom, initial_state)
from qiskit.aqua.utils import summarize_circuits
from qiskit.aqua.operators import (TPBGroupedWeightedPauliOperator, WeightedPauliOperator, MatrixOperator, op_converter)

# lib from Qiskit Aqua Chemistry
from qiskit.chemistry.drivers.pyscfd.integrals import compute_integrals
from qiskit.chemistry import FermionicOperator
from qiskit.chemistry.drivers import UnitsType
from qiskit.chemistry.drivers.pyscfd import PySCFDriver
from qiskit.chemistry.drivers.psi4d import PSI4Driver
from qiskit.chemistry.core import Hamiltonian
from qiskit.chemistry.core import TransformationType, QubitMappingType
from qiskit.chemistry.components.variational_forms import UCCSD
from qiskit.chemistry.components.initial_states import HartreeFock
from qiskit.chemistry.algorithms import QEomVQE
# lib from Qiskit tools
from qiskit.tools.visualization import plot_histogram, plot_state_city


# IBM mario code >
from active_space                    import overwrite_molecule
# IBM mario code <
import numpy as np
from pyscf import gto,scf,ao2mo,symm,mcscf

# do a "large" mean-field calculation

mol = gto.M(atom=[
                  ['C', (  -0.2316640,    1.1348450,    0.6956120)],
                  ['C', (  -0.8886300,    0.3253780,   -0.2344140)],
                  ['C', (  -0.1842470,   -0.1935670,   -1.3239330)],
                  ['C', (   1.1662930,    0.0801450,   -1.4737160)],
                  ['C', (   1.8089230,    0.8832220,   -0.5383540)],
                  ['C', (   1.1155860,    1.4218050,    0.5392780)],
                  ['S', (   3.5450920,    1.2449890,   -0.7349240)],
                  ['O', (   3.8606900,    1.0881590,   -2.1541690)],
                  ['C', (   4.3889120,   -0.0620730,    0.1436780)],
                  ['O', (   3.8088290,    2.4916780,   -0.0174650)],
                  ['C', (   4.6830900,    0.1064460,    1.4918230)],
                  ['C', (   5.3364470,   -0.9144080,    2.1705280)],
                  ['C', (   5.6895490,   -2.0818670,    1.5007820)],
                  ['C', (   5.4000540,   -2.2323130,    0.1481350)],
                  ['C', (   4.7467230,   -1.2180160,   -0.5404770)],
                  ['N', (  -2.2589180,    0.0399120,   -0.0793330)],
                  ['C', (  -2.8394600,   -1.2343990,   -0.1494160)],
                  ['C', (  -4.2635450,   -1.0769890,    0.0660760)],
                  ['C', (  -4.5212550,    0.2638010,    0.2662190)],
                  ['C', (  -3.2669630,    0.9823890,    0.1722720)],
                  ['C', (  -2.2678900,   -2.4598950,   -0.3287380)],
                  ['C', (  -3.1299420,   -3.6058560,   -0.3236210)],
                  ['C', (  -4.5179520,   -3.4797390,   -0.1395160)],
                  ['C', (  -5.1056310,   -2.2512990,    0.0536940)],
                  ['C', (  -5.7352450,    1.0074800,    0.5140960)],
                  ['C', (  -5.6563790,    2.3761270,    0.6274610)],
                  ['C', (  -4.4287740,    3.0501460,    0.5083650)],
                  ['C', (  -3.2040560,    2.3409470,    0.2746950)],
                  ['H', (  -0.7813570,    1.5286610,    1.5426490)],
                  ['H', (  -0.7079140,   -0.7911480,   -2.0611600)],
                  ['H', (   1.7161320,   -0.2933710,   -2.3302930)],
                  ['H', (   1.6308220,    2.0660550,    1.2427990)],
                  ['H', (   4.4214900,    1.0345500,    1.9875450)],
                  ['H', (   5.5773000,   -0.7951290,    3.2218590)],
                  ['H', (   6.2017810,   -2.8762260,    2.0345740)],
                  ['H', (   5.6906680,   -3.1381740,   -0.3739110)],
                  ['H', (   4.5337010,   -1.3031330,   -1.6001680)],
                  ['H', (  -1.1998460,   -2.5827750,   -0.4596910)],
                  ['H', (  -2.6937370,   -4.5881470,   -0.4657540)],
                  ['H', (  -5.1332290,   -4.3740010,   -0.1501080)],
                  ['H', (  -6.1752900,   -2.1516170,    0.1987120)],
                  ['H', (  -6.6812260,    0.4853900,    0.6017680)],
                  ['H', (  -6.5574610,    2.9529350,    0.8109620)],
                  ['H', (  -4.3980410,    4.1305040,    0.5929440)],
                  ['H', (  -2.2726630,    2.8838620,    0.1712760)]],
            charge=0,spin=0,basis='631g*',verbose=0)
mf = scf.RHF(mol)
mf = scf.newton(mf)
E  = mf.kernel()

print("SCF energy ",E)

# choose the list of orbitals to be used in the Qiskit calculation
# and the number of spin-up and spin-down electrons

orbitals = [99,100]
na,nb    = 1,1

# the Hamiltonian in the active space of interest is constructed below, and dumped on a file

norb   = len(orbitals)
freeze = [x for x in range(mol.nao_nr()) if mf.mo_occ[x]>0 and x not in orbitals]

mo_core   = mf.mo_coeff[:,freeze]
mo_active = mf.mo_coeff[:,orbitals]

core_dm = np.dot(mo_core,mo_core.T)*2
corevhf = mf.get_veff(mol,core_dm)

ecore  = mol.energy_nuc()
ecore += np.einsum('ij,ji',core_dm,mf.get_hcore())
ecore += np.einsum('ij,ji',core_dm,corevhf)*0.5

h1 = np.einsum('ab,ai,bk->ik',mf.get_hcore()+corevhf,mo_active,mo_active)
h2 = ao2mo.kernel(mol,mo_active)
h2 = ao2mo.restore(1,h2,norb)

fname = 'h.csv'

outf = open(fname,'w')
outf.write('%d,%d,%d\n' % (norb,na,nb))
outf.write('%f\n' % ecore)
for i in range(norb):
 for j in range(norb):
  outf.write('%d,%d,%.16f\n' % (i,j,h1[i,j]))
for p in range(norb):
 for r in range(norb):
  for q in range(norb):
   for s in range(norb):
    outf.write('%d,%d,%d,%d,%.16f \n' % (p,r,q,s,h2[p,r,q,s]))
outf.close()



# In[3]:
# calclulation results on classical computer
class cqr:

    def __init__(self, energy_shift, repulsion, num_AS_spin_orbitals, num_AS_particles):
        self._energy_shift = energy_shift
        self._repulsion = repulsion
        self._num_AS_spin_orbitals = num_AS_spin_orbitals
        self._num_AS_particles = num_AS_particles

# using driver to get fermionic Hamiltonian
# PySCF example
driver = PySCFDriver(atom='Li .0 .0 .0; H 1.530 .0 .0', unit=UnitsType.ANGSTROM, charge=0, spin=0, basis='sto3g')


#driver = PySCFDriver(atom=
#                     'Li 0.0000000  0.0000000  0.0000000;'
#                     'H  1.5300000  0.0000000  0.0000000',
#                     unit=UnitsType.ANGSTROM, charge=0, spin=0, basis='sto3g')
                     
qmolecule = driver.run()


# In[4]:


## CHECK 

vars(qmolecule)

# ----- changing matrix elements ----- #
overwrite_molecule('h.csv',qmolecule)


# In[5]:

# In[6]:

# In[7]:


hamiltonian = Hamiltonian(transformation=TransformationType.FULL,qubit_mapping=QubitMappingType.PARITY,two_qubit_reduction=True,freeze_core=False)
energy_shift = hamiltonian._energy_shift
repulsion = qmolecule.nuclear_repulsion_energy 
num_AS_spin_orbitals = norb + norb 
num_AS_particles = na + nb
classical_info = cqr(energy_shift, repulsion, num_AS_spin_orbitals, num_AS_particles) 
vars(hamiltonian)

# In[8]:


## CREATE WeightedPauliOperater (qubit_op)
## CAUTION: Take a long time

qubit_op, aux_ops = hamiltonian.run(qmolecule=qmolecule)


# In[9]:


vars(qubit_op)


# In[10]:


# FILE IO
wpOp = op_converter.to_weighted_pauli_operator(qubit_op)
for p in wpOp.paulis:
    p[0] = np.real(p[0])
with open('./PSPCz.hamiltonian', 'wb') as f:
    pickle.dump(classical_info, f)
wpOp.to_file('./PSPCz.wpOp')


# In[11]:


# Using exact eigensolver to get the smallest eigenvalue
energy_shift=hamiltonian._energy_shift
exact_eigensolver = ExactEigensolver(qubit_op, k=3)
ret = exact_eigensolver.run()


# In[12]:


# Using exact eigensolver to get the smallest eigenvalue
exact_eigensolver = ExactEigensolver(qubit_op, k=4)
ret = exact_eigensolver.run()
print('The computed E0 electronic energy is: %.12f' % ret['eigvals'][0].real)
print('The computed E1 electronic energy is: %.12f' % ret['eigvals'][1].real)
print('The computed E2 electronic energy is: %.12f' % ret['eigvals'][2].real)
print('The computed E3 electronic energy is: %.12f' % ret['eigvals'][3].real)
print('energy_shift: %.12f' % energy_shift)
print('nuclear_repulsion_energy: %.12f' % qmolecule.nuclear_repulsion_energy)
print('The total E0 state energy is: % .12f' % (ret['eigvals'][0].real + energy_shift + qmolecule.nuclear_repulsion_energy))
print('The total E1 state energy is: % .12f' % (ret['eigvals'][1].real + energy_shift + qmolecule.nuclear_repulsion_energy))
print('The total E2 state energy is: % .12f' % (ret['eigvals'][2].real + energy_shift + qmolecule.nuclear_repulsion_energy))
print('The total E3 state energy is: % .12f' % (ret['eigvals'][3].real + energy_shift + qmolecule.nuclear_repulsion_energy))
print('The excited energy E1-E0: % .12f a.u ' % (ret['eigvals'][1].real -ret['eigvals'][0].real))
print('The excited energy E2-E0: % .12f a.u ' % (ret['eigvals'][2].real -ret['eigvals'][0].real))
print('The excited energy E3-E0: % .12f a.u ' % (ret['eigvals'][3].real -ret['eigvals'][0].real))

print('The excited energy E2-E1: % .12f a.u ' % (ret['eigvals'][2].real -ret['eigvals'][1].real))


