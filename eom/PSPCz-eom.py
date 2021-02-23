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
import logging

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
from qiskit.aqua.components.optimizers import (Optimizer, COBYLA, SPSA, NELDER_MEAD, L_BFGS_B, TNC, CG, SLSQP)
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
from qiskit.aqua.components.variational_forms import RY, RYRZ
from qiskit.chemistry.components.initial_states import HartreeFock
from qiskit.chemistry.algorithms import QEomVQE
from qiskit.chemistry import QMolecule, set_qiskit_chemistry_logging

# lib from Qiskit tools
#from qiskit.tools.visualization import plot_histogram, plot_state_city


# In[3]:
# calclulation results on classical computer
class cqr:

    def __init__(self, energy_shift, repulsion, num_AS_spin_orbitals, num_AS_particles):
        self._energy_shift = energy_shift
        self._repulsion = repulsion
        self._num_AS_spin_orbitals = num_AS_spin_orbitals
        self._num_AS_particles = num_AS_particles


# In[4]:

# In[5]:

# In[6]:

with open('./PSPCz.hamiltonian', 'rb') as f:
    hamiltonian = pickle.load(f)
energy_shift = hamiltonian._energy_shift
repulsion = hamiltonian._repulsion
num_AS_spin_orbitals = hamiltonian._num_AS_spin_orbitals
num_AS_particles     = hamiltonian._num_AS_particles
print("number of spin orbitals in active space: {}".format(num_AS_spin_orbitals))
print("number of electrons in active space: {}".format(num_AS_particles))
vars(hamiltonian)

# In[7]:

# In[8]:

# In[9]:

# In[10]:


# FILE IO
qubit_op = WeightedPauliOperator.from_file('./PSPCz.wpOp')
vars(qubit_op)

# In[11]:


# Using exact eigensolver to get the smallest eigenvalue
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
print('nuclear_repulsion_energy: %.12f' % repulsion)
print('The total E0 state energy is: % .12f' % (ret['eigvals'][0].real + energy_shift + repulsion))
print('The total E1 state energy is: % .12f' % (ret['eigvals'][1].real + energy_shift + repulsion))
print('The total E2 state energy is: % .12f' % (ret['eigvals'][2].real + energy_shift + repulsion))
print('The total E3 state energy is: % .12f' % (ret['eigvals'][3].real + energy_shift + repulsion))
print('The excited energy E1-E0: % .12f a.u ' % (ret['eigvals'][1].real -ret['eigvals'][0].real))
print('The excited energy E2-E0: % .12f a.u ' % (ret['eigvals'][2].real -ret['eigvals'][0].real))
print('The excited energy E3-E0: % .12f a.u ' % (ret['eigvals'][3].real -ret['eigvals'][0].real))

print('The excited energy E2-E1: % .12f a.u ' % (ret['eigvals'][2].real -ret['eigvals'][1].real))


# In[25]:


# from qiskit import IBMQ
# IBMQ.load_accounts()
#Choose backend, qasm_simulator, qasm_simulator_py, state_vector_simulator, state_vector_simulator_py, unitary_simulator, clifford_simulator
backend = Aer.get_backend('qasm_simulator') 
#backend = Aer.get_backend('statevector_simulator') 


# In[26]:


#COBYLA, CG, L-BFGS-B, P-BFGS, Powell, TNC, Nealder-Mead, 

# COBYLA (constrained Optimization BY Linear Approximation) optimizer
#opt = COBYLA(maxiter=200)

# setup CG (conjugate gradient) optimizer
#opt = CG(maxiter=1000)

# setup SLSQP (conjugate gradient) optimizer
opt = SLSQP(maxiter=1000, ftol=1e-9 )

#opt = SPSA(max_trials=1000, save_steps=1, last_avg=10)


# In[27]:

# setup HartreeFock state
HF_state = HartreeFock(qubit_op.num_qubits, 
                       num_orbitals=num_AS_spin_orbitals, 
                       num_particles=num_AS_particles, 
                       qubit_mapping='parity',  
                       two_qubit_reduction=True)

# setup variational form
#var_form = UCCSD(qubit_op.num_qubits, depth=1, 
#                   num_orbitals=num_AS_spin_orbitals, num_particles=num_AS_particles, 
#                   active_occupied=[0], active_unoccupied=[0],
#                   initial_state=HF_state, qubit_mapping='parity', 
#                   two_qubit_reduction=True, num_time_slices=1)
var_form = RY(qubit_op.num_qubits, depth=1, entanglement="linear") 

params = np.zeros(var_form.num_parameters)

def store_intermediate_result(eval_count, params, mean, std):
    with open("intermediate.xvg", 'a') as fp:
        #paramstr = ' '.join([str(p) for p in params])
        paramstr=' '
        for p in params:
            paramstr += format(p, '10f')+ ' ' 
        content = "{:4d}  {:16.12f} {:8.5f} [{}]".format(eval_count, mean, std, paramstr)
        print(content, file=fp, flush=True)
        print(content)


# In[16]:


vars(var_form)


# In[28]:


# setup VQE
# vqe = VQE(qubitOp, var_form, cobyla, 'matrix')
operator_mode = 'grouped_paulis'
params=[ 1.53203614, -0.04493598,  1.55716405, -0.04073203] # VQD SPSA GS on qasm

#vqe = VQE(qubitOp, var_form, opt, operator_mode, initial_point=params, callback=store_intermediate_result) # redundant initial_point

vqe = QEomVQE(qubit_op, 
              var_form, 
              opt, 
              num_orbitals=num_AS_spin_orbitals, 
              num_particles=num_AS_particles, 
              initial_point=params, 
              max_evals_grouped=1, 
              callback=store_intermediate_result, 
              auto_conversion=True, 
              qubit_mapping='parity', 
              two_qubit_reduction=True, 
              is_eom_matrix_symmetric=True, 
              active_occupied=None, 
              active_unoccupied=None, 
              se_list=None, 
              de_list=None, 
              z2_symmetries=None, 
              untapered_op=None, 
              aux_operators=None)

quantum_instance = QuantumInstance(backend=backend, shots=8192)


# In[29]:


results = vqe.run(quantum_instance)


# In[30]:



print('The computed state energies: {}'.format(results['energies']))
print('The total E0 state energy is: % .12f' % (results['energies'][0] + energy_shift + repulsion))
print('The total E1 state energy is: % .12f' % (results['energies'][1] + energy_shift + repulsion))
print('The total E2 state energy is: % .12f' % (results['energies'][2] + energy_shift + repulsion))
print('The computed energy gaps from E0: {} Hartree'.format(results['energy_gap']))
print('The computed E1-E0 excited energy is: {} kcal/mol'.format(results['energy_gap'][0]*627.5096))


# In[20]:

