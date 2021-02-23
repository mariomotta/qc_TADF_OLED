
# coding: utf-8

# ## import

# In[1]:


import matplotlib.pyplot as plt
import numpy as np
from numpy import array
from scipy import optimize
from scipy.optimize import minimize
import time 
import pickle

from qiskit import QuantumRegister, QuantumCircuit, ClassicalRegister
from qiskit import Aer, IBMQ, execute

from qiskit.ignis.mitigation.measurement import CompleteMeasFitter
from qiskit.ignis.verification.tomography import state_tomography_circuits, StateTomographyFitter
from qiskit.quantum_info.operators import Operator, Pauli

import qiskit.aqua.operators.weighted_pauli_operator
from qiskit.aqua import aqua_globals
from qiskit.aqua import QuantumInstance
from qiskit.aqua.algorithms import ExactEigensolver
from qiskit.aqua.components.variational_forms import VariationalForm, RY, RYRZ
from qiskit.aqua.components.optimizers import (Optimizer, COBYLA, SPSA, NELDER_MEAD, L_BFGS_B, TNC, CG, SLSQP)
from qiskit.aqua.operators import Z2Symmetries, MatrixOperator
from qiskit.aqua.operators.common import evolution_instruction
from qiskit.chemistry import FermionicOperator
from qiskit.chemistry.drivers import PySCFDriver, UnitsType

# mario code >
from active_space import overwrite_molecule
# mario code <
from pyscf import gto, scf, ao2mo, symm, mcscf


# ## Load Backends / Make Quantum Instances

# In[7]:

backend_qasm = Aer.get_backend('qasm_simulator')


# In[8]:


def make_QI_unmitigated(backend):
    return QuantumInstance(backend=backend, shots=8192, skip_qobj_validation=False)
def make_QI_mitigated(backend):
    return QuantumInstance(backend=backend, shots=8192, skip_qobj_validation=False,
           measurement_error_mitigation_cls=CompleteMeasFitter, cals_matrix_refresh_period = 30)

# define Quantum Instances
QI_qasm_unmitigated      = make_QI_unmitigated(backend_qasm)
QI_qasm_mitigated        = make_QI_mitigated(  backend_qasm)


# ## Functions (originally from VQD_lib)

# In[9]:




### way of fermion-qubit mapping
qubit_reduction =  True # valid only for `parity` mapping
map_type = 'parity' 

### detail on whole simulation
basis = "631gs"
num_particles = 2

geometry = f"Li 0.0000000 0.0000000 0.0000000;H         1.5300000   0.0000000   0.0000000"

## run pyscf
thre = 1e-8
driver = PySCFDriver(atom=geometry, unit=UnitsType.ANGSTROM, charge=0, spin=0, basis=basis)
molecule = driver.run()

# ----- changing matrix elements ----- #
overwrite_molecule('h.csv',molecule)

nuclear_repulsion_energy = molecule.nuclear_repulsion_energy

###  ---- construct qubit Hamiltonian ----
h1 = molecule.one_body_integrals
h2 = molecule.two_body_integrals
ferOp = FermionicOperator(h1=h1, h2=h2)
core_energy = 0
    
# GLOBAL var : qubitOp
qubitOp = ferOp.mapping(map_type=map_type, threshold=thre)
if qubit_reduction and map_type == "parity":
    qubitOp = Z2Symmetries.two_qubit_reduction(qubitOp, num_particles)
qubitOp.chop(10**-10)
qubitOp.simplify() ## これがないとバグる
n_qubits = qubitOp.num_qubits

## prepare conserved quantities
tot_N  = ferOp.total_particle_number().mapping(map_type=map_type, threshold=thre)
tot_Ssq = ferOp.total_angular_momentum().mapping(map_type=map_type, threshold=thre)
if qubit_reduction and map_type == "parity":
    tot_N   = Z2Symmetries.two_qubit_reduction(tot_N, num_particles)
    tot_Ssq = Z2Symmetries.two_qubit_reduction(tot_Ssq, num_particles)
    
# for calculating <H> from psi obtained from state tomography 
pauli_table = vars(qubitOp)['_paulis_table']
coeffs = [ np.real(p[0]) for p in vars(qubitOp)['_paulis']]
h_matrix = sum([coeffs[val] * Pauli(label=key).to_matrix() for (key, val) in pauli_table.items()])
def eval_H(psi, h_matrix):
    return np.real(np.conjugate(psi).dot(h_matrix.dot(psi)))


# ## State tomography purification

# In[33]:


# 2qubit RY Ansatz => state vector
# Analytical expressions below are obtained by Mathematica
def ry00(th1, th2, th3, th4): 
    return  np.cos(th1/2)*np.cos(th2/2)*np.cos(th3/2)*np.cos(th4/2) -             np.cos(th2/2)*np.cos(th4/2)*np.sin(th1/2)*np.sin(th3/2) -             np.cos(th1/2)*np.cos(th3/2)*np.sin(th2/2)*np.sin(th4/2) -             np.sin(th1/2)*np.sin(th2/2)*np.sin(th3/2)*np.sin(th4/2)
def ry01(th1, th2, th3, th4): 
    return  np.cos(th2/2)*np.cos(th3/2)*np.cos(th4/2)*np.sin(th1/2) +             np.cos(th1/2)*np.cos(th2/2)*np.cos(th4/2)*np.sin(th3/2) +             np.cos(th3/2)*np.sin(th1/2)*np.sin(th2/2)*np.sin(th4/2) -             np.cos(th1/2)*np.sin(th2/2)*np.sin(th3/2)*np.sin(th4/2)
def ry10(th1, th2, th3, th4): 
    return  np.cos(th1/2)*np.cos(th3/2)*np.cos(th4/2)*np.sin(th2/2) +             np.cos(th4/2)*np.sin(th1/2)*np.sin(th2/2)*np.sin(th3/2) +             np.cos(th1/2)*np.cos(th2/2)*np.cos(th3/2)*np.sin(th4/2) -             np.cos(th2/2)*np.sin(th1/2)*np.sin(th3/2)*np.sin(th4/2)
def ry11(th1, th2, th3, th4): 
    return  -(np.cos(th3/2)*np.cos(th4/2)*np.sin(th1/2)*np.sin(th2/2)) +             np.cos(th1/2)*np.cos(th4/2)*np.sin(th2/2)*np.sin(th3/2) +             np.cos(th2/2)*np.cos(th3/2)*np.sin(th1/2)*np.sin(th4/2) +             np.cos(th1/2)*np.cos(th2/2)*np.sin(th3/2)*np.sin(th4/2)
def ry_state_vec(th1, th2, th3, th4): 
    return np.array([ry00(th1, th2, th3, th4), 
                     ry01(th1, th2, th3, th4), 
                     ry10(th1, th2, th3, th4), 
                     ry11(th1, th2, th3, th4)])

# rho(density matrix) => psi(statevector)
def get_psi(rho): 
    (evals, evecs) = np.linalg.eig(rho) # diagonalize rho
    evals = np.real(evals)              # eigenvalues
    evecs = np.transpose(evecs)         # eigenvectors
    psi = evecs[np.argmax(evals)]       # statevector with max-eigenvalue
    return psi

# 1 - |<Psi|RY(theta1,..,theta4)>|^2 :: minimize this to obtain (theta1,..,theta4)
def diff(psi):
    return lambda x: 1 - np.abs(psi.dot(ry_state_vec(x[0], x[1], x[2], x[3])))**2


# In[35]:


# RY-Ansatz parameters(theta1,..,theta4) obtained by VQD-iteration 
params = np.array([1.70075589,  1.72352203, -0.32227027, -1.34956145]) 

# circuit for 2-qubit RY Ansatz
var_form = RY(2, depth=1, entanglement="linear") 
circ = var_form.construct_circuit(parameters = params)
#circ.draw(output = 'mpl')


# In[36]:


# generate circuits for state tomography 
tomo_circuits = state_tomography_circuits(circ, circ.qregs)


# In[44]:


# exec state tomography experiments
job = execute(tomo_circuits, backend = backend_qasm, shots = 8192)


# In[10]:


job.status()


# In[51]:


# get rho(density matrix) from tomography data
tomo_data = StateTomographyFitter(job.result(), tomo_circuits)
rho = tomo_data.fit()
psi = get_psi(rho)
new_params = optimize.fmin(diff(psi), params) # minimize 1 - |<Psi|RY(theta1,..,theta4)>|^2
print("new params",new_params)
new_params


# In[59]:


new_Hamiltonian = eval_H(psi, h_matrix)
print("new_Hamiltonian",new_Hamiltonian)


# In[58]:


pauli_table = vars(qubitOp)['_paulis_table']
coeffs = [ np.real(p[0]) for p in vars(qubitOp)['_paulis']]
h_matrix = sum([coeffs[val] * Pauli(label=key).to_matrix() for (key, val) in pauli_table.items()])
def eval_H(psi, h_matrix):
    return np.real(np.conjugate(psi).dot(h_matrix.dot(psi)))


# ## Juncyard

# In[6]:


# operator
def to_matrix_operator_from_WeightedPauliOperator(operator):
    """
    Copy and paste from https://github.com/Qiskit/qiskit-aqua/blob/master/qiskit/aqua/operators/op_converter.py of version 0.6.1

    Converting a given operator to `MatrixOperator`
    Args:
        operator (WeightedPauliOperator):
            one of supported operator type
    Returns:
        MatrixOperator: the converted matrix operator
    Raises:
        AquaError: Unsupported type to convert
    """
    if operator.is_empty():
        return MatrixOperator(None)
    hamiltonian = 0
    for weight, pauli in operator.paulis:
        hamiltonian += weight * pauli.to_spmatrix()
    return MatrixOperator(matrix=hamiltonian, z2_symmetries=operator.z2_symmetries,
                              name=operator.name)

def overlap_operator(Op, q0, q1, shots=1024, backend=None, mode="shots", shallow_slicing=False):
    """ calculate |<0|q1^dagger*Op*q0|0>|^2
    Args:
        Op <qiskit.aqua.operators.weighted_pauli_operator.WeightedPauliOperator>: operator to evaluate
        q0 <qiskit.circuit.quantumcircuit.QuantumCircuit>: circuit to make state_0
        q1 <qiskit.circuit.quantumcircuit.QuantumCircuit>: circuit to make state_1
        shots <int>: the number of shots. used if and only if `mode`="shots"
        backend: Qiskit backend. If none and mode="shots", `qasm_simulator` of Qiskit Aqua will be used.
        mode ("shots" | "statevector): way to execute circuit
        shallow_slicing <bool>: option for implementing exp(i*P_i/4)
    Returns:
        <float>: |<0|q1^dagger*Op*q0|0>|^2

    DETAILS:
    Evaluate |<0|q1^dagger*A*q0|0>|^2 for A = \sum_i a_i P_i where a_i is real coefficient and P_i is Pauli operator P_i^2 = I.
    We assume
    - <0|q1^dagger*q0|0> = 0, i.e. q1 and q1 are orthogonal.
    - A is Hermite -> a_i is real.
    Our calculation is based on the following equations.
    |<0|q1^dagger*A*q0|0>|^2 = |<q1|A|q0>|^2
    = \sum_i |a_i|^2 |<q1|P_i|q0>|^2 + \sum_{i<j} a_i a_j 2Re(<q0|P_i|q1><q1|P_j|q0>)  

    U_{ij+} := 1/2*(I+iP_i)*(I+iP_j), U_{ij-} := 1/2*(I-iP_i)*(I-iP_j)
    2*(|<q1|U_{ij+}|q0>|^2 + |<q1|U_{ij-}|q0>|^2)
    =  |<q1|P_i|q0>|^2 + |<q1|P_j|q0>|^2 + |<q1|P_iP_j|q0>|^2
        + 2 * Re(<q0|P_i|q1><q1|P_j|q0>)
    """
    if not isinstance(Op, qiskit.aqua.operators.weighted_pauli_operator.WeightedPauliOperator):
        raise ValueError('operator must be type of qiskit.aqua.operators.weighted_pauli_operator.WeightedPauliOperator')

    ## remove Identity term
    n_qubits = Op.num_qubits
    identity_label = "".join(["I"]*n_qubits)
    n_pauli = 0
    coef_list = []; pauli_list = []
    for i in range(len(Op._paulis)):
        paulis = Op._paulis[i]
        if paulis[1].to_label() != identity_label:
            coef_list.append(paulis[0].real)
            pauli_list.append(paulis[1])
            n_pauli += 1

    q1_Pi_q0_sq = np.zeros(n_pauli, dtype=float)
    q1_PiPj_q0_sq = np.zeros((n_pauli,n_pauli), dtype=float)
    q1_Uij_plus_q0_sq  = np.zeros((n_pauli,n_pauli), dtype=float)
    q1_Uij_minus_q0_sq = np.zeros((n_pauli,n_pauli), dtype=float)
    Re_q0_Pi_q1_q1_Pj_q0 = np.zeros((n_pauli,n_pauli), dtype=float)

    ## determine diagonal part first
    for i in range(n_pauli):
        pauli_i = pauli_list[i]
        q2 = q0.copy()
        q2.append(pauli_i.to_instruction(), q2.qubits)
        q1_Pi_q0_sq[i] = overlap(q2, q1, shots=shots, backend=backend, mode=mode) 
    ## determine off-diagonal part
    for i in range(n_pauli):
        pauli_i = pauli_list[i]
        plus_i  = evolution_instruction([[1, pauli_i]],  np.pi/4., 1, shallow_slicing=shallow_slicing)
        minus_i = evolution_instruction([[1, pauli_i]], -np.pi/4., 1, shallow_slicing=shallow_slicing)
        for j in range(i+1, n_pauli):
            pauli_j = pauli_list[j]
            plus_j  = evolution_instruction([[1, pauli_j]],  np.pi/4., 1, shallow_slicing=shallow_slicing)
            minus_j = evolution_instruction([[1, pauli_j]], -np.pi/4., 1, shallow_slicing=shallow_slicing)
            q2 = q0.copy()
            q2.append(plus_j, q2.qubits)
            q2.append(plus_i, q2.qubits)
            q1_Uij_plus_q0_sq[i,j] = overlap(q1, q2, mode=mode, shots=shots)
            q2 = q0.copy()
            q2.append(minus_j, q2.qubits)
            q2.append(minus_i, q2.qubits)
            q1_Uij_minus_q0_sq[i,j] = overlap(q1, q2, mode=mode, shots=shots)
            q2 = q0.copy()
            q2.append(pauli_j.to_instruction(), q2.qubits)
            q2.append(pauli_i.to_instruction(), q2.qubits)
            q1_PiPj_q0_sq[i,j] =  overlap(q1, q2, mode=mode, shots=shots)
    ## post-process 
    for i in range(n_pauli):
        for j in range(i+1, n_pauli):
            Re_q0_Pi_q1_q1_Pj_q0[i,j] = q1_Uij_plus_q0_sq[i,j] + q1_Uij_minus_q0_sq[i,j] -0.5*(q1_Pi_q0_sq[i] + q1_Pi_q0_sq[j] + q1_PiPj_q0_sq[i,j])

    diag = np.sum(np.abs(coef_list)**2 * q1_Pi_q0_sq)
    off_diag = np.einsum( "i,ij,j", coef_list, Re_q0_Pi_q1_q1_Pj_q0, coef_list)
    return diag + 2*off_diag


# save_str = f"mode={mode}, shots={shots}, seed={my_seed}, shots={shots}, optimizer={optimizer}, options={options}, use_jac={use_jac}, weights={weights}"
# print(save_str)
