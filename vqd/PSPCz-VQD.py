
# coding: utf-8

# ## import

# In[1]:


import matplotlib.pyplot as plt
import numpy as np
from numpy import array
#from scipy import optimize
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

# mario code
from active_space import overwrite_molecule
from pyscf import gto, scf, ao2mo, symm, mcscf


# ## Load Backends / Make Quantum Instances

# In[7]:


#IBMQ.load_account()

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


def overlap(q0, q1):
    """ calculate |<0|q1^dagger*q0|0>|^2
    Args:
        q0 <qiskit.circuit.quantumcircuit.QuantumCircuit>: circuit to make state_0
        q1 <qiskit.circuit.quantumcircuit.QuantumCircuit>: circuit to make state_1
    Returns:
        <float>: |<0|q1^dagger*q0|0>|^2
    """
    global j_call

    # create circuit
    q1_inv = q1.inverse()
    qc = q0 + q1_inv   
    c = ClassicalRegister(qc.n_qubits)
    qc.add_register(c)
    qc.measure(qc.qubits, c)

    result = quantum_instance.execute(qc); j_call += 1 # job_call counter
    count = result.get_counts(0) # count = result.get_counts()

    zeros_key = "0" * qc.n_qubits
    if not(zeros_key in count): count.setdefault(zeros_key, 0)

    try:
        zero_count = float(count[zeros_key])
        return zero_count / sum([value for value in count.values()]) 
    except KeyError:
        return 0.

def expectation_value(Op, circuit):
    global j_call
    eval_circuits = Op.construct_evaluation_circuit(wave_function=circuit, statevector_mode=False)
    result = quantum_instance.execute(eval_circuits); j_call += 1 # job_call counter
    return Op.evaluate_with_result(result=result, statevector_mode=False)


# In[10]:


def expectation_value_H(params, n): # hamiltonian expectation value : <H>
    circ = init_circs[n] + var_form.construct_circuit(params)
    h_val = expectation_value(qubitOp, circ)[0].real
    return h_val 

def operator_penalty(params, n):
    circ = init_circs[n] + var_form.construct_circuit(params)
    op_penalty = 0.0
    if cost_ops is not None:
        for coef, Op in cost_ops:
            if coef != 0.0: op_penalty += coef * expectation_value(Op, circ)[0].real
    return op_penalty

def overlap_penalty(params, n):
    circ = init_circs[n] + var_form.construct_circuit(params)
    ov_penalty = 0.0
    for m in range(n):
        circ_m = init_circs[m] + var_form.construct_circuit(params_list[m]) 
        ov_penalty += weights[m] * overlap(circ, circ_m)
#        print("caoch ov_penalty",ov_penalty)
    return ov_penalty


# In[60]:


def cost_function(params, n, output=True):
    global f_call
    f_call += 1

    h_val = expectation_value_H(params, n)
    op_penalty = operator_penalty(params, n)
    ov_penalty = overlap_penalty(params, n)
    cost = h_val + op_penalty + ov_penalty

    if output:
        if (optimizer == "SPSA"): 
            print('------( function-call #%d, job_call #%d )--------------------------------------------' % (f_call, j_call))
            history = (f_call, cost, h_val, params.copy())
        else: 
            print('------( iteration #%d, function-call #%d, job_call #%d )--------------------------------------------' % (itr, f_call, j_call))
            history = (itr, cost, h_val, params.copy())
        print(' <H> = %.12f : operator_penalty = %.12f : overlap_penalty = %.12f' % (h_val, op_penalty, ov_penalty))
        print('cost = %.12f : params = %s' % (cost, repr(params)))
#        history = (itr, cost, h_val, params.copy())
        histories.append(history)
    return cost

def derivative_cost(params, n):
    shift = np.pi/2.
    shiftgrad = np.zeros(len(params))
    params_shifted = params.copy()
    for idx in range(len(params)):            
        params_shifted[idx] += shift
        plus = cost_function(params_shifted, n, output=False)
        params_shifted[idx] -= 2*shift
        minus = cost_function(params_shifted, n, output=False)
        params_shifted[idx] += shift
        shiftgrad[idx] = 0.5*(plus-minus)
    return shiftgrad


# In[46]:


def find_state(n):
    global params_list

    init_params = params_list[n].copy() # initial guess

    if (optimizer == "SLSQP"):
        opt_params = find_state_SLSQP(n, init_params, options_SLSQP)
    elif (optimizer == "SPSA"):
        opt_params = find_state_SPSA(n, init_params, options_SPSA)
    else:
        print(">>> No such optimizer !!!")
        return 
    
    params_list[n] = opt_params # overwrite :: initial guess => optimized values 


# In[44]:


def find_state_SLSQP(n, init_params, options):

    def callback(_): # just for iteration counter
        global itr
        itr += 1

    opt = minimize(cost_function, init_params, args=(n,), method=optimizer, options=options, jac=derivative_cost, callback=callback)
    print(opt)
    return opt.x


# In[45]:


def find_state_SPSA(n, init_params, options):
    
    fargs = (n,)

    spsa_max_trials = options["max_trials"]
    spsa_parameters = np.array([options["c0"], options["c1"], options["c2"], options["c3"], options["c4"]])

    print(">>> Enter SPSA parameter calibration...")
    num_steps_calibration = min(25, max(1, spsa_max_trials // 5))
    calibration(cost_function, init_params, num_steps_calibration, spsa_parameters, fargs)
    print("<<< Finish SPSA parameter calibration...", spsa_parameters)

    save_steps = 1
    last_avg=1
    (_, x, _, _, _, _) = spsa_optimization(cost_function, init_params, 
                         spsa_max_trials, save_steps, last_avg, spsa_parameters, fargs)
    return x


# ## Preparation

# In[13]:

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
qubitOp.simplify() 
n_qubits = qubitOp.num_qubits

## prepare conserved quantities
tot_N  = ferOp.total_particle_number().mapping(map_type=map_type, threshold=thre)
tot_Ssq = ferOp.total_angular_momentum().mapping(map_type=map_type, threshold=thre)
if qubit_reduction and map_type == "parity":
    tot_N   = Z2Symmetries.two_qubit_reduction(tot_N, num_particles)
    tot_Ssq = Z2Symmetries.two_qubit_reduction(tot_Ssq, num_particles)
    

# ## VQD execution

# In[29]:


def spsa_optimization(obj_fun, initial_theta, max_trials, save_steps, last_avg, spsa_c, spsa_fargs):
    """Minimizes obj_fun(theta) with a simultaneous perturbation stochastic
    approximation algorithm.

    Args:
        obj_fun (callable): the function to minimize
        initial_theta (numpy.array): initial value for the variables of
            obj_fun
        max_trials (int) : the maximum number of trial steps ( = function
            calls/2) in the optimization
        save_steps (int) : stores optimization outcomes each 'save_steps'
            trial steps
        last_avg (int) : number of last updates of the variables to average
            on for the final obj_fun
    Returns:
        list: a list with the following elements:
            cost_final : final optimized value for obj_fun
            theta_best : final values of the variables corresponding to
                cost_final
            cost_plus_save : array of stored values for obj_fun along the
            optimization in the + direction
            cost_minus_save : array of stored values for obj_fun along the
                optimization in the - direction
            theta_plus_save : array of stored variables of obj_fun along the
                optimization in the + direction
            theta_minus_save : array of stored variables of obj_fun along the
                optimization in the - direction
    """

    _parameters = spsa_c
    theta_plus_save = []
    theta_minus_save = []
    cost_plus_save = []
    cost_minus_save = []
    theta = initial_theta
    theta_best = np.zeros(initial_theta.shape)
    for k in range(max_trials):
        # SPSA Parameters
        a_spsa = float(_parameters[0]) / np.power(k + 1 + _parameters[4], _parameters[2])
        c_spsa = float(_parameters[1]) / np.power(k + 1, _parameters[3])
        delta = 2 * aqua_globals.random.randint(2, size=np.shape(initial_theta)[0]) - 1
        # plus and minus directions
        theta_plus = theta + c_spsa * delta
        theta_minus = theta - c_spsa * delta
        cost_plus  = obj_fun(theta_plus,  *spsa_fargs)
        cost_minus = obj_fun(theta_minus, *spsa_fargs)
        # derivative estimate
        g_spsa = (cost_plus - cost_minus) * delta / (2.0 * c_spsa)
        # updated theta
        theta = theta - a_spsa * g_spsa
        # saving
        if k % save_steps == 0:
#            print('Objective function at theta+ for step # %s: %1.7f', k, cost_plus)
#            print('Objective function at theta- for step # %s: %1.7f', k, cost_minus)
            theta_plus_save.append(theta_plus)
            theta_minus_save.append(theta_minus)
            cost_plus_save.append(cost_plus)
            cost_minus_save.append(cost_minus)

        if k >= max_trials - last_avg:
            theta_best += theta / last_avg
    # final cost update
    cost_final = obj_fun(theta_best, *spsa_fargs)

    return [cost_final, theta_best, cost_plus_save, cost_minus_save, theta_plus_save, theta_minus_save]

def calibration(obj_fun, initial_theta, stat, spsa_c, spsa_fargs):
    """Calibrates and stores the SPSA parameters back.

    SPSA parameters are c0 through c5 stored in parameters array

    c0 on input is target_update and is the aimed update of variables on the first trial step.
    Following calibration c0 will be updated.

    c1 is initial_c and is first perturbation of initial_theta.

    Args:
        obj_fun (callable): the function to minimize.
        initial_theta (numpy.array): initial value for the variables of
            obj_fun.
        stat (int) : number of random gradient directions to average on in
            the calibration.
    """

    _parameters = spsa_c
    target_update = _parameters[0]
    initial_c = _parameters[1]
    delta_obj = 0

    for i in range(stat):
        delta = 2 * aqua_globals.random.randint(2, size=np.shape(initial_theta)[0]) - 1
        theta_plus  = initial_theta + initial_c * delta
        theta_minus = initial_theta - initial_c * delta
        obj_plus  = obj_fun(theta_plus, *spsa_fargs)
        obj_minus = obj_fun(theta_minus,*spsa_fargs)
        delta_obj += np.absolute(obj_plus - obj_minus) / stat

    _parameters[0] = target_update * 2 / delta_obj * _parameters[1] * (_parameters[4] + 1)


# ### Settings

# In[22]:


##################################################
#
#    >>>> NOTICE :: GLOBAL variables referenced in VQD_lib functions
#    >>>>               These vars should be set before running VQD !!!
#
#  qubitOp
#  var_form
#  optimizer
#  quantum_instance
#  weights, cost_ops
#  params_list
#  init_circs


# In[40]:


# default options for Optimizers

#optimizer = "SLSQP" 
options_SLSQP = {"disp": True, "maxiter": 100, "ftol": 1e-4} #options = {"disp": True, "maxiter": 1000, "ftol": 1e-15}

#optimizer = "SPSA" 
options_SPSA = {"max_trials": 100, "c0": 0.06283185307179586, "c1": 0.1, "c2": 0.602, "c3": 0.101, "c4": 0}


# In[55]:


# Ansatz for VQD
var_form = RY(n_qubits, depth=1, entanglement="linear") 
n_params = var_form._num_parameters

# known / initial guess Ansatz parameters => params_list[i] for target state i will be 'overwritten' by optimized values
params_list = [ # for target = 1, 2
    array([ 2.02691544, -0.07453562,  0.97759513, -0.01523947]), # <= known Ground state (from VQE)
            0.2 * np.pi * np.random.rand(n_params),             # <= random initial guess for 1st Excited state
            0.2 * np.pi * np.random.rand(n_params),             # <= random initial guess for 2nd Excited state
    ]

k = len(params_list) # number of target states for VQD

# initial state for preparing Ansatz :
init_circs = [QuantumCircuit(n_qubits) for _ in range(0, k)] # k * [QuantumCircuit(n_qubits)] won't work!

weights = 2.*np.ones(k-1) # weights for overlap_penalty
coef = 0.  # coef for S^2 penalty -> triplet 
#coef = 1. # coef for S^2 penalty -> singlet

# penalty term
cost_ops = [[coef, tot_Ssq]]

# ### run

# In[61]:

# run & show
def run_VQD(target=0):
    global itr, f_call, j_call, histories
    print ('[ VQD-run : target = %d ]' % target)
    (itr, f_call, j_call, histories) = (0, 0, 0, [])
    find_state(n=target) 
    print('[ History ]')
    for h in unique(histories):print_history(h)

# utilities for output history
def unique(histories):
    def check_counter(first, second): (i, _, _, _) = first; (j, _, _, _) = second; return i != j
    histries_ = histories+[histories[0]]
    indice = [i for i in range(len(histries_)-1) if check_counter(*histries_[i:i+2])]
    return [histries_[index] for index in indice]

def print_history(history):
    (itr, cost, h_val, params) = history
    print('#%d : cost = %.12f : <H> = %.12f : params = %s' % (itr, cost, h_val, repr(params))) 



# backend and mitgation on/off
quantum_instance = QI_qasm_mitigated
optimizer = "SPSA"
params_list = [
    array([ 1.64362824, -0.11931101,  1.32144647, -0.04442336]), # <= given Ground state (from VQE)
    array([ 1.68937609,  1.6551028 , -0.39322725, -1.23607886]), # <= initial guess for 1st Excited state (from VQD)
    array([ 4.6238072 , 1.69426621, 3.72313986, 2.14560067]), # <= initial guess for 2st Excited state 
    ]
print("calculation for ES1")
run_VQD(target=1)


