import numpy as np

def print_molecule(m):

 print("===== Driver =====")
 print("origin_driver_name ",m.origin_driver_name)
 print("origin_driver_version ",m.origin_driver_version)
 print("origin_driver_config",m.origin_driver_config)
 
 print("===== Energies and orbits =====")
 print("hf_energy",m.hf_energy)
 print("nuclear_repulsion_energy",m.nuclear_repulsion_energy)
 print("num_orbitals",m.num_orbitals)
 print("num_alpha",m.num_alpha)
 print("num_beta",m.num_beta)
 print("mo_coeff",m.mo_coeff)
 print("mo_coeff_B",m.mo_coeff_B)
 print("orbital_energies",m.orbital_energies)
 print("orbital_energies_B",m.orbital_energies_B)

 print("===== Geometry, coords are in Bohr =====")
 print("molecular_charge",m.molecular_charge)
 print("multiplicity",m.multiplicity)
 print("num_atoms",m.num_atoms)
 print("atom_symbol",m.atom_symbol)
 print("atom_xyz",m.atom_xyz)

 print("===== 1 and 2 electron ints in AO basis =====")
 print("hcore",m.hcore)
 print("hcore_B",m.hcore_B)
 print("kinetic",m.kinetic)
 print("overlap",m.overlap)
 print("eri",m.eri)

 print("===== 1 and 2 electron ints in MO basis =====")
 print("mo_onee_ints",m.mo_onee_ints)
 print("mo_onee_ints_B",m.mo_onee_ints_B)
 print("mo_eri_ints",m.mo_eri_ints)
 print("mo_eri_ints_BB",m.mo_eri_ints_BB)
 print("mo_eri_ints_BA",m.mo_eri_ints_BA)

 print("===== dipole moment in AO basis =====")
 print("x_dip_ints",m.x_dip_ints)
 print("y_dip_ints",m.y_dip_ints)
 print("z_dip_ints",m.z_dip_ints)

 print("===== dipole moment in MO basis =====")
 print("x_dip_mo_ints",m.x_dip_mo_ints)
 print("x_dip_mo_ints_B",m.x_dip_mo_ints_B)
 print("y_dip_mo_ints",m.y_dip_mo_ints)
 print("y_dip_mo_ints_B",m.y_dip_mo_ints_B)
 print("z_dip_mo_ints",m.z_dip_mo_ints)
 print("z_dip_mo_ints_B",m.z_dip_mo_ints_B)
 print("nuclear_dipole_moment",m.nuclear_dipole_moment)
 print("reverse_dipole_sign",m.reverse_dipole_sign)

# ===================================================== #

def overwrite_molecule(filename,m):

 import numpy as np
 from pyscf import scf,gto,ao2mo,cc

 #============================================#
 # input file, must contain                   #
 # nbasis,nup,ndown                           #
 # E0                                         #
 # i,j,h[i,j]                                 #
 # p,r,q,s h[p,r,q,s] in Chemist's notation   #
 # in an orthonormal basis                    #
 #============================================#
 linf = open(filename,'r').read().split('\n')
 linf = [x.split(',') for x in linf]
 n,na,nb = linf[0]
 Enuc    = linf[1][0]
 n,na,nb = int(n),int(na),int(nb)
 Enuc    = float(Enuc)
 s1 = np.zeros((n,n))
 h1 = np.zeros((n,n))
 h2 = np.zeros((n,n,n,n))
 count = 2
 for mu in range(n**2):
  II,JJ = linf[count][0],linf[count][1]
  II,JJ = int(II),int(JJ)
  hij   = linf[count][2]
  hij   = float(hij)
  if(II==JJ): s1[II,JJ] = 1
  else:       s1[II,JJ] = 0
  h1[II,JJ] = hij
  count += 1
 for mu in range(n**4):
  PP,RR,QQ,SS = linf[count][0],linf[count][1],linf[count][2],linf[count][3]
  PP,RR,QQ,SS = int(PP),int(RR),int(QQ),int(SS)
  vprqs = linf[count][4]
  vprqs = float(vprqs)
  h2[PP,RR,QQ,SS] = vprqs
  count += 1

 mol = gto.M(verbose=2)
 mol.charge = 0
 mol.nelectron = na+nb
 mol.spin = na-nb
 mol.incore_anyway = True

 if(mol.spin==0): mf = scf.RHF(mol)
 if(mol.spin==1): mf = scf.ROHF(mol)

 mf.get_hcore  = lambda *args: h1
 mf.get_ovlp   = lambda *args: s1
 mf._eri       = ao2mo.restore(8,h2,n)
 mf.init_guess = '1e'
 E0 = mf.kernel()
 if(not mf.converged):
  mf = scf.newton(mf)
  E0 = mf.kernel()

 if(n<11):
  from pyscf import mcscf
  mycas = mcscf.CASSCF(mf,ncas=n,nelecas=mol.nelectron)
  E1 = mycas.kernel()[0]
  print("PySCF energies: SCF, FCI %f %f \n " % (E0+Enuc,E1+Enuc))
 else:
  print("PySCF energy: SCF %f " % (E0+Enuc))

 m.origin_driver_name = 'user-defined'
 m.origin_driver_version = None
 m.origin_driver_config = None

 m.hf_energy = E0+Enuc
 m.nuclear_repulsion_energy = Enuc
 m.num_orbitals = n
 m.num_alpha = na
 m.num_beta = nb
 m.mo_coeff = mf.mo_coeff
 m.mo_coeff_B = None
 m.orbital_energies = mf.mo_energy
 m.orbital_energies_B = None

 m.molecular_charge = None
 m.multiplicity = None
 m.num_atoms = None
 m.atom_symbol = None
 m.atom_xyz = None

 m.hcore = h1
 m.hcore_B = None
 m.kinetic = np.zeros((n,n))
 m.overlap = np.eye(n)
 m.eri = h2

 m.mo_onee_ints = h1
 m.mo_onee_ints_B = None
 m.mo_eri_ints = h2
 m.mo_eri_ints_BB = None
 m.mo_eri_ints_BA = None

 m.mo_onee_ints = np.einsum('pi,pr->ir',m.mo_coeff,m.mo_onee_ints)
 m.mo_onee_ints = np.einsum('ra,ir->ia',m.mo_coeff,m.mo_onee_ints)

 m.mo_eri_ints = np.einsum('pi,prqs->irqs',m.mo_coeff,m.mo_eri_ints)
 m.mo_eri_ints = np.einsum('ra,irqs->iaqs',m.mo_coeff,m.mo_eri_ints)
 m.mo_eri_ints = np.einsum('qj,iaqs->iajs',m.mo_coeff,m.mo_eri_ints)
 m.mo_eri_ints = np.einsum('sb,iajs->iajb',m.mo_coeff,m.mo_eri_ints)

 m.x_dip_ints = None
 m.y_dip_ints = None
 m.z_dip_ints = None

 m.x_dip_mo_ints = None
 m.x_dip_mo_ints_B = None
 m.y_dip_mo_ints = None
 m.y_dip_mo_ints_B = None
 m.z_dip_mo_ints = None
 m.z_dip_mo_ints_B = None
 m.nuclear_dipole_moment = None
 m.reverse_dipole_sign = None

