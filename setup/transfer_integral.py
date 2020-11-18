import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib import colors
C_dimer = np.loadtxt('dimer.ao_coeff.txt')
C_frag1 = np.loadtxt('frag1.ao_coeff.txt')
C_frag2 = np.loadtxt('frag2.ao_coeff.txt')
E_dimer = np.loadtxt('dimer.eigenvalue.txt')
C_loc = np.zeros([24, 24])
C_loc[0:0+C_frag1.shape[0],0:0+C_frag1.shape[1]]=C_frag1
C_loc[12:12+C_frag2.shape[0],12:12+C_frag2.shape[1]]=C_frag2
C_loc_t = C_loc.transpose()
C_loc_dym = C_loc_t.dot(C_dimer)
E_dimer_diag = np.diag(E_dimer)
F_int = C_loc_dym.dot(E_dimer_diag)
F_loc = F_int.dot(C_loc_dym)
print(F_loc[5,17])




