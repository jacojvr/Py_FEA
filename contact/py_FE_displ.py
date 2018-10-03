#%%
import numpy as np
import sys
sys.path.append('../0')
import os
sys.path.append(os.path.abspath('../0'))

from importlib import reload

from numpy import array, c_, matrix,  zeros
from numpy.linalg import norm


#%%

filename0 = 'contact_four_elements'
filename0 = 'contact_four_cpe4'

filename = 'ccx_examples/'+filename0


import setup_from_ccx_inp as FE_io
reload(FE_io)

#bvp_dict =  FE_io.bvp_dict(filename)
#bvp_dict.steps[1]['*nlgeom'] = False

bvp_dict =  FE_io.bvp_dict(filename)
#bvp_dict.steps[1]['*nlgeom'] = True

#%%


import FE_bvp
reload(FE_bvp)
bvp_dict =  FE_io.bvp_dict(filename)



#%%
soln = FE_bvp.solution(bvp_dict)


#%%


soln.solve_step()



#%%
# Write to output
import FE_io as FE_io2
reload(FE_io2)
FE_io2.write_to_vtu(soln,filename0,remove_old=True)




#%%

timef = 1
suse = soln
# test full residual tangent matrix 
residual,k_global,new_ISVs = suse.assemble(dU,time=timef)

pertval = 1e-6

k_global_FD = []
for i in range(residual.size):
    dU_FD = dU[suse.active_dofs]
    dU_FD[i] += pertval
    residualFw = suse.calc_residual(dU_FD,time=timef)
    
    dU_FD = dU[suse.active_dofs]
    dU_FD[i] -= pertval
    residualBw = suse.calc_residual(dU_FD,time=timef)
    
    k_global_FD += [0.5*array(residualFw-residualBw).flatten()/pertval]


KA = k_global.todense()
KB = matrix(k_global_FD)

compared_K = KA - KB

maxkb = np.max(KB-KB.T)
print('maxkb  = ',maxkb)

k_an_norm = norm(k_global.todense())
k_fd_norm = norm(k_global_FD)
k_cm_norm = norm(compared_K)

print('norm of AN : %.6e, FD : %.6e, compared : %.6e'%(k_an_norm,k_fd_norm,k_cm_norm))

#%%

from scipy import sparse


dU = soln.U.copy()
dU[soln.x_dof(7)] = -0.01*np.random.rand()
dU[soln.y_dof(7)] = -0.1*np.random.rand()
dU[soln.x_dof(8)] = -0.01*np.random.rand()
dU[soln.y_dof(8)] = -0.1*np.random.rand()
dU[soln.x_dof(18)] = -0.01*np.random.rand()
dU[soln.y_dof(18)] = -0.1*np.random.rand()
dU[soln.x_dof(22)] = -0.01*np.random.rand()
dU[soln.y_dof(22)] = -0.1*np.random.rand()
dU[soln.x_dof(24)] = -0.01*np.random.rand()
dU[soln.y_dof(24)] = -0.1*np.random.rand()

#%%

reload(FE_bvp)
soln2 = FE_bvp.solution(bvp_dict)
residual2,k_global2,new_ISVs2 = soln2.assemble(dU=soln.U)
C_stiff2 = soln2.check_contact(dU=soln.U,only_active=False,only_resid=False)



soln3 = FE_bvp.solution(bvp_dict)
soln3.U = soln.U.copy()
residual3,k_global3,new_ISVs3 = soln3.assemble()
C_stiff3 = soln3.check_contact(only_active=False,only_resid=False)

