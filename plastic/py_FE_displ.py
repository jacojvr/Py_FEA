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

#from time import time

filename0 = 'q8_patch_s11'

filename0 = 'Cook_8_CPE8R_ccx_plast'

filename = 'ccx_examples/'+filename0


import FE_io
reload(FE_io)


#bvp_dict =  FE_io.bvp_dict(filename)
#bvp_dict.steps[1]['*nlgeom'] = False

bvp_dict =  FE_io.bvp_dict(filename)
bvp_dict.steps[1]['*nlgeom'] = True

#%%


import FE_bvp
reload(FE_bvp)

# initialise solution
#soln = FE_bvp.solution(bvp_dict)

soln = FE_bvp.solution(bvp_dict)


#%%
timef=1

# Solve unknowns for full time step
#dU,new_ISVs,is_converged = soln.solve_increment(time=timef)
#
dU,new_ISVs,is_converged = soln.solve_increment(time=timef)
#

#%%
# again 
#over 10 steps


bvp_dict =  FE_io.bvp_dict(filename)
bvp_dict.steps[1]['*nlgeom'] = True# False
soln = FE_bvp.solution(bvp_dict)
dtime = 0.1


dU,new_ISVs,is_converged,n_iter = soln.solve_increment(time=0.001)
dU*=100
    
timef = 0
for i in range(10):
    timef += dtime
    dU,new_ISVs,is_converged,n_iter = soln.solve_increment(dU,time=timef)
    #
    if is_converged:
        print("*************** TIME %.2f converged\n\n"%timef)
        soln.statev = new_ISVs.copy()
        soln.U += dU
        soln.solved_time = timef
        soln.time_steps += [soln.solved_time]
        
        soln.archive_results(update_ISVs=True)
        
        
# write to output:
#FE_io.write_to_vtu(soln,'patch_lagrange_strain')
FE_io.write_to_vtu(soln,filename0+'_nlgeom')
        
        

#%%

import FE_bvp
reload(FE_bvp)
bvp_dict =  FE_io.bvp_dict(filename)
soln = FE_bvp.solution(bvp_dict)
soln.solve_step()
FE_io.write_to_vtu(soln,filename0)




#%%

suse = soln2
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
#
#def assemble_FD(dU,time):
#    residual = soln.calc_residual(dU,time=time)
#    #
#    k_global_FD = []
#    for i in range(residual.size):
#        dU_FD = dU[soln.active_dofs]
#        dU_FD[i] += 1e-8
#        residualFD = soln.calc_residual(dU_FD,time=timef)
#        k_global_FD += [array(residualFD-residual).flatten()/1e-8]
#        
#    return residual, matrix(k_global_FD)
        
#%%
    
#
## prescribed displacements and contribution to dU:
#U_i, U_val = soln.disp_ext(time,from_step)
#if U_i.size>0:
#    dU[U_i] = U_val-soln.U[U_i]
##
##free degrees of freedom
#fdof = [i for i in soln.active_dofs if i not in U_i]
##
#
#Ra,Ka = assemble_FD(dU,timef)    
#
#residual[fdof] = Ra
#
## add external foces to resdual:
#F_i, F_val = soln.force_ext(time,from_step)
#if F_i.size>0:
#    residual[F_i,0] = F_val





#%%
# try Newton with line search:
#
#from_step = 1
#do_sparse = True
#only_active = False
#
#from scipy import sparse
#from scipy.sparse import linalg
#
#
#stepTime = soln.steps[from_step]['*time'][1]
#time = timef
#dtime = time-soln.solved_time
##
## prescribed displacements and contribution to dU:
#U_i, U_val = soln.disp_ext(time,from_step)
#if U_i.size>0:
#    dU[U_i] = U_val-soln.U[U_i]
##
##free degrees of freedom
#fdof = [i for i in soln.active_dofs if i not in U_i]
##
#residual,k_global,new_ISVs = soln.assemble(dU,time,dtime,from_step,do_sparse,only_active)
##
## add external foces to resdual:
#F_i, F_val = soln.force_ext(time,from_step)
#if F_i.size>0:
#    residual[F_i,0] = F_val
##
## residual norm
## relative residual norm
#ResNrm = norm(residual[fdof,0]);
#print('Normalized residual at start of iteration  = %10.6e'%ResNrm)
## solve 
#Kaa = k_global[fdof,:].T[fdof,:].T;
#Pa  = residual[fdof];
#Kaa_LU = sparse.linalg.splu(Kaa.tocsc())
#deltaUf = -array(Kaa_LU.solve(Pa)).flatten()
#dUNrm = norm(deltaUf)
#print('Normalized update to Unknowns              = %10.6e'%dUNrm)
##
#def linefn(alpha):
#    dU_new = dU.copy()
#    dU_new[fdof] += alpha*deltaUf
#    residual_new = array(soln.calc_residual(dU_new,time,dtime,from_step,do_sparse,only_active)).flatten()
#    obj = sum(deltaUf*residual_new[fdof])
#    return obj
##lower
#a_L = 0    
#s_0 = sum(deltaUf*array(residual[fdof]).flatten())
#s_L = s_0
##upper
#a_U = 1
#s_N = linefn(a_U)
#s_U = s_N
#kount = 0
#while abs(s_N/s_0)>1e-15 and kount<20:
#    #
#    kount += 1
#    # next guess
#    a_N = a_U - s_U*(a_L - a_U)/(s_L - s_U)
#    # test and update bounds:
#    s_N = linefn(a_N)
#    if s_N*s_L<0:
#        a_U = a_N
#        s_U = s_N
#    else:
#        a_L = a_N
#        s_L = s_N
#        
## new guess for the update
#dU[fdof] += a_N*deltaUf

# full newton:
#dU[fdof] += deltaUf




#%%
#def assemble_system(bvp,U0,dU):
#    
#    
#    
#    return bvp
#
#
#
#
##outputall = True
#bvp = FE_solver.solve(filename,returnbvp=True,outputall=True)


#FE_solver.solve(filename,resultsfile='results/'+filename0,outputall=True)


### Create Boundary Value Problems (dictionar(y)/(ies))
#tic = time()

#bvp = FE_io.bvp_dict(filename)
##
##
#print('\n** BOUNDARY VALUE DICTIONARY CREATED FROM INPUT FILE << %s >> '%(filename.split('/')[-1].split('.')[0]))
#print('\t\t>> done in %.4s seconds'%(time()-tic))
##
##
##
##
## SOLVE step 1:

#from scipy import sparse
#from scipy.sparse import linalg
#import pickle


##
##
## CRATE A LIST TO LINK dof AND NODES:
##
## number of dofs per node:
#n_dofs = 2
## conversion from node_id to dof: (node_id's start at 1)
#def t_dof(node_id):
    #return 0#n_dofs*(node_id-1)
#def x_dof(node_id):
    #return n_dofs*(node_id-1)
#def y_dof(node_id):
    #return n_dofs*(node_id-1)+1
##
##
##





## list of nodes and elements
#nd_lst = list(bvp.nodes.keys())
#el_lst = list(bvp.elements.keys())
##
#n_nds = nd_lst.__len__()


##solve for step 1
#nstep = 1
##
## current step:
#cstep = bvp.steps[nstep]
#SolvedStep = 0 # start of the step

## maximum number of increments

#inc_count = 0
#try:
    #maxinc = cstep['*inc']
#except:
    #maxinc = 1000
    
#dTime0 = cstep['*time'][0]
#Timescale = 1
#stepTime = cstep['*time'][1]

#time_format = bvp.time_format

#while (inc_count<maxinc)&(SolvedStep<stepTime):
    
    #inc_count += 1
    
    ### Time stamp (also used for internal state variables)
    #PrevTime = bvp.solved_time
    #time_prevs = time_format%PrevTime
    ## Time increment:
    #dTime = dTime0*Timescale
    #if SolvedStep + dTime>stepTime:
        #print('*** Time increment reduced to reach total step time')
        #dTime = stepTime-SolvedStep        
    #sTime = SolvedStep + dTime
    #Time = bvp.solved_time+dTime
    #time_stamp = time_format%Time
    #print('** Attempting to solve for time at %s seconds'%time_stamp)
    ##
    ##
    ##
    ##
    ## All nodal Unknowns
    #U_all = np.zeros(n_nds*n_dofs)
    #dU_all = np.zeros(n_nds*n_dofs)
    ##
    ## initial values (from bvp.NodeVals)
    #if 'temperature' in bvp.NodeVals[time_prevs].keys():
        #prevTemp = bvp.NodeVals[time_prevs]['temperature']
        #for nd_Id in prevTemp.keys():
            #U_all[t_dof(nd_Id)] = prevTemp[nd_Id]
    
    
    #if 'displacement' in bvp.NodeVals[time_prevs].keys():
        #prevDisp = bvp.NodeVals[time_prevs]['displacement']
        #for nd_Id in prevDisp.keys():
            #nd_U = prevDisp[nd_Id]
            ## check if a list of values (i.e. x and y):
            #if type(nd_U) is list:
                #U_all[x_dof(nd_Id)],U_all[y_dof(nd_Id)] = nd_U
                
            #elif type(nd_U) is dict:
                #if 1 in prevDisp[nd_Id].keys():
                    #U_all[x_dof(nd_Id)] = prevDisp[nd_Id][1]
                #if 2 in prevDisp[nd_Id].keys():
                    #U_all[y_dof(nd_Id)] = prevDisp[nd_Id][2]
    ##
##
    ##
    ## at this time step, what is the known boundary condition:
    ##
    ##
    #U_known = []
    #U_k_idx = []
    ##
    ## first apply potential initial "fixed" / "homogeneous" case dofs:
    #if '*boundary' in bvp.steps[0].keys():
        #for nd_Id in bvp.steps[0]['*boundary'].keys():
            #for dof_nr in bvp.steps[0]['*boundary'][nd_Id].keys():
                ## values fixed at 0
                
                ##if dof_nr == 11: # Temperature:
                    ##U_k_idx += [t_dof(nd_Id)]
                    ##U_known += [scalef*fullvalue]
                #if dof_nr == 1: # X-dof:
                    #U_k_idx += [x_dof(nd_Id)]
                    #U_known += [0.]
                #if dof_nr == 2: # Y-dof:
                    #U_k_idx += [y_dof(nd_Id)]
                    #U_known += [0.]
                    
    ## now for the particular step    
    #if '*boundary' in cstep.keys():
        #for nd_Id in cstep['*boundary'].keys():
            #for dof_nr in cstep['*boundary'][nd_Id].keys():
                #fullvalue = cstep['*boundary'][nd_Id][dof_nr]['value']
                #ampname = cstep['*boundary'][nd_Id][dof_nr]['amplitude']
                ##check if there is an amplitude to take into account, else ramp
                #scalef = sTime/stepTime
                
                ##if dof_nr == 11: # Temperature:
                    ##U_k_idx += [t_dof(nd_Id)]
                    ##U_known += [scalef*fullvalue]
                #if dof_nr == 1: # X-dof:
                    #U_k_idx += [x_dof(nd_Id)]
                    #U_known += [scalef*fullvalue]
                #if dof_nr == 2: # Y-dof:
                    #U_k_idx += [y_dof(nd_Id)]
                    #U_known += [scalef*fullvalue]
    ##
    ##
    ##
    ## EXTERNAL FORCES:
    #F_known = []
    #F_k_idx = []
    #if '*cload' in cstep.keys():
        #for nd_Id in cstep['*cload'].keys():
            #for dof_nr in cstep['*cload'][nd_Id].keys():
                #fullvalue = cstep['*cload'][nd_Id][dof_nr]['value']
                #ampname = cstep['*cload'][nd_Id][dof_nr]['amplitude']
                ##check if there is an amplitude to take into account, else ramp
                #scalef = sTime/stepTime
                
                #if dof_nr == 1: # X-dof:
                    #F_k_idx += [x_dof(nd_Id)]
                    #F_known += [scalef*fullvalue]
                #if dof_nr == 2: # Y-dof:
                    #F_k_idx += [y_dof(nd_Id)]
                    #F_known += [scalef*fullvalue]
    ##









    ##
    ## initialise dictionary for ISVs 
    #bvp.ISVs[time_format%Time] = {}
    ##
    ## initialise dictinary for nodal values:
    #bvp.NodeVals[time_format%Time] = {}
    ##
    ## for a specific section definition (loop over sections)
    ##section_name = 'eall'
    ##current_section = bvp.sections[section_name]
    ##current_material = bvp.materials[current_section['*material']]
    ##
    ##
    ##
    ##

    ## check only active dofs and free degrees of freedom:
    ##adof = [] # active degrees of freedom
    ##fdof = [] # free degrees of freedom
    ##for nd_Id in nd_lst:
        ##cdof = t_dof(nd_Id)
        ##if cdof in row_vec:
            ##adof+=[cdof]
            ##if cdof not in U_k_idx: fdof +=[cdof]

    #adof = [] # active degrees of freedom
    #fdof = [] # free degrees of freedom
    ##
    ##
    ##
    #print('\nSOLVE STATIC DISPLACEMENT')
    #U0 = np.matrix(U_all.reshape(-1,1))
    #dU = np.matrix(dU_all.reshape(-1,1))
    ## assign known values
    #dU[U_k_idx] = np.matrix(U_known).reshape(-1,1) - U0[U_k_idx]
    ##
    ## Residuals and counters
    #tol    = 1.e-6;
    #ResNrm,dUNrm = 1.,1.
    #ii   = -1
    #while (ResNrm>tol)|(dUNrm>tol):
        #tic = time()
        #ii += 1;
        ##Main loop over elements. Compute k_elem and assemble
        ##Initialize global stiffness matrix vectors;
        #row_vec   = []# np.zeros((TotKvec*nelem,1));
        #col_vec   = []#np.zeros((TotKvec*nelem,1);
        #stiff_vec = []#np.zeros(TotKvec*nelem,1);
        #Residual  = np.zeros((n_dofs*n_nds,1));
        ##Initialize global load vector
        #F_ext = np.zeros((n_dofs*n_nds,1));
        ##
        ##pos_vec = 0;
        ## loop over sections:
        #for section_name in bvp.sections.keys():
            #current_section = bvp.sections[section_name]
            #current_material = bvp.materials[current_section['*material']]
            ## loop over each element in the set   
            #for el_nr in bvp.el_set[section_name]:
                ##
                ## element copnnectivity:
                #el_connect = bvp.elements[el_nr]
                ## Find reference coordinates of element nodes
                #X = [bvp.nodes[nd_Id][0] for nd_Id in el_connect]
                #Y = [bvp.nodes[nd_Id][1] for nd_Id in el_connect]
                #XY   = np.c_[X, Y];
                ##
                ## element internal state variables
                #try:
                    #EL_ISVs = bvp.ISVs[time_prevs][el_nr]
                #except:
                    #EL_ISVs = {}
                ## global dof's associated with element:
                ## for temperature:
                ##pg_temp = [t_dof(nd_Id) for nd_Id in el_connect]
                ## for displacement and temperature
                #pg_disp = []
                #Temp = []
                #for nd_Id in el_connect:
                    #pg_disp += [x_dof(nd_Id),y_dof(nd_Id)]
                    #try:
                        #Temp += [bvp.NodeVals[time_prevs]['temperature'][nd_Id]]
                    #except:
                        #Temp += [0.]
                ## 
                ## GET CURRENT GUESS FOR NODAL DISPLACEMENTS:
                #U0_el = U0[pg_disp]
                #dU_el = dU[pg_disp]
                ##            
                ## displacement residual:
                #[U_residual,U_tangent,EN_ISVs] = FE_element.Element_Disp(XY,U0_el,dU_el,Temp,dTime,EL_ISVs,bvp._el_type[el_nr],current_section,current_material,cstep['*nlgeom'])
                ##
                #bvp.ISVs[time_stamp][el_nr] = EN_ISVs
                ##
                ## Assemble residual
                #Residual[pg_disp,0] += U_residual;
                ## Assemble T_tangent into sparse k_global using vectors
                #col_vec += pg_disp*pg_disp.__len__()
                #row_vec += list(np.array([pg_disp]*pg_disp.__len__()).T.flatten())
                #stiff_vec += list(U_tangent)
            
        ##
        ##
        
        ## Assemble k_global from vectors
        #k_global = sparse.csr_matrix((stiff_vec,(row_vec,col_vec)),shape=(n_dofs*n_nds,n_dofs*n_nds));
        ##
        ## clear memory
        ##row_vec,col_vec,stiff_vec = [],[],[]
        
        #finish = time()-tic;
        #print(['Done assembling stiffness matrix: %.4f seconds.'%finish])
        
        #tic = time()
        
        ### check only active dofs:
        ##adof = [] # active degrees of freedom
        ##fdof = [] # free degrees of freedom
        #if fdof.__len__()==0:
            #adof = [] # active degrees of freedom
            #fdof = [] # free degrees of freedom
            #for nd_Id in nd_lst:
                #cdof = x_dof(nd_Id)
                #if cdof in row_vec:
                    #adof+=[cdof]
                    #if cdof not in U_k_idx: fdof +=[cdof]
                #cdof = y_dof(nd_Id)
                #if cdof in row_vec:
                    #adof+=[cdof]
                    #if cdof not in U_k_idx: fdof +=[cdof]
            #finish = time()-tic;
            #print(['Sorted through active and inactive degrees of freedom : %.4f seconds.'%finish])
            #tic = time()
        
        
        ##
        ## Add nodal loads to global load vector
        #F_ext[F_k_idx] = np.matrix(F_known).reshape(-1,1)
        
        
        ## Subtract internal nodal loads
        #F = Residual - F_ext;
        
        ## relative residual norm
        #ResNrm = np.linalg.norm(F[fdof,0]);
        #if ii == 0:
            #if ResNrm > 1e-4:
                #ResNrm0 = ResNrm;
            #else:
                #ResNrm0 = 1;
        #ResNrm = ResNrm/ResNrm0;
        
        #print('Normalized residual at start of iteration %i   = %10.6e'%(ii+1,ResNrm))
        
        ### Solve update of free dof's
        ### Solution for non-symmetric stiffness matrix
        
        #Kaa = k_global[fdof,:].T[fdof,:].T;
        #Pa  = F[fdof];
        
        #finish = time()-tic;
        #print('Done assembling stiffness matrix: %.4f seconds.'%finish)
        
        #tic = time();
        
        ## factorise: (incomplete LU):
        #Kaa_LU = sparse.linalg.spilu(Kaa.tocsc())
        ##
        ## solve using native solve?
        #deltaUf = -Kaa_LU.solve(Pa)
        ## gmres or something else?
        
        #finish = time()-tic;
        #print('Done solving system             : %.4f seconds.'%finish)
        
        
        #dUNrm = np.linalg.norm(deltaUf);
        #if ii == 0:
            #if dUNrm > 1e-4:
                #dUNrm0 = dUNrm;
            #else:
                #dUNrm0 = 1;
        #dUNrm = dUNrm/dUNrm0;
        #print('Normalized update to Unknowns            = %10.6e'%dUNrm)
        #print('                    --------------------')
        ## Sort Uf and Ub into A
        ##U0[fdof,0] += deltaUf[:,0]#;
        
        #dU[fdof] += deltaUf
        
        
        ##AllResNrm[ii] = ResNrm;
        ##AlldUNrm[ii]  = dUNrm;
        
    ### Get support reactions
    ###Fp = F[U_k_idx];

    ### update ISVs:
    ##ISVs = new_ISVs
    ##new_ISVs = {}


    ##print('Load increment %i converged after %i iterations.'%(iter_load+1,ii+1))
    ###All_iter(iter_load) = ii;
    ###All_soln(1+n_nds*(iter_load-1):n_nds*iter_load,:) = [U0(1:2:n_dofs*n_nds) U0(2:2:n_dofs*n_nds)];

    ##
    ## time step:
    #SolvedStep += dTime
    #bvp.solved_time += dTime
    #bvp.time_steps += [bvp.solved_time]

    ### save solution (node values to BVP dictionary
    ##nd_vals = bvp.NodeVals[time_stamp]['temperature'] = {}
    ##for nd_Id in bvp.nd_set['000']:
        ##dof = t_dof(nd_Id)
        ##if dof in adof:
            ##nd_vals[nd_Id] = U0[dof][0,0]
            
    ## add displacement to bvp nodal values
    #nd_vals = bvp.NodeVals[time_stamp]['displacement'] = {}
    #for nd_Id in bvp.nd_set['000']:
        #dof = x_dof(nd_Id)
        #x_disp = U0[dof][0,0]
        #if dof in adof: x_disp += dU[dof][0,0]
        #dof = y_dof(nd_Id)
        #y_disp = U0[dof][0,0]
        #if dof in adof: y_disp += dU[dof][0,0]
        ##
        #nd_vals[nd_Id] = [x_disp,y_disp]


    #if outputall:
        #FE_io.write_to_vtu(bvp,'results/'+filename0)

##if outputall:
    ##FE_io.write_to_vtu(bvp,filename0,bvp.NodeVals.keys())

##else:
    
#if not outputall:
    #FE_io.write_to_vtu(bvp,'results/'+filename0)

