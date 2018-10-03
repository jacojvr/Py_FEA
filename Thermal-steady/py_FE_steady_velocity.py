import numpy as np
import sys
sys.path.append('../0')

import FE_io
import FE_element

from importlib import reload
reload(FE_io)


from time import time


filename = 'ccx_examples/steady_heat_2d_slab_structured'
#filename = 'ccx_examples/steady_heat_1d_test_unstructured'

velocity = 0.1 # mm/s    = 1m/s
#velocity = 100



### Create Boundary Value Problems (dictionar(y)/(ies))
tic = time()

bvp = FE_io.bvp_dict(filename)

# ADD A SECTION VELOCITY 
bvp.sections['eall']['*velocity'] = [-velocity,0] # material velocity in the undeformed configuration
#
#
print('\n** BOUNDARY VALUE DICTIONARY CREATED FROM INPUT FILE << %s >> '%(filename.split('/')[-1].split('.')[0]))
print('\t\t>> done in %.4s seconds'%(time()-tic))
#
#
#
#
# SOLVE step 1:

from scipy import sparse
from scipy.sparse import linalg
import pickle


#
#
#
# CRATE A LIST TO LINK dof AND NODES:
#
# number of dofs per node:
n_dofs = 1
# conversion from node_id to dof:
def t_dof(node_id):
    return n_dofs*(node_id-1)
#
#
#





# list of nodes and elements
nd_lst = list(bvp.nodes.keys())
el_lst = list(bvp.elements.keys())
#
n_nds = nd_lst.__len__()

## Time stamp (also used for internal state variables)
PrevTime = bvp.solved_time
time_format = bvp.time_format
time_prevs = time_format%PrevTime


#solve for step 1
nstep = 1
#
# current step:
cstep = bvp.steps[nstep]

dtime = cstep['*time'][0]


Time=bvp.solved_time+dtime
time_stamp = time_format%Time
#
#
# All nodal Unknowns
U_all = np.zeros(n_nds*n_dofs)
#
# initial values (from bvp.NodeVals)
if 'temperature' in bvp.NodeVals[time_prevs].keys():
    prevTemp = bvp.NodeVals[time_prevs]['temperature']
    for nd_Id in prevTemp.keys():
        U_all[t_dof(nd_Id)] = prevTemp[nd_Id]
#
#
#
# at this time step, what is the known boundary condition:
#
#
U_known = []
U_k_idx = []
if '*boundary' in cstep.keys():
    for nd_Id in cstep['*boundary'].keys():
        for dof_nr in cstep['*boundary'][nd_Id].keys():
            fullvalue = cstep['*boundary'][nd_Id][dof_nr]['value']
            ampname = cstep['*boundary'][nd_Id][dof_nr]['amplitude']
            #check if there is an amplitude to take into account, else ramp
            scalef = Time/cstep['*time'][1]
            
            if dof_nr == 11: # Temperature:
                U_k_idx += [t_dof(nd_Id)]
                U_known += [scalef*fullvalue]

#
#
# known / prescribed dofs (pdof) and free degrees of freedom (fof):
pdof = U_k_idx
#fdof = [t_dof(nd_Id) for nd_Id in nd_lst if not t_dof(nd_Id) in pdof]



# steady
dt = 0.
#
# initialise dictionary for ISVs 
bvp.ISVs[time_format%Time] = {}
#
# initialise dictinary for nodal values:
bvp.NodeVals[time_format%Time] = {}
#
# for a specific section definition (loop over sections)
section_name = 'eall'
current_section = bvp.sections[section_name]
current_material = bvp.materials[current_section['*material']]
#
#
#
#

# check only active dofs and free degrees of freedom:
#adof = [] # active degrees of freedom
#fdof = [] # free degrees of freedom
#for nd_idx in nd_lst:
    #cdof = t_dof(nd_idx)
    #if cdof in row_vec:
        #adof+=[cdof]
        #if cdof not in U_k_idx: fdof +=[cdof]

adof = [] # active degrees of freedom
fdof = [] # free degrees of freedom
#
#
#
print('\nSOLVE STEADY STATE HEAT TRANSFER')
U = np.matrix(U_all.reshape(-1,1))
# assign known values
U[U_k_idx,0] = U_known
#
# Residuals and counters
tol    = 1.e-12#6;
ResNrm,dUNrm = 1.,1.
ii   = -1
while (ResNrm>tol)|(dUNrm>tol):
    tic = time()
    ii += 1;
    #Main loop over elements. Compute k_elem and assemble
    #Initialize global stiffness matrix vectors;
    row_vec   = []# np.zeros((TotKvec*nelem,1));
    col_vec   = []#np.zeros((TotKvec*nelem,1);
    stiff_vec = []#np.zeros(TotKvec*nelem,1);
    Residual  = np.zeros((n_dofs*n_nds,1));
    #Initialize global load vector
    F_ext = np.zeros((n_dofs*n_nds,1));
    #
    #pos_vec = 0;
    for el_nr in el_lst:
        #
        # element copnnectivity:
        el_connect = bvp.elements[el_nr]
        # Find reference coordinates of element nodes
        X = [bvp.nodes[nd_idx][0] for nd_idx in el_connect]
        Y = [bvp.nodes[nd_idx][1] for nd_idx in el_connect]
        XY   = np.c_[X, Y];
        #
        # element internal state variables
        try:
            EL_ISVs = bvp.ISVs[time_prevs][el_nr]
        except:
            EL_ISVs = {}
        # global dof's associated with element:
        # for temperature:
        pg_temp = [t_dof(nd_idx) for nd_idx in el_connect]
        ## for displacement
        ## pg_disp = [dx_dof(nd_idx) for nd_idx in el_connect]
        # 
        # GET CURRENT GUESS FOR NODAL TEMPERATURES:
        T_el = U[pg_temp]
        # get element temperature residuals and tangent:
        [T_residual,T_tangent,EN_ISVs] = FE_element.Element_Temp(XY,0.,T_el,dt,EL_ISVs,bvp._el_type[el_nr],current_section,current_material)
        #
        bvp.ISVs[time_stamp][el_nr] = EN_ISVs
        #
        # Assemble residual
        Residual[pg_temp,0] += T_residual;
        # Assemble T_tangent into sparse k_global using vectors
        col_vec += pg_temp*pg_temp.__len__()
        row_vec += list(np.array([pg_temp]*pg_temp.__len__()).T.flatten())
        stiff_vec += list(T_tangent)
        
    #
    #
    
    # Assemble k_global from vectors
    k_global = sparse.csr_matrix((stiff_vec,(row_vec,col_vec)),shape=(n_dofs*n_nds,n_dofs*n_nds));
    #
    # clear memory
    #row_vec,col_vec,stiff_vec = [],[],[]
    
    finish = time()-tic;
    print(['Done assembling stiffness matrix: %.4f seconds.'%finish])
    
    tic = time()
    
    ## check only active dofs:
    #adof = [] # active degrees of freedom
    #fdof = [] # free degrees of freedom
    if fdof.__len__()==0:
        adof = [] # active degrees of freedom
        fdof = [] # free degrees of freedom
        for nd_idx in nd_lst:
            cdof = t_dof(nd_idx)
            if cdof in row_vec:
                adof+=[cdof]
                if cdof not in U_k_idx: fdof +=[cdof]
        finish = time()-tic;
        print(['Sorted through active and inactive degrees of freedom : %.4f seconds.'%finish])
        tic = time()
    
    
    #
    #ResNrm = 0
    #dUNrm = 0
    # Add nodal loads to global load vector
    #loadpos = np.array((cload[:,0]-1)*2+cload[:,1]-1,int)
    #F_ext[loadpos,0] = LoadFac[iter_load]*cload[:,2]
    
    
    # Subtract internal nodal loads
    F = Residual - F_ext;
    
    # relative residual norm
    ResNrm = np.linalg.norm(F[fdof,0]);
    if ii == 0:
        if ResNrm > 1e-4:
            ResNrm0 = ResNrm;
        else:
            ResNrm0 = 1;
    ResNrm = ResNrm/ResNrm0;
    
    print('Normalized residual at start of iteration %i   = %10.6e'%(ii+1,ResNrm))
    
    ## Solve update of free dof's
    ## Solution for non-symmetric stiffness matrix
    
    Kaa = k_global[fdof,:].T[fdof,:].T;
    Pa  = F[fdof];
    
    finish = time()-tic;
    print('Done assembling stiffness matrix: %.4f seconds.'%finish)
    
    tic = time();
    
    # factorise: (incomplete LU):
    Kaa_LU = sparse.linalg.spilu(Kaa.tocsc())
    #
    # solve using native solve?
    deltaUf = -Kaa_LU.solve(Pa)
    # gmres or something else?
    
    finish = time()-tic;
    print('Done solving system             : %.4f seconds.'%finish)
    
    
    dUNrm = np.linalg.norm(deltaUf);
    if ii == 0:
        if dUNrm > 1e-4:
            dUNrm0 = dUNrm;
        else:
            dUNrm0 = 1;
    dUNrm = dUNrm/dUNrm0;
    print('Normalized update to Unknowns            = %10.6e'%dUNrm)
    print('                    --------------------')
    # Sort Uf and Ub into A
    #U[fdof,0] += deltaUf[:,0]#;
    
    U[fdof] += deltaUf
    
    
    #AllResNrm[ii] = ResNrm;
    #AlldUNrm[ii]  = dUNrm;
    
## Get support reactions
##Fp = F[pdof];

## update ISVs:
#ISVs = new_ISVs
#new_ISVs = {}


#print('Load increment %i converged after %i iterations.'%(iter_load+1,ii+1))
##All_iter(iter_load) = ii;
##All_soln(1+n_nds*(iter_load-1):n_nds*iter_load,:) = [U(1:2:n_dofs*n_nds) U(2:2:n_dofs*n_nds)];

#
# time step:
bvp.time_steps += [1.]

# save solution (node values to BVP dictionary
nd_vals = bvp.NodeVals[time_stamp]['temperature'] = {}
for n_id in bvp.nd_set['000']:
    dof = t_dof(n_id)
    if dof in adof:
        nd_vals[n_id] = U[dof][0,0]

    


FE_io.write_to_vtu(bvp,'steady_vel_%.4f'%velocity)



#tic = time();


#mesh = {'nodal_coord':coor,
        #'connectivity':elnodes[:,1:]-1,
        #'dimension':2}
##
## ISVs extrapolated to nodes = 
#ISVnodes = isv_to_nodes(mesh,ISVs)
##
##
#reslt = {'mesh':mesh,
         #'displ':U.reshape(-1,2),
         #'ISVs':ISVs,
         #'ISVnodes':ISVnodes}

##
## dump results to binary file:
#pickle.dump(reslt,open(filename+'.FEAresult','wb'))


##ISVnodes = isv_to_nodes(mesh,ISVs)
##reslt['ISVnodes'] = ISVnodes
#write_to_vtu(reslt,filename)


## write a vtk XML file with displacements, stresses and strains and other internal state variables:
##write_to_vtu(reslt,filename)

## Compute Von Mises and Tresca and write output to text based output file
##[StressNode,VonMises,Tresca] = write_output_file(file_out,U,displ,ndispl,Fp, ...
    ##plane,pois,elas,n_nds,nodes,nelem,elnodes,stress,strain);



