



def solve(bvp,resultsfile=None,outputall=False,returnbvp=False):
    
    return solve_displ(bvp,resultsfile,outputall,returnbvp)





def initialise_solution(bvp,from_step=1):
    
    # use bvp.steps[0] and bvp.steps[from_step]
    step0 = 
    
    cstep = bvp.steps[from_step]
        
    # number of dofs per node:
    n_dofs = 2
    #
    #
    bvp.t_dof = lambda node_id: 0
    bvp.x_dof = lambda node_id: n_dofs*(node_id-1)
    bvp.y_dof = lambda node_id: n_dofs*(node_id-1)+1
    
    return bvp


def get_active_dof(bvp,from_step=1):
    
    
    
    return bvp



#def assemble_system(bvp,U0,dU,):
    
    
    
    #return bvp








def solve_displ(bvp,resultsfile=None,outputall=False,returnbvp=False,solver='cholesky'):
    
    
    import FE_io
    import FE_element
    from time import time
    #
    
    calc_dR_u = True
    if solver =='aitken':
        calc_dR_u = False
    else:
        try:
            from sksparse.cholmod import cholesky
            use_cholesky = True
        except:
            if not solver=='LU':
                solver = 'default'
        
    #
    aitken_relax = 1.e-5
        
    #except:
        #use_cholesky = False
        #print("*"*20+"\n** Cholesky factorisation unavailable")
    from scipy import sparse
    from scipy.sparse import linalg 
    
    from numpy import array, c_, matrix, linalg, zeros
    #import pickle
    
    # bvp input can be either a ccx file name or a bvp dictionary:
    if type(bvp) is str:
        filename = bvp[:]
        tic = time()
        bvp = FE_io.bvp_dict(bvp)
        print('\n** BOUNDARY VALUE DICTIONARY CREATED FROM INPUT FILE << %s >> '%(filename.split('/')[-1].split('.')[0]))
        print('\t\t>> done in %.4s seconds'%(time()-tic))
        
    #
    if resultsfile is None:
        resultsfile = bvp._folder+bvp._name
    #
    # CRATE A Function TO LINK dof AND NODES:
    #
    # number of dofs per node:
    n_dofs = 2
    # conversion from node_id to dof: (node_id's start at 1)
    def t_dof(node_id):
        return 0#n_dofs*(node_id-1)
    def x_dof(node_id):
        return n_dofs*(node_id-1)
    def y_dof(node_id):
        return n_dofs*(node_id-1)+1
    #
    #
    #
    bvp.t_dof = lambda node_id: 0
    bvp.x_dof = lambda node_id: n_dofs*(node_id-1)
    bvp.y_dof = lambda node_id: n_dofs*(node_id-1)+1
    
    #
    #
    bvp.active_dofs = None
    bvp.free_dofs = None
    adof = [] # active degrees of freedom
    fdof = [] # free degrees of freedom



    # list of nodes and elements
    nd_lst = list(bvp.nodes.keys())
    el_lst = list(bvp.elements.keys())
    #
    n_nds = nd_lst.__len__()


    #solve for step 1
    nstep = 1
    #
    # current step:
    cstep = bvp.steps[nstep]
    SolvedStep = 0 # start of the step

    # maximum number of increments

    inc_count = 0
    try:
        maxinc = cstep['*inc']
    except:
        maxinc = 1000
        
    dTime0 = cstep['*time'][0]
    Timescale = 1
    stepTime = cstep['*time'][1]

    time_format = bvp.time_format
    
    U_all = zeros(n_nds*n_dofs)


    while (inc_count<maxinc)&(SolvedStep<stepTime):
        
        inc_count += 1
        
        ## Time stamp (also used for internal state variables)
        PrevTime = bvp.solved_time
        time_prevs = time_format%PrevTime
        # Time increment:
        dTime = dTime0*Timescale
        if SolvedStep + dTime>stepTime:
            print('*** Time increment reduced to reach total step time')
            dTime = stepTime-SolvedStep        
        sTime = SolvedStep + dTime
        Time = bvp.solved_time+dTime
        time_stamp = time_format%Time
        print('** Attempting to solve for time at %s seconds'%time_stamp)
        #
        #
        #
        #
        # All nodal Unknowns
        dU_all = zeros(n_nds*n_dofs)
        #
        # initial values (from bvp.NodeVals)
        if 'temperature' in bvp.NodeVals[time_prevs].keys():
            prevTemp = bvp.NodeVals[time_prevs]['temperature']
            for nd_Id in prevTemp.keys():
                U_all[t_dof(nd_Id)] = prevTemp[nd_Id]
        
        
        if 'displacement' in bvp.NodeVals[time_prevs].keys():
            prevDisp = bvp.NodeVals[time_prevs]['displacement']
            for nd_Id in prevDisp.keys():
                nd_U = prevDisp[nd_Id]
                # check if a list of values (i.e. x and y):
                if type(nd_U) is list:
                    U_all[x_dof(nd_Id)],U_all[y_dof(nd_Id)] = nd_U
                    
                elif type(nd_U) is dict:
                    if 1 in nd_U.keys():
                        U_all[x_dof(nd_Id)] = nd_U[1]
                    if 2 in nd_U.keys():
                        U_all[y_dof(nd_Id)] = nd_U[2]
        #
    #
        #
        # at this time step, what is the known boundary condition:
        #
        #
        U_known = []
        U_k_idx = []
        #
        # first apply potential initial "fixed" / "homogeneous" case dofs:
        if '*boundary' in bvp.steps[0].keys():
            for nd_Id in bvp.steps[0]['*boundary'].keys():
                for dof_nr in bvp.steps[0]['*boundary'][nd_Id].keys():
                    # values fixed at 0
                    
                    #if dof_nr == 11: # Temperature:
                        #U_k_idx += [t_dof(nd_Id)]
                        #U_known += [scalef*fullvalue]
                    if dof_nr == 1: # X-dof:
                        U_k_idx += [x_dof(nd_Id)]
                        U_known += [0.]
                    if dof_nr == 2: # Y-dof:
                        U_k_idx += [y_dof(nd_Id)]
                        U_known += [0.]
                        
        # now for the particular step    
        if '*boundary' in cstep.keys():
            for nd_Id in cstep['*boundary'].keys():
                for dof_nr in cstep['*boundary'][nd_Id].keys():
                    fullvalue = cstep['*boundary'][nd_Id][dof_nr]['value']
                    ampname = cstep['*boundary'][nd_Id][dof_nr]['amplitude']
                    #check if there is an amplitude to take into account, else ramp
                    scalef = sTime/stepTime
                    
                    #if dof_nr == 11: # Temperature:
                        #U_k_idx += [t_dof(nd_Id)]
                        #U_known += [scalef*fullvalue]
                    if dof_nr == 1: # X-dof:
                        U_k_idx += [x_dof(nd_Id)]
                        U_known += [scalef*fullvalue]
                    if dof_nr == 2: # Y-dof:
                        U_k_idx += [y_dof(nd_Id)]
                        U_known += [scalef*fullvalue]
        #
        #
        #
        # EXTERNAL FORCES:
        F_known = []
        F_k_idx = []
        if '*cload' in cstep.keys():
            for nd_Id in cstep['*cload'].keys():
                for dof_nr in cstep['*cload'][nd_Id].keys():
                    fullvalue = cstep['*cload'][nd_Id][dof_nr]['value']
                    ampname = cstep['*cload'][nd_Id][dof_nr]['amplitude']
                    #check if there is an amplitude to take into account, else ramp
                    scalef = sTime/stepTime
                    
                    if dof_nr == 1: # X-dof:
                        F_k_idx += [x_dof(nd_Id)]
                        F_known += [scalef*fullvalue]
                    if dof_nr == 2: # Y-dof:
                        F_k_idx += [y_dof(nd_Id)]
                        F_known += [scalef*fullvalue]
        #









        #
        # initialise dictionary for ISVs 
        bvp.ISVs[time_format%Time] = {}
        #
        # initialise dictinary for nodal values:
        bvp.NodeVals[time_format%Time] = {}
        #
        # for a specific section definition (loop over sections)
        #section_name = 'eall'
        #current_section = bvp.sections[section_name]
        #current_material = bvp.materials[current_section['*material']]
        #
        #
        #
        #

        # check only active dofs and free degrees of freedom:
        #adof = [] # active degrees of freedom
        #fdof = [] # free degrees of freedom
        #for nd_Id in nd_lst:
            #cdof = t_dof(nd_Id)
            #if cdof in row_vec:
                #adof+=[cdof]
                #if cdof not in U_k_idx: fdof +=[cdof]

        doflst = []
        #
        #
        #
        print('\nSOLVE STATIC DISPLACEMENT')
        U0 = matrix(U_all.reshape(-1,1))
        dU = matrix(dU_all.reshape(-1,1))
        # assign known values
        dU[U_k_idx] = matrix(U_known).reshape(-1,1) - U0[U_k_idx]
        #
        # previous residuals
        F_k  = zeros((n_dofs*n_nds,1));
        F_km1  = zeros((n_dofs*n_nds,1));
        #
        # Residuals and counters
        tol    = 1.e-6;
        ResNrm,dUNrm = 1.,1.
        ii   = -1
        while (ResNrm>tol)|(dUNrm>tol):
            tic = time()
            ii += 1;
            #Main loop over elements. Compute k_elem and assemble
            #Initialize global stiffness matrix vectors;
            row_vec   = []# zeros((TotKvec*nelem,1));
            col_vec   = []#zeros((TotKvec*nelem,1);
            stiff_vec = []#zeros(TotKvec*nelem,1);
            Residual  = zeros((n_dofs*n_nds,1));
            #Initialize global load vector
            F_ext = zeros((n_dofs*n_nds,1));
            #
            #pos_vec = 0;
            # loop over sections:
            for section_name in bvp.sections.keys():
                current_section = bvp.sections[section_name]
                current_material = bvp.materials[current_section['*material']]
                # loop over each element in the set   
                for el_nr in bvp.el_set[section_name]:
                    #
                    # element copnnectivity:
                    el_connect = bvp.elements[el_nr]
                    # Find reference coordinates of element nodes
                    X = [bvp.nodes[nd_Id][0] for nd_Id in el_connect]
                    Y = [bvp.nodes[nd_Id][1] for nd_Id in el_connect]
                    XY   = c_[X, Y];
                    #
                    # element internal state variables
                    try:
                        EL_ISVs = bvp.ISVs[time_prevs][el_nr]
                    except:
                        EL_ISVs = {}
                    # global dof's associated with element:
                    # for temperature:
                    #pg_temp = [t_dof(nd_Id) for nd_Id in el_connect]
                    # for displacement and temperature
                    pg_disp = []
                    Temp = []
                    for nd_Id in el_connect:
                        pg_disp += [x_dof(nd_Id),y_dof(nd_Id)]
                        try:
                            Temp += [bvp.NodeVals[time_prevs]['temperature'][nd_Id]]
                        except:
                            Temp += [0.]
                    # 
                    # GET CURRENT GUESS FOR NODAL DISPLACEMENTS:
                    U0_el = U0[pg_disp]
                    dU_el = dU[pg_disp]
                    #            
                    # displacement residual:
                    [U_residual,U_tangent,EN_ISVs] = FE_element.Element_Disp(XY,U0_el,dU_el,Temp,dTime,EL_ISVs,bvp._el_type[el_nr],current_section,current_material,cstep['*nlgeom'],el_nr,calc_dR_u)
                    #
                    bvp.ISVs[time_stamp][el_nr] = EN_ISVs
                    #
                    #
                    #print(U_residual)
                    # Assemble residual
                    Residual[pg_disp,0] += U_residual;
                    #
                    #
                    doflst += pg_disp
                    if calc_dR_u:
                        # Assemble T_tangent into sparse k_global using vectors
                        col_vec += pg_disp*pg_disp.__len__()
                        row_vec += list(array([pg_disp]*pg_disp.__len__()).T.flatten())
                        stiff_vec += list(U_tangent)
                
            #
            #
            #
            if calc_dR_u:
                
                # Assemble k_global from vectors
                k_global = sparse.csr_matrix((stiff_vec,(row_vec,col_vec)),shape=(n_dofs*n_nds,n_dofs*n_nds));
                #
                # clear memory
                #row_vec,col_vec,stiff_vec = [],[],[]
                
                finish = time()-tic;
                print(['Done assembling stiffness matrix: %.4f seconds.'%finish])
                
            
            
            ## check only active dofs:
            
            tic = time()
            #adof = [] # active degrees of freedom
            #fdof = [] # free degrees of freedom
            if fdof.__len__()==0:
                adof = [] # active degrees of freedom
                fdof = [] # free degrees of freedom
                for nd_Id in nd_lst:
                    cdof = x_dof(nd_Id)
                    if cdof in doflst:
                        adof+=[cdof]
                        if cdof not in U_k_idx: fdof +=[cdof]
                    cdof = y_dof(nd_Id)
                    if cdof in doflst:
                        adof+=[cdof]
                        if cdof not in U_k_idx: fdof +=[cdof]
                finish = time()-tic;
                print(['Sorted through active and inactive degrees of freedom : %.4f seconds.'%finish])
                
                
                bvp.active_dofs = adof
                bvp.free_dofs = fdof
                tic = time()
            
            
            #
            # Add nodal loads to global load vector
            F_ext[F_k_idx] = matrix(F_known).reshape(-1,1)
            
            
            # Subtract internal nodal loads
            F = Residual - F_ext;
            
            # relative residual norm
            ResNrm = linalg.norm(F[fdof,0]);
            if ii == 0:
                if ResNrm > 1e-4:
                    ResNrm0 = ResNrm;
                else:
                    ResNrm0 = 1;
            ResNrm = ResNrm/ResNrm0;
            
            print('Normalized residual at start of iteration %i   = %10.6e'%(ii+1,ResNrm))
            
            if calc_dR_u:
                ## Solve update of free dof's
                ## Solution for non-symmetric stiffness matrix
                
                Kaa = k_global[fdof,:].T[fdof,:].T;
                Pa  = F[fdof];
                
                finish = time()-tic;
                print('Done assembling stiffness matrix: %.4f seconds.'%finish)
                
                tic = time();
                
                ##
                #
                if solver == 'LU':
                # factorise: (incomplete LU):
                    Kaa_LU = sparse.linalg.splu(Kaa.tocsc())
                    deltaUf = -Kaa_LU.solve(Pa)
                #
                # factorise (Cholesky) #scikits.sparse.cholmod
                elif solver=='cholesky':
                # # conda install -c conda-forge scikit-sparse
                #if use_cholesky:
                    CFx = cholesky(Kaa.tocsc())
                    deltaUf = -matrix(CFx(Pa))
                    
                else:
                    
                    ## solve using native solve?
                    deltaUf = -matrix(sparse.linalg.spsolve(Kaa,Pa)).T
                
                #print(deltaUf)
                
                #deltaUf = -matrix(sparse.linalg.gmres(Kaa,Pa)[0]).T
                #deltaUf = -matrix(sparse.linalg.bicgstab(Kaa,Pa)[0]).T
                
                
                finish = time()-tic;
                print('Done solving system             : %.4f seconds.'%finish)
                
            else:
                # do Aitkens method:
                F_k_km1 = F_k-F_km1
                a_denominator = array(F_k_km1.T*F_k_km1).flatten()[0]
                a_numerator = array(F_km1.T*F_k_km1).flatten()[0]
                if a_denominator==0:
                    aitken_sf = -1.
                else:
                    aitken_sf = -a_numerator/a_denominator
                #
                print(aitken_sf)
                aitken_relax = aitken_relax*aitken_sf
                #deltaUf = aitken_relax*F[fdof]
                deltaUf = F[fdof]
                #
                F_km1 = F_k.copy()
                F_k = F.copy()
        
            
            dUNrm = linalg.norm(deltaUf);
            if ii == 0:
                if dUNrm > 1e-4:
                    dUNrm0 = dUNrm;
                else:
                    dUNrm0 = 1;
            dUNrm = dUNrm/dUNrm0;
            print('Normalized update to Unknowns            = %10.6e'%dUNrm)
            print('                    --------------------')
            # Sort Uf and Ub into A
            #U0[fdof,0] += deltaUf[:,0]#;
            
            dUsf = 1.
            
            dU[fdof] += dUsf*deltaUf
            
            
            #AllResNrm[ii] = ResNrm;
            #AlldUNrm[ii]  = dUNrm;
            
        ## Get support reactions
        ##Fp = F[U_k_idx];

        ## update ISVs:
        #ISVs = new_ISVs
        #new_ISVs = {}


        #print('Load increment %i converged after %i iterations.'%(iter_load+1,ii+1))
        ##All_iter(iter_load) = ii;
        ##All_soln(1+n_nds*(iter_load-1):n_nds*iter_load,:) = [U0(1:2:n_dofs*n_nds) U0(2:2:n_dofs*n_nds)];

        #
        # time step:
        SolvedStep += dTime
        bvp.solved_time += dTime
        bvp.time_steps += [bvp.solved_time]
        
        print("Solved To : %.4f of %.4f"%(SolvedStep,stepTime))
        
        
        #U_all += dU

        ## save solution (node values to BVP dictionary
        #nd_vals = bvp.NodeVals[time_stamp]['temperature'] = {}
        #for nd_Id in bvp.nd_set['000']:
            #dof = t_dof(nd_Id)
            #if dof in adof:
                #nd_vals[nd_Id] = U0[dof][0,0]
                
        # add displacement to bvp nodal values
        nd_vals = bvp.NodeVals[time_stamp]['displacement'] = {}
        for nd_Id in bvp.nd_set['000']:
            dof = x_dof(nd_Id)
            x_disp = U0[dof][0,0]
            if dof in adof: x_disp += dU[dof][0,0]
            dof = y_dof(nd_Id)
            y_disp = U0[dof][0,0]
            if dof in adof: y_disp += dU[dof][0,0]
            #
            nd_vals[nd_Id] = [x_disp,y_disp]



        if outputall:
            FE_io.write_to_vtu(bvp,resultsfile)

    #if outputall:
        #FE_io.write_to_vtu(bvp,filename0,bvp.NodeVals.keys())

    #else:
        
    if not outputall:
        FE_io.write_to_vtu(bvp,resultsfile)
        
        
    if returnbvp:
        return bvp
    
    return