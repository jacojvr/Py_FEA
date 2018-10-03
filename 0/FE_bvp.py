#
#
#
#
#    HERE WE CREATE A BOUNDARY VALUE PROBLEM CLASS OBJECT FROM A DICTIONARY OR CCX FILE AND SOLVE IT
#
#
#

from numpy import array, c_, r_, matrix, linalg, zeros
#
#

class solution():
    
    def __init__(self,bvp,from_step=1):
        
        #
        # transfer relevant information to solution class:
        for attr,val in bvp.__dict__.items():
            setattr(self,attr,val)
            
        #
        self.step_info = self.steps[from_step]
        #
        #
        # how many nodes are defined / maximum node
        self.n_nds = max(self.nodes.keys())
        # cycle through all sections to see which nodes should have dof's defined:
        node_list = []
        for section_name,section_info in self.sections.items():
            for el_nr in self.el_set[section_name]:
                for node_id in self.elements[el_nr]:
                    if not node_id in node_list: node_list += [node_id]
        #
        #
        self.active_nodes = node_list
        self.active_nodes.sort()
        #
        # active degrees of freedom:
        adof = []
        
        
        
        if self.step_info['*type'] is '*temperature':
            
            self.n_dofs = 1
            self.t_dof = lambda node_id: node_id - 1
            for node_id in self.active_nodes:
                adof += self.t_dof(node_id)
            
        else:
            # only x and y displacement degrees of freedom will be used but temperatures also saved in nodal values [will remain the initial condition]
            self.n_dofs = 3
            self.x_dof = lambda node_id: self.n_dofs*(node_id-1)
            self.y_dof = lambda node_id: self.n_dofs*(node_id-1)+1
            self.t_dof = lambda node_id: self.n_dofs*(node_id-1)+2 
            
            for node_id in self.active_nodes:
                adof += [self.x_dof(node_id)]
                adof += [self.y_dof(node_id)]
                
        self.active_dofs = adof
        #
        # initialise nodal values:
        U_all = zeros(self.n_nds*self.n_dofs)
        #
        # define fixed values:
        time_prevs = self.time_format%0
        prev_soln = self.NodeVals[time_prevs].copy() #previous solution
        #
        # initial values (from bvp.NodeVals)
        if 'temperature' in prev_soln.keys():
            prevTemp = prev_soln['temperature']
            for node_id in prevTemp.keys():
                U_all[self.t_dof(node_id)] = prevTemp[node_id]
        
        
        if 'displacement' in prev_soln.keys():
            prevDisp = prev_soln['displacement']
            for node_id in prevDisp.keys():
                nd_U = prevDisp[node_id]
                # check if a list of values (i.e. x and y):
                if type(nd_U) is list:
                    U_all[x_dof(node_id)],U_all[y_dof(node_id)] = nd_U
                    
                elif type(nd_U) is dict:
                    if 1 in nd_U.keys():
                        U_all[x_dof(node_id)] = nd_U[1]
                    if 2 in nd_U.keys():
                        U_all[y_dof(node_id)] = nd_U[2]
        
        #
        #
        #
        # fixed BCs: (using step 0 )
        
        U_known = []
        U_k_idx = []
        #
        # first apply potential initial "fixed" / "homogeneous" case dofs:
        if '*boundary' in self.steps[0].keys():
            fixedBC = self.steps[0]['*boundary']
            for node_id in fixedBC.keys():
                for dof_nr in fixedBC[node_id].keys():
                    # values fixed at 0
                    
                    if dof_nr == 11: # Temperature:
                        U_k_idx += [self.t_dof(node_id)]
                        U_known += [fixedBC[node_id][dof_nr]['value']]
                    if dof_nr == 1: # X-dof:
                        U_k_idx += [self.x_dof(node_id)]
                        U_known += [0.]
                    if dof_nr == 2: # Y-dof:
                        U_k_idx += [self.y_dof(node_id)]
                        U_known += [0.]
                        
            # add fixed values to U_all:
            U_all[U_k_idx] = U_known
            
        #
        #
        self.U = U_all
        self.U_fixed = array(U_known)
        self.U_fixed_dofs = array(U_k_idx)
        
        #
        #
        #initial ISVs
        self.statev = self.ISVs[time_prevs].copy()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    def force_ext(self,time=None,from_step=1):
        '''
        calculate the external forces at given time using the forces defined in self.steps[from_step]
        
        '''
        
        stepTime = self.steps[from_step]['*time'][1]
        if time is None: time = stepTime
        
        F_known = []
        F_k_idx = []
        if '*cload' in self.steps[from_step].keys():
            BC_info = self.steps[from_step]['*cload']
            for nd_Id in BC_info.keys():
                for dof_nr in BC_info[nd_Id].keys():
                    fullvalue = BC_info[nd_Id][dof_nr]['value']
                    ampname = BC_info[nd_Id][dof_nr]['amplitude']
                    #check if there is an amplitude to take into account, else ramp
                    scalef = time/stepTime
                    
                    if dof_nr == 1: # X-dof:
                        F_k_idx += [self.x_dof(nd_Id)]
                        F_known += [scalef*fullvalue]
                    if dof_nr == 2: # Y-dof:
                        F_k_idx += [self.y_dof(nd_Id)]
                        F_known += [scalef*fullvalue]
                        
                        
        return array(F_k_idx), array(F_known)
    
    
    #
    #
    #
    #
    #  
    
    
    
    
    
    
    def disp_ext(self,time=None,from_step=1):
        '''
        calculate the prescribed total displacement
        '''
        
        stepTime = self.steps[from_step]['*time'][1]
        if time is None: time = stepTime
        
        
        U_known = []
        U_k_idx = []
        
        if '*boundary' in self.steps[from_step].keys():
            BC_info = self.steps[from_step]['*boundary']
            for node_id in BC_info.keys():
                for dof_nr in BC_info[node_id].keys():
                    fullvalue = BC_info[node_id][dof_nr]['value']
                    ampname = BC_info[node_id][dof_nr]['amplitude']
                    #
                    #check if there is an amplitude to take into account, else ramp
                    #
                    #
                    scalef = time/stepTime
                    
                    if dof_nr == 11: # Temperature:
                        U_k_idx += [self.t_dof(node_id)]
                        U_known += [scalef*fullvalue]
                    if dof_nr == 1: # X-dof:
                        U_k_idx += [self.x_dof(node_id)]
                        U_known += [scalef*fullvalue]
                    if dof_nr == 2: # Y-dof:
                        U_k_idx += [self.y_dof(node_id)]
                        U_known += [scalef*fullvalue]
        
        return array(U_k_idx), array(U_known)
    
    
    
    
    
    
    
    #
    #
    # determine / set-up surface dictionaries:
    def _surf_dict_setup(self):
        
        if not hasattr(self,'_surf_dict'):
            self._surf_dict = {}
            
        for sname in self.surfaces.keys():
            if not sname in self._surf_dict.keys():
                surf_dict = self._surf_dict[sname] = {}
                
                #
                # if it is a segment type surface, loop through to collect all the nodes and segments:
                if self.surfaces[sname]['*type'] == 'segments':
                    
                    #print("\n\nSURFACE : [%s]"%sname)
                    
                    for es in self.surfaces[sname]['*definition']:
                        
                        el_nr,seg_nr = [int(val) for val in es.split('--s')]
                        el_nds = self.elements[el_nr]
                        nnds = el_nds.__len__()
                        
                        #print("\t element ",el_nr," segmnt ",seg_nr," nds = ",el_nds)
                    
                        
                        if seg_nr == 1:
                            nd_idx = [el_nds[0],el_nds[1]]
                            if nnds == 8:
                                nd_idx += [el_nds[4]]
                        elif seg_nr == 2:
                            nd_idx = [el_nds[1],el_nds[2]]
                            if nnds == 8:
                                nd_idx += [el_nds[5]]
                        elif seg_nr == 3:
                            nd_idx = [el_nds[2],el_nds[3]]
                            if nnds == 8:

                                nd_idx += [el_nds[6]]
                        else:
                            nd_idx = [el_nds[3],el_nds[0]]
                            if nnds == 8:
                                nd_idx += [el_nds[7]]
                                
                        #print("\t use nodes = ",nd_idx)
                        #
                        # segments defined in terms of node numbers:
                        #
                        #
                        for nn in nd_idx:
                            try:
                                surf_dict[nn] += [nd_idx]
                            except:
                                surf_dict[nn] = [nd_idx]
                                
                    surf_dict['nds'] = [i for i in surf_dict.keys()]
                                
                else:
                    # definition is a set of nodal coordinates
                    surf_dict['nds'] =  self.surfaces[sname]['*definition']
                           
    
    
    #
    #
    #
    
    
    #
    #
    #
    def check_contact(self,dU=None,do_sparse=True,only_active=True,only_resid=False,contct='penalty'):
        '''
        
         check contact boundary conditions for a specific contact pair
         
        
        '''
        
        from time import time as tt
        from scipy import sparse
        from scipy.spatial import KDTree
        import FE_contact
        
        tic = tt()
        
        
        if dU is None: dU = zeros(self.n_nds*self.n_dofs)
        #
        if dU.size == self.active_dofs.__len__():
            dU0 = dU.copy()
            dU = zeros(self.n_nds*self.n_dofs)
            dU[self.active_dofs] = dU0
            
        #
        # initialise all lists that should be returned:
        # contact residuals / energy and corresponding global degrees of freedom 
        residuals = []
        resid_dof = []
        # residual sensitivities / tangent (as sparse array lists)
        stiff_val = []
        stiff_col = []
        stiff_row = []
        
        # first row index of the lagrange contact:
        if 'lagrange' in contct:
            nr_contact_elems = 0
            row_indx_0 = 0
            #
            #
            if '*pairs' in self.contact.keys():
                if not hasattr(self,'_surf_dict'):
                    self._surf_dict_setup()
                
                for cnt_pair in self.contact['*pairs']:
                    
                    name_slave,name_mast = cnt_pair.split('--')
                    name_intct = self.contact['*pairs'][cnt_pair]['*interaction']
                    
                    # surface interaction
                    if name_intct in self.contact['*surface interaction'].keys():
                        cnt_interact = self.contact['*surface interaction'][name_intct]
                    else:
                        cnt_interact = {'params': [1e6, 1e-3, None], 'pressure-overclosure': 'linear'}
                    
                    #
                    #
                    # NOTE : CURRENTLY ONLY NODE TO SURFACE CONTACT IMPLEMENTED
                    # for each node on the slave surface, calculate the 
                    slav_nds = self._surf_dict[name_slave]['nds']

                    # create a kdtree of surface master nodes
                    mast_nds = self._surf_dict[name_mast]['nds']
                    # master nodes degrees of freedom:
                    XY = []
                    mast_indx = {}
                    mast_dof = {}
                    new_indx = 0
                    for node_id in mast_nds:
                        
                        mast_indx[node_id] = new_indx
                        new_indx+=1
                        
                        m_dof = [self.x_dof(node_id),self.y_dof(node_id)]
                        mast_dof[node_id] = m_dof
                        
                        U_mn = self.U[m_dof]+dU[m_dof]
                        #print('master: ',node_id,' ', U_mn)
                        
                        XY += [self.nodes[node_id][:2]+U_mn]
                    XY = array(XY)
                    
                    mKDT = KDTree(XY)
                    
                    # loop over slave nodes
                    for snd in slav_nds:
                                            
                        # slave node degrees of freedom:
                        s_dof = [self.x_dof(snd),self.y_dof(snd)]
                        U_sn = self.U[s_dof]+dU[s_dof]
                        # current coordinate:
                        S_coord = array(self.nodes[snd][:2]) + U_sn #xy coordinate
                        
                        #print('\n\nslave: ',snd,' ', S_coord)
                        
                        # closest master node:
                        m_close = mKDT.query(S_coord)[1]
                        # segments to test:
                        segmnts = self._surf_dict[name_mast][mast_nds[m_close]]
                        #
                        # check the segments:
                        for m_seg in segmnts:
                            mxii = [mast_indx[ii] for ii in m_seg]
                            M_coords = array([XY[i] for i in mxii])
                            
                            #print("m_seg :\n",m_seg)
                            #print("mxii :\n",mxii)
                            
                            Lslav,Lmast = FE_contact.lagrange_contact_weights(S_coord,M_coords)
                            
                            #print('Lslav: \n',Lslav )
                            #print('Lmast: \n',Lmast )
                            
                            xr,yr = Lslav
                            
                            pg_dof = s_dof.copy()
                            for ii in m_seg:
                                pg_dof += mast_dof[ii]      
                                
                            
                            #print('PG_DOF : \n',pg_dof)
                            
                            if abs(xr)>1e-20:
                                #print(Lmast)
                                nr_contact_elems += 1
                                stiff_val += [xr]
                                stiff_val += list(Lmast[:,0])
                                
                                stiff_col += pg_dof[::2]
                                stiff_row += [row_indx_0]*pg_dof[::2].__len__()
                                row_indx_0+=1
                                
                                
                            if abs(yr)>1e-20:
                                #print(Lmast)
                                nr_contact_elems += 1
                                stiff_val += [yr]
                                stiff_val += list(Lmast[:,1])
                                
                                stiff_col += pg_dof[1::2]
                                stiff_row += [row_indx_0]*pg_dof[1::2].__len__()
                                row_indx_0+=1
            
            
            #print(stiff_val)
            #print(stiff_row)
            #print(stiff_col)
            
            k_global = sparse.csr_matrix((stiff_val,(stiff_row,stiff_col)),shape=(row_indx_0,self.n_dofs*self.n_nds))
            
            
            if only_active:
                k_global = k_global.T[self.active_dofs].T
            
            
            if not do_sparse:
                k_global = k_global.todense()
            
            
            finish = tt()-tic;
            if nr_contact_elems>0:
                print('Done adding %i contact elements in : %.4f seconds.'%(nr_contact_elems,finish))
            
            # return :
            return k_global
        
        
        
        
        
        else:
            
            # PENALTY CONTACT:
            
            if '*pairs' in self.contact.keys():
                if not hasattr(self,'_surf_dict'):
                    self._surf_dict_setup()
                
                for cnt_pair in self.contact['*pairs']:
                    
                    name_slave,name_mast = cnt_pair.split('--')
                    name_intct = self.contact['*pairs'][cnt_pair]['*interaction']
                    
                    # surface interaction
                    if name_intct in self.contact['*surface interaction'].keys():
                        cnt_interact = self.contact['*surface interaction'][name_intct]
                    else:
                        cnt_interact = {'params': [1e6, 1e-3, None], 'pressure-overclosure': 'linear'}
                    
                    #
                    #
                    # NOTE : CURRENTLY ONLY NODE TO SURFACE CONTACT IMPLEMENTED
                    # for each node on the slave surface, calculate the 
                    slav_nds = self._surf_dict[name_slave]['nds']

                    # create a kdtree of surface master nodes
                    mast_nds = self._surf_dict[name_mast]['nds']
                    # master nodes degrees of freedom:
                    XY = []
                    mast_indx = {}
                    mast_dof = {}
                    new_indx = 0
                    for node_id in mast_nds:
                        
                        mast_indx[node_id] = new_indx
                        new_indx+=1
                        
                        m_dof = [self.x_dof(node_id),self.y_dof(node_id)]
                        mast_dof[node_id] = m_dof
                        
                        U_mn = self.U[m_dof]+dU[m_dof]
                        
                        XY += [self.nodes[node_id][:2]+U_mn]
                    XY = array(XY)
                    
                    mKDT = KDTree(XY)
                    
                    # loop over slave nodes
                    for snd in slav_nds:
                                            
                        # slave node degrees of freedom:
                        s_dof = [self.x_dof(snd),self.y_dof(snd)]
                        U_sn = self.U[s_dof]+dU[s_dof]
                        # current coordinate:
                        S_coord = array(self.nodes[snd][:2]) + U_sn #xy coordinate
                        
                        # closest master node:
                        m_close = mKDT.query(S_coord)[1]
                        # segments to test:
                        segmnts = self._surf_dict[name_mast][mast_nds[m_close]]
                        #
                        # check the segments:
                        for m_seg in segmnts:
                            mxii = [mast_indx[ii] for ii in m_seg]
                            M_coords = array([XY[i] for i in mxii])
                            
                            Cinfo = FE_contact.penalty_residual(S_coord,M_coords,cnt_interact,only_resid)
                            
                            
                            cresid = Cinfo[0]
                            pg_dof = s_dof.copy()
                            for ii in m_seg:
                                pg_dof += mast_dof[ii]
                            
                            # add to sparse residual vector 
                            residuals += list(cresid)
                            resid_dof += pg_dof
                                
                            if not only_resid:
                                stiff_col += pg_dof*pg_dof.__len__()
                                stiff_row += list(array([pg_dof]*pg_dof.__len__()).T.flatten())
                                stiff_val += list(Cinfo[1].flatten())
                            
            
                        
            residual = sparse.csr_matrix((residuals,(resid_dof,[0]*residuals.__len__())),shape=(self.n_dofs*self.n_nds,1)).todense()
            
            if only_resid:
                return residual
            
            
            k_global = sparse.csr_matrix((stiff_val,(stiff_row,stiff_col)),shape=(self.n_dofs*self.n_nds,self.n_dofs*self.n_nds))
            
            if only_active:
                k_global = k_global[self.active_dofs].T[self.active_dofs].T
                residual = residual[self.active_dofs]
                    
            if not do_sparse:
                k_global = k_global.todense()
            
            finish = tt()-tic;
            print('Done evaluating contact stiffness : %.4f seconds.'%finish)
            
            
            
            # return :
            return residual,k_global
            
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    def assemble(self,dU=None,time=None,dtime=None,from_step=1,do_sparse=True,only_active=True):
        '''
        assemble the system of equations
        
        '''
        from scipy import sparse
        from time import time as tt
        import FE_element
        tic = tt()
        #
        #
        if dU is None: dU = zeros(self.n_nds*self.n_dofs)
        #
        if dU.size == self.active_dofs.__len__():
            dU0 = dU.copy()
            dU = zeros(self.n_nds*self.n_dofs)
            dU[self.active_dofs] = dU0
        #
        #
        
        stepTime = self.steps[from_step]['*time'][1]
        if time is None: time = stepTime
        if dtime is None: dtime = time-self.solved_time
        
        #
        step_type = self.steps[from_step]['*type']
        nlgeom = self.steps[from_step]['*nlgeom']
        
        # internal state variables:
        new_ISVs = {}
        #
        # initialise the sparse global stiffness matrix
        row_vec   = []
        col_vec   = []
        stiff_vec = []
        #
        # initialise the residual vector
        residual  = zeros((self.n_nds*self.n_dofs,1))
        #
        #
        for section_name,section_info in self.sections.items():
            #
            # section material information
            material_info = self.materials[self.sections[section_name]['*material']]
            #
            # loop over elements in the section
            for el_nr in self.el_set[section_name]:
                #
                # element type:
                el_type = self._el_type[el_nr]
                # element connectivity
                el_connect = self.elements[el_nr]
                #
                # element nodal coordinates
                X = [self.nodes[node_id][0] for node_id in el_connect]
                Y = [self.nodes[node_id][1] for node_id in el_connect]
                XY   = c_[X, Y];
                #
                # element internal state variables
                try:
                    EL_ISVs = self.statev[el_nr]
                except:
                    EL_ISVs = {}
                #
                if step_type == '*static':
                    #
                    pg_disp = []
                    Temp = []
                    for node_id in el_connect:
                        pg_disp += [self.x_dof(node_id),self.y_dof(node_id)]
                        Temp += [self.U[self.t_dof(node_id)]+dU[self.t_dof(node_id)]]
                        
                    # GET CURRENT GUESS FOR NODAL DISPLACEMENTS:
                    U0_el = self.U[pg_disp]
                    dU_el = dU[pg_disp]
                    #            
                    # displacement residual:
                    [U_residual,U_tangent,EN_ISVs] = FE_element.Element_Disp(XY,U0_el,dU_el,Temp,dtime,EL_ISVs,el_type,section_info,material_info,nlgeom,el_nr,only_resid=False)
                    #
                    new_ISVs[el_nr] = EN_ISVs
                    #
                    residual[pg_disp,0] += U_residual;
                    #
                    col_vec += pg_disp*pg_disp.__len__()
                    row_vec += list(array([pg_disp]*pg_disp.__len__()).T.flatten())
                    stiff_vec += list(U_tangent)
                    
                    
                        
                        
                elif step_type == '*steady heat':
                    # element temperature degrees of freedom:
                    pg_temp = [self.t_dof(node_id) for node_id in el_connect]
                
                else:
                    print('only << static >> and << steady heat transfer >> step types currently implemented ')
               
               
        #
        #
        # CALCULATE THE CONTACT RESIDUALS AND TANGENTS: (IF CONTACT)
        #
        # return contact residuals and stiffness in [col, row, val] for sparse matrix assembly
          
        #if '*pairs' in self.contact.keys():
            
            #C_R,C_dof,C_stiff,C_col,C_row = self.check_contact(dU,only_resid=False)
            
            #CRlen = C_R.__len__()
            #if CRlen>0:
                #Cresid = residual*0
                #for knr in range(CRlen):
                    #Cresid[C_dof[knr]] += C_R[knr]
                    
                #C_nrm = linalg.norm(Cresid)
                #U_nrm = linalg.norm(residual)
                
                
                
                #print('Contact residual norm %.4e vs Internal %.4e'%(C_nrm,U_nrm))
                #residual += Cresid
                    
                #col_vec += C_col
                #row_vec += C_row
                #stiff_vec += C_stiff
                    
                
                

               
               
        k_global = sparse.csr_matrix((stiff_vec,(row_vec,col_vec)),shape=(self.n_dofs*self.n_nds,self.n_dofs*self.n_nds))
        
        if only_active:
            k_global = k_global[self.active_dofs].T[self.active_dofs].T
            residual = residual[self.active_dofs]
                   
        if not do_sparse:
            k_global = k_global.todense()
         
        finish = tt()-tic;
        print('Done assembling stiffness matrix: %.4f seconds.'%finish)
        
        
        return residual,k_global,new_ISVs
    






    def calc_residual(self,dU=None,time=None,dtime=None,from_step=1,do_sparse=True,only_active=True):
        '''
        get the residual of the system of equations
        
        '''
        from time import time as tt
        import FE_element
        tic = tt()
        #
        #
        if dU is None: dU = zeros(self.n_nds*self.n_dofs)
        #
        if dU.size == self.active_dofs.__len__():
            dU0 = dU.copy()
            dU = zeros(self.n_nds*self.n_dofs)
            dU[self.active_dofs] = dU0
        #
        #
        
        stepTime = self.steps[from_step]['*time'][1]
        if time is None: time = stepTime
        if dtime is None: dtime = time-self.solved_time
        
        #
        step_type = self.steps[from_step]['*type']
        nlgeom = self.steps[from_step]['*nlgeom']
        #
        #
        # initialise the residual vector
        residual  = zeros((self.n_nds*self.n_dofs,1))
        #
        #
        for section_name,section_info in self.sections.items():
            #
            # section material information
            material_info = self.materials[self.sections[section_name]['*material']]
            #
            # loop over elements in the section
            for el_nr in self.el_set[section_name]:
                #
                # element type:
                el_type = self._el_type[el_nr]
                # element connectivity
                el_connect = self.elements[el_nr]
                #
                # element nodal coordinates
                X = [self.nodes[node_id][0] for node_id in el_connect]
                Y = [self.nodes[node_id][1] for node_id in el_connect]
                XY   = c_[X, Y];
                #
                # element internal state variables
                try:
                    EL_ISVs = self.statev[el_nr]
                except:
                    EL_ISVs = {}
                #
                if step_type == '*static':
                    #
                    pg_disp = []
                    Temp = []
                    for node_id in el_connect:
                        pg_disp += [self.x_dof(node_id),self.y_dof(node_id)]
                        Temp += [self.U[self.t_dof(node_id)]+dU[self.t_dof(node_id)]]
                        
                    # GET CURRENT GUESS FOR NODAL DISPLACEMENTS:
                    U0_el = self.U[pg_disp]
                    dU_el = dU[pg_disp]
                    #            
                    # displacement residual:
                    U_residual = FE_element.Element_Disp(XY,U0_el,dU_el,Temp,dtime,EL_ISVs,el_type,section_info,material_info,nlgeom,el_nr,only_resid=True)
                    #
                    residual[pg_disp,0] += U_residual
                        
                elif step_type == '*steady heat':
                    # element temperature degrees of freedom:
                    pg_temp = [self.t_dof(node_id) for node_id in el_connect]
                
                else:
                    print('only << static >> and << steady heat transfer >> step types currently implemented ')
        #
        #
        # CALCULATE THE CONTACT RESIDUALS (IF CONTACT)
        
        #if '*pairs' in self.contact.keys():
            #C_R,C_dof,C_stiff,C_col,C_row = self.check_contact(dU,only_resid=True)
            
            #CRlen = C_R.__len__()
            #if CRlen>0:
                #Cresid = residual*0
                #for knr in range(CRlen):
                    #Cresid[C_dof[knr]] += C_R[knr]
                    
                ## displacement stress residual:
                #residual += Cresid
                    
               
                
        if only_active:
            residual = residual[self.active_dofs]
            
        return residual






    #def calc_residual_pp(self,dU=None,time=None,dtime=None,from_step=1,do_sparse=True,only_active=True):
        #'''
        #get the residual of the system of equations
        
        #'''
        
        #from time import time as tt
        #import FE_element
        #import multiprocessing
        
        #pool = multiprocessing.Pool(4)
        
        #tic = tt()
        ##
        ##
        #if dU is None: dU = zeros(self.n_nds*self.n_dofs)
        ##
        #if dU.size == self.active_dofs.__len__():
            #dU0 = dU.copy()
            #dU = zeros(self.n_nds*self.n_dofs)
            #dU[self.active_dofs] = dU0
        ##
        ##
        
        #stepTime = self.steps[from_step]['*time'][1]
        #if time is None: time = stepTime
        #if dtime is None: dtime = time-self.solved_time
        
        ##
        #step_type = self.steps[from_step]['*type']
        #nlgeom = self.steps[from_step]['*nlgeom']
        ##
        ##
        ## initialise the residual vector
        #residual  = zeros((self.n_nds*self.n_dofs,1))
        ##
        #def el_resid(el_nr):
            ## element type:
            #el_type = self._el_type[el_nr]
            ## element connectivity
            #el_connect = self.elements[el_nr]
            ##
            ## element nodal coordinates
            #X = [self.nodes[node_id][0] for node_id in el_connect]
            #Y = [self.nodes[node_id][1] for node_id in el_connect]
            #XY   = c_[X, Y];
            ##
            ## element internal state variables
            #try:
                #EL_ISVs = self.statev[el_nr]
            #except:
                #EL_ISVs = {}
            ##
            #if step_type == '*static':
                ##
                #pg_disp = []
                #Temp = []
                #for node_id in el_connect:
                    #pg_disp += [self.x_dof(node_id),self.y_dof(node_id)]
                    #Temp += [self.U[self.t_dof(node_id)]+dU[self.t_dof(node_id)]]
                    
                ## GET CURRENT GUESS FOR NODAL DISPLACEMENTS:
                #U0_el = self.U[pg_disp]
                #dU_el = dU[pg_disp]
                ##            
                ## displacement residual:
                #U_residual = FE_element.Element_Disp(XY,U0_el,dU_el,Temp,dtime,EL_ISVs,el_type,section_info,material_info,nlgeom,el_nr,only_resid=True)
                ##
                 
            #elif step_type == '*steady heat':
                ## element temperature degrees of freedom:
                #pg_temp = [self.t_dof(node_id) for node_id in el_connect]
            
            #else:
                #print('only << static >> and << steady heat transfer >> step types currently implemented ')
            #return pg_disp, U_residual
                    
            
            
            
        
        ##
        #for section_name,section_info in self.sections.items():
            ##
            ## section material information
            #material_info = self.materials[self.sections[section_name]['*material']]
            ##
            
            ## loop over elements in the section
            #for el_nr in self.el_set[section_name]:
                #pg_disp, U_residual = el_resid(el_nr)
                #residual[pg_disp,0] += U_residual
                   
            #pgs,uresids = zip(*pool.map(el_resid,self.el_set[section_name]))
            
            #print(pgs)
            #print(uresids)
                
        #if only_active:
            #residual = residual[self.active_dofs]
            
        #return residual




    def solve_increment(self,dU=None,time=None,from_step=1,do_sparse=True,only_active=False):
        '''
        solve a particular time increment
        
        '''
        from scipy import sparse
        from scipy.sparse import linalg
        from numpy.linalg import norm
        from time import time as tt
        import FE_element
        tic = tt()
        #
        # modify / construct initial increment
        if dU is None: dU = zeros(self.n_nds*self.n_dofs)
        #
        if dU.size == self.active_dofs.__len__():
            dU0 = dU.copy()
            dU = zeros(self.n_nds*self.n_dofs)
            dU[self.active_dofs] = dU0
        #
        #
        
        stepTime = self.steps[from_step]['*time'][1]
        if time is None: time = stepTime
        dtime = time-self.solved_time
        #
        # get fixed degrees of freedom:
        U_val = self.U_fixed
        U_i = self.U_fixed_dofs
        if U_i.size>0:
            dU[U_i] = U_val-self.U[U_i]
        #free degrees of freedom
        fdof = [i for i in self.active_dofs if i not in U_i]
        # prescribed displacements and contribution to dU:
        U_i, U_val = self.disp_ext(time,from_step)
        if U_i.size>0:
            dU[U_i] = U_val-self.U[U_i]
        #update free degrees of freedom
        fdof = [i for i in fdof if i not in U_i]
        # 
        # external foces to resdual:
        F_i, F_val = self.force_ext(time,from_step)
        #
        #
        # Residuals and counters
        tol    = 1.e-10;
        ResNrm,dUNrm = 1.,1.
        ResNrm0,dUNrm0 =0,0
        ii   = -1
        maxcount = 50
        while (ResNrm>tol)|(dUNrm>tol):
            # update counter
            ii+=1
            if ii>maxcount:
                print('*** WARNNG: solution did not converge in %i attempts'%maxcount)
                return dU,new_ISVs,False,ii
            #
            #
            residual,k_global,new_ISVs = self.assemble(dU,time,dtime,from_step,do_sparse,only_active)
            # add external forces
            if F_i.size>0:
                residual[F_i,0] = F_val
                
                
            resid_contact = 0*residual
            # contact
            contct = ['augment lagrange','penalty'][0]
            w_sf = 1
            
            if '*pairs' in self.contact.keys():
                C_stiff = self.check_contact(dU,do_sparse=True,only_active=only_active,only_resid=False,contct=contct)
                
                AtA = None
                
                if contct == 'augment lagrange':
                    
                    AtA = C_stiff.T*C_stiff
                    
                
                elif contct == 'penalty':
                    if type(C_stiff) == tuple:
                        resid_contact = C_stiff[0].reshape(-1,1)
                        AtA = C_stiff[1]
                
                
                if not AtA is None:
                    #
                    # largest global stiffness nonzeros = 
                    K_max = max(abs(array(k_global[k_global.nonzero()]).flatten()))
                    #
                    # largest contact stiffness nonzero:
                    C_nz = AtA.nonzero()
                    C_vals = AtA[C_nz]
                    C_max = max(abs(array(C_vals).flatten()))
                    
                    #scale the contact contribution:
                    try:
                        w_sf = max(K_max/C_max,1)
                                        
                    except:
                        w_sf = 1
                        
                    # add w*AtA to the stiffness matrix;
                    k_global = k_global + w_sf*AtA
                
            
            
            # relative residual norm
            if ResNrm0==0:
                ResNrm0 = max(norm(residual[fdof,0]),1)
            ResNrm = norm(residual[fdof,0])/ResNrm0
            print('Normalized residual at start of iteration  %i = %10.6e'%(ii+1,ResNrm))
            # solve 
            #
            #
            Kaa = k_global[fdof,:].T[fdof,:].T;
            Pa  = residual[fdof];
            
            if contct == 'lagrange':
                
                row_0 = Kaa.shape[0]
                rows,cols = Kaa.nonzero()
                rows = list(rows)
                cols = list(cols)
                vals = list(array(Kaa[rows,cols]).flatten())
                #
                # add stiffness values:
                c_r,c_c = C_stiff.nonzero()
                if c_r.size>0:
                    #print(vals.shape,c_r.shape)
                    c_v = list(w_sf*array(C_stiff[c_r,c_c]).flatten())
                    c_r += row_0
                    c_r = list(c_r)
                    c_c = list(c_c)
                    new_rows = rows+c_r+c_c
                    new_cols = cols+c_c+c_r
                    new_vals = vals+c_v+c_v
                else:
                    new_rows = rows
                    new_cols = cols
                    new_vals = vals
                
                
                shp0 = max(new_rows)+1
                
                Kaa_new = sparse.csr_matrix((new_vals,(new_rows,new_cols)),shape=(shp0,shp0))
                Pa_new = zeros((shp0,1))
                Pa_new[:row_0] = Pa
                
                Kaa_LU = sparse.linalg.splu(Kaa_new.tocsc())
                deltaUf = -array(Kaa_LU.solve(Pa_new)).flatten()
                print("LAGRANGE multipliers : ",deltaUf[row_0:])
                deltaUf = deltaUf[:row_0]               
                
                
            elif contct == 'augment lagrange':
                
                Kaa_LU = sparse.linalg.splu(Kaa.tocsc())
                deltaUf = -array(Kaa_LU.solve(Pa)).flatten()
                
                
                dUL = dU.copy()
                dUL[fdof] += deltaUf
                #
                d_lagrange_mult = -w_sf*C_stiff*dUL
                lagrange_mult = d_lagrange_mult.copy()
                #
                dlprev = 0
                dlcurr = norm(lagrange_mult)
                kcount = 0
                while (abs(dlcurr-dlprev)>1e-3)&(kcount<5):
                    print("LAGRANGE MULTIPLIER NORM [%i] = [%.4f]"%(kcount+1,dlcurr))
                    kcount += 1
                    resid_contact = C_stiff.T*lagrange_mult
                    residual2 = residual - resid_contact.reshape(-1,1)
                    #print(residual2)
                    #print(residual.shape,residual2.shape)
                    deltaUf = -array(Kaa_LU.solve(residual2[fdof])).flatten()
                        
                    dUL = dU.copy()
                    dUL[fdof] += deltaUf
                    #
                    d_lagrange_mult = -w_sf*C_stiff*dUL
                    lagrange_mult += d_lagrange_mult
                    dlprev = dlcurr
                    dlcurr = norm(lagrange_mult)
                
            else:
                # Penalty contact
                Pa  = residual[fdof] - w_sf*resid_contact[fdof]
                Kaa_LU = sparse.linalg.splu(Kaa.tocsc())
                deltaUf = -array(Kaa_LU.solve(Pa)).flatten()
                
                
                
                
                
            
            # displacement update norm
            if dUNrm0==0:
                dUNrm0 = max(norm(deltaUf),1)
                
            dUNrm_PREV = dUNrm
            dUNrm = norm(deltaUf)/dUNrm0
            #
            #
            a_N = 1.
            if dUNrm>dUNrm_PREV:
                print('Normalized update to Unknowns before linesearch   = %10.6e'%dUNrm)
                print('Doing Newton-Raphson with bisection line search')
                #try:
                if 1<2:
                    # do a line search :
                    def linefn(alpha):
                        dU_new = dU.copy()
                        dU_new[fdof] += alpha*deltaUf
                        residual_new = array(self.calc_residual(dU_new,time,dtime,from_step,do_sparse,only_active)).flatten()
                        #
                        # contact
                        if contct == 'penalty':
                            if '*pairs' in self.contact.keys():
                                resid_contact2 = self.check_contact(dU_new,only_active=only_active,only_resid=True,contct=contct)
                                residual_new = residual_new - w_sf*array(resid_contact2).flatten()
                        elif contct == 'augment lagrange':
                           residual_new = residual_new - array(resid_contact).flatten()
                        #
                        #
                        obj = sum(deltaUf*residual_new[fdof])
                        return obj
                    #lower bound
                    a_L = 0.    
                    s_0 = sum(deltaUf*array(residual[fdof]).flatten())
                    s_L = s_0
                    #upper bound
                    a_U = 1.
                    s_N = linefn(a_U)
                    s_U = s_N
                    kount = 0
                    a_N = 1.
                    while abs(s_N/s_0)>1e-5 and kount<5:
                        #
                        kount += 1
                        # next guess
                        a_N = a_U - s_U*(a_L - a_U)/(s_L - s_U)
                        # test and update bounds:
                        s_N = linefn(a_N)
                        if s_N*s_L<0:
                            a_U = a_N
                            s_U = s_N
                        else:
                            a_L = a_N
                            s_L = s_N
                        
                    # if less than 5%, try a smaller time step size 
                    if a_N<0.05:
                        print('Full bisection failed - [ 0.05 ] scale factor update')
                        a_N = 0.05
                        
                        #return dU,new_ISVs,False,ii
                        
                            
                
                    print('                              linesearch factor   = %10.6e'%a_N)
                #except:
                    #print('Bisection failed - full update')
                    
            
            deltaUf *= a_N
            #
            #
            dUNrm = norm(deltaUf)/dUNrm0
            print('Normalized update to Unknowns                     = %10.6e'%dUNrm)
            dU[fdof] += deltaUf
            
            
            
        print('***********   Solution Converged at time step %.6f\n\n'%time)
        return dU,new_ISVs,True,ii
        #
        
        
        
        
    def archive_results(self,U=None,for_time=None,update_ISVs=False):
        
        if U is None: U = self.U
        if for_time is None: for_time = self.solved_time
        
        time_stamp = self.time_format%for_time
        
        if not time_stamp in self.NodeVals.keys():
            self.NodeVals[time_stamp] = {}
        
        nd_dict = self.NodeVals[time_stamp]['displacement'] = {}
        for node_id in self.nd_set['000']:
            try:
                x_disp = U[self.x_dof(node_id)]
            except:
                x_disp = 0.
            
            try:
                y_disp = U[self.y_dof(node_id)]
            except:
                y_disp = 0.
                
            nd_dict[node_id] = [x_disp,y_disp]
            
            
        
        if update_ISVs:
            self.ISVs[time_stamp] = self.statev.copy()
            
        
        
        
        
    def solve_step(self, from_step=1, archive_results=True):
        
        
        
        nlgeom = self.steps[from_step]['*nlgeom']
        maxinc = self.steps[from_step]['*inc']
        dtime,steptime,mintime,maxtime = self.steps[from_step]['*time']
        
        timef = 0.
        dU0 = self.U*0
        dUsf = dtime
        
        time_up = 2.
        time_down = 0.5
        
        while timef<steptime:
            cdtime = dtime
            trytime = timef + dtime
            if trytime > steptime:
                trytime = steptime    
                dUsf *= (trytime-timef) / dtime
                
            dU,new_ISVs,is_converged,n_iter = self.solve_increment(dU0*dUsf,time=trytime)
            
            if is_converged:
                timef = trytime
                dU0 = dU
                dUsf = 1
                
                print("\n            >>> TIME %.6f converged\n"%timef)
                self.statev = new_ISVs.copy()
                self.U += dU
                self.solved_time = timef
                self.time_steps += [self.solved_time]
                
                if archive_results:
                    self.archive_results(update_ISVs=True)
            
            
                if n_iter<10:
                    print("\n            >>> TIME %.6f < 10 iterations, INCREASING automatic time step size \n"%trytime)
                    dtime *= time_up
                    dUsf *= time_up
                    
                
            if n_iter>20:
                print("\n            >>> TIME %.6f > 20 iterations, DECREASING automatic time step size \n"%trytime)
                dtime *= time_down
                dUsf *= time_down
                
            if not is_converged:
                print("\n            >>> TIME %.6f did not converge, DECREASING automatic time step size \n"%trytime)
                dtime *= time_down
                dUsf *= time_down
                
            if dtime<mintime:
                dtime = mintime
                
            if dtime>maxtime:
                dtime = maxtime
                
                
            if not cdtime == dtime:
                print("\n               >>> NEXT TRYING TIME INCREMENT  =  %.6f\n"%dtime)
            
            
            
            
            
            
    
    
    
    
                        
                        

                    
                
                
                
        
    
    #def output(self):
        #
        # output results to *.vtk and *.out files 
        #
        #
        
        
        
        
#    def solve_step(self,dU=None,from_step=1):
#        '''
#        solve a particular step
#        
#        '''
#        #
#        # modify / construct initial increment
#        if dU is None: dU = zeros(self.n_nds*self.n_dofs)
#        #
#        if dU.size == self.active_dofs.__len__():
#            dU0 = dU.copy()
#            dU = zeros(self.n_nds*self.n_dofs)
#            dU[self.active_dofs] = dU0
#        #
#        #
#        
#        stepTime = self.steps[from_step]['*time'][1]
#        dTime = self.steps[from_step]['*time'][0]
#        if time is None: time = stepTime
#        dtime = time-self.solved_time
        
        






#def solve_displ(bvp,resultsfile=None,outputall=False,returnbvp=False,solver='cholesky'):
    
    
    #import FE_io
    #import FE_element
    #from time import time
    ##
    
    #calc_dR_u = True
    #if solver =='aitken':
        #calc_dR_u = False
    #else:
        #try:
            #from sksparse.cholmod import cholesky
            #use_cholesky = True
        #except:
            #if not solver=='LU':
                #solver = 'default'
        
    ##
    #aitken_relax = 1.e-5
        
    ##except:
        ##use_cholesky = False
        ##print("*"*20+"\n** Cholesky factorisation unavailable")
    #from scipy import sparse
    #from scipy.sparse import linalg 
    
    #from numpy import array, c_, matrix, linalg, zeros
    ##import pickle
    
    ## bvp input can be either a ccx file name or a bvp dictionary:
    #if type(bvp) is str:
        #filename = bvp[:]
        #tic = time()
        #bvp = FE_io.bvp_dict(bvp)
        #print('\n** BOUNDARY VALUE DICTIONARY CREATED FROM INPUT FILE << %s >> '%(filename.split('/')[-1].split('.')[0]))
        #print('\t\t>> done in %.4s seconds'%(time()-tic))
        
    ##
    #if resultsfile is None:
        #resultsfile = bvp._folder+bvp._name
    ##
    ## CRATE A Function TO LINK dof AND NODES:
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
    #bvp.t_dof = lambda node_id: 0
    #bvp.x_dof = lambda node_id: n_dofs*(node_id-1)
    #bvp.y_dof = lambda node_id: n_dofs*(node_id-1)+1
    
    ##
    ##
    #bvp.active_dofs = None
    #bvp.free_dofs = None
    #adof = [] # active degrees of freedom
    #fdof = [] # free degrees of freedom



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
    
    #U_all = zeros(n_nds*n_dofs)


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
        #dU_all = zeros(n_nds*n_dofs)
        ##
        ## initial values (from bvp.NodeVals)
        #if 'temperature' in bvp.NodeVals[time_prevs].keys():
            #prevTemp = bvp.NodeVals[time_prevs]['temperature']
            #for node_id in prevTemp.keys():
                #U_all[t_dof(node_id)] = prevTemp[node_id]
        
        
        #if 'displacement' in bvp.NodeVals[time_prevs].keys():
            #prevDisp = bvp.NodeVals[time_prevs]['displacement']
            #for node_id in prevDisp.keys():
                #nd_U = prevDisp[node_id]
                ## check if a list of values (i.e. x and y):
                #if type(nd_U) is list:
                    #U_all[self.x_dof(node_id)],U_all[self.y_dof(node_id)] = nd_U
                    
                #elif type(nd_U) is dict:
                    #if 1 in nd_U.keys():
                        #U_all[x_dof(node_id)] = nd_U[1]
                    #if 2 in nd_U.keys():
                        #U_all[y_dof(node_id)] = nd_U[2]
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

        #doflst = []
        ##
        ##
        ##
        #print('\nSOLVE STATIC DISPLACEMENT')
        #U0 = matrix(U_all.reshape(-1,1))
        #dU = matrix(dU_all.reshape(-1,1))
        ## assign known values
        #dU[U_k_idx] = matrix(U_known).reshape(-1,1) - U0[U_k_idx]
        ##
        ## previous residuals
        #F_k  = zeros((n_dofs*n_nds,1));
        #F_km1  = zeros((n_dofs*n_nds,1));
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
            #row_vec   = []# zeros((TotKvec*nelem,1));
            #col_vec   = []#zeros((TotKvec*nelem,1);
            #stiff_vec = []#zeros(TotKvec*nelem,1);
            #Residual  = zeros((n_dofs*n_nds,1));
            ##Initialize global load vector
            #F_ext = zeros((n_dofs*n_nds,1));
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
                    #XY   = c_[X, Y];
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
                    #[U_residual,U_tangent,EN_ISVs] = FE_element.Element_Disp(XY,U0_el,dU_el,Temp,dTime,EL_ISVs,bvp._el_type[el_nr],current_section,current_material,cstep['*nlgeom'],el_nr,calc_dR_u)
                    ##
                    #bvp.ISVs[time_stamp][el_nr] = EN_ISVs
                    ##
                    ##
                    ##print(U_residual)
                    ## Assemble residual
                    #Residual[pg_disp,0] += U_residual;
                    ##
                    ##
                    #doflst += pg_disp
                    #if calc_dR_u:
                        ## Assemble T_tangent into sparse k_global using vectors
                        #col_vec += pg_disp*pg_disp.__len__()
                        #row_vec += list(array([pg_disp]*pg_disp.__len__()).T.flatten())
                        #stiff_vec += list(U_tangent)
                
            ##
            ##
            ##
            #if calc_dR_u:
                
                ## Assemble k_global from vectors
                #k_global = sparse.csr_matrix((stiff_vec,(row_vec,col_vec)),shape=(n_dofs*n_nds,n_dofs*n_nds));
                ##
                ## clear memory
                ##row_vec,col_vec,stiff_vec = [],[],[]
                
                #finish = time()-tic;
                #print(['Done assembling stiffness matrix: %.4f seconds.'%finish])
                
            
            
            ### check only active dofs:
            
            #tic = time()
            ##adof = [] # active degrees of freedom
            ##fdof = [] # free degrees of freedom
            #if fdof.__len__()==0:
                #adof = [] # active degrees of freedom
                #fdof = [] # free degrees of freedom
                #for nd_Id in nd_lst:
                    #cdof = x_dof(nd_Id)
                    #if cdof in doflst:
                        #adof+=[cdof]
                        #if cdof not in U_k_idx: fdof +=[cdof]
                    #cdof = y_dof(nd_Id)
                    #if cdof in doflst:
                        #adof+=[cdof]
                        #if cdof not in U_k_idx: fdof +=[cdof]
                #finish = time()-tic;
                #print(['Sorted through active and inactive degrees of freedom : %.4f seconds.'%finish])
                
                
                #bvp.active_dofs = adof
                #bvp.free_dofs = fdof
                #tic = time()
            
            
            ##
            ## Add nodal loads to global load vector
            #F_ext[F_k_idx] = matrix(F_known).reshape(-1,1)
            
            
            ## Subtract internal nodal loads
            #F = Residual - F_ext;
            
            ## relative residual norm
            #ResNrm = linalg.norm(F[fdof,0]);
            #if ii == 0:
                #if ResNrm > 1e-4:
                    #ResNrm0 = ResNrm;
                #else:
                    #ResNrm0 = 1;
            #ResNrm = ResNrm/ResNrm0;
            
            #print('Normalized residual at start of iteration %i   = %10.6e'%(ii+1,ResNrm))
            
            #if calc_dR_u:
                ### Solve update of free dof's
                ### Solution for non-symmetric stiffness matrix
                
                #Kaa = k_global[fdof,:].T[fdof,:].T;
                #Pa  = F[fdof];
                
                #finish = time()-tic;
                #print('Done assembling stiffness matrix: %.4f seconds.'%finish)
                
                #tic = time();
                
                ###
                ##
                #if solver == 'LU':
                ## factorise: (incomplete LU):
                    #Kaa_LU = sparse.linalg.splu(Kaa.tocsc())
                    #deltaUf = -Kaa_LU.solve(Pa)
                ##
                ## factorise (Cholesky) #scikits.sparse.cholmod
                #elif solver=='cholesky':
                ## # conda install -c conda-forge scikit-sparse
                ##if use_cholesky:
                    #CFx = cholesky(Kaa.tocsc())
                    #deltaUf = -matrix(CFx(Pa))
                    
                #else:
                    
                    ### solve using native solve?
                    #deltaUf = -matrix(sparse.linalg.spsolve(Kaa,Pa)).T
                
                ##print(deltaUf)
                
                ##deltaUf = -matrix(sparse.linalg.gmres(Kaa,Pa)[0]).T
                ##deltaUf = -matrix(sparse.linalg.bicgstab(Kaa,Pa)[0]).T
                
                
                #finish = time()-tic;
                #print('Done solving system             : %.4f seconds.'%finish)
                
            #else:
                ## do Aitkens method:
                #F_k_km1 = F_k-F_km1
                #a_denominator = array(F_k_km1.T*F_k_km1).flatten()[0]
                #a_numerator = array(F_km1.T*F_k_km1).flatten()[0]
                #if a_denominator==0:
                    #aitken_sf = -1.
                #else:
                    #aitken_sf = -a_numerator/a_denominator
                ##
                #print(aitken_sf)
                #aitken_relax = aitken_relax*aitken_sf
                ##deltaUf = aitken_relax*F[fdof]
                #deltaUf = F[fdof]
                ##
                #F_km1 = F_k.copy()
                #F_k = F.copy()
        
            
            #dUNrm = linalg.norm(deltaUf);
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
            
            #dUsf = 1.
            
            #dU[fdof] += dUsf*deltaUf
            
            
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
        
        #print("Solved To : %.4f of %.4f"%(SolvedStep,stepTime))
        
        
        ##U_all += dU

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
            #FE_io.write_to_vtu(bvp,resultsfile)

    ##if outputall:
        ##FE_io.write_to_vtu(bvp,filename0,bvp.NodeVals.keys())

    ##else:
        
    #if not outputall:
        #FE_io.write_to_vtu(bvp,resultsfile)
        
        
    #if returnbvp:
        #return bvp
    
    #return