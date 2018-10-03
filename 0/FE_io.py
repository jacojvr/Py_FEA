import numpy as np
from time import time
import os



from setup_from_ccx_inp import bvp_dict
import FE_shapefn


def write_to_vtu(bvp,filename,time_steps=None,remove_old=False):
    
    # requires vtk support:
    import vtk
    
    if remove_old:
        import os
        os.system('rm '+filename+'*.vtk');
    #
    # write 
    
    #fid = open(filename+'.vtk','w')
    #
    # mesh info to vtk objects:
    #XY = reslt['mesh']['nodal_coord']
    #bvp.nodes
    #EL = reslt['mesh']['connectivity']
    
    #DIM = reslt['mesh']['dimension']
    
    # if no time stamp, or list of time stamps, do the last
    time_format = bvp.time_format
    if time_steps is None:
        try:
            time_steps = bvp.time_steps
        except:
            print("No time step information in BVP dictionary \n\t>> EXITING")
            return
        
    #
    
    #
    # add points:
    points = vtk.vtkPoints()
    for pi,px in bvp.nodes.items():
        points.InsertPoint(pi,tuple(px))
        
    # number of points:
    nnodes = list(bvp.nodes.keys()).__len__()
    #
    # add cells:
    cellArray = vtk.vtkCellArray()
    cellTypes = []
    
    new_el_nrs = {}
    el_i = 0
    for old_el in bvp.elements.keys():
        new_el_nrs[old_el] = el_i
        el_i += 1
        
        if bvp._el_type[old_el] in ['cpe4','cpe4r','cps4','cps4r']:
            cell = vtk.vtkQuad()
            cellTypes += [vtk.VTK_QUAD]
            for i in range(4):
                cell.GetPointIds().SetId(i,bvp.elements[old_el][i])
                
        elif bvp._el_type[old_el] in ['cpe8','cpe8r','cps8','cps48']:
            cell = vtk.vtkQuadraticQuad()
            cellTypes += [vtk.VTK_QUADRATIC_QUAD]
            for i in range(8):
                cell.GetPointIds().SetId(i,bvp.elements[old_el][i])
                
                
                
                
            
        cellArray.InsertNextCell(cell)
        
    
    #for x,y in XY:
        #points.InsertNextPoint(x,y,0)
        
    
    
    #if DIM == 2:
        ## create VTK point array:
        #points = vtk.vtkPoints()
        #for x,y in XY:
            #points.InsertNextPoint(x,y,0)
            
        ## create VTK cell array:
        #cellArray = vtk.vtkCellArray()
        #cellTypes = []
        #for elnr,ndlst in enumerate(EL):
            ## only 2d so number of nodes dictate element type:
            ## 2 = line segment # VTK_LINE (cell type = 3)
            #if ndlst.size == 2:
                #cell = vtk.vtkLine()
                #cellTypes += [vtk.VTK_LINE]
            ## 3 = triangle # VTK_TRIANGLE (cell type = 5)
            #if ndlst.size == 3:
                #cell = vtk.vtkTriangle()
                #cellTypes += [vtk.VTK_TRIANGLE]
            ## 4 = quad # VTK_QUAD (cell type = 9)
            #if ndlst.size == 4:
                #cell = vtk.vtkQuad()
                #cellTypes += [vtk.VTK_QUAD]
            ## 6 = quadratic triangle # VTK_QUADRATIC_TRIANGLE (cell type = 22)
            #if ndlst.size == 6:
                #cell = vtk.vtkQuadraticTriangle()
                #cellTypes += [vtk.VTK_QUADRATIC_TRIANGLE]
            ## 8 = quadratic quad # VTK_QUADRATIC_QUAD (cell type = 23)
            #if ndlst.size == 8:
                #cell = vtk.vtkQuadraticQuad()
                #cellTypes += [vtk.VTK_QUADRATIC_QUAD]
            ## 9 = bi-quadratic quad # VTK_BIQUADRATIC_QUAD (cell type = 28)
            #if ndlst.size == 9:
                #cell = vtk.vtkBiQuadraticQuad()
                #cellTypes += [vtk.VTK_QUADRATIC_QUAD]
                
            #for i,v in enumerate(ndlst):
                #cell.GetPointIds().SetId(i,v)
            ##
            #cellArray.InsertNextCell(cell)
        
    #
        
        
    #grid0 = vtk.vtkUnstructuredGrid()
    #grid0.SetPoints(points)
    #grid0.SetCells(cellTypes,cellArray)
    
    
    grid = vtk.vtkUnstructuredGrid()
    grid.SetPoints(points)
    grid.SetCells(cellTypes,cellArray)
    
    
    #vtuWriter = vtk.vtkXMLPUnstructuredGridWriter()
    #vtuWriter = vtk.vtkXMLUnstructuredGridWriter()
    
    #vtuWriter.SetInputDataObject(grid)
    #vtuWriter.SetFileName(filename+'.vtu')
    #vtuWriter.SetNumberOfTimeSteps(time_steps.__len__())
    #vtuWriter.SetDataModeToAscii()
    #vtuWriter.SetCompressorTypeToNone()
    
    #vtuWriter.Start()
    
    
    nr_of_files = time_steps.__len__()
    
            
    file_suffix = ''
    
    if nr_of_files>1:
        # first find the time scale factor:
        time_sf = 1
        inctsf = any(np.array(time_steps)*10**time_sf%1)
        while inctsf:
            times2 = np.array(time_steps)*10**time_sf//1
            if np.unique(times2).size == nr_of_files:
                inctsf = False
            else:
                time_sf += 1
                inctsf = any(np.array(time_steps)*10**time_sf%1)
            
            
        file_suffix0 = '%'+'0%i'%(time_sf+1)+'i'
        
        
    
    
    for tt in time_steps:
        
        timef = time_format%tt        
            
        #grid = grid0.NewInstance()
        # element values:
        
        if timef in bvp.ISVs.keys():
            
            
            print('**** Writing output for time stamp <<%s>>'%timef)
            
            # add time to vtk file:
            TA = vtk.vtkDoubleArray()
            TA.SetName('TIME')
            TA.SetNumberOfComponents(1)
            TA.InsertNextTuple1(tt)
            FD = vtk.vtkFieldData()
            FD.AddArray(TA)
            grid.SetFieldData(FD)
            grid.Modified()
            
            
            
            
            El_vals = bvp.ISVs[timef]
            
            # element averaging of each possible element ISV:
            
            #isvnames = ['temp_el','temp_rate','strain','stress']
            isvnames = ['strain','stress','cauchy','temp_rate','plastic strain','equivalent plastic']
            
            isv4comp = ['strain','stress','cauchy','plastic strain']
            
            
            
            for isn in isvnames:
                
                do_isn = False
                
                Iavg = vtk.vtkFloatArray()
                Iavg.SetName('EA_'+isn)
                
                # also nodal Averaged:
                #nd_avg_name = 'NA_'+isn
                nd_avg_name = isn
                
                ND_I = bvp.NodeVals[timef][nd_avg_name] = {}
                ND_W = bvp.NodeVals[timef]['NA_weights'] = {}
                
                
                #number of components:
                ncomp = 1
                #if (isn == 'strain')|(isn=='stress'):
                if isn in isv4comp:
                    ncomp = 4
                    Imises = vtk.vtkFloatArray()
                    Imises.SetName('EA_'+isn+'_vM')
                    
                    
                # loop over elements
                for el_id in El_vals.keys():
                    
                    ind = new_el_nrs[el_id]
                    gp_vals = El_vals[el_id]
                    #
                    el_avg = [0]*4
                    
                    if isn in gp_vals[1].keys():
                        
                        # if any of the elements have a particular ISV = output that ISV (others are 0)
                        do_isn = True
                        
                        isvR = {}
                        for gp in gp_vals.keys():
                            isvR[gp] = gp_vals[gp][isn]
                        
                        # element connectivity and coordinates:
                        el_connect = bvp.elements[el_id]
                        #X = [bvp.nodes[nd_Id][0] for nd_Id in el_connect]
                        #Y = [bvp.nodes[nd_Id][1] for nd_Id in el_connect]
                        #XY   = np.c_[X, Y];
                        
                        
                        
                        if list(gp_vals.keys()).__len__() == 4: # 2x2 gauss (weighting = 1 per integration point)
                            el_avg = np.average([gp_vals[i][isn] for i in range(1,5)],0).flatten()
                            
                            #
                            NDvals = FE_shapefn.Extrapolate_2x2(isvR,el_connect.__len__())
                            for nd_nr in range(el_connect.__len__()):
                                nd_Id = el_connect[nd_nr]
                                try:
                                    ND_I[nd_Id] += NDvals[nd_nr]
                                    ND_W[nd_Id] += 1
                                except:
                                    ND_I[nd_Id] = NDvals[nd_nr]
                                    ND_W[nd_Id] = 1
                                    
                        
                        
                    # Add the ISV to the float array
                    if ncomp == 4:
                        Iavg.SetNumberOfComponents(4) # 11, 22, 33, 12
                        Iavg.InsertTuple4(ind,el_avg[0],el_avg[1],el_avg[2],el_avg[3])
                        
                        mises_val = np.sqrt(((el_avg[0]-el_avg[1])**2 + 
                                    (el_avg[1]-el_avg[2])**2 + 
                                    (el_avg[2]-el_avg[0])**2 +
                                    6*el_avg[3]**2)/2)
                        
                        Imises.InsertTuple1(ind,mises_val)
                        
                    else:
                        Iavg.InsertTuple1(ind,el_avg[0])
                
                    
                    
                if do_isn:
                    
                    grid.GetCellData().SetActiveScalars('EA_'+isn)
                    grid.GetCellData().SetScalars(Iavg)
                    if ncomp == 4:
                        grid.GetCellData().SetActiveScalars('EA_'+isn+'_vM')
                        grid.GetCellData().SetScalars(Imises)
                        
                        
                    # rescale the nodal averaged ISVs:
                    for nd_Id in ND_I.keys():
                        ND_I[nd_Id] /= ND_W[nd_Id]   
                
                else:
                    del(ND_I) # the specific internal state variable is unavailable
                #
                # 
            
            
        
        if timef in bvp.NodeVals.keys():
            
            Nd_vals = bvp.NodeVals[timef]
            
            if 'temperature' in Nd_vals.keys():
                # temeparature is a scalar value:
                scalars = vtk.vtkFloatArray()
                scalars.SetName('Temperature')
                for ind in bvp.nodes.keys():
                    val = 0.
                    if ind in Nd_vals['temperature'].keys():
                        val = Nd_vals['temperature'][ind]
                
                #for ind,val in Nd_vals['temp_nd'].items():
                    scalars.InsertTuple1(ind,val)
                    
                grid.GetPointData().SetActiveScalars('Temperature')                    
                grid.GetPointData().SetScalars(scalars)
                
                
                
            if 'displacement' in Nd_vals.keys():
                
                vector = vtk.vtkFloatArray()
                vector.SetNumberOfComponents(3)
                vector.SetName('Displacement')
                for ind in bvp.nodes.keys():
                    val = [0., 0., 0.] # displacement field
                    if ind in Nd_vals['displacement'].keys():
                        val = Nd_vals['displacement'][ind] + [0.]
                
                #for ind,val in Nd_vals['temp_nd'].items():
                    vector.InsertTuple3(ind,val[0],val[1],val[2])
                    
                    
                grid.GetPointData().SetActiveVectors('Displacement')
                grid.GetPointData().SetVectors(vector)
                
                
                
            if 'strain' in Nd_vals.keys():
                
                scalars = vtk.vtkFloatArray()
                scalars.SetNumberOfComponents(6)
                scalars.SetName('Strain')
                for ind in bvp.nodes.keys():
                    val = [0.]*4
                    if ind in Nd_vals['strain'].keys():
                        val = list(np.array(Nd_vals['strain'][ind]).flatten()) + [0.]*4
                
                #for ind,val in Nd_vals['temp_nd'].items():
                    scalars.InsertTuple6(ind,val[0],val[1],val[2],val[3],0.,0.)
                        #scalars.InsertTuple4(ind,val[0],val[1],val[2],val[3])
                    
                grid.GetPointData().SetActiveScalars('Strain')
                grid.GetPointData().SetScalars(scalars)
                
                
                
                
            if 'plastic strain' in Nd_vals.keys():
                
                #print("DOING PLASTIC")
                #print(Nd_vals.keys())
                #print(Nd_vals['plastic strain'])
                
                scalars = vtk.vtkFloatArray()
                scalars.SetNumberOfComponents(6)
                scalars.SetName('Plastic Strain')
                for ind in bvp.nodes.keys():
                    val = [0.]*4
                    if ind in Nd_vals['plastic strain'].keys():
                        val = list(np.array(Nd_vals['plastic strain'][ind]).flatten()) + [0.]*4
                
                #for ind,val in Nd_vals['temp_nd'].items():
                    scalars.InsertTuple6(ind,val[0],val[1],val[2],val[3],0.,0.)
                        #scalars.InsertTuple4(ind,val[0],val[1],val[2],val[3])
                        
                    
                grid.GetPointData().SetActiveScalars('Plastic Strain')
                grid.GetPointData().SetScalars(scalars)
                
                
            
            if 'equivalent plastic' in Nd_vals.keys():
                
                grid.GetPointData().SetActiveScalars('PEEQ')
                scalars = vtk.vtkFloatArray()
                scalars.SetName('PEEQ')
                for ind in bvp.nodes.keys():
                    val = 0.
                    if ind in Nd_vals['equivalent plastic'].keys():
                        val = Nd_vals['equivalent plastic'][ind]
                
                #for ind,val in Nd_vals['temp_nd'].items():
                    scalars.InsertTuple1(ind,val)
                    
                grid.GetPointData().SetScalars(scalars)
                
                
                
            if 'stress' in Nd_vals.keys():
                
                scalars = vtk.vtkFloatArray()
                scalars.SetNumberOfComponents(6)
                scalars.SetName('Stress')
                
                ssVM = vtk.vtkFloatArray()
                ssVM.SetName('von Mises Stress')
                
                for ind in bvp.nodes.keys():
                    val = [0.]*4
                    if ind in Nd_vals['stress'].keys():
                        val = list(np.array(Nd_vals['stress'][ind]).flatten())
                                    
                    scalars.InsertTuple6(ind,val[0],val[1],val[2],val[3],0.,0.)
                    
                    mises_val = np.sqrt(((val[0]-val[1])**2 + (val[1]-val[2])**2 + (val[2]-val[0])**2 + 6*val[3]**2)/2)
                    
                    ssVM.InsertTuple1(ind,mises_val)
                    
                
                grid.GetPointData().SetActiveScalars('Stress')
                grid.GetPointData().SetScalars(scalars)
                grid.GetPointData().SetActiveScalars('von Mises Stress')
                grid.GetPointData().SetScalars(ssVM)
            
                
                
                
            
            
        grid.Modified()
                
        
        if nr_of_files>1:
            file_suffix = file_suffix0%(int(round(tt*10**time_sf)))
        
        
                
        vtkWriter = vtk.vtkUnstructuredGridWriter()
        vtkWriter.SetInputData(grid)
        vtkWriter.SetFileTypeToASCII()
        #vtkWriter.SetFileName(filename+'_'+timef+'.vtk')
        vtkWriter.SetFileName(filename+file_suffix+'.vtk')
        vtkWriter.Write()
        
        
        
        #vtuWriter.WriteNextTime(tt)
        
        
    #vtuWriter.Stop()
    
    # return
    return #vtuWriter,grid
            
            
            
            
        
    
    #vector = vtk.vtkFloatArray()
    #vector.SetNumberOfComponents(3)
    #for ind,U in enumerate(reslt['displ']):
        #vector.InsertTuple3(ind,U[0],U[1],0.)
        
    #vector.SetName('U')
    
    #grid.GetPointData().SetVectors(vector)
    
    
    ##
    
    ## stresses and other ISVs at nodes:
    #if 'ISVnodes' in reslt.keys():
        ## do VM:
        #grid.GetPointData().SetActiveScalars('Mises')
        #vmscalar = vtk.vtkFloatArray()
        #vmscalar.SetName('Mises')
        #for ind,val in enumerate(reslt['ISVnodes']['Mises']):
            #vmscalar.InsertTuple1(ind,val)
        ##
        #grid.GetPointData().SetScalars(vmscalar)
        
        ##
        
        ## do VM:
        #grid.GetPointData().SetActiveScalars('S')
        #scalars = vtk.vtkFloatArray()
        #scalars.SetNumberOfComponents(4)
        #scalars.SetName('S')
        #for ind,val in enumerate(reslt['ISVnodes']['stress']):
            #scalars.InsertTuple4(ind,val[0],val[1],val[2],val[3])
        ##
        #grid.GetPointData().SetScalars(scalars)
        
        
        # do STRAINS:
        #grid.GetPointData().SetActiveScalars('S')
        #scalars = vtk.vtkFloatArray()
        #scalars.SetNumberOfComponents(4)
        #scalars.SetName('S')
        #for ind,val in enumerate(reslt['ISVnodes']['stress']):
            #scalars.InsertTuple4(ind,val[0],val[1],val[2],val[3])
        ##
        #grid.GetPointData().SetScalars(scalars)
            
        
        
    # use or calculate nodeISVs :nodeISVs = {}
    #for elnr,ndlst in enumerate(EL):
        
    
    #grid.GetPointData().SetActiveScalars('aa')

    ## scalar point data:
    #scalars = vtk.vtkFloatArray()
    #for ind,U in enumerate(reslt['displ']):
        #scalars.InsertTuple1(ind,U[1])#,U[1],0.)
        
    #scalars.SetName('aa')
    
    #grid.GetPointData().SetScalars(scalars)
    
    
    
    
    # vtk XML file writer
    #vtuWriter = vtk.vtkXMLPUnstructuredGridWriter()
    #vtuWriter.SetInputData(grid)
    #vtuWriter.SetFileName(filename+'.vtu')
    #vtuWriter.Update()
    
    #return #grid
    
    
    


#
#
#
#  POST_PROCESS
#
##  Interpolate stresses and strains or other internal state variables to nodes
##

def isv_to_nodes(mesh,ISVs,dofor=['stress','strain'],weight_by_centre_distance=True):#,GetMises=True,GetEig=True,):
    ##
    ## loop through elements and 
    
    
    #
    # only in 2d with either 3 or 4 stress and strain values (plane stress / strain?)   #
    #
    nodeISVs = {} # estrapolate all ISVs to nodes:
    
    XY = mesh['nodal_coord']
    EL = mesh['connectivity']
    
    nnodes = XY.shape[0]
    
    weights = np.zeros(nnodes)
    stresscomps = np.zeros((nnodes,5))
    straincomps = np.zeros((nnodes,5))
    
    
    for elnr,elISV in ISVs.items():
        cell = EL[elnr]
        # reorder ISVs in dict where arrays have name == statevariable
        isvR = {}
        isvkeys = elISV[1].keys()
        for kn in isvkeys:
            isvR[kn] = {}
            for gp in elISV.keys():
                isvR[kn][gp] = elISV[gp][kn]
        #
        nnds = cell.size
        # weight calculated as inverse distance from centre:
        el_xy = XY[cell,:]
        if weight_by_centre_distance:
            dist = el_xy - np.average(el_xy,0)
            el_weight = 1./np.sqrt(np.sum(dist**2,1))
        else: 
            el_weight = 1.
            
        weights[cell]+=el_weight
        #
        #for doisv in dofor:
        # stress at nodes
        if nnds == 8:
            nodestress = np.array(FE_shapefn.Extrapolate_2x2(isvR['stress'],8))
        # add weighted values to component array:
        stresscomps[cell,:nodestress.shape[1]] += (nodestress.T*el_weight).T
        
        #stresscomps[cell,:nodestress.shape[1]] = nodestress
        
    
    #allstresses = stresscomps#
    allstresses = (stresscomps.T/weights).T
    stresscomps = [] # clear memory
    nodeISVs['stress'] = allstresses
    # von Mises :
    sVM = (allstresses[:,0]-allstresses[:,1])**2
    sVM += (allstresses[:,1]-allstresses[:,2])**2
    sVM += (allstresses[:,2]-allstresses[:,0])**2
    sVM += 6*allstresses[:,3]**2
    nodeISVs['Mises'] = np.sqrt(sVM/2)
    
    
    return nodeISVs


