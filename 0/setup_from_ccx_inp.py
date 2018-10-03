import numpy as np
from time import time
import os

class bvp_dict:
    '''
    create a boundary value problem dictionary from a CalculiX input file
    
    '''    
    def __init__(self,filename):
        # make sure the filename ends in ".inp"
        if not filename[-4:] == '.inp': filename = filename+'.inp'
        #
        if os.path.exists(filename):
            fid = open(filename,'r')
        else:
            print('\n**\tNO CALCULIX INPUT FILE AT << %s >>'%filename)
            return
        
        # inherit name from file:
        self._name = filename.split('/')[-1].split('.')[0]
        self._folder = filename.split(self._name)[0]
        #
        # keywords recognised:
        self._keywrds = ['*node',
                        '*nset',
                        '*element',
                        '*elset',
                        '*boundary',
                        '*material',
                        '*solid section',
                        #'*shell section', ## shell elements aren't implemented
                        '*contact pair',
                        '*surface interaction',
                        '*surface',
                        '*initial conditions',
                        '*step']
        
        # known element types and nodes:
        self._element_nodes = {'cpe4':4,    # plane strain 4 node quadrilateral
                               'cpe4r':4,   #              " with reduced integration
                               'cpe8':8,    # plane strain 8 node quadrilateral
                               'cpe8r':8    #              " with reduced integration 
                               #'cps4':4,
                               #'cax4':4,
                               }
        
        # assign nodal coodinates as a dictionary of [x,y,z] lists with node index as key
        self.nodes = {}  # read_in  = 1
        # assign groups of nodes as dictionary of node lists
        self.nd_set = {}
        self.elements = {} # groups of elements
        self._el_type = {}
        self.el_set = {}
        #
        # section definitions: (solid / fluid / shell etc)
        self.sections = {}
        #
        # surfaces:
        self.surfaces = {}
        #
        # contact definitions
        self.contact = {}
        #self.element_material = {}
        
        # dictionary of material dictionaries
        self.materials = {}
        # dictionary of fixed boundary values
        #self.BC_fixed = {}
        # dictionary of initial conditions
        #self.initial = {}
        # dictionary of steps:
        self._nseps = 0 
        self.steps = {}
        # step 0 contains the boundary values defined ouside of step definitions
        # (i.e. always applied first and may be empty)
        self.steps[0] = {} 
        #
        # amplitudes:
        self.amplitudes = {}
        #
        # time stamp format and list of time steps:
        self.time_format = '%.8f'
        self.time_steps = []
        #
        # internal state variables:
        self.solved_time = 0
        self.ISVs = {}
        self.ISVs[self.time_format%self.solved_time]={}
        # nodal variables:
        self.NodeVals = {}
        self.NodeVals[self.time_format%self.solved_time]={}
        
        
        seekpos = 0
        newseek = 1
        while seekpos<newseek:
            seekpos = fid.tell()
            lntxt = fid.readline().lower()
            newseek = fid.tell()
            
            keyword = self.__contains_keyword(lntxt)
            
            if not keyword is None:
                #print('KEYWORD found : << %s >>'%keyword)
                newseek = self.__add_attribute(keyword,fid,seekpos)
                fid.seek(newseek)
        
        
        
     
     
     
     
     
     
     
    def __contains_keyword(self,linetxt):
        '''
        does the line text contain a keyword
        '''
        
        # first check if commented out
        if linetxt[:2] =='**':
            return None
        
        trykeys = [kn in linetxt.lower() for kn in self._keywrds]
        if any(trykeys):
            kwrd = [self._keywrds[ki] for ki in range(self._keywrds.__len__()) if trykeys[ki]][0]
            #print("************ %s "%kwrd)
            return kwrd
        
        else:
            return None
        
      
      
      
      
      
      
    def __add_attribute(self,keyword,fid,seekpos):
        
        if keyword == '*node':
            newseek = self.__add_nodes(fid,seekpos)
            #print(seekpos,newseek)
            
        elif keyword == '*nset':
            newseek = self.__add_nset(fid,seekpos)
            
        elif keyword == '*element':
            newseek = self.__add_elems(fid,seekpos)
            
        elif keyword == '*elset':
            newseek = self.__add_elset(fid,seekpos)
            
        elif keyword == '*boundary':
            newseek = self.__add_boundary(fid,seekpos)
            
        elif keyword == '*material':
            newseek = self.__add_material(fid,seekpos)
            
        elif keyword == '*solid section':
            newseek = self.__add_section(fid,seekpos)
            
        elif keyword == '*contact pair':
            newseek = self.__add_contact_pair(fid,seekpos)   
            
        elif keyword == '*surface interaction':
            newseek = self.__add_surf_interaction(fid,seekpos) 
            
        elif keyword == '*surface':
            newseek = self.__add_surface(fid,seekpos)              
            
        elif keyword == '*initial conditions':
            newseek = self.__add_ini(fid,seekpos)
            
        elif keyword == '*step':
            newseek = self.__add_step(fid,seekpos)
            
        #elif keyword == '*':
            #newseek = self.__add_(fid,seekpos)
            
            
            
        else: newseek = seekpos+1
        
        return newseek
    
        
        
        
        
#   ADD NODAL COORDINATES

        
        
        
    def __add_nodes(self,fid,seekpos=0):
        '''
        
        add nodes to list of node values, i.e.:
        
            self.nodes[node ID] = [x, y, z]
        
        '''        
        tic = time()
        nnodes = 0
        fid.seek(seekpos)
        lntxt = fid.readline().lower() # this line needs to contain the keyword "*nodes" and possibly a "*nset" flag
        newseek = fid.tell()
        if not "*node" in lntxt:
            print("WRONG ASSIGNMENT - no << *NODE >> flag in expected input file line")
            return newseek
        
        if 'file' in lntxt: return newseek
    
        # check if "nset" in the text
        nsetnames = ['000'] # default all nodes are added to nodeset "000"
        if "nset" in lntxt:
            nsetnames += [nn for nn in lntxt.split('nset')[1].split('=')[1].split('\n')[0].split(' ') if nn.__len__()>0]
        #        
        # while loop to add node info to repective dictionaries
        doseek = True
        while doseek:
            seekpos = fid.tell()
            lntxt = fid.readline().lower()
            newseek = fid.tell()
            if not newseek>seekpos:
                print('Assigned %i nodal coordinates to BVP in %.4f seconds'%(nnodes,time()-tic))
                return newseek
            
            #
            # try to see if linetext conforms to Id, X, Y, Z requirement (OR Z=0 assumed)
            lnarr = lntxt.split('\n')[0].split(',')
            if lnarr.__len__()==3: lnarr+=['0.'] # add z=0 if 2d
            
            # try to add node to dictionary
            try:
                nd_ID = int(lnarr[0])
                nd_XYZ = [float(lnarr[1]),float(lnarr[2]),float(lnarr[3])]
                # add node info
                self.__add_nd__(nd_ID,nd_XYZ,nsetnames)
                nnodes += 1
                
            except:
                # check if new keyword activated:
                keyword = self.__contains_keyword(lntxt)
                
                if not keyword is None:
                    #print ('keyword found')
                    #print('Keyword = << %s >>'%keyword)
                    print('Assigned %i nodal coordinates to BVP in %.4f seconds'%(nnodes,time()-tic))
                    return seekpos
                
        
        
    def __add_nd__(self,nd_ID,nd_XYZ,nsetnames=[]):
        self.nodes[nd_ID] = nd_XYZ
        for nset in nsetnames:
            if not nset in self.nd_set.keys():
                self.nd_set[nset] = [nd_ID]
            else:
                self.nd_set[nset] += [nd_ID]
        
        
        
        
##  ADD A NODE SET
        
    def __add_nset(self,fid,seekpos=0):
        '''
        
        add lists of nodes to existing nodeset or create new node set:
        
            self.nd_set[setname] = [i,i,i,...i]
        
        '''  
        print("--> ADDING NODE SET")
        tic = time()
        fid.seek(seekpos)
        lntxt = fid.readline().lower() # this line needs to contain the keyword "*nset" and possibly a "generate" flag
        newseek = fid.tell()
        if not "*nset" in lntxt:
            print("WRONG ASSIGNMENT - no << *NSET >> flag in expected input file line")
            return newseek
        
        #if 'file' in lntxt: return newseek
        
        # check if "nset" in the text
        nsetnames = [nn for nn in lntxt.split('nset')[2].split('=')[1].split('\n')[0].split(',')[0].split(' ') if nn.__len__()>0]
        #
        #
        # check if nset names are already defined (else create empty list):
        for nset in nsetnames:
            if not nset in self.nd_set.keys():
                self.nd_set[nset] = []
        #
        # if "generate" the next line may be 2 or 3 different integers (initial,final and interval)
        if 'generate' in lntxt:
            seekpos = fid.tell()
            lntxt = fid.readline().lower()
            newseek = fid.tell()
            # try or return:
            try:
                lnarr = [int(xx) for xx in lntxt.split('\n')[0].split(',')]
                intvl = 1
                if lnarr.__len__() == 3:
                    intvl = lnarr[2]
                self.nd_set[nset] = list(range(lnarr[0],lnarr[1],intvl))+[lnarr[1]]
                
            except:
                print("** 2 or 3 integer values expected in line after << GENERATE >> for << NSET = %s>>"%nset)
                
            return seekpos
            
        
        
        # while loop to add node info to repective dictionaries
        doseek = True
        while doseek:
            seekpos = fid.tell()
            lntxt = fid.readline().lower()
            newseek = fid.tell()
            if not newseek>seekpos:
                return newseek
            #
            # try to see if linetext conforms to list of id's or single node ID
            lnarr = lntxt.split('\n')[0].split(',')
            
            # try to add node to dictionary
            try:
                for kn in lnarr:
                    nd_ID = int(kn)
                    for nset in nsetnames:
                        self.nd_set[nset] += [nd_ID]
                
            except:
                # check if new keyword activated:
                keyword = self.__contains_keyword(lntxt)
                
                if not keyword is None:
                    print('Assigned node set %s'%nset)
                    return seekpos
        
        
        
        
    
#  ADD ELEMENTS
        
        
        
    def __add_elems(self,fid,seekpos=0):
        '''
        
        add elements to list of element definition , i.e.:
        
            self.elements[element ID] = [node_1, node_2, ...., node_n] #number of nodes depends on the element type
        
        '''        
        tic = time()
        nels = 0
        fid.seek(seekpos)
        lntxt = fid.readline().lower() # this line needs to contain the keyword "*element", a "type" and possibly an "*elset" flag
        newseek = fid.tell()
        
        if not "*element" in lntxt:
            print("WRONG ASSIGNMENT - no << *ELEMENT >> flag in expected input file line")
            return newseek


        ##what is the element type?
        eltype = None
        el_nodes = 0
        if 'type' in lntxt:
            try:
                eltype = [nn for nn in lntxt.split('type')[1].split('=')[1].split('\n')[0].split(',')[0].split(' ') if nn.__len__()>0][0]
            except:
                pass
        #
        if eltype is None:
            print("NO ELEMENT TYPE defined after << *ELEMENT >> flag in input file")
            return newseek
        
        if eltype in self._element_nodes.keys():
            el_nodes = self._element_nodes[eltype]
        else:
            print("UNKNOWN ELEMENT TYPE << %s >> defined after << *ELEMENT >> flag in input file"%eltype)
            return newseek
        
        
        
        
        # check if "elset" in the text
        elsetnames = ['000'] # default all nodes are added to nodeset "000"
        if "elset" in lntxt:
            elsetnames += [nn for nn in lntxt.split('elset')[1].split('=')[1].split('\n')[0].split(',')[0].split(' ') if nn.__len__()>0]
        
        
        # while loop to add element info to repective dictionaries
        doseek = True
        while doseek:
            seekpos = fid.tell()
            lntxt = fid.readline().lower()
            newseek = fid.tell()
            if not newseek>seekpos:
                print('Assigned %i, << %s >> type element definitions to BVP in %.4f seconds'%(nels,eltype,time()-tic))
                return newseek
            
            # try to add element to dictionary
            try:
                # try to see if linetext conforms to Id, nd_1, nd_2, .., nd_n (all integers)
                lnarr = [int(k) for k in lntxt.split('\n')[0].split(',')]
                #print(lnarr)
                el_ID = int(lnarr[0])
                el_connect = lnarr[1:el_nodes+1]
                # add element info
                self.__add_el__(el_ID,el_connect,eltype,elsetnames)
                nels += 1
                
            except:
                # check if new keyword activated:
                keyword = self.__contains_keyword(lntxt)
                
                if not keyword is None:
                    #print ('keyword found')
                    #print('Keyword = << %s >>'%keyword)
                    print('Assigned %i, << %s >> type element definitions to BVP in %.4f seconds'%(nels,eltype,time()-tic))
                    return seekpos
        

# ADD AN ELEMENT
        
    def __add_el__(self,el_ID,el_connect,eltype,elsetnames=[]):
        self.elements[el_ID] = el_connect
        self._el_type[el_ID]=eltype
        # add to list of element sets
        for elset in elsetnames:
            if not elset in self.el_set.keys():
                self.el_set[elset] = [el_ID]
            else:
                self.el_set[elset] += [el_ID]


#
#
#
#   ADD AN ELEMENT SET
        
    def __add_elset(self,fid,seekpos=0):
        '''
        
        add lists of elements to existing sets or create a new set:
        
            self.el_set[setname] = [i,i,i,...i]
        
        '''  
        print("--> ADDING ELEMENT SET")
        tic = time()
        fid.seek(seekpos)
        lntxt = fid.readline().lower() # this line needs to contain the keyword "*elset" and possibly a "generate" flag
        newseek = fid.tell()
        if not "*elset" in lntxt:
            print("WRONG ASSIGNMENT - no << *ELSET >> flag in expected input file line")
            return newseek
        
        
        # check if "elset" in the text
        elsetnames = [nn for nn in lntxt.split('elset')[2].split('=')[1].split('\n')[0].split(',')[0].split(' ') if nn.__len__()>0]
        # check if elset names are already defined (else create empty list):
        for elset in elsetnames:
            if not elset in self.el_set.keys():
                self.el_set[elset] = []
                
        # if "generate" the next line may be 2 or 3 different integers (initial,final and interval)
        if 'generate' in lntxt:
            seekpos = fid.tell()
            lntxt = fid.readline().lower()
            newseek = fid.tell()
            # try or return:
            try:
                lnarr = [int(xx) for xx in lntxt.split('\n')[0].split(',')]
                intvl = 1
                if lnarr.__len__() > 2:
                    intvl = lnarr[2]
                self.el_set[elset] = list(range(lnarr[0],lnarr[1],intvl))+[lnarr[1]]
                
            except:
                print("** 2 or 3 integer values expected in line after << GENERATE >> for << ELSET = %s>>"%elset)
                
            return seekpos
        
        # while loop to add node info to repective dictionaries
        doseek = True
        while doseek:
            seekpos = fid.tell()
            lntxt = fid.readline().lower()
            newseek = fid.tell()
            if not newseek>seekpos:
                return newseek
            #
            # try to see if linetext conforms to list of id's or single node ID
            lnarr = lntxt.split('\n')[0].split(',')
            
            # try to add node to dictionary
            try:
                for kn in lnarr:
                    nd_ID = int(kn)
                    for elset in elsetnames:
                        self.el_set[elset] += [nd_ID]
                
            except:
                # check if new keyword activated:
                keyword = self.__contains_keyword(lntxt)
                
                if not keyword is None:
                    print('Assigned element set %s'%elset)
                    return seekpos
            
        
        
        
        
        
        
        
        
##  ADD A MATERIAL DEFINITION
        
    def __add_material(self,fid,seekpos=0):
        '''
        
        add lists of elements to existing sets or create a new set:
        
            self.materials[matname] = {}
        
        '''  
        print("--> ADDING A MATERIAL DEFINITION")
        tic = time()
        fid.seek(seekpos)
        lntxt = fid.readline().lower() # this line needs to contain the keyword "*material" and a "NAME" flag
        newseek = fid.tell()
        
        
        if not "*material" in lntxt:
            print("WRONG ASSIGNMENT - no << *MATERIAL >> flag in expected input file line")
            return newseek
        
        
        # get name for "material" in the text
        try:
            matname = [nn for nn in lntxt.split('name')[1].split('=')[1].split('\n')[0].split(',')[0].split(' ') if nn.__len__()>0][0]
        except:
            print("WRONG ASSIGNMENT - << *MATERIAL >> requires a << NAME >> flag ")
            return newseek
        
        
        ##check if mat names are already defined (else create empty list):
        if not matname in self.materials.keys():
            self.materials[matname] = {}
            
        material = self.materials[matname]
        
        ## list of allowable material keywords:
        #matkeys = ['*elastic',
                   #'*plastic',
                   #'*density',
                   #'*conductivity',
                   #'*expansion',
                   #'*specific heat',
                   ##'*',
                   ##'*',
                   ##'*',
                   ##'*',
                   ##'*',
                   ##'*',
                   ##'*'
                   #]
                
                
        ## while loop to add material info to the material dictionary
        doseek = True
        while doseek:
            seekpos = fid.tell()
            lntxt = fid.readline().lower()
            newseek = fid.tell()
            if not newseek>seekpos:
                return newseek
                   
    
    
    
    
    
    
    
            # read ELASTIC PROPERTIES
            if '*elastic' in lntxt:
                material['*elastic'] = {}
                # (elastic also requires a type --> alternatively << Iso >> is assumed)
                etype = 'iso'
                if 'type' in lntxt:
                    for ae in ['iso','ortho','aniso','engineering constants']:
                        if ae in lntxt: etype=ae
                #
                material['*elastic']['type'] = etype
                do_key = True
                while do_key:
                    # if << ISO >>, read in the stiffness and poisson's ratio for each temperature:
                    if etype == 'iso':
                        
                        ctmp = 0 # if no temperature is specified, the << 0 >> value will apply to all temperatures, else its interpolated
                        seekpos = fid.tell()
                        lntxt = fid.readline().lower()
                        try:
                            lnarr = [float(x) for x in lntxt.split('\n')[0].split(',')]
                            if lnarr.__len__() == 2:
                                material['*elastic'][ctmp] = lnarr
                            elif lnarr.__len__() == 3:
                                material['*elastic'][lnarr[2]] = lnarr[:2]
                            else:
                                print(" ** INPUT FILE: << *ELASTIC >> of TYPE = ISO requires 2 or 3 floats in the following line(s)\n\t\tStiffness, Poisson, Temperature (optional)")
                                do_key = False
                        except:
                            do_key = False
                            fid.seek(seekpos)
                        # lntxt
                     
                    # if  << ORTHO >> read in a line pair 
                    #if etype == 'ortho':
                        
                        
                    else:
                        print(" ** INPUT FILE: !!! ONLY ELASTIC of TYPE = ISO currently supported !!")
                        do_key=False
                
                #
                
                
    
    
    
    
    
    
            # read conductivity
            elif '*conductivity' in lntxt:
                material['*conductivity'] = {}
                # (conductivity also requires a type --> alternatively << Iso >> is assumed)
                etype = 'iso'
                if 'type' in lntxt:
                    for ae in ['iso','ortho','aniso']:
                        if ae in lntxt: etype=ae
                #
                material['*conductivity']['type'] = etype
                do_key = True
                while do_key:
                    # if << ISO >>, read in the stiffness and poisson's ratio for each temperature:
                    if etype == 'iso':
                        
                        ctmp = 0 # if no temperature is specified, the << 0 >> value will apply to all temperatures, else its interpolated
                        seekpos = fid.tell()
                        lntxt = fid.readline().lower()
                        try:
                            lnarr = [float(x) for x in lntxt.split('\n')[0].split(',')]
                            if lnarr.__len__() == 1:
                                material['*conductivity'][ctmp] = lnarr[0]
                            elif lnarr.__len__() == 2:
                                material['*conductivity'][lnarr[1]] = lnarr[0]
                            else:
                                print(" ** INPUT FILE: << *CONDUCTIVITY >> of TYPE = ISO requires 1 or 2 floats in the following line(s)\n\t\tConductivity, Temperature (optional)")
                                do_key = False
                        except:
                            do_key = False
                            fid.seek(seekpos)
                            
                            
                            
                        # lntxt
                     
                    # if  << ORTHO >> read in a line pair 
                    #if etype == 'ortho':
                        
                        
                    else:
                        print(" ** INPUT FILE: !!! ONLY CONDUCTIVITY of TYPE = ISO currently supported !!")
                        do_key=False
                
                #
                
                
    
    
            # read expansion
            elif '*expansion' in lntxt:
                material['*expansion'] = {}
                # (expansion also requires a type --> alternatively << Iso >> is assumed)
                etype = 'iso'
                if 'type' in lntxt:
                    for ae in ['iso','ortho','aniso']:
                        if ae in lntxt: etype=ae
                #check if zero value is defined:
                zero = 0.
                if 'zero' in lntxt:
                    try:
                        zero = float(lntxt.split('\n')[0].split('zero')[1].split('=')[1].split(',')[0])
                    except:
                        print(" ** INPUT FILE: << EXPANSION >> Zero modification incorrect, using Zero = 0.0")
                #
                
                material['*expansion']['type'] = etype
                material['*expansion']['zero'] = zero
                do_key = True
                while do_key:
                    # if << ISO >>, read in the stiffness and poisson's ratio for each temperature:
                    if etype == 'iso':
                        
                        ctmp = 0 # if no temperature is specified, the << 0 >> value will apply to all temperatures, else its interpolated
                        seekpos = fid.tell()
                        lntxt = fid.readline().lower()
                        try:
                            lnarr = [float(x) for x in lntxt.split('\n')[0].split(',')]
                            if lnarr.__len__() == 1:
                                material['*expansion'][ctmp] = lnarr[0]
                            elif lnarr.__len__() == 2:
                                material['*expansion'][lnarr[1]] = lnarr[0]
                            else:
                                print(" ** INPUT FILE: << *EXPANSION >> of TYPE = ISO requires 1 or 2 floats in the following line(s)\n\t\tExpansion coefficient, Temperature (optional)")
                                do_key = False
                        except:
                            do_key = False
                            fid.seek(seekpos)
                            
                            
                            
                        # lntxt
                     
                    # if  << ORTHO >> read in a line pair 
                    #if etype == 'ortho':
                        
                        
                    else:
                        print(" ** INPUT FILE: !!! ONLY EXPANSION of TYPE = ISO currently supported !!")
                        do_key=False
                
                #
                
                
                
                
                
                
                
                '*specific heat'
            # read specific heat:
            elif '*specific heat' in lntxt:
                material['*specific heat'] = {}
                do_key = True
                while do_key:
                    ctmp = 0 # if no temperature is specified, the << 0 >> value will apply to all temperatures, else its interpolated
                    seekpos = fid.tell()
                    lntxt = fid.readline().lower()
                    try:
                        lnarr = [float(x) for x in lntxt.split('\n')[0].split(',')]
                        if lnarr.__len__() == 1:
                            material['*specific heat'][ctmp] = lnarr[0]
                        elif lnarr.__len__() == 2:
                            material['*specific heat'][lnarr[1]] = lnarr[0]
                        else:
                            print(" ** INPUT FILE: << *SPECIFIC HEAT >> requires 1 or 2 floats in the following line(s)\n\t\tSpecific heat, Temperature (optional)")
                            do_key = False
                    except:
                        do_key = False
                        fid.seek(seekpos)
                        
                        
                        
                        
                
                
                
            # read density:
            elif '*density' in lntxt:
                material['*density'] = {}
                do_key = True
                while do_key:
                    ctmp = 0 # if no temperature is specified, the << 0 >> value will apply to all temperatures, else its interpolated
                    seekpos = fid.tell()
                    lntxt = fid.readline().lower()
                    try:
                        lnarr = [float(x) for x in lntxt.split('\n')[0].split(',')]
                        if lnarr.__len__() == 1:
                            material['*density'][ctmp] = lnarr[0]
                        elif lnarr.__len__() == 2:
                            material['*density'][lnarr[1]] = lnarr[0]
                        else:
                            print(" ** INPUT FILE: << *DENSITY >> requires 1 or 2 floats in the following line(s)\n\t\tMass density, Temperature (optional)")
                            do_key = False
                    except:
                        do_key = False
                        fid.seek(seekpos)
                
                
            
            
            
            
            # read in plastic table:
            elif '*plastic' in lntxt:
                material['*plastic'] = {}
                # (elastic also requires a type --> alternatively << Iso >> is assumed)
                etype = 'isotropic'
                if 'hardening' in lntxt:
                    for ae in ['isotropic','kinematic','combined']:
                        if ae in lntxt: etype=ae
                #
                material['*plastic']['hardening'] = etype
                do_key = True
                while do_key:
                    # if << ISO >>, read in the stiffness and poisson's ratio for each temperature:
                    if etype == 'isotropic':
                        
                        ctmp = 0 # if no temperature is specified, the << 0 >> value will apply to all temperatures, else its interpolated
                        seekpos = fid.tell()
                        lntxt = fid.readline().lower()
                        try:
                            lnarr = [float(x) for x in lntxt.split('\n')[0].split(',')]
                            if lnarr.__len__() == 2:
                                try:
                                    material['*plastic'][ctmp] += lnarr
                                except:
                                    material['*plastic'][ctmp] = lnarr
                                    
                            elif lnarr.__len__() == 3:
                                try:
                                    material['*plastic'][lnarr[2]] += lnarr[:2]
                                except:
                                    material['*plastic'][lnarr[2]] = lnarr[:2]
                                    
                            else:
                                print(" ** INPUT FILE: << *ELASTIC >> of TYPE = ISO requires 2 or 3 floats in the following line(s)\n\t\tStiffness, Poisson, Temperature (optional)")
                                do_key = False
                        except:
                            do_key = False
                            fid.seek(seekpos)
                        # lntxt
                     
                    # if  << ORTHO >> read in a line pair 
                    #if etype == 'ortho':
                        
                        
                    else:
                        print(" ** INPUT FILE: !!! only ISOTROPIC PLASTICITY in table format is currently supported !!")
                        do_key=False
                
                #
                
                
                
                
                
                    
            else:
                ## check if new keyword activated:
                keyword = self.__contains_keyword(lntxt)
            
                if not keyword is None:
                    print('Assigned material << %s >> in %.4f seconds'%(matname,time()-tic))
                    return seekpos
                
                
                
                
                
                
                
            
            
##  ADD A SECTION DEFINITION
    
    
            
    def __add_section(self,fid,seekpos=0):
        '''
        
        add a solid section definition
        
            self.sections[elset] = {}
        
        '''  
        print("--> ADDING A SECTION DEFINITION")
        tic = time()
        fid.seek(seekpos)
        lntxt = fid.readline().lower().split('\n')[0] # this line needs to contain the keyword "*solid section"
        newseek = fid.tell()
        
        if not "*solid section" in lntxt:
            print("WRONG ASSIGNMENT - no << *SOLID SECTION >> flag in expected input file line")
            return newseek
        
        # assign to which element set?
        elset = lntxt.replace(" ","").split('elset=')[1].split(',')[0]
        matname = lntxt.replace(" ","").split('material=')[1].split(',')[0]
        
        self.sections[elset] = {'*type':'*solid section',
                                          '*material':matname}
        # add thickness (next line)
        try:
            self.sections[elset]['*thickness'] = float(fid.readline())
        except:
            pass
        
        return newseek
    
    
    
              
    def __add_contact_pair(self,fid,seekpos=0):
        '''
        
        add a contact pair definition
        
        '''  
        print("--> ADDING A CONTACT PAIR DEFINITION")
        tic = time()
        fid.seek(seekpos)
        lntxt = fid.readline().lower().split('\n')[0] # this line needs to contain the keyword "*contact pair"
        newseek = fid.tell()
        
        if not "*contact pair" in lntxt:
            print("WRONG ASSIGNMENT - no << *CONTACT PAIR >> flag in expected input file line")
            return newseek
        
        # 
        if 'interaction' in lntxt:
            int_name = lntxt.replace(" ","").split('interaction=')[1].split(',')[0]
        else:
            int_name = 'default'
            
        #
        #
        if 'type' in lntxt:
            int_type = lntxt.replace(" ","").split('type=')[1].split(',')[0]
        else:
            int_type = 'nodetosurface'
        
        #
        # contact pair
        lntxt = fid.readline().lower().split('\n')[0] # this line needs to contain the slave and master surface names 
        newseek = fid.tell()
        
        if not '*pairs' in self.contact.keys():
            self.contact['*pairs'] = {}
        
        try:
            Sslav,Smast = lntxt.replace(" ","").split(',')[:2]
            self.contact['*pairs'][Sslav+'--'+Smast] = {'*interaction':int_name,
                                             '*type':int_type
                }
            
        except:
            print("MAKE SURE << *CONTACT PAIR >>  FLAG IS FOLLOWED BY [slave, master] SURFACE NAMES")
            
        #
        # check if the interaction is there, else add the default values:
        if '*surface interaction' in self.contact.keys():
            if not int_name in self.contact['*surface interaction'].keys():
                self.contact['*surface interaction'][int_name] = {'pressure-overclosure':'linear',
                                                              'params':[1e6,None,None],
                                                             }
        else:
            self.contact['*surface interaction'] = {int_name:{'pressure-overclosure':'linear',
                                                              'params':[1e6,None,None],
                                                             }}
        
        return newseek
    
    
    
    
              
    def __add_surf_interaction(self,fid,seekpos=0):
        '''
        
        add a surface interaction definition
        
        '''  
        print("--> ADDING A SURFACE INTERACTION DEFINITION")
        tic = time()
        fid.seek(seekpos)
        lntxt = fid.readline().lower().split('\n')[0] # this line needs to contain the keyword '*surface interaction'
        newseek = fid.tell()
        
        if not '*surface interaction' in lntxt:
            print("WRONG ASSIGNMENT - no << *SURFACE INTERACTION >> flag in expected input file line")
            return newseek
        
        # requires a name
        try:
            int_name = [nn for nn in lntxt.split('name')[1].split('=')[1].split('\n')[0].split(',')[0].split(' ') if nn.__len__()>0][0]
        except:
            print("WRONG ASSIGNMENT - << *SURFACE INTERACTION >> requires a << NAME >> flag ")
            return newseek
        
        
        # surface behaviour
        lntxt = fid.readline().lower().split('\n')[0] # this line needs to contain the surface behaviour
        newseek = fid.tell()
        
        press_over = 'linear'
        press_params = [1e6,None,None]
        
        if 'pressure-overclosure' in lntxt:
            press_over = lntxt.replace(" ","").split('pressure-overclosure=')[1].split(',')[0]
            
        # parameters:
        param_lst = []
        doseek = True
        while doseek:
            newseek = fid.tell()
            lntxt = fid.readline().lower() # this line needs to contain the surface behaviour parameters
            seekpos = fid.tell()
            
            try:
                lnarr = [float(x) for x in lntxt.split('\n')[0].split(',') if x.__len__()>0]
            except:
                lnarr = []
                
            if not seekpos>newseek:
                doseek = False
                
            elif lnarr.__len__()>0:
                for pval in lnarr:
                    param_lst += [pval]
            
            else:
                ## check if new keyword activated:
                keyword = self.__contains_keyword(lntxt)
            
                if not keyword is None:
                    doseek = False
                   
        if param_lst.__len__()>0:
            press_params = param_lst + [None, None]
        
        
        
        
        # check if the interaction is there, else add the default values:
        if '*surface interaction' in self.contact.keys():
            self.contact['*surface interaction'][int_name] = {'pressure-overclosure':press_over,
                                                              'params':press_params,
                                                             }
        else:
            self.contact['*surface interaction'] = {int_name:{'pressure-overclosure':press_over,
                                                              'params':press_params,
                                                             }}
        
        
        
        return newseek
    
    
    
    
    
    
    
    
              
    def __add_surface(self,fid,seekpos=0):
        '''
        
        add a surface interaction definition
        
        '''  
        print("--> ADDING A SURFACE DEFINITION")
        tic = time()
        fid.seek(seekpos)
        lntxt = fid.readline().lower().split('\n')[0] # this line needs to contain the keyword '*surface interaction'
        newseek = fid.tell()
        
        if not '*surface' in lntxt:
            print("WRONG ASSIGNMENT - no << *SURFACE >> flag in expected input file line")
            return newseek
          
          
        # requires a name
        try:
            surf_name = [nn for nn in lntxt.split('name')[1].split('=')[1].split('\n')[0].split(',')[0].split(' ') if nn.__len__()>0][0]
        except:
            print("WRONG ASSIGNMENT - << *SURFACE >> requires a << NAME >> flag ")
            return newseek
        
        #
        # surface type: (defined as element faces or node set)
        surf_type = 'segments'
        if 'type=node' in lntxt.replace(' ',''):
            surf_type = 'node'
        
        
        # seek and add ot surface definition:
        surf = self.surfaces[surf_name] = {'*type':surf_type,
                                           '*definition':[]
                                           }
        
        doseek = True
        while doseek:
            newseek = fid.tell()
            lntxt = fid.readline().lower().split('\n')[0].replace(' ','')
            lnarr = [ss for ss in lntxt.split(',') if ss.__len__()>0]
            seekpos = fid.tell()
                
            if not seekpos>newseek:
                doseek = False
                
            elif lnarr.__len__()>0:
                # either node numbers or node sets
                
                if surf_type == 'node':
                    for nds in lnarr:
                        if nds in self.nd_set.keys():
                            nids = self.nd_set[nds]
                            
                        else:
                            try:
                                nids = [int(nds)]
                            except:
                                nids = []
                                doseek = False
                                pass
                            
                        for ii in nids:
                            surf['*definition'] += [ii]
                            
                            
                elif (lnarr.__len__()==2)&(lnarr[1][0]=='s'):
                    try:
                        int(lnarr[0])
                        surf['*definition'] += [lnarr[0]+'--'+lnarr[1]]
                    except:
                        doseek=False
                
                else:
                    doseek=False
                    
                    
                
            
            else:
                ## check if new keyword activated:
                keyword = self.__contains_keyword(lntxt)
            
                if not keyword is None:
                    doseek = False
                    
        
        return newseek
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ##
     # ADD AN INITIAL CONDITION
            
    def __add_ini(self,fid,seekpos=0):
        '''
        
        add initial conditions to a nodes or elements
        
            self.ISVs [self.time_format%0] [el_id] [gp_nr] [type] = gp_val
            self.NodeVal [self.time_format%0] [type] [nd_id] = nd_val
            
            
            
            ALLOWABLE TYPES
            
            for ISV values:
             >> PLASTIC STRAIN  ( el_id, gauss pt, EP11, EP22, etc ) 
             >> STRESS          ( el_id, gauss pt, S11, S22, S33, S12, S13, S23 )
             >> SOLUTION (ISVS) ( el_id, gauss pt, ISV1, ISV2, ... , ISVn ) 
            
            
            for Nodal values:
             >> TEMPERATURE     ( 2 values per row =   nd_id (or group), temperature value )
             >> DISPLACEMENT    ( 3 values per row =   nd_id (or group), dof value (1, 2, 3), displacement value )
        
        '''  
        print("--> ADDING INITIAL CONDITIONS")
        
        tic = time()
        fid.seek(seekpos)
        lntxt = fid.readline().lower().split('\n')[0] # this line needs to contain the keyword "*initial conditions"
        newseek = fid.tell()
        
        # the type dictates what to do:
        ini_types = ['plastic strain',
                     'stress',
                     'solution',
                     'temperature',
                     'displacement']
        
               
        
        try:
            ini_tp = lntxt.split('type')[1].split('=')[1].split(',')[0]
        except:
            ini_tp = 'unknown'
        
        
        if not "*initial conditions" in lntxt:
            print("WRONG ASSIGNMENT - no << *INITIAL CONDITIONS >> flag in expected input file line")
            return newseek
        
        if not ini_tp in ini_types:
            print(" ** WARNING: INITIAL CONDIDTIONS REQUIRE A VALID TYPE")
            return newseek
        
        
        #
        # 
        #
        # assign initial conditions (temperature)
        if ini_tp == 'temperature':
            NodeVals = self.NodeVals[self.time_format%0]
            if not 'temperature' in NodeVals:
                NodeVals['temperature'] = {}
            # following lines are either node or node_sets:
            doseek = True
            while doseek:
                seekpos = fid.tell()
                lntxt = fid.readline().lower()
                newseek = fid.tell()
            
                ln_arr = lntxt.replace(" ","").split(',')
                
                try:
                    try:
                        nd_set = [int(ln_arr[0])]
                    except:
                        nd_set = self.nd_set[ln_arr[0]]
                        
                    for nd_id in nd_set:
                        NodeVals['temperature'][nd_id] = float(ln_arr[1])
                    
                except:
                    doseek = False
                
            
        return newseek
        
        
        
        
        
        
        
        
        
# #  ADD A STEP

        
    def __add_step(self,fid,seekpos=0):
        '''
        
        add lists of elements to existing sets or create a new set:
        
            self.sections[elset] = {}
        
        '''  
        print("--> ADDING A STEP DEFINITION")
        tic = time()
        fid.seek(seekpos)
        lntxt = fid.readline().lower() # this line needs to contain the keyword "*step" 
        newseek = fid.tell()
        
        
        if not "*step" in lntxt:
            print("WRONG ASSIGNMENT - no << *STEP >> flag in expected input file line")
            return newseek
        
        # try the inc value:
        inc=1000
        #print(lntxt)
        if 'inc' in lntxt:
            try:
                inc = int(lntxt.split('inc')[1].split('=')[1].split(',')[0])
            except:
                pass
            
        nlgeom = False
        if 'nlgeom' in lntxt:
            nlgeom = True
            
        #
        self._nseps += 1
        cstep = self.steps[self._nseps] = {'*inc':inc,
                                           '*nlgeom':nlgeom,
                                           '*direct':False,
                                           '*explicit':False,
                                           '*steady state':False
                                           }
            
        #
        # next line should be step type
        step_types = ['*static',
                      '*dynamic',
                      '*heat transfer',
                      '*coupled temperature-displacement',
                      '*uncoupled temperature-displacement']
        
        # check step type
        lntxt = fid.readline().lower() # this line needs to contain the keyword "*step" 
        if '*dynamic' in lntxt:
            print(' >> STEP TYPE   =   << DYNAMIC >>')
            cstep['*type'] = '*dynamic'
        elif '*heat transfer' in lntxt:
            print(' >> STEP TYPE   =   << HEAT TRANSFER >>')
            cstep['*type'] = '*heat transfer'
        elif '*coupled temperature-displacement' in lntxt:
            print(' >> STEP TYPE   =   << COUPLED TEMPERATURE DISPLACEMENT >>')
            cstep['*type'] = '*coupled temperature-displacement'
        elif '*uncoupled temperature-displacement' in lntxt:
            print(' >> STEP TYPE   =   << UNCOUPLED TEMPERATURE DISPLACEMENT >>')
            cstep['*type'] = '*uncoupled temperature-displacement'
        else:
            print(' >> STEP TYPE   =   << STATIC >>')
            cstep['*type'] = '*static'
            
        
        if 'direct' in lntxt: cstep['*direct'] = True
        if 'explicit' in lntxt: cstep['*explicit'] = True
        if 'steady state' in lntxt: cstep['*steady state'] = True
        
        
        time_arr = [1.,1.,1.,1.]
        try:
            lntxt = fid.readline().lower()
            ln_arr = [float(lx) for lx in lntxt.split(',')]
            for i in range(ln_arr.__len__()):
                time_arr[i] = ln_arr[i]
        except:
            pass
        
        cstep['*time'] = time_arr
        
        #
        doseek = True
        while doseek:            
            seekpos = fid.tell()
            lntxt = fid.readline().lower().split('\n')[0]
            newseek = fid.tell()
            if not newseek>seekpos:
                return newseek
            if '*end step' in lntxt:
                return newseek
            
            #print('SEEKING')
            #print(lntxt)
            
            if '*boundary' in lntxt:
                if not '*boundary' in cstep.keys():
                    cstep['*boundary'] = {}
                
                ampname = None
                if 'amplitude' in lntxt:
                    ampname = lntxt.replace(" ","").split('amplitude=')[1].split(',')[0]
                    
                bc_nds = []
                bc_dof = []
                bc_vals = []
                
                do_key = True
                while do_key:
                    seekpos = fid.tell()
                    lntxt = fid.readline().lower().split('\n')[0].replace(" ","")
                    ln_arr = lntxt.split(',')
                    #print(ln_arr)
                    try:
                        try:
                            nd_set = [int(ln_arr[0])]
                        except:
                            nd_set = self.nd_set[ln_arr[0]]
                                                        
                        for nd_ID in nd_set:
                            #print(nd_ID)
                            #print(int(ln_arr[1]))
                            for dof_ID in range(int(ln_arr[1]),int(ln_arr[2])+1):
                                if nd_ID in cstep['*boundary'].keys():
                                    cstep['*boundary'][nd_ID][dof_ID] = {'value':float(ln_arr[3]),'amplitude':ampname}
                                else:
                                    cstep['*boundary'][nd_ID] = {dof_ID:{'value':float(ln_arr[3]),'amplitude':ampname}}
                    except:
                        do_key = False
                        fid.seek(seekpos)
                    
                
                    
                #
                
            
            
            if '*cload' in lntxt:
                if not '*cload' in cstep.keys():
                    cstep['*cload'] = {}
                
                ampname = None
                if 'amplitude' in lntxt:
                    ampname = lntxt.replace(" ","").split('amplitude=')[1].split(',')[0]
                    
                bc_nds = []
                bc_dof = []
                bc_vals = []
                
                do_key = True
                while do_key:
                    seekpos = fid.tell()
                    lntxt = fid.readline().lower().split('\n')[0].replace(" ","")
                    ln_arr = lntxt.split(',')
                    try:
                        try:
                            nd_set = [int(ln_arr[0])]
                        except:
                            nd_set = self.nd_set[ln_arr[0]]
                                                        
                        for nd_ID in nd_set:
                            if nd_ID in cstep['*cload'].keys():
                                cstep['*cload'][nd_ID][int(ln_arr[1])] = {'value':float(ln_arr[2]),'amplitude':ampname}
                            else:
                                cstep['*cload'][nd_ID] = {int(ln_arr[1]):{'value':float(ln_arr[2]),'amplitude':ampname}}
                    except:
                        do_key = False
                        fid.seek(seekpos)
                    
                
                    
            #
            # capture the amplitudes
            if '*amplitude' in lntxt:
                ampname = lntxt.replace(" ","").split('name=')[1].split(',')[0]
                usetime = 'step time'
                if 'total time' in lntxt:
                    usetime = 'total time'
                ampdict = self.amplitudes[ampname] = {'time':usetime,
                                                      'table':[]}
                
                table_ls = []
                
                do_key = True
                while do_key:
                    seekpos = fid.tell()
                    lntxt = fid.readline().lower().split('\n')[0]
                    try:
                        for lx in lntxt.split(','):
                            table_ls += [float(lx)]
                    except:
                        do_key = False
                        fid.seek(seekpos)
                    
                for k in range(table_ls.__len__()//2):
                    ampdict['table'] += [[table_ls[2*k],table_ls[2*k+1]]]
                    
                    
        return seekpos




# #  ADD AN INITIAL BOUNDARY

        
    def __add_boundary(self,fid,seekpos=0):
        '''
        
        add initial boundary conditions for fixed degrees of freedom
        
            self.steps[0] = {'*boundary':{}}
        
        '''  
        print("--> ADDING A FIXED BOUNDARY DEFINITION")
        tic = time()
        fid.seek(seekpos)
        lntxt = fid.readline().lower() # this line needs to contain the keyword "*step" 
        newseek = fid.tell()
        
        
        if not "*boundary" in lntxt:
            print("WRONG ASSIGNMENT - no << *BOUNDARY >> flag in expected input file line")
            return newseek
        
        #
        if not '*boundary' in self.steps[0].keys():
            self.steps[0]['*boundary'] = {}
            
        bdy_dict = self.steps[0]['*boundary']
                
        do_key = True
        while do_key:
            seekpos = fid.tell()
            lntxt = fid.readline().lower().split('\n')[0].replace(" ","")
            ln_arr = lntxt.split(',')
            #print(ln_arr)
            try:
                try:
                    nd_set = [int(ln_arr[0])]
                except:
                    nd_set = self.nd_set[ln_arr[0]]
                                                
                for nd_ID in nd_set:
                    try:
                        for dof_ID in range(int(ln_arr[1]),int(ln_arr[2])+1):
                            if nd_ID in bdy_dict.keys():
                                bdy_dict[nd_ID][dof_ID] = {'value':0.,'fixed':True}
                            else:
                                bdy_dict[nd_ID] = {dof_ID:{'value':0.,'fixed':True}}
                    except:
                        nd_ID = int(ln_arr[1])
                        if nd_ID in bdy_dict.keys():
                            bdy_dict[nd_ID][dof_ID] = {'value':0.,'fixed':True}
                        else:
                            bdy_dict[nd_ID] = {dof_ID:{'value':0.,'fixed':True}}
            except:
                do_key = False
                fid.seek(seekpos)         
                
            
                    
        return seekpos
        
        
        
        
        
##  ADD A SECTION DEFINITION
        
    #def __add_section(self,fid,seekpos=0):
        #'''
        
        #add lists of elements to existing sets or create a new set:
        
            #self.sections[elset] = {}
        
        #'''  
        #print("--> ADDING A SECTION DEFINITION")
        #tic = time()
        #fid.seek(seekpos)
        #lntxt = fid.readline().lower() # this line needs to contain the keyword "*material" and a "NAME" flag
        #newseek = fid.tell()
        
        
        #if not "*material" in lntxt:
            #print("WRONG ASSIGNMENT - no << *MATERIAL >> flag in expected input file line")
            #return newseek
        
        
        ## get name for "material" in the text
        #try:
            #matname = [nn for nn in lntxt.split('name')[1].split('=')[1].split('\n')[0].split(',')[0].split(' ') if nn.__len__()>0][0]
        #except:
            #print("WRONG ASSIGNMENT - << *MATERIAL >> requires a << NAME >> flag ")
            #return newseek
        
        
        ###check if mat names are already defined (else create empty list):
        #if not matname in self.materials.keys():
            #self.materials[matname] = {}
                
        ### while loop to add material info to the material dictionary
        #doseek = True
        #while doseek:
            #seekpos = fid.tell()
            #lntxt = fid.readline().lower()
            #newseek = fid.tell()
            #if not newseek>seekpos:
                #return newseek
            
            
            ###
            ### try to see if linetext conforms to list of id's or single node ID
            ##lnarr = lntxt.split('\n')[0].split(',')
            
            ### try to add node to dictionary
            ##try:
                ##for kn in lnarr:
                    ##nd_ID = int(kn)
                    ##for mat in matnames:
                        ##self.materials[mat] += [nd_ID]
                
            ##except:
                ### check if new keyword activated:
                ##keyword = self.__contains_keyword(lntxt)
                
                ##if not keyword is None:
                    ##print('Assigned element set %s'%mat)
                    ##return seekpos
                    
        #return seekpos
