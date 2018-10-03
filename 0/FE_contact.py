#
#
#  Contact residual and stiffness contribution to global degrees of freedom
#
#


#def map_to_isocoords(X,El_Coords):
    ### map global coordinates [x,y] to local "isoparamnetric coordinates" [xi, eta]
    
    #if El_Coords.shape[0]==4: # four node isoparametrix element
        
        
    #shapefn = [.25*(1-xi)*(1-eta),
               #.25*(1+xi)*(1-eta),
               #.25*(1+xi)*(1+eta),
               #.25*(1-xi)*(1+eta)]
    
    
def penalty_residual(Coord,Segment,Interaction,only_resid=True):
    
    from numpy import array, r_, zeros
    
    Cinfo = penalty_contact(Coord,Segment,Interaction)
    
    residual = r_[Cinfo[0],Cinfo[1].flatten()]
    
    
    #        
    # do finite difference 
    U0 = r_[Coord,Segment.flatten()]
    nr_dofs = U0.size
    
    if only_resid:
        return residual, zeros((nr_dofs,nr_dofs))
    
    if all(abs(residual)<1e-8):
        return residual, zeros((nr_dofs,nr_dofs))
    
    
    dRdU = []
    
    pert = 1e-8
    for k in range(nr_dofs):
        U = U0.copy()
        U[k] += pert
        resid2 = penalty_residual(U[:2],U[2:].reshape(-1,2),Interaction,True)[0]
        dRdU += [(resid2-residual)/pert]
        
    #
    # symmetric stiffness?:
    dRdU = 0.5*array(dRdU)
    stiff = dRdU.T+dRdU.T
    
    return residual, stiff
        
    
    
     
     
    
    
def penalty_contact(Coord,Segment,Interaction):
    '''
    take:
        Coord = [x,y]
    and
        Segment = [[x1, y1],
                   [x2, y2]] (if straight vs 3rd point if curved 
                    NOTE (line segment of quadratic element edge ie 3rd node in the middle of isoparamnetric definiton)
                   
    as well as the penalty contact iteraction dictionary and return contact residual 
    
    ALSO return contact senistivity w.r.t. Coord locations and Segment pointr locations for inclusion in global
           sensitivity matrix
    
    '''
    
    from numpy import array, sqrt, sum, exp, arctan, pi
    
    xy = None
    
    Fnd = 0*Coord
    Fdistr = 0*Segment

    if Segment.shape[0]==2: 
        # edge of linear element
        
        bv = 0.5*(-Segment[0]+Segment[1])
        cv = 0.5*(Segment[0]+Segment[1])
        
        xy     = lambda s: bv*s + cv
        dxy_ds = lambda s: bv
        
        shapefns = lambda s: 0.5*array([ 1-s , 1+s ])
        
    elif Segment.shape[0]==3:
        
        # edge of quadratic element
        av = 0.5*(Segment[0]+Segment[1])-Segment[2]
        bv = 0.5*(-Segment[0]+Segment[1])
        cv = Segment[2]
        
        xy     = lambda s: av*s**2 + bv*s + cv
        dxy_ds = lambda s: 2*av*s + bv
        
        shapefns = lambda s: array([ -0.5*s*(1-s) , 0.5*s*(1+s), 1-s**2 ])
        
        

    
    #
    if callable(xy):
        
        objf = lambda s : sqrt(sum((Coord-xy(s))**2))
        
        from scipy.optimize import golden
        
        sopt = golden(objf)
        
        if abs(sopt)>1:
            return Fnd,Fdistr
            
        xopt = xy(sopt)
        dist = xopt - Coord
            
        tang = dxy_ds(sopt)
        tang /= sqrt(sum(tang**2))
        
        # area
            
        norm = array([tang[1], -tang[0]])
            
        delta = sum(norm*dist)
                
        if True:#delta>-1e-3:
        
            # nodal force:
            F_pars = Interaction['params']
            if F_pars[0] is None: F_pars[0]=1e6
            if F_pars[1] is None: F_pars[1]=1e-3
            
            if Interaction['pressure-overclosure'] == 'linear':
                Fd = lambda d : F_pars[0]**2*d*(0.5 + arctan(d/F_pars[1])/pi)
            else:
                Fd = lambda d: F_pars[0]*exp(F_pars[1]*d)
                
            Fnd = Fd(delta)*norm
            # distributed force on 
            Fdistr = array([-w*Fnd for w in shapefns(sopt)])
            
            
    # return None if not a surface contact penalty
    return Fnd,Fdistr
        














#def lagrange_contact(Coord,Segment):
    
    #from numpy import array, r_, zeros
    
    #Cinfo = lagrange_contact_weights(Coord,Segment,Interaction)
    
    #residual = r_[Cinfo[0],Cinfo[1].flatten()]
    
    
    ##        
    ## do finite difference 
    #U0 = r_[Coord,Segment.flatten()]
    #nr_dofs = U0.size
    
    #if only_resid:
        #return -residual, zeros((nr_dofs,nr_dofs))
    
    #if all(abs(residual)<1e-8):
        #return -residual, zeros((nr_dofs,nr_dofs))
    
    
    ##dRdU = []
    
    ##pert = 1e-8
    ##for k in range(nr_dofs):
        ##U = U0.copy()
        ##U[k] += pert
        ##resid2 = penalty_residual(U[:2],U[2:].reshape(-1,2),Interaction,True)[0]
        ##dRdU += [(resid2-residual)/pert]
        
    ###
    ### symmetric stiffness?:
    ##dRdU = 0.5*array(dRdU)
    ##stiff = dRdU.T+dRdU.T
    
    #return -residual, stiff







    
def lagrange_contact_weights(Coord,Segment):
    '''
    take:
        Coord = [x,y]
    and
        Segment = [[x1, y1],
                   [x2, y2]] (if straight vs 3rd point if curved 
                    NOTE (line segment of quadratic element edge ie 3rd node in the middle of isoparamnetric definiton)
                   
    as well as the penalty contact iteraction dictionary and return contact residual 
    
    ALSO return contact senistivity w.r.t. Coord locations and Segment pointr locations for inclusion in global
           sensitivity matrix
    
    '''
    
    from numpy import array, sqrt, sum, exp, arctan, pi
    
    xy = None
    
    Fnd = 0*Coord
    Fdistr = 0*Segment

    if Segment.shape[0]==2: 
        # edge of linear element
        
        bv = 0.5*(-Segment[0]+Segment[1])
        cv = 0.5*(Segment[0]+Segment[1])
        
        xy     = lambda s: bv*s + cv
        dxy_ds = lambda s: bv
        
        shapefns = lambda s: 0.5*array([ 1-s , 1+s ])
        
    elif Segment.shape[0]==3:
        
        # edge of quadratic element
        av = 0.5*(Segment[0]+Segment[1])-Segment[2]
        bv = 0.5*(-Segment[0]+Segment[1])
        cv = Segment[2]
        
        xy     = lambda s: av*s**2 + bv*s + cv
        dxy_ds = lambda s: 2*av*s + bv
        
        shapefns = lambda s: array([ -0.5*s*(1-s) , 0.5*s*(1+s), 1-s**2 ])
        
        

    
    #
    if callable(xy):
        
        objf = lambda s : sqrt(sum((Coord-xy(s))**2))
        
        from scipy.optimize import golden
        
        sopt = golden(objf)
        
        if abs(sopt)>1:
            return Fnd,Fdistr
            
        xopt = xy(sopt)
        dist = xopt - Coord
            
        tang = dxy_ds(sopt)
        tang /= sqrt(sum(tang**2))
        
        # area
            
        norm = array([tang[1], -tang[0]])
            
        delta = sum(norm*dist)
                
        if delta>1e-8:
                        
            Fnd = norm
            # distributed force on 
            Fdistr = array([-w*Fnd for w in shapefns(sopt)])
            
            
    # return None if not a surface contact penalty
    return Fnd,Fdistr








       
if __name__ == '__main__':
        
    import numpy as np
    #

    from numpy import array, sqrt, sum, arctan, pi
    
    import matplotlib.pylab as pl
    
    pl.ion()
    
    pl.close('all')
    
    Coord  = array([0.15,0.23])
    Coord  = array([-1.1,-0.3])
    
    Segment = array([[-0.1,-0.2],
                       [0.5,0.5],
                       [0.1,0.15]])
    
    
    Interaction = {'params': [1000000000000.0, None, None], 'pressure-overclosure': 'linear'}
    
    
    av = 0.5*(Segment[0]+Segment[1])-Segment[2]
    bv = 0.5*(-Segment[0]+Segment[1])
    cv = Segment[2]
    xy     = lambda s: av*s**2 + bv*s + cv
    dxy_ds = lambda s: 2*av*s + bv
    
    shapefns = lambda s: array([ -0.5*s*(1-s) , 0.5*s*(1+s), 1-s**2 ])
    
    objf = lambda s : sqrt(sum((Coord-xy(s))**2))
            
    from scipy.optimize import golden
    
    sopt = golden(objf)
        
    xopt = xy(sopt)
    dist = xopt - Coord
        
    tang = dxy_ds(sopt)
    tang /= sqrt(sum(tang**2))
    
    # area
        
    norm = array([tang[1], -tang[0]])
        
    delta = sum(norm*dist)
    F_pars = Interaction['params']
    if F_pars[0] is None: F_pars[0]=1e6
    if F_pars[1] is None: F_pars[1]=1e-3
    Fd = lambda d : F_pars[0]*d*(0.5 + arctan(d/F_pars[1]/pi))
    
    #Fnd = Fd(delta)*norm
    Fnd = norm
    Fdistr = array([ -w*Fnd for w in shapefns(sopt)])
    
    
    mins = min(sopt,-1)
    maxs = max(sopt,1)
    
    line = array([xy(s) for s in np.linspace(mins,maxs,101)])
    
    
    fig = pl.figure(1)
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    ax.plot(line[:,0],line[:,1],'k')
    ax.plot(Segment[:,0],Segment[:,1],'kd')
    ax.scatter(Coord[0],Coord[1],marker='o',color='b')
    ax.scatter(xopt[0],xopt[1],marker='x',color='r')
    # tangent & normal:
    t_pts = np.c_[xopt,xopt+tang].T
    n_pts = np.c_[xopt,xopt+norm].T
    ax.plot(t_pts[:,0],t_pts[:,1],'b',lw=2)
    ax.plot(n_pts[:,0],n_pts[:,1],'r',lw=2)
    
    
    





    
#def assemble(U0,dU,contact_pair,):
    ### assemble the contact residual and stiffness contribution and return for use in a sparse matrix matrix assembly
    ## I.E. : residual,residual dof, tangent column, tangent rows, tangent stiffness value
    
    
    
    
    #el_type = self._el_type[el_nr]
    ## element connectivity
    #el_connect = self.elements[el_nr]
    ##
    ## element nodal coordinates
    #X = [self.nodes[node_id][0] for node_id in el_connect]
    #Y = [self.nodes[node_id][1] for node_id in el_connect]
    #XY   = c_[X, Y];

    
    
    #return resid,dof,col,row,stiff






    
    
    
    
    
#def residual():
    ### assemble the contact residual and relevant degrees of freedom of the residual contribution
    ##
    

    
    
    #return resid,dof
