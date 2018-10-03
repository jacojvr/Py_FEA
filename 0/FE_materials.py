
from numpy import array, matrix, zeros, sort, sqrt


#
# get all the material properties for a given temperature:
def interpolate_material_properties(T_gauss,material_info):
    
    T_gauss = array(T_gauss).flatten()[0]
    #
    other_keys = ['type','zero','hardening']
    matprops = {}
    for kn in material_info.keys():
        prop_dict = material_info[kn]
        prop_temps =  [float(k) for k in prop_dict.keys() if not k in other_keys]
        # if only one temperature -- >
        if prop_temps.__len__() == 1:
            matprops[kn] = prop_dict[prop_temps[0]]
          
        else:
            # if multiple temperatures:
            tmin = min(prop_temps)
            tmax = max(prop_temps)
            if T_gauss<=tmin:
                matprops[kn] = prop_dict[tmin]     
            elif T_gauss>=tmax:
                matprops[kn] = prop_dict[tmax]
            else:
                try:
                    # move tmin and tmax to the two closest (below and above)
                    for T in array(prop_temps):
                        if T>tmin and T<=T_gauss: tmin = T
                        if T<tmax and T>T_gauss: tmax = T
                    #
                    alpha = (T_gauss-tmin)/(tmax-tmin)                    
                    provals = alpha*array(prop_dict[tmax]) + (1-alpha)*array(prop_dict[tmin])
                    
                    matprops[kn] = list(provals)
                                        
                except:
                    print("WARNING: Problem interpolating << %s >> for Temperature << %.4f >>"%(kn,T_gauss))
                    print("\tplease check that you define the same number of data points per Temperature")
        
        
        for k2 in other_keys:
            if k2 in prop_dict.keys():
                matprops[kn+'_'+k2] = prop_dict[k2]
                
        #elif prop_temps.__len__() == 0:
            #matprops[kn]
            
    return matprops



def mechanical_response(stress_p,strain_p,dstrain,T_gauss,dt,material_info,isvP,nlgeom=False,only_resid=False):
    #
    mat_type = None
    #
    if '*elastic' in material_info.keys():
        mat_type = 'elastic'
        # add possible plastic deformation
        if '*plastic' in material_info.keys():
            mat_type = 'plastic'
            
    
    if mat_type == 'elastic':
        return simple_elastic(stress_p,strain_p,dstrain,T_gauss,dt,material_info,isvP,only_resid=only_resid)
    
    
    if mat_type == 'plastic':
        if nlgeom:
            return multiplicative_mises_plastic(stress_p,strain_p,dstrain,T_gauss,dt,material_info,isvP,only_resid=only_resid)
        else:
            return additive_mises_plastic(stress_p,strain_p,dstrain,T_gauss,dt,material_info,isvP,only_resid=only_resid)
            
    
    # add plastic deformation
    
    
    else:
        print('!!! IN MATERIAL << incremental strain >> NEEDS A VALID SET OF MATERIAL PARAMETERS FOR DEFORMATION RESIDUAL')
        return stress_p, None, isvP
    







def simple_elastic(stress_p,strain_p,dstrain,T_gauss,dt,material_info,isvP,only_resid=False):
    # 
    #print('simple elastic')
    matprops = interpolate_material_properties(T_gauss,material_info)
    
    # strain = [11,22,33,12]
    #
    #
    [E,nu] = matprops['*elastic']
    # Various elastic parameter values & lame params
    xk = E/(1-2*nu)
    eg2 = E/(1+nu)
    al = (xk-eg2)/3
    # bulk and shear modulus
    bulk_mod = xk/3
    eg = eg2/2
    #
    # d(stress)/d(strain)
    ddsdde = matrix([[al+eg2, al, al, 0.],
                        [al, al+eg2, al, 0.],
                        [al, al, al+eg2, 0.],
                        [0., 0., 0., eg]])
    
    strain = strain_p+dstrain
    stress = ddsdde*matrix(strain).reshape(-1,1)

    #ISN = {'stress':,
           #'strain':strain}
    
    return stress, ddsdde, {}
    
    
    





   
def additive_mises_plastic(stress_p,strain_in,dstrain,T_gauss,dt,material_info,isvP,only_resid=False):
    # 
    matprops = interpolate_material_properties(T_gauss,material_info)
    
    # strain = [11,22,33,12+21]
    #
    #
    [E,nu] = matprops['*elastic']
    # Various elastic parameter values & lame params
    xk = E/(1-2*nu)
    eg2 = E/(1+nu)
    al = (xk-eg2)/3
    # bulk and shear modulus
    bulk_mod = xk/3
    eg = eg2/2
    #
    # tolerance:
    toler = 1e-8
    # 2/3
    k23 = 2./3
    k32 = 3./2
    # sqrt(2./3)
    sq23 = sqrt(k23)
    #
    # Identity vector
    Ivec = matrix([[1],[1],[1],[0]])
    #
    # get previous plastic strains:
    try:
        strain_pl = matrix(isvP['plastic strain']).reshape(-1,1)
        eq_pl = isvP['equivalent plastic'][0]
    except:
        strain_pl = zeros((4,1))
        eq_pl = 0
    #
    #
    # elastic strain component
    strain_in = matrix(strain_in+dstrain).reshape(-1,1)
    strain = strain_in-strain_pl
    #
    # d(stress)/d(strain)
    ddsdde = matrix([[al+eg2, al, al, 0.],
                     [al, al+eg2, al, 0.],
                     [al, al, al+eg2, 0.],
                     [0., 0., 0., eg]])
    # stress
    stress = ddsdde*strain
    # von Mises equivalent stress
    smises=(stress[0,0]-stress[1,0])**2+(stress[1,0]-stress[2,0])**2+(stress[2,0]-stress[0,0])**2+6*stress[3,0]**2
    smises=sqrt(smises/2)
    #
    # call hardening function
    deqpl = 0.
    sy, dsy, isvN = isotropic_mises_hardening(eq_pl,deqpl,T_gauss,dt,matprops,isvP)
    #
    # trial yield function:
    ftrial = smises-sy
    #
    # check if necessary and do plastic
    if ftrial>toler:
        # plastic deformation:
        kount,kmax = 0,15 # newton loop initial and max
        eg3 = eg*3
        # newton loop to find new plastic strain
        while abs(ftrial)>toler and kount<kmax:
            kount+=1
            deqpl=abs(deqpl+ftrial/(eg3+dsy))
            sy, dsy, isvN = isotropic_mises_hardening(eq_pl,deqpl,T_gauss,dt,matprops,isvP)
            ftrial=smises-eg3*deqpl-sy
        #
        # update state
        #hydrostatic and deviatoric stress
        shydro = (stress[0,0]+stress[1,0]+stress[2,0])/3.
        stress_dev = stress - shydro*Ivec
        # direction of plastic flow
        flow_direction = stress_dev/smises
        stress = sy*flow_direction + shydro*Ivec
        
        if only_resid:
            return stress, 0, {}
        # update equivalent plastic strain
        eq_pl += deqpl
        gamma = deqpl*k32
        strain_pl = strain_pl + gamma*flow_direction
        #
        # LOCAL STIFFNESS
        effg=eg*sy/smises
        effg2=eg2*sy/smises
        effg3=k32*effg2
        efflam=(xk-effg2)/3
        effhrd=eg3*dsy/(eg3+dsy)-effg3
        
        ddsdde = matrix([[effg2+efflam, efflam, efflam, 0.],
                         [efflam, effg2+efflam, efflam, 0.],
                         [efflam, efflam, effg2+efflam, 0.],
                         [0., 0., 0., effg]])
        
        ddsdde += effhrd*flow_direction*flow_direction.T
    
    # store ISVs        
    isvN['equivalent plastic'] = array(eq_pl).flatten()
    isvN['plastic strain'] = array(strain_pl).flatten()
    isvN['stress'] = array(stress).flatten()
    
    return stress, ddsdde, isvN



import numba



def multiplicative_mises_plastic(stress_in,strain_in,dstrain,T_gauss,dt,material_info,isvP,only_resid=False):
    # 
    matprops = interpolate_material_properties(T_gauss,material_info)
    
    #print('multiplicative')
    
    # strain = [11,22,33,12+21]
    #
    #
    #  HARDCODED FOR PLANE STRAIN
    #
    #  e33=0 , i.e. C33 = 1 and inv(C)_33 = 1 etc.
    #
    #
    #
    #
    [E,nu] = matprops['*elastic']
    #if '*plastic' in matprops.keys()
        #plast =  array(matprops['*plastic']).reshape(-1,2)
    
    # Various elastic parameter values & lame params
    xk = E/(1-2*nu)               # xk
    eg2 = E/(1+nu)
    al = (xk-eg2)/3         # al
    # bulk and shear modulus
    bulk_mod = xk/3
    eg = eg2/2             # um
    #
    # tolerance:
    toler = 1e-8
    # 2/3
    k23 = 2./3
    # sqrt(2./3)
    sq23 = sqrt(k23)
    # 4/3
    k43 = 4./3
    
    #
    #
    #
    # get previous plastic strains:
    try:
        strain_pl = matrix(isvP['plastic strain']).reshape(-1,1)
        eq_pl = isvP['equivalent plastic'][0]
    except:
        strain_pl = zeros((4,1))
        eq_pl = 0
        
    # convert from Lagrange strains to right Cauchy Green:
    Ivec = matrix([[1],[1],[1],[0]])
    
    EtoC = matrix([[2,0,0,0],
                   [0,2,0,0],
                   [0,0,2,0],
                   [0,0,0,1]])
    
    C_pl = -EtoC*strain_pl + Ivec
    strain = strain_in+dstrain
    C = EtoC*strain + Ivec
    
    #C_pl = -2*strain_pl + Ivec
    #C_pl[3,0] *= 0.5  # given  strain[4,0] = E_12+E_21 = 0.5*(C_12+C_21)
    #C = 2*strain + Ivec
    #C[3,0] *= 0.5 # given strain[4,0] = E_12+E_21 = 0.5*(C_12+C_21)
    # determinant of Cauchy strain and deformation gradient:
    detC1 = C[0,0]*C[1,0] - C[3,0]*C[3,0] # if C_33 = 1
    vj2 = C[2,0]*detC1     # vj2
    vj = sqrt(vj2)       # vj
    # inverse right Cauchy strain:
    C_inv = matrix([[ C[1,0]*C[2,0] ],
                    [ C[0,0]*C[2,0] ],
                    [ detC1 ],
                    [ -C[3,0]*C[2,0]]] )/vj2
    #
    # for isochoric quantities:
    vj23 = vj**(2./3)  #vj23
    Cpl_bar = C_pl/vj23  
    Cpb_trace = (C[0,0]*Cpl_bar[0,0] +  #umb
                 C[1,0]*Cpl_bar[1,0] +
                 C[2,0]*Cpl_bar[2,0] + 
                 2*C[3,0]*Cpl_bar[3,0])/3.  
    Cpb_dev = Cpl_bar - Cpb_trace*C_inv #cplb
    #
    # deviatoric trial stress
    stress_dev = eg*Cpb_dev  #stril
    smises = sqrt(stress_dev[0,0]*stress_dev[0,0] +
                     stress_dev[1,0]*stress_dev[1,0] +
                     stress_dev[2,0]*stress_dev[2,0] +
                     2*stress_dev[3,0]*stress_dev[3,0]) 
    
    # call hardening function
    deqpl = 0.
    stress_y, d_stress_y, isvN = isotropic_mises_hardening(eq_pl,deqpl,T_gauss,dt,matprops,isvP)
    # convert yield stess (Cauchy) to Kirchhoff values:
    #  i.e. kirchhoff yield = yield stress * det (F)
    fiso0,dfiso0 = stress_y*vj, d_stress_y*vj # initial values
    sy = fiso0
    # trial yield function:
    ftrial = smises - fiso0*sq23
    #print(" TRIAL FLOW = ",ftrial)
    # if elastic: ftrial ~ 0
    if ftrial<toler:
        
        k0 = xk*vj2
        k1 = xk*(vj2-1)
        k2 = k1/2
        
        stress = k2*C_inv + stress_dev
        
        
        if only_resid:
            return stress, 0, {}
        
        #isvN['cauchy'] = array(stress).flatten(),
        #isvN['stress'] = array(stress).flatten()
    
        # LOCAL STIFFNESS MATRIX
        umb = Cpb_trace*eg
        umb2 = 2*umb
        
        xg1 = C_inv[0,0]
        xg2 = C_inv[1,0]
        xg3 = C_inv[2,0]
        xg4 = C_inv[3,0]
        
        xs1 = stress_dev[0,0]
        xs2 = stress_dev[1,0]
        xs3 = stress_dev[2,0]
        xs4 = stress_dev[3,0]
        
        #
        # d(stress)/d(strain)
        ddsdde = zeros((4,4))
        
        ddsdde[0,0] = k43*(umb*xg1*xg1-xs1*xg1)+k0*xg1*xg1-k1*xg1*xg1
        ddsdde[1,1] = k43*(umb*xg2*xg2-xs2*xg2)+k0*xg2*xg2-k1*xg2*xg2
        ddsdde[2,2] = k43*(umb*xg3*xg3-xs3*xg3)+k0*xg3*xg3-k1*xg3*xg3
        ddsdde[3,3] = umb*(xg1*xg2+xg4*xg4-2*xg4*xg4/3)-k43*xs4*xg4+k0*xg4*xg4 -k1*(xg1*xg2+xg4*xg4)/2
        ddsdde[0,1]=ddsdde[1,0] = umb2*(xg4*xg4-xg1*xg2/3) -2*(xs1*xg2+xg1*xs2)/3+k0*xg1*xg2-k1*xg4*xg4
        ddsdde[0,2]=ddsdde[2,0] = -umb2*xg1*xg3/3-2*(xs1*xg3+xg1*xs3)/3+k0*xg1*xg3
        ddsdde[0,3]=ddsdde[3,0] = k23*umb2*xg4*xg1-2*(xs1*xg4+xg1*xs4)/3+k0*xg1*xg4-k1*xg1*xg4
        ddsdde[1,2]=ddsdde[2,1] = -umb2*xg2*xg3/3-2*(xs2*xg3+xg2*xs3)/3+k0*xg2*xg3
        ddsdde[1,3]=ddsdde[3,1] = k23*umb2*xg2*xg4-2*(xs2*xg4+xg2*xs4)/3+k0*xg2*xg4-k1*(xg4*xg2+xg2*xg4)/2
        ddsdde[2,3]=ddsdde[3,2] = -umb2*xg3*xg4/3-2*(xs3*xg4+xg3*xs4)/3+k0*xg3*xg4
        
    else:
        #
        #
        # plastic deformation:
        gamma = 0. # rate of plastic deformation consistency parameter
        d_gamma = 1 # only to check and enter loop = thereafter 
        kount,kmax = 0,15 # newton loop initial and max
        #
        umb = eg*Cpb_trace
        umb2 = umb*2
        umb3 = umb*3
        # newton loop to find new plastic strain
        while abs(d_gamma)>toler and kount<kmax:
            kount+=1
            deqpl = sq23*gamma
            stress_y, d_stress_y, isvN = isotropic_mises_hardening(eq_pl,deqpl,T_gauss,dt,matprops,isvP)
            sy,dfiso = stress_y*vj, d_stress_y*vj
            d_gamma = (smises-umb2*gamma-sq23*sy)/(k23*dfiso+umb2)
            gamma = abs(gamma + d_gamma)
            
        # if caught in a loop / non-convegence
        # send back stress = None
        #
        #if not kount<kmax:
            #return None, None, isvP
        #
        # update state
        #
        k0 = xk*vj2
        k1 = xk*(vj2-1)
        k2 = k1/2 #c8
        #
        # direction of plastic flow
        flow_direction = stress_dev/smises
        #
        # update equivalent plastic strain
        eq_pl += deqpl
        # update right plastic cauchy strain
        k3 = 2*Cpb_trace*vj23*gamma
        C_pl = C_pl - k3*flow_direction
        # Plastic Lagrange strain
        strain_pl = -0.5*(C_pl - Ivec)
        strain_pl[3,0] *= 2  # given that strain[4,0] = E_12+E21 = 0.5*(C_12+C_21)
        k4 = umb2*gamma/smises  #c5
        stress = k2*C_inv + (1-k4)*stress_dev
        #
        
        if only_resid:
            return stress, 0, {}
        #
        # LOCAL STIFFNESS MATRIX
        #
        cc1 = C[0,0]
        cc2 = C[1,0]
        cc3 = C[2,0]
        cc4 = C[3,0]
        #
        xg1 = C_inv[0,0]
        xg2 = C_inv[1,0]
        xg3 = C_inv[2,0]
        xg4 = C_inv[3,0]
        
        xs1 = stress_dev[0,0]
        xs2 = stress_dev[1,0]
        xs3 = stress_dev[2,0]
        xs4 = stress_dev[3,0]
        #
        f0=umb2*gamma/smises
        d0=1+dfiso/umb3

        f1=1/d0-f0
        d1=2*f1*umb-(1/d0-1)*4*gamma*smises/3
        d2=2*smises*f1
        
        xn1=flow_direction[0,0]
        xn2=flow_direction[1,0]
        xn3=flow_direction[2,0]
        xn4=flow_direction[3,0]
        
        xd1=xn1*xn1*cc1+xn1*xn4*cc4+xn4*xn1*cc4+xn4*xn4*cc2
        xd2=xn4*xn4*cc1+xn4*xn2*cc4+xn2*xn4*cc4+xn2*xn2*cc2
        xd3=xn3*xn3*cc3
        xd4=xn1*xn4*cc1+xn1*xn2*cc4+xn4*xn4*cc4+xn4*xn2*cc2
        
        # trace and deviatoric part
        ctr = (xd1*cc1+xd2*cc2+xd3*cc3+2*xd4*cc4)/3
        # 
        xd1=xd1-ctr*xg1
        xd2=xd2-ctr*xg2
        xd3=xd3-ctr*xg3
        xd4=xd4-ctr*xg4
        
       
        elas1= ((umb-f0*umb)*(xg1*xg1+xg1*xg1-
                2*xg1*xg1/3)
                -2*(xs1*xg1+xg1*xs1)/3
                +f0*2*(xs1*xg1+xg1*xs1)/3
                -d1*xn1*xn1-d2*(xn1*xd1+
                xd1*xn1)/2+xk*vj2*xg1*xg1
                -xk*(vj2-1)*(xg1*xg1+xg1*xg1)/2)
        elas2=((umb-f0*umb)*(xg4*xg4+xg4*xg4-
                2*xg1*xg2/3)
                -2*(xs1*xg2+xg1*xs2)/3
                +f0*2*(xs1*xg2+xg1*xs2)/3
                -d1*xn1*xn2-d2*(xn1*xd2+
                xd1*xn2)/2+xk*vj2*xg1*xg2
                -xk*(vj2-1)*(xg4*xg4+xg4*xg4)/2)
        elas3=((umb-f0*umb)*(xg2*xg2+xg2*xg2-
                2*xg2*xg2/3)
                -2*(xs2*xg2+xg2*xs2)/3
                +f0*2*(xs2*xg2+xg2*xs2)/3
                -d1*xn2*xn2-d2*(xn2*xd2+
                xd2*xn2)/2+xk*vj2*xg2*xg2
                -xk*(vj2-1)*(xg2*xg2+xg2*xg2)/2)
        elas4=((umb-f0*umb)*(-2*xg1*xg3/3)
                -2*(xs1*xg3+xg1*xs3)/3
                +f0*2*(xs1*xg3+xg1*xs3)/3
                -d1*xn1*xn3-d2*(xn1*xd3+
                xd1*xn3)/2+xk*vj2*xg1*xg3)
        elas5=((umb-f0*umb)*(-2*xg2*xg3/3)
                -2*(xs2*xg3+xg2*xs3)/3
                +f0*2*(xs2*xg3+xg2*xs3)/3
                -d1*xn2*xn3-d2*(xn2*xd3+
                xd2*xn3)/2+xk*vj2*xg2*xg3)
        elas6=((umb-f0*umb)*(xg3*xg3+xg3*xg3-
                2*xg3*xg3/3)
                -2*(xs3*xg3+xg3*xs3)/3
                +f0*2*(xs3*xg3+xg3*xs3)/3
                -d1*xn3*xn3-d2*(xn3*xd3+
                xd3*xn3)/2+xk*vj2*xg3*xg3
                -xk*(vj2-1)*(xg3*xg3+xg3*xg3)/2)
        elas7=((umb-f0*umb)*(xg1*xg4+xg4*xg1-
                2*xg1*xg4/3)
                -2*(xs1*xg4+xg1*xs4)/3
                +f0*2*(xs1*xg4+xg1*xs4)/3
                -d1*xn1*xn4-d2*(xn1*xd4+
                xd1*xn4)/2+xk*vj2*xg1*xg4
                -xk*(vj2-1)*(xg1*xg4+xg4*xg1)/2)
        elas8=((umb-f0*umb)*(xg4*xg2+xg2*xg4-
                2*xg2*xg4/3)
                -2*(xs2*xg4+xg2*xs4)/3
                +f0*2*(xs2*xg4+xg2*xs4)/3
                -d1*xn2*xn4-d2*(xn2*xd4+
                xd2*xn4)/2+xk*vj2*xg2*xg4
                -xk*(vj2-1)*(xg4*xg2+xg2*xg4)/2)
        elas9=((umb-f0*umb)*(-2*xg3*xg4/3)
                -2*(xs3*xg4+xg3*xs4)/3
                +f0*2*(xs3*xg4+xg3*xs4)/3
                -d1*xn3*xn4-d2*(xn3*xd4+
                xd3*xn4)/2+xk*vj2*xg3*xg4)
        elas10=((umb-f0*umb)*(xg1*xg2+xg4*xg4-
                2*xg4*xg4/3)
                -2*(xs4*xg4+xg4*xs4)/3
                +f0*2*(xs4*xg4+xg4*xs4)/3
                -d1*xn4*xn4-d2*(xn4*xd4+
                xd4*xn4)/2+xk*vj2*xg4*xg4
                -xk*(vj2-1)*(xg1*xg2+xg4*xg4)/2)
        
   
        
        # d(stress)/d(strain)
        ddsdde = zeros((4,4))
         
        ddsdde[0,0] = elas1
        ddsdde[1,1] = elas3
        ddsdde[2,2] = elas6      
        ddsdde[3,3] = elas10
        ddsdde[0,1]=ddsdde[1,0] = elas2
        ddsdde[0,2]=ddsdde[2,0] = elas4
        ddsdde[0,3]=ddsdde[3,0] = elas7
        ddsdde[1,2]=ddsdde[2,1] = elas5
        ddsdde[1,3]=ddsdde[3,1] = elas8
        ddsdde[2,3]=ddsdde[3,2] = elas9
        
        
    # store ISVs        
    isvN['equivalent plastic'] = array(eq_pl).flatten()
    isvN['plastic strain'] = array(strain_pl).flatten()
    isvN['stress'] = array(stress).flatten()
    isvN['strain'] = array(strain).flatten()
    isvN['yield'] = array(sy).flatten()

    return stress, ddsdde, isvN
    
    
    
    
    
    
    
    
def isotropic_mises_hardening(eq_pl,deqpl,T_gauss,dt,matprops,isvP):
    
    # plastic properties:
    props = matprops['*plastic']
    # by default it is now assumed that the properties are tabulated yield vs plastic strain values in a single list,
    # i.e. props = [Sy0, Ep0, Sy1, Ep1, ... , SyN, EpN]
    # ADD AN IF statement here to check other keywords in the matprops dictionary
    
    isvN = isvP.copy() # new ISVs
    stress_y, d_stress_y = isotropic_mises_tabulated(eq_pl+deqpl,props[1::2],props[::2])
    
    return stress_y, d_stress_y, isvN


def isotropic_mises_tabulated(eq_pl,Pvals,Yvals):
    #
    # given an equivalent plastic strain value, return the yield stress and piecewise linear hardening rate:
    # 
    # NOTE perfect plastcity assumed after the final point (and before the first point if ep_0 >0
    #
    stress_y = Yvals[0] #initial yield
    d_stress_y = 0. #perfect plasticity
    #
    epnr = Pvals.__len__()
    for i in range(epnr):
        ep0,sy0 = Pvals[i],Yvals[i]
        if eq_pl>=ep0:
            if i+1<epnr:
                ep1,sy1 = Pvals[i+1],Yvals[i+1]
                d_stress_y = (sy1-sy0)/(ep1-ep0)
                alpha = (eq_pl-ep0)/(ep1-ep0)
                stress_y = sy0 + alpha*d_stress_y
            else:
                stress_y = sy0
                d_stress_y = 0.
    
    return stress_y, d_stress_y

    
