
from numpy import matrix,array,zeros,linalg

import FE_materials




# temperature residual:

def temperature_residual(N,B,Temp0,dTemp,dt,material_info,isvP):
    
    Temp = Temp0+dTemp
    # current temperature:
    T_gauss = N.T*Temp
    # material properties for current temperature
    matprops = FE_materials.interpolate_material_properties(T_gauss,material_info)
    #
    D = matprops['*conductivity']*matrix([[1,0],[0,1]])
    Cp = matprops['*specific heat']
    #
    # NOTE: if dt=0, while dTemp<>0 -->> a steady state heat transfer residual is returned
    dt = 0
    #
    if dt==0:
        # STEADY STATE:
        ### q:
        ## q = -D*B*Temp
        # Residual:
        R_t = B.T*D*(B*Temp)
        # Tangent:
        dR_t = B.T*D*B
    
    
    isvN = {'temp_el':array(T_gauss).flatten()}
    
    return R_t,dR_t,isvN
    



# include a constant velocity in the undeformed configuration

def temperature_residual_velocity(N,B,Temp0,dTemp,velocity,dt,material_info,isvP):
    
    Temp = Temp0+dTemp
    # velocity
    Vel = matrix(velocity).reshape(1,-1)
    # current temperature:
    T_gauss = N.T*Temp
    # material properties for current temperature
    matprops = FE_materials.interpolate_material_properties(T_gauss,material_info)
    #
    D = matprops['*conductivity']*matrix([[1,0],[0,1]])
    Cp = matprops['*specific heat']
    #
    # NOTE: if dt=0, a steady state heat transfer residual is returned
    dt = 0
    #
    #
    # current temperature:
    T_gauss = N.T*Temp
    #
    if dt==0:
        # STEADY STATE:
        ### q:
        ## q = -D*B*Temp
        #
        # rate of temperature change:
        T_dot = Vel*B*Temp
        
        # Residual:
        R_t = B.T*D*(B*Temp) + Cp*N*T_dot
        # Tangent:
        dR_t = B.T*D*B + Cp*N*(Vel*B)
        
        #
        # rate of temperature change 
    
    
    isvN = {'temp_el':array(T_gauss).flatten(),
            'temp_rate':array(T_dot).flatten(),}
    
    return R_t,dR_t,isvN










### DISPLACEMENT RESIDUAL


def small_displacement_2D_residual(N,B,Disp0,dDisp,Temp,dt,material_info,isvP,el_nr=0,el_gp=0,split='none',only_resid=False):
    #
    # split can have different values:
    # if split = 'isochoric', return only the pressure / volumetric component of the residual
    # if split = 'deviatoric', return only the deviatoric component of the residual
    # otherwise, return full residual
    #
    # make sure the displacements are in the correct form
    Disp0 = matrix(Disp0).reshape(-1,1)
    dDisp = matrix(dDisp).reshape(-1,1)
    Temp = matrix(Temp).reshape(-1,1)
    # total displacement
    Disp = Disp0+dDisp
    # current temperature:
    T_gauss = N.T*Temp
    #
    Ivec = matrix([1.,1.,1.,0.]).T
    #
    #
    # Small strain assumption
    B_eps = zeros((4,Disp.size))
    B_eps[0,0::2] = B[0]
    B_eps[1,1::2] = B[1]
    # B[2] = 0  since eps_33 is unaffected by displacement (eps_33 = 0 in plane strain)
    B_eps[3,0::2] = B[1]
    B_eps[3,1::2] = B[0]
    # strain = [e11, e22, e33, e12+e21]
    strain = B_eps*Disp
    #strain = array(B_eps*Disp).flatten()
    #
    # get old values of internal state:
    try:
        stress_p = isvP['stress']
        hydro_p = isvP['pressure']
        strain_p = isvP['strain']
    except:
        stress_p = zeros((4,1))
        hydro_p = 0.
        strain_p = zeros((4,1))
        
    #
    # if split = 'deviatoric', return only the pressure residual:
    if split=='deviatoric':
        
        hydro_strain = (strain[0,0]+strain[1,0]+strain[2,0])/3
        strain = strain - hydro_strain*Ivec
        # how is strain influenced by displacement
        dE_u = matrix([[2/3,-1/3,-1/3,0],[-1/3,2/3,-1/3,0],[-1/3,-1/3,2/3,0],[0,0,0,1]])*B_eps
        
    else:
        dE_u = B_eps
    
    
    #dstrain = strain-strain_p
    
    stress, ddsdde, isvN = FE_materials.mechanical_response(stress_p,strain,0.,T_gauss,dt,material_info,isvP,only_resid=only_resid)
    
    if not 'stress' in isvN.keys():
        isvN['stress'] = array(stress).flatten()
    if not 'strain' in isvN.keys():
        isvN['strain'] = array(strain).flatten()
    
    R_u = B_eps.T*stress
    #
    if only_resid:
        return R_u,0,{}
    #
    #
    dR_u = B_eps.T*ddsdde*dE_u
        
        #####FD tangent:
        #dR_u = zeros((R_u.size,R_u.size))
        #for i in range(Disp.size):
            #dDispFD = dDisp.copy()
            #dDispFD+=1e-10
            #R_u_FD = small_displacement_2D_residual(N,B,Disp0,dDispFD,Temp,dt,material_info,isvP,True)
            #dR_u[i,:] = array(R_u_FD-R_u).flatten()/1e-10
        
        # check determinant of dRdU
#        detRu = linalg.det(ddsdde)
        #print(' Element %i GP %i tangent determinant %e'%(el_nr,el_gp,detRu))
        #if detRu<=0:
            #if el_nr>0:
                #print(' Element %i GP %i tangent issues'%(el_nr,el_gp))
    
    return R_u,dR_u,isvN





def Svec_to_Smat(Svec):
    
    Smat = matrix([[Svec[0,0] ,  0      , Svec[2,0] ,  0],
            [0      , Svec[1,0] ,  0      , Svec[2,0]],
            [Svec[2,0] ,  0      , Svec[1,0] ,  0],
            [0      , Svec[2,0] ,  0      , Svec[0,0]]]);
    return Smat


def Fvec_to_Fmat(Fvec):
    Fmat = matrix([[Fvec[0,0],  0      , 0.5*Fvec[2,0] , 0.5*Fvec[2,0]],
            [0     , Fvec[1,0] , 0.5*Fvec[3,0] , 0.5*Fvec[3,0]],
            [0     , Fvec[2,0] , 0.5*Fvec[0,0] , 0.5*Fvec[0,0]],
            [Fvec[3,0],  0      , 0.5*Fvec[1,0] , 0.5*Fvec[1,0]]]);
    return Fmat







def lagrange_displacement_2D_residual(N,B,Disp0,dDisp,Temp,dt,material_info,isvP,only_resid=False):
    
    #
    # make sure the displacements are in the correct form
    Disp0 = matrix(Disp0).reshape(-1,1)
    dDisp = matrix(dDisp).reshape(-1,1)
    Temp = matrix(Temp).reshape(-1,1)
    # total displacement
    Disp = Disp0+dDisp
    #print("Lagrange Displacement")
    #print(Disp)
    # current temperature:
    T_gauss = N.T*Temp
    #
    #print(B)
    B_eps = zeros((4,Disp.size))
    B_eps[0,0::2] = B[0]
    B_eps[1,1::2] = B[1]
    # B[2] = 0  since eps_33 is unaffected by displacement (eps_33 = 0 in plane strain)
    B_eps[2,0::2] = B[1]  # eps_12
    B_eps[3,1::2] = B[0]  # eps_21
    #print(B_eps)
    # Identity
    Ivec = matrix([[1],[1],[0],[0]])
    #
    #print(B_eps*Disp)
    Fvec = Ivec + B_eps*Disp
    #print(Fvec)
    detF = Fvec[0,0]*Fvec[1,0]-Fvec[2,0]*Fvec[3,0]; # Determinant of F
    Fmat = Fvec_to_Fmat(Fvec)
    # right Cauchy strain:
    Cvec = Fmat.T*Fvec
    # Lagrange strain
    E = 0.5*(Cvec-Ivec)
    # matrix to convert lagrange vector E = [E11, E22, E12, E21] to strain = [e11, e22, e33, e12+e21]
    Econvert = matrix([[1,0,0,0],
                       [0,1,0,0],
                       [0,0,0,0],
                       [0,0,1,1]])
    
    strain = Econvert*E
    # get old values of internal state:
    try:
        stress_p = matrix(isvP['stress']).reshape(-1,1)
        strain_p = matrix(isvP['strain']).reshape(-1,1)
    except:
        stress_p = zeros((4,1))
        strain_p = zeros((4,1))
    #
    dstrain = strain-strain_p
    
    stress, ddsdde, isvN = FE_materials.mechanical_response(stress_p,strain_p,dstrain,T_gauss,dt,material_info,isvP,True,only_resid=only_resid)
    
    # convert stress from PK2 = [S11, S22, S33, S12] to [S11, S22, S12, S21] with S21 = S12
    Sconvert = matrix([[1, 0, 0, 0],
                       [0, 1, 0, 0],
                       [0, 0, 0, 1],
                       [0, 0, 0, 1]])
    #
    Svec = Sconvert*stress
    Smat = Svec_to_Smat(Svec)
    #
    # residual and tangent
    R_u = B_eps.T*(Fmat*Svec)
    
    if only_resid:
        return R_u,0,{}
    
    
    if not 'stress' in isvN.keys():
        isvN['stress'] = array(stress).flatten()
        #
        # stress returned is PK2
        PK2 = matrix([[stress[0,0],stress[3,0],0.],
                    [stress[3,0],stress[1,0],0.],
                    [0.,0.,stress[2,0]]])
        
        FF = matrix([[Fvec[0,0],Fvec[2,0],0.],
                    [Fvec[3,0],Fvec[1,0],0.],
                    [0.,0.,1.]])
        # Cauchy stress
        C0 = (FF*PK2*FF.T)/detF
        Carr = array([C0[0,0],C0[1,1],C0[2,2],(C0[0,1]+C0[1,0])/2])
        
        isvN['stress'] = array(Carr).flatten()
        
    if not 'strain' in isvN.keys():
        isvN['strain'] = array(strain).flatten()
    
    
    #
    dSdE = Sconvert*ddsdde*Econvert
    dR_u = B_eps.T*(Smat + Fmat*dSdE*Fmat.T)*B_eps
    
    return R_u,dR_u,isvN





##@numba.jit
#def material_residual(B,U0,dU,dt,material,isvP):
    #'''
    #return the material residual and tangent
    #'''
    #matcall  =  material['type']
    
    #if matcall == 'Mooney-Rivlin':
        #return mooney_rivlin(B,U0,dU,dt,material,isvP)
    #if matcall == 'Mult-Plast':
        #return mult_plasticity(B,U0,dU,dt,material,isvP)
    #else: # Default to an elastic 'Saint-Venant-Kirchoff' material:
        #return saint_venant(B,U0,dU,dt,material,isvP)
    
    
    

##@numba.jit
#def saint_venant(B,U0,dU,dt,material,isvP):
    ##
    ##
    ## total displacement
    #U = U0+dU
    ## identity
    #Ivec = matrix([1, 1, 0, 0]).T;
    ## deformation gradient
    #Fvec = Ivec + B*U;
    #detF = Fvec[0,0]*Fvec[1,0]-Fvec[2,0]*Fvec[3,0]; # Determinant of F
    #Fmat = Fvec_to_Fmat(Fvec);
##   right Cauchy-Green strain vector
    #Cvec = Fmat.T*Fvec;
##   Total Lagrange strain:
    #Evec = 0.5*(Cvec - Ivec);
    #Evec0 = array([Evec[0,0],Evec[1,0],0.,Evec[2,0],Evec[3,0]])
##
##   stiffness matrix: (plane stress)
    #e,nu = material['E'],material['nu']
    ## muduli
    #L=nu*e/(1.+nu)/(1.-2.*nu)
    #G=e/2./(1.+nu)
    #c1=L+2.*G
    #c2=2*G
    
    #Cmat = zeros((4,4));
    #Cmat[0,0] = c1;
    #Cmat[1,1] = c1;
    #Cmat[0,1] = L;
    #Cmat[1,0] = L;
    #Cmat[2,2] = c2;
    #Cmat[3,3] = c2;
##    
    #Svec = Cmat*Evec;
##   33 component
    #S33 = L*(Evec[0,0]+Evec[1,0])
    
    ##Svec0 = array([Svec[0,0],Svec[1,0],S33,Svec[2,0],Svec[3,0]])
    
    
    #Smat = Svec_to_Smat(Svec);
    ##
    ## residual R(U) and tangent dR(U)/dU
    #R_u  = Fmat*Svec
    #dR_u = (Smat + Fmat*Cmat*Fmat.T)*B
    ##
    ## transform PK2 to Cauchy stress:
    ## sigma = [F] x [S] x [F].T / detJ
    #PK2 = matrix([[Svec[0,0],Svec[2,0],0.],
                  #[Svec[3,0],Svec[1,0],0.],
                  #[0.,0.,S33]])
    
    #FF = matrix([[Fvec[0,0],Fvec[2,0],0.],
                 #[Fvec[3,0],Fvec[1,0],0.],
                 #[0.,0.,1.]])
    
    #C0 = (FF*PK2*FF.T)/detF
    #Carr = array([C0[0,0],C0[1,1],C0[2,2],C0[0,1],C0[1,0]])
    
    
    
    ##
    ## isvs:    
    ##FTF = matrix([[Fvec[0,0]**2,      Fvec[2,0]**2,    Fvec[0,0]*Fvec[2,0], Fvec[2,0]*Fvec[0,0]],
        ##[Fvec[3,0]**2,      Fvec[1,0]**2,    Fvec[3,0]*Fvec[1,0], Fvec[3,0]*Fvec[1,0]],
        ##[Fvec[0,0]*Fvec[3,0], Fvec[2,0]*Fvec[1,0], Fvec[0,0]*Fvec[1,0], Fvec[2,0]*Fvec[3,0]]]);
    ##C0                   = FTF*Svec;
    ##Carr = array([C0[0,0],C0[1,0],S33,C0[2,0]])/detF
    ##
    ##print("\nCauchy compare:")
    ##print(array(Cauchy).flatten())
    ##print(Carr)
    ##print("PK2:")
    ##print(Svec0/detF)
    #isvN = {
        #'stress':Carr,
        #'strain':array(Evec0).flatten()
        #}
    
    #return R_u,dR_u,isvN



##
##
##  multiplicative plasticity:
##
##

#def mult_plasticity(B,U0,dU,dt,material,isvP):
    ##
    ##
    ## total displacement
    #U = U0+dU
    ## identity
    #Ivec = matrix([1, 1, 0, 0]).T;
    ## deformation gradient
    #Fvec = Ivec + B*U;
    #detF = Fvec[0,0]*Fvec[1,0]-Fvec[2,0]*Fvec[3,0]; # Determinant of F
    ##
    #FF = matrix([[Fvec[0,0],Fvec[2,0],0.],
                 #[Fvec[3,0],Fvec[1,0],0.],
                 #[0.,0.,1.]])
    ## right cauchy strain
    #C = FT.T*FF
    ##
    ## previous plastic strain:
    #try:
        #Cp_vec = isvP['pl_strain']
        #peeq = isvP['pl_eq']
    #except:
        #Cp_vec = array([0.,0.,0.,0.])
        #peeq = 0.
    
    #Cp = matrix([[Cp_vec[0],Cp_vec[3],0.],
                 #[Cp_vec[3],Cp_vec[1],0.],
                 #[0.,0.,Cp_vec[2]]])
    ##isvP
    
    
    
    
    
    #Fmat = Fvec_to_Fmat(Fvec);
##   right Cauchy-Green strain vector
    #Cvec = Fmat.T*Fvec;
##   Total Lagrange strain:
    #Evec = 0.5*(Cvec - Ivec);
    #Evec0 = array([Evec[0,0],Evec[1,0],0.,Evec[2,0],Evec[3,0]])
##
##   stiffness matrix: (plane stress)
    #e,nu = material['E'],material['nu']
    ## muduli
    #L=nu*e/(1.+nu)/(1.-2.*nu)
    #G=e/2./(1.+nu)
    #c1=L+2.*G
    #c2=2*G
    
    #Cmat = zeros((4,4));
    #Cmat[0,0] = c1;
    #Cmat[1,1] = c1;
    #Cmat[0,1] = L;
    #Cmat[1,0] = L;
    #Cmat[2,2] = c2;
    #Cmat[3,3] = c2;
##    
    #Svec = Cmat*Evec;
##   33 component
    #S33 = L*(Evec[0,0]+Evec[1,0])
    
    ##Svec0 = array([Svec[0,0],Svec[1,0],S33,Svec[2,0],Svec[3,0]])
    
    
    #Smat = Svec_to_Smat(Svec);
    ##
    ## residual R(U) and tangent dR(U)/dU
    #R_u  = Fmat*Svec
    #dR_u = (Smat + Fmat*Cmat*Fmat.T)*B
    ##
    ## transform PK2 to Cauchy stress:
    ## sigma = [F] x [S] x [F].T / detJ
    #PK2 = matrix([[Svec[0,0],Svec[2,0],0.],
                  #[Svec[3,0],Svec[1,0],0.],
                  #[0.,0.,S33]])
    
    #FF = matrix([[Fvec[0,0],Fvec[2,0],0.],
                 #[Fvec[3,0],Fvec[1,0],0.],
                 #[0.,0.,1.]])
    
    #C0 = (FF*PK2*FF.T)/detF
    #Carr = array([C0[0,0],C0[1,1],C0[2,2],C0[0,1],C0[1,0]])
    
    
    
    ##
    ## isvs:    
    ##FTF = matrix([[Fvec[0,0]**2,      Fvec[2,0]**2,    Fvec[0,0]*Fvec[2,0], Fvec[2,0]*Fvec[0,0]],
        ##[Fvec[3,0]**2,      Fvec[1,0]**2,    Fvec[3,0]*Fvec[1,0], Fvec[3,0]*Fvec[1,0]],
        ##[Fvec[0,0]*Fvec[3,0], Fvec[2,0]*Fvec[1,0], Fvec[0,0]*Fvec[1,0], Fvec[2,0]*Fvec[3,0]]]);
    ##C0                   = FTF*Svec;
    ##Carr = array([C0[0,0],C0[1,0],S33,C0[2,0]])/detF
    ##
    ##print("\nCauchy compare:")
    ##print(array(Cauchy).flatten())
    ##print(Carr)
    ##print("PK2:")
    ##print(Svec0/detF)
    #isvN = {
        #'stress':Carr,
        #'strain':array(Evec0).flatten()
        #}
    
    #return R_u,dR_u,isvN