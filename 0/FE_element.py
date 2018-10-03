#
#
#  assemble and solve Finite element problem
#
#
#

from numpy import zeros, matrix, array, linalg, r_, sqrt

import FE_shapefn
import FE_residual
import FE_materials
from FE_residual import Svec_to_Smat, Fvec_to_Fmat

 
    
    
# return element temperature residual and tangent
def Element_Temp(XY,Temp0,dTemp,dt,isvP,el_type,section_info,material_info):
    
    if el_type == 'cpe4':
        return Residual_Temp_2x2(XY,Temp0,dTemp,dt,isvP,section_info,material_info)
    
    
    
    
    
    
def Element_Disp(XY,Disp0,dDisp,Temp,dt,isvP,el_type,section_info,material_info,nlgeom=False,el_nr=0,only_resid=False):
    
    if el_type == 'cpe4':
        if not '*plastic' in material_info.keys():
            return Residual_Disp_2x2(XY,Disp0,dDisp,Temp,dt,isvP,section_info,material_info,nlgeom,el_nr,only_resid=only_resid)
        else:
            return Residual_Disp_split_4_1(XY,Disp0,dDisp,Temp,dt,isvP,section_info,material_info,nlgeom,el_nr,only_resid=only_resid)
    
    if el_type == 'cpe8r':
        if not '*plastic' in material_info.keys():
            return Residual_Disp_2x2(XY,Disp0,dDisp,Temp,dt,isvP,section_info,material_info,nlgeom,el_nr,only_resid=only_resid)
        else:
           return Residual_Disp_2x2(XY,Disp0,dDisp,Temp,dt,isvP,section_info,material_info,nlgeom,el_nr,only_resid=only_resid)
            #return Residual_Disp_Fbar_2x2(XY,Disp0,dDisp,Temp,dt,isvP,section_info,material_info,nlgeom,el_nr,only_resid=only_resid)
         
    if el_type == 'cpe8':
        if not '*plastic' in material_info.keys():
            return Residual_Disp_3x3(XY,Disp0,dDisp,Temp,dt,isvP,section_info,material_info,nlgeom,el_nr,only_resid=only_resid)
        else:
            return Residual_Disp_3x3(XY,Disp0,dDisp,Temp,dt,isvP,section_info,material_info,nlgeom,el_nr,only_resid=only_resid)
    
    
    
    
    

def Residual_Temp_2x2(XY,Temp0,dTemp,dt,ISVs,section_info,material_info):
    #
    #
    # section thikness:
    thickness = section_info['*thickness']
    # section velocity in undeformed configuration
    if '*velocity' in section_info.keys():
        velocity = section_info['*velocity']
    else:
        velocity = None
    #
    # internal state variables:
    newISVs = {}
    # Gauss coordinates:
    Gv = 1/sqrt(3);             # Gauss point location
    GPs = [(-Gv,-Gv),(Gv,-Gv),(Gv,Gv),(-Gv,Gv)]  # 2 by 2 Gauss integration loops
    Gind = 0 #gauss point index
    for xi,eta in GPs:   
        #
        Gind += 1
        # Gauss point Internal State Variables:
        try:
            isvP = ISVs[Gind]
        except:
            isvP = {}
        #
        # Shape functions:
        N = matrix(FE_shapefn.shapefn(nnodes,xi,eta)).T
        # get B and detJ
        B,detJ = FE_shapefn.nabla_shapefn(nnodes,xi,eta,XY)
        #
        # material temperature residual, tangent and ISVs:
        if velocity is None:
            R_t,dR_t,isvN = FE_residual.temperature_residual(N,B,Temp0,dTemp,dt,material_info,isvP)
        else:
            R_t,dR_t,isvN = FE_residual.temperature_residual_velocity(N,B,Temp0,dTemp,velocity,dt,material_info,isvP)
        #
        # to integrate, sum over material volume and gauss point weight (2x2 : w=1)
        try:
            Res     += thickness*R_t*detJ
            Tangent += thickness*dR_t*detJ
        except:
            Res     = thickness*R_t*detJ
            Tangent = thickness*dR_t*detJ
        # state variables:
        newISVs[Gind] = isvN
    # Convert and return 1D arrays
    Res = array(Res).flatten()
    Tangent = array(Tangent).flatten()
    #
    return [Res,Tangent,newISVs]






#
# displacement residual using 2x2 gauss point integration
def Residual_Disp_2x2(XY,Disp0,dDisp,Temp,dt,ISVs,section_info,material_info,nlgeom=False,el_nr=0,only_resid=False):
    
    #
    #
    #if el_nr==1:
        #print('********** DISPLACEMENT RESIDUAL USING 2x2 GP INTEGRATION')
    #
    # section thikness:
    thickness = section_info['*thickness']
    # section velocity in undeformed configuration
    if '*velocity' in section_info.keys():
        velocity = section_info['*velocity']
    else:
        velocity = None
    #
    Disp0 = matrix(Disp0).reshape(-1,1)
    dDisp = matrix(dDisp).reshape(-1,1)
    Temp = matrix(Temp).reshape(-1,1)
    #
    # number of nodes:
    nnodes = XY.shape[0]
    # internal state variables:
    newISVs = {}
    # Gauss coordinates:
    Gv = 1/sqrt(3);             # Gauss point location
    GPs = [(-Gv,-Gv),(Gv,-Gv),(Gv,Gv),(-Gv,Gv)]  # 2 by 2 Gauss integration loops
    Gind = 0 #gauss point index
    for xi,eta in GPs:   
        #
        Gind += 1
        # Gauss point Internal State Variables:
        try:
            isvP = ISVs[Gind]
        except:
            isvP = {}
        #
        # Shape functions:
        N = matrix(FE_shapefn.shapefn(nnodes,xi,eta)).T
        B,detJ = FE_shapefn.nabla_shapefn(nnodes,xi,eta,XY)
        #
        if nlgeom:
            # lagrange strain E = (F.T*F - I)/2 = (U.T*B.T*B*U - I)/2  
            #if el_nr==1:
                #print('********** 2D Total Lagrange displacement')
            R_u,dR_u,isvN = FE_residual.lagrange_displacement_2D_residual(N,B,Disp0,dDisp,Temp,dt,material_info,isvP,only_resid=only_resid)
            
        else:
            #if el_nr==1:
                #print('********** 2D Small Strain displacement')
            # small strain formulation E = B*U
            R_u,dR_u,isvN = FE_residual.small_displacement_2D_residual(N,B,Disp0,dDisp,Temp,dt,material_info,isvP,el_nr,Gind,only_resid=only_resid)
        #
        # to integrate, sum over material volume and gauss point weight (2x2 : w=1)
        try:
            Res     += thickness*R_u*detJ
            Tangent += thickness*dR_u*detJ
        except:
            Res     = thickness*R_u*detJ
            Tangent = thickness*dR_u*detJ
        # state variables:
        newISVs[Gind] = isvN
    
    
    # Convert and return 1D arrays
    Res = array(Res).flatten()
    
    
    if only_resid:
        return Res
    #
        # 
#        Tangent *= 0.
#        for i in range(dDisp.size):
#            dDispFD = dDisp.copy()
#            dDispFD[i,0] += 1e-8
#            ResFD = Residual_Disp_split_4_1(XY,Disp0,dDispFD,Temp,dt,ISVs,section_info,material_info,nlgeom,0)
#            Tangent[i,:] = (ResFD - Res)/1e-8
#        Tangent = Tangent.T
#        
#        Tangent = array(Tangent).flatten()
        # Convert and return 1D arrays
    Tangent = array(Tangent).flatten()
    
    
    
    return [Res,Tangent,newISVs]




#
# displacement residual using 2x2 gauss point integration
def Residual_Disp_split_4_1(XY,Disp0,dDisp,Temp,dt,ISVs,section_info,material_info,nlgeom=False,el_nr=0,only_resid=False):
    #
    # section thikness:
    thickness = section_info['*thickness']
    # section velocity in undeformed configuration
    if '*velocity' in section_info.keys():
        velocity = section_info['*velocity']
    else:
        velocity = None
    #
    Disp0 = matrix(Disp0).reshape(-1,1)
    dDisp = matrix(dDisp).reshape(-1,1)
    Temp = matrix(Temp).reshape(-1,1)
    #
    # total displacement:
    Disp = Disp0+dDisp
    # number of nodes:
    nnodes = XY.shape[0]
    # internal state variables:
    newISVs = {}
    #
    #  4x4 integration on deviatoric part of deformtaion
    # first do the deviatoric part
    # Gauss coordinates:
    Gv = 1/sqrt(3);             # Gauss point location
    GPs = [(-Gv,-Gv),(Gv,-Gv),(Gv,Gv),(-Gv,Gv)]  # 2 by 2 Gauss integration loops
    Gind = 0 #gauss point index
    #
    # total volume:
    Vol = 0.
    for xi,eta in GPs:   
        #
        Gind += 1
        # Gauss point Internal State Variables:
        try:
            isvP = ISVs[Gind]
        except:
            isvP = {}
        #
        # Shape functions:
        N = matrix(FE_shapefn.shapefn(nnodes,xi,eta)).T
        B,detJ = FE_shapefn.nabla_shapefn(nnodes,xi,eta,XY)
        #
        Vol += detJ
#        print(B.shape)
        #
        # only small strain implemented:
        #
        R_u,dR_u,isvN = FE_residual.small_displacement_2D_residual(N,B,Disp0,dDisp,Temp,dt,material_info,isvP,el_nr,Gind,'deviatoric',only_resid)
        #
        # to integrate, sum over material volume and gauss point weight (2x2 : w=1)
        try:
            Res     += thickness*R_u*detJ
            Tangent += thickness*dR_u*detJ
            Ge      += detJ*B.T.flatten().T
        except:
            Res     = thickness*R_u*detJ
            Tangent = thickness*dR_u*detJ
            Ge      = detJ*B.T.flatten().T
#        # state variables:
        newISVs[Gind] = isvN
    #
    #
    # add the isochoric/volumetric/pressure term :
    
    # element pressure:
    N = matrix(FE_shapefn.shapefn(nnodes,0,0)).T
    T_gauss = N.T*Temp
    matprops = FE_materials.interpolate_material_properties(T_gauss,material_info)
    [E,nu] = matprops['*elastic']
    bulk_mod = E/(1-2*nu)/3
    Me = Vol/bulk_mod
    pressure = -1/Me * (Ge.T*Disp)
    Res += -Ge*pressure
    Tangent += Ge*Ge.T/Me
    
    # Convert and return 1D arrays
    Res = array(Res).flatten()
    
    
    if only_resid:
        return Res
    
#        Tangent *= 0.
##        Tangent = Tangent.reshape(int(sqrt(Tangent.size)),-1)
#        for i in range(dDisp.size):
#            dDispFD = dDisp.copy()
#            dDispFD[i,0] += 1e-6
#            ResFD = Residual_Disp_split_4_1(XY,Disp0,dDispFD,Temp,dt,ISVs,section_info,material_info,nlgeom,0)
#            Tangent[i,:] = array(ResFD - Res).flatten()/1e-6
#        Tangent = Tangent.T
        
    Tangent = array(Tangent).flatten()
    
    
    
    return [Res,Tangent,newISVs]






#
# displacement residual using 2x2 gauss point integration
def Residual_Disp_3x3(XY,Disp0,dDisp,Temp,dt,ISVs,section_info,material_info,nlgeom=False,el_nr=0,only_resid=False):
    #
    # section thikness:
    thickness = section_info['*thickness']
    # section velocity in undeformed configuration
    if '*velocity' in section_info.keys():
        velocity = section_info['*velocity']
    else:
        velocity = None
    #
    Disp0 = matrix(Disp0).reshape(-1,1)
    dDisp = matrix(dDisp).reshape(-1,1)
    Temp = matrix(Temp).reshape(-1,1)
    #
    # number of nodes:
    nnodes = XY.shape[0]
    # internal state variables:
    newISVs = {}
    # Gauss coordinates:
    Gv = sqrt(3/5);             # Gauss point location
    GPs = [(-Gv,-Gv),(0.,-Gv),(Gv,-Gv),
           (-Gv,0),(0.,0),(Gv,0),
           (-Gv,Gv),(0.,Gv),(Gv,Gv)]  # 3 by 3 Gauss integration loops
    w1,w2 = 5/9,8/9
    weights = [w1*w1, w2*w1, w1,w1,
               w1*w2, w2*w2, w1*w2,
               w1*w1, w2*w1, w1*w1]
    #
    #
    
    
    Gind = 0 #gauss point index
    for xi,eta in GPs:   
        #
        GPweight = weights[Gind] 
        #
        Gind += 1
        # Gauss point Internal State Variables:
        try:
            isvP = ISVs[Gind]
        except:
            isvP = {}
        #
        # Shape functions:
        N = matrix(FE_shapefn.shapefn(nnodes,xi,eta)).T
        B,detJ = FE_shapefn.nabla_shapefn(nnodes,xi,eta,XY)
        #
        if nlgeom:
            # lagrange strain E = (F.T*F - I)/2 = (U.T*B.T*B*U - I)/2  
            R_u,dR_u,isvN = FE_residual.lagrange_displacement_2D_residual(N,B,Disp0,dDisp,Temp,dt,material_info,isvP,only_resid=only_resid)
            
        else:
            # small strain formulation E = B*U
            R_u,dR_u,isvN = FE_residual.small_displacement_2D_residual(N,B,Disp0,dDisp,Temp,dt,material_info,isvP,el_nr,Gind,only_resid=only_resid)
        #
        # to integrate, sum over material volume and gauss point weight (2x2 : w=1)
        try:
            Res     += GPweight*thickness*R_u*detJ
            Tangent += GPweight*thickness*dR_u*detJ
        except:
            Res     = GPweight*thickness*R_u*detJ
            Tangent = GPweight*thickness*dR_u*detJ
        # state variables:
        newISVs[Gind] = isvN
    
    
    # Convert and return 1D arrays
    Res = array(Res).flatten()
    
    if only_resid:
        return Res
    #
    
        # 
#        Tangent *= 0.
#        for i in range(dDisp.size):
#            dDispFD = dDisp.copy()
#            dDispFD[i,0] += 1e-8
#            ResFD = Residual_Disp_split_4_1(XY,Disp0,dDispFD,Temp,dt,ISVs,section_info,material_info,nlgeom,0)
#            Tangent[i,:] = (ResFD - Res)/1e-8
#        Tangent = Tangent.T
#        
#        Tangent = array(Tangent).flatten()
        # Convert and return 1D arrays
    Tangent = array(Tangent).flatten()
    
    
    
    return [Res,Tangent,newISVs]




def Residual_Disp_Fbar_2x2(XY,Disp0,dDisp,Temp,dt,ISVs,section_info,material_info,nlgeom=False,el_nr=0,only_resid=False):
    #
    # section thikness:
    thickness = section_info['*thickness']
    # section velocity in undeformed configuration
    if '*velocity' in section_info.keys():
        velocity = section_info['*velocity']
    else:
        velocity = None
    #
    Disp0 = matrix(Disp0).reshape(-1,1)
    dDisp = matrix(dDisp).reshape(-1,1)
    Temp = matrix(Temp).reshape(-1,1)
    #
    # total displacement:
    Disp = Disp0+dDisp
    # number of nodes:
    nnodes = XY.shape[0]
    # internal state variables:
    newISVs = {}
    #
    # First , get the deformation gradient at the centre of the element:
    N = matrix(FE_shapefn.shapefn(nnodes,0,0)).T
    B,detJ = FE_shapefn.nabla_shapefn(nnodes,0,0,XY)
    #
    B_eps0 = zeros((4,Disp.size))
    B_eps0[0,0::2] = B[0]
    B_eps0[1,1::2] = B[1]
    # B[2] = 0  since eps_33 is unaffected by displacement (eps_33 = 0 in plane strain)
    B_eps0[2,0::2] = B[1]  # eps_12
    B_eps0[3,1::2] = B[0]  # eps_21
    # Identity
    Ivec = matrix([[1],[1],[0],[0]])
    #
    Fvec0 = Ivec + B_eps0*Disp
    detF0 = Fvec0[0,0]*Fvec0[1,0]-Fvec0[2,0]*Fvec0[3,0]; # Determinant of F
    #
    #
    #  4x4 integration on deviatoric part of deformtaion
    # first do the deviatoric part
    # Gauss coordinates:
    Gv = 1/sqrt(3);             # Gauss point location
    GPs = [(-Gv,-Gv),(Gv,-Gv),(Gv,Gv),(-Gv,Gv)]  # 2 by 2 Gauss integration loops
    Gind = 0 #gauss point index
    #
    # total volume:
#    Vol = 0.
    for xi,eta in GPs:   
        #
        Gind += 1
        # Gauss point Internal State Variables:
        try:
            isvP = ISVs[Gind]
        except:
            isvP = {}
        #
        # Shape functions:
        N = matrix(FE_shapefn.shapefn(nnodes,xi,eta)).T
        B,detJ = FE_shapefn.nabla_shapefn(nnodes,xi,eta,XY)
        #
        # temperature:
        T_gauss = N.T*Temp
        # Deformation gradient
        B_eps = zeros((4,Disp.size))
        B_eps[0,0::2] = B[0]
        B_eps[1,1::2] = B[1]
        # B[2] = 0  since eps_33 is unaffected by displacement (eps_33 = 0 in plane strain)
        B_eps[2,0::2] = B[1]  # eps_12
        B_eps[3,1::2] = B[0]  # eps_21
        # Identity
        Ivec = matrix([[1],[1],[0],[0]])
        #
        Fvec = Ivec + B_eps*Disp
        Fmat = Fvec_to_Fmat(Fvec)
        detF = Fvec[0,0]*Fvec[1,0]-Fvec[2,0]*Fvec[3,0]; # Determinant of F
        detSF = (detF0/detF)**(1./3) #isochoric 
        #
        # F-bar vectcor and matrix
        FBvec = detSF*Fvec
        FBmat = Fvec_to_Fmat(FBvec)
        
            
        # right Cauchy strain:
        Cvec = FBmat.T*FBvec
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
        
        stress, ddsdde, isvN = FE_materials.mechanical_response(stress_p,strain,dstrain,T_gauss,dt,material_info,isvP,True,only_resid=only_resid)
        
        
        if not 'stress' in isvN.keys():
            isvN['stress'] = array(stress).flatten()
            # stress returned is PK2
            PK2 = matrix([[stress[0,0],stress[3,0],0.],
                        [stress[3,0],stress[1,0],0.],
                        [0.,0.,stress[2,0]]])
            
            FF = matrix([[Fvec[0,0],Fvec[2,0],0.],
                        [Fvec[3,0],Fvec[1,0],0.],
                        [0.,0.,1.]])
        
            FFB = matrix([[FBvec[0,0],FBvec[2,0],0.],
                        [FBvec[3,0],FBvec[1,0],0.],
                        [0.,0.,1.]])
            # Cauchy stress
            C0 = (FF*PK2*FFB.T)/detF0
            Carr = array([C0[0,0],C0[1,1],C0[2,2],(C0[0,1]+C0[1,0])/2])
            
            isvN['stress'] = array(Carr).flatten()
            
            
            
        if not 'strain' in isvN.keys():
            isvN['strain'] = array(strain).flatten()
        
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
        #
        dSdE = Sconvert*ddsdde*Econvert
        dR_u = B_eps.T*(Smat + Fmat*dSdE*FBmat.T)*B_eps
    
        #
        #
        #
        #
        #
        #
        #
        # to integrate, sum over material volume and gauss point weight (2x2 : w=1)
        try:
            Res     += thickness*R_u*detJ
            Tangent += thickness*dR_u*detJ
            Ge      += detJ*B.T.flatten().T
        except:
            Res     = thickness*R_u*detJ
            Tangent = thickness*dR_u*detJ
            Ge      = detJ*B.T.flatten().T
#        # state variables:
        newISVs[Gind] = isvN
    #
    #
    #
    # Convert and return 1D arrays
    Res = array(Res).flatten()
    
    
    if only_resid:
        return Res
    
    # 
    #
#    Tangent *= 0.
##        Tangent = Tangent.reshape(int(sqrt(Tangent.size)),-1)
#    for i in range(dDisp.size):
#        dDispFD = dDisp.copy()
#        dDispFD[i,0] += 1e-8
#        ResFD = Residual_Disp_Fbar_2x2(XY,Disp0,dDispFD,Temp,dt,ISVs,section_info,material_info,nlgeom,0)
#        Tangent[i,:] = array(ResFD - Res).flatten()/1e-8
#    Tangent = Tangent.T
##        
#    Tangent = array(Tangent).flatten()
    
    
    
    return [Res,Tangent,newISVs]

