

from numpy import array, matrix, zeros, sort, sqrt, linalg
#from FE_materials import interpolate_material_properties, isotropic_mises_hardening
import FE_residual
import FE_shapefn
import FE_materials
from FE_residual import Svec_to_Smat, Fvec_to_Fmat
    
    




# displacement residual using 2x2 gauss point integration
def Residual_Disp_Fbar_2x2(XY,Disp0,dDisp,Temp,dt,ISVs,section_info,material_info,nlgeom=False,el_nr=0,calc_dR_u=False):
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
        
        stress, ddsdde, isvN = FE_materials.mechanical_response(stress_p,strain,dstrain,T_gauss,dt,material_info,isvP,True,calc_ddsdde=calc_dR_u)
        
        
        if not 'stress' in isvN.keys():
            isvN['stress'] = array(stress).flatten()
            # stress returned is PK2
            PK2 = matrix([[stress[0,0],stress[3,0],0.],
                        [stress[3,0],stress[1,0],0.],
                        [0.,0.,stress[2,0]]])
            
            FF = matrix([[FBvec[0,0],FBvec[2,0],0.],
                        [FBvec[3,0],FBvec[1,0],0.],
                        [0.,0.,1.]])
            # Cauchy stress
            C0 = (FF*PK2*FF.T)/detF
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
        R_u = B_eps.T*(FBmat*Svec)
        #
        dSdE = Sconvert*ddsdde*Econvert
        dR_u = B_eps.T*(Smat + FBmat*dSdE*FBmat.T)*B_eps
    
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
#            Ge      += detJ*B.T.flatten().T
        except:
            Res     = thickness*R_u*detJ
            Tangent = thickness*dR_u*detJ
#            Ge      = detJ*B.T.flatten().T
#        # state variables:
        newISVs[Gind] = isvN
    #
    #
    #
    # Convert and return 1D arrays
    Res = array(Res).flatten()
    
    if el_nr ==0:
        return Res
    
    # 
    #
    if calc_dR_u:
#        Tangent *= 0.
##        Tangent = Tangent.reshape(int(sqrt(Tangent.size)),-1)
#        for i in range(dDisp.size):
#            dDispFD = dDisp.copy()
#            dDispFD[i,0] += 1e-6
#            ResFD = Residual_Disp_Fbar_4_1(XY,Disp0,dDispFD,Temp,dt,ISVs,section_info,material_info,nlgeom,0)
#            Tangent[i,:] = array(ResFD - Res).flatten()/1e-6
#        Tangent = Tangent.T
#        
        Tangent = array(Tangent).flatten()
    
    
    
    return [Res,Tangent,newISVs]













# eigt noded F-bar element with small strain plasticity
# element coordinates:
XY = array([[0.2, 0],
            [0.8, 0.5],
            [1. , 2.],
            [0. , 1.],
            [0.5, 0.25],
            [0.9, 1.25],
            [0.5, 1.5],
            [0.1, 0.5]])

Disp0 = matrix([0]*16)

dDisp = matrix([0.,
                0.1,
                -0.01,
                0.04,
                -0.02,
                0.13,
                0.08,
                0.09,
                -0.04,
                0.1,
                0.,
                -0.11,
                -0.07,
                0.,
                0.,
                0.1]).T





material_info ={'*elastic':{0:[200000,0.3]},
                '*plastic':{0:[50,0,150,1]}
    }

section_info = {'*thickness':1.,
                '*material':'plmat'}

nlgeom = False
el_nr=1
calc_dR_u=True

ISVs = {}

dt = 1.
Temp = matrix([0]*8).T


Res,Tangent,newISVs =  Residual_Disp_Fbar_2x2(XY,Disp0,dDisp,Temp,dt,ISVs,section_info,material_info,nlgeom,el_nr,calc_dR_u)


#
# FINITE DIFFERENCE TANGENT
dR_u = zeros((Res.size,Res.size))
for i in range(dDisp.size):
    dDispFD = dDisp.copy()
    dDispFD[i,0]+=1e-6
    R_u_FD = Residual_Disp_Fbar_2x2(XY,Disp0,dDispFD,Temp,dt,ISVs,section_info,material_info,nlgeom,el_nr,calc_dR_u)[0]
    dR_u[i,:] = array(R_u_FD-Res).flatten()/1e-6
    
dR_u = dR_u.T.flatten()

print('Analytical vs FD:')
printstr = '%.7f\t%.7f\t\t%.7f\n'*dR_u.size
print(printstr%tuple(array([Tangent,dR_u,Tangent/dR_u]).T.flatten()))


#ddsdde_FD = zeros((8,8))
#for i in range(8):
#    strain_FD = strain.copy()
#    strain_FD[i,0] += 1e-5
#    stress_FD = array(additive_mises_plastic(0.,strain_FD,0.,0.,0.,material_info,isvP)[0]).flatten()
#    ddsdde_FD[i,:] = (stress_FD-array(stress).flatten())/1e-5
#    
#print("FINITE DIFFERENCE: ")
#for i in range(4):
#    print("%.4e, %.4e, %.4e, %.4e"%tuple(array(ddsdde_FD[i,:]).flatten()))
#    
#    
#print("ANALYTICAL: ")
#for i in range(4):
#    print("%.4e, %.4e, %.4e, %.4e"%tuple(array(ddsdde[i,:]).flatten()))
    

