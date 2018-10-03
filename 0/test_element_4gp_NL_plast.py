

from numpy import array, matrix, zeros, sort, sqrt, linalg
#from FE_materials import interpolate_material_properties, isotropic_mises_hardening
import FE_residual
import FE_shapefn
import FE_materials
    
    




# displacement residual using 2x2 gauss point integration
def Residual_Disp_2x2(XY,Disp0,dDisp,Temp,dt,ISVs,section_info,material_info,nlgeom=False,el_nr=0,calc_dR_u=False):
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
            R_u,dR_u,isvN = FE_residual.lagrange_displacement_2D_residual(N,B,Disp0,dDisp,Temp,dt,material_info,isvP,calc_dR_u=calc_dR_u)
            
        else:
            # small strain formulation E = B*U
            R_u,dR_u,isvN = FE_residual.small_displacement_2D_residual(N,B,Disp0,dDisp,Temp,dt,material_info,isvP,False,el_nr,Gind,calc_dR_u=calc_dR_u)
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
    
    
    if el_nr ==0:
        return Res
    #
    
    if calc_dR_u:
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

nlgeom = not False
el_nr=1
calc_dR_u=True

ISVs = {}

dt = 1.
Temp = matrix([0]*8).T


Res,Tangent,newISVs =  Residual_Disp_2x2(XY,Disp0,dDisp,Temp,dt,ISVs,section_info,material_info,nlgeom,el_nr,calc_dR_u)


#
# FINITE DIFFERENCE TANGENT
dR_u = zeros((Res.size,Res.size))
for i in range(dDisp.size):
    dDispFD = dDisp.copy()
    dDispFD[i,0]+=1e-6
    R_u_FD = Residual_Disp_2x2(XY,Disp0,dDispFD,Temp,dt,ISVs,section_info,material_info,nlgeom,el_nr,calc_dR_u)[0]
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
    

