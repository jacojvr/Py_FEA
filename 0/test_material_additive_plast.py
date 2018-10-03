

from numpy import array, matrix, zeros, sort, sqrt, linalg
from FE_materials import interpolate_material_properties, isotropic_mises_hardening
    
    
    
    
    
    
    
#def additive_mises_plastic(stress_p,strain_in,dstrain,T_gauss,dt,material_info,isvP):
    ## 
    #matprops = interpolate_material_properties(T_gauss,material_info)
    
    ## strain = [11,22,33,12]
    ##
    ##
    #[E,nu] = matprops['*elastic']
    ## Various elastic parameter values & lame params
    #xk = E/(1-2*nu)
    #eg2 = E/(1+nu)
    #al = (xk-eg2)/3
    ## bulk and shear modulus
    #bulk_mod = xk/3
    #eg = eg2/2
    ##
    ## tolerance:
    #toler = 1e-8
    ## 2/3
    #k23 = 2./3
    #k32 = 3./2
    ## sqrt(2./3)
    #sq23 = sqrt(k23)
    ##
    ## Identity vector
    #Ivec = matrix([[1],[1],[1],[0]])
    ##
    ## get previous plastic strains:
    #try:
        #strain_pl = matrix(isvP['plastic strain']).reshape(-1,1)
        #eq_pl = isvP['equivalent plastic'][0]
    #except:
        #strain_pl = zeros((4,1))
        #eq_pl = 0
    ##
    ##
    ## elastic strain component
    #strain_in = matrix(strain_in).reshape(-1,1)
    #strain = strain_in-strain_pl
    ##
    ## d(stress)/d(strain)
    #ddsdde = matrix([[al+eg2, al, al, 0.],
                     #[al, al+eg2, al, 0.],
                     #[al, al, al+eg2, 0.],
                     #[0., 0., 0., eg]])
    ## stress
    #stress = ddsdde*strain
    ## von Mises equivalent stress
    #smises=(stress[0,0]-stress[1,0])**2+(stress[1,0]-stress[2,0])**2+(stress[2,0]-stress[0,0])**2+6*stress[3,0]**2
    #smises=sqrt(smises/2)
    ##
    ## call hardening function
    #deqpl = 0.
    #sy, dsy, isvN = isotropic_mises_hardening(eq_pl,deqpl,T_gauss,dt,matprops,isvP)
    ##
    ##print('yield = ',sy)
    ## trial yield function:
    #ftrial = smises-sy
    ##print("ftrial = ",ftrial)
    ##
    ## check if necessary and do plastic
    #if ftrial>toler:
        ##print('PLASTIC')
        ## plastic deformation:
        #kount,kmax = 0,15 # newton loop initial and max
        #eg3 = eg*3
        ## newton loop to find new plastic strain
        #while abs(ftrial)>toler and kount<kmax:
            #kount+=1
            #deqpl=abs(deqpl+ftrial/(eg3+dsy))
            #sy, dsy, isvN = isotropic_mises_hardening(eq_pl,deqpl,T_gauss,dt,matprops,isvP)
            #ftrial=smises-eg3*deqpl-sy
            ##print("ftrial = ",ftrial)
        ##
        ## update state
        ##hydrostatic and deviatoric stress
        #shydro = (stress[0,0]+stress[1,0]+stress[2,0])/3.
        #stress_dev = stress - shydro*Ivec
        ##
        ## direction of plastic flow
        #flow_direction = stress_dev/smises
        ##flow_direction[3,0] *= 2 # e12+e21
        ###
        #stress = sy*flow_direction + shydro*Ivec
        ##smises=(stress[0,0]-stress[1,0])**2+(stress[1,0]-stress[2,0])**2+(stress[2,0]-stress[0,0])**2+6*stress[3,0]**2
        ##smises=sqrt(smises/2)
        ##print("S1 = %.4f, %.4f, %.4f, %.4f"%tuple(array(stress).flatten()))
        ##print("Smises 1 vs Syield = ",smises,' vs ',sy)
        ## update equivalent plastic strain
        #eq_pl += deqpl
        ## update plastic strain tensor
        ##what should plastic strain be:
        ##strain_el = linalg.inv(ddsdde)*stress 
        ##strain_pl0 = strain_in - strain_el
        ##print("plastic strain 1 = ",strain_pl0)
        ##
        #gamma = deqpl*3/2
        #strain_pl = strain_pl + gamma*flow_direction
        ##print("plastic strain 2 = ",strain_pl)
        ###
        ##print(" SF plastic = ")
        ##print(array(strain_pl)/array(strain_pl0))
        ## update stress: (ddsdde still the elastic one):
        ## elastic strain component
        ##strain = strain_in-strain_pl
        ##stress = ddsdde*strain
        ##smises=(stress[0,0]-stress[1,0])**2+(stress[1,0]-stress[2,0])**2+(stress[2,0]-stress[0,0])**2+6*stress[3,0]**2
        ##smises=sqrt(smises/2)
        ##print("S2 = %.4f, %.4f, %.4f, %.4f"%tuple(array(stress).flatten()))
        ##print("Smises 2 = ",smises)
        ##
        ##
        #effg=eg*sy/smises
        #effg2=eg2*sy/smises
        #effg3=k32*effg2
        #efflam=(xk-effg2)/3
        #effhrd=eg3*dsy/(eg3+dsy)-effg3
        
        #ddsdde = matrix([[effg2+efflam, efflam, efflam, 0.],
                         #[efflam, effg2+efflam, efflam, 0.],
                         #[efflam, efflam, effg2+efflam, 0.],
                         #[0., 0., 0., effg]])
        
        #ddsdde += effhrd*flow_direction*flow_direction.T
      
      
    
    
    ## store ISVs        
    #isvN['equivalent plastic'] = array(eq_pl).flatten()
    #isvN['plastic strain'] = array(strain_pl).flatten()
    #isvN['stress'] = array(stress).flatten()
    
    #return stress, ddsdde, isvN


def additive_mises_plastic(stress_p,strain_in,dstrain,T_gauss,dt,material_info,isvP):
    # 
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
    strain_in = matrix(strain_in).reshape(-1,1)
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
        # update equivalent plastic strain
        eq_pl += deqpl
        gamma = deqpl*3/2
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


material_info ={'*elastic':{0:[200000,0.3]},
                '*plastic':{0:[50,0,150,1]}
    }

isvP = {}

strain = 10*matrix([0.01,0.025,0.,0.09]).T

stress, ddsdde, isvN = additive_mises_plastic(0.,strain,0.,0.,0.,material_info,isvP)

ddsdde_FD = zeros((4,4))
for i in range(4):
    strain_FD = strain.copy()
    strain_FD[i,0] += 1e-5
    stress_FD = array(additive_mises_plastic(0.,strain_FD,0.,0.,0.,material_info,isvP)[0]).flatten()
    ddsdde_FD[i,:] = (stress_FD-array(stress).flatten())/1e-5
    
print("FINITE DIFFERENCE: ")
for i in range(4):
    print("%.4e, %.4e, %.4e, %.4e"%tuple(array(ddsdde_FD[i,:]).flatten()))
    
    
print("ANALYTICAL: ")
for i in range(4):
    print("%.4e, %.4e, %.4e, %.4e"%tuple(array(ddsdde[i,:]).flatten()))
    

