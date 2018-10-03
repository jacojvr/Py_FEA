import numpy as np


def shapefn(nnodes,xi,eta):
    
    if nnodes==4: return q4_shapefn(xi,eta)
    if nnodes==8: return q8_shapefn(xi,eta)



def nabla_shapefn(nnodes,xi,eta,XY):
    
    if nnodes==4: return q4_nabla(xi,eta,XY)
    if nnodes==8: return q8_nabla(xi,eta,XY)





def q4_nabla(xi,eta,X):
    #
    # Derivatives of shape functions wrt xi & eta
    dNdxi = np.matrix([[eta-1,1-eta,1+eta,-1-eta],
                       [xi-1,-1-xi,1+xi,1-xi]])*.25
    #
    J      = dNdxi*X;                # dXdxi
    detJ   = np.linalg.det(J);       # Determinant of Jacobian J
    invJ   = np.linalg.inv(J);       # Inverse of Jacobian
    #
    ## for a plain B matrix where B_ij= dN_i / dX_j:
    B = invJ*dNdxi #Shape function derivatives wrt x and y
    
    # B values may be inseted into various other places in so that 
    #   >>  eps = B*U     OR    >>  F = I+B*U
    
    return [B,detJ]




def q4_shapefn(xi,eta):
    
    shapefn = [.25*(1-xi)*(1-eta),
               .25*(1+xi)*(1-eta),
               .25*(1+xi)*(1+eta),
               .25*(1-xi)*(1+eta)]
    
    return shapefn




def q8_nabla(xi,eta,X):
    #
    # Derivatives of shape functions wrt xi & eta
    dNdxi = np.matrix([[1/4*eta+1/2*xi*(1-eta)-1/4*eta**2,
                       -1/4*eta+1/2*xi*(1-eta)+1/4*eta**2,
                        1/4*eta+1/4*eta**2+1/2*xi*(1+eta),
                       -1/4*eta+1/2*xi*(1+eta)-1/4*eta**2,
                       -xi*(1-eta),
                        1/2-1/2*eta**2,
                       -xi*(1+eta),
                       -1/2+1/2*eta**2],
    
                       [1/4*xi-1/4*xi**2+(1/2-1/2*xi)*eta,
                       -1/4*xi-1/4*xi**2+(1/2+1/2*xi)*eta,
                        1/4*xi+(1/2+1/2*xi)*eta+1/4*xi**2,
                       -1/4*xi+1/4*xi**2+(1/2-1/2*xi)*eta,
                       -1/2+1/2*xi**2,
                       -2*(1/2+1/2*xi)*eta,
                        1/2-1/2*xi**2,
                       -2*(1/2-1/2*xi)*eta]]);
    #
    J      = dNdxi*X
    detJ   = np.linalg.det(J)       # Determinant of Jacobian J
    invJ   = np.linalg.inv(J)       # Inverse of Jacobian
    #
    ## for a plain B matrix where B_ij= dN_i / dX_j:
    B = invJ*dNdxi #Shape function derivatives wrt x and y
    
    # B values may be inseted into various other places in so that 
    #   >>  eps = B*U     OR    >>  F = I+B*U
    return [B,detJ]



def q8_shapefn(xi,eta):
    
    shapefn = [-.25*(1-xi)*(1-eta)*(1+xi+eta),
               -.25*(1+xi)*(1-eta)*(1-xi+eta),
               -.25*(1+xi)*(1+eta)*(1-xi-eta),
               -.25*(1-xi)*(1+eta)*(1+xi-eta),
               0.5*(1-xi**2)*(1-eta),
               0.5*(1+xi)*(1-eta**2),
               0.5*(1-xi**2)*(1+eta),
               0.5*(1-xi)*(1-eta**2)]
    
    return shapefn



def Extrapolate_2x2(gaussvals,nnodes = 4):
    # extrapolate to 4, 8 or 9 nodes using (2x2 Gauss interpolation)
    Gmat = np.matrix([gaussvals[i] for i in range(1,5)])
    
    sq32 = 0.8660254037844386#np.sqrt(3)/2.
    
    Wmat = np.array([[1+sq32, -.5, 1-sq32, -.5],
                     [-.5, 1+sq32, -.5, 1-sq32],
                     [1-sq32, -.5, 1+sq32, -.5],
                     [-.5, 1-sq32, -.5, 1+sq32]])
    
    #
    if nnodes >4 :
        sq34 = 0.4330127018922193 #sqrt(3)/4

        Wmat = np.r_[Wmat,np.array([[.25+sq34, .25+sq34,.25-sq34, .25-sq34 ],
                                    [.25-sq34, .25+sq34,.25+sq34, .25-sq34 ],
                                    [.25-sq34, .25-sq34,.25+sq34, .25+sq34 ],
                                    [.25+sq34, .25-sq34,.25-sq34, .25+sq34]])]
        
    if nnodes == 9:
        Wmat = np.r_[Wmat,np.array([[.25, .25, .25, .25 ]])]
        
    #
    return Wmat*Gmat

    
    
    