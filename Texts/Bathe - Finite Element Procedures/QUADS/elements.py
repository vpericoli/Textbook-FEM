import numpy as np


def QUADS(nel, itype, nint, thic, ym, pr, xx, s, iout):
    #
    #   Input Variables:
    #       nel     =   number of element 
    #       itype   =   element type   
    #                       0 = axisymmetric 
    #                       1 = plane strain 
    #                       2 = plane stress 
    #       nint    =   gauss numerical integration order 
    #       thic    =   thickness of element
    #       ym      =   young's modulus 
    #       pr      =   poisson's ratio
    #       xx[2,4] =   element node coords
    #       iout    =   unit number used for output
    #
    #   Output Variables:
    #       s[8,8]  =   calculated stiffness matrix
    #
    
    s = np.zeros(8, 8) 
    
    #
    #   Obtain Stress-Strain Law
    #
    
    #
    #   Plane Strain Analysis
    #
    
    #  
    #   Axisymmetric Analysis 
    #
    
    #
    #   Plane Stress 
    #
    
    #
    #   Calculate element stiffness 
    #
    
    #
    #   Evaluate derivative operator B and the Jacobian determinant DET
    #
    
    stdm(xx, B, det, ri, si, xbar, nel, itype, out)
    
    #
    #  Add contribution to element stiffness 
    #
    
    return s
    
def stdm(xx, B, det, R, S, xbar, nel, itype, iout):
    #
    #   interpolation functions 
    #
    
    #
    #   natural coordinate derivatives of the interpolation functions 
    #
    #       1. with respect to R 
    #       2. with respect to S 
    #
    
    #
    #   Evaluate the Jacobian matrix at point (r, s)
    #
    
    #
    #   Compute the determinant of the jacobian matrix at point (r, s)
    #
    
    #
    #   Compute inverse of the jacobian matrix 
    #
    
    #
    #   Evaluate global derivative operator S 
    #
    
    #
    #   In case of plane strain or plane stress analysis do not include
    #   the normal strain component.
    #
    
    #
    #   Compute the radius at point (r,s)
    #
    
    #
    #   Evaluate the hoop strain-displacement relation
    #
    
    #
    #   For the case of zero radius requate radial to hoop strain. 
    #
    
    #
    #   Non-zero radius
    #
    
    return 0
    
    
    
    
    
    
    
    