import numpy as np


def pointPotential(x,y,q,posx,posy):
    """Returns the potential for a point particle [meters][coulombs]"""
    k=8.987551787*10**9
    Vxy=k*q/np.sqrt((x-posx)**2+(y-posy)**2)
    return Vxy

def dipolePotential(x,y,q,d_x,d_y):
    """Returns the electric potential for a dipole [meters][coulombs]"""
    Vxy=pointPotential(x,y,q,d_x/2.,d_y/2.)-pointPotential(x,y,q,-d_x/2.,-d_y/2.)
    return Vxy
