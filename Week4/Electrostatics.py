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

def pointField(x,y,q,Xq,Yq):
    """Returns a tuple of the electric field components from a charge q[C] at a position (Xq,Yq)"""
    k=8.987551787*10**9
    Ex=k*q*(x-Xq)/((x-Xq)**2+(y-Yq)**2)**(1./2.)
    Ey=k*q*(y-Yq)/((x-Xq)**2+(y-Yq)**2)**(1./2.)
    return (Ex,Ey)
