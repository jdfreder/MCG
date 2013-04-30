import numpy as np

def finiteDifference(x,y):
    """Returns the centered finite difference derivative of y with respect to x"""
    dydx=np.zeros(y.shape,float)
    dydx[1:-1]=(y[2:]-y[:-2])/(x[2:]-x[:-2])
    dydx[0]=(y[1]-y[0])/(x[1]-x[0])
    dydx[-1]=(y[-1]-y[-2])/(x[-1]-x[-2])
    return dydx

def fourPtFiniteDiff(x,y):
    """Returns the four point finite difference derivative of y with respect to x"""
    dydx=np.zeros(y.shape,float)
    dydx[0]=(y[1]-y[0])/(x[1]-x[0])
    dydx[1]=(y[2]-y[1])/((x[2] - x[1]))
    dydx[-1]=(y[-1]-y[-2])/(x[-1]-x[-2])
    dydx[-2]=(y[-2]-y[-3])/((x[-2] - x[-3]))
    for i, l in enumerate(y[2:-2]):
        i+=2
        dydx[i] = (y[i-2]-8*y[i-1]+8*y[i+1]-y[i+2])/(12*(x[i+1]-x[i]))
    return dydx
