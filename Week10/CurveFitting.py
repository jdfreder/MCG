import numpy as np  

def LinearLeastSquaresFit(x,y):
    """Take in arrays represting (x,y) values for a set of linearly varying data and perform a linear least squares regression. Return the resulting slope and intercept parameters of the best fit line with their uncertainties."""
    avX=1/np.float(len(x))*np.sum(x)
    avY=1/np.float(len(y))*np.sum(y)
    avX2=1/np.float(len(x))*np.sum(x**2)
    avXY=1/np.float(len(x))*np.sum(x*y)
    delI=np.zeros(x.shape,float)
    slope=(avXY-avX*avY)/(avX2-avX**2)
    intercept=(avX2*avY-avX*avXY)/(avX2-avX**2)
    delI[:]=y[:]-(slope*x[:]+intercept)
    Del2=sum(delI**2)
    slerr=(1/(float(len(x))-2)*(Del2/float(len(x)))/(avX2-avX**2))**.5
    interr=(1/(float(len(x))-2)*(Del2/float(len(x)))*avX2/(avX2-avX**2))**.5
    return slope,intercept,slerr,interr
