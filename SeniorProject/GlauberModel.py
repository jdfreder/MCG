import IPython.core.pylabtools as pyt
import numpy as np
import matplotlib.pyplot as plt
from datetime import date
from urllib2 import urlopen
from scipy.optimize import curve_fit
from scipy.special import j0
import time

#Importing data from particle data group
#Attempts to use data from current year
#If that data is not available, drops down a year until data is found or defaults to 2013 data
Y=date.today().year
i=0
while i<=Y-2013:
    try:
        TotalData = urlopen('http://pdg.lbl.gov/'+str(Y-i)+'/hadronic-xsections/rpp'+str(Y-i)+'-pp_total.dat')
        print "Using "+str(Y-i)+" data for total cross section."
        DataFound1=True
        break
    except:
        print str(Y-i)+" total cross section data is unavailable. The Particle Data Group website may not have the latest data or may have changed format."
        i+=1
        if i>Y-2013:
            print "---\nData not found. Please check your internet connection to http://pdg.lbl.gov/2013/html/computer_read.html\n---"
            DataFound1=False
l=0
while l<=Y-2013:
    try:
        ElasticData = urlopen('http://pdg.lbl.gov/'+str(Y-l)+'/hadronic-xsections/rpp'+str(Y-l)+'-pp_elastic.dat')
        print "Using "+str(Y-l)+" data for elastic cross section."
        DataFound2=True
        break
    except:
        l+=1
        if l>Y-2013:
            print "---\nData not found. Please check your internet connection to http://pdg.lbl.gov/2013/html/computer_read.html\n---"
            DataFound2=False

if DataFound1==True:
    data=np.loadtxt(TotalData,float,usecols=(0,1,2,3,4,5,6,7,8),skiprows=11)
    Point=data[:,0]
    Plab=data[:,1] #GeV/c
    Plab_min=data[:,2]
    Plab_max=data[:,3]
    Sig=data[:,4]
    StEr_H=data[:,5]
    StEr_L=data[:,6]
    SyEr_H=data[:,7]
    SyEr_L=data[:,8]
if DataFound2==True:
    Edata=np.loadtxt(ElasticData,float,usecols=(0,1,2,3,4,5,6,7,8),skiprows=11)
    EPoint=Edata[:,0]
    EPlab=Edata[:,1] #GeV/c
    EPlab_min=Edata[:,2]
    EPlab_max=Edata[:,3]
    ESig=Edata[:,4]
    EStEr_H=Edata[:,5]
    EStEr_L=Edata[:,6]
    ESyEr_H=Edata[:,7]
    ESyEr_L=Edata[:,8]

def Ecm(Plab):
    """Converts Plab momenta to center of mass energy [GeV]."""
    E=(((Plab**2+.938**2)**(1/2.)+.938)**2-(Plab**2))**(1/2.)
    return E
if DataFound1==True and DataFound2==True:
    #Automatically converts all P_lab momenta to corresponding center-of-mass energy [GeV]
    E_cm=Ecm(Plab)
    eE_cm=Ecm(EPlab)
    cm_min=Ecm(Plab_min)
    cm_max=Ecm(Plab_max)
    ecm_min=Ecm(EPlab_min)
    ecm_max=Ecm(EPlab_max)

#Define best fit curve given by the particle data group
def func(s,P,H,M,R1,R2,n1,n2):
    m=.93827 #Proton mass GeV/c^2
    sM=(2*m+M)**2 #Mass^2 (GeV/c^2)^2
    hbar=6.58211928*10**-25 #GeV*s
    c=2.99792458*10**8 #m/s
    sigma=H*(np.log(s**2/sM))**2+P+R1*(s**2/sM)**(-n1)-R2*(s**2/sM)**(-n2)
    return sigma

#Apply best fit curve to the elastic cross-section data
s=eE_cm[:]
y=ESig[:]
p0=[4.45,.0965,2.127,11,4,.55,.55]
popt,pcov=curve_fit(func,s,y,p0)

#Apply best fit curve to total cross-section data
s2=E_cm[90:]
y2=Sig[90:]
p0=[34.49,.2704,2.127,12.98,7.38,.451,.549]
popt2,pcov2=curve_fit(func,s2,y2,p0)

def SigI(BE):
    """Returns the proton-proton cross-sectional area [fm^2] for given beam energy [GeV]"""
    return .1*(func(BE,popt2[0],popt2[1],popt2[2],popt2[3],popt2[4],popt2[5],popt2[6])-func(BE,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5],popt[6]))

def DisplayData():
    """Displays the Proton-Proton Cross Section Data."""
    pyt.figsize(12,7)
    plt.loglog(E_cm,Sig,ls=' ',marker='.',markersize=3,color='black',label='PDG Data')
    plt.loglog(eE_cm,ESig,ls=' ',marker='.',markersize=3,color='black')
    plt.loglog(E_cm[90:],func(E_cm[90:],popt2[0],popt2[1],popt2[2],popt2[3],popt2[4],popt2[5],popt2[6]),color='blue')
    plt.loglog(E_cm,func(E_cm,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5],popt[6]),color='blue',label='Fit')
    plt.scatter(7000,70.5,label='TOTEM EPL 101 21003',color='red')
    plt.scatter([2760,7000],[62.1,72.7],label='ALICE 2011',color='blue')
    plt.scatter([7000,8000],[72.9,74.7],label='TOTEM 2013',color='green')
    plt.errorbar([2760,7000,7000,7000,8000],[62.1,70.5,72.7,72.9,74.7],yerr=[5.9,3.4,6.2,1.5,1.7],fmt=' ',color='black')
    plt.loglog(E_cm[90:],10*SigI(E_cm[90:]))
    plt.errorbar(E_cm,Sig,xerr=[E_cm-cm_min,cm_max-E_cm],yerr=[StEr_L,StEr_H],ms=.5,mew=0,fmt=None,ecolor='black')
    plt.errorbar(eE_cm,ESig,xerr=[eE_cm-ecm_min,ecm_max-eE_cm],yerr=[EStEr_L,EStEr_H],ms=.5,mew=0,fmt=None,ecolor='black')
    plt.annotate("Total",fontsize=11,xy=(7,46),xytext=(7,46))
    plt.annotate("Elastic",fontsize=11,xy=(1000,10),xytext=(1000,10))
    plt.annotate("Inelastic",fontsize=11,xy=(35,25),xytext=(35,25))
    plt.title("pp Cross Section Data",fontsize=16)
    plt.ylabel("Cross Section [mb]",fontsize=12)
    plt.xlabel("$\sqrt{s}\,[GeV]$",fontsize=16)
    plt.ylim(1,400)
    plt.grid(which='minor',axis='y')
    plt.grid(which='major',axis='x')
    plt.legend(loc=4)
    plt.show()

#Reads in parameters to calculate nuclear charge densities (NCD)
parameters=np.loadtxt("WoodSaxonParameters.txt",dtype='string',delimiter='\t')
FBdata=np.loadtxt("FourierBesselParameters.txt",str,delimiter='\t')
pNucleus=parameters[:,0]
pModel=parameters[:,1]
pr2=parameters[:,2]
pC_A=parameters[:,3]
pZ_Alpha=parameters[:,4]
pw=parameters[:,5]
FBnucleus=FBdata[:,0]
FBrms=FBdata[:,1]
FBR=FBdata[:,2]
FBa=np.zeros((len(FBnucleus)-1,17),float)
for i in range(len(FBnucleus)-1):
    FBa[i,:]=FBdata[i+1,3:]
FBa = FBa.astype(np.float)

def NCD(Nucleus,Model,Range=2,Bins=100):
    """Returns the Nuclear Charge Distribution for a given Nucleus with specified
    model. Creates radial distribution from 0 to Range*nuclear radius with n 
    number of bins. If no values are set, defaults to 197Au using two-parameter
    Fermi model up to twice the nuclear radius with 100 bins."""
    #For multiple models of the same nucleus takes the first set of parameters and notifies the user which parameters are used.
    j=[]
    for index in range(len(pNucleus)):
        if pNucleus[index]==Nucleus and pModel[index]==Model:
            j.append(index)
            i=index
    j=np.array(j,dtype=int)
    if len(j)>1:
        #print "Multiple parameters detected for given model. Using primary values."
        i=j[0]
    r=np.linspace(0,Range*float(pr2[i]),Bins)
    if Model=='HO':
        return (1+float(pZ_Alpha[i])*(r/float(pC_A[i]))**2)*np.exp(-1*(r/float(pC_A[i]))**2)
    elif Model=='MHO':
        print "Warning: Model not yet supported\nPlease use a different model."
        return None
    elif Model=='Mi':
        print "Warning: Model not yet supported\nPlease use a different model."
        return None
    elif Model=='FB':
        #print "Warning: Fourier-Bessel Model currently contains support for He-3 only. If not simulating He-3 collisions, please choose another model."
        for FBindex in range(len(FBnucleus)):
            if FBnucleus[FBindex]==Nucleus:
                iFB=FBindex
        r=np.linspace(0,float(FBR[iFB]),Bins)
        p=np.zeros(np.size(r),float)
        v=np.arange(0,17,1)
        for i in range(len(r)):
            p[i]=abs(sum(FBa[iFB-1,v]*j0((v+1)*np.pi*r[i]/float(FBR[iFB]))))
        return p
    elif Model=='SOG':
        print "Warning: Model not yet supported\nPlease use a different model."
        return None
    elif Model=='2pF':
        return 1/(1+np.exp((r-float(pC_A[i]))/float(pZ_Alpha[i])))
    elif Model=='3pF':
        return (1+float(pw[i])*r**2/float(pC_A[i])**2)/(1+np.exp((r-float(pC_A[i]))/float(pZ_Alpha[i])))
    elif Model=='3pG':
        return (1+float(pw[i])*r**2/float(pC_A[i])**2)/(1+np.exp((r**2-float(pC_A[i])**2)/float(pZ_Alpha[i])**2))
    elif Model=='UG':
        print "Warning: Model not yet supported\nPlease use a different model."
        return None
    else:
        print 'Error: Model not found\nPlease check that the model was typed in correctly. (Case Sensitive)'
        return None

def distribute1D(x,prob,N):
    """Takes any numerical distribution probability, on 
    the interval defined by the array x, and returns an 
    array of N sampled values that are statistically the 
    same as the input data."""
    y=prob*4*np.pi*x**2
    A=np.cumsum(y)/(np.cumsum(y)[(len(x)-1)])
    z=np.random.random_sample(N)
    B=np.searchsorted(A,z)
    return x[B]

def Collider(N,Particle1,A1,Particle2,A2,model1,model2,Energy,bRange=1.1,Range=2,Bins=100):
    """
    Simulates N collisions between specified Elements (with number of nucleons A)
    and using Center of Mass Energy [GeV]. 
    Returns the matrices corresponding to center-to-center seperation distance 
    (Random value between 0 and bRange*(radius of nucleus1 + radius of nucleus2) [fm]), 
    Nuclei 1 and 2, number of participating nucleons, and the number of binary 
    collisions. Additionally returns the interaction distance of the nucleons 
    given the choosen beam energy and the radii of both nuclei.
    """
    #Set Rp1 and Rp2 equal to the radii of the nuclei choosen
    j1=[]
    j2=[]
    i1="Unassigned"
    i2="Unassigned"
    for index in range(len(pNucleus)):
        if pNucleus[index]==Particle1 and pModel[index]==model1:
            j1.append(index)
            i1=index
        if pNucleus[index]==Particle2 and pModel[index]==model2:
            j2.append(index)
            i2=index
    j1=np.array(j1,dtype=int)
    j2=np.array(j2,dtype=int)
    if len(j1)>1 or len(j2)>1:
        print "Multiple parameters detected for specified model. Using primary values."
        i1=j1[0]
        i2=j2[0]
    if i1=="Unassigned" or i2=="Unassigned":
        print 'Error: Model not found\nPlease check that the model was typed in correctly. (Case Sensitive)'
        return None,None,None,None,None,None,None,None
    if model1 != 'FB':
        Rp1=float(pr2[i1])
    else:
        for FBindex in range(len(FBnucleus)):
            if FBnucleus[FBindex]==Particle1:
                iFB=FBindex
        Rp1=float(FBR[iFB])
    if model2 != 'FB':
        Rp2=float(pr2[i2])
    else:
        for FBindex in range(len(FBnucleus)):
            if FBnucleus[FBindex]==Particle2:
                iFB=FBindex
        Rp2=float(FBR[iFB])
    b=(Rp1+Rp2)*bRange*np.random.random_sample(N) #Create array of random impact parameters
    if model1 != 'FB':
        r1=np.linspace(0,Range*Rp1,Bins) #Array of radial data used for plotting
    else:
        r1=np.linspace(0,Rp1,Bins)
    if model2 != 'FB':
        r2=np.linspace(0,Range*Rp2,Bins)
    else:
        r2=np.linspace(0,Rp2,Bins)
    Npart=np.zeros(N,float)
    Ncoll=np.zeros(N,float)
    Maxr=np.sqrt(SigI(Energy)/np.pi) #Radius within which two nucleons will interact
    #Runs N number of times; each run creates a new array containing the parameters of interest: Ncoll, Npart, etc.
    for L in range(N):
        Nucleus1=np.zeros((A1,7),float)
        Nucleus2=np.zeros((A2,7),float)
        #Gives each nucleon is own radial distance from the center of the nucleus
        #Sampled from the NCD function and distributed with a factor 4*pi*r^2*p(r)
        Nucleus1[:,0]=distribute1D(r1,NCD(Particle1,model1,Range,Bins),A1)[:]
        Nucleus2[:,0]=distribute1D(r2,NCD(Particle2,model2,Range,Bins),A2)[:]
        #Nucleons are then given random azimuthal distances such that no particular point is more likely to be populated than another
        #Cartesian coordinates (x,y,z) are determined from spherical coordinates and passed to the nuclei arrays
        for i in range(A1):
            Nucleus1[i,1]=np.arccos(2*np.random.random_sample(1)-1)
            Nucleus1[i,2]=2*np.pi*np.random.random_sample(1)
            Nucleus1[i,3]=Nucleus1[i,0]*np.sin(Nucleus1[i,1])*np.cos(Nucleus1[i,2])
            Nucleus1[i,4]=Nucleus1[i,0]*np.sin(Nucleus1[i,1])*np.sin(Nucleus1[i,2])
            Nucleus1[i,5]=Nucleus1[i,0]*np.cos(Nucleus1[i,1])
        for i in range(A2):
            Nucleus2[i,1]=np.arccos(2*np.random.random_sample(1)-1)
            Nucleus2[i,2]=2*np.pi*np.random.random_sample(1)
            Nucleus2[i,3]=Nucleus2[i,0]*np.sin(Nucleus2[i,1])*np.cos(Nucleus2[i,2])
            Nucleus2[i,4]=Nucleus2[i,0]*np.sin(Nucleus2[i,1])*np.sin(Nucleus2[i,2])
            Nucleus2[i,5]=Nucleus2[i,0]*np.cos(Nucleus2[i,1])
        for p1 in range(A1):
            for p1x in range(A1):
                FailSafe=0 #Prevents program from running indefinitely in some cases
                if p1x==p1:
                    pass
                else:
                     while np.sqrt((Nucleus1[p1,3]-Nucleus1[p1x,3])**2+(Nucleus1[p1,4]+Nucleus1[p1x,4])**2+(Nucleus1[p1,5]+Nucleus1[p1x,5])**2)<Maxr:
                        Nucleus1[p1x,1]=np.arccos(2*np.random.random_sample(1)-1)
                        Nucleus1[p1x,2]=2*np.pi*np.random.random_sample(1)
                        Nucleus1[p1x,3]=Nucleus1[p1x,0]*np.sin(Nucleus1[p1x,1])*np.cos(Nucleus1[p1x,2])
                        Nucleus1[p1x,4]=Nucleus1[p1x,0]*np.sin(Nucleus1[p1x,1])*np.sin(Nucleus1[p1x,2])
                        Nucleus1[p1x,5]=Nucleus1[p1x,0]*np.cos(Nucleus1[p1x,1])
                        FailSafe+=1
                        if FailSafe>10:
                            Nucleus1[p1x,0]=distribute1D(r1,NCD(Particle1,model1,Range,Bins),A1)[p1x]
        for p2 in range(A2):
            for p2x in range(A2):
                FailSafe=0
                if p2x==p2:
                    pass
                else:
                    while np.sqrt((Nucleus2[p2,3]-Nucleus2[p2x,3])**2+(Nucleus2[p2,4]-Nucleus2[p2x,4])**2+(Nucleus2[p2,5]-Nucleus2[p2x,5])**2)<Maxr:
                        Nucleus2[p2x,1]=np.arccos(2*np.random.random_sample(1)-1)
                        Nucleus2[p2x,2]=2*np.pi*np.random.random_sample(1)
                        Nucleus2[p2x,3]=Nucleus2[p2x,0]*np.sin(Nucleus2[p2x,1])*np.cos(Nucleus2[p2x,2])
                        Nucleus2[p2x,4]=Nucleus2[p2x,0]*np.sin(Nucleus2[p2x,1])*np.sin(Nucleus2[p2x,2])
                        Nucleus2[p2x,5]=Nucleus2[p2x,0]*np.cos(Nucleus2[p2x,1])
                        FailSafe+=1
                        if FailSafe>10:
                            Nucleus2[p2x,0]=distribute1D(r2,NCD(Particle2,model2,Range,Bins),A2)[p2x]
        for p1 in range(A1): #Loops through all nucleons in Nucleus 1
            count=0 #Arbitrary count variable used to determine whether a nucleon escaped the collision unwounded
            for p2 in range(A2): #Loops through all nucleons in Nucleus 2
                if np.sqrt((Nucleus1[p1,3]-Nucleus2[p2,3])**2+(b[L]+Nucleus2[p2,4]-Nucleus1[p1,4])**2) <= Maxr:
                    Ncoll[L]+=1
                    #If distance between nucleons is smaller than maxr, nucleons interact and 1 collision is counted.
                if np.sqrt((Nucleus1[p1,3]-Nucleus2[p2,3])**2+(b[L]+Nucleus2[p2,4]-Nucleus1[p1,4])**2) > Maxr:
                    count+=1
                    #If distance is greater than maxr, +1 is added to count
                if count==A2:
                    #If count reaches the number of total nucleons in nucleus 2, we can infer that
                    #the nucleon in nucleus 1 has not hit any other nucleons in the other nucleus. Therefore
                    #this nucleon is flagged as an unwounded particle
                    Nucleus1[p1,6]=1
        for p2 in range(A2):
            count=0
            for p1 in range(A1):
                if np.sqrt((Nucleus1[p1,3]-Nucleus2[p2,3])**2+(b[L]+Nucleus2[p2,4]-Nucleus1[p1,4])**2) > Maxr:
                    count+=1
                if count==A1:
                    Nucleus2[p2,6]=1
        #Number of participating particles is the total number of particles minus all the flagged unwounded particles
        Npart[L]=A1+A2-(sum(Nucleus1[:,6])+sum(Nucleus2[:,6]))
    return b,Nucleus1,Nucleus2,Npart,Ncoll,Maxr,Rp1,Rp2

def PlotNuclei(Nucleus1,Nucleus2,Particle1,Particle2,model1,model2,Rp1,Rp2,Range,Bins):
    """
    Plots the nuclear charge density for each nucleus and the 
    cummulative distribution of the nucleons in each nucleus.
    Blue corresponds to nucleus 1 and green to nucleus 2.
    """
    if model1 != 'FB':
        Range1=Range
    else:
        Range1=1
    if model2 != 'FB':
        Range2=Range
    else:
        Range2=1
    r1=np.linspace(0,Range1*Rp1,Bins)
    r2=np.linspace(0,Range2*Rp2,Bins)
    #n1,bins,patches=plt.hist(Nucleus1[:,0],40,normed=True,alpha=.75,label=str(Particle1)+ " nucleons")
    #n2,bins,patches=plt.hist(Nucleus2[:,0],40,normed=True,alpha=.75,label=str(Particle2)+" nucleons")
    plt.plot(r1,NCD(Particle1,model1,Range,Bins)/max(NCD(Particle1,model1,Range,Bins)),lw=2.5,label=str(Particle1)+" Radial Density")
    plt.plot(r2,NCD(Particle2,model2,Range,Bins)/max(NCD(Particle2,model2,Range,Bins)),lw=2.5,label=str(Particle2)+" Radial Density")
    #if max(n1)>max(n2):
        #plt.plot(r1,NCD(Particle1,model1,Range,Bins)*max(n1),lw=2.5,label=str(Particle1)+" Radial Density")
        #plt.plot(r2,NCD(Particle2,model2,Range,Bins)*max(n1),lw=2.5,label=str(Particle2)+" Radial Density")
    #else:
        #plt.plot(r1,NCD(Particle1,model1,Range,Bins)*max(n2),lw=2.5,label=str(Particle1)+" Radial Density")
        #plt.plot(r2,NCD(Particle2,model2,Range,Bins)*max(n2),lw=2.5,label=str(Particle2)+" Radial Density")
    plt.xlabel("Radial Distance [fm]",fontsize=14)
    plt.ylabel("Density",fontsize=14)
    plt.legend()
    plt.show()

def ShowCollision(N,Particle1,A1,Particle2,A2,Rp1,Rp2,Nucleus1,Nucleus2,b,Npart,Ncoll,Maxr):
    """Plots a cross-sectional and horizontal view of the last collision."""
    fig,ax=plt.subplots()
    N1=plt.Circle((Rp1,Rp1),Rp1,color='b',fill=False,lw=2)
    N2=plt.Circle((Rp1+b[N-1],Rp1),Rp2,color='g',fill=False,lw=2)
    fig.gca().add_artist(N1)
    fig.gca().add_artist(N2)
    for i in range(A1):
        if Nucleus1[i,6]==1:
            ax.plot(Rp1+Nucleus1[i,4],Rp1+Nucleus1[i,3],'b.',ms=26,alpha=.6,mew=0,mec='blue')
        if Nucleus1[i,6]==0:
            ax.plot(Rp1+Nucleus1[i,4],Rp1+Nucleus1[i,3],'r.',ms=26,alpha=.6,mew=0,mec='red')
    for i in range(A2):
        if Nucleus2[i,6]==1:
            ax.plot(b[N-1]+Rp1+Nucleus2[i,4],Rp1+Nucleus2[i,3],'g.',ms=26,alpha=.6,mew=0,mec='green')
        if Nucleus2[i,6]==0:
            ax.plot(b[N-1]+Rp1+Nucleus2[i,4],Rp1+Nucleus2[i,3],'y.',ms=26,alpha=.6,mew=0,mec='yellow')
    zed=1.2*(Rp1+Rp2)+b[N-1]
    ax.annotate('Npart='+str(Npart[N-1])+'\nNcoll='+str(Ncoll[N-1]),xy=(1,0),xytext=(0,1.015*zed),fontsize=16)
    ax.annotate('Maxr: '+str(Maxr)[:5]+' fm',xy=(0,2*Rp1),xytext=(.01*zed,.95*zed),fontsize=12)
    ax.plot([(.01*zed),(.01*zed)+Maxr],[zed*.93,zed*.93],color='r',ls='-',lw=3)
    plt.xlim((0,zed))
    plt.ylim((0,zed))
    plt.xlabel('Horizontal Cross Section [fm]',fontsize=15)
    plt.ylabel('Vertical Position [fm]',fontsize=15)
    fig.set_size_inches(6,6)
    fig1,ax1=plt.subplots()
    N3=plt.Circle((Rp1,Rp1),Rp1,color='b',fill=False,lw=2)
    N4=plt.Circle(((Rp1+Rp2)*2,Rp1+b[N-1]),Rp2,color='g',fill=False,lw=2)
    for i in range(A1):
        if Nucleus1[i,6]==1:
            ax1.plot(Rp1+Nucleus1[i,5],Rp1+Nucleus1[i,4],'b.',ms=26,alpha=.6,mew=0,mec='blue')
        if Nucleus1[i,6]==0:
            ax1.plot(Rp1+Nucleus1[i,5],Rp1+Nucleus1[i,4],'r.',ms=26,alpha=.6,mew=0,mec='red')
    for i in range(A2):
        if Nucleus2[i,6]==1:
            ax1.plot(2*(Rp1+Rp2)+Nucleus2[i,5],b[N-1]+Rp1+Nucleus2[i,4],'g.',ms=26,alpha=.6,mew=0,mec='green')
        if Nucleus2[i,6]==0:
            ax1.plot(2*(Rp1+Rp2)+Nucleus2[i,5],b[N-1]+Rp1+Nucleus2[i,4],'y.',ms=26,alpha=.6,mew=0,mec='yellow')
    ax1.annotate("",xy=(2*Rp1+Rp2,Rp1), xycoords='data',xytext=(2*Rp1,Rp1), textcoords='data',arrowprops=dict(arrowstyle='-|>',connectionstyle="arc"))
    ax1.annotate("",xy=(2*Rp1,b[N-1]+Rp1), xycoords='data',xytext=(2*(Rp1+Rp2)-Rp2, b[N-1]+Rp1), textcoords='data',arrowprops=dict(arrowstyle='-|>',connectionstyle="arc"))
    fig1.gca().add_artist(N3)
    fig1.gca().add_artist(N4)
    zed=2.5*(Rp1+Rp2)
    plt.xlim((0,zed))
    plt.ylim((0,zed))
    plt.ylabel('Vertical Position [fm]',fontsize=15)
    plt.xlabel('Horizontal Position [fm]',fontsize=15)
    fig1.set_size_inches(6,6)
    plt.show()

def PlotResults(b,Npart,Ncoll,Particle1,Particle2,N,Energy,bins=10):
    """Plots number of collisions and participants as a function of impact parameter. 
    Shows average trend over data using specified number of bins."""
    xmin=0
    xmax=max(b)
    x=b
    y=Ncoll
    j = bins
    E = np.zeros(j)
    H = np.zeros(j)
    L = np.zeros(j)#,int64)
    newx = np.linspace(xmin,xmax,j)
    binwidth = (xmax-xmin)/float(j)
    #Shift by half a bin so the values plot at the right location?  
    #If the bins are small enough or the function not too steep, this doesn't matter
    plotx = newx - 0.5*binwidth 
    for i in range(len(x)):
        #find the array element in newx where the value from x belongs
        val = x[i]
        bin = newx.searchsorted(val)
        L[bin] += 1
        E[bin] += y[i]**2
        H[bin] += y[i]
    h = H/L
    spr = np.sqrt(E/L - h**2)
    err = spr/np.sqrt(L)
    y2=Npart
    E2 = np.zeros(j)
    H2 = np.zeros(j)
    L2 = np.zeros(j)#,int64)
    newx2 = np.linspace(xmin,xmax,j)
    binwidth2 = (xmax-xmin)/float(j)
    plotx2 = newx2 - 0.5*binwidth2
    for i in range(len(x)):
        val2 = x[i]
        bin = newx2.searchsorted(val2)
        L2[bin] += 1
        E2[bin] += y2[i]**2
        H2[bin] += y2[i]
    h2 = H2/L2
    spr2 = np.sqrt(E2/L2 - h2**2)
    err2 = spr2/np.sqrt(L2)
    plt.plot(x,y,"ro",alpha=.9,label='Ncoll')
    plt.plot(x,y2,"bo",alpha=.5,label='Npart')
    plt.plot(plotx2,h2,"g-",linewidth=4,label='Avg Npart')
    plt.plot(plotx,h,"y-",linewidth=4,label='Avg Ncoll')
    plt.xlim(0,max(x))
    plt.ylim(0,1.1*max(y))
    plt.legend()
    plt.ylabel('Npart / Ncoll')
    plt.xlabel('Impact parameter [fm]')
    plt.title(str(Particle1)+' + '+str(Particle2)+'. '+str(N)+' iterations. '+str(Energy)+' center-of-mass energy [GeV].')
    plt.show()