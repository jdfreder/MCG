import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from datetime import date
from urllib2 import urlopen

#Importing data from particle data group
Y=date.today().year
i=0
while i<=Y-2013:
    try:
        TotalData = urlopen('http://pdg.lbl.gov/'+str(Y-i)+'/hadronic-xsections/rpp'+str(Y-i)+'-pp_total.dat')
        print "Using "+str(Y-i)+" data for total cross section."
        DataFound1=True
        break
    except:
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

#Fitting a function to PDG's total data
def func(s,Z,M,Y1,Y2,n1,n2):
    m=.938
    sM=(2*m+M)**2
    hbar=6.58211928*10**-16
    c=3e8
    sigma=Z+np.pi*(hbar*c)**2/M**2*(np.log((s)**2/sM))**2+Y1*(sM/(s)**2)**n1-Y2*(sM/(s)**2)**n2
    return sigma

#Func2-4 define fits for PDG's elastic data
def func2(s,a,b):
    return a*np.e**(b*s)+.42

def func3(s,a,b):
    return np.log(a*s)+b

def func4(s,a,b):
    return np.e**(a*np.log(s)+b)

#Set up parameters to plot PDG data with fits
s=EPlab[:]
y=ESig[:]
p0=[5,2.15,12.72,7.35,.462,.55]
popt,pcov=curve_fit(func,s,y,p0)
sp=np.zeros(np.size(s)+1,float)
sp[:-1]=np.copy(s)
sp[-1]=200
sp2=np.arange(2e3,9e3,1000)
sp3=np.arange(8.6e3,6e8,1000)

s2=Plab[:]
y2=Sig[:]
p0=[34.71,2.15,12.72,7.35,.462,.550]
popt2,pcov2=curve_fit(func,s2,y2,p0)

g=np.arange(5,6*10**8,1000)
SigInel=func(g,popt2[0],popt2[1],popt2[2],popt2[3],popt2[4],popt2[5])-func(g,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5])

def SigI(BE):
    """Returns the cross-sectional area [fm^2] for given beam energy [GeV]"""
    if BE<4*10**3:
        return .1*(func(BE,popt2[0],popt2[1],popt2[2],popt2[3],popt2[4],popt2[5])-func(BE,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5]))
    else:
        return .1*(func(BE,popt2[0],popt2[1],popt2[2],popt2[3],popt2[4],popt2[5])-np.e**(9*10*-2*np.log(BE)+1.21))

def DisplayData():
    """Displays the Particle Data Group Data."""
    plt.loglog(Plab,Sig,ls=' ',marker='.',markersize=3,color='black')
    plt.loglog(EPlab,ESig,ls=' ',marker='.',markersize=3,color='black')
    plt.loglog(sp,func(sp,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5]),color='b',lw=1.5)
    plt.loglog(sp2,func2(sp2,6.566,.00001),lw=1.5)
    plt.loglog(sp3,func4(sp3,9*10**-2,1.21),color='b',lw=1.5)
    plt.loglog(s2[100:],func(s2[100:],popt2[0],popt2[1],popt2[2],popt2[3],popt2[4],popt2[5]),color='b',lw=1.5)
    plt.errorbar(Plab,Sig,xerr=[Plab-Plab_min,Plab_max-Plab],yerr=[StEr_L,StEr_H],ms=.5,mew=0,fmt=None,ecolor='black')
    plt.errorbar(EPlab,ESig,xerr=[EPlab-EPlab_min,EPlab_max-EPlab],yerr=[EStEr_L,EStEr_H],ms=.5,mew=0,fmt=None,ecolor='black')
    plt.annotate("$P_{lab}(GeV/c)$",fontsize=16,xy=(.1,1),xytext=(10e5,1.5))
    plt.annotate("total",fontsize=11,xy=(300,60),xytext=(300,60))
    plt.annotate("elastic",fontsize=11,xy=(300,10),xytext=(300,10))
    plt.annotate("pp",fontsize=11,xy=(70,15),xytext=(68,18))
    plt.ylabel("Cross Section (mb)",fontsize=12)
    plt.ylim(1,400)
    plt.grid(which='minor',axis='y')
    plt.grid(which='major',axis='x')

def WoodsSaxon(y,v):
    """Returns the Woods-Saxon Density profile for given element or atomic number"""
    A={'C':12,'O':16,'Al':27,'S':32,'Ca':40,'Ni':58,'Cu':63,'W':186,'Au':197,'Pb':208,'U':238}
    if type(y)==str:
        x=y
    else:
        for An,Val in A.iteritems():
            if Val==y:
                x=An
    a={'C':0,'O':.513,'Al':.519,'S':.61,'Ca':.586,'Ni':.516,'Cu':.596,'W':.535,'Au':.535,'Pb':.546,'U':.6}
    w={'C':0,'O':-.051,'Al':0,'S':0,'Ca':-.161,'Ni':-.1308,'Cu':0,'W':0,'Au':0,'Pb':0,'U':0}
    R={'C':2.47,'O':2.608,'Al':3.07,'S':3.458,'Ca':3.76,'Ni':4.309,'Cu':4.2,'W':6.51,'Au':6.38,'Pb':6.68,'U':6.68}
    r=np.arange(0,1.5*R[x]+.01,.01)
    if v==True:
        Rho=4*np.pi*r**2*(1+w[x]*(r**2)/(R[x]**2))/(1+np.e**((r-R[x])/a[x]))
    else:
        Rho=(1+w[x]*(r**2)/(R[x]**2))/(1+np.e**((r-R[x])/a[x]))
    return Rho

def distribute1D(x,prob,N):
    """Takes any numerical distribution prob, on the interval defined by the array x, and returns an array of N sampled values that are statistically the same as the input data."""
    y=prob
    A=np.cumsum(y)/(np.cumsum(y)[(len(x)-1)])
    z=np.random.random_sample(N)
    B=np.searchsorted(A,z)
    return x[B]

def Collider(N,Particle,Energy,v):
    """
    Simulates N collisions given specific Element and Center of Mass Energy [GeV]. 
    Returns the matrices that correspond to center-to-center seperation distance 
    (Random value between 0 and 150% of the Nucleus diameter [fm]), 
    Nuclei 1 and 2, number participating nucleons, and the number of binary 
    collisions. Additionally returns the interactions distance of the nucleons 
    given the choosen beam energy.
    """
    Plab=((Energy**2/(2*.938)-.938)**2-.938**2)**(1/2.) #Convert CM energy to corresponding Plab
    
    #wood saxon parameters
    w={'C':0,'O':-.051,'Al':0,'S':0,'Ca':-.161,'Ni':-.1308,'Cu':0,'W':0,'Au':0,'Pb':0,'U':0}
    A={'C':12,'O':16,'Al':27,'S':32,'Ca':40,'Ni':58,'Cu':63,'W':186,'Au':197,'Pb':208,'U':238}
    R={'C':2.47,'O':2.608,'Al':3.07,'S':3.458,'Ca':3.76,'Ni':4.309,'Cu':4.2,'W':6.51,'Au':6.38,'Pb':6.68,'U':6.68}
    a={'C':0,'O':.513,'Al':.519,'S':.61,'Ca':.586,'Ni':.516,'Cu':.596,'W':.535,'Au':.535,'Pb':.546,'U':.6}
    #Set Rp equal to the number of nucleons in a nucleus
    Rp=R[Particle] 
    b=2*Rp*np.random.random_sample(N) #Create array of random impact parameter
    r=np.arange(0,1.5*Rp+.01,.01) #Array for radial data used for plotting
    Npart=np.zeros(N,float)
    Ncoll=np.zeros(N,float)
    Sig1=SigI(Plab) #Find cross-sectional area from Plab from PDG fit functions. See function SigI.
    Maxr=(float(SigI(Plab))/np.pi)**.5 #Radius within which 2 nucleons will interact
    #Runs N number of times; each run creates a new array containing useful parameters, Ncoll, Npart, etc.
    for L in range(N):
        Nucleus1=np.zeros((A[Particle],2),float)
        Nucleus2=np.zeros((A[Particle],2),float)
        Nucleus1[:,0]=distribute1D(r,WoodsSaxon(Particle,v),A[Particle])[:] #Create nuclei from woods saxon function
        Nucleus2[:,0]=distribute1D(r,WoodsSaxon(Particle,v),A[Particle])[:]
        for i in range(A[Particle]):
            Nucleus1[i,1]=2*np.pi*np.random.random_sample(1) #give each nucleon in both nuclei a random azimuthal angle
            Nucleus2[i,1]=2*np.pi*np.random.random_sample(1) #from 0 to 2*pi
        Colide1=np.copy(Nucleus1) #Creates a copy of the original nuclei
        Colide2=np.copy(Nucleus2) #Colide arrays are used for plotting only
        for p1 in range(A[Particle]): #Loops through all nucleons in Nucleus 1
            count=0 #Arbitrary count variable used to determine whether a nucleon escaped the collision unwounded
            for p2 in range(A[Particle]): #Loops through all nucleons in Nucleus 2
                #x distance between two nucleons is: b+r1*cos(a1)-r2*cos(a2) where b is the random impact parameter, r is the radial distance from the center of the nucleus, and a is the azimuthal angle
                #y distance between two nucleons is r1*sin(a1)-r2*sin(a2)
                #distance between two nucleons is sqrt(x^2+y^2). If statements check to see if this value is smaller than Maxr
                if ((b[L]+Nucleus2[p2,0]*np.cos(Nucleus2[p2,1])-Nucleus1[p1,0]*np.cos(Nucleus1[p1,1]))**2+(Nucleus2[p2,0]*np.sin(Nucleus2[p2,1])-Nucleus1[p1,0]*np.sin(Nucleus1[p1,1]))**2)**.5 <= Maxr/2.:
                    Ncoll[L]+=1
                    #If distance between nucleons is smaller than maxr, nucleons interact and 1 collision is counted.
                if ((b[L]+Nucleus2[p2,0]*np.cos(Nucleus2[p2,1])-Nucleus1[p1,0]*np.cos(Nucleus1[p1,1]))**2+(Nucleus2[p2,0]*np.sin(Nucleus2[p2,1])-Nucleus1[p1,0]*np.sin(Nucleus1[p1,1]))**2)**.5 > Maxr/2.:
                    count+=1
                    #If distance is greater than maxr, +1 is added to count
                if count==A[Particle]:
                    #If count reaches the number of total nucleons in a nucleus, we infer that
                    #that nucleon has not hit any other nucleons in the other nucleus. Therefore
                    #this nucleon is tagged in the colide array.
                    Colide1[p1]=1000
        for p2 in range(A[Particle]):
            count=0
            for p1 in range(A[Particle]):
                if ((b[L]+Nucleus2[p2,0]*np.cos(Nucleus2[p2,1])-Nucleus1[p1,0]*np.cos(Nucleus1[p1,1]))**2+(Nucleus2[p2,0]*np.sin(Nucleus2[p2,1])-Nucleus1[p1,0]*np.sin(Nucleus1[p1,1]))**2)**.5 > Maxr:
                    count+=1
                if count==A[Particle]:
                    Colide2[p2]=1000
        #any nucleons not tagged in the colide array have therefore participated in the collisoin and are added to Npart
        for i in Colide1[:,0]:
            if i<1000:
                Npart[L]+=1
        for i in Colide2[:,0]:
            if i<1000:
                Npart[L]+=1
    return b,Nucleus1,Nucleus2,Npart,Ncoll,Maxr,Colide1,Colide2

def PlotNuclei(Nucleus1,Nucleus2,Particle,v):
    """
    Plots a histogram showing the relation between radial distance and 
    the number of nucleons from each nucleus. 
    Blue corresponds to nucleus 1 and green to nucleus 2.
    """
    w={'C':0,'O':-.051,'Al':0,'S':0,'Ca':-.161,'Ni':-.1308,'Cu':0,'W':0,'Au':0,'Pb':0,'U':0}
    A={'C':12,'O':16,'Al':27,'S':32,'Ca':40,'Ni':58,'Cu':63,'W':186,'Au':197,'Pb':208,'U':238}
    R={'C':2.47,'O':2.608,'Al':3.07,'S':3.458,'Ca':3.76,'Ni':4.309,'Cu':4.2,'W':6.51,'Au':6.38,'Pb':6.68,'U':6.68}
    a={'C':0,'O':.513,'Al':.519,'S':.61,'Ca':.586,'Ni':.516,'Cu':.596,'W':.535,'Au':.535,'Pb':.546,'U':.6}
    Rp=R[Particle]
    r=np.arange(0,1.5*Rp+.01,.01)
    n1,bins,patches=plt.hist(Nucleus1[:,0],15,normed=True,alpha=1)
    n2,bins,patches=plt.hist(Nucleus2[:,0],15,normed=True,alpha=.7)
    if max(n1)<=max(n2):
        plt.plot(r,WoodsSaxon(Particle,v)*max(n2)/max(WoodsSaxon(Particle,v)),lw=2.5)
    else:
        plt.plot(r,WoodsSaxon(Particle,v)*max(n1)/max(WoodsSaxon(Particle,v)),lw=2.5)
    plt.xlabel("Radial Distance [fm]",fontsize=14)
    plt.ylabel("Density",fontsize=14)

def PlotCollision(N,Particle,ImpactDistance,Nucleus1,Nucleus2,Participants,BinaryCollisions,Maxr,Col1,Col2):
    """Plots a cross-sectional view of the colliding particles."""
    b=ImpactDistance
    Npart=Participants
    Ncoll=BinaryCollisions
    w={'C':0,'O':-.051,'Al':0,'S':0,'Ca':-.161,'Ni':-.1308,'Cu':0,'W':0,'Au':0,'Pb':0,'U':0}
    A={'C':12,'O':16,'Al':27,'S':32,'Ca':40,'Ni':58,'Cu':63,'W':186,'Au':197,'Pb':208,'U':238}
    R={'C':2.47,'O':2.608,'Al':3.07,'S':3.458,'Ca':3.76,'Ni':4.309,'Cu':4.2,'W':6.51,'Au':6.38,'Pb':6.68,'U':6.68}
    a={'C':0,'O':.513,'Al':.519,'S':.61,'Ca':.586,'Ni':.516,'Cu':.596,'W':.535,'Au':.535,'Pb':.546,'U':.6}
    Rp=R[Particle]
    N1=plt.Circle((Rp,Rp),Rp,color='b',fill=False,lw=2)
    N2=plt.Circle((Rp+b[N-1],Rp),Rp,color='g',fill=False,lw=2)
    fig = plt.gcf()
    ax = plt.gca()
    #Any nucleons tagged in the colide array of the COLLIDE FUNCTION are not plotted in red because they have not interacted with any other nucleons
    ax.plot(Rp+Nucleus1[:,0]*np.cos(Nucleus1[:,1]),Rp+Nucleus1[:,0]*np.sin(Nucleus1[:,1]),'b.',ms=26,alpha=.8,mew=1,mec='blue')
    ax.plot(Rp+b[N-1]+Nucleus2[:,0]*np.cos(Nucleus2[:,1]),Rp+Nucleus2[:,0]*np.sin(Nucleus2[:,1]),'g.',ms=26,alpha=.8,mew=1,mec='green')
    ax.plot(Rp+Col1[:,0]*np.cos(Col1[:,1]),Rp+Col1[:,0]*np.sin(Col1[:,1]),'r.',ms=26,alpha=.4,mew=1,mec='red')
    ax.plot(Rp+b[N-1]+Col2[:,0]*np.cos(Col2[:,1]),Rp+Col2[:,0]*np.sin(Col2[:,1]),'r.',ms=26,alpha=.4,mew=1,mec='red')
    ax.annotate('Npart='+str(Npart[N-1])+'\nNcoll='+str(Ncoll[N-1]),xy=(1,0),xytext=(-3,26.5),fontsize=16)
    ax.annotate('Maxr: '+str(Maxr)+' fm',xy=(0,2*Rp),xytext=(-3,24.5),fontsize=12)
    fig.gca().add_artist(N1)
    fig.gca().add_artist(N2)
    plt.plot([-2.5,Maxr-2.5],[24,24],color='r',ls='-',lw=3)
    plt.xlim((-4,26))
    plt.ylim((-4,26))
    plt.xlabel('Horizontal Position [fm]',fontsize=15)
    plt.ylabel('Vertical Position [fm]',fontsize=15)
    fig.set_size_inches(6,6)

def PlotResults(b,Npart,Ncoll):
    """Plots the number or wounded nucleons and binary collisions as a function of impact parameter."""
    plt.scatter(b,Npart,label='Participants',color='b',alpha=.6)
    plt.scatter(b,Ncoll,label='Collisions',color='r',alpha=.6)
    plt.xlabel('Separation distance [fm]',fontsize=14)
    plt.xlim((0,max(b)))
    plt.ylim((0,max(Ncoll)*1.1))
    plt.legend()