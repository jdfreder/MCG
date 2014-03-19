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

def Ecm(Plab):
    """Converts Plab momenta to center of mass energy [GeV]."""
    E=(((Plab**2+.938**2)**(1/2.)+.938)**2-(Plab**2))**(1/2.)
    return E
if DataFound1==True and DataFound2==True:
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

#Apply best fit curve to the elastic data
s=eE_cm[:]
y=ESig[:]
p0=[4.45,.0965,2.127,11,4,.55,.55]
popt,pcov=curve_fit(func,s,y,p0)

#Apply best fit curve to total data
s2=E_cm[90:]
y2=Sig[90:]
p0=[34.49,.2704,2.127,12.98,7.38,.451,.549]
popt2,pcov2=curve_fit(func,s2,y2,p0)

def SigI(BE):
    """Returns the cross-sectional area [fm^2] for given beam energy [GeV]"""
    return .1*(func(BE,popt2[0],popt2[1],popt2[2],popt2[3],popt2[4],popt2[5],popt2[6])-func(BE,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5],popt[6]))

def DisplayData():
    """Displays the Particle Data Group Data."""
    #figsize(10,7)
    plt.loglog(E_cm,Sig,ls=' ',marker='.',markersize=3,color='black')
    plt.loglog(eE_cm,ESig,ls=' ',marker='.',markersize=3,color='black')
    plt.loglog(E_cm[90:],func(E_cm[90:],popt2[0],popt2[1],popt2[2],popt2[3],popt2[4],popt2[5],popt2[6]))
    plt.loglog(E_cm,func(E_cm,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5],popt[6]))
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

parameters=np.loadtxt("WoodSaxonParameters.txt",dtype='string',delimiter='\t')
pNucleus=parameters[:,0]
pModel=parameters[:,1]
pr2=parameters[:,2]
pC_A=parameters[:,3]
pZ_Alpha=parameters[:,4]
pw=parameters[:,5]

def NCD(Nucleus='197Au',Model='2pF',Range=2,Bins=100):
    """Returns the Nuclear Charge Distribution for a given Nucleus with specified
    model. Creates radial distribution from 0 to Range*nuclear radius with n 
    number of bins. If no values are set, defaults to 197Au using two-parameter
    Fermi model up to twice the nuclear radius with 100 bins."""
    #For multiple models of the same nucleus takes the first set of parameters and notifies the user which parameters are used.
    for index in range(len(pNucleus)):
        if pNucleus[index]==Nucleus and pModel[index]==Model:
            i=index
    r=np.linspace(0,Range*float(pr2[i]),Bins)
    if Model=='HO':
        return (1+float(pZ_Alpha[i])*(r/float(pC_A[i]))**2)*np.exp(-1*(r/float(pC_A[i]))**2)
    elif Model=='MHO':
        print "Model not yet supported"
    elif Model=='Mi':
        print "Model not yet supported"
    elif Model=='FB':
        print "Model not yet supported"
    elif Model=='SOG':
        print "Model not yet supported"
    elif Model=='2pF':
        return 1/(1+np.exp((r-float(pC_A[i]))/float(pZ_Alpha[i])))
    elif Model=='3pF':
        return (1+float(pw[i])*r**2/float(pC_A[i])**2)/(1+np.exp((r-float(pC_A[i]))/float(pZ_Alpha[i])))
    elif Model=='3pG':
        return (1+float(pw[i])*r**2/float(pC_A[i])**2)/(1+np.exp((r**2-float(pC_A[i])**2)/float(pZ_Alpha[i])**2))
    elif Model=='UG':
        print "Model not yet supported"
    else:
        print 'Model not found'

def distribute1D(x,prob,N):
    """Takes any numerical distribution prob, on the interval defined by the array x, and returns an array of N sampled values that are statistically the same as the input data."""
    y=prob*4*np.pi*x**2
    A=np.cumsum(y)/(np.cumsum(y)[(len(x)-1)])
    z=np.random.random_sample(N)
    B=np.searchsorted(A,z)
    return x[B]

def Collider(N,Particle,A,Energy,model,Range,Bins):
    """
    Simulates N collisions given specific Element and Center of Mass Energy [GeV]. 
    Returns the matrices that correspond to center-to-center seperation distance 
    (Random value between 0 and 200% of the Nucleus diameter [fm]), 
    Nuclei 1 and 2, number participating nucleons, and the number of binary 
    collisions. Additionally returns the interactions distance of the nucleons 
    given the choosen beam energy.
    """
    
    #Set Rp equal to the radius of a nucleus
    for index in range(len(pNucleus)):
        if pNucleus[index]==Particle and pModel[index]==model:
            i=index
    Rp=float(pr2[i])
    b=2*Rp*np.random.random_sample(N) #Create array of random impact parameter
    r=np.linspace(0,2*Rp,Bins) #Array for radial data used for plotting
    Npart=np.zeros(N,float)
    Ncoll=np.zeros(N,float)
    Maxr=(float(SigI(Energy))/np.pi)**.5 #Radius within which 2 nucleons will interact
    #Runs N number of times; each run creates a new array containing useful parameters, Ncoll, Npart, etc.
    for L in range(N):
        Nucleus1=np.zeros((A,3),float)
        Nucleus2=np.zeros((A,3),float)
        Nucleus1[:,0]=distribute1D(r,NCD(Particle,model,Range,Bins),A)[:]
        Nucleus2[:,0]=distribute1D(r,NCD(Particle,model,Range,Bins),A)[:]
        for i in range(A):
            Nucleus1[i,1]=2*np.pi*np.random.random_sample(1) #give each nucleon in both nuclei a random azimuthal angle
            Nucleus2[i,1]=2*np.pi*np.random.random_sample(1) #from 0 to 2*pi
            Nucleus1[i,2]=2*np.pi*np.random.random_sample(1)-np.pi
            Nucleus2[i,2]=2*np.pi*np.random.random_sample(1)-np.pi

    FailSafe=0
    for p1 in range(A):
        for p1x in range(A):
            if p1x==p1:
                skip1=1
            else:
                while np.sqrt((Nucleus1[p1,0]*np.cos(Nucleus1[p1,2])*np.cos(Nucleus1[p1,1])-Nucleus1[p1x,0]*np.cos(Nucleus1[p1x,2])*np.cos(Nucleus1[p1x,1]))**2+(Nucleus1[p1,0]*np.cos(Nucleus1[p1,2])*np.sin(Nucleus1[p1,1])-Nucleus1[p1x,0]*np.cos(Nucleus1[p1x,2])*np.sin(Nucleus1[p1x,1]))**2+(Nucleus1[p1,0]*np.sin(Nucleus1[p1,2])-Nucleus1[p1x,0]*np.sin(Nucleus1[p1x,2]))**2)<Maxr:
                    Nucleus1[p1x,0]=distribute1D(r,NCD(Particle,model,Range,Bins),A)[p1x]
                    Nucleus1[p1x,1]=2*np.pi*np.random.random_sample(1)
                    Nucleus1[p1x,2]=2*np.pi*np.random.random_sample(1)-np.pi
                    FailSafe+=1
                    if FailSafe>=1000:
                        break

    FailSafe=0
    for p2 in range(A):
        for p2x in range(A):
            if p2x==p2:
                skip2=1
            else:
                while np.sqrt((Nucleus2[p2,0]*np.cos(Nucleus2[p2,2])*np.cos(Nucleus2[p2,1])-Nucleus2[p2x,0]*np.cos(Nucleus2[p2x,2])*np.cos(Nucleus2[p2x,1]))**2+(Nucleus2[p2,0]*np.cos(Nucleus2[p2,2])*np.sin(Nucleus2[p2,1])-Nucleus2[p2x,0]*np.cos(Nucleus2[p2x,2])*np.sin(Nucleus2[p2x,1]))**2+(Nucleus2[p2,0]*np.sin(Nucleus2[p2,2])-Nucleus2[p2x,0]*np.sin(Nucleus2[p2x,2]))**2)<Maxr:
                    Nucleus2[p2x,0]=distribute1D(r,NCD(Particle,model,Range,Bins),A)[p2x]
                    Nucleus2[p2x,1]=2*np.pi*np.random.random_sample(1)
                    Nucleus2[p2x,2]=2*np.pi*np.random.random_sample(1)-np.pi
                    FailSafe+=1
                    if FailSafe>=1000:
                        break

        Colide1=np.copy(Nucleus1) #Creates a copy of the original nuclei
        Colide2=np.copy(Nucleus2) #Colide arrays are used for plotting only
        for p1 in range(A): #Loops through all nucleons in Nucleus 1
            count=0 #Arbitrary count variable used to determine whether a nucleon escaped the collision unwounded
            for p2 in range(A): #Loops through all nucleons in Nucleus 2
                if ((b[L]+Nucleus2[p2,0]*np.cos(Nucleus2[p2,1])*np.cos(Nucleus2[p2,2])-Nucleus1[p1,0]*np.cos(Nucleus1[p1,1])*np.cos(Nucleus1[p1,2]))**2+(Nucleus2[p2,0]*np.sin(Nucleus2[p2,1])*np.cos(Nucleus2[p2,2])-Nucleus1[p1,0]*np.sin(Nucleus1[p1,1])*np.cos(Nucleus1[p1,2]))**2)**.5 <= Maxr:
                    Ncoll[L]+=1
                    #If distance between nucleons is smaller than maxr, nucleons interact and 1 collision is counted.
                if ((b[L]+Nucleus2[p2,0]*np.cos(Nucleus2[p2,1])*np.cos(Nucleus2[p2,2])-Nucleus1[p1,0]*np.cos(Nucleus1[p1,1])*np.cos(Nucleus1[p1,2]))**2+(Nucleus2[p2,0]*np.sin(Nucleus2[p2,1])*np.cos(Nucleus2[p2,2])-Nucleus1[p1,0]*np.sin(Nucleus1[p1,1])*np.cos(Nucleus1[p1,2]))**2)**.5 > Maxr:
                    count+=1
                    #If distance is greater than maxr, +1 is added to count
                if count==A:
                    #If count reaches the number of total nucleons in a nucleus, we infer that
                    #that nucleon has not hit any other nucleons in the other nucleus. Therefore
                    #this nucleon is tagged in the colide array.
                    Colide1[p1,0]=100
                    Colide1[p1,1]=0
                    Colide1[p1,2]=0
        for p2 in range(A):
            count=0
            for p1 in range(A):
                if ((b[L]+Nucleus2[p2,0]*np.cos(Nucleus2[p2,1])*np.cos(Nucleus2[p2,2])-Nucleus1[p1,0]*np.cos(Nucleus1[p1,1])*np.cos(Nucleus1[p1,2]))**2+(Nucleus2[p2,0]*np.sin(Nucleus2[p2,1])*np.cos(Nucleus2[p2,2])-Nucleus1[p1,0]*np.sin(Nucleus1[p1,1])*np.cos(Nucleus1[p1,2]))**2)**.5 > Maxr:
                    count+=1
                if count==A:
                    Colide2[p2,0]=100
                    Colide2[p2,1]=0
                    Colide2[p2,2]=0
        #any nucleons not tagged in the colide array have therefore participated in the collisoin and are added to Npart
        for i in Colide1[:,0]:
            if i<100:
                Npart[L]+=1
        for i in Colide2[:,0]:
            if i<100:
                Npart[L]+=1
    return b,Nucleus1,Nucleus2,Npart,Ncoll,Maxr,Colide1,Colide2

def PlotNuclei(Nucleus1,Nucleus2,Particle,model,Range,Bins):
    """
    Plots a histogram showing the relation between radial distance and 
    the number of nucleons from each nucleus. 
    Blue corresponds to nucleus 1 and green to nucleus 2.
    """
    for index in range(len(pNucleus)):
        if pNucleus[index]==Particle and pModel[index]==model:
            i=index
    Rp=float(pr2[i])
    r=np.linspace(0,Range*Rp,Bins)
    n1,bins,patches=plt.hist(Nucleus1[:,0],40,normed=True,alpha=1)
    n2,bins,patches=plt.hist(Nucleus2[:,0],40,normed=True,alpha=.75)
    if max(n1)<=max(n2):
        plt.plot(r,NCD(Particle,model,Range,Bins)*max(n2)/max(NCD(Particle,model,Range,Bins)),lw=2.5)
    else:
        plt.plot(r,NCD(Particle,model,Range,Bins)*max(n1)/max(NCD(Particle,model,Range,Bins)),lw=2.5)
    plt.xlabel("Radial Distance [fm]",fontsize=14)
    plt.ylabel("Density",fontsize=14)

def PlotCollisionSlice(N,Particle,model,ImpactDistance,Nucleus1,Nucleus2,Participants,BinaryCollisions,Maxr,Col1,Col2):
    """Plots a cross-sectional view of the colliding particles."""
    b=ImpactDistance
    Npart=Participants
    Ncoll=BinaryCollisions
    for index in range(len(pNucleus)):
        if pNucleus[index]==Particle and pModel[index]==model:
            i=index
    Rp=float(pr2[i])
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

def PlotCollisionHorizontal(N,Particle,ImpactDistance,Nucleus1,Nucleus2,Participants,BinaryCollisions,Colide1,Colide2):
    """Plots a horizontal view of the colliding particles."""
    b=ImpactDistance
    Npart=Participants
    Ncoll=BinaryCollisions
    for index in range(len(pNucleus)):
        if pNucleus[index]==Particle and pModel[index]==model:
            i=index
    Rp=float(pr2[i])
    N1=plt.Circle((Rp,Rp),Rp,color='b',fill=False,lw=2)
    N2=plt.Circle((Rp*4,Rp+b[N-1]),Rp,color='g',fill=False,lw=2)
    fig = plt.gcf()
    ax = plt.gca()
    ax.plot(Rp+Nucleus1[:,0]*np.sin(Nucleus1[:,1]),Rp+Nucleus1[:,0]*np.cos(Nucleus1[:,1])*np.cos(Nucleus1[:,2]),'b.',ms=26,alpha=.8,mew=1,mec='blue')
    ax.plot(4*Rp+Nucleus2[:,0]*np.sin(Nucleus2[:,1]),b[N-1]+Rp+Nucleus2[:,0]*np.cos(Nucleus2[:,1])*np.cos(Nucleus2[:,2]),'g.',ms=26,alpha=.8,mew=1,mec='green')
    ax.plot(Rp+Colide1[:,0]*np.sin(Colide1[:,1]),Rp+Colide1[:,0]*np.cos(Colide1[:,1])*np.cos(Colide1[:,2]),'r.',ms=26,alpha=.4,mew=1,mec='red')
    ax.plot(4*Rp+Colide2[:,0]*np.sin(Colide2[:,1]),b[N-1]+Rp+Colide2[:,0]*np.cos(Colide2[:,1])*np.cos(Colide2[:,2]),'y.',ms=26,alpha=.4,mew=1,mec='yellow')
    ax.annotate("",xy=(2.5*Rp,Rp), xycoords='data',xytext=(2*Rp,Rp), textcoords='data',arrowprops=dict(arrowstyle='-|>',connectionstyle="arc"))
    ax.annotate("",xy=(2.5*Rp,b[N-1]+Rp), xycoords='data',xytext=(3*Rp, b[N-1]+Rp), textcoords='data',arrowprops=dict(arrowstyle='-|>',connectionstyle="arc"))
    fig.gca().add_artist(N1)
    fig.gca().add_artist(N2)
    zed=5*Rp
    plt.xlim((0,zed))
    plt.ylim((0,zed))
    plt.xlabel('Horizontal Position [fm]',fontsize=15)
    plt.ylabel('Vertical Position [fm]',fontsize=15)
    fig.set_size_inches(6,6)

def PlotResults(b,Npart,Ncoll):
    """Plots the number or wounded nucleons and binary collisions as a function of impact parameter."""
    HData=zeros((len(b),2))
    PData=zeros((len(b),2))
    HData[:,0]=b[:]
    PData[:,0]=b[:]
    HData[:,1]=Ncoll[:]
    PData[:,1]=Npart[:]
    SortHData=HData[HData[:,0].argsort()]
    SortPData=PData[PData[:,0].argsort()]
    AvgNColl=zeros(floor(max(SortHData[:,0]))+1)
    AvgNPart=zeros(floor(max(SortPData[:,0]))+1)
    CBins=zeros(len(AvgNColl))
    PBins=zeros(len(AvgNPart))
    xColl=linspace(0,len(AvgNColl)-1,len(AvgNColl))
    xPart=linspace(0,len(AvgNPart)-1,len(AvgNPart))
    lc=0
    lp=0
    for n in SortHData[:,0]:
        AvgNColl[floor(n)]+=SortHData[lc,1]
        CBins[floor(n)]+=1
        lc+=1
    for i in SortPData[:,0]:
        AvgNPart[floor(i)]+=SortPData[lp,1]
        PBins[floor(i)]+=1
        lp+=1
    for p in range(len(CBins)):
        AvgNColl[p]=AvgNColl[p]/CBins[p]
        AvgNPart[p]=AvgNPart[p]/PBins[p]
    plt.bar(xColl,AvgNColl,width=1,alpha=.25,color='red')
    plt.bar(xPart,AvgNPart,width=1,alpha=.5,color='blue')
    plt.scatter(b,Npart,label='Participants',color='b',alpha=.6)
    plt.scatter(b,Ncoll,label='Collisions',color='r',alpha=.6)
    plt.xlabel('Separation distance [fm]',fontsize=14)
    plt.xlim((0,max(b)))
    plt.ylim((0,max(Ncoll)*1.1))
    plt.legend()