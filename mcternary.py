import numpy as np
import matplotlib.pyplot as plt

def getLightMass(avemL, sigmL, A):
    '''
    The mass of the light fragment is sampled 
    from a Gaussian distribution with an 
    average input value of avemL sigma of sigL. 
    AvemL, sigmL are obtained from reference x TBD.
    '''
    #Eventual code goes here
    mL= np.random.normal(avemL, sigmL, 1)
    #print np.mean(mL)
    #print np.std(mL, ddof=1)
    
    '''
    count, bins, ignored = plt.hist(mL, 50, normed=True)
    plt.plot(bins, 1/(sigmL * np.sqrt(2 * np.pi)) *
         np.exp( - (bins - avemL)**2 / (2 * sigmL**2) ),
         linewidth=2, color='r')
    plt.show()
    '''
    #mandZandA[0]=mL
    
    return mL#mandZandA


def getfragmentZandA(mL, actA):
    '''
     Let the charge of the light and heavy 
     fragment be expressed as variable name 
     ZL and ZH. Using the mass from step 1, 
     the AL is computed by converting MeV to 
     amu. Since the alpha particle has Zalpha 
     of 2, use the value of AL to compute ZL, 
     ZH, and AH assuming an equal mass charge 
     distribution, that is, the ratio of ZL/AL 
     is equal to ZH/AH.
    '''
     #Eventual code goes here
    ZL=37 #e C
    ZH=90-ZL #e C
    Zalpha=2 #e C
    AL=95 #amu
    AH=232-AL #amu
    Aalpha=4 #amu
    malpha=4 #amu
    
    temp=np.zeros(7)
    temp[0]=malpha
    temp[1]=ZL
    temp[2]=ZH
    temp[3]=Zalpha
    temp[4]=AL
    temp[5]=AH
    temp[6]=Aalpha
    return temp


def getmH(mZA):
    #Semi-empirical mass formula constant terms
    av = 15.5 #MeV
    as1 = 17.2 #MeV
    ac = 0.72 #MeV
    asym = 23.2 #MeV
    ap = 12 #MeV
    delta = 0.
    if mZA[7]%2==0 and mZA[4]%2==0:
        delta = ap/np.sqrt(mZA[7])
    if mZA[7]%2==0 and mZA[4]%2!=0:
        delta = -ap/np.sqrt(mZA[7])
    #print("delta:", delta, "pairing term")

    B = av*mZA[7] - as1*mZA[7]**(2./3.) - (ac*mZA[4]*(mZA[4] - 1.0))/(mZA[7]**(1./3.)) - (asym*(mZA[7] - 2.0*mZA[4])**2)/(mZA[7]) + delta
    
    return B


def getE_alpha_scissandE_L_scissandE_H_sciss(EKB_act_infinity, EKB_act_sciss, EKT_act_infinity, deltaE_act, EKT_act_sciss, E_alpha_infinity):
    '''
    To compute the energy available for the three 
    particles coming out of the ternary fission, 
    we use the published numbers for the TKE of 
    binary fission minus the typical energy of the 
    LRA to figure out how much each fragment can have. 
    The variables and their initial values for U236 are: 
    EKB_act_infty = 168 +/- 4.5 MeV,
    EKB_act_sciss = 18.1 MeV, 
    EKT_act_infty = 155.5 +/- 0.8 MeV. 
    At scission, the alpha-accompanied fission
    has the initial values of: DeltaE_act = 4-5 MeV, 
    and EKT_act_sciss = 13 MeV. Furthermore, average
    kinetic energy of alpha particle at infinity 
    measured to be: E_alpha_infty = 16 MeV.
    
    These values have been obtained from Negele et al. 1978,
    and Asghar et al. 1970
    
    '''
    #Eventual code goes here
    E_alpha_sciss=128 #MeV
    E_L_sciss=256 #MeV
    E_H_sciss=512 #MeV
    
    temp=np.zeros(3)
    temp[0]=E_alpha_sciss
    temp[1]=E_L_sciss
    temp[2]=E_H_sciss
    return temp


def getd0andD0(mZA):
    '''
    Using initial value of ZL and ZH, calculate 
    distance d0 of the electrostatic saddle point 
    by using equation 21, relating it to separation 
    distance D0.
    '''
    #Eventual code goes here
    D0=21.9 #fm
    d0=D0/((mZA[4]/mZA[3])+1)
    
    temp=np.zeros(2)
    temp[0]=d0
    temp[1]=D0
    return temp


def getAvedandAveD(EKT_act_infinity, EKT_act_sciss, E_alpha_infinity,E_alpha_sciss,mandZandA):
    '''
    Use a guessed value for mean alpha KE 
    and equation 20 in order to calculate 
    average of d and D
    '''
    #Eventual code goes here
    aved=10 #fm
    aveD=22 #fm
    
    temp=np.zeros(2)
    temp[0]=aved
    temp[1]=aveD
    return temp


def getdandD(sigma_xp, dD, sigma_D):
    '''
    To account for uncertainty in separation 
    distance D, choose random values of D 
    from a Gaussian, equation 23 (with a width 
    of sigma_D = 1 fm). The uncertainty will 
    thus be given by summing a random variable
    (defined by equation 24b) with the value 
    of the separation distance obtained from 
    equation 20.  Likewise, as for the distance 
    d of the electrostatic saddle point, its 
    uncertainty is given by sampling d from a
    Gaussian distribution (equation 25), and 
    then by summing a random variable 
    (equation 26b) to d value from equation 20. 
    Where sigma_xp=0.93.
    '''
    
    d=np.random.normal(dD[2],sigma_xp, 1) #fm
    D=np.random.normal(dD[3], sigma_D, 1) #fm
    
    temp=np.zeros(2)
    temp[0]=d
    temp[1]=D
    return temp


def getxpandypandzpandRc(d, xpin, ypin, zpin,sigma_xp, sigma_yp, sigma_zp):
    '''
    To compute the uncertainty of the alpha 
    particle position along the primed axis, 
    use Gaussian distribution around d0, given 
    by equation 25, where xp=d-d0. One can 
    also use the same approach for figuring out 
    the uncertainty along yp and zp, while 
    assuming that sigma_yp=sigma_zp. Using 
    equation 28, one can then compute Rc from 
    equation 27.
    '''
    #Eventual code goes here
    xpout=1.2 #fm
    ypout=2.3 #fm
    zpout=2.2 #fm
    Rc=ypout #fm
    
    temp=np.zeros(3)
    temp[0]=xpout
    temp[1]=ypout
    temp[2]=zpout
    
    return temp


def getRc(paxisandpP):
    
    Rc=paxisandpP[1]
    
    return Rc


def getPxpandPypandPzp(Pxpin, Pypin, Pzpin, sigma_xp, sigma_yp, sigma_zp):
    '''
    Use the uncertainty principle to select 
    the initial momentum of the alpha particle 
    using a Gaussian distribution, employing 
    equations 30 and 31, where 
    sigma_xp = 0.93 fm 
    and sigma_yp = sigma_zp = 1.3 fm.
    
    The choice
    of sigma_xp is obtained from the assumption that an estimate
    of sigma_xp, may be of an order of magnitude of the
    root mean square of the radius of the alpha particle.
    The values of sigma_yp =sigma_zp, are chosen to be sqrt(2) sigma_p rather
    arbitrarily (larger and smaller values of  were
    used, and they seemed to give poorer agreement between
    theory and experiment).
    '''
    #Eventual code goes here
    Pxpout=2.4 #kg*m/s
    Pypout=4.6 #kg*m/s
    Pzpout=4.4 #kg*m/s
    
    temp=np.zeros(3)
    temp[0]=Pxpout
    temp[1]=Pypout
    temp[2]=Pzpout
    return temp


def getVL0andVH0(mandZandA):
    '''
    From conservation of momentum, explain what is 
    used to calculate them
    
    two default values given the numbers in the 
    paper come out to be
    VL0 = 0.013c
    and VH0 = 0.009c.
    '''
   
    M=236
    VH0 = 2*mandZandA[0]*M/(mandZandA[1]*(mandZandA[1] + 2*mandZandA[2])) #c m/s
    VL0 = (M-mandZandA[1]*VH0**2)/mandZandA[0] #c m/s
    
    temp=np.zeros(2)
    temp[0]=VL0
    temp[1]=VH0
    return temp


def getR0andV0andT0(mandZandA, DEpars):
    '''
    Depending on your fission fragments, 
    select the radius of each one from reference 
    x TBD and compute the average radii using 
    equation 11. Compute the average speed using 
    equation 12 and find T0 using equation 13.
    '''
    #Eventual code goes here
    R0=4.86 #fm
    V0=(DEpars[15]+DEpars[16])/2 #c m/s
    T0=V0/R0 #1/c m/s
    
    temp=np.zeros(3)
    temp[0]=R0
    temp[1]=V0
    temp[2]=T0
    return temp


def getalphapositionandvelocity(mandZandA, dD, DEpars):
    '''
     Using equations in Table 1, compute initial coordinates 
     and velocities of light, heavy and alpha 
     particles in the center of mass frame. The 
     velocity of alpha fragment is to be chosen 
     at random following a Gaussian distribution - with what parameter values? 
     In addition, velocities of the heavy and light 
     fragments are to be modified as to keep the 
     total momentum of the system as zero in the 
     cm frame. In other words, the values which the 
     L and H fragment can attain will be constricted 
     to whatever keeps total CM momentum zero
    '''
    #Eventual code goes here
    Vxalpha=0.6 #fm/s
    Vyalpha=-0.2 #fm/s
    Vzalpha=2.4 #fm/s
    
    mLH=mandZandA[0]+mandZandA[1]
    M=mandZandA[0]+mandZandA[1]+mandZandA[2]
    deltax=dD[4]-(mandZandA[1]/mLH)*dD[5]
    deltaV=(mandZandA[0]*DEpars[15]-mandZandA[1]*DEpars[16]+mandZandA[2]*Vxalpha)/mLH
    
    xalpha= -(mLH/M)*deltax #fm
    yalpha=(mLH/M)*dD[6] #fm
    zalpha=20. #fm  this is supposed to be 0
    
    
    temp=np.zeros(6)
    temp[0]=xalpha
    temp[1]=yalpha
    temp[2]=zalpha
    temp[3]=Vxalpha
    temp[4]=Vyalpha
    temp[5]=Vzalpha
    return temp


def getLpositionandvelocity(mandZandA, dD, DEpars, alphapars):
    '''
     Using equations in Table 1, compute initial coordinates 
     and velocities of light, heavy and alpha 
     particles in the center of mass frame. The 
     velocity of alpha fragment is to be chosen 
     at random following a Gaussian distribution - with what parameter values? 
     In addition, velocities of the heavy and light 
     fragments are to be modified as to keep the 
     total momentum of the system as zero in the 
     cm frame. In other words, the values which the 
     L and H fragment can attain will be constricted 
     to whatever keeps total CM momentum zero
    '''
    #Eventual code goes here
    
    mLH=mandZandA[0]+mandZandA[1]
    M=mandZandA[0]+mandZandA[1]+mandZandA[2]
    deltax=dD[4]-(mandZandA[1]/mLH)*dD[5]
    deltaV=(mandZandA[0]*DEpars[15]-mandZandA[1]*DEpars[16]+mandZandA[2]*alphapars[3])/mLH
    
    xL=(mandZandA[2]/M)*deltax+(mandZandA[1]/M)*dD[5] #fm
    yL= -(mandZandA[2]/M)*dD[6] #fm
    zL=0. #fm
    
    VxL=0.5 #fm/s
    VyL=-(mandZandA[2]/mLH)*alphapars[4] #fm/s
    VzL=-(mandZandA[2]/mLH)*alphapars[5] #fm/s
    
    
    temp=np.zeros(6)
    temp[0]=xL
    temp[1]=yL
    temp[2]=zL
    temp[3]=VxL
    temp[4]=VyL
    temp[5]=VzL
    return temp


def getHpositionandvelocity(mandZandA, dD, DEpars, alphapars, Lpars):
    '''
     Using equations in Table 1, compute initial coordinates 
     and velocities of light, heavy and alpha 
     particles in the center of mass frame. The 
     velocity of alpha fragment is to be chosen 
     at random following a Gaussian distribution - with what parameter values? 
     In addition, velocities of the heavy and light 
     fragments are to be modified as to keep the 
     total momentum of the system as zero in the 
     cm frame. In other words, the values which the 
     L and H fragment can attain will be constricted 
     to whatever keeps total CM momentum zero
    '''
    #Eventual code goes here
    
    mLH=mandZandA[0]+mandZandA[1]
    M=mandZandA[0]+mandZandA[1]+mandZandA[2]
    deltax=dD[4]-(mandZandA[1]/mLH)*dD[5]
    deltaV=(mandZandA[0]*DEpars[15]-mandZandA[1]*DEpars[16]+mandZandA[2]*alphapars[3])/mLH
    
    xH= -(mLH/M)*deltax-dD[5]+dD[4] #fm
    yH=Lpars[1] #fm
    zH=0. #fm 
    
    VxH=-0.3 #fm/s
    VyH=Lpars[4] #fm/s
    VzH=Lpars[5] #fm/s
    
    temp=np.zeros(6)
    temp[0]=xH
    temp[1]=yH
    temp[2]=zH
    temp[3]=VxH
    temp[4]=VyH
    temp[5]=VzH
    return temp


def gethatpars(mandZandA, DEpars,  Lpars, alphapars, t):
    '''
    Using results from equations 11-13, and 
    the initial positions/velocities found in 
    previous step, compute differential
    equation parameters found in equations 7-10.
    '''
    c=299792458 #m/s
    
    xhL=Lpars[0]/DEpars[12] #no unit
    xhalpha=alphapars[0]/DEpars[12] #no unit
    vhxL=Lpars[3]/DEpars[13] #no unit
    vhxalpha=alphapars[3]/DEpars[13] #no unit
    
    yhL=Lpars[1]/DEpars[12] #no unit
    yhalpha=alphapars[1]/DEpars[12] #no unit
    vhyL=Lpars[4]/DEpars[13] #no unit
    vhyalpha=alphapars[4]/DEpars[13] #no unit
    
    zhL=Lpars[2]/DEpars[12] #no unit
    zhalpha=alphapars[2]/DEpars[12] #no unit
    vhzL=Lpars[5]/DEpars[13] #no unit
    vhzalpha=alphapars[5]/DEpars[13] #no unit
    
    th=t/DEpars[14] #no unit
    beta0=DEpars[13]/c #no unit
    
    
    temp=np.zeros(14)
    temp[0]=xhL
    temp[1]=yhL
    temp[2]=zhL
    temp[3]=vhxL
    temp[4]=vhyL
    temp[5]=vhzL
    temp[6]=xhalpha
    temp[7]=yhalpha
    temp[8]=zhalpha
    temp[9]=vhxalpha
    temp[10]=vhyalpha
    temp[11]=vhzalpha
    temp[12]=th
    temp[13]=beta0
    return temp


def getDEpars(mandZandA, DEpars,  Lpars, alphapars, t, hatpars):
    '''
    Using results from equations 11-13, and 
    the initial positions/velocities found in 
    previous step, compute differential
    equation parameters found in equations 7-10.
    '''
  

    Ax=hatpars[0]-hatpars[6] #no unit
    Bx=(1+mandZandA[0]/mandZandA[1])*hatpars[0]+(mandZandA[2]/mandZandA[1])*hatpars[6] #no unit
    Cx=(mandZandA[0]/mandZandA[1])*hatpars[0]+(1+mandZandA[2]/mandZandA[0])*hatpars[6] #no unit
    
    Ay=hatpars[1]-hatpars[7] #no unit
    By=(1+mandZandA[0]/mandZandA[1])*hatpars[1]+(mandZandA[2]/mandZandA[1])*hatpars[7] #no unit
    Cy=(mandZandA[0]/mandZandA[1])*hatpars[1]+(1+mandZandA[2]/mandZandA[0])*hatpars[7] #no unit
    
    Az=hatpars[2]-hatpars[7] #no unit
    Bz=(1+mandZandA[0]/mandZandA[1])*hatpars[2]+(mandZandA[2]/mandZandA[1])*hatpars[7] #no unit
    Cz=(mandZandA[0]/mandZandA[1])*hatpars[2]+(1+mandZandA[2]/mandZandA[0])*hatpars[7] #no unit
    
    
    A=(Ax**2+Ay**2+Az**2)**(3/2) #no unit
    B=(Bx**2+By**2+Bz**2)**(3/2) #no unit
    C=(Cx**2+Cy**2+Cz**2)**(3/2) #no unit
    
    
    temp=np.zeros(12)
    temp[0]=A
    temp[1]=Ax
    temp[2]=Ay
    temp[3]=Az
    temp[4]=B
    temp[5]=Bx
    temp[6]=By
    temp[7]=Bz
    temp[8]=C
    temp[9]=Cx
    temp[10]=Cy
    temp[11]=Cz
    return temp


def getPalpha(alphapars, mandZandA):
    from math import pow
    pxalpha=mandZandA[2]*alphapars[3]
    pyalpha=mandZandA[2]*alphapars[4]
    pzalpha=mandZandA[2]*alphapars[5]
    Palpha=pow((pow(pxalpha, 2)+pow(pyalpha, 2)+pow(pzalpha, 2)), 0.5)
        
    temp=np.zeros(4)
    temp[0]=pxalpha
    temp[1]=pyalpha
    temp[2]=pzalpha
    temp[3]=Palpha
    return temp


def getdistanceL(alphapars, Lpars):
    from math import pow
    distanceL=pow((pow(alphapars[0]-Lpars[0], 2)+pow(alphapars[1]-Lpars[1], 2)+pow(alphapars[2]-Lpars[2], 2)), 0.5)
    return distanceL


def getdistanceH(alphapars, Hpars):
    from math import pow
    distanceH=pow((pow(alphapars[0]-Hpars[0], 2)+pow(alphapars[1]-Hpars[1], 2)+pow(alphapars[2]-Hpars[2], 2)), 0.5)
    return distanceH


def initparameters(alphapars, Lpars, Hpars, mandZandA, paxisandpP, dD, DEpars, hatpars, energies):
    avemL=27.
    sigmL=3.
    M=232
    actA=239.
    EKB_act_infinity=124
    EKB_act_sciss=114
    EKT_act_infinity=154
    deltaE_act=186
    EKT_act_sciss=287
    EKT_act_sciss=268
    E_alpha_infinity=298
    Zalpha=2
    sigma_xp=1.3
    sigma_D=1
    xpin=25
    ypin=12
    zpin=14
    sigma_yp=0.93
    sigma_zp=0.93
    Pxpin=13
    Pypin=23
    Pzpin=33
    inputs=50
    mH=137
    malpha=4
    zalpha=2
    A=236
    t=4
    
    mandZandA[0]=getLightMass(avemL, sigmL, A)
    mandZandA[2:9]=getfragmentZandA(mandZandA, actA)
    mandZandA[1]=getmH(mandZandA)
    energies=getE_alpha_scissandE_L_scissandE_H_sciss(EKB_act_infinity, EKB_act_sciss, EKT_act_infinity, deltaE_act, EKT_act_sciss, E_alpha_infinity)
    dD[2:4]=getd0andD0(mandZandA)
    dD[0:2]=getAvedandAveD(EKT_act_infinity, EKT_act_sciss, E_alpha_infinity, energies,mandZandA)
    dD[4:6]=getdandD(sigma_xp, dD, sigma_D)
    paxisandpP[0:3]=getxpandypandzpandRc(dD, xpin, ypin, zpin,sigma_xp, sigma_yp, sigma_zp)
    paxisandpP[3:6]=getPxpandPypandPzp(Pxpin, Pypin, Pzpin,sigma_xp, sigma_yp, sigma_zp)
    dD[6]=getRc(paxisandpP)
    DEpars[15:17]=getVL0andVH0(mandZandA)
    DEpars[12:15]=getR0andV0andT0(mandZandA, DEpars)
    alphapars[0:6]=getalphapositionandvelocity(mandZandA, dD, DEpars)
    Lpars[0:6]=getLpositionandvelocity(mandZandA, dD, DEpars, alphapars)
    Hpars[0:6]=getHpositionandvelocity(mandZandA, dD, DEpars, alphapars, Lpars)
    hatpars=gethatpars(mandZandA, DEpars,  Lpars, alphapars, t)
    DEpars[0:12]=getDEpars(mandZandA, DEpars,  Lpars, alphapars, t, hatpars)
    alphapars[6:10]=getPalpha(alphapars, mandZandA)
    distanceL=getdistanceL(alphapars, Lpars)
    distanceH=getdistanceH(alphapars, Hpars)

    
    return alphapars, Lpars, Hpars, mandZandA, paxisandpP, dD, DEpars, hatpars, energies


def writeresults(filename, results, label, mandZandA, alphapars, Lpars, Hpars):
    '''
    Once numerical solution for each equation 
    of motion has been completed, write file 
    to filename with particle trajectory data.
    '''
    
    #matrix=np.array([alphapars, Lpars, Hpars, mandZandA])
    writer=file('test.txt', 'a')
    np.savetxt(writer, (mandZandA, alphapars, Lpars, Hpars), fmt='%.3f', header=label)
    writer.close()
    #!cat test.txt
    #return matrix
    #Does this return anything?
    
    
def solveequationsofmotion(DEpars,hatpars,mandZandA, momentum, distanceL, distanceH):
    '''
    All of the preceding parameters are to 
    be passed in order to set up differential 
    equations of motion for the light and alpha 
    particles, equations 5 and 6.
    The coupled first-order differential 
    equations are to be solved numerically by 
    employing the Adams-Moulton predictor-corrector 
    method. This will compute the solutions to the 
    equations of motion (while the equation of motion 
    for the heavy fragment can be obtained via 
    conservation of momentum).
    '''
    mome=momentum-1
    distanceL=distanceL-1
    distanceH=distanceH-1
    print mome
    print distanceL
    print distanceH
            
    #Eventual code goes here
    #Return results  -- WHAT FORM ARE THE RESULTS IN?
    
    return mome, distanceL, distanceH
    #return xalpha, yalpha, zalpha, pxalpha, pyalpha, pzalpha, xL, yL, zL, xH, yH, zH
    
    
def m(A, E, alphapars, Lpars, Hpars, mandZandA, paxisandpP, dD, DEpars, hatpars, energies):
    
    e=0
    events=np.array([e,E])
    while events[0]<E:
        events[0]=events[0]+1
        alphapars, Lpars, Hpars, mandZandA, paxisandpP, dD, DEpars, hatpars, energies=initparameters(alphapars, Lpars, Hpars, mandZandA, paxisandpP, dD, DEpars, hatpars, energies)
        
        distanceL=getdistanceL(alphapars, Lpars)
        distanceH=getdistanceH(alphapars, Hpars)
        momentum=alphapars[9]
        
        filename="test.txt"
        writer=file('test.txt', 'a')
        np.savetxt(writer, events, fmt='%d',header='event', newline=' ')
        writer.close()
        writeresults(filename, 1, "test datas", mandZandA, alphapars, Lpars, Hpars)
        
        while True:
            #numerical solver
            
                momentum, distanceL, distanceH=solveequationsofmotion(DEpars,hatpars,mandZandA, momentum, distanceL, distanceH)
    
                if momentum<=0 or distanceL<=0.8 or distanceH<=0.8:
                #write end results
                        print "Results"
                        alphapars[9]=momentum
                    
                        writeresults(filename, 1, "test datae", mandZandA, alphapars, Lpars, Hpars)
                    
                        break