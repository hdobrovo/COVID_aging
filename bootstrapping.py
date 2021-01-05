import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy import integrate
from scipy.optimize import minimize

#time[0], old_a[1], young_a[2], old_b[3], young_b[4]
global par
data_file = "rockxdata.dat" 
tfile = np.loadtxt(data_file)
t = tfile[:,0]
ageda = tfile[:,1]
younga = tfile[:,2]
agedb = tfile[:,3]
youngb = tfile[:,4]

#=====================================================
# ODE function describing the virus and cell equations
def virus_dt(Z,t):
    global par
    b = par[0]
    p = par[1]
    c = 10 #fixed
    k = 3 #fixed
    d = par[2]

    [T, I1, I2, V] = [0, 1, 2, 3]
    dT = -b*Z[T]*Z[V]  #target
    dI1 = b*Z[T]*Z[V]- k*Z[I1] #infected
    dI2 = k*Z[I1] - d*Z[I2] #productive
    dV = p*Z[I2] - c*Z[V] #virus

    return [dT, dI1, dI2, dV]

#=====================================================
#inside of ssr
#3.Score Fit of System
#=========================================================
def SSR(pp):
    global par
    global sync
    global dead

    fpar=np.power(10,pp)
    par = []
    par.extend(fpar[:])
    upd = []
    sample_time=np.arange(25.)
 #  [dT, dI1, dI2, dV]
    y = [1.0, 0.0, 0.0, fpar[-1]]
    soln = odeint(virus_dt, y, sample_time)
    Vmodel = soln[:,-1]
    
    #Vmodel = np.where(Vmodel < 0.1, 0.1, Vmodel)
    findindex=lambda x:np.where(sample_time==x)[0][0]
    mindex=list(map(findindex,t))
    Vm=Vmodel[mindex]

    def ss(dat1, mod1):
        summ=0
        for i in range(len(dat1)): 
           if dat1[i]>-0.96:
               summ=summ+(dat1[i]-mod1[i])**2
        return summ	
    return ss(ndata[0], np.log10(Vm))
#========================================================
 #2.best fittings
# parameter models: b,p,d,v0,
aged_a=[4.55172351e-03, 4.04141857e+04, 8.17014160e-01, 2.31425556e-01] 
#ssra=1.9058684044693173
young_a=[1.89419364e-01, 3.59753428e+03, 3.46312442e+00, 4.78525884e-02] 
#ssray=3.2352004014225497
aged_b=[3.16227766e-01, 1.50096339e+03, 8.32120867e-01, 3.16227764e-01] 
#ssryob=2.5716744761052923
young_b=[6.30957344e-01, 6.72875572e+02, 1.29900649e+00, 1.00000000e+00] 
#2.6298538752335903

#### get best fitting data
bfpar=aged_b                                                       #####change parameter value
fpar=bfpar
par = []
par.extend(fpar[:])
upd = []
sample_time=np.arange(25.)
 #  [dT, dI1, dI2, dV]
y = [1.0, 0.0, 0.0, fpar[-1]]
soln = odeint(virus_dt, y, sample_time)
Vmodel = soln[:,-1]
Vmodel = np.where(Vmodel < 0.1, 0.1, Vmodel)
findindex=lambda x:np.where(sample_time==x)[0][0]
mindex=list(map(findindex,t))
bfdata=Vmodel[mindex]
bfdata = np.log10(bfdata)


fd=open('boots_aged_b2.dat','w')                                     ##change file
# Bootstrapping
res = bfdata - agedb                                                #####change data(aged,young\a,b)
nres = np.shape(res)[0]
 
j=0
while j<1000:
    fpar = np.log10(bfpar)
    idx = np.floor(nres*np.random.rand(1,nres))
    ndata = bfdata + res[idx.astype(int)]
    ndata = np.where(ndata < -1, -1, ndata)
    plt.plot(agedb,'ro',label='V data')                             #####change data
    plt.plot(ndata[0],'r-',label='V fit')
    plt.show()
    bnds=((-5,1),(0,7),(-4,1),(-4,2))
    answ = minimize(SSR,fpar, method='L-BFGS-B', bounds=bnds)
    if int(answ.success)==1:
        j=j+1
        print('ans',j,answ.fun)
        fd.write(str([answ.fun, 10**(answ.x),"\n"]))
    #print(i)
    
fd.close()
