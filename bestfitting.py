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
olda = tfile[:,1]
younga = tfile[:,2]
oldb = tfile[:,3]
youngb = tfile[:,4]
data = olda                                        ##change

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

def virus(t):
    global par2
    vp = par2[0]
    tp = par2[1]
    lg = par2[2]
    ld = par2[3]
    v=2*vp/(np.exp(-lg*(t-tp))+np.exp(ld*(t-tp)))
    return v
#=====================================================
# parameter models: b,p,d,v0,
old_a=[4.55172351e-03, 4.04141857e+04, 8.17014160e-01, 2.31425556e-01] 
ssra=1.9058684044693173
young_a=[1.89419364e-01, 3.59753428e+03, 3.46312442e+00, 4.78525884e-02] 
ssray=3.2352004014225497

old_b=[3.16227766e-01, 1.50096339e+03, 8.32120867e-01, 3.16227764e-01] 
ssryob=2.5716744761052923
young_b=[6.30957344e-01, 6.72875572e+02, 1.29900649e+00, 1.00000000e+00] 
#2.6298538752335903

# parameters equations: vp, tp, lu, ld
taged_a=[2.51890561e+03, 3.48946018e+00, 2.80808611e+00, 9.29798458e-01] #1.2836154349305864
tyoung_a=[368.16723669,  1.03555333, 103.19043954,   2.54146864] #3.233827925227035 
taged_b=[116.60483706,   1.18495742,   6.52690098,   1.12504036] #2.4184731017296444 
tyoung_b=[20.73555724,  0.96125494,  6.32663129,  2.44462398] #0.46365883053489826  

#solve model
#=========================================================
fpar=old_a                                           ###############change
par = []
par.extend(fpar[:])
sample_time=np.arange(0,25.,0.25)

#odeint model
y = [1.0, 0.0, 0.0, fpar[-1]]
soln = odeint(virus_dt, y, sample_time)
Vmodel = soln[:,-1]
#Vmodel = np.where(Vmodel < 0.1, 0.1, Vmodel)
Vmodel = np.log10(Vmodel)

#triangle function
fpar2=tyoung_b#np.power(10,pp)
par2 = []
par2.extend(fpar2[:])

Vtriang = virus(sample_time)

#plotting
plt.plot(sample_time,Vmodel,'k-',label='Predicted Model')
#plt.plot(sample_time,np.log10(Vtriang),'b--',label='Predicted Equation')
plt.plot(t,data,'ro',label='Experimental')
plt.legend( prop={'size': 14})
plt.xlabel('Days Post Infection', fontsize=14)
plt.ylabel('SARS-CoV-2 Viral Titer $Log_{10}$', fontsize=14)
plt.axis([0, 25, -2, 4.2])
plt.title('Nasal Swab Aged Animal', fontsize=18)
#plt.title('Throat Swab Aged Animal', fontsize=18)               #####change
plt.gca().spines['right'].set_color('none')
plt.gca().spines['top'].set_color('none')
plt.savefig('nasalaged2.pdf')                                    ####change
plt.show()

