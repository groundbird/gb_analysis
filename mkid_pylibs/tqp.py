import numpy as np

# Al typical values
Tc   = 1.24       #K
kb   = 8.617e-5   #eV/K
#kbt  = kb*T       #eV
#kbtc = kb*Tc      #eV
delta= 1.76*kb*Tc #eV
N0   = 1.74e10    #eV**-1 microm**-3
tau0 = 450e-9     #sec
ita  = 0.57

def t2tqp(T):#macro m
    kbt = kb*T
    kbtc = kb*Tc
    return tau0*((kbtc/(2*delta))**2.5)*np.sqrt(Tc/T)*np.exp(delta/kbt)/np.sqrt(np.pi)

def t2nqp(T):
    kbt = kb*T
    kbtc = kb*Tc
    return 2*N0*np.sqrt(2*np.pi*kbt*delta)*np.exp(-delta/kbt)

def tqp2p(tqp, T):# 1/sec
    kbt = kb*T
    kbtc = kb*Tc
    V = 86.3584#microm**-3
#    Nqp = 2*V*1.74e10*np.sqrt(2*np.pi*kbt*delta)*np.exp(-delta/kbt)
    Nqp = V*tau0*N0*(kbtc)**3/(tqp*2*(delta**2))
    retev = Nqp*1.76*kbt/(tqp*ita)#eV/s
    retw = ev2w(retev)#W
    return retw

def ev2w(ev):
    return 1.602e-19*ev

def Pdiff(tqp_all, T):
    tqp_thermal = t2tqp(T)*1e-6
    print(tqp_thermal)
    return tqp2p(tqp_all,T) - tqp2p(tqp_thermal,T)