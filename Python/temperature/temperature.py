"""
Check magnetic field profile. What I do here is to find the proper 
conditions on internal density
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as spline
plt.rcParams.update({'font.size': 12})

r0, c, rho, p, G0, T = np.loadtxt("cptrho.l5bi.d.15c",unpack=True)

RSun = 696.34  # Mm
r0 = RSun*(r0-1)

K  = 1.2e15      # Polytropic constant
c  = 3e10        # Light speed
Bc = 1.0*8.5e-9  # inverse of Ideal EOS constant
C  = 1/Bc        # ideal gas equation constant

alpha = 1/1.5
beta  = 1

Bext = 0.1   # Gauss
Bint = 1000

# --- B constant ---

rhoin = rho/alpha - (Bint**2 - Bext**2)/(8*np.pi*alpha*T)

plt.plot(r0,T,label='Temperature (Model S)')
plt.plot(r0,p/(rho*C),label='Temperature (Ideal gas)')
#plt.plot(r0,p,label='P (Model S)')
#plt.plot(r0,rho,'k',label='P (Model S)')
plt.axvline(x=0*RSun,c='k',alpha=0.5)
plt.axvline(x=-10,c='k',alpha=0.5)
plt.xlim(-11,1.001)
plt.ylim(-1e-4,1e5)
#plt.ylim(-2e8,6.5e9)
#plt.ylim(-1e-4,0.0008)
plt.grid(ls='--')
plt.legend()
plt.show()

exit()
rhonew = p/(T*C)
Bin = np.sqrt(8*np.pi*C*(1-alpha*beta)*rho*T)


#for i in [1]:#np.linspace(0.5,1.5,10):
#  plt.plot(r0,(p/T*(B*i)),label='Ideal EOS '+str(round(i,2)))

#plt.plot(r0,rho,label=r'$\rho$ (Model S)',color='k')
#plt.plot(r0,rhonew,label='Ideal EOS')
#plt.plot(r0,(1e6*(p/K)**(3/5)),'-',label='Polytropic (constant)',color='r')
plt.xlim(-11,1.001)
plt.ylim(-1e-4,600)
#plt.ylim(-0,2)
plt.axhline(y=1,c='k',alpha=0.5)
plt.axvline(x=0*RSun,c='k',alpha=0.5)
plt.legend()

plt.show()

