"""
Check magnetic field profile. What I do here is to find the proper 
conditions on internal density
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as spline

r0, c, rho, p, G0, T = np.loadtxt("cptrho.l5bi.d.15c",unpack=True)

RSun = 696.34  # Mm
r0 = RSun*(r0-1)
rho *= 1e6

K  = 1.2e15      # Polytropic constant
c  = 3e10        # Light speed
Bc = 1.0*8.5e-3  # inverse of Ideal EOS constant
C  = 1/Bc        # ideal gas equation constant

alpha_m = 1/1.5
beta  = 1

Bext = 0.1   # Gauss
Bext = np.linspace(0.1,10000,len(rho))
Bint = 1000


alpha  = 1.25
alpha2 = 18.4
B0 = 3000
h0 = 10

def L(z):
  return 1/(1+np.exp(-z))
def Bz(z):
  return B0*np.exp(-(z/alpha)*L(-z/alpha)*L(z/alpha2))
def hz(z):
  return h0*np.sqrt(B0/Bz(z))
def Brz(r,z):
  return Bz(z)*np.exp(-np.log(2)*r**2/hz(z)**2)

z = np.linspace(-2.5,0,1000)
r = np.linspace(-20,20,1000)
R, Z = np.meshgrid(r,z)
plt.imshow(np.log(Brz(R,Z)),aspect='auto',extent=[np.min(R),np.max(R),np.min(Z),np.max(Z)],cmap='viridis', origin='lower')
plt.colorbar(label='Magnetic Field [Gauss]')
plt.grid()
plt.show()
exit()


plt.plot(z,hz(z))
plt.show()


exit()
rr = np.linspace(0,30,1000)
#plt.semilogy(r0,Bz(r0))
#plt.semilogy(rr,Brz(rr,0))
plt.plot(rr,(180/np.pi)*np.arccos(Brz(rr,0)/B0))
plt.xlim(0,30)
#plt.ylim(1e3,1e6)
plt.grid()
plt.show()
exit()







# ------- TESTS ----------

# --- B = B(y) ---

for beta in np.arange(0.1,1/alpha_m+0.1,0.1):
  Bint2 = np.sqrt(8*np.pi*C*(1-alpha_m*beta)*rho*T+Bext**2)
  plt.semilogy(r0,Bint2,label=r'$B_y$ $\beta=$'+str(round(beta,2)))
  
plt.semilogy(r0,Bext,'--',color='mediumblue',label=r'$B_{ext}$')
plt.axvline(x=0*RSun,c='k',alpha_m=0.5)
plt.axhline(y=0,c='k',alpha_m=0.5)
plt.xlim(-11,1.001)
plt.ylim(-1e4,5e5)
plt.legend()
plt.grid()
plt.show()
exit()





# --- B constant ---

rhoin = rho/alpha_m - (Bint**2 - Bext**2)/(8*np.pi*alpha_m*T)
prsin = p - (Bint**2 - Bext**2)/(8*np.pi)

plt.plot(r0,rho,'-',label=r'$\rho$ (Model S)',color='k')
plt.semilogy(r0,rhoin,'mediumblue',label=r'$\rho_{int}$')
#plt.plot(r0,p,'-',label=r'$\rho$ (Model S)',color='k')
#plt.plot(r0,prsin,'mediumblue',label=r'$\rho_{int}$')
plt.axvline(x=0*RSun,c='k',alpha_m=0.5)
plt.axhline(y=0,c='k',alpha_m=0.5)
plt.xlim(-11,1.001)
plt.ylim(-1e-4,600)
#plt.ylim(-1e-4,8e9)
plt.grid()
plt.legend()
plt.show()

exit()
rhonew = p/(T*C)
Bin = np.sqrt(8*np.pi*C*(1-alpha_m*beta)*rho*T)


#for i in [1]:#np.linspace(0.5,1.5,10):
#  plt.plot(r0,(p/T*(B*i)),label='Ideal EOS '+str(round(i,2)))

#plt.plot(r0,rho,label=r'$\rho$ (Model S)',color='k')
#plt.plot(r0,rhonew,label='Ideal EOS')
#plt.plot(r0,(1e6*(p/K)**(3/5)),'-',label='Polytropic (constant)',color='r')
plt.xlim(-11,1.001)
plt.ylim(-1e-4,600)
#plt.ylim(-0,2)
plt.axhline(y=1,c='k',alpha_m=0.5)
plt.axvline(x=0*RSun,c='k',alpha_m=0.5)
plt.legend()

plt.show()

