"""
 Check the sound speed (adiabatic) with Gamma*P/rho
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as spline

r0, c, rho, p, G0, T = np.loadtxt("cptrho.l5bi.d.15c",unpack=True)
#rho2, G12, gamma1, c2 = np.loadtxt("variable",unpack=True)

RSun = 696.34e3  # km
r1 = r0*RSun
#x2 = x2[::-1]
unit_pressure = (4*np.pi*1e-6)*8e8


#plt.plot(r0,c)
#plt.plot(r0,np.sqrt(G0*p/rho))
plt.plot(r0,1-c/np.sqrt(G0*p/rho))
plt.axhline(y=0,c='k',alpha=0.5)
plt.axvline(x=1,c='k',alpha=0.5)
plt.show()


exit()


K = 5.5e13
c = 3e10


#plt.plot(r0,G0,label=r'$\Gamma_1$ (Model S)')
#plt.plot(r0,G12,label=r'$\Gamma_1$ (thermodynam)')
#plt.plot(r0,gamma1,label=r'$\gamma$ (polytropic index)')
#plt.plot(r0,p)
plt.plot(r0,rho*c**2-p)

#plt.plot(r0,rho,label=r'$\rho$ (Model S)')
#plt.plot(r0,rho2,label=r'$\rho$ (thermodynam)')
#plt.plot(r0,(p/K)**(1/gamma1))
#plt.plot(r0,rho-rho2,label=r'$\Delta \rho$')
#plt.axhline(y=4/3,c='k',alpha=0.5)
#plt.axhline(y=5/3,c='k',alpha=0.5)
plt.legend()
plt.show()



exit()

DP, DR = [], []
for i in range(len(r1)):
  if i==0:
    Dp = p[i+1] - p[i]
    Dr = r0[i+1] - r0[i]
  else:
    Dp = p[i] - p[i-1]
    Dr = r0[i] - r0[i-1]
  DP.append(Dp)
  DR.append(Dr*RSun)
DR = np.array(DR)
DP = np.array(DP)

# Interpolation
Pnew = p[::-1]
Rnew = r1[::-1]
rhonew = rho[::-1]
PRS = spline(Rnew,Pnew)
RHON = spline(Rnew,rhonew)
XX2 = -1*ri+RSun
XX = np.linspace(np.min(XX2),np.max(XX2),1000)
unit_density = 1e-6
unit_velocity = 1e5
unit_pressure = unit_density*unit_velocity**2


DPRESS = []
DRAD = []
RHO = []
for i in range(len(XX)):
  if i==0:
    Dp = PRS(XX[i+1]) - PRS(XX[i])
    Dr = XX[i+1] - XX[i]
  else:
    Dp = PRS(XX[i]) - PRS(XX[i-1])
    Dr = XX[i] - XX[i-1]
  DPRESS.append(Dp)
  DRAD.append(Dr)
  RHO.append(RHON(XX[i]))
DPRESS = np.array(DPRESS)
DRAD = np.array(DRAD)
RHO = np.array(RHO)


"""
plt.figure()
plt.plot(XX-RSun,DPRESS/DRAD,'o',color='orange',label="interpolation")
plt.plot(-1*ri,-dpc*unit_pressure/dxc,'.',color='mediumblue',label="PLUTO")
plt.plot(r1-RSun,DP/DR,'--',color='lime',label="Model S")
plt.axvline(x=0,color='k')
plt.xlim(XX[0]-5e2-RSun,XX[-1]+5e2-RSun), plt.ylim(-2.1e6,0.1e6)
plt.ylabel(r"$\nabla P$ [dyne cm$^{-3}$]"), plt.xlabel(r"Depth [km]")
plt.legend(), plt.grid()

plt.figure()
plt.semilogy(XX-RSun,PRS(XX),'o',color='orange',label="interpolation")
plt.semilogy(-1*ri,pres_new*unit_pressure,'.',color='mediumblue',label="PLUTO")
plt.semilogy(r1-RSun,p,'--',color='lime',label="Model S")
plt.axvline(x=0,color='k')
plt.xlim(XX[0]-5e2-RSun,XX[-1]+5e2-RSun), plt.ylim(1.0e4,1.9e10)
plt.ylabel(r"Pressure [dyne cm$^{-2}$]"), plt.xlabel(r"Depth [km]")
plt.legend(), plt.grid()
"""

plt.figure()
plt.semilogy(XX-RSun,RHON(XX),'o',color='orange',label="interpolation")
plt.semilogy(-1*ri,rhoc*unit_density,'.',color='mediumblue',label="PLUTO")
plt.semilogy(r1-RSun,rho,'--',color='lime',label="Model S")
plt.axvline(x=0,color='k')
plt.xlim(XX[0]-5e2-RSun,XX[-1]+5e2-RSun), plt.ylim(1.0e-7,2e-3)
plt.ylabel(r"Density [g cm$^{-3}$]"), plt.xlabel(r"Depth [km]")
plt.legend(), plt.grid()

plt.figure()
NRES0 = 1e-3 # Conversion m to km
NRES1 = 1e-2 # Conversion cm to m
NRES = NRES0*NRES1
plt.plot(XX-RSun,DPRESS/(DRAD*RHO)*NRES,'o',color='orange',
    label="interpolation")
plt.plot(-1*ri,-dpc*unit_pressure/(dxc*rhoc*unit_density)*NRES,'.', \
    color='mediumblue',label="PLUTO")
plt.plot(r1-RSun,DP/(DR*rho)*NRES,'--',color='lime',label="Model S")
plt.axvline(x=0,color='k')
plt.xlim(XX[0]-5e2-RSun,XX[-1]+5e2-RSun)
#plt.ylim(-2.1e6,0.1e6)
plt.ylabel(r"Gravity [m s$^{-1}$]"), plt.xlabel(r"Depth [km]")
plt.legend(), plt.grid()
plt.show()
exit()


DER = (DP/DR)
plt.plot(r1-RSun,DER/rho*1e18,'.',color='orange')
plt.plot(-1*ri,-dpdx*1e8,'.',color='mediumblue')
plt.axvline(x=1,color='k')
plt.xlim(-1.2e4,5e2)
plt.show()
exit()

#plt.plot(r1-RSun,rho,'o-',color='green')
##plt.plot(r2,rho2/1e6,color='mediumblue')
##plt.plot(r3,rho3/1e6,color='red')
plt.xlim(-1.2e4,5e2)
#plt.ylim(-1e-3,5e-3)
plt.ylim(-0.01e12,0.1e12)
plt.plot(r1-RSun,p,'o-',color='green')
plt.plot(x2-1e4,PRS*1e4,'o')
plt.show()
exit()

plt.plot(r1-RSun,DER/rho*1e-10,'.')
plt.plot(x-1e4,grav*1e3,'.',color='orange')
plt.plot(x2-1e4,grav2*1e3,'.',color='mediumblue')
plt.xlim(-1.2e4,5e2)
#plt.ylim(-3.5e-1,-2.5e-1)
plt.axvline(x=1,color='k')
plt.show()



exit()


G = G0[::-1]
def delta():
  lnT = np.log(T[::-1])
  lnP = np.log(p[::-1])
  r = r0[::-1]
  dlnT, dlnP = [], []
  for i in range(len(r)):
    if i==0:
      dlnT.append(lnT[i+1]-lnT[i])
      dlnP.append(lnP[i+1]-lnP[i])
    else:
      dlnT.append(lnT[i]-lnT[i-1])
      dlnP.append(lnP[i]-lnP[i-1])
  return np.array(r), np.array(dlnT), np.array(dlnP)


r, dx, dy = delta()

Delad = 1-1/G
Del   = dx/dy
Rrel  = RSun*(r-1)
plt.figure(figsize=(11,8))
plt.plot(r,Del,'-',color='mediumblue',label=r'$\nabla$')
plt.plot(r,Delad,'-',color='red',label=r'$\nabla_{ad}$')
#plt.plot(Rrel,Del-Delad)
plt.grid()
plt.axhline(y=0,c='k',lw=1)
#plt.xscale('log')
plt.legend()
plt.show()
