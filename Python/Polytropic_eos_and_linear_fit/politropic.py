import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from numpy.polynomial import polynomial
plt.rcParams.update({'font.size': 12})

r0, c, rho, p, G0, T = np.loadtxt("cptrho.l5bi.d.15c",unpack=True)
rho2, G12, gamma1, c2 = np.loadtxt("variable",unpack=True)

RSun = 696.34 # km
r1 = RSun*(r0-1)
#x2 = x2[::-1]
unit_pressure = (4*np.pi*1e-6)*8e8

K = 1.2e15
K2= 1.2e15
c = 3e10
A = 2.3e4


# -- Polynomial fitting
x0 = 0*np.ones(len(r1))
x1 = -10*np.ones(len(r1))

s0 = np.argwhere(np.diff(np.sign(x0-r1)))
s1 = np.argwhere(np.diff(np.sign(x1-r1)))
idx0, idx1 = s0[0][0]+1, s1[0][0]+2

r2 = np.linspace(r1[idx0],r1[idx1],len(r1[idx0:idx1]))

xx = r1[idx0:idx1]   # Radius in range [-10,0]
yy = rho[idx0:idx1]  # Rho in range
yy = p[idx0:idx1]  # PRS in range
yyr= rho[idx0:idx1]  # Rho in range
yyd = p[idx0-1:idx1-1]  # PRS diff in range
xxd = r1[idx0-1:idx1-1]  # PRS diff in range
#xx = r1[0:idx1]     # Radius in range [-10,:]
#yy = rho[0:idx1]    # Rho in range of r in [-10,:]
xx2 = np.linspace(np.min(xx),np.max(xx),300)
#xx2 = np.linspace(-10.05,0.5,300)

print(r0[0],r0[1])
print(xx2[0],xx2[-1])
poly = 25
colors = plt.cm.jet_r(np.linspace(0,1,len(range(0,poly))))

# Example
RHO = np.polyfit(xx, yy, 17)
derivative = np.polyder(RHO)
dpdr = np.poly1d(derivative)

RHOr= np.polyfit(xx, yyr, 17)

rhonew2 = np.poly1d(RHO)
rhonewr = np.poly1d(RHOr)

plt.plot(xx,(dpdr(xx)/rhonewr(xx))*-1)
plt.semilogy(xx,(8.11224e9*xx-2.7379e12)*-1,label='linear')
plt.xlabel('Depth [Mm]')
plt.ylabel(r'|$\nabla p/\rho$| [AU]')
plt.grid()
plt.ylim(2.6e12,2.87e12)
#plt.plot(xx,rhonew2(xx))
plt.show()
exit()

plt.plot(xx,yy)
plt.plot(xx2,rhonew2(xx2))
plt.grid()
plt.show()


exit()


"""
  Plots for density, pressure and pressure gradient
"""
for i in range(1,poly+1):
  RHO = np.polyfit(xx, yy, i)
  rhonew2 = np.poly1d(RHO)
  fig, axs = plt.subplots(1, 2, figsize=(13,7))
  axs[0].set_xlabel('Depth [Mm]')
  axs[1].set_xlabel('Depth [Mm]')
  axs[1].axhline(y=0,color='k',alpha=0.8)

  
  DIFF0 = abs(np.diff(yy)/np.diff(xx))
  DIFF1 = np.abs(np.diff(rhonew2(xx2))/np.diff(xx2))
  DIFF2 = np.abs(np.diff(rhonew2(xx))/np.diff(xx))


  axs[0].semilogy(xx,yy,color='k',alpha=0.2,label='Model S',lw=10)
  #axs[0].semilogy(xx[0:-1],DIFF0,color='k',alpha=0.2,
  #    label='Model S',lw=10)
  axs[0].grid()
  axs[1].grid()

  #axs[0].semilogy(xx2[0:-1],DIFF1,'-', label='Polyfit n=%2d'%i, 
  #    color=colors[i-1])
  axs[0].semilogy(xx2,(rhonew2(xx2)),'-',label='Polyfit n=%2d'%i,
      color=colors[i-1])
  axs[0].set_ylim(1e-7,1e-3)  # limits for rho
  #axs[0].set_ylim(1e4,1e10)   # limits for p
  #axs[0].set_ylim(1e5,4e9)  # limits for dp/dr log
  #axs[0].set_ylim(0,2.2e9)  # limits for dp/dr lin
  axs[0].set_ylabel(r'Density [g cm$^{-3}$]')
  #axs[0].set_ylabel(r'Pressure [dyne]',labelpad=-4)
  #axs[0].set_ylabel(r"$\nabla p$ [dyne cm$^{-1}$]")
  
  axs[1].plot(xx,100*(rhonew2(xx)-yy)/yy,label='Polyfit n='+str(i),
      color=colors[i-1])
  #axs[1].plot(xx[0:-1],100*(DIFF2-DIFF0)/DIFF0,label='Polyfit n='+str(i),
  #    color=colors[i-1])
  axs[1].set_ylim(-10,10)
  #axs[1].set_ylim(-2.5,2.5) # limits for p
  #axs[1].set_ylim(-7.5,7.5) # limits for dp/dr
  axs[1].set_ylabel(r'% Error Density',labelpad=-4)
  #axs[1].set_ylabel(r'% Error Pressure',labelpad=-4)
  #axs[1].set_ylabel(r"% Error $\nabla p$",labelpad=-4)

  axs[0].legend(loc='lower left')
  axs[1].legend(loc='lower center')
  plt.savefig('pic.%02d.png'%i)
  plt.close('all')

plt.show()
exit()










#plt.plot(r0,G0,label=r'$\Gamma_1$ (Model S)')
#plt.plot(r0,G12,label=r'$\Gamma_1$ (thermodynam)')
#plt.plot(r0,gamma1,label=r'$\gamma$ (polytropic index)')
#plt.plot(r0,p)
#plt.plot(r0,rho*c**2-p,label='Difference')
#plt.plot(r0,rho2,label=r'$\rho$ (thermodynam)')
#plt.plot(r1,rho*T/A,label='Ideal EOS')
#plt.xlim(0.99,1.001)
#plt.ylim(-1e-7,1.0e-3)


#plt.plot(r1,rho,'-',label=r'$\rho$ (Model S)',color='mediumblue')
#plt.plot(r1,(p/K2)**(1/gamma1),label='Polytropic')
#plt.plot(r1,(p/K)**(3/5),label='Polytropic (constant)')
#plt.plot(r1[idx0:idx1],rhonew,label='Polyfit n='+str(i))
#plt.ylim(-1e-4,8e-4)
#plt.xlim(-11,1.001)



#plt.plot(r0,rho-rho2,label=r'$\Delta \rho$')
#plt.axhline(y=4/3,c='k',alpha=0.5)
#plt.axhline(y=5/3,c='k',alpha=0.5)
plt.axhline(y=0,c='k',alpha=0.5)
#plt.axvline(x=0,c='k',alpha=0.5)
plt.axvline(x=r1[idx1],c='k',alpha=0.5)
plt.axvline(x=r1[idx0],c='k',alpha=0.5)
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
derivativeplt.plot(Rrel,Del-Delad)
plt.grid()
plt.axhline(y=0,c='k',lw=1)
#plt.xscale('log')
plt.legend()
plt.show()
