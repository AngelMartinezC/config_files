import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from numpy.polynomial import polynomial
plt.rcParams.update({'font.size': 12})

#r0, c, rho, p, G0, T = np.loadtxt("cptrho.l5bi.d.15c",unpack=True)
#rc,mc,pc,rhoc,Tc,cs,gc,Ac,Hc,Nc = np.loadtxt("CD_ATMOSPHERE",unpack=True)
r0,mc,p,rho,T,cs,gc,Ac,Hc,Nc = np.loadtxt("CD_ATMOSPHERE",unpack=True)

RSun = 696.34e3 # km
#r1 = RSun*(r0-1)
#r1c = RSun*(-rc)
r1 = RSun*(-r0)

# -- Polynomial fitting

#x0 =  0.496e3*np.ones(len(r1))     # Upper and lower limits in km
x0 =  0.000e3*np.ones(len(r1))     # Upper and lower limits in km
x1 = -1.500e3*np.ones(len(r1))

s0 = np.argwhere(np.diff(np.sign(x0-r1)))
s1 = np.argwhere(np.diff(np.sign(x1-r1)))
idx0, idx1 = s0[0][0], s1[0][0]+2  # To get more in depth (+2)

r2 = np.linspace(r1[idx0],r1[idx1],len(r1[idx0:idx1]))
xx = r1[idx0:idx1+1]   # Radius in stablished range [-1.5e3,0]  #km
yy = rho[idx0:idx1+1]  # Rho in range
yy = p[idx0:idx1+1]    # PRS in range
#yyr= rho[idx0:idx1+1]  # Rho in range
yyd = p[idx0-1:idx1]  # PRS diff in range
xxd = r1[idx0-1:idx1]  # PRS diff in range
print(r2[0],r2[-1])
exit()

poly = 25
colors = plt.cm.jet_r(np.linspace(0,1,len(range(0,poly))))
for i in range(1,poly+1):
  #i = poly
  var = np.polyfit(xx, yy, i)
  #derivative = np.polyder(RHO)
  #dpdr = np.poly1d(derivative)
  
  #RHOr= np.polyfit(xx, yyr, 17)
  interp_var = np.poly1d(var)
  #rhonewr = np.poly1d(RHOr)
  
  fig, axs = plt.subplots(1, 2, figsize=(13,7))
  
  axs[0].semilogy(xx,yy,color='k',alpha=0.2,label='Model S',lw=10)
  axs[0].semilogy(xx,interp_var(xx),'-',label='Polyfit n=%2d'%i,
      color=colors[i-1])
  axs[0].set_ylim(1e-7,1e-5)  # limits for rho
  axs[0].set_ylim(4e4,1e7)   # limits for p
  #axs[0].set_ylim(1e5,4e9)  # limits for dp/dr log
  #axs[0].set_ylim(0,2.2e9)  # limits for dp/dr lin
  axs[0].set_ylabel(r'Density [g cm$^{-3}$]')
  axs[0].set_ylabel(r'Pressure [dyne]',labelpad=-4)
  #axs[0].set_ylabel(r"$\nabla p$ [dyne cm$^{-1}$]")
  
  
  axs[1].plot(xx,100*(interp_var(xx)-yy)/yy,label='Polyfit n='+str(i),
      color=colors[i-1])
  #axs[1].plot(xx[0:-1],100*(DIFF2-DIFF0)/DIFF0,label='Polyfit n='+str(i),
  #    color=colors[i-1])
  axs[1].set_ylim(-10,10)
  axs[1].set_ylim(-1.0,1.0) # limits for p
  #axs[1].set_ylim(-7.5,7.5) # limits for dp/dr
  axs[1].set_ylabel(r'% Error Density',labelpad=-4)
  axs[1].set_ylabel(r'% Error Pressure',labelpad=-4)
  #axs[1].set_ylabel(r"% Error $\nabla p$",labelpad=-4)
  
  axs[0].set_xlabel('Depth [Mm]')
  axs[1].set_xlabel('Depth [Mm]')
  axs[1].axhline(y=0,color='k',alpha=0.8)
  axs[0].grid()
  axs[1].grid()
  
  axs[0].legend(loc='lower left')
  axs[1].legend(loc='lower center')
  plt.savefig('pic.%02d.png'%i)
  if i==poly:
    pass
  else:
    plt.close('all')

plt.show()

