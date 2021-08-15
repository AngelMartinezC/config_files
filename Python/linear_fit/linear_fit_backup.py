import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from numpy.polynomial import polynomial
from scipy.interpolate import interp1d
import matplotlib as mpl
#mpl.rcParams['text.usetex'] = True
#mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command
plt.style.use(['science','std-colors'])#,'nature'])#,'ieee'])
#plt.rcParams['font.size'] = 20
#plt.rcParams['legend.fontsize'] = 16
#plt.rcParams['xtick.direction'] = 'in'
#plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.size'] = 5.0
plt.rcParams['xtick.minor.size'] = 3.0
plt.rcParams['ytick.major.size'] = 5.0
plt.rcParams['ytick.minor.size'] = 3.0

r0, c, rho, p, G0, T = np.loadtxt("cptrho.l5bi.d.15c",unpack=True)


G_const = 6.67e-11  # SI
RSun = 695.99036e3 # [km]
MSun = 1.989e33   # [g]
r1 = RSun*(r0-1)
#x2 = x2[::-1]
#x0 =  0.496e3*np.ones(len(r1))     # Upper and lower limits in km
x0 =  0.000e3*np.ones(len(r1))     # Upper and lower limits in km
x1 = -1.500e3*np.ones(len(r1))

s0 = np.argwhere(np.diff(np.sign(x0-r1)))
s1 = np.argwhere(np.diff(np.sign(x1-r1)))
idx0, idx1 = s0[0][0], s1[0][0]+2  # To get more in depth (+2)

UNIT_LENGTH   = 1e5  # conversion from cm to km
UNIT_DENSITY  = 1e-6 # conversion to 10x^{-6} g/cm3
UNIT_VELOCITY = 1e5  # conversion to leave TIME in s
UNIT_PRESSURE = UNIT_DENSITY*UNIT_VELOCITY**2

r2 = np.flip(r0[idx0+1::])

r = np.linspace(0,1,len(r2))
frho = interp1d(np.flip(r0), np.flip(rho), kind='cubic')


Mr, r_new, Mass = np.zeros(len(r2)), np.zeros(len(r2)), 0
for i in range(0,len(r2)):
  if i < len(r2)-1:
    Mass += 4*np.pi*r2[i]**2*frho(r2[i])*(r2[i+1]-r2[i])*RSun**3*1e5**3 #cgs
    Mr[i], r_new[i] = Mass, r2[i]
  else:
    Mr[i], r_new[i] = Mass, r2[i]
r2g = r2
fMr = interp1d(r2, Mr/MSun, kind='cubic')
g_gravity = G_const*MSun*fMr(r2)/(r2*RSun*1e3)**2/1e3


# -- Polynomial fitting

r2 = np.linspace(r1[idx0],r1[idx1],len(r1[idx0:idx1]))
xx1  = r1[idx0:idx1+1]   # Radius in stablished range [-1.5e3,0]  #km
xx2 = r1[idx0-1:idx1]  # PRS diff in range
rho1 = rho[idx0:idx1+1]  # Rho in range
prs1 = p[idx0:idx1+1]    # PRS in range
prs2 = p[idx0-1:idx1]  # PRS diff in range

irho = 16
iprs = 24
var_rho = np.polyfit(xx1, rho1, irho)
var_prs = np.polyfit(xx1, prs1, iprs)
rho_int = np.poly1d(var_rho)
prs_int = np.poly1d(var_prs)

dprs = prs1 - prs2
dxx  = xx1 - xx2
der0 = dprs/dxx

# Gravity interplation in range [-1.5e3,0]
fgravity = interp1d(-1*RSun*(1-r2g),g_gravity)
xx3 = np.linspace(0,min(xx1),len(xx1))

SI1 = 1e3*1e3*1e1
der1 = np.gradient(prs1,xx1)
der2 = np.gradient(prs_int(xx1),xx1)




#-- PLOTS
r0c,mc,p,rho,T,cs,gc,Ac,Hc,Nc = np.loadtxt("CD_ATMOSPHERE",unpack=True)
r1c = RSun*(-r0c)
xx2, rrho2 = np.loadtxt("GRAVITY0",unpack=True)

xx2 *=-1
rrho2 *= -1
nr =     xx2 #np.array(sorted(xx2))
nrho = rrho2 #np.array(sorted(rrho2))

print(max(nrho),min(nrho))


#print(min(xx3),max(xx3))
#fgravityc = interp1d(np.linspace(min(r1c),max(r1c),len(gc)),gc/1e2)
fgravityc = interp1d(r1c,gc/1e2)

plt.figure(figsize=(5,4))
#plt.figure(figsize=(10,4))
plt.plot(xx1/1e3, abs(der1/rho1)/SI1, label='Model S (np.gradient)',lw=1.5,alpha=0.8)
plt.plot(xx1/1e3, abs(der2/rho1)/SI1, lw=1.5, c='orange',
    label=r'Polyfit  $\rho: \mathcal{O}(x^{16})$, $p: \mathcal{O}(x^{24})$')
plt.plot(xx3/1e3, fgravity(xx3), label=r'$g=\frac{GM(r)}{r^2}$',lw=2,c='limegreen')
plt.plot(xx3/1e3,fgravityc(xx3), '--', label=r'$g$ Charlie Model l4b.14a', c='r')
plt.plot(nr/1e3, (nrho+0.02)*UNIT_PRESSURE/(UNIT_LENGTH*UNIT_DENSITY)/1e2, 'o-.',
    ms=2,label='PLUTO simulation $+0.02$',color='m')
plt.ylabel(r'gravity $(-\nabla p/\rho)~[\text{m s}^{-2}]$',labelpad=10)
plt.xlabel(r'$\text{Depth} ~[\text{Mm}]$')
plt.xlim(-1.6,0.1)
plt.ylim(240,290)
plt.grid(ls='--')
plt.legend(frameon=True)#,fancybox=True,framealpha=0.7)
plt.savefig('gravity_2.png',dpi=300)
plt.savefig('gravity_2.pdf')
plt.savefig('gravity_2.eps')
plt.show()


exit()

poly = 1
colors = plt.cm.jet_r(np.linspace(0,1,len(range(0,poly))))
for i in range(1,poly+1):
  #i = poly
  var = np.polyfit(xx, yy, i)
  derivative = np.polyder(var)
  dpdr = np.poly1d(derivative)
  
  RHOr= np.polyfit(xx, yyr, 17)
  interp_var = np.poly1d(var)
  rhonewr = np.poly1d(RHOr)
  
  fig, axs = plt.subplots(1, 2, figsize=(13,7))
  
  DIFF0 = abs(np.diff(yy)/np.diff(xx))
  DIFF2 = np.abs(np.diff(rhonew2(xx))/np.diff(xx))

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


exit()


"""
  Plots for density, pressure and pressure gradient
"""
for i in range(1,poly+1):
  RHO = np.polyfit(xx, yy, i)
  rhonew2 = np.poly1d(RHO)

  #axs[0].semilogy(xx[0:-1],DIFF0,color='k',alpha=0.2,
  #    label='Model S',lw=10)

  #axs[0].semilogy(xx2[0:-1],DIFF1,'-', label='Polyfit n=%2d'%i, 
  #    color=colors[i-1])
  axs[0].set_ylim(1e-7,1e-3)  # limits for rho
  #axs[0].set_ylim(1e4,1e10)   # limits for p
  #axs[0].set_ylim(1e5,4e9)  # limits for dp/dr log
  #axs[0].set_ylim(0,2.2e9)  # limits for dp/dr lin
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




