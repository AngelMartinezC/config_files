import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from numpy.polynomial import polynomial
from scipy.interpolate import interp1d
import matplotlib as mpl
import warnings
warnings.filterwarnings("ignore")
#mpl.rcParams['text.usetex'] = True
#mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command

#plt.style.use(['science','std-colors'])#,'nature'])#,'ieee'])
#plt.rcParams['font.size'] = 20
#plt.rcParams['legend.fontsize'] = 16
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.size'] = 5.0
plt.rcParams['xtick.minor.size'] = 3.0
plt.rcParams['ytick.major.size'] = 5.0
plt.rcParams['ytick.minor.size'] = 3.0


# -- Read functions containing variables
r0, c, rho, p, G0, T = np.loadtxt("cptrho.l5bi.d.15c",unpack=True)
Dm,h,T,Vt,ne,nH,tau = np.loadtxt("Ding_Fang.txt",unpack=True)
r0c,mc,pc,rhoc,Tc,cs,gc,Ac,Hc,Nc = np.loadtxt("CD_ATMOSPHERE",unpack=True)

u = 1.66054e-24  # uma
mH = 1.00784*u
me = 9.10938e-28
aHe = 0.1   # Helium abundance 
kB = 1.3807e-16  # Boltzmann constant

# -- Constants for Model S
G_const = 6.67e-11  # SI
RSun = 695.99036e3 # [km]
MSun = 1.989e33   # [g]
r1  = RSun*(r0-1)
r1c = RSun*(-r0c)

# -- Limits of interpolation
x0 =  0.490e3*np.ones(len(r1))     # Upper and lower limits in km
x1 = -1.500e3*np.ones(len(r1))
s0 = np.argwhere(np.diff(np.sign(x0-r1)))
s1 = np.argwhere(np.diff(np.sign(x1-r1)))
idx0, idx1 = s0[0][0], s1[0][0]  # To get more in depth (+2)
#UNIT_LENGTH   = 1e5  # conversion from cm to km
#UNIT_DENSITY  = 1e-6 # conversion to 10x^{-6} g/cm3
#UNIT_VELOCITY = 1e5  # conversion to leave TIME in s
#UNIT_PRESSURE = UNIT_DENSITY*UNIT_VELOCITY**2
r_range = r1[idx0:idx1]
rho_range = rho[idx0:idx1]
prs_range = p[idx0:idx1]
T_range = T[idx0:idx1]

# Calculate Pressure and Density
rhoDF = me*ne + mH*nH
prsDF = ((1+aHe)*nH+ne)*kB*T + 0.5*rhoDF*Vt**2

# -- Plots Thermodynamical variables from Model S and Maltby
def plot_p():
  with plt.style.context(['science']):
    #plt.figure(figsize=(6.04,4.3))
    plt.semilogy(r_range/1e3,prs_range,'k',label='Model S')
    #plt.semilogy((rM-550-70)/1e3,pM*0.66,'g',ls='--')
    plt.axvline(x=-70/1e3,ls='--',color='k')
    plt.axvline(x=(-70-550)/1e3,ls='--',color='b')
    plt.semilogy((h-300)/1e3,prsDF*1,'r',ls='--')
    plt.ylim(1e2,1e7)
    plt.xlim(-1.5,0.5)
    plt.xlabel('z [Mm]')
    plt.legend(frameon=True)
    plt.ylabel(r'Density $[\text{g cm}^{-3}]$')
    plt.grid()
    plt.savefig("IMAGES/penumbral_density.png",dpi=300)
    plt.show()
#plot_p()

# -- Weight Function --
z0 = -70
zp = -240# -300 #220
def wp(z,zp0=220):
  I2 = lambda z: np.cos(np.pi*(zp-80-z)/(120))/2 + 1/2
  I4 = lambda z: np.cos(np.pi*(zp+20-z)/(1100-zp))/2 + 1/2
  wu = np.piecewise(z,[z>-1e6,z>zp-zp0,z>zp-80,z>zp+20,z>800],[0,I2,1,I4,0])
  return wu
def plot_weight():
  with plt.style.context(['science']):
    z = np.linspace(-1500,1000, 100)
    #plt.style.use(['science','std-colors'])#,'nature'])#,'ieee'])
    plt.axvline(x=(zp-280)/1e3,ls='--',color='b'     ,label='$z_p-280$ km')
    plt.axvline(x=(zp- 80)/1e3,ls='--',color='orange',label='$z_p-80$ km')
    plt.axvline(x=(zp+ 20)/1e3,ls='--',color='g'     ,label='$z_p+20$ km')
    plt.axvline(x=(   800)/1e3,ls='--',color='m'     ,label='$800$ km')
    plt.plot(z/1e3,wp(z),'-',color='k')
    plt.axvline(x=(z0)/1e3,ls='--',color='k')
    plt.axvline(x=(zp)/1e3,ls='--',color='r')
    plt.xlabel("Depth [Mm]")
    plt.ylabel("$w_p(z)$")
    plt.legend(frameon=True,fontsize=8)
    plt.savefig("IMAGES/Weight_Function_p.png",dpi=300)
    #plt.savefig("Weight_Function.pdf")
    #plt.savefig("Weight_Function.eps")
    plt.show()
#plot_weight()
#exit()


z = np.linspace(-1500,500,1000)
x2int   = h+z0+zp
rho2int, prs2int       = rhoDF, prsDF
rho_int_DF, prs_int_DF = spline(-x2int,rho2int,k=2), spline(-x2int,prs2int,k=2)
P_int_S   = spline(-r1,p)
rho_int_S = spline(-r1,rho)

# -- Termodynamical variables as a function of z
def Pqs(z):
  return P_int_S(-z)
def rhoqs(z):
  return rho_int_S(-z)
def PDF(z):
  return prs_int_DF(-z)
def rhoDF(z):
  return rho_int_DF(-z)

# --------------------------------------------------------------------------
# -- Umbral thermodynamical variables as function of depth and plots
def Pp(z):
  I1 = lambda z: Pqs(z)
  I2 = lambda z: np.exp( (1-wp(z))*np.log(Pqs(z)) + wp(z)*np.log(PDF(z)) )
  ss = np.piecewise(z,[z<-1000,z>=-1000],[I1,I2])
  return ss
def rhop(z,zz=-675):
  I1 = lambda z: rhoqs(z)
  #I2 = lambda z: np.exp( (1-wp(z))*np.log(rhoqs(z)) + wp(z)*np.log(0.66*rhoDF(z)) )
  I2 = lambda z: np.exp( (1-wp(z,110))*np.log(rhoqs(z)) + wp(z,110)*np.log(rhoDF(z)) )
  ss = np.piecewise(z,[z<zz,z>=zz],[I1,I2])
  return ss

def plot_density_u():
  with plt.style.context(['science']):
    plt.axvline(x=-70/1e3,ls='--',color='k')
    plt.axvline(x=(zp)/1e3,ls='--',color='r')
    plt.semilogy(z/1e3,rhoqs(z),'k',label='Model S')
    plt.semilogy((h+z0+zp)/1e3,rho2int,'o-',color='orange',label=r'Ding \& Fang')
    plt.semilogy(z/1e3,rhop(z),'b',ls='-.',label=r'$\rho_p(z)$')
    plt.ylim(1e-9,1e-5),plt.xlim(-1.5,0.5),plt.grid()
    plt.ylabel("Density [g cm$^{-3}$]"),plt.xlabel("Depth [Mm]")
    plt.legend(frameon=True,fontsize=8)
    plt.savefig("Density.png",dpi=300)
    plt.show()
#plot_density_u()
def plot_pressure_u():
  with plt.style.context(['science']):
    #plt.figure(figsize=(6,4.3))
    plt.semilogy(z/1e3,Pqs(z),'k',label='Model S')
    plt.semilogy(x2int/1e3,prs2int,'o-',color='orange',label=r'Ding \& Fang')
    #plt.semilogy(z/1e3,PMaltby(z),'b')
    plt.semilogy(z/1e3,Pp(z),'b',ls='-.',label='$P_p(z)$')
    plt.axvline(x=-70/1e3,ls='--',color='k')
    plt.axvline(x=(zp)/1e3,ls='--',color='r')
    #plt.ylim(1e-9,1e-5)
    plt.ylim(1e2,1e7)
    plt.xlim(-1.5,0.5)
    plt.grid()
    plt.xlabel("Depth [Mm]")
    plt.ylabel("Pressure [dyn]")
    plt.legend(frameon=True,fontsize=8)
    plt.savefig("Pressure.png",dpi=300)
    plt.show()
#plot_pressure_u()
# --------------------------------------------------------------------------



# -- Interpolation for umbral pressure ---------------------------------
def interp_prs_penumbra():
  zi = np.linspace(-1500,-650,1000)
  zf = np.linspace(-650,-300,1000)
  zg = np.linspace(-300,0,1000)
  zh = np.linspace(0,500,1000)
  P2intpluto = np.polyfit(zi,Pp(zi),25)
  P2intpluto2= np.polyfit(zf,Pp(zf),22)
  P2intpluto3= np.polyfit(zg,Pp(zg),10)
  P2intpluto4= np.polyfit(zh,Pp(zh),15)
  P_int_pluto = np.poly1d(P2intpluto)
  P_int_pluto2= np.poly1d(P2intpluto2)
  P_int_pluto3= np.poly1d(P2intpluto3)
  P_int_pluto4= np.poly1d(P2intpluto4)

  var = P2intpluto4
  print(" ")
  print(" ")
  for i in range(len(var)):
    if (i==len(var)-1):
      print("%.9e*pow(z,%d);"%(np.flip(var)[i],i),end=" ")
      break
    elif (i==0):
      print("%.9e +"%(np.flip(var)[i]),end=" ")
    elif (i==1):
      print("%.9e*z +"%(np.flip(var)[i]))
    elif (i%2 == 1):
      print("%.9e*pow(z,%d) +"%(np.flip(var)[i],i))
    else:
      print("%.9e*pow(z,%d) +"%(np.flip(var)[i],i),end=" ")
  print(" ")
  print(" ")
  
  plt.subplot(121)
  #plt.semilogy(zi/1e3,rhoqs(zi))
  plt.semilogy(zi/1e3,Pp(zi),'k')
  plt.semilogy(zf/1e3,Pp(zf),'k')
  plt.semilogy(zg/1e3,Pp(zg),'k')
  plt.semilogy(zh/1e3,Pp(zh),'k')
  plt.semilogy(zi/1e3,P_int_pluto(zi))
  plt.semilogy(zf/1e3,P_int_pluto2(zf))
  plt.semilogy(zg/1e3,P_int_pluto3(zg))
  plt.semilogy(zh/1e3,P_int_pluto4(zh))
  plt.ylim(1e2,1e7)
  plt.xlim(-1.6,0.5),plt.grid()
  
  plt.subplot(122)
  plt.plot(zi/1e3,100*(Pp(zi)-P_int_pluto(zi)) /Pp(zi))
  plt.plot(zf/1e3,100*(Pp(zf)-P_int_pluto2(zf))/Pp(zf))
  plt.plot(zg/1e3,100*(Pp(zg)-P_int_pluto3(zg))/Pp(zg))
  plt.plot(zh/1e3,100*(Pp(zh)-P_int_pluto4(zh))/Pp(zh))
  plt.ylim(-1.5,1.5), plt.xlim(-1.6,0.5), plt.grid()
  plt.show()
interp_prs_penumbra()
exit()

# -- Interpolation for umbral density ---------------------------------
def interp_rho_penumbra():
  zi = np.linspace(-1500,-650,1000)
  zf = np.linspace(-650,-300,1000)
  zg = np.linspace(-300,0,1000)
  zh = np.linspace(0,500,1000)
  rho2intpluto = np.polyfit(zi,rhop(zi),19)
  rho2intpluto2= np.polyfit(zf,rhop(zf),19)
  rho2intpluto3= np.polyfit(zg,rhop(zg),10)
  rho2intpluto4= np.polyfit(zh,rhop(zh),10)
  rho_int_pluto = np.poly1d(rho2intpluto)
  rho_int_pluto2= np.poly1d(rho2intpluto2)
  rho_int_pluto3= np.poly1d(rho2intpluto3)
  rho_int_pluto4= np.poly1d(rho2intpluto4)
  
  var = rho2intpluto4
  print(" ")
  print(" ")
  for i in range(len(var)):
    if (i==len(var)-1):
      print("%.9e*pow(z,%d);"%(np.flip(var)[i],i),end=" ")
      break
    elif (i==0):
      print("%.9e +"%(np.flip(var)[i]),end=" ")
    elif (i==1):
      print("%.9e*z +"%(np.flip(var)[i]))
    elif (i%2 == 1):
      print("%.9e*pow(z,%d) +"%(np.flip(var)[i],i))
    else:
      print("%.9e*pow(z,%d) +"%(np.flip(var)[i],i),end=" ")
  print(" ")
  print(" ")
  #exit()
  
  plt.subplot(121)
  #plt.semilogy(zi/1e3,rhoqs(zi))
  plt.semilogy(zi/1e3,rhop(zi),'k')
  plt.semilogy(zf/1e3,rhop(zf),'k')
  plt.semilogy(zg/1e3,rhop(zg),'k')
  plt.semilogy(zh/1e3,rhop(zh),'k')
  plt.semilogy(zi/1e3,rho_int_pluto(zi))
  plt.semilogy(zf/1e3,rho_int_pluto2(zf))
  plt.semilogy(zg/1e3,rho_int_pluto3(zg))
  plt.semilogy(zh/1e3,rho_int_pluto4(zh))
  plt.ylim(1e-9,1e-5)
  plt.xlim(-1.6,0.5),plt.grid()
  
  plt.subplot(122)
  plt.plot(zi/1e3,100*(rhop(zi)- rho_int_pluto(zi))/rhop(zi))
  plt.plot(zf/1e3,100*(rhop(zf)-rho_int_pluto2(zf))/rhop(zf))
  plt.plot(zg/1e3,100*(rhop(zg)-rho_int_pluto3(zg))/rhop(zg))
  plt.plot(zh/1e3,100*(rhop(zh)-rho_int_pluto4(zh))/rhop(zh))
  plt.ylim(-1.5,1.5), plt.xlim(-1.6,0.5), plt.grid()
  plt.show()
interp_rho_penumbra() 

























exit()
# Plot Density
plt.figure(figsize=(5.2,3.7))
plt.semilogy((h-300)/1e3,rhoDF)
plt.xlim(-1.5,0.5)
plt.ylim(1e-9,1e-5)
plt.show()
exit()

# Plot Temperature
plt.figure(figsize=(5.2,3.7))
plt.plot((h-300)/1e3,T/1e3)
plt.xlim(-1.5,0.5)
plt.ylim(0,20)
plt.show()



exit()

# -- Constants for Model S
G_const = 6.67e-11  # SI
RSun = 695.99036e3 # [km]
MSun = 1.989e33   # [g]
r1  = RSun*(r0-1)
r1c = RSun*(-r0c)


# -- Limits of interpolation
x0 =  0.490e3*np.ones(len(r1))     # Upper and lower limits in km
x1 = -1.500e3*np.ones(len(r1))
s0 = np.argwhere(np.diff(np.sign(x0-r1)))
s1 = np.argwhere(np.diff(np.sign(x1-r1)))
idx0, idx1 = s0[0][0], s1[0][0]  # To get more in depth (+2)
#UNIT_LENGTH   = 1e5  # conversion from cm to km
#UNIT_DENSITY  = 1e-6 # conversion to 10x^{-6} g/cm3
#UNIT_VELOCITY = 1e5  # conversion to leave TIME in s
#UNIT_PRESSURE = UNIT_DENSITY*UNIT_VELOCITY**2
r_range = r1[idx0:idx1]
rho_range = rho[idx0:idx1]
prs_range = p[idx0:idx1]
T_range = T[idx0:idx1]


# -- Plots Thermodynamical variables from Model S and Maltby
def plot_p():
  plt.figure(figsize=(6.04,4.3))
  plt.semilogy(r_range/1e3,prs_range,'k',label='Model S')
  #plt.semilogy((rM-550-70)/1e3,pM*0.66,'g',ls='--')
  plt.semilogy((rM-550-70)/1e3,pM*1,'r',ls='--')
  plt.semilogy((rM-550-70)/1e3,pMM*1,'g',ls='--')
  plt.semilogy((rM-550-70)/1e3,pML*1,'m',ls='--')
  plt.axvline(x=-70/1e3,ls='--',color='k')
  plt.axvline(x=(-70-550)/1e3,ls='--',color='r')
  plt.ylim(1e2,1e7)
  plt.xlim(-1.5,0.5)
  plt.xlabel('z [Mm]')
  plt.legend()
  plt.ylabel(r'Density $[\text{g cm}^{-3}]$')
  plt.grid()
  plt.savefig("Model.png",dpi=300)
  plt.show()
def plot_rho():
  x = rM-550-70
  y = rhoM
  xs = np.linspace(-1500,0.0,100)
  spl = UnivariateSpline(-1*x, y)
  # -- Density
  plt.figure(figsize=(6.04,4.3))
  plt.semilogy(r_range/1e3,rho_range,'k',label='Model S')
  plt.semilogy((rM-550-70)/1e3,rhoM*0.66,'tab:orange',ls='--')
  plt.semilogy((rM-550-70)/1e3,rhoM*1,ls='--',c='b')
  plt.semilogy(r1c/1e3,rhoc,'g',ls='--',c='g')
  spl.set_smoothing_factor(0.5)
  #plt.plot(xs/1e3,spl(-xs),lw=4,color='limegreen')
  #plt.semilogy((rM-550-70)/1e3,rhoML*1,'g',ls='--')
  #plt.semilogy((rM-550-70)/1e3,rhoML*1,'m',ls='--')
  plt.axvline(x=-70/1e3,ls='--',color='k')
  plt.axvline(x=(-70-550)/1e3,ls='--',color='r')
  plt.ylim(1e-9,1e-5)
  plt.xlim(-1.5,0.5)
  plt.ylabel(r'Density $[\text{g cm}^{-3}]$')
  plt.xlabel('z [Mm]')
  #plt.grid()
  plt.legend()

# -- Weight Function --
z0 = -70
zu = z0 - 550
def wu(z):
  I2 = lambda z: np.cos(np.pi*(zu-80-z)/(300))/2 + 1/2
  I4 = lambda z: np.cos(np.pi*(zu+20-z)/(800-zu))/2 + 1/2
  wu = np.piecewise(z,[z>-1e6,z>zu-380,z>zu-80,z>zu+20,z>900],[0,I2,1,I4,0])
  return wu
def plot_weight():
  with plt.style.context(['science']):
    z = np.linspace(-1500,1000, 100)
    #plt.style.use(['science','std-colors'])#,'nature'])#,'ieee'])
    plt.axvline(x=(zu-380)/1e3,ls='--',color='b'     ,label='$z_u-380$ km')
    plt.axvline(x=(zu- 80)/1e3,ls='--',color='orange'     ,label='$z_u-80$ km')
    plt.axvline(x=(zu+ 20)/1e3,ls='--',color='g',label='$z_u+20$ km')
    plt.axvline(x=(   900)/1e3,ls='--',color='m'     ,label='$900$ km')
    plt.plot(z/1e3,wu(z),'-',color='k')
    plt.axvline(x=(z0)/1e3,ls='--',color='k')
    plt.axvline(x=(zu)/1e3,ls='--',color='r')
    plt.xlabel("Depth [Mm]")
    plt.ylabel("$w_u(z)$")
    plt.legend(frameon=True,fontsize=8)
    plt.savefig("Weight_Function.png",dpi=300)
    plt.savefig("Weight_Function.pdf")
    plt.savefig("Weight_Function.eps")
    plt.show()
#plot_weight()


z = np.linspace(-1500,500,1000)
x2int   = rM-550-70
rho2int, prs2int       = rhoM, pM
rho_int_ME, prs_int_ME = spline(-x2int,rho2int,k=2), spline(-x2int,prs2int,k=2)
P_int_S   = spline(-r1,p)
rho_int_S = spline(-r1,rho)

# -- Termodynamical variables as a function of z
def Pqs(z):
  return P_int_S(-z)
def rhoqs(z):
  return rho_int_S(-z)
def PMaltby(z):
  return prs_int_ME(-z)
def rhoMaltby(z):
  return rho_int_ME(-z)

# --------------------------------------------------------------------------
# -- Umbrel thermodynamical variables as function of depth and plots
def Pu(z):
  I1 = lambda z: Pqs(z)
  I2 = lambda z: np.exp( (1-wu(z))*np.log(Pqs(z)) + wu(z)*np.log(PMaltby(z)) )
  ss = np.piecewise(z,[z<-1000,z>=-1000],[I1,I2])
  return ss
def rhou(z,zz=-675):
  I1 = lambda z: rhoqs(z)
  I2 = lambda z: np.exp( (1-wu(z))*np.log(rhoqs(z)) + wu(z)*np.log(0.66*rhoMaltby(z)) )
  ss = np.piecewise(z,[z<zz,z>=zz],[I1,I2])
  return ss

def plot_density_u():
  with plt.style.context(['science']):
    plt.axvline(x=-70/1e3,ls='--',color='k')
    plt.axvline(x=(-70-550)/1e3,ls='--',color='r')
    plt.semilogy(z/1e3,rhoqs(z),'k',label='Model S')
    plt.semilogy(x2int/1e3,rho2int,'o-',color='orange',label='Maltby E')
    plt.semilogy(z/1e3,rhou(z),'b',ls='-.',label=r'$\rho_u(z)$')
    plt.ylim(1e-9,1e-5),plt.xlim(-1.5,0.5),plt.grid()
    plt.ylabel("Density [g cm$^{-3}$]"),plt.xlabel("Depth [Mm]")
    plt.legend(frameon=True,fontsize=8)
    plt.savefig("Density.png",dpi=300)
    plt.show()
plot_density_u()
def plot_pressure_u():
  with plt.style.context(['science']):
    plt.semilogy(z/1e3,Pqs(z),'k',label='Model S')
    plt.semilogy(x2int/1e3,prs2int,'o-',color='orange',label='Maltby E')
    #plt.semilogy(z/1e3,PMaltby(z),'b')
    plt.semilogy(z/1e3,Pu(z),'b',ls='-.',label='$P_u(z)$')
    plt.axvline(x=-70/1e3,ls='--',color='k')
    plt.axvline(x=(-70-550)/1e3,ls='--',color='r')
    #plt.ylim(1e-9,1e-5)
    plt.ylim(1e2,1e7)
    plt.xlim(-1.5,0.5)
    plt.grid()
    plt.xlabel("Depth [Mm]")
    plt.ylabel("Pressure [dyn]")
    plt.legend(frameon=True,fontsize=8)
    plt.savefig("Pressure.png",dpi=300)
    plt.show()
#plot_pressure_u()
# --------------------------------------------------------------------------



# -- Interpolation for umbral pressure ---------------------------------
def interp_prs_umbra():
  zi = np.linspace(-1500,-650,1000)
  zf = np.linspace(-650,-300,1000)
  zg = np.linspace(-300,0,1000)
  zh = np.linspace(0,500,1000)
  P2intpluto = np.polyfit(zi,Pu(zi),25)
  P2intpluto2= np.polyfit(zf,Pu(zf),22)
  P2intpluto3= np.polyfit(zg,Pu(zg),10)
  P2intpluto4= np.polyfit(zh,Pu(zh),15)
  P_int_pluto = np.poly1d(P2intpluto)
  P_int_pluto2= np.poly1d(P2intpluto2)
  P_int_pluto3= np.poly1d(P2intpluto3)
  P_int_pluto4= np.poly1d(P2intpluto4)
  
  plt.subplot(121)
  #plt.semilogy(zi/1e3,rhoqs(zi))
  plt.semilogy(zi/1e3,Pu(zi),'k')
  plt.semilogy(zf/1e3,Pu(zf),'k')
  plt.semilogy(zg/1e3,Pu(zg),'k')
  plt.semilogy(zh/1e3,Pu(zh),'k')
  plt.semilogy(zi/1e3,P_int_pluto(zi))
  plt.semilogy(zf/1e3,P_int_pluto2(zf))
  plt.semilogy(zg/1e3,P_int_pluto3(zg))
  plt.semilogy(zh/1e3,P_int_pluto4(zh))
  plt.ylim(1e2,1e7)
  plt.xlim(-1.6,0.5),plt.grid()
  
  plt.subplot(122)
  plt.plot(zi/1e3,100*(Pu(zi)-P_int_pluto(zi)) /Pu(zi))
  plt.plot(zf/1e3,100*(Pu(zf)-P_int_pluto2(zf))/Pu(zf))
  plt.plot(zg/1e3,100*(Pu(zg)-P_int_pluto3(zg))/Pu(zg))
  plt.plot(zh/1e3,100*(Pu(zh)-P_int_pluto4(zh))/Pu(zh))
  plt.ylim(-1.5,1.5), plt.xlim(-1.6,0.5), plt.grid()
  plt.show()
  

# -- Interpolation for umbral density ---------------------------------
def interp_rho_umbra():
  zi = np.linspace(-1500,-650,1000)
  zf = np.linspace(-650,-300,1000)
  zg = np.linspace(-300,0,1000)
  zh = np.linspace(0,500,1000)
  rho2intpluto = np.polyfit(zi,rhou(zi),19)
  rho2intpluto2= np.polyfit(zf,rhou(zf),19)
  rho2intpluto3= np.polyfit(zg,rhou(zg),10)
  rho2intpluto4= np.polyfit(zh,rhou(zh),10)
  rho_int_pluto = np.poly1d(rho2intpluto)
  rho_int_pluto2= np.poly1d(rho2intpluto2)
  rho_int_pluto3= np.poly1d(rho2intpluto3)
  rho_int_pluto4= np.poly1d(rho2intpluto4)
  
  var = rho2intpluto4
  for i in range(len(var)):
    if (i==len(var)-1):
      print("%.9e*pow(y,%d);"%(np.flip(var)[i],i),end=" ")
      break
    elif (i==0):
      print("%.9e +"%(np.flip(var)[i]),end=" ")
    elif (i==1):
      print("%.9e*y +"%(np.flip(var)[i]))
    elif (i%2 == 1):
      print("%.9e*pow(y,%d) +"%(np.flip(var)[i],i))
    else:
      print("%.9e*pow(y,%d) +"%(np.flip(var)[i],i),end=" ")
  print(" ")
  #exit()
  
  plt.subplot(121)
  #plt.semilogy(zi/1e3,rhoqs(zi))
  plt.semilogy(zi/1e3,rhou(zi),'k')
  plt.semilogy(zf/1e3,rhou(zf),'k')
  plt.semilogy(zg/1e3,rhou(zg),'k')
  plt.semilogy(zi/1e3,rho_int_pluto(zi))
  plt.semilogy(zf/1e3,rho_int_pluto2(zf))
  plt.semilogy(zg/1e3,rho_int_pluto3(zg))
  plt.semilogy(zh/1e3,rho_int_pluto4(zh))
  plt.ylim(1e-9,1e-5)
  plt.xlim(-1.6,0.5),plt.grid()
  
  plt.subplot(122)
  plt.plot(zi/1e3,100*(rhou(zi)-rho_int_pluto(zi))/rhou(zi))
  plt.plot(zf/1e3,100*(rhou(zf)-rho_int_pluto2(zf))/rhou(zf))
  plt.plot(zg/1e3,100*(rhou(zg)-rho_int_pluto3(zg))/rhou(zg))
  plt.plot(zh/1e3,100*(rhou(zh)-rho_int_pluto4(zh))/rhou(zh))
  plt.ylim(-1.5,1.5), plt.xlim(-1.6,0.5), plt.grid()
  plt.show()
  








# -- Temperature
#plt.figure(figsize=(6.04,4.3))
#plt.semilogy(r_range/1e3,T_range,'k',label='Model S')
##plt.semilogy((rM-550-70)/1e3,TM*0.66,'tab:cyan',ls='--')
#plt.semilogy((rM-550-70)/1e3,TM*1,'r',ls='--')
##plt.semilogy((rM-550-70)/1e3,rhoML*1,'g',ls='--')
##plt.semilogy((rM-550-70)/1e3,rhoML*1,'m',ls='--')
#plt.axvline(x=-70/1e3,ls='--',color='k')
#plt.axvline(x=(-70-550)/1e3,ls='--',color='r')
##plt.ylim(1e-9,1e-5)
#plt.xlim(-1.5,0.5)
#plt.xlabel('z [Mm]')
#plt.legend()
##plt.ylabel(r'Density $[\text{g cm}^{-3}]$')
#plt.grid()
#plt.show()

