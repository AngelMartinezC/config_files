import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os
import matplotlib
import warnings   # ------ librería que esconde los warnings
import multiprocessing
from joblib import Parallel, delayed
import time as tt
warnings.filterwarnings("ignore")
plt.rcParams.update({'font.size': 14})
os.system('clear')

# -- Constante
RSun    = 6.9634e10  # [cm]
MSun    = 1.989e33   # [g]
G_const = 6.67e-8    # cm^3 g^{-1} s^{-2}

# -- Read variables from model S Christensen-Dalsgaard and interpolate
r0,c0,rho0,p0,G0,T0 = np.loadtxt('cptrho.l5bi.d.15c', unpack=True)
r = np.linspace(0,1,1000)
rho = interp1d(r0, rho0, kind='cubic')
c = interp1d(r0, c0, kind='cubic')
p = interp1d(r0, p0, kind='cubic')
G1 = interp1d(r0, G0, kind='cubic')

# -- Calculate mass from standard model
Mr, r_new  = [], []
Mass = 0
for i in range(0,len(r)):
  if i==len(r)-1:
    Mr.append(Mass)
    r_new.append(r[i])
  else:  
    Mass += 4*np.pi*r[i]**2*rho(r[i])*(r[i+1]-r[i])*RSun**3
    Mr.append(Mass)
    r_new.append(r[i])
M = interp1d(r_new,Mr, kind='cubic')

# -- Calculate gravity
def g(r):
  return G_const*M(r)/(r**2*RSun**2)

# -- Calculate derivatives of ln(p), ln(rho) for bouyancy and interpolate
dp_dr, drho_dr = [], []
for i in range(0,len(r)):
  if i==0:
    h = r[i+1]-r[i]
    up_p   = p(r[i]+h)   - p(r[i])
    up_rho = rho(r[i]+h) - rho(r[i])
    dp_dr.append(up_p/(RSun*h))
    drho_dr.append(up_rho/(RSun*h))
  elif i==len(r)-1:
    h = r[i]-r[i-1]
    up_p   = p(r[i])   - p(r[i]-h)
    up_rho = rho(r[i]) - rho(r[i]-h)
    dp_dr.append(up_p/(RSun*h))
    drho_dr.append(up_rho/(RSun*h))
  else:
    h = r[i+1]-r[i]
    up_p = p(r[i]+h)     - p(r[i]-h)
    up_rho = rho(r[i]+h) - rho(r[i]-h)
    dp_dr.append(up_p/(RSun*2*h))
    drho_dr.append(up_rho/(RSun*2*h))
dpdr   = interp1d(r_new,dp_dr,kind='cubic')
drhodr = interp1d(r_new,drho_dr,kind='cubic')

# -- Lamb and Boyuancy frequencies
def S_l(l):
  return np.sqrt(l*(l+1))*(c(r)/(r*RSun))
N = np.sqrt(g(r) * (dpdr(r)/(p(r)*G1(r))-drhodr(r)/rho(r)) )
def L(l):
  return np.sqrt(l*(l+1))
Hp = abs((-1*dpdr(r)/p(r))**(-1))

# -- Plots of frequencies
def Pressure_Scale_Height():
  plt.figure(figsize=(8,6))
  plt.plot(r,Hp,'b-')
  plt.yscale('log')
  plt.xscale('log')
  plt.grid()
  plt.xlabel(r'$r/R\odot$')
  plt.ylabel(r'Pressure scale height,  $H_p=\left(\frac{d\ln p}{dr}\right)^{-1}$')
  plt.xlim([1e-4,1.3])
  plt.show()

def Lamb_Bouyancy():
  plt.figure(figsize=(8,6))
  ax1 = plt.subplot()
  l = [1, 2, 5, 10, 50, 100, 200]
  #l = np.arange(5,100,1)
  for i in l:
    plt.plot(r,1e3*S_l(i)/(2*np.pi),'-',label=r'$l=$'+str(i),lw=2)
  plt.plot(r,1e3*N/(2*np.pi),'--',lw=2)
  plt.legend(loc=8,fontsize=11)
  plt.yscale('log')
  plt.yticks([0.1,0.5,1,5])
  plt.yticks([1e-3,0.01,0.1,0.5,1,5],[1e-3,0.01,0.1,0.5,1.0,5.0])
  plt.xlabel(r'$r/R_\odot$')
  plt.ylabel(r'Frequency $\nu$ [mHz]')
  ax1.get_xaxis().set_tick_params(which='minor', size=10)
  ax1.get_yaxis().set_tick_params(which='minor', size=4, direction='in',
      right=True,labelright=True)
  ax1.get_yaxis().set_tick_params(which='major', size=8, direction='in',
      right=True,labelright=True)
  ax1.tick_params(labelright=False, axis='y', grid_linewidth=1)
  ax1.get_xaxis().set_tick_params(which='minor', width=0)
  plt.grid()
  plt.ylim([4e-2,5])
  plt.show()
  return None

#Pressure_Scale_Height()
#Lamb_Bouyancy()

# -- End of frequency plots -----------------------------------------------


# -- Turning Point
def rt(omega,l,rmin=0.8):
  def r1(omega,l,r0):
    y = c(r0)*L(l)/omega
    return y
  rad = np.linspace(rmin,1,1000)
  inter = []
  for i in rad:
    rad0 = r1(omega,l,i)
    rad1 = i*RSun
    y = rad0-rad1
    inter.append(y)
  intersection = np.array(inter)
  idx = np.argwhere(np.diff(np.sign(intersection - rad))).flatten()
  if len(idx)==0:
    return np.nan
  else:
    return rad[idx][0]


# -- Plots turning point vs (frequency and degree number)
def plot_depth_vs_l():
  l = np.arange(1,1101,1)
  def depth_vs_l(omega=0.005):
    depth = []
    for i in l:
      depth.append(rt(omega,l=i,rmin=0))
    return np.array(depth)
  plt.figure(figsize=(8,6))
  for i in np.arange(0.001,0.012,0.001):
    ff = i*1000
    print('Making frequency \u03BD=%.2f mHz'%ff)
    plt.plot(l, depth_vs_l(i*2*np.pi),label=r'$\nu$= %.1f mHz'%(ff))
  plt.legend(fontsize=12)
  plt.xlabel(r'Degree number $l$')
  plt.ylabel(r'Turning point: (1-$r_t$) $r/R_\odot$')
  plt.xscale('log')
  plt.grid()
  plt.show()

def plot_depth_vs_w():
  w = np.arange(0.001,0.012,0.001)
  def depth_vs_w(l=100):
    print("l number",l)
    depth = []
    for i in w:
      print("now is the turn of %.3f mHz"%i)
      depth.append(rt(i*2*np.pi,l,rmin=0.1))
    return np.array(depth)
  plt.figure(figsize=(8,6))
  for i in np.arange(10,110,10):
    print('\nMaking l=%f'%i)
    plt.plot(w*1000, depth_vs_w(i),'.-',label=r'$l=$ %d '%(i))
  plt.legend(fontsize=12)
  plt.xlabel(r'Frequency $\nu$ [mHz]')
  plt.ylabel(r'Turning point: (1-$r_t$) $r/R_\odot$')
  plt.grid()
  plt.show()

#plot_depth_vs_l()
#plot_depth_vs_w()





# -- RAY PATH APPROXIMATION -----
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -- Calculate integrals for radius, theta and time (ray-path approximation)
def path_inside(w=0.003,l=10,n=12000,rmin=0.8):
  """
    Calculate integrals of radius and theta as a function of time according
    to the ray-path approximation.
    This function returns the arrays of the path traveled by the p-mode
    inside the Sun, i.e., plots can be made of these paths
  """
  turn = rt(w,l,rmin=rmin)
  rint0 = np.linspace(turn,1,int(n/2))
  rint1 = np.linspace(1,2-1.0000*turn,int(n/2))
  #rint1 = np.linspace(1,2-turn,int(n/2))
  rint = np.array(list(rint0)+list(rint1))
  theta, radius, time = [], [], []
  th,t = 0, 0
  for i in range(len(rint)):
    # First data (do not edit)
    if i==0:
      th += 0
      t += 0
      theta.append(th)
      radius.append(rint[i])
      time.append(t)
    # rest of the data
    else:
      if rint[i]<=1: # If radius is lesser than the solar
        fraction = (w*rint[i]*RSun/(L(l)*c(rint[i])))**2
      else: # Otherwise, make calculation downwards
        fraction = (w*(2-rint[i])*RSun/(L(l)*c(2-rint[i])))**2
  
      if fraction < 1:  # Ensure negative values do not contribute
        th += 0.0
        t  += 0.1
        theta.append(th)
        time.append(t)
        if rint[i] <= 1:
          radius.append((rint[i]))
        else:
          radius.append((2-rint[i]))
        pass
  
      else: # Append values given by the integration
        if i==len(rint)-1:
          dr = rint[i]-rint[i-1]
        else:
          dr = rint[i+1]-rint[i]
        if rint[i]<=1:
          root  = np.sqrt(fraction-1)
          roott = np.sqrt(1-(1/fraction))
          th += (1/root)*dr/rint[i]
          t  += (1/roott)*dr*RSun/c(rint[i])
          radius.append(rint[i])
        else:
          root = np.sqrt(fraction-1)
          roott = np.sqrt(1-(1/fraction))
          th += (1/root)*dr/(2-rint[i])
          t  += (1/roott)*dr*RSun/c(2-rint[i])
          radius.append((2-rint[i]))
        theta.append(th) 
        time.append(t) 
  return np.array(radius), np.array(theta), np.array(time)

def plot_TD(single=True, bounces=False):
  """
    Plot of the TD relation with ray path approximation. For the first part,
    I only show the reflection of the second bounce given the 5 h interval.
    -- single is for a single TD relation.
    -- bounces has 6 skips on the solar surface
  """
  count = 0
  if single:
    plt.figure(figsize=(10,7))
    TIME, DIST = [], []
    for f in [0.003,0.005,0.007,0.009,0.011]:
      SD = np.arange(500,1201,5)
      for l in SD:
        w = 2*np.pi*f
        r, th, t = path_inside(w,l,n=6000,rmin=0.95)
        f0 = f*1000
        tot = len(SD)*5
        part = count/tot*100
        print('Frequency \u03BD=%.1f mHz. Degree l: %d. (%.2f %s)'
            %(f0,l,part,chr(37)))
        TIME.append(t[-1]/60)
        DIST.append(th[-1]*RSun/100/1e6)
        count += 1

    TIME = np.array(TIME)
    DIST = np.array(DIST)
    plt.plot(DIST,TIME,'.',color='yellow',markersize=7)
    RES = np.array([TIME])
    np.savetxt('ray_path.txt',np.transpose(RES))
    plt.savefig('ray_path_transparent.png',transparent=True)
    plt.xlabel('Distance [Mm]')
    plt.ylabel('Time [min]')
    #plt.legend()
    #plt.xlim([0,180])
    plt.ylim([0,60])
    plt.grid()
    plt.tight_layout()
    plt.savefig('ray_path_approximation.png',transparent=False)
    plt.savefig('ray_path_approximation.pdf',transparent=False)
    plt.show()


  if bounces:
    XX, YY = [], []
    for f in [0.003,0.005,0.007,0.009,0.011]:
      TIME, DIST, ARRAY = [], [], []
      TIME2, DIST2      = [], []
      TIME3, DIST3      = [], []
      TIME4, DIST4      = [], []
      TIME5, DIST5      = [], []
      TIME6, DIST6      = [], []
      SD = np.arange(10,1201,5)
      for l in SD:
        w = 2*np.pi*f
        f0 = f*1000
        tot = len(SD)*5
        part = count/tot*100
        print('Frequency \u03BD=%.1f mHz. Degree l: %d. (%.2f %s)'
            %(f0,l,part,chr(37)),end=" ")
        if f==0.011:
          r, th, t = path_inside(w,l,n=6000,rmin=0.1)
        elif f==0.009:
          r, th, t = path_inside(w,l,n=6000,rmin=0.12)
        elif f==0.007:
          r, th, t = path_inside(w,l,n=6000,rmin=0.15)
        elif f==0.005:
          r, th, t = path_inside(w,l,n=6000,rmin=0.2)
        elif f==0.003:
          r, th, t = path_inside(w,l,n=6000,rmin=0.3)
        if 2*th[-1]>np.pi:
          DIST2.append((2*np.pi-2*th[-1])*180/np.pi)
        else:
          DIST2.append(2*th[-1]*180/np.pi)#*RSun/100/1e6)
        if np.isnan(th[-1]):
          print(' -->  is nan')
        else:
          print(" ")
        TIME.append(t[-1]/60)
        DIST.append(th[-1]*180/np.pi)#*RSun/100/1e6)
        TIME2.append(2*t[-1]/60)
        TIME3.append(3*t[-1]/60)
        DIST3.append(3*th[-1]*180/np.pi)#*RSun/100/1e6)
        TIME4.append(4*t[-1]/60)
        DIST4.append(4*th[-1]*180/np.pi)#*RSun/100/1e6)
        TIME5.append(5*t[-1]/60)
        DIST5.append(5*th[-1]*180/np.pi)#*RSun/100/1e6)
        TIME6.append(6*t[-1]/60)
        DIST6.append(6*th[-1]*180/np.pi)#*RSun/100/1e6)
        count += 1
      plt.plot(DIST,TIME,'b.',markersize=4,label=r'$\nu=$%.1f mHz'%f0)
      plt.plot(DIST2,TIME2,'.',color='darkorange',markersize=4,label=r'$\nu=$%.1f mHz'%f0)
      plt.plot(DIST3,TIME3,'.',color='green',markersize=4,label=r'$\nu=$%.1f mHz'%f0)
      plt.plot(DIST4,TIME4,'.',color='r',markersize=4,label=r'$\nu=$%.1f mHz'%f0)
      plt.plot(DIST5,TIME5,'.',color='lime',markersize=4,label=r'$\nu=$%.1f mHz'%f0)
      plt.plot(DIST6,TIME6,'.',color='gold',markersize=4,label=r'$\nu=$%.1f mHz'%f0)
    plt.xlabel(r'Distance [degree, $(^\circ)$]')
    plt.ylabel('Time [min]')
    plt.xlim([0,180])
    plt.ylim([0,400])
    plt.tight_layout()
    plt.show()
  return None

#plot_TD()

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------





# -- Plots paths of rays in a polar plot (like cover)
def plot_many_paths():
  l = 10
  w = 0.011
  fig, ax = plt.subplots(subplot_kw={'projection': 'polar'},figsize=(9,9))
  colors = ['b','r','g','purple','olive','orange','cyan','lime']
  cc = 0
  for l in [5,10,20,30,50,70,80,90]:
    radius, theta, time = path_inside(w=w, l=l)
    print("Making plot of l=%d"%l)
    for i in range(6):
      if i==0:
        plt.plot(theta+i*np.max(theta),radius,color=colors[cc],label=r'$l=$%d'%l)
      else:
        plt.plot(theta+i*np.max(theta),radius,color=colors[cc])
    plt.tight_layout()
    name = np.arange(0,2*np.pi,np.pi/4)
    plt.xticks(name,['']*len(name))
    lines, labels = plt.rgrids( [0,0.25,0.5,0.75,1], ['','','','','R'])
    cc += 1
  plt.legend()
  plt.show()

#plot_many_paths()

# Function to make one skip for a given frequency and degree number.
def path_for_video(w=0.005,l=30,rmin=0.8,skips=1):
  """
    Function to make one skip for a given frequency and degree number.
    Here integrals of time, theta and time are made in order to return a
    specific number of skips. By default, just one skip (both right and right)
    is plotted. 
  """
  radius, theta, time = path_inside(w=w, l=l, rmin=rmin)
  idx = np.argmax(radius)
  tiempo00 = np.array(time[0:idx])
  nn = 0
  #for k in range(len(tiempo00)-1):
  #  if tiempo00[k] == tiempo00[k+1]:
  #    nn += 1
  #  else:
  #    break
  tiempo0 = np.array(time[0:idx])
  radio0 = np.array(radius[0:idx])
  radio1 = np.array(radius[idx::])
  theta0 = np.array(theta[0:idx])-max(theta[0:idx])+min(theta[0:idx])
  theta1 = np.array(theta[idx::])-max(theta[idx::])+min(theta[idx::])
  theta2 = -np.array(theta[idx::])+max(theta[idx::])-min(theta[idx::])
  theta3 = theta0 + 2*max(theta1)
  theta4 = -theta0 - 2*max(theta1) 
  _ = np.array(time[idx::])
  tiempo1 = _ - min(_)
  tiempo2 = tiempo0 + max(tiempo1)
  t_t   = np.array(list(tiempo1)+list(tiempo2))
  r_t   = np.array(list(radio1)+list(radio0))
  th_iz = np.array(list(theta1)+list(theta3))
  th_de = np.array(list(theta2)+list(theta4))

  thetad, thetai, radio, tiempo = [], [], [], []
  for i in range(skips):
    thetad += list(th_de-i*max(abs(th_de)))
    thetai += list(th_iz+i*max(th_iz))
    radio  += list(r_t)
    tiempo += list(t_t+i*max(t_t))
  #print(tiempo)
  return np.array(radio),np.array(thetai),np.array(thetad),np.array(tiempo,dtype=float)

#rr,thi,thd,tt=path_for_video(0.005, 440,rmin=0.6,skips=8)
#plt.plot(rr)
#rr,thi,thd,tt=path_for_video(0.005, 440,rmin=0.9,skips=7)
#plt.plot(rr) plt.show()
#print(len(rr))
#print(rr)
#exit()

# temp is the array with the smallest maximum time, tmin is a seed
# This will not check whether the minimun time corresponds to the 
# array with the minimum size
freq = 0.005
radio, thiz, thde, tiempo, tmin, ctemp = [], [], [], [], 1e8, 0
skips_arr = [2,3,3,4,5,7,8,12,13,14]
depth_l = [0.6,0.95,0.97,0.97,0.98,0.98,0.99,0.994,0.997]
l_array = [100, 300, 400, 500, 700, 900,1000,1200, 1700]#,1000,1200]
ccc = ['orange','g','r','b','purple','gold','lime','deeppink','cyan']
modes = len(l_array)
T = 3600
tt = np.arange(0,T+1,1)
for l in l_array:
  #if l==440:
  #  r, thi, thd, t = path_for_video(2*np.pi*freq,l,rmin=0.95,skips=skips_arr[ctemp])
  #elif l==300:
  #  r, thi, thd, t = path_for_video(2*np.pi*freq,l,rmin=0.77,skips=skips_arr[ctemp])
  #elif l==200:
  #  r, thi, thd, t = path_for_video(2*np.pi*freq,l,rmin=0.6,skips=skips_arr[ctemp])
  #else:
  #  r, thi, thd, t = path_for_video(2*np.pi*freq,l,rmin=0.6,skips=skips_arr[ctemp])
  r, thi, thd, t = path_for_video(2*np.pi*freq,l,rmin=depth_l[ctemp],skips=skips_arr[ctemp])
  conta = 1
  for i in range(len(t)-1): # if duplicate
    if t[i+1] == t[i]:
      t[i+1] = t[i+1]+0.00000001*conta
      conta += 1
    else: pass
  print(l,max(t))
  fr   = interp1d(t,r,kind='linear')
  fthd = interp1d(t,thd,kind='linear')
  fthi = interp1d(t,thi,kind='linear')
  radio.append(fr(tt))
  thiz.append(fthi(tt)), thde.append(fthd(tt))
  tiempo.append(tt)
  if np.max(t)<tmin:
    tmin = np.max(t)
    temp = ctemp
  else: pass
  ctemp += 1

var_dir = {}
a = 81*np.pi/180
t0 = np.linspace(a,np.pi-a,1000)
nn = 10

for i in range(1,modes+1):
  var_dir['p_r%d'%i] = []
  var_dir['p_th%diz'%i], var_dir['p_th%dde'%i] = [], []


# -- Begin of plot --
counter = 0
ttt = tiempo[temp]
for i in range(0,len(ttt),1):
  if i==len(ttt)-1:
    dt = ttt[i]-ttt[i-1]
  else:
    dt = ttt[i+1]-ttt[i]
  f=plt.figure(figsize=(16,5.5))
  f.subplots_adjust(bottom=0.07, left=0.05, right=0.96, top=0.95, \
      wspace=0.0, hspace=0.0)
  stamp = int(ttt[i]//nn)

  for j in range(1,modes+1):
    p_r = var_dir['p_r%d'%j]
    p_thiz = var_dir['p_th%diz'%j]
    p_thde = var_dir['p_th%dde'%j]

    p_r.append(radio[j-1][i])
    p_thiz.append(thiz[j-1][i]), p_thde.append(thde[j-1][i])
    ap_r    = np.array(p_r)
    ap_thiz, ap_thde = np.array(p_thiz), np.array(p_thde)
    x0, y0 = ap_r*np.cos(ap_thiz+np.pi/2), ap_r*np.sin(ap_thiz+np.pi/2)
    x2, y2 = ap_r*np.cos(ap_thde+np.pi/2), ap_r*np.sin(ap_thde+np.pi/2)
    if j==1:
      plt.plot([0],[0],'white',label='Time: %4.0f s'%ttt[i])
      plt.plot(x0,y0,'-',linewidth=1.5, color=ccc[j-1])
      plt.plot(x2,y2,'-',linewidth=1.5, color=ccc[j-1])
    else:
      plt.plot(x0,y0,'-',linewidth=1.5, color=ccc[j-1])
      plt.plot(x2,y2,'-',linewidth=1.5, color=ccc[j-1])
  
  for k in range(modes):
    plt.plot([0],[0],color=ccc[k],label='$l=$%d'%l_array[k])

  plt.plot(np.cos(t0),np.sin(t0),'k')
  plt.plot(0.95*np.cos(t0),0.95*np.sin(t0),ls='dotted',color='gray',lw=1.0)
  plt.plot(0.96*np.cos(t0),0.96*np.sin(t0),'--',color='gray',lw=1)
  plt.plot(0.97*np.cos(t0),0.97*np.sin(t0),ls='dotted',color='gray',lw=1.0)
  plt.plot(0.98*np.cos(t0),0.98*np.sin(t0),'--',color='gray',lw=1)
  plt.plot(0.99*np.cos(t0),0.99*np.sin(t0),ls='dotted',color='gray',lw=1.0)
  plt.legend(loc=3)
  freq_mHz = freq*1000
  plt.text(1.007*np.cos(a+0.03),1.007*np.sin(a+0.03),r'$\nu$=%.1f mHz'%freq_mHz)
  plt.axes().set_aspect('equal')
  plt.ylim([0.94,1.01])
  plt.xlim([-1.1*np.cos(a),1.1*np.cos(a)])
  if i==0:
    plt.savefig('IMAGES/pic.%04ds.png'%int(ttt[i]))
  else:
    if stamp == counter:
      pass
    else:
      plt.savefig('IMAGES/pic.%04ds.png'%int(counter+1))
      counter += 1
  porc = i/len(ttt)*100
  print('Time: %.2f s.  Dt: %.4f s.  Frame %d of %d (%.2f %s)  Time stamp %d. TT %f'
         %(ttt[i],dt,i,len(ttt),porc,chr(37),stamp,tiempo[0][i]))
  plt.close('all')





