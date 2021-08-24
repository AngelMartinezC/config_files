import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pyPLUTO as pypl
import pyPLUTO.pload as pp
import pyPLUTO.Image as img
import pyPLUTO.Tools as tl
import multiprocessing
from joblib import Parallel, delayed
import os
import time
from matplotlib.ticker import StrMethodFormatter
#plt.style.use(['science','std-colors'])#,'nature'])#,'ieee'])


# -- To save in folder
if os.path.exists('IMAGES'): print('Folder exists.')
else: os.system('mkdir IMAGES')

# Normalization units
unit_v = 1e5
unit_rho = 1e-6
unit_b = np.sqrt(4*np.pi*unit_rho*unit_v**2)

wdir = 'dbl_files/'
X0 = float(os.popen("grep 'X1' '"+wdir+"/grid.out' | tr -d '[a-z]' | \
    awk '{print $3}' | sed 's/,//g'").read())
X1 = float(os.popen("grep 'X1' '"+wdir+"/grid.out' | tr -d '[a-z]' | \
    awk '{print $4}' | sed 's/,//g'").read())
Y0 = float(os.popen("grep 'X2' '"+wdir+"/grid.out' | tr -d '[a-z]' | \
    awk '{print $3}' | sed 's/,//g'").read())
Y1 = float(os.popen("grep 'X2' '"+wdir+"/grid.out' | tr -d '[a-z]' | \
    awk '{print $4}' | sed 's/,//g'").read())
resX = 5
resY = 3
XX = np.around(np.linspace(X0,X1,resX)/1e3,2).astype(int)
YY = np.around(np.linspace(Y0,Y1,resY)/1e3,2)#.astype(int)
nlinf = pypl.nlast_info(w_dir=wdir)#, datatype='flt')
number = nlinf['nlast']
SSS = []*number

def image(num,plot=False): #i=0, cmap='jet', diff=False,unit=1,res=25,plot=False):
  D, I = pp.pload(num, w_dir=wdir), img.Image()
  vmin =  None #np.log10(0.2)  #np.log10(20)   #None  #0 #None #0.2
  vmax =  None #np.log10(5.0)  #np.log10(730.0)   #None #13 #None #5.0
  VAR  =  (D.rho).T #D.Bx2*unit_b #np.log10(D.rho)  #D.Bx2*unit_b
  VX1  =  (D.vx1).T
  VX2  =  (D.vx2).T
  VX3  =  (D.vx3).T
  BX1  =  (D.Bx1).T
  BX2  =  (D.Bx2).T
  BX3  =  (D.Bx3).T
  
  DIVB1 = np.copy(VAR)*0
  DIVB2 = np.copy(VAR)*0
  DIVB3 = np.copy(VAR)*0
  DIVV1 = np.copy(VAR)*0
  DIVV2 = np.copy(VAR)*0
  DIVV3 = np.copy(VAR)*0
  
  # in y
  i = 40
  for i in range(len(VAR[0,:])):
    for j in range(len(VAR[:,0])):
      if j==0:
        DX2 = (D.x2[j+1] - D.x2[j])*1e5
        DIVV2[j,i] = (VX2[j+1,i] - VX2[j,i])/DX2
        DIVB2[j,i] = (BX2[j+1,i] - BX2[j,i])/DX2
      else:
        DX2 = (D.x2[j] - D.x2[j-1])*1e5
        DIVV2[j,i] = (VX2[j,i] - VX2[j-1,i])/DX2
        DIVB2[j,i] = (BX2[j,i] - BX2[j-1,i])/DX2
  
  for j in range(len(VAR[:,0])):
    for i in range(len(VAR[0,:])):
      if i<len(VAR[0,:])-1:
        DX1 = (D.x1[i+1] - D.x1[i])*1e5
        DIVV1[j,i] = (VX1[j,i+1] - VX1[j,i])/DX1
        DIVB1[j,i] = (BX1[j,i+1] - BX1[j,i])/DX1
      else:
        DX1 = (D.x1[j] - D.x1[i-1])*1e5
        DIVV1[j,i] = (VX1[j,i] - VX1[j,i-1])/DX1
        DIVB1[j,i] = (BX1[j,i] - BX1[j,i-1])/DX1
      
  fig = plt.figure(figsize=(10,9))
  fig.subplots_adjust(bottom=0.08, left=0.12, right=0.96, top=0.95, \
      wspace=0.0, hspace=0.1)
  
  unit_mag = 3.545e+02
  ax1 = plt.subplot(311)
  minute, second = divmod(D.SimTime,60)
  hour, minute   = divmod(minute,60)
  day, hour      = divmod(hour,24)
  plt.title('Time: {} days  {} h  {} min {} s'.format(int(day),int(hour),
      int(minute),round(second,2)), fontsize=20)
  plt.imshow(BX1*unit_mag,origin='lower',cmap='jet',aspect='auto', \
      extent=[D.x1.min(),D.x1.max(),D.x2.min(),D.x2.max()])
  plt.colorbar(label=r'$B_x$ [Gauss]')
  plt.xticks(np.linspace(D.x1[0],D.x1[-1],resX),[]*len(XX))
  plt.yticks(np.linspace(D.x2[0],D.x2[-1],resY),YY)
  plt.grid(ls='--')
  
  ax2 = plt.subplot(312)
  plt.imshow(BX2*unit_mag,origin='lower',cmap='jet',aspect='auto', \
      extent=[D.x1.min(),D.x1.max(),D.x2.min(),D.x2.max()],vmin=0,vmax=400)
  plt.colorbar(label=r'$B_y$ [Gauss]')
  plt.xticks(np.linspace(D.x1[0],D.x1[-1],resX),[]*len(XX))
  plt.yticks(np.linspace(D.x2[0],D.x2[-1],resY),YY)
  plt.ylabel(r"Depth [Mm]")
  plt.grid(ls='--')
  
  ax3 = plt.subplot(313)
  plt.imshow((DIVB1+DIVB2)*unit_mag,origin='lower',cmap='jet',aspect='auto', \
      extent=[D.x1.min(),D.x1.max(),D.x2.min(),D.x2.max()])
  plt.xticks(np.linspace(D.x1[0],D.x1[-1],resX),XX)
  plt.yticks(np.linspace(D.x2[0],D.x2[-1],resY),YY)
  #plt.imshow(DIVV1+DIVV2,origin='lower',cmap='jet',aspect='auto')#,vmin=-0.0075,vmax=0.0075)
  plt.colorbar(label=r'$\nabla \cdot \vec{B}$ [G cm$^{-1}$]')
  plt.xlabel(r"$x$ [Mm]")
  plt.grid(ls='--')
  
  #plt.tight_layout()
  string = '%0.4d' % (num)  # To save CUSTOM FRAME With name pic.000#.png
  print("Image %04d"%num)
  plt.savefig('IMAGES_DIV_MAG/pic.'+string+'.png')#,dpi=200)


  fig = plt.figure(figsize=(10,9))
  fig.subplots_adjust(bottom=0.08, left=0.12, right=0.96, top=0.95, \
      wspace=0.0, hspace=0.1)
  
  norm_v = 1.000e+03
  ax1 = plt.subplot(311)
  plt.title('Time: {} days  {} h  {} min {} s'.format(int(day),int(hour),
      int(minute),round(second,2)), fontsize=20)
  plt.imshow(VX1*norm_v,origin='lower',cmap='jet',aspect='auto', \
      extent=[D.x1.min(),D.x1.max(),D.x2.min(),D.x2.max()])
  plt.colorbar(label=r'$v_x$ [m s$^{-1}$]')
  plt.xticks(np.linspace(D.x1[0],D.x1[-1],resX),[]*len(XX))
  plt.yticks(np.linspace(D.x2[0],D.x2[-1],resY),YY)
  plt.grid(ls='--')
  
  ax2 = plt.subplot(312)
  plt.imshow(VX2*norm_v,origin='lower',cmap='jet',aspect='auto', \
      extent=[D.x1.min(),D.x1.max(),D.x2.min(),D.x2.max()])
  plt.colorbar(label=r'$v_y$ [m s$^{-1}$]')
  plt.xticks(np.linspace(D.x1[0],D.x1[-1],resX),[]*len(XX))
  plt.yticks(np.linspace(D.x2[0],D.x2[-1],resY),YY)
  plt.ylabel(r"Depth [Mm]")
  plt.grid(ls='--')
  
  ax3 = plt.subplot(313)
  plt.imshow((DIVV1+DIVV2)*norm_v*100,origin='lower',cmap='jet',aspect='auto', \
      extent=[D.x1.min(),D.x1.max(),D.x2.min(),D.x2.max()])
  plt.xticks(np.linspace(D.x1[0],D.x1[-1],resX),XX)
  plt.yticks(np.linspace(D.x2[0],D.x2[-1],resY),YY)
  #plt.imshow(DIVV1+DIVV2,origin='lower',cmap='jet',aspect='auto')#,vmin=-0.0075,vmax=0.0075)
  plt.colorbar(label=r'$\nabla \cdot \vec{v}$ [s$^{-1}$]')
  plt.xlabel(r"$x$ [Mm]")
  plt.grid(ls='--')
  
  #plt.tight_layout()
  string = '%0.4d' % (num)  # To save CUSTOM FRAME With name pic.000#.png
  print("Image %04d"%num)
  plt.savefig('IMAGES_DIV_VEL/pic.'+string+'.png')#,dpi=200)

  if plot: plt.show()
  else: plt.close("all")

nlinf = pypl.nlast_info(w_dir=wdir)#, datatype='flt')
number = nlinf['nlast']
#image(1381,plot=True)#,unit=354.490770,res=20)
#exit()

s = time.time()
print("Number of files is %d"%number)
num_cores = multiprocessing.cpu_count()
processed_list = Parallel(n_jobs=num_cores)(delayed(image)(i)
    for i in range(1769,number,1))
#    for i in range(1499,number,1))
print("Time spent {:.2f} s".format(time.time()-s))

