import os
import matplotlib.pyplot as plt
import pyPLUTO.pload as pp
import pyPLUTO.Image as img
import pyPLUTO.Tools as tl
import numpy as np
import multiprocessing
from joblib import Parallel, delayed
import time

I = img.Image()
wdir = 'dbl_files/'


def image(num, plot=False):
  plt.rcParams.update({'font.size': 13})
  D  = pp.pload(num,datatype='dbl',w_dir=wdir) 
  D0 = pp.pload(num+1,datatype='dbl',w_dir=wdir) 
  
  #xy
  cut = 224 #224
  VAR = D.rho[:,:,cut]
  I.pldisplay(D, (VAR), x1=D.x1, x2=D.x2,
      vmin=None, vmax=(np.mean(VAR)+np.max(VAR))/2,
      polar=[True, True], subplot=131, cmap='gnuplot',
      cbar=[False, 'vertical'],figsize=[17,7],
      label1='x', label2='y',title=r'Density [$\rho$]  $z=0$')
  T, newdims = tl.Tools(), (20,7)
  s0 = 15
  rr = np.linspace(D.x1[s0],D.x1[-1],len(D.x1[s0::]))
  tt = np.linspace(D.x2[0],D.x2[-1],len(D.x2))
  RRR, TTT = np.meshgrid(rr,tt)
  Xmesh, Ymesh = RRR*np.cos(TTT), RRR*np.sin(TTT)
  xcong = T.congrid(Xmesh,newdims,method='linear')
  ycong = T.congrid(Ymesh,newdims,method='linear')

  VX = D.vx1[s0::,:,cut].T*np.cos(TTT) - D.vx2[s0::,:,cut].T*np.sin(TTT)
  VY = D.vx1[s0::,:,cut].T*np.sin(TTT) + D.vx2[s0::,:,cut].T*np.cos(TTT)
  velxcong = T.congrid(VX,newdims,method='linear')
  velycong = T.congrid(VY,newdims,method='linear')
  plt.gca().quiver(xcong, ycong, velxcong, velycong,ec='w',fc='k',lw=0.8,\
      width=0.0030, headwidth=6, alpha=0.5)


  VAR = D.vx3[:,:,cut] - D0.vx3[:,:,cut]
  I.pldisplay(D, (VAR), x1=D.x1, x2=D.x2,
      #vmin=None, vmax=(np.mean(VAR)+np.max(VAR))/2,
      #vmin=-0.0002, vmax=0.0002,
      polar=[True, True], subplot=132, cmap='gnuplot',
      cbar=[False, 'vertical'],figsize=[15,8],
      label1='x', label2='y',title=r'Density [$\rho$]  $z=0$')
  T, newdims = tl.Tools(), (20,7)
  s0 = 15
  rr = np.linspace(D.x1[s0],D.x1[-1],len(D.x1[s0::]))
  tt = np.linspace(D.x2[0],D.x2[-1],len(D.x2))
  RRR, TTT = np.meshgrid(rr,tt)
  Xmesh, Ymesh = RRR*np.cos(TTT), RRR*np.sin(TTT)
  xcong = T.congrid(Xmesh,newdims,method='linear')
  ycong = T.congrid(Ymesh,newdims,method='linear')

  VX = D.vx1[s0::,:,cut].T*np.cos(TTT) - D.vx2[s0::,:,cut].T*np.sin(TTT)
  VY = D.vx1[s0::,:,cut].T*np.sin(TTT) + D.vx2[s0::,:,cut].T*np.cos(TTT)
  velxcong = T.congrid(VX,newdims,method='linear')
  velycong = T.congrid(VY,newdims,method='linear')
  plt.gca().quiver(xcong, ycong, velxcong, velycong,ec='w',fc='k',lw=0.8,\
      width=0.0030, headwidth=6, alpha=0.5)

  minute, second = divmod(D.SimTime,60)
  hour, minute   = divmod(minute,60)
  day, hour      = divmod(hour,24)
  plt.suptitle('Time: {} days  {} h  {} min {} s'.format(int(day),int(hour),
      int(minute),round(second,2)), fontsize=16, y=0.93, c='k')


  VAR = D.Bx3[:,0,:] - D0.Bx3[:,0,:] #D.rho[:,0,:] - D0.rho[:,0,:]
  #I.pldisplay(D, np.log10(D.rho[:,0,:]), x1=D.x1, x2=D.x3,
  I.pldisplay(D, VAR, x1=D.x1, x2=D.x3,
      polar=[False, False], subplot=133, cmap='gnuplot', 
      #vmin=-2.5, vmax=2.5,
      #vmin=-0.005, vmax=0.005,
      cbar=[False, 'vertical'],
      label1='r', label2='z',title=r'Density [$\rho$]  $\theta=0$')
  T, newdims = tl.Tools(), (20,10)
  idx0, idx1 = 10, 5
  Xmesh, Ymesh = np.meshgrid(D.x1.T,D.x3[idx0:-idx1].T)
  xcong = T.congrid(Xmesh,newdims,method='linear')
  ycong = T.congrid(Ymesh,newdims,method='linear')
  velxcong = T.congrid(D.vx1[:,0,idx0:-idx1].T,newdims,method='linear')
  velycong = T.congrid(D.vx3[:,0,idx0:-idx1].T,newdims,method='linear')
  plt.gca().quiver(xcong, ycong, velxcong, velycong,ec='w',fc='k',lw=0.8,\
      width=0.0030, headwidth=6, alpha=0.7)

  
  string = '%0.4d' % (num)  # To save CUSTOM FRAME With name pic.000#.png
  print("Image %04d"%num)
  plt.tight_layout()
  plt.savefig('IMAGES/pic.'+string+'.png',dpi=100)
  if plot:
    plt.show()
  else:
    plt.close('all')
#plt.savefig('amr_flowcyc.png')



number = 951
s = time.time()
print("Number of files is %d\n"%number)
num_cores = multiprocessing.cpu_count()
processed_list = Parallel(n_jobs=num_cores)(delayed(image)(i)
    for i in range(715,746+1,1))
print("Time spent {:.2f} s".format(time.time()-s))
