"""
  Example of usage and calling pluto_image.
"""

import matplotlib.pyplot as plt
import pluto_image as pi
import numpy as np
import os
import multiprocessing
from joblib import Parallel, delayed
import time

# -- To save in folder
if os.path.exists('IMAGES'):
  print('Folder exists.')
else:
  os.system('mkdir IMAGES')

figsize = (11,10)
unit = 1e-2
def image(i):
  xr, yr, array = pi.image(i, cmap='jet', aspect='equal', figsize=figsize,\
      cbarlabel=r'$\Delta$Density ($\times 10^{-4}$) [g cm$^{-3}$]', \
      vmin=0, vmax=7.0, \
      #vmin=-0.004, vmax=0.004, diff=True, step=1, \
      unit=unit, wdir='dbl_files/',var='density')
  xran = [np.min(xr),np.max(xr)]
  yran = [np.min(yr),np.max(yr)]
  XX0 = np.linspace(xran[0],xran[1],9)
  XX1 = np.linspace(-5,5,9,dtype=int)
  YY0 = np.linspace(yran[0],yran[1],9)
  YY1 = np.linspace(10,0,9,dtype=int)
  plt.xticks(XX0,XX1)
  plt.yticks(YY0,YY1)
  plt.xlabel(r'x [Mm]')
  plt.ylabel(r'z [Mm]')
  minute = i/2
  days = int(minute // 1440)
  hours = int((minute - (days*1440)) // 60)
  minutes = int((minute - (days*1440) - (hours*60)) // 1)
  plt.title('Time: {} days  {} h  {} min'.format(days,hours,minutes), 
      fontsize=18)
  
  print("image {}\n".format(i))
  string = '%0.4d' % (i)  # To save custom frame with name pic.000#.png
  plt.savefig('IMAGES/pic.'+string+'.png')
  plt.close("all")

s = time.time()

#num_cores = multiprocessing.cpu_count()
num_cores = 2
processed_list = Parallel(n_jobs=num_cores)(delayed(image)(i) 
    for i in range(0,401,1))

print("Time spent {:.2f} s".format(time.time()-s))
print("END")
