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

#- For a 3D image (for a 2D, dim=2):
counter = 0
#for i in range(0,10):
def image(i):
  xr, yr, array = pi.image(i, cmap='jet', aspect='equal', figsize=figsize) 
  xran = [np.min(xr),np.max(xr)]
  yran = [np.min(yr),np.max(yr)]
  XX0 = np.linspace(xran[0],xran[1],9)
  XX1 = np.linspace(0,10,9,dtype=int)
  YY0 = np.linspace(yran[0],yran[1],9)
  YY1 = np.linspace(0,10,9,dtype=int)
  plt.xticks(XX0,XX1)
  plt.yticks(YY0,YY1)
  #plt.xlabel('x [km]')
  #plt.ylabel('Depth [km]')
  
  print("image {}\n".format(i))
  string = '%0.4d' % (i)  # To save custom frame with name pic.000#.png
  plt.savefig('IMAGES/pic.'+string+'.png')
  #counter += 1
  plt.close("all")

s = time.time()
#for i in range(195,226):
#  image(i)
#print("Time spent {:.2f} s".format(time.time()-s))
#exit()

num_cores = multiprocessing.cpu_count()
processed_list = Parallel(n_jobs=num_cores)(delayed(image)(i) 
    for i in range(215,350))
print("Time spent {:.2f} s".format(time.time()-s))

#plt.show()
print("END")
