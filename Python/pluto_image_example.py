"""
  Example of usage and calling pluto_image.
"""

import matplotlib.pyplot as plt
import pluto_image as pi
import numpy as np
import os

# -- To save in folder
if os.path.exists('IMAGES'):
  print('Folder exists.')
else:
  os.system('mkdir IMAGES')

#- For a 3D image (for a 2D, dim=2):
xr, yr, array = pi.image(2, cmap='jet', aspect='auto', dim=3, dslice='13',n=55) 
xran = [np.min(xr),np.max(xr)]
yran = [np.min(yr),np.max(yr)]
XX0 = np.linspace(xran[0],xran[1],9)
XX1 = np.linspace(-5,5,9,dtype=int)
YY0 = np.linspace(yran[0],yran[1],9)
YY1 = np.linspace(100,0,9,dtype=int)
plt.xticks(XX0,XX1)
plt.yticks(YY0,YY1)
#plt.xlabel('x [km]')
#plt.ylabel('Depth [km]')

string = '%0.4d' % (2)  # To save custom frame with name pic.000#.png
#plt.savefig('IMAGES/pic.'+string+'.png')

plt.show()
