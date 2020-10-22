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
pi.image(2, cmap='jet', aspect='auto', dim=3, dslice='13',n=55) 
#plt.xticks(np.linspace(0,1e4,5), np.linspace(-5,5,5,dtype=int))
#plt.yticks(np.linspace(0,2e5,5), np.linspace(0,200,5,dtype=int))
plt.xlabel('x [km]')
plt.ylabel('Depth [km]')

string = '%0.4d' % (2)  # To save custom frame with name pic.000#.png
#plt.savefig('IMAGES/pic.'+string+'.png')

plt.show()
