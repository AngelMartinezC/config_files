"""
IMAGE CUBE
"""
import numpy as np
import matplotlib.pyplot as plt
import pluto_image as pi
import warnings   # ------ librería que esconde los warnings
import os
from mpl_toolkits.mplot3d import Axes3D
from rebin import rebin
warnings.filterwarnings("ignore")

os.system('clear')

def checkList(lst):
  #Function to check if all elements of 2D array are the same
  ele = lst[0,0]
  chk = True
  for y in range(lst.shape[0]):
    for x in range(lst.shape[1]):
      if ele != lst[y,x]:
        chk = False
        return False
        break;
  return True


unit = 1
def image(i,var='rho',vectorial=True,image=True,title='',vmin=None,vmax=None,
    dslice='12'):
  x, y, arr = pi.image(i,var=var,dim=3,dslice=dslice,n=0,figsize=(10,9),
      vmin=vmin, vmax=vmax, 
      cbarlabel=' ', xlabel=' ',cbar=False,
      unit=unit, image=image, vectorial=vectorial, aspect='equal',title=title)
  string = '%0.4d' %(i)
  #plt.savefig('IMAGES_2/pic.'+string+'.pdf')
  #plt.close("all")
  return arr


def cube(NUM=0,var='rho'):
  S0a = image(NUM,var=var,dslice='12',vectorial=True,image=False)
  S1a = image(NUM,var=var,dslice='23',vectorial=True,image=False)
  S2a = image(NUM,var=var,dslice='13',vectorial=True,image=False)
  
  Z0 = checkList(S0a)
  Z1 = checkList(S1a)
  Z2 = checkList(S2a)
  #print(Z0,Z1,Z2)
  
  var0 = np.zeros(S0a.shape)
  var1 = np.zeros(S1a.shape)
  var2 = np.zeros(S2a.shape)
  
  # if all elements are the same:
  # this lines set the same "colorbar" for every slice
  if Z0 and Z1 and Z2:
    maximo = [1.25]
    minimo = [0.75]
  elif Z0 and Z1:
    maximo = [np.max(S2a)]
    minimo = [np.min(S2a)]
  elif Z1 and Z2:
    maximo = [np.max(S0a)]
    minimo = [np.min(S0a)]
  elif Z2 and Z0:
    maximo = [np.max(S1a)]
    minimo = [np.min(S1a)]
  elif Z0:
    maximo = [np.max(S1a),np.max(S2a)]
    minimo = [np.min(S1a),np.min(S2a)]
  elif Z1:
    maximo = [np.max(S0a),np.max(S2a)]
    minimo = [np.min(S0a),np.min(S2a)]
  elif Z2:
    maximo = [np.max(S0a),np.max(S1a)]
    minimo = [np.min(S0a),np.min(S1a)]
  else:
    maximo = minimo = [0]
  
  #plt.subplot(131)
  #plt.imshow(S0a,origin='lower',vmin=0.94,vmax=1.06)
  #plt.colorbar()
  #plt.subplot(132)
  #plt.imshow(S1a,origin='lower',vmin=0.94,vmax=1.06)
  #plt.colorbar()
  #plt.subplot(133)
  #plt.imshow(S2a,origin='lower',vmin=0.94,vmax=1.06)
  #plt.colorbar()
  #plt.show()
  #exit()
  
  var0[0,0]   = np.max(maximo)
  var0[-1,-1] = np.min(minimo)-S0a[-1,-1]
  var0[-1,-1] = np.min(minimo)-S0a[-1,-1]-1000
  var1[0,0]   = np.max(maximo)
  var1[-1,-1] = np.min(minimo)-S1a[-1,-1]
  var2[0,0]   = np.max(maximo)
  var2[-1,-1] = np.min(minimo)-S2a[-1,-1]
  
  S0 = rebin(S0a+var0,S0a.shape)
  S1 = rebin(S1a+var1,S0a.shape)
  S2 = rebin(S2a+var2,S0a.shape)
  
  x0, y0 = S0.shape
  x, y = np.arange(0,x0,1), np.arange(0,y0,1)
  X, Y = np.meshgrid(x, y)
  
  f = fig = plt.figure(figsize=(9,8))
  f.subplots_adjust(bottom=0.07, left=0.02, right=0.96, top=0.95, \
    wspace=0.2, hspace=0.0)
  ax = fig.gca(projection='3d')
  cset = [[],[],[]]

  if var=='rho':
    vmin, vmax = 0.94, 1.04
  
  # this is the example that worked for you:
  cset[0] = ax.contourf(X, Y, S0, zdir='z', offset=y[-1],
      levels=np.linspace(np.min(minimo),np.max(maximo),100),
      vmin=vmin,vmax=vmax, cmap='viridis',alpha=0.8)
  
  # now, for the x-constant face, assign the contour to the x-plot-variable:
  cset[1] = ax.contourf(np.transpose(S1), Y, X, zdir='x', offset=x[-1],
      levels=np.linspace(np.min(minimo),np.max(maximo),100),
      vmin=vmin,vmax=vmax, cmap='viridis',alpha=0.8)
  
  # likewise, for the y-constant face, assign the contour to the y-plot-variable:
  cset[2] = ax.contourf(X, S2, Y, zdir='y', offset=0,
      levels=np.linspace(np.min(minimo),np.max(maximo),100),
      vmin=vmin,vmax=vmax, cmap='viridis',alpha=0.8)
  
  # setting 3D-axis-limits:    
  ax.set_xlabel('X1',labelpad=20)
  ax.set_ylabel('X2',labelpad=20)
  ax.set_zlabel('X3',labelpad=20)
  ax.set_xlim3d(x[0],x[-1])
  ax.set_ylim3d(y[0],y[-1])
  ax.set_zlim3d(0,x[-1])

  string = '%0.4d' %(NUM)
  plt.savefig('IMAGES_2/pic.'+string+'.png')

cube(15)
plt.show()

exit()

for i in range(0,101):
  print('image %d\n' %i)
  cube(i)
  plt.close('all')



exit()
#for i in range(0, 11):
#  f = plt.figure(figsize=(6,6))
#  f.subplots_adjust(bottom=0.00, left=0.00, right=1.0 , top=1.0 , \
#    wspace=0.2, hspace=0.0)
#  image(i,vmin=0.96,vmax=1.04)
#exit()

#"rho.0013.dbl
IM = 1
f = plt.figure(figsize=(19,6))
f.subplots_adjust(bottom=0.07, left=0.02, right=0.96, top=0.95, \
   wspace=0.2, hspace=0.0)

plt.subplot(131)
S=image(IM,var='rho',image=True)#,vmin=0.999, vmax=1.001)
#plt.imshow(S,origin='lower')
#plt.colorbar()

plt.subplot(132)
S=image(IM,var='vx1',image=True)
#plt.imshow(S,origin='lower')
#plt.colorbar()

plt.subplot(133)
S=image(IM,var='vx2',image=True)
#plt.imshow(S,origin='lower')
#plt.colorbar()

plt.show()

exit()
#"dbl_files/rho.0020.dbl"
for i in range(0, 13):
  image(i)
