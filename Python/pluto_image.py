"""
  Python file to be called for plotting .dbl files in 2 and 3 dimensions.
  
  Add:
    bash> export PYTHONPATH=$PYTHONPATH:$HOME/config_files/Python
  at the end of .bashrc file.

"""

import numpy as np
import matplotlib.pyplot as plt
import pyPLUTO as pp
import os
import sys 
plt.rcParams.update({'font.size': 18})

wdir = os.popen('echo `pwd`/').read()
wdir = wdir.replace('\n','')

#p = sys.argv[1]
def image(frame, var='density', aspect='auto', xlabel='x', ylabel='y', \
    cbar=True, title=r'Density $\rho$', vmin=None, vmax=None, Mm=True, \
    Mmx=True, Mmy=True, cmap='jet', figsize=(8,10), labelpad=10.0, \
    cbarlabel=r'Density ($\times$10$^{10}$) [gr cm$^{-3}$]', pad=0.05, 
    dim=2, n=0, dslice='12', unit=None, diff=None, step=1, **kwargs):

  """
    If the file has 3 dimensions, select a dslice (transversal o profile 
    image)
  """
  
  if diff is None:
    D = pp.pload(frame,w_dir=wdir)
    if (var == 'rho') or (var == 'density'):
      variable = D.rho
    elif (var == 'pressure') or (var == 'prs'):
      variable = D.prs
    elif (var == 'velocity1') or (var == 'vx1'):
      variable = D.vx1
    elif (var == 'velocity2') or (var == 'vx2'):
      variable = D.vx2
    elif (var == 'velocity3') or (var == 'vx3'):
      variable = D.vx3
    elif (var == 'mag_b1') or (var == 'Bx1'):
      variable = D.Bx1
    elif (var == 'mag_b2') or (var == 'Bx2'):
      variable = D.Bx2
    elif (var == 'mag_b3') or (var == 'Bx3'):
      variable = D.Bx3
    elif (var == 'divV'):# or (var == 'Bx3'):
      variable = D.divV
    elif (var == 'divB'):# or (var == 'Bx3'):
      variable = D.divB
    elif (var == 'Temp') or (var == 'temperature'):
      variable = D.Temp
    elif (var == 'PTOT') or (var == 'P_TOT'):
      variable = D.PTOT
  else:
    D2= pp.pload(frame+step,w_dir=wdir)
    if (var == 'rho') or (var == 'density'):
      variable2= D2.rho
    elif (var == 'pressure') or (var == 'prs'):
      variable2= D2.prs
    elif (var == 'velocity1') or (var == 'vx1'):
      variable2= D2.vx1
    elif (var == 'velocity2') or (var == 'vx2'):
      variable2= D2.vx2
    elif (var == 'velocity3') or (var == 'vx3'):
      variable2= D2.vx3
    elif (var == 'mag_b1') or (var == 'Bx1'):
      variable2= D2.Bx1
    elif (var == 'mag_b2') or (var == 'Bx2'):
      variable2= D2.Bx2
    elif (var == 'mag_b3') or (var == 'Bx3'):
      variable2= D2.Bx3
    elif (var == 'divV'):# or (var == 'Bx3'):
      variable2= D2.divV
    elif (var == 'divB'):# or (var == 'Bx3'):
      variable2= D2.divB
    elif (var == 'Temp') or (var == 'temperature'):
      variable2= D2.Temp
    elif (var == 'PTOT') or (var == 'P_TOT'):
      variable2= D2.PTOT

  if dim == 3:
    if dslice == '12':
      variable = variable[:,:,n]
      xran = D.x1
      yran = D.x2
    elif dslice == '13':
      variable = variable[:,n,:]
      xran = D.x1
      yran = D.x3
    elif dslice == '23':
      variable = variable[n,:,:]
      xran = D.x2
      yran = D.x3
  else:
      xran = D.x1
      yran = D.x2

  if unit is None: unit = 1
    
  if vmin is None: vmin = np.min(variable)*unit
  if vmax is None: vmax = np.max(variable)*unit

  I = pp.Image()
  if diff is None:
    I.pldisplay(D, variable, x1=xran,x2=yran,label1=xlabel,label2=ylabel,\
      title=title, cbar=(cbar,'vertical'), vmin=vmin, vmax=vmax, pad=pad, \
      aspect=aspect, figsize=figsize, cmap=cmap, cbarlabel=cbarlabel, \
      labelpad=labelpad, unit=unit, **kwargs)
  else:
    I.pldisplay(D, variable2-variable, x1=xran, x2=yran, label1=xlabel, \
      label2=ylabel, title=title, cbar=(cbar,'vertical'), vmin=vmin, \
      vmax=vmax, pad=pad, aspect=aspect, figsize=figsize, cmap=cmap, \
      cbarlabel=cbarlabel, labelpad=labelpad, unit=unit, **kwargs)
  dx = np.ones(10)
  dy = np.ones(10)
  #I.field_line(D.vx1,D.vx2,D.x1,D.x2,D.dx1,D.dx2,dx,dy)
  #I.myfieldlines(D,np.linspace(-5,5,10),np.linspace(0,10,0))

  #x0arr = np.linspace(0.0,np.max(xran)*1e3,20) 
  #y0arr = np.linspace(0.0,np.max(yran)*1e3,20) 
  #I.myfieldlines(D,x0arr,y0arr,colors='k',ls='--',lw=1.0) 
  #print(I)

  return xran, yran

