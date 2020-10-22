import numpy as np
import matplotlib.pyplot as plt
import pyPLUTO as pp
import os
import sys 
plt.rcParams.update({'font.size': 18})

wdir = os.popen('echo `pwd`/').read()
wdir = wdir.replace('\n','')
print(wdir)

#p = sys.argv[1]
def image(frame, var='density', aspect='auto', xlabel='x', ylabel='y', \
    cbar=True, title=r'Density $\rho$', vmin=None, vmax=None, Mm=True, \
    Mmx=True, Mmy=True, cmap='viridis', figsize=(8,10), labelpad=10.0, \
    cbarlabel=r'Density ($\times$10$^{10}$) [gr cm$^{-3}$]', pad=0.05, 
    dim=2, n=0, dslice='12', **kwargs):
  
  D = pp.pload(frame,w_dir=wdir)
  if (var == 'rho') or (var == 'density'):
    variable = D.rho
  elif (var == 'pressure') or (var == 'prs'):
    variable = D.prs

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
    
  if vmin == None: vmin = np.min(variable)
  if vmax == None: vmax = np.max(variable)

  I = pp.Image()
  I.pldisplay(D, variable, x1=xran, x2=yran, label1=xlabel, label2=ylabel,\
      title=title, cbar=(cbar,'vertical'), vmin=vmin, vmax=vmax, pad=pad, \
      aspect=aspect, figsize=figsize, cmap=cmap, cbarlabel=cbarlabel, \
      labelpad=labelpad, **kwargs)
  
  #x0arr = np.linspace(0.0,np.max(xran)*1e3,20) 
  #y0arr = np.linspace(0.0,np.max(yran)*1e3,20) 
  #I.myfieldlines(D,x0arr,y0arr,colors='k',ls='--',lw=1.0) 
  #print(I)

