import numpy as np
from sunpy.map import Map
import matplotlib.pyplot as plt
from scipy.fft import fft2, ifft2, ifftshift
from numpy import unravel_index
import glob
import os

# -- To save in folder
if os.path.exists('IMAGES2'):
    print('Folder exists.')
else:
  os.system('mkdir IMAGES2')

# -- Shift data
def shifted(img0, img1, hx=1641, hy=2478,width=1000):
  # -- Get data
  def correlation(img0, img1, hx, hy,width):
    a = img0.data
    b = img1.data
    
    # -- Correlation box
    bx = width
    by = width
    
    # -- Center of figure correlation
    hx = hx
    hy = hy
    
    # -- Get half box
    mx = int(bx/2)
    my = int(by/2)
    
    ssz = a.shape
    nx = ssz[1]
    ny = ssz[0]
    
    # -- Cut in FoV and Fourier transform
    na  = a[hy-my:hy+my-1,hx-mx:hx+mx-1]
    nb  = b[hy-my:hy+my-1,hx-mx:hx+mx-1]
    fna = fft2(na*1)
    fnb = fft2(nb*1)
    
    # -- Cross correlation function and normalize
    ccf0 = np.real(ifft2(fna*np.conj(fnb)))
    ccf  = ccf0/np.max(ccf0)
    newx, newy = ccf.shape
    
    # -- Create arrays to be shifted with the information of the correlation
    CCF = np.zeros((newx,newy))
    CCF0 = np.zeros((newx,newy))
    CCF0[:,0:mx]  = ccf[:,-mx::]
    CCF0[:,-mx::] = ccf[:,0:mx]
    CCF[0:my,:]  = CCF0[-my::,:]
    CCF[-my::,:] = CCF0[0:my,:]
    #plt.subplot(121)
    #plt.imshow(ccf)
    #plt.subplot(122)
    #plt.imshow(CCF)
    #plt.show()
  
    # -- Find maximum index
    s01, s00 = unravel_index(CCF.argmax(), CCF.shape)
    sx = s00 - mx
    sy = s01 - my
    return sx, sy

  sx, sy = correlation(a0,b0,hx,hy,width)
  mx = int(width/2)
  my = int(width/2)
  
  na  = img0.data[hy-my:hy+my-1,hx-mx:hx+mx-1]
  nb  = img1.data[hy-my:hy+my-1,hx-mx:hx+mx-1]
  nb2 = img1.data[hy-my-sy:hy+my-1-sy,hx-mx-sx:hx+mx-1-sx]
  return na, nb2

LIST = sorted(glob.glob("hmi*"))

a0 = Map(LIST[0])

hx = 1641
hy = 2478
width = 1000
mx = my = int(width/2)

plt.figure(figsize=(9,7))
plt.subplot(111)
plt.imshow(a0.data[hy-my:hy+my-1,hx-mx:hx+mx-1])#,vmin=-5000,vmax=5000)
plt.colorbar()
#plt.show()
plt.savefig("IMAGES2/pic.0000.png")
plt.close("all")

for i in range(1,len(LIST)):
  b0 = Map(LIST[i])
  na, sb = shifted(a0, b0, hx, hy, width)

  plt.figure(figsize=(9,7))
  plt.subplot(111)
  plt.imshow(sb)#,vmin=-5000,vmax=5000)
  plt.colorbar()
  string = '%0.4d' % (i)  # To save CUSTOM FRAME With name pic.000#.png
  #plt.tight_layout()
  plt.savefig('IMAGES2/pic.'+string+'.png')
  #plt.show()
  plt.close("all")
  




