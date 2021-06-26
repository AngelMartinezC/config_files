"""
  Plot an HMI map into its Carrington coordinates
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from sunpy.map import Map
from astropy import units as u
import glob
import multiprocessing
from joblib import Parallel, delayed
import os
import time


x0 = 2545
y0 = 1948
width = 500


# -- To save in folder
if os.path.exists('IMAGES'):
    print('Folder exists.')
else:
  os.system('mkdir IMAGES')


# -- Constants
RSun = 6.9634e8   # Sun radius im m

# -- Convert pix to arcsec
def pix2xy(x0,y0,X0,Y0,RES):
  t1 = (X0 - x0)*RES
  t2 = (Y0 - y0)*RES
  return t1, t2

# -- Cut map
def cut_map(data,t1,t2,width=400):
  t1 *= u.arcsec
  t2 *= u.arcsec
  width *= u.arcsec
  top_right = SkyCoord(t1+width/2, t2+width/2, 
      frame=data.coordinate_frame)
  bottom_left = SkyCoord(t1-width/2, t2-width/2, 
      frame=data.coordinate_frame)
  submap = data.submap(bottom_left, top_right=top_right)
  return submap

# -- Get CRLN CRLT from data
def xy2lonlat(x0,y0,X0,Y0,RES,B0,L0,d):
  # -- Convert pix to arcsec
  t1 = (X0 - x0)*RES
  t2 = (Y0 - y0)*RES
  # -- Convert helioprojective to heliocentric
  x = (d-RSun)*np.pi/180*t1/3600
  y = (d-RSun)*np.pi/180*t2/3600
  # -- Convert to Carrington coordinates
  r = RSun
  z = np.sqrt(RSun**2-x**2-y**2)
  CRLT = 180/np.pi*np.arcsin((y*np.cos(B0)+z*np.sin(B0))/r)
  CRLN = 180/np.pi*np.arctan(x/(z*np.cos(B0)-y*np.sin(B0))) + L0
  #print(t1,t2, CRLT, CRLN)
  return CRLT, CRLN

# -- Get variables to conversion (Carrington to Stonyhurst
def lonlat2xy(CRLT, CRLN, L0, B0, d):
  HGLT  = CRLT*np.pi/180
  HGLN2 = (CRLN-L0)*np.pi/180
  B0 *= np.pi/180
  # -- Convert coordinates
  r = RSun
  x = r*np.cos(HGLT)*np.sin(HGLN2)
  y = r*(np.sin(HGLT)*np.cos(B0) - np.cos(HGLT)*np.cos(HGLN2)*np.sin(B0))
  z = r*(np.sin(HGLT)*np.sin(B0) + np.cos(HGLT)*np.cos(HGLN2)*np.cos(B0))
  # -- Find arcsec from conversion
  t1 = x*3600/(d-RSun)*180/np.pi
  t2 = y*3600/(d-RSun)*180/np.pi
  #print("L0 is: ",L0)
  #print(t1,t2)
  return t1, t2


# -- main function

LIST = np.array(sorted(glob.glob("*.fits")))

data = Map(LIST[0])
header = data.fits_header
X0 = header["CRPIX1"]
Y0 = header["CRPIX2"]
B0 = header["CRLT_OBS"]*np.pi/180
L0 = header["CRLN_OBS"]
d  = header["DSUN_OBS"]
RES= header["CDELT1"]

CRLT, CRLN = xy2lonlat(x0,y0,X0,Y0,RES,B0,L0,d)
t10, t20 = pix2xy(x0,y0,X0,Y0,RES)
#t1, t2 = lonlat2xy(CRLT,CRLN,L0, B0)

#LIST = ["hmi.ic_45s.20110215_013945_TAI.continuum.fits","hmi.ic_45s.20110215_020945_TAI.continuum.fits"]
#for i in len(LIST):
def image(i,count=0):
  frame = LIST[i]
  data = Map(frame)
  header = data.fits_header
  X0 = header["CRPIX1"]
  Y0 = header["CRPIX2"]
  B0 = header["CRLT_OBS"]#*np.pi/180
  L0 = header["CRLN_OBS"]
  d  = header["DSUN_OBS"]
  RES= header["CDELT1"]
  if count == 0:
    frame0 = LIST[0]
    data0 = Map(frame0)
    header = data0.fits_header
    X0 = header["CRPIX1"]
    Y0 = header["CRPIX2"]
    RES= header["CDELT1"]
    t1, t2 = pix2xy(x0,y0,X0,Y0,RES)
  else:
    t1, t2 = lonlat2xy(CRLT,CRLN, L0, B0, d)
    #print("Here")
  submap = cut_map(data,t1,t2,width=width)
  f = plt.figure(figsize=(6,6))
  f.subplots_adjust(bottom=0.05, left=0.05, right=0.93, top=0.93, \
      wspace=0.0, hspace=0.0)
  ax = plt.subplot(111,projection=submap)
  plt.imshow(submap.data)
  submap.plot()
  ax.tick_params(labelleft=False, axis='y', grid_linewidth=1)
  ax.tick_params(labelbottom=False, axis='x', grid_linewidth=1)
  #plt.title(" ")
  plt.xlabel(" ")
  plt.ylabel(" ")
  #submap.draw_grid()
  #print("Image %04d"%i)
  print("Image %04d  t1: %.5f  t2: %.5f  CRLT %.4f  LEN1 %d  LEN2 %d"
      %(i,t1,t2,CRLT,len(submap.data[:,0]),len(submap.data[0])))
  string = '%0.4d' % (i)  # To save CUSTOM FRAME With name pic.000#.png
  #plt.tight_layout()
  plt.savefig('IMAGES/pic.'+string+'.png')
  #plt.show()
  plt.close('all')

#image(0)
#exit()


s = time.time()
number = int(os.popen('ls hmi* | grep .fits | wc -l').read())
print("Number of files is %d"%number)

num_cores = multiprocessing.cpu_count()
processed_list = Parallel(n_jobs=num_cores)(delayed(image)(i,i)
    for i in range(0,number,1))
print("Time spent {:.2f} s".format(time.time()-s))




exit()
t1, t2 = pix2xy(x0,y0)
submap = cut_map(data,t1,t2)
submap.plot()
plt.show()

exit()

INT  = Map("INT.fits")

plt.figure(figsize=(11,7))
plt.subplot(121,projection=data)
data.plot()
data.draw_grid(grid_spacing=10*u.deg)
plt.grid()

plt.subplot(122,projection=INT)
INT.plot()
INT.draw_grid()
plt.show()



