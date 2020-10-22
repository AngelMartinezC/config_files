#-*- coding:utf-8 -*-

"""
	Program to make the time-distance plot of a given array.

	It is necessary to have the location x0 and y0 of the suspected sunquake,
	as well as two angles (in degrees) in which it is located.
	If a .sav file exists, it is possible to make the distribution of times
	and distance automatically, otherwise it will plot only the time 
	(in every 45 seconds) and the distance in pixels
	
	Package polarTransform is needed. This package can be downloaded from 
		https://github.com/addisonElliott/polarTransform
		-----
			pip install PolarTransform
		-----
	Installing process is detailed on the page
	
	Date created: July 10 2018
	Date last modified: August 19 2018

	
"""


from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sp
import time
import polarTransform
from datetime import datetime
import matplotlib.patches as mpatches




class cubo:
	
	def getdata(self,array,savfile):
		
		if isinstance(array,(str)):
			cube = fits.getdata(array,0)
			if savfile:
				head = sp.readsav(savfile)
				header = head.aiatimes
				return cube, header
			else:
				return cube, None
		elif isinstance(array,(int,float)):
			return "Not a valid array"
		else: 
			return NotImplemented




class td:

	rsun_m = 695.5e6
	rsun_pix = 1884 # Value for HMI
	
	def __init__(self,array, x0=148, y0=140, theta0=192, theta1=260, \
		savfile=False, path=None, radius=None, radius0=0, time1=False, \
		time0=False, cartesian=None, xticks=5, yticks=6):
	
		self.array = array
		self.time1 = time1
		self.time0 = time0
		self.x0 = x0
		self.y0 = y0
		self.theta0 = theta0
		self.theta1 = theta1
		self.radius = radius
		self.radius0 = radius0
		self.cartesian = cartesian
		self.xticks = xticks
		self.yticks = yticks
		self.sav = savfile
		self.path = path


	################################################################################################
	
	
	# Obtain the information of the header as a readable file
	# If the sav file is not given, the output labels of the tdplot graph
	# will be time in seconds and distance in pixels
	
	def values(self):
		
		if self.sav:
			headi = cubo().getdata(self.array,self.sav)
			head = headi[1]
			frame = int(len(head)/2)
			datef = str(head[frame][0:10])
			date = datef[2:12]
			d = datetime.strptime(date, '%Y-%m-%d')
			data = d.strftime(r"%B %d$^{\mathrm{th}}$, %Y")
			value = []
			for i in range(len(head)):
				value.append(head[i][11:16].decode("utf-8")) 
			return data, value
		else: 
			pass

	################################################################################################

	
	# Function wich makes the time-distance plot
	
	def tdplot(self, colorbar=False, plot=None, save=None, noprint=None, \
		notime=None, nolabel=False,title=' ',notitle=False,slide=False,magneto=False,intensity=False,cmap=None,**kwargs):
		
		self.colorbar = colorbar
		self.plot = plot
		self.save = save
		self.noprint = noprint
		self.notime = notime
		self.nolabel = nolabel
		self.title = title
		self.notitle = notitle
		self.slide = slide
		self.magneto = magneto
		self.intensity = intensity
		self.cmap = cmap

		# Call the array
		cube = cubo().getdata(self.array,self.sav)
		s = time.time()
		
		# Convert the angles to radians
		t0 = self.theta0*np.pi/180
		t1 = self.theta1*np.pi/180
		
		# To ensure the minimum radius in the plot
		# This is to set finalRadius in converting to polar coordinates
		if self.radius == None:
			mx = cube[0][0].shape[0]/2
			my = cube[0][0].shape[1]/2
			if self.x0<=mx and self.y0<=my:
				rad = min(self.x0,self.y0)
			elif self.y0>my and self.x0<=mx:
				rad = min(self.x0,my*2 -self.y0)
			elif self.x0>mx and self.y0<=my:
				rad = min(self.y0,mx*2 -self.x0)
			elif self.x0>mx and self.y0>my:
				rad = min(mx*2 -self.x0,my*2 -self.y0)
		else: 
			rad = self.radius
		
		if self.noprint == None:
			print("\n Radius in pixels is {}".format(rad))
		else: 
			pass

		# If put, the final image will be up the selected one, 
		#  otherwise it uses the entire range-1
		if self.time0:
			image0 = self.time0
		else: image0 = 0

		if self.time1:
			image1 = self.time1
		else: image1 = len(cube[0])-1
		
		# td graph
		image_td = []
		for ij in range(image0,image1):
			polarImage, ptSettings = polarTransform.convertToPolarImage(cube[0][ij],\
								center=[self.x0, self.y0], initialRadius = self.radius0, \
								finalRadius = rad, initialAngle = t0, finalAngle = t1)
			slices = []
			for i in range(len(polarImage[0,:])):
				slices.append(np.mean(polarImage[:,i]))
			image_td.append(slices)
		e = time.time()
		
		"""-----------------------------------------------------------------------------
		   Variables for the graph. The final graph wil be shown if is
		   selected the keyword "graph", and will be saved if the keyword
		   "save" is selected
		"""
		if slide:
			pass
		else:
			plt.rcParams.update({'font.size': 12})
			plt.figure(figsize=(8,6))
			plt.subplots_adjust(right=0.99,bottom=0.1,top=0.95,left=0.11)
		
		if self.sav:
			
			# Call the function which has information of the flare
			main = self.values()

			# Labels for x range
			xdiv = self.xticks
			xrangex = np.linspace(0,len(image_td[0]),xdiv)
			x2Mm1 = rad*695.51/1884 # Rsun[Mm]=695.5 & Rsun[pixels]=1884
			x2Mm0 = self.radius0*695.51/1884
			xlabel = np.linspace(x2Mm0,x2Mm1,xdiv,dtype=int)
			xlabel = list(map(str,xlabel))

			# Labels for y range
			ydiv = self.yticks
			yrange = np.linspace(0,int(image1-image0),ydiv)
			
			ysec = np.array(yrange,dtype=int)
			ylabel = [main[1][i+image0] for i in ysec]
			
			if self.notitle:
				plt.title(self.title)
			else: 
				plt.title(self.title+main[0])
			plt.xticks(xrangex,xlabel)
			plt.yticks(yrange,ylabel)
			plt.xlabel("Distance (Mm)")
			plt.ylabel("Time [UT]")
			# Image of the variables with all labels defined
	
		else:
			if self.nolabel:
				pass
			else:
				plt.xlabel("Distance [pixels]")
				plt.ylabel("Time [45 seconds]")
		
		if self.cmap is None or self.cmap == False:
			plt.imshow(image_td,origin='lower',cmap='Greys_r',interpolation = 'spline36',aspect='auto',**kwargs)
		else:
			plt.imshow(image_td,origin='lower',interpolation = 'spline36',aspect='auto',cmap=self.cmap,**kwargs)
		
		
		if self.colorbar:
			if self.intensity:
				plt.colorbar(label=r'$\Delta$ Counts',pad=0.03,aspect=18)
			elif self.magneto:
				plt.colorbar(label='Magnetic field (Gauss)',pad=0.03,aspect=18)
			else:
				plt.colorbar(label='Velocity (m/s)',pad=0.03,aspect=18)
		else: 
			pass
		
		"""
		   Final of the graph 
		   -----------------------------------------------------------------------------
		"""
		
		# Save the graph
		if self.save:
			if self.path == None:
				plt.savefig("Sunquake_x:"+str(self.x0)+"_y:"+str(self.y0)+\
					"_t0:"+ str(self.theta0)+"_t1:"+str(self.theta1)+"_r:"+\
					str(rad)+".png",format="png",dpi=450)
				plt.savefig("Sunquake_x:"+str(self.x0)+"_y:"+str(self.y0)+\
					"_t0:"+str(self.theta0)+"_t1:"+str(self.theta1)+"_r:"+\
					str(rad)+".pdf",format="pdf",dpi=450)
			else:	
				plt.savefig(self.path+"Sunquake_x:"+str(self.x0)+"_y:"+\
					str(self.y0)+"_t0:"+str(self.theta0)+"_t1:"+\
					str(self.theta1)+"_r:"+str(rad)+".png",format="png")
				plt.savefig(self.path+"Sunquake_x:"+str(self.x0)+"_y:"+\
					str(self.y0)+"_t0:"+str(self.theta0)+"_t1:"+\
					str(self.theta1)+"_r:"+str(rad)+".pdf",format="pdf")
		else:
			pass
		
		# Show the time spent making the graph
		if self.notime == None:
			print(" The time spent is {:0.4f} s".format(e-s))
		else: 
			pass
		
		# Show the graph
		if self.plot:
			plt.show()
		else: 
			pass

		# Return a cartesian version of the set 
		if self.cartesian == None:
			pass
		else:
			cartesianImage = ptSettings.convertToCartesianImage(polarImage)
			return cartesianImage
		
		
		# At the end the function returns teh array of teh time-distance plot

		return image_td
	
	
	
	################################################################################################
	
	"""
		The next functions ar not required to show the td plot.
		
		Functions:
		----------
		
		
		-- test: --
		
			Given the number of columns and rows, it plot as many time-distance plots as are in 
			columns*rows. 
			
			
		-- slider: --
		
			Useful function to slide between pixels and angles without need to compile by hand new 
			values of pixels and angles.
			
			
		-- cbar_slider: --
				
			To handle color contrast (min and max of colorbar)
		
	"""

	################################################################################################


	# Makes a test of many pixels around (5 columns and 4 rows)
	# Reference to IDL procedure (Martinez-Oliveros)
	
	def test(self,columns=5, rows=4):
		
		s = time.time()
		self.x = columns
		self.y = rows
		fig = plt.figure(figsize=(38, 38))
		n=0
		print("\n")
		for j in range(self.y,0,-1):
			for i in range(0,self.x):
				n+=1
				plt.subplot(self.y,self.x,n)
				print("   Plotting image {} of {}".format(n,self.x*self.y))
				Y0 = self.y0-int(self.y/2)+j
				X0 = self.x0-int(self.x/2)+i
				rad = self.radius
				
				image = td(self.array, x0=X0, y0=Y0, theta0=self.theta0, theta1=self.theta1,time0 = self.time0,time1 = self.time1, 
				radius=rad, savfile=self.sav, path=self.path,\
				radius0=self.radius0)
				
				final = image.tdplot(plot=False, colorbar=False,noprint=False,notime=False)
				plt.title(" ")
				plt.xlabel(" ")
				plt.ylabel(" ")
				plt.imshow(final,cmap='Greys_r',origin='lower',interpolation='spline36')
				red_patch = mpatches.Patch(color='red', label="x="+str(X0)+"\ny="+str(Y0))
				plt.legend(handles=[red_patch],bbox_to_anchor=(1, .33))
		e = time.time()
		print("\n Time spent in the test {:0.4f} s".format(e-s))	
		plt.show()		
		
	################################################################################################
	
	
	def slider(self):
		
		from matplotlib.widgets import Slider, Button, RadioButtons

		# Image configuration and first call
		plt.figure(figsize=(10,8))
		fig, ax = plt.subplots()
		fig.subplots_adjust(bottom=0.15)
		
		def data(x0,y0,t0,t1,radius0):
			image = td(self.array, x0=int(x0), y0=int(y0), theta0=int(t0), theta1=int(t1),radius0=int(radius0), radius=self.radius, savfile=self.sav, path=self.path,time0 = self.time0,time1 = self.time1)
			im = image.tdplot(plot=False,colorbar=False,slide=True)
			return im
		
		im = data(self.x0,self.y0,self.theta0,self.theta1,self.radius0)
		image = ax.imshow(im, cmap='Greys_r', interpolation='spline36',origin='lower')

		# Remove actual image to not be shown
		#plt.show(block=False)
		#plt.pause(0.0001)
		plt.close('all')

		# Set up sliders
		axt1 = plt.axes([0.1, 0.06, 0.34, 0.04], facecolor='lightcyan')
		axt0 = plt.axes([0.1, 0.01, 0.34, 0.04], facecolor='pink')
		axy0 = plt.axes([0.55, 0.06, 0.35, 0.04], facecolor='lightcyan')
		axx0 = plt.axes([0.55, 0.01, 0.35, 0.04], facecolor='pink')
		axrd = plt.axes([0.1,0.86,0.34,0.04], facecolor='pink')
		
		slider_r0 = Slider(axrd,r'radio',0,50,valinit = self.radius0)
		slider_t1 = Slider(axt1,r'$\theta _1(^\circ)$',0, 359, valinit = self.theta1)
		slider_t0 = Slider(axt0,r'$\theta _0(^\circ)$',0, 359, valinit = self.theta0)
		slider_x0 = Slider(axy0,r'x',0,400, valinit = self.x0) 
		slider_y0 = Slider(axx0,r'y',0,400, valinit = self.y0)

		plt.axes([0.25, 0.25, 0.6, 0.6])

		# Update the image 
		def update(val):
			theta0 = slider_t0.val
			theta1 = slider_t1.val
			x0 = slider_x0.val
			y0 = slider_y0.val
			radius0 = slider_r0.val
			image.set_data(data(int(x0),int(y0),int(theta0),int(theta1),int(radius0)))
			fig.canvas.draw()
			ax.clear()
			plt.draw()
		
		# Change sliders
		slider_r0.on_changed(update)
		slider_t0.on_changed(update)
		slider_t1.on_changed(update)
		slider_x0.on_changed(update)
		slider_y0.on_changed(update)

		plt.show()
	
	################################################################################################

	
	# Visualization of the time distance plot cahnging (not necessary)
	def cbar_slider(self,**kwargs):
		
		from matplotlib.widgets import Slider, Button, RadioButtons
		
		fig = plt.figure()
		ax = fig.add_subplot(111)
		fig.subplots_adjust(left=0.15, bottom=0.25)
		
		image = td(self.array, x0=self.x0, y0=self.y0, theta0=self.theta0, theta1=self.theta1,time0=self.time0,time1=self.time1,\
		radius = self.radius, savfile=self.sav, path=self.path, radius0=self.radius0)
		im = np.array(image.tdplot(plot=False, colorbar=False, noprint=False,slide=True)) 
		
		im1 = ax.imshow(im,cmap='Greys_r',origin='lower',interpolation='spline36',aspect='auto',**kwargs)
		fig.colorbar(im1)
		axcolor = 'gainsboro'
		axmin = fig.add_axes([0.1, 0.05, 0.65, 0.03], facecolor=axcolor)
		axmax  = fig.add_axes([0.1, 0.1, 0.65, 0.03], facecolor=axcolor)
		
		min0 = im.min()
		max0 = im.max()
		smin = Slider(axmin, 'Min', -1*abs(im.min()*1.5), abs(im.max()*1.5), valinit=min0, color='dimgray')
		smax = Slider(axmax, 'Max', -1*abs(im.min()*1.5), abs(im.max()*1.5), valinit=max0, color='dimgray')
		
		def update(val):
			im1.set_clim([smin.val,smax.val])
			fig.canvas.draw()
		smin.on_changed(update)
		smax.on_changed(update)
		
		plt.show()
	
	
		
	################################################################################################
	
	### ----- end of td class
	
	################################################################################################


if __name__ == "__main__":

	path = "/home/angel/IDLWorkspace/Flare/F29/HMIDoppler/Work/"
	flare = "/home/angel/IDLWorkspace/Flare/F29/HMIDoppler/Work/HMIDoppler_29_diff.fits"
	savfile = '/home/angel/IDLWorkspace/Flare/F29/HMIDoppler/Work/HMIDoppler_coord_29.sav'



	#  Initial parameters


	x0 = 149
	y0 = 149
	x0 = 148-4
	y0 = 169+2
	theta0 = 130
	theta1 = 150
	theta0 = 140
	theta1 = 160
	radius = 220
	radius = 120

	time0 = 180
	time1 = 390
	time0 = 200
	time1 = 300


	data = td(flare, x0=x0, y0=y0, theta0=theta0, theta1=theta1,\
	radius=radius,radius0=0, savfile=savfile, path=path, time0=time0,time1=time1,xticks=7)

	plot = data.tdplot(save=False,colorbar=True,nolabel=False,plot=True,notitle=True,title='SOL2014-02-04T03:00-M1.1')#,vmin=0)


	exit()


	# Paths where the files are located
	plt.figure(figsize=(10,8))
	#savfile = "/home/angel/IDLWorkspace/Flare/Cycle24/F1/HMIDoppler/Work/HMIDoppler_coord_1.sav"
	flare = "/home/angel/IDLWorkspace/Flare/F29/HMIDoppler/Work/HMIDoppler_29_diff.fits"
	savfile = '/home/angel/IDLWorkspace/Flare/F29/HMIDoppler/Work/HMIDoppler_coord_29.sav'
	#flare = '/home/angel/Newcodes/fourier/5_5-7_diff.fits'

	x0 = 134
	y0 = 155
	theta0 = 130
	theta1 = 150
	radius = None
	radius = 130

	# Values (parameters) which are set by default
	image = td(flare,  x0=x0, y0=y0, theta0=theta0, theta1=theta1,radius0=0,radius=radius,xticks=7,time0=280,time1=400) 
	
	# Plot the td image

	final = image.tdplot(save=None,plot=True)
	#Points=plt.ginput(4)
	#print(Points)

	#testing = image.test(rows=3,columns=3)

	#colorbar = image.cbar_slider()

	#slide = image.slider()










	exit()
	# Test to make many plots around the given center (x0, y0)
	testing = image.test(rows=3,columns=3)
	
	# Just for change the visualization
	colorbar = image.cbar_slider()
	
	# Slider across the angles and pixels
	slide = image.slider()


	exit()
	savfile = "/home/angel/HMIDoppler_coord_29.sav"
	flare = "/home/angel/HMIDoppler_diff_29.fits"

	# Values (parameters) which are set by default
	image = td(flare, x0=190, y0=145, theta0=110, theta1=140, \
		savfile=None,radius0=0,radius=130) 
	
	# Plot the td image

	final = image.tdplot(save=None,plot=True)
	Points=plt.ginput(4)
	print(Points)

	exit()
	# Test to make many plots around the given center (x0, y0)
	testing = image.test(rows=3,columns=3)
	
	# Just for change the visualization
	colorbar = image.cbar_slider()
	
	# Slider across the angles and pixels
	slide = image.slider()

