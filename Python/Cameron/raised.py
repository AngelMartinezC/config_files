import numpy as np
import matplotlib.pyplot as plt

T = 0.1
b = 0.5


plt.figure(figsize=(6.65,4.74))
def raised(T,b,z0=0):
  z = np.linspace(0,30,1000)
  I = lambda z: (1/2)*(1+np.cos((np.pi*T/b)*(abs(z-z0)-(1-b)/(2*T))))
  H = np.piecewise(z,[z>-1e9+z0,z>-(1+b)/(2*T)+z0,z>-(1-b)/(2*T)+z0,z>(1-b)/(2*T)+z0,z>(1+b)/(2*T)+z0],[0,I,1,I,0])
  
  plt.plot(z,H)

raised(0.05,0.3,0)
raised(0.1,0.6,15)
raised(0.05,0.3,30)


plt.ylim(0,1.2)
plt.xlim(0,30)
plt.grid()
plt.show()















