import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
from pint import UnitRegistry
ureg = UnitRegistry()
ureg.define('electron_volt = 1.60217657*10**(-19)*joules')

dk_1d = 0.1

a_cc = 0.142 #* ureg.nanometer
a_1 = a_cc/2 * np.array([3, np.sqrt(3)]) 
a_2 = a_cc/2 * np.array([3, -np.sqrt(3)]) 

b_1 = 2*np.pi/(3*a_cc) * np.array([1, np.sqrt(3)]) 
b_2 = 2*np.pi/(3*a_cc) * np.array([1, -np.sqrt(3)]) 

K = 2*np.pi/(3*a_cc) * np.array([1, 1/np.sqrt(3)])
K_prime = 2*np.pi/(3*a_cc) * np.array([1, -1/np.sqrt(3)]) 

def E(k_x, k_y):
	t = 2.8 #* ureg.electron_volt
	return np.array([t*np.sqrt(3+f(k_x, k_y)), -t*np.sqrt(3+f(k_x, k_y))])

def f(k_x, k_y):
	return 2*np.cos(np.sqrt(3)*k_y*a_cc)+ \
               4*np.cos(np.sqrt(3)/2*k_y*a_cc)*np.cos(3/2*k_x*a_cc)

fig = plt.figure()

im = plt.imread('Arpes3crop.png')

k_x = np.arange(-2.0, 2.0, dk_1d)
k_y = np.ones(len(k_x)) * np.sqrt(K[0]**2+K[1]**2)
plt.imshow(im, extent=[-2.0, 2.0,-1.2+0.2, 0.2+0.2])
plt.plot(k_x, E(k_x, k_y)[0,:], 'r', k_x, E(k_x, k_y)[1,:], 'k')
plt.title('Graphene Tight Binding Brillouin Zone Band Structure')
plt.xlabel(r'$k_x$ ($\frac{1}{nm}$)')
plt.ylabel(r'Energy (eV)')
#plt.show()
plt.savefig('TBvsARPES46Bicent.png')


