import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
from pint import UnitRegistry
ureg = UnitRegistry()
ureg.define('electron_volt = 1.60217657*10**(-19)*joules')

dk_3d = 0.5
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

k_x = np.arange(-27, 27, dk_3d)
k_y = np.arange(-24, 24, dk_3d)
k_x, k_y = np.meshgrid(k_x, k_y)
Z_p = E(k_x, k_y)[0,:,:]
Z_n = E(k_x, k_y)[1,:,:]

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(k_x, k_y, Z_p, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
surf2 = ax.plot_surface(k_x, k_y, Z_n, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
plt.title('Graphene Tight Binding Brillouin Zone Band Structure')
plt.xlabel(r'k_x ($\frac{1}{nm}$)')
plt.ylabel(r'k_y ($\frac{1}{nm}$)')
ax.set_zlabel('Energy (eV)')

plt.savefig('BZ_TightBinding.png')

fig = plt.figure()

k_x = np.arange(-K[0], K[0], dk_1d)
k_y = np.ones(len(k_x)) * np.sqrt(K[0]**2+K[1]**2)
plt.plot(k_x, E(k_x, k_y)[0,:], 'r', k_x, E(k_x, k_y)[1,:], 'k')
plt.title('Graphene Tight Binding Brillouin Zone Band Structure')
plt.xlabel(r'k_x ($\frac{1}{nm}$)')
plt.ylabel(r'Energy (eV)')
plt.savefig('Kpoint_TightBinding.png')


