import matplotlib.pyplot as plt
import numpy as np

# Define the band structure equation
def E(u, m_1, m_2, m_3):
	return 1./3.*( (u-2*(-m_1+m_2+m_3))**2 + \
			(u-2*(m_1-m_2+m_3))**2 + \
			(u-2*(m_1+m_2-m_3))**2 )

# values which k can exist over normalized to [-1,1]
u = np.linspace(-1,1,50)

# Plot each band
fig = plt.figure()
plt.xlabel('u')
plt.ylabel('E')
plt.title('FCC Lattice Band Structure in the Reduced Zone Scheme')
plt.plot(u, E(u,0,0,0))
plt.plot(u, E(u,1,0,0))
plt.plot(u, E(u,-1,0,0))
plt.plot(u, E(u,1,1,0))
plt.plot(u, E(u,-1,-1,0))
plt.plot(u, E(u,1,1,1))
plt.plot(u, E(u,-1,-1,-1))
fig.savefig('FCCBandStructure.png')
plt.show()


