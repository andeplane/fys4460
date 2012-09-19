from numpy import *
import matplotlib.pyplot as plt

arr = genfromtxt('energy.dat', dtype=None)

plt.plot(arr[:,0],arr[:,1])
plt.show()