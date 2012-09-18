from numpy import *
import matplotlib.pyplot as plt
# from sh import md
arr = genfromtxt('temperature.dat', dtype=None)

plt.plot(arr[:,0],arr[:,1])
plt.show()