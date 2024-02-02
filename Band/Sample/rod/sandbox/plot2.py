import matplotlib.pyplot as plt
import numpy as np

te_freqs = np.loadtxt('te.dat', delimiter=',', skiprows=1, usecols=(1,6,7,8,9,10,11,12,13))

print(te_freqs)

