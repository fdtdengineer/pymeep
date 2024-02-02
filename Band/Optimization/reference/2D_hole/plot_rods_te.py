import matplotlib.pyplot as plt
import numpy as np

te_freqs = np.loadtxt('te.dat', delimiter=',', skiprows=1, usecols=(1,6,7,8,9,10,11,12,13))
# tm_freqs = np.loadtxt('tm.dat', delimiter=',', skiprows=1, usecols=(1,6,7,8,9,10,11,12,13))


fig, ax = plt.subplots()
x = range(len(te_freqs))
# Plot bands
# Scatter plot for multiple y values, see https://stackoverflow.com/a/34280815/2261298
for xz, tez in zip(x,te_freqs):
    ax.scatter([xz]*len(tez), tez, color='red', facecolors='none')
ax.plot(te_freqs, color='red')
ax.set_ylim([0, 1])
ax.set_xlim([x[0], x[-1]])

# Plot labels
ax.text(13.05, 0.235, 'TE bands', color='red', size=15)

points_in_between = (len(te_freqs) - 4) / 3
tick_locs = [i*points_in_between+i for i in range(4)]
tick_labs = ['Γ', 'K', 'M', 'Γ']
ax.set_xticks(tick_locs)
ax.set_xticklabels(tick_labs, size=16)
ax.set_ylabel('frequency (c/a)', size=16)
ax.grid(True)

plt.show()