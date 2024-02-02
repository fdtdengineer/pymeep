from __future__ import division

import math
import meep as mp
from meep import mpb

# A honeycomb lattice of dielectric rods in air.  (This structure has
# a complete (overlapping TE/TM) band gap.)  A honeycomb lattice is really
# just a triangular lattice with two rods per unit cell, so we just
# take the lattice, k_points, etcetera from mpb_tri_rods.py.

r = 0.2  # the rod radius
eps = 6.25  # the rod dielectric constant

# triangular lattice:
geometry_lattice = mp.Lattice(size=mp.Vector3(1, 1),
                              basis1=mp.Vector3(math.sqrt(3) / 2, 0.5),
                              basis2=mp.Vector3(math.sqrt(3) / 2, -0.5))

# Two rods per unit cell, at the correct positions to form a honeycomb
# lattice, and arranged to have inversion symmetry:
geometry = [mp.Cylinder(r, center=mp.Vector3(1 / 6, 1 / 6), height=mp.inf,
                        material=mp.Medium(epsilon=eps)),
            mp.Cylinder(r, center=mp.Vector3(1 / -6, 1 / -6), height=mp.inf,
                        material=mp.Medium(epsilon=eps))]

# The k_points list, for the Brillouin zone of a triangular lattice:
k_points = [
    mp.Vector3(),               # Gamma
    mp.Vector3(1 / -3, 1 / 3),  # K
    mp.Vector3(y=0.5),          # M
    mp.Vector3()                # Gamma
]

k_interp = 4  # number of k_points to interpolate
k_points = mp.interpolate(k_interp, k_points)

resolution = 32
num_bands = 8

ms = mpb.ModeSolver(
    geometry_lattice=geometry_lattice,
    geometry=geometry,
    k_points=k_points,
    resolution=resolution,
    num_bands=num_bands
)


def main():
    ms.run_tm()
    ms.run_te()

# Since there is a complete gap, we could instead see it just by using:
# run()
# The gap is between bands 12 and 13 in this case.  (Note that there is
# a false gap between bands 2 and 3, which disappears as you increase the
# k_point resolution.)

if __name__ == '__main__':
    main()


tm_freqs = ms.all_freqs
tm_gaps = ms.gap_list
te_freqs = ms.all_freqs
te_gaps = ms.gap_list

###################

import matplotlib.pyplot as plt

md = mpb.MPBData(rectify=True, periods=3, resolution=32)
eps = ms.get_epsilon()
converted_eps = md.convert(eps)
plt.imshow(converted_eps.T, interpolation='spline36', cmap='binary')
plt.axis('off')
plt.show()



# Display the band dispersion

fig, ax = plt.subplots()
x = range(len(tm_freqs))
# Plot bands
# Scatter plot for multiple y values, see https://stackoverflow.com/a/34280815/2261298
for xz, tmz, tez in zip(x, tm_freqs, te_freqs):
    ax.scatter([xz]*len(tmz), tmz, color='blue')
    ax.scatter([xz]*len(tez), tez, color='red', facecolors='none')
ax.plot(tm_freqs, color='blue')
ax.plot(te_freqs, color='red')
ax.set_ylim([0, 1])
ax.set_xlim([x[0], x[-1]])

# Plot gaps
for gap in tm_gaps:
    if gap[0] > 1:
        ax.fill_between(x, gap[1], gap[2], color='blue', alpha=0.2)

for gap in te_gaps:
    if gap[0] > 1:
        ax.fill_between(x, gap[1], gap[2], color='red', alpha=0.2)


# Plot labels
ax.text(12, 0.04, 'TM bands', color='blue', size=15)
ax.text(13.05, 0.235, 'TE bands', color='red', size=15)

points_in_between = (len(tm_freqs) - 4) / 3
tick_locs = [i*points_in_between+i for i in range(4)]
tick_labs = ['Γ', 'K', 'M', 'Γ']
ax.set_xticks(tick_locs)
ax.set_xticklabels(tick_labs, size=16)
ax.set_ylabel('frequency (c/a)', size=16)
ax.grid(True)

plt.show()