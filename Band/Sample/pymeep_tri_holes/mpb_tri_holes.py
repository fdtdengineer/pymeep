# Band dispersion of tri-rods, circle hole

from __future__ import division

import math
import meep as mp
from meep import mpb

# 2d system: triangular lattice of air holes in dielectric
# This structure has a complete band gap (i.e. a gap in both TE and TM
# simultaneously) for a hole radius of 0.45a and a dielectric constant of
# 12.   (See, e.g., the book "Photonic Crystals" by Joannopoulos et al.)

# first, define the lattice vectors and k-points for a triangular lattice:

geometry_lattice = mp.Lattice(size=mp.Vector3(1, 1),
                              basis1=mp.Vector3(math.sqrt(3) / 2, 0.5),
                              basis2=mp.Vector3(math.sqrt(3) / 2, -0.5))

kz = 0  # use non-zero kz to consider vertical propagation

k_points = [
    mp.Vector3(z=kz),               # Gamma
    mp.Vector3(1 / -3, 1 / 3, kz),  # K
    mp.Vector3(0, 0.5, kz),         # M
    mp.Vector3(z=kz)                # Gamma
]

k_interp = 4
k_points = mp.interpolate(k_interp, k_points)

# Now, define the geometry, etcetera:

eps = 11.7  # the dielectric constant of the background
r = 0.3  # the hole radius

default_material = mp.Medium(epsilon=eps)
geometry = [mp.Cylinder(r, material=mp.air)]

resolution = 32
num_bands = 8

ms = mpb.ModeSolver(
    geometry_lattice=geometry_lattice,
    geometry=geometry,
    k_points=k_points,
    default_material=default_material,
    resolution=resolution,
    num_bands=num_bands
)

ms.run_tm(mpb.output_at_kpoint(mp.Vector3(-1./3, 1./3), mpb.fix_efield_phase,
          mpb.output_efield_z))
tm_freqs = ms.all_freqs
tm_gaps = ms.gap_list
ms.run_te()
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
ax.set_ylim([0, 0.6])
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