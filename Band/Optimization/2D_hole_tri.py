#%%
import math
import meep as mp
from meep import mpb
import numpy as np
import matplotlib.pyplot as plt

#屈折率
n_Air = 1
n_Si = 2.6
radius = 0.25 # a
s0 = 0.5 # one side of triangular hole
h = 20

#計算する固有周波数の数
num_bands = 4

#メッシュの細かさ
resolution = 32

#Γ-K, K-M, M-Γ間の点の個数
num_of_kpoint = 30

#単位格子
geometry_lattice = mp.Lattice(size=mp.Vector3(1, 1),
                              basis1=mp.Vector3(1./2, np.sqrt(3)/2),
                              basis2=mp.Vector3(1./2, -np.sqrt(3)/2))

#構造
#default_material = mp.Medium(epsilon=n_Si**2)
hole_triangular_up = [
        mp.Vector3(-1,-1)*s0,
        mp.Vector3(1,0)*s0,
        mp.Vector3(0,1)*s0
      ]
hole_triangular2 = [
        mp.Vector3(0.5, -0.5) * s0,
        mp.Vector3(0.5,    1) * s0,
        mp.Vector3( -1, -0.5) * s0
      ]

radius_pole = 0.15
s1 = s0*0.7
geometry = [
    mp.Block(material=mp.Medium(epsilon=n_Si**2),
             size=mp.Vector3(mp.inf, mp.inf)),
    #mp.Cylinder(radius, material=mp.Medium(epsilon=n_Air**2)),
    mp.Prism(hole_triangular2, center=mp.Vector3(0), height=h,
            material=mp.Medium(epsilon=n_Air**2)),
    mp.Cylinder(radius_pole, center=mp.Vector3(0.5, -0.5)*s1, 
                material=mp.Medium(epsilon=n_Air**2)),
    mp.Cylinder(radius_pole, center=mp.Vector3(0.5,    1)*s1, 
                material=mp.Medium(epsilon=n_Air**2)),
    mp.Cylinder(radius_pole, center=mp.Vector3( -1, -0.5)*s1, 
                material=mp.Medium(epsilon=n_Air**2)),
]


# ▷


#ブリルアンゾーン
k_points = [
    mp.Vector3(),               # Gamma
    mp.Vector3(1./3, 1./3),    # K
    mp.Vector3(1./2, 0),          # M
    mp.Vector3(),               # Gamma
]
k_points = mp.interpolate(num_of_kpoint, k_points)


#計算
ms = mpb.ModeSolver(
    geometry=geometry,
    geometry_lattice=geometry_lattice,
    k_points=k_points,
    resolution=resolution,
    num_bands=num_bands
)

ms.run_tm(mpb.output_at_kpoint(mp.Vector3(1./3, 1./3),
                               mpb.fix_efield_phase,
                               mpb.output_efield_z))
tm_freqs = ms.all_freqs
tm_gaps = ms.gap_list
ms.run_te(mpb.output_at_kpoint(mp.Vector3(1./3, 1./3),
                               mpb.fix_efield_phase,
                               mpb.output_efield_z))
te_freqs = ms.all_freqs
te_gaps = ms.gap_list


# show permittivity
md = mpb.MPBData(rectify=True, periods=3, resolution=64)
eps = ms.get_epsilon()
converted_eps = md.convert(eps)
plt.imshow(converted_eps.T, interpolation='spline36', cmap="binary")
plt.axis('off')
plt.show()


# plot
fs = 20
fig, ax = plt.subplots(figsize=(6,6))
x = range(len(te_freqs))

ax.plot(tm_freqs, color='lightseagreen', linestyle="--")
ax.plot(te_freqs, color='black')
ax.set_ylim([te_freqs.min(), te_freqs.max()])
ax.set_xlim([x[0], x[-1]])

# Plot gaps
for gap in tm_gaps:
    if gap[0] > 1:
        ax.fill_between(x, gap[1], gap[2], color='turquoise', alpha=0.2)

for gap in te_gaps:
    if gap[0] > 1:
        ax.fill_between(x, gap[1], gap[2], color='black', alpha=0.2)



# Plot labels
ax.text(13.05, 0.235, 'TM bands', color='lightseagreen', size=fs)
ax.text(12, 0.04, 'TE bands', color='black', size=fs)

points_in_between = (len(te_freqs) - 4) / 3
tick_locs = [i*points_in_between+i for i in range(4)]
tick_labs = ['Γ', 'K', 'M', 'Γ']
ax.set_xticks(tick_locs)
ax.set_xticklabels(tick_labs, size=fs)
ax.set_ylabel('frequency (c/a)', size=fs)
ax.grid(True)
plt.tick_params(labelsize=fs)
#plt.savefig("2D_pillar_honeycomb_π_6_normFreq.png", bbox_inches="tight")

plt.show()



# %%
