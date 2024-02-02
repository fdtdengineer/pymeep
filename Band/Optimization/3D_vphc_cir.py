#%%
import math
import meep as mp
from meep import mpb
import numpy as np
import matplotlib.pyplot as plt

n_Air = 1
n_Si = 3.48

a = 1
#h = 1.2*a
#r1 = 0.14*a
#r2 = 0.38*a

h, r1, r2 = 1.2, 0.38, 0.14
#計算する固有周波数の数
num_bands = 3

#メッシュの細かさ
resolution = 16

#Γ-K, K-M, M-Γ間の点の個数
num_k = 5#10

#単位格子
geometry_lattice = mp.Lattice(size=mp.Vector3(1, 1, 10*h),
                              basis1=mp.Vector3(1./2, np.sqrt(3)/2),
                              basis2=mp.Vector3(1./2, -np.sqrt(3)/2))


geometry = [
    mp.Block(material=mp.Medium(epsilon=1),
             size=mp.Vector3(mp.inf, mp.inf, 10*h)),
    mp.Block(material=mp.Medium(epsilon=n_Si**2),
             size=mp.Vector3(mp.inf, mp.inf, h)),
    mp.Cylinder(r1, center=mp.Vector3(2, 1)/3, 
                material=mp.Medium(epsilon=n_Air**2)),
    mp.Cylinder(r2, center=mp.Vector3(1, 2)/3, 
                material=mp.Medium(epsilon=n_Air**2)),
]

#ブリルアンゾーン
k_points = [
    mp.Vector3(),               # Gamma
    mp.Vector3(1./3, 1./3),    # K
    mp.Vector3(1./2, 0),          # M
    mp.Vector3(),               # Gamma
]
k_points = mp.interpolate(num_k, k_points)


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
md = mpb.MPBData(rectify=True, periods=3, resolution=16)
eps = ms.get_epsilon()
converted_eps = md.convert(eps)
plt.imshow(converted_eps[:,:,converted_eps.shape[2]//2].T, interpolation='spline36', cmap="binary")
plt.axis('off')
plt.show()

#%%
# plot
fs = 20
fig, ax = plt.subplots(figsize=(5,5))
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

# lightline
x1 = [0,num_k+1]
y11 = [0, 2/3]
y12 = np.ones_like(y11)
x2 = [2*(num_k+1), 3*(num_k+1),]
y21 = [0.577350269189626, 0]
y22 = np.ones_like(y21)

ax.plot(x1, y11, color="black")
plt.fill_between(x1, y11, y12, where=y11<y12, facecolor="gray", alpha=1)
ax.plot(x2,y21, color="black")
plt.fill_between(x2, y21, y22, where=y21<y22, facecolor="gray", alpha=1)

# Plot labels
ax.text(num_k, 0.3, 'TM bands', color='lightseagreen', size=fs)
ax.text(num_k, 0.2, 'TE bands', color='black', size=fs)

points_in_between = (len(te_freqs) - 4) / 3
tick_locs = [i*points_in_between+i for i in range(4)]
tick_labs = ['Γ', 'K', 'M', 'Γ']
ax.set_xticks(tick_locs)
ax.set_xticklabels(tick_labs, size=fs)
ax.set_ylabel(r'$\omega a / 2\pi c$', size=fs)
ax.grid(True)
plt.tick_params(labelsize=fs)

plt.show()


idx_K = num_k+1
idx_M = 2*(num_k+1)

f_te1_K = te_freqs[idx_K,0]
f_te2_K = te_freqs[idx_K,1]
f_te1_M = te_freqs[idx_M,0]
f_te2_M = te_freqs[idx_M,1]
f_tm1_K = tm_freqs[idx_K,0]
f_tm2_K = tm_freqs[idx_K,1]
f_tm1_M = tm_freqs[idx_M,0]
f_tm2_M = tm_freqs[idx_M,1]

lap_1 = max(f_tm2_K, f_tm2_M) - min(f_te1_K, f_te1_M)
lap_2 = max(f_te2_K, f_te2_M) - min(f_tm1_K, f_tm1_M)
lap = min(lap_1, lap_2)

print("Overlap bandgap between band 1 and 2:", lap)

# %%
