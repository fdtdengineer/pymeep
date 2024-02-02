#%%
import meep as mp
import numpy as np
from matplotlib import pyplot as plt
from IPython.display import Video

# https://github.com/NanoComp/meep/blob/master/python/examples/holey-wvg-cavity.ipynb

resolution = 20  # pixels/um

eps = 13  # dielectric constant of waveguide
w = 1.2  # width of waveguide
r = 0.36  # radius of holes
d = 1.4  # defect spacing (ordinary spacing = 1)
N = 3  # number of holes on either side of defect

sy = 6  # size of cell in y direction (perpendicular to wvg.)
pad = 2  # padding between last hole and PML edge
dpml = 1  # PML thickness

sx = 2 * (pad + dpml + N) + d - 1  # size of cell in x directionsx = 2 * (pad + dpml + N) + d - 1  # size of cell in x direction


def sim_cavity(N=3, sy=6, fcen=0.25, df=0.2):
    sx = 2 * (pad + dpml + N) + d - 1  # size of cell in x direction
    cell = mp.Vector3(sx, sy, 0)
    blk = mp.Block(size=mp.Vector3(mp.inf, w, mp.inf), material=mp.Medium(epsilon=eps))
    geometry = [blk]

    for i in range(N):
        geometry.append(mp.Cylinder(r, center=mp.Vector3(d / 2 + i)))
        geometry.append(mp.Cylinder(r, center=mp.Vector3(-(d / 2 + i))))

    src = [mp.Source(mp.GaussianSource(fcen, fwidth=df), mp.Hz, mp.Vector3(0))]
    sym = [mp.Mirror(mp.Y, phase=-1), mp.Mirror(mp.X, phase=-1)]
    pml_layers = [mp.PML(1.0)]


    sim = mp.Simulation(
        cell_size=cell,
        geometry=geometry,
        boundary_layers=pml_layers,
        sources=src,
        symmetries=sym,
        resolution=resolution,
    )

    return sim


sim = sim_cavity()
f = plt.figure(dpi=100)
sim.plot2D(ax=f.gca())
plt.show()

# Get eigen-modes
fcen = 0.25  # pulse center frequency
df = 0.2  # pulse frequency width
h = mp.Harminv(mp.Hz, mp.Vector3(), fcen, df)
sim.run(mp.after_sources(h), until_after_sources=400)

#f = plt.figure(dpi=150)
#animate = mp.Animate2D(f=f, fields=mp.Hz, realtime=False, normalize=True)
#sim.run(mp.at_every(1 / fcen / 20, animate), until=1 / fcen)
#plt.close()
f = [m.freq for m in h.modes]
Q = [m.Q for m in h.modes]

for fiter, qiter in zip(f, Q):
    print(f"Resonant frequency: {fiter}, Q: {qiter}")

# %%
