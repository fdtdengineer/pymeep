#%%
# Calculating 2d ring-resonator modes, from the Meep tutorial.
import meep as mp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

n = 3.4  # index of waveguide
w = 1  # width of waveguide
r = 1  # inner radius of ring
pad = 4  # padding between waveguide and edge of PML
dpml = 2  # thickness of PML
sxy = 2 * (r + w + pad + dpml)  # cell size
resolution = 16 # mesh size

# Create a ring waveguide by two overlapping cylinders - later objects
# take precedence over earlier objects, so we put the outer cylinder first.
# and the inner (air) cylinder second.

c1 = mp.Cylinder(radius=r + w, material=mp.Medium(index=n))
c2 = mp.Cylinder(radius=r)

# If we don't want to excite a specific mode symmetry, we can just
# put a single point source at some arbitrary place, pointing in some
# arbitrary direction.  We will only look for Ez-polarized modes.

fcen = 0.15  # pulse center frequency
df = 0.1  # pulse width (in frequency)

src = mp.Source(mp.GaussianSource(fcen, fwidth=df), mp.Ez, mp.Vector3(r + 0.1))

sim = mp.Simulation(
    cell_size=mp.Vector3(sxy, sxy),
    geometry=[c1, c2],
    sources=[src],
    resolution=resolution,
    symmetries=[mp.Mirror(mp.Y)],
    boundary_layers=[mp.PML(dpml)],
)

f = plt.figure(dpi=100)
sim.plot2D(ax=f.gca())
plt.show()

t_after_sources = 300
h = mp.Harminv(mp.Ez, mp.Vector3(r + 0.1), fcen, df)
sim.reset_meep()
sim.run(mp.after_sources(h), until_after_sources=t_after_sources)


f = [m.freq for m in h.modes]
Q = [m.Q for m in h.modes]

##### Plot eigen-modes #####
for fiter, qiter in zip(f, Q):
    print(f"Resonant frequency: {fiter}, Q: {qiter}")

df_results = pd.DataFrame(np.array([f,Q]).T)
print(df_results)
print("Done")


    

# %%
fcen_cavity=0.118139 # 0.27685639922705274
df_cavity=0.01
sim.reset_meep()
h_cavity = mp.Harminv(mp.Ez, mp.Vector3(r + 0.1), fcen_cavity, df_cavity)
sim.run(mp.after_sources(h_cavity), until_after_sources=t_after_sources*10)

f_cavity = [m.freq for m in h_cavity.modes]
Q_cavity = [m.Q for m in h_cavity.modes]
df_cavity_results = pd.DataFrame(np.array([f_cavity,Q_cavity]).T)
print(df_cavity_results)

f = plt.figure(dpi=150)
sim.plot2D(ax=f.gca(), fields=mp.Ez)
plt.show()

print("Done")

# %%
