#%%
if True:
    import geometry
    import parse_to_meep    
    import numpy as np
    import pandas as pd
    import meep as mp
    import matplotlib.pyplot as plt
    #from matplotlib import rc
    #rc('text', usetex=False)
    plt.rcParams['font.family']= 'sans-serif'
    #plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams["font.size"] = 15 # 全体のフォントサイズが変更されます。
    plt.rcParams['xtick.direction'] = 'in'#x軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')
    plt.rcParams['ytick.direction'] = 'in'#y軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')


# mpirun -np 4 python Cavity2D\harminv_cavity.py
##### settings of geometry #####
a = 1 # 0.4
nx = 8 # 14
ny = 30 # 50
offset_x = 0
offset_y = 0
n_cavity = 5 # 3
barrier = 4
wgi = 0 #1.1
holeshift = 0.2
hslab=0.5 #mp.inf
hall=8*2

eps_r=3.48#2.6**2
fcen=0.27
df=0.05
resolution=8
pml_buffer=4
t_after_sources=3000

##### settings of simulation #####
# 2D cavity geometry
cell = mp.Vector3(ny+pml_buffer, (nx+pml_buffer/2)*np.sqrt(3),hall)
blk = mp.Block(size=mp.Vector3(ny+pml_buffer, (nx+pml_buffer/2)*np.sqrt(3), hslab), material=mp.Medium(epsilon=eps_r))
ld = geometry.LineDefect(a, nx, ny, offset_x, offset_y, n_cavity, barrier, wgi, holeshift)
arr_obj = parse_to_meep.parse_geometry(ld, thick_slab=hslab)
arr_geometry = [blk] + arr_obj


src = [
    mp.Source(mp.GaussianSource(fcen, fwidth=df), mp.Ey, mp.Vector3(0.), amplitude=1),
    #mp.Source(mp.GaussianSource(fcen, fwidth=df), mp.Ey, mp.Vector3(1,0.1), amplitude=1),
    #mp.Source(mp.GaussianSource(fcen, fwidth=df), mp.Ey, mp.Vector3(-1), amplitude=1),
    #mp.Source(mp.GaussianSource(fcen, fwidth=df), mp.Hz, mp.Vector3(-0.5),amplitude=-1),
       ]

#src = [
#    mp.Source(mp.GaussianSource(fcen, fwidth=df), mp.Hz, mp.Vector3(0), amplitude=1),
#       ]

#sym = [mp.Mirror(mp.Y, phase=-1), mp.Mirror(mp.X, phase=-1)]
#sym = [mp.Mirror(mp.Y, phase=-1), mp.Mirror(mp.X, phase=1)]
#sym = [mp.Mirror(mp.Y, phase=-1), mp.Mirror(mp.X, phase=1)]
sym = [mp.Mirror(mp.Y, phase=-1), mp.Mirror(mp.X, phase=1), mp.Mirror(mp.Z, phase=1)]

#sym = []

pml_layers = [mp.PML(1.0)]


sim = mp.Simulation(
    cell_size=cell,
    geometry=arr_geometry,
    boundary_layers=pml_layers,
    sources=src,
    symmetries=sym,
    resolution=resolution,
    dimensions=3
)

# check
#f = plt.figure(dpi=100)
#sim.plot2D(ax=f.gca())
#plt.show()


##### Get eigen-modes #####
h = mp.Harminv(mp.Ey, mp.Vector3(1), fcen, df)
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

#%%
##### Get the field distributions of eigen-modes #####
# To see the mode, we will simply run the simulation again with a narrow-band source.

fcen_cavity=0.267 #0.27111168664110774 # 0.27685639922705274
df_cavity=0.01
src = [
    #mp.Source(mp.GaussianSource(fcen, fwidth=df), mp.Ey, mp.Vector3(0.05,0.05), amplitude=1),
    mp.Source(mp.GaussianSource(fcen, fwidth=df), mp.Ey, mp.Vector3(0), amplitude=1),
       ]
sim_cavity = mp.Simulation(
    cell_size=cell,
    geometry=arr_geometry,
    boundary_layers=pml_layers,
    sources=src,
    symmetries=sym,
    resolution=resolution,
)

sim_cavity.reset_meep()
h_cavity = mp.Harminv(mp.Ey, mp.Vector3(0.5,0.,0), fcen_cavity, df_cavity)
sim_cavity.run(mp.after_sources(h_cavity), until_after_sources=t_after_sources)

f_cavity = [m.freq for m in h_cavity.modes]
Q_cavity = [m.Q for m in h_cavity.modes]
df_cavity_results = pd.DataFrame(np.array([f_cavity,Q_cavity]).T)

f = plt.figure(dpi=150)
sim_cavity.plot2D(ax=f.gca(), fields=mp.Hz)
plt.show()
print(df_cavity_results)
print("Done")



# %%
Ey=sim.get_efield_y()
arr = Ey[:,:,64].T
vmax = np.max(arr) /100
plt.imshow(arr,cmap="bwr", vmax=vmax, vmin=-vmax)
