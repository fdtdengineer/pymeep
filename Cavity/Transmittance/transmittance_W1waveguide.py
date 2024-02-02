#%%
if True:
    import geometry
    import parse_to_meep
    import transmittance
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
nx = 10 # 14
ny = 30 # 50
offset_x = 0
offset_y = 0
barrier = 6
wgo = 1
wgi = 1
holeshift = 0
dpml=1
# connection waveguide
len_siwg = 10 # length of straight input/output waveguide
eps_r=2.6**2
fcen=0.27
df=0.05
resolution=16
pml_buffer=2
t_after_sources=600
nfreq = 500
thres_conv = 1e-7 # 1e-7 # convergence threshold


##### settings of simulation #####
# 2D cavity geometry
# Area of simulation cell
cell = mp.Vector3(a*(ny+pml_buffer+len_siwg*2+dpml-1), a*(nx+2)*np.sqrt(3),0)

# Bulk material
arr_blk = [mp.Block(
    size=mp.Vector3(a*ny, a*(nx+2)*np.sqrt(3), mp.inf), 
    material=mp.Medium(epsilon=eps_r)
)]

# Waveguide
arr_wg = [
    mp.Block(
        mp.Vector3(a*(len_siwg+pml_buffer+dpml-1), wgi*np.sqrt(3),mp.inf),
        center=mp.Vector3(-a*(ny/2+(len_siwg+pml_buffer+dpml-1)/2), 0),
        material=mp.Medium(epsilon=eps_r)
    ),
    mp.Block(
        mp.Vector3(a*(len_siwg+pml_buffer+dpml-1), wgi*np.sqrt(3),mp.inf),
        center=mp.Vector3(a*(ny/2+(len_siwg+pml_buffer+dpml-1)/2), 0),
        material=mp.Medium(epsilon=eps_r)
    ),
]

# Photonics crystal cavity
wm = geometry.WidthModulated(a, nx, ny, offset_x, offset_y, barrier, wgo, wgi, holeshift)
arr_obj = parse_to_meep.parse_geometry(wm)

# All Geometry
arr_geometry = arr_blk + arr_wg + arr_obj

# Source
src = [
    mp.Source(
        mp.GaussianSource(fcen, fwidth=df),
        component=mp.Hz,
        center=mp.Vector3(-a*(ny/2+len_siwg-(dpml-1)/2), 0),
        size=mp.Vector3(0, wgi*np.sqrt(3))
    )
]

# Symmetry
sym = [mp.Mirror(mp.Y, phase=-1)]

# PML
pml_layers = [mp.PML(dpml)]

# Build simulation object
sim = mp.Simulation(
    cell_size=cell,
    geometry=arr_geometry,
    boundary_layers=pml_layers,
    sources=src,
    symmetries=sym,
    resolution=resolution,
)

# Output flux
vec_out = mp.Vector3(a*(ny/2+len_siwg-1-(dpml-1)/2), 0)
tran_out = mp.FluxRegion(center=vec_out,size=mp.Vector3(0, 2*wgi))
trans_out = sim.add_flux(fcen, df, nfreq, tran_out)

# Plot geometry
f = plt.figure(dpi=100)
sim.plot2D(ax=f.gca())
plt.show()


# Run simulation
sim.run(
    until_after_sources=mp.stop_when_fields_decayed(
        50, mp.Hz, vec_out, thres_conv
    )
)

freqs = np.array(mp.get_flux_freqs(trans_out))
psd_out = np.array(mp.get_fluxes(trans_out))


# reference
cls_ref = transmittance.MeepTransmittance(dim=2,fcen=fcen, df=df, nfreq=nfreq, resolution=resolution,  thres_conv=thres_conv)
freqs_ref, psd_ref = cls_ref.get_reference_transmittance()


lattice_const = 400
wl = lattice_const / freqs_ref
df = pd.DataFrame(np.array([freqs_ref, wl, psd_out / psd_ref]).T, columns=["freq","wl","transmittance"])
df.to_csv(r"outputs/transmittance_W1waveguide.csv")

lattice_const = 400
wl = lattice_const / freqs_ref

# plot
plt.figure(figsize=(8,6))
plt.plot(freqs, psd_out / psd_ref, label="W1", marker="o", markersize=3)
plt.xlabel("Frequency")
plt.ylabel("Intensity")
plt.yscale("log")
plt.legend()
plt.tight_layout()
plt.savefig(r"outputs/transmittance_W1waveguide.svg")
plt.show()


plt.figure(figsize=(8,6))
plt.plot(wl, psd_out / psd_ref, label="W1", marker="o", markersize=3, color="black")
plt.xlabel("Wavelength [nm]")
plt.ylabel("Intensity")
plt.yscale("log")
plt.legend()
plt.tight_layout()
plt.savefig(r"outputs/transmittance_wl_W1waveguide.svg")
plt.show()




#print("freqs:", freqs)


# %%
