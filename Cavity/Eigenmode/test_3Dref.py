#%%
if True:
    import geometry
    import parse_to_meep
    import transmittance
    import numpy as np
    import pandas as pd
    import meep as mp
    import matplotlib.pyplot as plt
    #from mayavi import mlab
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


# 3D component
hslab = 0.5
h=8

# connection waveguide
len_siwg = 10 # length of straight input/output waveguide
eps_r=3.48**2
fcen=0.25
df=0.1
resolution=8#16
pml_buffer=2
t_after_sources=600
nfreq = 500
thres_conv = 1e-7 # 1e-7 # convergence threshold



# reference
cls_ref = transmittance.MeepTransmittance(
    dim=3,
    area_z=h, 
    thick_slab=hslab,
    fcen=fcen,
    df=df,
    thres_conv = thres_conv,
    nfreq = nfreq,
    resolution=resolution,
    dpml=1
    )
freqs_ref, psd_ref = cls_ref.get_reference_transmittance(
    a=a, 
    nx=nx, 
    ny=ny, 
    len_siwg=len_siwg,
    eps_r=eps_r,
    wg=wgi, 
    pml_buffer=pml_buffer,
    )

lattice_const = 400
wl = lattice_const / freqs_ref
df_ref = pd.DataFrame(np.array([freqs_ref, wl, psd_ref]).T, columns=["freq","wl","transmittance"])
df_ref.to_csv(r"outputs/test_3Dref_230729.csv")

# plot
plt.figure(figsize=(8,6))
plt.plot(freqs_ref, psd_ref, label="W1", marker="o", markersize=3)
plt.xlabel("Frequency")
plt.ylabel("Intensity")
plt.yscale("log")
plt.legend()
plt.tight_layout()
plt.savefig(r"outputs/test_3Dref_230729.svg")
plt.show()




#print("freqs:", freqs)


# %%
