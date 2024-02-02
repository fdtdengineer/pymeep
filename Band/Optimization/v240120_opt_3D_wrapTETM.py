#%%
import meep as mp
from meep import mpb
import numpy as np

# configuration
n_Air = 1
n_Si = 3.48
a = 1
num_bands = 3 #計算する固有周波数の数
resolution = 16 #メッシュの細かさ
num_k = 0#5#10 #Γ-K, K-M, M-Γ間の点の個数

"""
s0 ... one side of triangular hole
h = 0.5*a ... slab thickness
radius_pole ... radius of pole circles
s1 = s0 * 0.5 ... position of pole circles from center triangle
"""

def gap_tri(s0, h, radius_pole, s1):

    #単位格子
    geometry_lattice = mp.Lattice(size=mp.Vector3(1, 1, 10*h),
                                basis1=mp.Vector3(1./2, np.sqrt(3)/2),
                                basis2=mp.Vector3(1./2, -np.sqrt(3)/2))

    #構造
    hole_triangular2 = [
            mp.Vector3(0.5, -0.5) * s0 / np.sqrt(3),
            mp.Vector3(0.5,    1) * s0 / np.sqrt(3),
            mp.Vector3( -1, -0.5) * s0 / np.sqrt(3),
        ]

    geometry = [
        mp.Block(material=mp.Medium(epsilon=1),
                size=mp.Vector3(mp.inf, mp.inf, 10*h)),
        mp.Block(material=mp.Medium(epsilon=n_Si**2),
                size=mp.Vector3(mp.inf, mp.inf, h)),
        mp.Prism(hole_triangular2, center=mp.Vector3(0), height=float(h),
                material=mp.Medium(epsilon=n_Air**2)),
        mp.Cylinder(radius_pole, center=mp.Vector3(0.5, -0.5)*s1, 
                    material=mp.Medium(epsilon=n_Air**2)),
        mp.Cylinder(radius_pole, center=mp.Vector3(0.5,    1)*s1, 
                    material=mp.Medium(epsilon=n_Air**2)),
        mp.Cylinder(radius_pole, center=mp.Vector3( -1, -0.5)*s1, 
                    material=mp.Medium(epsilon=n_Air**2)),
    ]

    #ブリルアンゾーン
    k_points = [
        #mp.Vector3(),               # Gamma
        mp.Vector3(1./3, 1./3),    # K
        mp.Vector3(1./2, 0),          # M
        #mp.Vector3(),               # Gamma
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
    ms.run_te(mpb.output_at_kpoint(mp.Vector3(1./3, 1./3),
                                mpb.fix_efield_phase,
                                mpb.output_efield_z))
    te_freqs = ms.all_freqs

    idx_K = 0#num_k+1
    idx_M = 1#2*(num_k+1)

    f_te1_K = te_freqs[idx_K,0]
    f_te2_K = te_freqs[idx_K,1]
    f_te1_M = te_freqs[idx_M,0]
    f_te2_M = te_freqs[idx_M,1]
    f_tm1_K = tm_freqs[idx_K,0]
    f_tm2_K = tm_freqs[idx_K,1]
    f_tm1_M = tm_freqs[idx_M,0]
    f_tm2_M = tm_freqs[idx_M,1]

    lap = [0]*6
    # Bandgap of TE, TM
    lap[0] = min(f_tm2_K, f_tm2_M) - max(f_tm1_K, f_tm1_M)
    lap[1] = min(f_te2_K, f_te2_M) - max(f_te1_K, f_te1_M)
    
    # Overlap bandgap between band 1 and 2
    lap[2] = min(f_tm2_K, f_tm2_M) - max(f_te1_K, f_te1_M)
    lap[3] = min(f_te2_K, f_te2_M) - max(f_tm1_K, f_tm1_M)
    
    lap[4] = max(f_tm2_K, f_tm2_M) - min(f_te1_K, f_te1_M)
    lap[5] = max(f_te2_K, f_te2_M) - min(f_tm1_K, f_tm1_M)
    
    minlap = min(lap)
    print("\n\n\n\n\n\n\n")
    print("Overlap bandgap between band 1 and 2:", minlap)
    print("\n\n\n\n\n\n\n")

    return minlap


if __name__ == "__main__":
    #s0 = 0.85*a 
    #h = 0.5*a 
    #radius_pole = 0.1 
    #s1 = s0 * 0.5 
    s0, h, radius_pole, s1= (0.8, 0.7584089890588084, 0.1, 0.45)

    #lap = gap_tri(s0, h, radius_pole, s1)

    # optimization
    import scipy.optimize as opt

    def gap_tri_opt(x):
        s0 = x[0]
        h = x[1]
        radius_pole = x[2]
        s1 = x[3]
        return -gap_tri(s0, h, radius_pole, s1)

    x0 = np.array([s0, h, radius_pole, s1])
    bounds = [(0, 0.9), (0, 1), (0, 1), (0, 1)] 
    res = opt.minimize(gap_tri_opt, x0, bounds=bounds, method='Nelder-Mead')  # "SLSQP", options={"maxiter":1})
    print(res)
    print("s0 = ", res.x[0])
    print("h = ", res.x[1])
    print("radius_pole = ", res.x[2])
    print("s1 = ", res.x[3])


    import pandas as pd
    pd.DataFrame(res.x).to_csv("opt_result.csv")

    print("Done!")



# %%
