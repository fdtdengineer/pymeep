import time
import datetime
import math
import meep as mp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mayavi import mlab


def phc_trans(PhC = True, lengthPhC = 20, resolution = 8, checkpoint=0):
    #checkpoint = 0 #lengthPhC/2-3/2
    widthPhC = 10

    ConnectionWaveguide = 5
    length = lengthPhC + 2*ConnectionWaveguide
    width = widthPhC

    # Si厚さ
    hslab = 0.5
    # 空気層
    h=15
    
    wgi = 1
    Nx = int(lengthPhC)
    Ny = int(widthPhC)
    r = 1/4
    eps = 11.7

    cell = mp.Vector3(length,width*np.sqrt(3),h)
    
    if PhC == False:
        waveguide = mp.Block(mp.Vector3(mp.inf,wgi*np.sqrt(3),hslab),
                                 center=mp.Vector3(0, 0),
                                 material=mp.Medium(epsilon=eps))
        geometry = [waveguide]
    
    if PhC == True:
        blk = mp.Block(mp.Vector3(lengthPhC,widthPhC*np.sqrt(3),hslab),
                             center=mp.Vector3(),
                             material=mp.Medium(epsilon=eps))

        waveguide = mp.Block(mp.Vector3(mp.inf,wgi*np.sqrt(3),hslab),
                             center=mp.Vector3(),
                             material=mp.Medium(epsilon=eps))
        geometry = [blk]
        geometry.append(waveguide)

        for j in range(Ny):
            for i in range(Nx+1):
                shift_y = np.sqrt(3)
                geometry.append(mp.Cylinder(r, center=mp.Vector3(i-Nx/2, wgi*np.sqrt(3)/2 + shift_y*j), height = hslab))
                geometry.append(mp.Cylinder(r, center=mp.Vector3(i-Nx/2, -(wgi*np.sqrt(3)/2 + shift_y*j)), height = hslab))

                geometry.append(mp.Cylinder(r, center=mp.Vector3(i-(Nx+1)/2, wgi*np.sqrt(3)/2 + shift_y*(j+1/2)), height = hslab))
                geometry.append(mp.Cylinder(r, center=mp.Vector3(i-(Nx+1)/2, -(wgi*np.sqrt(3)/2 + shift_y*(j+1/2))), height = hslab))
                #geometry.append(mp.Cylinder(r, center=mp.Vector3(i-N/2,-wgi*np.sqrt(3)/2)))

    fcen = 0.3   # pulse center frequency
    df = 0.1    # pulse frequency width
    sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df),
                         component=mp.Hz,
                         center=mp.Vector3(-length/2 +1, 0),
                         size=mp.Vector3(0,wgi*np.sqrt(3), hslab))
              ]

    pml_layers = [mp.PML(1.0)]
    symmetries = [mp.Mirror(mp.Z,+1)]

    sim = mp.Simulation(cell_size=cell,
                        boundary_layers=pml_layers,
                        geometry=geometry,
                        sources=sources,
                        dimensions=3,
                        symmetries=symmetries,
                        resolution=resolution)

    # show geometry (xmingの起動必須)
    sim.init_sim()
    eps_data = sim.get_epsilon()
    
    #from mayavi import mlab
    s = mlab.contour3d(eps_data, colormap="YlGnBu")
    #mlab.show()
    if PhC:
        filename = "geometry_PhC"
    else:
        filename = "geometry_Siwaveguide"
    mlab.savefig(filename + ".png")

    tran_out = mp.FluxRegion(center=mp.Vector3(length/2-3/2, 0, 0),size=mp.Vector3(0, 2*wgi, hslab))
    nfreq = 500 # number of frequencies at which to compute flux
    
    #trans_in = sim.add_flux(fcen, df, nfreq, tran_in)
    trans_out = sim.add_flux(fcen, df, nfreq, tran_out)

    sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Hz, mp.Vector3(checkpoint, 0, 0), 1e-3))

    freqs = mp.get_flux_freqs(trans_out)
    #psd_in = mp.get_fluxes(trans_in)
    psd_out = mp.get_fluxes(trans_out)

    return freqs, psd_out



if __name__ == "__main__":  
    time_start = time.perf_counter()
    dt_now = datetime.datetime.now()
    name_dt = dt_now.strftime('%Y%m%d_%H%M%S')

    a = 400
    c_const = 299792458
    freqs_wo, psd_out_wo = phc_trans(PhC = False, lengthPhC = 20)#, decay_check=10, T_decay=500, name_time=name_dt)
    freqs_w,  psd_out_w  = phc_trans(PhC = True, lengthPhC = 50)#, decay_check=20, T_decay=10000, name_time=name_dt)

    freqs = a / np.array(freqs_w)

    #print(freqs)
    #print(np.array(psd_out_w)/np.array(psd_out_wo))

    df = pd.DataFrame()
    df["normalized_frequency"] = np.array(freqs_w)
    df["wavelength"] = freqs
    df["transmittance"] = np.array(psd_out_w)/np.array(psd_out_wo)
    df.to_csv("transmittance_" + name_dt + ".csv", index=False)

    time_end = time.perf_counter()
    time = time_end - time_start
    print("The necessary time: {:.3f}s".format(time))

    f = open(name_dt + 'totaltime.txt', 'w')
    f.write("The necessary time: {:.3f}s".format(time))
    f.close()