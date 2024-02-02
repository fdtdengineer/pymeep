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


class MeepTransmittance:
    def __init__(
            self, 
            sim=[],
            endpoint = mp.Vector3(),
            tran_out = [],
            dim=2,
            area_z=8, 
            thick_slab=0,
            fcen=0.25,
            df=0.05,
            thres_conv = 1e-3,
            nfreq = 500,
            resolution=16,
            dpml=1
            ) -> None:
        
        self.sim = sim
        self.endpoint = endpoint
        self.fcen = fcen
        self.df = df
        self.dict_thick_slab = {2:mp.inf, 3:thick_slab}
        self.dict_area_z = {2:0, 3:area_z} 
        self.dict_beamz = {2:0, 3:thick_slab}
        self.area_z = self.dict_area_z[dim]
        self.beamz = self.dict_beamz[dim]
        self.dim = dim
        self.thick_slab = self.dict_thick_slab[dim]
        self.thres_conv = thres_conv
        self.nfreq = nfreq
        self.resolution = resolution
        self.dpml = dpml
        pass

    # get the transmittance of the reference wabveguide
    def get_reference_transmittance(
            self, 
            a=1, 
            nx=5, 
            ny=30, 
            len_siwg=10, 
            eps_r=2.6**2,
            wg=1, 
            pml_buffer=2,
            ):
        cell = mp.Vector3(a*(ny+pml_buffer+len_siwg*2+self.dpml-1), a*(nx+pml_buffer/2)*np.sqrt(3), self.area_z)

        # Waveguide
        arr_wg = [
            mp.Block(
                mp.Vector3(a*(ny+2*(len_siwg+pml_buffer+self.dpml-1)), wg*np.sqrt(3), self.thick_slab),
                center=mp.Vector3(0, 0),
                material=mp.Medium(epsilon=eps_r)
            ),
        ]

        # All Geometry
        arr_geometry = arr_wg

        # Source
        src = [
            mp.Source(
                mp.GaussianSource(self.fcen, fwidth=self.df),
                component=mp.Hz,
                center=mp.Vector3(-a*(ny/2+len_siwg-(self.dpml-1)/2), 0),
                size=mp.Vector3(0, wg*np.sqrt(3),self.beamz)
            )
        ]

        # Symmetry
        sym = [mp.Mirror(mp.Y, phase=-1)]

        if self.dim==3:
            sym.append(mp.Mirror(mp.Z, phase=1))

        # PML
        pml_layers = [mp.PML(self.dpml)]

        # Build simulation object
        self.sim_ref = mp.Simulation(
            cell_size=cell,
            geometry=arr_geometry,
            boundary_layers=pml_layers,
            sources=src,
            dimensions=self.dim,
            symmetries=sym,
            resolution=self.resolution,
        )

        # Output flux
        vec_out = mp.Vector3(a*(ny/2+len_siwg-1+(self.dpml-2)/2), 0)
        tran_out = mp.FluxRegion(center=vec_out,size=mp.Vector3(0, 2*wg,self.beamz))
        trans_out = self.sim_ref.add_flux(self.fcen, self.df, self.nfreq, tran_out)

        # Plot geometry
        if self.dim==2:
            f = plt.figure(dpi=100)
            self.sim_ref.plot2D(ax=f.gca())
            plt.show()

        # Run simulation
        self.sim_ref.run(
            until_after_sources=mp.stop_when_fields_decayed(
                50, mp.Hz, vec_out, self.thres_conv
            )
        )

        self.freqs_ref = np.array(mp.get_flux_freqs(trans_out))
        self.psd_out = np.array(mp.get_fluxes(trans_out))

        return self.freqs_ref, self.psd_out

    def plot_reference_transmittance(self):
        plt.figure(figsize=(6,4))
        plt.plot(self.freqs_ref, self.psd_out, label="Reference")
        plt.xlabel("Frequency")
        plt.ylabel("Power Spectral Density")
        plt.yscale("log")
        plt.legend()
        plt.tight_layout()
        plt.show()




if __name__ == "__main__":
    cls_ref = MeepTransmittance(dpml=1)
    cls_ref.get_reference_transmittance()
    cls_ref.plot_reference_transmittance()


# %%
