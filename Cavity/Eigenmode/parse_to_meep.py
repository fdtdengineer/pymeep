import numpy as np
import meep as mp

dict_geometry_mp = {
    "circle":mp.Cylinder,
    "triangle":mp.Prism,
    "block":mp.Block,
    }

def parse_geometry(obj_phc, thick_slab=0):
    # obj_phc: object of class geometry
    # thick_slab: thickness of PhC slab
    arr_obj = []
    dict_geometry = obj_phc.get_geometry()
    dict_kargs = obj_phc.get_kargs()

    for key, value in dict_geometry.items():
        kargs = dict_kargs[key]
        ptr_mpobj = dict_geometry_mp[key]
        npr_geometry = dict_geometry[key]
        list_obj = [0] * npr_geometry.shape[0]
        
        kargs = dict_kargs[key]
        for i, locs in enumerate(dict_geometry[key]):
            x, y = locs
            list_obj[i] = ptr_mpobj(center=mp.Vector3(x,y,thick_slab), **kargs)

        arr_obj.append(list_obj)

    # flatten
    arr_obj = np.array(arr_obj).flatten().tolist()
    return arr_obj


