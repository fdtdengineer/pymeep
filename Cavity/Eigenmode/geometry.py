#%%
if True:
    import numpy as np
    
    import matplotlib.pyplot as plt
    #from matplotlib import rc
    #rc('text', usetex=False)
    plt.rcParams['font.family']= 'sans-serif'
    #plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams["font.size"] = 15 # 全体のフォントサイズが変更されます。
    plt.rcParams['xtick.direction'] = 'in'#x軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')
    plt.rcParams['ytick.direction'] = 'in'#y軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')


class GeometryPhC:
    dict_geometry_mp = {
        "circle":["radius"],
        }
    dict_kargs = {}
    dict_geometry = {}

    def get_geometry(self):
        return self.dict_geometry

    def get_kargs(self):
        return self.dict_kargs


class LineDefect(GeometryPhC):
    def __init__(self, a, nx, ny, offset_x, offset_y, n_cavity, len_barrier, wgi, holeshift, radius=0.25):
        # a: lattice constant
        # nx, ny: number of holes in x and y directions
        # offset_x, offset_y: offset of the center of the unit cell
        # len_barrier: barrier length
        # wgi: waveguide width
        # holeshift: hole shift in both sides of holes of the cavity
        # radius: radius of holes (for parser to meep)
        self.a = a
        self.nx = nx
        self.ny = ny
        self.offset_x = offset_x
        self.offset_y = offset_y
        self.n_cavity = n_cavity
        self.len_barrier = len_barrier
        self.wgi = wgi
        self.holeshift = holeshift
        self.radius = radius

        # Build geometry
        self.gen_geometry()


    def gen_geometry(self, holetype="circle",roration=True):
        # holetype: "circle" as default
        # roration: whether apply 90 degree rotatation the geometry or not
        arr_lochole = []
        dx = self.a * np.sqrt(3)/2 # 0.3464
        dy = self.a
        dwgi = (self.wgi-1)*self.a*np.sqrt(3)
        x_center = self.offset_x - self.nx * dx
        y_center = self.offset_y - self.ny * dy / 2

        loc_waveguide_l = self.ny - self.n_cavity - self.len_barrier - 6
        loc_waveguide_r = self.ny + self.n_cavity + self.len_barrier + 6


        for j in range(self.nx+1):
            for i in range(self.ny+1):
                if 2*j == self.nx and abs(2*i - self.ny) <= self.n_cavity - 1: # cavity部分に穴を配置しない
                    continue
                if self.wgi!=0: # waveguide部分に穴を配置しない (wgi=0のときはwaveguideを無視)
                    if 2*j == self.nx - (self.len_barrier + 1) and 2*i <= loc_waveguide_l: #結合導波路 左（バリア領域の長さが奇数のとき）
                        continue
                    elif 2*j == self.nx + (self.len_barrier + 1) and 2*i >= loc_waveguide_r: #結合導波路 右（バリア領域の長さが奇数のとき）
                        continue
                #else:
                if 2*j == self.nx and abs(2*i - self.ny) == self.n_cavity + 1: #cavityの穴シフト
                    y = y_center + dy*i + self.holeshift*self.a*(2*i - self.ny)/abs(2*i - self.ny)
                else:
                    y = y_center + dy*i
                
                if self.wgi!=0:
                    if 2*j < self.nx - (self.len_barrier + 1) and 2*i <= loc_waveguide_l + 1: #結合導波路 左 を (wgi-1)だけ下にシフト
                        x = x_center + 2*dx *j - dwgi
                    elif 2*j > self.nx + (self.len_barrier + 1) and 2*i >= loc_waveguide_r - 1:#結合導波路 右 を (wgi-1)だけ上にシフト
                        x = x_center + 2*dx *j + dwgi
                    else:
                        x = x_center + 2*dx *j
                else:
                    x = x_center + 2*dx *j

                arr_lochole.append([x,y])

        for j in range(self.nx):
            for i in range(self.ny):
                if self.wgi!=0:
                    if 2*j == self.nx - (self.len_barrier + 2) and 2*i <= loc_waveguide_l: #結合導波路 左（バリア領域の長さが偶数のとき）
                        continue
                    elif 2*j == self.nx + self.len_barrier and 2*i >= loc_waveguide_r - 2: #結合導波路 右（バリア領域の長さが偶数のとき）
                        continue
                #else:
                y = y_center + dy/2 + dy*i
                
                if self.wgi!=0:
                    if 2*j < self.nx - (self.len_barrier + 2) and 2*i <= loc_waveguide_l + 1: #結合導波路 左 を (wgi-1)だけ下にシフト
                        x = x_center + dx + 2*dx *j - dwgi
                    elif 2*j > self.nx + self.len_barrier and 2*i >= loc_waveguide_r - 3:
                        x = x_center + dx + 2*dx *j + dwgi
                    else:
                        x = x_center + dx + 2*dx *j
                else:
                    x = x_center + dx + 2*dx *j

                arr_lochole.append([x,y])

        npr_geometry = np.array(arr_lochole)
        
        if roration:
            npr_geometry = np.array([npr_geometry[:,1],-npr_geometry[:,0]]).T
        self.dict_geometry[holetype] = npr_geometry

        # kargs (for parser to meep)
        self.dict_kargs[holetype] = {self.dict_geometry_mp[holetype][0]:self.radius}


    def plot_geometry(self, scale=1./3.):
        dict_geometry = self.get_geometry()
        plt.figure(figsize=(self.ny*scale, self.nx*scale*np.sqrt(3)))
        plt.scatter(self.dict_geometry[dict_geometry.keys[0]][:,0], dict_geometry[self.dict_geometry.keys[0]][:,1])
        plt.xlabel("x")
        plt.ylabel("y")
        plt.show()
        

class WidthModulated(GeometryPhC):
    def __init__(self, a, nx, ny, offset_x, offset_y, len_barrier, wgo, wgi, holeshift, radius=0.25):
        # a: lattice constant
        # self.nx, self.ny: number of holes in x and y directions
        # offset_x, offset_y: offset of the center of the unit cell
        # cavity: cavity type (ex. L3: 3, L5:5)
        # len_barrier: barrier length
        # wgo: outer waveguide width
        # wgi: internal waveguide width
        # holeshift: hole shift in both sides of holes of the cavity
        # radius: radius of holes (for parser to meep)
        self.a = a
        self.nx = nx
        self.ny = ny
        self.offset_x = offset_x
        self.offset_y = offset_y
        self.len_barrier = len_barrier
        self.wgo = wgo
        self.wgi = wgi
        self.holeshift = holeshift
        self.radius = radius

        # Build geometry
        self.gen_geometry()


    def gen_geometry(self, holetype="circle",roration=True):
        # holetype: "circle" as default
        # roration: whether apply 90 degree rotatation the geometry or not
        arr_lochole = []
        dx = self.a * np.sqrt(3)/2
        dy = self.a
        dwgo = (self.wgo-1)*self.a*np.sqrt(3)
        dwgi = (self.wgi-1)*self.a*np.sqrt(3)
        x_center = self.offset_x - self.nx * dx
        y_center = self.offset_y - self.ny * dy / 2


        for j in range(self.nx+1):
            for i in range(self.ny+1):
                if 2*j == self.nx: #中央を空ける
                    continue
                else:
                    y = y_center + dy*i

                    if (abs(2*j - self.nx) == 2 and abs(2*i - self.ny) <= 2): 
                        x = x_center + 2*dx *j + dwgo*(2*j - self.nx)/abs(2*j - self.nx)  + 2*self.holeshift*(2*j - self.nx)/abs(2*j - self.nx)
                    elif (abs(2*j - self.nx) == 2 and abs(2*i - self.ny) == 4):
                        x = x_center + 2*dx *j + dwgo*(2*j - self.nx)/abs(2*j - self.nx)  + 1*self.holeshift*(2*j - self.nx)/abs(2*j - self.nx)
                    elif abs(2*i - self.ny) <= 2*self.len_barrier :
                        x = x_center + 2*dx *j + dwgo*(2*j - self.nx)/abs(2*j - self.nx) 
                    else:
                        x = x_center + 2*dx *j + dwgi*(2*j - self.nx)/abs(2*j - self.nx)

                arr_lochole.append([x,y])

        for j in range(self.nx):
            for i in range(self.ny):
                y = y_center + dy/2 + dy*i

                if abs(2*j - self.nx + 1) == 1 and abs(2*i - self.ny + 1) <= 1:
                    x = x_center +  dx + 2*dx *j + dwgo*(2*j - self.nx + 1)/abs(2*j - self.nx + 1)  + 3*self.holeshift*(2*j - self.nx + 1)/abs(2*j - self.nx + 1)
                elif abs(2*j - self.nx + 1) == 1 and abs(2*i - self.ny + 1) == 3:
                    x = x_center + dx + 2*dx *j + dwgo*(2*j - self.nx + 1)/abs(2*j - self.nx + 1)  + 2*self.holeshift*(2*j - self.nx + 1)/abs(2*j - self.nx + 1)
                elif (abs(2*j - self.nx + 1) == 3 and abs(2*i - self.ny + 1) <= 3) or (abs(2*j - self.nx + 1) == 1 and abs(2*i - self.ny + 1) == 5):
                    x = x_center + dx + 2*dx *j + dwgo*(2*j - self.nx + 1)/abs(2*j - self.nx + 1)  + 1*self.holeshift*(2*j - self.nx + 1)/abs(2*j - self.nx + 1)
                elif abs(2*i - self.ny + 1) <= 2*self.len_barrier:
                    x = x_center + dx + 2*dx *j + dwgo*(2*j - self.nx + 1)/abs(2*j - self.nx + 1) 
                else:
                    x = x_center + dx + 2*dx *j + dwgi*(2*j - self.nx + 1)/abs(2*j - self.nx + 1) 

                arr_lochole.append([x,y])


        npr_geometry = np.array(arr_lochole)
        
        if roration:
            npr_geometry = np.array([npr_geometry[:,1],-npr_geometry[:,0]]).T
        self.dict_geometry[holetype] = npr_geometry

        # kargs (for parser to meep)
        self.dict_kargs[holetype] = {self.dict_geometry_mp[holetype][0]:self.radius}


    def plot_geometry(self, scale=1./3.):
        dict_geometry = self.get_geometry()
        cols = list(dict_geometry.keys())
        plt.figure(figsize=(self.ny*scale, self.nx*scale*np.sqrt(3)))
        plt.scatter(dict_geometry[cols[0]][:,0], dict_geometry[cols[0]][:,1])
        plt.xlabel("x")
        plt.ylabel("y")
        plt.show()


if __name__ == "__main__":
    """
    a = 1 # 0.4
    nx = 8 # 14
    ny = 30 # 50
    offset_x = 0
    offset_y = 0
    n_cavity = 5 # 3
    barrier = 4
    
    wgi = 0 #1.1
    holeshift = 0.2

    ld = LineDefect(a, nx, ny, offset_x, offset_y, n_cavity, barrier, wgi, holeshift)
    ld.plot_geometry()
    """
    a = 1 # 0.4
    nx = 8 # 14
    ny = 30 # 50
    offset_x = 0
    offset_y = 0
    barrier = 6
    wgo = 0.98    
    wgi = 1.1
    holeshift = 0.003 / 0.4

    wm = WidthModulated(a, nx, ny, offset_x, offset_y, barrier, wgo, wgi, holeshift)
    wm.plot_geometry()
    
    print("Done")

# %%
