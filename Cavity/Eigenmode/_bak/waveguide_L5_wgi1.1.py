import math
import insert
insert = insert.insert
import waveguide

flag410 = 1
a = 0.4 + 0.01*flag410
nx = 14
ny = 50
position_x = 100*1 +75*flag410
position_y = 12900*1 -200*flag410
shift_x = 150
shift_y = 0
layer = 1
Nx = 1
cavity = 5
barrier = 4
wgi = 1.1
obj_b1_l = "C:\\Users\\fkh\\CAD\\waveguide_bottom1_left.dwg"
obj_b1_r = "C:\\Users\\fkh\\CAD\\waveguide_bottom1_right.dwg"
obj_t1_l = "C:\\Users\\fkh\\CAD\\waveguide_top1_left.dwg"
obj_t1_r = "C:\\Users\\fkh\\CAD\\waveguide_top1_right.dwg"

for barrier_i in range(1,5):
	position_xi = position_x +2500*(barrier_i - 1)
	waveguide.waveguide(a, nx, ny, position_xi, position_y, shift_x, shift_y, layer, Nx, cavity, barrier_i, wgi, obj_b1_l, obj_b1_r, obj_t1_l, obj_t1_r)