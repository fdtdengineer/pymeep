import math
import insert
insert = insert.insert

a = 0.4
nx = 14
ny = 50
position_x = 100*1 +2500*3
position_y = 12900*1
shift_x = 100
shift_y = 0
layer = 1
Nx = 1
cavity = 3
barrier = 4
wgi = 1.05
obj_b1_l = "C:\\Users\\fkh\\CAD\\waveguide_bottom1_left.dwg"
obj_b1_r = "C:\\Users\\fkh\\CAD\\waveguide_bottom1_right.dwg"
obj_t1_l = "C:\\Users\\fkh\\CAD\\waveguide_top1_left.dwg"
obj_t1_r = "C:\\Users\\fkh\\CAD\\waveguide_top1_right.dwg"

def waveguide(a, nx, ny, position_x, position_y, shift_x, shift_y, layer, Nx, cavity, barrier, wgi, obj_b1_l, obj_b1_r, obj_t1_l, obj_t1_r):
	output = list()

#PhC接続部分
	for i in range(Nx):
		position_xi = position_x + shift_x*(i)
		position_yi = position_y + shift_y*(i)

		x_b1_l = position_xi - a*(barrier + 2 + (wgi-1)*2 )*math.sqrt(3)/2
		y_b1_l = position_yi - a*ny/2
		x_b1_r = position_xi - a*(barrier)*math.sqrt(3)/2
		y_b1_r = position_yi - a*ny/2
		x_t1_l = position_xi + a*(barrier + 2 + (wgi-1)*2 )*math.sqrt(3)/2
		y_t1_l = position_yi + a*ny/2
		x_t1_r = position_xi + a*(barrier)*math.sqrt(3)/2
		y_t1_r = position_yi + a*ny/2

		#PhCに直接繋がる部分
		output.append("-LAYER") 
		output.append("S")
		output.append("33")
		output.append("\n\n")
		output.append("RECTANG")
		output.append("{:.4f},{:.4f}".format(x_b1_l, y_b1_l))
		output.append("{:.4f},{:.4f}".format(x_b1_l - 1, y_b1_l - 10))
		output.append("RECTANG")
		output.append("{:.4f},{:.4f}".format(x_b1_r, y_b1_r))
		output.append("{:.4f},{:.4f}".format(x_b1_r + 1, y_b1_r - 10))
		output.append("RECTANG")
		output.append("{:.4f},{:.4f}".format(x_t1_l, y_t1_l))
		output.append("{:.4f},{:.4f}".format(x_t1_l + 1, y_t1_l + 10))
		output.append("RECTANG")
		output.append("{:.4f},{:.4f}".format(x_t1_r, y_t1_r))
		output.append("{:.4f},{:.4f}".format(x_t1_r - 1, y_t1_r + 10))

		#PhCに直接繋がる部分の両脇
		output.append("-LAYER")
		output.append("S")
		output.append("36")
		output.append("\n\n")
		output.append("RECTANG")
		output.append("{:.4f},{:.4f}".format(x_b1_l - 1, y_b1_l))
		output.append("{:.4f},{:.4f}".format(position_xi - a*(nx+1)*math.sqrt(3)/2 , y_b1_l - 2))
		output.append("RECTANG")
		output.append("{:.4f},{:.4f}".format(x_b1_r + 1, y_b1_l))
		output.append("{:.4f},{:.4f}".format(position_xi + a*(nx+1)*math.sqrt(3)/2 , y_b1_l - 2))
		output.append("RECTANG")
		output.append("{:.4f},{:.4f}".format(x_t1_l + 1, y_t1_l))
		output.append("{:.4f},{:.4f}".format(position_xi + a*(nx+1)*math.sqrt(3)/2 , y_t1_l + 2))
		output.append("RECTANG")
		output.append("{:.4f},{:.4f}".format(x_t1_r - 1, y_t1_l))
		output.append("{:.4f},{:.4f}".format(position_xi - a*(nx+1)*math.sqrt(3)/2 , y_t1_l + 2))

		#曲がったポリライン
		output.append("-insert")
		output.append(obj_b1_l)
		output.append("{:.4f},{:.4f}".format(x_b1_l, y_b1_l - 10))
		output.append("\n\n")
		output.append("-insert")
		output.append(obj_b1_r)
		output.append("{:.4f},{:.4f}".format(x_b1_r, y_b1_r - 10))
		output.append("\n\n")
		output.append("-insert")
		output.append(obj_t1_l)
		output.append("{:.4f},{:.4f}".format(x_t1_l, y_t1_l + 10))
		output.append("\n\n")
		output.append("-insert")
		output.append(obj_t1_r)
		output.append("{:.4f},{:.4f}".format(x_t1_r, y_t1_l + 10))
		output.append("\n\n")

		#曲がったポリライン にくっ付いている長方形
		output.append("-LAYER") 
		output.append("S")
		output.append("34")
		output.append("\n\n")
		output.append("RECTANG")
		output.append("{:.4f},{:.4f}".format(x_b1_l - 1.118, y_b1_l - 25.48))
		output.append("{:.4f},{:.4f}".format(x_b1_l - 3.118, y_b1_l - 30))
		output.append("RECTANG")
		output.append("{:.4f},{:.4f}".format(x_b1_r + 1.118, y_b1_r - 25.48))
		output.append("{:.4f},{:.4f}".format(x_b1_r + 3.118, y_b1_r - 30))
		output.append("RECTANG")
		output.append("{:.4f},{:.4f}".format(x_t1_l + 1.118, y_t1_l + 25.48))
		output.append("{:.4f},{:.4f}".format(x_t1_l + 3.118, y_t1_l + 30))
		output.append("RECTANG")
		output.append("{:.4f},{:.4f}".format(x_t1_r - 1.118, y_t1_r + 25.48))
		output.append("{:.4f},{:.4f}".format(x_t1_r - 3.118, y_t1_r + 30))

		#入力部分との接続
		output.append("PLINE")
		output.append("{:.4f},{:.4f}".format(x_b1_l - 1.118, y_b1_l - 30))
		output.append("{:.4f},{:.4f}".format(x_b1_l - 3.118, y_b1_l - 30))
		output.append("{:.4f},{:.4f}".format(position_xi - 6, y_b1_l - 85.48))
		output.append("{:.4f},{:.4f}".format(position_xi - 4, y_b1_l - 85.48))
		output.append("c")
		output.append("PLINE")
		output.append("{:.4f},{:.4f}".format(x_b1_r + 1.118, y_b1_r - 30))
		output.append("{:.4f},{:.4f}".format(x_b1_r + 3.118, y_b1_r - 30))
		output.append("{:.4f},{:.4f}".format(position_xi + 6, y_b1_r - 85.48))
		output.append("{:.4f},{:.4f}".format(position_xi + 4, y_b1_r - 85.48))
		output.append("c")
		output.append("PLINE")
		output.append("{:.4f},{:.4f}".format(x_t1_l + 1.118, y_t1_l + 30))
		output.append("{:.4f},{:.4f}".format(x_t1_l + 3.118, y_t1_l + 30))
		output.append("{:.4f},{:.4f}".format(position_xi + 6, y_t1_l + 85.48))
		output.append("{:.4f},{:.4f}".format(position_xi + 4, y_t1_l + 85.48))
		output.append("c")
		output.append("PLINE")
		output.append("{:.4f},{:.4f}".format(x_t1_r - 1.118, y_t1_r + 30))
		output.append("{:.4f},{:.4f}".format(x_t1_r - 3.118, y_t1_r + 30))
		output.append("{:.4f},{:.4f}".format(position_xi - 6, y_t1_r + 85.48))
		output.append("{:.4f},{:.4f}".format(position_xi - 4, y_t1_r + 85.48))
		output.append("c")

		output.append("RECTANG")
		output.append("{:.4f},{:.4f}".format(position_xi - 4, y_b1_l - 85.48))
		output.append("{:.4f},{:.4f}".format(position_xi - 6, position_yi - 100))
		output.append("RECTANG")
		output.append("{:.4f},{:.4f}".format(position_xi + 4, y_b1_l - 85.48))
		output.append("{:.4f},{:.4f}".format(position_xi + 6, position_yi - 100))
		output.append("RECTANG")
		output.append("{:.4f},{:.4f}".format(position_xi + 4, y_t1_l + 85.48))
		output.append("{:.4f},{:.4f}".format(position_xi + 6, position_yi + 100))
		output.append("RECTANG")
		output.append("{:.4f},{:.4f}".format(position_xi - 4, y_t1_l + 85.48))
		output.append("{:.4f},{:.4f}".format(position_xi - 6, position_yi + 100))


	lattice_const = str(int(a*1000))
	barrier_str = str(barrier)
	wgi_str = str(wgi)
	cavity_name = str(cavity)
	X_offset_str = str(position_x)
	Y_offset_str = str(position_y)

	filename = "Waveguide_L" + cavity_name  + "_a" + lattice_const + "_c" + barrier_str + "_wgi" + wgi_str + "_spoint[" + X_offset_str + "," + Y_offset_str + "].scr"

	f = open(filename, 'w')	#ファイル書き出し
	for i in range(0,len(output)):
		f.write(output[i] + "\n")
	f.close()
	print("Done.")

#waveguide(a, nx, ny, position_x, position_y, shift_x, shift_y, layer, Nx, cavity, barrier, wgi, obj_b1_l, obj_b1_r, obj_t1_l, obj_t1_r)