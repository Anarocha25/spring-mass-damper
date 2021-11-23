import vedo as vd  
import cv2  
import os
import ast
import shutil
import numpy as np 
from z_system import MassSpringDamperSys1, MassSpringDamperSys2

folder_name = 'temp'
directory = os.path.join(os.path.dirname(__file__), folder_name)
dir_file = os.path.join(directory, 'param1_file.txt')
# READ
file = open(dir_file, "r")
contents = file.read()
param1 = ast.literal_eval(contents)
file.close()   

if param1["two_sys"]:
    dir_file = os.path.join(directory, 'param2_file.txt')
    # READ
    file = open(dir_file, "r")
    contents = file.read()
    param2 = ast.literal_eval(contents)
    file.close()  

vp = vd.Plotter(axes=4, offscreen=True)
vp.show(azimuth=90, elevation=80, roll=90, resetcam=True)

# Create mass spring damper system
sys1 = MassSpringDamperSys1(param1["mass_par"], param1["damp_cons_par"], param1["spring_cons_par"], \
                            param1["spring_leng_par"], param1["side_cube_par"], param1["end_time_par"], \
                            param1["fram_res_par"], param1["amp_par"], param1["freq_par"])

yc_max = np.max(sys1.yc)
side_cube_max = sys1.side_cube

if param1["two_sys"]:
    sys2 = MassSpringDamperSys2(sys1.yc, param2["mass_par2"], param2["damp_cons_par2"], param2["spring_cons_par2"], \
                            param2["spring_leng_par2"], param2["side_cube_par2"], param2["end_time_par"], \
                            param2["fram_res_par"], param2["amp_par2"], param2["freq_par2"])
    
    base_cyli2 = np.array([-sys1.side_cube/4, 0, sys2.init_point])
    top_cyli2 = np.array([-sys1.side_cube/4, 0, sys2.init_point + 0.6])
    
    yc_max = np.maximum(np.max(sys2.yc), np.max(sys1.yc))
    side_cube_max = np.maximum(sys2.side_cube, sys1.side_cube)
    
wall = vd.shapes.Box(pos=(0, 0, -0.8), length=1.5 * sys1.side_cube, width=1.5 * sys1.side_cube, \
                    height=0.05 * sys1.side_cube, c=(128,128,128))  

const = 2 if param1["two_sys"] else 1

height_surf = const*(0.8 + yc_max + side_cube_max)

surface = vd.shapes.Box(pos=(0, -side_cube_max, -0.8 + (height_surf/2)), length=2 * side_cube_max, width=0.02,
                        height=height_surf, alpha=0)

vp += [wall, surface, sys1.block, sys1.spring, sys1.tube, sys1.cylinder, sys1.redpt,
            sys1.line_cos, sys1.line_x, sys1.line_y, sys1.txt_U, sys1.txt_time, sys1.redpt_block]

vp += [wall, surface, sys1.block, sys1.spring, sys1.tube]

if param1["two_sys"]:
    vp += [sys2.block, sys2.spring, sys2.tube, sys2.cylinder, sys2.redpt,
            sys2.line_cos, sys2.line_x, sys2.line_y, sys2.txt_U, sys2.txt_time, sys2.redpt_block]

vp.resetCamera()
video = vd.io.Video(param1["name_save"]+".mp4", fps=param1["fram_res_par"], backend='opencv') # backend='opencv'
video.addFrame()

for i, yc_i in enumerate(sys1.yc):
    vp.actors[6].pos(sys1.pos_points[i]) #redpoint
    vp.actors[2].pos([0, 0, yc_i + sys1.block_ratio])  #block
    vp.actors[12].pos(0, sys1.block_ratio, yc_i + sys1.block_ratio) #redpt_block
    vp.actors[3].stretch([sys1.side_cube/4, 0,-0.8],[sys1.side_cube/4, 0, yc_i]) #spring
    vp.actors[4].stretch([-sys1.side_cube/4, 0, -0.2], [-sys1.side_cube/4, 0, yc_i]) #tube
    
    if param1["two_sys"]:
        init_point = yc_i + sys1.side_cube
        yc2_i = sys2.yc[i]
        
        vp.actors[17].pos(sys2.pos_points[i]) #redpt2
        vp.actors[13].pos([0, 0, yc2_i + sys2.block_ratio + 0.8 + init_point]) #block
        vp.actors[23].pos(0, sys2.block_ratio, yc2_i + sys2.block_ratio + 0.8 + init_point) #redpt_block2
        vp.actors[14].stretch([sys1.side_cube/4, 0,init_point], [sys1.side_cube/4, 0, yc2_i + init_point + 0.8]) #spring
        vp.actors[15].stretch([-sys1.side_cube/4, 0,init_point + 0.6], [-sys1.side_cube/4, 0, yc2_i + init_point + 0.8]) #tube

        base_cyli2[2] = init_point
        top_cyli2[2]  = init_point + 0.6
        pos = (base_cyli2 + top_cyli2)/2
        vp.actors[16].pos(pos) #cylinder
            
    video.addFrame() 
 
video.close()

if os.path.exists(directory):
    shutil.rmtree(directory)