import sys
import numpy as np
from scipy import signal
import vedo as vd

def map_value(x, out_min, out_max):
    return (x - min(x)) * (out_max - out_min) / (max(x) - min(x)) + out_min

def simulation(mass, damp_cons, spring_cons, u_vet, time):
    # Transfer function model
    G = signal.TransferFunction([damp_cons, spring_cons], [mass, damp_cons, spring_cons])

    # Integration
    t, y, x = signal.lsim(G, u_vet, time)
    return y

def block(pos, side, color):
    return vd.shapes.Cube(pos=pos, side=side, c=color)

def spring(start_point, end_point, ratio, color):
    return vd.shapes.Spring(start_point, end_point, r=ratio, thickness=0.01, c=color)

def cylinder(pos, color, ratio):
    return vd.shapes.Cylinder(pos=pos, c=color, r=ratio)

def tube(points, ratio, color):
    return vd.shapes.Tube(points=points, r=ratio, c=color)

def pos_points(side_cube, x_vals, z):
    
    points = np.empty((len(x_vals), 3))
    points[:, 0] = x_vals
    block_ratio = side_cube/2
    points[:, 1] = block_ratio
    points[:, 2] = z
    return points

def line(points1, points2=None):
    return vd.shapes.Line(points1, points2)

def point(pos, color, ratio):
    return vd.pointcloud.Point(r=ratio).c(color).pos(pos)

def text(txt, size, pos, angle):
    txt_3d = vd.shapes.Text3D(txt=txt, s=size, pos=pos) 
    txt_3d.rotateX(angle, locally=True)
    return txt_3d

def mass_spring_damper(mass, damp_cons, spring_cons, spring_leng, side_cube, end_time, fram_res, amp, freq, two_sys=False, mass2=None, damp_cons2=None, spring_cons2=None, spring_leng2=None, side_cube2=None, amp2=None, freq2=None):
    dt   = 1/fram_res                                # Time resolution        [s]
    time = np.linspace(0,end_time,end_time*fram_res) # Time 

    w_freq = 2*np.pi*freq                 # Frequency        [rad/s]
    u_vet = amp*np.cos(w_freq*time)       # Displacement  

    yc = simulation(mass, damp_cons, spring_cons, u_vet, time)
    yc = yc + spring_leng  

    plt = vd.Plotter(interactive=0, axes=4)
  
    block_ratio = side_cube/2

    redpt_ratio = 5*side_cube

    # System 1
    wall = vd.shapes.Box(pos=(0, 0, -0.8), length=1.5*side_cube, width=1.5*side_cube, height=0.05*side_cube, c=(128,128,128))  
    const = 2 if two_sys else 1
    height_surf = const*(0.8+np.max(yc)+side_cube)
    surface = vd.shapes.Box(pos=(0, -side_cube, -0.8+(height_surf/2)), length=2*side_cube, width=0.02, height=height_surf, alpha=0)

    block1 = block([0,0,yc[0]+block_ratio], side_cube, (128,128,128))

    spring1 = spring([side_cube/4,0,-0.8], [side_cube/4,0,yc[0]], 0.15*side_cube, 'gray5')

    cylinder1 = cylinder([[-side_cube/4,0,-0.8], [-side_cube/4,0,-0.2]], 'gray6', side_cube*0.2)

    tube1 = tube([[-side_cube/4,0,-0.2], [-side_cube/4,0,yc[0]]], side_cube*0.15, 'gray5')

    x_vals = map_value(time, side_cube, 2.5*side_cube)

    pos_points1 = pos_points(side_cube, x_vals, yc+block_ratio)

    line_cos1 = line(pos_points1)

    redpt1 = point(pos_points1[0,:], 'red', redpt_ratio) 

    redpt_block1 = point([0,block_ratio,yc[0]+block_ratio], 'red', redpt_ratio)

    line_x1 = line([x_vals[0], block_ratio, np.median(yc)+block_ratio], [x_vals[-1]+2*dt, block_ratio, np.median(yc)+block_ratio])
    
    txt_time1 = text('time', 0.05*side_cube, (x_vals[-1]+2*dt, block_ratio, np.median(yc)+block_ratio), 270)

    line_y1 = line([x_vals[0], block_ratio, max(yc)+block_ratio], [x_vals[0], block_ratio, min(yc)+block_ratio])
    
    txt_U1 = text('U', 0.05*side_cube, (x_vals[0], block_ratio, min(yc)+block_ratio), 270)
  
    plt += [wall, surface, block1, spring1, tube1, cylinder1, redpt1, line_cos1, line_x1, line_y1, txt_U1, txt_time1, redpt_block1]

    # System 2
    if two_sys:

        if amp2 is None:
            amp2 = amp
        
        if freq2 is None:
            freq2 = freq
        
        if (amp2 is not None) or (freq2 is not None):
            w_freq = 2*np.pi*freq2                  # Frequency        [rad/s]
            u_vet2 = amp2*np.cos(w_freq*time)       # Displacement 
        else:
            u_vet2 = u_vet

        if mass2 is None:
            mass2 = mass

        if damp_cons2 is None:
            damp_cons2 = damp_cons  

        if spring_cons2 is None:
            spring_cons2 = spring_cons
        
        if (mass2 is not None) or (damp_cons2 is not None) or (spring_cons2 is not None):
            yc2 = simulation(mass2, damp_cons2, spring_cons2, u_vet2, time)
        else:
            yc2 = yc       

        if (spring_leng2 is not None):
            yc2 = yc2 + spring_leng2
        else:
            yc2 = yc2 + spring_leng

        if (side_cube2 is None):
            side_cube2 = side_cube
            block_ratio2 = block_ratio
            redpt_ratio2 = redpt_ratio
        else:
            if side_cube2 < side_cube:
                sys.exit('side_cube2 < side_cube')
            block_ratio2 = side_cube2/2
            redpt_ratio2 = 5

        init_point = yc[0] + side_cube

        block2 = block([0,0, init_point+block_ratio2+yc2[0]+0.8], side_cube2, (128,128,128))

        spring2 = spring([side_cube/4,0,init_point], [side_cube/4,0,yc2[0]+init_point+0.8], 0.15*side_cube, 'gray5')

        cylinder2 = cylinder([[-side_cube/4,0,init_point], [-side_cube/4,0,init_point+0.6]], 'gray6', side_cube*0.2)

        tube2 = tube([[-side_cube/4,0,init_point+0.6], [-side_cube/4,0,yc2[0]+init_point+0.8]], side_cube*0.15, 'gray5')

        x_vals2 = map_value(time, side_cube2, 2.5*side_cube2)

        pos_points2 = pos_points(side_cube2, x_vals2, yc2+block_ratio2+0.8+init_point)

        line_cos2 = line(pos_points2)

        redpt2 = point(pos_points2[0,:], 'red', redpt_ratio2)

        redpt_block2 = point([0,block_ratio2,yc2[0]+init_point+0.8+block_ratio2], 'red', redpt_ratio2)

        line_x2 = line([x_vals2[0], block_ratio2, np.median(yc2)+init_point+0.8+block_ratio2], [x_vals2[-1]+2*dt, block_ratio2, np.median(yc2)+init_point+0.8+block_ratio2])
        
        txt_time2 = text('time', 0.05*side_cube2, (x_vals2[-1]+2*dt, block_ratio2, np.median(yc2)+init_point+0.8+block_ratio2), 270)

        line_y2 = line([x_vals2[0], block_ratio2, max(yc2)+block_ratio2+init_point+0.8], [x_vals2[0], block_ratio2, min(yc2)+block_ratio2+init_point+0.8])
        
        txt_U2 = text('U', 0.05*side_cube2, (x_vals2[0], block_ratio2, min(yc2)+block_ratio2+init_point+0.8), 270)
    
        plt += [block2, spring2, tube2, cylinder2, redpt2, line_cos2, line_x2, line_y2, txt_U2, txt_time2, redpt_block2]

        base_cyli2 = np.array([-side_cube/4,0,init_point])
        top_cyli2 = np.array([-side_cube/4,0,init_point+0.6])

        yc_max = np.maximum(np.max(yc2), np.max(yc))
        side_cube_max = np.maximum(side_cube, side_cube2)
    else:
        yc_max = np.max(yc)
        side_cube_max = side_cube

    wall = vd.shapes.Box(pos=(0, 0, -0.8), length=1.5*side_cube, width=1.5*side_cube, height=0.05*side_cube, c=(128,128,128))  
    
    const = 2 if two_sys else 1
    height_surf = const*(0.8+yc_max+side_cube_max)
    surface = vd.shapes.Box(pos=(0, -side_cube_max, -0.8+(height_surf/2)), length=2*side_cube_max, width=0.02, height=height_surf, alpha=0)
    
    plt += [wall, surface]
    
    plt.show(azimuth=90, elevation=80, roll=90)
    for i, yc_i in enumerate(yc):
        redpt1.pos(pos_points1[i])
        block1.pos([0,0,yc_i+block_ratio])  
        redpt_block1.pos(0,block_ratio,yc_i+block_ratio)
        spring1.stretch([side_cube/4,0,-0.8],[side_cube/4,0,yc_i])
        tube1.stretch([-side_cube/4,0,-0.2], [-side_cube/4,0,yc_i])

        if two_sys:
            init_point = yc_i + side_cube
            yc2_i = yc2[i]
            
            redpt2.pos(pos_points2[i])
            block2.pos([0,0,yc2_i+block_ratio2+0.8+init_point])  
            redpt_block2.pos(0,block_ratio2,yc2_i+block_ratio2+0.8+init_point)
            spring2.stretch([side_cube/4,0,init_point], [side_cube/4,0,yc2_i+init_point+0.8])
            tube2.stretch([-side_cube/4,0,init_point+0.6], [-side_cube/4,0,yc2_i+init_point+0.8])

            base_cyli2[2] = init_point
            top_cyli2[2]  = init_point+0.6
            pos = (base_cyli2 + top_cyli2) / 2
            cylinder2.pos(pos)

        plt.show()
        if plt.escaped: break # if ESC is hit during the loop
    plt.show(interactive=1).close()

if __name__ == "__main__":

    # Video
    end_time   = 30                # Final time             [s]
    fram_res   = 30                # Frame rate             [fps]

    # System
    mass = 1000                    # Mass                   [kg]
    spring_cons = 48000            # Spring constant        [N/m]
    damp_cons = 4000               # Damping constant  
    
    # Animation model
    spring_leng  = 0.9             # Spring relaxed length  [m]
    side_cube = 1                  # Size cube              [m]

    # Base input
    amp = 0.2                      # Amplitude              [m]
    freq = 1.0                     # Frequency             [Hz]

    two_sys = True

    """If None use the same values as system 1"""
    mass2 = None                   # Mass                 [kg]
    spring_cons2 = None            # Spring constant        [N/m]
    damp_cons2 = None              # Damping constant  
    
    # Animation model
    spring_leng2  = None           # Spring relaxed length  [m]
    side_cube2 = 2                 # Size cube              [m]

    # Base input
    amp2 = 0.5                     # Amplitude              [m]
    freq2 = None                   # Frequency             [Hz]

    mass_spring_damper(mass, damp_cons, spring_cons, spring_leng, side_cube, end_time, fram_res, amp, freq, two_sys, mass2, spring_cons2, damp_cons2, spring_leng2, side_cube2, amp2, freq2)