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

def main(mass, damp_cons, spring_cons, spring_leng, side_cube, end_time, fram_res, amp, freq, two_sys=False):
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

    height_surf = 2*(0.8+np.max(yc)+side_cube)
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

    line_x1 = line([x_vals[0], block_ratio, yc[0]+block_ratio], [x_vals[-1]+2*dt, block_ratio, yc[0]+block_ratio])
    
    txt_time1 = text('time', 0.05*side_cube, (x_vals[-1]+2*dt, block_ratio, yc[0]+block_ratio), 270)

    line_y1 = line([x_vals[0], block_ratio, max(yc)+block_ratio], [x_vals[0], block_ratio, min(yc)+block_ratio])
    
    txt_U1 = text('U', 0.05*side_cube, (x_vals[0], block_ratio, min(yc)+block_ratio), 270)
  
    plt += [wall, surface, block1, spring1, tube1, cylinder1, redpt1, line_cos1, line_x1, line_y1, txt_U1, txt_time1, redpt_block1]

    # System 2
    if two_sys:
        init_point = yc[0] + side_cube

        block2 = block([0,0, init_point+block_ratio+yc[0]+0.8], side_cube, (128,128,128))

        spring2 = spring([side_cube/4,0,init_point], [side_cube/4,0,yc[0]+init_point+0.8], 0.15*side_cube, 'gray5')

        cylinder2 = cylinder([[-side_cube/4,0,init_point], [-side_cube/4,0,init_point+0.6]], 'gray6', side_cube*0.2)

        tube2 = tube([[-side_cube/4,0,init_point+0.6], [-side_cube/4,0,yc[0]+init_point+0.8]], side_cube*0.15, 'gray5')

        pos_points2 = pos_points(side_cube, x_vals, yc+block_ratio+0.8+init_point)

        line_cos2 = line(pos_points2)

        redpt2 = point(pos_points2[0,:], 'red', redpt_ratio)

        redpt_block2 = point([0,block_ratio,yc[0]+init_point+0.8+block_ratio], 'red', redpt_ratio)

        line_x2 = line([x_vals[0], block_ratio, yc[0]+init_point+0.8+block_ratio], [x_vals[-1]+2*dt, block_ratio, yc[0]+init_point+0.8+block_ratio])
        
        txt_time2 = text('time', 0.05*side_cube, (x_vals[-1]+2*dt, block_ratio, yc[0]+init_point+0.8+block_ratio), 270)

        line_y2 = line([x_vals[0], block_ratio, max(yc)+block_ratio+init_point+0.8], [x_vals[0], block_ratio, min(yc)+block_ratio+init_point+0.8])
        
        txt_U2 = text('U', 0.05*side_cube, (x_vals[0], block_ratio, min(yc)+block_ratio+init_point+0.8), 270)
    
        plt += [block2, spring2, tube2, cylinder2, redpt2, line_cos2, line_x2, line_y2, txt_U2, txt_time2, redpt_block2]

        base_cyli2 = np.array([-side_cube/4,0,init_point])
        top_cyli2 = np.array([-side_cube/4,0,init_point+0.6])
    
    plt.show(azimuth=90, elevation=80, roll=90)
    for i, yc_i in enumerate(yc):
        redpt1.pos(pos_points1[i])
        block1.pos([0,0,yc_i+block_ratio])  
        redpt_block1.pos(0,block_ratio,yc_i+block_ratio)
        spring1.stretch([side_cube/4,0,-0.8], [side_cube/4,0,yc_i])
        tube1.stretch([-side_cube/4,0,-0.2], [-side_cube/4,0,yc_i])

        if two_sys:
            init_point = yc_i + side_cube
            
            redpt2.pos(pos_points2[i])
            block2.pos([0,0,yc_i+block_ratio+0.8+init_point])  
            redpt_block2.pos(0,block_ratio,yc_i+block_ratio+0.8+init_point)
            spring2.stretch([side_cube/4,0,init_point], [side_cube/4,0,yc_i+init_point+0.8])
            tube2.stretch([-side_cube/4,0,init_point+0.6], [-side_cube/4,0,yc_i+init_point+0.8])

            base_cyli2[2] = init_point
            top_cyli2[2]  = init_point+0.6
            pos = (base_cyli2 + top_cyli2) / 2
            cylinder2.pos(pos)

        plt.show()
        if plt.escaped: break # if ESC is hit during the loop
    plt.show(interactive=1).close()

if __name__ == "__main__":
    # System
    mass = 1000                    # Mass                   [kg]
    spring_cons = 48000            # Spring constant        [N/m]
    damp_cons = 4000               # Damping constant  
    
    # Animation model
    spring_leng  = 0.4             # Spring relaxed length  [m]
    side_cube = 1                  # Size cube              [m]

    # Video
    end_time   = 30                # Final time             [s]
    fram_res   = 30                # Frame rate             [fps]

    # Base input
    amp = 0.1                      # Amplitude              [m]
    freq = 1.0                     # Frequency             [Hz]

    main(mass, damp_cons, spring_cons, spring_leng, side_cube, end_time, fram_res, amp, freq, two_sys=True)
        

    # ## Parameters

    # # System
    # m = 1000                       # Mass                          [kg]
    # k = 480000                      # Spring constant               [N/m]
    # c = 400000                       # Damping constant  

    # # Animation model
    # L0  = 0.9                      # Spring relaxed length  [m]
    # side_cube = 1                  # Size cube             [m]

    # # Video
    # tF   = 30                      # Final time                    [s]
    # fR   = 30                      # Frame rate                    [fps]
    # dt   = 1/fR                    # Time resolution               [s]
    # time = np.linspace(0,tF,tF*fR) # Time     

    # # Base input
    # A = 0.5                        # Amplitude                     [m]
    # f = 1.0                        # Frequency                     [Hz]
    # w = 2*np.pi*f                  # Frequency                     [rad/s]
    # u_vet = A*np.cos(w*time)       # Displacement   

    # ## Simulation

    # # Transfer function model
    # G = signal.TransferFunction([c, k], [m, c, k])

    # # Integration
    # t, y, x = signal.lsim(G, u_vet, time)

    # ## Animation

    # # Mass absolute vertical position
    # yc = y + L0 

    # plt = vd.Plotter(interactive=0, axes=1)

    # block_ratio = side_cube/2
    # block = vd.shapes.Cube(pos=[0,0,yc[0]+block_ratio], side=side_cube, c=(128,128,128), alpha=0.7)

    # wall = vd.shapes.Box(pos=(0, 0, -0.8), length=1.5*side_cube, width=1.5*side_cube, height=0.05*side_cube, c=(128,128,128))  # wall

    # height_surf=2*(0.8+yc[0]+side_cube)+A*(0.8+yc[0]+side_cube)
    # surface = vd.shapes.Box(pos=(0, -side_cube, -0.8+(height_surf/2)), length=2*side_cube, width=0.02, height=height_surf, alpha=0)  # surface

    # spring = vd.shapes.Spring([side_cube/4,0,-0.8], [side_cube/4,0,yc[0]], r=side_cube*0.15, thickness=0.01, c='gray5')

    # cylinder = vd.shapes.Cylinder(pos=[[-side_cube/4,0,-0.8], [-side_cube/4,0,-0.2]], c='gray6', r=side_cube*0.2)

    # tube = vd.shapes.Tube(points=[[-side_cube/4,0,-0.2], [-side_cube/4,0,yc[0]]], r=side_cube*0.15, c='gray5')

    # x_vals = map_value(time, side_cube, 2.5*side_cube)
    # points = np.empty((len(x_vals), 3))
    # points[:, 0] = x_vals
    # points[:, 1] = block_ratio
    # points[:, 2] = yc+block_ratio
    # line_cos = vd.shapes.Line(points) #seno

    # redpt = vd.pointcloud.Point(r=5).c('red').pos(points[0,:])

    # redpt_cube = vd.pointcloud.Point(r=5).c('red').pos([0,block_ratio,yc[0]+block_ratio])

    # line_x = vd.shapes.Line((x_vals[0], block_ratio, yc[0]+block_ratio), (x_vals[-1]+2*dt, block_ratio, yc[0]+block_ratio))
    # txt_time = vd.shapes.Text3D(txt='time', s=0.05*side_cube, pos=(x_vals[-1]+2*dt, block_ratio, yc[0]+block_ratio)) 
    # txt_time.rotateX(270, locally=True)

    # line_y = vd.shapes.Line((x_vals[0], block_ratio, max(yc)+block_ratio), (x_vals[0], block_ratio, min(yc)+block_ratio),res=3)
    # txt_U = vd.shapes.Text3D(txt='U', s=0.05*side_cube, pos=(x_vals[0], block_ratio, min(yc)+block_ratio)) 
    # txt_U.rotateX(270, locally=True)

    # plt += [wall, surface, block, spring, tube, cylinder, redpt, line_cos, line_x, line_y, txt_U, txt_time, redpt_cube]

    # ##
    # init_point = yc[0] + side_cube
    # block2 = vd.shapes.Cube(pos=[0,0, init_point+block_ratio+yc[0]+0.8], side=side_cube, c=(128,128,128), alpha=0.7)

    # spring2 = vd.shapes.Spring([side_cube/4,0, init_point], [side_cube/4,0,yc[0]+init_point+0.8], r=side_cube*0.15, thickness=0.01, c='gray5')

    # cylinder2 = vd.shapes.Cylinder(pos=[[-side_cube/4,0,init_point], [-side_cube/4,0,init_point+0.6]], c='gray6', r=side_cube*0.2)

    # tube2 = vd.shapes.Tube(points=[[-side_cube/4,0,init_point+0.6], [-side_cube/4,0,yc[0]+init_point+0.8]], r=side_cube*0.15, c='gray5')

    # points2 = np.empty((len(x_vals), 3))
    # points2[:, 0] = x_vals
    # points2[:, 1] = block_ratio 
    # points2[:, 2] = yc+block_ratio+0.8+init_point
    # line_cos2 = vd.shapes.Line(points2) #seno

    # redpt2 = vd.pointcloud.Point(r=5).c('red').pos(points2[0,:])

    # redpt_block2 = vd.pointcloud.Point(r=5).c('red').pos([0,block_ratio,yc[0]+init_point+0.8+block_ratio])

    # line_x2 = vd.shapes.Line((x_vals[0], block_ratio, yc[0]+init_point+0.8+block_ratio), (x_vals[-1]+2*dt, block_ratio, yc[0]+init_point+0.8+block_ratio))
    # txt_time2 = vd.shapes.Text3D(txt='time', s=0.05*side_cube, pos=(x_vals[-1]+2*dt, block_ratio, yc[0]+init_point+0.8+block_ratio)) 
    # txt_time2.rotateX(270, locally=True)

    # line_y2 = vd.shapes.Line((x_vals[0], block_ratio, max(yc)+block_ratio+init_point+0.8), (x_vals[0], block_ratio, min(yc)+block_ratio+init_point+0.8),res=3)
    # txt_U2 = vd.shapes.Text3D(txt='U', s=0.05*side_cube, pos=(x_vals[0], block_ratio, min(yc)+block_ratio+init_point+0.8)) 
    # txt_U2.rotateX(270, locally=True)


    # plt += [block2,spring2,cylinder2,tube2,redpt2, line_cos2, line_x2, line_y2, txt_U2, txt_time2, redpt_block2]
    # plt.show(azimuth=90, elevation=80, roll=90)

    # ##
    # base = np.array([-side_cube/4,0,init_point])
    # top = np.array([-side_cube/4,0,init_point+0.6])
    # for i, yc_i in enumerate(yc):

    #     block.pos([0,0,yc_i+block_ratio])  
    #     redpt_cube.pos(0,block_ratio,yc_i+block_ratio)
    #     redpt.pos(points[i])
    #     spring.stretch([side_cube/4,0,-0.8], [side_cube/4,0,yc_i])
    #     tube.stretch([-side_cube/4,0,-0.2], [-side_cube/4,0,yc_i])

    #     init_point = yc_i + side_cube

    #     block2.pos([0,0,yc_i+block_ratio+0.8+init_point])  
    #     redpt_block2.pos(0,block_ratio,yc_i+block_ratio+0.8+init_point)
    #     redpt2.pos(points2[i])
    #     spring2.stretch([side_cube/4,0,init_point], [side_cube/4,0,yc_i+init_point+0.8])
    #     tube2.stretch([-side_cube/4,0,init_point+0.6], [-side_cube/4,0,yc_i+init_point+0.8])

    #     base[2] = init_point
    #     top[2]  = init_point+0.6
    #     pos = (base + top) / 2

    #     cylinder2.pos(pos)

    #     plt.show()
    #     if plt.escaped: break # if ESC is hit during the loop

    # plt.show(interactive=1).close()