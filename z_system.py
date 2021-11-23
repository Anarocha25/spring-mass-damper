import numpy as np
from scipy import signal
import vedo as vd

class MassSpringDamper():
    def __init__(self, mass, damp_cons, spring_cons, spring_leng, side_cube, end_time, fram_res, amp, freq):
        self.mass = mass
        self.damp_cons = damp_cons
        self.spring_cons = spring_cons
        
        self.spring_len = spring_leng 
        self.side_cube = side_cube

        self.dt   = 1/fram_res                                # Time resolution        [s]
        self.time = np.linspace(0, end_time, end_time * fram_res) # Time 

        w_freq = 2 * np.pi * freq                           # Frequency        [rad/s]
        self.u_vet = amp * np.cos(w_freq * self.time)       # Displacement  

        self.yc = self.__simulation()
        self.yc += spring_leng

        self.block_ratio = self.side_cube/2

    def __simulation(self):
        # Transfer function model
        G = signal.TransferFunction([self.damp_cons, self.spring_cons], [self.mass, self.damp_cons, self.spring_cons])

        # Integration
        t, y, x = signal.lsim(G, self.u_vet, self.time)
        return y

    def map_value(self, x, out_min, out_max):
        return (x - min(x)) * (out_max - out_min) / (max(x) - min(x)) + out_min

    def block(self, pos, side, color):
        return vd.shapes.Cube(pos=pos, side=side, c=color)

    def spring(self, start_point, end_point, ratio, color):
        return vd.shapes.Spring(start_point, end_point, r=ratio, thickness=0.01, c=color)

    def cylinder(self, pos, color, ratio):
        return vd.shapes.Cylinder(pos=pos, c=color, r=ratio)

    def tube(self, points, ratio, color):
        return vd.shapes.Tube(points=points, r=ratio, c=color)

    def pos_points(self, side_cube, x_vals, z):
        points = np.empty((len(x_vals), 3))
        points[:, 0] = x_vals
        block_ratio = side_cube/2
        points[:, 1] = block_ratio
        points[:, 2] = z
        return points

    def line(self, points1, points2=None):
        return vd.shapes.Line(points1, points2)

    def point(self, pos, color, ratio):
        return vd.pointcloud.Point(r=ratio).c(color).pos(pos)

    def text(self, txt, size, pos, angle):
        txt_3d = vd.shapes.Text3D(txt=txt, s=size, pos=pos) 
        txt_3d.rotateX(angle, locally=True)
        return txt_3d

class MassSpringDamperSys1(MassSpringDamper):
        def __init__(self, mass, damp_cons, spring_cons, spring_leng, side_cube, end_time, fram_res, amp, freq):
            super().__init__(mass, damp_cons, spring_cons, spring_leng, side_cube, end_time, fram_res, amp, freq)

            self.redpt_ratio = 5 * self.side_cube

            self.block = self.block([0, 0, self.yc[0]+self.block_ratio], self.side_cube, (128,128,128))

            self.spring = self.spring([self.side_cube/4, 0, -0.8], [self.side_cube/4, 0, self.yc[0]], 0.15 * self.side_cube, 'gray5')

            self.cylinder = self.cylinder([[-self.side_cube/4, 0, -0.8], [-self.side_cube/4, 0, -0.2]], 'gray6', self.side_cube * 0.2)

            self.tube = self.tube([[-self.side_cube/4, 0, -0.2], [-self.side_cube/4, 0, self.yc[0]]], self.side_cube * 0.15, 'gray5')

            x_vals = self.map_value(self.time, self.side_cube, 2.5 * self.side_cube)

            self.pos_points = self.pos_points(side_cube, x_vals, self.yc + self.block_ratio)

            self.line_cos = self.line(self.pos_points)

            self.redpt = self.point(self.pos_points[0,:], 'red', self.redpt_ratio) 

            self.redpt_block = self.point([0, self.block_ratio, self.yc[0] + self.block_ratio], 'red', self.redpt_ratio)

            self.line_x = self.line([x_vals[0], self.block_ratio, np.median(self.yc) + self.block_ratio], [x_vals[-1] + 2 * self.dt, self.block_ratio, np.median(self.yc) + self.block_ratio])
            
            self.txt_time = self.text('time', 0.05*side_cube, (x_vals[-1]+2 * self.dt, self.block_ratio, np.median(self.yc) + self.block_ratio), 270)

            self.line_y = self.line([x_vals[0], self.block_ratio, max(self.yc) + self.block_ratio], [x_vals[0], self.block_ratio, min(self.yc) + self.block_ratio])
            
            self.txt_U = self.text('U', 0.05 * self.side_cube, (x_vals[0], self.block_ratio, min(self.yc) + self.block_ratio), 270)

class MassSpringDamperSys2(MassSpringDamper):
        def __init__(self, yc, mass, damp_cons, spring_cons, spring_leng, side_cube, end_time, fram_res, amp, freq):
            super().__init__(mass, damp_cons, spring_cons, spring_leng, side_cube, end_time, fram_res, amp, freq)

            self.init_point = yc[0] + self.side_cube

            self.redpt_ratio = 5

            self.block = self.block([0, 0,self.init_point + self.block_ratio + self.yc[0] + 0.8], self.side_cube, (128,128,128))

            self.spring = self.spring([side_cube/4, 0,self.init_point], [side_cube/4, 0, self.yc[0] +self.init_point + 0.8], 0.15 * side_cube, 'gray5')

            self.cylinder = self.cylinder([[-side_cube/4, 0,self.init_point], [-side_cube/4, 0,self.init_point + 0.6]], 'gray6', side_cube * 0.2)

            self.tube = self.tube([[-side_cube/4, 0,self.init_point + 0.6], [-side_cube/4, 0, self.yc[0] +self.init_point + 0.8]], side_cube * 0.15, 'gray5')

            x_vals2 = self.map_value(self.time, self.side_cube, 2.5 * self.side_cube)

            self.pos_points = self.pos_points(self.side_cube, x_vals2, self.yc + self.block_ratio + 0.8 +self.init_point)

            self.line_cos = self.line(self.pos_points)

            self.redpt = self.point(self.pos_points[0,:], 'red', self.redpt_ratio)

            self.redpt_block = self.point([0, self.block_ratio, self.yc[0] + self.init_point + 0.8 + self.block_ratio], 'red', self.redpt_ratio)

            self.line_x = self.line([x_vals2[0], self.block_ratio, np.median(self.yc) +self.init_point + 0.8 + self.block_ratio], \
                           [x_vals2[-1] + 2 * self.dt, self.block_ratio, np.median(self.yc) + self.init_point + 0.8 + self.block_ratio])
            
            self.txt_time = self.text('time', 0.05 * self.side_cube, (x_vals2[-1] + 2 * self.dt, self.block_ratio, \
                                    np.median(self.yc) + self.init_point + 0.8 + self.block_ratio), 270)

            self.line_y = self.line([x_vals2[0], self.block_ratio, max(self.yc)+self.block_ratio + self.init_point+0.8], \
                                [x_vals2[0], self.block_ratio, min(self.yc) + self.block_ratio + self.init_point + 0.8])
            
            self.txt_U = self.text('U', 0.05 * self.side_cube, (x_vals2[0], self.block_ratio, min(self.yc) + self.block_ratio + self.init_point + 0.8), 270)
