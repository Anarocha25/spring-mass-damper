import numpy as np
from scipy import signal
import vedo as vd

import sys
import ast
from PyQt5 import Qt, QtWidgets, QtCore
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

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

        self.redpt_ratio = 5 * self.side_cube

        self.block = self.block([0, 0, self.yc[0]+self.block_ratio], self.side_cube, (128,128,128))

        self.spring = self.spring([self.side_cube/4, 0, -0.8], [self.side_cube/4, 0, self.yc[0]], 0.15 * self.side_cube, 'gray5')

        self.cylinder = self.cylinder([[-self.side_cube/4, 0, -0.8], [-self.side_cube/4, 0, -0.2]], 'gray6', self.side_cube * 0.2)

        self.tube = self.tube([[-self.side_cube/4, 0, -0.2], [-self.side_cube/4, 0, self.yc[0]]], self.side_cube * 0.15, 'gray5')

        x_vals = self.__map_value(self.time, self.side_cube, 2.5 * self.side_cube)

        self.pos_points = self.pos_points(side_cube, x_vals, self.yc + self.block_ratio)

        self.line_cos = self.line(self.pos_points)

        self.redpt = self.point(self.pos_points[0,:], 'red', self.redpt_ratio) 

        self.redpt_block = self.point([0, self.block_ratio, self.yc[0] + self.block_ratio], 'red', self.redpt_ratio)

        self.line_x = self.line([x_vals[0], self.block_ratio, np.median(self.yc) + self.block_ratio], [x_vals[-1] + 2 * self.dt, self.block_ratio, np.median(self.yc) + self.block_ratio])
        
        self.txt_time = self.text('time', 0.05*side_cube, (x_vals[-1]+2 * self.dt, self.block_ratio, np.median(self.yc) + self.block_ratio), 270)

        self.line_y = self.line([x_vals[0], self.block_ratio, max(self.yc) + self.block_ratio], [x_vals[0], self.block_ratio, min(self.yc) + self.block_ratio])
        
        self.txt_U = self.text('U', 0.05 * self.side_cube, (x_vals[0], self.block_ratio, min(self.yc) + self.block_ratio), 270)

    def __simulation(self):
        # Transfer function model
        G = signal.TransferFunction([self.damp_cons, self.spring_cons], [self.mass, self.damp_cons, self.spring_cons])

        # Integration
        t, y, x = signal.lsim(G, self.u_vet, self.time)
        return y

    def __map_value(self, x, out_min, out_max):
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

class PopUp(QtWidgets.QDialog):
    def __init__(self, warnings):
        super().__init__()

        self.warnings = warnings

        self.setWindowTitle("Warning")
        self.layout = QtWidgets.QVBoxLayout()

        for warning in self.warnings:
            self.layout.addWidget(warning)
        self.setLayout(self.layout)

class Parameters():
    def __init__(self):
        # QLineEdit
        self.mass = QtWidgets.QLineEdit()
        self.damp_cons = QtWidgets.QLineEdit()
        self.spring_cons = QtWidgets.QLineEdit()
        self.spring_leng = QtWidgets.QLineEdit()
        self.side_cube = QtWidgets.QLineEdit()
        
        self.end_time = QtWidgets.QLineEdit()
        self.fram_res = QtWidgets.QLineEdit()
        self.amp = QtWidgets.QLineEdit()
        self.freq = QtWidgets.QLineEdit()
    
        self.warnings = []

        self._set_default()
        self.update_params()

        self.text = """ vazzio por enquanto """

    def add_widgtes(self, layout):
        layout.addWidget(QtWidgets.QLabel('Final time'))
        layout.addWidget(self.end_time)

        layout.addWidget(QtWidgets.QLabel('Frame rate'))
        layout.addWidget(self.fram_res)

        layout.addWidget(QtWidgets.QLabel('Mass'))
        layout.addWidget(self.mass)

        layout.addWidget(QtWidgets.QLabel('Spring constant'))
        layout.addWidget(self.spring_cons)

        layout.addWidget(QtWidgets.QLabel('Damping constant'))
        layout.addWidget(self.damp_cons)

        layout.addWidget(QtWidgets.QLabel('Spring relaxed length'))
        layout.addWidget(self.spring_leng)

        layout.addWidget(QtWidgets.QLabel('Side cube'))
        layout.addWidget(self.side_cube)

        layout.addWidget(QtWidgets.QLabel('Amplitude'))
        layout.addWidget(self.amp)

        layout.addWidget(QtWidgets.QLabel('Frequency'))
        layout.addWidget(self.freq)
    
    def _set_default(self):       
        self.end_time.setText('30')
        self.fram_res.setText('30')

        self.mass.setText('1000')
        self.spring_cons.setText('48000')
        self.damp_cons.setText('4000')

        self.spring_leng.setText('0.9')
        self.side_cube.setText('1')

        self.amp.setText('0.2')
        self.freq.setText('1')

    def update_params(self): 
        self.end_time_par = ast.literal_eval(self.end_time.text())
        self.fram_res_par = ast.literal_eval(self.fram_res.text())
        
        self.mass_par   = ast.literal_eval(self.mass.text())
        self.spring_cons_par   = ast.literal_eval(self.spring_cons.text())
        self.damp_cons_par = ast.literal_eval(self.damp_cons.text())
        
        self.spring_leng_par = ast.literal_eval(self.spring_leng.text())
        self.side_cube_par   = ast.literal_eval(self.side_cube.text())
        
        self.amp_par = ast.literal_eval(self.amp.text())
        self.freq_par  = ast.literal_eval(self.freq.text())

    def check_params(self):
        pass

class ParametersText():
    def __init__(self):
        self.editor = QtWidgets.QLabel()
        self.editor.setTextInteractionFlags(QtCore.Qt.LinksAccessibleByMouse)
        self.editor.setTextFormat(QtCore.Qt.RichText)
        self.editor.setTextInteractionFlags(QtCore.Qt.TextBrowserInteraction)
        self.editor.setOpenExternalLinks(True)

    def set_text(self, text):
        self.editor.setText(text)

class MainWindow(QtWidgets.QDialog):
    def __init__(self):
        super(MainWindow, self).__init__()
        # Progress bar
        self.pbar = QtWidgets.QProgressBar()
        
        # lateral menu
        left_layout = QtWidgets.QVBoxLayout()

        self.param = Parameters()
        self.param.add_widgtes(left_layout)

        self.button_run = QtWidgets.QPushButton('Run')
        self.button_stop = QtWidgets.QPushButton('Stop')
        left_layout.addWidget(self.button_run)
        self.counter_button_run = 0
        self.button_run.clicked.connect(lambda: self.run())

        left_layout.addStretch(5)
        left_layout.setSpacing(20)
        self.left_widget = QtWidgets.QWidget()
        self.left_widget.setLayout(left_layout)
        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(self.left_widget)
        scroll.setWidgetResizable(True)

        # window 
        self.right_widget = QtWidgets.QTabWidget()
        param_text = ParametersText()
        param_text.set_text(self.param.text)
        self.tab1 = self.ui1(param_text.editor)
        self.right_widget.addTab(self.tab1, '')

        self.tab2 = self.ui2()
        self.right_widget.addTab(self.tab2, '')
        
        self.right_widget.setCurrentIndex(0)
        self.right_widget.setStyleSheet('''QTabBar::tab{width: 0; \
            height: 0; margin: 0; padding: 0; border: none;}''')

        scroll_right = QtWidgets.QScrollArea()
        scroll_right.setWidget(self.right_widget)
        scroll_right.setWidgetResizable(True)

        # main window 
        main_layout = QtWidgets.QHBoxLayout()
        main_layout.addWidget(scroll)
        main_layout.addWidget(scroll_right)
        main_layout.setStretch(0, 56)
        main_layout.setStretch(1, 300)
        self.setLayout(main_layout) 

    # ------- pages -------
    def ui1(self, param_text):
        layout_ui1 = QtWidgets.QVBoxLayout()
        layout_ui1.addWidget(param_text)
        layout_ui1.addStretch(20)
        
        main = QtWidgets.QWidget()
        main.setLayout(layout_ui1)
        return main
    
    def ui2(self):
        self.init_ui2_layout = QtWidgets.QVBoxLayout()
        main = QtWidgets.QWidget()
        main.setLayout(self.init_ui2_layout)
        return main

    def add_sys_ui2(self):
        self.frame = Qt.QFrame()
        new_layout = Qt.QVBoxLayout()
        self.vtkWidget = QVTKRenderWindowInteractor(self.frame)

        # Create renderer and add the vedo objects
        self.vp = vd.Plotter(axes=4, qtWidget=self.vtkWidget)
        self.vp.show(azimuth=90, elevation=80, roll=90, resetcam=True)
        

        # Set-up the rest of the Qt window
        self.pbar = QtWidgets.QProgressBar()
        new_layout.addWidget(self.pbar)
        new_layout.addWidget(self.vtkWidget)
        self.update_ui2_layout(new_layout)

    def update_ui2_layout(self, new_layout):
        if self.tab2.layout() is not None:
            QtWidgets.QWidget().setLayout(self.tab2.layout())
        self.tab2.setLayout(new_layout)   

    # ------- buttons ------- 
    @QtCore.pyqtSlot()
    def run(self):
        self.param.check_params()
        
        if len(self.param.warnings) != 0:
            dlg = PopUp(self.param.warnings)      
            dlg.exec_()
        else:
            self.counter_button_run += 1
            if self.counter_button_run == 1:
                self.right_widget.setCurrentIndex(1)
                
            self.param.update_params()
            self.add_sys_ui2()
            
            self.worker = WorkerThread(parent=self)
            self.worker.start()
            self.worker.complete_worker.connect(lambda: self.worker.quit())
            self.worker.complete_worker.connect(lambda: self.worker.deleteLater())
            self.worker.update_plot.connect(self.evt_update)
            self.worker.update_progess.connect(self.evt_update_progress)
            print('passou um aviao')
            
    def evt_update(self):
        self.vp.interactor.Render()

    def evt_update_progress(self, val):
        self.pbar.setValue(val)

class WorkerThread(QtCore.QThread):
    
    update_plot = QtCore.pyqtSignal(vd.plotter.Plotter)

    update_progess = QtCore.pyqtSignal(int)

    complete_worker = QtCore.pyqtSignal(bool)

    def __init__(self, parent):    
        super().__init__()
        self.parent = parent
            
    @QtCore.pyqtSlot()
    def run(self):
        # Create mass spring damper system
        self.sys1 = MassSpringDamper(self.parent.param.mass_par, self.parent.param.damp_cons_par, self.parent.param.spring_cons_par, \
                                  self.parent.param.spring_leng_par, self.parent.param.side_cube_par, self.parent.param.end_time_par, \
                                  self.parent.param.fram_res_par, self.parent.param.amp_par, self.parent.param.freq_par)

        wall = vd.shapes.Box(pos=(0, 0, -0.8), length=1.5 * self.sys1.side_cube, width=1.5 * self.sys1.side_cube, \
                            height=0.05 * self.sys1.side_cube, c=(128,128,128))  
        
        const = 1#2 if self.parent.param.two_sys_par else 1
        
        height_surf = const*(0.8 + np.max(self.sys1.yc) + self.sys1.side_cube)
        
        surface = vd.shapes.Box(pos=(0, -self.sys1.side_cube, -0.8 + (height_surf/2)), length=2 * self.sys1.side_cube, width=0.02, \
                                height=height_surf, alpha=0)

        self.parent.vp += [wall, surface, self.sys1.block, self.sys1.spring, self.sys1.tube, self.sys1.cylinder, self.sys1.redpt, \
                self.sys1.line_cos, self.sys1.line_x, self.sys1.line_y, self.sys1.txt_U, self.sys1.txt_time, self.sys1.redpt_block]

        self.update_plot.emit(self.parent.vp)

        QtWidgets.QApplication.processEvents()

        for i, yc_i in enumerate(self.sys1.yc):
            self.parent.vp.actors[6].pos(self.sys1.pos_points[i])
            self.parent.vp.actors[2].pos([0, 0, yc_i + self.sys1.block_ratio])  
            self.parent.vp.actors[-1].pos(0, self.sys1.block_ratio, yc_i + self.sys1.block_ratio)
            self.parent.vp.actors[3].stretch([self.sys1.side_cube/4, 0,-0.8],[self.sys1.side_cube/4, 0, yc_i])
            self.parent.vp.actors[4].stretch([-self.sys1.side_cube/4, 0, -0.2], [-self.sys1.side_cube/4, 0, yc_i])
            self.update_plot.emit(self.parent.vp)
            self.update_progess.emit(i)
            QtWidgets.QApplication.processEvents() 
        self.complete_worker.emit(True)


if __name__ == "__main__":

    app = QtWidgets.QApplication(sys.argv)
    widget = QtWidgets.QStackedWidget()
    mainwindow = MainWindow()
    widget.addWidget(mainwindow)
    widget.showMaximized()
    sys.exit(app.exec_())

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

    # sys1 = MassSpringDamper(mass, damp_cons, spring_cons, spring_leng, side_cube, end_time, fram_res, amp, freq)

    # sys2 = sys1.copy()

    # plt = vd.Plotter(interactive=0, axes=4)

