import numpy as np
import vedo as vd

import sys
import os
import ast
import cv2
from PyQt5 import Qt, QtWidgets, QtCore
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

from z_system import MassSpringDamperSys1, MassSpringDamperSys2

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

        self.two_sys_check = QtWidgets.QCheckBox("Two systems")

        self.mass2 = QtWidgets.QLineEdit()
        self.mass2.setDisabled(True)
        self.damp_cons2 = QtWidgets.QLineEdit()
        self.damp_cons2.setDisabled(True)
        self.spring_cons2 = QtWidgets.QLineEdit()
        self.spring_cons2.setDisabled(True)

        self.spring_leng2 = QtWidgets.QLineEdit()
        self.spring_leng2.setDisabled(True)
        self.side_cube2 = QtWidgets.QLineEdit()
        self.side_cube2.setDisabled(True)
        self.amp2 = QtWidgets.QLineEdit()
        self.amp2.setDisabled(True)
        self.freq2 = QtWidgets.QLineEdit()
        self.freq2.setDisabled(True)

        self.save_check = QtWidgets.QCheckBox("Save")
        self.name_save_line = QtWidgets.QLineEdit()
        self.name_save_line.setDisabled(True)

        self.warnings = []

        self._set_default()
        self.update_params()

        self.text =  """ 
                    <p style="text-align:center;"> <b> <font size="+8"> Spring mass damper system </font> </b>
                   <hr>
                    
                    <p> <b> <font size="+5"> Parameters </font> </b>
                    
                    <p> <b> <font size="+3"> Video </font> </b>
                    
                    <p style="margin-left:2em"> <font size="+3"> end_time (<i>int</i>):  Final time [s]  </font> </p>
                    <p style="margin-left:2em"> <font size="+3"> fram_res (<i>int</i>): Frame rate [fps] </font> </p>
                    
                    <p> <b> <font size="+3"> System </font> </b>
                    
                    <p style="margin-left:2em"> <font size="+3"> mass (<i>float</i>): Mass [kg] </font> </p>
                    <p style="margin-left:2em"> <font size="+3"> spring_cons (<i>float</i>): Spring constant [N/m] </font> </p>
                    <p style="margin-left:2em"> <font size="+3"> damp_cons (<i>float</i>): Damping constant </font> </p>
                    
                    <p> <b> <font size="+3"> Animation Model </font> </b>

                    <p style="margin-left:2em"> <font size="+3"> spring_leng  (<i>float</i>): Spring relaxed length [m] </font> </p>
                    <p style="margin-left:2em"> <font size="+3"> side_cube  (<i>int</i>): Size cube [m] </font> </p>
                    
                    <p> <b> <font size="+3"> Base input </font> </b>

                    <p style="margin-left:2em"> <font size="+3">  amp (<i>float</i>): Amplitude [m] </font> </p>
                    <p style="margin-left:2em"> <font size="+3">  freq (<i>int</i>): Frequency [Hz] </font> </p>
                    
                    <p> <b> <font size="+3"> Two systems </font> </b>

                    <p style="margin-left:2em"> <font size="+3"> two_sys (<i>bool</i>): If True is considered two mass spring damper systems. </font> </p>
                    <p style="margin-left:14em"> <font size="+3"> The parameters are the same from the first system. </p>
                    
                    """

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

        layout.addWidget(self.two_sys_check)

        layout.addWidget(QtWidgets.QLabel('Second system parameters:'))

        layout.addWidget(QtWidgets.QLabel('Mass'))
        layout.addWidget(self.mass2)

        layout.addWidget(QtWidgets.QLabel('Spring constant'))
        layout.addWidget(self.spring_cons2)

        layout.addWidget(QtWidgets.QLabel('Damping constant'))
        layout.addWidget(self.damp_cons2)

        layout.addWidget(QtWidgets.QLabel('Spring relaxed length'))
        layout.addWidget(self.spring_leng2)

        layout.addWidget(QtWidgets.QLabel('Side cube'))
        layout.addWidget(self.side_cube)

        layout.addWidget(QtWidgets.QLabel('Amplitude'))
        layout.addWidget(self.amp2)

        layout.addWidget(QtWidgets.QLabel('Frequency'))
        layout.addWidget(self.freq2)

        layout.addWidget(self.save_check)

        layout.addWidget(QtWidgets.QLabel('Name to save'))
        layout.addWidget(self.name_save_line)
    
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

        self.mass2.setText('1000')
        self.spring_cons2.setText('48000')
        self.damp_cons2.setText('4000')

        self.spring_leng2.setText('0.9')
        self.side_cube.setText('1')

        self.amp2.setText('0.2')
        self.freq2.setText('1')

        self.name_save_line.setText('system')

    def update_params(self): 
        self.end_time_par = ast.literal_eval(self.end_time.text())
        self.fram_res_par = ast.literal_eval(self.fram_res.text())
        
        self.mass_par = ast.literal_eval(self.mass.text())
        self.spring_cons_par = ast.literal_eval(self.spring_cons.text())
        self.damp_cons_par = ast.literal_eval(self.damp_cons.text())
        
        self.spring_leng_par = ast.literal_eval(self.spring_leng.text())
        self.side_cube_par   = ast.literal_eval(self.side_cube.text())
        
        self.amp_par = ast.literal_eval(self.amp.text())
        self.freq_par = ast.literal_eval(self.freq.text())

        # Second system
        self.two_sys = self.two_sys_check.checkState()
        if self.two_sys:
            self.mass_par2 = ast.literal_eval(self.mass2.text())
            self.spring_cons_par2 = ast.literal_eval(self.spring_cons2.text())
            self.damp_cons_par2 = ast.literal_eval(self.damp_cons2.text())
            
            self.spring_leng_par2 = ast.literal_eval(self.spring_leng2.text())
            self.side_cube_par2   = ast.literal_eval(self.side_cube.text())
            
            self.amp_par2 = ast.literal_eval(self.amp2.text())
            self.freq_par2 = ast.literal_eval(self.freq2.text())
        
        self.save = self.save_check.checkState()
        self.name_save = self.name_save_line.text()

    def export_param(self):
        param1 = {"fram_res_par":self.fram_res_par, "mass_par":self.mass_par, "damp_cons_par":self.damp_cons_par, "spring_cons_par":self.spring_cons_par, \
        "spring_leng_par":self.spring_leng_par, "side_cube_par":self.side_cube_par, "end_time_par":self.end_time_par, "fram_res_par":self.fram_res_par, \
        "amp_par":self.amp_par, "freq_par":self.freq_par, "two_sys":self.two_sys, "name_save":self.name_save}

        if self.two_sys:
            param2 = {"mass_par2":self.mass_par2, "damp_cons_par2":self.damp_cons_par2, \
            "spring_cons_par2":self.spring_cons_par2, "spring_leng_par2":self.spring_leng_par2, "side_cube_par2":self.side_cube_par2, \
            "amp_par2":self.amp_par2, "freq_par2":self.freq_par2}
        
        folder_name = 'temp'
        directory = os.path.join(os.path.dirname(__file__), folder_name)
        os.makedirs(directory, exist_ok=True)

        try: 
            # WRITE
            dir_file = os.path.join(directory, 'param1_file.txt')
            # open file for writing
            f = open(dir_file,"w")
            # write file
            f.write(str(param1))
            # close file
            f.close()

            if self.two_sys:
                # WRITE
                dir_file = os.path.join(directory, 'param2_file.txt')
                # open file for writing
                f = open(dir_file,"w")
                # write file
                f.write(str(param2))
                # close file
                f.close()

        except: 
            print("Something went wrong")

    def check_param(self, input_val, types, war):
        try:
            for type in types:
                isinstance(ast.literal_eval(input_val), type)
        except ValueError:
            warning = QtWidgets.QLabel(war)
            self.warnings.append(warning)
        
    def check_params(self):
        self.warnings = []

        self.check_param(self.mass.text(), [int, float], 'Mass must be an int or float')
        self.check_param(self.damp_cons.text(), [int, float], 'Spring constant must be an int or a float')
        self.check_param(self.spring_cons.text(), [int, float], 'Mass must be an int or a float')

        self.check_param(self.end_time.text(), [int], 'Final time must be an int')
        self.check_param(self.fram_res.text(), [int], 'Frame rate must be an int')

        self.check_param(self.spring_leng.text(), [int, float], 'Spring relaxed length must be an int or float')
        self.check_param(self.side_cube.text(), [int, float], 'Size cube must be an int or float')
        
        self.check_param(self.amp.text(), [int, float], 'Amplitude must be an int or float')
        self.check_param(self.freq.text(), [int], 'Frequency must be an int')

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
        # lateral menu
        left_layout = QtWidgets.QVBoxLayout()

        self.param = Parameters()
        self.param.add_widgtes(left_layout)

        self.param.two_sys_check.toggled.connect(self.param.mass2.setEnabled)
        self.param.two_sys_check.toggled.connect(self.param.spring_cons2.setEnabled)
        self.param.two_sys_check.toggled.connect(self.param.damp_cons2.setEnabled)
        self.param.two_sys_check.toggled.connect(self.param.spring_leng2.setEnabled)
        self.param.two_sys_check.toggled.connect(self.param.side_cube2.setEnabled)
        self.param.two_sys_check.toggled.connect(self.param.amp2.setEnabled)
        self.param.two_sys_check.toggled.connect(self.param.freq2.setEnabled)

        self.param.save_check.toggled.connect(self.param.name_save_line.setEnabled)

        self.button_run = QtWidgets.QPushButton('Run')
        left_layout.addWidget(self.button_run)
        self.counter_button_run = 0
        self.button_run.clicked.connect(lambda: self.run())
    
        self.button_stop = QtWidgets.QPushButton('Stop')
        self.button_stop.setEnabled(False)
        left_layout.addWidget(self.button_stop)

        self.stop_thread = False
        self.process_save = None
        self.worker = None

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
            
        new_layout.addWidget(self.vtkWidget)

        self.text = QtWidgets.QPlainTextEdit()
        self.text.setFixedSize(150, 50)
        self.text.setReadOnly(True)
        new_layout.addWidget(self.text)
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
            self.button_stop.setEnabled(True)
            self.button_run.setEnabled(False)

            self.counter_button_run += 1
            if self.counter_button_run == 1:
                self.right_widget.setCurrentIndex(1)
                
            self.param.update_params()
            self.add_sys_ui2()

            self.message("Running")

            if self.param.save:
                self.param.export_param() 

                self.process_save = QtCore.QProcess()
                #self.process_save.readyReadStandardOutput.connect(self.handle_stdout)
                #self.process_save.readyReadStandardError.connect(self.handle_stderr)
                #self.process_save.stateChanged.connect(self.handle_state)
                self.process_save.start("python", ['z_process.py'])
                self.process_save.finished.connect(lambda: self.finish_process())
                

            self.worker = WorkerThread(parent=self)
            self.worker.start()
            self.worker.update_plot.connect(self.evt_update)
            self.worker.complete_worker.connect(lambda: self.worker.terminate())
            self.worker.complete_worker.connect(lambda: self.worker.deleteLater())
            self.worker.complete_worker.connect(lambda: self.finish_worker())        
            
            self.button_stop.clicked.connect(lambda: self.stop_execution())
            #self.button_stop.clicked.connect(lambda: self.button_run.setEnabled(True))
            #self.button_stop.clicked.connect(lambda: self.button_stop.setEnabled(False))
          

    def message(self, s):
        self.text.setPlainText(s)

    # def handle_stderr(self):
    #     data = self.process_save.readAllStandardError()
    #     stderr = bytes(data).decode("utf8")
    #     # Extract progress if it is in the data.
    #     self.message(stderr)

    # def handle_stdout(self):
    #     data = self.process_save.readAllStandardOutput()
    #     stdout = bytes(data).decode("utf8")
    #     self.message(stdout)

    # def handle_state(self, state):
    #     states = {
    #         QtCore.QProcess.NotRunning: 'Not running',
    #         QtCore.QProcess.Starting: 'Starting',
    #         QtCore.QProcess.Running: 'Running',
    #     }
    #     state_name = states[state]
    #     print(state_name)

    def stop_execution(self):
        self.message("Stopping")
        if self.worker is not None:
            self.stop_thread = True
        
        if self.process_save is not None:
            self.process_save.kill()
            self.del_process()

        if (self.worker is None) and (self.process_save is None):
            self.message("Finished")

    def del_process(self):
        self.process_save = None

    def finish_worker(self):
        self.worker = None
        if self.process_save is None:
            self.message("Finished")
            self.button_stop.setEnabled(False)
            self.button_run.setEnabled(True)
        else:
            self.message("Saving video")

    def finish_process(self):
        if self.worker is None:
            self.message("Finished")
            self.button_stop.setEnabled(False)
            self.button_run.setEnabled(True)
        self.del_process()

    def evt_update(self):
        self.vp.interactor.Render()            

class WorkerThread(QtCore.QThread):
    
    update_plot = QtCore.pyqtSignal(vd.plotter.Plotter)

    complete_worker = QtCore.pyqtSignal(bool)

    def __init__(self, parent):    
        super().__init__()
        self.parent = parent
            
    @QtCore.pyqtSlot()
    def run(self):
        # Create mass spring damper system
        self.sys1 = MassSpringDamperSys1(self.parent.param.mass_par, self.parent.param.damp_cons_par, self.parent.param.spring_cons_par, \
                                  self.parent.param.spring_leng_par, self.parent.param.side_cube_par, self.parent.param.end_time_par, \
                                  self.parent.param.fram_res_par, self.parent.param.amp_par, self.parent.param.freq_par)
        
        yc_max = np.max(self.sys1.yc)
        side_cube_max = self.sys1.side_cube
  
        if self.parent.param.two_sys:
            self.sys2 = MassSpringDamperSys2(self.sys1.yc, self.parent.param.mass_par2, self.parent.param.damp_cons_par2, self.parent.param.spring_cons_par2, \
                                    self.parent.param.spring_leng_par2, self.parent.param.side_cube_par2, self.parent.param.end_time_par, \
                                    self.parent.param.fram_res_par, self.parent.param.amp_par2, self.parent.param.freq_par2)
            
            base_cyli2 = np.array([-self.sys1.side_cube/4, 0, self.sys2.init_point])
            top_cyli2 = np.array([-self.sys1.side_cube/4, 0, self.sys2.init_point + 0.6])
            
            yc_max = np.maximum(np.max(self.sys2.yc), np.max(self.sys1.yc))
            side_cube_max = np.maximum(self.sys2.side_cube, self.sys1.side_cube)
            
        wall = vd.shapes.Box(pos=(0, 0, -0.8), length=1.5 * self.sys1.side_cube, width=1.5 * self.sys1.side_cube, \
                            height=0.05 * self.sys1.side_cube, c=(128,128,128))  
        
        const = 2 if self.parent.param.two_sys else 1
        
        height_surf = const*(0.8 + yc_max + side_cube_max)
        
        surface = vd.shapes.Box(pos=(0, -side_cube_max, -0.8 + (height_surf/2)), length=2 * side_cube_max, width=0.02, \
                                height=height_surf, alpha=0)

        self.parent.vp += [wall, surface, self.sys1.block, self.sys1.spring, self.sys1.tube, self.sys1.cylinder, self.sys1.redpt, \
                self.sys1.line_cos, self.sys1.line_x, self.sys1.line_y, self.sys1.txt_U, self.sys1.txt_time, self.sys1.redpt_block]
        
        if self.parent.param.two_sys:
            self.parent.vp += [self.sys2.block, self.sys2.spring, self.sys2.tube, self.sys2.cylinder, self.sys2.redpt, \
                    self.sys2.line_cos, self.sys2.line_x, self.sys2.line_y, self.sys2.txt_U, self.sys2.txt_time, self.sys2.redpt_block]
        
        self.parent.vp.resetCamera()

        self.update_plot.emit(self.parent.vp)

        QtWidgets.QApplication.processEvents()

        for i, yc_i in enumerate(self.sys1.yc):
            self.parent.vp.actors[6].pos(self.sys1.pos_points[i]) #redpoint
            self.parent.vp.actors[2].pos([0, 0, yc_i + self.sys1.block_ratio])  #block
            self.parent.vp.actors[12].pos(0, self.sys1.block_ratio, yc_i + self.sys1.block_ratio) #redpt_block
            self.parent.vp.actors[3].stretch([self.sys1.side_cube/4, 0,-0.8],[self.sys1.side_cube/4, 0, yc_i]) #spring
            self.parent.vp.actors[4].stretch([-self.sys1.side_cube/4, 0, -0.2], [-self.sys1.side_cube/4, 0, yc_i]) #tube
            
            if self.parent.param.two_sys:
                init_point = yc_i + self.sys1.side_cube
                yc2_i = self.sys2.yc[i]
                
                self.parent.vp.actors[17].pos(self.sys2.pos_points[i]) #redpt2
                self.parent.vp.actors[13].pos([0, 0, yc2_i + self.sys2.block_ratio + 0.8 + init_point]) #block
                self.parent.vp.actors[23].pos(0, self.sys2.block_ratio, yc2_i + self.sys2.block_ratio + 0.8 + init_point) #redpt_block2
                self.parent.vp.actors[14].stretch([self.sys1.side_cube/4, 0,init_point], [self.sys1.side_cube/4, 0, yc2_i + init_point + 0.8]) #spring
                self.parent.vp.actors[15].stretch([-self.sys1.side_cube/4, 0,init_point + 0.6], [-self.sys1.side_cube/4, 0, yc2_i + init_point + 0.8]) #tube

                base_cyli2[2] = init_point
                top_cyli2[2]  = init_point + 0.6
                pos = (base_cyli2 + top_cyli2)/2
                self.parent.vp.actors[16].pos(pos) #cylinder
            
            self.update_plot.emit(self.parent.vp)
            QtWidgets.QApplication.processEvents()
            if self.parent.stop_thread:
                self.parent.stop_thread = False
                break
        self.complete_worker.emit(True)

if __name__ == "__main__":

    app = QtWidgets.QApplication(sys.argv)
    widget = QtWidgets.QStackedWidget()
    mainwindow = MainWindow()
    widget.addWidget(mainwindow)
    widget.showMaximized()
    sys.exit(app.exec_())