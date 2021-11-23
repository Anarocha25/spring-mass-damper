from PyQt5 import QtWidgets, QtGui
from PyQt5.QtCore import Qt, QObject, QThread, pyqtSignal, pyqtSlot, QRunnable, QThreadPool, QProcess

import sys

import matplotlib.pyplot as plt
import numpy as np
import pyqtgraph as pg

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from class_fem_solution import SolutionFEM2D
from class_plots_fem2d import PlotsFem2d
from class_parameters import ParametersFEM2D, ParametersOpt, ParametersText
from class_beam import Beam
from class_optimization import MmaOpt, GcmmaOpt
from class_plot_opt import PlotsOpt

class PopUp(QtWidgets.QDialog):
    def __init__(self, warnings):
        super().__init__()

        self.warnings = warnings

        self.setWindowTitle("Warning")
        self.layout = QtWidgets.QVBoxLayout()

        for warning in self.warnings:
            self.layout.addWidget(warning)
        self.setLayout(self.layout)

class MainWindow(QtWidgets.QDialog):
    def __init__(self):
        super(MainWindow, self).__init__()
        layout = QtWidgets.QHBoxLayout()

        button_fem2d = QtWidgets.QPushButton('FEM2D', self)
        button_fem2d.clicked.connect(self.gotoFem2d)
        layout.addWidget(button_fem2d)

        button_opt = QtWidgets.QPushButton('Optimization', self)
        button_opt.clicked.connect(self.gotoOpt)
        layout.addWidget(button_opt)

        button_fem3d = QtWidgets.QPushButton('FEM3D', self)
        button_fem3d.clicked.connect(self.gotoFem2d)
        layout.addWidget(button_fem3d)
        self.setLayout(layout)
       
    def gotoFem2d(self):
        param_2d = ParametersFEM2D()
        fem2d = WindowsFem2d(param_2d)
        widget.addWidget(fem2d)
        widget.setCurrentIndex(widget.currentIndex()+1)

    def gotoOpt(self):
        param_opt = ParametersOpt()
        opt = WindowsOptimization(param_opt)
        widget.addWidget(opt)
        widget.setCurrentIndex(widget.currentIndex()+1)

class SecondWindow(QtWidgets.QDialog):
    def __init__(self, parameters):
        super(SecondWindow, self).__init__()

        self.pbar = QtWidgets.QProgressBar()
        
        # lateral menu
        left_layout = QtWidgets.QVBoxLayout()
        self.button_back = QtWidgets.QPushButton('Back')
        self.button_back.clicked.connect(self.gotoScreen1)
        left_layout.addWidget(self.button_back)

        self.parameters = parameters
        self.parameters.add_widgtes(left_layout)

        self.button_run = QtWidgets.QPushButton('Run')
        left_layout.addWidget(self.button_run)
        self._counter_button_run = 0
        #self.button_run.clicked.connect(lambda: self.run())

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
        param_text.set_text(self.parameters.text)
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
        main_layout.setStretch(0, 56) #aqui aumenta a barra lateral
        main_layout.setStretch(1, 300)
        self.setLayout(main_layout) 

    # ------- pages -------
    def gotoScreen1(self):
        mainwindow = MainWindow()
        widget.addWidget(mainwindow)
        widget.setCurrentIndex(widget.currentIndex()+1)

    def ui1(self, parameters_text):
        layout_ui1 = QtWidgets.QVBoxLayout()
        layout_ui1.addWidget(parameters_text)
        layout_ui1.addStretch(20)
        
        main = QtWidgets.QWidget()
        main.setLayout(layout_ui1)
        return main
    
    def ui2(self):
        self.init_ui2_lay = QtWidgets.QVBoxLayout()
        main = QtWidgets.QWidget()
        main.setLayout(self.init_ui2_lay)
        return main

    # ------- buttons ------- 
    def run(self):
        self.parameters.check_params()
        if len(self.parameters.warnings) != 0:
            dlg = PopUp(self.parameters.warnings)      
            dlg.exec_()
        else:
            self._counter_button_run += 1
            if self._counter_button_run == 1:
                self.right_widget.setCurrentIndex(1)
                print("entrou")
            self.parameters.update_params()

    # ------- functions ------- 
    def update_pbar(self, current_idx):
        self.pbar.setValue(current_idx)

class WindowsFem2d(SecondWindow):
    def __init__(self, parameters):
        super().__init__(parameters)
        self.add_plot()
        self.button_run.clicked.connect(lambda: self.run())
        self.parameters.freqrsp_check.toggled.connect(self.parameters.freq_range_spin.setEnabled)
        self.parameters.freqrsp_check.toggled.connect(self.parameters.node_plot_spin.setEnabled)

    # ------- pages -------
    def update_ui2(self, parameters, canvas): 
        self.update_pbar(2)
        fem_2d = SolutionFEM2D(parameters.nelx_par, parameters.nely_par, parameters.lx_par, parameters.ly_par, parameters.E_par, parameters.v_par, parameters.rho_par, parameters.alpha_par, parameters.beta_par, parameters.eta_par, parameters.freq_par, parameters.freq_range_par, parameters.matrix_F_par, parameters.constrain_par, parameters.freqrsp, parameters.node_plot_par)
        self.update_pbar(23)
        
        plots_2d = PlotsFem2d(fem_2d.coord, fem_2d.connect, fem_2d.lx, fem_2d.ly, fem_2d.load_matrix, fem_2d.restri_matrix, fem_2d.freq_rsp, fem_2d.delta, parameters.factor_par, fem_2d.disp_vector)
        deform_mesh = plots_2d.build_collection()
        self.update_pbar(38)
        canvas.figure.clear()

        if parameters.freqrsp:
            ax = canvas.figure.subplots(2, 1)
            plots_2d.plot_collection(ax[0], deform_mesh)
            self.update_pbar(50)
            vector_U = fem_2d.freqresponse(self.update_pbar)
            plots_2d.plot_freqresponse(ax[1], vector_U)
            ax[1].set_xlim(fem_2d.freq_rsp[0] + fem_2d.delta, fem_2d.freq_rsp[1])
        else:
            ax = canvas.figure.subplots()
            plots_2d.plot_collection(ax, deform_mesh)
            self.update_pbar(100)
        canvas.draw()

        if parameters.save:
            self.canvas.figure.savefig('plots.png')
            # TODO:https://stackoverflow.com/questions/4325733/save-a-subplot-in-matplotlib salvar separado

    # ------- buttons -------
    def run(self):
        super().run()
        if len(self.parameters.warnings) == 0:
            self.update_ui2(self.parameters, self.canvas)
        
    # ------- functions -------
    def add_plot(self):
        ui2_layout = QtWidgets.QVBoxLayout()
        ui2_layout.addWidget(self.pbar)

        self.canvas = FigureCanvas(plt.Figure())
        self.canvas.figure.tight_layout()
        ui2_layout.addWidget(self.canvas)

        QtWidgets.QWidget().setLayout(self.init_ui2_lay)
        self.tab2.setLayout(ui2_layout)