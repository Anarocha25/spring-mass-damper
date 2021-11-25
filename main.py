import sys
from PyQt5 import QtWidgets
from z_window import MainWindow
if __name__ == "__main__":

    app = QtWidgets.QApplication(sys.argv)
    widget = QtWidgets.QStackedWidget()
    mainwindow = MainWindow()
    widget.addWidget(mainwindow)
    widget.showMaximized()
    sys.exit(app.exec_())