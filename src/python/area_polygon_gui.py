#
# Programmer: Jan Mooiman
# email: jan.mooiman@outlook.com
#
import numpy as np
from PyQt6.QtWidgets import *
from PyQt6.QtCore import Qt
from PyQt6 import QtGui
from PyQt6 import QtCore

class DFTgui(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.setWindowTitle('Area polygon')

        # Main Widget
        widget = QWidget()
        self.setCentralWidget(widget)

        # Layout
        layout = QVBoxLayout()
        widget.setLayout(layout)

        # Horizontal layouts
        edit_layout = QGridLayout()

        # Quit button
        quit_button = QPushButton('Quit')
        quit_button.clicked.connect(app.exit)

        # Simulate System Button
        simulate_button = QPushButton('Run')
        simulate_button.clicked.connect(self.run)

        # Labels
        self.groupbox_1 = QGroupBox('Regularization parameters')

        # Line edit fields
        nrow = -1


        self.lx_labelx1 = QLabel('x1')
        self.lx_labely1 = QLabel('y1')
        self.lx_labelx2 = QLabel('x2')
        self.lx_labely2 = QLabel('y2')
        self.lx_labelx3 = QLabel('x3')
        self.lx_labely3 = QLabel('y3')
        self.lx_labelx4 = QLabel('x4')
        self.lx_labely4 = QLabel('y4')
        #self.lx_labely.setText('Y')
        self.lx_editx1 = QLineEdit(self)
        self.lx_editx2 = QLineEdit(self)
        self.lx_editx3 = QLineEdit(self)
        self.lx_editx4 = QLineEdit(self)
        self.lx_edity1 = QLineEdit(self)
        self.lx_edity2 = QLineEdit(self)
        self.lx_edity3 = QLineEdit(self)
        self.lx_edity4 = QLineEdit(self)

        nrow += 1
        edit_layout.addWidget(self.lx_labelx1, nrow, 0)
        edit_layout.addWidget(self.lx_editx1, nrow, 1)
        edit_layout.addWidget(self.lx_labely1, nrow, 2)
        edit_layout.addWidget(self.lx_edity1, nrow, 3)

        nrow += 1
        edit_layout.addWidget(self.lx_labelx2, nrow, 0)
        edit_layout.addWidget(self.lx_editx2, nrow, 1)
        edit_layout.addWidget(self.lx_labely2, nrow, 2)
        edit_layout.addWidget(self.lx_edity2, nrow, 3)

        nrow += 1
        edit_layout.addWidget(self.lx_labelx3, nrow, 0)
        edit_layout.addWidget(self.lx_editx3, nrow, 1)
        edit_layout.addWidget(self.lx_labely3, nrow, 2)
        edit_layout.addWidget(self.lx_edity3, nrow, 3)

        nrow += 1
        edit_layout.addWidget(self.lx_labelx4, nrow, 0)
        edit_layout.addWidget(self.lx_editx4, nrow, 1)
        edit_layout.addWidget(self.lx_labely4, nrow, 2)
        edit_layout.addWidget(self.lx_edity4, nrow, 3)

        self.lx_editx1.setText("50.")
        self.lx_editx2.setText("100.")
        self.lx_editx3.setText("114.")
        self.lx_editx4.setText("60.")
        self.lx_edity1.setText("10.")
        self.lx_edity2.setText("5.")
        self.lx_edity3.setText("50.")
        self.lx_edity4.setText("42.5")

        # layout.addStretch()
        layout.addLayout(edit_layout)
        layout.addWidget(simulate_button)
        layout.addWidget(quit_button)

    def run(self):
        nx = 4
        x = np.zeros(nx, dtype=np.float64)
        y = np.zeros(nx, dtype=np.float64)
        edge_x = np.zeros(nx, dtype=np.float64)
        edge_y = np.zeros(nx, dtype=np.float64)
        scv_x = np.zeros(nx, dtype=np.float64)
        scv_y = np.zeros(nx, dtype=np.float64)
        scv = np.zeros(nx, dtype=np.float64)
        x[0] = self.lx_editx1.text()
        x[1] = self.lx_editx2.text()
        x[2] = self.lx_editx3.text()
        x[3] = self.lx_editx4.text()
        y[0] = self.lx_edity1.text()
        y[1] = self.lx_edity2.text()
        y[2] = self.lx_edity3.text()
        y[3] = self.lx_edity4.text()

        cv = self.compute_area(x, y)

        mc_x = 0.0
        mc_y = 0.0
        for i in range(0, len(x)):
            mc_x += x[i]
        for i in range(0, len(x)):
            mc_y += y[i]
        mc_x = mc_x/len(x)
        mc_y = mc_y/len(x)

        for i in range(0, len(x)-1):
            edge_x[i] = 0.5 *(x[i] + x[i+1])
        i = len(x)-1
        edge_x[i] = 0.5 * (x[i] + x[0])
        for i in range(0, len(x)-1):
            edge_y[i] = 0.5 *(y[i] + y[i+1])
        i = len(x)-1
        edge_y[i] = 0.5 * (y[i] + y[0])

        scv_x[0] = x[0]
        scv_x[1] = edge_x[0]
        scv_x[2] = mc_x
        scv_x[3] = edge_x[3]
        scv_y[0] = y[0]
        scv_y[1] = edge_y[0]
        scv_y[2] = mc_y
        scv_y[3] = edge_y[3]
        scv[0] = self.compute_area(scv_x, scv_y)

        scv_x[0] = x[1]
        scv_x[1] = edge_x[1]
        scv_x[2] = mc_x
        scv_x[3] = edge_x[0]
        scv_y[0] = y[1]
        scv_y[1] = edge_y[1]
        scv_y[2] = mc_y
        scv_y[3] = edge_y[0]
        scv[1] = self.compute_area(scv_x, scv_y)

        scv_x[0] = x[2]
        scv_x[1] = edge_x[2]
        scv_x[2] = mc_x
        scv_x[3] = edge_x[1]
        scv_y[0] = y[2]
        scv_y[1] = edge_y[2]
        scv_y[2] = mc_y
        scv_y[3] = edge_y[1]
        scv[2] = self.compute_area(scv_x, scv_y)

        scv_x[0] = x[3]
        scv_x[1] = edge_x[3]
        scv_x[2] = mc_x
        scv_x[3] = edge_x[2]
        scv_y[0] = y[3]
        scv_y[1] = edge_y[3]
        scv_y[2] = mc_y
        scv_y[3] = edge_y[2]
        scv[3] = self.compute_area(scv_x, scv_y)

        print(' Control volume: %f' % cv)
        print(' Sub-CV        : %f, %f, %f, %f' % (scv[0], scv[1], scv[2], scv[3]))
        print(' Difference    : %f' % (cv - scv[0] - scv[1] - scv[2] - scv[3]))

        return

    def compute_area(self, x, y):
        area = 0.0
        for i in range(0, len(x)-1):
            area += x[i] * y[i+1]
        area += x[len(x)-1] * y[0]
        for i in range(0, len(x)-1):
            area -= x[i+1] * y[i]
        area -= x[0] * y[len(x)-1]
        area *= 0.5
        return area


if __name__ == '__main__':
    app = QApplication([])
    app.setStyle('Fusion')

    interface = DFTgui()
    interface.setMinimumWidth(325)
    interface.move(5, 5)
    interface.show()
    app.exec()
