#
# Programmer: Jan Mooiman
# email: jan.mooiman@outlook.com
#
from PyQt6.QtWidgets import *
from PyQt6.QtCore import Qt
from PyQt6 import QtGui
from PyQt6 import QtCore

import reg_1d_scalar

class DFTgui(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.setWindowTitle('Regularization')

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

        nrow += 1
        self.bath_label = QLabel(self)
        self.bath_label.setText('Bathymetry:')
        self.bath_combobox = QComboBox(self)
        self.bath_combobox.addItem("tanh + step")
        self.bath_combobox.addItem("tanh + step (deeper)")
        self.bath_combobox.addItem("Frank Platzek")
        self.bath_combobox.addItem("shoal: -10 [m] to -2.5 [m]")
        self.bath_combobox.addItem("Weir")
        self.bath_combobox.addItem("Step function")
        self.bath_combobox.addItem("Constant")
        self.bath_combobox.addItem("Step at right boundary")
        self.bath_combobox.setToolTip("Several scalar profiles")
        edit_layout.addWidget(self.bath_label, nrow, 0)
        edit_layout.addWidget(self.bath_combobox, nrow, 1)
        self.bath_combobox.setCurrentIndex(5)

        nrow += 1
        self.lx_label = QLabel(self)
        self.lx_label.setText('Length')
        self.lx_edit = QLineEdit(self)
        self.lx_edit.setToolTip("")
        edit_layout.addWidget(self.lx_label, nrow, 0)
        edit_layout.addWidget(self.lx_edit, nrow, 1)
        self.lx_edit.setText("1000.")

        nrow += 1
        self.dx_label = QLabel(self)
        self.dx_label.setText('Dx:')
        self.dx_edit = QLineEdit(self)
        self.dx_edit.setToolTip("Dx > 0")
        edit_layout.addWidget(self.dx_label, nrow, 0)
        edit_layout.addWidget(self.dx_edit, nrow, 1)
        self.dx_edit.setText("20.")

        nrow += 1
        self.cpsi_label = QLabel(self)
        self.cpsi_label.setText('c_psi:')
        self.cpsi_edit = QLineEdit(self)
        self.cpsi_edit.setToolTip("c_psi")
        edit_layout.addWidget(self.cpsi_label, nrow, 0)
        edit_layout.addWidget(self.cpsi_edit, nrow, 1)
        self.cpsi_edit.setText("4.")

        nrow += 1
        self.step_left_label = QLabel(self)
        self.step_left_label.setText('step left:')
        self.step_left_edit = QLineEdit(self)
        self.step_left_edit.setToolTip("step right = 100. * step left: ")
        edit_layout.addWidget(self.step_left_label, nrow, 0)
        edit_layout.addWidget(self.step_left_edit, nrow, 1)
        self.step_left_edit.setText("0.1")

        # layout.addStretch()
        layout.addLayout(edit_layout)
        layout.addWidget(simulate_button)
        layout.addWidget(quit_button)


    def run(self):
        lx = self.lx_edit.text()
        dx = self.dx_edit.text()
        c_psi = self.cpsi_edit.text()
        bath = self.bath_combobox.currentIndex()
        step_left = self.step_left_edit.text()

        reg_1d_scalar.main(bath, lx, dx, c_psi, step_left)


if __name__ == '__main__':
    app = QApplication([])
    app.setStyle('Fusion')

    interface = DFTgui()
    interface.setMinimumWidth(325)
    interface.move(5, 5)
    interface.show()
    app.exec()
