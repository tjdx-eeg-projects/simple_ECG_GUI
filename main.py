# -*- coding: utf-8 -*-
import sys
from PyQt4.QtGui import QApplication, QDialog
import ui_hrv

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = QDialog()
    ui = ui_hrv.Ui_Dialog()
    ui.setupUi(window)
    window.show()
    sys.exit(app.exec_())