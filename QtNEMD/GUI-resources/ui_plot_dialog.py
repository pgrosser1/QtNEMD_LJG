# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'GUI-resources/plot_dialog.ui'
#
# Created by: PyQt5 UI code generator 5.12
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_PlotDialog(object):
    def setupUi(self, PlotDialog):
        PlotDialog.setObjectName("PlotDialog")
        PlotDialog.resize(400, 300)
        self.verticalLayoutWidget = QtWidgets.QWidget(PlotDialog)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(-1, 9, 401, 291))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.ttcf_checkbox = QtWidgets.QCheckBox(self.verticalLayoutWidget)
        self.ttcf_checkbox.setObjectName("ttcf_checkbox")
        self.verticalLayout.addWidget(self.ttcf_checkbox)
        self.msd_checkbox = QtWidgets.QCheckBox(self.verticalLayoutWidget)
        self.msd_checkbox.setObjectName("msd_checkbox")
        self.verticalLayout.addWidget(self.msd_checkbox)
        self.buttonBox = QtWidgets.QDialogButtonBox(self.verticalLayoutWidget)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setCenterButtons(True)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout.addWidget(self.buttonBox)

        self.retranslateUi(PlotDialog)
        QtCore.QMetaObject.connectSlotsByName(PlotDialog)

    def retranslateUi(self, PlotDialog):
        _translate = QtCore.QCoreApplication.translate
        PlotDialog.setWindowTitle(_translate("PlotDialog", "Form"))
        self.ttcf_checkbox.setText(_translate("PlotDialog", "TTCF"))
        self.msd_checkbox.setText(_translate("PlotDialog", "Mean-square displacement"))


