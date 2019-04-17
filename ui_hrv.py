# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\Zhang\HRV.ui'
#
# Created by: PyQt4 UI code generator 4.12.1
#
# WARNING! All changes made in this file will be lost!

import sys
from PyQt4 import QtCore, QtGui
from PyQt4.QtGui import QApplication, QDialog
#from matplotlib.backends.qt_compat import QtCore, QtGui
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
import do_wfdb
import matplotlib.pyplot as plt

import numpy as np
import wfdb
from wfdb.io._signal import downround, upround

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_Dialog(QtGui.QDialog):
    def setupUi(self, Dialog):
        Dialog.setObjectName(_fromUtf8("Dialog"))
        Dialog.resize(1071, 600)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Dialog.sizePolicy().hasHeightForWidth())
        Dialog.setSizePolicy(sizePolicy)
        self.horizontalLayoutWidget_3 = QtGui.QWidget(Dialog)
        self.horizontalLayoutWidget_3.setGeometry(QtCore.QRect(0, 0, 1071, 600))
        self.horizontalLayoutWidget_3.setObjectName(_fromUtf8("horizontalLayoutWidget_3"))
        self.horizontalLayout_3 = QtGui.QHBoxLayout(self.horizontalLayoutWidget_3)
        self.horizontalLayout_3.setMargin(0)
        self.horizontalLayout_3.setObjectName(_fromUtf8("horizontalLayout_3"))
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.pushButton = QtGui.QPushButton(self.horizontalLayoutWidget_3)
        self.pushButton.setObjectName(_fromUtf8("pushButton"))
        self.pushButton.clicked.connect(self.read_file)
        self.verticalLayout.addWidget(self.pushButton)
        '''
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.label_3 = QtGui.QLabel(self.horizontalLayoutWidget_3)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.horizontalLayout_2.addWidget(self.label_3)
        self.comboBox = QtGui.QComboBox(self.horizontalLayoutWidget_3)
        self.comboBox.setObjectName(_fromUtf8("comboBox"))
        self.horizontalLayout_2.addWidget(self.comboBox)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        '''
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.label_2 = QtGui.QLabel(self.horizontalLayoutWidget_3)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.horizontalLayout.addWidget(self.label_2)
        self.spinBox = QtGui.QSpinBox(self.horizontalLayoutWidget_3)
        self.spinBox.setMaximum(1000000)
        self.spinBox.setSingleStep(1000)
        self.spinBox.setObjectName(_fromUtf8("spinBox"))
        self.horizontalLayout.addWidget(self.spinBox)
        self.label = QtGui.QLabel(self.horizontalLayoutWidget_3)
        self.label.setObjectName(_fromUtf8("label"))
        self.horizontalLayout.addWidget(self.label)
        self.spinBox_2 = QtGui.QSpinBox(self.horizontalLayoutWidget_3)
        self.spinBox_2.setMaximum(1000000)
        self.spinBox_2.setSingleStep(1000)
        self.spinBox_2.setProperty("value", 10000)
        self.spinBox_2.setObjectName(_fromUtf8("spinBox_2"))
        self.horizontalLayout.addWidget(self.spinBox_2)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.pushButton_2 = QtGui.QPushButton(self.horizontalLayoutWidget_3)
        self.pushButton_2.setObjectName(_fromUtf8("pushButton_2"))
        self.verticalLayout.addWidget(self.pushButton_2)
        self.label_4 = QtGui.QLabel(self.horizontalLayoutWidget_3)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.verticalLayout.addWidget(self.label_4)
        self.fig = plt.figure()
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.verticalLayout.addWidget(self.toolbar)
        self.verticalLayout.addWidget(self.canvas)
        '''
        self.graphicsView = QtGui.QGraphicsView(self.horizontalLayoutWidget_3)
        self.graphicsView.setObjectName(_fromUtf8("graphicsView"))
        self.verticalLayout.addWidget(self.graphicsView)
        '''
        self.horizontalLayout_3.addLayout(self.verticalLayout)
        self.verticalLayout_2 = QtGui.QVBoxLayout()
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.pushButton_3 = QtGui.QPushButton(self.horizontalLayoutWidget_3)
        self.pushButton_3.setObjectName(_fromUtf8("pushButton_3"))
        self.verticalLayout_2.addWidget(self.pushButton_3)
        self.textBrowser = QtGui.QTextBrowser(self.horizontalLayoutWidget_3)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(12)
        self.textBrowser.setFont(font)
        self.textBrowser.setObjectName(_fromUtf8("textBrowser"))
        self.verticalLayout_2.addWidget(self.textBrowser)
        self.horizontalLayout_3.addLayout(self.verticalLayout_2)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(_translate("Dialog", "HRV", None))
        #Dialog.setToolTip(_translate("Dialog", "<html><head/><body><p><br/></p></body></html>", None))
        #Dialog.setWhatsThis(_translate("Dialog", "<html><head/><body><p><br/></p></body></html>", None))
        self.pushButton.setText(_translate("Dialog", "读取文件", None))
        #self.label_3.setText(_translate("Dialog", "频道", None))
        self.label_2.setText(_translate("Dialog", "样本数", None))
        self.label.setText(_translate("Dialog", "到", None))
        self.pushButton_2.setText(_translate("Dialog", "绘图", None))
        self.label_4.setText(_translate("Dialog", "注释：", None))
        self.pushButton_3.setText(_translate("Dialog", "计算特征", None))
        self.pushButton_2.clicked.connect(self.ui_plot)
        self.pushButton_3.clicked.connect(self.calculate)
    
    def read_file(self):
        file_name = QtGui.QFileDialog.getOpenFileName(caption = "open file dialog", filter = "Header files(*.hea)")
        self.file_name = file_name[:-4]
        self.record = do_wfdb.read_ecg(file_path=self.file_name)
        self.label_4.setText(_translate("Dialog", ''.join(self.record.comments), None))
        
        
    
    def ui_plot(self):
        self.fig.clear()
        record_to_plot = wfdb.rdrecord(self.file_name, channels = [0], sampfrom=self.spinBox.value(), sampto=self.spinBox_2.value())
        self.plot_wfdb(record=record_to_plot,
                    time_units='seconds', return_fig = True)
        self.canvas.draw()  
        
    def calculate(self):
        self.rr_intervals = do_wfdb.compute_rr_interval(self.record)
        self.mean = do_wfdb.compute_mean(self.rr_intervals)
        self.rmssd = do_wfdb.compute_rmssd(self.rr_intervals)
        self.sdnn = do_wfdb.compute_sdnn(self.rr_intervals)
        self.phh50 = do_wfdb.compute_phh50(self.rr_intervals)
        self.cv = do_wfdb.compute_cv(self.rr_intervals)
        self.textBrowser.append("mean: %.4f\nrmssd: %.4f\nsdnn: %.4f\nphh50: %.4f\ncv: %.4f" %(self.mean, self.rmssd, self.sdnn, self.phh50, self.cv))
        

    def plot_items(self, signal=None, ann_samp=None, ann_sym=None, fs=None,
               time_units='samples', sig_name=None, sig_units=None,
               ylabel=None, title=None, sig_style=[''], ann_style=['r*'],
               ecg_grids=[], figsize=None, return_fig=False):
        # Figure out number of subplots required
        sig_len, n_sig, n_annot, n_subplots = self.get_plot_dims(signal, ann_samp)
    
        # Create figure
        axes = self.create_figure(n_subplots, figsize)
    
        if signal is not None:
            self.plot_signal(signal, sig_len, n_sig, fs, time_units, sig_style, axes)

        if ecg_grids:
            self.plot_ecg_grids(ecg_grids, fs, sig_units, time_units, axes)
    
        # Add title and axis labels.
        self.label_figure(axes, n_subplots, time_units, sig_name, sig_units, ylabel, title)
    
    
    def get_plot_dims(self, signal, ann_samp):
        "Figure out the number of plot channels"
        if signal is not None:
            if signal.ndim == 1:
                sig_len = len(signal)
                n_sig = 1
            else:
                sig_len = signal.shape[0]
                n_sig = signal.shape[1]
        else:
            sig_len = 0
            n_sig = 0

        if ann_samp is not None:
            n_annot = len(ann_samp)
        else:
            n_annot = 0

        return sig_len, n_sig, n_annot, max(n_sig, n_annot)
    
    
    def create_figure(self, n_subplots, figsize):
        axes = []

        for i in range(n_subplots):
            axes.append(self.fig.add_subplot(n_subplots, 1, i+1))

        return axes
        
        
    def plot_signal(self, signal, sig_len, n_sig, fs, time_units, sig_style, axes):
        "Plot signal channels"

        # Extend signal style if necesary
        if len(sig_style) == 1:
            sig_style = n_sig * sig_style

        # Figure out time indices
        if time_units == 'samples':
            t = np.linspace(0, sig_len-1, sig_len)
        else:
            downsample_factor = {'seconds':fs, 'minutes':fs * 60,
                                'hours':fs * 3600}
            t = np.linspace(0, sig_len-1, sig_len) / downsample_factor[time_units]

        # Plot the signals
        if signal.ndim == 1:
            axes[0].plot(t, signal, sig_style[0], zorder=3)
        else:
            for ch in range(n_sig):
                axes[ch].plot(t, signal[:,ch], sig_style[ch], zorder=3)
                
                    
    def plot_ecg_grids(self, ecg_grids, fs, units, time_units, axes):
        "Add ecg grids to the axes"
        if ecg_grids == 'all':
            ecg_grids = range(0, len(axes))

    
        for ch in ecg_grids:
            # Get the initial plot limits
            auto_xlims = axes[ch].get_xlim()
            auto_ylims= axes[ch].get_ylim()
    
            (major_ticks_x, minor_ticks_x, major_ticks_y,
                minor_ticks_y) = self.calc_ecg_grids(auto_ylims[0], auto_ylims[1],
                                            units[ch], fs, auto_xlims[1],
                                                    time_units)

            min_x, max_x = np.min(minor_ticks_x), np.max(minor_ticks_x)
            min_y, max_y = np.min(minor_ticks_y), np.max(minor_ticks_y)

            for tick in minor_ticks_x:
                axes[ch].plot([tick, tick], [min_y,  max_y], c='#ededed',
                            marker='|', zorder=1)
            for tick in major_ticks_x:
                axes[ch].plot([tick, tick], [min_y, max_y], c='#bababa',
                            marker='|', zorder=2)
            for tick in minor_ticks_y:
                axes[ch].plot([min_x, max_x], [tick, tick], c='#ededed',
                            marker='_', zorder=1)
            for tick in major_ticks_y:
                axes[ch].plot([min_x, max_x], [tick, tick], c='#bababa',
                            marker='_', zorder=2)

            # Plotting the lines changes the graph. Set the limits back
            axes[ch].set_xlim(auto_xlims)
            axes[ch].set_ylim(auto_ylims)
            
    def calc_ecg_grids(self, minsig, maxsig, sig_units, fs, maxt, time_units):
        # Get the grid interval of the x axis
        if time_units == 'samples':
            majorx = 0.2 * fs
            minorx = 0.04 * fs
        elif time_units == 'seconds':
            majorx = 0.2
            minorx = 0.04
        elif time_units == 'minutes':
            majorx = 0.2 / 60
            minorx = 0.04/60
        elif time_units == 'hours':
            majorx = 0.2 / 3600
            minorx = 0.04 / 3600

        # Get the grid interval of the y axis
        if sig_units.lower()=='uv':
            majory = 500
            minory = 125
        elif sig_units.lower()=='mv':
            majory = 0.5
            minory = 0.125
        elif sig_units.lower()=='v':
            majory = 0.0005
            minory = 0.000125
        else:
            raise ValueError('Signal units must be uV, mV, or V to plot ECG grids.')

        major_ticks_x = np.arange(0, upround(maxt, majorx) + 0.0001, majorx)
        minor_ticks_x = np.arange(0, upround(maxt, majorx) + 0.0001, minorx)

        major_ticks_y = np.arange(downround(minsig, majory),
                                upround(maxsig, majory) + 0.0001, majory)
        minor_ticks_y = np.arange(downround(minsig, majory),
                              upround(maxsig, majory) + 0.0001, minory)

        return (major_ticks_x, minor_ticks_x, major_ticks_y, minor_ticks_y)
        
    
    def label_figure(self, axes, n_subplots, time_units, sig_name, sig_units, ylabel,
                 title):
        if title:
            axes[0].set_title(title)

        # Determine y label
        # Explicit labels take precedence if present. Otherwise, construct labels
        # using signal names and units
        if not ylabel:
            ylabel = []
            # Set default channel and signal names if needed
            if not sig_name:
                sig_name = ['ch_'+str(i) for i in range(n_subplots)]
            if not sig_units:
                sig_units = n_subplots * ['NU']

            ylabel = ['/'.join(pair) for pair in zip(sig_name, sig_units)]

            # If there are annotations with channels outside of signal range
            # put placeholders
            n_missing_labels = n_subplots - len(ylabel)
            if n_missing_labels:
                ylabel = ylabel + ['ch_%d/NU' % i for i in range(len(ylabel),
                                                             n_subplots)]

        for ch in range(n_subplots):
            axes[ch].set_ylabel(ylabel[ch])

        axes[-1].set_xlabel('/'.join(['time', time_units[:-1]]))
        
        
    def plot_wfdb(self, record=None, annotation=None, plot_sym=False,
            time_units='samples', title=None, sig_style=[''],
            ann_style=['r*'], ecg_grids=[], figsize=None, return_fig=False):
        (signal, ann_samp, ann_sym, fs,
        ylabel, record_name) = self.get_wfdb_plot_items(record=record,
                                                   annotation=annotation,
                                                   plot_sym=plot_sym)

        return self.plot_items(signal=signal, ann_samp=ann_samp, ann_sym=ann_sym, fs=fs,
                        time_units=time_units, ylabel=ylabel,
                        title=(title or record_name),
                        sig_style=sig_style,
                        ann_style=ann_style, ecg_grids=ecg_grids,
                        figsize=figsize, return_fig=return_fig)
                        
                        
    def get_wfdb_plot_items(self, record, annotation, plot_sym):
        """
        Get items to plot from wfdb objects
        """
        # Get record attributes
        if record:
            if record.p_signal is not None:
                signal = record.p_signal
            elif record.d_signal is not None:
                signal = record.d_signal
            else:
                raise ValueError('The record has no signal to plot')

            fs = record.fs
            sig_name = record.sig_name
            sig_units = record.units
            record_name = 'Record: %s' % record.record_name
            ylabel = ['/'.join(pair) for pair in zip(sig_name, sig_units)]
        else:
            signal = fs = ylabel = record_name = None

        # Get annotation attributes
        if annotation:
            # Get channels
            ann_chans = set(annotation.chan)
            n_ann_chans = max(ann_chans) + 1

            # Indices for each channel
            chan_inds = n_ann_chans * [np.empty(0, dtype='int')]

            for chan in ann_chans:
                chan_inds[chan] = np.where(annotation.chan == chan)[0]

            ann_samp = [annotation.sample[ci] for ci in chan_inds]

            if plot_sym:
                ann_sym = n_ann_chans * [None]
                for ch in ann_chans:
                    ann_sym[ch] = [annotation.symbol[ci] for ci in chan_inds[ch]]
            else:
                ann_sym = None
    
            # Try to get fs from annotation if not already in record
            if fs is None:
                fs = annotation.fs

            record_name = record_name or annotation.record_name
        else:
            ann_samp = None
            ann_sym = None
        
        # Cleaning: remove empty channels and set labels and styles.

        return signal, ann_samp, ann_sym, fs, ylabel, record_name