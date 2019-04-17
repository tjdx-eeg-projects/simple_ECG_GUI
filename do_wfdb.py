# -*- coding: utf-8 -*-
import numpy as np
from wfdb import io
from wfdb import plot
from wfdb import processing


def read_ecg(file_path):
    #读取ecg信号
    record = io.rdrecord(file_path)
    return record

'''
def plot_ecg(file_path, channels = 0, sampfrom=0, sampto=None):
    #画出指定channel的一段样本从sampfrom到sampto的ecg图像
    record_to_plot = io.rdrecord(file_path, channels = [channels], sampfrom=sampfrom, sampto=sampto)
    fig = plot.plot_wfdb(record=record_to_plot,
                time_units='seconds', return_fig = True)
    return fig
'''
               
def compute_rr_interval(record, channels = 0):
    #计算rr间期，返回包含所有间期的一维numpy数组
    qrs_inds = processing.xqrs_detect(sig=record.p_signal[:,channels], fs=record.fs)
    '''
    i = 1
    rr_interval = []
    while i < qrs_inds.size:
        rr_interval.append((qrs_inds[i] - qrs_inds[i-1]) / record.fs)
        i = i+1
    '''
    rr_intervals = processing.calc_rr(qrs_inds, fs=record.fs, rr_units='seconds')
    return rr_intervals

def compute_mean(rr_intervals):
    return np.mean(rr_intervals)

def compute_rmssd(rr_intervals):
    differences = np.diff(rr_intervals)
    return np.std(differences)

def compute_sdnn(rr_intervals):
    return np.std(rr_intervals)
    
def compute_phh50(rr_intervals):
    differences = np.diff(rr_intervals)
    count = 0
    for i in differences:
        if(abs(i) > 0.05):
            count += 1
    return count/rr_intervals.size

def compute_cv(rr_intervals):
    sdnn = np.std(rr_intervals)
    mean = np.mean(rr_intervals)
    return sdnn/mean