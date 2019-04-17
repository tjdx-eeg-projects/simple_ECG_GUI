# -*- coding: utf-8 -*-
import numpy as np
from wfdb import io
from wfdb import plot
from wfdb import processing


def read_ecg(file_path):
    #read ecg
    record = io.rdrecord(file_path)
    return record
               
def compute_rr_interval(record, channels = 0):
    qrs_inds = processing.xqrs_detect(sig=record.p_signal[:,channels], fs=record.fs)
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
