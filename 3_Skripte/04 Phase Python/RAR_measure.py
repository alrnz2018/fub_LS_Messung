# -*- coding: utf-8 -*-
"""
Created on Wed May 22 11:50:29 2019

@author: flau
"""

import numpy as np
from scipy.io import wavfile
from scipy.signal import deconvolve
import matplotlib.pyplot as plt



fs, measure = wavfile.read('Messung.wav')

fs, test = wavfile.read('Testsignal.wav')

difference=len(measure)-len(test)
correlation=np.correlate(test,measure,"full")

maxsample=len(measure)-(np.ceil(np.argmax(correlation)))
zeros=np.zeros(int(maxsample))
measure=np.concatenate([zeros,measure])


test=test[:144000]
measure=measure[:144000]
spectrum_m = np.fft.fft(measure)
magnitude_m = np.abs(spectrum_m)
phase_m = np.angle(spectrum_m,deg=True)

spectrum_t = np.fft.fft(test)
magnitude_t = np.abs(spectrum_t)
phase_t = np.angle(spectrum_t,deg=True)




spectrum_f=spectrum_m/spectrum_t
len2=np.ceil(len(spectrum_f)/2)
spectrum_f=spectrum_f[:int(len2)]
magnitude_f = np.abs(spectrum_f)
phase_f = np.angle(spectrum_f,deg=True)

impulsantwort=np.fft.ifft(spectrum_f)


plt.plot(phase_f)
plt.xscale('log')
#plt.show