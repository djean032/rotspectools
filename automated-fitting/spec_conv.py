from file_stitcher import load_spectrum, file_stitcher
from scipy import signal
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import scienceplots


mpl.rc("text", usetex=True)
plt.style.use(["science", "ieee"])

path1 = "/Users/dairenjean/OneDrive/Projects/test_ROT/2022-09-30-12_08_(3-cyano)methylenecyclopropane_235-345GHz_12mtorr.spe"
path2 = "/Users/dairenjean/OneDrive/Projects/test_ROT/2022-09-27-10_10_(3-cyano)methylenecyclopropane_340-500GHz_12mtorr.spe"
x_low, y_low = load_spectrum(path1)
x_up, y_up = load_spectrum(path2)
x_stitched, y_stitched = file_stitcher(x_low, y_low, x_up, y_up, 345000)

path3 = "/Users/dairenjean/OneDrive/Projects/automated-fitting/ketene.spe"
path4 = "/Users/dairenjean/OneDrive/Projects/automated-fitting/ketene_d2.spe"
path5 = "/Users/dairenjean/OneDrive/Projects/automated-fitting/2022-12-02-10_48_HD_ketene_SW_e037_4_mtorr_350-500_GHz.spe"
x, y = load_spectrum(path3)
x2, y2 = load_spectrum(path4)
x3, y3 = load_spectrum(path5)

low_range = 350000
high_range = 500000

def prepare_spectrum(y):
    y_baseline = signal.savgol_filter(y, 211, 7)
    y -= y_baseline
    return y

y = y[np.where((x > low_range) & (x < high_range))[0]]
y2 = y2[np.where((x2 > low_range) & (x2 < high_range))[0]]
y3 = y3[np.where((x3 > low_range) & (x3 < high_range))[0]]

y = prepare_spectrum(y)
y2 = prepare_spectrum(y2)
y3 = prepare_spectrum(y3)

y_norm = np.linalg.norm(y)
y2_norm = np.linalg.norm(y2)
y3_norm = np.linalg.norm(y3)

y = y / y_norm
y2 = y2 / y2_norm
y3 = y3 / y3_norm

y_stitched = y_stitched[np.where((x_stitched > low_range) & (x_stitched < high_range))[0]]
y_stitched = prepare_spectrum(y_stitched)
y_stitched_norm = np.linalg.norm(y_stitched)
y_stitched = y_stitched / y_stitched_norm

autocorr = signal.fftconvolve(y2, y3[::-1], mode='full')
autocorr2 = signal.fftconvolve(y, y2[::-1], mode='full')
autocorr3 = signal.fftconvolve(y, y3[::-1], mode='full')
hh = signal.fftconvolve(y, y[::-1], mode='full')
hd = signal.fftconvolve(y3, y3[::-1], mode='full')
dd = signal.fftconvolve(y2, y2[::-1], mode='full')
ref = signal.fftconvolve(y, y_stitched[::-1], mode='full')
ref2 = signal.fftconvolve(y2, y_stitched[::-1], mode='full')
ref3 = signal.fftconvolve(y3, y_stitched[::-1], mode='full')

print("D2 and HD")
print(np.max(autocorr))
print("H2 and D2")
print(np.max(autocorr2))
print("H2 and HD")
print(np.max(autocorr3))

print("Correlation w/ other signal")
print(np.max(ref), np.max(ref2), np.max(ref3))
print("Autocorrelations")
print(np.max(hh), np.max(hd), np.max(dd))

#plt.plot(np.arange(-len(autocorr)/2, len(autocorr)/2), autocorr, linewidth=0.4)
#plt.plot(np.arange(-len(autocorr2)/2, len(autocorr2)/2), autocorr2, linewidth=0.4)
plt.plot(np.arange(-len(autocorr3)/2, len(autocorr3)/2), autocorr3, linewidth=0.4, alpha=0.5)
plt.plot(np.arange(-len(ref)/2, len(ref)/2), ref, linewidth=0.4, alpha=0.5)
plt.savefig("ketene_correlations_part_deuterated.png")


plt.cla()

plt.plot(np.arange(-len(hh)/2, len(hh)/2), hh, linewidth=0.5)
plt.savefig("autocorr.png")
