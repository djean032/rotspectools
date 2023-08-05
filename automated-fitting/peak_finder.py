import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import multiprocessing as mp
import time
import pickle

path = "/Users/dairenjean/OneDrive/Projects/test_ROT/test_ascii.spe"


def load_spectrum(path):
    x, y = np.loadtxt(
        path,
        usecols=(0, 1),
        skiprows=3,
        unpack=True,
    )
    return x, y


def file_stitcher(x_lower, y_lower, x_upper, y_upper, stitch_freq):
    l_ind = np.where(x_lower < stitch_freq)[0]
    u_ind = np.where(x_upper >= stitch_freq)[0]
    x_l = x_lower[l_ind]
    y_l = y_lower[l_ind]
    x_u = x_upper[u_ind]
    y_u = y_upper[u_ind]
    x = np.concatenate((x_l, x_u))
    y = np.concatenate((y_l, y_u))
    return x, y


def prepare_spectrum(y):
    y_baseline = sp.signal.savgol_filter(y, 211, 7)
    y -= y_baseline
    y = sp.signal.savgol_filter(y, 15, 7)
    return y


def find_peaks(chunks):
    window = 20
    chunk_x, chunk_y = chunks
    chunk_y = -1 * chunk_y
    tent_peaks_index = sp.signal.find_peaks_cwt(chunk_y, np.arange(1, 25), min_snr=2.5)
    actual_peaks_index = []
    actual_peaks_x = []
    heights = []
    for i in tent_peaks_index:
        if i < window:
            min_index = 0
        else:
            min_index = i - window
        if i > len(chunk_y) - window:
            max_index = len(chunk_y) - 1
        else:
            max_index = i + window
        peak_range = chunk_y[min_index:max_index]
        rel_max = peak_range.max()
        rel_min = peak_range.min()
        height = rel_max - rel_min
        rel_max_index = np.where(chunk_y == rel_max)[0][0]
        actual_peaks_index.append(rel_max_index)
        actual_peaks_x = chunk_x[actual_peaks_index]
        heights.append(height)
    heights = np.array(heights)
    return actual_peaks_x, heights


def chunk_spectrum(x, y):
    chunk_size = len(x) / 2900
    chunk_x = np.array_split(x, chunk_size)
    chunk_y = np.array_split(y, chunk_size)
    chunks = zip(chunk_x, chunk_y)
    return chunks


def find_peaks_mp(chunks):
    with mp.Pool() as p:
        peaks = np.array([])
        heights = np.array([])
        for peaks1, heights1 in p.map(find_peaks, chunks):
            peaks = np.append(peaks, peaks1)
            heights = np.append(heights, heights1)
    return peaks, heights


if __name__ == "__main__":
    begin = time.time()
    path1 = '2022-09-30-12_08_(3-cyano)methylenecyclopropane_235-345GHz_12mtorr.spe'
    path2 = '2022-09-27-10_10_(3-cyano)methylenecyclopropane_340-500GHz_12mtorr.spe'
    x_low, y_low = load_spectrum(path1)
    x_up, y_up = load_spectrum(path2)
    x, y = file_stitcher(x_low, y_low, x_up, y_up, 340000)
    y = prepare_spectrum(y)
    chunks = chunk_spectrum(x, y)
    peaks, heights = find_peaks_mp(chunks)
    end = time.time()
    print(end - begin)
    max_ind = np.where(heights == heights.max())[0][0]

    x_plot = np.where((x >= 429000) & (x <= 429908))[0][0]
    window = 450
    x_range = x[x_plot - window : x_plot + window]
    y_range = y[x_plot - window : x_plot + window]
    plt.plot(x_range, y_range, linewidth=0.7, color="black")
    plot_peaks_mask = (peaks >= x_range[0]) & (peaks <= x_range[-1])
    plot_peaks = peaks[plot_peaks_mask]
    plot_heights = heights[plot_peaks_mask]
    print(plot_peaks, plot_heights)
    plt.vlines(
        plot_peaks,
        ymin=min(y_range),
        ymax=max(y_range),
        color="r",
        linestyles="dashed",
        linewidth=0.6,
    )
    plt.show()

    peaks = peaks.tolist()
    heights = heights.tolist()

    with open("peaks.pkl", "wb") as f:
        pickle.dump(peaks, f)

    with open("heights.pkl", "wb") as f:
        pickle.dump(heights, f)
