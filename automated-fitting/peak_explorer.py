# %%
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

with open("peaks.pkl", "rb") as f:
    peaks = pickle.load(f)

with open("heights.pkl", "rb") as f:
    heights = pickle.load(f)

peaks = np.array(peaks, dtype=np.float64)
heights = np.array(heights, dtype=np.float64)
# %%
print(len(peaks))
print(len(heights))
# %%
df = pd.DataFrame(heights)
df.describe()

# %%
intensity_mask = np.where(heights >= 0.04)[0]
# %%
print(len(intensity_mask))
# %%
new_peaks = peaks[intensity_mask]
new_heights = heights[intensity_mask]

df2 = pd.DataFrame(new_heights)
df2.describe()
# %%
fig, ax = plt.subplots()
ax.hist(new_heights, bins=1000)
ax.set_xlim(0, 0.9)
# %%
10 ** (-7.0)
# %%
10 ** (-5)
# %%
