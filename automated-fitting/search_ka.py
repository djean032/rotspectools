import itertools
import pickle
from tkinter import Tk
from tkinter.filedialog import askopenfilename

import numpy as np
import pandas as pd
from pandas import DataFrame

import peak_finder


def read_cat() -> DataFrame:
    Tk().withdraw()
    filename = askopenfilename()
    dataframe: DataFrame = pd.read_fwf(  # type: ignore
        filename,
        header=None,
        names=[
            "Frequency",
            "Error",
            "Integrated Intensity",
            "Deg of Freedom",
            "Lower State Energy (cm-1)",
            "Upper State Deg",
            "Tag",
            "QNFMT",
            "Upper N",
            "Upper Ka",
            "Upper Kc",
            "Upper V",
            "Upper J",
            "Upper F",
            "Lower N",
            "Lower Ka",
            "Lower Kc",
            "Lower V",
            "Lower J",
            "Lower F",
        ],
        widths=[13, 8, 8, 2, 10, 3, 7, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    )
    return dataframe


def map_strings_to_numeric(df: DataFrame):
    df = df.astype(str)
    df.drop(
        columns=["Upper J", "Upper F", "Lower J", "Lower F"],
        inplace=True,
    )
    quantum_list = [
        "Upper N",
        "Upper Ka",
        "Upper Kc",
        "Upper V",
        "Lower N",
        "Lower Ka",
        "Lower Kc",
        "Lower V",
    ]
    for column_name in quantum_list:
        column = df[column_name]
        numeric_values = []

        for value in column:
            if value.isdigit():
                numeric_values.append(value)
                continue

            if pd.isnull(value):  # Skip empty values
                numeric_values.append(value)
                continue

            value = str(value)
            letter = value[0]  # Extract the first character (letter)
            number = int(value[1:])  # Extract the remaining characters as the number
            letter_value = (
                ord(letter.upper()) - 65
            )  # Convert letter to value (A=0, B=1, ...)
            if letter.islower():  # Check if the letter is lowercase
                numeric_value = (
                    -10 * (letter_value + 1) - number
                )  # Calculate the final numeric value for lowercase letter
            else:
                numeric_value = (
                    100 + letter_value * 10 + number
                )  # Calculate the final numeric value for uppercase letter
            numeric_values.append(numeric_value)

        df[column_name] = numeric_values
    df = df.apply(pd.to_numeric)
    return df


def categorize_trans_type(df: DataFrame, column_name1: str, column_name2: str):
    column1 = df[column_name1]
    column2 = df[column_name2]
    transition_types = []

    for Ka, Kc in zip(column1, column2):
        if Kc % 2 == 0:
            transition_types.append("c-type")
        elif Ka % 2 == 0:
            transition_types.append("a-type")
        else:
            transition_types.append("b-type")
    df["Transition Type"] = transition_types
    return df


def categorize_branch(df: DataFrame, column_name: str):
    column = df[column_name]
    branches = []

    for value in column:
        if value == 1:
            branches.append("R-Branch")
        elif value == -1:
            branches.append("P-Branch")
        else:
            branches.append("Q-Branch")
    df["Branch"] = branches
    return df


def qn_delta(df: DataFrame) -> DataFrame:
    df["Delta N"] = df["Upper N"] - df["Lower N"]
    df["Delta Ka"] = df["Upper Ka"] - df["Lower Ka"]
    df["Delta Kc"] = df["Upper Kc"] - df["Lower Kc"]
    return df


def find_Ka_poss(peaks, heights, upper, lower, range):
    fil_mask_peaks = np.where((peaks >= lower) & (peaks <= upper))
    filt_x_peaks = peaks[fil_mask_peaks]
    filt_y_peaks = heights[fil_mask_peaks]

    max_loc = np.argsort(filt_y_peaks.astype(float))[::-1]
    filt_x_peaks = filt_x_peaks[max_loc]
    peaks_1 = filt_x_peaks[0:5]
    filt_y_peaks = filt_y_peaks[max_loc]

    fil_mask_peaks2 = np.where((peaks >= upper) & (peaks <= (upper + range)))
    filt_x_peaks2 = peaks[fil_mask_peaks2]
    filt_y_peaks2 = heights[fil_mask_peaks2]

    max_loc2 = np.argsort(filt_y_peaks2.astype(float))[::-1]
    filt_x_peaks2 = filt_x_peaks2[max_loc2]
    peaks_2 = filt_x_peaks2[0:5]
    filt_y_peaks2 = filt_y_peaks2[max_loc2]

    fil_mask_peaks3 = np.where(
        (peaks >= (upper + range)) & (peaks <= (upper + 2 * range))
    )
    filt_x_peaks3 = peaks[fil_mask_peaks3]
    filt_y_peaks3 = heights[fil_mask_peaks3]

    max_loc3 = np.argsort(filt_y_peaks3.astype(float))[::-1]
    filt_x_peaks3 = filt_x_peaks3[max_loc3]
    peaks_3 = filt_x_peaks3[0:5]
    filt_y_peaks3 = filt_y_peaks3[max_loc3]

    combination = [p for p in itertools.product(peaks_1, peaks_2, peaks_3)]

    possibility = []
    for i in combination:
        i = list(i)
        difference = i[1] - i[0]
        i.append(difference)
        difference2 = i[2] - i[1]
        i.append(difference2)
        difference3 = difference2 - difference
        i.append(difference3)
        ind = np.where(i[0] == filt_x_peaks)
        ind2 = np.where(i[1] == filt_x_peaks2)
        ind3 = np.where(i[2] == filt_x_peaks3)
        int_sum = filt_y_peaks[ind][0] + filt_y_peaks2[ind2][0] + filt_y_peaks3[ind3][0]
        i.append(int_sum)
        if difference < (range + 120) and difference > (range - 120):
            if difference2 < (range + 120) and difference2 > (range - 120):
                possibility.append(i)
    possibility = np.array(possibility)
    possibility = possibility[possibility[:, 6].argsort()[::-1]]
    strong_possibility = possibility[0:5]
    return strong_possibility


def find_strong_pred(pred_spectrum: DataFrame, upper: int, lower: int) -> DataFrame:
    filt_mask = np.where(
        (pred_spectrum["Frequency"] >= lower) & (pred_spectrum["Frequency"] <= upper)
    )
    filt_pred_spectrum = pred_spectrum.loc[filt_mask]
    filt_pred_spectrum = filt_pred_spectrum.sort_values(
        by=["Integrated Intensity"], ascending=False
    )
    return filt_pred_spectrum


def construct_line_df(lin_dataframe):
    complete_lin_dataframe = lin_dataframe.copy()
    complete_lin_dataframe.drop(["Frequency"], axis=1, inplace=True)

    complete_lin_dataframe["Filler 1"] = 0
    complete_lin_dataframe["Filler 2"] = 0
    complete_lin_dataframe["Filler 3"] = 0
    complete_lin_dataframe["Filler 4"] = 0

    complete_lin_dataframe["Error"] = 0.050000
    complete_lin_dataframe["Weight"] = 1.000000

    complete_lin_dataframe.sort_values(by=["Assigned Frequency"], inplace=True)
    complete_lin_dataframe[
        [
            "Upper N",
            "Upper Ka",
            "Upper Kc",
            "Upper V",
            "Lower N",
            "Lower Ka",
            "Lower Kc",
            "Lower V",
            "Filler 1",
            "Filler 2",
            "Filler 3",
            "Filler 4",
            "Assigned Frequency",
            "Error",
            "Weight",
        ]
    ]
    return complete_lin_dataframe


def convert_to_fwf(df, fname):
    line_list = []
    for index, row in df.iterrows():
        upper_n = f'{int(row["Upper N"]):3}'
        upper_ka = f'{int(row["Upper Ka"]):3}'
        upper_kc = f'{int(row["Upper Kc"]):3}'
        upper_v = f'{int(row["Upper V"]):3}'
        lower_n = f'{int(row["Lower N"]):3}'
        lower_ka = f'{int(row["Lower Ka"]):3}'
        lower_kc = f'{int(row["Lower Kc"]):3}'
        lower_v = f'{int(row["Lower V"]):3}'
        filler_1 = f'{int(row["Filler 1"]):3}'
        filler_2 = f'{int(row["Filler 2"]):3}'
        filler_3 = f'{int(row["Filler 3"]):3}'
        filler_4 = f'{int(row["Filler 4"]):3}'
        assigned_freq = f'{row["Assigned Frequency"]:16.6f}'
        error = f'{row["Error"]:13.6f}'
        weight = f'{row["Weight"]:10.6f}'
        line = (
            upper_n
            + upper_ka
            + upper_kc
            + upper_v
            + lower_n
            + lower_ka
            + lower_kc
            + lower_v
            + filler_1
            + filler_2
            + filler_3
            + filler_4
            + assigned_freq
            + error
            + weight
        )
        line_list.append(line)
    with open(fname, "w") as f:
        f.write("\n".join(line_list))


if __name__ == "__main__":
    path1 = "/users/dairenjean/onedrive/projects/test_rot/2022-09-30-12_08_(3-cyano)methylenecyclopropane_235-345ghz_12mtorr.spe"
    path2 = "/users/dairenjean/onedrive/projects/test_rot/2022-09-27-10_10_(3-cyano)methylenecyclopropane_340-500ghz_12mtorr.spe"
    x_low, y_low = peak_finder.load_spectrum(path1)
    x_up, y_up = peak_finder.load_spectrum(path2)
    x_spec, y_spec = peak_finder.file_stitcher(x_low, y_low, x_up, y_up, 340000)
    y_spec = peak_finder.prepare_spectrum(y_spec)
    #    chunks = peak_finder.chunk_spectrum(x_spec, y_spec)
    #    peaks, heights = peak_finder.find_peaks_mp(chunks)

    with open("peaks.pkl", "rb") as f:
        peaks = pickle.load(f)

    with open("heights.pkl", "rb") as f:
        heights = pickle.load(f)

    peaks = np.array(peaks, dtype=np.float64)
    heights = np.array(heights, dtype=np.float64)

    exp_spectrum = pd.DataFrame(x_spec, y_spec)
    pred_spectrum = read_cat()
    pred_spectrum = map_strings_to_numeric(pred_spectrum)
    qn_delta(pred_spectrum)
    categorize_branch(pred_spectrum, "Delta N")
    categorize_trans_type(pred_spectrum, "Delta Ka", "Delta Kc")

    filt_pred_spectrum = find_strong_pred(pred_spectrum, 399000, 395000)
    filt_pred_spectrum2 = find_strong_pred(pred_spectrum, 404000, 399000)
    filt_pred_spectrum3 = find_strong_pred(pred_spectrum, 408000, 404000)
    ka0_1 = filt_pred_spectrum.loc[
        (filt_pred_spectrum["Upper Ka"] == 0) & (filt_pred_spectrum["Delta Ka"] == 0)
    ]
    ka0_2 = filt_pred_spectrum2.loc[
        (filt_pred_spectrum2["Upper Ka"] == 0) & (filt_pred_spectrum2["Delta Ka"] == 0)
    ]
    ka0_3 = filt_pred_spectrum3.loc[
        (filt_pred_spectrum3["Upper Ka"] == 0) & (filt_pred_spectrum3["Delta Ka"] == 0)
    ]
    ka0_series_high = pd.concat([ka0_1, ka0_2, ka0_3])

    strong_possibility = find_Ka_poss(peaks, heights, 399000, 395000, 4200)
    print(strong_possibility)
    ka0_series_high["Assigned Frequency"] = strong_possibility[0][0:3]
    high_lin_dataframe = construct_line_df(ka0_series_high)

    filt_pred_spectrum_low = find_strong_pred(pred_spectrum, 275000, 271000)
    filt_pred_spectrum2_low = find_strong_pred(pred_spectrum, 279000, 275000)
    filt_pred_spectrum3_low = find_strong_pred(pred_spectrum, 283000, 279000)

    ka0_1_low = filt_pred_spectrum_low.loc[
        (filt_pred_spectrum_low["Upper Ka"] == 0)
        & (filt_pred_spectrum_low["Delta Ka"] == 0)
    ]

    ka0_2_low = filt_pred_spectrum2_low.loc[
        (filt_pred_spectrum2_low["Upper Ka"] == 0)
        & (filt_pred_spectrum2_low["Delta Ka"] == 0)
    ]

    ka0_3_low = filt_pred_spectrum3_low.loc[
        (filt_pred_spectrum3_low["Upper Ka"] == 0)
        & (filt_pred_spectrum3_low["Delta Ka"] == 0)
    ]

    ka0_series_low = pd.concat([ka0_1_low, ka0_2_low, ka0_3_low])

    strong_possibility_low = find_Ka_poss(peaks, heights, 275000, 271000, 4200)
    print(strong_possibility_low)
    ka0_series_low["Assigned Frequency"] = strong_possibility_low[0][0:3]
    low_lin_dataframe = construct_line_df(ka0_series_low)

    final_df = pd.concat([high_lin_dataframe, low_lin_dataframe])

    convert_to_fwf(final_df, "ka0_series.txt")
