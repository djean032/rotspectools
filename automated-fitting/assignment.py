import search_ka as sk
import numpy as np
import pandas as pd
import scipy as sp
import pickle


def assign_ka0(pred, peaks, heights, upper_j, low_freq, high_freq, search_range):
    filt_pred_spectrum = pred.loc[
        (pred["Upper Ka"] == 0)
        & (pred["Lower Ka"] == 0)
        & (pred["Delta Kc"] == 1)
        & (pred["Delta N"] == 1)
        & (pred["Upper N"] <= upper_j)
        & (pred["Frequency"] >= low_freq)
        & (pred["Frequency"] <= high_freq)
    ]

    assigned_list = []
    previous_assigned = 0
    average_diff = [0, 0]
    for index, row in filt_pred_spectrum.iterrows():
        pred_freq = row["Frequency"]
        upper_freq = pred_freq + (search_range / 2)
        lower_freq = pred_freq - (search_range / 2)
        freq_mask = np.where((peaks <= upper_freq) & (peaks >= lower_freq))
        test_peaks = peaks[freq_mask]
        test_heights = heights[freq_mask]
        arrange_heights_mask = np.argsort(test_heights)[::-1]
        arrange_peaks = test_peaks[arrange_heights_mask]
        if filt_pred_spectrum["Frequency"].iloc[0] == row["Frequency"]:
            assigned_freq = arrange_peaks[0]
            previous_assigned = assigned_freq
            assigned_list.append(assigned_freq)
        elif filt_pred_spectrum["Frequency"].iloc[1] == row["Frequency"]:
            assigned_freq = arrange_peaks[0]
            average_diff[0] = assigned_freq - previous_assigned
            previous_assigned = assigned_freq
            assigned_list.append(assigned_freq)
        elif filt_pred_spectrum["Frequency"].iloc[2] == row["Frequency"]:
            assigned_freq = arrange_peaks[0]
            average_diff[1] = assigned_freq - previous_assigned
            previous_assigned = assigned_freq
            assigned_list.append(assigned_freq)
        else:
            m, b, *_ = sp.stats.linregress(assigned_list[1:3], average_diff[0:2])
            assigned_freq = arrange_peaks[0]
            i = 0
            tent_assigned_freq = 0
            difference = assigned_freq - previous_assigned
            extrapolated_diff = m * assigned_freq + b
            lower_limit = extrapolated_diff * 0.99
            upper_limit = extrapolated_diff * 1.01
            while difference < lower_limit or difference > upper_limit:
                if i >= len(arrange_peaks) - 1:
                    break
                tent_assigned_freq = arrange_peaks[i + 1]
                i += 1
                extrapolated_diff = m * tent_assigned_freq + b
                difference = tent_assigned_freq - previous_assigned
                lower_limit = extrapolated_diff * 0.99
                upper_limit = extrapolated_diff * 1.01
            if tent_assigned_freq == 0:
                assigned_freq = arrange_peaks[0]
            else:
                assigned_freq = tent_assigned_freq
            previous_assigned = assigned_freq
            average_diff.append(assigned_freq - previous_assigned)
            assigned_list.append(assigned_freq)

    filt_pred_spectrum["Assigned Frequency"] = assigned_list

    lin_dataframe = pd.DataFrame(
        filt_pred_spectrum,
        columns=[
            "Assigned Frequency",
            "Upper N",
            "Upper Ka",
            "Upper Kc",
            "Upper V",
            "Lower N",
            "Lower Ka",
            "Lower Kc",
            "Lower V",
            "Frequency",
            "Integrated Intensity",
        ],
    )

    return assigned_list, lin_dataframe


def assign_lines_not_ka0(
    pred, peaks, heights, upper_ka, lower_ka, upper_j, search_range, series_variance
):
    filt_pred_spectrum = pred.loc[
        (pred["Upper Ka"] == upper_ka)
        & (pred["Lower Ka"] == lower_ka)
        & (pred["Delta Kc"] == 1)
        & (pred["Delta N"] == 1)
        & (pred["Upper N"] <= upper_j)
        & (pred["Lower Ka"] + pred["Lower Kc"] == pred["Lower N"])
        & (pred["Frequency"] >= 235000)
        & (pred["Frequency"] <= 500000)
    ]

    assigned_list = []
    previous_assigned = 0
    average_diff = [0, 0]
    for index, row in filt_pred_spectrum.iterrows():
        pred_freq = row["Frequency"]
        upper_freq = pred_freq + search_range / 2
        lower_freq = pred_freq - search_range / 2
        freq_mask = np.where((peaks <= upper_freq) & (peaks >= lower_freq))
        test_peaks = peaks[freq_mask]
        test_heights = heights[freq_mask]
        arrange_heights_mask = np.argsort(test_heights)[::-1]
        arrange_peaks = test_peaks[arrange_heights_mask]
        if filt_pred_spectrum["Frequency"].iloc[0] == row["Frequency"]:
            assigned_freq = arrange_peaks[0]
            previous_assigned = assigned_freq
            assigned_list.append(assigned_freq)
        elif filt_pred_spectrum["Frequency"].iloc[1] == row["Frequency"]:
            assigned_freq = arrange_peaks[0]
            average_diff[0] = assigned_freq - previous_assigned
            previous_assigned = assigned_freq
            assigned_list.append(assigned_freq)
        elif filt_pred_spectrum["Frequency"].iloc[2] == row["Frequency"]:
            assigned_freq = arrange_peaks[0]
            average_diff[1] = assigned_freq - previous_assigned
            previous_assigned = assigned_freq
            assigned_list.append(assigned_freq)
        else:
            m, b, *_ = sp.stats.linregress(assigned_list[1:3], average_diff[0:2])
            low_variance = 1 - (series_variance / 2)
            upper_variance = 1 + (series_variance / 2)
            if len(arrange_peaks) != 0:
                assigned_freq = arrange_peaks[0]
                i = 0
                tent_assigned_freq = 0
                difference = assigned_freq - previous_assigned
                extrapolated_diff = m * assigned_freq + b
                lower_limit = extrapolated_diff * low_variance
                upper_limit = extrapolated_diff * upper_variance
                while difference < lower_limit or difference > upper_limit:
                    if i >= len(arrange_peaks) - 1:
                        break
                    tent_assigned_freq = arrange_peaks[i + 1]
                    i += 1
                    extrapolated_diff = m * tent_assigned_freq + b
                    difference = tent_assigned_freq - previous_assigned
                    lower_limit = extrapolated_diff * low_variance
                    upper_limit = extrapolated_diff * upper_variance
                if tent_assigned_freq == 0:
                    assigned_freq = arrange_peaks[0]
                else:
                    assigned_freq = tent_assigned_freq
                previous_assigned = assigned_freq
                average_diff.append(assigned_freq - previous_assigned)
                assigned_list.append(assigned_freq)
            else:
                assigned_freq = 0
                assigned_list.append(assigned_freq)

    filt_pred_spectrum["Assigned Frequency"] = assigned_list

    lin_dataframe = pd.DataFrame(
        filt_pred_spectrum,
        columns=[
            "Assigned Frequency",
            "Upper N",
            "Upper Ka",
            "Upper Kc",
            "Upper V",
            "Lower N",
            "Lower Ka",
            "Lower Kc",
            "Lower V",
            "Frequency",
            "Integrated Intensity",
        ],
    )

    return assigned_list, lin_dataframe


def assign_lines_after_pred(pred, peaks, heights, upper_ka, lower_ka):
    filt_pred_spectrum = pred.loc[
        (pred["Frequency"] >= 235000)
        & (pred["Frequency"] <= 500000)
        # (pred["Upper Ka"] == upper_ka)
        & (pred["Lower Ka"] == lower_ka)
        # & (pred["Delta Ka"] == 1)
        # & (pred["Delta Kc"] == -1)
        & (pred["Delta N"] == 0)
        # & (pred["Upper N"] <= 100)
        # & (pred["Lower Ka"] + pred["Lower Kc"] == pred["Lower N"])
        & (pred["Integrated Intensity"] >= -6.8)
    ]

    assigned_list = []
    for index, row in filt_pred_spectrum.iterrows():
        pred_freq = row["Frequency"]
        upper_freq = pred_freq + 0.5
        lower_freq = pred_freq - 0.5
        freq_mask = np.where((peaks <= upper_freq) & (peaks >= lower_freq))
        test_peaks = peaks[freq_mask]
        test_heights = heights[freq_mask]
        arrange_heights_mask = np.argsort(test_heights)[::-1]
        arrange_peaks = test_peaks[arrange_heights_mask]
        if len(arrange_peaks) != 0:
            assigned_freq = arrange_peaks[0]
        else:
            assigned_freq = 0
        assigned_list.append(assigned_freq)

    filt_pred_spectrum["Assigned Frequency"] = assigned_list

    lin_dataframe = pd.DataFrame(
        filt_pred_spectrum,
        columns=[
            "Assigned Frequency",
            "Upper N",
            "Upper Ka",
            "Upper Kc",
            "Upper V",
            "Lower N",
            "Lower Ka",
            "Lower Kc",
            "Lower V",
            "Frequency",
            "Integrated Intensity",
        ],
    )

    return assigned_list, lin_dataframe


def find_degenerates(prediction, lin_dataframe):
    additional_lines = pd.Series.isin(
        prediction["Frequency"], lin_dataframe["Frequency"]
    )
    additional_lines = prediction.loc[additional_lines]
    additional_lines_dataframe = pd.DataFrame(
        additional_lines,
        columns=[
            "Upper N",
            "Upper Ka",
            "Upper Kc",
            "Upper V",
            "Lower N",
            "Lower Ka",
            "Lower Kc",
            "Lower V",
            "Frequency",
            "Integrated Intensity",
        ],
    )
    return additional_lines_dataframe


def construct_line_df(assigned_list, lin_dataframe, additional_lines_dataframe):
    complete_lin_dataframe = pd.concat(
        [lin_dataframe, additional_lines_dataframe],
        axis=0,
    )
    freqs = lin_dataframe["Frequency"].to_list()
    freq_assign_dict = dict(zip(freqs, assigned_list))

    complete_lin_dataframe["Assigned Frequency"] = complete_lin_dataframe[
        "Frequency"
    ].map(freq_assign_dict)
    complete_lin_dataframe.drop(["Frequency"], axis=1, inplace=True)
    counts = complete_lin_dataframe.value_counts("Assigned Frequency")
    count_assign_dict = counts.to_dict()
    complete_lin_dataframe["Count"] = complete_lin_dataframe["Assigned Frequency"].map(
        count_assign_dict
    )
    test = complete_lin_dataframe.copy()
    test.set_index("Assigned Frequency", inplace=True)
    test_dict = (
        test.groupby(["Assigned Frequency"])["Integrated Intensity"].min().to_dict()
    )
    complete_lin_dataframe["Max Int"] = complete_lin_dataframe[
        "Assigned Frequency"
    ].map(test_dict)
    complete_lin_dataframe["Integrated Intensity"] = (
        complete_lin_dataframe["Integrated Intensity"].abs()
        / complete_lin_dataframe["Max Int"].abs()
    )

    complete_lin_dataframe["Filler 1"] = 0
    complete_lin_dataframe["Filler 2"] = 0
    complete_lin_dataframe["Filler 3"] = 0
    complete_lin_dataframe["Filler 4"] = 0

    complete_lin_dataframe["Error"] = 0.050000
    complete_lin_dataframe["Weight"] = (
        1.0 / complete_lin_dataframe["Count"]
    ) * complete_lin_dataframe["Integrated Intensity"]
    complete_lin_dataframe.sort_values(by=["Assigned Frequency"], inplace=True)
    complete_lin_dataframe.drop(
        complete_lin_dataframe[
            complete_lin_dataframe["Assigned Frequency"] == 0.0
        ].index,
        inplace=True,
    )
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
    complete_lin_dataframe.drop_duplicates(
        subset=[
            "Upper N",
            "Upper Ka",
            "Upper Kc",
            "Upper V",
            "Lower N",
            "Lower Ka",
            "Lower Kc",
            "Lower V",
        ],
        keep="last",
        inplace=True,
    )
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


pred = sk.read_cat()
pred = sk.map_strings_to_numeric(pred)
sk.qn_delta(pred)
sk.categorize_branch(pred, "Delta N")
sk.categorize_trans_type(pred, "Delta Ka", "Delta Kc")

with open("peaks.pkl", "rb") as f:
    peaks = pickle.load(f)

with open("heights.pkl", "rb") as f:
    heights = pickle.load(f)

peaks = np.array(peaks, dtype=np.float64)
heights = np.array(heights, dtype=np.float64)

intensity_mask = np.where(heights > 0.09)[0]

peaks = peaks[intensity_mask]
heights = heights[intensity_mask]

remaining_heights = heights
remaining_peaks = peaks

list_of_lines = []
assigned_list = []

for i in range(0, 70, 1):
    remaining_peaks = remaining_peaks[
        np.isin(remaining_peaks, assigned_list, invert=True, assume_unique=True)
    ]
    remaining_heights = remaining_heights[
        np.isin(remaining_heights, assigned_list, invert=True, assume_unique=True)
    ]
    assigned_list = []

    assigned_list, lin_dataframe = assign_lines_after_pred(
        pred, remaining_peaks, remaining_heights, i, i
    )

    additional_lines_dataframe = find_degenerates(pred, lin_dataframe)

    complete_lin_dataframe = construct_line_df(
        assigned_list, lin_dataframe, additional_lines_dataframe
    )
    print(len(complete_lin_dataframe))

    list_of_lines.append(complete_lin_dataframe)
    print("ka" + str(i) + " complete")

complete_lin_dataframe = pd.concat(list_of_lines)

convert_to_fwf(complete_lin_dataframe, "0_65.fwf")
