import rotspectools.rotspectools as rt
import pandas as pd
import search_ka

test = rt.Experiment(5000, "asymmetric_rotor", "cyanomethcycloprop_gs")

test.clean_lin(2.0)


def read_lin(fname):
    lines = pd.read_fwf(
        fname,
        header=None,
        widths=[3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 16, 13, 10],
    )

    lines.rename(
        columns={
            0: "Upper N",
            1: "Upper Ka",
            2: "Upper Kc",
            3: "Upper V",
            4: "Lower N",
            5: "Lower Ka",
            6: "Lower Kc",
            7: "Lower V",
            8: "Filler 1",
            9: "Filler 2",
            10: "Filler 3",
            11: "Filler 4",
            12: "Assigned Frequency",
            13: "Error",
            14: "Weight",
        },
        inplace=True,
    )
    lines.sort_values(by="Assigned Frequency", ascending=True, inplace=True)
    lines.drop(lines[lines["Assigned Frequency"] == 0.0].index, inplace=True)
    lines.drop_duplicates(
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
    return lines


lines = read_lin("cyanomethcycloprop_gs.lin")
search_ka.convert_to_fwf(lines, "cyanomethcycloprop_gs.lin")

# test.plot_data_dist(51, 120)
