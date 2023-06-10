import math
from pathlib import Path
from typing import Tuple

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scienceplots
from matplotlib.legend_handler import HandlerPatch
from pandas import DataFrame

mpl.rc("text", usetex=True)
plt.style.use(["science", "ieee"])


class Experiment:
    def __init__(
        self,
        num_lines: int,
        molecule_type: str,
        file_name: str,
        file_path_res: str = "", 
        file_path_lin: str = "",
        res_dataframe: DataFrame = pd.DataFrame(),
        lin_dataframe: DataFrame = pd.DataFrame(),
    ):
        self.num_lines = num_lines
        self.molecule_type = molecule_type
        self.file_name = file_name
        self.file_path_res = file_path_res
        self.file_path_lin = file_path_lin
        self.res_dataframe = res_dataframe
        self.lin_dataframe = lin_dataframe
        self.rel_file_path()
        self.read_res()
        self.df_error = self.res_dataframe.copy(deep=True)
        self.read_lin()

    def rel_file_path(self) -> str:
        base_path = Path(__file__).parent
        file_name_res = self.file_name + ".res"
        file_name_lin = self.file_name + ".lin"
        self.file_path_res = (base_path / file_name_res).resolve()
        self.file_path_lin = (base_path / file_name_lin).resolve()

    def read_lin(self) -> DataFrame:
        line_width = len(str(self.num_lines)) + 1
        self.lin_dataframe = pd.read_fwf(
            self.file_path_lin,
            header=None,
            skiprows=7,
            names=[
                "J",
                "Ka",
                "Kc",
                "v",
                "J'",
                "Ka'",
                "Kc'",
                "v'",
                "01",
                "02",
                "03",
                "04",
                "Measured Frequency",
                "Max Error",
                "Relative Weight",
            ],
            widths=[line_width, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 16, 13, 10],
        )

    def clean_lin(self, threshold: float) -> DataFrame:
        self.df_error.drop(["Line Number"], axis=1, inplace=True)
        self.df_error[self.df_error["Error"] != "UNFITTD"]
        self.lin_dataframe.apply(pd.to_numeric, errors="ignore")
        self.lin_dataframe.drop_duplicates(
            ["J", "Ka", "Kc", "J'", "Ka'", "Kc'"], keep="last", inplace=True
        )
        self.df_error["blend O-C"] = self.df_error["blend O-C"].fillna(
            self.df_error["O-C"]
        )
        self.df_error["blend weight"] = self.df_error["blend weight"].fillna(
            self.df_error["Error"]
        )
        self.df_error["Relative Error Blend"] = abs(self.df_error["blend O-C"]) / abs(
            self.df_error["Error"]
        )
        df_lin_res = self.lin_dataframe.join(
            self.df_error["Relative Error Blend"], how="right"
        )
        df_remove = df_lin_res[df_lin_res["Relative Error Blend"] > threshold]
        removal_list = df_remove.index.tolist()
        with open(self.file_path_lin) as f:
            lines = f.readlines()
        keys = np.arange(0, len(lines), 1)
        line_dict = dict(zip(keys, lines))
        for i in removal_list:
            line_dict.pop(i)
        with open(self.file_path_lin, "w") as output:
            for i in line_dict.values():
                output.write(str(i))

    def read_res(self) -> DataFrame:
        line_width = len(str(self.num_lines)) + 1
        self.res_dataframe = pd.read_fwf(
            self.file_path_res,
            header=None,
            skiprows=7,
            names=[
                "Line Number",
                "J",
                "Ka",
                "Kc",
                "v",
                "J'",
                "Ka'",
                "Kc'",
                "v'",
                "Obs-Freq",
                "O-C",
                "Error",
                "blend O-C",
                "blend weight",
            ],
            widths=[
                line_width,
                3,
                3,
                3,
                3,
                5,
                3,
                3,
                3,
                24,
                9,
                7,
                9,
                5,
            ],
        )
        last_row = self.res_dataframe[
            self.res_dataframe["Line Number"].str.contains("---", na=False)
        ].index[0]
        self.res_dataframe = self.res_dataframe.iloc[:last_row]
        self.res_dataframe[self.res_dataframe["Error"] != "UNFITTD"]
        self.delta_values()
        self.categorize_trans_type("Delta Ka", "Delta Kc")
        self.categorize_branch("Delta J")

    def categorize_trans_type(self, column_name1: str, column_name2: str):
        column1 = self.res_dataframe[column_name1].astype(int)
        column2 = self.res_dataframe[column_name2].astype(int)
        transition_types = []

        for Ka, Kc in zip(column1, column2):
            if Kc % 2 == 0:
                transition_types.append("c-type")
            elif Ka % 2 == 0:
                transition_types.append("a-type")
            else:
                transition_types.append("b-type")
        self.res_dataframe["Transition Type"] = transition_types

    def categorize_branch(self, column_name: str):
        column = self.res_dataframe[column_name]
        branches = []

        for value in column:
            if value == 1:
                branches.append("R-Branch")
            elif value == -1:
                branches.append("P-Branch")
            else:
                branches.append("Q-Branch")
        self.res_dataframe["Branch"] = branches

    def delta_values(self) -> DataFrame:
        self.res_dataframe = self.res_dataframe.apply(pd.to_numeric, errors="ignore")
        self.res_dataframe["Delta J"] = (
            self.res_dataframe["J"] - self.res_dataframe["J'"]
        )
        self.res_dataframe["Delta Ka"] = (
            self.res_dataframe["Ka"] - self.res_dataframe["Ka'"]
        )
        self.res_dataframe["Delta Kc"] = (
            self.res_dataframe["Kc"] - self.res_dataframe["Kc'"]
        )

    def categorize_error(self) -> DataFrame:
        self.res_dataframe["blend O-C"].fillna(
            value=self.res_dataframe["O-C"], inplace=True
        )
        self.res_dataframe["blend O-C"] = self.res_dataframe["blend O-C"].abs()
        self.res_dataframe["blend O-C/error"] = (
            self.res_dataframe["blend O-C"].abs() / self.res_dataframe["Error"]
        )
        self.res_dataframe["blend O-C/error"] = self.res_dataframe[
            "blend O-C/error"
        ].apply(lambda x: math.ceil(x))
        self.res_dataframe["blend O-C/error"].replace(0, 1, inplace=True)
        self.res_dataframe = self.res_dataframe.sort_values(
            by=["blend O-C/error", "Obs-Freq"], ascending=[True, False]
        )

    # Do I need this???
    def split_branches(self) -> Tuple[DataFrame, DataFrame, DataFrame]:
        df_r = self.res_dataframe.loc[self.res_dataframe["Branch"] == "R-Branch"]
        df_p = self.res_dataframe.loc[self.res_dataframe["Branch"] == "P-Branch"]
        df_q = self.res_dataframe.loc[self.res_dataframe["Branch"] == "Q-Branch"]
        return (df_r, df_p, df_q)

    def split_lines(self) -> Tuple[DataFrame, DataFrame, DataFrame]:
        df_ar = self.res_dataframe.loc[
            (self.res_dataframe["Transition Type"] == "a-type")
            & (self.res_dataframe["Branch"] == "R-Branch")
        ]
        df_br = self.res_dataframe.loc[
            (self.res_dataframe["Transition Type"] == "b-type")
            & (self.res_dataframe["Branch"] == "R-Branch")
        ]
        df_cr = self.res_dataframe.loc[
            (self.res_dataframe["Transition Type"] == "c-type")
            & (self.res_dataframe["Branch"] == "R-Branch")
        ]

        df_ap = self.res_dataframe.loc[
            (self.res_dataframe["Transition Type"] == "a-type")
            & (self.res_dataframe["Branch"] == "P-Branch")
        ]
        df_bp = self.res_dataframe.loc[
            (self.res_dataframe["Transition Type"] == "b-type")
            & (self.res_dataframe["Branch"] == "P-Branch")
        ]
        df_cp = self.res_dataframe.loc[
            (self.res_dataframe["Transition Type"] == "c-type")
            & (self.res_dataframe["Branch"] == "P-Branch")
        ]

        df_aq = self.res_dataframe.loc[
            (self.res_dataframe["Transition Type"] == "a-type")
            & (self.res_dataframe["Branch"] == "Q-Branch")
        ]
        df_bq = self.res_dataframe.loc[
            (self.res_dataframe["Transition Type"] == "b-type")
            & (self.res_dataframe["Branch"] == "Q-Branch")
        ]
        df_cq = self.res_dataframe.loc[
            (self.res_dataframe["Transition Type"] == "c-type")
            & (self.res_dataframe["Branch"] == "Q-Branch")
        ]

        return (df_ar, df_br, df_cr, df_ap, df_bp, df_cp, df_aq, df_bq, df_cq)

    def plot_data_dist(self, max_Ka: int, max_J: int):
        mpl.rc("text", usetex=True)
        plt.style.use(["science", "ieee"])
        self.categorize_error()
        self.res_dataframe = self.res_dataframe.sort_values(
            by=["Transition Type"], ascending=True
        )
        dataframes = self.split_branches()
        fig, ax = plt.subplots(1, 2, figsize=(8, 4))
        for i in range(0, 2):
            condition = dataframes[i]["blend O-C/error"] >= 4
            df_high_error = dataframes[i][condition]
            df_normal = dataframes[i][~condition]
            ax[0].scatter(
                df_normal["Ka'"],
                df_normal["J'"],
                s=4 * df_normal["blend O-C/error"],
                facecolors=df_normal["Transition Type"].map(
                    {
                        "a-type": "#0000FF65",
                        "b-type": "#00800045",
                        "c-type": "#FFFFFF00",
                    }
                ),
                edgecolors=df_normal["Transition Type"].map(
                    {"a-type": "blue", "b-type": "green", "c-type": "orange"}
                ),
                linewidth=0.8,
            )
            ax[0].scatter(
                df_high_error["Ka'"],
                df_high_error["J'"],
                s=4 * df_high_error["blend O-C/error"],
                facecolors="#FF000065",
                edgecolors="red",
                linewidth=0.8,
            )
        condition_q = dataframes[2]["blend O-C/error"] >= 4
        df_high_error_q = dataframes[2][condition_q]
        df_normal_q = dataframes[2][~condition_q]
        ax[1].scatter(
            df_high_error_q["Ka'"],
            df_high_error_q["J'"],
            s=4 * df_high_error_q["blend O-C/error"],
            facecolors="#FF000065",
            edgecolors="red",
            linewidth=0.8,
        )
        ax[1].scatter(
            df_normal_q["Ka'"],
            df_normal_q["J'"],
            s=4 * df_normal_q["blend O-C/error"],
            facecolors=df_normal_q["Transition Type"].map(
                {"a-type": "#0000FF65", "b-type": "#00800045", "c-type": "#FFFFFF00"}
            ),
            edgecolors=df_normal_q["Transition Type"].map(
                {"a-type": "blue", "b-type": "green", "c-type": "orange"}
            ),
            linewidth=0.8,
        )

        ax[0].axline((0, 0), slope=1, color="lightgrey", linewidth=1.5)
        ax[0].set_xlabel("$K_a''$", fontsize=14)
        ax[0].set_ylabel("$J''$", fontsize=14)
        ax[0].set_ylim(0, max_J)
        ax[0].set_xlim(-0.5, max_Ka)
        ax[0].set_xticks((0, 10, 20, 30, 40, 50, 60))
        ax[0].set_yticks((0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120))
        ax[0].text(38, 110, "R-Branch", fontsize=18)

        ax[1].axline((0, 0), slope=1, color="lightgrey", linewidth=1.5)
        ax[1].set_xlabel("$K_a''$", fontsize=14)
        ax[1].set_ylabel("$J''$", fontsize=14)
        ax[1].set_ylim(0, max_J)
        ax[1].set_xlim(-0.5, max_Ka)
        ax[1].set_xticks((0, 10, 20, 30, 40, 50, 60))
        ax[1].set_yticks((0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120))
        ax[1].text(38, 110, "Q-Branch", fontsize=18)

        c = [
            mpatches.Circle(
                (0.5, 0.5),
                radius=0.25,
                facecolor="#0000FF65",
                edgecolor="blue",
                label="a-type",
            ),
            mpatches.Circle(
                (0.5, 0.5),
                radius=0.25,
                facecolor="#00800045",
                edgecolor="green",
                label="b-type",
            ),
            mpatches.Circle(
                (0.5, 0.5),
                radius=0.25,
                facecolor="#FFA50065",
                edgecolor="orange",
                label="c-type",
            ),
            mpatches.Circle(
                (0.5, 0.5),
                radius=0.25,
                facecolor="#FF000065",
                edgecolor="red",
                label=r"$>$ $3\sigma$",
            ),
        ]
        texts = ["a-type", "b-type", "c-type", r"$>$ $3\sigma$"]
        legend = ax[0].legend(
            c,
            texts,
            loc="lower right",
            fontsize=8,
            frameon=True,
            edgecolor="black",
            framealpha=1,
            borderaxespad=0,
            fancybox=False,
            handler_map={mpatches.Circle: HandlerEllipse()},
        )
        frame = legend.get_frame()
        frame.set_linewidth(0.5)

        legend1 = ax[1].legend(
            c,
            texts,
            loc="lower right",
            fontsize=8,
            frameon=True,
            edgecolor="black",
            framealpha=1,
            borderaxespad=0,
            fancybox=False,
            handler_map={mpatches.Circle: HandlerEllipse()},
        )
        frame1 = legend1.get_frame()
        frame1.set_linewidth(0.5)
        file_path = self.file_name + ".jpg"
        fig.savefig(file_path, dpi=800)

    def plot_data_dist_notebook(self, max_Ka: int, max_J: int):
        mpl.rcParams.update(mpl.rcParamsDefault)
        self.categorize_error()
        dataframes = self.split_lines()
        fig, ax = plt.subplots(1, 2, figsize=(8, 4))
        for i in range(0, 6):
            condition = dataframes[i]["blend O-C/error"] >= 4
            df_high_error = dataframes[i][condition]
            df_normal = dataframes[i][~condition]
            ax[0].scatter(
                df_normal["Ka'"],
                df_normal["J'"],
                s=4 * df_normal["blend O-C/error"],
                facecolors="white",
                edgecolors=df_normal["Transition Type"].map(
                    {"a-type": "#0000FF65", "b-type": "#00800045", "c-type": "orange"}
                ),
                linewidth=0.8,
                alpha=0.5,
            )
            ax[0].scatter(
                df_high_error["Ka'"],
                df_high_error["J'"],
                s=4 * df_high_error["blend O-C/error"],
                facecolors="white",
                edgecolors="red",
                linewidth=0.8,
                alpha=0.5,
            )
        for i in range(8, 5):
            condition_q = dataframes[i]["blend O-C/error"] >= 4
            df_high_error_q = dataframes[i][condition_q]
            df_normal_q = dataframes[i][~condition_q]
            ax[1].scatter(
                df_high_error_q["Ka'"],
                df_high_error_q["J'"],
                s=4 * df_high_error_q["blend O-C/error"],
                facecolors="white",
                edgecolors="red",
                linewidth=0.8,
                alpha=0.5,
            )
            ax[1].scatter(
                df_normal_q["Ka'"],
                df_normal_q["J'"],
                s=4 * df_normal_q["blend O-C/error"],
                facecolors=df_normal_q["Transition Type"].map(
                    {"a-type": "#0000FF65", "b-type": "#00800045", "c-type": "orange"}
                ),
                edgecolors=df_normal_q["Transition Type"].map(
                    {"a-type": "#0000FF65", "b-type": "#00800045", "c-type": "orange"}
                ),
                alpha=0.5,
                linewidth=0.8,
            )

        ax[0].axline((0, 0), slope=1, color="lightgrey", linewidth=1.5)
        ax[0].set_xlabel("$K_a''$", fontsize=14)
        ax[0].set_ylabel("$J''$", fontsize=14)
        ax[0].set_ylim(0, max_J)
        ax[0].set_xlim(-0.5, max_Ka)
        ax[0].set_xticks((0, 10, 20, 30, 40, 50, 60))
        ax[0].set_yticks((0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120))
        ax[0].text(38, 110, "R-Branch", fontsize=18)

        ax[1].axline((0, 0), slope=1, color="lightgrey", linewidth=1.5)
        ax[1].set_xlabel("$K_a''$", fontsize=14)
        ax[1].set_ylabel("$J''$", fontsize=14)
        ax[1].set_ylim(0, max_J)
        ax[1].set_xlim(-0.5, max_Ka)
        ax[1].set_xticks((0, 10, 20, 30, 40, 50, 60))
        ax[1].set_yticks((0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120))
        ax[1].text(38, 110, "Q-Branch", fontsize=18)

        c = [
            mpatches.Circle(
                (0.5, 0.5),
                radius=0.25,
                facecolor="#0000FF65",
                edgecolor="none",
                label="a-type",
            ),
            mpatches.Circle(
                (0.5, 0.5),
                radius=0.25,
                facecolor="#00800045",
                edgecolor="none",
                label="b-type",
            ),
            mpatches.Circle(
                (0.5, 0.5),
                radius=0.25,
                facecolor="orange",
                edgecolor="none",
                label="c-type",
            ),
            mpatches.Circle(
                (0.5, 0.5),
                radius=0.25,
                facecolor="red",
                edgecolor="none",
                label=r"$>$ $3\sigma$",
            ),
        ]
        texts = ["a-type", "b-type", "c-type", r"$>$ $3\sigma$"]
        legend = ax[0].legend(
            c,
            texts,
            loc="lower right",
            fontsize=8,
            frameon=True,
            edgecolor="black",
            framealpha=1,
            borderaxespad=0,
            fancybox=False,
            handler_map={mpatches.Circle: HandlerEllipse()},
        )
        frame = legend.get_frame()
        frame.set_linewidth(0.5)

        legend1 = ax[1].legend(
            c,
            texts,
            loc="lower right",
            fontsize=8,
            frameon=True,
            edgecolor="black",
            framealpha=1,
            borderaxespad=0,
            fancybox=False,
            handler_map={mpatches.Circle: HandlerEllipse()},
        )
        frame1 = legend1.get_frame()
        frame1.set_linewidth(0.5)
        plt.show()


class HandlerEllipse(HandlerPatch):
    def create_artists(
        self,
        legend,
        orig_handle,
        xdescent,
        ydescent,
        width,
        height,
        fontsize,
        trans,
    ):
        center = 0.5 * width - 0.5 * xdescent, 0.5 * height - 0.5 * ydescent
        p = mpatches.Ellipse(
            xy=center, width=height + xdescent, height=height + ydescent
        )
        self.update_prop(p, orig_handle, legend)
        p.set_transform(trans)
        return [p]
