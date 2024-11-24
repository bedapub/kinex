import numpy as np
import pandas as pd
import plotly.express as px
from statsmodels.stats.multitest import multipletests

from kinex.resources import (
    get_ser_thr_family,
    get_ser_thr_family_colors,
    get_tyr_family,
    get_tyr_family_colors,
)
from kinex.sequence import SequenceType
from kinex.table2x2 import Table2x2


class Enrichment:
    def __init__(self, sequence_type: SequenceType, all_kinases: set):
        """
        Initialize the Enrichment class.

        :param sequence_type: Type of sequence (SER_THR or TYR).
        :param all_kinases: Set of all kinases.
        """
        self.sequence_type = sequence_type
        self.total_upregulated = 0
        self.total_downregulated = 0
        self.total_unregulated = 0
        self.regulation_list = []
        self.top15_kinases_list = []
        self.enrichment_table = pd.DataFrame(
            columns=["kinase", "upregulated", "downregulated", "unregulated"]
        )
        self.all_kinases = all_kinases

    def adjust_background_sites(self):
        """
        Adjust the background sites if the total unregulated is zero.
        """
        if self.total_unregulated == 0:
            self.total_unregulated = (
                np.min([self.total_upregulated, self.total_downregulated]) / 2
            )

    def fisher_statistics(self):
        """
        Perform Fisher's exact test to analyze kinase enrichment in regulation states.

        This method calculates the statistical enrichment of kinases associated with
        upregulated or downregulated sites using Fisher's exact test. The test determines
        whether there is a significant association between two categorical variables:
        - Presence or absence of a specific kinase.
        - Regulation state (upregulated or downregulated).

        Outputs:
        - p-value: Probability of observing the given enrichment by chance.
        - Enrichment values: Ratio of regulated sites (up or down) to the total number of sites.

        Updates the enrichment table with the following columns:
        - `upregulated_enrichment_value`: Enrichment value for upregulated sites.
        - `downregulated_enrichment_value`: Enrichment value for downregulated sites.
        - `dominant_direction`: Regulation state with the highest enrichment value.
        """

        new_columns = [
            "upregulated_enrichment_value",
            "upregulated_enrichment_value_log2",
            "upregulated_p_value",
            "upregulated_p_value_log10_abs",
            "upregulated_adjusted_p_value",
            "upregulated_adjusted_p_value_log10_abs",
            "downregulated_enrichment_value",
            "downregulated_enrichment_value_log2",
            "downregulated_p_value",
            "downregulated_p_value_log10_abs",
            "downregulated_adjusted_p_value",
            "downregulated_adjusted_p_value_log10_abs",
            "dominant_direction",
            "dominant_enrichment_value_log2",
            "dominant_p_value_log10_abs",
            "dominant_adjusted_p_value_log10_abs",
        ]

        # Create new columns in the enrichment table
        self.enrichment_table = self.enrichment_table.assign(
            **{col: None for col in new_columns}
        )

        for i in range(len(self.enrichment_table)):
            self._calculate_enrichment_for_row(i)

        self._calculate_adjusted_p_values()

        self.enrichment_table = self.enrichment_table.set_index("kinase")

        self._reindex_missing_kinases()

    def _calculate_enrichment_for_row(self, index: int):
        """
        Calculate enrichment values for a specific row.

        :param index: Index of the row in the enrichment table.
        """
        # Add background
        # TODO Look into adding unregulated_hit backgroud
        if self.enrichment_table.loc[index, "unregulated"] == 0:
            self.enrichment_table.loc[index, "unregulated"] += 1

        # Save hits to the vars
        upregulated_hit = int(self.enrichment_table.iloc[index]["upregulated"])
        downregulated_hit = int(self.enrichment_table.iloc[index]["downregulated"])
        unregulated_hit = int(self.enrichment_table.iloc[index]["unregulated"])

        # Upregulated enrichment and p values
        if self.total_upregulated == 0 or upregulated_hit == 0:
            upregulated_enrichment_value = 0
            upregulated_p_value = 1
            self.enrichment_table.loc[index, "upregulated_enrichment_value_log2"] = 0
            self.enrichment_table.loc[index, "upregulated_p_value_log10_abs"] = 0
        else:
            upregulated_contingency_table = Table2x2(
                [
                    [upregulated_hit, self.total_upregulated - upregulated_hit],
                    [unregulated_hit, self.total_unregulated - unregulated_hit],
                ],
                shift_zeros=True,
            )
            upregulated_enrichment_value = upregulated_contingency_table.odds_ratio()
            upregulated_p_value = upregulated_contingency_table.p_val(mode="greater")
            self.enrichment_table.loc[index, "upregulated_enrichment_value_log2"] = (
                np.abs(np.log2(upregulated_enrichment_value))
            )
            self.enrichment_table.loc[index, "upregulated_p_value_log10_abs"] = (
                np.absolute(np.log10(upregulated_p_value))
            )

        # Downregulated enrichment and p values
        if self.total_downregulated == 0 or downregulated_hit == 0:
            downregulated_enrichment_value = 0
            downregulated_p_value = 1
            self.enrichment_table.loc[index, "downregulated_enrichment_value_log2"] = 0
            self.enrichment_table.loc[index, "downregulated_p_value_log10_abs"] = 0
        else:
            downregulated_contingency_table = Table2x2(
                [
                    [downregulated_hit, self.total_downregulated - downregulated_hit],
                    [unregulated_hit, self.total_unregulated - unregulated_hit],
                ],
                shift_zeros=True,
            )
            downregulated_enrichment_value = (
                downregulated_contingency_table.odds_ratio()
            )
            downregulated_p_value = downregulated_contingency_table.p_val(
                mode="greater"
            )

            self.enrichment_table.loc[index, "downregulated_enrichment_value_log2"] = (
                -np.abs(np.log2(downregulated_enrichment_value))
            )
            self.enrichment_table.loc[index, "downregulated_p_value_log10_abs"] = (
                np.absolute(np.log10(downregulated_p_value))
            )

        # Set the enrichment and p values
        self.enrichment_table.loc[index, "upregulated_enrichment_value"] = (
            upregulated_enrichment_value
        )
        self.enrichment_table.loc[index, "upregulated_p_value"] = upregulated_p_value

        self.enrichment_table.loc[index, "downregulated_enrichment_value"] = (
            downregulated_enrichment_value
        )
        self.enrichment_table.loc[index, "downregulated_p_value"] = (
            downregulated_p_value
        )

        self._determine_dominant_direction(index)

    def _determine_dominant_direction(self, index: int):
        """
        Determine the dominant direction (upregulated or downregulated) for a specific row.

        :param index: Index of the row in the enrichment table.
        """
        if (
            self.enrichment_table.loc[index, "upregulated_enrichment_value"]
            > self.enrichment_table.loc[index, "downregulated_enrichment_value"]
        ):
            self.enrichment_table.loc[index, "dominant_direction"] = "upregulated set"
            self.enrichment_table.loc[index, "dominant_enrichment_value_log2"] = (
                self.enrichment_table.loc[index, "upregulated_enrichment_value_log2"]
            )
            self.enrichment_table.loc[index, "dominant_p_value_log10_abs"] = (
                self.enrichment_table.loc[index, "upregulated_p_value_log10_abs"]
            )
        else:
            self.enrichment_table.loc[index, "dominant_direction"] = "downregulated set"
            self.enrichment_table.loc[index, "dominant_enrichment_value_log2"] = (
                self.enrichment_table.loc[index, "downregulated_enrichment_value_log2"]
            )
            self.enrichment_table.loc[index, "dominant_p_value_log10_abs"] = (
                self.enrichment_table.loc[index, "downregulated_p_value_log10_abs"]
            )

    def _calculate_adjusted_p_values(self):
        """
        Calculate adjusted p-values using the Benjamini-Hochberg method.
        """
        upregulated_adjusted_p_value = multipletests(
            self.enrichment_table["upregulated_p_value"], method="fdr_bh"
        )
        self.enrichment_table["upregulated_adjusted_p_value"] = (
            upregulated_adjusted_p_value[1]
        )
        downregulated_adjusted_p_value = multipletests(
            self.enrichment_table["downregulated_p_value"], method="fdr_bh"
        )
        self.enrichment_table["downregulated_adjusted_p_value"] = (
            downregulated_adjusted_p_value[1]
        )

        for i in range(len(self.enrichment_table)):
            # Adjusted p values log10 abs and dominant adjusted p values log10 abs
            self.enrichment_table.loc[i, "upregulated_adjusted_p_value_log10_abs"] = (
                np.absolute(
                    np.log10(
                        self.enrichment_table.loc[i, "upregulated_adjusted_p_value"]
                    )
                )
            )
            self.enrichment_table.loc[i, "downregulated_adjusted_p_value_log10_abs"] = (
                np.absolute(
                    np.log10(
                        self.enrichment_table.loc[i, "downregulated_adjusted_p_value"]
                    )
                )
            )

            if (
                self.enrichment_table.loc[i, "dominant_direction"]
                == "downregulated set"
            ):
                self.enrichment_table.loc[i, "dominant_adjusted_p_value_log10_abs"] = (
                    self.enrichment_table.loc[
                        i, "downregulated_adjusted_p_value_log10_abs"
                    ]
                )
            elif (
                self.enrichment_table.loc[i, "dominant_direction"] == "upregulated set"
            ):
                self.enrichment_table.loc[i, "dominant_adjusted_p_value_log10_abs"] = (
                    self.enrichment_table.loc[
                        i, "upregulated_adjusted_p_value_log10_abs"
                    ]
                )

    def _reindex_missing_kinases(self):
        """
        Reindex the enrichment table to include missing kinases.
        """
        missing_kinases = list(self.all_kinases - set(self.enrichment_table.index))
        self.enrichment_table = self.enrichment_table.reindex(
            self.enrichment_table.index.union(missing_kinases), fill_value=0
        )

    def plot(self, use_adjusted_pval: bool = False):
        match self.sequence_type:
            case SequenceType.SER_THR:
                kinase_family = get_ser_thr_family()
                family_colors = get_ser_thr_family_colors()
            case SequenceType.TYR:
                kinase_family = get_tyr_family()
                family_colors = get_tyr_family_colors()
            case _:
                raise ValueError("Invalid sequence type")

        family = []
        for kinase in self.enrichment_table.index:
            family.append(kinase_family[kinase])

        # Plot for 2 different cases: with and without adjusted p-val
        if use_adjusted_pval:
            y = "dominant_adjusted_p_value_log10_abs"
            y_axis_title = "-Log\u2081\u2080(adjusted p-value)"
        else:
            y = "dominant_p_value_log10_abs"
            y_axis_title = "-Log\u2081\u2080(p-value)"

        fig = px.scatter(
            self.enrichment_table,
            x="dominant_enrichment_value_log2",
            y=y,
            hover_name=self.enrichment_table.index,
            color=family,
            text=self.enrichment_table.index,
            color_discrete_map=family_colors,  # Use the defined color mapping
            # category_orders={"family": category_orders},
            template="none",
        )

        # Add horizontal line at y = 1.3 which is absolute val of log10(0.05)
        fig.add_hline(
            y=1.3,
            line_dash="dash",
            line_width=1,
            line_color="black",
            annotation_text="p ≤ 0.05",
            annotation_position="top right",
            annotation_font_size=11,
            annotation_font_color="black",
            opacity=0.7,
        )
        # Add vertical line at position x = 0
        fig.add_vline(
            x=0,
            annotation_text="Inhibited | Activated",
            annotation_position="top",
            annotation_font_size=11,
            line_width=1,
            line_color="black",
            opacity=0.5,
        )
        # Format the text that appears on the points. (the kinases names)
        fig.update_traces(
            textfont_size=7,
            textposition="middle right",
            marker=dict(size=6),
        )
        # Format the legend and axis
        fig.update_layout(
            legend=dict(font=dict(size=10)),
            legend_title=dict(font=dict(size=14, color="black")),
            legend_title_text="Family",
            xaxis_title_font=dict(size=12),  # Set x-axis label font size
            yaxis_title_font=dict(size=12),  # Set y-axis label font size
            xaxis_title="Log\u2082(EOR)",
            yaxis_title=y_axis_title,
            title={
                "text": "Kinase inference",
                "x": 0.5,
                "xanchor": "center",
                "font_size": 15,
            },
            xaxis=dict(ticks="outside", mirror=True, showline=True),
            yaxis=dict(ticks="outside", mirror=True, showline=True),
            width=600,
            height=600,
        )
        return fig
