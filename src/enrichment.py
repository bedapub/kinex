import numpy as np
import pandas as pd
import plotly.express as px
from statsmodels.stats.multitest import multipletests

from kinex.table2x2 import Table2x2
from kinex.data import get_groups


class Enrichment:
    """
    Enrichment results including all the necessary information for the input sites

    Attributes
    ----------
    enrichment_table: pandas.DataFrame
        Table containing enrichment analysis results. 
    input_sites: pandas.DataFrame
        A DataFrame containing the phosphosite sequences in the first column, logarithmised Fold Change in the second column, 
        regulation in the third column, and top15 kinases most likely to target each sequence in the fourth column.
    failed_sites: list
        A list of the sites labeled invalid by the check_sequence function.
    total_upregulated: int
        Total number of upregulated phosphosites
    total_downregulated: int
        Total number of downregulated phosphosites
    total_unregulated: int
        Total number of unregulated phosphosites

    Methods
    -------
    plot(self, use_adjusted_pval: bool = False) -> None
        Vulcano plot of the side displaying significant enrichment for each kinase vs the corresponding p-value.

    """

    def __init__(self, enrichment_table: pd.DataFrame, input_sites: pd.DataFrame, failed_sites: list, total_upregulated: int, total_downregulated: int, total_unregulated: int, all_kinases: set) -> None:

        self.enrichment_table = enrichment_table
        self.total_upregulated = total_upregulated
        self.total_downregulated = total_downregulated
        self.total_unregulated = total_unregulated
        self.input_sites = input_sites
        self.failed_sites = failed_sites

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

        # Create new columns
        self.enrichment_table = self.enrichment_table.assign(
            **{col: None for col in new_columns})

        for i in range(len(self.enrichment_table)):
            # Add background
            # TODO Look into adding unregulated_hit backgroud
            if self.enrichment_table.loc[i, "unregulated"] == 0:
                self.enrichment_table.loc[i, "unregulated"] += 1

            # Save hits to the vars
            upregulated_hit = int(self.enrichment_table.iloc[i]["upregulated"])
            downregulated_hit = int(
                self.enrichment_table.iloc[i]["downregulated"])
            unregulated_hit = int(self.enrichment_table.iloc[i]["unregulated"])

            # Upregulated enrichment and p values
            if self.total_upregulated == 0 or upregulated_hit == 0:
                upregulated_enrichment_value = 0
                upregulated_p_value = 1
                self.enrichment_table.loc[i,
                                          "upregulated_enrichment_value_log2"] = 0
                self.enrichment_table.loc[i,
                                          "upregulated_p_value_log10_abs"] = 0
            else:
                upregulated_contingency_table = Table2x2([[upregulated_hit, self.total_upregulated - upregulated_hit],
                                                          [unregulated_hit, self.total_unregulated - unregulated_hit]],
                                                         shift_zeros=True
                                                         )
                upregulated_enrichment_value = upregulated_contingency_table.odds_ratio()
                upregulated_p_value = upregulated_contingency_table.p_val(
                    mode="greater")
                self.enrichment_table.loc[i, "upregulated_enrichment_value_log2"] = np.abs(np.log2(
                    upregulated_enrichment_value))
                self.enrichment_table.loc[i, "upregulated_p_value_log10_abs"] = np.absolute(
                    np.log10(upregulated_p_value))

            # Downregulated enrichment and p values
            if self.total_downregulated == 0 or downregulated_hit == 0:
                downregulated_enrichment_value = 0
                downregulated_p_value = 1
                self.enrichment_table.loc[i,
                                          "downregulated_enrichment_value_log2"] = 0
                self.enrichment_table.loc[i,
                                          "downregulated_p_value_log10_abs"] = 0
            else:
                downregulated_contingency_table = Table2x2([[downregulated_hit, total_downregulated - downregulated_hit],
                                                            [unregulated_hit, total_unregulated - unregulated_hit]],
                                                           shift_zeros=True
                                                           )
                downregulated_enrichment_value = downregulated_contingency_table.odds_ratio()
                downregulated_p_value = downregulated_contingency_table.p_val(
                    mode="greater")

                self.enrichment_table.loc[i, "downregulated_enrichment_value_log2"] = -np.abs(np.log2(
                    downregulated_enrichment_value))
                self.enrichment_table.loc[i, "downregulated_p_value_log10_abs"] = np.absolute(
                    np.log10(downregulated_p_value))

            # Set the enrichment and p values
            self.enrichment_table.loc[i,
                                      "upregulated_enrichment_value"] = upregulated_enrichment_value
            self.enrichment_table.loc[i,
                                      "upregulated_p_value"] = upregulated_p_value

            self.enrichment_table.loc[i,
                                      "downregulated_enrichment_value"] = downregulated_enrichment_value
            self.enrichment_table.loc[i,
                                      "downregulated_p_value"] = downregulated_p_value

            # Determine the dominant direction, either upregulated or downregulated.
            if upregulated_enrichment_value > downregulated_enrichment_value:
                self.enrichment_table.loc[i,
                                          "dominant_direction"] = "upregulated set"
                self.enrichment_table.loc[i, "dominant_enrichment_value_log2"] = self.enrichment_table.loc[i,
                                                                                                           "upregulated_enrichment_value_log2"]
                self.enrichment_table.loc[i, "dominant_p_value_log10_abs"] = self.enrichment_table.loc[i,
                                                                                                       "upregulated_p_value_log10_abs"]
            else:
                self.enrichment_table.loc[i,
                                          "dominant_direction"] = "downregulated set"
                self.enrichment_table.loc[i, "dominant_enrichment_value_log2"] = self.enrichment_table.loc[i,
                                                                                                           "downregulated_enrichment_value_log2"]
                self.enrichment_table.loc[i, "dominant_p_value_log10_abs"] = self.enrichment_table.loc[i,
                                                                                                       "downregulated_p_value_log10_abs"]

        # Calculate adjusted p values
        upregulated_adjusted_p_value = multipletests(
            self.enrichment_table["upregulated_p_value"], method="fdr_bh")
        self.enrichment_table["upregulated_adjusted_p_value"] = upregulated_adjusted_p_value[1]
        downregulated_adjusted_p_value = multipletests(
            self.enrichment_table["downregulated_p_value"], method="fdr_bh")
        self.enrichment_table["downregulated_adjusted_p_value"] = downregulated_adjusted_p_value[1]

        for i in range(len(self.enrichment_table)):

            # adjusted p values log10 abs and dominant adjusted p values log10 abs
            self.enrichment_table.loc[i, "upregulated_adjusted_p_value_log10_abs"] = np.absolute(
                np.log10(self.enrichment_table.loc[i, "upregulated_adjusted_p_value"]))
            self.enrichment_table.loc[i, "downregulated_adjusted_p_value_log10_abs"] = np.absolute(
                np.log10(self.enrichment_table.loc[i, "downregulated_adjusted_p_value"]))

            if self.enrichment_table.loc[i, "dominant_direction"] == "downregulated set":
                self.enrichment_table.loc[i, "dominant_adjusted_p_value_log10_abs"] = self.enrichment_table.loc[i,
                                                                                                                "downregulated_adjusted_p_value_log10_abs"]
            elif self.enrichment_table.loc[i, "dominant_direction"] == "upregulated set":
                self.enrichment_table.loc[i, "dominant_adjusted_p_value_log10_abs"] = self.enrichment_table.loc[i,
                                                                                                                "upregulated_adjusted_p_value_log10_abs"]
        self.enrichment_table = self.enrichment_table.set_index("kinase")
        missing_kinases = list(all_kinases - set(self.enrichment_table.index))
        self.enrichment_table = self.enrichment_table.reindex(
            self.enrichment_table.index.union(missing_kinases), fill_value=0)

    def __repr__(self) -> str:
        return f"Total number of upregulated phospho-sequences is: {self.total_upregulated}\nTotal number of downregulated phospho-sequences is: {self.total_downregulated}\nTotal number of unregulated phospho-sequences is: {self.total_unregulated}"

    def plot(self, use_adjusted_pval: bool = False):

        kinase_family = get_groups()

        family = []
        plotting_table = self.enrichment_table[self.enrichment_table.dominant_p_value_log10_abs > 0.3]

        for kinase in plotting_table.index:
            family.append(kinase_family[kinase])

        # Colors of the kinase groups which correspond to the ones from Johnson et al. https://doi.org/10.1038/s41586-022-05575-3
        family_colors = {
            "TKL": "#0400A9",
            "STE": "#754391",
            "CK1": "#499FD2",
            "AGC": "#992B2C",
            "CAMK": "#DE643A",
            "Other": "#957C64",
            "CMGC": "#2C643C",
            "Alpha": "#C88AB5",
            "PDHK": "#316D79",
            "PIKK": "#020101",
            "FAM20C": "#C22F7F"
        }
        # Plot for 2 different cases: with and without adjusted p-val
        if use_adjusted_pval:
            fig = px.scatter(
                plotting_table,
                x="dominant_enrichment_value_log2",
                y="dominant_adjusted_p_value_log10_abs",
                hover_name=plotting_table.index,
                color=family,
                text=plotting_table.index,
                color_discrete_map=family_colors,  # Use the defined color mapping
                category_orders={"family": [
                    'TKL', 'STE', 'CK1', 'AGC', 'CAMK', 'Other', 'CMGC', 'Alpha', 'PDHK', 'PIKK', 'FAM20C']},
                template="none"
            )
            y_axis_title = "-Log\u2081\u2080(adjusted p-value)"

        else:
            fig = px.scatter(
                plotting_table,
                x="dominant_enrichment_value_log2",
                y="dominant_p_value_log10_abs",
                hover_name=plotting_table.index,
                color=family,
                text=plotting_table.index,
                color_discrete_map=family_colors,  # Use the defined color mapping
                category_orders={"family": [
                    'TKL', 'STE', 'CK1', 'AGC', 'CAMK', 'Other', 'CMGC', 'Alpha', 'PDHK', 'PIKK', 'FAM20C']},
                template="none"
            )
            y_axis_title = "-Log\u2081\u2080(p-value)"

        # add horizontal line at y = 1.3 which is absolute val of log10(0.05)
        fig.add_hline(
            y=1.3,
            line_dash="dash",
            line_width=1,
            line_color="black",
            annotation_text="p ≤ 0.05",
            annotation_position="top right",
            annotation_font_size=11,
            annotation_font_color="black",
            opacity=0.7
        )
        # add vertical line at position x = 0
        fig.add_vline(
            x=0,
            annotation_text="Inhibited | Activated",
            annotation_position="top",
            annotation_font_size=11,
            line_width=1,
            line_color="black",
            opacity=0.5
        )
        # format the text that appears on the points. (the kinases names)
        fig.update_traces(
            textfont_size=7,
            textposition="middle right",
            marker=dict(size=6),
        )
        # format the legend and axis
        fig.update_layout(
            legend=dict(font=dict(size=10)),
            legend_title=dict(font=dict(size=14, color="black")),
            legend_title_text="Family",
            xaxis_title_font=dict(size=12),  # Set x-axis label font size
            yaxis_title_font=dict(size=12),  # Set y-axis label font size
            xaxis_title='Log\u2082(EOR)',
            yaxis_title=y_axis_title,
            title={
                'text': "Kinase inference",
                'x': 0.5,
                'xanchor': 'center',
                "font_size": 15},
            xaxis=dict(ticks="outside", mirror=True, showline=True),
            yaxis=dict(ticks="outside", mirror=True, showline=True),
            width=600,
            height=600
        )
#         fig.show()
        return fig
