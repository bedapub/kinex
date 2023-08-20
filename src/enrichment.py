import numpy as np
import pandas as pd
import plotly.express as px
from statsmodels.stats.multitest import multipletests

from table2x2 import Table2x2
from data import get_groups


class Enrichment:
    """
    Enrichment results including all the necessary information for the input sites

    Attributes
    ----------
    enrichment_table: pandas.DataFrame
        Table containing enrichment analysis results. 
    input_sites: pandas.DataFrame
        A DataFrame containing the phosphosite sequences in the first column, logarithmised Fold Change in the second column, and regulation in the third column.
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

    def __init__(self, enrichment_table: pd.DataFrame, input_sites: pd.DataFrame, failed_sites: list, total_upregulated: int, total_downregulated: int, total_unregulated: int) -> None:

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
            # Save hits to the vars
            upregulated_hit = int(self.enrichment_table.iloc[i]["upregulated"])
            downregulated_hit = int(
                self.enrichment_table.iloc[i]["downregulated"])
            unregulated_hit = int(self.enrichment_table.iloc[i]["unregulated"])

            # Contingency tables for upregulated / downregulated counts
            upregulated_contingency_table = Table2x2([[upregulated_hit, total_upregulated - upregulated_hit],
                                                      [unregulated_hit, total_unregulated - unregulated_hit]],
                                                     shift_zeros=True
                                                     )
            downregulated_contingency_table = Table2x2([[downregulated_hit, total_downregulated - downregulated_hit],
                                                        [unregulated_hit, total_unregulated - unregulated_hit]],
                                                       shift_zeros=True
                                                       )

            # Calculate enrichment values
            upregulated_enrichment_value = upregulated_contingency_table.odds_ratio()
            downregulated_enrichment_value = downregulated_contingency_table.odds_ratio()

            self.enrichment_table.loc[i,
                                      "upregulated_enrichment_value"] = upregulated_enrichment_value
            self.enrichment_table.loc[i,
                                      "downregulated_enrichment_value"] = downregulated_enrichment_value

            # log2(enrichment_value)
            self.enrichment_table.loc[i, "upregulated_enrichment_value_log2"] = np.log2(
                upregulated_enrichment_value)
            self.enrichment_table.loc[i, "downregulated_enrichment_value_log2"] = np.log2(
                downregulated_enrichment_value)

            # p values
            upregulated_p_value = upregulated_contingency_table.p_val(
                mode="greater")
            downregulated_p_value = downregulated_contingency_table.p_val(
                mode="greater")

            self.enrichment_table.loc[i,
                                      "upregulated_p_value"] = upregulated_p_value
            self.enrichment_table.loc[i,
                                      "downregulated_p_value"] = downregulated_p_value

            self.enrichment_table.loc[i, "upregulated_p_value_log10_abs"] = np.absolute(
                np.log10(upregulated_p_value))
            self.enrichment_table.loc[i, "downregulated_p_value_log10_abs"] = np.absolute(
                np.log10(downregulated_p_value))
            
            # Determine the dominant direction, either upregulated or downregulated. 

            if upregulated_enrichment_value > downregulated_enrichment_value:
                self.enrichment_table.loc[i,
                                          "dominant_direction"] = "upregulated set"
                self.enrichment_table.loc[i, "dominant_enrichment_value_log2"] = np.absolute(
                    np.log2(upregulated_enrichment_value))
                self.enrichment_table.loc[i, "dominant_p_value_log10_abs"] = np.absolute(
                    np.log10(upregulated_p_value))
            else:
                self.enrichment_table.loc[i,
                                          "dominant_direction"] = "downregulated set"
                self.enrichment_table.loc[i, "dominant_enrichment_value_log2"] = -np.absolute(
                    np.log2(downregulated_enrichment_value))
                self.enrichment_table.loc[i, "dominant_p_value_log10_abs"] = np.absolute(
                    np.log10(downregulated_p_value))

        # Calculate adjusted p values
        upregulated_adjusted_p_value = multipletests(
            self.enrichment_table["upregulated_p_value"], method="fdr_bh")
        self.enrichment_table["upregulated_adjusted_p_value"] = upregulated_adjusted_p_value[1]
        downregulated_adjusted_p_value = multipletests(
            self.enrichment_table["downregulated_p_value"], method="fdr_bh")
        self.enrichment_table["downregulated_adjusted_p_value"] = downregulated_adjusted_p_value[1]

        
        for i in range(len(self.enrichment_table)):
            # Add background
            if self.enrichment_table.loc[i, "unregulated"] == 0:
                self.enrichment_table.loc[i, "unregulated"] += 1
            
            # adjusted p values log10 abs and dominant adjusted p values log10 abs 
            self.enrichment_table.loc[i, "upregulated_adjusted_p_value_log10_abs"] = np.absolute(np.log10(self.enrichment_table.loc[i, "upregulated_adjusted_p_value"]))
            self.enrichment_table.loc[i, "downregulated_adjusted_p_value_log10_abs"] = np.absolute(np.log10(self.enrichment_table.loc[i, "downregulated_adjusted_p_value"]))


            if self.enrichment_table.loc[i, "dominant_direction"] == "downregulated set":
                self.enrichment_table.loc[i, "dominant_adjusted_p_value_log10_abs"] = self.enrichment_table.loc[i, "downregulated_adjusted_p_value_log10_abs"]
            elif self.enrichment_table.loc[i, "dominant_direction"] == "upregulated set":
                self.enrichment_table.loc[i, "dominant_adjusted_p_value_log10_abs"] = self.enrichment_table.loc[i, "upregulated_adjusted_p_value_log10_abs"]

        self.enrichment_table = self.enrichment_table.set_index("kinase")

    def __repr__(self) -> str:
        return f"Total number of upregulated sites is: {self.total_upregulated}\nTotal number of downregulated sites is: {self.total_downregulated}\nTotal number of unregulated sites is: {self.total_unregulated}"

    def plot(self, use_adjusted_pval: bool = False):

        kinase_group = get_groups()

        group = []
        for kinase in self.enrichment_table.index:
            group.append(kinase_group[kinase])
        self.enrichment_table["group"] = group

        if use_adjusted_pval:
            fig = px.scatter(
                self.enrichment_table,
                x="dominant_enrichment_value_log2",
                y="dominant_adjusted_p_value_log10_abs",
                hover_name=self.enrichment_table.index,
                color=self.enrichment_table["group"],
                text=self.enrichment_table.index,
            )
            fig.add_hline(
                y=1.3,
                line_dash="dash",
                line_width=1,
                line_color="black",
                annotation_text="p ≤ 0.05",
                annotation_position="top right",
                annotation_font_size=15,
                annotation_font_color="black",
                opacity=0.7
            )
        else:
            fig = px.scatter(
                self.enrichment_table,
                x="dominant_enrichment_value_log2",
                y="dominant_p_value_log10_abs",
                hover_name=self.enrichment_table.index,
                color=self.enrichment_table["group"],
                text=self.enrichment_table.index,
            )
            fig.add_hline(
                y=1.3,
                line_dash="dash",
                line_width=1,
                line_color="black",
                annotation_text="p ≤ 0.05",
                annotation_position="top right",
                annotation_font_size=15,
                annotation_font_color="black",
                opacity=0.7
            )
        fig.update_traces(textfont_size=7, textposition="middle right") 
        fig.show()
        return fig