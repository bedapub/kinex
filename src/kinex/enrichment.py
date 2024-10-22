import numpy as np
import pandas as pd
import plotly.express as px
from statsmodels.stats.multitest import multipletests

from kinex.sequence import SequenceType
from kinex.table2x2 import Table2x2
from kinex.resources import get_ser_thr_family, get_ser_thr_family_colors, get_tyr_family, get_tyr_family_colors


class Enrichment:
    def __init__(self, sequence_type: SequenceType, all_kinases: set):
        self.sequence_type = sequence_type
        self.total_upregulated = 0
        self.total_downregulated = 0
        self.total_unregulated = 0
        self.regulation_list = []
        self.top15_kinases_list = []
        self.enrichment_table = pd.DataFrame(columns=['kinase', 'upregulated', 'downregulated', 'unregulated'])
        self.all_kinases = all_kinases

    def adjust_background_sites(self):
        if self.total_unregulated == 0:
            self.total_unregulated = np.min(
                [self.total_upregulated, self.total_downregulated]) / 2

    def fisher_statistics(self):
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
                downregulated_contingency_table = Table2x2(
                    [[downregulated_hit, self.total_downregulated - downregulated_hit],
                     [unregulated_hit, self.total_unregulated - unregulated_hit]],
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

        missing_kinases = list(self.all_kinases - set(self.enrichment_table.index))
        self.enrichment_table = self.enrichment_table.reindex(
            self.enrichment_table.index.union(missing_kinases), fill_value=0)

    def plot(self, use_adjusted_pval: bool = False):
        match self.sequence_type:
            case SequenceType.SER_THR:
                kinase_family = get_ser_thr_family()
                family_colors = get_ser_thr_family_colors()
            case SequenceType.TYR:
                kinase_family = get_tyr_family()
                family_colors = get_tyr_family_colors()
            case _:
                raise ValueError(f"Invalid sequence type")

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
            template="none"
        )

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
        return fig