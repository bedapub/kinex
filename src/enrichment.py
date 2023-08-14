import numpy as np
import pandas as pd
import plotly.express as px
from statsmodels.stats.multitest import multipletests

from table2x2 import Table2x2


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

        kinase_group = {
            'AAK1': 'Other',
            'ACVR2A': 'TKL',
            'ACVR2B': 'TKL',
            'AKT1': 'AGC',
            'AKT2': 'AGC',
            'AKT3': 'AGC',
            'ALK2': 'TKL',
            'ALK4': 'TKL',
            'ALPHAK3': 'Alpha',
            'AMPKA1': 'CAMK',
            'AMPKA2': 'CAMK',
            'ANKRD3': 'TKL',
            'ASK1': 'STE',
            'ATM': 'PIKK',
            'ATR': 'PIKK',
            'AURA': 'Other',
            'AURB': 'Other',
            'AURC': 'Other',
            'BCKDK': 'PDHK',
            'BIKE': 'Other',
            'BMPR1A': 'TKL',
            'BMPR1B': 'TKL',
            'BMPR2': 'TKL',
            'BRAF': 'TKL',
            'BRSK1': 'CAMK',
            'BRSK2': 'CAMK',
            'BUB1': 'Other',
            'CAMK1A': 'CAMK',
            'CAMK1B': 'CAMK',
            'CAMK1D': 'CAMK',
            'CAMK1G': 'CAMK',
            'CAMK2A': 'CAMK',
            'CAMK2B': 'CAMK',
            'CAMK2D': 'CAMK',
            'CAMK2G': 'CAMK',
            'CAMK4': 'CAMK',
            'CAMKK1': 'Other',
            'CAMKK2': 'Other',
            'CAMLCK': 'CAMK',
            'CDC7': 'Other',
            'CDK1': 'CMGC',
            'CDK10': 'CMGC',
            'CDK12': 'CMGC',
            'CDK13': 'CMGC',
            'CDK14': 'CMGC',
            'CDK16': 'CMGC',
            'CDK17': 'CMGC',
            'CDK18': 'CMGC',
            'CDK19': 'CMGC',
            'CDK2': 'CMGC',
            'CDK3': 'CMGC',
            'CDK4': 'CMGC',
            'CDK5': 'CMGC',
            'CDK6': 'CMGC',
            'CDK7': 'CMGC',
            'CDK8': 'CMGC',
            'CDK9': 'CMGC',
            'CDKL1': 'CMGC',
            'CDKL5': 'CMGC',
            'CHAK1': 'Alpha',
            'CHAK2': 'Alpha',
            'CHK1': 'CAMK',
            'CHK2': 'CAMK',
            'CK1A': 'CK1',
            'CK1A2': 'CK1',
            'CK1D': 'CK1',
            'CK1E': 'CK1',
            'CK1G1': 'CK1',
            'CK1G2': 'CK1',
            'CK1G3': 'CK1',
            'CK2A1': 'Other',
            'CK2A2': 'Other',
            'CLK1': 'CMGC',
            'CLK2': 'CMGC',
            'CLK3': 'CMGC',
            'CLK4': 'CMGC',
            'COT': 'STE',
            'CRIK': 'AGC',
            'DAPK1': 'CAMK',
            'DAPK2': 'CAMK',
            'DAPK3': 'CAMK',
            'DCAMKL1': 'CAMK',
            'DCAMKL2': 'CAMK',
            'DLK': 'TKL',
            'DMPK1': 'AGC',
            'DNAPK': 'PIKK',
            'DRAK1': 'CAMK',
            'DSTYK': 'TKL',
            'DYRK1A': 'CMGC',
            'DYRK1B': 'CMGC',
            'DYRK2': 'CMGC',
            'DYRK3': 'CMGC',
            'DYRK4': 'CMGC',
            'EEF2K': 'Alpha',
            'ERK1': 'CMGC',
            'ERK2': 'CMGC',
            'ERK5': 'CMGC',
            'ERK7': 'CMGC',
            'FAM20C': 'FAM20',
            'GAK': 'Other',
            'GCK': 'STE',
            'GCN2': 'Other',
            'GRK1': 'AGC',
            'GRK2': 'AGC',
            'GRK3': 'AGC',
            'GRK4': 'AGC',
            'GRK5': 'AGC',
            'GRK6': 'AGC',
            'GRK7': 'AGC',
            'GSK3A': 'CMGC',
            'GSK3B': 'CMGC',
            'HASPIN': 'Other',
            'HGK': 'STE',
            'HIPK1': 'CMGC',
            'HIPK2': 'CMGC',
            'HIPK3': 'CMGC',
            'HIPK4': 'CMGC',
            'HPK1': 'STE',
            'HRI': 'Other',
            'HUNK': 'CMGC',
            'ICK': 'CMGC',
            'IKKA': 'Other',
            'IKKB': 'Other',
            'IKKE': 'Other',
            'IRAK1': 'TKL',
            'IRAK4': 'TKL',
            'IRE1': 'Other',
            'IRE2': 'Other',
            'JNK1': 'CMGC',
            'JNK2': 'CMGC',
            'JNK3': 'CMGC',
            'KHS1': 'STE',
            'KHS2': 'STE',
            'KIS': 'Other',
            'LATS1': 'AGC',
            'LATS2': 'AGC',
            'LKB1': 'CAMK',
            'LOK': 'STE',
            'LRRK2': 'TKL',
            'MAK': 'CMGC',
            'MAP3K15': 'STE',
            'MAPKAPK2': 'CAMK',
            'MAPKAPK3': 'CAMK',
            'MAPKAPK5': 'CAMK',
            'MARK1': 'CAMK',
            'MARK2': 'CAMK',
            'MARK3': 'CAMK',
            'MARK4': 'CAMK',
            'MASTL': 'AGC',
            'MEK1': 'STE',
            'MEK2': 'STE',
            'MEK5': 'STE',
            'MEKK1': 'STE',
            'MEKK2': 'STE',
            'MEKK3': 'STE',
            'MEKK6': 'STE',
            'MELK': 'CAMK',
            'MINK': 'STE',
            'MLK1': 'TKL',
            'MLK2': 'TKL',
            'MLK3': 'TKL',
            'MLK4': 'STE',
            'MNK1': 'CAMK',
            'MNK2': 'CAMK',
            'MOK': 'CMGC',
            'MOS': 'Other',
            'MPSK1': 'Other',
            'MRCKA': 'AGC',
            'MRCKB': 'AGC',
            'MSK1': 'AGC',
            'MSK2': 'AGC',
            'MST1': 'STE',
            'MST2': 'STE',
            'MST3': 'STE',
            'MST4': 'STE',
            'MTOR': 'PIKK',
            'MYLK4': 'CAMK',
            'MYO3A': 'STE',
            'MYO3B': 'STE',
            'NDR1': 'AGC',
            'NDR2': 'AGC',
            'NEK1': 'Other',
            'NEK11': 'Other',
            'NEK2': 'Other',
            'NEK3': 'Other',
            'NEK4': 'Other',
            'NEK5': 'Other',
            'NEK6': 'Other',
            'NEK7': 'Other',
            'NEK8': 'Other',
            'NEK9': 'Other',
            'NIK': 'STE',
            'NIM1': 'CAMK',
            'NLK': 'CMGC',
            'NUAK1': 'CAMK',
            'NUAK2': 'CAMK',
            'OSR1': 'STE',
            'P38A': 'CMGC',
            'P38B': 'CMGC',
            'P38D': 'CMGC',
            'P38G': 'CMGC',
            'P70S6K': 'AGC',
            'P70S6KB': 'AGC',
            'P90RSK': 'AGC',
            'PAK1': 'STE',
            'PAK2': 'STE',
            'PAK3': 'STE',
            'PAK4': 'STE',
            'PAK5': 'STE',
            'PAK6': 'STE',
            'PASK': 'CAMK',
            'PBK': 'Other',
            'PDHK1': 'PDHK',
            'PDHK4': 'PDHK',
            'PDK1': 'AGC',
            'PERK': 'Other',
            'PHKG1': 'CAMK',
            'PHKG2': 'CAMK',
            'PIM1': 'CAMK',
            'PIM2': 'CAMK',
            'PIM3': 'CAMK',
            'PINK1': 'Other',
            'PKACA': 'AGC',
            'PKACB': 'AGC',
            'PKACG': 'AGC',
            'PKCA': 'AGC',
            'PKCB': 'AGC',
            'PKCD': 'AGC',
            'PKCE': 'AGC',
            'PKCG': 'AGC',
            'PKCH': 'AGC',
            'PKCI': 'AGC',
            'PKCT': 'AGC',
            'PKCZ': 'AGC',
            'PKG1': 'AGC',
            'PKG2': 'AGC',
            'PKN1': 'AGC',
            'PKN2': 'AGC',
            'PKN3': 'AGC',
            'PKR': 'Other',
            'PLK1': 'Other',
            'PLK2': 'Other',
            'PLK3': 'Other',
            'PLK4': 'Other',
            'PRKD1': 'CAMK',
            'PRKD2': 'CAMK',
            'PRKD3': 'CAMK',
            'PRKX': 'AGC',
            'PRP4': 'CMGC',
            'PRPK': 'Other',
            'QIK': 'CAMK',
            'QSK': 'CAMK',
            'RAF1': 'TKL',
            'RIPK1': 'TKL',
            'RIPK2': 'TKL',
            'RIPK3': 'TKL',
            'ROCK1': 'AGC',
            'ROCK2': 'AGC',
            'RSK2': 'AGC',
            'RSK3': 'AGC',
            'RSK4': 'AGC',
            'SBK': 'Other',
            'SGK1': 'AGC',
            'SGK3': 'AGC',
            'SIK': 'CAMK',
            'SKMLCK': 'CAMK',
            'SLK': 'STE',
            'SMG1': 'PIKK',
            'SMMLCK': 'CAMK',
            'SNRK': 'CAMK',
            'SRPK1': 'CMGC',
            'SRPK2': 'CMGC',
            'SRPK3': 'CMGC',
            'SSTK': 'CAMK',
            'STK33': 'CAMK',
            'STLK3': 'STE',
            'TAK1': 'TKL',
            'TAO1': 'STE',
            'TAO2': 'STE',
            'TAO3': 'STE',
            'TBK1': 'Other',
            'TGFBR1': 'TKL',
            'TGFBR2': 'TKL',
            'TLK1': 'Other',
            'TLK2': 'Other',
            'TNIK': 'STE',
            'TSSK1': 'CAMK',
            'TSSK2': 'CAMK',
            'TTBK1': 'CK1',
            'TTBK2': 'CK1',
            'TTK': 'Other',
            'ULK1': 'Other',
            'ULK2': 'Other',
            'VRK1': 'CK1',
            'VRK2': 'CK1',
            'WNK1': 'Other',
            'WNK3': 'Other',
            'WNK4': 'Other',
            'YANK2': 'AGC',
            'YANK3': 'AGC',
            'YSK1': 'STE',
            'YSK4': 'STE',
            'ZAK': 'STE'
        }

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