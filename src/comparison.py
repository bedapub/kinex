import pandas as pd
import numpy as np
from sklearn.manifold import MDS, TSNE
from sklearn import preprocessing
        
import plotly.express as px
import plotly.graph_objects as go
from functions import get_distances
from data import get_experiments


class Comparison:

    def __init__(self) -> None:
        self.experiments = get_experiments()
    
    def get_comparison(self, input_data: pd.DataFrame, method: str = "tSNE", center: bool = True):

        method = method.upper()
        supportedMethods = ["TSNE", "MDS"]

        if method not in supportedMethods:
            raise ValueError(
                f"The method {method} is not supported. Supported methods: " + supportedMethods)

        # Add Input and DMSO
        nExperiments = len(self.experiments)

        self.experiments.append({"id": nExperiments, "experiment_name": "DMSO", "dose": 0,
                            "dominant_enrichment_value_log2": np.zeros(303), "dominant_p_value_log10_abs": np.zeros(303)})

        self.experiments.append({"id": nExperiments + 1, "experiment_name": input_data.experiment_name.unique()[0], "dose": input_data.dose.unique()[
            0], "dominant_enrichment_value_log2": np.array(input_data["dominant_enrichment_value_log2"]), "dominant_p_value_log10_abs": np.array(input_data["dominant_p_value_log10_abs"])})

        # Calculate how many distances need to be computed for experiments + DMS0 + Input
        nExperiments += 2
        nDistances = (int)((nExperiments)*(nExperiments - 1) / 2)

        # Compute the distances
        distancesIndex = np.empty(nDistances, dtype=tuple)
        distances = np.empty((nDistances, 303))

        # Color and text indices
        colorIndex = np.empty(nExperiments, dtype=object)
        textIndex = np.empty(nExperiments, dtype=object)

        counter = 0
        for id1, experiment1 in enumerate(self.experiments):
            # Save the values for plot indexing
            colorIndex[id1] = experiment1["experiment_name"]
            textIndex[id1] = experiment1["dose"]

            for experiment2 in self.experiments[id1+1:]:
                distances[counter][:] = get_distances(
                    experiment1=experiment1, experiment2=experiment2)
                distancesIndex[counter] = (experiment1["id"], experiment2["id"])
                counter += 1

        # Scale the data
        distances = preprocessing.scale(distances)

        # Shift values to remove negative values
        minValue = distances.min().min()
        if minValue < 0:
            distances -= minValue

        # Create a dissimilarity matrix
        dissimilarityMatrix = np.zeros((nExperiments, nExperiments))

        for id, distance in enumerate(distances):
            ids = distancesIndex[id]

            dist = distance.sum()
            dissimilarityMatrix[ids[0]][ids[1]] = dist
            dissimilarityMatrix[ids[1]][ids[0]] = dist

        if method == 'TSNE':
            X_transform = TSNE(n_components=2, learning_rate='auto', init="random",
                            perplexity=50, metric='precomputed', random_state=0).fit_transform(dissimilarityMatrix)
        elif method == 'MDS':
            X_transform = MDS(n_components=2, dissimilarity='precomputed', normalized_stress="auto",
                            random_state=0).fit_transform(dissimilarityMatrix)

        if center:
            X_transform -= X_transform[-2]

        fig1 = px.scatter(X_transform[:-1], x=0, y=1, color=colorIndex[:-1],
                        hover_name=textIndex[:-1], opacity=0.6)

        fig1.update_traces(marker=dict(size=5, line=dict(width=0.5)),
                        textposition='top right')

        fig2 = px.scatter(X_transform[-1:], x=0, y=1, color_discrete_sequence=[
            'red'], text=["Input"])
        fig2.update_traces(marker=dict(size=10, line=dict(width=2)),
                        textposition='top right')

        fig = go.Figure(data=fig1.data + fig2.data)
        fig.update_layout(title=f"{method} ALL DRUGS", xaxis_title=f"{method}1", yaxis_title=f"{method}2", template="none", showlegend=False, xaxis=dict(ticks="outside",
                                                                                                                                                        mirror=True, showline=True), yaxis=dict(ticks="outside", mirror=True, showline=True), legend=dict(title="Drug"), width=800, height=800)
        fig.update_yaxes(scaleanchor="x", scaleratio=1)

        return fig