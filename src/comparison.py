import pandas as pd
import numpy as np
from sklearn.manifold import MDS, TSNE
from sklearn import preprocessing
from os import listdir
from umap import UMAP

import plotly.express as px
import plotly.graph_objects as go
from kinex.functions import get_distances
from kinex.data import get_experiments

class Comparison:
    """
    Methods
    -------
    get_comparison(self, input_data: pd.DataFrame = None, data_path: str = None, method: str = "tSNE")
        Compares samples
    """

    def __init__(self) -> None:
        self.experiments = get_experiments()

    def get_comparison(self, input_data: pd.DataFrame = None, data_path: str = None, method: str = "MDS", center: bool = True):
        """
        Method to either compare the input with a collection of drugs from Zecha et al. or compare multiple samples between eachother

        Attributes
        ----------
        data_path:
            Can be eitehr the path to the collection of drugs profiles (a json file given with the paper) 
            Or can be a path to the user generated Kinex enrichment tables. Must contain enrichment values and p values. The name of the samples will be given by the filename.     
        input_data: pandas.DataFrame  
            Optional
            A DataFrame representing the Kinex enrichment table that will be used to compare it to the collection of drugs from Zecha et al.
        method: str
            Dimensionality redcution method which can be either "umap", "tsne" or "mds". By default: MDS
        """

        method = method.upper()
        supportedMethods = ["TSNE", "MDS", "UMAP"]

        if method not in supportedMethods:
            raise ValueError(
                f"The method {method} is not supported. Supported methods: " + supportedMethods)

        if input_data is None and data_path is None:
            raise ValueError(
                "Please give at least the input or the path to the experiments")

        if data_path is None:
            # Add Input and DMSO
            nExperiments = len(self.experiments)

            self.experiments.append({"id": nExperiments, "experiment_name": "DMSO", "dose": 0,
                                    "dominant_enrichment_value_log2": np.zeros(303), "dominant_p_value_log10_abs": np.zeros(303)})

            self.experiments.append({"id": nExperiments + 1, "experiment_name": "Input", "dose": "Input", "dominant_enrichment_value_log2": np.array(
                input_data["dominant_enrichment_value_log2"]), "dominant_p_value_log10_abs": np.array(input_data["dominant_p_value_log10_abs"])})

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
                    distancesIndex[counter] = (
                        experiment1["id"], experiment2["id"])
                    counter += 1

            # Scale the data
            distances = preprocessing.scale(distances, with_mean=False)

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
            elif method == 'UMAP':
                reducer = UMAP(n_neighbors=30)
                X_transform = reducer.fit_transform(dissimilarityMatrix)

            scaler = preprocessing.MinMaxScaler(feature_range=(-1, 1))
            X_transform = scaler.fit_transform(X_transform)

            X_transform -= X_transform[-2]

            fig1 = px.scatter(X_transform[:-1], x=0, y=1, color=colorIndex[:-1],
                              hover_name=textIndex[:-1], opacity=0.75)

            fig1.update_traces(marker=dict(size=5, line=dict(width=0.5)),
                               textposition='top right')

            fig2 = px.scatter(X_transform[-1:], x=0, y=1, color_discrete_sequence=[
                'red'], text=["Input"])
            fig2.update_traces(marker=dict(size=10, line=dict(width=2)),
                               textposition='top right')

            fig = go.Figure(data=fig1.data + fig2.data)
            fig.update_layout(title=f"{method} ALL DRUGS", xaxis_title=f"{method}1", yaxis_title=f"{method}2", template="none", showlegend=False, xaxis=dict(ticks="outside",
                                                                                                                                                             mirror=True, showline=True), yaxis=dict(ticks="outside", mirror=True, showline=True), legend=dict(title="Drug"), width=600, height=600)
            fig.update_yaxes(scaleanchor="x", scaleratio=1)

            return fig

        else:
            table_paths = []
            index = []

            for file in listdir(data_path):
                if (file != "my-directory-list.csv"):
                    table_paths.append(file)
                    index.append(file.split(".csv")[0])

            index.append("DMSO")
            # Empty array for tables
            tables = np.empty(len(table_paths) + 1, dtype=pd.DataFrame)
            # Add DMSO to the tables
            tables[len(tables) - 1] = pd.DataFrame({"dominant_enrichment_value_log2": np.zeros(
                303), "dominant_p_value_log10_abs": np.zeros(303)})

            for id, path in enumerate(table_paths):
                tables[id] = pd.read_csv(f"{data_path}/{path}")

            nExperiments = len(tables)
            nDistances = (int)((nExperiments)*(nExperiments - 1) / 2)

            # Compute the distances
            distancesIndex = np.empty(nDistances, dtype=tuple)
            distances = np.empty((nDistances, 303))

            counter = 0
            for id1, table1 in enumerate(tables):
                for id2, table2 in enumerate(tables[id1+1:]):
                    distances[counter][:] = get_distances(
                        experiment1=table1, experiment2=table2)
                    distancesIndex[counter] = (id1, id1+1 + id2)
                    counter += 1

            # Scale the distances
            distances = preprocessing.scale(distances, with_mean=False)

            # Create a dissimilarity matrix
            dissimilarityMatrix = np.zeros((nExperiments, nExperiments))

            for id, distance in enumerate(distances):
                ids = distancesIndex[id]
                dist = distance.sum()
                dissimilarityMatrix[ids[0]][ids[1]] = dist
                dissimilarityMatrix[ids[1]][ids[0]] = dist

            if method == 'TSNE':
                X_transform = TSNE(n_components=2, learning_rate='auto', init="random",
                                   perplexity=3, metric='precomputed', random_state=0).fit_transform(dissimilarityMatrix)
            elif method == 'MDS':
                X_transform = MDS(n_components=2, dissimilarity='precomputed', normalized_stress="auto",
                                random_state=0).fit_transform(dissimilarityMatrix)
            elif method == 'UMAP':
                reducer = UMAP()
                X_transform = reducer.fit_transform(dissimilarityMatrix)
            
            # Set the points range
            scaler = preprocessing.MinMaxScaler(feature_range=(-1, 1))
            X_transform = scaler.fit_transform(X_transform)

            # Center to DMSO
            X_transform -= X_transform[-1]

            fig = px.scatter(X_transform, x=0, y=1, color=index,
                             hover_name=index, opacity=1)
            fig.update_layout(title=f"{method} ALL DRUGS", xaxis_title=f"{method}1", yaxis_title=f"{method}2", template="none", showlegend=False, xaxis=dict(ticks="outside",
                                                                                                                                                             mirror=True, showline=True), yaxis=dict(ticks="outside", mirror=True, showline=True), legend=dict(title="Drug"), width=600, height=600)
            fig.update_yaxes(scaleanchor="x", scaleratio=1)

            return fig
