import pandas as pd
import numpy as np
from os import listdir
from sklearn.manifold import MDS, TSNE
from sklearn import preprocessing
        
import plotly.express as px
import plotly.graph_objects as go
from functions import get_distances, get_distance


class Comparison:

    def __init__(self, DATA_PATH: str, precomputed_distance_matrix: pd.DataFrame,) -> None:

        # List with path to all files
        table_paths = [f for f in listdir(DATA_PATH) if f != 'my-directory-list.csv']
        # Read all the enrichment tables
        tables = [pd.read_csv(f'{DATA_PATH}{file}') for file in table_paths]

        # Split table into individual tables based on the experiment
        self.coordinates = {f'{experiment}_{dose}': table[table.experiment_name == experiment][table.dose == dose]
                    for table in tables for experiment in table.experiment_name.unique() for dose in table[table.experiment_name == experiment].dose.unique()}

        # Index for the plot
        self.color_index = [f'{experiment}_{dose}' for table in tables for experiment in table.experiment_name.unique() for dose in table[table.experiment_name == experiment].dose.unique()]
        self.text_index = [dose for table in tables for experiment in table.experiment_name.unique() for dose in table[table.experiment_name == experiment].dose.unique()]
        self.precomputed_distance_matrix = precomputed_distance_matrix
        self.precomputed_distance_matrix.columns = range(precomputed_distance_matrix.columns.size)


    def get_comparison(self, input_sample: pd.DataFrame):
        coordinates_keys = list(self.coordinates.keys())
        # calculate the distances from input to every exepriment
        input_distances = pd.DataFrame([get_distances(input_sample[['dominant_enrichment_value_log2', 'dominant_p_value_log10_abs']], self.coordinates[sample][[
                                    'dominant_enrichment_value_log2', 'dominant_p_value_log10_abs']]) for sample in coordinates_keys])
        # Set the index
        input_distances.index = [
            f"{row};{len(coordinates_keys)}" for row, _ in enumerate(coordinates_keys)]

        input_distances.columns = range(input_distances.columns.size)
        
        # Add distances to input to precomputed matrix
        dist_matrix = pd.concat([input_distances, self.precomputed_distance_matrix])
        # Scale the distances
        dist_matrix = pd.DataFrame(preprocessing.scale(dist_matrix), index=dist_matrix.index)

        # Sum the columns
        len_keys = len(coordinates_keys) + 1
        dist_matrix_scaled = np.zeros([len_keys, len_keys])

        for row in range(len_keys):
            for column in range(row + 1, len_keys):
                dist_matrix_scaled[row, column] = np.sum(
                    dist_matrix.loc[f"{row};{column}", :])
                dist_matrix_scaled[column, row] = np.sum(
                    dist_matrix.loc[f"{row};{column}", :])

        if dist_matrix_scaled.min().min() < 0:
            dist_matrix_scaled += abs(dist_matrix_scaled.min().min())

        X_transform = TSNE(n_components=2, learning_rate='auto', init="random",
                        perplexity=50, metric='precomputed', random_state=0).fit_transform(dist_matrix_scaled)

        fig1 = px.scatter(X_transform[:-1], x=0, y=1, color=self.color_index,
                        hover_name=self.text_index, opacity=0.6)

        fig1.update_traces(marker=dict(size=5, line=dict(width=0.5)),
                        textposition='top right')

        fig2 = px.scatter(X_transform[-1:], x=0, y=1, color_discrete_sequence=[
                        'red'], text=["Input"])
        fig2.update_traces(marker=dict(size=10, line=dict(width=2)),
                        textposition='top right')

        fig = go.Figure(data=fig1.data + fig2.data)
        fig.update_layout(title='t-SNE ALL DRUGS', xaxis_title=f't-SNE1', yaxis_title=f't-SNE2', template="none", showlegend=False, xaxis=dict(ticks="outside",
                        mirror=True, showline=True), yaxis=dict(ticks="outside", mirror=True, showline=True), legend=dict(title="Drug"), width=800, height=800)
        fig.update_yaxes(scaleanchor="x", scaleratio=1)
        return fig