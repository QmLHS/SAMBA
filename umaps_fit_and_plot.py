import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from umap import umap_
from sklearn.feature_extraction.text import TfidfTransformer



def compute_tfidf_remove_0s_compute_UMAP_embedding(_df, name: str):
    transformer = TfidfTransformer()
    tfidf_matrix = transformer.fit_transform(_df.values.T).T
    tfidf_df = pd.DataFrame(
        tfidf_matrix.toarray(),
        index=_df.index,
        columns=_df.columns
    )
    # Now remove columns with all 0s
    non_zero_cols = (tfidf_df != 0).any(axis=0)
    cleaned_df = tfidf_df.loc[:, non_zero_cols]

    removed_count = len(tfidf_df.columns) - len(cleaned_df.columns)
    print(f"Removed {removed_count} columns that contained all 0s.")

    reducer = umap_.UMAP(n_neighbors=200, min_dist=0.1, metric='euclidean', random_state=42)
    embedding = reducer.fit_transform(cleaned_df)
    print(embedding.shape)

    df_ = pd.DataFrame(embedding, columns=['x', 'y'])
    df_.to_csv(f'Umap_{name}.csv', index=False)
    print(f'Umap_{name}.csv written')



###############
#### NCBI #####
###############

ncbi = pd.read_csv('ncbi_filtd_UMAP.csv', low_memory=False)
print('Data loaded.')

ncbi = ncbi.fillna(0)
# All Kingdoms
compute_tfidf_remove_0s_compute_UMAP_embedding(ncbi.iloc[:, 10:], 'ncbi_all')

# Only Metazoa
ncbi_ = ncbi[ncbi['Kingdom']=='Metazoa']
compute_tfidf_remove_0s_compute_UMAP_embedding(ncbi_.iloc[:, 10:], 'ncbi_metazoa')

# Only Pseudomonadati
ncbi_ = ncbi[ncbi['Kingdom']=='Pseudomonadati']
compute_tfidf_remove_0s_compute_UMAP_embedding(ncbi_.iloc[:, 10:], 'ncbi_pseudomonadati')

# Only Viridiplantae
ncbi_ = ncbi[ncbi['Kingdom']=='Viridiplantae']
compute_tfidf_remove_0s_compute_UMAP_embedding(ncbi_.iloc[:, 10:], 'ncbi_viridiplantae')

# Only Fungi
ncbi_ = ncbi[ncbi['Kingdom']=='Fungi']
compute_tfidf_remove_0s_compute_UMAP_embedding(ncbi_.iloc[:, 10:], 'ncbi_fungi')

###############
#### Bold #####
###############


bold = pd.read_csv('bold_filtd_UMAP.csv', low_memory=False)
print('Data loaded.')

bold = bold.fillna(0)

# All Kingdoms
compute_tfidf_remove_0s_compute_UMAP_embedding(bold.iloc[:, 10:], 'bold_all')

# Only Animalia
bold_ = bold[bold['Kingdom']=='Animalia']
compute_tfidf_remove_0s_compute_UMAP_embedding(bold_.iloc[:, 10:], 'bold_animalia')

# Only Protozoa
bold_ = bold[bold['Kingdom']=='Protozoa']
compute_tfidf_remove_0s_compute_UMAP_embedding(bold_.iloc[:, 10:], 'bold_protozoa')

# Only Fungi
bold_ = bold[bold['Kingdom']=='Fungi']
compute_tfidf_remove_0s_compute_UMAP_embedding(bold_.iloc[:, 10:], 'bold_fungi')

