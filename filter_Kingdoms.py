import pandas as pd
import numpy as np
from umap import umap_
from sklearn.feature_extraction.text import TfidfTransformer

print('load ncbi...')
ncbi = pd.read_csv('ncbi_filtd_UMAP.csv', low_memory=False)
print('ncbi loaded.')

mask = (ncbi['Kingdom'] == 'Metazoa') & ((ncbi['Phylum'] == 'Arthropoda') | (ncbi['Phylum'] == 'Chordata') | (ncbi['Phylum'] == 'Mollusca') | (ncbi['Phylum'] == 'Annelida'))
ncbi = ncbi[mask]
ncbi.to_csv('ncbi_metazoa_ACMA.csv', index=False)

print('ncbi done.')

# -------------------------------------------

print('load bold...')
bold = pd.read_csv('bold_filtd_UMAP.csv', low_memory=False)
print('bold loaded.')

mask = (bold['Kingdom'] == 'Metazoa') & ((bold['Phylum'] == 'Arthropoda') | (bold['Phylum'] == 'Chordata') | (bold['Phylum'] == 'Mollusca') | (bold['Phylum'] == 'Annelida'))
bold = bold[mask]
bold.to_csv('bold_metazoa_ACMA.csv', index=False)

print('bold done.')

# -------------------------------------------

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

# -------------------------------------------


ncbi = ncbi.fillna(0)
bold = bold.fillna(0)
compute_tfidf_remove_0s_compute_UMAP_embedding(ncbi.iloc[:, 10:], 'ncbi_metazoa_ACMA')
compute_tfidf_remove_0s_compute_UMAP_embedding(bold.iloc[:, 10:], 'bold_metazoa_ACMA')
