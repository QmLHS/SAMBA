import pandas as pd
import matplotlib.pyplot as plt
import umap

###############
#### NCBI #####
###############

ncbi = pd.read_csv('ncbi_filtd_UMAP.csv')
print('Data loaded.')

ncbi = ncbi.fillna(0)

reducer = umap.UMAP()

embedding = reducer.fit_transform(ncbi.iloc[:, 10:])
print(embedding.shape)

# Plots:
labels = ncbi['Kingdom']
unique_species = labels.unique()
species_to_num = {species: i for i, species in enumerate(unique_species)}
numeric_labels = labels.map(species_to_num)

plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=numeric_labels,
    cmap='viridis',
    s=50,
    alpha=0.8
)
cbar = plt.colorbar()
cbar.set_label('Label')

plt.gca().set_aspect('equal', 'datalim')
plt.title('UMAP Projection of the Dataset (Kingdom)')
plt.savefig('ncbi_Kingdom_UMAP.png')



labels = ncbi['Phylum']
unique_species = labels.unique()
species_to_num = {species: i for i, species in enumerate(unique_species)}
numeric_labels = labels.map(species_to_num)

plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=numeric_labels,
    cmap='viridis',
    s=50, 
    alpha=0.8  
)
cbar = plt.colorbar()
cbar.set_label('Label')

plt.gca().set_aspect('equal', 'datalim')
plt.title('UMAP Projection of the Dataset (Phylum)')
plt.savefig('ncbi_Phylum_UMAP.png')



labels = ncbi['Class']
unique_species = labels.unique()
species_to_num = {species: i for i, species in enumerate(unique_species)}
numeric_labels = labels.map(species_to_num)

plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=numeric_labels,
    cmap='viridis',
    s=50, 
    alpha=0.8 
)
cbar = plt.colorbar()
cbar.set_label('Label')  

plt.gca().set_aspect('equal', 'datalim')
plt.title('UMAP Projection of the Dataset (Class)')
plt.savefig('ncbi_Class_UMAP.png')



labels = ncbi['Order']
unique_species = labels.unique()
species_to_num = {species: i for i, species in enumerate(unique_species)}
numeric_labels = labels.map(species_to_num)

plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=numeric_labels,
    cmap='viridis',
    s=50,
    alpha=0.8 
)
cbar = plt.colorbar()
cbar.set_label('Label')  # Set label for the colorbar

plt.gca().set_aspect('equal', 'datalim')
plt.title('UMAP Projection of the Dataset (Order)')
plt.savefig('ncbi_Order_UMAP.png')



labels = ncbi['Family']
unique_species = labels.unique()
species_to_num = {species: i for i, species in enumerate(unique_species)}
numeric_labels = labels.map(species_to_num)

plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=numeric_labels,
    cmap='viridis',
    s=50,
    alpha=0.8
)
cbar = plt.colorbar()
cbar.set_label('Label')

plt.gca().set_aspect('equal', 'datalim')
plt.title('UMAP Projection of the Dataset (Family)', fontsize=24)
plt.savefig('ncbi_Family_UMAP.png')



labels = ncbi['Genus']
unique_species = labels.unique()
species_to_num = {species: i for i, species in enumerate(unique_species)}
numeric_labels = labels.map(species_to_num)

plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=numeric_labels,
    cmap='viridis',
    s=50,
    alpha=0.8
)
cbar = plt.colorbar()
cbar.set_label('Label')

plt.gca().set_aspect('equal', 'datalim')
plt.title('UMAP Projection of the Dataset (Genus)')
plt.savefig('ncbi_Genus_UMAP.png')



labels = ncbi['Species']
unique_species = labels.unique()
species_to_num = {species: i for i, species in enumerate(unique_species)}
numeric_labels = labels.map(species_to_num)

plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=numeric_labels,
    cmap='viridis',
    s=50,
    alpha=0.8
)

cbar = plt.colorbar()
cbar.set_label('Label')

plt.gca().set_aspect('equal', 'datalim')
plt.title('UMAP Projection of the Dataset (Species)')
plt.savefig('ncbi_Species_UMAP.png')

del ncbi
print('NCBI plots saved.')


###############
#### Bold #####
###############


bold = pd.read_csv('bold_filtd_UMAP.csv')
print('Data loaded.')

bold = bold.fillna(0)

reducer = umap.UMAP()

embedding = reducer.fit_transform(bold.iloc[:, 10:])
print(embedding.shape)

# Plots:
labels = bold['Kingdom']
unique_species = labels.unique()
species_to_num = {species: i for i, species in enumerate(unique_species)}
numeric_labels = labels.map(species_to_num)

plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=numeric_labels,
    cmap='viridis',
    s=50,
    alpha=0.8
)
cbar = plt.colorbar()
cbar.set_label('Label')

plt.gca().set_aspect('equal', 'datalim')
plt.title('UMAP Projection of the Dataset (Kingdom)')
plt.savefig('bold_Kingdom_UMAP.png')



labels = bold['Phylum']
unique_species = labels.unique()
species_to_num = {species: i for i, species in enumerate(unique_species)}
numeric_labels = labels.map(species_to_num)

plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=numeric_labels,
    cmap='viridis',
    s=50, 
    alpha=0.8  
)
cbar = plt.colorbar()
cbar.set_label('Label')

plt.gca().set_aspect('equal', 'datalim')
plt.title('UMAP Projection of the Dataset (Phylum)')
plt.savefig('bold_Phylum_UMAP.png')



labels = bold['Class']
unique_species = labels.unique()
species_to_num = {species: i for i, species in enumerate(unique_species)}
numeric_labels = labels.map(species_to_num)

plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=numeric_labels,
    cmap='viridis',
    s=50, 
    alpha=0.8 
)
cbar = plt.colorbar()
cbar.set_label('Label')  

plt.gca().set_aspect('equal', 'datalim')
plt.title('UMAP Projection of the Dataset (Class)')
plt.savefig('bold_Class_UMAP.png')



labels = bold['Order']
unique_species = labels.unique()
species_to_num = {species: i for i, species in enumerate(unique_species)}
numeric_labels = labels.map(species_to_num)

plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=numeric_labels,
    cmap='viridis',
    s=50,
    alpha=0.8 
)
cbar = plt.colorbar()
cbar.set_label('Label')  # Set label for the colorbar

plt.gca().set_aspect('equal', 'datalim')
plt.title('UMAP Projection of the Dataset (Order)')
plt.savefig('bold_Order_UMAP.png')



labels = bold['Family']
unique_species = labels.unique()
species_to_num = {species: i for i, species in enumerate(unique_species)}
numeric_labels = labels.map(species_to_num)

plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=numeric_labels,
    cmap='viridis',
    s=50,
    alpha=0.8
)
cbar = plt.colorbar()
cbar.set_label('Label')

plt.gca().set_aspect('equal', 'datalim')
plt.title('UMAP Projection of the Dataset (Family)', fontsize=24)
plt.savefig('bold_Family_UMAP.png')



labels = bold['Genus']
unique_species = labels.unique()
species_to_num = {species: i for i, species in enumerate(unique_species)}
numeric_labels = labels.map(species_to_num)

plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=numeric_labels,
    cmap='viridis',
    s=50,
    alpha=0.8
)
cbar = plt.colorbar()
cbar.set_label('Label')

plt.gca().set_aspect('equal', 'datalim')
plt.title('UMAP Projection of the Dataset (Genus)')
plt.savefig('bold_Genus_UMAP.png')



labels = bold['Species']
unique_species = labels.unique()
species_to_num = {species: i for i, species in enumerate(unique_species)}
numeric_labels = labels.map(species_to_num)

plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=numeric_labels,
    cmap='viridis',
    s=50,
    alpha=0.8
)

cbar = plt.colorbar()
cbar.set_label('Label')

plt.gca().set_aspect('equal', 'datalim')
plt.title('UMAP Projection of the Dataset (Species)')
plt.savefig('bold_Species_UMAP.png')

print('Bold plots saved.')
