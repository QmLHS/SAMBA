import pandas as pd
import sys
from pathlib import Path



## Importing of the necessary files: sequences.fasta and taxonomy.tsv
if len(sys.argv) != 3:
  print("Usage: script.py <fasta_file> <taxonomy_file>")
  sys.exit(1)

try:
  fasta_path = Path(sys.argv[1])
  taxa_path = Path(sys.argv[2])
  if not (fasta_path.exists() and taxa_path.exists()):
    raise FileNotFoundError("Input files not found")
except Exception as e:
  print(f"Error with input files: {e}")
  sys.exit(1)


# First thing we want to read the FASTA file
def read_fasta(fasta_file):
  ids = []
  sequences = []
  current_sequence = []
  with open(fasta_file, 'r') as f:
    for line in f:
      line = line.strip()
      if line.startswith('>'):
        if current_sequence:
          sequences.append(''.join(current_sequence))
        ids.append(line[1:])
        current_sequence = []
      else:
        current_sequence.append(line)                
  if current_sequence:
    sequences.append(''.join(current_sequence))
  df = pd.DataFrame({
    'Feature ID': ids,
    'Sequence': sequences
  })  
  return df

my_fasta = sys.argv[1]
seqs = read_fasta(my_fasta)


# Then, we want to read and standardize the taxonomy information so that ML 
# doesn't output "non-consistency" errors
def standardize_taxonomy(taxonomy_string):
  standard_ranks = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
  required_length = len(standard_ranks)
  # Extract the actual taxonomy part
  taxonomy_part = taxonomy_string.replace('tax=', '')
  ranks = taxonomy_part.split(';')
  ranks = [r.strip() for r in ranks]
  standardized_ranks = []
  # Process existing ranks
  for i in range(required_length):
    prefix = standard_ranks[i]
    if i < len(ranks) and ranks[i]:
      rank = ranks[i]
      # Add prefix if it's missing
      if not rank.startswith(prefix):
        rank = prefix + rank
      standardized_ranks.append(rank)
    else:
      # Always append empty rank with prefix
      standardized_ranks.append(prefix)
  # Ensure we have exactly 7 ranks
  assert len(standardized_ranks) == required_length, f"Expected {required_length} ranks, got {len(standardized_ranks)}"
  # Add back the 'tax=' prefix and join with semicolons
  return 'tax=' + ';'.join(standardized_ranks)

def read_taxonomy(tsv_file):
  df = pd.read_csv(tsv_file, sep='\t')
  df.columns = ['Feature ID', 'Taxon']
  df['Feature ID'] = df['Feature ID'].astype(str)
  # Only standardize if the taxonomy doesn't already start with 'tax='
  df['Taxon'] = df['Taxon'].apply(standardize_taxonomy)
  return df

my_taxa = sys.argv[2]
taxa = read_taxonomy(my_taxa)


## Creation of the DataFrame objects + Output in files compatible with QIIME2

# The first thing is to make a combined df so that we can remove sequences with
# a taxa label that only appear once.
# Process that needs to be repeated even after the splitting in train-test-val
# because a class with only 2 instances can end up split in train and test

combined_df = pd.merge(seqs, taxa, on='Feature ID', how='inner')

del seqs, taxa

# Inside the taxonomy file, we have to be sure that no class is present in
# only one record, so we retrieve the indices of records that appear only
# one time
def analyze_rare_taxa(df, min_count=2):
  taxa_counts = df['Taxon'].value_counts()
  rare_taxa = taxa_counts[taxa_counts < min_count]
  rare_taxa_df = df[df['Taxon'].isin(rare_taxa.index)]
  return rare_taxa_df['Feature ID'], len(rare_taxa)
 
taxa_id_to_remove, num_ids_to_remove = analyze_rare_taxa(combined_df)
print(f"Number of taxa labels removed from combined df: {num_ids_to_remove}")
combined_df_filtd = combined_df[~combined_df['Feature ID'].isin(taxa_id_to_remove)].copy()

del combined_df

# Now this df needs to be divided in train-test for classical ML and 
# train-val-test for Bert fine-tuning

sampled_ids_train_val = combined_df_filtd['Feature ID'].sample(frac=0.85, random_state=42)
combined_df_filtd_train_val = combined_df_filtd[combined_df_filtd['Feature ID'].isin(sampled_ids_train_val)].copy()
combined_df_filtd_test = combined_df_filtd[~combined_df_filtd['Feature ID'].isin(sampled_ids_train_val)].copy()

taxa_id_to_remove, num_ids_to_remove = analyze_rare_taxa(combined_df_filtd_train_val)
print(f"Number of taxa labels removed from train-val: {num_ids_to_remove}")
combined_df_filtd_train_val_filtd = combined_df_filtd_train_val[~combined_df_filtd_train_val['Feature ID'].isin(taxa_id_to_remove)].copy()

taxa_id_to_remove, num_ids_to_remove = analyze_rare_taxa(combined_df_filtd_test)
print(f"Number of taxa labels removed from test: {num_ids_to_remove}")
combined_df_filtd_test_filtd = combined_df_filtd_test[~combined_df_filtd_test['Feature ID'].isin(taxa_id_to_remove)].copy()

del combined_df_filtd_train_val, combined_df_filtd_test

## Output these df objects as files that can be read by the qiime artifacts

def format_sequence(sequence, line_length=80):
  """Set new line char after each 80 bp as QIIME2 would output it."""
  return '\n'.join(sequence[i:i+line_length] for i in range(0, len(sequence), line_length))

seqs_filtd_train = combined_df_filtd_train_val_filtd[['Feature ID', 'Sequence']].copy()
seqs_filtd_test = combined_df_filtd_test_filtd[['Feature ID', 'Sequence']].copy()
taxa_filtd_train = combined_df_filtd_train_val_filtd[['Feature ID', 'Taxon']].copy()
taxa_filtd_test = combined_df_filtd_test_filtd[['Feature ID', 'Taxon']].copy()

with open(f'{sys.argv[1][:-6]}-clean-train.fasta', 'w') as f:
  for _, row in seqs_filtd_train.iterrows():
    f.write(f">{row['Feature ID']}\n")
    f.write(format_sequence(row['Sequence']) + '\n')
with open(f'{sys.argv[1][:-6]}-clean-test.fasta', 'w') as f:
  for _, row in seqs_filtd_test.iterrows():
    f.write(f">{row['Feature ID']}\n")
    f.write(format_sequence(row['Sequence']) + '\n')

taxa_filtd_train.to_csv(f'{sys.argv[2][:-4]}-clean-train.tsv', sep='\t', index=False)
taxa_filtd_test.to_csv(f'{sys.argv[2][:-4]}-clean-test.tsv', sep='\t', index=False)

del seqs_filtd_train, seqs_filtd_test, taxa_filtd_train, taxa_filtd_test

## Bert needs validation set too
# The validation set that Bert needs is sampled from the training set,
# so the test remain the same

sampled_ids = combined_df_filtd_train_val_filtd['Feature ID'].sample(frac=0.15, random_state=42)
combined_df_filtd_val_filtd = combined_df_filtd_train_val_filtd[combined_df_filtd_train_val_filtd['Feature ID'].isin(sampled_ids)].copy()
combined_df_filtd_train_filtd = combined_df_filtd_train_val_filtd[~combined_df_filtd_train_val_filtd['Feature ID'].isin(sampled_ids)].copy()

taxa_id_to_remove, num_ids_to_remove = analyze_rare_taxa(combined_df_filtd_val_filtd)
print(f"Number of taxa labels removed from val: {num_ids_to_remove}")
combined_df_filtd_val_filtd_filtd = combined_df_filtd_val_filtd[~combined_df_filtd_val_filtd['Feature ID'].isin(taxa_id_to_remove)].copy()

taxa_id_to_remove, num_ids_to_remove = analyze_rare_taxa(combined_df_filtd_train_filtd)
print(f"Number of taxa labels removed from train: {num_ids_to_remove}")
combined_df_filtd_train_filtd_filtd = combined_df_filtd_train_filtd[~combined_df_filtd_train_filtd['Feature ID'].isin(taxa_id_to_remove)].copy()

del combined_df_filtd_train_val_filtd, combined_df_filtd_train_filtd, combined_df_filtd_val_filtd

# Also for the Bert fine-tuning process, we need labels as integer numbers,
# and a map file so that we can later bring the predictions to their
# original format

unique_taxa = combined_df_filtd['Taxon'].unique()
taxa_to_label = {taxon: i for i, taxon in enumerate(sorted(unique_taxa))}

combined_df_filtd_train_filtd_filtd['label'] = combined_df_filtd_train_filtd_filtd['Taxon'].map(taxa_to_label)
combined_df_filtd_val_filtd_filtd['label'] = combined_df_filtd_val_filtd_filtd['Taxon'].map(taxa_to_label)
combined_df_filtd_test_filtd['label'] = combined_df_filtd_test_filtd['Taxon'].map(taxa_to_label)

combined_df_filtd_train_filtd_filtd = combined_df_filtd_train_filtd_filtd.drop('Taxon', axis=1)
combined_df_filtd_val_filtd_filtd = combined_df_filtd_val_filtd_filtd.drop('Taxon', axis=1)
combined_df_filtd_test_filtd = combined_df_filtd_test_filtd.drop('Taxon', axis=1)

# Make files for the training
combined_df_filtd_train_filtd_filtd.to_csv(f'train-{sys.argv[2][9:13]}.csv', index=False)
combined_df_filtd_val_filtd_filtd.to_csv(f'val-{sys.argv[2][9:13]}.csv', index=False)
combined_df_filtd_test_filtd.to_csv(f'test-{sys.argv[2][9:13]}.csv', index=False)

# Write the mapping file for later
with open(f'taxa_labels_mapping.txt', 'a') as f:
  for taxon, label in taxa_to_label.items():
    f.write(f"{label}, {taxon}\n")
    
