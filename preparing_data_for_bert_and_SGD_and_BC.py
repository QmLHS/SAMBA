import pandas as pd


def read_fasta(fasta_file):
  sequences = {}
  current_id = None
  current_sequence = []
  with open(fasta_file, 'r') as f:
    for line in f:
      line = line.strip()
      if line.startswith('>'):
        if current_id is not None:
          sequences[current_id] = ''.join(current_sequence)
        current_id = line[1:]
        current_sequence = []
      else:
        current_sequence.append(line)
  if current_id is not None:
    sequences[current_id] = ''.join(current_sequence)
  return sequences


def read_taxonomy(tsv_file):
  df = pd.read_csv(tsv_file, sep='\t')
  df.columns = ['Feature ID', 'Taxon']
  df['Feature ID'] = df['Feature ID'].astype(str)
  df['Taxon'] = df['Taxon'].str.replace('tax=', '')
  return df


def create_combined_dataframe(fasta_file, tsv_file):
  sequences = read_fasta(fasta_file)
  taxonomy_df = read_taxonomy(tsv_file)
  sequences_df = pd.DataFrame(list(sequences.items()), columns=['Feature ID', 'sequence'])
  sequences_df['Feature ID'] = sequences_df['Feature ID'].astype(str)
  combined_df = pd.merge(sequences_df, taxonomy_df, on='Feature ID', how='inner')
  return combined_df


def analyze_unique_taxa(df):
  taxa_counts = df['Taxon'].value_counts()
  unique_taxa = taxa_counts[taxa_counts == 1]
  unique_taxa_df = df[df['Taxon'].isin(unique_taxa.index)]
  non_unique_taxa_df = df[~df['Taxon'].isin(unique_taxa.index)]
  return unique_taxa_df['Feature ID'], len(unique_taxa)


def create_labeled_datasets(df, suffix, train_size=0.7, test_size=0.1, dev_size=0.2, seed=42):
  unique_taxa = df['Taxon'].unique()
  taxa_to_label = {taxon: i for i, taxon in enumerate(sorted(unique_taxa))}

  full_df = pd.DataFrame()
  full_df['sequence'] = df['sequence']
  full_df['label'] = df['Taxon'].map(taxa_to_label)

  full_df.to_csv(f'complete_{suffix}.csv', index=False)

  assert abs(train_size + test_size + dev_size - 1.0) < 1e-10

  shuffled_df = full_df.sample(frac=1, random_state=seed)
  total_size = len(shuffled_df)
  train_end = int(total_size * train_size)
  test_end = train_end + int(total_size * test_size)

  train_df = shuffled_df[:train_end]
  test_df = shuffled_df[train_end:test_end]
  dev_df = shuffled_df[test_end:]

  train_df.to_csv(f'train_{suffix}.csv', index=False)
  test_df.to_csv(f'test_{suffix}.csv', index=False)
  dev_df.to_csv(f'dev_{suffix}.csv', index=False)

  with open(f'taxa_labels_mapping_{suffix}.txt', 'w') as f:
    for taxon, label in taxa_to_label.items():
      f.write(f"{label}, {taxon}\n")

  print("\nDataset sizes:")
  print(f"Full dataset: {len(full_df)}")
  print(f"Train set: {len(train_df)} ({len(train_df)/total_size:.1%})")
  print(f"Test set: {len(test_df)} ({len(test_df)/total_size:.1%})")
  print(f"Dev set: {len(dev_df)} ({len(dev_df)/total_size:.1%})")


def format_sequence(sequence, line_length=80):
  return '\n'.join(sequence[i:i+line_length] for i in range(0, len(sequence), line_length))


def clean_dataset(fasta_file, tsv_file, ids_file, output_fasta, output_tsv):
  with open(ids_file, 'r') as f:
    ids_to_remove = set(line.strip() for line in f)

  sequences = read_fasta(fasta_file)
  filtered_sequences = {
    seq_id: seq
    for seq_id, seq in sequences.items()
    if seq_id.split()[0] not in ids_to_remove
  }
  with open(output_fasta, 'w') as f:
    for seq_id, sequence in filtered_sequences.items():
      f.write(f">{seq_id}\n")
      f.write(format_sequence(sequence) + '\n')

  df = pd.read_csv(tsv_file, sep='\t')
  df.columns = ['Feature ID', 'Taxon']
  filtered_df = df[~df['Feature ID'].astype(str).isin(ids_to_remove)]

  filtered_df.to_csv(output_tsv, sep='\t', index=False)

  print(f"\nResults for {fasta_file}:")
  print(f"Original sequences: {len(sequences)}")
  print(f"Filtered sequences: {len(filtered_sequences)}")
  print(f"Removed sequences: {len(sequences) - len(filtered_sequences)}")
  
 
df_ncbi = create_combined_dataframe('dna-sequences-ncbi.fasta', 'taxonomy-ncbi.tsv')
df_bold = create_combined_dataframe('dna-sequences-bold.fasta', 'taxonomy-bold.tsv')

uniques_ncbi, n_uniques_ncbi = analyze_unique_taxa(df_ncbi)
uniques_bold, n_uniques_bold = analyze_unique_taxa(df_bold)
print(f"NCBI dataset has {n_uniques_ncbi} sigle used labels for {df_ncbi.shape[0]} rows.")
print(f"Bold dataset has {n_uniques_bold} sigle used labels for {df_bold.shape[0]} rows.")
uniques_bold.to_csv('bold_ids_to_remove.txt', header=False, index=False)
uniques_ncbi.to_csv('ncbi_ids_to_remove.txt', header=False, index=False)


create_labeled_datasets(df_ncbi, 'ncbi')
create_labeled_datasets(df_bold, 'bold')

clean_dataset(
  fasta_file='dna-sequences-bold.fasta',
  tsv_file='taxonomy-bold.tsv',
  ids_file='bold_ids_to_remove.txt',
  output_fasta='dna-sequences-bold-clean.fasta',
  output_tsv='taxonomy-bold-clean.tsv'
)
clean_dataset(
  fasta_file='dna-sequences-ncbi.fasta',
  tsv_file='taxonomy-ncbi.tsv',
  ids_file='ncbi_ids_to_remove.txt',
  output_fasta='dna-sequences-ncbi-clean.fasta',
  output_tsv='taxonomy-ncbi-clean.tsv'
)
