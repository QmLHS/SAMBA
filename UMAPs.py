import pandas as pd
import numpy as np
import umap
from itertools import product

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

def parse_taxonomy_to_columns(taxonomy_string):
  standard_ranks = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
  rank_names = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
  taxonomy_part = taxonomy_string.replace('tax=', '')
  ranks = taxonomy_part.split(';')
  ranks = [r.strip() for r in ranks]

  result = {rank_name: "" for rank_name in rank_names}

  for i, rank_prefix in enumerate(standard_ranks):
    if i < len(ranks) and ranks[i]:
      rank_value = ranks[i]
      if not rank_value.startswith(rank_prefix):
        rank_value = rank_prefix + rank_value

      name_part = rank_value[len(rank_prefix):]
      if name_part:  # Only add non-empty values
        result[rank_names[i]] = name_part
      else:
        result[rank_names[i]] = np.nan
  return result

def read_taxonomy(tsv_file):
  df = pd.read_csv(tsv_file, sep='\t')
  df.columns = ['Feature ID', 'Taxon']
  df['Feature ID'] = df['Feature ID'].astype(str)
  taxonomy_data = df['Taxon'].apply(parse_taxonomy_to_columns)
  taxonomy_df = pd.DataFrame(taxonomy_data.tolist())
  result_df = pd.concat([df[['Feature ID']], taxonomy_df], axis=1)
  return result_df

def generate_kmers(sequence, k):
  return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

def count_kmers(sequence, k_max=5):
  kmer_counts = {}
  for k in range(1, k_max + 1):
    kmers = generate_kmers(sequence, k)
    for kmer in kmers:
      if kmer in kmer_counts:
        kmer_counts[kmer] += 1
      else:
        kmer_counts[kmer] = 1
  return kmer_counts

def create_kmer_matrix(sequences, k_max=5):
    all_kmers = set()
    for seq in sequences:
        for k in range(1, k_max + 1):
            kmers = generate_kmers(seq, k)
            all_kmers.update(kmers)
    kmer_list = sorted(all_kmers)
    kmer_matrix = pd.DataFrame(0, index=range(len(sequences)), columns=kmer_list)
    for i, seq in enumerate(sequences):
        kmer_counts = count_kmers(seq, k_max)
        kmer_matrix.loc[i] = kmer_counts
    return kmer_matrix


seqs_bold_train = read_fasta('dna-sequences-bold-clean-train.fasta')
seqs_bold_test = read_fasta('dna-sequences-bold-clean-test.fasta')
seqs_ncbi_train = read_fasta('dna-sequences-ncbi-clean-train.fasta')
seqs_ncbi_test = read_fasta('dna-sequences-ncbi-clean-test.fasta')
taxa_bold_train = read_taxonomy('taxonomy-bold-clean-train.tsv')
taxa_bold_test = read_taxonomy('taxonomy-bold-clean-test.tsv')
taxa_ncbi_train = read_taxonomy('taxonomy-ncbi-clean-train.tsv')
taxa_ncbi_test = read_taxonomy('taxonomy-ncbi-clean-test.tsv')

print('All file loaded.')

concatenated_seqs_bold = pd.concat([seqs_bold_train, seqs_bold_test], axis=0, ignore_index=True)
concatenated_seqs_ncbi = pd.concat([seqs_ncbi_train, seqs_ncbi_test], axis=0, ignore_index=True)
concatenated_taxa_bold = pd.concat([taxa_bold_train, taxa_bold_test], axis=0, ignore_index=True)
concatenated_taxa_ncbi = pd.concat([taxa_ncbi_train, taxa_ncbi_test], axis=0, ignore_index=True)

print('All concatenation done.')

del seqs_bold_train, seqs_bold_test, seqs_ncbi_train, seqs_ncbi_test
del taxa_bold_train, taxa_bold_test, taxa_ncbi_train, taxa_ncbi_test



bold = pd.merge(concatenated_seqs_bold, concatenated_taxa_bold, on='Feature ID', how='inner')
ncbi = pd.merge(concatenated_seqs_ncbi, concatenated_taxa_ncbi, on='Feature ID', how='inner')

print('All merging done.')

del concatenated_seqs_bold, concatenated_seqs_ncbi
del concatenated_taxa_bold, concatenated_taxa_ncbi


# Additional length check
bold['seq_len'] = bold['Sequence'].str.len()
ncbi['seq_len'] = ncbi['Sequence'].str.len()

bold_filtd = bold[bold['seq_len'] < 200]
ncbi_filtd = ncbi[ncbi['seq_len'] < 200]

bold_filtd = bold_filtd.drop(columns=['seq_len'])
ncbi_filtd = ncbi_filtd.drop(columns=['seq_len'])

print('Filtered for length < 190.')

del bold, ncbi



bold_kmer_matrix = create_kmer_matrix(bold_filtd['Sequence'])
bold_filtd = pd.concat([bold_filtd, bold_kmer_matrix], axis=1)

del bold_kmer_matrix

bold_filtd = bold_filtd.drop(columns=['Sequence'])
bold_filtd.to_csv('bold_filtd_UMAP.csv')

print('Bold file ready.')



ncbi_kmer_matrix = create_kmer_matrix(ncbi_filtd['Sequence'])
ncbi_filtd = pd.concat([ncbi_filtd, ncbi_kmer_matrix], axis=1)

del ncbi_kmer_matrix

ncbi_filtd = ncbi_filtd.drop(columns=['Sequence'])
ncbi_filtd.to_csv('ncbi_filtd_UMAP.csv')

print('NCBI file ready.')
