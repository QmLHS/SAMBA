# Plot the results for the homopolymers:

import pandas as pd
import matplotlib.pyplot as plt

df_h = pd.read_csv('homopolymer_frequencies.txt', sep='\s+', header=None, names=['homopolymer', 'N'])

print(df_h)

# Remove lines with N=0 because that example homopolymer is not present in our Seqs
df_h = df_h.loc[df_h['N']>0]

print(df_h)

df_h['length'] = df_h['homopolymer'].str.len()

grouped = df_h.groupby(['length', 'homopolymer'])['N'].sum().unstack(fill_value=0)

grouped.plot(kind='bar', stacked=True, figsize=(12, 8))
plt.title('Stacked Bar Plot of Homopolymers by Length')
plt.xlabel('Length of Homopolymer')
plt.ylabel('Frequency (N)')
plt.legend(title='Homopolymers', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig('homopolymers_visual.png', dpi=300)
