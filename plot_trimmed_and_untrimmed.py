import pandas as pd
import matplotlib.pyplot as plt

# import data

# n.bp: is the number of bases
# N: is the number of sequences with that number of bases

df_untrimmed = pd.read_csv('seqlength_untrimmed_hist.txt', sep='\s+', header=None, names=['N', 'n.bp'])
df_Ntrimmed = pd.read_csv('seqlength_Ntrimmed_hist.txt', sep='\s+', header=None, names=['N', 'n.bp'])

bins = list(range(0, 1100, 100)) + [float('inf')]
labels = [f'{i}-{i+100}' for i in range(0, 1000, 100)] + ['1075+']

df_untrimmed['n.bp class'] = pd.cut(df_untrimmed['n.bp'], bins=bins, labels=labels, right=False)
df_Ntrimmed['n.bp class'] = pd.cut(df_Ntrimmed['n.bp'], bins=bins, labels=labels, right=False)

untrimmed_hist = df_untrimmed.groupby('n.bp class', observed=False)['N'].sum()
Ntrimmed_hist = df_Ntrimmed.groupby('n.bp class', observed=False)['N'].sum()


# make plots
plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
plt.bar(untrimmed_hist.index, untrimmed_hist.values, color='skyblue', edgecolor='black')
plt.title('Histogram for untrimmed Seq')
plt.xlabel('bp')
plt.ylabel('Seq Count')
plt.xticks(rotation=45)

plt.subplot(1, 2, 2)
plt.bar(Ntrimmed_hist.index, Ntrimmed_hist.values, color='green', edgecolor='black')
plt.title('Histogram for Ntrimmed Seq')
plt.xlabel('bp')
plt.ylabel('Seq Count')
plt.xticks(rotation=45)

plt.tight_layout()
plt.savefig('Ntrimmed_untrimmed_comparison.png', dpi=300)
