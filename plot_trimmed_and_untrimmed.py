import pandas as pd
import matplotlib.pyplot as plt

# import data
df_Ntrimmed = pd.read_csv('seqlength_Ntrimmed_hist.txt', sep='\s+', header=None, names=['Frequency', 'ID'])
df_untrimmed = pd.read_csv('seqlength_untrimmed_hist.txt', sep='\s+', header=None, names=['Frequency', 'ID'])


# make plots
xmin, xmax, ymin, ymax = 0, 2000, 0, max(df_Ntrimmed['Frequency'].max(), df_untrimmed['Frequency'].max()) + 10

plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
plt.bar(df_Ntrimmed['ID'], df_Ntrimmed['Frequency'], color='blue', alpha=0.7)
plt.title('Histogram for Ntrimmed Seq')
plt.xlabel('ID')
plt.ylabel('Seq Count')
plt.xlim([xmin, xmax])
plt.ylim([ymin, ymax])

plt.subplot(1, 2, 2)
plt.bar(df_untrimmed['ID'], df_untrimmed['Frequency'], color='blue', alpha=0.7)
plt.title('Histogram for untrimmed Seq')
plt.xlabel('ID')
plt.ylabel('Seq Count')
plt.xlim([xmin, xmax])
plt.ylim([ymin, ymax])

plt.tight_layout()
plt.savefig('Ntrimmed_untrimmed_comparison.png', dpi=300)
