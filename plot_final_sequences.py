import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('anml_seq_lengths.freq.table', sep='\s+', header=None, names=['N', 'bp'])

plt.figure(figsize=(10, 6))
plt.bar(df['bp'], df['N'], color='blue', label='Final sequences')
plt.title("Fineal sequences before last processing")
plt.xlabel("bp")
plt.ylabel("Count")
plt.grid(True)
plt.legend()
plt.savefig('final_Seqs_before_preproc.png', dpi=300)
