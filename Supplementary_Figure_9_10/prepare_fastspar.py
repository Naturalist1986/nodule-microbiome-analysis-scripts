import pandas as pd

df = pd.read_csv('/mnt/c/Temp/combined_abundance_long_renormalized.tsv', sep='\t')

wide = df.pivot_table(index='Bin_Prefixed', columns='Sample', values='Relative_Abundance', fill_value=0)
wide.index.name = '#OTU ID'
wide = wide.reset_index()

wide.to_csv('/mnt/c/Temp/fastspar_input.tsv', sep='\t', index=False)
print(f"Output shape: {wide.shape[0]} rows x {wide.shape[1]} cols")
print(wide.iloc[:3, :5])
