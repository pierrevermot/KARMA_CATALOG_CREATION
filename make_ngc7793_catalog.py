import pandas as pd
import numpy as np

df = pd.read_csv('./NGC7793/table.dat', sep='\t')

# Keep only the first 4 columns and rename them
df = df.iloc[:, :4]
df.columns = ["ID", "RA", "DE", "F(Ha)"]

# Replace "−" with "-" in the dataframe
df = df.replace('−', '-', regex=True)

# Split RA and DE into subcolumns
df[['RAh', 'RAm', 'RAs']] = df['RA'].str.split(':', expand=True)
df[['DEd', 'DEm', 'DEs']] = df['DE'].str.split(':', expand=True)

# Keep only the part of F(Ha) before "±"
df['F(Ha)'] = df['F(Ha)'].str.split('±').str[0].str.strip()

# Convert F(Ha) to numeric
df['F(Ha)'] = pd.to_numeric(df['F(Ha)'], errors='coerce')

# Sort by decreasing value of F(Ha)
df = df.sort_values(by='F(Ha)', ascending=False)

# Drop the original RA and DE columns
df = df.drop(columns=['RA', 'DE'])

output_path = './NGC7793/merged_sorted_table.csv'
df.to_csv(output_path, index=False)
