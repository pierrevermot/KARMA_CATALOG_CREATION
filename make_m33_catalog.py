import pandas as pd
import numpy as np

# Define the column specifications based on byte ranges for table1
colspecs_table1 = [(0, 8),   # ID
                   (9, 11),  # RAh
                   (12, 14), # RAm
                   (15, 20), # RAs
                   (21, 23), # DEd
                   (24, 26), # DEm
                   (27, 31), # DEs
                   (32, 35), # D
                   (36, 39), # rho
                   (40, 47), # L10
                   (48, 56), # LL14
                   (57, 64), # XMM
                   (65, 68), # X-ray
                   (69, 73), # Prev
                   (74, 77), # New
                   (78, 81)] # Confirm

# Define column names for table1
colnames_table1 = ["ID", "RAh", "RAm", "RAs", "DEd", "DEm", "DEs", "D", "rho", "L10", "LL14", "XMM", "X-ray", "Prev", "New", "Confirm"]

# Define the column specifications based on byte ranges for table2
colspecs_table2 = [(0, 8),  # ID
                   (9, 13), # FHa
                   (14, 15), # f_Hb
                   (15, 19), # Hb
                   (20, 21), # f_[OIII]
                   (21, 25), # [OIII]
                   (26, 27), # f_[OI]
                   (27, 31), # [OI]
                   (32, 35), # Ha
                   (36, 37), # f_[NII]
                   (37, 41), # [NII]
                   (42, 43), # f_[SII]6717
                   (43, 47), # [SII]6717
                   (48, 49), # f_[SII]6731
                   (49, 53), # [SII]6731
                   (54, 55), # f_[SII]/Ha
                   (55, 60), # [SII]/Ha
                   (61, 64)] # FWHM

# Define column names for table2
colnames_table2 = ["ID", "FHa", "f_Hb", "Hb", "f_[OIII]", "[OIII]", "f_[OI]", "[OI]", "Ha", "f_[NII]", "[NII]", "f_[SII]6717", "[SII]6717", "f_[SII]6731", "[SII]6731", "f_[SII]/Ha", "[SII]/Ha", "FWHM"]

# Read the file into a pandas DataFrame
def read_table_to_dataframe(file_path, colspecs, colnames):
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()

        # Find the starting line of the table data (after description)
        start_idx = 0
        for i, line in enumerate(lines):
            if line.startswith('L10-'):  # Table data starts here
                start_idx = i
                break

        # Read only the relevant lines into a DataFrame
        df = pd.read_fwf(file_path, colspecs=colspecs, names=colnames, dtype=str, skiprows=start_idx)
        return df
    except Exception as e:
        print(f"Error reading file: {e}")
        return None

# Function to filter objects with the most references
def filter_objects_with_references(df):
    try:
        # Check for required columns containing valid identifiers
        filtered_df = df[(df["L10"] != "-") &
                         (df["LL14"] != "-") &
                         (df["XMM"] != "-") &
                         (df["X-ray"] == "yes") &
                         (df["Prev"] != "-") &
                         (df["Confirm"] == "yes")]
        return filtered_df
    except Exception as e:
        print(f"Error filtering objects: {e}")
        return None

# Merge two tables on the source identifier
def merge_tables(filtered_df, table2_df):
    try:
        # Perform inner merge to get objects present in both DataFrames
        merged_df = pd.merge(filtered_df, table2_df, on="ID", how="inner")
        return merged_df
    except Exception as e:
        print(f"Error merging tables: {e}")
        return None

# Sort by FHa and save the merged DataFrame
def sort_and_save_dataframe(df, output_path):
    try:
        # Convert FHa to numeric for sorting
        df["FHa"] = pd.to_numeric(df["FHa"], errors="coerce")
        # Sort the DataFrame by FHa in descending order
        df_sorted = df.sort_values(by="FHa", ascending=False)
        # Save to file
        df_sorted.to_csv(output_path, index=False)
        print(f"DataFrame saved to {output_path}")
    except Exception as e:
        print(f"Error sorting and saving DataFrame: {e}")

# Example usage
file_path_table1 = './M33/table1.txt'
file_path_table2 = './M33/table2.txt'
output_path = './M33/merged_sorted_table.csv'

dataframe_table1 = read_table_to_dataframe(file_path_table1, colspecs_table1, colnames_table1)
dataframe_table2 = read_table_to_dataframe(file_path_table2, colspecs_table2, colnames_table2)

if dataframe_table1 is not None and dataframe_table2 is not None:
    filtered_dataframe_table1 = filter_objects_with_references(dataframe_table1)
    if filtered_dataframe_table1 is not None:
        # Merge filtered table1 and table2 on "ID"
        merged_dataframe = merge_tables(filtered_dataframe_table1, dataframe_table2)
        if merged_dataframe is not None:
            # Sort and save the merged DataFrame
            sort_and_save_dataframe(merged_dataframe, output_path)