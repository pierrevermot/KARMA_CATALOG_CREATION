import pandas as pd

# Read table to DataFrame
def read_table_to_dataframe(file_path, colspecs, colnames):
    try:
        df = pd.read_fwf(file_path, colspecs=colspecs, names=colnames)
        return df
    except Exception as e:
        print(f"Error reading table {file_path}: {e}")
        return None

# Merge DataFrames
def merge_dataframes(df1, df2):
    try:
        merged_df = pd.concat([df1, df2], ignore_index=True)
        return merged_df
    except Exception as e:
        print(f"Error merging tables: {e}")
        return None

# Filter the DataFrame
def filter_dataframe(df):
    try:
        filtered_df = df[(df["F(Ha)"] > 50) | (df["OName"].str.strip() != '') | (df["X"] == 'y')]
        return filtered_df
    except Exception as e:
        print(f"Error filtering DataFrame: {e}")
        return None

# Sort by F(Ha) and save the merged DataFrame
def sort_and_save_dataframe(df, output_path):
    try:
        # Convert F(Ha) to numeric for sorting
        df["F(Ha)"] = pd.to_numeric(df["F(Ha)"], errors="coerce")
        # Sort the DataFrame by F(Ha) in descending order
        df_sorted = df.sort_values(by="F(Ha)", ascending=False)
        # Save to file
        df_sorted.to_csv(output_path, index=False)
        print(f"DataFrame saved to {output_path}")
    except Exception as e:
        print(f"Error sorting and saving DataFrame: {e}")

# Column specifications and names for table2 and table3
colspecs_table2 = [(0, 3), (4, 5), (6, 8), (9, 11), (12, 17), (18, 21), (22, 24), (25, 29), (30, 34), (35, 38), (39, 44), (45, 49), (50, 55), (56, 60), (61, 65), (66, 67), (68, 90)]
colnames_table2 = ["ID", "f_Seq", "RAh", "RAm", "RAs", "DEd", "DEm", "DEs", "Rad", "Diam", "F(Ha)", "F([SII])", "F([OIII])", "[SII]/Ha", "[OIII]/Ha", "X", "OName"]

colspecs_table3 = colspecs_table2
colnames_table3 = colnames_table2

# Example usage
file_path_table2 = './M83/table2.dat'
file_path_table3 = './M83/table3.dat'
output_path = './M83/merged_sorted_table.csv'

dataframe_table2 = read_table_to_dataframe(file_path_table2, colspecs_table2, colnames_table2)
dataframe_table3 = read_table_to_dataframe(file_path_table3, colspecs_table3, colnames_table3)

if dataframe_table2 is not None and dataframe_table3 is not None:
    merged_dataframe = merge_dataframes(dataframe_table2, dataframe_table3)
    if merged_dataframe is not None:
        filtered_dataframe = filter_dataframe(merged_dataframe)
        if filtered_dataframe is not None:
            sort_and_save_dataframe(filtered_dataframe, output_path)