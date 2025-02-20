"""
This script is used to prepare a catalog for a KARMA/VLT observation. It takes properly formatted lists of sources as input,
computes the optimal pointing position to include as many sources as possible in the field of view, then queries the SIMBAD
database for reference and guide stars with the proper properties.

Steps:
1. Load the input DataFrame containing the list of sources.
2. Compute the optimal pointing position (barycenter) to include as many sources as possible in the field of view.
3. Remove entries that are further away from the barycenter than a specified distance (FoV of the arms).
4. Find empty positions in the field of view that are further apart than a minimum distance from all other positions to use as "Sky Positions".
5. Query the SIMBAD database for guide stars and reference stars.
6. Add the guide stars and reference stars to the DataFrame.
7. Assign target priorities to the sources.
8. Save the final catalog to a file.
"""

import pandas as pd
import random
import numpy as np
from datetime import datetime
import kmos_functions as kfunc

R_fov = 0.06
R_ifu = 1e-3
# Calculate the time delta in years between 1st of January 2000 and 1st of August 2025
date_j2000 = datetime(2000, 1, 1)
date_target = datetime(2025, 8, 1)
time_delta = (date_target - date_j2000).days / 365.25

# Convert the DataFrame to the KARMA catalogue format
def convert_to_karma_format(df, output_path, ref_bands=['J', 'H'], ranked_priority=False, rad_ref = 3.6):
    
    random.seed(25)
    
    # Prepare the KARMA catalogue DataFrame
    karma_df = pd.DataFrame()

    # Use the ID column for the name
    karma_df['ID'] = df['ID']

    # Combine RA and Dec into required formats
    karma_df['RA'] = df['RAh'].astype(str) + ":" + df['RAm'].astype(str) + ":" + df['RAs'].astype(str)
    karma_df['Dec'] = df['DEd'].astype(str) + ":" + df['DEm'].astype(str) + ":" + df['DEs'].astype(str)

    # Set the Type column to 'O' for scientific targets
    karma_df['Type'] = 'O'

    # Calculate the barycenter
    barycenter_ra_h, barycenter_ra_m, barycenter_ra_s, barycenter_dec_d, barycenter_dec_m, barycenter_dec_s = kfunc.calculate_optimal_center(df, R_fov, 100000)
    # Add the barycenter to the DataFrame
    barycenter_row = {
        'ID': 'Barycenter',
        'RA': f"{barycenter_ra_h}:{barycenter_ra_m}:{barycenter_ra_s:.2f}",
        'Dec': f"{barycenter_dec_d}:{barycenter_dec_m}:{barycenter_dec_s:.2f}",
        'Type': 'C'
    }
    karma_df = pd.concat([karma_df, pd.DataFrame([barycenter_row])], ignore_index=True)
    karma_df = kfunc.remove_distant_entries(karma_df, R_fov)
    karma_df = kfunc.find_empty_positions(karma_df, R_fov, R_ifu, 1000, 20)
    # Query SIMBAD for guide stars
    barycenter_ra, barycenter_dec = kfunc.hms_to_decimal(barycenter_ra_h, barycenter_ra_m, barycenter_ra_s, barycenter_dec_d, barycenter_dec_m, barycenter_dec_s)
    guide_stars = kfunc.query_simbad_for_guide_stars(barycenter_ra, barycenter_dec, 0, 30, 8, 14, time_delta=time_delta)

    # Add guide stars to the DataFrame
    karma_df = kfunc.add_stars_to_dataframe(karma_df, guide_stars)

    for ref_band in ref_bands:
        
        karma_df_t = karma_df.copy()
        # Query SIMBAD for reference stars
        reference_stars = kfunc.query_simbad_for_reference_stars(barycenter_ra, barycenter_dec, rad_ref, ref_band, time_delta=time_delta)

        # Add reference stars to the DataFrame
        karma_df_t = kfunc.add_stars_to_dataframe(karma_df_t, reference_stars)

        if ranked_priority:
            karma_df_t = kfunc.add_ranked_target_priority(karma_df_t)
        else:
            karma_df_t = kfunc.add_target_priority(karma_df_t)
        karma_df_t = kfunc.add_default_magnitude_and_band(karma_df_t)
        # Save the KARMA catalogue
        karma_df_t.to_csv(output_path+'_'+ref_band+'.cat', sep=' ', index=False, header=False)
        print(f"KARMA catalogue saved to {output_path}")

def main():
    input_path = './M33/merged_sorted_table.csv'
    output_path = './M33/karma_catalogue'

    merged_dataframe = kfunc.load_dataframe(input_path)
    if merged_dataframe is not None:
        convert_to_karma_format(merged_dataframe, output_path, ref_bands=["J", "H"])
        
    input_path = './M83/merged_sorted_table.csv'
    output_path = './M83/karma_catalogue'

    merged_dataframe = kfunc.load_dataframe(input_path)
    if merged_dataframe is not None:
        convert_to_karma_format(merged_dataframe, output_path, ref_bands=['J','H'], ranked_priority=True)

    input_path = './NGC7793/merged_sorted_table.csv'
    output_path = './NGC7793/karma_catalogue'

    merged_dataframe = kfunc.load_dataframe(input_path)
    if merged_dataframe is not None:
        convert_to_karma_format(merged_dataframe, output_path, ref_bands=['J','H'], ranked_priority=True, rad_ref=6.)

if __name__ == "__main__":
    main()

