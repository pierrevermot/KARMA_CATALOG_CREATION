import pandas as pd
import numpy as np
import random
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u

# Load the merged DataFrame
def load_dataframe(file_path):
    """
    Loads a DataFrame from a CSV file.

    Parameters:
    file_path (str): The path to the CSV file to be loaded.

    Returns:
    pd.DataFrame: The loaded DataFrame if successful, None otherwise.

    Raises:
    Exception: If there is an error loading the DataFrame, it prints the error message.
    """
    try:
        df = pd.read_csv(file_path)
        return df
    except Exception as e:
        print(f"Error loading DataFrame: {e}")
        return None

# Convert HMS to decimal degrees
def hms_to_decimal(ra_h, ra_m, ra_s, dec_d, dec_m, dec_s):
    """
    Convert right ascension (RA) and declination (Dec) from hours, minutes, and seconds to decimal degrees.

    Parameters:
    ra_h (float): Hours of right ascension.
    ra_m (float): Minutes of right ascension.
    ra_s (float): Seconds of right ascension.
    dec_d (float): Degrees of declination.
    dec_m (float): Minutes of declination.
    dec_s (float): Seconds of declination.

    Returns:
    tuple: A tuple containing:
        - ra (float): Right ascension in decimal degrees.
        - dec (float): Declination in decimal degrees.
    """
    ra = 15 * (ra_h + ra_m / 60 + ra_s / 3600)
    dec = np.sign(dec_d)*(abs(dec_d) + dec_m / 60 + dec_s / 3600)
    return ra, dec

# Convert decimal degrees to HMS
def decimal_to_hms(ra, dec):
    """
    Convert right ascension (RA) and declination (Dec) from decimal degrees to hours, minutes, and seconds.

    Parameters:
    ra (float): Right ascension in decimal degrees.
    dec (float): Declination in decimal degrees.

    Returns:
    tuple: A tuple containing:
        - ra_h (int): Right ascension hours.
        - ra_m (int): Right ascension minutes.
        - ra_s (float): Right ascension seconds.
        - dec_d (int): Declination degrees (with sign).
        - dec_m (int): Declination minutes.
        - dec_s (float): Declination seconds.
    """
    ra_h = int(ra // 15)
    ra_m = int((ra % 15) * 4)
    ra_s = ((ra % 15) * 4 - ra_m) * 60

    sign = int(np.sign(dec))
    dec = abs(dec)
    dec_d = int(dec)
    dec_m = int((dec - dec_d) * 60)
    dec_s = ((dec - dec_d) * 60 - dec_m) * 60

    return ra_h, ra_m, ra_s, sign*dec_d, dec_m, dec_s

# Calculate the barycenter of the positions
def calculate_barycenter(df):
    """
    Calculate the barycenter (average position) of a set of celestial coordinates.

    This function takes a DataFrame containing celestial coordinates in hours, minutes, and seconds
    for right ascension (RA) and degrees, minutes, and seconds for declination (DEC), converts them
    to decimal degrees, and then calculates the average RA and DEC to determine the barycenter.

    Parameters:
    df (pandas.DataFrame): A DataFrame with columns 'RAh', 'RAm', 'RAs', 'DEd', 'DEm', and 'DEs' representing
                           the right ascension and declination in hours, minutes, seconds and degrees, minutes, seconds.

    Returns:
    tuple: A tuple containing the barycenter RA and DEC in hours, minutes, seconds and degrees, minutes, seconds.
    """
    ra_list = []
    dec_list = []

    for index, row in df.iterrows():
        ra, dec = hms_to_decimal(row['RAh'], row['RAm'], row['RAs'], row['DEd'], row['DEm'], row['DEs'])
        ra_list.append(ra)
        dec_list.append(dec)

    barycenter_ra = sum(ra_list) / len(ra_list)
    barycenter_dec = sum(dec_list) / len(dec_list)

    return decimal_to_hms(barycenter_ra, barycenter_dec)

def calculate_optimal_center(df, R, n_it):
    """
    Calculate the optimal center coordinates (RA, Dec) that maximizes the number of targets within a given radius.
    Parameters:
    df (pandas.DataFrame): DataFrame containing the target coordinates with columns 'RAh', 'RAm', 'RAs', 'DEd', 'DEm', 'DEs'.
    R (float): Radius within which to count the number of targets.
    n_it (int): Number of iterations for the random search.
    Returns:
    tuple: Optimal center coordinates in the format (RA, Dec) as strings in HMS and DMS format respectively.
    """
    ra_list = []
    dec_list = []

    for index, row in df.iterrows():
        ra, dec = hms_to_decimal(row['RAh'], row['RAm'], row['RAs'], row['DEd'], row['DEm'], row['DEs'])
        ra_list.append(ra)
        dec_list.append(dec)

    ra_min, ra_max = min(ra_list), max(ra_list)
    dec_min, dec_max = min(dec_list), max(dec_list)

    max_count = 0
    optimal_center = (0, 0)
    
    # First guess
    center_ra = np.random.uniform(ra_min, ra_max, n_it)
    center_dec = np.random.uniform(dec_min, dec_max, n_it)
    ras = np.array(ra_list)
    decs = np.array(dec_list)
    distances = np.sqrt((ras[None, :] - center_ra[:, None])**2 + (decs[None, :] - center_dec[:, None])**2)
    counts = np.sum(distances<=R, 1)
    max_count = np.max(counts)
    arg = np.argmax(counts)
    optimal_center = (center_ra[arg], center_dec[arg])

    # Iterate on increasingly smaller regions
    for K in np.arange(1, 51, 5):
        center_ra = np.random.uniform(optimal_center[0]-R/K, optimal_center[0]+R/K, n_it)
        center_dec = np.random.uniform(optimal_center[1]-R/K, optimal_center[1]+R/K, n_it)
        distances = np.sqrt((ras[None, :] - center_ra[:, None])**2 + (decs[None, :] - center_dec[:, None])**2)
        counts = np.sum(distances<=R, 1)
        max_count = np.max(counts)
        arg = np.argmax(counts)
        optimal_center = (center_ra[arg], center_dec[arg])

    print("Number of targets within FOV : ", max_count)
    return decimal_to_hms(optimal_center[0], optimal_center[1])

def angular_distance(ra1, dec1, ra2, dec2):
    """
    Calculate the angular distance between two points on the celestial sphere given their right ascension and declination.

    Parameters:
    ra1 (float): Right ascension of the first point in degrees.
    dec1 (float): Declination of the first point in degrees.
    ra2 (float): Right ascension of the second point in degrees.
    dec2 (float): Declination of the second point in degrees.

    Returns:
    float: Angular distance between the two points in degrees.
    """
    # Convert degrees to radians
    ra1, dec1, ra2, dec2 = map(np.radians, [ra1, dec1, ra2, dec2])
    # Haversine formula
    d_ra = ra2 - ra1
    d_dec = dec2 - dec1
    a = np.sin(d_dec / 2)**2 + np.cos(dec1) * np.cos(dec2) * np.sin(d_ra / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    return np.degrees(c)
    
def remove_distant_entries(df, R_exclude):
    """
    Remove entries from the DataFrame that are farther than a specified angular distance from the barycenter.
    Parameters:
    df (pandas.DataFrame): The input DataFrame containing celestial object data. 
                           It must have columns 'ID', 'RA', and 'Dec'.
    R_exclude (float): The angular distance threshold. Entries farther than this distance from the barycenter will be removed.
    Returns:
    pandas.DataFrame: A DataFrame with entries farther than the specified angular distance from the barycenter removed.
    Notes:
    - The 'RA' and 'Dec' columns should contain coordinates in the format 'HH:MM:SS' for RA and 'DD:MM:SS' for Dec.
    - The DataFrame must contain an entry with 'ID' equal to 'Barycenter' to determine the barycenter coordinates.
    - If the barycenter entry is not found, the function will print a warning and return the original DataFrame.
    """
    # Look for the barycenter entry in the DataFrame
    barycenter_row = df[df['ID'] == 'Barycenter']
    
    if barycenter_row.empty:
        print("Warning: Barycenter entry not found in the DataFrame.")
        return df

    # Extract the barycenter coordinates
    barycenter_ra_h, barycenter_ra_m, barycenter_ra_s = map(float, barycenter_row['RA'].values[0].split(':'))
    barycenter_dec_d, barycenter_dec_m, barycenter_dec_s = map(float, barycenter_row['Dec'].values[0].split(':'))
    barycenter_ra, barycenter_dec = hms_to_decimal(barycenter_ra_h, barycenter_ra_m, barycenter_ra_s, barycenter_dec_d, barycenter_dec_m, barycenter_dec_s)

    filtered_df = df.copy()
    for index, row in df.iterrows():
        ra_h, ra_m, ra_s = map(float, row['RA'].split(':'))
        dec_d, dec_m, dec_s = map(float, row['Dec'].split(':'))
        ra, dec = hms_to_decimal(ra_h, ra_m, ra_s, dec_d, dec_m, dec_s)
        distance = angular_distance(ra, dec, barycenter_ra, barycenter_dec)
        if distance > R_exclude:
            filtered_df.drop(index, inplace=True)

    return filtered_df

def find_empty_positions(df, R, r_min, N1, N_it, N_lim = 20):
    def find_empty_positions(df, R, r_min, N1, N_it, N_lim=20):
        """
        Find empty Sky positions within a specified radius from the barycenter in a DataFrame.
        This function searches for empty positions around a barycenter within a given radius R,
        ensuring that the positions are at least r_min distance apart from each other and from
        existing positions in the DataFrame. It iterates N_it times to find the maximum number
        of valid positions, with a limit of N_lim positions.
        Parameters:
        df (pd.DataFrame): DataFrame containing the existing positions with columns 'ID', 'RA', and 'Dec'.
        R (float): Maximum radius within which to search for empty positions.
        r_min (float): Minimum distance between any two positions.
        N1 (int): Number of positions to attempt to generate in each iteration.
        N_it (int): Number of iterations to perform.
        N_lim (int, optional): Maximum number of positions to find. Default is 20.
        Returns:
        pd.DataFrame: Updated DataFrame with the new positions added.
        """
    # Look for the barycenter entry in the DataFrame
    barycenter_row = df[df['ID'] == 'Barycenter']
    
    if barycenter_row.empty:
        print("Warning: Barycenter entry not found in the DataFrame.")
        return df

    # Extract the barycenter coordinates
    barycenter_ra_h, barycenter_ra_m, barycenter_ra_s = map(float, barycenter_row['RA'].values[0].split(':'))
    barycenter_dec_d, barycenter_dec_m, barycenter_dec_s = map(float, barycenter_row['Dec'].values[0].split(':'))
    barycenter_ra, barycenter_dec = hms_to_decimal(barycenter_ra_h, barycenter_ra_m, barycenter_ra_s, barycenter_dec_d, barycenter_dec_m, barycenter_dec_s)

    max_positions = []
    
    for _ in range(N_it):
        temp_positions = []
        
        for _ in range(N1):
            if len(temp_positions) < N_lim:
                # Draw a random position within a distance R of the barycenter
                random_angle = random.uniform(0, 2 * np.pi)
                random_radius = random.uniform(0, R)
                random_ra = barycenter_ra + random_radius * np.cos(random_angle)
                random_dec = barycenter_dec + random_radius * np.sin(random_angle)
                
                # Check the distance with all previous positions
                valid_position = True
                for index, row in df.iterrows():
                    ra_h, ra_m, ra_s = map(float, row['RA'].split(':'))
                    dec_d, dec_m, dec_s = map(float, row['Dec'].split(':'))
                    ra, dec = hms_to_decimal(ra_h, ra_m, ra_s, dec_d, dec_m, dec_s)
                    ang_dist = angular_distance(random_ra, random_dec, ra, dec)
                    if ang_dist < r_min:
                        valid_position = False
                        break
                
                if valid_position:
                    for pos in temp_positions:
                        ang_dist = angular_distance(random_ra, random_dec, pos[0], pos[1])
                        if ang_dist < r_min:
                            valid_position = False
                            break
                
                if valid_position:
                    temp_positions.append((random_ra, random_dec))
                    
        if (len(temp_positions) > len(max_positions)) and (len(temp_positions)<=N_lim):
            max_positions = temp_positions
        elif (len(temp_positions) >= len(max_positions)) and (len(temp_positions)>=N_lim):
            min_distance_temp = min(
            list(angular_distance(ra1, dec1, ra2, dec2)
            for i, (ra1, dec1) in enumerate(temp_positions)
            for j, (ra2, dec2) in enumerate(temp_positions)
            if i != j
            ))
            min_distance_max = min(
            list(angular_distance(ra1, dec1, ra2, dec2)
            for i, (ra1, dec1) in enumerate(max_positions)
            for j, (ra2, dec2) in enumerate(max_positions)
            if i != j
            ))
            if min_distance_temp > min_distance_max:
                max_positions = temp_positions
    print(len(max_positions))
    # Add the targets within the max_positions list to the dataframe
    for i, (ra, dec) in enumerate(max_positions):
        ra_h, ra_m, ra_s, dec_d, dec_m, dec_s = decimal_to_hms(ra, dec)
        sky_position_row = {
            'ID': f"SkyPos{i+1}",
            'RA': f"{ra_h}:{ra_m}:{ra_s:.2f}",
            'Dec': f"{dec_d}:{dec_m}:{dec_s:.2f}",
            'Type': 'S'
        }
        df = pd.concat([df, pd.DataFrame([sky_position_row])], ignore_index=True)
    
    return df

# Query SIMBAD for guide stars
def query_simbad_for_guide_stars(ra, dec, radius_min, radius_max, mag_min, mag_max, time_delta=0):
    """
    Queries the SIMBAD database for guide stars within a specified region and magnitude range.
    Parameters:
    ra (float): Right Ascension of the center of the search region in degrees.
    dec (float): Declination of the center of the search region in degrees.
    radius_min (float): Minimum radius of the search region in arcminutes.
    radius_max (float): Maximum radius of the search region in arcminutes.
    mag_min (float): Minimum magnitude of the stars to be included in the search.
    mag_max (float): Maximum magnitude of the stars to be included in the search.
    Returns:
    list: A list of dictionaries, each containing information about a guide star. Each dictionary has the following keys:
        - 'name': The identifier of the star.
        - 'RA': The Right Ascension of the star in HH:MM:SS.ss format.
        - 'Dec': The Declination of the star in DD:MM:SS.ss format.
        - 'Type': The type of the star (always 'G' for guide stars).
        - 'magnitude': The magnitude of the star in the R band.
        - 'wavelength_band': The wavelength band of the magnitude (always 'R').
    """
    custom_simbad = Simbad()
    custom_simbad.add_votable_fields('ra', 'dec', 'R')
    custom_simbad.ROW_LIMIT = -1
    custom_simbad.TIMEOUT = 120
    custom_simbad.add_votable_fields('otype')
    custom_simbad.add_votable_fields('pmra')
    custom_simbad.add_votable_fields('pmdec')
    center = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')
    result = custom_simbad.query_region(center, radius=radius_max*u.arcmin)
    
    guide_stars = []
    for row in result:
        if '*' == row['otype']:
            pmra, pmdec = 0., 0.
            if not np.ma.is_masked(row['pmra']) and not np.ma.is_masked(row['pmdec']):
                pmra, pmdec = 1e-3 * row['pmra'] * time_delta, 1e-3 * row['pmdec'] * time_delta
                pmra, pmdec = pmra / 3600, pmdec / 3600  # Convert from milliarcseconds to degrees
            star_coord = SkyCoord(ra=row['ra']+pmra, dec=row['dec']+pmdec, unit=(u.deg, u.deg), frame='icrs')
            separation = center.separation(star_coord).arcmin
            if '4547' in row['main_id']:
                star_coord = SkyCoord(ra=row['ra']+pmra*0, dec=row['dec']+pmdec*0, unit=(u.deg, u.deg), frame='icrs')
                separation_t = center.separation(star_coord).arcmin
                print(row['main_id'],pmra, pmdec, separation, separation_t)
            if radius_min <= separation <= radius_max and mag_min <= row['R'] <= mag_max:
                ra_h, ra_m, ra_s, dec_d, dec_m, dec_s = decimal_to_hms(row['ra'], row['dec'])
                guide_stars.append({
                    'name': row['main_id'],
                    'RA': f"{ra_h}:{ra_m}:{ra_s:.2f}",
                    'Dec': f"{dec_d}:{dec_m}:{dec_s:.2f}",
                    'Type': 'G',
                    'magnitude': row['R'],
                    'wavelength_band': 'R'
                })
    return guide_stars

# Query SIMBAD for reference stars
def query_simbad_for_reference_stars(ra, dec, radius, band='J', mag_min=8, mag_max=15, time_delta=0):
    """
    Query the SIMBAD database for reference stars within a specified region and magnitude range.
    Parameters:
    ra (float): Right Ascension of the center of the search region in degrees.
    dec (float): Declination of the center of the search region in degrees.
    radius (float): Radius of the search region in arcminutes.
    band (str, optional): Photometric band to use for magnitude filtering (default is 'J').
    mag_min (float, optional): Minimum magnitude for filtering stars (default is 8).
    mag_max (float, optional): Maximum magnitude for filtering stars (default is 15).
    Returns:
    list: A list of dictionaries, each containing information about a reference star:
        - 'name': The identifier of the star.
        - 'RA': Right Ascension in HH:MM:SS.ss format.
        - 'Dec': Declination in DD:MM:SS.ss format.
        - 'Type': Type of the star (always 'R' for reference stars).
        - 'magnitude': Magnitude of the star in the specified band.
        - 'wavelength_band': The photometric band used for the magnitude.
    """
    custom_simbad = Simbad()
    custom_simbad.add_votable_fields(f'{band}')
    custom_simbad.ROW_LIMIT = -1
    custom_simbad.TIMEOUT = 120
    custom_simbad.add_votable_fields('otype')
    custom_simbad.add_votable_fields('pmra')
    custom_simbad.add_votable_fields('pmdec')
    center = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')
    result = custom_simbad.query_region(center, radius=radius*u.arcmin)
    
    reference_stars = []
    for row in result:
        if '*' in row['otype']:
            pmra, pmdec = 0., 0.
            if not np.ma.is_masked(row['pmra']) and not np.ma.is_masked(row['pmdec']):
                pmra, pmdec = 1e-3 * row['pmra'] * time_delta, 1e-3 * row['pmdec'] * time_delta
                pmra, pmdec = pmra / 3600, pmdec / 3600  # Convert from milliarcseconds to degrees
            star_coord = SkyCoord(ra=row['ra']+pmra, dec=row['dec']+pmdec, unit=(u.deg, u.deg), frame='icrs')
            ra_h, ra_m, ra_s, dec_d, dec_m, dec_s = decimal_to_hms(row['ra'], row['dec'])
            separation = center.separation(star_coord).arcmin
            if separation <= radius and mag_min <= row[band.upper()] <= mag_max:
                reference_stars.append({
                    'name': row['main_id'],
                    'RA': f"{ra_h}:{ra_m}:{ra_s:.2f}",
                    'Dec': f"{dec_d}:{dec_m}:{dec_s:.2f}",
                    'Type': 'R',
                    'magnitude': row[f'{band.upper()}'],
                    'wavelength_band': band.upper()
                })
    return reference_stars

def add_target_priority(df):
    """
    Adds a "Target Priority" column to the DataFrame with default priority 3.
    Sets the priority to 1 for rows where the "Type" column is 'O'.
    Parameters:
    df (pandas.DataFrame): The input DataFrame containing a "Type" column.
    Returns:
    pandas.DataFrame: The DataFrame with the added "Target Priority" column.
    """
    # Initialize the "Target Priority" column with default priority 3
    df['Target Priority'] = 3
    
    # Set priority 1 for scientific targets ("Type = 'O'")
    df.loc[df['Type'] == 'O', 'Target Priority'] = 1
    
    return df

def add_ranked_target_priority(df):
    """
    Adds a "Target Priority" column to the DataFrame based on the type of target.
    The function initializes the "Target Priority" column with a default priority of 3.
    It then assigns priorities to scientific targets (Type 'O') and sky positions (Type 'S')
    based on their order in the DataFrame.
    Parameters:
    df (pandas.DataFrame): The input DataFrame containing a 'Type' column.
    Returns:
    pandas.DataFrame: The DataFrame with an added "Target Priority" column.
    """
    # Initialize the "Target Priority" column with default priority 3
    df['Target Priority'] = 3
    
    # Set priority for scientific targets ("Type = 'O'")
    o_targets = df[df['Type'] == 'O']
    df.loc[o_targets.index[:20], 'Target Priority'] = 1
    df.loc[o_targets.index[20:40], 'Target Priority'] = 2
    
    # Set priority for sky positions ("Type = 'S'")
    s_targets = df[df['Type'] == 'S']
    df.loc[s_targets.index[:5], 'Target Priority'] = 1
    df.loc[s_targets.index[5:40], 'Target Priority'] = 2
    
    return df

def add_default_magnitude_and_band(df):
    """
    Add default values for magnitude and wavelength band in a DataFrame.
    This function modifies the input DataFrame by filling missing values in the 
    'magnitude' column with 0 and missing values in the 'wavelength_band' column 
    with 'V'.
    Parameters:
    df (pandas.DataFrame): The input DataFrame containing 'magnitude' and 
                           'wavelength_band' columns.
    Returns:
    pandas.DataFrame: The modified DataFrame with default values filled in.
    """
    # Add default magnitude 0 for targets without a magnitude entry
    df['magnitude'] = df['magnitude'].fillna(0)
    
    # Add default wavelength band 'V' for targets without a wavelength band entry
    df['wavelength_band'] = df['wavelength_band'].fillna('V')
    
    return df

def add_stars_to_dataframe(df, stars):
    """
    Adds a list of stars to a given DataFrame.

    Parameters:
    df (pandas.DataFrame): The DataFrame to which the stars will be added.
    stars (list of dict): A list of dictionaries, where each dictionary represents a star with the following keys:
        - 'name' (str): The name of the star.
        - 'RA' (float): The right ascension of the star.
        - 'Dec' (float): The declination of the star.
        - 'Type' (str): The type of the star.
        - 'magnitude' (float): The magnitude of the star.
        - 'wavelength_band' (str): The wavelength band of the star.

    Returns:
    pandas.DataFrame: The updated DataFrame with the new stars added.
    """
    for star in stars:
        star_row = {
            'ID': star['name'].replace(' ', '-').replace('+', 'p').replace('-', 'm').replace('.', '').replace('[', '').replace(']', ''),
            'RA': star['RA'],
            'Dec': star['Dec'],
            'Type': star['Type'],
            'magnitude': star['magnitude'],
            'wavelength_band': star['wavelength_band']
        }
        df = pd.concat([df, pd.DataFrame([star_row])], ignore_index=True)
    return df