# KARMA/VLT Observation Catalog Preparation

This project contains scripts to prepare a catalog for a KMOS/VLT observation. It takes properly formatted lists of sources as input, computes the optimal pointing position to include as many sources as possible in the field of view, then queries the SIMBAD database for reference and guide stars with the proper properties and positions.

## Steps

0. Create DataFrame from published data tables.
1. Load the input DataFrame containing the list of sources.
2. Compute the optimal pointing position (barycenter) to include as many sources as possible in the field of view.
3. Remove entries that are further away from the barycenter than a specified distance.
4. Find empty positions in the field of view that are further apart than a minimum distance from all other positions to use as "Sky Positions".
5. Query the SIMBAD database for guide stars and reference stars.
6. Add the guide stars and reference stars to the DataFrame.
7. Assign target priorities to the sources.
8. Save the final catalog to a file.

## Files

- `make_m33_catalog.py`: Script to prepare the catalog for M33.
- `make_m83_catalog.py`: Script to prepare the catalog for M83.
- `make_ngc7793_catalog.py`: Script to prepare the catalog for NGC7793.
- `make_karma_catalogs.py`: Script to prepare the final KARMA catalog.
- `kmos_functions.py`: Contains utility functions used by the scripts.
- `run_all_scripts.sh`: Bash script to sequentially run all the Python scripts.

## Usage

### Running the Python Scripts

To run the individual Python scripts, use the following commands:

```bash
python make_m33_catalog.py
python make_m83_catalog.py
python make_ngc7793_catalog.py
python make_karma_catalogs.py
```

### Running all Scripts Sequentially

To run all the scripts sequentially, use the provided bash script:

```bash
chmod +x run_all_scripts.sh
./run_all_scripts.sh
```

## Dependencies

Python 3.x
pandas
numpy
astroquery
astropy

```bash
pip install pandas numpy astroquery astropy
```

## License

This project is licensed under the MIT License. See the LICENSE file for details.

## Acknowledgments

This project uses data from the SIMBAD database, operated at CDS, Strasbourg, France.
