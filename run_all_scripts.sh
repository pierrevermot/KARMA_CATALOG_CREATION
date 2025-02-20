#!/bin/bash

# Activating conda environment
echo "Activating conda environment..."
conda init
conda activate /home/pierre/anaconda3/

# Run make_m33_catalog.py
echo "Running make_m33_catalog.py..."
python /home/pierre/Documents/KMOS/make_m33_catalog.py

# Run make_m83_catalog.py
echo "Running make_m83_catalog.py..."
python /home/pierre/Documents/KMOS/make_m83_catalog.py

# Run make_ngc7793_catalog.py
echo "Running make_ngc7793_catalog.py..."
python /home/pierre/Documents/KMOS/make_ngc7793_catalog.py

# Run make_karma_catalogs.py
echo "Running make_karma_catalogs.py..."
python /home/pierre/Documents/KMOS/make_karma_catalogs.py

echo "All scripts have been run successfully."
