#!/bin/sh

# Download chembl_19.sdf.gz
curl -o chembl_19.sdf.gz ftp://ftp.ebi.ac.uk//pub/databases/chembl/ChEMBLdb/releases/chembl_19/chembl_19.sdf.gz

# Load chembl_19.sdf.gz into a MongoDB collection called chembl
mchem load chembl_19.sdf.gz --collection chembl --idfield chembl_id

# Add fingerprints
mchem addfp --collection chembl --fp morgan --radius 2
mchem addfp --collection chembl --fp morgan --radius 3
mchem addfp --collection chembl --fp morgan --radius 4
mchem addfp --collection chembl --fp morgan --radius 2 --length 512
mchem addfp --collection chembl --fp morgan --radius 2 --length 1024
mchem addfp --collection chembl --fp morgan --radius 2 --length 2048

# Count fingerprints
mchem countfp --collection chembl --fp morgan --radius 2
mchem countfp --collection chembl --fp morgan --radius 3
mchem countfp --collection chembl --fp morgan --radius 4
mchem countfp --collection chembl --fp morgan --radius 2 --length 512
mchem countfp --collection chembl --fp morgan --radius 2 --length 1024
mchem countfp --collection chembl --fp morgan --radius 2 --length 2048

# Create random sample of 1000 molecules
mchem sample --collection chembl --size 1000 --seed 201405291515 data/sample_chembl_1000.txt

# Test the different screening methods
mchem analyse screening --collection chembl --sample data/sample_chembl_1000.txt

# Test different fingerprint folding
mchem analyse fingerprint --collection chembl --sample data/sample_chembl_1000.txt --length 2048
mchem analyse fingerprint --collection chembl --sample data/sample_chembl_1000.txt --length 1024
mchem analyse fingerprint --collection chembl --sample data/sample_chembl_1000.txt --length 512

# Test different fingerprint radius
mchem analyse fingerprint --collection chembl --sample data/sample_chembl_1000.txt --radius 3
mchem analyse fingerprint --collection chembl --sample data/sample_chembl_1000.txt --radius 4
