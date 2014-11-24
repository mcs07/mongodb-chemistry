#!/bin/sh

# TODO: Download PubChem SDF files

# Load PubChem SDF files into a MongoDB collection called pubchem
mchem load pubchem/*.sdf.gz --collection pubchem --id PUBCHEM_COMPOUND_CID

# Add fingerprints
mchem addfp --collection pubchem --fp morgan --radius 2
mchem addfp --collection pubchem --fp morgan --radius 3
mchem addfp --collection pubchem --fp morgan --radius 4
mchem addfp --collection pubchem --fp morgan --radius 2 --length 512
mchem addfp --collection pubchem --fp morgan --radius 2 --length 1024
mchem addfp --collection pubchem --fp morgan --radius 2 --length 2048

# Count fingerprints
mchem countfp --collection pubchem --fp morgan --radius 2
mchem countfp --collection pubchem --fp morgan --radius 3
mchem countfp --collection pubchem --fp morgan --radius 4
mchem countfp --collection pubchem --fp morgan --radius 2 --length 512
mchem countfp --collection pubchem --fp morgan --radius 2 --length 1024
mchem countfp --collection pubchem --fp morgan --radius 2 --length 2048
