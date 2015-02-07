#!/bin/sh

sudo apt-get update
sudo apt-get install -y git build-essential libboost-all-dev flex bison g++ cmake make
sudo apt-get install -y postgresql postgresql-server-dev-all postgresql-doc postgresql-contrib
sudo apt-get install -y python-dev python-pip python-numpy python-scipy python-matplotlib python-psycopg2

echo "Creating swapfile"
# This is for if we are on low RAM server - need plenty of memory to compile RDKit
sudo fallocate -l 4G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile

echo "Installing RDKit"
curl -L -o Release_2014_09_2.tar.gz https://github.com/rdkit/rdkit/archive/Release_2014_09_2.tar.gz
tar -zxvf Release_2014_09_2.tar.gz
cd rdkit-Release_2014_09_2
mkdir build
cd build
export PREFIX=/usr/local
export RDBASE=$HOME/rdkit-Release_2014_09_2
export LD_LIBRARY_PATH=$RDBASE/lib:$PREFIX/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$RDBASE:$PYTHONPATH
cmake -DRDK_INSTALL_INTREE=OFF -DCMAKE_INSTALL_PREFIX=$PREFIX ..
make install
#ctest
cd ../Code/PgSQL/rdkit
make
sudo make install
#make installcheck

echo "Removing swapfile"
sudo swapoff -a -v
sudo rm /swapfile

echo "Setting up PostgreSQL database"
sudo -u postgres createuser -dsr root
createdb chembl_19
curl -o chembl_19_postgresql.tar.gz ftp://ftp.ebi.ac.uk//pub/databases/chembl/ChEMBLdb/releases/chembl_19/chembl_19_postgresql.tar.gz
tar -zxvf chembl_19_postgresql.tar.gz
psql chembl_19 < chembl_19_postgresql/chembl_19.pgdump.sql


# Backup postgresql.conf
NOW=$(date +"%Y-%m-%d-%H-%M-%S")
cp "/etc/postgresql/9.3/main/postgresql.conf" "/etc/postgresql/9.3/main/postgresql.conf-$NOW"

# Update shared_buffers and work_mem in postgresql.conf
sed -i "s/^#*shared_buffers \=.*/shared_buffers = 256MB/" "/etc/postgresql/9.3/main/postgresql.conf"
sed -i "s/^#*work_mem \=.*/work_mem = 128MB/" "/etc/postgresql/9.3/main/postgresql.conf"

# Restart PostgreSQL
sudo service postgresql restart

# Set up fingerprints etc.
pgchem -v --db chembl_19 load

# Get sample mols
curl -o sample_chembl_1000.txt https://raw.githubusercontent.com/mcs07/mongodb-chemistry/master/data/sample_chembl_1000.txt

# Run benchmarks
pgchem -v --db chembl_19 profile --fp m2l2048 --sample sample_chembl_1000.txt --threshold 1.0
pgchem -v --db chembl_19 profile --fp m2l2048 --sample sample_chembl_1000.txt --threshold 0.95
pgchem -v --db chembl_19 profile --fp m2l2048 --sample sample_chembl_1000.txt --threshold 0.9
pgchem -v --db chembl_19 profile --fp m2l2048 --sample sample_chembl_1000.txt --threshold 0.85
pgchem -v --db chembl_19 profile --fp m2l2048 --sample sample_chembl_1000.txt --threshold 0.8
pgchem -v --db chembl_19 profile --fp m2l2048 --sample sample_chembl_1000.txt --threshold 0.75
pgchem -v --db chembl_19 profile --fp m2l2048 --sample sample_chembl_1000.txt --threshold 0.7
pgchem -v --db chembl_19 profile --fp m2l2048 --sample sample_chembl_1000.txt --threshold 0.65
pgchem -v --db chembl_19 profile --fp m2l2048 --sample sample_chembl_1000.txt --threshold 0.6
pgchem -v --db chembl_19 profile --fp m2l2048 --sample sample_chembl_1000.txt --threshold 0.55
pgchem -v --db chembl_19 profile --fp m2l2048 --sample sample_chembl_1000.txt --threshold 0.5
