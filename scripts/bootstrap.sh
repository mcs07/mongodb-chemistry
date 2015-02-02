#!/bin/sh

sudo apt-get -y update

# Install MongoDB
sudo apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv 7F0CEB10
echo "deb http://downloads-distro.mongodb.org/repo/ubuntu-upstart dist 10gen" | tee -a /etc/apt/sources.list.d/mongodb.list
sudo apt-get -y update
sudo apt-get -y install mongodb-org

sudo apt-get -y install git

# Install Python stuff
sudo apt-get -y install build-essential python-dev python-pip python-numpy python-scipy python-matplotlib
sudo pip install -U pip
sudo pip install -U setuptools
sudo pip install pymongo

# Install RDKit
sudo apt-get -y install python-rdkit librdkit1 rdkit-data

# Install mongodb-chemistry
pip install git+https://github.com/mcs07/mongodb-chemistry.git
