# MongoDB Chemistry

Some ideas for running chemical similarity searches in MongoDB.

See:

http://blog.matt-swain.com/post/87093745652/chemical-similarity-search-in-mongodb

Also:

http://blog.rguha.net/?p=1300  
http://datablend.be/?p=254  
http://datablend.be/?p=256  
http://datablend.be/?p=265


## Installation

    brew install rdkit
    pip install -r requirements.txt
    python setup.py install
    

## Usage

A simple example:

    mchem load mymols.sdf
    mchem addfp
    mchem countfp
    mchem similar O=C(Oc1ccccc1C(=O)O)C
