# MongoDB Chemistry

Chemical similarity search implementation in MongoDB, with performance analysis.

See [this blog post](http://blog.matt-swain.com/post/87093745652/chemical-similarity-search-in-mongodb) for more information.

See also:

- [Presentation slides](https://speakerdeck.com/mcs07/chemical-structure-handling-in-mongodb)
- [Blog post by Rajarshi Guha](http://blog.rguha.net/?p=1261)
- [Blog post by Davy Suvee](http://datablend.be/?p=254) ([part 2](http://datablend.be/?p=256), [part 3](http://datablend.be/?p=265))


## Installation

### Mac OS X

    brew install rdkit mongodb
    pip install git+https://github.com/mcs07/mongodb-chemistry.git
    
### Ubuntu

See [scripts/bootstrap.sh](https://github.com/mcs07/mongodb-chemistry/blob/master/scripts/bootstrap.sh).

## Usage

A simple example:

    mchem load mymols.sdf
    mchem addfp
    mchem countfp
    mchem similar O=C(Oc1ccccc1C(=O)O)C

See [scripts/chembl.sh](https://github.com/mcs07/mongodb-chemistry/blob/master/scripts/chembl.sh) for a more detailed example.
