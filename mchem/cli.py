# -*- coding: utf-8 -*-
"""
mchem.cli
~~~~~~~~~

Command line interface for the mongodb-chemistry project.

:copyright: Copyright 2014 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""
import logging
import random

import click
import pymongo
from rdkit import Chem

from . import __version__
from . import build, fps, similarity, screening, plot


MONGODB_URI = 'mongodb://localhost:27017'
MONGODB_DB = 'mchem'
MONGODB_COLL = 'mols'


@click.group()
@click.option('--uri', '-u', default=MONGODB_URI, envvar='MCHEM_MONGODB_URI', help='MongoDB URI (default: mongodb://localhost:27017).')
@click.option('--db', '-d', default=MONGODB_DB, envvar='MCHEM_MONGODB_DB', help='MongoDB database name (default: mchem).')
@click.option('--verbose', '-v', is_flag=True, help='Verbose debug logging.')
@click.version_option(__version__, '--version', '-V')
@click.help_option('--help', '-h')
@click.pass_context
def cli(ctx, uri, db, verbose):
    """mongodb-chemistry command line interface."""
    click.echo('mchem v%s' % __version__)
    click.echo('Connecting to %s/%s' % (uri, db))
    logging.basicConfig(level=logging.DEBUG if verbose else logging.INFO, format='%(levelname)s: %(message)s')
    ctx.obj = pymongo.MongoClient(uri)[db]


@cli.command()
@click.option('--collection', '-c', default=MONGODB_COLL, envvar='MCHEM_MONGODB_COLL', help='MongoDB collection name (default: mols).')
@click.option('--idfield', '-i', help='SDF property field to use for molecule ID.')
@click.argument('sdf', type=click.Path(exists=True, resolve_path=True), nargs=-1)
@click.pass_obj
def load(db, sdf, collection, idfield):
    """Load molecules from SDF file(s) into a MongoDB collection."""
    click.echo('mchem.load')
    if db[collection].count() > 0:
        click.confirm('Add to existing %s collection?' % collection, abort=True)
    build.load_sdfs(sdf, db[collection], idfield)


@cli.command()
@click.argument('collection')
@click.pass_obj
def drop(db, collection):
    """Drop a MongoDB collection."""
    click.echo('mchem.drop')
    if click.confirm('Are you sure you want to drop the collection %s?' % collection.name):
        db.drop_collection(collection)
        click.echo('Dropped the collection %s!' % collection.name)


def get_fingerprinter(name, radius, length=None):
    fingerprinter = {
        'morgan': fps.MorganFingerprinter(radius=radius, length=length)
        # Add other fingerprinters here in future
    }[name]
    return fingerprinter


@cli.command()
@click.option('--collection', '-c', default=MONGODB_COLL, envvar='MCHEM_MONGODB_COLL', help='Molecule collection (default: mols).')
@click.option('--fp', default='morgan', type=click.Choice(['morgan']), help='Fingerprint type (default: morgan).')
@click.option('--radius', default=2, help='Radius for morgan fingerprint (default: 2).')
@click.option('--length', type=click.IntRange(0), help='Fold length for fingerprint (default: no folding).')
@click.pass_obj
def addfp(db, collection, fp, radius, length):
    """Generate a fingerprint for every molecule in a MongoDB collection."""
    click.echo('mchem.addfp')
    fingerprinter = get_fingerprinter(fp, radius, length)
    fp_collection = db['%s.%s' % (collection, fingerprinter.name)]
    if fp_collection.count() > 0:
        click.confirm('Add %s fingerprints to existing %s collection?' % (collection, fp_collection.name), abort=True)
    fps.generate(db[collection], fp_collection, fingerprinter)


@cli.command()
@click.option('--collection', '-c', default=MONGODB_COLL, envvar='MCHEM_MONGODB_COLL', help='Molecule collection (default: mols).')
@click.option('--fp', default='morgan', type=click.Choice(['morgan']), help='Fingerprint type (default: morgan).')
@click.option('--radius', default=2, help='Radius for morgan fingerprint (default: 2).')
@click.option('--length', type=click.IntRange(0), help='Fold length for fingerprint (default: no folding).')
@click.pass_obj
def countfp(db, collection, fp, radius, length):
    """Pre-calculate count of each fingerprint bit for a specific fingerprint in a specific molecule collection."""
    click.echo('mchem.countfp')
    fingerprinter = get_fingerprinter(fp, radius, length)
    fp_collection = db['%s.%s' % (collection, fingerprinter.name)]
    count_collection = db['%s.counts' % fp_collection.name]
    fps.count(fp_collection, count_collection)


class FloatRange(click.ParamType):
    """A parameter that works similar to :data:`click.FLOAT` but restricts the value to fit into a range."""
    name = 'float range'

    def __init__(self, min=None, max=None, clamp=False):
        self.min = min
        self.max = max
        self.clamp = clamp

    def convert(self, value, param, ctx):
        try:
            rv = float(value)
        except ValueError:
            self.fail('%s is not a valid floating point value' % value, param, ctx)
        if self.clamp:
            if self.min is not None and rv < self.min:
                return self.min
            if self.max is not None and rv > self.max:
                return self.max
        if self.min is not None and rv < self.min or self.max is not None and rv > self.max:
            if self.min is None:
                self.fail('%s is bigger than the maximum valid value %s.' % (rv, self.max), param, ctx)
            elif self.max is None:
                self.fail('%s is smaller than the minimum valid value %s.' % (rv, self.min), param, ctx)
            else:
                self.fail('%s is not in the valid range of %s to %s.' % (rv, self.min, self.max), param, ctx)
        return rv

    def __repr__(self):
        return 'FloatRange(%r, %r)' % (self.min, self.max)


@cli.command()
@click.argument('smiles', required=True)
@click.option('--collection', '-c', default=MONGODB_COLL, envvar='MCHEM_MONGODB_COLL', help='Molecule collection (default: mols).')
@click.option('--threshold', default=0.8, type=FloatRange(0, 1), help='Similarity search threshold (default 0.8).')
@click.option('--fp', default='morgan', type=click.Choice(['morgan']), help='Fingerprint type (default: morgan).')
@click.option('--radius', default=2, help='Radius for morgan fingerprint (default: 2).')
@click.option('--length', type=click.IntRange(0), help='Fold length for fingerprint (default: no folding).')
@click.pass_obj
def similar(db, smiles, collection, threshold, fp, radius, length):
    """Perform a similarity search."""
    click.echo('mchem.similar')
    fingerprinter = get_fingerprinter(fp, radius, length)
    fp_collection = db['%s.%s' % (collection, fingerprinter.name)]
    count_collection = db['%s.counts' % fp_collection.name]
    mol = Chem.MolFromSmiles(smiles.encode())
    results = similarity.similarity_search(mol, fingerprinter, fp_collection, threshold, count_collection)
    for result in results:
        print result['_id']


@cli.command()
@click.option('--collection', '-c', default=MONGODB_COLL, envvar='MCHEM_MONGODB_COLL', help='Molecule collection (default: mols).')
@click.option('--size', type=click.IntRange(0), default=1000, help='Sample size.')
@click.option('--seed', default=201405291515, help='Random seed.')
@click.argument('output', type=click.File('w'), required=True)
@click.pass_obj
def sample(db, collection, size, seed, output):
    """Choose a random sample of molecules from a collection."""
    click.echo('mchem.sample')
    # TODO: I think this will break down for larger collections.
    # Instead assign random number field to molecules when creating?
    ids = [m['_id'] for m in db[collection].find().sort('_id')]
    random.seed(seed)
    rands = random.sample(ids, size)
    for rand in rands:
        output.write('%s\n' % rand)


def load_sample_mols(db, collection, f):
    ids = f.read().strip().split('\n')
    return [Chem.Mol(mol['rdmol']) for mol in db[collection].find({'_id': {'$in': ids}})]


@cli.command()
@click.argument('test', type=click.Choice(['screening', 'fingerprint']), required=True)
@click.option('--collection', '-c', default=MONGODB_COLL, envvar='MCHEM_MONGODB_COLL', help='Molecule collection (default: mols).')
@click.option('--sample', type=click.File('r'), help='File containing sample ids.')
@click.option('--fp', default='morgan', type=click.Choice(['morgan']), help='Fingerprint type (default: morgan).')
@click.option('--radius', default=2, help='Radius for morgan fingerprint (default: 2).')
@click.option('--length', type=click.IntRange(0), help='Fold length for fingerprint (default: no folding).')
@click.pass_obj
def analyse(db, test, collection, sample, fp, radius, length):
    """Analyse different screening methods."""
    click.echo('mchem.test')
    mols = load_sample_mols(db, collection, sample)
    fingerprinter = get_fingerprinter(fp, radius, length)
    fp_collection = db['%s.%s' % (collection, fingerprinter.name)]
    count_collection = db['%s.counts' % fp_collection.name]
    result_collection = db['%s.test' % collection]
    thresholds = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]

    for threshold in thresholds:
        screening.test_screening(mols, fingerprinter, fp_collection, result_collection, threshold, count_collection)
        # Test alternative screening methods
        if test == 'screening':
            screening.test_screening(mols, fingerprinter, fp_collection, result_collection, threshold, reqbits=False)
            screening.test_screening(mols, fingerprinter, fp_collection, result_collection, threshold, counts=False)
            screening.test_screening(mols, fingerprinter, fp_collection, result_collection, threshold, count_collection, counts=False)




@cli.command()
@click.argument('test', type=click.Choice(['constraints', 'rhist', 'folding', 'radius']), required=True)
@click.pass_obj
def results(db, test):
    """Run various profiling tests. Lots of hardcoded stuff here that needs fixing."""
    click.echo('mchem.results')
    if test == 'constraints':
        plot.plot_constraints(db)
    elif test == 'rhist':
        plot.plot_radius_hist(db)
    elif test == 'folding':
        plot.plot_folding(db)
    elif test == 'radius':
        plot.plot_radius(db)


@cli.command()
@click.argument('test', type=click.Choice(['build', 'profile']), required=True)
@click.pass_obj
def pg(db, test):
    """Build and benchmark PostgreSQL."""
    from . import postgres
    click.echo('mchem.postgres')
    if test == 'build':
        postgres.build_postgres()
    elif test == 'profile':
        postgres.profile_postgres(db)

