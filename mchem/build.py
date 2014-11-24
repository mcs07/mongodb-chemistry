# -*- coding: utf-8 -*-
"""
mchem.build
~~~~~~~~~~~

Functions for building the database.

:copyright: Copyright 2014 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import gzip
import logging
import os
import urllib2

from bson import Binary
import pymongo
from rdkit import Chem

from . import errors


log = logging.getLogger(__name__)


def download_chembl():
    """Download chembl_19.sdf.gz to data subdirectory."""
    url = 'ftp://ftp.ebi.ac.uk//pub/databases/chembl/ChEMBLdb/releases/chembl_19/chembl_19.sdf.gz'
    dest = os.path.realpath(os.path.join(os.path.dirname(__file__), '../data/chembl_19.sdf.gz'))
    log.info('Downloading %s' % url)
    if not os.path.isfile(dest):
        r = urllib2.urlopen(url)
        with open(dest, 'wb') as f:
            f.write(r.read())
    else:
        log.info('Already downloaded: %s' % dest)


def load_sdfs(paths, collection, idfield=None):
    """Insert molecules from list of SDF file paths into database collection.

    :param string paths: Paths to SDF files.
    :param collection: MongoDB collection.
    :param string idfield: SDF property field to use for molecule ID.
    """
    log.info('Inserting molecules into MongoDB %s' % collection.name)
    for path in paths:
        log.info('Loading %s' % path)
        sdf = gzip.open(path) if path[-3:] == '.gz' else open(path)
        load_sdf(sdf, collection, idfield)


def load_sdf(sdf, collection, idfield):
    """Insert molecules from SDF file into database collection.

    :param file sdf: SDF file object.
    :param collection: MongoDB collection.
    :param string idfield: SDF property field to use for molecule ID.
    """

    success, fail, skip = 0, 0, 0
    for rdmol in Chem.ForwardSDMolSupplier(sdf):
        if rdmol is None:
            log.debug('Failed to read molecule')
            fail += 1
            continue
        mol = {
            'smiles': Chem.MolToSmiles(rdmol, isomericSmiles=True),
            'rdmol': Binary(rdmol.ToBinary()),
        }
        if idfield:
            mol['_id'] = rdmol.GetProp(idfield.encode())
        try:
            collection.insert(mol)
            log.debug('Inserted %s' % mol['_id'])
            success += 1
        except pymongo.errors.DuplicateKeyError:
            log.debug('Skipped %s: Already exists' % mol['_id'])
            skip += 1
    log.info('%s successes, %s failures, %s skipped' % (success, fail, skip))


