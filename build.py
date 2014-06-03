#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Build database for chemical similarity searches in MongoDB."""

from bson import Binary
import psycopg2
import pymongo
from rdkit import Chem
from rdkit.Chem import AllChem


db = pymongo.MongoClient().chem


def build_db():
    """Import molecules into database from ChEMBL SDF.

    Each molecule document contains a SMILES, CHEMBL ID, and a pickle of the RDKit Mol object.
    """
    db.mols.drop()
    for rdmol in Chem.ForwardSDMolSupplier('data/chembl_18.sdf'):
        if rdmol is None: continue
        db.molecules.insert({
            'smiles': Chem.MolToSmiles(rdmol, isomericSmiles=True),
            'chembl_id': rdmol.GetProp('chembl_id'),
            'rdmol': Binary(rdmol.ToBinary()),
        })


def add_fps():
    """Generate fingerprints for every molecule in the database.

    Each fingerprint is stored as an array of "on" bits, along with the total count of this array.
    """
    for molecule in db.molecules.find({'mfp3': {'$exists': False}}, timeout=False):
        rdmol = Chem.Mol(molecule['rdmol'])
        print molecule['chembl_id']
        #mfp = list(AllChem.GetMorganFingerprintAsBitVect(rdmol, 2, nBits=2048).GetOnBits())
        #molecule['mfp'] = {'bits': mfp, 'count': len(mfp)}
        # mfp2 = AllChem.GetMorganFingerprint(rdmol, 2).GetNonzeroElements().keys()
        # molecule['mfp2'] = {'bits': mfp2, 'count': len(mfp2)}
        mfp3 = list(AllChem.GetMorganFingerprintAsBitVect(rdmol, 2, nBits=512).GetOnBits())
        molecule['mfp3'] = {'bits': mfp3, 'count': len(mfp3)}
        db.molecules.save(molecule)


def count_fps():
    """Build collection containing total counts of all occurrences of each fingerprint bit.

    This creates a collection that stores counts of the number of times each individual bit occurs in the molecules
    collection. This is used for an efficiency shortcut. The resulting documents have an _id that corresponds to the
    fingerprint bit and a single count field. (e.g. { "_id" : 511, "count" : 148 })
    """
    for fp in ['mfp3']:  # 'mfp', 'mfp2'
        db['{}_counts'.format(fp)].drop()
        counts = {}
        for molecule in db.molecules.find({fp: {'$exists': True}}, timeout=False):
            for bit in molecule[fp]['bits']:
                counts[bit] = counts.get(bit, 0) + 1
        for k, v in counts.items():
            db['{}_counts'.format(fp)].insert({'_id': k, 'count': v})


def ensure_indices():
    """Build index on relevant fields."""
    # db.molecules.ensure_index('mfp.bits')
    # db.molecules.ensure_index('mfp.count')
    # db.molecules.ensure_index('mfp2.bits')
    # db.molecules.ensure_index('mfp2.count')
    db.molecules.ensure_index('mfp3.bits')
    db.molecules.ensure_index('mfp3.count')


def build_postgres():
    """Build PostgreSQL database."""
    # Start by creating the database and loading the chembl dump via the command line:
    # createdb chembl
    # psql chembl < chembl_18.pgdump.sql
    conn = psycopg2.connect("dbname=chembl user=matt")
    cur = conn.cursor()
    cur.execute('create extension if not exists rdkit;')
    cur.execute('create schema rdk;')
    cur.execute('select * into rdk.mols from (select molregno,mol_from_ctab(molfile::cstring) m  from compound_structures) tmp where m is not null;')
    cur.execute('create index molidx on rdk.mols using gist(m);')
    cur.execute('alter table rdk.mols add primary key (molregno);')
    cur.execute('select molregno, m into rdk.fps from rdk.mols;')
    cur.execute('alter table rdk.fps add column mfp bfp;')
    cur.execute('alter table rdk.fps add column mfp2 sfp;')
    cur.execute('alter table rdk.fps add column mfp3 bfp;')
    cur.execute('set rdkit.morgan_fp_size=2048;')
    cur.execute('update rdk.fps set mfp = morganbv_fp(m);')
    cur.execute('update rdk.fps set mfp2 = morgan_fp(m);')
    cur.execute('set rdkit.morgan_fp_size=512;')
    cur.execute('update rdk.fps set mfp3 = morganbv_fp(m);')
    cur.execute('alter table rdk.fps drop column m;')
    cur.execute('create index fps_mfp_idx on rdk.fps using gist(mfp);')
    cur.execute('create index fps_mfp2_idx on rdk.fps using gist(mfp2);')
    cur.execute('create index fps_mfp3_idx on rdk.fps using gist(mfp3);')
    cur.execute('alter table rdk.fps add primary key (molregno);')
    cur.close()
    conn.close()


if __name__ == '__main__':
    #build_db()
    #add_fps()
    #count_fps()
    #ensure_indices()
    build_postgres()
