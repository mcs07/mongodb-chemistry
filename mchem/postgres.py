# -*- coding: utf-8 -*-
"""
mchem.postgres
~~~~~~~~~~~~~~

Functions to build and benchmark PostgreSQL database for comparison.

:copyright: Copyright 2014 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import logging
import time

import numpy as np
import psycopg2


log = logging.getLogger(__name__)


def build_postgres():
    """Build PostgreSQL database."""
    # Start by creating the database and loading the chembl dump via the command line:
    # createdb chembl
    # psql chembl < chembl_19.pgdump.sql
    conn = psycopg2.connect("host=127.0.0.1 dbname=chembl user=matt")
    cur = conn.cursor()
    cur.execute('create extension if not exists rdkit;')
    cur.execute('create schema rdk;')
    cur.execute('select * into rdk.mols from (select molregno,mol_from_ctab(molfile::cstring) m  from compound_structures) tmp where m is not null;')
    cur.execute('create index molidx on rdk.mols using gist(m);')
    cur.execute('alter table rdk.mols add primary key (molregno);')
    cur.execute('select molregno, m into rdk.fps from rdk.mols;')
    cur.execute('alter table rdk.fps add column m2l512 bfp;')
    cur.execute('alter table rdk.fps add column m2l2048 bfp;')
    cur.execute('alter table rdk.fps add column m2 sfp;')
    cur.execute('alter table rdk.fps add column m3 sfp;')
    cur.execute('update rdk.fps set m2 = morgan_fp(m);')
    cur.execute('update rdk.fps set m3 = morgan_fp(m, 3);')
    cur.execute('set rdkit.morgan_fp_size=2048;')
    cur.execute('update rdk.fps set m2l2048 = morganbv_fp(m);')
    cur.execute('set rdkit.morgan_fp_size=512;')
    cur.execute('update rdk.fps set m2l512 = morganbv_fp(m);')
    cur.execute('alter table rdk.fps drop column m;')
    cur.execute('create index fps_m2_idx on rdk.fps using gist(m2);')
    cur.execute('create index fps_m3_idx on rdk.fps using gist(m3);')
    cur.execute('create index fps_m2l2048_idx on rdk.fps using gist(m2l2048);')
    cur.execute('create index fps_m2l512_idx on rdk.fps using gist(m2l512);')
    cur.execute('alter table rdk.fps add primary key (molregno);')
    cur.close()
    conn.close()


def profile_postgres(db):
    conn = psycopg2.connect("host=127.0.0.1 dbname=chembl user=matt")
    cur = conn.cursor()
    with open('../data/sample_chembl_1000.txt') as f:
        chembl_ids = f.read().strip().split('\n')
    for fp in ['m2l512', 'm2', 'm3', 'm2l2048']:
        for threshold in [1, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]:
            times = []
            cur.execute("set rdkit.tanimoto_threshold=%s;", (threshold,))
            for i, chembl_id in enumerate(chembl_ids[:100]):
                log.debug('Query molecule %s of %s: %s' % (i+1, len(chembl_ids), chembl_id))
                # ARGH! The CHEMBL ID vs. molregno thing is a nightmare
                cur.execute("select entity_id from chembl_id_lookup where chembl_id = %s", (chembl_id,))
                molregno = cur.fetchone()[0]
                cur.execute("select m from rdk.mols where molregno = %s", (molregno,))
                smiles = cur.fetchone()[0]
                cur.execute("select %s from rdk.fps where molregno = %s", (fp, molregno,))
                qfp = cur.fetchone()[0]
                log.debug(chembl_id)
                start = time.time()
                cur.execute("select molregno from rdk.fps where %s%%%s", (fp, qfp,))
                #cur.execute("select molregno from rdk.fps where %s%%morganbv_fp(%s)", (fp, smiles,))  # using smiles
                results = cur.fetchall()
                end = time.time()
                times.append(end - start)
            # Save results
            result = {
                'median_time': np.median(times),
                'mean_time': np.mean(times),
                'fp': fp,
                'threshold': threshold,
                'sample': 'chembl_1000'
            }
            db.profile.postgres.insert(result)
    cur.close()
    conn.close()
