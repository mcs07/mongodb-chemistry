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

import datetime
import time

import numpy as np
import psycopg2


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


#
# def profile_postgres():
#     conn = psycopg2.connect("dbname=chembl user=matt")
#     cur = conn.cursor()
#
#     for threshold in [0.95, 0.9, 0.85, 0.8, 0.75, 0.7]:
#         report = []
#         times = []
#         counts = []
#         cur.execute("set rdkit.tanimoto_threshold=%s;", (threshold,))
#         for chembl_id in chembl_ids:
#             # ARGH! The CHEMBL ID vs. molregno thing is a nightmare
#             cur.execute("select entity_id from chembl_id_lookup where chembl_id = %s", (chembl_id,))
#             molregno = cur.fetchone()[0]
#             cur.execute("select m from rdk.mols where molregno = %s", (molregno,))
#             smiles = cur.fetchone()[0]
#             cur.execute("select mfp3 from rdk.fps where molregno = %s", (molregno,))
#             qfp = cur.fetchone()[0]
#             report.append('Query: {} - {}'.format(chembl_id, smiles))
#             print(chembl_id)
#             start = time.time()
#             cur.execute("select molregno from rdk.fps where mfp3%%%s", (qfp,))
#             #cur.execute("select molregno from rdk.fps where mfp%%morganbv_fp(%s)", (smiles,))
#             results = cur.fetchall()
#             end = time.time()
#             counts.append(len(results))
#             report.append('Results ({})'.format(len(results)))
#             print(report[-1])
#             for r in results:
#                 cur.execute("select chembl_id from chembl_id_lookup where entity_type = 'COMPOUND' and entity_id = %s", (r,))
#                 report.append(cur.fetchone()[0])
#                 print(report[-1])
#             times.append(end - start)
#
#         # Produce a report of the results
#         report.append('Counts: {}'.format(' '.join(str(c) for c in counts)))
#         report.append('Mean: {}'.format(np.mean(counts)))
#         report.append('Median: {}'.format(np.median(counts)))
#         report.append('Times: {}'.format(' '.join(str(t) for t in times)))
#         report.append('Mean: {}'.format(np.mean(times)))
#         report.append('Median: {}'.format(np.median(times)))
#         report.append('95th Percentile: {}'.format(np.percentile(times, 95)))
#         report = '\n'.join(report)
#         print(report)
#         with open('benchmarks/{}-{}-{}-{}.txt'.format('unfolded', 'postgres', threshold, datetime.datetime.utcnow()), 'w') as f:
#             f.write(report)
#     cur.close()
#     conn.close()


if __name__ == '__main__':
    build_postgres()
