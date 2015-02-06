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

import click
import numpy as np
import psycopg2
from psycopg2.extensions import AsIs


log = logging.getLogger(__name__)


# Start by creating the database and loading the chembl dump via the command line:
# createdb chembl
# psql chembl < chembl_19.pgdump.sql


@click.group()
@click.option('--db', '-d', default='mchem', envvar='MCHEM_POSTGRES_DB', help='PostgreSQL database name (default: mchem).')
@click.option('--user', '-u', default='root', envvar='MCHEM_POSTGRES_USER', help='PostgreSQL username (default: root).')
@click.option('--password', '-p', default=None, envvar='MCHEM_POSTGRES_PASSWORD', help='PostgreSQL password.')
@click.option('--verbose', '-v', is_flag=True, help='Verbose debug logging.')
@click.help_option('--help', '-h')
@click.pass_context
def cli(ctx, db, user, password, verbose):
    """PostgreSQL command line interface."""
    click.echo('Connecting %s@%s' % (user, db))
    logging.basicConfig(level=logging.DEBUG if verbose else logging.INFO, format='%(levelname)s: %(message)s')
    ctx.obj = psycopg2.connect(database=db, user=user, password=password)


@cli.command()
@click.pass_obj
def load(conn):
    """Build PostgreSQL database."""
    cur = conn.cursor()
    cur.execute('create extension if not exists rdkit;')
    cur.execute('create schema rdk;')
    cur.execute('drop table if exists biotherapeutics, drug_mechanism, activities, assays, assay_parameters, compound_records, compound_properties, molecule_hierarchy, ligand_eff, predicted_binding_domains, molecule_synonyms, docs, formulations, molecule_atc_classification cascade;')
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
    conn.commit()
    cur.close()
    conn.close()


@cli.command()
@click.option('--sample', type=click.File('r'), help='File containing sample ids.')
@click.option('--fp', default='m2', type=click.Choice(['m2', 'm3', 'm2l2048', 'm2l512', 'm3l2048', 'm3l512']), help='Fingerprint type (default: m2).')
@click.pass_obj
def profile(conn, sample, fp):
    cur = conn.cursor()
    mol_ids = sample.read().strip().split('\n')
    for threshold in [1, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.5]:
        times = []
        cur.execute("set rdkit.tanimoto_threshold=%s;", (threshold,))
        for i, mol_id in enumerate(mol_ids[:100]):
            log.debug('Query molecule %s of %s: %s' % (i+1, len(mol_ids), mol_id))
            # ARGH! The CHEMBL ID vs. molregno thing is a nightmare
            cur.execute("select entity_id from chembl_id_lookup where chembl_id = %s", (mol_id,))
            molregno = cur.fetchone()[0]
            #cur.execute("select m from rdk.mols where molregno = %s", (molregno,))
            #smiles = cur.fetchone()[0]
            cur.execute("select %s from rdk.fps where molregno = %s", (AsIs(fp), molregno,))
            qfp = cur.fetchone()[0]
            log.debug(mol_id)
            start = time.time()
            cur.execute("select molregno from rdk.fps where %s%%%s", (AsIs(fp), qfp,))
            #cur.execute("select molregno from rdk.fps where %s%%morganbv_fp(%s)", (fp, smiles,))  # using smiles
            results = cur.fetchall()
            end = time.time()
            times.append(end - start)
        # Save results
        result = {
            'median_time': np.median(times),
            'mean_time': np.mean(times),
            'fp': fp,
            'threshold': threshold
        }
        log.info(result)
    cur.close()
    conn.close()
