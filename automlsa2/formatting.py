#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import logging
import os
import json
import glob
from rich.progress import track

# from .pbar import track_wide
from .helper_functions import (
    generate_hash,
    json_writer,
    remove_intermediates,
    end_program,
    sanitize_path,
)
from .blast_functions import make_blast_database
from Bio import SeqIO  # type: ignore
from collections import defaultdict
from typing import List, Dict, Any
from shutil import copy


def convert_fasta(
    rundir: str,
    fastas: List[str],
    labels: List[str],
    makeblastdb: str,
    runid: str,
    protect: bool,
) -> List[str]:
    """Converts fasta file with placeholder label names

    input - unformatted fasta files
    return - Set of formatted fasta files
    also generates BLAST dbs per fasta
    """
    logger = logging.getLogger(__name__)
    fastadir: str = os.path.join(rundir, 'fasta')
    new_fastas: List[str] = []
    renamef: str = os.path.join(rundir, '.autoMLSA', 'rename.json')
    deleted: bool = False
    if os.path.exists(renamef):
        with open(renamef, 'r') as rf:
            rename_info: Dict[Any, Any] = json.load(rf)
    else:
        rename_info = defaultdict(dict)
    if not os.path.exists(fastadir):
        os.mkdir(fastadir)

    for fasta in track(fastas, 'Working...'):
        base = os.path.basename(fasta)
        label = labels.index(base)
        labelfastaf = os.path.join(fastadir, base)
        if base not in rename_info:
            try:
                rename_info[base]['index'] = label
            except KeyError:
                rename_info[base] = {}
                rename_info[base]['index'] = label
        else:
            if check_hash(fasta, base, rename_info[base]['hash']):
                deleted = True
        new_fastas.append(labelfastaf)
        if not os.path.exists(labelfastaf):
            msg = 'Writing renamed fasta for {}'
            logger.debug(msg.format(base))
            seq = ''
            with open(fasta, 'r') as ff, open(labelfastaf, 'w') as lff:
                for rec in SeqIO.parse(ff, 'fasta'):
                    lff.write('>{} {}\n'.format(label, rec.id))
                    lff.write('{}\n'.format(rec.seq))
                    seq += str(rec.seq)
            seqhash = generate_hash(seq)
            rename_info[base]['hash'] = seqhash
    logger.info('Generating and/or validating BLAST DBs.')
    for labelfastaf in track(new_fastas, 'makeblastdb...'):
        make_blast_database(makeblastdb, labelfastaf)
    json_writer(renamef, rename_info)

    expected_fastas_fn: str = os.path.join(rundir, '.autoMLSA', 'expected_fastas.json')
    if os.path.exists(expected_fastas_fn):
        with open(expected_fastas_fn, 'r') as eff:
            expected_fastas: List[str] = json.load(eff)
    else:
        expected_fastas = []

    updated: bool = set(expected_fastas) != set(new_fastas)
    if updated or deleted:
        logger.debug('Found new genome files')
        remove_intermediates(runid, ['genome'], protect)
    else:
        logger.debug('Found no new genome files')
    json_writer(expected_fastas_fn, new_fastas)

    return new_fastas


def check_hash(fasta: str, base: str, seqhash: str) -> bool:
    """
    Checks hash, deletes files related to genome if necessary
    return - value if something is deleted
    """
    logger = logging.getLogger(__name__)
    deleted: bool = False
    with open(fasta, 'r') as ff:
        newseqhash = generate_hash(
            ''.join([str(rec.seq) for rec in SeqIO.parse(ff, 'fasta')])
        )
    if newseqhash != seqhash:
        deleted = True
        msg = 'Old genome file {} sequence has changed. Updating.'
        logger.info(msg.format(base))
        for fa in glob.iglob(os.path.join('fasta', base + '*')):
            os.remove(fa)
        for res in glob.iglob(os.path.join('blast', '*_' + base + '.tab')):
            os.remove(res)
    else:
        logger.debug('Genome file {} sequence has not changed.'.format(base))
    return deleted


def get_queries(
    runid: str, rundir: str, dups: bool, queries: List[str], protect: bool
) -> List[str]:
    """Converts query fasta file(s) into individual fastas for BLAST

    input - FASTA files with queries as individual entries
    return - individual fasta files
    """
    logger = logging.getLogger(__name__)
    querydir: str = os.path.join(rundir, 'queries')
    backups: str = os.path.join(rundir, '.autoMLSA', 'backups')

    new_queries: List[str] = []
    seen: Dict[str, Dict[str, str]] = {}
    hashes: Dict[str, str] = {}
    if not os.path.exists(querydir):
        os.mkdir(querydir)
    for query_file in queries:
        i = 1
        query_base: str = os.path.basename(query_file)
        query_backup: str = os.path.join(backups, query_base)
        if not os.path.exists(query_file):
            if os.path.exists(query_backup):
                msg = 'Using query backup file as {} is missing.'
                logger.warning(msg.format(query_file))
                query_file = query_backup
            else:
                msg = 'Query file {} seems to have been removed and/or lost.'
                logger.critical(msg.format(query_file))
                msg = 'No backup was found in the backups directory.'
                logger.critical(msg)
                msg = (
                    'Unable to continue without the queries. Either replace'
                    ' the file or start again from a new analysis.'
                )
                logger.critical(msg)
                end_program(66)
        with open(query_file, 'r') as qf:
            logger.debug('Reading query file ({})'.format(query_file))
            for rec in SeqIO.parse(qf, 'fasta'):
                safeid: str = sanitize_path(rec.id)
                seqhash: str = generate_hash(str(rec.seq))

                if seqhash in hashes and hashes[seqhash] == safeid:
                    # Skip duplicate header + seq
                    continue
                elif seqhash in hashes and hashes[seqhash] != safeid:
                    # Same seq but different names, throw error
                    mismatchid: str = hashes[seqhash]
                    msg = (
                        'Same sequence found in two query inputs with '
                        'different sequence IDs:'
                    )
                    logger.error(msg)
                    msg = '+++++++++++++++Offending queries+++++++++++++++'
                    logger.error(msg)
                    if query_file == seen[mismatchid]['path']:
                        msg = '{} - TWICE - header ID {} & {}'
                        logger.error(
                            msg.format(
                                os.path.basename(query_file),
                                seen[mismatchid]['unsafeid'],
                                rec.id,
                            )
                        )
                    else:
                        msg = '{} header ID {}'
                        logger.error(msg.format(os.path.basename(query_file), rec.id))
                        logger.error(
                            msg.format(
                                os.path.basename(seen[mismatchid]['path']),
                                seen[mismatchid]['unsafeid'],
                            )
                        )
                    msg = '+++++++++++++++++++++++++++++++++++++++++++++++'
                    logger.error(msg)
                    msg = (
                        'Check your sequences to make sure they aren\'t '
                        'misnamed, fix the problem, and try again.'
                    )
                    logger.error(msg)
                    end_program(65)
                else:
                    # Expected outcome; new id, new seq
                    hashes[seqhash] = safeid

                # Check for same id, allow if --dups flag, otherwise error
                if safeid in seen and not dups:
                    msg = 'Same query name ({}) found in two query inputs:'
                    logger.error(msg.format(safeid))
                    msg = '+++++++++++++++Offending queries+++++++++++++++'
                    logger.error(msg)
                    if query_file == seen[safeid]['path']:
                        msg = '{} - TWICE - sequences {} & {}'
                        logger.error(
                            msg.format(
                                os.path.basename(query_file), seen[safeid]['index'], i
                            )
                        )
                    else:
                        msg = '{} sequence {}'
                        logger.error(msg.format(os.path.basename(query_file), i))
                        logger.error(
                            msg.format(
                                os.path.basename(seen[safeid]['path']),
                                seen[safeid]['index'],
                            )
                        )
                    msg = '+++++++++++++++++++++++++++++++++++++++++++++++'
                    logger.error(msg)
                    msg = 'Remove or rename one of these to continue.'
                    logger.error(msg)
                    msg = (
                        'Alternatively, if this is intended, use the --dups'
                        ' flag to include both copies.'
                    )
                    logger.error(msg)
                    end_program(65)
                else:
                    if safeid in seen and dups:
                        msg = (
                            'Keeping additional query {} with duplicate ID '
                            'from file {}.'
                        )
                        logger.info(msg.format(safeid, os.path.basename(query_file)))
                    seen[safeid] = {}
                    seen[safeid]['path'] = query_file
                    seen[safeid]['index'] = str(i)
                    seen[safeid]['unsafeid'] = rec.id

                fn = os.path.join(querydir, safeid + '_' + seqhash + '.fas')
                new_queries.append(fn)
                if not os.path.exists(fn):
                    msg = 'Writing {} for seq {} in {}'
                    logger.debug(
                        msg.format(
                            os.path.basename(fn), i, os.path.basename(query_file)
                        )
                    )
                    with open(fn, 'w') as fh:
                        fh.write('>{}\n'.format(safeid))
                        fh.write('{}\n'.format(rec.seq))
                i += 1
        if query_file != query_backup:
            copy(query_file, backups)

    # Check if expected equals new_queries
    expected_queries_fn = os.path.join(rundir, '.autoMLSA', 'expected_queries.json')
    if os.path.exists(expected_queries_fn):
        logger.debug('Reading query file {}'.format(expected_queries_fn))
        with open(expected_queries_fn, 'r') as eqf:
            expected_queries: List[str] = json.load(eqf)
    else:
        expected_queries = []

    updated: bool = set(expected_queries) != set(new_queries)
    if updated:
        if all(query in new_queries for query in expected_queries):
            logger.debug('Found new query sequences')
            remove_intermediates(runid, ['query'], protect)
        else:
            logger.debug('Query sequences have been removed')
            remove_intermediates(runid, ['query', 'genome'], protect)
    else:
        logger.debug('Found no new query sequences')
    json_writer(expected_queries_fn, new_queries)

    if not new_queries:
        msg = 'No query sequences found. Check your query file and try again.'
        logger.error(msg)
        end_program(-1)

    open(os.path.join('.autoMLSA', 'checkpoint', 'get_queries'), 'w').close()

    return new_queries
