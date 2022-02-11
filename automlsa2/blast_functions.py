#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import print_function
import subprocess
import numpy as np  # type: ignore
import pandas as pd  # type: ignore
import json
import os
import logging
import shlex
import csv
import sys
import concurrent.futures as cf
from rich.progress import track
from typing import List, Dict, Any, DefaultDict
from .helper_functions import (
    json_writer,
    checkpoint_reached,
    checkpoint_tracker,
    remove_intermediates,
    end_program,
    SUFFIXES,
    worker_init,
)
from collections import defaultdict

# from .pbar import track_wide

SINGLE_COPY_ESTIMATE = 0.90
PERCENT_CUTOFF = 50.0


def make_blast_database(makeblastdb: str, fasta: str) -> None:
    """Generates blast DB if necessary per fasta

    input  - fasta file
    return - null
    """
    exist: List[bool] = [os.path.exists(fasta + x) for x in SUFFIXES]
    if not all(exist):
        with open(os.path.join('fasta', 'makeblastdb.log'), 'a') as textlog:
            subprocess.run(
                [makeblastdb, '-dbtype', 'nucl', '-in', fasta],
                stdout=textlog,
                stderr=subprocess.STDOUT,
            )


def reader(fn: str) -> List[List[str]]:
    """
    Read tsv file into list
    input  - filename to read
    return - list of results
    """
    dat: List[List[str]] = []
    with open(fn, 'r') as tsv:
        tsvreader = csv.reader(
            [row for row in tsv if not row.startswith('#')], delimiter='\t'
        )
        for row in tsvreader:
            dat.append(row)
    return dat


def generate_blast_list(
    rundir: str,
    exe: str,
    queries: List[str],
    targets: List[str],
    evalue: float,
    threads: int,
    checkpoint: bool,
) -> List[str]:
    """
    Generates list of BLAST searches that need to be run

    input  - blast type, exes, queries, and targets
    output - list of blast output names
    """
    logger = logging.getLogger(__name__)
    cmds: List[List[str]] = []
    # blastout = defaultdict(list)
    blastout: List[str] = []
    blastdir: str = os.path.join(rundir, 'blast')
    outfmt: str = (
        '7 qseqid sseqid saccver pident qlen length bitscore qcovhsp stitle sseq'
    )
    base_cmd: List[str] = [exe, '-evalue', str(evalue), '-outfmt', outfmt]
    cmd: List[str] = []
    if not os.path.exists(blastdir):
        os.mkdir(blastdir)
    for target in targets:
        target_base: str = os.path.basename(target)
        for query in queries:
            query_base: str = os.path.splitext(os.path.basename(query))[0]
            outname: str = '{}_vs_{}.tab'.format(query_base, target_base)
            outpath: str = os.path.join(blastdir, outname)
            blastout.append(os.path.join(outpath))
            cmd = base_cmd + ['-db', target, '-query', query, '-out', outpath]
            if not os.path.exists(outpath) or os.path.getsize(outpath) == 0:
                cmds.append(cmd)

    if cmds:
        if checkpoint:
            blastfile: str = os.path.join(rundir, 'blastcmds.txt')
            with open(blastfile, 'w') as fh:
                for cmd in cmds:
                    fh.write(' '.join([shlex.quote(x) for x in cmd]) + '\n')
            msg = 'BLAST commands written to {}. Exiting.'
            logger.info(msg.format(blastfile))
            checkpoint_reached('prior to BLAST searches')
        else:
            msg = 'Running {} BLAST searches using {} CPUs.'
            ncmds: int = len(cmds)
            logger.info(msg.format(ncmds, threads))
            p: cf.ProcessPoolExecutor = cf.ProcessPoolExecutor(
                threads, worker_init(os.getpid())
            )
            with p:
                for _ in track(
                    p.map(subprocess.run, cmds),
                    'Blast Search...',
                    ncmds,
                ):
                    pass
    else:
        logger.info('No BLAST searches remaining. Moving to parse.')

    open(os.path.join('.autoMLSA', 'checkpoint', 'generate_blast_list'), 'w').close()
    return blastout


def read_blast_results(
    blastfiles: List[str], coverage: int, identity: int, threads: int
) -> pd.DataFrame:
    """
    Parses set of BLAST results. Expects list of files.

    input  - list of BLAST output filenames.
    return - pandas data frame of BLAST output
    """
    logger: logging.Logger = logging.getLogger(__name__)
    headers: List[str] = [
        'qseqid',
        'sseqid',
        'saccver',
        'pident',
        'qlen',
        'length',
        'bitscore',
        'qcovhsp',
        'stitle',
        'sseq',
    ]

    logger.info('Reading BLAST results.')
    dtypes: Dict[str, Any] = {
        'qseqid': 'category',
        'sseqid': 'category',
        'saccver': 'string',
        'pident': np.float,
        'qlen': np.int,
        'length': np.int,
        'bitscore': np.float,
        'qcovhsp': np.int,
        'stitle': 'string',
        'sseq': 'string',
    }

    results_fn: str = os.path.join('.autoMLSA', 'blast_results.tsv')
    if os.path.exists(results_fn):
        logger.debug('Reading from existing BLAST results.')
        results: pd.DataFrame = pd.read_csv(results_fn, dtype=dtypes, sep='\t')
    else:
        # results = pd.DataFrame(columns=headers)
        blastrows: List[List[str]] = []
        nfiles = len(blastfiles)
        p: cf.ProcessPoolExecutor = cf.ProcessPoolExecutor(
            threads, worker_init(os.getpid())
        )
        with p:
            for dat in track(
                p.map(reader, blastfiles),
                'Reading Blast...',
                nfiles,
            ):
                for row in dat:
                    blastrows.append(row)
        results = pd.DataFrame(data=blastrows, columns=headers)
        results = results.astype(dtypes)
        results.query('(pident >= @identity) & (qcovhsp >= @coverage)', inplace=True)
        results.to_csv(results_fn, sep='\t')
    checkpoint_tracker('read_blast_results')
    return results


def print_blast_summary(
    runid: str,
    blastout: pd.DataFrame,
    labels: List[str],
    nallowed: int,
    missing_check: bool,
    checkpoint: bool,
    protect: bool,
) -> pd.DataFrame:
    """
    Generates summary of BLAST results.

    input  - pandas data frame with BLAST results
    return - pandas data frame with BLAST results to keep
    """
    logger: logging.Logger = logging.getLogger(__name__)
    logger.info('Summarizing and filtering BLAST hits.')
    summary: DefaultDict[str, Dict[str, Any]] = defaultdict(dict)
    missing_count: DefaultDict[int, List[str]] = defaultdict(list)
    genome_missing_dict: Dict = {}
    keeps: List[str] = []
    keepsidx: List[str] = []
    genome_map: Dict[str, int] = {k: i for i, k in enumerate(labels)}

    grouped: pd.DataFrame = blastout.groupby(['sseqid', 'qseqid'])
    presence_matrix: pd.DataFrame = grouped.size().unstack('qseqid', fill_value=0)
    label_indexes = list(map(int, presence_matrix.index.values))
    indexes: List[int] = presence_matrix.index.tolist()
    presence_matrix.index = [labels[i] for i in label_indexes]
    presence_matrix.to_csv('presence_matrix.tsv', sep='\t')
    summary['queries']['names'] = presence_matrix.columns.tolist()
    summary['queries']['count'] = len(presence_matrix.columns)
    summary['genomes']['names'] = presence_matrix.index.tolist()  # type: ignore
    summary['genomes']['indexes'] = indexes
    summary['genomes']['count'] = len(presence_matrix.index)

    # Get zero counts for filtering
    # Queries
    removal_candidates: Dict[str, str] = {}
    summary['queries']['missing'] = defaultdict(dict)
    zero_counts_query: pd.DataFrame = presence_matrix.apply(
        lambda x: x[x == 0].index.tolist(), axis=0
    )
    zero_counts_query_dict: Dict[str, List[str]] = zero_counts_query.to_dict()
    for query, genomes in zero_counts_query_dict.items():
        ngenomes = len(genomes)
        percent = ngenomes / summary['genomes']['count'] * 100
        percentf = '{:.2f}'.format(percent)
        if percent > PERCENT_CUTOFF:
            removal_candidates[query] = percentf
        summary['queries']['missing'][query]['genomes'] = genomes
        summary['queries']['missing'][query]['count'] = ngenomes
        summary['queries']['missing'][query]['percent'] = percentf
    summary['queries']['missing'] = dict(summary['queries']['missing'])

    # Genomes
    summary['genomes']['missing'] = defaultdict(dict)
    zero_counts_genome: pd.DataFrame = presence_matrix.apply(
        lambda x: x[x == 0].index.tolist(), axis=1
    )
    zero_counts_genome_dict: Dict[str, List[str]] = zero_counts_genome.to_dict()
    for genome, queries in zero_counts_genome_dict.items():
        nmissing: int = len(queries)
        summary['genomes']['missing'][genome]['queries'] = queries
        summary['genomes']['missing'][genome]['count'] = nmissing
        summary['genomes']['missing'][genome]['percent'] = '{:.2f}'.format(
            nmissing / summary['queries']['count'] * 100
        )
        missing_count[nmissing].append(genome)
        if nmissing == 0:
            keeps.append(genome)
        else:
            genome_missing_dict[genome] = queries
    summary['genomes']['missing'] = dict(summary['genomes']['missing'])

    for genome in genome_missing_dict:
        nmissing = len(genome_missing_dict[genome])
        if nmissing > nallowed:
            msg = 'Genome {} is going to be removed due to missing queries.'
            logger.warning(msg.format(genome))
            msg = 'Increase --allow_missing to {} from {} to keep this genome.'
            logger.warning(msg.format(nmissing, nallowed))
        else:
            msg = 'Keeping genome {} missing {} queries due to --allow_missing'
            logger.debug(msg.format(genome, nmissing))
            keeps.append(genome)
    msg = 'Keeping these genomes ({}):\n         - {}'
    logger.info(msg.format(len(keeps), '\n         - '.join(keeps)))
    keepsidx = [str(genome_map[x]) for x in keeps]
    json_writer(os.path.join('.autoMLSA', 'keepsidx.json'), keepsidx)
    blastout_keeps: pd.DataFrame = (
        blastout.query('sseqid in @keepsidx')
        .sort_values(['bitscore', 'qcovhsp'], ascending=False)
        .drop_duplicates(['qseqid', 'sseqid'])
        .sort_index()
    )

    blastout_keeps.to_csv(os.path.join('.autoMLSA', 'blast_filtered.tsv'), sep='\t')

    # Estimate single copy
    # Currently UNUSED #
    single_copy: pd.DefaultDict = presence_matrix.apply(
        lambda x: ((x.values == 1).sum() / len(presence_matrix)) > SINGLE_COPY_ESTIMATE,
        axis=0,
    )
    json_writer(os.path.join('.autoMLSA', 'single_copy.json'), single_copy.to_dict())
    json_writer('missing_counts.json', missing_count)
    json_writer('blast_summary.json', summary)
    json_writer('missing_by_genome.json', genome_missing_dict)
    if len(keeps) < summary['genomes']['count']:
        logger.warning(
            'Some genomes ({}) are missing one or more genes'.format(
                len(genome_missing_dict)
            )
        )
        logger.warning(
            'Check out the presence_matrix.tsv, and '
            'missing_by_genome.json files to review.'
        )
        if removal_candidates:
            logger.warning('Consider removing these genes from the analysis:')
            for query, percent in removal_candidates.items():
                sys.stderr.write(f'         - {query} -> missing {percent}%\n')
        if not missing_check:
            logger.warning('Use --missing_check to continue.')
            end_program(1)
    expected_filt_fn: str = os.path.join('.autoMLSA', 'expected_filt.json')
    if os.path.exists(expected_filt_fn):
        with open(expected_filt_fn, 'r') as eqf:
            expected_filt: List[str] = json.load(eqf)
    else:
        expected_filt = []

    updated = set(expected_filt) != set(keeps)
    if updated:
        logger.debug('Found new filtered sequences')
        logger.debug('Removing downstream files, if present')
        remove_intermediates(runid, ['genome'], protect)
    if not updated:
        logger.debug('Found no new filtered sequences')
    json_writer(expected_filt_fn, keeps)

    if checkpoint:
        checkpoint_reached('after BLAST result filtering')

    checkpoint_tracker('print_blast_summary')

    return blastout_keeps


def print_fasta_files(blastout: pd.DataFrame, labels: List[str]) -> List[str]:
    """
    Writes unaligned FASTA files as output from BLAST results

    input  - pandas DataFrame with blast output, list of labels
    return - List of unaligned FASTA files
    """
    logger: logging.Logger = logging.getLogger(__name__)
    fastdir: str = 'unaligned'
    labels = [os.path.splitext(x)[0] for x in labels]
    unaligned: List[str] = []
    if not os.path.exists(fastdir):
        os.mkdir(fastdir)
    msg = 'Writing unaligned FASTA sequences, if necessary.'
    logger.info(msg)
    for name, group in track(blastout.groupby('qseqid'), 'Writing FASTA...'):
        fasta: str = os.path.join(fastdir, '{}.fas'.format(name))
        unaligned.append(fasta)
        if os.path.exists(fasta):
            logger.debug('Unaligned file {} already found, skipping.'.format(fasta))
        else:
            if os.path.exists(fasta):
                text = 'Overwriting'
            else:
                text = 'Writing'
            logger.debug('{} {} FASTA file'.format(text, fasta))
            with open(fasta, 'w') as fh:
                for row in group.itertuples():
                    fh.write('>{}\n'.format(labels[int(row.sseqid)]))
                    fh.write('{}\n'.format(row.sseq.replace('-', '')))

    checkpoint_tracker('print_fasta_files')

    return unaligned
