#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import logging
import os
import glob
import json
import argparse
import subprocess
import pprint
from typing import List, Dict, Any
from .helper_functions import (
    check_if_fasta,
    end_program,
    checkpoint_reached,
    json_writer,
    remove_intermediates,
)
from ._settings import MAFFT, IQTREE


def get_fasta_list(dirpath: str, rundir: str) -> List[str]:
    """
    Finds FASTA files from directory
    """
    logger = logging.getLogger(__name__)
    if not os.path.exists(dirpath):
        fastadir = os.path.join(rundir, 'fasta')
        if os.path.exists(fastadir):
            msg = 'Original genome directory ({}) does not exist.'
            logger.warning(msg.format(dirpath))
            msg = 'Attempting to run from already renamed FASTA files.'
            logger.warning(msg)
            dirpath = fastadir
    fasta_list: List[str] = []
    logger.debug('Checking for FASTA format files in {}'.format(dirpath))
    for fa in glob.iglob(os.path.join(dirpath, '*')):
        # sys.stderr.write('Checking if {} is FASTA.\n'.format(fa))
        if os.path.isdir(fa):
            continue
        if check_if_fasta(fa):
            # sys.stderr.write('{} is FASTA.\n'.format(fa))
            fasta_list.append(os.path.abspath(fa))
        else:
            msg = '{} does not appear to be FASTA file, skipping.'
            logger.debug(msg.format(os.path.basename(fa)))
    if fasta_list:
        msg = 'Found FASTA files in {}'
        logger.info(msg.format(dirpath))
    return fasta_list


def read_config(configfile: str) -> Dict[str, Any]:
    """Finds settings in config file (json format). path/rundir/config.json

    input  - rundir path, logger
    return - dict with configuration settings
    """
    logger = logging.getLogger(__name__)

    logger.info('Reading from configuration file: {}'.format(configfile))

    with open(configfile, 'r') as cf:
        config: Dict[str, Any] = json.load(cf)
    return config


def validate_arguments(
    args: argparse.Namespace, config: dict, checkpoint: bool
) -> argparse.Namespace:
    """Checks for executables as well as makes sure files exist.

    input  - args from argparse; config from config file
    return - validated args from argparse
    """
    # config dict is from file
    # Determine if defaults need to be set
    # e = 1e-5
    # m = 500
    # c = 50
    # p = tblastn
    # t = 1

    logger = logging.getLogger(__name__)
    logger.debug('Validating & reconciling arguments.')
    error: bool = False

    autoMLSA = os.path.join(args.rundir, '.autoMLSA')
    checkpoint_dir = os.path.join(autoMLSA, 'checkpoint')

    if not os.path.exists(autoMLSA):
        os.mkdir(autoMLSA)

    if not os.path.exists(checkpoint_dir):
        os.mkdir(checkpoint_dir)

    # Reconcile config file and command line options
    changed: bool = False
    if args.evalue is None:
        if 'evalue' in config.keys() and config['evalue'] is not None:
            args.evalue = float(config['evalue'])
        else:
            args.evalue = 1e-5
    else:
        try:
            check = args.evalue != float(config['evalue'])
        except KeyError:
            pass
        else:
            if check:
                changed = True
    if args.coverage is None:
        if 'coverage' in config.keys() and config['coverage'] is not None:
            args.coverage = int(config['coverage'])
        else:
            args.coverage = 50
    else:
        try:
            check = args.coverage != int(config['coverage'])
        except KeyError:
            pass
        else:
            if check:
                changed = True
    if args.identity is None:
        if 'identity' in config.keys() and config['identity'] is not None:
            args.identity = int(config['identity'])
        else:
            args.identity = 30
    else:
        try:
            check = args.identity != int(config['identity'])
        except KeyError:
            pass
        else:
            if check:
                changed = True
    if args.allow_missing is None:
        if 'allow_missing' in config.keys() and config['allow_missing'] is not None:
            args.allow_missing = int(config['allow_missing'])
        else:
            args.allow_missing = 0
    else:
        try:
            check = args.allow_missing != int(config['allow_missing'])
        except KeyError:
            pass
        else:
            if check:
                changed = True
    if not args.missing_check:
        if 'missing_check' in config.keys() and config['missing_check']:
            args.allow_missing = True
        else:
            args.missing_check = False
    if args.program is None:
        if 'program' in config.keys():
            if config['program'] in ('tblastn', 'blastn'):
                args.program = config['program']
            else:
                msg = 'Program specified in config file "{}" is not valid.'
                logger.error(msg.format(config['program']))
                msg = 'Please give either "tblastn" or "blastn" and try again.'
                logger.error(msg)
                end_program(78)
        else:
            args.program = 'tblastn'
    else:
        try:
            check = args.program != config['program']
        except KeyError:
            pass
        else:
            if check:
                changed = True
    if args.threads is None:
        if 'threads' in config.keys() and config['threads'] is not None:
            args.threads = int(config['threads'])
        else:
            args.threads = 1
    if args.outgroup is None:
        if 'outgroup' in config.keys() and config['outgroup'] is not None:
            args.outgroup = config['outgroup']
        else:
            args.outgroup = ''
    else:
        try:
            check = args.outgroup != config['outgroup']
        except KeyError:
            pass
        else:
            if check:
                pass

    if args.iqtree is None:
        if 'iqtree' in config.keys() and config['iqtree'] is not None:
            args.iqtree = config['iqtree']
        else:
            args.iqtree = IQTREE
    else:
        try:
            check = args.iqtree != config['iqtree']
        except KeyError:
            pass
        else:
            if check:
                changed = True

    if args.mafft is None:
        if 'mafft' in config.keys() and config['mafft'] is not None:
            args.mafft = config['mafft']
        else:
            args.mafft = MAFFT
    else:
        try:
            check = args.mafft != config['mafft']
        except KeyError:
            pass
        else:
            if check:
                changed = True

    if changed:
        logger.debug('One or more settings was changed from config.')
        remove_intermediates(args.runid, ['genome'], args.protect)

    if args.checkpoint is None:
        if 'checkpoint' in config.keys() and config['checkpoint'] is not None:
            args.checkpoint = config['checkpoint']
        else:
            args.checkpoint = 'none'
        if args.checkpoint == 'validate':
            checkpoint = True
    if args.protect is None:
        if 'protect' in config.keys() and config['protect'] is not None:
            args.protect = config['protect']
        else:
            args.protect = False

    if args.evalue > 10:
        msg = 'Specified evalue "{}" is greater than 10.'
        logger.error(msg.format(args.evalue))
        msg = 'Please specify an evalue < 10 and try again.'
        logger.error(msg)
        end_program(78)
    if args.coverage not in range(0, 100):
        msg = 'Coverage value is not between 0 and 100.'
        logger.error(msg.format(args.coverage))
        msg = 'Please specify a coverage between 0 and 100 and try again.'
        logger.error(msg)
        end_program(78)

    cmd = 'lscpu | grep -G "^CPU(s):" | grep -o -E "[0-9]+"'
    try:
        maxthreads = (
            subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
            .decode()
            .strip()
        )
    except subprocess.CalledProcessError:
        msg = 'Unable to check number of available threads.'
        logger.warning(msg)
        msg = 'Make sure the number of specified threads is correct.'
        logger.warning(msg)
    else:
        logger.debug('Max threads found {}'.format(maxthreads))
        if args.threads > int(maxthreads):
            msg = 'Threads specified {} greater than number of available ' 'threads {}'
            logger.error(msg.format(args.threads, maxthreads))
            msg = 'Specify threads less than or equal to {} and try again.'
            logger.error(msg.format(maxthreads))
            end_program(78)

    # Set up lists for files, dirs, and queries
    if args.files is None:
        args.files = []
    if 'files' in config.keys() and config['files'] is not None:
        for f in config['files']:
            if f:
                args.files.append(f)

    if args.dir is None:
        args.dir = []
    if 'dir' in config.keys() and config['dir'] is not None:
        for d in config['dir']:
            if d:
                args.dir.append(d)

    if args.query is None:
        args.query = []
    if 'query' in config.keys() and config['query'] is not None:
        for q in config['query']:
            if q:
                args.query.append(q)
    for query in args.query.copy():
        if not check_if_fasta(query):
            msg = 'Specified query file {} does not appear to be FASTA file.'
            logger.error(msg.format(query))
            logger.debug('Removing {} from query list'.format(query))
            args.query.remove(query)
            error = True

    args.fasta = []
    # Combine lists
    for fastadir in args.dir:
        fastas = get_fasta_list(fastadir, args.rundir)
        args.fasta.extend(fastas)
    args.fasta.extend(args.files)

    # Remove duplicates
    args.fasta = list(dict.fromkeys(args.fasta))
    args.dir = list(dict.fromkeys(args.dir))
    args.query = list(dict.fromkeys(args.query))

    seen: Dict[str, str] = {}
    for fasta in args.fasta:
        base = os.path.basename(fasta)
        if base in seen:
            msg = 'Same genome name found in two locations'
            logger.error(msg)
            msg = '++++++++++++Offending files++++++++++++'
            logger.error(msg)
            logger.error('{}'.format(fasta))
            logger.error('{}'.format(seen[base]))
            msg = '+++++++++++++++++++++++++++++++++++++++'
            logger.error(msg)
            msg = 'Remove or rename one of these to continue.'
            logger.error(msg)
            end_program(65)
        else:
            seen[base] = fasta

    if 'dups' in config.keys() and config['dups']:
        args.dups = True

    if not args.fasta:
        msg = 'No valid FASTA files provided as --dir or --files.'
        logger.error(msg)
        error = True

    logger.debug('Validated arguments:')
    logger.debug('\n' + pprint.pformat(vars(args), indent=4))

    write_config(args)
    if error:
        end_program(-1)

    if checkpoint:
        checkpoint_reached('at validate arguments')

    return args


def write_config(args: argparse.Namespace) -> None:
    """Writes new config file, overwriting old if necessary

    input  - valiated & compiled arguments
    return - NULL
    """
    logger = logging.getLogger(__name__)
    logger.debug('Writing config file {}'.format(args.configfile))
    configdict = vars(args).copy()
    configdict.pop('config')
    configdict.pop('debug')
    configdict.pop('quiet')
    configdict.pop('configfile')
    configdict.pop('rundir')
    configdict.pop('checkpoint')
    json_writer(args.configfile, configdict)


def get_labels(rundir: str, fastas: List[str]) -> List[str]:
    """Determines the arbitrary labels per fasta file
    labels are index of list

    input - fasta files
    return - label list
    """
    logger = logging.getLogger(__name__)
    labelsf: str = os.path.join(rundir, '.autoMLSA', 'labels.json')
    bases: List[str] = [os.path.basename(x) for x in fastas]
    if os.path.exists(labelsf):
        with open(labelsf, 'r') as lf:
            labels: List[str] = json.load(lf)
        for base in bases:
            if base not in labels:
                logger.debug('Adding {} to list of labels'.format(base))
                labels.append(base)
    else:
        logger.debug('Generating list of labels {}'.format(labelsf))
        labels = bases
    json_writer(labelsf, labels)
    return labels
