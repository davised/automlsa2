#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import argparse
import logging
import pprint
import sys
from rich.logging import RichHandler
from rich.console import Console
from .helper_functions import check_if_fasta
from .validate_requirements import (
    check_program,
    validate_requirements,
    install_blast,
    install_mafft,
    install_iqtree,
)
from .__version__ import __version__
from ._settings import MAFFT, IQTREE


class InstallDeps(argparse.Action):
    def __init__(self, nargs='?', **kw):
        super().__init__(nargs=nargs, **kw)

    def __call__(self, parser, namespace, values, option_string=None):
        if values:
            external = values
        else:
            external = os.path.join('~', '.local', 'external')
        if '~' in external:
            external = os.path.expanduser(external)
        external = os.path.abspath(external)
        if not os.path.basename(external) == 'external':
            external = os.path.join(external, 'external')

        init_logger(debug=False, quiet=False, rundir='', runid='')
        logger = logging.getLogger(__name__)
        logger.info('Checking for dependencies and installing.')
        install = False
        if check_program('tblastn', 'tblastn', '-version'):
            logger.info('tblastn already found, skipping install.')
        else:
            install_blast(external)
            install = True
        if check_program('mafft', 'mafft', '--version'):
            logger.info('mafft already found, skipping install.')
        else:
            install_mafft(external)
            install = True
        if check_program('iqtree2', 'iqtree', '--version'):
            logger.info('iqtree2 already found, skipping install.')
        else:
            install_iqtree(external)
            install = True
        validate_requirements(external)
        if install:
            msg = (
                'Please provide given path {} as the --external flag to '
                'use the installed programs.'
            )
            logger.info(msg.format(external))
            msg = '~/.local/external will be found automatically.'
            logger.info(msg)
        parser.exit()


def config_rundir(runid: str) -> str:
    """Identifies existing (or not) run directory from input runid

    input  - runid as given as input
    output - full path to the run directory
    """
    rundir: str = ''
    if os.path.exists('../{}'.format(runid)):
        rundir = os.path.abspath('../{}'.format(runid))
    else:
        rundir = os.path.abspath('./{}'.format(runid))
        if not os.path.exists(rundir):
            try:
                os.mkdir(rundir)
                os.mkdir(os.path.join(rundir, 'backup'))
                os.mkdir(os.path.join(rundir, '.autoMLSA'))
                os.mkdir(os.path.join(rundir, '.autoMLSA', 'backups'))
                os.mkdir(os.path.join(rundir, '.autoMLSA', 'checkpoint'))
            except OSError as ose:
                msg = 'Unable to generate rundir "{}" : {}'
                sys.stderr.write(msg.format(rundir, ose) + '\n')
                sys.stderr.write('Check your path and try again.\n')
                exit(71)
    os.chdir(rundir)
    return rundir


def extant_file(x: str) -> str:
    """
    'Type' for argparse - checks that file exists
    input  - path to file
    return - absolute path to file
    """
    if not os.path.exists(x):
        raise argparse.ArgumentTypeError("{} does not exist".format(x))
    return os.path.abspath(x)


def get_fasta_files(x: str) -> str:
    """
    'Type' for argparse - checks if file is FASTA format
    input  - path to file
    return - absolute path to FASTA file
    """
    if not os.path.exists(x):
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    elif check_if_fasta(x):
        return os.path.abspath(x)
    else:
        raise argparse.ArgumentTypeError("{0} is not a FASTA file".format(x))


def init_logger(debug: bool, quiet: bool, rundir: str, runid: str) -> None:
    """Sets up logging system to file and stderr

    input  - parsed arguments
    return - None; use logger = logging.getLogger(__name__)
    """

    # must set logging level to DEBUG to print DEBUG to file
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    stderr_handler = RichHandler(rich_tracebacks=True)
    if not sys.stderr.isatty():
        stderr_handler = RichHandler(rich_tracebacks=True, console=Console(width=119))
    if debug:
        stderr_handler.setLevel(logging.DEBUG)
    elif quiet:
        stderr_handler.setLevel(logging.WARNING)
    else:
        stderr_handler.setLevel(logging.INFO)
    logger.addHandler(stderr_handler)

    # Print to file as well, if provided.
    if rundir and runid:
        logfile = os.path.join(rundir, '{}.log'.format(runid))
        file_handler = RichHandler(console=Console(file=open(logfile, 'a'), width=119))
        logger.addHandler(file_handler)


def run_argparse() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='This is a rewrite of autoMLSA.pl. Generates automated '
        'multi-locus sequence analyses.',
        epilog='Quick usage: automlsa2 --dir genomedir --query queries.fas '
        '-t # -- runID',
    )
    parser.add_argument('runid', help='Name of the run directory.', type=str)
    parser.add_argument(
        '--query', help='Path to file with input seq(s).', type=extant_file, nargs='+'
    )
    parser.add_argument(
        '--files',
        help='Path to the target genome FASTA files.',
        type=get_fasta_files,
        nargs='+',
    )
    parser.add_argument(
        '--dir',
        help='Path to the target genome directory with FASTA files.',
        type=extant_file,
        nargs='+',
    )
    parser.add_argument(
        '-e', '--evalue', help='E-value cutoff for BLAST searches. [1e-5]', type=float
    )
    parser.add_argument(
        '-c', '--coverage', help='Sets the coverage cut-off threshold. [50]', type=int
    )
    parser.add_argument(
        '-i', '--identity', help='Sets the identity cut-off threshold. [30]', type=int
    )
    parser.add_argument(
        '-p',
        '--program',
        help='Which BLAST program to run. [tblastn]',
        choices=['blastn', 'tblastn'],
        type=str,
    )
    parser.add_argument(
        '--config', help='Path to configuration json file to copy.', type=extant_file
    )
    parser.add_argument(
        '--missing_check',
        help='Use this to confirm that settings have been '
        'checked when genes are missing.',
        action='store_true',
        default=False,
    )
    parser.add_argument(
        '-t', '--threads', help='Number of threads to use. [1]', type=int
    )
    parser.add_argument(
        '--dups',
        help='Allow for duplicate query names for more sequence '
        'coverage across disparate organisms.',
        action='store_true',
        default=False,
    )
    parser.add_argument(
        '--allow_missing', help='Allow for N missing genes per genome. [0]', type=int
    )
    parser.add_argument(
        '--outgroup', help='Name of outgroup file or strain to root on.', type=str
    )
    parser.add_argument(
        '--protect',
        help='Save files from getting overwritten. By default, '
        'as input files update, older alignments and trees are deleted.',
        action='store_true',
    )
    parser.add_argument(
        '--checkpoint',
        help='Name of stage to stop computing on. [none]',
        type=str,
        choices=[
            'validate',
            'preblast',
            'filtering',
            'prealign',
            'postalign',
            'nexus',
            'none',
        ],
    )
    parser.add_argument(
        '--install_deps',
        help='Install dependencies into given directory. ' '[~/.local/external]',
        action=InstallDeps,
    )
    parser.add_argument(
        '--external',
        help='Path to installed external programs. ' '[~/.local/external]',
        type=extant_file,
    )
    parser.add_argument('--mafft', help=f'mafft settings [{MAFFT}]', type=str)
    parser.add_argument('--iqtree', help=f'iqtree2 settings [{IQTREE}]', type=str)
    parser.add_argument(
        '--debug',
        help='Turn on debugging messages.',
        action='store_true',
        default=False,
    )
    parser.add_argument('--version', action='version', version=__version__)
    parser.add_argument(
        '--quiet',
        help='Turn off progress messages.',
        action='store_true',
        default=False,
    )
    args = parser.parse_args()
    args.rundir = config_rundir(args.runid)
    args.configfile = os.path.join(args.rundir, 'config.json')
    init_logger(args.debug, args.quiet, args.rundir, args.runid)
    logger = logging.getLogger(__name__)
    logger.debug('Started autoMLSA.py for run: {}'.format(args.runid))
    logger.debug(pprint.pformat(vars(args)))
    if not args.external:
        args.external = ''
    return args
