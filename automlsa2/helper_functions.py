#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import print_function
import os
import json
import shutil
import glob
import logging
import signal
import sys
import time
import psutil  # type: ignore
from typing import List
from hashlib import blake2b

SUFFIXES = ['.nsq', '.nin', '.nhr', '.nto', '.not', '.ndb', '.ntf', '.njs']


# Forked worker processes inherit this value, so they can tell they are not
# the main process.
MAIN_PID = os.getpid()


def kill_child_processes() -> None:
    """
    Terminates all descendant processes, e.g. BLAST workers and their BLAST
    searches, mafft, or iqtree.
    """
    try:
        children = psutil.Process().children(recursive=True)
    except psutil.Error:
        return
    for child in children:
        try:
            child.terminate()
        except psutil.NoSuchProcess:
            pass
    _, alive = psutil.wait_procs(children, timeout=5)
    for child in alive:
        try:
            child.kill()
        except psutil.NoSuchProcess:
            pass


def siginthandler(sig, frame):
    if os.getpid() != MAIN_PID:
        return
    # Ignore repeated Ctrl-C while cleaning up
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    kill_child_processes()
    end_program(2)


def _ignore_sigint(signal_num, frame):
    pass


def worker_init():
    """
    Initializer for worker processes. Interrupts are handled by the main
    process, which terminates the workers. A Python-level handler is used
    rather than SIG_IGN so external programs started by the workers get the
    default SIGINT behavior back when they exec.
    """
    signal.signal(signal.SIGINT, _ignore_sigint)


def remove_intermediates(runid: str, intermediates: List[str], protect: bool) -> None:
    logger = logging.getLogger(__name__)
    if protect:
        msg = 'One or more intermediates of {} type need to be removed'
        logger.info(msg.format(','.join(intermediates)))
        msg = (
            'Delete the "protect" flag in the config.json file or restart '
            'your analysis with another runid to continue.'
        )
        logger.info(msg)
        end_program(1)

    genome_tmpfiles = [
        'blast_results.tsv',
        'keepsidx.json',
        'single_copy.json',
        'blast_filtered.tsv',
        'expected_filt.json',
    ]
    if 'blast' in intermediates:
        for blastfile in glob.iglob(os.path.join('blast', '*.tab')):
            os.remove(blastfile)
    if 'genome' in intermediates:
        for deldir in ['unaligned', 'aligned']:
            if os.path.exists(deldir):
                shutil.rmtree(deldir)
        for fullname in [os.path.join('.autoMLSA', fn) for fn in genome_tmpfiles]:
            if os.path.exists(fullname):
                os.remove(fullname)
    delfiles: List[str] = sorted(glob.glob(glob.escape(runid) + '.nex*'))
    if delfiles:
        # Each backup gets its own directory so earlier backups are kept
        backupdir = os.path.join('backup', time.strftime('%Y%m%d-%H%M%S'))
        i = 1
        while os.path.exists(backupdir):
            backupdir = os.path.join('backup', time.strftime('%Y%m%d-%H%M%S') + f'-{i}')
            i += 1
        os.makedirs(backupdir)
        logger.info('Backing up nexus and tree files to {}.'.format(backupdir))
        for delfile in delfiles:
            shutil.move(delfile, backupdir)


def end_program(code: int) -> None:
    """
    Program message including success or failure of the program.
    """
    logger = logging.getLogger(__name__)
    if code == 1:
        msg = 'Program was stopped at an intermediate stage.'
        logger.info(msg)
    elif code == 0:
        msg = 'Program was a success! Congratulations!'
        logger.info(msg)
    elif code == 2:
        print('\033[?25h')
        msg = 'Program was interrupted using the keyboard.'
        logger.info(msg)
        msg = 'Partial files may exist. Re-running the analysis is suggested.'
        logger.info(msg)
    else:
        msg = 'Program exiting with code ({}) indicating failure.'
        logger.info(msg.format(code))
        msg = 'Check error messages to resolve the problem.'
        logger.info(msg)
    sys.exit(code)


def check_if_fasta(fa):
    """
    Checks first line of file for '>' to determine if FASTA, returns bool
    """
    for suffix in SUFFIXES:
        if suffix in fa:
            return False
    with open(fa, 'r') as fah:
        try:
            test = fah.readline().startswith('>')
        except UnicodeDecodeError:
            return False
        else:
            if test:
                return True
            else:
                return False


def json_writer(fn, x):
    """
    Writes object to file as json format.
    """
    with open(fn, 'w') as fh:
        json.dump(x, fh, indent=4)
        fh.write('\n')


def sanitize_path(s):
    """
    Turns string into filename without symbols
    """
    s = s.replace(' ', '_')
    filename = ''.join(x for x in s if (x.isalnum() or x in '._-'))
    return filename


def generate_hash(s):
    """
    Generates blake2b hash for sequence
    input  - nucleotide sequence
    return - blake2b hash digest
    """
    seqhash = blake2b(s.encode(), digest_size=16)
    return seqhash.hexdigest()


def checkpoint_reached(stage):
    logger = logging.getLogger(__name__)
    msg = 'Checkpoint reached {}. Stopping...'
    logger.info(msg.format(stage))
    end_program(1)


def checkpoint_tracker(fn_name):
    checkpath = os.path.join('.autoMLSA', 'checkpoint', fn_name)
    if not os.path.exists(checkpath):
        open(os.path.join(checkpath), 'w').close()


def exit_successfully(rundir: str, treefile: str) -> None:
    """Temporary command

    input  - rundir and treefile to print to log
    return - None
    """
    logger = logging.getLogger(__name__)

    msg = 'Your treefile is ready: {}/{}'
    logger.info(msg.format(rundir, treefile))
    end_program(0)
