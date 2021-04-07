#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import logging
import shlex
import subprocess
from .helper_functions import checkpoint_reached, checkpoint_tracker
from rich.progress import track

# from .pbar import track_wide
from typing import List


def run_mafft(
    threads: int,
    mafft: str,
    unaligned: List[str],
    checkpoint: bool,
    opts: List[str],
) -> List[str]:
    """
    input  - unaligned fasta files per query
    return - list of aligned files per query
    """
    logger = logging.getLogger(__name__)
    base_cmd: List[str] = [mafft, '--thread', str(threads)]
    base_cmd.extend(opts)
    aligneddir: str = 'aligned'
    aligned: List[str] = []
    cmdstrs: List[str] = []
    logger.info('Aligning FASTA sequences, if necessary.')
    if not os.path.exists(aligneddir):
        os.mkdir(aligneddir)
    for unalign in track(unaligned, 'Running mafft...'):
        outname: str = '{}/{}.aln'.format(
            aligneddir, os.path.basename(os.path.splitext(unalign)[0])
        )
        logname: str = '{}.log'.format(outname)
        logger.debug('Aligning and writing to {}'.format(outname))
        aligned.append(outname)
        if os.path.exists(outname):
            msg = 'Aligned file {} already found, skipping.'
            logger.debug(msg.format(outname))
        else:
            cmd: List[str] = base_cmd + [unalign]
            cmdstr = ' '.join([shlex.quote(x) for x in cmd]) + ' > {}'.format(
                shlex.quote(outname)
            )
            logger.debug(cmdstr)
            if checkpoint == 'prealign':
                cmdstrs.append(cmdstr)
            else:
                with open(outname, 'w') as fh, open(logname, 'w') as logfh:
                    subprocess.run(cmd, stdout=fh, stderr=logfh, text=True)
    if cmdstrs:
        with open('mafftcmds.txt', 'w') as fh:
            fh.write('\n'.join(cmdstrs))
        logger.info('MAFFT alignment commands written to mafftcmds.txt.')
        logger.info('Run these commands and resubmit to continue.')
        checkpoint_reached('prior to MAFFT alignments')
    elif checkpoint == 'postalign':
        checkpoint_reached('after MAFFT alignments')

    checkpoint_tracker('run_mafft')

    return aligned
