#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import logging
import shlex
import subprocess
from .helper_functions import checkpoint_reached, checkpoint_tracker, end_program
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
            cmdstr = shlex.join(cmd) + ' > {}'.format(shlex.quote(outname))
            logger.debug(cmdstr)
            if checkpoint == 'prealign':
                cmdstrs.append(cmdstr)
            else:
                # Write to a temporary file so a failed run is not mistaken
                # for a finished alignment when the analysis is resumed.
                tmpname: str = '{}.tmp'.format(outname)
                with open(tmpname, 'w') as fh, open(logname, 'w') as logfh:
                    result = subprocess.run(
                        cmd, stdout=fh, stderr=logfh, stdin=subprocess.DEVNULL
                    )
                if result.returncode != 0 or os.path.getsize(tmpname) == 0:
                    os.remove(tmpname)
                    msg = 'mafft failed on {} (exit code {}). Check {} for details.'
                    logger.critical(msg.format(unalign, result.returncode, logname))
                    end_program(73)
                os.replace(tmpname, outname)
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
