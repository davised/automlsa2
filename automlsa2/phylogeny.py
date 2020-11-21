#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import logging
import subprocess
from .helper_functions import (
    checkpoint_reached,
    checkpoint_tracker,
    end_program,
)
from typing import List


def generate_nexus(runid: str, aligned: List[str], checkpoint: bool) -> str:
    """
    Generate nexus file containing partitions per aligned gene

    input  - list of aligned FASTA filenames
    return - path to nexus file
    """
    logger = logging.getLogger(__name__)
    # Example nexus file
    # #nexus
    # begin sets;
    #         charset part1 = rpoB.all.aln: *;
    #         charset part2 = ppsA.all.aln: *;
    #         charset part3 = gyrB.all.aln: *;
    #         charset part4 = dnaE.all.aln: *;
    #         charset part5 = pyrC.all.aln: *;
    #         charset part6 = aroE.all.aln: *;
    #         charset part7 = acsA.all.aln: *;
    #         charset part8 = recA.all.aln: *;
    # end;

    nexus: str = '{}.nex'.format(runid)
    if os.path.exists(nexus):
        msg = 'Nexus file {} already found, skipping nexus file generation.'
        logger.info(msg.format(nexus))
    else:
        with open(nexus, 'w') as nex:
            nex.write('#nexus\n')
            nex.write('begin sets;\n')
            for fn in aligned:
                base = os.path.basename(os.path.splitext(fn)[0])
                nex.write('\tcharset {} = {}: *;\n'.format(base, fn))
            nex.write('end;\n')
    if checkpoint:
        checkpoint_reached('after generating nexus file {}'.format(nexus))

    checkpoint_tracker('generate_nexus')
    return nexus


def run_iqtree(
    threads: int, iqtree: str, nexus: str, outgroup: str, opts: List[str]
) -> str:
    """
    Runs external iqtree command to generate phylogeny

    input  - nexus file
    return - path to output file
    """
    logger = logging.getLogger(__name__)
    out_tree: str = '{}.treefile'.format(nexus)
    logger.info('Generating phylogenetic tree ({}).'.format(out_tree))
    # cmd: List[str] = [
    #     iqtree,
    #     '-p', nexus,
    #     '-m', 'MFP',
    #     '-B', '1000',
    #     '-alrt', '1000',
    #     '--msub', 'nuclear',
    #     '--merge', 'rclusterf',
    #     '-nt', str(threads)]
    cmd: List[str] = [iqtree, '-p', nexus, '-nt', str(threads)]
    cmd.extend(opts)
    if outgroup:
        cmd.extend(['-o', outgroup])
    if os.path.exists(out_tree):
        logger.info('Treefile {} already found. Skipping iqtree.'.format(out_tree))
    else:
        logger.info('Running: {}'.format(' '.join(cmd)))
        subprocess.run(cmd)
        if not os.path.exists(out_tree):
            msg = 'iqtree2 seems to have failed.'
            logger.critical(msg)
            msg = (
                'Check the log files for error messages to see if they can '
                'be resolved.'
            )
            logger.critical(msg)
            end_program(73)

    checkpoint_tracker('run_iqtree')

    return out_tree
