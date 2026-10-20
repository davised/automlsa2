#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import re
import shutil
import subprocess
import logging
import packaging.version as pv
from typing import Dict, List, Optional, Tuple
from .helper_functions import end_program
from ._exceptions import ProgramMismatchError
from ._versions import (
    BLASTVERMIN,
    MAFFTVERMIN,
    IQTREEVERMIN,
    BLASTVERTESTED,
    MAFFTVERTESTED,
    IQTREEVERTESTED,
)

INSTALL_HINT = (
    'Install the external dependencies with conda: '
    'conda install -c conda-forge -c bioconda blast mafft iqtree'
)

# program key -> (candidate executable names, version flag, min ver, tested ver)
PROGRAMS: Dict[str, Tuple[List[str], str, str, str]] = {
    'tblastn': (['tblastn'], '-version', BLASTVERMIN, BLASTVERTESTED),
    'blastn': (['blastn'], '-version', BLASTVERMIN, BLASTVERTESTED),
    'makeblastdb': (['makeblastdb'], '-version', BLASTVERMIN, BLASTVERTESTED),
    'mafft': (['mafft'], '--version', MAFFTVERMIN, MAFFTVERTESTED),
    'iqtree': (['iqtree3', 'iqtree2', 'iqtree'], '--version', IQTREEVERMIN, IQTREEVERTESTED),
}


def get_program_version(verstr: str, program_name: str) -> pv.Version:
    """
    Parses the version string printed by an external program
    input  - output of the version command, program key
    return - parsed version
    """
    # Refs from October 2026
    # mafft: v7.526 (2024/Apr/26)
    # tblastn: tblastn: 2.17.0+
    #           Package: blast 2.17.0, build Aug 22 2026 11:52:50
    # iqtree3: IQ-TREE version 3.1.4 for Linux x86 64-bit built Sep 29 2026
    # iqtree2: IQ-TREE multicore version 2.1.1 COVID-edition for Linux 64-bit ...
    match: Optional[re.Match]
    if program_name in ('tblastn', 'blastn', 'makeblastdb'):
        match = re.search(r'^\S+:\s+(\d+(?:\.\d+)+)\+?', verstr, re.MULTILINE)
    elif program_name == 'mafft':
        match = re.search(r'v(\d+\.\d+)', verstr)
    elif program_name == 'iqtree':
        match = re.search(r'IQ-TREE.*?version (\d+(?:\.\d+)+)', verstr)
    else:
        raise ValueError('Unknown program: {}'.format(program_name))
    if match is None:
        raise pv.InvalidVersion(
            'Unable to determine version for {}: {}'.format(program_name, verstr)
        )
    return pv.Version(match.group(1))


def find_program(names: List[str], external: str = '') -> Optional[str]:
    """
    Finds the first available executable, searching the --external directory
    (if given) before $PATH.
    input  - candidate executable names, optional directory
    return - absolute path to executable, or None
    """
    for name in names:
        if external:
            exe = shutil.which(name, path=external)
            if exe:
                return os.path.abspath(exe)
    for name in names:
        exe = shutil.which(name)
        if exe:
            return os.path.abspath(exe)
    return None


def check_program(exe: str, program_name: str, version_flag: str) -> bool:
    """
    Runs the version command for an executable and compares it to the
    supported versions.
    input  - path to executable, program key, version flag
    return - True if usable; raises ProgramMismatchError if too old
    """
    logger = logging.getLogger(__name__)
    _, _, minver, testedver = PROGRAMS[program_name]
    try:
        verstr: str = subprocess.run(
            [exe, version_flag],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            stdin=subprocess.DEVNULL,
            text=True,
            timeout=60,
        ).stdout
    except FileNotFoundError:
        logger.error('Unable to find {} executable ({})'.format(program_name, exe))
        return False
    except PermissionError:
        msg = (
            'Unable to run {} executable ({}). Make sure it is installed properly '
            'and has the correct permissions to run.'
        )
        logger.error(msg.format(program_name, exe))
        return False
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as err:
        logger.error('Unable to check {} version ({}): {}'.format(program_name, exe, err))
        return False

    ver = get_program_version(verstr, program_name)
    refver = pv.Version(minver)
    logger.debug('Program version found: {} -> {}; min -> {}'.format(program_name, ver, refver))
    if ver < refver:
        raise ProgramMismatchError(
            '{} version {} needs to be updated to {} at minimum. {}'.format(
                program_name, ver, refver, INSTALL_HINT
            )
        )
    if ver.major > pv.Version(testedver).major:
        msg = '{} version {} is newer than tested software ({}). Proceed at your own risk.'
        logger.warning(msg.format(program_name, ver, testedver))
    return True


def validate_requirements(external: str = '') -> Dict[str, str]:
    """
    Identifies installed software and sets paths for running

    input  - optional directory to search before $PATH
    output - exes dict with paths to executables
    """
    logger = logging.getLogger(__name__)
    exes: Dict[str, str] = {}
    for program_name, (names, version_flag, _, _) in PROGRAMS.items():
        exe = find_program(names, external)
        if exe is None or not check_program(exe, program_name, version_flag):
            msg = 'Unable to find a working {} ({}) in $PATH{}.'
            where = ' or {}'.format(external) if external else ''
            logger.error(msg.format(program_name, ' or '.join(names), where))
            logger.error(INSTALL_HINT)
            end_program(78)
        logger.debug('Using {}: {}'.format(program_name, exe))
        exes[program_name] = exe  # type: ignore
    return exes
