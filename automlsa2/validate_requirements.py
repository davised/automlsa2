#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import print_function
import os.path
import subprocess
import logging
import platform
import urllib.request
import shutil
import re
import packaging.version as pv
from pathlib import Path
from hashlib import md5
from typing import Dict, Tuple, Union, cast
from .helper_functions import end_program
from ._exceptions import ProgramMismatchError
from ._versions import (
    BLASTVERMIN,
    MAFFTVERMIN,
    IQTREEVERMIN,
    BLASTVERCURR,
    MAFFTVERCURR,
    IQTREEVERCURR,
)

EXTERNALDIR = os.path.join(os.path.expanduser('~'), '.local', 'external')


def get_external_path(external: str = '') -> str:
    """
    Gets path to external directory
    input  - None
    return - Path to external dir
    """

    if not external:
        external = EXTERNALDIR

    if not os.path.exists(external):
        os.makedirs(external)
    return external


def get_program_version(
    verstr: str, program_name: str
) -> Tuple[pv.Version, pv.Version]:
    # Refs from December 2020
    # mafft: v7.471 (2020/Jul/3)
    # tblastn: tblastn: 2.10.0+
    #           Package: blast 2.10.0, build Mar 15 2020 22:36:43
    # iqtree2: IQ-TREE multicore version 2.1.1 COVID-edition for Linux 64-bit built Aug 22 2020
    #          Developed by Bui Quang Minh, James Barbetti, Nguyen Lam Tung,
    #          Olga Chernomor, Heiko Schmidt, Dominik Schrempf, Michael Woodhams.
    logger = logging.getLogger(__name__)
    ver: Union[pv.Version, pv.LegacyVersion]
    refver: Union[pv.Version, pv.LegacyVersion]

    if program_name in ('tblastn', 'blastn', 'makeblastdb'):
        ver_line = verstr.splitlines()[0]
        ver = pv.parse(ver_line.split()[1].rstrip('+'))
        refver = pv.parse(BLASTVERMIN)
    elif program_name == 'mafft':
        ver_line = verstr
        ver = pv.parse(ver_line.split()[0])
        refver = pv.parse(MAFFTVERMIN)
    elif program_name == 'iqtree':
        ver_line = verstr.splitlines()[0]
        ver_search: str = re.search('version (.+?) ', ver_line).group(1)  # type: ignore
        ver = pv.parse(ver_search)
        refver = pv.parse(IQTREEVERMIN)
    if isinstance(ver, pv.Version):
        logger.debug(
            'Program version found: {} -> {}; ref -> {}'.format(
                program_name, ver, refver
            )
        )
        return (cast(pv.Version, ver), cast(pv.Version, refver))
    else:
        raise pv.InvalidVersion(
            'Unable to determine accurate version for {}: {}'.format(program_name, ver)
        )


def check_program(exe: str, program_name: str, version_flag: str) -> bool:
    logger = logging.getLogger(__name__)
    try:
        verstr: str = subprocess.run(
            [exe, version_flag],
            check=True,
            stdout=subprocess.PIPE,
            text=True,
            stderr=subprocess.STDOUT,
        ).stdout
    except FileNotFoundError:
        msg = 'Unable to find {} executable in path or in provided dir ({})'
        logger.error(msg.format(program_name, exe))
        success = False
    except PermissionError:
        msg = (
            'Unable to run {} executable in ({}). Make sure it is installed properly '
            'and has the correct permissions to run.'
        )
        logger.error(msg.format(program_name, exe))
        success = False
    else:
        ver, refver = get_program_version(verstr, program_name)
        if ver.major == refver.major:
            if ver >= refver:
                logger.debug('{} found: {}'.format(program_name, ver))
                success = True
            else:
                raise ProgramMismatchError(
                    '{} version {} needs to be updated to {} at minimum.'.format(
                        program_name, ver, refver
                    )
                )
        elif ver.major < refver.major:
            raise ProgramMismatchError(
                '{} version {} needs to be updated to {} at minimum.'.format(
                    program_name, ver, refver
                )
            )
        else:
            msg = (
                '{} version {} is newer than tested software. Proceed at your own risk.'
            )
            logger.warning(msg.format(program_name, ver))

    return success


def validate_requirements(external: str = '') -> Dict[str, str]:
    """
    Identifies installed software and sets paths for running

    input  - logger
    output - exes dict with paths to executables
    """
    logger = logging.getLogger(__name__)

    # Set local path to blast executable DIR
    BLASTPATH: str = ''
    externalpath = get_external_path(external)
    blastcheck = os.path.join(
        externalpath, 'ncbi-blast-{}+'.format(BLASTVERCURR), 'bin'
    )
    if os.path.exists(blastcheck):
        logger.debug('BLAST found in external dir.')
        BLASTPATH = blastcheck

    if platform.system() == 'Linux':
        mafftsuffix = 'linux64'
        iqtreesuffix = os.path.join('{}-Linux'.format(IQTREEVERCURR), 'bin')
    elif platform.system() == 'Darwin':
        mafftsuffix = 'mac'
        iqtreesuffix = os.path.join('{}-MacOSX'.format(IQTREEVERCURR), 'bin')
    elif platform.system() == 'Windows':
        mafftsuffix = ''
        iqtreesuffix = os.path.join('{}-Windows'.format(IQTREEVERCURR), 'bin')
    else:
        logger.warning('Unable to determine os version.')
        mafftsuffix = ''
        iqtreesuffix = ''

    mafftcheck = os.path.join(externalpath, 'mafft-{}'.format(mafftsuffix))
    iqtreecheck = os.path.join(externalpath, 'iqtree-{}'.format(iqtreesuffix))

    # Expects in $PATH, can hard-code to full path here
    tblastn: str = 'tblastn'
    blastn: str = 'blastn'
    makeblastdb: str = 'makeblastdb'
    mafft: str = 'mafft'
    if os.path.exists(mafftcheck):
        logger.debug('MAFFT found in external dir.')
        mafft = os.path.join(mafftcheck, 'mafft.bat')
    iqtree: str = 'iqtree2'
    if os.path.exists(iqtreecheck):
        logger.debug('iqtree found in external dir.')
        iqtree = os.path.join(iqtreecheck, 'iqtree2')

    if os.path.exists(os.path.join(BLASTPATH, tblastn)):
        tblastn = os.path.join(BLASTPATH, tblastn)
    if os.path.exists(os.path.join(BLASTPATH, blastn)):
        blastn = os.path.join(BLASTPATH, blastn)
    if os.path.exists(os.path.join(BLASTPATH, makeblastdb)):
        makeblastdb = os.path.join(BLASTPATH, makeblastdb)

    logger.debug('Checking tblastn: {}'.format(tblastn))
    if not check_program(tblastn, 'tblastn', '-version'):
        msg = 'Please add tblastn to your $PATH or supply the path as --external'
        logger.error(msg)
        end_program(78)

    logger.debug('Checking blastn: {}'.format(blastn))
    if not check_program(blastn, 'blastn', '-version'):
        msg = 'Please add blastn to your $PATH or supply the path as --external'
        logger.error(msg)
        end_program(78)

    logger.debug('Checking makeblastdb: {}'.format(makeblastdb))
    if not check_program(makeblastdb, 'makeblastdb', '-version'):
        msg = 'Please add makeblastdb to your $PATH or supply the path as --external'
        logger.error(msg)
        end_program(78)

    # MAFFT testing
    logger.debug('Checking mafft: {}'.format(mafft))
    if not check_program(mafft, 'mafft', '--version'):
        msg = 'Please add mafft to your $PATH or supply the path as --external'
        logger.error(msg)
        end_program(78)

    # iqtree testing
    logger.debug('Checking iqtree: {}'.format(iqtree))
    if not check_program(iqtree, 'iqtree', '--version'):
        msg = 'Please add iqtree to your $PATH or supply the path as --external'
        logger.error(msg)
        end_program(78)

    exes: Dict[str, str] = {}
    exes['tblastn'] = tblastn
    exes['blastn'] = blastn
    exes['makeblastdb'] = makeblastdb
    exes['mafft'] = mafft
    exes['iqtree'] = iqtree
    return exes


def download_file(out_file: str, url: str, program: str) -> None:
    """
    Downloads given file from url
    """
    logger = logging.getLogger(__name__)
    logger.info('Downloading {} from {}.'.format(program, url))
    if not os.path.exists(out_file) or os.path.getsize(out_file) == 0:
        with urllib.request.urlopen(url) as response, open(out_file, 'wb') as tarfh:
            shutil.copyfileobj(response, tarfh)
    else:
        logger.info('{} already downloaded.'.format(program))


def install_blast(external: str = '') -> None:
    """
    Installs BLAST executables to external dir
    """
    logger = logging.getLogger(__name__)
    externalpath = get_external_path(external)
    base_url = 'https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/{}'.format(
        BLASTVERCURR
    )
    base_file = 'ncbi-blast-{}+-x64-'.format(BLASTVERCURR)

    # Assume Linux, update if necessary
    base_out_file = base_file + 'linux.tar.gz'
    if platform.system() == 'Darwin':
        base_out_file = base_file + 'macosx.tar.gz'
    elif platform.system() == 'Windows':
        base_out_file = base_file + 'win64.tar.gz'

    out_file = os.path.join(externalpath, base_out_file)
    blast_tar_url = '/'.join([base_url, base_out_file])
    md5_url = '/'.join([base_url, f'{base_out_file}.md5'])
    download_file(out_file, blast_tar_url, 'BLAST')

    blast_md5_digest = md5(Path(out_file).read_bytes()).hexdigest()
    logger.debug('calculated md5 digest: {}'.format(blast_md5_digest))

    download_file(out_file + '.md5', md5_url, 'BLAST md5')
    with open(out_file + '.md5', 'r') as md5fh:
        md5_hash = md5fh.readline().split(' ', 1)[0]
    logger.debug('given md5 digest: {}'.format(md5_hash))

    if blast_md5_digest != md5_hash:
        logger.info('BLAST download was unsuccessful. Delete file and try again.')
        end_program(-1)
    tblastn = os.path.join(
        externalpath, 'ncbi-blast-{}+'.format(BLASTVERCURR), 'bin', 'tblastn'
    )
    if not os.path.exists(tblastn):
        logger.info('Unpacking BLAST tar.gz file.')
        shutil.unpack_archive(out_file, externalpath, 'gztar')
    logger.info('Checking tblastn version.')
    try:
        tblastn_ver: bytes = subprocess.check_output(
            [tblastn, '-version'], stderr=subprocess.STDOUT
        )
    except FileNotFoundError:
        msg = 'Unable to find tblastn executable in path or in provided dir ({})'
        logger.error(msg.format(tblastn))
        msg = 'Please add tblastn to your $PATH or supply a BLASTPATH in {}'
        logger.error(msg.format(__file__))
        end_program(78)
    else:
        logger.info('tblastn found: {}'.format(tblastn_ver.strip().decode()))
        logger.info(
            'Add {} to your $PATH to have access to this program outside of automlsa2'.format(
                os.path.dirname(tblastn)
            )
        )


def install_mafft(external: str = '') -> None:
    logger = logging.getLogger(__name__)

    base_url = 'https://mafft.cbrc.jp/alignment/software'
    base_file = 'mafft-{}'.format(MAFFTVERCURR)
    if platform.system() == 'Linux':
        base_out_file = base_file + '-linux.tgz'
        md5_hash = '01ff001442faefdd5d1065973976ea3e'
    elif platform.system() == 'Darwin':
        base_out_file = base_file + '-mac.zip'
        md5_hash = '8206c78533286a4cd6a38821902fbb7a'
    elif platform.system() == 'Windows':
        logger.error('Windows is not currently supported for MAFFT install.')
        logger.error(
            'Look at install options here: '
            'https://mafft.cbrc.jp/alignment/software/windows.html'
        )
        return

    externalpath = get_external_path(external)
    out_file = os.path.join(externalpath, base_out_file)
    tar_url = '/'.join([base_url, base_out_file])
    download_file(out_file, tar_url, 'MAFFT')

    mafft_md5_digest = md5(Path(out_file).read_bytes()).hexdigest()
    logger.debug('calculated md5 digest: {}'.format(mafft_md5_digest))
    logger.debug('given md5 digest: {}'.format(md5_hash))
    if mafft_md5_digest != md5_hash:
        logger.info('MAFFT download was unsuccessful. Delete file and try again.')
        exit(-1)
    if platform.system() == 'Linux':
        mafft = os.path.join(externalpath, 'mafft-linux64', 'mafft.bat')
        if not os.path.exists(mafft):
            logger.info('Unpacking MAFFT tar.gz file.')
            shutil.unpack_archive(out_file, externalpath, 'gztar')
    elif platform.system() == 'Darwin':
        mafft = os.path.join(externalpath, 'mafft-mac', 'mafft.bat')
        if not os.path.exists(mafft):
            logger.info('Unpacking MAFFT zip file.')
            shutil.unpack_archive(out_file, externalpath, 'zip')
        logger.debug('Changing permissions.')
        for root, _, files in os.walk(os.path.join(externalpath, 'mafft-mac')):
            for exe in files:
                os.chmod(os.path.join(root, exe), 0o755)
    logger.info('Checking mafft version.')
    try:
        mafft_ver: bytes = subprocess.check_output(
            [mafft, '--version'], stderr=subprocess.STDOUT
        )
    except FileNotFoundError:
        msg = 'Unable to find mafft executable in path or in provided dir ({})'
        logger.error(msg.format(mafft))
        msg = 'Please add mafft to your $PATH or supply a BLASTPATH in {}'
        logger.error(msg.format(__file__))
        end_program(78)
    else:
        logger.info('mafft found: {}'.format(mafft_ver.strip().decode()))
        logger.info(
            'Add {} to your $PATH to have access to this program outside of automlsa2'.format(
                os.path.dirname(mafft)
            )
        )


def install_iqtree(external: str = '') -> None:
    logger = logging.getLogger(__name__)

    base_url = 'https://github.com/iqtree/iqtree2/releases/download/v{}'.format(
        IQTREEVERCURR
    )
    base_file = 'iqtree-{}'.format(IQTREEVERCURR)
    if platform.system() == 'Linux':
        base_out_file = base_file + '-Linux.tar.gz'
        md5_hash = 'fcb2d06c547e597a70c27aff06e63e38'
    elif platform.system() == 'Darwin':
        base_out_file = base_file + '-MacOSX.zip'
        md5_hash = '889709744247f079e5e79ad8885f78d0'
    elif platform.system() == 'Windows':
        base_out_file = base_file + '-Windows.zip'
        md5_hash = '99af9d49bbf4b90bb6a1739a305f5fac'

    externalpath = get_external_path(external)
    out_file = os.path.join(externalpath, base_out_file)
    tar_url = '/'.join([base_url, base_out_file])
    download_file(out_file, tar_url, 'iqtree')

    iqtree_md5_digest = md5(Path(out_file).read_bytes()).hexdigest()
    logger.debug('calculated md5 digest: {}'.format(iqtree_md5_digest))
    logger.debug('given md5 digest: {}'.format(md5_hash))
    if iqtree_md5_digest != md5_hash:
        logger.info('iqtree2 download was unsuccessful. Delete file and try again.')
        exit(-1)
    if platform.system() == 'Linux':
        iqtree = os.path.join(
            externalpath, 'iqtree-{}-Linux'.format(IQTREEVERCURR), 'bin', 'iqtree2'
        )
        if not os.path.exists(iqtree):
            logger.info('Unpacking iqtree tar.gz file.')
            shutil.unpack_archive(out_file, externalpath, 'gztar')
    elif platform.system() == 'Darwin':
        iqtree = os.path.join(
            externalpath, 'iqtree-{}-MacOSX'.format(IQTREEVERCURR), 'bin', 'iqtree2'
        )
        if not os.path.exists(iqtree):
            logger.info('Unpacking iqtree zip file.')
            shutil.unpack_archive(out_file, externalpath, 'zip')
        os.chmod(iqtree, 0o755)
    elif platform.system() == 'Windows':
        iqtree = os.path.join(
            externalpath,
            'iqtree-{}-Windows'.format(IQTREEVERCURR),
            'bin',
            'iqtree2.exe',
        )
        if not os.path.exists(iqtree):
            logger.info('Unpacking iqtree zip file.')
            shutil.unpack_archive(out_file, externalpath, 'zip')
    logger.info('Checking iqtree version.')
    try:
        iqtree_ver: bytes = subprocess.check_output(
            [iqtree, '--version'], stderr=subprocess.STDOUT
        )
    except FileNotFoundError:
        msg = 'Unable to find iqtree2 executable in path or in provided dir ({})'
        logger.error(msg.format(iqtree))
        msg = 'Please add iqtree2 to your $PATH or supply a full iqtree path in {}'
        logger.error(msg.format(__file__))
        end_program(78)
    else:
        logger.info('iqtree2 found: {}'.format(iqtree_ver.strip().decode()))
        logger.info(
            'Add {} to your $PATH to have access to this program outside of automlsa2'.format(
                os.path.dirname(iqtree)
            )
        )
