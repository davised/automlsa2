# -*- coding: utf-8 -*-
import os
import stat

import pytest

from automlsa2._exceptions import ProgramMismatchError
from automlsa2.validate_requirements import (
    check_program,
    find_program,
    get_program_version,
    validate_requirements,
)

VERSION_OUTPUT = {
    'tblastn': 'tblastn: 2.17.0+\n Package: blast 2.17.0, build Aug 22 2026 11:52:50\n',
    'blastn': 'blastn: 2.17.0+\n Package: blast 2.17.0, build Aug 22 2026 11:52:50\n',
    'makeblastdb': 'makeblastdb: 2.17.0+\n Package: blast 2.17.0, build Aug 22 2026\n',
    'mafft': 'v7.526 (2024/Apr/26)\n',
    'iqtree3': 'IQ-TREE version 3.1.4 for Linux x86 64-bit built Sep 29 2026\n',
}


def make_exe(directory, name, output):
    path = os.path.join(directory, name)
    with open(path, 'w') as fh:
        fh.write("#!/bin/sh\nprintf '%s' '{}'\n".format(output))
    os.chmod(path, os.stat(path).st_mode | stat.S_IEXEC)
    return path


@pytest.mark.parametrize(
    'verstr,program,expected',
    [
        (VERSION_OUTPUT['tblastn'], 'tblastn', '2.17.0'),
        ('makeblastdb: 2.10.0+\n Package: blast 2.10.0\n', 'makeblastdb', '2.10.0'),
        (VERSION_OUTPUT['mafft'], 'mafft', '7.526'),
        (VERSION_OUTPUT['iqtree3'], 'iqtree', '3.1.4'),
        (
            'IQ-TREE multicore version 2.1.1 COVID-edition for Linux 64-bit built '
            'Aug 22 2020\n',
            'iqtree',
            '2.1.1',
        ),
        ('IQ-TREE multicore version 2.2.2.2 for Linux 64-bit\n', 'iqtree', '2.2.2.2'),
    ],
)
def test_get_program_version(verstr, program, expected):
    assert str(get_program_version(verstr, program)) == expected


def test_get_program_version_invalid():
    with pytest.raises(Exception):
        get_program_version('garbage', 'mafft')


def test_find_program_prefers_external(tmp_path, monkeypatch):
    pathdir = tmp_path / 'path'
    extdir = tmp_path / 'external'
    pathdir.mkdir()
    extdir.mkdir()
    make_exe(str(pathdir), 'iqtree3', '')
    ext_iqtree2 = make_exe(str(extdir), 'iqtree2', '')
    monkeypatch.setenv('PATH', str(pathdir))
    assert find_program(['iqtree3', 'iqtree2'], str(extdir)) == ext_iqtree2
    assert find_program(['iqtree3', 'iqtree2']) == str(pathdir / 'iqtree3')
    assert find_program(['missing']) is None


def test_check_program_too_old(tmp_path):
    exe = make_exe(str(tmp_path), 'mafft', 'v7.310 (2017/Mar/17)\n')
    with pytest.raises(ProgramMismatchError):
        check_program(exe, 'mafft', '--version')


def test_check_program_missing(tmp_path):
    assert not check_program(str(tmp_path / 'nope'), 'mafft', '--version')


def test_validate_requirements(tmp_path, monkeypatch):
    for name, output in VERSION_OUTPUT.items():
        make_exe(str(tmp_path), name, output)
    monkeypatch.setenv('PATH', str(tmp_path))
    exes = validate_requirements()
    assert exes['iqtree'] == str(tmp_path / 'iqtree3')
    assert exes['mafft'] == str(tmp_path / 'mafft')
    assert set(exes) == {'tblastn', 'blastn', 'makeblastdb', 'mafft', 'iqtree'}


def test_validate_requirements_missing(tmp_path, monkeypatch):
    monkeypatch.setenv('PATH', str(tmp_path))
    with pytest.raises(SystemExit) as err:
        validate_requirements()
    assert err.value.code == 78
