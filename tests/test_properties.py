# -*- coding: utf-8 -*-
"""Property-based tests using hypothesis."""
import argparse
import contextlib
import json
import os
import string
import tempfile

import numpy as np
import pandas as pd
import packaging.version as pv
import pytest
from Bio import SeqIO
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from automlsa2.blast_functions import (
    print_blast_summary,
    print_fasta_files,
    query_name,
    reader,
)
from automlsa2.configuration import validate_arguments
from automlsa2.formatting import get_queries
from automlsa2.helper_functions import generate_hash, sanitize_path
from automlsa2.phylogeny import generate_nexus
from automlsa2.validate_requirements import get_program_version

SETTINGS = settings(
    max_examples=50,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow],
)

SAFE_CHARS = string.ascii_letters + string.digits + '._-'
seq_ids = st.text(SAFE_CHARS, min_size=1, max_size=20).filter(
    lambda s: not s.startswith('.')
)
aa_seqs = st.text('ACDEFGHIKLMNPQRSTVWY', min_size=1, max_size=60)


@contextlib.contextmanager
def rundir():
    """Temporary run directory laid out like config_rundir creates it."""
    old = os.getcwd()
    with tempfile.TemporaryDirectory() as tmp:
        for sub in ('backup', '.autoMLSA/backups', '.autoMLSA/checkpoint'):
            os.makedirs(os.path.join(tmp, sub))
        os.chdir(tmp)
        try:
            yield tmp
        finally:
            os.chdir(old)


# ----------------------------------------------------------------------------
# Version parsing
# ----------------------------------------------------------------------------
version_parts = st.lists(st.integers(0, 999), min_size=2, max_size=4)


@given(parts=version_parts, program=st.sampled_from(['tblastn', 'blastn', 'makeblastdb']))
def test_blast_version_roundtrip(parts, program):
    ver = '.'.join(map(str, parts))
    out = '{}: {}+\n Package: blast {}, build Jan 1 2026 00:00:00\n'.format(
        program, ver, ver
    )
    assert get_program_version(out, program) == pv.Version(ver)


@given(major=st.integers(0, 99), minor=st.integers(0, 999))
def test_mafft_version_roundtrip(major, minor):
    ver = '{}.{}'.format(major, minor)
    out = 'v{} (2024/Apr/26)\n'.format(ver)
    assert get_program_version(out, 'mafft') == pv.Version(ver)


@given(
    parts=version_parts,
    flavor=st.sampled_from(['IQ-TREE', 'IQ-TREE multicore', 'IQ-TREE multicore']),
    suffix=st.sampled_from(['for Linux x86 64-bit', 'COVID-edition for Linux 64-bit']),
)
def test_iqtree_version_roundtrip(parts, flavor, suffix):
    ver = '.'.join(map(str, parts))
    out = '{} version {} {} built Sep 29 2026\nDeveloped by ...\n'.format(
        flavor, ver, suffix
    )
    assert get_program_version(out, 'iqtree') == pv.Version(ver)


# ----------------------------------------------------------------------------
# Naming helpers
# ----------------------------------------------------------------------------
@given(st.text())
def test_sanitize_path_is_safe_and_idempotent(s):
    safe = sanitize_path(s)
    assert all(c.isalnum() or c in '._-' for c in safe)
    assert os.sep not in safe
    assert sanitize_path(safe) == safe


@given(seqid=seq_ids, seq=aa_seqs)
def test_query_name_roundtrip(seqid, seq):
    # get_queries names files {seqid}_{hash}.fas; BLAST reports seqid
    fn = os.path.join('queries', '{}_{}.fas'.format(seqid, generate_hash(seq)))
    assert query_name(fn) == seqid


@given(st.text())
def test_generate_hash(s):
    h = generate_hash(s)
    assert len(h) == 32
    assert all(c in string.hexdigits for c in h)
    assert generate_hash(s) == h


# ----------------------------------------------------------------------------
# BLAST tabular reader
# ----------------------------------------------------------------------------
field = st.text(
    st.characters(blacklist_categories=['Cs', 'Cc'], blacklist_characters='\t\n\r"'),
    max_size=10,
)


@SETTINGS
@given(
    rows=st.lists(st.tuples(seq_ids, st.lists(field, min_size=1, max_size=5))),
    comments=st.lists(field.map(lambda s: '# ' + s)),
)
def test_reader_skips_comments(rows, comments):
    with tempfile.TemporaryDirectory() as tmp:
        fn = os.path.join(tmp, 'out.tab')
        with open(fn, 'w') as fh:
            for c in comments:
                fh.write(c + '\n')
            for first, rest in rows:
                fh.write('\t'.join([first] + rest) + '\n')
        assert reader(fn) == [[first] + rest for first, rest in rows]


# ----------------------------------------------------------------------------
# Query extraction
# ----------------------------------------------------------------------------
@SETTINGS
@given(
    records=st.lists(
        st.tuples(seq_ids, aa_seqs),
        min_size=1,
        max_size=8,
        unique_by=(lambda r: r[0].lower(), lambda r: r[1]),
    )
)
def test_get_queries_writes_one_file_per_query(records):
    with rundir() as tmp:
        qfile = os.path.join(tmp, 'input.fas')
        with open(qfile, 'w') as fh:
            for seqid, seq in records:
                fh.write('>{} description\n{}\n'.format(seqid, seq))
        queries = get_queries('run', tmp, False, [qfile], False)

        assert len(queries) == len(records)
        for fn, (seqid, seq) in zip(queries, records):
            assert query_name(fn) == seqid
            recs = list(SeqIO.parse(fn, 'fasta'))
            assert [(r.id, str(r.seq)) for r in recs] == [(seqid, seq)]
        # Re-running with the same input is stable
        assert get_queries('run', tmp, False, [qfile], False) == queries


# ----------------------------------------------------------------------------
# BLAST summary and filtering
# ----------------------------------------------------------------------------
HEADERS = [
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
DTYPES = {
    'qseqid': 'category',
    'sseqid': 'category',
    'saccver': 'string',
    'pident': np.float64,
    'qlen': np.int64,
    'length': np.int64,
    'bitscore': np.float64,
    'qcovhsp': np.int64,
    'stitle': 'string',
    'sseq': 'string',
}


@st.composite
def blast_hits(draw):
    """Hit counts per (genome, query), plus a matching BLAST results table."""
    ngenomes = draw(st.integers(1, 5))
    nqueries = draw(st.integers(1, 4))
    queries = ['gene{}'.format(i) for i in range(nqueries)]
    counts = draw(
        st.lists(
            st.lists(st.integers(0, 3), min_size=nqueries, max_size=nqueries),
            min_size=ngenomes,
            max_size=ngenomes,
        )
    )
    rows = []
    for g, row in enumerate(counts):
        for q, n in enumerate(row):
            for _ in range(n):
                bitscore = draw(st.floats(10, 1000, allow_nan=False))
                seq = draw(st.text('ACDEFGHIKLMNPQRSTVWY-', min_size=1, max_size=20))
                rows.append(
                    [queries[q], str(g), 'acc', 90.0, 100, 90, bitscore, 90, 't', seq]
                )
    allow_missing = draw(st.integers(0, nqueries))
    return ngenomes, queries, counts, rows, allow_missing


@SETTINGS
@given(data=blast_hits())
def test_print_blast_summary_filters_by_missing(data):
    ngenomes, queries, counts, rows, allow_missing = data
    labels = ['genome{}.fna'.format(i) for i in range(ngenomes)]
    query_files = ['queries/{}_{}.fas'.format(q, generate_hash(q)) for q in queries]
    blastout = pd.DataFrame(rows, columns=HEADERS).astype(DTYPES)
    expected_keeps = [
        labels[g] for g, row in enumerate(counts) if row.count(0) <= allow_missing
    ]

    with rundir():
        args = (
            'run',
            blastout,
            labels,
            [str(i) for i in range(ngenomes)],
            query_files,
            allow_missing,
            True,
            False,
            False,
        )
        if not expected_keeps:
            with pytest.raises(SystemExit):
                print_blast_summary(*args)
            return
        kept = print_blast_summary(*args)

        # The presence matrix covers every genome and query, hits or not
        matrix = pd.read_csv('presence_matrix.tsv', sep='\t', index_col=0)
        assert list(matrix.index) == labels
        assert list(matrix.columns) == queries
        assert matrix.values.tolist() == counts

        with open(os.path.join('.autoMLSA', 'expected_filt.json')) as fh:
            assert sorted(json.load(fh)) == sorted(expected_keeps)

        # One best hit per kept (query, genome) pair
        kept_idx = {str(labels.index(k)) for k in expected_keeps}
        assert set(kept['sseqid'].astype(str)) <= kept_idx
        pairs = list(zip(kept['qseqid'].astype(str), kept['sseqid'].astype(str)))
        assert len(pairs) == len(set(pairs))
        expected_pairs = {
            (queries[q], str(g))
            for g, row in enumerate(counts)
            for q, n in enumerate(row)
            if n and str(g) in kept_idx
        }
        assert set(pairs) == expected_pairs
        for q, g in expected_pairs:
            best = blastout[
                (blastout.qseqid == q) & (blastout.sseqid == g)
            ].bitscore.max()
            got = kept[(kept.qseqid == q) & (kept.sseqid == g)].bitscore.iloc[0]
            assert got == best

        # Unaligned FASTA output has one ungapped sequence per kept genome
        unaligned = print_fasta_files(kept, labels)
        for fn in unaligned:
            q = os.path.splitext(os.path.basename(fn))[0]
            recs = list(SeqIO.parse(fn, 'fasta'))
            ids = [r.id for r in recs]
            assert len(ids) == len(set(ids))
            for r in recs:
                assert '-' not in str(r.seq)
                g = str(labels.index(r.id + '.fna'))
                assert (q, g) in expected_pairs


# ----------------------------------------------------------------------------
# Nexus partition file
# ----------------------------------------------------------------------------
@SETTINGS
@given(names=st.lists(seq_ids, min_size=1, max_size=10, unique=True))
def test_generate_nexus_lists_each_alignment(names):
    aligned = [os.path.join('aligned', n + '.aln') for n in names]
    with rundir():
        nexus = generate_nexus('run', aligned, False)
        with open(nexus) as fh:
            lines = fh.read().splitlines()
    assert lines[0] == '#nexus'
    charsets = [ln.strip() for ln in lines if ln.strip().startswith('charset')]
    assert charsets == [
        'charset {} = {}: *;'.format(n, a) for n, a in zip(names, aligned)
    ]


# ----------------------------------------------------------------------------
# Argument validation
# ----------------------------------------------------------------------------
def make_args(tmp, **kw):
    fasta = os.path.join(tmp, 'genome.fna')
    query = os.path.join(tmp, 'query.fas')
    for fn, seq in ((fasta, 'ACGT'), (query, 'MKV')):
        with open(fn, 'w') as fh:
            fh.write('>x\n{}\n'.format(seq))
    args = dict(
        runid='run',
        rundir=tmp,
        configfile=os.path.join(tmp, 'config.json'),
        config=None,
        debug=False,
        quiet=False,
        external='',
        evalue=1e-5,
        coverage=50,
        identity=30,
        allow_missing=0,
        missing_check=False,
        program='tblastn',
        threads=1,
        outgroup=None,
        iqtree=None,
        mafft=None,
        checkpoint='none',
        protect=None,
        files=[fasta],
        dir=None,
        query=[query],
        dups=False,
    )
    args.update(kw)
    return argparse.Namespace(**args)


def accepted(**kw):
    with rundir() as tmp:
        try:
            validate_arguments(make_args(tmp, **kw), {}, False)
        except SystemExit as err:
            assert err.code == 78
            return False
    return True


@SETTINGS
@given(st.integers(-50, 150))
def test_coverage_range(coverage):
    assert accepted(coverage=coverage) == (0 <= coverage <= 100)


@SETTINGS
@given(st.integers(-50, 150))
def test_identity_range(identity):
    assert accepted(identity=identity) == (0 <= identity <= 100)


@SETTINGS
@given(st.floats(-1, 20, allow_nan=False))
def test_evalue_range(evalue):
    assert accepted(evalue=evalue) == (0 < evalue <= 10)


@SETTINGS
@given(st.integers(-4, len(os.sched_getaffinity(0)) + 4))
def test_threads_range(threads):
    maxthreads = len(os.sched_getaffinity(0))
    assert accepted(threads=threads) == (1 <= threads <= maxthreads)


@SETTINGS
@given(
    cli=st.one_of(st.none(), st.just(True)),
    saved=st.one_of(st.none(), st.booleans()),
    allow_missing=st.integers(0, 5),
)
def test_config_reconciliation(cli, saved, allow_missing):
    """Saved config values are used unless overridden on the command line."""
    with rundir() as tmp:
        config = {'allow_missing': allow_missing}
        if saved is not None:
            config['protect'] = saved
            config['missing_check'] = saved
        args = make_args(
            tmp, protect=cli, missing_check=bool(cli), allow_missing=None
        )
        args = validate_arguments(args, config, False)
    expected = bool(cli) or bool(saved)
    assert args.protect is expected
    assert args.missing_check is expected
    assert args.allow_missing == allow_missing
