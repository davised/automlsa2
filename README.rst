automlsa2
=========

.. image:: https://forthebadge.com/images/badges/built-with-science.svg
    :alt: Made with science
    :target: https://cgrb.oregonstate.edu
.. image:: https://img.shields.io/github/v/release/davised/automlsa2
    :alt: GitHub version
    :target: https://github.com/davised/automlsa2/releases/latest
.. image:: https://img.shields.io/pypi/v/automlsa2
    :alt: PyPI package version
    :target: https://pypi.org/project/automlsa2
.. image:: https://img.shields.io/conda/vn/bioconda/automlsa2
    :alt: Bioconda package version
    :target: https://bioconda.github.io/recipes/automlsa2/README.html
.. image:: https://img.shields.io/pypi/pyversions/automlsa2
    :alt: Supported Python versions
    :target: https://pypi.org/project/automlsa2
.. image:: https://img.shields.io/github/release-date/davised/automlsa2
    :alt: GitHub Release Date
    :target: https://github.com/davised/automlsa2/releases
.. image:: https://img.shields.io/badge/Maintained%3F-yes-green.svg
    :alt: Maintained? Yes
    :target: https://github.com/davised/automlsa2/graphs/commit-activity
.. image:: https://img.shields.io/badge/Made%20with-Python-1f425f.svg
    :alt: Made with python
    :target: https://python.org
.. image:: https://github.com/davised/automlsa2/actions/workflows/tests.yml/badge.svg
    :alt: Tests
    :target: https://github.com/davised/automlsa2/actions/workflows/tests.yml
.. image:: https://img.shields.io/badge/License-BSD_3--Clause-blue.svg
    :alt: License: BSD 3-Clause
    :target: https://github.com/davised/automlsa2/blob/main/LICENSE
.. image:: https://img.shields.io/github/last-commit/davised/automlsa2
    :alt: Latest commit
    :target: https://github.com/davised/automlsa2/commits
.. image:: https://img.shields.io/github/commits-since/davised/automlsa2/latest
    :alt: GitHub commits since latest release (by date)
    :target: https://github.com/davised/automlsa2/commits/main
.. image:: https://img.shields.io/pypi/dm/automlsa2
    :alt: Downloads per month
    :target: https://pypistats.org/packages/automlsa2
.. image:: https://img.shields.io/conda/dn/bioconda/automlsa2
    :alt: Bioconda downloads
    :target: https://anaconda.org/bioconda/automlsa2
.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :alt: Black format
    :target: https://github.com/psf/black
.. image:: https://img.shields.io/badge/mypy-checked-blue
    :alt: Checked with mypy
    :target: http://mypy-lang.org/
.. image:: https://img.shields.io/github/issues/davised/automlsa2
    :alt: GitHub issues
    :target: https://github.com/davised/automlsa2/issues
.. image:: https://img.shields.io/github/forks/davised/automlsa2
    :alt: GitHub forks
    :target: https://github.com/davised/automlsa2/network
.. image:: https://img.shields.io/github/stars/davised/automlsa2
    :alt: GitHub stars
    :target: https://github.com/davised/automlsa2/stargazers

.. contents:: **Table of Contents**
    :backlinks: none

Who might want this software?
-----------------------------

The intended audience is scientific researchers, computational biologists, and
bioinformaticians who are interested in exploring phylogenetic and phylogenomic
relationships between genes and organisms. The general idea is to allow the
input of sequence data along with marker genes and output a robust phylogenetic
tree. All dependencies can be installed with conda, and I hope the software is
easy to use.

If you have feature requests or otherwise have ideas about how to make this
software better, please submit an issue with your ideas!

PhD/Masters students and undergraduates are especially encouraged to submit
issues if they are having trouble using this software.

Windows users, please install WSL to make use of this software. Using a Linux
distribution will make your life as a computational researcher significantly
easier.

Installation
------------

automlsa2 runs on Linux and macOS with Python 3.10 or newer. Windows users can
run it under WSL.

Bioconda (recommended)
~~~~~~~~~~~~~~~~~~~~~~

The easiest way to install automlsa2 along with all of its dependencies,
including BLAST+, MAFFT, and IQ-TREE, is from `bioconda
<https://bioconda.github.io>`_ using `conda <https://docs.conda.io>`_ or
`mamba <https://mamba.readthedocs.io>`_:

.. code-block:: bash

    $ conda create -n automlsa2 -c conda-forge -c bioconda automlsa2
    $ conda activate automlsa2
    $ automlsa2 --version

PyPI
~~~~

automlsa2 is also distributed on `PyPI <https://pypi.org/project/automlsa2/>`_,
but pip only installs the Python dependencies. The external programs need to be
installed separately, e.g. with conda:

.. code-block:: bash

    $ conda create -n automlsa2 -c conda-forge -c bioconda python blast mafft iqtree
    $ conda activate automlsa2
    $ python3 -m pip install automlsa2

From source
~~~~~~~~~~~

.. code-block:: bash

    $ git clone https://github.com/davised/automlsa2.git
    $ cd automlsa2
    $ conda env create -f environment.yml
    $ conda activate automlsa2
    $ python3 -m pip install --no-deps .

For development, install in editable mode instead and run the tests:

.. code-block:: bash

    $ python3 -m pip install --no-deps -e .
    $ python3 -m pytest

Upgrading from 0.9 or earlier
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``--install_deps`` option has been removed; install the external
programs with conda as shown above. Programs previously installed by
``--install_deps`` in ``~/.local/external`` are no longer used automatically
and that directory can be deleted. To keep using programs installed in a
specific directory, pass it with ``--external``.

Dependencies
------------

Python modules:

1. pandas
2. numpy
3. biopython
4. rich
5. packaging
6. psutil

See ``pyproject.toml`` or ``environment.yml`` for more info.

External programs:

1. `NCBI BLAST+ >= 2.12.0 <https://blast.ncbi.nlm.nih.gov>`_
2. `mafft >= 7.487 <https://mafft.cbrc.jp/alignment/software/>`_
3. `IQ-TREE >= 2.1.1 <https://iqtree.github.io>`_ (IQ-TREE 3 is preferred;
   ``iqtree3``, ``iqtree2``, and ``iqtree`` are searched for in that order)

Tested with BLAST+ 2.17.0, mafft 7.526, and IQ-TREE 3.1.4.

External programs are found in your ``$PATH``. To use programs installed
elsewhere, give the directory containing them with ``--external``; it is
searched before ``$PATH``.

Just tell me how to run it
--------------------------

.. code-block:: bash

    $ automlsa2 --files Genus_species_1.fna Genus_species_2.fna ... \
      Genus_species_N.fna --query queries.fasta -t THREADS -- runID

Alternatively:

.. code-block:: bash

    $ automlsa2 --dir path/to/genomes --query queries.fasta -t THREADS \
      -- runID


Overview
--------

automlsa2 is a re-imagination of `autoMLSA.pl
<https://github.com/osuchanglab/autoMLSA>`_

The entire codebase has been re-written in python. While the general algorithm
produces similar output, and several steps are shared, there are many
updates and differences between the two programs, which will be covered later.

The general overview can be summarized here:

0. Input is a set of marker genes as queries, and a set of target genome FASTA
   files.
1. BLAST databases are generated for each target genome, and each query gene
   is extracted from the input query FASTA files.
2. BLAST searches are done with the extracted sequences and genomes.
3. Per genome hits are calculated pending the cut-offs, and genomes are
   filtered from the analysis.
4. Sequences are extracted from the BLAST results as unaligned multi-FASTAs.
5. Unaligned sequences are aligned using mafft.
6. A nexus file is generated pointing to all aligned sequences.
7. A phylogenetic tree is generated using the nexus file as input.

BLAST searches are threaded, or, optionally, written to a file to be submitted
to a compute cluster. mafft alignment commands can also be written to a file
for submission to a compute cluster.

Input query files and genome directories are scanned for updates - if
sequences are added, removed, or changed, the analysis is re-done.

Multiple queries targeting the same gene sequence can be used to improve
coverage of disparate gene sequences, e.g. attempting to cover an entire
phylum with multiple reference genomes being used.

Usage
-----

.. code-block:: text

    $ automlsa2 -h
    usage: automlsa2 [-h] [--query QUERY [QUERY ...]] [--files FILES [FILES ...]]
                     [--dir DIR [DIR ...]] [-e EVALUE] [-c COVERAGE] [-i IDENTITY]
                     [-p {blastn,tblastn}] [--config CONFIG] [--missing_check]
                     [-t THREADS] [--dups] [--allow_missing ALLOW_MISSING]
                     [--outgroup OUTGROUP] [--protect]
                     [--checkpoint {validate,preblast,filtering,prealign,postalign,nexus,none}]
                     [--external EXTERNAL] [--mafft MAFFT] [--iqtree IQTREE] [--debug]
                     [--version] [--quiet]
                     runid

    This is a rewrite of autoMLSA.pl. Generates automated multi-locus sequence analyses.

    positional arguments:
      runid                 Name of the run directory.

    options:
      -h, --help            show this help message and exit
      --query QUERY [QUERY ...]
                            Path to file with input seq(s).
      --files FILES [FILES ...]
                            Path to the target genome FASTA files.
      --dir DIR [DIR ...]   Path to the target genome directory with FASTA files.
      -e, --evalue EVALUE   E-value cutoff for BLAST searches. [1e-5]
      -c, --coverage COVERAGE
                            Sets the coverage cut-off threshold. [50]
      -i, --identity IDENTITY
                            Sets the identity cut-off threshold. [30]
      -p, --program {blastn,tblastn}
                            Which BLAST program to run. [tblastn]
      --config CONFIG       Path to configuration json file to copy.
      --missing_check       Use this to confirm that settings have been checked when genes
                            are missing.
      -t, --threads THREADS
                            Number of threads to use. [1]
      --dups                Allow for duplicate query names for more sequence coverage
                            across disparate organisms.
      --allow_missing ALLOW_MISSING
                            Allow for N missing genes per genome. [0]
      --outgroup OUTGROUP   Name of outgroup genome (file name or strain) to root on.
                            Separate multiple genomes with commas.
      --protect             Save files from getting overwritten. By default, as input
                            files update, older alignments and trees are deleted.
      --checkpoint {validate,preblast,filtering,prealign,postalign,nexus,none}
                            Name of stage to stop computing on. [none]
      --external EXTERNAL   Directory with external programs to use before those in $PATH.
      --mafft MAFFT         mafft settings [--localpair --maxiterate 1000 --reorder]
      --iqtree IQTREE       iqtree settings [-m MFP -B 1000 --alrt 1000 --msub nuclear
                            --merge rclusterf --keep-ident]
      --debug               Turn on debugging messages.
      --version             show program's version number and exit
      --quiet               Turn off progress messages.

One or more input target genome FASTA files is required, either using
``--files`` or ``--dir``. Genome files must be uncompressed FASTA. Additionally, one or more query FASTA files
containing one or more query gene sequences is necessary for analysis.

By default, protein queries are expected, and nucleotide FASTA sequence is
required for the target genomes. ``tblastn`` is used to target the genome
sequences using the amino acid queries. ``blastn`` is also available, targeting
the genome sequences using nucleotide queries.

Threads will speed things up significantly. BLAST searches are threaded in
python; submitting multiple threads to the blast executable often does not
result in much speed up, so each BLAST search is run with one CPU given.

Query marker genes often come from a well-studied representative of, at most,
the same genus. Intergenera phylogenies should have a representative sequence
from each genus. This can be accomplished by giving all examples of a
particular gene the same name in the reference FASTA file. e.g.

.. code-block:: bash

  >Gene1 Refgenus1 refspecies ABC
  <AA sequence>
  >Gene1 Refgenus2 refspecies DEF
  <AA sequence>
  >Gene1 Refgenus3 refspecies GHI
  <AA sequence>

This ^ FASTA ^ file would have three representatives of Gene1 in the analysis.
The resulting alignments would have one copy of the gene, with the best hits
from each target genome included.

Target genome files will be named based on the filename in the final output.
Generally, one will want to have Genus_species_strain.fasta or
G_species_strain.fasta as the filenames prior to analysis.

Genomes can be downloaded using my ``get_assemblies`` program, here:
https://pypi.org/project/get-assemblies/. Locally produced genomes can be
renamed as required.

To root the tree, give ``--outgroup`` the name of one of the genomes, either
as the strain name used in the tree (the file name without its extension) or
as the genome file name or path, e.g. ``--outgroup Escherichia_coli_K-12``
or ``--outgroup genomes/Escherichia_coli_K-12.fna``. Several genomes can be
given separated by commas.

Missing genes and multiple copies
---------------------------------

After the BLAST searches, hits below the ``--identity`` and ``--coverage``
cut-offs are discarded and a presence matrix of hits per genome and gene is
written to ``presence_matrix.tsv``.

* A genome missing more genes than ``--allow_missing`` (default 0) is removed
  from the analysis. Genomes with no hits at all are reported and removed the
  same way.
* Genes missing from more than half of the genomes are listed as candidates
  for removal from the query file.
* If any genome is missing genes, the run stops so you can review the
  results. Re-run with ``--missing_check`` to continue with the genomes that
  were kept.

When a genome has more than one hit for a gene (e.g. multi-copy genes such as
16S rRNA, or several query sequences sharing one gene name), only the best hit
is used: the one with the highest bitscore, then the highest query coverage.

Output files
------------

All files are written to the run directory (``runID``). The main results are:

* ``runID.nex.treefile`` - the final maximum-likelihood tree (Newick), with
  SH-aLRT and ultrafast bootstrap support values on each branch
  (``SH-aLRT/UFBoot``).
* ``runID.nex.contree`` - the bootstrap consensus tree.
* ``runID.nex.iqtree`` - the IQ-TREE report, including the selected models
  and partitioning scheme.
* ``runID.nex.splits.nex`` - bootstrap support for splits, e.g. for viewing
  as a network in SplitsTree.
* ``runID.nex.log`` - the IQ-TREE log.
* ``presence_matrix.tsv`` - number of BLAST hits per genome and gene.
* ``blast_summary.json``, ``missing_by_genome.json``, ``missing_counts.json``
  - summaries of missing genes per genome and per gene.
* ``aligned/`` - the per-gene alignments, and ``runID.nex``, the partition
  file pointing to them.
* ``runID.log`` - the automlsa2 log.

Intermediate files are kept in ``fasta/`` (renamed genomes and BLAST
databases), ``queries/``, ``blast/``, ``unaligned/``, and ``.autoMLSA/``.

Re-running and resuming
-----------------------

Running the same command again with the same ``runID`` resumes the analysis,
skipping steps that are already done. This also applies after stopping a run
with Ctrl-C or with ``--checkpoint``. Settings are saved in
``runID/config.json`` and reused, so later runs only need the options you want
to change.

When inputs or settings change, the affected steps are redone automatically:

* Adding, removing, or changing genomes or query sequences redoes the
  affected BLAST searches and everything after them.
* Changing ``--evalue`` or ``--program`` redoes all BLAST searches.
* Changing ``--identity``, ``--coverage``, ``--allow_missing``, ``--mafft``,
  or ``--iqtree`` redoes the filtering, alignments, and tree.

Previous nexus and tree files are moved to a timestamped folder under
``backup/``. Use ``--protect`` to stop instead of deleting or moving any
results.

Running Tests
-------------

Unit and property-based tests (using pytest and hypothesis) do not need the
external programs and run in a few seconds:

.. code-block:: bash

    $ python3 -m pip install -e '.[test]'
    $ python3 -m pytest

An end-to-end test runs the full pipeline on example data, which is included
in a separate repository - https://github.com/davised/automlsa2-examples

To download this test data, you can run ``git submodule update --init``.

Alternatively, you can run the ``run_tests.sh`` script and it will download
the git submodule for you and run the test command.

The end-to-end test takes about 2-3 minutes with 4 CPUs.

TODO
----

☑ Write list of output files.

☐ Compare functionality of this version to prior autoMLSA.pl version.

☑ Check for version numbers for external programs.

Contributing
------------

Bug reports are encouraged! Submit a github issue and I'll be happy to take
a look. Also, feel free to fork and submit pull requests.

Author Contact
--------------

`Ed Davis <mailto:ed@cgrb.oregonstate.edu>`_

Acknowledgments
----------------

Special thanks for helping me test the software and get the python code packaged:

* `Alex Weisberg <https://github.com/alexweisberg>`_
* `Shawn O'Neil <https://github.com/oneilsh>`_

Also, thanks to these groups for supporting me through my scientific career:

* `OSU Chang Lab <https://github.com/osuchanglab>`_
* `Center for Genome Research and Biocomputing @ OSU <https://cgrb.oregonstate.edu>`_

License
-------

automlsa2 is distributed under the BSD 3-Clause License. See the ``LICENSE``
file for details.

Copyrights
----------

Copyright (c) 2026, Edward Davis, Oregon State University
