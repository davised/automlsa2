automlsa2
=========

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black
.. image:: https://badge.fury.io/py/automlsa2.svg
    :target: https://badge.fury.io/py/automlsa2
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
tree. I've implemented commands to help with the installation of external
dependencies, and I hope the software is easy to use.

If you have feature requests or otherwise have ideas about how to make this
software better, please submit an issue with your ideas!

PhD/Masters students and undergraduates are especially encouraged to submit
issues if they are having trouble using this software.

Windows users, please install WSL to make use of this software. Using a Linux
distribution will make your life as a computational researcher significantly
easier.

Installation
------------

automlsa2 is distributed on `PyPI <https://pypi.org/project/automlsa2/>`_ as
a universal wheel and is available on Linux/macOS and Windows (untested) and
supports Python 3.7+ and PyPy.

.. code-block:: bash

    $ python3 -m pip install -U automlsa2

git version install:

.. code-block:: bash

    $ git clone https://github.com/davised/automlsa2.git
    $ cd automlsa2
    $ python3 -m pip install -r requirements.txt
    $ python3 -m pip install -U .

Dependencies
------------

Python modules:

1. pandas
2. numpy
3. biopython
4. tqdm

See ``requirements.txt`` for more info.

External programs:

1. `NCBI BLAST+ >= 2.10.1 <https://blast.ncbi.nlm.nih.gov>`_
2. `mafft >= 7.471 <https://mafft.cbrc.jp/alignment/software/>`_
3. `IQ-TREE COVID-19 release >= 2.1.1 <http://www.iqtree.org>`_

You can install external programs using the ``automlsa2 --install_deps``
command. These will be installed to ``${HOME}/.local/external`` unless
otherwise specified.

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

.. code-block:: bash

    $ automlsa2 -h
    usage: automlsa2 [-h] [--query QUERY [QUERY ...]] [--files FILES [FILES ...]]
                 [--dir DIR [DIR ...]] [-e EVALUE] [-c COVERAGE] [-i IDENTITY]
                 [-p {blastn,tblastn}] [--config CONFIG] [--missing_check]
                 [-t THREADS] [--dups] [--allow_missing ALLOW_MISSING]
                 [--outgroup OUTGROUP] [--protect]
                 [--checkpoint {validate,preblast,filtering,prealign,postalign,nexus,none}]
                 [--install_deps [INSTALL_DEPS]] [--external EXTERNAL]
                 [--debug] [--version] [--quiet]
                 runid

    This is a rewrite of autoMLSA.pl. Generates automated multi-locus sequence analyses.

    positional arguments:
      runid                 Name of the run directory.

optional arguments:

-h, --help            show this help message and exit
--query <QUERY [QUERY ...]>
                      Path to file with input seq(s).
--files <FILES [FILES ...]>
                      Path to the target genome FASTA files.
--dir <DIR [DIR ...]>
                      Path to the target genome directory with FASTA files.
-e EVALUE, --evalue EVALUE
                      E-value cutoff for BLAST searches. [1e-5]
-c COVERAGE, --coverage COVERAGE
                      Sets the coverage cut-off threshold. [50]
-i IDENTITY, --identity IDENTITY
                      Sets the identity cut-off threshold. [30]
-p PROGRAM, --program PROGRAM
                      Which BLAST program to run. [tblastn]
                      {tblastn, blastn}
--config CONFIG       Path to configuration json file to copy.
--missing_check       Use this to confirm that settings have been checked when
                      genes are missing.
-t THREADS, --threads THREADS
                      Number of threads to use. [1]
--dups                Allow for duplicate query names for more sequence
                      coverage across disparate organisms.
--allow_missing ALLOW_MISSING
                      Allow for N missing genes per genome. [0]
--outgroup OUTGROUP   Name of outgroup file or strain to root on.
--protect             Save files from getting overwritten. By default, as input
                      files update, older alignments and trees are deleted.
--checkpoint CHECKPOINT
                      Name of stage to stop computing on. [none]
                      {validate,preblast,filtering,prealign,postalign,nexus,none}
--install_deps <[INSTALL_DEPS]>
                      Install dependencies into given directory. [~/.local/external]
--external EXTERNAL   Path to installed external programs. [~/.local/external]
--debug               Turn on debugging messages.
--version             show program's version number and exit
--quiet               Turn off progress messages.

One or more input target genome FASTA files is required, either using
``--files`` or ``--dir``. Additionally, one or more query FASTA files
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

TODO
----

1. Write detailed list of intermediate files.
2. Compare functionality of this version to prior autoMLSA.pl version.
3. Check for version numbers for external programs.

Contributing
------------

Bug reports are encouraged! Submit a github issue and I'll be happy to take
a look. Also, feel free to clone and submit merge requests.

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

automlsa2 is distributed under the terms listed in the ``LICENSE`` file. The
software is free for non-commercial use.

Copyrights
----------

Copyright (c) 2020 Oregon State University

All Rights Reserved.
