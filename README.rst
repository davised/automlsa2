automlsa2
=========

.. contents:: **Table of Contents**
    :backlinks: none

Installation
------------

automlsa2 is distributed on `PyPI <https://pypi.org>`_ as a universal
wheel and is available on Linux/macOS and Windows and supports
Python 3.7+ and PyPy.

.. code-block:: bash

    $ pip install --upgrade automlsa2

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
command. These will be installed to ``${HOME}/.local/external``

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

License
-------

automlsa2 is distributed under the terms listed in the ``LICENSE`` file.
