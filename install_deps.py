#!/usr/bin/env python
from automlsa2.validate_requirements import install_blast, install_mafft,\
    install_iqtree, validate_requirements
from automlsa2.parse_args import init_logger


def main():
    init_logger(debug=False, quiet=False, rundir='', runid='')
    install_blast()
    install_mafft()
    install_iqtree()
    validate_requirements()


if __name__ == '__main__':
    main()
