#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import print_function
import argparse
import os.path
import logging
import ete3  # type: ignore

# from collections import defaultdict
from signal import signal, SIGPIPE, SIGINT, SIG_DFL

signal(SIGPIPE, SIG_DFL)
signal(SIGINT, SIG_DFL)


def extant_file(x):
    """
    'Type' for argparse - checks that file exists but does not open.
    """
    if not os.path.exists(x):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x


def init_logger(args):
    logger = logging.getLogger(__name__)
    ch = logging.StreamHandler()
    if args.debug:
        logger.setLevel(logging.DEBUG)
        ch.setLevel(logging.DEBUG)
    elif args.verbose:
        logger.setLevel(logging.INFO)
        ch.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)
        ch.setLevel(logging.WARNING)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return logger


def run_argparse():
    parser = argparse.ArgumentParser(
        description='Tool to print tree images to file. Handles multiple '
        'bootstrap support values.'
    )
    parser.add_argument('treefile', help='Input newick format file.', type=extant_file)
    parser.add_argument(
        '--fsize', help='Default fontsize. [%(default)i]', type=int, default=10
    )
    parser.add_argument(
        '--ftype', help='Default fonttype. [%(default)s]', type=str, default='Calibri'
    )
    parser.add_argument(
        '--dpi',
        help='dots-per-inch value for output. [%(default)i]',
        type=int,
        default=300,
    )
    parser.add_argument(
        '--width',
        help='width (inches) value for output. [%(default)i]',
        type=int,
        default=5,
    )
    parser.add_argument(
        '--debug', help='Turn on debugging messages.', action='store_true'
    )
    parser.add_argument(
        '--verbose', help='Print verbose progress messages.', action='store_true'
    )
    args = parser.parse_args()
    args.logger = init_logger(args)
    return (args.logger, args.treefile, args.dpi, args.width, args.fsize, args.ftype)


def main():
    logger, treefile, dpi, width, fsize, ftype = run_argparse()
    # Testing
    t = ete3.Tree(treefile, format=1)
    ts = ete3.TreeStyle()
    ts.show_leaf_name = False
    ts.legend.add_face(
        ete3.TextFace('Legend', fsize=int(fsize * 0.9), ftype=ftype, penwidth=6),
        column=1,
    )
    ts.legend.add_face(ete3.TextFace(''), column=0)
    ts.legend.add_face(ete3.CircleFace(3, 'FireBrick'), column=0)
    ts.legend.add_face(
        ete3.TextFace('Well Supported', fsize=int(fsize * 0.8), ftype=ftype), column=1
    )
    ts.legend_position = 4
    # ts.scale = 120
    ts.branch_vertical_margin = 10
    # Convert to +/- support
    for i, node in t.traverse():
        if not node.is_leaf():
            if node.name:
                alrt, boot = node.name.split("/", maxsplit=1)
                node.add_feature('alrt', float(alrt))
                node.add_feature('bootstrap', int(boot))
            else:
                node.add_feature('alrt', float(0.0))
                node.add_feature('bootstrap', int(0))
            if node.alrt > 80 and node.bootstrap > 95:
                node.support = 100
            else:
                node.support = 0

    print('# Original')
    print(t)
    outgroup = t.get_midpoint_outgroup()
    t.set_outgroup(outgroup)
    print('# Midpoint Rooted')
    print(t)

    for i, node in enumerate(t.traverse()):
        nstyle = ete3.NodeStyle()
        if i == 0:
            nstyle['size'] = 0
        elif not node.is_leaf():
            nstyle['size'] = 6
            if node.support == 100:
                nstyle['fgcolor'] = 'FireBrick'
            else:
                nstyle['fgcolor'] = 'black'
        else:
            node.add_face(ete3.AttrFace('name', ftype=ftype, fsize=int(fsize)), 0)
            nstyle['size'] = 0
        node.set_style(nstyle)

    t.render(treefile + '.png', dpi=dpi, w=width, units='in', tree_style=ts)


if __name__ == '__main__':
    main()
