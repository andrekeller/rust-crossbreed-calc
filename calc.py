#!/usr/bin/env python

"""Tool to determine possible crossbreed candidates for RUST"""

from collections import defaultdict
import itertools
import operator
from pprint import pprint
import click


@click.command()
@click.option('-x', '--exclude', multiple=True)
@click.option('-X', '--default_exclude', multiple=True, default=['w', 'x'])
@click.option('-p', '--possibilities', default=1, type=click.IntRange(1, 2))
@click.argument('geneseq_file', type=click.File('r'))
def plants(possibilities, geneseq_file, exclude, default_exclude):
    """get resulting plant possibilities for list of genesequences"""

    # read list of genes from file and duplicate them, as each gene can be
    # used up to two times to crossbreed plants
    geneseqs_raw = list(set(geneseq_file.read().upper().splitlines()))
    geneseqs = list(itertools.chain(*zip(geneseqs_raw, geneseqs_raw)))

    if exclude is not None:
        default_exclude = default_exclude + exclude

    results = {}
    for geneseq_combination in itertools.combinations(geneseqs, 4):
        newplants = crossbreed(geneseq_combination)
        if len(newplants) <= possibilities:
            for newplant in newplants:
                if not any(x.upper() in newplant for x in default_exclude):
                    results["".join(newplant)] = sorted(geneseq_combination)

    pprint(sorted(results.items(), key=operator.itemgetter(0)))


def crossbreed(geneseqs):
    """calculate crossbreed result of 4 genesequences"""
    results = []

    genespots = defaultdict(list)
    for index in range(0, 6):
        for gene in geneseqs:
            genespots[index].append(gene[index])

    for spot in genespots.values():
        max_score = 0
        possibilities = []
        for gene, score in weight(spot):
            # find the most likely genes per spot. If two genes have the
            # heighest weight, report them both.
            if max_score in (0, score):
                max_score = score
                possibilities.append(gene)
        results.append(possibilities)

    # generate cartisian product of all possible combinations given the gene
    # spot possibilities
    return list(itertools.product(*results))


def weight(geneseq):
    """calculate weight of individual genes"""
    values = {
        'X': 0.8,
        'W': 0.8,
        'Y': 0.5,
        'G': 0.5,
        'H': 0.5,
    }
    result = defaultdict(float)

    for gene in geneseq:
        result[gene] += values[gene]

    return sorted(result.items(), key=operator.itemgetter(1), reverse=True)

if __name__ == '__main__':
    # parameters are added by the click.command decorator
    plants()  # pylint: disable=no-value-for-parameter
