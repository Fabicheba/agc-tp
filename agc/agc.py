#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""
# Modules

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
#import nwalign3 as nw

__author__ = "Fatoumata Binta BARRY"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["fatoumata Binta BARRY"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Fatoumata Binta BARRY"
__email__ = "fatoumatabinta.barry@outlook.fr"
__status__ = "Developpement"


# Fonctions


def isfile(path):
    """verifie que l'existence du fichier
      :Parameters:
          path: chemin du fichier
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Recupère les arguments du programme
      Returns: un objet contenant les arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()


# 1 Dé-duplication en séquence "complète"

def read_fasta(amplicon_file, minseqlen):
    """ prend deux arguments correspondant au fichier fasta et à
    la longueur minimale des séquences et retourne un générateur 
    de séquences de longueur l >= minseqlen: yield sequence """
    prot_dict = {}
    with gzip.open(amplicon_file, "rt") as amplicon:
        prot_dict = {}
        prot_id = ""
        for line in amplicon:
            if line.startswith(">"):
                prot_id = line[1:].split()[0]
                prot_dict[prot_id] = ""
            else:
                prot_dict[prot_id] += line.strip()
        for i in prot_dict:
            sequence = prot_dict[i]
            if len(sequence) >= minseqlen:
                yield sequence


def dereplication_fulllength(amplicon, minseqlen, mincount):
    """ Prend trois arguments correspondant au fichier fasta, 
        la longueur minimale des séquences et leur comptage minimum.
        Elle fait appel au générateur fourni par read_fasta et retourne
        un générateur des séquences uniques ayant une occurrence O>=mincount
        ainsi que leur occurrence.
    """
        
    amplicon_unique = {}
    for seq in read_fasta(amplicon, minseqlen):
        if seq in amplicon_unique:
            amplicon_unique[seq] +=1
        else:
            amplicon_unique[seq] =1
    amplicon_filtre = []
    for seq in amplicon_unique:
        occ = amplicon_unique[seq]
        if occ >= mincount:
            seq_count = [seq, occ]
            amplicon_filtre.append(seq_count)
    amplicon_filtre.sort(key= lambda i : i[1], reverse = True)
    for i in amplicon_filtre:
        yield i







#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()


if __name__ == '__main__':
    #main()
    gi=dereplication_fulllength("essai.gz", 400, 2)
    #print(gi)
     
    print(next(gi), "\n")
    print(next(gi), "\n")
