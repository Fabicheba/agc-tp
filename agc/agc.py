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
import nwalign3 as nw

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
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile,
                        required=True, help="Amplicon is a compressed fasta file (.fasta.gz)")
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


def get_chunks(sequence, chunk_size):
    """ prend une séquence et une longueur de segment l: chunk_size et
        retourne une liste de sous-séquences de taille l non chevauchantes.
        A minima 4 segments par séquence. """
    sous_sequences = []
    for i in range(0,len(sequence),chunk_size):
        sous_sequence = sequence[i:i+chunk_size]
        sous_sequences.append(sous_sequence)
    if len(sous_sequences) >= 4:
        return sous_sequences


def cut_kmer(sequence, kmer_size):
    """ prend une séquence et une longueur de séquence k et retourne un
        générateur de tous les mots de longueur k présents dans cette séquence """
    for i in range(len(sequence) - kmer_size + 1):
        yield sequence[i:i+kmer_size]


def get_unique_kmer(kmer_dict, sequence, id_sequence, kmer_size):
    """ prend un dictionnaire ayant pour clé un index de kmer et pour valeur
        une liste d’identifiant des séquences dont ils proviennen """
    kmer_dict = {}
    for kmer in cut_kmer(sequence, kmer_size):
        if kmer in kmer_dict:
            kmer_dict[kmer].append(id_sequence)
        else:
            kmer_dict[kmer]=[id_sequence]
    return kmer_dict


def search_mates(kmer_dict, sequence, kmer_size):
    """ prend un dictionnaire ayant pour clé un index de kmer et pour valeur une liste
        d’identifiant des séquences dont ils proviennent, une séquence et une longueur
        de kmer: kmer_size. """
    return [i[0] for i in Counter([ids for kmer in cut_kmer(sequence, kmer_size)
            if kmer in kmer_dict for ids in kmer_dict[kmer]]).most_common(8)]


def get_identity(alignment):
    """ prend un alignement (sous forme de liste) et calcule le pourcentage
        d’identité entre les deux séquences"""
    taille = len(alignment[0]) # longueur de l'alignement
    seq_1 = alignment[0]
    seq_2 = alignment[1]
    ident = 0
    for i in range(taille):
        if seq_1[1] == seq_2[i]:
            ident +=1
    pourcent_ident = round(ident *100/ taille, 2)
    return pourcent_ident


def get_unique(ids):
    """ """
    return {}.fromkeys(ids).keys()


def detect_chimera(perc_identity_matrix):
    """ prend une matrice donnant par segment le taux d’identité entre
        la séquence candidate et deux séquences parentes et retourne un
        booléen indiquant si la séquence candidate est une chimère (True)
        ou ne l’est pas (False) """
    sd = []
    bool_1 = False
    bool_2 = False

    for i in range(len(perc_identity_matrix)):
        perc_1 = perc_identity_matrix[i][0]
        perc_2 = perc_identity_matrix[i][1]
        sd.append(statistics.stdev(perc_identity_matrix[i]))
        if perc_1 > perc_2:
            bool_1 = True
        if perc_1 < perc_2:
            bool_2 = True

    if statistics.mean(sd) > 5 and bool_1 and bool_2:
        return True

    return False

def common(lst1, lst2):
    return list(set(lst1) & set(lst2))


def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """ Fait appel au générateur fourni par dereplication_fulllength et retourne
        un générateur des séquences non chimérique au format:
        yield [sequence, count] """
    no_chimere = []
    kmer_dict = {}

    for i, amplicon in enumerate(dereplication_fulllength(amplicon_file, minseqlen, mincount)):
        chunks = get_chunks(amplicon[0], chunk_size)

        chunk_mates = []
        for sequence in chunks:
            mates = search_mates(kmer_dict, sequence, kmer_size)
            chunk_mates.append(mates)
        ens = []
        for j in range(len(chunk_mates)):
            ens = common(ens, chunk_mates[j])
        if len(ens) > 1:
            for one in ens[0:2]:
                seq = get_chunks(no_chimere[one], chunk_size)
                mat_identitque = [[]]
                for k, chunk in enumerate(chunks):
                    align = nw.global_align(chunk, seq[k])
                    identite =  get_identity(align)
                    mat_identitque[k].append(identite)
            chimera = detect_chimera(mat_identitque)
        if not chimera:
            kmer_dict = get_unique_kmer(kmer_dict, amplicon[0], i, kmer_size)
            yield amplicon


def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """ Fait appel à chimera removal et retourne un liste d’OTU, cette liste indiquera pour
        chaque séquence son occurrence """
    OTU = []
    amplicons = chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size)
    for amplicon in amplicons:
        OTU.append(amplicon)

    return OTU


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width]
            for i in range(0, len(text), width))


def write_OTU(OTU_list, output_file):
    """ fichier de sortie des OTU"""
    with open(output_file, "w") as filout:
        for i in range(0,len(OTU_list)):
            filout.write(">OTU_{} occurrence: {}\n{}\n".format(i,
                            OTU_list[i][1], fill(OTU_list[i][0])))


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    OTU_cluster = abundance_greedy_clustering(args.amplicon_file, args.minseqlen,
                  args.mincount, args.chunk_size, args.kmer_size)
    write_OTU(OTU_cluster, args.output_file)


if __name__ == '__main__':
    main()
