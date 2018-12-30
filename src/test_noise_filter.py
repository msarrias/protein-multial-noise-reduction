#!/usr/bin/env python
# coding: utf-8

from Bio.Seq import Seq, transcribe, translate
from Bio import SeqIO, AlignIO
from collections import Counter

def noise_filter(column):
    indel = '-'
    amino_counter = Counter(column)
    if column.count(indel) > len(column) / 2:
        return True
    if len([key for key, value in amino_counter.items() if value == 1 ]) >= len(column) / 2:
        return True
    if max(amino_counter.values()) <= 2:
        return True
    return False     

column_case1 = "A---------" 

column_case2 = "ARTLSFFFFF"

column_case3 = "AARRTTLLSS"

column_case4 = "AFFFFFFFFF" 

noise_filter(column_case1)

noise_filter(column_case2)


noise_filter(column_case3)


noise_filter(column_case4)



