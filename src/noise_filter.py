#!/usr/bin/env python
# coding: utf-8

from Bio.Seq import Seq, transcribe, translate
from Bio import SeqIO, AlignIO
from collections import Counter

def noise_filter(column):
    """
    noise_filter is a function that determines whether an alignment column is noisy or not.
    :param column: multiple sequence alignment column.
    :return: boolean.
    """
    indel = '-'
    amino_counter = Counter(column)
    if column.count(indel) > len(column) / 2:
        return True
    if len([key for key, value in amino_counter.items() if value == 1 ]) >= len(column) / 2:
        return True
    if max(amino_counter.values()) <= 2:
        return True
    return False

def reduce_noise(multiple_seq_alignment):
    """
    reduce_noise is a function that filters all noisy columns.
    :param multiple_seq_alignment: MultipleSeqAlignment.
    :return: a MultipleSeqAlignment.
    """
    noise_free_column_list = []
    new_aligned_list = []
    for i in range(my_sequence_recorded.get_alignment_length()):
        if not noise_filter(my_sequence_recorded[:, i]):
            noise_free_column_list.append(my_sequence_recorded[:, i])
    for i in range(len(my_sequence_recorded)):
        new_aligned_list.append(''.join([x[i] for x in noise_free_column_list]))
    record_list = Align.MultipleSeqAlignment([
        SeqRecord.SeqRecord(Seq(new_seq), id=record.id, name=record.name, description=record.description)
        for new_seq, record in zip(new_aligned_list, my_sequence_recorded)
     ])


if __name__ == '__main__':

    file_in = sys.argv[1]