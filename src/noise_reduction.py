
import sys
import glob
import os
from Bio.Seq import Seq
from Bio import SeqIO, AlignIO, Align, SeqRecord
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


def reduce_noise(multiple_seq_alignment, _alignment_path):
    """
    reduce_noise is a function that filters all noisy columns.
    :param multiple_seq_alignment: MultipleSeqAlignment.
    :param _alignment_path: msa path.
    :return: a MultipleSeqAlignment.
    """
    noise_free_column_list = []
    new_aligned_list = []
    alignment_filename_ = _alignment_path.split('/')[-1]

    for i in range(multiple_seq_alignment.get_alignment_length()):
        if not noise_filter(multiple_seq_alignment[:, i]):
            noise_free_column_list.append(multiple_seq_alignment[:, i])
    for i in range(len(multiple_seq_alignment)):
        new_aligned_list.append(''.join([x[i] for x in noise_free_column_list]))
    record_list_ = Align.MultipleSeqAlignment([
        SeqRecord.SeqRecord(Seq(new_seq), id=record.id, name=record.name, description=record.description)
        for new_seq, record in zip(new_aligned_list, multiple_seq_alignment)
     ])

    if multiple_seq_alignment.get_alignment_length() == record_list_.get_alignment_length():
        print(f'{alignment_filename_}'+ '>> WARNING: no noise reduction')
    if 0 < record_list_.get_alignment_length()  < multiple_seq_alignment.get_alignment_length() / 2:
        print(f'{alignment_filename_}'+ '>> WARNING: more than 0.5 of the sequence is noise')
    return record_list_

def perform_noise_reduction(alignment_path_):
    """
    perform_noise_reduction is a function that takes the parameter alignment path and applies
    the reduce_noise function.
    :param alignment_path_: msa path.
    :return: a MultipleSeqAlignment.
    """
    with open(alignment_path_, mode='r') as aligned_file:
        my_sequence_recorded = AlignIO.read(aligned_file, 'fasta')
    return reduce_noise(my_sequence_recorded, alignment_path_)


if __name__ == '__main__':
    args_in = sys.argv[1:]
    if len(args_in) > 0:
        print('         ')
        sys.exit('WARNING: run_program do not require any input.')
    print('processing data...')

    original_dir = './results/raw_test_data'
    reduced_dir = './results/reduced_test_data'
    directories = []

    for folder in glob.glob(original_dir +'/*'):
        sub_folder_name = folder.split('/')[-1]
        directories.append(sub_folder_name)

    for directory in directories:
        original_dir_path = os.path.join(original_dir, directory)
        new_folder_path = os.path.join(reduced_dir, directory)
        os.makedirs(new_folder_path, exist_ok=True)

        for alignment_path in glob.glob(original_dir_path + '/*.msl'):
            alignment_filename = alignment_path.split('/')[-1]
            if os.stat(alignment_path).st_size == 0:
                sys.exit(f'{alignment_filename}'+ '>> ERROR: empty msl file.')

            record_list = perform_noise_reduction(alignment_path)
            if record_list.get_alignment_length() == 0:
                sys.exit(f'{alignment_filename}'+ '>> WARNING: all columns are noisy.')

            reduced_filename_out = os.path.join(new_folder_path, alignment_filename)
            AlignIO.write(record_list, reduced_filename_out, "fasta")

    
    

