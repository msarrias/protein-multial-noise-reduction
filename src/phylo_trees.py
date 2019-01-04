import sys
import glob
import os
import json
from Bio.Seq import Seq
from Bio import SeqIO, AlignIO, Align, SeqRecord
from collections import Counter
import dendropy
from dendropy.calculate import treecompare
from Bio import Phylo
import subprocess
from subprocess import PIPE, Popen

def computing_and_writing_alignment_tree(msa_filename, tree_outfile_):
    """
    computing_and_writing_alignment_tree is a function that pipes two processes,
    it first creates a distance matrix for each msa using fastprot,
    an then infers an alignment tree using fnj. The function writes the inferred
    tree in a file per msa.
    :param msa_filename: msa path.
    :param tree_outfile_: inferred tree path.
    :return: file containing the inferred tree.
    """
    args= "cat " + msa_filename + " | fastprot | fnj -O newick -o " + tree_outfile_
    child_noise_reduced = subprocess.Popen(args, shell=True)
    child_noise_reduced.wait()



if __name__ == '__main__':
    args_in = sys.argv[1:]
    if len(args_in) > 0:
        print('         ')
        sys.exit('WARNING: run_program do not require any input.')
    print('processing data...')

    original_dir = './results/raw_test_data'
    reduced_dir = './results/reduced_test_data'
    directories = []
    tns = dendropy.TaxonNamespace()

    for folder in glob.glob(original_dir +'/*'):
        sub_folder_name = folder.split('/')[-1]
        directories.append(sub_folder_name)

    for directory in directories:
        original_dir_path = os.path.join(original_dir, directory)
        new_folder_path = os.path.join(reduced_dir, directory)

        for alignment_path in glob.glob(original_dir_path + '/*.msl'):
            reduced_alignment_name = alignment_path.split('/')[-1]
            reduced_filename_out = os.path.join(new_folder_path, reduced_alignment_name)
            # computing and writing noise reduced alignment trees
            tree_outfile_reduced = reduced_filename_out[:-3] + 'tree'
            computing_and_writing_alignment_tree(reduced_filename_out, tree_outfile_reduced)

            # computing and writing original alignment trees
            alignment_name = alignment_path.split('/')[-1]
            filename_out = os.path.join(original_dir_path, alignment_name)
            tree_outfile = filename_out[:-3] + 'tree'
            computing_and_writing_alignment_tree(alignment_path, tree_outfile)

        for tree_path in glob.glob(original_dir_path + '/*.tree'):
            if os.stat(tree_path).st_size == 0:
                tree_file_name = tree_path.split('/')[-3:]
                sys.exit(f'{tree_path}'+ '>> ERROR: empty tree file.')

    print('successfully inferred alignments tree.')