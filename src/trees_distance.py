import sys
import glob
import os
import json
import dendropy
from dendropy.calculate import treecompare


def compare_trees(ref_tree, inferred_tree, is_bipartitions_updated=False):
    """
    compare_trees is a function that computes the simmetric difference between
    the reference tree and the inferred tree.
    :param ref_tree: reference tree read using dendropy.
    :param inferred_tree: inferred tree read using dendropy.
    :param is_bipartitions_updated: recalculates bipartitions.
    :return: (int) symmetric distance between the reference and inferred tree.
    """
    sd = treecompare.symmetric_difference(ref_tree, inferred_tree)
    return sd

def compute_distance_between_trees(inferred_tree_file_path, reference_tree_):
    """
    compute_distance_between_trees is a function that reads inferred trees using
    dendropy and computes the distance using the function compare_trees
    :param inferred_tree_file_path: path
    :param reference_tree_: reference tree read using dendropy.
    :return: (int) symmetric distance between the reference and inferred tree.
    """
    with open(inferred_tree_file_path, mode='r') as reduced_noise_tree_file:
            reduced_noise_tree_str = ''.join(list(reduced_noise_tree_file))
            reduced_noise_tree = dendropy.Tree.get_from_string(
                reduced_noise_tree_str,
                schema="newick",
                taxon_namespace=tns
            )
    return compare_trees(reference_tree_, reduced_noise_tree)

if __name__ == '__main__':
    args_in = sys.argv[1:]
    if len(args_in) > 0:
        print('         ')
        sys.exit('WARNING: run_program do not require any input.')
    print('processing data...')

    original_dir = './results/raw_test_data'
    reduced_dir = './results/reduced_test_data'
    result_dir = './results/test_results'
    directories = []
    tns = dendropy.TaxonNamespace()

    for folder in glob.glob(original_dir +'/*'):
        sub_folder_name = folder.split('/')[-1]
        directories.append(sub_folder_name)

    compare_trees_dictionary = dict()

    for directory in directories:
        original_dir_path = os.path.join(original_dir, directory)
        new_folder_path = os.path.join(reduced_dir, directory)

        reference_tree_path = glob.glob(original_dir_path + '/*.tree')[0]
        with open(reference_tree_path, mode='r') as ref_tree_file:
            ref_tree_str = ''.join(list(ref_tree_file))
            reference_tree = dendropy.Tree.get_from_string(
                ref_tree_str,
                schema="newick",
                taxon_namespace=tns
            )
        folder_dict_key_name = original_dir_path.split('/')[-1]

        compare_trees_dictionary[folder_dict_key_name] = dict()

        for alignment_path in glob.glob(original_dir_path + '/*.msl'):
            # computing distance between ref tree and infered trees
            reduced_alignment_name = alignment_path.split('/')[-1]
            reduced_filename_out = os.path.join(new_folder_path, reduced_alignment_name)
            tree_outfile_reduced = reduced_filename_out[:-3] + 'tree'

            alignment_name = alignment_path.split('/')[-1]
            filename_out = os.path.join(original_dir_path, alignment_name)
            tree_outfile = filename_out[:-3] + 'tree'

            noise_reduced_distance = compute_distance_between_trees(tree_outfile_reduced, reference_tree)
            original_distance = compute_distance_between_trees(tree_outfile, reference_tree)

            alignment_key_name = alignment_name[:-3]
            compare_trees_dictionary[folder_dict_key_name][alignment_key_name] = (
                original_distance, noise_reduced_distance
            )
    distance_results_path = os.path.join(result_dir, 'distance_result_dict')
    with open(distance_results_path, 'w') as result_dir_file:
        json.dump(compare_trees_dictionary, result_dir_file)

    print('distances file created.')