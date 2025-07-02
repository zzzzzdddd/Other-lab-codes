import nibabel as nib
import numpy as np
from scipy import ndimage
import pandas as pd
import os

def process_nifti_file(nifti_file_path, binary_label_path):
    img = nib.load(nifti_file_path)
    data = img.get_fdata()

    label_img = nib.load(binary_label_path)
    label_data = label_img.get_fdata()

    binary_data = data > 0.9

    # Remove clusters of size less than 5
    labeled_array, num_features = ndimage.label(binary_data, structure=ndimage.generate_binary_structure(3, 1))
    cluster_sizes = np.bincount(labeled_array.ravel())
    small_clusters = cluster_sizes < 30
    remove_small_clusters = small_clusters[labeled_array]
    binary_data_filtered = binary_data.copy()
    binary_data_filtered[remove_small_clusters] = 0

    # Count non-overlapping clusters
    overlap = binary_data_filtered & label_data.astype(bool)
    non_overlapping_clusters = binary_data_filtered & ~overlap
    labeled_array, num_features = ndimage.label(non_overlapping_clusters, structure=ndimage.generate_binary_structure(3, 1))
    unique_clusters = np.unique(labeled_array)
    num_non_overlapping_clusters = len(unique_clusters) - 1  # Exclude background

    has_overlap = np.any(overlap)

    return has_overlap, num_non_overlapping_clusters, binary_data_filtered

def list_files(directory):
    file_paths = []  # List to store the file paths
    for root, dirs, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            file_paths.append(file_path)
    return file_paths


def process_based_on_prob(nifti_file_path, binary_label_path):
    img = nib.load(nifti_file_path)
    data = img.get_fdata()

    label_img = nib.load(binary_label_path)
    label_data = label_img.get_fdata()

    structure = np.array([[[0, 0, 0], [0, 1, 0], [0, 0, 0]],
                          [[0, 1, 0], [1, 1, 1], [0, 1, 0]],
                          [[0, 0, 0], [0, 1, 0], [0, 0, 0]]], dtype=bool)
    structure = ndimage.generate_binary_structure(3, 1)
    labeled_array, num_features = ndimage.label(data > 0, structure=structure)
    cluster_means = ndimage.mean(data, labels=labeled_array, index=range(1, num_features + 1))
    mask = np.isin(labeled_array, np.where(cluster_means < 0.95, 0, range(1, num_features + 1)))
    filtered_data = data * mask

    # labeled_array, num_features = ndimage.label(filtered_data , structure=ndimage.generate_binary_structure(3, 1))
    # cluster_sizes = np.bincount(labeled_array.ravel())
    # small_clusters = cluster_sizes < 10
    # remove_small_clusters = small_clusters[labeled_array]
    # binary_data_filtered = filtered_data.copy()
    # binary_data_filtered[remove_small_clusters] = 0

    binary_data = filtered_data > 0.1

    overlap = binary_data & label_data.astype(bool)
    non_overlapping_clusters = binary_data & ~overlap
    labeled_array, num_features = ndimage.label(non_overlapping_clusters, structure=ndimage.generate_binary_structure(3, 1))
    unique_clusters = np.unique(labeled_array)
    num_non_overlapping_clusters = len(unique_clusters) - 1  # Exclude background
    has_overlap = np.any(overlap)

    return has_overlap, num_non_overlapping_clusters, binary_data


nifti_file_paths = list_files('C:/voxeloutput/voxeloutput')
binary_label_paths = list_files('C:/voxeloutput/voxellabel')
img = nib.load('C:/voxeloutput/header.nii')

results = []
ol = 0
FPperP = 0
for nifti_file_path, binary_label_path in zip(nifti_file_paths, binary_label_paths):
    has_overlap, num_non_overlapping_clusters, newimg = process_nifti_file(nifti_file_path, binary_label_path)
    # has_overlap, num_non_overlapping_clusters, newimg = process_based_on_prob(nifti_file_path, binary_label_path)
    results.append([nifti_file_path, int(has_overlap), num_non_overlapping_clusters])
    if has_overlap:
        ol = ol+1
    FPperP += num_non_overlapping_clusters

    newimg = nib.Nifti1Image(newimg, img.affine, img.header)
    fname = 'C:/voxeloutput/postprocessed/' + nifti_file_path[-18:-9] + '_processed.nii'
    nib.save(newimg, fname)

summary_df = pd.DataFrame(results, columns=["File", "Overlap", "# of FP Clusters"])
summary_df.to_excel('C:/voxeloutput/summary.xlsx', index=False)
print(summary_df)
print("Total overlapped:\n", ol)
print("FP per patient:\n", FPperP/43)

