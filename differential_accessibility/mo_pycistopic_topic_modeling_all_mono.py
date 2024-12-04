"""
This script is used to create the pycistopic object

authors: Martijn van der Werf

"""

from pycisTopic.cistopic_class import CistopicObject
from scipy.sparse import csr_matrix
from scipy.io import mmread
import gzip
from sklearn.preprocessing import binarize
import pandas as pd
import os
from pycisTopic.lda_models import run_cgs_models_mallet
from numpy.random import choice
from pycisTopic.lda_models import evaluate_models
import pickle
import numpy as np
import glob

os.chdir('/groups/umcg-franke-scrna/tmp03/users/umcg-mwvanderwerff/Functional_genomics_PhD/Multiome/scenicplus_scripts/')

"""
To build a CistopicObject:

class CistopicObject:
    
    cisTopic data class.

    :class:`CistopicObject` contains the cell by fragment matrices (stored as counts :attr:`fragment_matrix` and as binary accessibility :attr:`binary_matrix`),
    cell metadata :attr:`cell_data`, region metadata :attr:`region_data` and path/s to the fragments file/s :attr:`path_to_fragments`.

    LDA models from :class:`CisTopicLDAModel` can be stored :attr:`selected_model` as well as cell/region projections :attr:`projections` as a dictionary.

    Attributes
    ----------
    fragment_matrix: sparse.csr_matrix
        A matrix containing cell names as column names, regions as row names and fragment counts as values.
    binary_matrix: sparse.csr_matrix
        A matrix containing cell names as column names, regions as row names and whether regions as accessible (0: Not accessible; 1: Accessible) as values.
    cell_names: list
        A list containing cell names.
    region_names: list
        A list containing region names.
    cell_data: pd.DataFrame
        A data frame containing cell information, with cells as indexes and attributes as columns.
    region_data: pd.DataFrame
        A data frame containing region information, with region as indexes and attributes as columns.
    path_to_fragments: str or dict
        A list containing the paths to the fragments files used to generate the :class:`CistopicObject`.
    project: str
        Name of the cisTopic project.
"""

# Load and read matrix
print("Loading sparse matrix..")
fragment_matrix_gzip = open('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/matrices/monocyte/matrix.mtx')
coo_fragment_matrix = mmread(fragment_matrix_gzip)
fragment_matrix = csr_matrix(coo_fragment_matrix)
fragment_matrix_gzip.close()
print("Loaded sparse matrix!")

# Keep regions if 0.001 cells have a region
# Get number of non-zero values per columns
non_zero_per_row = fragment_matrix.getnnz(axis=1)
# Calculate fraction compared to number of cells
non_zero_per_row_fraction = non_zero_per_row / fragment_matrix.shape[1]
# Filter for cells >= 0.001
regions_pass_filter = np.argwhere(non_zero_per_row_fraction >= 0.001).ravel()
# Filter chrom 6 matrix on regions where 0.001 cells have that region
fragment_matrix_filtered = fragment_matrix[regions_pass_filter,:]

# Binarization sparse matrix 
print("Binarization of sparse matrix..")
fragment_matrix_binarized = binarize(fragment_matrix_filtered, threshold=0)
print("Binarization done!")

print("Loading files for cisTopic object..")
# Read region names
region_names_file = gzip.open('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/matrices/monocyte/features.tsv.gz', 'rb')
region_names = pd.read_csv(region_names_file, sep='\t', header=None).iloc[:,0].to_list()
region_names = [i.replace('-', ':', 1) for i in region_names]
# Filter region names that are matching index position with regions_pass_filter
region_names_filtered = [region_names[i] for i in regions_pass_filter]
# Read cell names
cell_names_file = gzip.open("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/matrices/monocyte/barcodes.tsv.gz", 'rb')
cell_names = pd.read_csv(cell_names_file, sep='\t', header=None).iloc[:,0].to_list()

# Read region data
region_data = pd.read_csv("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/signac_peaks/output/mo_peaks_lane1to80_monocyte.bed", sep="\t")
region_data.name = [i.replace('-', ':', 1) for i in region_data.name]
# Filter region_data for regions that are in region_names, dims are not equal
region_data = region_data[region_data['name'].isin(region_names)]
region_data.index = region_data.name
region_data.index.name = None
# Filter region data that are matching index position with regions_pass_filter
region_data_filtered = region_data[region_data.index.isin(region_names_filtered)]

# Fragment count cPeaks path
fragment_path = "/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/matrices/monocyte/"

# Cell metadata
# Change dtypes based on this warning: 
# <stdin>:1: DtypeWarning: Columns (15,20,36,51) have mixed types. Specify dtype option on import or set low_memory=False. 
cell_metadata = pd.read_csv("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/metadata/mo_celllevel_metadata.tsv.gz", sep="\t", dtype={'confined_condition': 'category', 
                                                                                                                                                'unconfined_condition': 'category',
                                                                                                                                                'final_condition': 'category',
                                                                                                                                  'LONG_COVID': 'category'})

# Subset cell metadata 
cell_metadata = cell_metadata[cell_metadata.barcode_lane.isin(cell_names)]
# Filter cell names
cell_names = cell_metadata.barcode_lane.to_list()
# Names as rownames 
cell_metadata.index = cell_metadata.barcode_lane

# Initialize cisTopic object
cistopic_obj = CistopicObject(
    fragment_matrix=fragment_matrix_filtered,
    binary_matrix=fragment_matrix_binarized,
    cell_names=cell_names,
    region_names=region_names_filtered,
    region_data=region_data_filtered,
    cell_data=cell_metadata,
    path_to_fragments=fragment_path,
    project="Monocytes"
)
print("Finished cisTopic object!")
print(cistopic_obj)

with open("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_all_cells/cistopic_objects/cistopic_obj.pkl", 'wb') as f:
   pickle.dump(cistopic_obj, f)

with open("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_all_cells/cistopic_objects/cistopic_obj.pkl", 'rb') as f:
   cistopic_obj = pickle.load(f)
   
# Run mallet for LDA
os.environ['MALLET_MEMORY'] = '300G'

# Mallet path
mallet_path = "/groups/umcg-franke-scrna/tmp03/users/umcg-mwvanderwerff/scenicplus/Mallet-202108/bin/mallet"

print("Runnign Mallet LDA..")
models=run_cgs_models_mallet(
    cistopic_obj,
    n_topics=[10,20],
    n_cpu=10,
    n_iter=500,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    tmp_path=os.environ["TMPDIR"],
    save_path="/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_all_cells/topic_models",
    mallet_path=mallet_path,
)

print("Finished Mallet LDA")

file_list = glob.glob(os.path.join("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_all_cells/topic_models", "Topic*.pkl"))

models = []

for file_path in file_list:
    with open(file_path, "rb") as f_input:
        models.append(pickle.load(f_input))


model = evaluate_models(
    models,
    return_model = True,
    select_model = 10,
    #save='figures_all_mono/model_evaluation.pdf'
)

# Add model to cistopic object 
cistopic_obj.add_LDA_model(model)

# Dump model with 1 topics
with open("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_all_cells/cistopic_objects/cistopic_obj_10_topics_model.pkl", 'wb') as f:
   pickle.dump(cistopic_obj, f)

with open("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_all_cells/cistopic_objects/cistopic_obj_10_topics_model.pkl", 'rb') as f:
   cistopic_obj = pickle.load(f)

# Add correct metadata (if this code is re-ran not needed, because when saving the cistopic object the barcodes were scrambled)
# cistopic_obj.cell_data = cell_metadata

# Clustering and visualization
from pycisTopic.clust_vis import (
    find_clusters,
    run_umap,
    run_tsne,
    plot_metadata,
    plot_topic,
    cell_topic_heatmap
)

# Identify clusters
find_clusters(
    cistopic_obj,
    target  = 'cell',
    k = 10,
    res = [0.6, 1.2, 3],
    prefix = 'pycisTopic_',
    scale = True,
    split_pattern = '_'
)

# Run UMAP
run_umap(
    cistopic_obj,
    target  = 'cell', scale=True)

# Dump object with clusters
with open("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_all_cells/cistopic_objects/cistopic_obj_10_topics_model.pkl", 'wb') as f:
   pickle.dump(cistopic_obj, f)

with open("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_all_cells/cistopic_objects/cistopic_obj_10_topics_model.pkl", 'rb') as f:
   cistopic_obj = pickle.load(f)

cell_ut = cistopic_obj.cell_data.index[cistopic_obj.cell_data['condition_final'] == 'UT'].tolist()
cell_stim = cistopic_obj.cell_data.index[cistopic_obj.cell_data['condition_final'] == '24hCA'].tolist()

# Plot condition
plot_metadata(
    cistopic_obj,
    reduction_name='UMAP',
    variables=['condition_final'],
    target='cell', num_columns=2,
    text_size=10,
    dot_size=1,
    save='figures_all_mono/umap_condition.pdf')

cell_sub_celltypes = cistopic_obj.cell_data.index[cistopic_obj.cell_data['celltype_imputed'] != 'CD14 Mono'].tolist()

# Plot only certain days the day
cistopic_obj.cell_data['day'] = cistopic_obj.cell_data['batch'].str[:6]

# Plot celltype
plot_metadata(
    cistopic_obj,
    reduction_name='UMAP',
    variables=['batch'],
    target='cell', num_columns=2,
    text_size=3,
    dot_size=5,
    selected_features=cell_sub_celltypes,
    save='figures_all_mono/umap_cellsubtypes_batch.pdf')

plot_metadata(
  cistopic_obj,
  reduction_name='UMAP',
  variables=['day'],
  target='cell', num_columns=2,
  text_size=10,
  dot_size=5,
  selected_features=cistopic_obj.cell_data.index[cistopic_obj.cell_data['day'].isin(['230105', '230316'])].tolist(),
  save=f'figures_all_mono/umap_.pdf')

# Plot metadata 
# Plots for visual inspection topics
cell_topics = cistopic_obj.selected_model.cell_topic.T
# Highest contributing topic per cells
highest_topic = cell_topics.apply(lambda x: cell_topics.columns[np.argmax(x)], axis = 1)
cistopic_obj.cell_data['most_contributing_topic'] = highest_topic
# second highest condtributing topic
second_highest_topic = cell_topics.apply(lambda row: row.nlargest(2).index[-1],axis=1)
# Add to cell_data \
cistopic_obj.cell_data['second_most_contributing_topic'] = second_highest_topic

# most contributing topics plots 
data_barplot_stim = cistopic_obj.cell_data[['most_contributing_topic', 'condition_final']]
cross_tab_prop = pd.crosstab(index=data_barplot_stim['most_contributing_topic'], columns=data_barplot_stim['condition_final'], normalize='index')
cross_tab_prop.plot(kind='bar', stacked=True, title='Stimulation status cells per topic proportion').get_figure().savefig('/groups/umcg-franke-scrna/tmp03/users/umcg-mwvanderwerff/Functional_genomics_PhD/Multiome/scenicplus_scripts/figures_all_mono/stim_status_10_topics_all_mono_most_contributing.pdf')

import matplotlib.pyplot as plt

# Plot topics on x-axis and bars in bars
data_barplot_day = cistopic_obj.cell_data[['most_contributing_topic', 'batch']]
data_barplot_day['batch'] = data_barplot_day['batch'].str[:6]
data_barplot_day = data_barplot_day.groupby(by=['most_contributing_topic', 'batch']).size().reset_index(name='counts')
data_barplot_day_crosstab = pd.crosstab(index=data_barplot_day['most_contributing_topic'], columns=data_barplot_day['batch'], normalize='index')
data_barplot_day_crosstab.plot(kind='bar', stacked=True, title='Batch day cells per topic', figsize=(5,5))
plt.legend(bbox_to_anchor=(1.0,1.0))
plt.tight_layout()
plt.savefig('/groups/umcg-franke-scrna/tmp03/users/umcg-mwvanderwerff/Functional_genomics_PhD/Multiome/scenicplus_scripts/figures_all_mono/day_batch_10_topics_all_mono_most_contributing.pdf')

# Plot bar on x-axis and topics in bars
data_barplot_day = cistopic_obj.cell_data[['most_contributing_topic', 'batch']]
data_barplot_day = data_barplot_day.groupby(by=['most_contributing_topic', 'batch']).size().reset_index(name='counts')
data_barplot_day_crosstab = pd.crosstab(index=data_barplot_day['batch'], columns=data_barplot_day['most_contributing_topic'], normalize='index')
data_barplot_day_crosstab.plot(kind='bar', stacked=True, title='Frequency cells of each lane in topics', figsize=(12,6), width=1)
plt.ylabel('Frequency')
plt.legend(bbox_to_anchor=(1.0,1.0))
plt.tight_layout()
plt.savefig('/groups/umcg-franke-scrna/tmp03/users/umcg-mwvanderwerff/Functional_genomics_PhD/Multiome/scenicplus_scripts/figures_all_mono/topic_contribution_each_batch_all_mono_most_contributing.pdf')

# Plot monocytes and color for topic contribution
plot_topic(
    cistopic_obj,
    reduction_name = 'UMAP',
    target = 'cell',
    num_columns=4,
    save='figures_all_mono/umap_monocytes_colored_by_10_topics.pdf'
)

############################################################
# VALIDATION TOPICS 
############################################################

from scipy.stats import mannwhitneyu
# Select topic for each cell where that cell contributes most to
cell_topics = cistopic_obj.selected_model.cell_topic.T
# Highest contributing topic per cells
highest_topic = cell_topics.apply(lambda x: cell_topics.columns[np.argmax(x)], axis = 1)
cistopic_obj.cell_data['most_contributing_topic'] = highest_topic

# Find regions that have the highest score for that topic
region_topic = cistopic_obj.selected_model.region_topic
highest_region = region_topic.apply(lambda x: cell_topics.columns[np.argmax(x)], axis=1)

# Get indices stim/unstim
ut_cells_indices = cistopic_obj.cell_data.reset_index()[cistopic_obj.cell_data.reset_index()['condition_final'] == 'UT'].index.tolist()
stim_cells_indices = cistopic_obj.cell_data.reset_index()[cistopic_obj.cell_data.reset_index()['condition_final'] == '24hCA'].index.tolist()

printed = False
for day in cistopic_obj.cell_data['day'].unique():
  # Get indices day
  day_data_indices = cistopic_obj.cell_data.reset_index()[cistopic_obj.cell_data.reset_index()['day'] == day].index.tolist()
  # Topic 6 most associated with UT cells, so we pick that one for the test. Here we extract the indices of the cells
  topic_6_cell_indices = cistopic_obj.cell_data.reset_index()[cistopic_obj.cell_data.reset_index()['most_contributing_topic'] == 'Topic6'].index.tolist()
  # Get topic 6 cells data 
  topic_6_cell_data = cistopic_obj.cell_data.iloc[topic_6_cell_indices]
  # Get regions contributing most to topic 6
  topic_6_regions_indices = highest_region.reset_index()[highest_region.reset_index().iloc[:,1] == 'Topic6'].index.tolist()
  # Get indices matching for day, ut/stim, topic6
  topic_6_cells_ut_day = list(set(day_data_indices).intersection(ut_cells_indices, topic_6_cell_indices))
  topic_6_cells_stim_day = list(set(day_data_indices).intersection(stim_cells_indices, topic_6_cell_indices))
  if len(topic_6_cells_ut_day)==0 or len(topic_6_cells_stim_day)==0:
     continue
  # Split filtered (based on regions) region-by-cell matrix in stim and unstim
  region_by_cell_topic6_ut_day = fragment_matrix_filtered[:,topic_6_cells_ut_day]
  region_by_cell_topic6_ut_day = region_by_cell_topic6_ut_day[topic_6_regions_indices,:]
  if not printed: 
    print(f'Total number of regions for topic 6: {len(topic_6_regions_indices)}')
    printed = True
  print(f'{day} has {region_by_cell_topic6_ut_day.shape[1]} UT cells for topic 6')
  region_by_cell_topic6_stim_day = fragment_matrix_filtered[:,topic_6_cells_stim_day]
  region_by_cell_topic6_stim_day = region_by_cell_topic6_stim_day[topic_6_regions_indices,:]
  print(f'{day} has {region_by_cell_topic6_stim_day.shape[1]} 24hCa cells for topic 6')
  # Wilcoxon signed rank test between stimulated and unstimulated
  U1, p = mannwhitneyu(x=np.asarray(region_by_cell_topic6_ut_day.todense()).sum(axis=1), y=np.asarray(region_by_cell_topic6_stim_day.todense()).sum(axis=1), alternative='less')
  print(f'p-value: {p}')

# Sum the region values for stim and unstim
summed_region_scores_ut = np.asarray(region_by_cell_ut.sum(axis=1)).ravel()
summed_region_scores_stim = np.asarray(region_by_cell_stim.sum(axis=1)).ravel()

from pycisTopic.topic_binarization import binarize_topics

# Region binarization using Otsu, Li and nTopics models
binarized_otsu = binarize_topics(
    cistopic_obj, method='otsu',
    plot=True, num_columns=4#, save='quality_contol_topics/topics_20/binarized_otsu.pdf'
)

binarized_li = binarize_topics(
    cistopic_obj,
    target='cell',
    method='li',
    plot=True,
    num_columns=5, nbins=100
    #save='quality_contol_topics/topics_20/binarized_li.pdf'
    )

# run ntop for downstream analysis
binarized_3k_ntop = binarize_topics(
    cistopic_obj, method='ntop', ntop = 3_000,
    plot=True, num_columns=5
    #save='quality_contol_topics/topics_20/binarized_ntop_3k.pdf'
)

with open("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_all_cells/cistopic_objects/cistopic_obj_10_topics_model.pkl", 'wb') as f:
   pickle.dump(cistopic_obj, f)

"""
Following, we can compute the topic quality control metrics. These include:

    Number of assignments

    Topic coherence (Mimno et al., 2011): Measures to which extent high scoring regions in the topic are actually co-accessible in the original data. If it is low it indicates that the topic is rather random. The higher, the better is a topic.

    The marginal topic distribution: Indicates how much each topic contributes to the model. The higher, the better is a topic.

    The gini index: Value between 0 and 1, that indicates the specificity of topics (0: General, 1:Specific)

    If topics have been binarized, the number of regions/cells per topic will be added.

"""

from pycisTopic.topic_qc import compute_topic_metrics, plot_topic_qc, topic_annotation
import matplotlib.pyplot as plt
from pycisTopic.utils import fig2img

topic_qc_metrics = compute_topic_metrics(cistopic_obj)

# Create dictionary with all QC figures
fig_dict={}
fig_dict['CoherenceVSAssignments']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Log10_Assignments', var_color='Gini_index', plot=False, return_fig=True, save='figures_all_mono/quality_contol_topics/topics_10/coherence_vs_assignments_monocytes.pdf')
fig_dict['AssignmentsVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Log10_Assignments', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True, save='figures_all_mono/quality_contol_topics/topics_10/assignments_vs_cell_in_bin_monocytes.pdf')
fig_dict['CoherenceVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True, save='figures_all_mono/quality_contol_topics/topics_10/coherence_vs_cell_in_bin_monocytes.pdf')
fig_dict['CoherenceVSRegions_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Regions_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True, save='figures_all_mono/quality_contol_topics/topics_10/coherence_vs_regions_monocytes.pdf')
fig_dict['CoherenceVSMarginal_dist']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Marginal_topic_dist', var_color='Gini_index', plot=False, return_fig=True, save='figures_all_mono/quality_contol_topics/topics_10/coherence_vs_marginal_dist_monocytes.pdf')
fig_dict['CoherenceVSGini_index']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Gini_index', var_color='Gini_index', plot=False, return_fig=True, save='figures_all_mono/quality_contol_topics/topics_10/coherence_vs_gini_index_monocytes.pdf')


topic_annot = topic_annotation(
    cistopic_obj,
    annot_var='celltype_imputed_lowerres',
    binarized_cell_topic=binarized_li,
    general_topic_thr = 0.2
)

# DAR identification
from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features
)

# Impute regions
imputed_acc_obj = impute_accessibility(
    cistopic_obj,
    selected_cells=None,
    selected_regions=None,
    scale_factor=10**6
)

# Normalize object
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)

variable_regions = find_highly_variable_features(
    normalized_imputed_acc_obj,
    min_disp = 0.05,
    min_mean = 0.0125,
    max_mean = 3,
    max_disp = np.inf,
    n_bins=20,
    n_top_features=None,
    plot=True,
    #save='figures_all_mono/variable_regions_monocytes.pdf'
)

# Get the amount of variable regions
len(variable_regions)

# Save variable regions
with open('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/scenic_workflow/monocytes_all/variable_regions_mono.txt', 'w') as f:
    for region in variable_regions:
        f.write(f"{region}\n")

from pycisTopic.clust_vis import plot_imputed_features

plot_imputed_features(
    cistopic_obj,
    reduction_name='UMAP',
    imputed_data=imputed_acc_obj,
    features=[markers_dict[x].index.tolist()[0] for x in ['24hCa', 'UT']],
    scale=False,
    num_columns=4,
    save='figures/DARS_24hCa_UT_20_topics_5000_monocytes.pdf'
)

# Write cistopic object with model to a pickle file
with open("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_10000_cells_20_topics/cistopic_obj.pkl", 'wb') as f:
   pickle.dump(cistopic_obj, f, protocol=5)

with open("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_10000_cells_20_topics/cistopic_obj.pkl", "rb") as f:
    cistopic_obj = pickle.load(f)

# Topic 5 seems to show the most information of CD16 monocytes 

# Get list of all topics (second most contributing)
#topics_dars_compare_with_topic_5 = cistopic_obj.cell_data.second_most_contributing_topic.unique().tolist()
# Remove topic 5 which we want to compare 
#topics_dars_compare_with_topic_5.pop(3)

# Removed topic7 beause of issues in SCENIC+ workflow
#input_query = [[['Topic11'], ['Topic9', 'Topic7', 'Topic4']]]

# Run DAR analysis
markers_dict= find_diff_features(
    cistopic_obj,
    imputed_acc_obj,
    variable='most_contributing_topic',
    var_features=variable_regions,
 #   contrasts=input_query,
    adjpval_thr=0.05,
    log2fc_thr=np.log2(1.5),
    n_cpu=4,
    _temp_dir=os.environ["TMPDIR"],
    split_pattern = '_'
)

from pycisTopic.clust_vis import plot_imputed_features

# Number of DARs
print("Number of DARs found:")
print("---------------------")
for x in markers_dict:
    print(f"  {x}: {len(markers_dict[x])}")


# Get non-empty results (some dicts were emtpy --> no DARs)
markers_dict = {k:v for (k,v) in markers_dict.items() if not v.empty}

from pycisTopic.utils import region_names_to_coordinates

# Save DARs monocytes 
for topic in markers_dict:
    region_names_to_coordinates(
        markers_dict[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_all_cells/", "region_sets/DARS_topics_vs_others", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )

# Save 3k topics, removed topic 18 because of issues
for topic in binarized_3k_ntop:
    region_names_to_coordinates(
        binarized_3k_ntop[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_all_cells/', "region_sets", "Topics_top_3k", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )

for topic in binarized_otsu:
    region_names_to_coordinates(
        binarized_otsu[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_all_cells/', "region_sets", "Topics_otsu", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )

with open("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_10000_cells_20_topics/cistopic_obj.pkl", 'wb') as f:
   pickle.dump(cistopic_obj, f)

# Write sampled cells to file so we can use those for subsetting the Seurat object
with open('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_10000_cells_20_topics/cell_barcodes.txt', 'w') as f:
    for barcode in cistopic_obj.cell_data.index.to_list():
        f.write(f"{barcode}\n")