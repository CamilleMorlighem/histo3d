import pandas as pd
import geopandas as gpd
import pickle


"""
This script contains some functions to write the results of the maps alignment process to files (text file and shapefiles). 
It is used by the script one_to_many_alignment.py. 

Credit: implementation based on the open source code of 
Kai Sun , Yingjie Hu , Jia Song & Yunqiang Zhu (2020): Aligning geographic
entities from historical maps for building knowledge graphs, International Journal of Geographical
Information Science, DOI: 10.1080/13658816.2020.1845702. 
Source code available at https://figshare.com/articles/dataset/Aligning_geographic_entities_from_historical_maps_for_building_knowledge_graphs/13158098
"""

def write_result_file(df_result, method, basepath_plots_to_fp):
    """ Write the matches between building plots and footprints in a text file named by name of matching method

    Keyword arguments:
    df_result -- geodf with matches 
    method -- name of the similarity metric (approximately within relationship)
    basepath_plots_to_fp -- path of the folder in which to save the file 

    Returns:
    path of the text file generated 
    """

    file = basepath_plots_to_fp + method + '.txt'
    df_result.to_csv(file, header=0, index=0, sep='\t')

    return file
    

#DEPRECATED FUNCTION  
def eval_perf(matched_result_path, ground_truth_path):
    """ Deprecated function (not used in the existing implementation)
    Can be used if ground truth about the matches is available to compute some metrics and assess the maps alignment 
    """

    df_matched = pd.read_csv(matched_result_path, header=None, sep='\t')
    df_ground_truth = pd.read_csv(ground_truth_path, header=None, sep='\t')
    intersection = pd.merge(df_matched, df_ground_truth, how='inner')
    num_positive = len(intersection)

    precision = num_positive/len(df_matched)
    recall = num_positive/(len(df_ground_truth))
    f1 = 2.0 * ((precision * recall) / (precision + recall))

    return precision, recall, f1


def write_result_shapefile(matched_result, epsg_code, basepath_plots_to_fp, year):
    """ Write the matches to two shapefiles; one with the matched building plot and one with the matched building footprint 

    Keyword arguments:
    matched_result -- geodf with matches 
    epsg_code -- epsg code of the crs of the entities (int or str) 
    basepath_plots_to_fp -- path of the folder in which to save the files

    Returns:
    nothing
    """

    matched_result.columns = ['sou_id', 'tar_id', "sou_feature", "tar_feature", "bouwjaar", "identifica"]
    df_results=matched_result

    gpd_map1=gpd.GeoDataFrame(df_results[['sou_id', 'tar_id','bouwjaar', "identifica"]], geometry=df_results.sou_feature)
    gpd_map1.columns = ['gid', 'match', 'bouwjaar', "identifica", 'geometry']
    df_results[['tar_id', 'sou_id']]=df_results[['sou_id', 'tar_id']]
    gpd_map2=gpd.GeoDataFrame(df_results[['sou_id', 'tar_id']], geometry=df_results.tar_feature)
    gpd_map2.columns = ['id', 'match', 'geometry']

    gpd_map1.set_crs(epsg=int(epsg_code))
    gpd_map2.set_crs(epsg=int(epsg_code))

    date_map2=year
    date_map1="current_dataset"

    gpd_map1.to_file(basepath_plots_to_fp + 'aligned_entities_' + date_map1 + '.shp')
    gpd_map2.to_file(basepath_plots_to_fp + 'aligned_entities_' + date_map2 + '.shp')



 

