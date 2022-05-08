import geopandas as gpd
import evaluate_performance_alignment
import pandas as pd
import pickle
from shapely.geometry import MultiPoint
import shapely
import numpy as np
import os
from shapely.geometry import Point
from shapely.geometry import MultiPoint
from shapely.validation import explain_validity



"""
This script performs the one-to-many alignment step between building plots and building footprints.
It is used by the script main_whole_workflow.py. 

Credit: implementation based on the open source code of 
Kai Sun , Yingjie Hu , Jia Song & Yunqiang Zhu (2020): Aligning geographic
entities from historical maps for building knowledge graphs, International Journal of Geographical
Information Science, DOI: 10.1080/13658816.2020.1845702. 
Source code available at https://figshare.com/articles/dataset/Aligning_geographic_entities_from_historical_maps_for_building_knowledge_graphs/13158098
"""


def atr_within(row, radius=10):
    """Compute main metric on which the 1:n alignment is based; approximately within relation between building footprint and buffered plot. It is actually an area ratio. 

    Keyword arguments:
    row -- row of a geodf containing two geometries, the geometry of the building footprint and the one of the building plot 
    radius -- buffer size for the building plot (default 10 m but can be changed)

    Returns:
    area ratio (float)
    """

    geometry_sou = row['sou_feature']
    geometry_tar = row['tar_feature']
    buf_2 = geometry_tar.buffer(radius,cap_style=3, join_style=2)
    buf_1=geometry_sou

    if explain_validity(geometry_sou) == 'Valid Geometry' and explain_validity(geometry_tar) == 'Valid Geometry':
        area = buf_1.intersection(buf_2).area
        area_ratio = area/(min(buf_1.area, buf_2.area))
        return area_ratio


def approx(df_similarity, approx_value=0.6):
    """Determines whether a building plot and a building footprint are a match comparing their approximately within relationship (area ratio) with a threshold 

    Keyword arguments:
    df_similarity -- geodf with building plot, building footprint and their area ratio 
    approx_value -- threshold value to compare the ratio with 

    Returns:
    geodf enriched with matches (one column with building plot, one with building footprint and additional attributes)
    """

    df_result = pd.DataFrame(columns=('sou_id', 'tar_id', "sou_feature", 'tar_feature', 'bouwjaar', 'identifica'))
    # checks if atr_withing is smaller than approx value 
    df_dist_approx = df_similarity[df_similarity['atr_within'] >= approx_value]
    # for each group of plots that match with a footprint, keep the one with highest approx value
    # a footprint can match with several plots because approx_value is computed within the buffered plot 
    for sou_id, group in df_dist_approx.groupby(['sou_id']):
        idx_max=group["atr_within"].argmax()
        group_max=group.iloc[idx_max]
        
        df_result = df_result.append(group_max[['sou_id', 'tar_id', 'sou_feature', 'tar_feature', 'bouwjaar', "identifica"]], ignore_index=True)
    return df_result


def transformation_crs(entity_set2, entity_set_crs):
    """Change crs of geodf to the crs of another geodf 

    Keyword arguments:
    entity_set2 -- geodf of which we want to change the crs 
    entity_set_crs -- geodf with crs we want to transform to 

    Returns:
    transformed geodf 
    """

    entity_set_transformed= entity_set2.to_crs(entity_set_crs)
    return entity_set_transformed



def overlapping_entity_pairs(entity_set1_overlaid, entity_set2_overlaid):
    """Search the entities which are within the overlapping area of two entity sets. Only entities within the overlapping area will be processed further.

    Keyword arguments:
    entity_set1_overlaid -- geodf 1 
    entity_set2_overlaid -- geodf 2 

    Returns:
    subset of geodf 1 and subset of geodf 2 which are located in the intersection area of the two input geodf 
    """
    
    all_vertices1 = []
    all_vertices2 = []

    # Obtain vertices of all entities to be matched.
    for index, row in entity_set1_overlaid.iterrows():
        if row.geometry.geom_type == 'Polygon':
            all_vertices1.append(list(row.geometry.exterior.coords))
        else:
            all_vertices1.append(list(row.geometry.coords))
    all_vertices1 = MultiPoint(sum(all_vertices1, []))

    for index, row in entity_set2_overlaid.iterrows():
        if row.geometry.geom_type == 'Polygon':
            all_vertices2.append(list(row.geometry.exterior.coords))
        else:
            all_vertices2.append(list(row.geometry.coords))
    all_vertices2 = MultiPoint(sum(all_vertices2, []))

    # Compute overlapping area by computing the intersection area of convex_hulls which are generated with all vertices.
    all_vertices1_convexhull = all_vertices1.convex_hull
    all_vertices2_convexhull = all_vertices2.convex_hull
    overlapping_area = all_vertices1_convexhull.intersection(all_vertices2_convexhull)
    links_gdf = gpd.GeoDataFrame(geometry=[overlapping_area])

    # Those entities which do not intersect with the overlapping area will be removed.
    entity_set1_overlaid = entity_set1_overlaid.reset_index(drop = True)
    entity_set1_overlaid_copy = entity_set1_overlaid.copy()
    for index, row in entity_set1_overlaid.iterrows():
        if explain_validity(row.geometry) == 'Valid Geometry':
            if row.geometry.intersection(overlapping_area).is_empty:
                delete_index = entity_set1_overlaid_copy[entity_set1_overlaid_copy['gid'] == row.gid].index.tolist()
                entity_set1_overlaid_copy.drop(entity_set1_overlaid_copy.index[delete_index], inplace=True)
                entity_set1_overlaid_copy.index = range(len(entity_set1_overlaid_copy))

    entity_set2_overlaid_copy = entity_set2_overlaid.copy()
    for index, row in entity_set2_overlaid.iterrows():
        if explain_validity(row.geometry) == 'Valid Geometry':
            if row.geometry.intersection(overlapping_area).is_empty:
                delete_index = entity_set2_overlaid_copy[entity_set2_overlaid_copy['FeaID'] == row.FeaID].index.tolist()
                entity_set2_overlaid_copy.drop(entity_set2_overlaid_copy.index[delete_index], inplace=True)
                entity_set2_overlaid_copy.index = range(len(entity_set2_overlaid_copy))

    return entity_set1_overlaid_copy, entity_set2_overlaid_copy


def execute_alignment(pickl_file, shapefile1, shapefile2, year):
    """Main fction of the maps alignment which (i) read the input shapefiles and transform them to the same crs, (ii) select the entities that are in 
    the intersection of the extent of the shapefiles and (iii) execute the maps alignment by calling the functions computing the main metric (area ratio)

    Keyword arguments:
    pickl_file -- path to pickl_file that will contain the building plots, building footprints and similarity score 
    shapefile 1 -- shapefile with building footprints 
    shapefile 2 -- shapefile with building plots 
    year -- date of the historical map (str or int)

    Returns:
    nothing 
    """

    entity_set1 = gpd.read_file(shapefile1).explode()
    #select footprints based on construction year 
    entity_set1["bouwjaar"]=entity_set1['bouwjaar'].apply(lambda x: int(x[0:4]))
    entity_set1 = entity_set1.loc[entity_set1["bouwjaar"]<= int(year)]

    entity_set2 = gpd.read_file(shapefile2).explode()
    
    if not entity_set1.empty:
        entity_set_crs1 = entity_set1.crs
    
    if not entity_set2.empty:
        entity_set_crs2 = entity_set2.crs

    # Transform the CRS of entity_set1 to the CRS of entity_set2.
    if entity_set_crs1 != entity_set_crs2:
        entity_set1_overlaid = transformation_crs(entity_set1, entity_set_crs2)
    else: 
        entity_set1_overlaid=entity_set1
    entity_set1_overlapping, entity_set2_overlapping = overlapping_entity_pairs(entity_set1_overlaid, entity_set2)
   
    similarity_calculation(pickl_file, entity_set1_overlapping, entity_set2_overlapping)

        
def similarity_calculation(pickl_file, entity_set1_processed, entity_set2_processed):
    """Compute similarity score (area ratio metric) between building plots and building footprints located nearby and write the results to a pickl file 

    Keyword arguments:
    pickl_file -- path to pickl_file that will contain the building plots, building footprints and similarity score 
    entity_set1_processed -- geodf with building footprints 
    entity_set2_processed -- geodf with building plots 

    Returns:
    nothing 
    """
    
    df_similarity = pd.DataFrame(columns=('sou_id', 'tar_id', 'sou_feature', 'tar_feature', 'bouwjaar', "identifica"))
    num = 0

    # for one-to-many mapping: 
    # a plot and a footprint are a POTENTIAL match if the footprint is included into the buffered plot (still need to further check the similarity score)
    for index1, item1 in entity_set1_processed.iterrows():
        for index2, item2 in entity_set2_processed.iterrows():
            buffered_item2=item2.geometry.buffer(3, cap_style=3, join_style=2) 
            if not item1.geometry.intersection(buffered_item2).is_empty : 
                df_similarity.loc[num] = [item1['gid'], item2['FeaID'], item1.geometry, item2.geometry, item1.bouwjaar, item1.identifica]
                num = num + 1
     
    df_similarity['atr_within'] = df_similarity.apply(atr_within, axis=1) #args=(radius, ), axis=1)

    # The computed dataframe of similarity will be written in a pkl file.
    with open(pickl_file, 'wb') as pickle_file:
        pickle.dump(df_similarity, pickle_file)



def alignment_classification(df_similarity_path, epsg_code, basepath_plots_to_fp, year, ground_truth=None):
    """From the pkl file containing the computed similarity scores, this function compares the similarity score with threshold and identifies the matches 

    Keyword arguments:
    df_similarity_path -- path to pickl_file that will contain the building plots, building footprints and similarity score 
    epsg_code -- epsg code of the crs of the entities (int or str) 
    basepath_plots_to_fp -- path of the folder where the ouput files should be created 
    year -- date of the historical map 
    ground_truth -- do not consider; used in deprecated version of this function 

    Returns:
    nothing 
    """

    # Read the pkl file of similarity.
    with open(df_similarity_path, 'rb') as pickle_file:
        df_similarity = pickle.load(pickle_file)
    df_result = approx(df_similarity)
    result_file = evaluate_performance_alignment.write_result_file(df_result, "approx", basepath_plots_to_fp)
    evaluate_performance_alignment.write_result_shapefile(df_result, epsg_code, basepath_plots_to_fp, year)

