from time import time
import numpy as np
import pandas as pd 
from laspy.file import File
from fiona import open as fopen
from fiona.crs import from_epsg
from shapely.geometry import Polygon, Point, point
import geopandas as gpd
from scipy.spatial import cKDTree
import random 


"""
This script assigns reference ground and roof heights to building footprints 
It is used by the script main_whole_workflow.py 
"""


##### The first 4 functions are used to compute a GROUND reference height #####
def read_pointcloud(file_name):
    """Read a point cloud file and store the x,y,z locations in a numpy array

    Keyword arguments:
    file_name -- Filename of the .las file

    Returns:
    array of the points [vector] [[x,y,z], ..., [x,y,z]]
    """

    if file_name[-3:]=="las": 
        with File(file_name, mode = "r") as file:
            # numpy matrix contains: x, y, z
            pts = np.vstack((file.x, file.y, file.z)).transpose()  

    elif file_name[-3:]=="shp": 
        points_shp = gpd.read_file(file_name)
        pts = np.vstack((points_shp["x"], points_shp["y"], points_shp["z"])).transpose()  

    return pts #, header


def cluster_pointcloud(pts, cluster_size):
    """clusters points based on x,y location 

    Keyword arguments:
    pts -- all ground points from .las file in a [n,3] 2D array
    cluster_size -- size of the cluster (nber of points in a cluster)

    Returns:
    dictionnary with items as clusters; key as the grid location and value as a list of all ground points in the cluster 
    """
    
    pts_rounded = pts[:,:2] - (pts[:,:2] % cluster_size)

    pts_dict = {}
    # go over all the points
    for i in range(len(pts)):
        # get the key value (the rounded position)
        key = tuple(pts_rounded[i])
        
        # if the key exist, append the item to the list, if not, create one.
        if(key in pts_dict):
            pts_dict[key].append(pts[i])
        else:
            pts_dict[key] = [pts[i]]

    return pts_dict


def find_point_in_polygon(pts, polygons, cluster_size):
    """for each polygon, check which points are inside, compute a serie of percentiles height from the points (z value) and add this info to the df where one row = one poygon 

    Keyword arguments:
    pts --  all ground points from .las file clustered in squares
    polygons -- df with polygons geometries 
    cluster_size -- size of the cluster (nber of points in a cluster)

    Returns:
    none (geodf is expanded)
    """

    #for each polygon 
    for index, row in polygons.iterrows(): 
        if pd.isna(row["ground-0.0"]): 
            geo = row['geometry']
            points_in_polygon_z = []
            buffer_value=1

            #while we don't have any ground points lying into the buffered polygon 
            while points_in_polygon_z==[]: 

                buffered_geo=geo.buffer(buffer_value, cap_style=3, join_style=2)
                # get the bounding box of the polygon
                bbox = np.array(buffered_geo.bounds)
                # convert the bbox to grid locations
                bmin = (bbox[:2] - (bbox[:2] % cluster_size)).astype(np.int)
                bmax = (bbox[2:] - (bbox[2:] % cluster_size) + 2 * cluster_size).astype(np.int)

                # go over the (x,y) grid clusters to find the points on top of the building
                for x in range(bmin[0], bmax[0], cluster_size):
                    for y in range(bmin[1], bmax[1], cluster_size):
                        # verify that there are points in 
                        if((x,y) in pts):
                            for i, point in enumerate(pts[(x,y)]):
                                if(buffered_geo.intersects(Point(point[0], point[1]))):
                                    points_in_polygon_z.append(point[2])
                                    # remove item from pts list to prevent checking it again
                                    del(pts[(x,y)][i])
                
                #if we found ground points at the polygon location 
                if(len(points_in_polygon_z) > 0): 
                    percentile_0= np.percentile(points_in_polygon_z, 0)
                    percentile_10= np.percentile(points_in_polygon_z, 10)
                    percentile_20= np.percentile(points_in_polygon_z, 20)
                    percentile_30= np.percentile(points_in_polygon_z, 30)
                    percentile_40= np.percentile(points_in_polygon_z, 40)
                    percentile_50= np.percentile(points_in_polygon_z, 50)
                    
                    polygons.at[index, "ground-0.0"]=percentile_0
                    polygons.at[index, "ground-0.1"]=percentile_10
                    polygons.at[index, "ground-0.2"]=percentile_20
                    polygons.at[index, "ground-0.3"]=percentile_30
                    polygons.at[index, "ground-0.4"]=percentile_40
                    polygons.at[index, "ground-0.5"]=percentile_50
                    polygons.at[index, "buffer_value"]=buffer_value
                    polygons.at[index, "nr_ground_pts"]=len(points_in_polygon_z)
                buffer_value+=0.5


def get_ground_reference_heights(input_pc, input_fp_df, cluster_size=10):
    """compute ground reference heights for the polygons from a point in polygon procedure with a point cloud 

    Keyword arguments:
    input_pc -- las file containing the pc 
    input_fp_df -- input geodf with the building footprints 
    cluster_size -- size of the cluster (nber of points in a cluster)

    Returns:
    input geodf expanded with reference ground heights 
    """

    pts = read_pointcloud(input_pc)
    clustered_pts = cluster_pointcloud(pts, cluster_size)
    find_point_in_polygon(clustered_pts, input_fp_df, cluster_size)

    return input_fp_df


##### The following functions are used to compute a ROOF reference height #####
def assign_height_using_neighbours(polygons): 
    """ for input polygons, find their neighbours and compute their roof height using the height of their neighbours 

    Keyword arguments:
    polygons -- df with polygons geometries, some of which already have roof height attributes (aligned building footprints)

    Returns:
    input geodf expanded with reference roof heights 
    """

    #create a dict that maps from the polygon centroid to the height attributes of the polygon 
    centroid_roof_height_correspondance=dict()
    centroid_list=[]
    for index, row in polygons.iterrows(): 
        centroid=row["geometry"].centroid
        centroid_roof_height_correspondance[(centroid.x, centroid.y)]=(row["plot_roof-0.25"], row["plot_roof-0.75"], row["plot_roof-0.99"])
        centroid_list.append((centroid.x, centroid.y))
    
    #create a kdtree with the centroids 
    pts=np.array(centroid_list)
    kdtree=cKDTree(pts)
    still_na_values=True 

    #while there are still some polygons with Na values for the height attributes 
    k = 0
    kneighbours = 5
    #to avoid continuous loop for ever, at max iter (100), increase the number of neighbours to consider for roof height
    while still_na_values: 
        still_na_values=False 
        if k > 100: 
            kneighbours = kneighbours + 1

        #for each polygon 
        for index, row in polygons.iterrows(): 
            #if it has Na values for any height attribute 
            if (not pd.isna(row["plot_roof-0.25"])) and (not pd.isna(row["plot_roof-0.75"])) and (not pd.isna(row["plot_roof-0.99"])): 
                continue 

            #find the 5 closest neighbour 
            centro=row["geometry"].centroid
            centro=(centro.x, centro.y)
            dist, indices=kdtree.query(centro, k=kneighbours)
            
            roof_025=[]
            roof_075=[]
            roof_099=[]
            
            #skip the first neighbour which is the centroid itself 
            #gather the height attribute values of the neighbours 
            for i in range(1,len(indices)): 
                neighbour_idx=indices[i]
                neighbour_centroid=pts[neighbour_idx]
                roof_025_h= centroid_roof_height_correspondance[tuple(neighbour_centroid)][0]
                if not pd.isna(roof_025_h): 
                    roof_025.append(roof_025_h)

                roof_075_h=centroid_roof_height_correspondance[tuple(neighbour_centroid)][1]
                if not pd.isna(roof_075_h): 
                    roof_075.append(roof_075_h)

                roof_099_h= centroid_roof_height_correspondance[tuple(neighbour_centroid)][2]
                if not pd.isna(roof_099_h): 
                    roof_099.append(roof_099_h)

            #if the polygon has at least one neighbour with non Na values for the height attributes 
            if roof_025!=[] and roof_075!=[] and roof_099!=[]: 
                #compute the median heights 
                med_roof_025=np.median(roof_025)
                med_roof_075=np.median(roof_075)
                med_roof_099=np.median(roof_099)
                
                #update the polygon height attributes (na values) with the median heights computed
                if pd.isna(row["plot_roof-0.25"]): 
                    polygons.at[index, "plot_roof-0.25"]=med_roof_025
                else: 
                    med_roof_025=row["plot_roof-0.25"]
                if pd.isna(row["plot_roof-0.75"]): 
                    polygons.at[index, "plot_roof-0.75"]=med_roof_075
                else: 
                    med_roof_075=row["plot_roof-0.75"]
                if pd.isna(row["plot_roof-0.99"]): 
                    polygons.at[index, "plot_roof-0.99"]=med_roof_099
                else: 
                    med_roof_099=row["plot_roof-0.99"]
                
                #update the height attributes of the polygon in the centroid-based dictionary 
                centroid_roof_height_correspondance[(centro[0], centro[1])]=(med_roof_025, med_roof_075, med_roof_099)
            else: 
                still_na_values=True
                k+=1
    
    #return the polygon dataset enriched with height attributes values 
    return polygons


def assign_random_roof_heights(polygons, min_height = 6.7, max_height = 17, min_mansard = 2.6, max_mansard = 4):
    """ for input polygons, assign they random roof heights based on threshold values 

    Keyword arguments:
    polygons -- df with polygons geometries, some of which already have roof height attributes (aligned building footprints)
    min_height -- minimum height betw ground and building at the base of the roof (ie does not include roof)
    max_height -- max height betw ground and building at the base of the roof (ie does not include roof)
    min_mansard -- minimum height of the roof (between roof base and top)
    max_mansard -- max height of the roof (between roof base and top)
    Note the default values are regulated values for Brussels (ref: de Pange, I. (2004). Histoire de l’architecture `a Saint-Gilles. 
                                                                Inventaire du Patrimoine architectural de la R´egion de Bruxelles-Capitale.)

    Returns:
    input geodf expanded with reference roof heights 
    """

    polygons["h_building"]=np.nan
    polygons["h_attic"]=np.nan
    polygons["roof-0.25"]=np.nan
    polygons["roof-0.99"]=np.nan
    for index, row in polygons.iterrows(): 
        polygons.at[index, "h_building"]=round(random.uniform(min_height, max_height), 1)
        polygons.at[index, "h_attic"]=round(random.uniform(min_mansard, max_mansard),1)
        polygons.at[index, "roof-0.25"] = polygons.at[index, "ground-0.2"] + polygons.at[index, "h_building"]
        polygons.at[index, "roof-0.99"] = polygons.at[index, "ground-0.2"] + polygons.at[index, "h_attic"] + polygons.at[index, "h_building"]

    return polygons
       
