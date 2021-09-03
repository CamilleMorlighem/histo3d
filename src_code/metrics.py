import shapely 
import math 
import geopandas as gpd
import matplotlib.pyplot as plt
import stat
import os.path 
import pandas as pd 


"""
This script is used to compute all the quality metrics of the building plots after they were (i) extracted from the historical maps, (ii) vectorised and (iii) generalised 
It is used by the script building_plots_extraction.py 
"""

def compute_area(row): 
    return round(row.geometry.area,3) 

def n_vertices(row): 
    multi = row.geometry.type.startswith("Multi")
    if multi:
        n = 0
        # iterate over all parts of multigeometry
        for part in row.geometry:
            n += len(part.exterior.coords)
    else:
        n = len(row.geometry.exterior.coords)
    return n 

def perimeter(row): 
    return round(row.geometry.length) 

def compactness(row): 
    return round((4*math.pi*row.geometry.area)/(row.geometry.length**2), 3) 

def elongation(row): 
    mbr=row.geometry.minimum_rotated_rectangle
    side1=shapely.geometry.LineString([mbr.exterior.coords[0], mbr.exterior.coords[1]])
    side2= shapely.geometry.LineString([mbr.exterior.coords[1], mbr.exterior.coords[2]])
    if side1.length < side2.length: 
        return side1.length/side2.length
    else: 
        return side2.length/side1.length 

def roundness(row): 
    convex_perim=row.geometry.convex_hull.length
    return (4*math.pi*row.geometry.area)/(convex_perim**2)

def convexity(row): 
    convex_perim=row.geometry.convex_hull.length
    return convex_perim/row.geometry.length

def solidity(row): 
    convex_area= row.geometry.convex_hull.area
    return row.geometry.area/convex_area

def rectangularity(row): 
    mbr_area=row.geometry.minimum_rotated_rectangle.area 
    return row.geometry.area/mbr_area

def has_hole(row): 
    if list(row.geometry.interiors) == []: 
        return 0 
    else: 
        return 1 
    

def add_metrics(input_df, prefix, csv_file): 
    """Compute a serie of quality metrics for input polygons and create new (or add to existing) csv file these metrics and the polygons ID 

    Keyword arguments:
    input_df --- input dataframe containing the polygons (geometries) for which to compute metrics 
    prefix --- prefix to add to the metric (same metrics computed for polygon after different step, prefix is usually the name of the step)
    csv_file --- path of csv file to be created (to add to)

    Returns:
    The input df enriched with the metrics 
    """
   
    #shorten col name for fitting in output csv file 
    col_prefix=prefix[0:3]
    input_df["area_"+str(col_prefix)]= input_df.apply(compute_area, axis=1)
    input_df["n_vert_"+str(col_prefix)]=input_df.apply(n_vertices, axis=1)
    input_df["perim_"+str(col_prefix)]=input_df.apply(perimeter, axis=1)
    input_df["compac_"+str(col_prefix)]=input_df.apply(compactness, axis=1)
    input_df["elong_"+str(col_prefix)]=input_df.apply(elongation, axis=1)
    input_df["round_"+str(col_prefix)]=input_df.apply(roundness, axis=1)
    input_df["convex_"+str(col_prefix)]=input_df.apply(convexity, axis=1)
    input_df["solid_"+str(col_prefix)]=input_df.apply(solidity, axis=1)
    input_df["rect_"+str(col_prefix)]=input_df.apply(rectangularity, axis=1)
    input_df["hole_"+str(col_prefix)]=input_df.apply(has_hole, axis=1)

    #compute statistics for whole dataset
    d={} 
    d["file"]=[str(prefix)+"_plots"]
    d["avg_area"]=[round(input_df["area_"+str(col_prefix)].mean(), 3)]
    d["med_area"]=[round(input_df["area_"+str(col_prefix)].median(), 3)]
    d["std_area"]=[round(input_df["area_"+str(col_prefix)].std(axis=0, skipna=True), 3)]
    d["min_area"]=[round(input_df["area_"+str(col_prefix)].min(), 3)]
    d["max_area"]=[round(input_df["area_"+str(col_prefix)].max(), 3)]

    d["avg_n_vert"]=[round(input_df["n_vert_"+str(col_prefix)].mean(), 3)]
    d["med_n_vert"]=[round(input_df["n_vert_"+str(col_prefix)].median(), 3)]
    d["std_n_vert"]=[round(input_df["n_vert_"+str(col_prefix)].std(axis=0, skipna=True), 3)]

    d["avg_perim"]=[round(input_df["perim_"+str(col_prefix)].mean(), 3)]
    d["med_perim"]=[round(input_df["perim_"+str(col_prefix)].median(), 3)]
    d["std_perim"]=[round(input_df["perim_"+str(col_prefix)].std(axis=0, skipna=True), 3)]
    
    d["avg_compac"]=[round(input_df["compac_"+str(col_prefix)].mean(), 3)]
    d["med_compac"]=[round(input_df["compac_"+str(col_prefix)].median(), 3)]
    d["std_compac"]=[round(input_df["compac_"+str(col_prefix)].std(axis=0, skipna=True), 3)]
   
    d["avg_elong"]=[round(input_df["elong_"+str(col_prefix)].mean(), 3)] 
    d["med_elong"]=[round(input_df["elong_"+str(col_prefix)].median(), 3)]
    d["std_elong"]=[round(input_df["elong_"+str(col_prefix)].std(axis=0, skipna=True), 3)]

    d["avg_round"]=[round(input_df["round_"+str(col_prefix)].mean(), 3)]
    d["med_round"]=[round(input_df["round_"+str(col_prefix)].median(), 3)]
    d["std_round"]=[round(input_df["round_"+str(col_prefix)].std(axis=0, skipna=True), 3)]

    d["avg_convex"]=[round(input_df["convex_"+str(col_prefix)].mean(), 3)]
    d["med_convex"]=[round(input_df["convex_"+str(col_prefix)].median(), 3)]
    d["std_convex"]=[round(input_df["convex_"+str(col_prefix)].std(axis=0, skipna=True), 3)]
    d["nber_concave_poly"]=[len(input_df.loc[input_df["convex_"+str(col_prefix)]!=1.0])]

    d["avg_solid"]=[round(input_df["solid_"+str(col_prefix)].mean(), 3)]
    d["med_solid"]=[round(input_df["solid_"+str(col_prefix)].median(), 3)]
    d["std_solid"]=[round(input_df["solid_"+str(col_prefix)].std(axis=0, skipna=True), 3)]

    d["avg_rect"]=[round(input_df["rect_"+str(col_prefix)].mean(), 3)]
    d["med_rect"]=[round(input_df["rect_"+str(col_prefix)].median(), 3)]
    d["std_rect"]=[round(input_df["rect_"+str(col_prefix)].std(axis=0, skipna=True), 3)]

    d["nber_poly_with_hole"]=[len(input_df.loc[input_df["hole_"+str(col_prefix)]==1])]
    #count unique FeaID (sometimes two identical ID when we had a multipolygon which was split with function explode
    d["nber_plots"]=[len(input_df["FeaID"].value_counts())]

    statistics_df = pd.DataFrame(data=d)

    if os.path.exists(csv_file): 
        input_csv=pd.read_csv(csv_file)
        output_csv=input_csv.append(statistics_df)
        output_csv.to_csv(csv_file, index=False)

    else: 
        statistics_df.to_csv(csv_file, index=False)
    
    return input_df

