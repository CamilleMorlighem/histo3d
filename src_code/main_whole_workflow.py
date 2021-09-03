## import needed modules 
import time
import os, sys
import csv
import math
import geopandas as gpd
import statistics as stat 
import shapely
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import one_to_many_alignment
import cjio 
from cjio import cityjson
from cjio import subset 

import building_plots_extraction as bfe
import assign_reference_heights
import plots_to_footprints as ptf

"""
input_map= "delft_1961"
grassbin='C:\\\"Program Files\"\\\"GRASS GIS 7.8\"\\grass78.bat'
full_epsg = ['EPSG:28992', "urn:ogc:def:crs:EPSG::7415"] 
#full_epsg = ['EPSG:31370',"urn:ogc:def:crs:EPSG::6190"]  
basepath= "C:\\Users\\Camille\\Documents\\mem_docs\\3_ANALYSIS\\final_workflow\\"
blender_path = '\"C:\\Program Files\\Blender Foundation\\Blender 2.83\\blender.exe\"'
val3dity_path = "C:\\Users\\Camille\\Documents\\val3dity-windows-x64-v2.2.0\\val3dity.exe" 

#require to place this files at the right place so that they can be accessed automatically: 
myfile=basepath + input_map + "\\input_data\\" + input_map +  ".tif"


training_pts = basepath + input_map + "\\input_data\\training_points_" + input_map + ".shp" 
ground_truth = basepath + input_map + "\\input_data\\ground_truth_" + input_map + '.shp'
ground_truth_generalized= basepath + input_map + "\\input_data\\ground_truth_generalized_" + input_map  + ".shp"
bag_footprints= None #basepath + input_map + "\\input_data\\bag3d_delft.shp"
input_pc= basepath + input_map + "\\input_data\\delf_ground.las"
current_lod2_model = basepath + input_map + "\\input_data\\all_merged_delft_cleaned.json"
all_footprints = basepath + input_map + "\\plots_to_fp\\final_all_footprints.shp" 

#input parameters 
facade_depth_range = [12, 16]
facade_len_range = [5, 7]
avg_area = 300
min_area = 10
roof_height_range = [6.7, 17]
mansard_height_range = [2.6, 4]

alpha = 12
distances = 5
angles = 180
"""

#############################################################################################################


def building_footprints_extraction(basepath, grassbin, myepsg, input_map, myfile, ground_truth, training_pts, ground_truth_generalized, alpha, distances, angles, final_plots_file): 
    """Perform building plots extraction from input historical map 

    Keyword arguments:
    basepath -- path of the folder where scripts and input/output files are located 
    grassbin -- GRASS bin environment variable 
    myepsg -- EPSG code of the coordinate reference system 
    input_map -- string representing the name of the historical map (eg delft_1880)
    myfile -- input historical map (.tif file)
    ground_truth -- points shapefile with the centroids of the ground truth building plots 
    training_pts -- shapefile with training points 
    ground_truth_generalized -- shapefile with ground truth building plots manually digitalized 
    alpha -- integer representing the alpha value 
    distances -- distance threshold value in commandeur generalisation 
    angles -- angle threshold value in commandeur generalisation 
    final_plots_file -- output shapefile that will contain building plots after both generalisation methods are implemented

    Returns:
    nothing 
    """

    basepath_bfe = basepath + input_map + "\\building_plots_extraction\\"
    if not os.path.exists(basepath_bfe):
        os.makedirs(basepath_bfe)
    raw_building_plots = basepath_bfe + "raw_building_plots_" + input_map + '.shp'
    raw_building_plots_area_filter = basepath_bfe + "area_filter_raw_building_plots_" + input_map + '.shp'
    metrics_txt = basepath_bfe + "metrics_" + input_map + '.txt'
    metrics_csv = basepath_bfe + "metrics_" + input_map + '.csv'
    alpha_shape_script = basepath + 'src_code\\simplify_map_alpha_shape.R'
    commandeur_script= basepath + "src_code\\commandeur_method\\bin\\Cell_Decomposition.exe"
    out_alpha_shape=basepath_bfe + input_map + '_building_plots_alpha_shape_'
    out_generalized= basepath_bfe + 'generalized_building_plots_' + input_map 
    csv_gen_results = basepath_bfe + "generalization_metrics_" + input_map + '.csv'

    bfe.OBIA_segmentation_grass(grassbin, myepsg, myfile, basepath, input_map, training_pts, raw_building_plots)
    bfe.area_filter_and_assess_classif(raw_building_plots, metrics_csv, metrics_txt, raw_building_plots_area_filter, ground_truth, area_threshold = 150)
    bfe.final_generalization_workflow(raw_building_plots_area_filter, csv_gen_results, alpha_shape_script, input_map, out_alpha_shape, alpha, distances, angles, commandeur_script, out_generalized, final_plots_file, ground_truth_generalized, metrics_csv, metrics_txt) 



def perform_maps_alignment(basepath, bag_footprints, input_map, myepsg, final_plots):
    """Perform maps alignment 

    Keyword arguments:
    basepath -- path of the folder where scripts and input/output files are located 
    bag_footprints -- shapefile with BAG building footprints (if study area not Netherlands, current building footprints dataset)
    input_map -- string representing the name of the historical map (eg delft_1880)
    myepsg -- EPSG code of the coordinate reference system 
    final_plots -- input shapefile that contains the extracted building plots 

    Returns:
    nothing 
    """

    basepath_plots_to_fp = basepath + input_map + "\\plots_to_fp\\"
    if not os.path.exists(basepath_plots_to_fp):
        os.makedirs(basepath_plots_to_fp)

    pickl_file = basepath_plots_to_fp + input_map + ".pkl"
    year=input_map[-4:]
    one_to_many_alignment.execute_alignment(pickl_file, bag_footprints, final_plots, year)
    one_to_many_alignment.alignment_classification(pickl_file, myepsg.split(":")[1], basepath_plots_to_fp, year)



def generate_final_fp_with_height(basepath, bag_footprints_shp, aligned_fp, final_plots, input_pc, facade_len_min, facade_len_max, facade_depth, max_facade_depth, avg_area, min_area, max_perim=math.inf): 
    """This function is used in the complete methodology version (ie including the maps alignment step). 
    It does the following: (i) proceeds to 2D procedural modelling, (ii) creates one shapefile from the maps alignment and 2PM results 
    and (iii) assigns roof and ground height attributes to the building footprints 

    Keyword arguments:
    basepath -- path of the folder where scripts and input/output files are located 
    bag_footprints_shp -- shapefile with BAG building footprints (if study area not Netherlands, current building footprints dataset)
    aligned_fp -- shapefile that contains the aligned building footprints
    final_plots -- input shapefile that contains the extracted building plots
    input_pc --  input las file containing a 3D point cloud 
    #Procedural modelling parameters: 
    avg_area -- average building footprint area 
    min_area -- minimum area of footprint to be considered as valid footprint 
    max_perim -- max perimeter allowed for footprint to be valid 
    facade_len_min -- minimum facade lenght allowed 
    facade_len_max -- max facade length allowed 
    facade_depth -- min facade depth allowed 
    max_facade_depth -- max facade depth allowed 

    Returns:
    nothing 
    """

    #open shapefiles 
    aligned_footprints=gpd.read_file(aligned_fp)
    plots= gpd.read_file(final_plots).explode()
    bag_footprints=gpd.read_file(bag_footprints_shp)
    
    #add to aligned footprints the bag attributes 
    bag_footprints_pd= pd.DataFrame(bag_footprints.drop(columns=['geometry', "bouwjaar", "identifica"] )) 
    bag_footprints["gid"]=pd.to_numeric(bag_footprints["gid"]) 
    aligned_footprints["gid"]=pd.to_numeric(aligned_footprints["gid"])
    aligned_footprints=pd.merge(aligned_footprints, bag_footprints_pd, how="left", on="gid")

    #replace height attributes by nan values when height is not valid 
    aligned_footprints.loc[aligned_footprints['height_val'] == "false", ['roof-0.25', 'roof-0.75', 'roof-0.99' ]] = np.nan 

    #initialize all variables 
    df_splitted_footprints = pd.DataFrame(columns=('FeaID', 'bouwjaar', "pand_id", 'roof-0.25', 'roof-0.75', 'roof-0.99', 'geometry'))
    index_r=0
    #avg area is avg footprint area of the bag footprints 
    bag_footprints["area"]=bag_footprints["geometry"].area
    #avg_area = bag_footprints["area"].mean()
    """
    sys.stderr.write(str(avg_area) + "\n")
    sys.stderr.write(str(min_area) + "\n")
    sys.stderr.write(str(facade_len_max) + "\n")
    sys.stderr.write(str(facade_len_min) + "\n")
    sys.stderr.write(str(facade_depth) + "\n")
    sys.stderr.write(str(max_facade_depth) + "\n")
    """
    plots["plot_roof-0.25"]=np.nan 
    plots["plot_roof-0.75"]=np.nan 
    plots["plot_roof-0.99"]=np.nan 

    #for each plot: 
    for index, row in plots.iterrows(): 
        FeaID=row["FeaID"]
        median_roof_025=np.nan
        median_roof_075=np.nan 
        median_roof_099=np.nan 

        #if there are bag footprints aligned with it: 
        FeaID_aligned_fp=aligned_footprints[aligned_footprints["match"]==str(FeaID)]
        if FeaID_aligned_fp.empty==False:
            #compute avg area as mean area of the footprints in the plot 
            FeaID_aligned_fp["area"]=FeaID_aligned_fp["geometry"].area
            #avg_area = FeaID_aligned_fp["area"].mean()
            #get bag fp with valid height 
            fp_valid_height=FeaID_aligned_fp[FeaID_aligned_fp["height_val"] == "true"]
            flat_roofs=fp_valid_height[fp_valid_height["roof_flat"] == "true"]
            not_flat_roofs=fp_valid_height[fp_valid_height["roof_flat"] == "false"]
            #compute median roof height for flat and not flat roofs 
            if flat_roofs.empty==False: 
                median_roof_075=flat_roofs["roof-0.75"].median()
                plots.at[index, "plot_roof-0.75"]=median_roof_075
            if not_flat_roofs.empty==False: 
                median_roof_099=not_flat_roofs["roof-0.99"].median()
                median_roof_025=not_flat_roofs["roof-0.25"].median()
                plots.at[index, "plot_roof-0.25"]=median_roof_025
                plots.at[index, "plot_roof-0.99"]=median_roof_099

        #split plot into footprints 
        input_poly = row["geometry"]
        #make sure plot is oriented ccw 
        input_poly=shapely.geometry.polygon.orient(input_poly, sign=1.0)
        footprints=ptf.split_into_footprints(input_poly, avg_area, min_area, max_perim, facade_len_min, facade_len_max, facade_depth, max_facade_depth)
        
        #for each fp we generated: 
        if FeaID_aligned_fp.empty==False: 
            #check if overlap with any existing bag footprint
            for fp in footprints: 
                overlap=False
                for idx, rw in FeaID_aligned_fp.iterrows(): 
                    if fp.intersects(rw["geometry"])==True: 
                        overlap=True
                        break 
                #if does not overlap, add the fp to output gpd
                if overlap==False: 
                    df_splitted_footprints.loc[index_r] = [row['FeaID'], 9999, "nan", median_roof_025, median_roof_075, median_roof_099, fp]
                    index_r+=1
                else: 
                    pass 
        else: 
            for fp in footprints: 
                df_splitted_footprints.loc[index_r] = [row['FeaID'], 9999, "nan", median_roof_025, median_roof_075, median_roof_099, fp]
                index_r+=1

    #merge aligned and self-generated footprints 
    #nber of columns, col names and types should be the same 
    df_splitted_footprints["ground-0.0"]=np.nan
    df_splitted_footprints["ground-0.1"]=np.nan
    df_splitted_footprints["ground-0.2"]=np.nan
    df_splitted_footprints["ground-0.3"]=np.nan
    df_splitted_footprints["ground-0.4"]=np.nan
    df_splitted_footprints["ground-0.5"]=np.nan
    df_splitted_footprints["nr_ground_pts"]=np.nan 
    df_splitted_footprints["pand_id"]=np.nan 
    aligned_footprints=aligned_footprints[["match", "bouwjaar", "identifica", "roof-0.25", "roof-0.75", "roof-0.99", "geometry", "ground-0.0", "ground-0.1", "ground-0.2", "ground-0.3", "ground-0.4", "ground-0.5", "nr_ground_"]]
    aligned_footprints.columns=["FeaID", "bouwjaar", "pand_id", "roof-0.25", "roof-0.75", "roof-0.99", "geometry", "ground-0.0", "ground-0.1", "ground-0.2", "ground-0.3", "ground-0.4", "ground-0.5", "nr_ground_pts"]
    aligned_footprints["FeaID"]=pd.to_numeric(aligned_footprints["FeaID"])
    df_splitted_footprints["FeaID"]=pd.to_numeric(df_splitted_footprints["FeaID"])
    all_fp=aligned_footprints.append(df_splitted_footprints)
    all_fp["buffer_value"]=0.5
    #assign ground height 
    all_fp=assign_reference_heights.get_ground_reference_heights(input_pc, all_fp)
    #assign roof height to each plot based on neighbours 
    plots_with_height=assign_reference_heights.assign_height_using_neighbours(plots)
    plots_with_height=plots_with_height.drop(columns=["geometry"])
    plots_with_height["FeaID"]=pd.to_numeric(plots_with_height["FeaID"])

    fp_merged_with_plots=pd.merge(all_fp, plots_with_height, how="left", on="FeaID")
    #where the roof height is null for a fp, replace it with the median height of the plot in which the footprint is located 
    fp_merged_with_plots['roof-0.25'] = np.where(pd.isna(fp_merged_with_plots['roof-0.25']), fp_merged_with_plots["plot_roof-0.25"], fp_merged_with_plots['roof-0.25'])
    fp_merged_with_plots['roof-0.75'] = np.where(pd.isna(fp_merged_with_plots['roof-0.75']), fp_merged_with_plots["plot_roof-0.75"], fp_merged_with_plots['roof-0.75'])
    fp_merged_with_plots['roof-0.99'] = np.where(pd.isna(fp_merged_with_plots['roof-0.99']), fp_merged_with_plots["plot_roof-0.99"], fp_merged_with_plots['roof-0.99'])
    print("write output file")
    #prepare dataframe before exporting to shapefile and export 
    all_fp=fp_merged_with_plots[["FeaID", "bouwjaar", "pand_id", "roof-0.25", "roof-0.75", "roof-0.99", "geometry", "ground-0.0", "ground-0.1", "ground-0.2", "ground-0.3", "ground-0.4", "ground-0.5", "nr_ground_pts"]]
    all_fp=all_fp.reset_index(drop=True)
    all_fp=all_fp.reset_index()
    all_fp=all_fp.astype({"nr_ground_pts": int, "bouwjaar": int})
    all_fp.to_file(all_footprints)#, float_format='%.3g')



def generate_final_fp_with_height_2(basepath, final_plots, input_pc,  facade_len_min, facade_len_max, facade_depth, max_facade_depth, avg_area, min_area, min_height, max_height, min_mansard, max_mansard, max_perim=math.inf):  
    """This function is used in the shortened methodology version (ie not including the maps alignment step). 
    It does the following: (i) proceeds to 2D procedural modelling and (ii) assigns roof and ground height attributes to the building footprints. 

    Keyword arguments:
    basepath -- path of the folder where scripts and input/output files are located 
    final_plots -- input shapefile that contains the extracted building plots
    input_pc --  input las file containing a 3D point cloud 
    #Procedural modelling parameters: 
    avg_area -- average building footprint area 
    min_area -- minimum area of footprint to be considered as valid footprint 
    max_perim -- max perimeter allowed for footprint to be valid 
    facade_len_min -- minimum facade lenght allowed 
    facade_len_max -- max facade length allowed 
    facade_depth -- min facade depth allowed 
    max_facade_depth -- max facade depth allowed 
    #Roof height assignment parameters: 
    min_height -- minimum height betw ground and building at the base of the roof (ie does not include roof)
    max_height -- max height betw ground and building at the base of the roof (ie does not include roof)
    min_mansard -- minimum height of the roof (between roof base and top)
    max_mansard -- max height of the roof (between roof base and top)

    Returns:
    nothing 
    """
    
    basepath_plots_to_fp = basepath + input_map + "\\plots_to_fp\\"
    if not os.path.exists(basepath_plots_to_fp):
        os.makedirs(basepath_plots_to_fp)

    plots= gpd.read_file(final_plots).explode()

    ### split plots into footprints 
    df_splitted_footprints = pd.DataFrame(columns=('FeaID', 'geometry'))
    idx=0
    for index, row in plots.iterrows(): 
            input_poly = row["geometry"]
            input_poly=shapely.geometry.polygon.orient(input_poly, sign=1.0)
            footprints=ptf.split_into_footprints(input_poly, avg_area, min_area, max_perim, facade_len_min, facade_len_max, facade_depth, max_facade_depth)
            for fp in footprints: 
                df_splitted_footprints.loc[idx] = [row['FeaID'], fp]
                idx+=1
            idx+=1
    myepsg = full_epsg[0]
    df_splitted_footprints=gpd.GeoDataFrame(df_splitted_footprints, crs=myepsg, geometry="geometry")

    
    df_splitted_footprints["bouwjaar"]=9999
    df_splitted_footprints["ground-0.0"]=np.nan
    df_splitted_footprints["ground-0.1"]=np.nan
    df_splitted_footprints["ground-0.2"]=np.nan
    df_splitted_footprints["ground-0.3"]=np.nan
    df_splitted_footprints["ground-0.4"]=np.nan
    df_splitted_footprints["ground-0.5"]=np.nan

    df_splitted_footprints["ground-0.0"]=pd.to_numeric(df_splitted_footprints["ground-0.0"])
    df_splitted_footprints["ground-0.1"]=pd.to_numeric(df_splitted_footprints["ground-0.1"])
    df_splitted_footprints["ground-0.2"]=pd.to_numeric(df_splitted_footprints["ground-0.2"])
    df_splitted_footprints["ground-0.3"]=pd.to_numeric(df_splitted_footprints["ground-0.3"])
    df_splitted_footprints["ground-0.4"]=pd.to_numeric(df_splitted_footprints["ground-0.4"])
    df_splitted_footprints["ground-0.5"]=pd.to_numeric(df_splitted_footprints["ground-0.5"])

    footprints_splitted = basepath + input_map + "\\plots_to_fp\\" + "splitted_fp_" + input_map + ".shp"
    df_splitted_footprints.to_file(footprints_splitted)

    #footprints_splitted = basepath + input_map + "\\plots_to_fp\\" + "footprints_subset_31370.shp"

    df_splitted_footprints = gpd.read_file(footprints_splitted)
    #assign ground height 
    all_fp=assign_reference_heights.get_ground_reference_heights(input_pc, df_splitted_footprints)
    #assign roof height   
    all_fp = assign_reference_heights.assign_random_roof_heights(all_fp, min_height, max_height, min_mansard, max_mansard)
    all_fp=all_fp.reset_index(drop=True)
    all_fp=all_fp.reset_index()
    all_fp=all_fp.astype({"nr_ground_pts": int})
    all_fp = all_fp.astype({"ground-0.0":'float64', "ground-0.1":'float64', "ground-0.2":'float64', "ground-0.3":'float64', "ground-0.4":'float64', "ground-0.5":'float64'})
    all_fp.to_file(all_footprints)#, float_format='%.3g')
    


def procedural_modelling(all_footprints, blender_path, current_lod2_model, full_epsg, input_map, use_current_city_model = True): 
    """Implements the procedural modelling step to reconstruct 3D buildings from the building footprints dataset. Calls Blender

    Keyword arguments:
    all_footprints -- input shapefile with building footprints 
    blender_path -- path to blender program 
    current_lod2_model -- current 3D city model of the city of interest 
    full_epsg -- tuple containing the EPSG code of the 2D CRS and the EPSG code of the 3D CRS 
    input_map -- string representing the name of the historical map (eg delft_1880)
    use_current_city_model -- binary (T or F) telling whether current 3D city model is used or not to recover some historic buildings (T if maps alignment was implemented, otherwise F)

    Returns:
    nothing 
    """

 
    basepath_pm = basepath + input_map + "\\procedural_modelling\\"
    if not os.path.exists(basepath_pm):
        os.makedirs(basepath_pm)

    myepsg = full_epsg[0]

    sg_footprints = basepath_pm + "final_2PM_footprints.shp"
    output_cityjson_sg = basepath_pm + "2PM_3D_model.json"
    cjio_report = basepath_pm + "cjio_report.txt"
    val3dity_report = basepath_pm + "val3dity_report"
    aligned_3d_model = basepath_pm + "aligned_3D_model.json" 
    cleaned_cityjson_sg = basepath_pm + "2PM_3D_model_cleaned.json" 
    blender_script = basepath + "src_code\\blender_procedural_mod.py"
    final_model = basepath_pm + "histo_model_" + input_map + ".json"
    rulefiles = basepath + "src_code\\rulefiles\\"

    #separate footprints into bag and sg 
    all_footprints_gpd=gpd.read_file(all_footprints)
    sg_footprints_gpd=all_footprints_gpd[all_footprints_gpd["bouwjaar"] == 9999]
    bag_footprints = all_footprints_gpd[all_footprints_gpd["bouwjaar"] != 9999]
    sg_footprints_gpd["bouwjaar"] = int(input_map.split('_')[1])
    sg_footprints_gpd.to_file(sg_footprints)

    #run blender script for sg footprints 
    blender_command = '\"' + blender_path + '\"' " --background --python " + blender_script + " -- " + sg_footprints + " " + output_cityjson_sg + " " + full_epsg[0] +" " + full_epsg[1] + " final_2PM_footprints " + rulefiles
    os.system(blender_command)

    if use_current_city_model: 
        #validates 
        clear_command = "cjio " + output_cityjson_sg + " remove_duplicate_vertices 3 save " + cleaned_cityjson_sg
        os.system(clear_command)
        validate_command = "cjio " + cleaned_cityjson_sg + " validate > " + cjio_report
        os.system(validate_command)
        val3dity_command = val3dity_path + " " + cleaned_cityjson_sg + " --report " + val3dity_report
        os.system(val3dity_command)

        #extract the bag buildings from the bag Lod2 
        cm = cityjson.load(current_lod2_model)
        cm_df = cm.to_dataframe()

        cm_df["identificatie"] = cm_df["identificatie"].apply(lambda x: x[14:])
        cm_df["identificatie"]=pd.to_numeric(cm_df["identificatie"])
        cm_df = cm_df.reset_index()
        cm_df = cm_df[["index", "identificatie"]]
        cm_df.columns = ["id", "pand_id"]
        bag_footprints = bag_footprints[["FeaID", "pand_id"]]
        bag_footprints["pand_id"]= pd.to_numeric(bag_footprints["pand_id"])
        bag_3d_aligned=pd.merge(bag_footprints, cm_df, how="left", on="pand_id")
        bag_3d_aligned["id"]= bag_3d_aligned["id"].apply(lambda x: str(x))
        bag_needed_ids = bag_3d_aligned["id"].tolist()

        co_ids = cm.get_cityobjects(id=bag_needed_ids)
        cm.cityobjects = co_ids
        cityjson.save(cm, aligned_3d_model)

        #merge with SG buildings 
        merge_command = "cjio " + aligned_3d_model + " merge " +  cleaned_cityjson_sg +  " save " + final_model
        os.system(merge_command)
    
    else: 
        #validates 
        clear_command = "cjio " + output_cityjson_sg + " remove_duplicate_vertices 3 save " + final_model
        os.system(clear_command)
        validate_command = "cjio " + final_model + " validate > " + cjio_report
        os.system(validate_command)
        val3dity_command = val3dity_path + " " + final_model + " --report " + val3dity_report
        os.system(val3dity_command)
        


def h3dcm(): 
    """Main function implementing automatic reconstruction of 3D city models from historical maps. Calls the main steps 
    of the methodology: (i) building plots extraction, (ii) reconstruction of individual building footprints and 
    (iii) reconstruction of 3D buildings with 3D procedural modelling. 

    Returns:
    Dictionnary with the time it took to run the different steps 
    """

    #initialise time dictionnary 
    time_dict = dict()
    start = time.time()

    sys.stderr.write("########## RUNNING STEP1: BUILDING PLOTS EXTRACTION ##########\n")
    building_footprints_extraction(basepath, grassbin_path, myepsg, input_map, myfile, ground_truth, training_pts, ground_truth_generalized, alpha, distances, angles, final_plots)
    end = time.time()
    time_dict["Building plots extraction"]= (end-start)/60

    #if current building footprints dataset available, implement complete methodo version 
    if bag_footprints != None: 
        sys.stderr.write("########## RUNNING STEP2: MAPS ALIGNMENT PROCESS ##########\n")
        perform_maps_alignment(basepath, bag_footprints, input_map, myepsg, final_plots)
        end = time.time()
        time_dict["Maps alignment"]= (end-start)/60

        sys.stderr.write("########## RUNNING STEP2: 2D PROCEDURAL MODELLING ##########\n")
        generate_final_fp_with_height(basepath, bag_footprints, aligned_fp, final_plots, input_pc, facade_len_min, facade_len_max, facade_depth, max_facade_depth, avg_area, min_area)
        end = time.time()
        time_dict["2D procedural modelling and height assignment"]= (end-start)/60

        sys.stderr.write("########## RUNNING STEP3: 3D PROCEDURAL MODELLING ##########\n")
        procedural_modelling(all_footprints, blender_path, current_lod2_model, full_epsg, input_map, True) 

    #shortened methodo version 
    else: 
        sys.stderr.write("########## RUNNING STEP2: 2D PROCEDURAL MODELLING ##########\n")
        generate_final_fp_with_height_2(basepath, final_plots, input_pc, facade_len_min, facade_len_max, facade_depth, max_facade_depth, avg_area, min_area, min_height, max_height, min_mansard, max_mansard)
        end = time.time()
        time_dict["2D procedural modelling and height assignment"]= (end-start)/60

        sys.stderr.write("########## RUNNING STEP3: 3D PROCEDURAL MODELLING ##########\n")
        procedural_modelling(all_footprints, blender_path, current_lod2_model, full_epsg, input_map, False) 
   
    end = time.time()
    time_dict["3D procedural modelling"]= (end-start)/60

    sys.stderr.write("Processing time in minutes: \n")
    for step in time_dict: 
        sys.stderr.write(step + " : " + str(round(time_dict[step], 2)) +"\n")

    return time_dict



if __name__ == "__main__":
    basepath = sys.argv[0].partition("src_code")[0] #+ "study_areas\\"
    input_map = sys.argv[1]
    full_epsg = [sys.argv[2], sys.argv[3]]
    grassbin = sys.argv[4]
    blender_path = sys.argv[5]
    val3dity_path = sys.argv[6].replace("\\", "\\\\")

    grassbin_elements = grassbin.split("\\")
    grassbin_path = grassbin_elements[0]
    for el in grassbin_elements[1:-1]:
        grassbin_path = grassbin_path + "\\\\\\\"" + el + "\\\"" 

    grassbin_path = grassbin_path + "\\\\" + grassbin_elements[-1]

    gene_param = input("Enter generalisation parameters as alpha,distanceThreshold,angleThreshold: ")
    gene_param = gene_param.split(",")
    facade_ranges = input("Enter facade depth and facade length ranges as depthMin,depthMax,lengthMin,lengthMax: ")
    facade_ranges = facade_ranges.split(",")
    area_fp = input("Enter average and minimum building footprint area as avgArea,minArea: ")
    area_fp = area_fp.split(",")
    avg_area, min_area = float(area_fp[0]), float(area_fp[1])
    alpha, distances, angles = float(gene_param[0]), float(gene_param[1]), float(gene_param[2])
    facade_depth, max_facade_depth = float(facade_ranges[0]), float(facade_ranges[1])
    facade_len_min, facade_len_max = float(facade_ranges[2]), float(facade_ranges[3])
    
    #get input files 
    myfile = basepath + input_map + "\\input_data\\" + input_map + ".tif"
    training_pts = basepath + input_map + "\\input_data\\training_points_" + input_map + ".shp" 
    ground_truth = basepath + input_map + "\\input_data\\ground_truth_" + input_map + '.shp'
    ground_truth_generalized= basepath + input_map + "\\input_data\\ground_truth_generalized_" + input_map  + ".shp"
    input_pc= basepath + input_map + "\\input_data\\point_cloud_" + input_map.split("_")[0] + ".las"

    bag_footprints= basepath + input_map + "\\input_data\\building_footprints_" + input_map.split("_")[0] + ".shp"
    if os.path.exists(bag_footprints): 
        sys.stderr.write("########## FULL METHODOLOGY VERSION IS IMPLEMENTED ##########\n")
        current_lod2_model = basepath + input_map + "\\input_data\\3D_buildings_" + input_map.split("_")[0] + ".json"

    else: 
        sys.stderr.write("########## SHORTENED METHODOLOGY VERSION IS IMPLEMENTED ##########\n")
        roof_heights = input("Enter roof height and mansard height ranges as minRoof,maxRoof,minMansard,maxMansard: ")
        roof_heights = roof_heights.split(",")
        min_height, max_height = float(roof_heights[0]), float(roof_heights[1]) 
        min_mansard, max_mansard = float(roof_heights[2]), float(roof_heights[3]) 
        bag_footprints = None 
        current_lod2_model = None 

    #define output files 
    all_footprints = basepath + input_map + "\\plots_to_fp\\final_all_footprints.shp" 
    final_plots = basepath + input_map + "\\building_plots_extraction\\" + "final_plots_" + input_map + ".shp"
    aligned_fp = basepath + input_map + "\\plots_to_fp\\aligned_entities_current_dataset.shp" 

    #derive some needed parameters from input 
    myepsg= full_epsg[0]  

    h3dcm()
    