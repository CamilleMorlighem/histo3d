#!/usr/bin/env python3
import skgeom
import os, sys
import subprocess
import shutil
import binascii
import tempfile
import csv
import math
import logging
import geopandas as gpd
import statistics as stat 
import shapely
import eval_generalisation
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import metrics


"""
This script enables the building plots extraction step of the general methodology worklow. It contains (i) the OBIA segmentation and classification method, 
(ii) functions implementing generalisation methods and (iii) functions for computing metrics and thus assessing the quality of the building plots extraction step. 
It is used by the script main_whole_workflow.py. 
"""

def OBIA_segmentation_grass(grassbin, myepsg, myfile, basepath, input_map, training_pts, raw_building_plots): 
    """Performs OBIA segmentation in GRASS GIS 

    Keyword arguments:
    grassbin -- GRASS bin environment variable 
    myepsg -- EPSG code of the coordinate reference system 
    myfile -- input historical map (.tif file)
    basepath -- path of the folder where scripts and input/output files are located 
    input_map -- string representing the name of the historical map (eg delft_1880)
    training_pts -- shapefile with training points 
    raw_building_plots -- output shapefile that will contain the extracted building plots 

    Returns:
    nothing 
    """
    os.environ['GRASSBIN']= grassbin
    from grass_session import Session
    import grass.script as gscript
    import grass.script.setup as gsetup
    from grass.pygrass.modules.shortcuts import general as g
    from grass.pygrass.modules.shortcuts import raster as r
    from grass.pygrass.modules.shortcuts import vector as v
    from grass.pygrass.modules.shortcuts import temporal as t
    from grass.pygrass.modules import Module

    ## get gisbase path 
    cmd = "{grassbin} --config path".format(grassbin=grassbin)
    try:
        p = subprocess.Popen(cmd, shell=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
    except OSError as error:
        sys.exit("ERROR: Cannot find GRASS GIS start script"
                " {cmd}: {error}".format(cmd=startcmd[0], error=error))
    if p.returncode != 0:
        sys.exit("ERROR: Issues running GRASS GIS start script"
                " {cmd}: {error}"
                .format(cmd=' '.join(startcmd), error=err))
    gisbase = out.strip(os.linesep.encode())

    ## set environment variables
    init_proj_lib_path = os.environ["PROJ_LIB"] 
    os.environ["PROJ_LIB"] = gisbase.decode() + '\\share\\proj'
    os.environ['GISBASE'] = gisbase.decode()
    os.environ['PYTHONPATH']=os.path.join(os.environ['GISBASE'], 'etc', 'python')
    sys.path.append(os.path.join(os.environ['GISBASE'], 'etc', 'python'))
    os.environ['PATH']=os.path.join(os.environ['GISBASE'], 'etc', 'python;')+gisbase.decode()+"\\lib;" + gisbase.decode()+"\\bin;"+ gisbase.decode()+"\\extrabin;"+ gisbase.decode()+"\\scripts;"+ os.environ['GRASS_ADDON_BASE']+"\\bin;"+os.environ['PATH']

    ## create new location if does not exist already: starts by creating a permanent mapset (ALWAYS) and then create the new mapset
    #with statement both opens and closes the session 
    mygisdb = basepath + "GRASS_GIS_DB"
    if not os.path.exists(mygisdb):
        os.makedirs(mygisdb)

    #mylocation = input_map.split('_')[0]
    mylocation = input_map
    mymapset = input_map

    if not os.path.isdir(mygisdb+"/" + mylocation): 
        with Session(gisdb=mygisdb, location=mylocation, create_opts=myepsg):
            print("Created new location with permanent mapset")

    ## Launch a new session 
    #open session with the permanent mapset and then with g.mapset, creates a new one and switch to that one 
    gsetup.init(gisbase.decode(), mygisdb, mylocation, "PERMANENT")
    gscript.run_command("g.mapset", location = mylocation, flags="c",mapset=mymapset)

    ## install needed extensions 
    if not gscript.find_program('r.neighborhoodmatrix', '--help'):
        gscript.run_command("g.extension", extension="r.neighborhoodmatrix")
    if not gscript.find_program('i.segment.uspo', '--help'):
        gscript.run_command("g.extension", extension="i.segment.uspo")
    if not gscript.find_program('i.segment.stats', '--help'):
        gscript.run_command("g.extension", extension="i.segment.stats")
    if not gscript.find_program('v.class.mlR', '--help'):
        gscript.run_command("g.extension", extension="v.class.mlR")
    if not gscript.find_program('r.object.thickness', '--help'):
        gscript.run_command("g.extension", extension="r.object.thickness")
    if not gscript.find_program('r.fill.category', '--help'):
        gscript.run_command("g.extension", extension="r.fill.category")
    
    gscript.run_command("db.connect", database = '$GISDBASE/$LOCATION_NAME/$MAPSET/sqlite/sqlite.db')
    
    gscript.run_command("r.in.gdal", input=myfile, output=input_map, overwrite=True) #creates a group called input_map
    gscript.run_command("i.group", group=input_map+"_group", input= [input_map + ".blue@" + mymapset, input_map + ".red@" + mymapset, input_map + ".green@" + mymapset], overwrite=True)
    gscript.run_command("g.region", raster=input_map + ".blue@"+mymapset, save="reg", overwrite=True)
    csvfile=gscript.gisenv()['GISDBASE'] + '/' + gscript.gisenv()['LOCATION_NAME'] + '/' + gscript.gisenv()['MAPSET'] + '/ortho_uspo.csv' 

    gscript.run_command("i.segment.uspo", regions ="reg", group=input_map+"_group", output=csvfile, segment_map="segment_optimized", 
    threshold_start=0.07, threshold_stop=0.1, threshold_step=0.01, minsizes=(7,8,10), number_best=3, processes=4, overwrite=True)#, memory=2000) 
    
    #gscript.run_command("i.segment", group=input_map+"_group", output="segment_optimized_reg_rank1", threshold="0.1", minsize=7, overwrite=True)

    gscript.run_command("v.in.ogr", input=training_pts, output="training_pts", overwrite=True) 

    text_label=gscript.read_command("db.select", sql="select distinct label from training_pts where class='text'")
    text_label=int(text_label.strip("label\r\n"))

    symb_label=gscript.read_command("db.select", sql="select distinct label from training_pts where class='symbol'")
    if symb_label.strip("label\r\n")!="": 
        symb_label=int(symb_label.strip("label\r\n"))
    
    gscript.run_command("i.segment.stats", map="segment_optimized_reg_rank1", 
    rasters=(input_map + ".blue", input_map + ".red", input_map + ".green"), raster_statistics=("median", "mean", "first_quart", "third_quart"), 
    area_measures=("compact_circle","fd","compact_square"), vectormap="segment_vector_map", overwrite=True)
    
    gscript.run_command("v.vect.stats", points="training_pts", areas="segment_vector_map", method="maximum", count_column="nber_pts", points_column="label", stats_column="label", overwrite=True) 
    gscript.run_command("v.db.dropcolumn", map="segment_vector_map", columns="nber_pts")
    gscript.run_command("v.extract", input="segment_vector_map", output="training_map", where="label not NULL", overwrite=True)
    gscript.run_command("v.db.dropcolumn", map="segment_vector_map", columns="label")

    gscript.run_command("v.class.mlR", segments_map="segment_vector_map", training_map="training_map",
    weighting_modes=("smv"), train_class_column="label", classifiers=[ "rf", "knn", "rpart", "svmRadial", "svmLinear"] , classified_map="classified_raster_map", 
    raster_segments_map="segment_optimized_reg_rank1", overwrite=True) #, r_script_file="script2.r")
    
    values_text=gscript.read_command("r.object.thickness", input="classified_raster_map_smv", category=text_label, tspace="1", tsize="15", vmedian="medians_text", transects="transects_text", overwrite=True)
    values_text=values_text.strip("()\r\n")
    values_text=values_text.split(",")
    values_text=list(map(float, values_text))

    values_symb=gscript.read_command("r.object.thickness", input="classified_raster_map_smv", category=symb_label, tspace="3", tsize="20", vmedian="medians_symb", transects="transects_symb", overwrite=True)
    values_symb=values_symb.strip("()\r\n")
    values_symb=values_symb.split(",")
    values_symb=list(map(float, values_symb))
    
    nsize_symb=int(round(values_symb[1]/2+1))
    if (nsize_symb % 2) == 0:
        nsize_symb=nsize_symb+1

    nsize_text=int(round(values_text[1]/2+1))
    if (nsize_text % 2) == 0:
        nsize_text=nsize_text+1
    
    gscript.run_command("r.fill.category", input="classified_raster_map_smv", output="classified_without_txt", category=text_label, maxiter="25", nsize=nsize_text, overwrite=True)
    gscript.run_command("r.fill.category", input="classified_without_txt", output="classified_without_symbols", category=symb_label, maxiter="25", nsize=nsize_symb, overwrite=True)
    
    gscript.run_command("r.to.vect", input="classified_without_symbols", output="classified_vector_map", type="area", overwrite=True)
    #gscript.run_command("r.to.vect", input="classified_raster_map_smv", output="classified_vector_map", type="area", overwrite=True)

    build_label=gscript.read_command("db.select", sql="select distinct label from training_pts where class='building'")
    build_label=build_label.strip("label\r\n")
    #build_label_2=gscript.read_command("db.select", sql="select distinct label from training_pts where class='building_2'")
    #build_label_2=build_label_2.strip("label\r\n")
    sql_select_buildings="value="+ build_label #+ " or value=" + build_label_2
    gscript.run_command("v.extract", input="classified_vector_map", output="raw_building_plots", where=sql_select_buildings, overwrite=True)
    gscript.run_command("v.db.dropcolumn", map="raw_building_plots", columns="label,value")

    gscript.run_command("v.out.ogr", input="raw_building_plots", output=raw_building_plots, format="ESRI_Shapefile", overwrite=True)

    ## end the session
    gsetup.finish()

    #restore initial proj lib path 
    os.environ["PROJ_LIB"] = init_proj_lib_path



def area_filter_and_assess_classif(raw_building_plots, metrics_csv, metrics_txt, raw_building_plots_area_filter, ground_truth, area_threshold = 150): 
    """Remove tiny building plots (come from missclassifications) and compute a serie of metrics to assess the building plots extraction step 

    Keyword arguments:
    raw_building_plots -- input shapefile with extracted raw building plots 
    metrics_csv -- csv file where to write all the computed metrics (detailed)
    metrics_txt -- text file that will contain global metric values (summarized)
    raw_building_plots_area_filter -- output shapefile that will contain the filtered building plots 
    ground_truth -- points shapefile with the centroids of the ground truth building plots 
    area_treshold -- area threshold to filter out missclassifications (salt and pepper effect)

    Returns:
    nothing 
    """

    OBIA_output=gpd.read_file(raw_building_plots).explode()

    # filter on area (remove small missclassification errors)
    OBIA_output['area_raw']=OBIA_output.area
    OBIA_output.area_raw=OBIA_output.area_raw.round(3)
    OBIA_output=OBIA_output.loc[OBIA_output["area_raw"] > area_threshold]

    # add some metrics
    OBIA_output=OBIA_output.reset_index(drop=True)
    OBIA_output['cat']= OBIA_output.index
    OBIA_output=OBIA_output.rename(columns={'cat': "FeaID"})
    OBIA_output=metrics.add_metrics(OBIA_output, "raw", metrics_csv)

    # assess obia segmentation results
    try:  
        ground_truth= gpd.read_file(ground_truth)
        join=gpd.sjoin(ground_truth, OBIA_output[['geometry']], how='left', op='intersects')

        # drop duplicates in case the same ground truth point was localised several time in the same plot (overlap) (improve falsely the number of true positive)
        join=join.drop_duplicates(subset=['index_right'])
        # count() function does not count NA (point has no corresponding plot)
        num_positive = join["index_right"].count()
        precision = round(num_positive/len(OBIA_output),4)
        recall = round(num_positive/len(ground_truth),4)
        fscore = round(2.0 * ((precision * recall) / (precision + recall)),4)
        
        f= open(metrics_txt,"w+")
        f.write("RESULTS OBIA SEGMENTATION: \r\n")
        f.write("\t * Precision: %f\r\n" % (precision))
        f.write("\t * Recall: %f\r\n" % (recall))
        f.write("\t * F-measure: %f\r\n" % (fscore))
        f.close()
    except: 
        pass 

    OBIA_output.to_file(raw_building_plots_area_filter)
    


def final_generalization_workflow(raw_building_plots_area_filter, csv_gen_results, alpha_shape_script, input_map, out_alpha_shape, alpha, distances, angles, commandeur_script, out_generalized, final_plots_file, ground_truth_generalized, metrics_csv, metrics_txt, perpendicularity = False, parallelism_and_perpend = False): 
    """Generalizes building plots with (1) alpha-shape based generalisation and (2) Commandeur's generalisation method. In addition, compute metrics to assess the generalisation implemented 

    Keyword arguments:
    raw_building_plots_area_filter -- input shapefile with building plots to be generalized 
    csv_gen_results -- output csv file that will contain the generalisation metrics (shape similarity and turning function) for all the x building plots taken into account to compute these values 
    alpha_shape_scrit -- path to R file with alpha shape generalisation code 
    input_map -- string representing the name of the historical map (eg delft_1880)
    out_alpha_shape -- output shapefile that will contain building plots generalised with alpha shape generalisation 
    alpha -- integer representing the alpha value 
    distances -- distance threshold value in commandeur generalisation 
    angles -- angle threshold value in commandeur generalisation 
    commandeur_script -- path to Commandeur code 
    out_generalized -- output shapefile that will contain building plots generalised with commandeur generalisation 
    final_plots_file -- output shapefile that will contain building plots after both generalisation methods are implemented
    ground_truth_generalized -- shapefile with ground truth building plots manually digitalized 
    perpendicularity -- binary (T or F) telling whether perpendiculary should be ensured in commandeur code 
    parallelism_and_perpend -- binary (T or F) telling whether both perpendiculary and // should be ensured in commandeur code 
    metrics_csv -- csv file where to write all the computed metrics (detailed)
    metrics_txt -- text file that will contain global metric values (summarized)

    Returns:
    nothing 
    """
    
    alpha = str(alpha)
    # alpha shape generalization
    
    alpha_shape_simplification(raw_building_plots_area_filter, alpha_shape_script, input_map, metrics_csv, out_alpha_shape, (alpha,))
        
    # add metrics 
    alpha_shape_output=gpd.read_file(out_alpha_shape + str(int(float(alpha))) +'.shp')
    alpha_shape_output=metrics.add_metrics(alpha_shape_output, "alpha", metrics_csv)
    alpha_shape_output.to_file(out_alpha_shape + str(int(float(alpha))) + ".shp")
    
    generalized_plots_path = tom_commandeur_generalization((distances,), (angles,), raw_building_plots_area_filter, int(float(alpha)), commandeur_script, out_alpha_shape, out_generalized)

    #add_metrics 
    generalized_output=gpd.read_file(generalized_plots_path).explode()

    #whole_df=generalized_output.merge(alpha_shape_output, how='left', on="FeaID")
    generalized_output=metrics.add_metrics(generalized_output, "generalized", metrics_csv)
    generalized_output.to_file(final_plots_file, index=False)
    
    # assess the generalisation results: 
    try: 
        generalization_results=eval_generalisation.eval_generalisation(ground_truth_generalized, final_plots_file) 
        # write results to csv file 
        pd.DataFrame(generalization_results).to_csv(csv_gen_results, sep=";")
        avg_shape_similarity= round(stat.mean(generalization_results["shape_similarity"]), 4)
        avg_diff_turning_fct = round(stat.mean(generalization_results["diff_turning_f"]), 4)
        #write statistics results to text file 
        f= open(metrics_txt,"a+")
        f.write("AFTER COMMANDEUR' S GENERALISATION: \r\n")
        f.write("\t * Average difference between shapes: %f\r\n" % (avg_shape_similarity))
        f.write("\t * Average difference between shape turning functions: %f\r\n" % (avg_diff_turning_fct))
        f.close()
    except: 
        pass 
    
    
def alpha_shape_simplification(raw_building_plots_area_filter, alpha_shape_script, input_map, metrics_csv, out_alpha_shape, alpha_parameters):
    """Performs alpha shape generalisation 

    Keyword arguments:
    raw_building_plots_area_filter -- input shapefile with building plots to be generalized 
    alpha_shape_scrit -- path to R file with alpha shape generalisation code 
    input_map -- string representing the name of the historical map (eg delft_1880)
    out_alpha_shape -- output shapefile that will contain building plots generalised with alpha shape generalisation 
    alpha_parameters -- list/tuple with alpha values 
    metrics_csv -- do not consider; used in deprecated version of the function 

    Returns:
    nothing 
    """

    # alpha shape generalization
    for alpha in alpha_parameters: 
        alpha = str(alpha)
        command = 'R --vanilla --silent --slave -f ' + alpha_shape_script + ' --args ' + alpha + " " + raw_building_plots_area_filter + ' ' + 'area_filter_raw_building_plots_' + input_map +' ' + out_alpha_shape 
        logging.debug(command)
        os.system(command)
        

def tom_commandeur_generalization(distances, angles, raw_building_plots, alpha, commandeur_script, out_alpha_shape, out_generalized, perpendicularity = False, parallelism_and_perpend = False): 
    """Performs Commandeur's generalisation 

    Keyword arguments:
    distances -- distance threshold value in commandeur generalisation 
    angles -- angle threshold value in commandeur generalisation
    raw_building_plots -- input shapefile with the extracted (raw) building plots (before first generalisation method), used only for the name of the shapefile 
    alpha -- integer representing the alpha value that was used in prior generalisation method 
    commandeur_script -- path to Commandeur code 
    out_alpha_shape -- input shapefile with building plots to be generalised 
    out_generalized -- output shapefile that will contain building plots generalised with commandeur generalisation 
    perpendicularity -- binary (T or F) telling whether perpendiculary should be ensured in commandeur code 
    parallelism_and_perpend -- binary (T or F) telling whether both perpendiculary and // should be ensured in commandeur code 
    
    Returns:
    Path of shapefile with generalised building plots 
    """

    alpha = str(alpha)
    for dist in distances: 
        for angle in angles: 
            generalized_plots_path = out_generalized + '_' + alpha + '_' + str(dist) + '_' + str(angle) 
            command = commandeur_script + ' -i ' + out_alpha_shape + alpha +'.shp' + ' -odir ' + generalized_plots_path + ' -distance ' + str(dist) + ' -angle ' + str(angle) + ' -overwrite'
            os.system(command)
            os.system("copy " + raw_building_plots[:-4] + ".prj " + generalized_plots_path + "\\out_Generalized_Polygon.prj")
            if parallelism_and_perpend == True: 
                command = commandeur_script + ' -i ' + out_alpha_shape + alpha +'.shp' + ' -odir ' + generalized_plots_path + "_pra" + ' -distance ' + str(dist) + ' -angle ' + str(angle) + ' -parallel -rightAngles' + ' -overwrite'
                os.system(command)
                logging.debug(command)
                os.system("copy " + raw_building_plots[:-4] + ".prj " + generalized_plots_path +  "_pra" + "\\out_Generalized_Polygon.prj")
            if perpendicularity ==True: 
                command = commandeur_script + ' -i ' + out_alpha_shape + alpha +'.shp' + ' -odir ' + generalized_plots_path + "_ra" + ' -distance ' + str(dist) + ' -angle ' + str(angle) + ' -rightAngles' + ' -overwrite'
                os.system(command)
                logging.debug(command)
                os.system("copy " + raw_building_plots[:-4] + ".prj " + generalized_plots_path +  "_ra" + "\\out_Generalized_Polygon.prj")

    return generalized_plots_path + "\\out_Generalized_Polygon.shp"


def shapely_generalisation(out_alpha_shape, raw_building_plots, output_file, tolerance, preserve_topo=True, alpha=None): 
    """Performs shapely generalisation 

    Keyword arguments:
    out_alpha_shape -- input shapefile with raw building plots to be generalised 
    raw_building_plots -- input shapefile with building plots, already generalised a first time with alpha-shape method, to be generalised 
    #Only one of the two files is used as input for generalisation, depending on the value of the alpha parameter
    output_file -- output shapefile that will contain generalised building plots 
    tolerance -- tolerance value for generalisation 
    preserve_topo -- binary (T or F) indicating whether the topology should be preserved or not 
    alpha -- alpha parameter value of the alpha shape generalisation method implemented, IF implemented. If not implemented, default value is chosen, which is none 

    Returns:
    nothing 
    """
    alpha = str(alpha)
    if alpha != None: 
        input_geom=gpd.read_file(out_alpha_shape + str(int(float(alpha))) +'.shp')
        output= output_file + alpha + "_" + str(tolerance) +".shp"

    else: 
        input_geom=gpd.read_file(raw_building_plots)
        output= output_file + "_raw" + "_" + str(tolerance) +".shp"

    simplified_df = input_geom.copy()

    for index, row in simplified_df.iterrows(): 
        simplified=row["geometry"].simplify(tolerance, preserve_topology=preserve_topo)
        simplified_df.at[index,"geometry"]=simplified

    simplified_df.to_file(output)


def scikit_geom_generalisation(out_alpha_shape, raw_building_plots, output_file, tolerance, mode, preserve_topo=True, alpha=None): 
    """Performs scikit geom generalisation 

    Keyword arguments:
    out_alpha_shape -- input shapefile with raw building plots to be generalised 
    raw_building_plots -- input shapefile with building plots, already generalised a first time with alpha-shape method, to be generalised 
    #Only one of the two files is used as input for generalisation, depending on the value of the alpha parameter
    output_file -- output shapefile that will contain generalised building plots 
    tolerance -- tolerance value for generalisation 
    mode -- string indicating the mode to be used for generalisation (ratio, count or sqeuclidean)
    preserve_topo -- binary (T or F) indicating whether the topology should be preserved or not 
    alpha -- alpha parameter value of the alpha shape generalisation method implemented, IF implemented. If not implemented, default value is chosen, which is none 

    Returns:
    nothing 
    """

    alpha = str(alpha)
    if alpha != None: 
        input_geom=gpd.read_file(out_alpha_shape + str(int(float(alpha))) +'.shp')
        output=output_file + mode + "_" + alpha + "_" + str(tolerance) +".shp"

    else: 
        input_geom=gpd.read_file(raw_building_plots)
        output=output_file + mode + "_" + str(tolerance) +".shp"

    simplified_df = input_geom.copy()

    for index, row in simplified_df.iterrows(): 
        x1,y1 = row["geometry"].exterior.xy
    
        p=list(zip(list(x1), list(y1)))
        poly=skgeom.Polygon(np.array(p))

        simplified=skgeom.simplify(poly,tolerance, mode, preserve_topology=preserve_topo)
        simplified_df.at[index,"geometry"]=shapely.geometry.Polygon(simplified.coords)

    simplified_df.to_file(output)
