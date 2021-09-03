@echo off & setLocal EnableDelayedExpansion
set output_dec_dir=Breda\out_decomposition
set output_point_dir=Breda\out_points

set input_building_file=Breda\BAG_Breda_clean_clipped.shp
set input_point_file=Breda\Breda_50cm_points.shp
set decomposed_file=!output_dec_dir!\out_Decomposed_Polygon.shp

set output_segmented_point=!output_point_dir!\out_points_Breda_segmented.shp
set output_3d=Breda\Breda_Noord_3D.shp

set log=Breda\Breda_Noord_complete_log.txt
set err_log=Breda\Breda_Noord_complete_error_log.txt

@echo on
rem Cell_Decomposition.exe -i !input_building_file! -odir !output_dec_dir! 1>>!log! 2>>!err_log!
Create_segmented_clipped_shp.exe -i !input_point_file! -shp -spat 12 -clipsrc !decomposed_file! -o !output_segmented_point! 1>>!log! 2>>!err_log!
Create_3D_buildings.exe -ip !decomposed_file! -is !output_segmented_point! -o !output_3d! -overwrite 1>>!log! 2>>!err_log!