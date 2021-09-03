import bpy, bmesh
import mathutils 
import math 
import random 
import os 
from mathutils.bvhtree import BVHTree
from mathutils.kdtree import KDTree
from mathutils import Matrix
from mathutils import Vector
from math import radians, sqrt
import json
import sys

argv = sys.argv
argv = argv[argv.index("--") + 1:]
myepsg= argv[2]
epsg_3d = argv[3]
coll_name = argv[4]
rulefiles_path = argv[5]


def prepare_blender_scene(): 
    """Prepare blender scene parameters and load shapefile 

    Returns:
    nothing  
    """

    ### Clear scene from any object 
    for obj in bpy.context.scene.objects:
        obj.select_set(True)
    bpy.ops.object.delete()
    
    ### Enable needed add-ons
    bpy.ops.preferences.addon_enable(module="bcga-master")
    bpy.ops.preferences.addon_enable(module="BlenderGIS-master")
    bpy.ops.preferences.addon_enable(module="Up3date-main")
    
    ### Prepare the scene crs 
    #create new default crs 
    prefs = bpy.context.preferences.addons["BlenderGIS-master"].preferences
    #bpy.ops.bgis.add_predef_crs("EXEC_DEFAULT", crs= myepsg, name= "New RD")
    data = json.loads(prefs.predefCrsJson)
    data.append([myepsg, "New crs", "none"])
    prefs.predefCrsJson = json.dumps(data)
    bpy.ops.wm.save_userpref()
    
    #set it as the predefCrs
    bpy.context.preferences.addons["BlenderGIS-master"].preferences.predefCrs=myepsg 
    #set predefCrs as scene crs 
    #bpy.ops.geoscene.set_crs("EXEC_DEFAULT") 

    ### Load shapefile using blender-gis addon 
    bpy.ops.importgis.shapefile_props_dialog("EXEC_DEFAULT", filepath=r"{}".format(argv[0]), 
    vertsElevSource='NONE', separateObjects=True, shpCRS=myepsg)


def find_facade_edge_idx(obj, bm, world_matrix, kd, obj_to_edge_dict, obj_to_facade): 
    """Find edge where facade will be built

    Keyword arguments:
    obj -- object being processed 
    bm -- mesh of the object being processed 
    world_matrix -- world matrix of the blender scene 
    kd -- kd tree containing the mid points (and their ID) of all edges of all building footprints 
    obj_to_edge_dict -- dictionnary mapping from the id of the edges mid points (key) to the name of the object to which they belong (value)
    obj_to_facade -- dict mapping from the name of the object (key) to the street edge of the object (value) stored as a tuple(end point 1, end point 2)

    Returns:
    index of the street edge 
    """

    #first step consists in finding out which edge of the mesh are shared with other meshes ==> no facade on those ones 
    non_overlapping_e= []
    neighbours = []
    #for each edge of the mesh 
    for e in bm.edges: 
        #get mid point of the edge 
        pt0_global = world_matrix @ e.verts[0].co
        pt1_global = world_matrix @ e.verts[1].co
        half_pt_x = (pt0_global.x + pt1_global.x)/2
        half_pt_y = (pt0_global.y + pt1_global.y)/2
        half_pt_z = (pt0_global.z + pt1_global.z)/2
        
        #find the 5 nearest midpts  
        nearest_midpts = kd.find_n(mathutils.Vector((half_pt_x, half_pt_y, 0)), 5)
        edge_shared = False
       
        #if one of that midpt is on the edge segment, it means that the edge is shared by two meshes 
        for midpt_tup in nearest_midpts: 
            midpt = midpt_tup[0]
            neighbour_name = obj_to_edge_dict[(midpt_tup[1],)] 
            if neighbour_name != obj.name: 
                if collinear(pt0_global, pt1_global, midpt) and is_on(pt0_global, pt1_global, midpt): 
                    #the edge is overlapping with another edge 
                    edge_shared = True  
                    #keep list with neighbours of the mesh 
                    neighbours.append(neighbour_name)
                if edge_shared: 
                    break 

        #if the edge is not shared, it is a potential facade edge 
        if edge_shared == False: 
            non_overlapping_e.append(e.index) 
   
    #CASE 1: if there is only one edge not shared, make the facade at this one
    edge_idx = "not_found_yet" 
    if len(non_overlapping_e) == 1: 
        edge_idx = non_overlapping_e[0]
        
    #CASE 2: if there are more than one edge without neighbour, check if neighbouring buildings are already made 3D to put the facade on the same side 
    else: 
        shortest_dist = math.inf 
        looking_for_facade_edge = True 
        #for each neighbour 
        for neighbour in neighbours: 
            #if the neighbour has been extruded (already have a facade)
            if neighbour in obj_to_facade: 
                vtx_neighbour = obj_to_facade[neighbour]
                dir_neighb = vtx_neighbour[0] - vtx_neighbour[1]
                half_x_neighb = (vtx_neighbour[0].x + vtx_neighbour[1].x)/2
                half_y_neighb = (vtx_neighbour[0].y + vtx_neighbour[1].y)/2
                half_coord_neighb = mathutils.Vector((half_x_neighb, half_y_neighb))
                #for each potential facade edge 
                for e in non_overlapping_e: 
                    v0_global = world_matrix @ bm.edges[e].verts[0].co 
                    v1_global = world_matrix @ bm.edges[e].verts[1].co
                    dir_e = v0_global - v1_global 
                    length_prod = bm.edges[e].calc_length() * (vtx_neighbour[0]-vtx_neighbour[1]).length 
                    dot_prod = dir_e.x * dir_neighb.x + dir_e.y  * dir_neighb.y 
                    #if the neighbour facade edge and this edge are parallel (dot product = product of their lenght)
                    if (dot_prod - length_prod) >= -0.5 and (dot_prod - length_prod) <= 0.5: 
                        #get the distance between the half point of the two edges 
                        half_x_e = (v0_global.x + v1_global.x)/2
                        half_y_e = (v0_global.y + v1_global.y)/2
                        half_coord_e = mathutils.Vector((half_x_e, half_y_e))
                        dist = (half_coord_e - half_coord_neighb).length 
                        #the facade edge will be the closest to the neighbour facade edge 
                        if dist < shortest_dist: 
                            edge_idx = e 
                            looking_for_facade_edge = False 
                            
            if looking_for_facade_edge == False: 
                break 
        
        #CASE 3: if the mesh has neighbours but none of them have a facade already, check alignment with neighbouring edges 
        if looking_for_facade_edge: 
            #if the polygon has neighbours 
            if neighbours!= []: 
                counts= []
                #for each potential facade edge, we get is midpoint and count the number of edge of the neibghbours with which it is collinear 
                for e in non_overlapping_e: 
                    count_collinear_edges=0
                    v0_global = world_matrix @ bm.edges[e].verts[0].co 
                    v1_global = world_matrix @ bm.edges[e].verts[1].co
                    half_x_e = (v0_global.x + v1_global.x)/2
                    half_y_e = (v0_global.y + v1_global.y)/2
                    half_coord_e = mathutils.Vector((half_x_e, half_y_e))
                    
                    for neighbour in neighbours: 
                        neigh_obj = bpy.context.scene.objects[neighbour]
                        neigh_mesh = neigh_obj.data
                        
                        edge_numb = 0
                        for edge in neigh_mesh.edges: 
                            c1 = (neigh_mesh.vertices[edge.vertices[0]].co.z - neigh_obj["ground-0.2"]) >= -0.2 and (neigh_mesh.vertices[edge.vertices[0]].co.z - neigh_obj["ground-0.2"]) <= 0.2 
                            c2 = (neigh_mesh.vertices[edge.vertices[1]].co.z - neigh_obj["ground-0.2"]) >= -0.2 and (neigh_mesh.vertices[edge.vertices[1]].co.z - neigh_obj["ground-0.2"]) <= 0.2 
                            if c1 and c2:
                                pt0_global = neigh_obj.matrix_world @ neigh_mesh.vertices[edge.vertices[0]].co
                                pt1_global = neigh_obj.matrix_world @ neigh_mesh.vertices[edge.vertices[1]].co
                                edge_numb +=1
                                if collinear(pt0_global, pt1_global, half_coord_e): 
                                    count_collinear_edges += 1 
    
                    counts.append((count_collinear_edges, e)) 
                
                #the potential facade edge with the highest count is chosen as the facade edge  (from the observation that facade edges are aligned with neighbours)
                
                counts.sort(reverse = True) 
                if counts != []: 
                    if counts[0][0] != counts[1][0]: 
                        edge_idx = counts[0][1]
    #CASE 4: if we still havent find any edge facade, it will be the shortest edge of the mesh    
    if edge_idx == "not_found_yet": 
        shortest_len = math.inf 
        for index in non_overlapping_e: 
            edge_len = bm.edges[index].calc_length()
            if edge_len <= shortest_len: 
                shortest_len = edge_len 
                edge_idx = index   
    if edge_idx == "not_found_yet": 
        edge_idx = 0 
        
    #add facade edge to dictionary keeping track of facade edges 
    start_pt = bm.edges[edge_idx].verts[0]
    end_pt = bm.edges[edge_idx].verts[1]    
    s_pt_global = world_matrix @ start_pt.co 
    e_pt_global = world_matrix @ end_pt.co 
    obj_to_facade[obj.name] = (s_pt_global, e_pt_global) 

    return edge_idx  



def procedural_modelling(kd, obj_to_edge_dict): 
    """Main function implementing the procedural modelling with the BCGA addon. Reconstruct the 3D buildings from the building footprints by (i) choosing 
    a roof type, (ii) identifying the street edge and (iii) generating and executing a rule file per building footprint 

    Keyword arguments:
    kd -- kd tree containing the mid points (and their ID) of all edges of all building footprints 
    obj_to_edge_dict -- dictionnary mapping from the id of the edges mid points (key) to the name of the object to which they belong (value)
 
    Returns:
    nothing
    """

    #make kd tree for keeping track of the edges where special facade is built 
    #size = 2*len(bpy.data.collections['test_footprints_blender_2'].objects) 
    #kd_facade_vertices = KDTree(size)
    obj_to_facade = dict()
    
    h=0
    #for each footprint 
    for obj in bpy.data.collections[coll_name].objects: 
        #make object active object 
        bpy.context.view_layer.objects.active = obj
        mesh=obj.data 
        n_vertices=len(mesh.vertices)
       
        #randomly choose a type of roofs which depend on the year of construction and the nber of vertices 
        if (obj["roof-0.99"] -  obj["roof-0.25"]) < 1: 
            roof_type = "flat" 
        elif n_vertices == 4: 
            if obj["bouwjaar"] > 1930: 
                roofs_lst=["flat", "hip", "mansard", "gable", "dutch"]
                weight=[25, 25, 10, 25, 15]
            else: 
                roofs_lst=["flat", "hip", "mansard", "gable", "dutch"]
                weight=[10, 20, 10, 20, 40]
            roof_type = random.choices(roofs_lst, weights=weight, k=1)[0] 
        else: 
            #complex footprint limits the types of roofs available 
            if obj["bouwjaar"] > 1930: 
                roofs_lst=["flat", "complex_hip", "mansard"]
                weight=[40, 40, 20]
                roofs_lst=["flat", "mansard"]
                weight=[60, 40]
            else: 
                roofs_lst=["flat", "complex_hip", "mansard"]
                weight=[15, 65, 20]
                roofs_lst=["flat",  "mansard"]
                weight=[60,40]
            roof_type = random.choices(roofs_lst, weights=weight, k=1)[0] 
                
        #choose random parameters values depending on the type of roof previously selected 
        #store the parameters in a dictionary dict_param() 
      
        param_dict= dict()
        if roof_type == "flat": 
            obj["is_flat"] = 1
            rulepath= rulefiles_path + "flat_roof" 
            param_dict["height"]= obj["roof-0.25"] - obj["ground-0.2"]
        
        elif roof_type == "hip": 
            obj["is_flat"] = 0 
            rulepath=rulefiles_path + 'hip_roof2'
            param_dict["height"] = obj["roof-0.25"] - obj["ground-0.2"]
            param_dict["top_height"] = obj["roof-0.99"] - obj["ground-0.2"]
            
            world_matrix = obj.matrix_world 
            bpy.ops.object.mode_set(mode="EDIT")
            bm=bmesh.from_edit_mesh(mesh)
            bm.edges.ensure_lookup_table()

            #find an edge for building the facade based on neighbouring polygons 
            edge_idx = find_facade_edge_idx(obj, bm, world_matrix, kd, obj_to_edge_dict, obj_to_facade) 
            
            start_pt = bm.edges[edge_idx].verts[0]
            end_pt = bm.edges[edge_idx].verts[1] 
            current_edge = bm.edges[edge_idx]
            
            #add the mid point of the facade edge to the dictionary 
            half_x = (start_pt.co.x + end_pt.co.x)/2
            half_y = (start_pt.co.y + end_pt.co.y)/2
            half_z = (start_pt.co.z + end_pt.co.z)/2
           
            #compute min pitch value 
            v = mathutils.Vector((start_pt.co.x, start_pt.co.y, obj["roof-0.25"])) 
            v2 = mathutils.Vector((end_pt.co.x, end_pt.co.y, obj["roof-0.25"])) 
            half_pt = mathutils.Vector((half_x, half_y, obj["roof-0.99"])) 
            #find opposite edge 
            opp_edge = bm.edges[edge_idx - len(bm.edges) + 2]
            v3 = opp_edge.verts[0].co
            v4 = opp_edge.verts[1].co 
            v3 = mathutils.Vector((v3.x, v3.y, obj["roof-0.25"])) 
            v4 = mathutils.Vector((v4.x, v4.y, obj["roof-0.25"])) 
            half_x_op = (opp_edge.verts[0].co.x + opp_edge.verts[1].co.x)/2
            half_y_op = (opp_edge.verts[0].co.y + opp_edge.verts[1].co.y)/2
            half_pt_op = mathutils.Vector((half_x_op, half_y_op, obj["roof-0.99"])) 
            
            center_pol_x= (half_x_op + half_x)/2
            center_pol_y = (half_y_op + half_y)/2
            center_pol = mathutils.Vector((center_pol_x, center_pol_y, obj["roof-0.99"])) 
            vec1_in_vplane = v - v2
            vec2_in_vplane = half_pt - v2 
            normal_vplane = vec1_in_vplane.cross(vec2_in_vplane)
            
            vec1_in_splane = v - v2
            vec2_in_splane = center_pol - v2 
            normal_splane = vec1_in_splane.cross(vec2_in_splane)

            a1, b1, c1 = normal_vplane.x, normal_vplane.y, normal_vplane.z 
            a2, b2, c2 = normal_splane.x, normal_splane.y, normal_splane.z 
        
            d = ( a1 * a2 + b1 * b2 + c1 * c2 )
            e1 = math.sqrt( a1 * a1 + b1 * b1 + c1 * c1)
            e2 = math.sqrt( a2 * a2 + b2 * b2 + c2 * c2)
            d = d / (e1 * e2)
            A = math.degrees(math.acos(d))
            if A > 180: 
                A = A - 180
                
            max_angle = 90 - A
            
            vec1_in_vplane = v3 - v4 
            vec2_in_vplane = half_pt_op - v4 
            normal_vplane = vec1_in_vplane.cross(vec2_in_vplane)
            
            vec1_in_splane = v3 - v4
            vec2_in_splane = center_pol - v4
            normal_splane = vec1_in_splane.cross(vec2_in_splane)
          
            a1, b1, c1 = normal_vplane.x, normal_vplane.y, normal_vplane.z 
            a2, b2, c2 = normal_splane.x, normal_splane.y, normal_splane.z 
        
            d = ( a1 * a2 + b1 * b2 + c1 * c2 )
            e1 = math.sqrt( a1 * a1 + b1 * b1 + c1 * c1)
            e2 = math.sqrt( a2 * a2 + b2 * b2 + c2 * c2)
            d = d / (e1 * e2)
            B = math.degrees(math.acos(d))
            if B > 180: 
                B = B - 180
        
            max_angle_2 = 90 - B
           
            if max_angle > max_angle_2: 
                the_one = int(max_angle) +1
            else: 
                the_one = int(max_angle_2) +1
            
            try: 
                if the_one < 60: 
                    the_one = 60
                param_dict["roof_pitch1"] = random.choice(range(the_one,90,1))
                param_dict["roof_pitch2"] = param_dict["roof_pitch1"]
            except: 
                param_dict["roof_pitch1"] = 90
                param_dict["roof_pitch2"] = param_dict["roof_pitch1"]
            
            param_dict["edge_midpt"]= (half_x, half_y, half_z) 
            mesh.update()
            bpy.ops.object.mode_set(mode="OBJECT")
           
        
        elif roof_type == "complex_hip": 
            obj["is_flat"] = 0 
            rulepath=rulefiles_path + 'complex_hip_roof'
            param_dict["height"] = obj["roof-0.25"] - obj["ground-0.2"]
            param_dict["top_height"] = obj["roof-0.99"] - obj["ground-0.2"]
        
        elif roof_type == "mansard": 
            obj["is_flat"] = 0 
            rulepath=rulefiles_path + 'mansard_roof'
            param_dict["height"] = obj["roof-0.25"] - obj["ground-0.2"]
            param_dict["roof_height"] = obj["roof-0.99"] - obj["roof-0.25"]
            
            #for inset we need to check that 2*inset < shortest_edge_len otherwise there will be some self-intersection in the output roof 
            max_inset = 2 
            min_inset = 1
            #find shortest edge 
            bpy.ops.object.mode_set(mode="EDIT")
            bm=bmesh.from_edit_mesh(mesh)
            shortest_len = math.inf 
            for e in bm.edges:
                if e.calc_length() < shortest_len: 
                    shortest_len = e.calc_length()
            bpy.ops.object.mode_set(mode="OBJECT")
            
            if shortest_len < 2*max_inset: 
                max_inset = (round(shortest_len,1) * 10)/2
                if shortest_len < 2*min_inset or int(max_inset) == int(min_inset*10): 
                    min_inset = 0.1 *10 
                else:
                    min_inset = min_inset *10
            else: 
                max_inset = max_inset *10
                min_inset = min_inset *10
            
            try: 
                inset1 = random.choice(range(int(min_inset), int(max_inset), 1))
            except: 
                inset1 = max_inset 
           
            param_dict["inset1"] = inset1/10
            
        elif roof_type == "dutch": 
            obj["is_flat"] = 0 
            rulepath=rulefiles_path + 'dutch_roof'
            #building height is difference between z coordinates of the top and of the ground
            param_dict["height"] = obj["roof-0.25"] - obj["ground-0.2"]
            param_dict["top_height"] = obj["roof-0.99"] - obj["ground-0.2"]
            param_dict["roof_pitch"] = 70 #random.choice(range(80,90,1))
            step_height = random.choice(range(5, 7, 1))
            if (obj["roof-0.99"] -  obj["roof-0.25"])/ 3 < step_height/10: 
                step_height = (obj["roof-0.99"] -  obj["roof-0.25"])/ 3 - 0.1
            else: 
                step_height = step_height/10
            param_dict["step_height"] = step_height
            
            param_dict["facade_thickness"] = 0.5
            world_matrix = obj.matrix_world 
            
            bpy.ops.object.mode_set(mode="EDIT")
            bm=bmesh.from_edit_mesh(mesh)
            bm.edges.ensure_lookup_table()

            #find an edge for building the facade based on neighbouring polygons 
            edge_idx = find_facade_edge_idx(obj, bm, world_matrix, kd, obj_to_edge_dict, obj_to_facade) 
            
            start_pt = bm.edges[edge_idx].verts[0]
            end_pt = bm.edges[edge_idx].verts[1] 
            current_edge = bm.edges[edge_idx]
            half_x = (start_pt.co.x + end_pt.co.x)/2
            half_y = (start_pt.co.y + end_pt.co.y)/2
            half_z = (start_pt.co.z + end_pt.co.z)/2
        
            #compute min pitch value 
            v = mathutils.Vector((start_pt.co.x, start_pt.co.y, obj["roof-0.25"])) 
            v2 = mathutils.Vector((end_pt.co.x, end_pt.co.y, obj["roof-0.25"])) 
            half_pt = mathutils.Vector((half_x, half_y, obj["roof-0.99"])) 
            
            #find opposite edge 
            opp_edge = bm.edges[edge_idx - len(bm.edges) + 2]
            v3 = opp_edge.verts[0].co
            v4 = opp_edge.verts[1].co 
            v3 = mathutils.Vector((v3.x, v3.y, obj["roof-0.25"])) 
            v4 = mathutils.Vector((v4.x, v4.y, obj["roof-0.25"])) 
            half_x_op = (opp_edge.verts[0].co.x + opp_edge.verts[1].co.x)/2
            half_y_op = (opp_edge.verts[0].co.y + opp_edge.verts[1].co.y)/2
            half_pt_op = mathutils.Vector((half_x_op, half_y_op, obj["roof-0.99"])) 
            
            center_pol_x= (half_x_op + half_x)/2
            center_pol_y = (half_y_op + half_y)/2
            center_pol = mathutils.Vector((center_pol_x, center_pol_y, obj["roof-0.99"])) 
            
            vec1_in_vplane = v3 - v4 
            vec2_in_vplane = half_pt_op - v4 
            normal_vplane = vec1_in_vplane.cross(vec2_in_vplane)
            
            vec1_in_splane = v3 - v4
            vec2_in_splane = center_pol - v4
            normal_splane = vec1_in_splane.cross(vec2_in_splane)
          
            a1, b1, c1 = normal_vplane.x, normal_vplane.y, normal_vplane.z 
            a2, b2, c2 = normal_splane.x, normal_splane.y, normal_splane.z 
        
            d = ( a1 * a2 + b1 * b2 + c1 * c2 )
            e1 = math.sqrt( a1 * a1 + b1 * b1 + c1 * c1)
            e2 = math.sqrt( a2 * a2 + b2 * b2 + c2 * c2)
            d = d / (e1 * e2)
            B = math.degrees(math.acos(d))
            if B > 180: 
                B = B - 180
        
            max_angle_2 = int(90 - B)+1
            
            try: 
                if max_angle_2 < 60: 
                    max_angle_2 = 60
                param_dict["roof_pitch"] = random.choice(range(max_angle_2,90,1))
            except: 
                param_dict["roof_pitch"] = 90 
                
            
            #shorten the input fp with the facade width value 
            #otherwise output 3D model footprint is bigger than input fp 
            pts=[start_pt, end_pt]
            directions_vec=[]

            #find adjacent edges to facade edge and translate the vertices towards the interior of the mesh 
            for pt in pts: 
                #get edges ajdacent to the vertex 
                for e in pt.link_edges: 
                    if e != current_edge: 
                        if e.verts[0] != pt: 
                            dir = (e.verts[0].co - pt.co).normalized()
                        else: 
                            dir = (e.verts[1].co - pt.co).normalized()
                        directions_vec.append(dir)     
        
            j=0
            for dir in directions_vec:
                translate_direction = dir * param_dict["facade_thickness"]
                #bmesh.ops.translate(bm, vec = translate_direction, verts=[pts[j]])
                j=j+1
            
            #add the mid point of the facade edge to the dictionary 
            half_x = (start_pt.co.x + end_pt.co.x)/2
            half_y = (start_pt.co.y + end_pt.co.y)/2
            half_z = (start_pt.co.z + end_pt.co.z)/2
            param_dict["edge_midpt"]= (half_x, half_y, half_z) 
           
            mesh.update()
            bpy.ops.object.mode_set(mode="OBJECT")
            
        elif roof_type == "gable": 
            obj["is_flat"] = 0 
            rulepath=rulefiles_path + 'gable_roof'
            param_dict["height"] = obj["roof-0.25"] - obj["ground-0.2"]
            param_dict["top_height"] = obj["roof-0.99"] - obj["ground-0.2"]
          
            world_matrix = obj.matrix_world 
            bpy.ops.object.mode_set(mode="EDIT")
            bm=bmesh.from_edit_mesh(mesh)
            bm.edges.ensure_lookup_table()

            #find an edge for building the facade based on neighbouring polygons 
            edge_idx = find_facade_edge_idx(obj, bm, world_matrix, kd, obj_to_edge_dict, obj_to_facade) 
            
            start_pt = bm.edges[edge_idx].verts[0]
            end_pt = bm.edges[edge_idx].verts[1] 
            current_edge = bm.edges[edge_idx]
            
            #add the mid point of the facade edge to the dictionary 
            half_x = (start_pt.co.x + end_pt.co.x)/2
            half_y = (start_pt.co.y + end_pt.co.y)/2
            half_z = (start_pt.co.z + end_pt.co.z)/2
            param_dict["edge_midpt"]= (half_x, half_y, half_z) 
           
            mesh.update()
            bpy.ops.object.mode_set(mode="OBJECT")

        #iterate over the dictionary and add the parameters to the template rulefile
        a_file = open(rulepath + '.py' , "r")
        list_of_lines = a_file.readlines()
        a_file.close()
        #start in line 3 
        j=2
        for key in param_dict: 
            if key == "edge_midpt": 
                list_of_lines[j] = key + "= " + str(param_dict[key]) + "\n" 
            else: 
                list_of_lines[j] = key + "= param(" + str(param_dict[key]) + ")\n" 
            j=j+1

        #a new file is created each time because parameters are attached to a given python rulefile (module). So even if we change the parameters, 
        #the output building all have the same paramaters than the first ones which were used in that file. 
        rulefile_path = rulepath + str(h) + '.py'
        a_file = open(rulefile_path, "w")
        a_file.writelines(list_of_lines)
        a_file.close()

        bpy.context.scene.bcgaScript=rulefile_path
        #execute rulefile 
        bpy.ops.object.apply_pro_script()
        os.remove(rulefile_path)
        h=h+1
        
def is_on(a, b, c):
    "Return true iff point c intersects the line segment from a to b."
    # (or the degenerate case that all 3 points are coincident)
    return (collinear(a, b, c)
            and (within(a.x, c.x, b.x) if a.x != b.x else 
                 within(a.y, c.y, b.y)))

def collinear(a, b, c):
    "Return true iff a, b, and c all lie on the same line."
    dif=((b.x - a.x) * (c.y - a.y)) - ((c.x - a.x) * (b.y - a.y)) 
  
    return dif >=-0.5 and dif <= 0.5

def within(p, q, r):
    "Return true iff q is between p and r (inclusive)."
    return p <= q <= r or r <= q <= p


def prepare_for_cityjson_export(): 
    """Prepare the objects for exporting to cityjson and assign semantics to their faces 
    
    Returns:
    nothing 
    """

    # clear all materials (so that code can be rerun) 
    for material in bpy.data.materials:
        material.user_clear()
        bpy.data.materials.remove(material)

    # create one new material for each type of semantic surfaces 
    mat_wall=bpy.data.materials.new(name="WallSurface01") 
    mat_wall["type"]="WallSurface" 
    mat_ground=bpy.data.materials.new(name="GroundSurface01") 
    mat_ground["type"]="GroundSurface"
    mat_roof=bpy.data.materials.new(name="RoofSurface01") 
    mat_roof["type"]="RoofSurface" 

    # Get materials of semantics surfaces 
    #mat_wall = bpy.data.materials.get("WallSurface01") 
    idx=0
    for obj in bpy.data.collections[coll_name].objects:
        
        obj.name= "2: [LoD2] " + str(idx)
        obj["lod"]=2
        obj["type"]="Solid" 
        
        bpy.ops.object.empty_add(type='CUBE') 
        bpy.context.active_object.name = str(idx)
        bpy.context.active_object["type"]="Building" 
        idx+=1
       
        """
        bpy.context.active_object["bouwjaar"] = obj["bouwjaar"]
        bpy.context.active_object["flat_roof"] = obj["is_flat"]
        bpy.context.active_object["ground_height"] = obj["ground-0.2"]
        bpy.context.active_object["roof-0.25"] = obj["roof-0.25"]
        bpy.context.active_object["roof-0.99"] = obj["roof-0.99"]
        bpy.context.active_object["origin"] = "historical map"
        """
        # assign the materials to the object 
        for mat in [mat_wall, mat_ground, mat_roof]:
            obj.data.materials.append(mat) 
            
        mesh = obj.data 
        faces = mesh.polygons 
        ground_normal = mathutils.Vector((0,0,-1))
       
        # assign a semantic to each surface 
        for face in faces: 
            try: 
                # if the surface has vertices with an height equals to the roof height, it is a roof surface 
                # need to change this for facade stepped 
                # problem with ground surface angle 180 degrees 
                semantic_assigned=False
                
                #if obj["is_flat"] == 0: 
                roof_height = obj["roof-0.25"] 
                #else: 
                    #roof_height = obj["roof-0.75"] 
                
                for v in face.vertices: 
                    v= mesh.vertices[v]
                    if v.co.z >= roof_height - 0.2: 
                        semantic_assigned=True 
                    else: 
                        semantic_assigned= False 
                        break 
                
                if semantic_assigned == True: 
                    face.material_index = 2 #roof material has index 2 in the data.materials list 
                    continue 
                
                # if the angle betw the surface normal and the ground normal is betw -5 and 5 degrees, it is a ground surface 
                angle_with_ground = math.degrees(face.normal.angle(ground_normal))
                if angle_with_ground >= 180-0.2 : 
                    angle_with_ground = angle_with_ground - 180
                if angle_with_ground <= 5 and angle_with_ground >= -5:
                    face.material_index = 1 #ground material has index 1 
                    continue 
                
                #if the angle betw surface normal and horizontal vector with same x and y is betw -5 and 5 degrees, it is a wall surface 
                wall_normal = mathutils.Vector((face.normal.x, face.normal.y, 0))
                angle_with_wall = math.degrees(face.normal.angle(wall_normal))
                if angle_with_wall >= 180-0.2 : 
                    angle_with_wall= angle_with_wall - 180
                if angle_with_wall <= 5 and angle_with_wall >= -5:
                    face.material_index = 0 #wall material has index 0
                    continue 
            
            except: 
                pass 
            
        mesh.update()

    #add world properties about the translation betw blender crs and original crs 
    #needed when exporting to cityjson so that it can translates the coordinates back to the original position 
    bpy.context.scene.world['Axis_Origin_X_translation']=-bpy.data.scenes["Scene"]["crs x"]
    bpy.context.scene.world['Axis_Origin_Y_translation']=-bpy.data.scenes["Scene"]["crs y"]
    bpy.context.scene.world['Axis_Origin_Z_translation']=0
    #note that I need a 3D crs, i can't just use EPSG:28992!! So i can use EPSG:7415 which is Amersfoort / RD New + NAP height 
    bpy.context.scene.world['CRS']= epsg_3d #"urn:ogc:def:crs:EPSG::6190" #"urn:ogc:def:crs:EPSG::7415"



def build_kd_tree_with_footprints():
    """Build a kd tree from the mid point edges of all building footprints. Later used to identify the neighbours of a given building footprint. 

    Returns: 
    kd -- kd tree containing the mid points (and their ID) of all edges of all building footprints 
    obj_to_edge_dict -- dictionnary mapping from the id of the edges mid points (key) to the name of the object to which they belong (value)
    """

    edge_mid_points=[]
    for obj in bpy.data.collections[coll_name].objects:
        bpy.context.view_layer.objects.active = obj
        mesh = obj.data
        # Get a BMesh representation
        bm = bmesh.new()   # create an empty BMesh
        bm.from_mesh(mesh)   # fill it in from a Mesh
        bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=0.1)
        bm.to_mesh(mesh)
        bm.free()  # free and prevent further access
    
        for vert in mesh.vertices:
            new_location = vert.co
            new_location[2] = obj["ground-0.2"]
            vert.co = new_location
        mesh.update()

        for e in mesh.edges:
            pt0_global = obj.matrix_world @ mesh.vertices[e.vertices[0]].co
            pt1_global = obj.matrix_world @ mesh.vertices[e.vertices[1]].co
            half_pt_x = (pt0_global.x + pt1_global.x)/2
            half_pt_y = (pt0_global.y + pt1_global.y)/2
            half_pt_z = (pt0_global.z + pt1_global.z)/2
            edge_mid_points.append((obj.name, mathutils.Vector((half_pt_x, half_pt_y, 0)))) 
           
    obj_to_edge_midpt = dict()
    size = len(edge_mid_points)

    kd = KDTree(size)
    for i in range(0, len(edge_mid_points)): 
        kd.insert(edge_mid_points[i][1], i) 
        kd.balance()
        obj_to_edge_midpt[(i,)] = edge_mid_points[i][0]
    
    return kd, obj_to_edge_midpt
         
        
prepare_blender_scene()
kd, obj_to_edge_midpt = build_kd_tree_with_footprints() 
procedural_modelling(kd, obj_to_edge_midpt)
prepare_for_cityjson_export()
bpy.ops.cityjson.export_file("EXEC_DEFAULT", filepath = r"{}".format(argv[1]))
