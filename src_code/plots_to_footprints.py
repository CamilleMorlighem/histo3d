from skgeom.draw import draw
import functools
from shapely.geometry import LineString, MultiPoint, Point, Polygon
from shapely.ops import split
from numpy import ones,vstack
from numpy.linalg import lstsq
import shapely 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import skgeom
import math 
import poly_decomp as pdd
import geopandas as gpd
import random 
from shapely.ops import unary_union

"""
This script is used split building plots into individual building footprints 
It is used by the script obia_segmentation_full_worklow.py
"""

###################  HELPERS ################### 

def draw_skeleton(polygon, skeleton, show_time=False):
    """draw straight skeleton of given polygon 

    Keyword arguments:
    skeleton -- straight skeleton of given polygon 
    polygon -- input polygon 

    Returns:
    none
    """
    draw(polygon)
    for h in skeleton.halfedges:
        if h.is_bisector:
            p1 = h.vertex.point
            p2 = h.opposite.vertex.point
            plt.plot([p1.x(), p2.x()], [p1.y(), p2.y()], 'r-', lw=2)

    if show_time:
        for v in skeleton.vertices:
            plt.gcf().gca().add_artist(plt.Circle(
                (v.point.x(), v.point.y()),
                v.time, color='blue', fill=False))

def get_line_param_eq(line): 
    """get parametric equation of shapely linestring under the form: 
    X = oax + abx*t 
    Y = oay + aby*t

    Keyword arguments:
    line -- shapely linestring 
    
    Returns:
    oa, abx, aby -- members of line parametric equation; (abx, aby) is directional vector of the line 
    """
    oa=line.coords[0]
    abx=line.coords[-1][0]-line.coords[0][0]
    aby=line.coords[-1][1]-line.coords[0][1]
    return oa, abx, aby
    
def generate_segment_from_line(input_polygon, pt, abx, aby):
    """create a perpendicular segment to a input line, passing through a point 

    Keyword arguments:
    input_polygon -- input_polygon in which the lines and segments are created (used to make sure the created segment intersects with the input_polygon boundary) 
    pt -- point in which the perpendicular segment passes 
    abx, aby -- directional vector of the input line 
    
    Returns:
    perpendicular segment under the form of a shapely linestring 
    """ 

    #uses the extreme polygon coordinates to make sure the created segment intersects with the input_polygon boundary 
    #get the input polygon BBOX 
    minx, miny, maxx, maxy = input_polygon.bounds
    bounding_box = shapely.geometry.box(minx, miny, maxx, maxy)
    #get coordinates of point i through which the segment must pass 
    ix=pt[0]
    iy=pt[1]
    
    #if segment to be created is vertical, x is cst and will always be ix 
    if abx==0: 
        points_on_boundary_lines = [Point(ix, maxy), Point(ix, miny)]
    #elif segment is horizontal, y is cst and will always be iy 
    elif aby==0: 
        points_on_boundary_lines = [Point(minx, iy), Point(maxx, iy)]
    #segment is neither vertical nor horizontal 
    else: 
        t=(minx-ix)/abx
        y0=iy+(aby)*t
        t=(maxx-ix)/abx
        y1=iy+(aby)*t
        t=(miny-iy)/(aby)
        x0=ix+abx*t
        t=(maxy-iy)/(aby)
        x1=ix+abx*t
        points_on_boundary_lines = [Point(minx, y0), Point(maxx, y1), 
                                        Point(x0, miny), Point(x1, maxy)]
        
    #select the points which are the closest to the bounding box and create a segment with them 
    points_sorted_by_distance = sorted(points_on_boundary_lines, key=bounding_box.distance)
    extended_line = LineString(points_sorted_by_distance[:2])
    return extended_line


def add_footprint(footprint, min_area, max_perim, footprints): 
    """add previously created footprint to list of footprints already created to subdivide the input building plot. Check if current footprint meets area requirement to be added to footprints list

    Keyword arguments:
    footprint -- newly created footprint to be added in the form of a shapely polygon 
    min_area -- minimum area of footprint to be considered as valid footprint 
    max_perim -- max perimeter allowed for footprint to be valid 
    footprints -- list of footprints already created (and thus valid)
   
    Returns:
    none (expand footprints list)
    """ 

    if footprint.is_valid and footprint.area >= min_area and footprint.length <= max_perim: 
        #dont include triangles 
        if len(footprint.exterior.coords)-1 != 3: 
            #remove duplicate vertices 
            footprint = footprint.simplify(0.01)
            footprint=shapely.geometry.polygon.orient(footprint, sign=1.0)
            footprints.append(footprint)
           

def add_end_footprints(input_poly, extended_line, endpt, s_lon_side, e_lon_side, reverse, min_area, max_perim, footprints): 
    """add footprint at the end/beginning of a building row (used in case 2 where the building plot is considered to be a unique row of building footprints)

    Keyword arguments:
    input_poly -- input building plot (unique row of building footprints) under the form of shapely polygon  
    extended_line -- last perpendicular segment created to generate the last footprint in the row 
    endpt -- shapely point either equal to s_lon_side or e_lon_side, depending whether we are at the first building of the row or last building of the row 
    s_lon_side -- starting shapely point of the line segment of the building plot along which we are walking and creating footprints 
    e_lon_side -- end shapely point of the line segment of the building plot along which we are walking and creating footprints 
    reverse -- binary (T or F) indicating on which side of the building footprint we shorten it for adding randomness in the building depth 
    min_area -- minimum area of footprint to be considered as valid footprint 
    max_perim -- max perimeter allowed for footprint to be valid 
    footprints -- list of footprints already created (and thus valid)
   
    Returns:
    none (expand footprints list)
    """ 

    #split the input_polygon with the perpendicular segment with endpoints pt and intersection 
    splitted_poly=shapely.ops.split(input_poly, extended_line)
    #the first footprint/last fp will be the polygon that is the closest to the first/last endpoint of the longest_side segment 
    shortest_d = math.inf
    for poli in splitted_poly.geoms:
        d = poli.distance(endpt)
        if d < shortest_d: 
            shortest_d = d
            pol = poli 

    if shortest_d < math.inf: 
        footprint=shapely.geometry.polygon.orient(pol, sign=1.0)
        #change facade depth of the footprint IF it has 4 vertices 
        if len(footprint.exterior.coords)>4: 
            fp_ext_coords = footprint.exterior.coords[0:-1]
            for coord in range(0, len(footprint.exterior.coords)-1):
                bdary_line=LineString([footprint.exterior.coords[coord],footprint.exterior.coords[coord+1]])
                s_bd_line, e_bd_line = bdary_line.boundary
                #once we found the linestring of the footprint which lies on the splitting segment 
                if (s_bd_line.distance(e_lon_side) < 0.1  or s_bd_line.distance(s_lon_side) < 0.1 or is_on(s_lon_side, e_lon_side, s_bd_line)) and (e_bd_line.distance(s_lon_side) < 0.1 or is_on(s_lon_side, e_lon_side, e_bd_line) or e_bd_line.distance(e_lon_side) < 0.1 ):  
                    #we get the two adjacents linestrings 
                    prev_vtx = fp_ext_coords[coord-1]
                    i_pt_1 = coord 
                    i_r_pt_1 = coord-1

                    if (coord+1) == len(fp_ext_coords): 
                        next_vtx = fp_ext_coords[1]
                        i_pt_2 = 0
                        i_r_pt_2 = 1
                    elif (coord+1) == len(fp_ext_coords)-1: 
                        next_vtx = fp_ext_coords[0]
                        i_pt_2 = -1
                        i_r_pt_2 = 0
                    else: 
                        next_vtx = fp_ext_coords[coord + 2]
                        i_pt_2 = coord + 1
                        i_r_pt_2 = coord + 2
                    if not reverse: 
                        depth_edge_1 = LineString([s_bd_line, prev_vtx])
                        depth_edge_2 = LineString([e_bd_line, next_vtx])
                    else: 
                        depth_edge_1 = LineString([prev_vtx, s_bd_line])
                        depth_edge_2 = LineString([next_vtx, e_bd_line])
                    if depth_edge_1.length < depth_edge_2.length: 
                        shortest_len = depth_edge_1.length
                    else: 
                        shortest_len = depth_edge_2.length
                    facade_depth_random = round(random.uniform(0, shortest_len/3),1)
                    #we move the two points of the linestring on the splitting side to add randomness 
                    pt_1 = depth_edge_1.interpolate(facade_depth_random)
                    pt_2 = depth_edge_2.interpolate(facade_depth_random)
                    if not reverse: 
                        fp_ext_coords[i_pt_1] = pt_1 
                        fp_ext_coords[i_pt_2] = pt_2 
                    else: 
                        fp_ext_coords[i_r_pt_1] = pt_1 
                        fp_ext_coords[i_r_pt_2] = pt_2 
                    
                    fp_ext_coords.append(fp_ext_coords[0])
                    footprint = Polygon(fp_ext_coords) 
        
        add_footprint(footprint, min_area, max_perim, footprints)


def generate_facades_len(line, facade_len_min, facade_len_max, add_end_pts = True): 
    """split a line into different segments of various lengths representing the facade edge of building footprints

    Keyword arguments:
    line -- shapely linestring belonging to a input building plot along which we are going to walk to create different building footprints 
    facade_len_min -- minimum facade lenght allowed 
    facade_len_max -- max facade length allowed 
    add_end_pts -- binary (T or F) telling whether we are adding the end points of the input line to the list of points returned 
   
    Returns:
    a list of points [a,b,c,d,e] where ab is facade 1, bc is facade 2, cd is facade 3, de is facade 4 
    """ 

    splitter = []
    accumulative_line_length = 0 
    if line.length > 2*facade_len_min: 
        while True: #(line.length - accumulative_line_length) > facade_len_min: 
            facade_len_random = round(random.uniform(facade_len_min, facade_len_max), 1)
            pt = line.interpolate(accumulative_line_length + facade_len_random)

            if (line.length - accumulative_line_length) > facade_len_random:  
                splitter.append(pt)
            elif (line.length - accumulative_line_length) >= facade_len_min and (line.length - accumulative_line_length) <= facade_len_max: 
                break 
            elif prev_len + (line.length - accumulative_line_length) > 2 * facade_len_min:
                splitter.pop()
                pt = line.interpolate(accumulative_line_length - prev_len + (prev_len + (line.length - accumulative_line_length))/2)
                splitter.append(pt)
                break 
            elif line.length - accumulative_line_length < facade_len_min: 
                splitter.pop() 
                break 

            prev_len = facade_len_random
            accumulative_line_length += facade_len_random

    splitter = MultiPoint(splitter)
    s,e = line.boundary

    if add_end_pts: 
        segments_pts=[(s.x, s.y)]
        segments_pts=segments_pts+[(p.x, p.y) for p in splitter]
        segments_pts.append((e.x, e.y))
    else: 
        segments_pts=[(p.x, p.y) for p in splitter]

    return segments_pts

def is_on(a, b, c):
    "Return true iff point c intersects the line segment from a to b."
    # (or the degenerate case that all 3 points are coincident)
    return (collinear(a, b, c)
            and (within(a.x, c.x, b.x) if a.x != b.x else 
                 within(a.y, c.y, b.y)))

def collinear(a, b, c):
    "Return true iff a, b, and c all lie on the same line."
    dif=((b.x - a.x) * (c.y - a.y)) - ((c.x - a.x) * (b.y - a.y)) 
    return dif >=-0.002 and dif <= 0.002

def within(p, q, r):
    "Return true iff q is between p and r (inclusive)."
    return p <= q <= r or r <= q <= p


################### CASES ################### 

def case_2(input_poly, min_area, max_perim, facade_depth, facade_len_min, facade_len_max, splitting_side = None, reverse = False): 
    """Case 2 subdivides a building plot into individual building footprints, considering that the input building plot is only made of one row of building footprints 

    Keyword arguments:
    input_poly -- input building plot (shapely polygon)
    min_area -- minimum area of footprint to be considered as valid footprint 
    max_perim -- max perimeter allowed for footprint to be valid 
    facade_len_min -- minimum facade lenght allowed 
    facade_len_max -- max facade length allowed 
    facade_depth -- do not consider; not used here, used in deprecated version of this function 
    splitting_side -- side of the building plot along which to walk and creates bf (default is none, in this case the longest side of the building plot is chosen as splitting side)
    reverse -- binary (T or F) indicating on which side of the building footprint we shorten it for adding randomness in the building depth 
   
    Returns:
    list of individual footprints as shapely polygons 
    """ 

    #input poly must be oriented ccw 
    input_poly=shapely.geometry.polygon.orient(input_poly, sign=1.0)
    if splitting_side == None: 
        #take longest side of the polygon 
        longest_len= 0 
        for coord in range(0, len(input_poly.exterior.coords)-1):
            ext_linestr=LineString([input_poly.exterior.coords[coord],input_poly.exterior.coords[coord+1] ])
            if ext_linestr.length>longest_len: 
                longest_len=ext_linestr.length 
                longest_side=ext_linestr
        splitting_side = longest_side
    #splitting side must have the same direction as the one it has in the input_polygon 
    else: 
        
        for coord in range(0, len(input_poly.exterior.coords)-1):
            ext_linestr=LineString([input_poly.exterior.coords[coord],input_poly.exterior.coords[coord+1]])
            s_lon_side, e_lon_side = splitting_side.boundary
            s_op_side, e_op_side = ext_linestr.boundary
            if ext_linestr.centroid.distance(splitting_side.centroid) < 0.1 and s_lon_side.distance(s_op_side) > 0.1: 
                splitting_side = LineString([e_lon_side, s_lon_side])
               
    s_lon_side, e_lon_side = splitting_side.boundary
    segments_pts=generate_facades_len(splitting_side, facade_len_min, facade_len_max, add_end_pts = False)

    #get longest_side parametric equation 
    oa, abx, aby = get_line_param_eq(splitting_side)
    
    #for each "facade" along the longest_side: 
    footprints=[]
    shortest_len = None 
    for pt in segments_pts: 
        #generate perpendicular line segment that crosses the polygon exterior boundary 
        extended_line=generate_segment_from_line(input_poly, pt, aby, -abx)
        #for each linestring of the input polygon, check the intersection betw him and the perpendicular segment 
        #you will get two intersections: the one on the longest_side and the one on the opposide side of the longest_side 
        for coord in range(0, len(input_poly.exterior.coords)-1):
            ext_linestr=LineString([input_poly.exterior.coords[coord],input_poly.exterior.coords[coord+1]])
            #makes sure that the linestring that intersects is not the longest_side itself 
            if ext_linestr.intersects(extended_line)==True and ext_linestr.centroid.distance(splitting_side.centroid) > 0.01: #ext_linestr!=splitting_side and ext_linestr != reverse_splitting_side : 
                intersection=extended_line.intersection(ext_linestr)
                opposite_side=ext_linestr
                
        try: 
            #get end points of the longest_side and opp side 
            s_lon_side, e_lon_side = splitting_side.boundary
            s_op_side, e_op_side = opposite_side.boundary
        except: 
            add_footprint(input_poly, min_area, max_perim, footprints)
            return footprints 
       
        #if we are working on the first footprint of the plot 
        if pt==segments_pts[0]: 
            add_end_footprints(input_poly, extended_line, s_lon_side, s_lon_side, e_lon_side, reverse, min_area, max_perim, footprints)
            #store needed coordinates for next footprint: 
            pt_minus_1=pt
            intersection_minus_1=intersection
            prev_op_side = opposite_side
        else: 
            if not reverse: 
                depth_edge_1 = LineString([pt, intersection])
                depth_edge_2 = LineString([pt_minus_1, intersection_minus_1])
            else: 
                depth_edge_1 = LineString([intersection, pt])
                depth_edge_2 = LineString([intersection_minus_1, pt_minus_1])
           
            if shortest_len==None or shortest_len/3 > depth_edge_1.length/2 or shortest_len/3 > depth_edge_2.length/2: 
                if depth_edge_1.length < depth_edge_2.length: 
                    shortest_len = depth_edge_1.length
                else: 
                    shortest_len = depth_edge_2.length
            facade_depth_random = round(random.uniform(0, shortest_len/3),1)
            pt_1 = depth_edge_1.interpolate(facade_depth_random, normalized=False)
            pt_2 = depth_edge_2.interpolate(facade_depth_random, normalized = False)
            if not reverse: 
                pt_11 = pt_1 
                pt_22 = pt_2 
                int_1 = intersection 
                int_2 = intersection_minus_1 
            else: 
                pt_11 = pt 
                pt_22 = pt_minus_1 
                int_1 = pt_1 
                int_2 = pt_2 

            #if the opposite side and the previous ones are not the same, we have to include at least one corner of the input_polygon in the footprint 
            if opposite_side!=prev_op_side: 
                s_prev_op_side, e_prev_op_side = prev_op_side.boundary 
                #if endpoints of the exterior linestrings are the same, there is one corner to add 
                if e_op_side==s_prev_op_side: 
                    footprint= Polygon([pt_11, intersection, e_op_side, intersection_minus_1, pt_22])
                #if endpoints are not the same, we add two corners (note that they could be more but we don't take all of them into account)
                else: 
                    footprint= Polygon([pt_11, intersection, e_op_side, s_prev_op_side, intersection_minus_1, pt_22])

                add_footprint(footprint, min_area, max_perim, footprints)

            else: 
                footprint= Polygon([pt_11, int_1, int_2, pt_22])
                add_footprint(footprint, min_area, max_perim, footprints)

            pt_minus_1=pt
            intersection_minus_1=intersection
            prev_op_side = opposite_side 

        #elif we are working on the second to last and on the last footprints of the plot 
        if pt==segments_pts[-1]: 
            add_end_footprints(input_poly, extended_line, e_lon_side, s_lon_side, e_lon_side, reverse, min_area, max_perim, footprints)

    return footprints 
            

def case_3(input_poly, min_area, max_perim,facade_depth, facade_len_min, facade_len_max): 
    """Case 3 subdivides a building plot into individual building footprints, considering that the input building plot is made of TWO rows of building footprints 

    Keyword arguments:
    input_poly -- input building plot (shapely polygon)
    min_area -- minimum area of footprint to be considered as valid footprint 
    max_perim -- max perimeter allowed for footprint to be valid 
    facade_len_min -- minimum facade lenght allowed 
    facade_len_max -- max facade length allowed 
    facade_depth -- min facade depth allowed 

    Returns:
    list of individual footprints as shapely polygons 
    """ 

    #take longest_median
    longest_len=0

    #if polygon has an even number of sides 
    if (len(input_poly.exterior.coords)-1)%2==0: 
        n=int((len(input_poly.exterior.coords)-1)/2)
            
        #for each linestring in the exterior bdary of the input polygon 
        for coord in range(0,len(input_poly.exterior.coords)-1): 
            #get half point 
            if coord==len(input_poly.exterior.coords)-2: 
                half_pt_x=(input_poly.exterior.coords[0][0]+input_poly.exterior.coords[coord][0])/2
                half_pt_y=(input_poly.exterior.coords[0][1]+input_poly.exterior.coords[coord][1])/2

            else: 
                half_pt_x=(input_poly.exterior.coords[coord+1][0]+input_poly.exterior.coords[coord][0])/2
                half_pt_y=(input_poly.exterior.coords[coord+1][1]+input_poly.exterior.coords[coord][1])/2
            
            #find opposite half pt 
            if coord+n>=len(input_poly.exterior.coords)-1: 
                cn=coord+n-(len(input_poly.exterior.coords)-1)
                cn1=coord+n+1-(len(input_poly.exterior.coords)-1)
            else: 
                cn=coord+n
                cn1=coord+n+1
            
            opp_half_x=(input_poly.exterior.coords[cn1][0]+input_poly.exterior.coords[cn][0])/2
            opp_half_y=(input_poly.exterior.coords[cn1][1]+input_poly.exterior.coords[cn][1])/2
   
            #create median line 
            median_line=LineString([(half_pt_x, half_pt_y), (opp_half_x, opp_half_y)])
           
            if median_line.length>longest_len: 
                longest_len=median_line.length
                longest_median=median_line 

    else: 
        n=(1+(len(input_poly.exterior.coords))-1)//2
        #for each linestring in the exterior bdary of the input polygon 
        for coord in range(0,len(input_poly.exterior.coords)-1): 
            #get half point 
            if coord==len(input_poly.exterior.coords)-2: 
                half_pt_x=(input_poly.exterior.coords[0][0]+input_poly.exterior.coords[coord][0])/2
                half_pt_y=(input_poly.exterior.coords[0][1]+input_poly.exterior.coords[coord][1])/2

            else: 
                half_pt_x=(input_poly.exterior.coords[coord+1][0]+input_poly.exterior.coords[coord][0])/2
                half_pt_y=(input_poly.exterior.coords[coord+1][1]+input_poly.exterior.coords[coord][1])/2
            
            #find opposite polygon vertex 
            if coord+n>=len(input_poly.exterior.coords)-1: 
                cn=coord+n-(len(input_poly.exterior.coords)-1)
                cn1=coord+n+1-(len(input_poly.exterior.coords)-1)
            else: 
                cn=coord+n
            opp_vertex=input_poly.exterior.coords[cn]

            #create median line 
            median_line=LineString([(half_pt_x, half_pt_y), opp_vertex])
           
            if median_line.length>longest_len: 
                longest_len=median_line.length
                longest_median=median_line 
    
    #get longest_median parametric equation 
    oa, abx, aby = get_line_param_eq(longest_median)
    #extend the median 
    median_extended=generate_segment_from_line(input_poly, oa, abx, aby)

    #split the input_poly with the extended median 
    splitted_poly=shapely.ops.split(input_poly, median_extended)
    footprints=[]
    for pol in splitted_poly.geoms: 
        footprints=footprints+case_2(pol, min_area, max_perim, facade_depth, facade_len_min, facade_len_max, longest_median)
    return footprints 


def case_5(input_poly, avg_area, min_area, max_perim, facade_len_min, facade_len_max, facade_depth, max_facade_depth):
    """Case 5 subdivides a CONCAVE building plot into individual building footprints, using convex decomposition 

    Keyword arguments:
    input_poly -- input building plot (shapely polygon)
    avg_area -- average building footprint area 
    min_area -- minimum area of footprint to be considered as valid footprint 
    max_perim -- max perimeter allowed for footprint to be valid 
    facade_len_min -- minimum facade lenght allowed 
    facade_len_max -- max facade length allowed 
    facade_depth -- min facade depth allowed 
    max_facade_depth -- max facade depth allowed 
    
    Returns:
    list of individual footprints as shapely polygons 
    """ 

    #translate shapely polygon to list of vertices 
    polygon = []
    x1,y1 = input_poly.exterior.xy
    p=list(zip(list(x1), list(y1)))[:-1]
    for x in p: 
        polygon.append(list(x)) 

    #decompose concave polygon into convex polygons 
    convex_poly=pdd.polygonQuickDecomp(polygon)
    footprints=[]
    
    #for each convex sub polygon, call the split_into_footprints function 
    summed_area = 0
    to_merge = []
    for cp in convex_poly: 
        pol=Polygon(cp)
        summed_area = pol.area + summed_area
        to_merge.append(pol)
    merged = unary_union(to_merge)
    if summed_area < input_poly.area - 0.1:
        polpol = input_poly.difference(merged).simplify(0.05, True)
        to_merge.append(polpol)

    for pol in to_merge: #convex_poly: 
        x,y = pol.exterior.coords.xy
       
        longest_len= 0 
        for coord in range(0, len(pol.exterior.coords)-1):
            ext_linestr=LineString([pol.exterior.coords[coord],pol.exterior.coords[coord+1] ])
            if ext_linestr.length>longest_len: 
                longest_len=ext_linestr.length 
                longest_side=ext_linestr
        reverse = False 
        for coord in range(0, len(input_poly.exterior.coords)-1):
            ext_linestr=LineString([input_poly.exterior.coords[coord],input_poly.exterior.coords[coord+1] ])
            if ext_linestr.centroid.distance(longest_side.centroid)<0.1: 
                reverse = True 
                break 
        
        fp=split_into_footprints(pol, avg_area, min_area, max_perim, facade_len_min, facade_len_max, facade_depth, max_facade_depth, reverse, False)
        footprints= footprints+fp 
    
    footprints_filtered=[]
    
    #now we need to check that all footprints created have an access from the interior (not surrounded by other footprints on all sides)
    #for each footprint 
    for fp in footprints: 
        #for each linestring of the footprint 
        searching_intersection = True 
        for cd in range(0, len(fp.exterior.coords)-1): 
            #get half pt 
            half_pt_x=(fp.exterior.coords[cd+1][0]+ fp.exterior.coords[cd][0])/2
            half_pt_y=(fp.exterior.coords[cd+1][1]+ fp.exterior.coords[cd][1])/2
            
            #check if halp point is on the exterior bdary of the input poly 
            #if it is not the case, don't keep it 
            for coord in range(0, len(input_poly.exterior.coords)-1):
                ext_linestr=LineString([input_poly.exterior.coords[coord],input_poly.exterior.coords[coord+1]])
                s,e=ext_linestr.boundary
                hp=Point([(half_pt_x,half_pt_y)])
                if collinear(s,e,hp): 
                    if is_on(s,e,hp): 
                        footprints_filtered.append(fp)
                        searching_intersection=False 
                        break 
            if searching_intersection==False: 
                break 
    
    return footprints_filtered


def case_4(input_poly, min_area, max_perim, facade_len_min, facade_len_max, max_facade_depth, min_facade_depth, max_area = 300):
    """Case 4 subdivides a building plot into individual bf, considering the building plot contains an interior courtyard around the bf are distributed 

    Keyword arguments:
    input_poly -- input building plot (shapely polygon)
    min_area -- minimum area of footprint to be considered as valid footprint 
    max_perim -- max perimeter allowed for footprint to be valid 
    facade_len_min -- minimum facade lenght allowed 
    facade_len_max -- max facade length allowed 
    facade_depth -- min facade depth allowed 
    max_facade_depth -- max facade depth allowed 
    max_area -- max bf area 
    
    Returns:
    list of individual footprints as shapely polygons 
    """ 
    input_poly=shapely.geometry.polygon.orient(input_poly, sign=1.0)
    x1,y1 = input_poly.exterior.xy

    p=list(zip(list(x1), list(y1)))[:-1]
    sg_input_poly=skgeom.Polygon(np.array(p))

    #get interior straight skeleton 
    skel=skgeom.skeleton.create_interior_straight_skeleton(sg_input_poly)
    #get offset polygon (negative buffer)
    #poly_with_holes = functools.reduce(lambda a, b: skgeom.boolean_set.difference(a, b)[0], skel.offset_polygons(facade_depth), sg_input_poly)
    try: 
        pts=np.array(skel.offset_polygons(max_facade_depth)[0].coords) 
        facade_depth = max_facade_depth
    except: 
        try: 
            pts=np.array(skel.offset_polygons(min_facade_depth)[0].coords) 
            facade_depth = min_facade_depth
        except: 
            footprints = case_2(input_poly, min_area, max_perim, min_facade_depth, facade_len_min, facade_len_max)
            return footprints 

 
    pts_lst=pts.tolist()
    pts_lst=[tuple(l) for l in pts_lst]
    pts_lst.append(pts_lst[0])

    #prepare empty list to store the coordinates of the footprints on the corners 
    pts_on_corners=[]
    footprints=[]

    #walk along the interior polygon (offset) boundary 
    #for each linestring in the boundary 
    for i in range(0, len(pts_lst)-1): 
        line = LineString([pts_lst[i], pts_lst[i+1]])
        oa, abx, aby= get_line_param_eq(line)

        segments_pts = generate_facades_len(line, facade_len_min, facade_len_max, add_end_pts = True)
        segments_pts_copy=segments_pts.copy()

        #i think that u need to change it here maybe we can put one facade alors que ici ça devrait etre qd on c mm pas mettre une facade dc qd line lenght < facade len min 
        if line.length > facade_len_min: #len(segments_pts_copy) > 2: #n>1: 
            #prepare empty list to store the intersection points between the exterior boundary and perpendicular transects  
            pts_on_ext_boundary=[]
            index=0
            #for each point along the interior linestring:  
            for pt in segments_pts: 
                #get point coordinates 
                ix=pt[0]
                iy=pt[1]
                #get line equation of the perpendicular line to the interior linestring passing through pt 
                #line equation of perpendicular in starting point i of the segment is X=ix + aby*t and Y = iy-abx*t 
                extended_line=generate_segment_from_line(input_poly, pt, aby, -abx)
               
                #checking the intersection between the input polygon and the extended line segment (which is thus a line)
                intersection_pts=[]
                ext_linestr_lst={}
                
                for coord in range(0, len(input_poly.exterior.coords)-1):
                    ext_linestr=LineString([input_poly.exterior.coords[coord],input_poly.exterior.coords[coord+1] ])
                    if ext_linestr.intersects(extended_line)==True: 
                        intersection=extended_line.intersection(ext_linestr)
                        intersection_pts.append(intersection)
                        ext_linestr_lst[(intersection.x, intersection.y)]= ext_linestr
                
                #the line will intersect twice the polygon, keep the closest intersection 
                closest_pt=sorted(intersection_pts, key=line.distance)[0]
                
                #get d1, directional vector of the exterior line segment 
                s,e=ext_linestr_lst[(closest_pt.x, closest_pt.y)].boundary
                d1x=(e.x-s.x)
                d1y=(e.y-s.y)
                #get d2, directional vector of the interior line segment 
                s,e=line.boundary
                d2x=(e.x-s.x)
                d2y=(e.y-s.y)
                product_of_the_length=math.sqrt(d1x**2+d1y**2)*math.sqrt(d2x**2+d2y**2)
                dot_product=d1x*d2x+d2y*d1y
                #if the exterior line segment and the interior line segment are parallel, the dot product of the directional vector
                #equals the product of their length 
                #if they are parallel, continue. Otherwise we don't create the footprint at this point because it will likely intersects with other footprints 
                if dot_product-product_of_the_length<=0.002 and dot_product-product_of_the_length>=-0.002: 
                    pts_on_ext_boundary.append((closest_pt.x, closest_pt.y))

                    #if we are working on the first footprint 
                    #add one corner point 
                    if index == 0: 
                        pts_on_corners.append("s")
                        pts_on_corners.append((closest_pt.x, closest_pt.y))
                        s,e=ext_linestr_lst[(closest_pt.x, closest_pt.y)].boundary
                        pts_on_corners.append((s.x, s.y))

                    #if we working on the last footprint 
                    #add 3 corner points 
                    if index == len(segments_pts)-1:
                        pts_on_corners.append("e")
                        s,e=ext_linestr_lst[(closest_pt.x, closest_pt.y)].boundary
                        closest_corner=e
                        pts_on_corners.append((closest_corner.x, closest_corner.y))
                        pts_on_corners.append((closest_pt.x, closest_pt.y))
                        pts_on_corners.append((ix, iy))
                
                #if interior and exterior linestr are not paralell, don't create the footprints 
                else: 
                    segments_pts_copy.remove((ix,iy)) 
                index=index+1

            #generate the footprints (except the corner footprints)
            for i in range(0,len(pts_on_ext_boundary)-1): 
                depth_edge_1 = LineString([pts_on_ext_boundary[i], segments_pts_copy[i]])
                depth_edge_2 = LineString([pts_on_ext_boundary[i+1], segments_pts_copy[i+1]])
                facade_depth_random = round(random.uniform(2* facade_depth/3, facade_depth), 1)
                depth_pt_1 = depth_edge_1.interpolate(facade_depth_random)
                depth_pt_2 = depth_edge_2.interpolate(facade_depth_random)
              
                footprint_pts=[depth_pt_1, depth_pt_2, pts_on_ext_boundary[i+1], pts_on_ext_boundary[i]]
                footprint_poly=Polygon(footprint_pts)
                add_footprint(footprint_poly, min_area, max_perim, footprints)
                x,y = footprint_poly.exterior.coords.xy
               
        #if we can't put a footprint with facade len along the interior linestring, we'll merge it with one of its adjacent footprint 
        else: #n<=1
            #take mid point of the interior segment so that i can create a perpendicular to the segment at the mid point and by getting the intersection 
            #with the exterior footprint i can access the exterior linestring 
            ix= (segments_pts[-1][0]+segments_pts[0][0])/2
            iy= (segments_pts[-1][1]+segments_pts[0][1])/2

            #generate perpendicular segment 
            extended_line=generate_segment_from_line(input_poly, (ix, iy), aby, -abx)
            
            #get intersection betw perpendicular segment and exterior linestrings 
            intersection_pts=[]
            ext_linestr_lst={}
            for coord in range(0, len(input_poly.exterior.coords)-1):
                ext_linestr=LineString([input_poly.exterior.coords[coord],input_poly.exterior.coords[coord+1] ])
                if ext_linestr.intersects(extended_line)==True: 
                    intersection=extended_line.intersection(ext_linestr)
                    intersection_pts.append(intersection)
                    ext_linestr_lst[(intersection.x, intersection.y)]= ext_linestr
            
            #closest intersection pt is on the exterior linestring we are looking for 
            closest_pt=sorted(intersection_pts, key=line.distance)[0]
            s,e=ext_linestr_lst[(closest_pt.x, closest_pt.y)].boundary
            
            #add the end points of the interior and exterior linestrings to the corner points list 
            pts_on_corners.append("m")
            pts_on_corners.append(segments_pts[-1])
            pts_on_corners.append((e.x, e.y))
                
    #first point that was added is part of the last corner foorprint so we move it and put it at the end of the points list 
    for i in range(0,len(pts_on_corners)): 
        if pts_on_corners[i]=="e": 
            pts_on_corners=pts_on_corners[i:]+pts_on_corners[:i]
            break 
    #to build the footprint on the corners we used flags s-e-m
    #normal corner footprints are made by the coordinates e123 - s1
    #footprints merged with adjacent footprints are made by the coordinates  e123 - m12 - m12 - ... - s1
    #a footprint always end with s and always starts with e 
    #iniatialize the first footprint 
    footprint_corner_pts=None 
    
    for i in range(0,len(pts_on_corners)): 
        #alsmost always starts with e
        if pts_on_corners[i]=="e":  
            footprint_corner_pts=pts_on_corners[i+1:i+4]
            e1 = Point(pts_on_corners[i+1])
 
        elif pts_on_corners[i]=="s":  
            #when we reach a flag s, we have the final coordinate of the current footprint 
            #so we create the shapely polygon 
            footprint_corner_pts=footprint_corner_pts+[pts_on_corners[i+1]] 
            if len(footprint_corner_pts) > 3: 
                p2 = Point(pts_on_corners[i+2])
                p0= Point(pts_on_corners[i+1])
                p1=Point(footprint_corner_pts[0])
                split_side_1 = footprint_corner_pts[-1]
                split_side_2 = footprint_corner_pts[-2]
                e2 = None 
                if p2.distance(p1) > 0.1: 
                    footprint_corner_pts=footprint_corner_pts +[pts_on_corners[i+2]]
                    e2=p2
                
                if len(footprint_corner_pts) >=4: 
                    footprint_corner=Polygon(footprint_corner_pts)
                    splitting_side = LineString([split_side_1, split_side_2])
              
                    #footprints = footprints + case_2(footprint_corner, min_area, max_perim, facade_depth, facade_len, splitting_side)
                    corner_fp = case_2(footprint_corner, min_area, max_perim, facade_depth, facade_len_min, facade_len_max, splitting_side)
                    for fp in corner_fp: 
                        #if the corner is too big ! identify the longest side for splitting along it. Determine the revese value depeneding if the longest 
                        #is on the exterir of te input polygon or not 
                        if fp.area > max_area: 
                            longest_len= 0 
                            for coord in range(0, len(fp.exterior.coords)-1):
                                ext_linestr=LineString([fp.exterior.coords[coord], fp.exterior.coords[coord+1] ])
                                if ext_linestr.length>longest_len: 
                                    longest_len=ext_linestr.length 
                                    longest_side=ext_linestr
                                    a,b = longest_side.boundary
                            if collinear(a,b,e1) or (e2!= None and collinear(a,b,e2)): 
                            #if a.distance(e1) <0.1 or b.distance(e1) < 0.1 or b.distance(e2) <0.1 or a.distance(e2) < 0.1 : 
                                reverse = True 
                            else: 
                                reverse = False 
                            footprints = footprints + case_2(fp, min_area, max_perim, facade_depth, facade_len_min, facade_len_max, None, reverse)
                        else: 
                            footprints = footprints + [fp]

            footprint_corner_pts = []
           
        elif pts_on_corners[i]=="m": 
            #if the footprint is not None, we just add the points to the current footprint being created 
            if footprint_corner_pts!=None: 
                footprint_corner_pts=[pts_on_corners[i+2]]+footprint_corner_pts+[pts_on_corners[i+1]]
                #if we are at the last flag on the list of corner points, it means that the footprint ends here
                #so we create a shapely polygon 
                if i==len(pts_on_corners)-3:
                    if len(footprint_corner_pts) > 3:  
                        footprint_corner=Polygon(footprint_corner_pts)
                        footprints.append(footprint_corner)    
            #if the footprint is None, it means that the whole list is made of flags m (all footprints are merged) 
            else: 
                #call case 3 to try another subdivision case that would suit better than this one AND returns the footprints here.
                #footprint_corner_pts=[pts_on_corners[i+2]]+[pts_on_corners[i+1]]
                footprints = case_3(input_poly, min_area, max_perim, facade_depth, facade_len_min, facade_len_max)
                return footprints 
      
    return footprints

        
def split_into_footprints(input_poly, avg_area, min_area, max_perim, facade_len_min, facade_len_max, min_facade_depth, max_facade_depth, reverse = False, concave= True): 
    """subdivides a building plot into individual building footprints considering 5 cases 

    Keyword arguments:
    input_poly -- input building plot (shapely polygon)
    avg_area -- average building footprint area 
    min_area -- minimum area of footprint to be considered as valid footprint 
    max_perim -- max perimeter allowed for footprint to be valid 
    facade_len_min -- minimum facade lenght allowed 
    facade_len_max -- max facade length allowed 
    facade_depth -- min facade depth allowed 
    max_facade_depth -- max facade depth allowed 
    concave -- binary (T or F) telling whether input_poly is concave or not 
    reverse -- binary (T or F) indicating on which side of the building footprint we shorten it for adding randomness in the building depth 

    
    Returns:
    list of individual footprints as shapely polygons 
    """ 

    facade_depth = min_facade_depth
    try: 
        #CASE 1 : plot is already a footprint. Keep original plot OR take MBR to correct for weird building shapes. 
        footprints=[]
        if input_poly.area <= avg_area and len(input_poly.exterior.coords)-1 > 3: 
            if input_poly.is_valid and input_poly.area >= min_area : # and input_poly.length <= max_perim: 
                building_footprint=input_poly
                footprints.append(building_footprint)

        #other CASES 
        else: 

            #CASE 5 polygon is concave 
            if concave and ((input_poly.convex_hull.area - input_poly.area) <=-0.02 or (input_poly.convex_hull.area - input_poly.area) >=0.02): 
                footprints=footprints+case_5(input_poly, avg_area, min_area, max_perim, facade_len_min, facade_len_max, facade_depth, max_facade_depth)
            else: 
                #take input polygon MBR 
                mbr=input_poly.minimum_rotated_rectangle
                #extract length of shortest MBR side 
                shortest_len=math.inf 
                for coord in range(0, 2):
                    ext_linestr=LineString([mbr.exterior.coords[coord],mbr.exterior.coords[coord+1] ])
                    if ext_linestr.length<shortest_len: 
                        shortest_len=ext_linestr.length
                #CASE 2: there is only one row of houses in the plot (no courtyard)
                if shortest_len // facade_depth <= 1:
                    footprints=footprints+case_2(input_poly, min_area, max_perim, facade_depth, facade_len_min, facade_len_max, None, reverse)
                
                #CASE 3: there is two rows of houses in the plot (no courtyard)
                elif shortest_len // facade_depth == 2: 
                    footprints=footprints+case_3(input_poly, min_area, max_perim, facade_depth, facade_len_min, facade_len_max)

                #CASE 4: the footprint is big enough to make an interior courtyard 
                else: 
                    footprints=footprints+case_4(input_poly, min_area, max_perim, facade_len_min, facade_len_max, max_facade_depth, min_facade_depth)
        return footprints  

    except: 
        return [input_poly]

