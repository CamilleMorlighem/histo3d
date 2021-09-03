import numpy as np
import matplotlib.pyplot as plt
from numpy.lib.function_base import diff
import shapely
import math
import shapely.geometry
import pandas as pd
import geopandas as gpd


"""
This script is used to evaluate the generalisation and vectorisation of the building plots computing two metrics: difference in turning functions and shape similarity 
It is used by the script building_plots_extraction.py 

Credit: the turning function implementation is based on the matlab implementation of Matan WEISSBUCH (https://github.com/matan2050/TFMatcher)
"""


def tf_value(turningFunction, xValue): 
    """Compute the value of the theta function

    Keyword arguments:
    turningFunction -- turning function 
    xValue -- value in the x axis in which we want the value of the turning function (xValue between 0 and 1)

    Returns:
    The value of the theta function
    """

    for i in range(0, len(turningFunction), 1): 
        if xValue < turningFunction[i+1, 0] and xValue>=turningFunction[i,0]: 
            fValue = turningFunction[i,1]
            break 

        if xValue == turningFunction[i+1, 0] or abs(turningFunction[i+1, 0] - xValue) < 0.000001: 
            fValue= turningFunction[i+1, 1]
            break 
    return fValue 


def create_turning_function_degrees(tp): 
    """Find the turning function of a given polygon in DEGREES 

    Keyword arguments:
    tp -- input polygon 
    
    Returns:
    Matrix with x and y values forming the turning function 
    """

    tf_function_y=[]
    tf_function_x=[]
    tf_inv_function_x=[]
    tf_inv_function_y=[]
    n_turning_pts=len(tp.exterior.coords)-1
    lineVec = np.empty((n_turning_pts,2))
    lineVec[:] = np.NaN

    invLineVec = np.empty((n_turning_pts,2))
    invLineVec[:] = np.NaN
    perimeter = 0

    for i in range(0,n_turning_pts): 
        #calculating current element in turning function
        lineVec[i,0] = tp.exterior.coords[i+1][0]-tp.exterior.coords[i][0]
        lineVec[i,1] = tp.exterior.coords[i+1][1]-tp.exterior.coords[i][1]

        #summing perimeter of the polygon
        perimeter= perimeter+np.linalg.norm(lineVec[i,:])
        
        #calculating current element in reverse turning function
        invLineVec[i,0] = tp.exterior.coords[n_turning_pts-i-1][0] - tp.exterior.coords[n_turning_pts+1-i-1][0]
        invLineVec[i,1]= tp.exterior.coords[n_turning_pts-i-1][1] - tp.exterior.coords[n_turning_pts+1-i-1][1]

    unitLineVec = np.empty((n_turning_pts,2))
    unitLineVec[:] = np.NaN
    unitInvLineVec = np.empty((n_turning_pts,2))
    unitInvLineVec[:] = np.NaN

    for i in range(0,len(lineVec)): 
        unitLineVec[i,:] = lineVec[i,:] / perimeter
        unitInvLineVec[i,:] = invLineVec[i,:] / perimeter

    #using dot products between unit vectors to create the turning functions y
    #values
    angles = np.empty((n_turning_pts, 1))
    angles[:] = np.NaN
    invAngles = np.empty((n_turning_pts, 1))
    invAngles[:] = np.NaN

    #calculating values of turning angles
    angles[0]= math.degrees(np.arccos(np.dot(unitLineVec[0,:], np.array([1, 0]))))
    invAngles[0] = math.degrees(np.arccos(np.dot(unitInvLineVec[0,:], np.array([1, 0]))))

    for i in range (1, n_turning_pts-1): 
        angles[i] = math.degrees(np.arccos(np.dot(unitLineVec[i,:], unitLineVec[i+1,:])))
        invAngles[i] =math.degrees(np.arccos(np.dot(unitInvLineVec[i,:], unitInvLineVec[i+1,:])))

    #finally, adding the last angle between the last and first lines
    angles[-1] =math.degrees(np.arccos(np.dot(unitLineVec[n_turning_pts-1, :], unitLineVec[0, :])))
    invAngles[-1] =math.degrees(np.arccos(np.dot(unitInvLineVec[n_turning_pts-1, :], unitInvLineVec[0,:])))

    #summing angles to create the turning function's y values
    nAngles = len(angles)
    tf_function_y=  np.empty((nAngles+1, 1))
    tf_function_y[:] = np.NaN
    tf_inv_function_y= np.empty((nAngles+1,1))
    tf_inv_function_y[:] = np.NaN
    summ = 0
    invSum = 0

    for i in range(0,nAngles): 
        summ = summ+ angles[i]
        invSum = invSum + invAngles[i]
        tf_function_y[i] = summ 
        tf_inv_function_y[i] = invSum

    tf_function_y[-1] = summ 
    tf_inv_function_y[-1] = invSum

    #summing perimeter lengths for function's x values
    tf_function_x=  np.empty((nAngles+1, 1))
    tf_function_x[:] = np.NaN
    tf_inv_function_x=  np.empty((nAngles+1, 1))
    tf_inv_function_x[:] = np.NaN
    summ = 0
    invSum = 0

    tf_inv_function_x[0] = 0
    tf_function_x[0] = 0

    for i in range(1, nAngles): 
        summ = summ + np.linalg.norm(unitLineVec[i,:]) 
        invSum = invSum + np.linalg.norm(unitInvLineVec[i,:]) 
        tf_function_x[i] = summ 
        tf_inv_function_x[i] = invSum 

    #closing the function on total length 1
    tf_function_x[-1]= 1
    tf_inv_function_x[-1]= 1

    #for interfacing purposes, matrix that hold xs and ys
    tfMat=np.asmatrix(np.column_stack((tf_function_x,tf_function_y)))
    itfMat  = np.asmatrix(np.column_stack((tf_inv_function_x,tf_function_y)))

    return tfMat


def create_turning_function(tp): 
    """Find the turning function of a given polygon in RADIANS (default)

    Keyword arguments:
    tp -- input polygon 
    
    Returns:
    Matrix with x and y values forming the turning function 
    """

    tf_function_y=[]
    tf_function_x=[]
    tf_inv_function_x=[]
    tf_inv_function_y=[]
    n_turning_pts=len(tp.exterior.coords)-1
    lineVec = np.empty((n_turning_pts,2))
    lineVec[:] = np.NaN

    invLineVec = np.empty((n_turning_pts,2))
    invLineVec[:] = np.NaN
    perimeter = 0

    for i in range(0,n_turning_pts): 
        #calculating current element in turning function
        lineVec[i,0] = tp.exterior.coords[i+1][0]-tp.exterior.coords[i][0]
        lineVec[i,1] = tp.exterior.coords[i+1][1]-tp.exterior.coords[i][1]

        #summing perimeter of the polygon
        perimeter= perimeter+np.linalg.norm(lineVec[i,:])
        
        #calculating current element in reverse turning function
        invLineVec[i,0] = tp.exterior.coords[n_turning_pts-i-1][0] - tp.exterior.coords[n_turning_pts+1-i-1][0]
        invLineVec[i,1]= tp.exterior.coords[n_turning_pts-i-1][1] - tp.exterior.coords[n_turning_pts+1-i-1][1]

    unitLineVec = np.empty((n_turning_pts,2))
    unitLineVec[:] = np.NaN
    unitInvLineVec = np.empty((n_turning_pts,2))
    unitInvLineVec[:] = np.NaN

    for i in range(0,len(lineVec)): 
        unitLineVec[i,:] = lineVec[i,:] / perimeter
        unitInvLineVec[i,:] = invLineVec[i,:] / perimeter

    #using dot products between unit vectors to create the turning functions y
    #values
    angles = np.empty((n_turning_pts, 1))
    angles[:] = np.NaN
    invAngles = np.empty((n_turning_pts, 1))
    invAngles[:] = np.NaN

    #calculating values of turning angles
    angles[0]= (np.arccos(np.dot(unitLineVec[0,:], np.array([1, 0]))))
    invAngles[0] = (np.arccos(np.dot(unitInvLineVec[0,:], np.array([1, 0]))))

    for i in range (1, n_turning_pts-1): 
        angles[i] = (np.arccos(np.dot(unitLineVec[i,:], unitLineVec[i+1,:])))
        invAngles[i] = (np.arccos(np.dot(unitInvLineVec[i,:], unitInvLineVec[i+1,:])))

    #finally, adding the last angle between the last and first lines
    angles[-1] = (np.arccos(np.dot(unitLineVec[n_turning_pts-1, :], unitLineVec[0, :])))
    invAngles[-1] = (np.arccos(np.dot(unitInvLineVec[n_turning_pts-1, :], unitInvLineVec[0,:])))

    #summing angles to create the turning function's y values
    nAngles = len(angles)
    tf_function_y=  np.empty((nAngles+1, 1))
    tf_function_y[:] = np.NaN
    tf_inv_function_y= np.empty((nAngles+1,1))
    tf_inv_function_y[:] = np.NaN
    summ = 0
    invSum = 0

    for i in range(0,nAngles): 
        summ = summ+ angles[i]
        invSum = invSum + invAngles[i]
        tf_function_y[i] = summ 
        tf_inv_function_y[i] = invSum

    tf_function_y[-1] = summ 
    tf_inv_function_y[-1] = invSum

    #summing perimeter lengths for function's x values
    tf_function_x=  np.empty((nAngles+1, 1))
    tf_function_x[:] = np.NaN
    tf_inv_function_x=  np.empty((nAngles+1, 1))
    tf_inv_function_x[:] = np.NaN
    summ = 0
    invSum = 0

    tf_inv_function_x[0] = 0
    tf_function_x[0] = 0

    for i in range(1, nAngles): 
        summ = summ + np.linalg.norm(unitLineVec[i,:]) 
        invSum = invSum + np.linalg.norm(unitInvLineVec[i,:]) 
        tf_function_x[i] = summ 
        tf_inv_function_x[i] = invSum 

    #closing the function on total length 1
    tf_function_x[-1]= 1
    tf_inv_function_x[-1]= 1

    #for interfacing purposes, matrix that hold xs and ys
    tfMat=np.asmatrix(np.column_stack((tf_function_x,tf_function_y)))
    itfMat  = np.asmatrix(np.column_stack((tf_inv_function_x,tf_function_y)))

    return tfMat


def compute_difference_tf(row): 
    """Compute difference between turning functions 

    Keyword arguments:
    row -- row of the dataframe with 2 polygons between which we want to compute the diff in turning function
    
    Returns:
    Difference of turning functions -- difference in areas under the curves made by the functions in cartesian coordinate system 
    """

    poly1=row["geometry_x"]
    poly2=row["geometry_y"]
    area_lst= []
    pts = list(poly1.exterior.coords)[:-1]
    for i in range(0, len(pts)):
        poly1 = shapely.geometry.Polygon(list(map(list,pts[i:] + pts[0:i] + [pts[i]])))
        if poly1!= None and poly2!=None: 
            fction1 = create_turning_function(poly1)
            fction2= create_turning_function(poly2)
    
            from shapely.geometry import Polygon, LineString
            from shapely.ops import unary_union, polygonize
            x_y_curve1 = fction1.tolist()
            x_y_curve2= fction2.tolist()
            lst1= []
            for i in x_y_curve1: 
                lst1.append(tuple(i))

            lst2= []
            for i in x_y_curve2: 
                lst2.append(tuple(i))

            x_y_curve1 = lst1
            x_y_curve2=lst2
            polygon_points = [] #creates a empty list where we will append the points to create the polygon

            for xyvalue in x_y_curve1:
                polygon_points.append([xyvalue[0],xyvalue[1]]) #append all xy points for curve 1

            for xyvalue in x_y_curve2[::-1]:
                polygon_points.append([xyvalue[0],xyvalue[1]]) #append all xy points for curve 2 in the reverse order (from last point to first point)

            for xyvalue in x_y_curve1[0:1]:
                polygon_points.append([xyvalue[0],xyvalue[1]]) #append the first point in curve 1 again, to it "closes" the polygon

            polygon = Polygon(polygon_points)
            area = polygon.area

            x,y = polygon.exterior.xy
                # original data
            ls = LineString(np.c_[x, y])
                # closed, non-simple
            lr = LineString(ls.coords[:] + ls.coords[0:1])
            lr.is_simple  # False
            mls = unary_union(lr)
            mls.geom_type  # MultiLineString'

            Area_cal =[]

            for polygon in polygonize(mls):
                Area_cal.append(polygon.area)
                Area_poly = (np.asarray(Area_cal).sum())
            area_lst.append(Area_poly)
        
        return min(area_lst)


"""
###DEPRECATED FUNCTION###
def compute_difference(poly1, poly2): 
    if poly1!= None and poly2!=None: 
        tfL = create_turning_function(poly1)
        tfR = create_turning_function(poly2)
        m=len(tfL)
        n=len(tfR)

        #building the xEventOrder vector that holds the x values of each point
        #where there is a change in the turning function y value
        xEventOrder = tfL[:,0]
        xEventOrder=np.append(xEventOrder, [[tfL[len(tfL)-1,1-1]]], axis=0)
        xEventOrder = np.append(xEventOrder, tfR[:,0], axis=0)
        xEventOrder=np.append(xEventOrder, [[tfR[len(tfR)-1,1-1]]], axis=0)
        xEventOrder= np.sort(xEventOrder, axis=0)
        xEventOrder=np.unique(xEventOrder, axis=0)

        difference = 0

        for i in range(1,len(xEventOrder), 1):  
            stripLength = xEventOrder[i] - xEventOrder[i-1] 
            difference = difference + ((tf_value(tfL,xEventOrder[i]) - tf_value(tfR,xEventOrder[i]))**2)*stripLength
    else: 
        difference=9999
    return float(difference)
"""

def plot_turning_functions(tfL, tfR):
    plt.subplot(1,2, 1)
    plt.step(tfL[:,0],tfL[:,1], where='post')
    plt.subplot(1,2, 2)
    plt.step(tfR[:,0],tfR[:,1], where='post')
    plt.show()

    
def compute_shape_similarity(row): 
    """Compute shape similarity between two polygons, based on the buffer-based measure of Samal2004

    Keyword arguments:
    row -- row of the dataframe with 2 polygons between which we want to compute the diff in turning function
    
    Returns:
    Shape similarity 
    """

    poly1=row["geometry_x"]
    poly2=row["geometry_y"]
    if poly1 !=None and poly2!=None: 
        buffered1=poly1.buffer(1, cap_style=3, join_style=2)
        buffered2=poly2.buffer(1, cap_style=3, join_style=2)
       
        p1_in_b2=100*(poly1.intersection(buffered2).area)/poly1.area
        p2_in_b1= 100*(poly2.intersection(buffered1).area)/poly2.area
        shape_similarity=(p1_in_b2+p2_in_b1)/2
    else: 
        shape_similarity=9999
    return shape_similarity


def eval_generalisation(shapefile_list1, shapefile_list2):
    """Compare a list of polygon in a shapefile with another list in another shapefile by computing two metrics (shape similarity and difference between turning functions)

    Keyword arguments:
    shapefile_list1 -- shapefile 1 
    shapefile_list2 -- shapefile 2 
    
    Returns:
    DF with four columns: polygon1, polygon2, their shape similarity and their diff in turning fction 
    """
    entity_set1 = gpd.read_file(shapefile_list1)
    entity_set2 = gpd.read_file(shapefile_list2)
  
    joined=entity_set1.merge(entity_set2, how='left', on="FeaID")
    joined["shape_similarity"] = joined.apply(compute_shape_similarity, axis=1)
    joined["diff_turning_f"] = joined.apply(compute_difference_tf, axis=1)
    joined=joined[joined["shape_similarity"]!=9999]
    return joined 


def graph_thesis(poly2): 
    fction2= create_turning_function_degrees(poly2)
    plt.step(fction2[:,0],fction2[:,1], where='post', color = "darkblue", label="$\it{T_{A}}$")
    plt.xlabel("Normalised accumulated length", labelpad=8, fontsize=14)
    plt.ylabel("Accumulated tangent angle (°)", labelpad=8, fontsize=14)
    #plt.title("Turning function of a rectangle", fontsize=12, fontweight="bold")
    plt.subplots_adjust(wspace=0.3, hspace=0.6)
    plt.ylim(1, 360)
    plt.legend(fontsize=14)
    plt.subplots_adjust(bottom=0.2)
    plt.show()


"""
#test
poly2=shapely.geometry.Polygon([[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]])
graph_thesis(poly2)
"""


