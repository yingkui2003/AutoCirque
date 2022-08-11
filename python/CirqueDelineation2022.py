#-------------------------------------------------------------------------------
# Name: CirqueDelineation.py
# Purpose: This tool automatically delineates cirques from DEM. The outlet
# feature will be automatically identified based on the turning points of
# the first order stream lines. However, the results may not good enough.
# The best practice is to manually draw points (can also determine the cirque
# outlet points using the thid tool in this toolbox) or cross sections for
# potential cirques and then run the cirque delineation use the second tool.
# 
# Author: Dr. Yingkui Li
# Created:     09/21-12/29/2020
# Updated:     05/10-05/15/2021
# revised:     09/01/2021 (based on review comments for the paper in Geomorphology)
# Department of Geography, University of Tennessee
# Knoxville, TN 37996
#-------------------------------------------------------------------------------

# Import arcpy module
import arcpy, sys
from arcpy import env
from arcpy.sa import *
#import numpy
import numpy as np
import time

arcpy.env.overwriteOutput = True
arcpy.env.XYTolerance= "0.01 Meters"

arcpy.Delete_management("in_memory") ### Empty the in_memory
ArcGISPro = 0
arcpy.AddMessage("The current python version is: " + str(sys.version_info[0]))
if sys.version_info[0] == 2:  ##For ArcGIS 10, need to check the 3D and Spatial Extensions
    try:
        if arcpy.CheckExtension("Spatial")=="Available":
            arcpy.CheckOutExtension("Spatial")
        else:
            raise Exception ("not extension available")
            #print "not extension available"
    except:
        raise Exception ("unable to check out extension")
        #print "unable to check out extension"

    try:
        if arcpy.CheckExtension("3D")=="Available":
            arcpy.CheckOutExtension("3D")
        else:
            raise Exception ("not extension available")
            #print "not extension available"
    except:
        raise Exception ("unable to check out extension")
        #print "unable to check out extension"
elif sys.version_info[0] == 3:  ##For ArcGIS Pro
    ArcGISPro = 1
    #pass ##No need to Check
else:
    raise Exception("Must be using Python 2.x or 3.x")
    exit()   
#---------------------------------------------------------------------------------------------------------------
## The function to clean extrlines based on from and to nodes
## if only one to node and no corresponding from node, except for the highest facc section, marking for deletion
## The same processes are iterated to remove all extra lines
##This is more efficient than another clean line method based on the intersect of the to node points
#---------------------------------------------------------------------------------------------------------------
def cleanextralineswithtopology(inline,outline, field, faccThreshold):
    bflag = 1
    while bflag:
        bflag = 0
        lineArray = arcpy.da.FeatureClassToNumPyArray(inline,['OID@','from_node','to_node', field])
        fromnode = np.array([item[1] for item in lineArray])
        tonode = np.array([item[2] for item in lineArray])
        facc = np.array([item[3] for item in lineArray])
        uniquetonode = np.unique(tonode)
        maxfacc = max(facc)
        #maxfacc = min(faccThreshold,maxfacc)
        if faccThreshold > 0:
            maxfacc = faccThreshold
            
        lineid = [] ##Record the id for deletion
        for i in range(len(lineArray)):
            linetonode = lineArray[i][2]
            nodecount = np.count_nonzero(uniquetonode == linetonode)
            if nodecount == 1 and not (linetonode in fromnode) and lineArray[i][3] < maxfacc: ###only one tonode except for the highest facc section
                ##print "mark for deletion"
                lineid.append(lineArray[i][0])
                bflag = 1

        ##Delete the line marked for deletion
        with arcpy.da.UpdateCursor(inline, "OID@") as cursor:
            for row in cursor:
                if int(row[0]) in lineid:
                    cursor.deleteRow()     
        del cursor, row

    arcpy.CopyFeatures_management(inline, outline)
    return outline

#---------------------------------------------------------------------------------------------------------------
# This function calculates the SL ratio of a set of XY points
#---------------------------------------------------------------------------------------------------------------
def SLRatio(x, y):
    x_diff = x[1:] - x[0:-1]
    y_diff = y[1:] - y[0:-1]
    upxdiff = x_diff[0:-1]
    dnxdiff = x_diff[1:]
    upydiff = y_diff[0:-1]
    dnydiff = y_diff[1:]
    upsl = upydiff/upxdiff #* x[1:-1]
    dnsl = dnydiff/dnxdiff #* x[1:-1]
    xbothside = x[2:] - x[0:-2]
    ybothside = y[2:] - y[0:-2]
    slbothside = ybothside / xbothside #* x[1:-1]
    slbothside[slbothside == 0] = -0.001 ##try to get rid of the micro changes in flat areas
    Ratio = (dnsl - upsl) / slbothside
    ##try to get rid of the micro changes in flat areas
    Ratio [dnsl > -0.087] = 0.0 ##if the downslope is less than 5 degree, set the ratio to zero
    Ratio [slbothside > -0.087] = 0.0 ##if the both is less than 5 degree, set the ratio to zero

    ##using the NDVI format to derive the ratio
    #Ratio = (dnsl - upsl)/ (dnsl + upsl)
    
    return Ratio

#---------------------------------------------------------------------------------------------------------------
# This function calculates the RDE ratio of a set of XY points
#---------------------------------------------------------------------------------------------------------------
def RDERatio(x, y):
    deltax = x[2:] - x[0:-2]
    deltay = y[2:] - y[0:-2]
    RDEs = -1.0* deltay / deltax * x[1:-1]
    H = y[0] - y[-1]
    L = x[-1]
    RDEt = H / max(0.0001, math.log(L))
    Ratio = RDEs / RDEt
    return Ratio

#---------------------------------------------------------------------------------------------------------------
# This function determines the turning points along a line
#---------------------------------------------------------------------------------------------------------------    
def turning_points(streamLength, streamZ, turning_points = 10, cluster_radius = 200):

    new_x, new_y = np.array(streamLength), np.array(streamZ)
    SL = SLRatio(new_x, new_y)
    turn_point_idx = np.argsort(SL)[::-1]
    t_points = []
    t_ratios = []
    #determine the maximum number of turning points based on the normal distrbution of one standard dieviation (100% - 68%)/2 = 16%
    SL_positive = [i for i in SL if i >= 0]
    ##print len(SL_positive)
    turning_points = min(turning_points, int(len(SL_positive)*0.16+0.5))
    arcpy.AddMessage("turing_points is adjusted to: " + str(turning_points)) 

    while len(t_points) < turning_points and len(turn_point_idx) > 0:
        Ratio = SL[turn_point_idx[0]]
        if Ratio < 0.0:
            break
        else:
            t_points += [turn_point_idx[0]]
            t_ratios.append(SL[turn_point_idx[0]])
            cumLength = new_x[turn_point_idx[0]]
            trueidx = np.where(np.abs(new_x-cumLength) < cluster_radius)
            if len(trueidx[0])> 0:
                for i in range(len(trueidx[0])):
                    index = trueidx[0][i]
                    turn_point_idx = np.delete(turn_point_idx, np.where(turn_point_idx == index))
    
    return t_points, t_ratios

#---------------------------------------------------------------------------------------------------------------
# This function determines the turning points along a line based on RDE ratios
#---------------------------------------------------------------------------------------------------------------      
def turning_points_RDE(streamLength, streamZ, turning_points = 10, cluster_radius = 200):

    #new_x, new_y = streamLength, streamZ
    new_x, new_y = np.array(streamLength), np.array(streamZ)
    SL = RDERatio(new_x, new_y)
    turn_point_idx = np.argsort(SL)[::-1]
    t_points = []
    t_ratios = []
    #for i in range(len(SL)):
    #    if SL[i] > 2.0:
    #        t_points += [i+1]  ##use i +1 becasue the SL starts from 1 to n-1

    while len(t_points) < turning_points and len(turn_point_idx) > 0:
        SLRatio = SL[turn_point_idx[0]]
        if SLRatio < 1.5:
            break
        else:
            t_points += [turn_point_idx[0]]
            t_ratios.append(SL[turn_point_idx[0]])
            cumLength = new_x[turn_point_idx[0]]
            trueidx = np.where(np.abs(new_x-cumLength) < cluster_radius)
            if len(trueidx[0])> 0:
                for i in range(len(trueidx[0])):
                    index = trueidx[0][i]
                    turn_point_idx = np.delete(turn_point_idx, np.where(turn_point_idx == index))

    return t_points, t_ratios


#---------------------------------------------------------------------------------------------------------------
# This function determines the turning points along a line based on RDP Algorithm
#---------------------------------------------------------------------------------------------------------------      
###rdp only positive distance!!! for turning point detection

def Knickpoints_rdp(points, epsilon, turn_points, angles):
    # get the start and end points
    start = np.tile(np.expand_dims(points[0], axis=0), (points.shape[0], 1))
    end = np.tile(np.expand_dims(points[-1], axis=0), (points.shape[0], 1))
    linedist = Dist(start[0][0],start[0][1],end[0][0],end[0][1])
    ##print linedist
    dist_point_to_line = np.cross(end - start, points - start, axis=-1) / np.linalg.norm(end - start, axis=-1)
    ##print dist_point_to_line
    #max_idx = np.argmax(np.abs(dist_point_to_line))
    #max_value = dist_point_to_line[max_idx]/linedist
    max_idx = np.argmax(dist_point_to_line)
    max_value = dist_point_to_line[max_idx]##/linedist
    ##print max_value
    ##print points[max_idx]
    
    if abs(max_value) > epsilon:  ##the distance is at least 1 m from the line
        if max_value > 0:
             ##Calculate the angle of this maximum point to the both side
             #centerpnt =  points[max_idx]
             #angle = math.atan2(end[0][1] - centerpnt[1], end[0][0] - centerpnt[0]) - math.atan2(start[0][1] - centerpnt[1], start[0][0] - centerpnt[0]);
             turn_points.append(points[max_idx])
             #ba = start[0] - centerpnt
             #bc = end[0] - centerpnt
             #cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
             #angle = np.arccos(cosine_angle)*180/np.pi
             #angles.append(180 - angle)
             angles.append(max_value)
             #angles.append(180 - (360 + angle/3.14159*180))

        partial_results_left = Knickpoints_rdp(points[:max_idx+1], epsilon, turn_points,angles)
        partial_results_right = Knickpoints_rdp(points[max_idx:], epsilon, turn_points,angles)
        
def turning_points_RDP(streamLength, streamZ, turning_points = 10, cluster_radius = 200):
    stream_points = np.concatenate([np.expand_dims(streamLength, axis=-1), np.expand_dims(streamZ, axis=-1)], axis=-1)

    epsilon = 0.01
    turn_points = []
    turn_angles = []

    Knickpoints_rdp(stream_points, epsilon, turn_points, turn_angles)
    ##print "turn_points"
    ##print turn_points
    ##print turn_angles

    if len(turn_points) < 1:
        return [], []
    
    new_x = np.exp(np.array(turn_points)[:,0]) ##The value for the calculation is ln value; Need to convert the length to the reallength value, 
    ##print new_x
    turn_point_idx = np.argsort(turn_angles)[::-1]

    t_pointsID = []
    t_angles = []

    while len(t_pointsID) < turning_points and len(turn_point_idx) > 0:
        angle = turn_angles[turn_point_idx[0]]
        ##print angle
        if angle < 0.01: ##this is the angle based on the log value
            break
        else:
            t_pointsID += [turn_point_idx[0]]
            t_angles.append(turn_angles[turn_point_idx[0]])
            cumLength = new_x[turn_point_idx[0]]
            ##print cumLength
            ##print np.abs(new_x - cumLength)
            
            trueidx = np.where(np.abs(new_x - cumLength) < cluster_radius)
            if len(trueidx[0])> 0:
                for i in range(len(trueidx[0])):
                    index = trueidx[0][i]
                    turn_point_idx = np.delete(turn_point_idx, np.where(turn_point_idx == index))

    t_points = []
    
    ##Find the original  point-ID infomation
    for i in range(len(t_pointsID)):
        points = turn_points[t_pointsID[i]]
        t_point_idx = np.where(PointZ == points[1])
        ##print t_point_idx[0]
        ##print t_point_idx[0][-1] ##take the last point if it multiple points have the same elevation
        t_points.append(t_point_idx[0][0])
    ##print t_points
    ##print t_angles
    return t_points, t_angles



#---------------------------------------------------------------------------------------------------------------
# This function derives the convex angles for each points along the stream profile
#---------------------------------------------------------------------------------------------------------------    
def convex_angle(x, y):
    #Derive the the upslope grident along the stream profile and use the <5 degree to determine the potential cirque threholds
    dh = y[:-2] - y[1:-1]
    dx = x[1:-1] - x[:-2]

    upslopeAngle = np.degrees(np.arctan(dh/dx))

    #Derive the the downslope grident along the stream profile and use the >5 degree to determine the potential cirque threholds
    dh = y[1:-1] - y[2:]
    dx = x[2:] - x[1:-1]

    downslopeAngle = np.degrees(np.arctan(dh/dx))

    angleDiff = downslopeAngle - upslopeAngle

    sign1 = np.where(np.array(angleDiff) > 10 , 1, 0) ##for the big turn with angle difference of > 10
    sign2 = np.where(np.array(upslopeAngle) < 5 , 1, 0)
    sign3 = np.where(np.array(downslopeAngle) > 5 , 1, 0)
    sign4 = np.where(np.array(angleDiff) > 3 , 1, 0)

    sign = np.where((sign1 + sign2*sign3*sign4) > 0, 1, 0)
    
    sign_angle = sign*angleDiff
    
    return sign_angle

#---------------------------------------------------------------------------------------------------------------
# This function determines the turning points along a line based on RDE ratios
#---------------------------------------------------------------------------------------------------------------    
def turning_points_ConvexAngle(streamLength, streamZ, turning_points = 10, cluster_radius = 200):
    
    #new_x, new_y = streamLength, streamZ
    new_x, new_y = np.array(streamLength), np.array(streamZ)
    SL = convex_angle(new_x, new_y)

    turn_point_idx = np.argsort(SL)[::-1]

    t_points = []
    t_ratios = []

    while len(t_points) < turning_points and len(turn_point_idx) > 0:
        SLRatio = SL[turn_point_idx[0]]
        ##print SLRatio
        if SLRatio < 3:
            break
        else:
            t_points += [turn_point_idx[0]]
            t_ratios.append(SL[turn_point_idx[0]])
            cumLength = new_x[turn_point_idx[0]]
            trueidx = np.where(np.abs(new_x-cumLength) < cluster_radius)
            if len(trueidx[0])> 0:
                for i in range(len(trueidx[0])):
                    index = trueidx[0][i]
                    turn_point_idx = np.delete(turn_point_idx, np.where(turn_point_idx == index))

    return t_points, t_ratios

#---------------------------------------------------------------------------------------------------------------
# This function calculates the distance between two points
#--------------------------------------------------------------------------------------------------------------- 
def Dist(x1,y1,x2,y2):
    return math.sqrt(math.pow(math.fabs(x1-x2),2)+math.pow(math.fabs(y1-y2),2))

'''
#---------------------------------------------------------------------------------------------------------------
# This function derives the boundary based on the skyline analysis
#--------------------------------------------------------------------------------------------------------------- 
def SkylineBoundaryExtract(inDEM, instreams, wspoly, Elevationshift, outpoly):
    arcpy.env.extent = inDEM
    ##Simplify stream
    simplestreams = "in_memory\\simplestreams"
    arcpy.cartography.SimplifyLine(instreams, simplestreams,"POINT_REMOVE", 5)

    ##Stream to points
    streampoints = "in_memory\\streampoints"
    arcpy.FeatureVerticesToPoints_management(simplestreams, streampoints, "ALL")
    try:
        ##DEM adds Elevationshift
        outDEM = Plus(inDEM, Elevationshift)

        streampoints3D = "in_memory\\streampoints3D"
        arcpy.InterpolateShape_3d(outDEM, streampoints, streampoints3D)

        #Skyline
        skylineFc = "in_memory\\skylineFc"
        arcpy.Skyline_3d(streampoints3D, skylineFc, inDEM)

        ##Skyline to points
        SkylinePoints = "in_memory\\SkylinePoints"
        arcpy.FeatureVerticesToPoints_management(skylineFc, SkylinePoints, "ALL")

        #Point density
        pdensOut = PointDensity(SkylinePoints, "NONE", 30, NbrCircle(30, "MAP"))

        outCon = Con(pdensOut > 2, 1,0)

        wsline = "in_memory\\wsline"
        arcpy.PolygonToLine_management(wspoly, wsline)

        #Thin
        thinOut = Thin(outCon, "ZERO", "FILTER", "ROUND", "#")
        thinployline = "in_memory\\thinployline"
        arcpy.RasterToPolyline_conversion (thinOut, thinployline)

        ##Feature to polygon
        thinpolygon = "in_memory\\thinpolygon"
        arcpy.FeatureToPolygon_management([thinployline,wsline], thinpolygon)
        
        ##Select the polygons that intersect with the turning points
        #poly_layer = arcpy.MakeFeatureLayer_management(thinpolygon, "in_memory\\poly_layer")
        #stream_layer = arcpy.MakeFeatureLayer_management(instreams, "in_memory\\stream_layer")
        #arcpy.SelectLayerByLocation_management(poly_layer,"INTERSECT", stream_layer,"1 METERS","NEW_SELECTION","")
        sel_poly = "in_memory\\sel_poly"
        arcpy.SpatialJoin_analysis(thinpolygon, instreams, sel_poly, "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "1 Meters", "#")

        selcount = int(arcpy.GetCount_management(sel_poly)[0])
        if selcount > 0:
            arcpy.Dissolve_management(sel_poly, outpoly)
        else:
            arcpy.CopyFeatures_management(wspoly, outpoly)

    except arcpy.ExecuteError:
        arcpy.CopyFeatures_management(wspoly, outpoly)

    return outpoly

def BoundaryRefineBySlopeAndSkyline(dem, slopeRaster, stream,Elevationshift, max_slope, min_slope, max_angle):
    ##Test the process to just obtain the slope of >27 and the gentle area close to the stream
    sloperange = max_slope - min_slope + 1  ##the smallest angle is 20
    bndcleanpoly = "in_memory\\bndcleanpoly"
    sel_zero_poly = "in_memory\\sel_zero_poly"
    sel_zero_poly_over_stream = "in_memory\\sel_zero_poly_over_stream"
    sel_slope_poly = "in_memory\\sel_slope_poly"

    ##Skyline analysis to get the raster with the skyline boundary
    ##Simplify stream
    simplestreams = "in_memory\\simplestreams"
    arcpy.cartography.SimplifyLine(stream, simplestreams,"POINT_REMOVE", 5)

    ##Stream to points
    streampoints = "in_memory\\streampoints"
    arcpy.FeatureVerticesToPoints_management(simplestreams, streampoints, "ALL")

    ##DEM adds Elevationshift
    outDEM = Plus(dem, Elevationshift)

    streampoints3D = "in_memory\\streampoints3D"
    arcpy.InterpolateShape_3d(outDEM, streampoints, streampoints3D)

    #Skyline
    skylineFc = "in_memory\\skylineFc"
    arcpy.Skyline_3d(streampoints3D, skylineFc, dem)

    ##Skyline to points
    SkylinePoints = "in_memory\\SkylinePoints"
    arcpy.FeatureVerticesToPoints_management(skylineFc, SkylinePoints, "ALL")

    #Point density
    pdensOut = PointDensity(SkylinePoints, "NONE", 30, NbrCircle(30, "MAP"))

    outskylineCon = Con(pdensOut > 2, 1,0)

    #arcpy.CopyRaster_management(outskylineCon, "c:\\test\\outskylineCon.tif")
    
    ##Run the slope analysis
    for i in range (sloperange):
        angle = max_slope - i
        #arcpy.AddMessage( "The angle is: " + str(angle))
        ##print angle
        slp_gt = Con(slopeRaster > angle, 1, 0) ##make a big difference is not use 0
        
        ##plus the outskylineConraster
        slp_gt_skyline = BooleanOr(slp_gt, outskylineCon)

        #arcpy.CopyRaster_management(slp_gt_skyline, "c:\\test\\slp_gt_skyline.tif")
        
        ##Boundary clean?? check if it is necessaryy???
        OutBndCln = BoundaryClean(slp_gt_skyline, "ASCEND", "TWO_WAY")
        
        arr = arcpy.RasterToNumPyArray(OutBndCln,nodata_to_value=0)
        if arr.sum() > 0:
            arcpy.RasterToPolygon_conversion(OutBndCln, bndcleanpoly, "#","VALUE")
            polyarray = arcpy.da.FeatureClassToNumPyArray(bndcleanpoly,['SHAPE@AREA', 'gridcode'])
            polyarea = np.array([item[0] for item in polyarray])
            con_code = np.array([item[1] for item in polyarray])
            slopearea = polyarea[con_code > 0]
            if len(slopearea) > 1: ##if there are multiple polygons
                area_cutoff = max(slopearea)*0.5
                num_big_poly = (polyarea > area_cutoff).sum()
                if num_big_poly < 2:
                    flag = 1
                else:
                    flag = 0
            else:
                flag = 1
        else:
            slp_gt = Con(slopeRaster > 0, 1)
            arcpy.RasterToPolygon_conversion(slp_gt, bndcleanpoly, "#","VALUE")
            flag = 0
            
        if flag > 0: ##check if the polygon is bigger enough to make the stream within the slopped valley
            #Spatial join to select the polygon interested with the stream
            arcpy.Select_analysis(bndcleanpoly, sel_zero_poly, "gridcode < 1")
            arcpy.SpatialJoin_analysis(sel_zero_poly, stream, sel_zero_poly_over_stream, "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "1 Meters", "#")
            ##get the highest elevation of the dem
            polycountResult = arcpy.GetCount_management(sel_zero_poly_over_stream)
            polycount = int(polycountResult.getOutput(0))

            if polycount > 0:
                zeroPolyDEM = ExtractByMask(dem, sel_zero_poly_over_stream)
                #maxEleZeroPoly = arcpy.GetRasterProperties_management(zeroPolyDEM,"MAXIMUM")
                #maxEleZeroPoly_float = float(maxEleZeroPoly.getOutput(0))
                maxEleZeroPoly_float = zeroPolyDEM.maximum
            else:
                maxEleZeroPoly_float = 0
            arcpy.Select_analysis(bndcleanpoly, sel_slope_poly, "gridcode > 0")
            slopePolyDEM = ExtractByMask(dem, sel_slope_poly)
            #maxEleslopePoly = arcpy.GetRasterProperties_management(slopePolyDEM,"MAXIMUM")
            #maxEleslopePoly_float = float(maxEleslopePoly.getOutput(0))
            maxEleslopePoly_float = slopePolyDEM.maximum

            if maxEleslopePoly_float > maxEleZeroPoly_float:
                break

    try:
        arcpy.Append_management(sel_zero_poly_over_stream, sel_slope_poly, "NO_TEST")
        merge_poly = "in_memory\\merge_poly"
        arcpy.Dissolve_management(sel_slope_poly, merge_poly,"", "","SINGLE_PART", "DISSOLVE_LINES")

        polyline = "in_memory\\polyline"
        arcpy.FeatureToLine_management(merge_poly, polyline)

        ##choose only the longest line
        lineArray = arcpy.da.FeatureClassToNumPyArray(polyline,['SHAPE@LENGTH'])
        lineLength = np.array([item[0] for item in lineArray])
        count = len(lineLength)
        if count > 1:
            maxLength = max(lineLength)
            with arcpy.da.UpdateCursor(polyline, ['SHAPE@LENGTH']) as cursor:
                for row in cursor:
                    length = float(row[0])
                    if length < maxLength: ##This step reomve the small spurious ploygons as well
                        #arcpy.AddMessage("delete spurious lines..." )
                        cursor.deleteRow()
            del cursor, row

        newline = remove_bigturn(polyline, max_angle)
        newpoly = "in_memory\\newpoly"
        arcpy.FeatureToPolygon_management(newline, newpoly)
        
        return newpoly, True

    except:
        return bndcleanpoly, False

'''

#------------------------------------------------------------------------------------------------------------
# This function smooths the line by removing the big turns.
#------------------------------------------------------------------------------------------------------------
def remove_bigturn(line, max_angle):
    
    arcpy.FeatureVerticesToPoints_management(line, "in_memory\\line_points", 'All')

    ###Create the new line after removing the outlier points
    spatialref=arcpy.Describe(line).spatialReference
    field = 'ORIG_FID' ##The first two fields are FID and Geometry

    new_line = arcpy.CreateFeatureclass_management("in_memory", "new_line","POLYLINE", line,"","", spatialref)
    arcpy.AddField_management(new_line, field, "LONG")

    pointarray = arcpy.da.FeatureClassToNumPyArray("in_memory\\line_points", ('SHAPE@X', 'SHAPE@Y',field))
    line_ids = np.array([item[2] for item in pointarray])
    unique_line_ids = np.unique(line_ids)
    
    for fid in unique_line_ids:
        arr = pointarray[line_ids == fid]
        #arcpy.AddMessage(str(len(arr)))
        pntx = []
        pnty = []
        for i in range(len(arr)):
            pntx.append(arr[i][0])
            pnty.append(arr[i][1])
        #b = zip(*arr) ##This only works in ArcGIS not ArcGIS Pro

        #arcpy.AddMessage("pnty")
        #arcpy.AddMessage(pnty)
        
        #b3 = b2.T
        #arcpy.AddMessage("b2")
        #arcpy.AddMessage(b3)
        
        points = np.array(list(zip(pntx,pnty)))
        
        #points = np.array(list(zip(b3[0],b3[1])))
        #arcpy.AddMessage("points")
        #arcpy.AddMessage(points)

        turn_count = 1
        while turn_count > 0:
            turn_count = 0
            ##derive the updated angles
            end = points[2:]
            start = points[0:-2]
            center = points[1:-1]
            ba = start-center
            bc = end-center
            ca = end - start
            dot = []
            for i in range(len(bc)):
                dot.append(np.dot(ba[i], bc[i]))

            cosine_angle = np.array(dot) / (np.linalg.norm(ba, axis=-1) * np.linalg.norm(bc, axis=-1))
            #angles = np.degrees(np.arccos(cosine_angle))
            angles = np.degrees(np.arccos(np.maximum(np.minimum(1,cosine_angle), -1)))
            #Derive the distance to the line, so that to determine the convex slope
            dist = np.cross(ca, -ba, axis=-1) / np.linalg.norm(ca, axis=-1)
            
            for row in range(len(points)):
                if row <(len(points)-1) and row > 0:#if it is not first or last point of all
                    pntangle = angles[row-1] ##get the angle
                    pntdist = dist[row-1]    ##get the direction of the angle
                    if pntangle < max_angle and pntdist < 0: ##adjust to the mid point
                        turn_count += 1
                        midpoint = (points[row-1] + points[row+1])/2
                        points[row][0] = midpoint[0]
                        points[row][1] = midpoint[1]
                        arr[row][0] = midpoint[0]
                        arr[row][1] = midpoint[1]
            ##print turn_count
        ##Make the new feature class
        numpy_array_to_features(new_line, arr, ['SHAPE@X', 'SHAPE@Y'], field)

    ##Delete the created in-memory dataset
    arcpy.Delete_management ("in_memory\\line_points")

    return new_line
'''
def BoundaryRefineBySlopeBak(dem, slopeRaster, stream, max_slope, min_slope, max_angle):
    ##Test the process to just obtain the slope of >27 and the gentle area close to the stream
    sloperange = max_slope - min_slope  ##the smallest angle is 20
    bndcleanpoly = "in_memory\\bndcleanpoly"
    sel_zero_poly = "in_memory\\sel_zero_poly"
    sel_zero_poly_over_stream = "in_memory\\sel_zero_poly_over_stream"
    sel_slope_poly = "in_memory\\sel_slope_poly"

    for i in range (sloperange):
        angle = max_slope - i
        arcpy.AddMessage( "Processing headwall angle: " + str(angle))
        #print "Processing headwall angle: " + str(angle)
        ##print angle
        slp_gt = Con(slopeRaster > angle, 1, 0) ##make a big difference is not use 0
        OutBndCln = BoundaryClean(slp_gt)##, "ASCEND", "TWO_WAY")
        arr = arcpy.RasterToNumPyArray(OutBndCln,nodata_to_value=0)
        if arr.sum() > 0:
            arcpy.RasterToPolygon_conversion(OutBndCln, bndcleanpoly, "#","VALUE")
            polyarray = arcpy.da.FeatureClassToNumPyArray(bndcleanpoly,['SHAPE@AREA', 'gridcode'])
            polyarea = np.array([item[0] for item in polyarray])
            con_code = np.array([item[1] for item in polyarray])
            slopearea = polyarea[con_code > 0]
            if len(slopearea) > 1: ##if there are multiple polygons
                area_cutoff = max(slopearea)*0.5
                num_big_poly = (polyarea > area_cutoff).sum()
                if num_big_poly < 2:
                    flag = 1
                    #arcpy.AddMessage( "The flag 01 is: " + str(flag))
                else:
                    flag = 0
                    #arcpy.AddMessage( "The flag 02 is: " + str(flag))
            else:
                flag = 1
                #arcpy.AddMessage( "The flag 03 is: " + str(flag))
        else:
            slp_gt = Con(slopeRaster > 0, 1)
            #arr = arcpy.RasterToNumPyArray(slp_gt,nodata_to_value=0)
            #arcpy.AddMessage("the slope raster length is: " + str(len(arr)))            
            arcpy.RasterToPolygon_conversion(slp_gt, bndcleanpoly, "#","VALUE")
            flag = 0
            #arcpy.AddMessage( "The flag 04 is: " + str(flag))

        #arcpy.AddMessage( "The flag is: " + str(flag))
        
        if flag > 0: ##check if the polygon is bigger enough to make the stream within the slopped valley
            #Spatial join to select the polygon interested with the stream
            arcpy.Select_analysis(bndcleanpoly, sel_zero_poly, "gridcode < 1")
            arcpy.SpatialJoin_analysis(sel_zero_poly, stream, sel_zero_poly_over_stream, "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "1 Meters", "#")
            ##get the highest elevation of the dem
            polycountResult = arcpy.GetCount_management(sel_zero_poly_over_stream)
            polycount = int(polycountResult.getOutput(0))

            if polycount > 0:
                zeroPolyDEM = ExtractByMask(dem, sel_zero_poly_over_stream)
                #maxEleZeroPoly = arcpy.GetRasterProperties_management(zeroPolyDEM,"MAXIMUM")
                #maxEleZeroPoly_float = float(maxEleZeroPoly.getOutput(0))
                maxEleZeroPoly_float = zeroPolyDEM.maximum
            else:
                maxEleZeroPoly_float = 0

            arcpy.Select_analysis(bndcleanpoly, sel_slope_poly, "gridcode > 0")

            #polycountResult = arcpy.GetCount_management(sel_slope_poly)
            #polycount = int(polycountResult.getOutput(0))
            #arcpy.AddMessage( "The polycount is: " + str(polycount))
            #if polycount > 0:
            slopePolyDEM = ExtractByMask(dem, sel_slope_poly)
            #maxEleslopePoly = arcpy.GetRasterProperties_management(slopePolyDEM,"MAXIMUM")
            #maxEleslopePoly_float = float(maxEleslopePoly.getOutput(0))
            maxEleslopePoly_float = slopePolyDEM.maximum
            #else:
            #    maxEleZeroPoly_float = maxEleZeroPoly_float + 1 ##make sure that the value is higher

            if maxEleslopePoly_float > maxEleZeroPoly_float:
                break
            #else:
            #    #print "reduce more angles to make the stream within the slope polygon"
        #else: ##if flag == 0
        #    arcpy.CopyFeatures_management(bndcleanpoly, sel_slope_poly)
    #try:
    #polycountResult = arcpy.GetCount_management(sel_slope_poly)
    #polycount = int(polycountResult.getOutput(0))
    #arcpy.AddMessage( "The polycount is: " + str(polycount))
    merge_poly = "in_memory\\merge_poly"
    if arcpy.Exists(sel_slope_poly):
        arcpy.Append_management(sel_zero_poly_over_stream, sel_slope_poly, "NO_TEST")
        arcpy.Dissolve_management(sel_slope_poly, merge_poly,"", "","SINGLE_PART", "DISSOLVE_LINES")
    else:
        arcpy.Dissolve_management(bndcleanpoly, merge_poly,"", "","SINGLE_PART", "DISSOLVE_LINES")
        
    ##Need to delete spurious polgon 
    polyArray = arcpy.da.FeatureClassToNumPyArray(merge_poly,['SHAPE@AREA'])
    PolyNum = np.array([item[0] for item in polyArray])
    if len(PolyNum) > 1:
        maxArea = max(PolyNum)
        ##print maxLength
        with arcpy.da.UpdateCursor(merge_poly, ['SHAPE@AREA']) as cursor:
            for row in cursor:
                area = float(row[0])
                if area < maxArea: ##This step reomve the small spurious ploygons as well
                    arcpy.AddMessage("delete spurious polygons..." )
                    #print "delete spurious polygons..."
                    cursor.deleteRow()
        del cursor, row

    #arcpy.CopyFeatures_management(merge_poly, "c:\\test\\merge_poly_lyk2.shp")

    ##Try to get rid of the islands within the polygon
    polyline = "in_memory\\polyline"
    arcpy.FeatureToLine_management(merge_poly, polyline)
    #arcpy.CopyFeatures_management(polyline, "c:\\test\\merge_polyline.shp")

    ##choose only the longest line
    lineArray = arcpy.da.FeatureClassToNumPyArray(polyline,['SHAPE@LENGTH'])
    lineLength = np.array([item[0] for item in lineArray])
    count = len(lineLength)
    #arcpy.AddMessage( "The line count is: " + str(count))
    if count > 1:
        maxLength = max(lineLength)
        ##print maxLength
        with arcpy.da.UpdateCursor(polyline, ['SHAPE@LENGTH']) as cursor:
            for row in cursor:
                length = float(row[0])
                if length < maxLength: ##This step reomve the small spurious ploygons as well
                    arcpy.AddMessage("delete spurious lines..." )
                    #print "delete spurious lines..."
                    cursor.deleteRow()
        del cursor, row

    newpoly = "in_memory\\newpoly"
    #if (ArcGISPro == 1):
    #    arcpy.FeatureToPolygon_management(polyline, newpoly)    
    #else:
    newline = remove_bigturn(polyline, max_angle)
    arcpy.FeatureToPolygon_management(newline, newpoly)

    #arcpy.FeatureToPolygon_management(newline, newpoly)
    #arcpy.FeatureToPolygon_management(polyline, newpoly)
    
    return newpoly, True

    #except:
    #    return bndcleanpoly, False

def BoundaryRefineBySlopebak2(dem, slopeRaster, stream, max_slope, min_slope, max_angle):
    ##Test the process to just obtain the slope of >27 and the gentle area close to the stream
    sloperange = max_slope - min_slope  ##the smallest angle is 20
    bndcleanpoly = "in_memory\\bndcleanpoly"
    sel_zero_poly = "in_memory\\sel_zero_poly"
    sel_zero_poly_over_stream = "in_memory\\sel_zero_poly_over_stream"
    sel_slope_poly = "in_memory\\sel_slope_poly"

    for i in range (sloperange):
        angle = max_slope - i
        arcpy.AddMessage( "Processing headwall angle: " + str(angle))
        #print angle
        slp_gt = Con(slopeRaster > angle, 1, 0) ##make a big difference is not use 0
        #slp_gt.save("c:\\testdata\\slpgt")
        OutBndCln = BoundaryClean(slp_gt, "NO_SORT", "ONE_WAY")
        #OutBndCln.save("c:\\testdata\\outbndcln")
        arr = arcpy.RasterToNumPyArray(OutBndCln,nodata_to_value=0)
        if arr.sum() > 0:
            arcpy.RasterToPolygon_conversion(OutBndCln, bndcleanpoly, "#","VALUE")
            polyarray = arcpy.da.FeatureClassToNumPyArray(bndcleanpoly,['SHAPE@AREA', 'gridcode'])
            polyarea = np.array([item[0] for item in polyarray])
            con_code = np.array([item[1] for item in polyarray])
            slopearea = polyarea[con_code > 0]
            if len(slopearea) > 1: ##if there are multiple polygons
                area_cutoff = max(slopearea)*0.5
                num_big_poly = (polyarea > area_cutoff).sum()
                if num_big_poly < 2:
                    flag = 1
                    #arcpy.AddMessage( "The flag 01 is: " + str(flag))
                else:
                    flag = 0
                    #arcpy.AddMessage( "The flag 02 is: " + str(flag))
            else:
                flag = 1
                #arcpy.AddMessage( "The flag 03 is: " + str(flag))
        else:
            slp_gt = Con(slopeRaster > 0, 1)
            #arr = arcpy.RasterToNumPyArray(slp_gt,nodata_to_value=0)
            #arcpy.AddMessage("the slope raster length is: " + str(len(arr)))            
            arcpy.RasterToPolygon_conversion(slp_gt, bndcleanpoly, "#","VALUE")
            flag = 0
            #arcpy.AddMessage( "The flag 04 is: " + str(flag))

        #arcpy.AddMessage( "The flag is: " + str(flag))
        
        if flag > 0: ##check if the polygon is bigger enough to make the stream within the slopped valley
            #Spatial join to select the polygon interested with the stream
            arcpy.Select_analysis(bndcleanpoly, sel_zero_poly, "gridcode < 1")
            arcpy.SpatialJoin_analysis(sel_zero_poly, stream, sel_zero_poly_over_stream, "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "1 Meters", "#")
            ##get the highest elevation of the dem
            polycountResult = arcpy.GetCount_management(sel_zero_poly_over_stream)
            polycount = int(polycountResult.getOutput(0))

            if polycount > 0:
                zeroPolyDEM = ExtractByMask(dem, sel_zero_poly_over_stream)
                maxEleZeroPoly = arcpy.GetRasterProperties_management(zeroPolyDEM,"MAXIMUM")
                maxEleZeroPoly_float = float(maxEleZeroPoly.getOutput(0))
            else:
                maxEleZeroPoly_float = 0

            arcpy.Select_analysis(bndcleanpoly, sel_slope_poly, "gridcode > 0")

            #polycountResult = arcpy.GetCount_management(sel_slope_poly)
            #polycount = int(polycountResult.getOutput(0))
            #arcpy.AddMessage( "The polycount is: " + str(polycount))
            #if polycount > 0:
            slopePolyDEM = ExtractByMask(dem, sel_slope_poly)
            maxEleslopePoly = arcpy.GetRasterProperties_management(slopePolyDEM,"MAXIMUM")
            maxEleslopePoly_float = float(maxEleslopePoly.getOutput(0))
            #else:
            #    maxEleZeroPoly_float = maxEleZeroPoly_float + 1 ##make sure that the value is higher

            if maxEleslopePoly_float > maxEleZeroPoly_float:
                break
            #else:
            #    print "reduce more angles to make the stream within the slope polygon"
        #else: ##if flag == 0
        #    arcpy.CopyFeatures_management(bndcleanpoly, sel_slope_poly)
    #try:
    #polycountResult = arcpy.GetCount_management(sel_slope_poly)
    #polycount = int(polycountResult.getOutput(0))
    #arcpy.AddMessage( "The polycount is: " + str(polycount))
    merge_poly = "in_memory\\merge_poly"
    if arcpy.Exists(sel_slope_poly):
        arcpy.Append_management(sel_zero_poly_over_stream, sel_slope_poly, "NO_TEST")
        arcpy.Dissolve_management(sel_slope_poly, merge_poly,"", "","SINGLE_PART", "DISSOLVE_LINES")
    else:
        arcpy.Dissolve_management(bndcleanpoly, merge_poly,"", "","SINGLE_PART", "DISSOLVE_LINES")
        
    ##Need to delete spurious polgon 
    polyArray = arcpy.da.FeatureClassToNumPyArray(merge_poly,['SHAPE@AREA'])
    PolyNum = np.array([item[0] for item in polyArray])
    if len(PolyNum) > 1:
        maxArea = max(PolyNum)
        #print maxLength
        with arcpy.da.UpdateCursor(merge_poly, ['SHAPE@AREA']) as cursor:
            for row in cursor:
                area = float(row[0])
                if area < maxArea: ##This step reomve the small spurious ploygons as well
                    arcpy.AddMessage("delete spurious polygons..." )
                    cursor.deleteRow()
        del cursor, row

    #arcpy.CopyFeatures_management(merge_poly, "c:\\test\\merge_poly_lyk2.shp")

    ##Try to get rid of the islands within the polygon
    polyline = "in_memory\\polyline"
    arcpy.FeatureToLine_management(merge_poly, polyline)
    #arcpy.CopyFeatures_management(polyline, "c:\\testdata\\merge_polyline.shp")

    ##choose only the longest line
    lineArray = arcpy.da.FeatureClassToNumPyArray(polyline,['SHAPE@LENGTH'])
    lineLength = np.array([item[0] for item in lineArray])
    count = len(lineLength)
    #arcpy.AddMessage( "The line count is: " + str(count))
    if count > 1:
        maxLength = max(lineLength)
        #print maxLength
        with arcpy.da.UpdateCursor(polyline, ['SHAPE@LENGTH']) as cursor:
            for row in cursor:
                length = float(row[0])
                if length < maxLength: ##This step reomve the small spurious ploygons as well
                    arcpy.AddMessage("delete spurious lines..." )
                    cursor.deleteRow()
        del cursor, row

    newline = remove_bigturn(polyline, max_angle)
    #arcpy.CopyFeatures_management(newline, "c:\\testdata\\newline.shp")
    newpoly = "in_memory\\newpoly"
    arcpy.FeatureToPolygon_management(newline, newpoly)
    #arcpy.CopyFeatures_management(newpoly, "c:\\testdata\\newpoly.shp")
    return newpoly, True

    #except:
    #    return bndcleanpoly, False
'''    


#------------------------------------------------------------------------------------------------------------
# This function insert new features into an existing feature class (polygon, polyline or multipoint) based on
# a NumPy array. This function is from by online free code.
#------------------------------------------------------------------------------------------------------------
def numpy_array_to_features(in_fc, in_array, geom_fields, id_field):
    """
    Insert new features into an existing feature class (polygon,
    polyline or multipoint) based on a NumPy array.
 
    Parameters
    ----------
    in_fc : string
        An existing feature class to which new features will be added.
 
    in_array : structured NumPy array
        Array must include fields representing x and y coordinates, and
        an ID field.
 
    geom_fields: list of strings | string
        Field(s) representing x- and y-coordinates.
        If only a single numpy field is required (such as a field that
        has x,y coordinates included in a tuple) the field name can be
        passed in within a list or as a string.
 
    id_field: string
        The field that identifies how coordinates are grouped.  All
        coordinates with a common id value will be combined (in order
        of occurrence) into an output feature.
        The id_field is used in both the array and the feature class
        (i.e., the field name must exist in the feature class)
 
    """
    # Establish minimum number of x,y pairs to create proper geometry
    min_xy_dict = {'Polygon': 3, 'Polyline': 2, 'Multipoint': 1}
    min_xy_pairs = min_xy_dict[arcpy.Describe(in_fc).shapeType]
 
    if isinstance(geom_fields, list) and len(geom_fields) == 1:
        # Can't access a single field via a list later, extract the
        # only value
        geom_fields = geom_fields[0]
 
    with arcpy.da.InsertCursor(in_fc, ['SHAPE@', id_field]) as cursor:
        unique_array = np.unique(in_array[id_field])  # unique ids
 
        # Iterate through unique sets, get array that matches unique
        # value, convert coordinates to a list and insert via cursor.
        for unique_value in unique_array:
            a = in_array[in_array[id_field] == unique_value]
            if len(a) >= min_xy_pairs:  # skip if not enough x,y pairs
                cursor.insertRow([a[geom_fields].tolist(), unique_value])
            else:
                pass  # skip if not enough x,y pairs
    del cursor ##delete cursor
    
    return      


#------------------------------------------------------------------------------------------------------------
# This function calculates the plan closure for a cirque. The codes are modified from ACME code
#------------------------------------------------------------------------------------------------------------
def plan_clos(cirqueDEM):
    #tool to get the plan_closure
    meanH = int(cirqueDEM.mean)
    ##print meanH
    mid_height=(meanH)
    midHcontour = Contour(cirqueDEM, "in_memory/cont", 10000, mid_height)
    geometry = arcpy.CopyFeatures_management(midHcontour, arcpy.Geometry())
    #get total length of geometry to select only the longest contour if at that elevation there are more than one line
    lenghtlist=[]
    for x in geometry:
        #arcpy.AddMessage("length is:" + str(x.length))
        if x.length > 100: ##contour needs to > 100 m
            lenghtlist.append(x.length)
    if len(lenghtlist) > 0:
        maxlength=max(lenghtlist)
        index=lenghtlist.index(maxlength)
        goodgeometry=geometry[index]
        maximum=int(geometry[index].length)
        #find the coordinates of first, mid and last points along the mid_height contour
        point_coord_list=[]
        for m in range(0, maximum+1, int(maximum/2)):
            points = goodgeometry.positionAlongLine (m)
            centroid = points.centroid
            point_coord_list.append(centroid.X)
            point_coord_list.append(centroid.Y)
        #define the coordinates
        x_start,y_start=point_coord_list[0],point_coord_list[1]
        x_end,y_end=point_coord_list[4],point_coord_list[5]
        x_mid,y_mid=point_coord_list[2],point_coord_list[3]

        #to avoid dividing by 0 in the next section, i.e. when getting s1 and s2, x_mid and x_end and x1 cannot have the same value
        if x_end==x_start:
            x_end+=0.00001
        elif x_mid==x_end:
            x_end+=0.00001
        else:
            pass

        #end_start and mid_end lines midpoints coordinates

        m1x,m1y= (x_end+x_start)/2, (y_end+y_start)/2
        m2x,m2y= (x_mid+x_end)/2, (y_mid+y_end)/2

        # slope of the end_start and mid_end lines
        s1=(y_end-y_start)/(x_end-x_start)
        if s1 == 0:
            s1 = 0.00001
        ##print s1
        s2=(y_mid-y_end)/(x_mid-x_end)
        ##print s2
        if s2 == 0:
            s2 = 0.00001

        #inverse slope
        is1=-1*(1/s1)
        is2=-1*(1/s2)

        #equations that enable to define the point of intersection between the two lines passing by the midpoints and perpendicular to the end_start and mid_end segments
        a=np.array([[-is1,1],[-is2,1]])
        b=np.array([(-m1x*is1+m1y),(-m2x*is2+m2y)])
        try:
            centre=np.linalg.solve(a,b)
        except:
            arcpy.AddMessage("three points are probably colinear")
            #print "three points are probably colinear"
            return 0

        #measure distances between key points
        dist_centre_start = math.sqrt((math.pow((x_start-centre[0]),2))+(math.pow((y_start-centre[1]),2)))
        dist_centre_end = math.sqrt((math.pow((x_end-centre[0]),2))+(math.pow((y_end-centre[1]),2)))
        dist_centre_mid = math.sqrt((math.pow((x_mid-centre[0]),2))+(math.pow((y_mid-centre[1]),2)))
        dist_start_end = math.sqrt((math.pow((x_start-x_end),2))+(math.pow((y_start-y_end),2)))
        #define end_start and mid_centre segments as polylines
        array_centre_mid=arcpy.Array([arcpy.Point(x_mid, y_mid),arcpy.Point(centre[0], centre[1])])
        segment_centre_mid=arcpy.Polyline(array_centre_mid)
        array_start_end=arcpy.Array([arcpy.Point(x_start, y_start),arcpy.Point(x_end, y_end)])
        segment_start_end=arcpy.Polyline(array_start_end)
        #verify whether the mid_centre segment intersect end_start segment
        if segment_centre_mid.crosses(segment_start_end)==False:
            #calculate 360 degrees - the angle between centre, start and end points
            Angle = ((2*math.pi - (math.acos(((math.pow(dist_centre_end,2)+math.pow(dist_centre_start,2)-math.pow(dist_start_end,2))/(2*dist_centre_end*dist_centre_start)))))*180/math.pi)
        else:
            #calculate the angle between centre, start and end points
            Angle = (((math.acos(((math.pow(dist_centre_end,2)+math.pow(dist_centre_start,2)-math.pow(dist_start_end,2))/(2*dist_centre_end*dist_centre_start)))))*180/math.pi)
    else:
        Angle = 0

    ##delete the temp dataset
    arcpy.Delete_management ("in_memory/cont")

    return Angle

#------------------------------------------------------------------------------------------------------------
# This function refine the boundary based on the slope analysis.
#------------------------------------------------------------------------------------------------------------
def BoundaryRefineBySlope(dem, slopeRaster, outpnt, max_slope, min_slope, max_angle):
    ##Test the process to just obtain the slope of >27 and the gentle area close to the stream
    sloperange = max_slope - min_slope  ##the smallest angle is 20
    bndcleanpoly = "in_memory\\bndcleanpoly"
    #sel_zero_poly = "in_memory\\sel_zero_poly"
    sel_poly_over_outpnt = "in_memory\\sel_poly_over_outpnt"
    sel_slope_poly = "in_memory\\sel_slope_poly"
    refinedPoly = "in_memory\\refinedPoly"
    
    for i in range (sloperange):
        angle = max_slope - i
        arcpy.AddMessage( "Processing headwall angle: " + str(angle))
        #print angle
        slp_gt = Con(slopeRaster > angle, 1, 0) ##make a big difference is not use 0
        OutBndCln = BoundaryClean(slp_gt)#, "DESCEND", "ONE_WAY")
        arr = arcpy.RasterToNumPyArray(OutBndCln,nodata_to_value=0)
        if arr.sum() > 0:
            arcpy.RasterToPolygon_conversion(OutBndCln, bndcleanpoly, "NO_SIMPLIFY","VALUE")
            polyarray = arcpy.da.FeatureClassToNumPyArray(bndcleanpoly,['SHAPE@AREA', 'gridcode'])
            polyarea = np.array([item[0] for item in polyarray])
            con_code = np.array([item[1] for item in polyarray])
            slopearea = polyarea[con_code > 0]
            if len(slopearea) > 1: ##if there are multiple polygons
                area_cutoff = max(slopearea)*0.5
                num_big_poly = (polyarea > area_cutoff).sum()
                if num_big_poly < 2:
                    flag = 1
                else:
                    flag = 0
            else:
                flag = 1 ##why?? set this as 1
        else:
            slp_gt = Con(slopeRaster > 0, 1)
            arcpy.RasterToPolygon_conversion(slp_gt, bndcleanpoly, "NO_SIMPLIFY","VALUE")
            flag = 0

        
        if flag > 0: ##check if the polygon is bigger enough to make the stream within the slopped valley
            ##It seems that the best way to do it is to use the elevation of the flat zone, to check if the flat zone are higher than the slope zone
            ##Need to rewrite the program to loop and extract the elevation from each zone (centrid points maybe simple)
            ##09/28-2021
            ##using zonal statistics to get the mean elevation for each zone
            PolyID = arcpy.Describe(bndcleanpoly).OIDFieldName
            polyMeanEle = ZonalStatistics(bndcleanpoly, PolyID, dem, "MEAN")

            ##Get the meanElevation from the steep slope zone
            arcpy.Select_analysis(bndcleanpoly, sel_slope_poly, "gridcode > 0")
            #delete the small polygons
            polyArray = arcpy.da.FeatureClassToNumPyArray(sel_slope_poly,['SHAPE@AREA'])
            PolyNum = np.array([item[0] for item in polyArray])
            maxArea = max(PolyNum)
            with arcpy.da.UpdateCursor(sel_slope_poly, ['SHAPE@AREA']) as cursor:
                for row in cursor:
                    area = float(row[0])
                    if area < maxArea: ##This step reomve the small spurious ploygons as well
                        cursor.deleteRow()
            del cursor, row
            #arcpy.CopyFeatures_management(sel_slope_poly, "c:\\testdata\\sel_slope_poly.shp")

            slopePolyDEM = ExtractByMask(dem, sel_slope_poly)

            try: ##if the slopePoly is empty
                steepMeanEle = slopePolyDEM.mean + 0.5 ##Add 0.5 meter to make sure the selection 
            except:
                arcpy.AddMessage("Error! cannot delineate for this point")
                #arcpy.CopyFeatures_management(merge_poly, "c:\\testdata\\merge_poly.shp")
                return sel_slope_poly, False

            ## Find the raster of <= steepMeanEle
            refinedRst = Con(polyMeanEle < steepMeanEle, 1)
            arcpy.RasterToPolygon_conversion(refinedRst, refinedPoly, "NO_SIMPLIFY","VALUE")

            ##How to determine the angle is enough???
            ##point within polygon??
            arcpy.SpatialJoin_analysis(refinedPoly, outpnt, sel_poly_over_outpnt, "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "10 Meters", "#")
            ##get the highest elevation of the dem
            polycountResult = arcpy.GetCount_management(sel_poly_over_outpnt)
            polycount = int(polycountResult.getOutput(0))

            ##Delete tempory datasets
            try:
                arcpy.Delete_management (polyMeanEle)
            except:
                pass
            try:
                arcpy.Delete_management (slopePolyDEM)
            except:
                pass
            try:
                arcpy.Delete_management (refinedRst)
            except:
                pass
            if polycount > 0:
                break  ##out of loop             

    merge_poly = "in_memory\\merge_poly"
    if arcpy.Exists(refinedPoly):
        arcpy.Append_management(sel_poly_over_outpnt, refinedPoly, "NO_TEST")
        arcpy.Dissolve_management(refinedPoly, merge_poly,"", "","SINGLE_PART", "DISSOLVE_LINES")
    else:
        arcpy.Dissolve_management(bndcleanpoly, merge_poly,"", "","SINGLE_PART", "DISSOLVE_LINES")
        
    ##Need to delete spurious polgon 
    polyArray = arcpy.da.FeatureClassToNumPyArray(merge_poly,['SHAPE@AREA'])
    PolyNum = np.array([item[0] for item in polyArray])
    if len(PolyNum) > 1:
        maxArea = max(PolyNum)
        #print maxLength
        with arcpy.da.UpdateCursor(merge_poly, ['SHAPE@AREA']) as cursor:
            for row in cursor:
                area = float(row[0])
                if area < maxArea: ##This step reomve the small spurious ploygons as well
                    arcpy.AddMessage("delete spurious polygons..." )
                    cursor.deleteRow()
        del cursor, row

    #arcpy.CopyFeatures_management(merge_poly, "c:\\test\\merge_poly_lyk2.shp")

    ##Try to get rid of the islands within the polygon
    polyline = "in_memory\\polyline"
    try:
        arcpy.FeatureToLine_management(merge_poly, polyline)
        #arcpy.CopyFeatures_management(polyline, "c:\\test\\merge_polyline.shp")
    except:
        arcpy.AddMessage("Error! cannot delineate for this point")
        #arcpy.CopyFeatures_management(merge_poly, "c:\\testdata\\merge_poly.shp")
        return merge_poly, False
    
    ##choose only the longest line
    lineArray = arcpy.da.FeatureClassToNumPyArray(polyline,['SHAPE@LENGTH'])
    lineLength = np.array([item[0] for item in lineArray])
    count = len(lineLength)
    #arcpy.AddMessage( "The line count is: " + str(count))
    if count > 1:
        maxLength = max(lineLength)
        #print maxLength
        with arcpy.da.UpdateCursor(polyline, ['SHAPE@LENGTH']) as cursor:
            for row in cursor:
                length = float(row[0])
                if length < maxLength: ##This step reomve the small spurious ploygons as well
                    arcpy.AddMessage("delete spurious lines..." )
                    cursor.deleteRow()
        del cursor, row

    newline = remove_bigturn(polyline, max_angle)
    newpoly = "in_memory\\newpoly"
    arcpy.FeatureToPolygon_management(newline, newpoly)

    ##Delete the created in-memory dataset
    try:
        arcpy.Delete_management ("in_memory\\line_points")
        arcpy.Delete_management (bndcleanpoly)
        arcpy.Delete_management (sel_poly_over_outpnt)
        arcpy.Delete_management (sel_slope_poly)
        arcpy.Delete_management (refinedPoly)
        arcpy.Delete_management (polyline)
        arcpy.Delete_management (merge_poly)
        ##Delete in-memory Rasters
        arcpy.Delete_management (OutBndCln)
        arcpy.Delete_management (slp_gt)
    except:
        pass
    ##Return    
    return newpoly, True

#------------------------------------------------------------------------------------------------------------
# This function Create cross sections for a set points along the lines.
#------------------------------------------------------------------------------------------------------------
def cross_sections(points, pointfield, line, linefield, window, distance):
    ##Create a point buffer first
    Str_buffer_dis = str(window) + " Meters"
    pointsbuf = "in_memory\\pointsbuf"
    arcpy.Buffer_analysis(points, pointsbuf, Str_buffer_dis)
    cliplines = "in_memory\\cliplines"
    arcpy.Clip_analysis(line, pointsbuf, cliplines)
    singlecliplines = "in_memory\\singlecliplines"
    arcpy.MultipartToSinglepart_management(cliplines, singlecliplines)
    
    perpendicular_line = arcpy.CreateFeatureclass_management("in_memory", "perpsline","POLYLINE","","","",line)
    arcpy.AddField_management(perpendicular_line, pointfield, "Long")
    new_line_cursor = arcpy.da.InsertCursor(perpendicular_line, ('SHAPE@', pointfield))
    
    ##Loop for each point
    pointID = arcpy.Describe(points).OIDFieldName

    pointarray = arcpy.da.FeatureClassToNumPyArray(points,pointfield)
    pointfieldArr = np.array([item[0] for item in pointarray])
    pntcount = len(pointfieldArr)
    #arcpy.AddMessage(pointfieldArr)
    #pntcount_result = arcpy.GetCount_management(points)
    #pntcount = int(pntcount_result.getOutput(0))
    select_point = "in_memory\\select_point"
    select_line = "in_memory\\select_line"
    for i in range(pntcount):
        #spatialjoin to get the line section corresponding to the point
        query = pointID+" = "+str(i+1)
        #query = pointID+" = "+str(i)
        arcpy.Select_analysis(points, select_point, query) ###Just do a simply select analysis
        arcpy.SpatialJoin_analysis(singlecliplines, select_point, select_line, "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "30 Meters", "#")

        linearray = arcpy.da.FeatureClassToNumPyArray(select_line,[linefield])
        linefieldArr = np.array([item[0] for item in linearray])
        if len(linefieldArr > 0):  
            maxValue = max(linefieldArr)
            with arcpy.da.UpdateCursor(select_line, [linefield]) as cursor:
                for row in cursor:
                    fieldvalue = float(row[0])
                    if fieldvalue < maxValue: ##This step reomve the small spurious ploygons as well
                        cursor.deleteRow()
            del cursor, row

            ##Get the point_x and point_y
            for rowpoint in arcpy.da.SearchCursor(select_point, ["SHAPE@XY"]):
                pointx, pointy = rowpoint[0]      
            del rowpoint
        
            with arcpy.da.SearchCursor(select_line, "SHAPE@") as cursor:
                for row in cursor:
                    firstPnt = row[0].firstPoint
                    startx = firstPnt.X
                    starty = firstPnt.Y

                    endPnt = row[0].lastPoint
                    endx = endPnt.X
                    endy = endPnt.Y

                    if starty==endy or startx==endx:
                        if starty == endy:
                            y1 = pointy + distance
                            y2 = pointy - distance
                            x1 = pointx
                            x2 = pointx
                        if startx == endx:
                            y1 = pointy
                            y2 = pointy 
                            x1 = pointx + distance
                            x2 = pointx - distance     
                    else:
                        m = ((starty - endy)/(startx - endx)) #get the slope of the line
                        negativereciprocal = -1*((startx - endx)/(starty - endy))    #get the negative reciprocal
                        if m > 0:
                            if m >= 1:
                                y1 = negativereciprocal*(distance)+ pointy
                                y2 = negativereciprocal*(-distance) + pointy
                                x1 = pointx + distance
                                x2 = pointx - distance
                            if m < 1:
                                y1 = pointy + distance
                                y2 = pointy - distance
                                x1 = (distance/negativereciprocal) + pointx
                                x2 = (-distance/negativereciprocal)+ pointx           
                        if m < 0:
                            if m >= -1:
                                y1 = pointy + distance
                                y2 = pointy - distance
                                x1 = (distance/negativereciprocal) + pointx
                                x2 = (-distance/negativereciprocal)+ pointx     
                            if m < -1:
                                y1 = negativereciprocal*(distance)+ pointy
                                y2 = negativereciprocal*(-distance) + pointy
                                x1 = pointx + distance
                                x2 = pointx - distance
                    array = arcpy.Array([arcpy.Point(x1,y1),arcpy.Point(x2, y2)])
                    polyline = arcpy.Polyline(array)
                    pntID = pointfieldArr[i]
                    #pntID = arr[row][2]
                    #segID = arr[row][4]
                    new_line_cursor.insertRow([polyline, pntID])

            del cursor, row  

    del new_line_cursor


    ##Delete the created in-memory dataset
    try:
        arcpy.Delete_management (pointsbuf)
        arcpy.Delete_management (cliplines)
        arcpy.Delete_management (singlecliplines)
        arcpy.Delete_management (select_point)
        arcpy.Delete_management (select_line)
    except:
        pass
    return perpendicular_line

 
##Main program
# Script arguments

##To run the tool in ArcGIS as a geoprocessing tool. Obtain the parameters from the interface
InputDEM = arcpy.GetParameterAsText(0)
EleThreshold = arcpy.GetParameter(1)
StreamThresholdKM2 = arcpy.GetParameter(2)
TributaryThresholdKM2 = arcpy.GetParameter(3)
#MinCirqueAreaKm2 = arcpy.GetParameter(4)
#MaxCirqueAreaKm2 = arcpy.GetParameter(5)
Min_Cir_Index = arcpy.GetParameter(4)
#Min_Plan_Clos = arcpy.GetParameter(5)
Min_steep_slope_percent = arcpy.GetParameter(5)
Min_gentle_slope_percent = arcpy.GetParameter(6)
Max_backwall_slope =  arcpy.GetParameter(7)
Flexible_backwall_slope =  arcpy.GetParameter(8)

Out3DProfiles = arcpy.GetParameterAsText(9)
OutKnickpoints = arcpy.GetParameterAsText(10) ##the turning points for cirques; 
OutCrossSections = arcpy.GetParameterAsText(11) ##cross sections
OutCirques = arcpy.GetParameterAsText(12)

Min_backwall_slope = int(Max_backwall_slope) - int(Flexible_backwall_slope)

'''
##To run the tool in python by assign parameters in the file, not from the ArcGIS interface
##This can speed up the runing time
InputDEM = "c:\\testdata\\tsDEM1.tif"
EleThreshold = 3400
StreamThresholdKM2 = 0.03
TributaryThresholdKM2 = 0.5
#MinCirqueAreaKm2 = 0.4
#MaxCirqueAreaKm2 = 2.0
Min_Cir_Index = 0.7
Min_steep_slope_percent = 10
Min_gentle_slope_percent = 5
Max_backwall_slope =  27
Min_backwall_slope =  23

Out3DProfiles = "c:\\testdata\\ts3dprofile265up.shp"
OutKnickpoints = "c:\\testdata\\tsknickpoints265up.shp" ##the turning points for cirques
OutCirques = "c:\\testdata\\tsCirques265up.shp"
'''
#Record the start time
start_time = time.time()

arcpy.Delete_management("in_memory")


cellsize = arcpy.GetRasterProperties_management(InputDEM,"CELLSIZEX")
cellsize_int = int(float(cellsize.getOutput(0)))

StreamThreshold = int(float(StreamThresholdKM2) * 1e6 / (cellsize_int * cellsize_int))
TributaryThreshold = int(float(TributaryThresholdKM2) * 1e6 / (cellsize_int * cellsize_int))

#minCirqueArea = MinCirqueAreaKm2 *1e6
#maxCirqueArea = MaxCirqueAreaKm2 *1e6
#set the minCirque area as the min contrbuting area for a stream
minCirqueArea = float(StreamThresholdKM2) * 1e6

###temporay files
#outPolygon =  "in_memory\\polygon.shp"
MaxFccTable = "in_memory\\MaxFccTable"
StreamOrderTable = "in_memory\\StreamOrderTable"
CleanStream = "in_memory\\CleanStream"
tmpoutStream = "in_memory\\tmpoutStream"
TmpStream = "in_memory\\TmpStream"

###Step 1: Stream network
arcpy.AddMessage("Step 1: Stream extraction...")
#print "Step 1: Stream extraction..."
## Extract the DEM only higher than the EleThreshold
OutBnd = Con(Raster(InputDEM)> EleThreshold,1)
outExtractDEM = ExtractByMask(InputDEM, OutBnd)

##create a slope raster for the original DEM
outslope = Slope(outExtractDEM)

#Hydro analysis
fillDEM =Fill(outExtractDEM)  ##Fill the sink first
fdir = FlowDirection(fillDEM,"NORMAL") ##Flow direction
facc = FlowAccumulation(fdir) ##Flow accmulation
outGreaterThan = Con(facc > StreamThreshold, 1,0)  ##Determine the highest flowaccumuation part

# Process: Stream analysis
outStreamLink = StreamLink(outGreaterThan, fdir)
StreamToFeature(outStreamLink, fdir, TmpStream, "SIMPLIFY")

# Process: Zonal Statistics as Table
ZonalStatisticsAsTable(outStreamLink, "VALUE", facc, MaxFccTable, "DATA", "MAXIMUM")

# Process: Stream Order ##It is not necessary
outStreamOrder = StreamOrder(outStreamLink, fdir, "STRAHLER")

# Process: Zonal Statistics as Table (2)
ZonalStatisticsAsTable(outStreamLink, "VALUE", outStreamOrder, StreamOrderTable, "DATA", "MAXIMUM")

# Process: Join Field
arcpy.JoinField_management(TmpStream, "grid_code", MaxFccTable, "Value", "MAX")  ##Join to get the flow accumulation value
arcpy.JoinField_management(TmpStream, "grid_code", StreamOrderTable, "Value", "MAX") ##Join to get the stream order value  THe field name is MAX_1

#arcpy.CopyFeatures_management(TmpStream, "c:\\testdata\\outStreamAll.shp")

##Step 2: clean streams                
arcpy.AddMessage("Step 2: Filtering streams...")
#print "Step 2: Clean streams..."
###This TmpStream already have a to_node in the attibute table, so that it can be used to make the decision
lineArray = arcpy.da.FeatureClassToNumPyArray(TmpStream,['OID@','to_node','MAX'])
tonode = np.array([item[1] for item in lineArray])
uniquenode = np.unique(tonode)
lineid = [] ##Record the id for deletion
for i in range(len(uniquenode)):
    selArr = lineArray[tonode == uniquenode[i]]
    fcclist = []
    if len(selArr) > 1: ##Sometimes having more than two end points
        for j in range(len(selArr)):
            fcclist.append(selArr[j][2])

        numselected = len(fcclist)
        while numselected > 1:
            minfcc = min(fcclist)
            if minfcc < TributaryThreshold:
                for j in range(len(selArr)):
                    if selArr[j][2] == minfcc: ##Remove the smaller fcc one
                        lineid.append(selArr[j][0])
                        fcclist.pop(j)
                        selArr = np.delete(selArr, j)##Remove this one and loop to delete others
                        numselected = len(selArr)
                        break ##Only remove one each time

            else: ##quit the loop if all minfcc are larger than the theshold?? Try to remove more based on the ratio
                break

##Delete the line marked for deletion
with arcpy.da.UpdateCursor(TmpStream, "OID@") as cursor:
    for row in cursor:
        if int(row[0]) in lineid:
            cursor.deleteRow()     
del cursor, row             

##Clean extralines based on the end points intersection 09/24/2020
cleanextralineswithtopology(TmpStream,tmpoutStream, 'MAX', TributaryThreshold)  ## clean the extra lines before dissolving
arcpy.Dissolve_management(tmpoutStream, CleanStream, '#', 'MAX MAX;MAX_1 MIN', 'SINGLE_PART', 'UNSPLIT_LINES') 

#arcpy.CopyFeatures_management(CleanStream, "c:\\testdata\\outCleanStreams.shp")

###Step 3: Find the knickpoints for the first order streams only
arcpy.AddMessage("Step 3: Find the potential threshold points for the first order streams...")
#print "Step 3: Find the knickpoints for the first order streams..."

FirstOrderStream = "in_memory\\FirstOrderStream"
arcpy.Select_analysis(CleanStream, FirstOrderStream, '"MIN_MAX_1" < 2')

FirstOrderStream3D = "in_memory\\FirstOrderStream3D"
##use filled DEM and 3*cellsize as spacing; save the 3d feature as one output: out3DProfiles
arcpy.InterpolateShape_3d(fillDEM, FirstOrderStream, FirstOrderStream3D, cellsize_int*3) 

Knickpoints = arcpy.CreateFeatureclass_management("in_memory", "Knickpoints","POINT","", "DISABLED", "DISABLED", FirstOrderStream3D)
arcpy.AddField_management(Knickpoints, "LineID", "Long")
arcpy.AddField_management(Knickpoints, "Ratio", "DOUBLE")
arcpy.AddField_management(Knickpoints, "Ele", "DOUBLE")

LineIDlist = []
with arcpy.da.SearchCursor(FirstOrderStream3D, ["OID@","SHAPE@", "SHAPE@Length"]) as cursor:
    for row in cursor: ##Loop for each line
        PointX = []
        PointY = []
        LengthfromStart = []
        PointZ = []
        lineId = row[0]
        lineLength = float(row[2])
        cumLength = 0
        for part in row[1]:
            # Step through each vertex in the feature
            #start = 0
            #startx = 0
            pntCount = 0
            cumLength = 0
            segmentdis = 0
            for pnt in part:
                if pnt:
                    if pntCount > 0:
                        cumLength += Dist(startx, starty, pnt.X, pnt.Y) 

                    PointX.append(pnt.X)
                    PointY.append(pnt.Y)
                    PointZ.append(pnt.Z)
                    LengthfromStart.append(cumLength)

                    startx = pnt.X
                    starty = pnt.Y
                    pntCount += 1

        #t_points, t_ratios = turning_points_ConvexAngle(LengthfromStart, PointZ, 10, 200)
        #t_points, t_ratios = turning_points_ConvexAngle(LengthfromStart, PointZ, 3, cellsize_int*3)
        logstreamlength = np.log(LengthfromStart[1:])
        t_points, t_ratios = turning_points_RDP(logstreamlength, PointZ[1:], 3, cellsize_int*5) ##only top 3 turing points should be enough
        #t_points, t_ratios = turning_points_RDE(LengthfromStart, PointZ, 10, 200)
        #t_points, t_ratios = turning_points(LengthfromStart, PointZ, turning_points = 5, cluster_radius = cellsize_int*5)
        #t_points, t_ratios = turning_points_ConvexAngle(LengthfromStart, PointZ, 10, 200)
        
        #print t_points
        #print t_ratios
        
        ##Make the turn point layer
        if len(t_points) > 0:
            CursorPnts = arcpy.da.InsertCursor(Knickpoints, ['SHAPE@', 'LineID', 'Ratio', 'Ele'])
            LineIDlist.append(lineId)
            for i in range(len(t_points)):
                #idx = t_points[i] + 1 ##the idx should plus 1 because the t_points are only indexed except for the start and end points
                #arcpy.AddMessage("idx is: " + str(idx))
                idx = t_points[i] ##+ 1 ##do not add 1 for the RDP method
                Ratio = t_ratios[i]
                Pnt = arcpy.Point(PointX[idx], PointY[idx])
                CursorPnts.insertRow((Pnt, lineId, Ratio, PointZ[idx])) ##Add Z as a point attribute here
            del CursorPnts ##Delete cursor object
del cursor, row  ##Delete cursor object

#arcpy.CopyFeatures_management(FirstOrderStream3D, "c:\\test\\FirstOrderStream3D.shp")
#arcpy.CopyFeatures_management(Knickpoints, "c:\\test\\Knickpoints.shp")
##create the cross sections for the threshold points along the first order stream
CrossSections = cross_sections(Knickpoints, "LineID", FirstOrderStream, "MAX_MAX", 50, 150)

###Step 4: Select the turning points that likely represents the cirque mouth
arcpy.AddMessage("Step 4: Delineating potential cirque outlines..." )
#print "Step 4: Delineating potential cirque outlines..."
##Need to change to the delineation for cross sections
##regenerate the streamnetwork (a more detailed stream network)
Streams = "in_memory\\Streams"
FCselected = "in_memory\\FCselected"  ##Set a in_memory file for each moraine feature
Singlepoint = "in_memory\\Singlepoint"
Singlestream = "in_memory\\Singlestream"
SingleWs = "in_memory\\SingleWs"
singleBND = "in_memory\\singleBND"
singleLine = "in_memory\\singleLine"
Linepoints = "in_memory\\Linepoints"
pntWsBnd = "in_memory\\pntWsBnd"
smoothBnd = "in_memory\\smoothBnd"

outGreaterThan = Con(facc > 10, 1,0)  ##Determine the highest flowaccumuation part
# Process: Stream Link
outStreamLink = StreamLink(outGreaterThan, fdir)
# Process: Stream to Feature
StreamToFeature(outStreamLink, fdir, Streams, "NO_SIMPLIFY")

##Here are the new codes 9/29/2021 Yingkui Li
selectedTCrossSections = arcpy.CreateFeatureclass_management("in_memory", "selectedTCrossSections","POLYLINE","","","",Knickpoints)

slopeBnds = arcpy.CreateFeatureclass_management("in_memory", "slopeBnds","POLYGON","","","",Knickpoints)

arcpy.AddField_management(slopeBnds, "LineID", "Long")
#arcpy.AddField_management(slopeBnds, "PointID", "Long")
arcpy.AddField_management(slopeBnds, "CirIndex", "DOUBLE")
arcpy.AddField_management(slopeBnds, "SteepPct", "DOUBLE")
arcpy.AddField_management(slopeBnds, "GentlePct", "DOUBLE")
#remove plan_clos becasue it takes a long time and this value can be derived later using ACME tool; 07/29/2022 Yingkui Li
#arcpy.AddField_management(slopeBnds, "Plan_Clos", "DOUBLE")
arcpy.AddField_management(slopeBnds, "Slp_Range", "DOUBLE")

# Set the progressor
arcpy.SetProgressor("step", "Select the turning points likely representing the cirque outlets", 0, len(LineIDlist), 1)

##Need to use the extent of the fdir/facc/filldem to make sure the watershed is correct!!!!
arcpy.env.extent = fdir
arcpy.env.snapRaster = fdir ##Need to make sure the generated raster align with the cell alignment of the fdir
smoothtolerance = cellsize_int * 5


#TotalBND = arcpy.CreateFeatureclass_management("in_memory", "TotalBND","POLYGON","","","",InputFC)
linearray = arcpy.da.FeatureClassToNumPyArray(CrossSections,"LineID")
lineIDArr = np.array([item[0] for item in linearray])
count = len(lineIDArr)

FcID = arcpy.Describe(CrossSections).OIDFieldName

for ifeature in range (count):
    arcpy.AddMessage("Generating cirque "+str(ifeature + 1)+" of "+str(count))
    query = FcID +" = "+str(ifeature+1)
    arcpy.Select_analysis(CrossSections, FCselected, query)

    arcpy.Intersect_analysis([Streams, FCselected], Singlepoint, "", "0", "point")
    pntcountResult = arcpy.GetCount_management(Singlepoint)
    pntcount = int(pntcountResult.getOutput(0))

    if (pntcount == 0): ##if no intersect points, try one more time: use the buffer to get the intersection points for another time
        tmpbuf = "in_memory\\tmpbuf"
        arcpy.Buffer_analysis(FCselected, tmpbuf, "5 Meters")
        arcpy.Intersect_analysis([Streams, tmpbuf], Singlepoint, "#", "#", "POINT")
        pntcountResult = arcpy.GetCount_management(Singlepoint)
        pntcount = int(pntcountResult.getOutput(0))
        #arcpy.AddMessage("the number of point is:" + str(pntcount))
    
    if (pntcount > 0):
        #Calculate Watershed
        outPour = SnapPourPoint(Singlepoint, facc, 0)
        outWs1 = Watershed(fdir, outPour)
        outWs = Con(outWs1 >= 0, 1)  ##Determine the highest flowaccumuation part and set others as Nodata
        
        arcpy.RasterToPolygon_conversion(outWs, SingleWs, "NO_SIMPLIFY", "VALUE")
        
        arcpy.Clip_analysis(Streams, SingleWs, Singlestream)
        
        singleDEM = ExtractByMask(fillDEM, outWs)
        ext_slope = ExtractByMask(outslope, outWs)

        SlopeBnd, status = BoundaryRefineBySlope(singleDEM, ext_slope, Singlepoint, int(Max_backwall_slope), int(Min_backwall_slope), 120)
        slopepolycountResult = arcpy.GetCount_management(SlopeBnd)
        slopepolycount = int(slopepolycountResult.getOutput(0))
        if (slopepolycount>0 and status): ##meet the cirque standards
            arcpy.cartography.SmoothPolygon(SlopeBnd, smoothBnd, "PAEK", smoothtolerance)##Do not smooth the polygon
            #arcpy.cartography.SmoothPolygon(SlopeBnd, smoothBnd, "BEZIER_INTERPOLATION", 0)##Do not smooth the polygon
            ##Do the slope percentage analysis
            ext_slope2 = ExtractByMask(ext_slope, smoothBnd)
            arr = arcpy.RasterToNumPyArray(ext_slope2)
            totalPixels = (arr > 0).sum()
            if totalPixels > 0:
                steepPixels = (arr > 31).sum() ##Calculate the pixels of the slope > 31 based on Evans and Cox (1974)
                steep_percent = float(steepPixels)/ float(totalPixels) * 100
                gentlePixels = (arr < 20).sum() - (arr <= 0).sum() ##Calculate the pixels of the slope < 20 based on Evans and Cox (1974)
                gentle_percent = float(gentlePixels)/ float(totalPixels) * 100
                max_slp = min(np.nanmax(arr),90)
                #arcpy.AddMessage(str(max_slp))
                min_slp = max(np.nanmin(arr),0)
                #arcpy.AddMessage(str(min_slp))
                slprange = max_slp-min_slp
                #arcpy.AddMessage(str(slprange))

                if (steep_percent > float(Min_steep_slope_percent) and gentle_percent > float(Min_gentle_slope_percent)): ##float(Min_Plan_Clos)): ##Only count the steep percent of >30% within the watershed, need to make a parameter for this value
                    #arcpy.CopyFeatures_management(SlopeBnd, smoothBnd)
                    ##Remove the potential spurious polygons and recalculate the CirIndex
                    arcpy.AddField_management(smoothBnd, "LineID", "Long")
                    #arcpy.AddField_management(smoothBnd, "PointID", "Long")
                    arcpy.AddField_management(smoothBnd, "CirIndex", "DOUBLE")
                    arcpy.AddField_management(smoothBnd, "SteepPct", "DOUBLE")
                    arcpy.AddField_management(smoothBnd, "GentlePct", "DOUBLE")
                    #arcpy.AddField_management(smoothBnd, "Plan_Clos", "DOUBLE")
                    arcpy.AddField_management(smoothBnd, "Slp_Range", "DOUBLE")

                    ##Calculate Plan closure
                    #cirqueDEM = ExtractByMask(singleDEM, smoothBnd)
                    #planClos = plan_clos(cirqueDEM)
                    #arcpy.AddMessage("Plan Closure is " + str(planClos))
                    #print planClos
                    
                    meetcount2  = 0
                    #with arcpy.da.UpdateCursor(smoothBnd, ['SHAPE@LENGTH', 'SHAPE@AREA', 'LineID', 'CirIndex', 'SteepPct', 'GentlePct', 'Plan_Clos', 'Slp_Range']) as cursor:
                    with arcpy.da.UpdateCursor(smoothBnd, ['SHAPE@LENGTH', 'SHAPE@AREA', 'LineID', 'CirIndex', 'SteepPct', 'GentlePct', 'Slp_Range']) as cursor:
                        for row in cursor:
                            length = float(row[0])
                            area = float(row[1])
                            Circularindex = (4.0 * 3.14 * area) / (length * length)
                            #if ((area > minCirqueArea) and (area < maxCirqueArea) and (Circularindex > float(Min_Cir_Index))):
                            #if ((area > minCirqueArea) and (Circularindex > float(Min_Cir_Index)) and planClos > 0):
                            if ((area > minCirqueArea) and (Circularindex > float(Min_Cir_Index)) ):
                                meetcount2  += 1
                                row[2] = lineIDArr[ifeature] #LineID
                                #row[3] = i+1  ##PointID
                                row[3] = Circularindex
                                #row[4] = planClos
                                row[4] = steep_percent
                                row[5] = gentle_percent
                                #row[6] = planClos
                                row[6] = slprange
                                #print row[7]
                                cursor.updateRow(row)
                            else:
                                cursor.deleteRow()
                    del cursor, row
                    if meetcount2 > 0:                
                        arcpy.AddMessage("Add one cirque outline")
                        #print "Add one cirque outline"
                        arcpy.Append_management(smoothBnd, slopeBnds, "NO_TEST") ##save the ws boundary that meet the criteria
                        arcpy.Append_management(FCselected, selectedTCrossSections, "NO_TEST") ##save the single point that meet the criteria
                
        
    arcpy.SetProgressorPosition()

##Smooth polygons
#arcpy.CopyFeatures_(TotalBND, OutCirques, "PAEK", smoothtolerance)
#arcpy.CopyFeatures_management(TotalBND, OutCirques)
#polycountResult = arcpy.GetCount_management(slopeBnds)
#polycount = int(polycountResult.getOutput(0))
#if polycount > 0:
#
#else:
#    arcpy.AddMessage("Fail to delineate the polygons!")

arcpy.SpatialJoin_analysis(FirstOrderStream3D, selectedTCrossSections, Out3DProfiles, "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "1 Meters", "#")
arcpy.CopyFeatures_management(slopeBnds, OutCirques)
arcpy.CopyFeatures_management(selectedTCrossSections, OutCrossSections)
##Get the threshold points by spatial join
arcpy.SpatialJoin_analysis(Knickpoints, selectedTCrossSections, OutKnickpoints, "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "1 Meters", "#")

##Finish the process and display the processing time
arcpy.AddMessage("Complete the cirque delineation!!!")
#print "Complete the cirque delineation!!!"

#process_time = time.time() - start_time
#minutes = int(process_time / 60)
#seconds = process_time - minutes * 60
#print("---%s minutes and %s seconds ---" % (minutes, seconds))


##Delete intermidiate data
arcpy.Delete_management("in_memory") ### Empty the in_memory

'''
##The following are the old codes

singleLine = "in_memory\\singleLine"
Linepoints = "in_memory\\Linepoints"
Singlepoint = "in_memory\\Singlepoints"
pntWsBnd = "in_memory\\pntWsBnd"
smoothBnd = "in_memory\\smoothBnd"


selectedTpoints = arcpy.CreateFeatureclass_management("in_memory", "selectedTpoints","POINT","","","",Knickpoints)

slopeBnds = arcpy.CreateFeatureclass_management("in_memory", "slopeBnds","POLYGON","","","",Knickpoints)
arcpy.AddField_management(slopeBnds, "LineID", "Long")
arcpy.AddField_management(slopeBnds, "PointID", "Long")
arcpy.AddField_management(slopeBnds, "CirIndex", "DOUBLE")
arcpy.AddField_management(slopeBnds, "SteepPct", "DOUBLE")
arcpy.AddField_management(slopeBnds, "GentlePct", "DOUBLE")
arcpy.AddField_management(slopeBnds, "Plan_Clos", "DOUBLE")
arcpy.AddField_management(slopeBnds, "Slp_Range", "DOUBLE")
# Set the progressor
arcpy.SetProgressor("step", "Select the turning points likely representing the cirque outlets", 0, len(LineIDlist), 1)

##Need to use the extent of the fdir/facc/filldem to make sure the watershed is correct!!!!
arcpy.env.extent = fdir
arcpy.env.snapRaster = fdir ##Need to make sure the generated raster align with the cell alignment of the fdir
smoothtolerance = cellsize_int * 5
linecount = len(LineIDlist)
ClipBndCount = 0
##print linecount
for id in range(linecount):
#for id in range(50): ##just do 200 lines
    # Update the progressor label for current shapefile
    arcpy.SetProgressorLabel("Loading {0}...".format(id))
    
    ##Select all turn points of each line
    query = "LineID = " + str(LineIDlist[id])
    arcpy.AddMessage("Processing first order stream: " + str(id+1)+ " out of " + str(linecount))
    #print "Processing first order stream: " + str(id)+ " out of " + str(linecount)
    arcpy.Select_analysis(Knickpoints, Linepoints, query)
    countResult = arcpy.GetCount_management(Linepoints)
    count = int(countResult.getOutput(0))
    #arcpy.AddMessage("The number of points: " + str(count))
    PointID = arcpy.Describe(Linepoints).OIDFieldName

    StreamID = arcpy.Describe(FirstOrderStream3D).OIDFieldName
    query = StreamID +" = " + str(LineIDlist[id])
    arcpy.Select_analysis(FirstOrderStream3D, singleLine, query)
    

    ###First to select potential cirque theshold points using the area and circular index based on the watershed boundary
    ### This step don't need to just pickup the maximum circuar one, rather save all for later refinement. Also, need to save the boundary as well
    ### After this step, run the slope based analysis (include the steep slope and gentle slope percentages) for the refined points and watershed boundaries
    for i in range (count):
    #for i in range (2):
        query = PointID +" = "+ str(i+1)  ##THe FID starts from 1 for the in_memory data; Which is different from shape (FID from 0) and feature class in the geodatabase (objectID from 1)
        arcpy.Select_analysis(Linepoints, Singlepoint, query)
        outSnapPour = SnapPourPoint(Singlepoint, facc, cellsize_int)##single point watershed
        outpntWs = Watershed(fdir, outSnapPour)
        arcpy.RasterToPolygon_conversion(outpntWs, pntWsBnd, "#","VALUE")
        ##Simply select potential cirque theshold points using the area and circular index based on the watershed boundary
        polycountResult = arcpy.GetCount_management(pntWsBnd)
        polycount = int(polycountResult.getOutput(0))
        if (polycount>0): ##meet the cirque standards
            meetcount  = 0
            with arcpy.da.UpdateCursor(pntWsBnd, ['SHAPE@LENGTH', 'SHAPE@AREA']) as cursor:
                for row in cursor:
                    length = float(row[0])
                    area = float(row[1])
                    #Circularindex = (4.0 * 3.14 * area) / (length * length)
                    ##remove the constrain for the max cirque area
                    if ((area > minCirqueArea)): ##and (area < maxCirqueArea*1.2)): ##and (Circularindex > float(Min_Cir_Index))): ##This step reomve the small spurious ploygons as well
                        meetcount  += 1
                    else:
                        cursor.deleteRow()
            del cursor, row
            if meetcount > 0:                
                ##refine the boundary by slope analysis
                ext_slope = ExtractByMask(outslope, pntWsBnd)
                ext_DEM = ExtractByMask(fillDEM, pntWsBnd)
                #arr = arcpy.RasterToNumPyArray(ext_slope)
                #totalPixels = (arr > 0).sum()
                #steepPixels = (arr > 31).sum() ##Calculate the pixels of the slope > 31 based on Evans and Cox (1974)
                #steep_percent = float(steepPixels)/ float(totalPixels) * 100
                #gentlePixels = (arr < 20).sum() - (arr <= 0).sum() ##Calculate the pixels of the slope < 20 based on Evans and Cox (1974)
                #gentle_percent = float(gentlePixels)/ float(totalPixels) * 100

                #if (steep_percent > float(Min_steep_slope_percent) and gentle_percent > float(Min_gentle_slope_percent)): ##Only count the steep percent of >30% within the watershed, need to make a parameter for this value
                ##Run the slope analysis to get rid of the gentle area above the headwall
                SlopeBnd, status = BoundaryRefineBySlope(ext_DEM, ext_slope, singleLine, int(Max_backwall_slope), int(Min_backwall_slope), 120)
                #SlopeBnd, status = BoundaryRefineBySlopeAndSkyline(ext_DEM, ext_slope, singleLine, 1, int(Max_backwall_slope), int(Min_backwall_slope), 120)
                
                slopepolycountResult = arcpy.GetCount_management(SlopeBnd)
                slopepolycount = int(slopepolycountResult.getOutput(0))
                #print "slopepolycount is: " + str(slopepolycount)

                if (slopepolycount>0 and status): ##meet the cirque standards
                    arcpy.cartography.SmoothPolygon(SlopeBnd, smoothBnd, "PAEK", smoothtolerance)##Do not smooth the polygon
                    #arcpy.cartography.SmoothPolygon(SlopeBnd, smoothBnd, "BEZIER_INTERPOLATION", 0)##Do not smooth the polygon
                    ##Do the slope percentage analysis
                    ext_slope2 = ExtractByMask(ext_slope, smoothBnd)
                    arr = arcpy.RasterToNumPyArray(ext_slope2)
                    totalPixels = (arr > 0).sum()
                    if totalPixels > 0:
                        steepPixels = (arr > 31).sum() ##Calculate the pixels of the slope > 31 based on Evans and Cox (1974)
                        steep_percent = float(steepPixels)/ float(totalPixels) * 100
                        gentlePixels = (arr < 20).sum() - (arr <= 0).sum() ##Calculate the pixels of the slope < 20 based on Evans and Cox (1974)
                        gentle_percent = float(gentlePixels)/ float(totalPixels) * 100
                        max_slp = min(np.nanmax(arr),90)
                        #arcpy.AddMessage(str(max_slp))
                        min_slp = max(np.nanmin(arr),0)
                        #arcpy.AddMessage(str(min_slp))
                        slprange = max_slp-min_slp
                        #arcpy.AddMessage(str(slprange))

                        ##Calculate Plan closure
                        cirqueDEM = ExtractByMask(ext_DEM, smoothBnd)
                        planClos = plan_clos(cirqueDEM)
                        arcpy.AddMessage("Plan Closure is" + str(planClos))
                        #print planClos
                        if (steep_percent > float(Min_steep_slope_percent) and gentle_percent > float(Min_gentle_slope_percent) and planClos > 0): ##float(Min_Plan_Clos)): ##Only count the steep percent of >30% within the watershed, need to make a parameter for this value
                            #arcpy.CopyFeatures_management(SlopeBnd, smoothBnd)
                            ##Remove the potential spurious polygons and recalculate the CirIndex
                            arcpy.AddField_management(smoothBnd, "LineID", "Long")
                            arcpy.AddField_management(smoothBnd, "PointID", "Long")
                            arcpy.AddField_management(smoothBnd, "CirIndex", "DOUBLE")
                            arcpy.AddField_management(smoothBnd, "SteepPct", "DOUBLE")
                            arcpy.AddField_management(smoothBnd, "GentlePct", "DOUBLE")
                            arcpy.AddField_management(smoothBnd, "Plan_Clos", "DOUBLE")
                            arcpy.AddField_management(smoothBnd, "Slp_Range", "DOUBLE")
                            
                            meetcount2  = 0
                            with arcpy.da.UpdateCursor(smoothBnd, ['SHAPE@LENGTH', 'SHAPE@AREA', 'LineID', 'PointID','CirIndex', 'SteepPct', 'GentlePct', 'Plan_Clos', 'Slp_Range']) as cursor:
                                for row in cursor:
                                    length = float(row[0])
                                    area = float(row[1])
                                    Circularindex = (4.0 * 3.14 * area) / (length * length)
                                    #if ((area > minCirqueArea) and (area < maxCirqueArea) and (Circularindex > float(Min_Cir_Index))):
                                    if ((area > minCirqueArea) and (Circularindex > float(Min_Cir_Index))):
                                        meetcount2  += 1
                                        row[2] = LineIDlist[id] #LineID
                                        row[3] = i+1  ##PointID
                                        row[4] = Circularindex
                                        #row[4] = planClos
                                        row[5] = steep_percent
                                        row[6] = gentle_percent
                                        row[7] = planClos
                                        row[8] = slprange
                                        #print row[7]
                                        cursor.updateRow(row)
                                    else:
                                        cursor.deleteRow()
                            del cursor, row
                            if meetcount2 > 0:                
                                arcpy.AddMessage("Add one cirque outline")
                                #print "Add one cirque outline"
                                arcpy.Append_management(smoothBnd, slopeBnds, "NO_TEST") ##save the ws boundary that meet the criteria
                                arcpy.Append_management(Singlepoint, selectedTpoints, "NO_TEST") ##save the single point that meet the criteria
                
        
    arcpy.SetProgressorPosition()

arcpy.ResetProgressor()

#arcpy.CopyFeatures_management(slopeBnds, "c:\\test\\tmppntslopeBnds.shp")
#arcpy.CopyFeatures_management(selectedTpoints, "c:\\test\\tmpselectpoints.shp")

arcpy.SpatialJoin_analysis(FirstOrderStream3D, selectedTpoints, Out3DProfiles, "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "1 Meters", "#")
arcpy.CopyFeatures_management(slopeBnds, OutCirques)
arcpy.CopyFeatures_management(selectedTpoints, OutKnickpoints)

##Finish the process and display the processing time
arcpy.AddMessage("Complete the cirque delineation!!!")
#print "Complete the cirque delineation!!!"

process_time = time.time() - start_time
minutes = int(process_time / 60)
seconds = process_time - minutes * 60
#print("---%s minutes and %s seconds ---" % (minutes, seconds))


##Delete intermidiate data
arcpy.Delete_management("in_memory") ### Empty the in_memory
'''




