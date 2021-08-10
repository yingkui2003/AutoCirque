#-------------------------------------------------------------------------------
# Name: CirquePotentialOutletPoints.py
# Purpose: This tool generates the potential cirque outlet points based on DEM
# and the min elevation. The user can then select the outlet points to help
# delineate the cirques.
# 
# Author: Dr. Yingkui Li
# Created:     09/21-12/29/2020
# Department of Geography, University of Tennessee
# Knoxville, TN 37996
#-------------------------------------------------------------------------------

# Import arcpy module
import arcpy, sys
from arcpy import env
from arcpy.sa import *
#import numpy
import numpy as np

arcpy.env.overwriteOutput = True
arcpy.env.XYTolerance= "0.01 Meters"

arcpy.Delete_management("in_memory") ### Empty the in_memory

try:
    if arcpy.CheckExtension("Spatial")=="Available":
        arcpy.CheckOutExtension("Spatial")
    else:
        print "not extension available"
except:
    print "unable to check out extension"

try:
    if arcpy.CheckExtension("3D")=="Available":
        arcpy.CheckOutExtension("3D")
    else:
        print "not extension available"
except:
    print "unable to check out extension"

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
        if faccThreshold > 0:
            maxfacc = faccThreshold
            
        lineid = [] ##Record the id for deletion
        for i in range(len(lineArray)):
            linetonode = lineArray[i][2]
            nodecount = np.count_nonzero(uniquetonode == linetonode)
            if nodecount == 1 and not (linetonode in fromnode) and lineArray[i][3] < maxfacc: ###only one tonode except for the highest facc section
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
    upsl = upydiff/upxdiff * x[1:-1]
    dnsl = dnydiff/dnxdiff * x[1:-1]
    xbothside = x[2:] - x[0:-2]
    ybothside = y[2:] - y[0:-2]
    slbothside = ybothside / xbothside * x[1:-1]
    slbothside[slbothside == 0] = -0.001
    Ratio = (dnsl - upsl) / slbothside
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

    while len(t_points) < turning_points and len(turn_point_idx) > 0:
        Ratio = SL[turn_point_idx[0]]
        if Ratio < 1.0:
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


###rdp only positive distance!!! for turning point detection

def Knickpoints_rdp(points, epsilon, turn_points, angles):
    # get the start and end points
    start = np.tile(np.expand_dims(points[0], axis=0), (points.shape[0], 1))
    end = np.tile(np.expand_dims(points[-1], axis=0), (points.shape[0], 1))
    linedist = Dist(start[0][0],start[0][1],end[0][0],end[0][1])
    dist_point_to_line = np.cross(end - start, points - start, axis=-1) / np.linalg.norm(end - start, axis=-1)
    max_idx = np.argmax(np.abs(dist_point_to_line))
    max_value = dist_point_to_line[max_idx]/linedist
    
    if abs(max_value) > epsilon:  ##the distance is at least 1 m from the line
        if max_value > 0:
             ##Calculate the angle of this maximum point to the both side
             centerpnt =  points[max_idx]
             angle = math.atan2(end[0][1] - centerpnt[1], end[0][0] - centerpnt[0]) - math.atan2(start[0][1] - centerpnt[1], start[0][0] - centerpnt[0]);
             turn_points.append(points[max_idx])
             angles.append(180 - (360 + angle/3.14159*180))

        partial_results_left = Knickpoints_rdp(points[:max_idx+1], epsilon, turn_points,angles)
        partial_results_right = Knickpoints_rdp(points[max_idx:], epsilon, turn_points,angles)
        
def turning_points_RDP(streamLength, streamZ, turning_points = 10, cluster_radius = 200):
    stream_points = np.concatenate([np.expand_dims(streamLength, axis=-1), np.expand_dims(streamZ, axis=-1)], axis=-1)

    epsilon = 0.01
    turn_points = []
    turn_angles = []

    Knickpoints_rdp(stream_points, epsilon, turn_points, turn_angles)

    if len(turn_points) < 1:
        return [], []
    
    new_x = np.array(turn_points)[:,0]
    turn_point_idx = np.argsort(turn_angles)[::-1]

    t_pointsID = []
    t_angles = []

    while len(t_pointsID) < turning_points and len(turn_point_idx) > 0:
        angle = turn_angles[turn_point_idx[0]]
        if angle < 3:
            break
        else:
            t_pointsID += [turn_point_idx[0]]
            t_angles.append(turn_angles[turn_point_idx[0]])
            cumLength = new_x[turn_point_idx[0]]
            trueidx = np.where(np.abs(new_x-cumLength) < cluster_radius)
            if len(trueidx[0])> 0:
                for i in range(len(trueidx[0])):
                    index = trueidx[0][i]
                    turn_point_idx = np.delete(turn_point_idx, np.where(turn_point_idx == index))

    t_points = []
    
    ##Find the original  point-ID infomation
    for i in range(len(t_pointsID)):
        points = turn_points[t_pointsID[i]]
        t_point_idx = np.where(PointZ == points[1])
        t_points.append(t_point_idx[0])

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
        print SLRatio
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

#---------------------------------------------------------------------------------------------------------------
# This function derives the boundary based on the skyline analysis
#--------------------------------------------------------------------------------------------------------------- 
def BoundaryExtractbySkyline(inDEM, instreams, wspoly, t_points, Elevationshift, outpoly):
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
        poly_layer = arcpy.MakeFeatureLayer_management(thinpolygon, "in_memory\\poly_layer")
        point_layer = arcpy.MakeFeatureLayer_management(t_points, "in_memory\\point_layer")
        arcpy.SelectLayerByLocation_management(poly_layer,"INTERSECT", point_layer,"1 METERS","NEW_SELECTION","")

        selcount = int(arcpy.GetCount_management(poly_layer)[0])
        if selcount > 0:
            arcpy.CopyFeatures_management(poly_layer, outpoly)
        else:
            arcpy.CopyFeatures_management(wspoly, outpoly)


    except arcpy.ExecuteError:
        #arcpy.AddMessage("Skip this part and continue the processing...")
        arcpy.CopyFeatures_management(wspoly, outpoly)

    return outpoly
        
    
##Main program
# Script arguments
InputDEM = arcpy.GetParameterAsText(0)
EleThreshold = arcpy.GetParameter(1)
StreamThresholdKM2 = arcpy.GetParameter(2)
TributaryThresholdKM2 = arcpy.GetParameter(3)
Out3DProfiles = arcpy.GetParameterAsText(4)
OutKnickpoints = arcpy.GetParameterAsText(5) ##the turning points for cirques

cellsize = arcpy.GetRasterProperties_management(InputDEM,"CELLSIZEX")
cellsize_int = int(float(cellsize.getOutput(0)))

StreamThreshold = int(float(StreamThresholdKM2) * 1e6 / (cellsize_int * cellsize_int))
TributaryThreshold = int(float(TributaryThresholdKM2) * 1e6 / (cellsize_int * cellsize_int))
minCirqueAreaThresholdKM2 = 0.5 * TributaryThreshold
minCirqueArea = int(float(minCirqueAreaThresholdKM2) * 1e6 / (cellsize_int * cellsize_int))

###temporay files
MaxFccTable = "in_memory\\MaxFccTable"
StreamOrderTable = "in_memory\\StreamOrderTable"
CleanStream = "in_memory\\CleanStream"
tmpoutStream = "in_memory\\tmpoutStream"
TmpStream = "in_memory\\TmpStream"

###Step 1: Stream network
arcpy.AddMessage("Step 1: Stream extraction...")
## Extract the DEM only higher than the EleThreshold
OutBnd = Con(Raster(InputDEM)> EleThreshold,1)
outExtractDEM = ExtractByMask(InputDEM, OutBnd)

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

##Step 2: clean streams                
arcpy.AddMessage("Step 2: Clean streams...")
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

##Copy the data out for test
#arcpy.CopyFeatures_management(CleanStream, "c:\\test\\cleanstream.shp")

###Step 3: Find the knickpoints for the first order streams only
arcpy.AddMessage("Step 3: Find the knickpoints for the first order streams...")

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

        #arcpy.AddMessage("Start determine the turning points")
        #t_points, t_ratios = turning_points_RDP(LengthfromStart, PointZ, 10, 200)
        #t_points, t_ratios = turning_points_RDE(LengthfromStart, PointZ, 10, 200)
        #t_points, t_ratios = turning_points(LengthfromStart, PointZ) ##, turning_points = 10, cluster_radius = 200)
        t_points, t_ratios = turning_points_ConvexAngle(LengthfromStart, PointZ, 10, 200)

        ##Make the turn point layer
        if len(t_points) > 0:
            CursorPnts = arcpy.da.InsertCursor(Knickpoints, ['SHAPE@', 'LineID', 'Ratio', 'Ele'])
            LineIDlist.append(lineId)
            for i in range(len(t_points)):
                idx = t_points[i] + 1 ##the idx should plus 1 because the t_points are only indexed except for the start and end points
                Ratio = t_ratios[i]
                Pnt = arcpy.Point(PointX[idx], PointY[idx])
                CursorPnts.insertRow((Pnt, lineId, Ratio, PointZ[idx])) ##Add Z as a point attribute here

del cursor, row             

##Copy to the outputs
arcpy.CopyFeatures_management(FirstOrderStream3D, Out3DProfiles)
arcpy.CopyFeatures_management(Knickpoints, OutKnickpoints)

##Delete intermidiate data
arcpy.Delete_management("in_memory") ### Empty the in_memory
