#-------------------------------------------------------------------------------
# Tool Name: CirqueDelineationUsingAssignedFeatures.py
# Purpose: This tool delineates cirques based on input outlet features. The outlet
# feature can be points (outlet points) or polyline features, such as cross sections
# the cirque outlet or curved boundaries (or moraines) for the cirque.
# 
# Author: Dr. Yingkui Li
# Created:     09/21-12/29/2020
# modified:    09/07-10/01/2021
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
 
    return      

       
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
InputDEM = arcpy.GetParameterAsText(0)
InputFC  = arcpy.GetParameterAsText(1) ##Input turning points or cross sections around the outlet points
StreamThresholdKM2 = arcpy.GetParameter(2)
Max_backwall_slope =  arcpy.GetParameter(3)
Flexible_backwall_slope =  arcpy.GetParameter(4)

OutCirques = arcpy.GetParameterAsText(5)

Min_backwall_slope = int(Max_backwall_slope) - int(Flexible_backwall_slope)

FcType = arcpy.Describe(InputFC).shapeType
if FcType == "Point" or FcType == "Multipoint":
    bFCLine = 0
elif FcType == "Polyline":
    bFCLine = 1
else:
    arcpy.AddError("The Input Feature should be points and line features!")
    arcpy.GetMessage(0)


cellsize = arcpy.GetRasterProperties_management(InputDEM,"CELLSIZEX")
cellsize_int = int(float(cellsize.getOutput(0)))

#StreamThreshold = int(float(StreamThresholdKM2) * 1e6 / (cellsize_int * cellsize_int))
StreamThreshold = 10

###temporay files
smoothtolerance = cellsize_int * 5
FCselected = "in_memory\\FCselected"  ##Set a in_memory file for each moraine feature
Singlepoint = "in_memory\\Singlepoint"
Streams = "in_memory\\Streams"
Singlestream = "in_memory\\Singlestream"
SingleWs = "in_memory\\SingleWs"
singleBND = "in_memory\\singleBND"

###Get the count of cirque features
countResult = arcpy.GetCount_management(InputFC)
count = int(countResult.getOutput(0))

if count < 1:
    arcpy.AddMessage("There is no features to identify cirques")
    sys.exit()

###Step 1: Stream network
arcpy.AddMessage("Flow direction and accumulation analysis...")

#Hydro analysis
outslope = Slope(InputDEM)
fillDEM =Fill(InputDEM)  ##Fill the sink first
fdir = FlowDirection(fillDEM,"NORMAL") ##Flow direction

facc = FlowAccumulation(fdir) ##Flow accmulation

outGreaterThan = Con(facc > StreamThreshold, 1,0)  ##Determine the highest flowaccumuation part

# Process: Stream Link
outStreamLink = StreamLink(outGreaterThan, fdir)
    
# Process: Stream to Feature
StreamToFeature(outStreamLink, fdir, Streams, "NO_SIMPLIFY")
#arcpy.CopyFeatures_management(Streams, "c:\\testdata\\fstreams3.shp")


inputFCcopy = "in_memory\\inputFCcopy"
inputFCcs = "in_memory\\inputFCcs"
arcpy.CopyFeatures_management(InputFC, inputFCcopy)
FcID = arcpy.Describe(inputFCcopy).OIDFieldName


if bFCLine == 0: ##if the inputFc is point feature, need to create the cross sections for each point
    MaxFccTable = "in_memory\\MaxFccTable"
    TmpStream = "in_memory\\TmpStream"

    #arcpy.CopyFeatures_management(Singlestream, "c:\\testdata\\fstreams2.shp")
    StreamToFeature(outStreamLink, fdir, TmpStream, "SIMPLIFY")
    ZonalStatisticsAsTable(outStreamLink, "VALUE", facc, MaxFccTable, "DATA", "MAXIMUM")
    # Process: Join Field
    arcpy.JoinField_management(TmpStream, "grid_code", MaxFccTable, "Value", "MAX")  ##Join to get the flow accumulation value

    ##test the cross section function here!!!!
    crosssections = cross_sections(inputFCcopy, FcID, TmpStream, "MAX", 50, 120)
    #arcpy.CopyFeatures_management(crosssections, "c:\\testdata\\crosssection1.shp")
    arcpy.CopyFeatures_management(crosssections, inputFCcs)
    try:
        arcpy.Delete_management (MaxFccTable)
        arcpy.Delete_management (TmpStream)
    except:
        pass
else:
    arcpy.CopyFeatures_management(inputFCcopy, inputFCcs)
    
TotalBND = arcpy.CreateFeatureclass_management("in_memory", "TotalBND","POLYGON","","","",InputFC)
countResult = arcpy.GetCount_management(inputFCcs)
count = int(countResult.getOutput(0))
FcID = arcpy.Describe(inputFCcs).OIDFieldName
#FcID = arcpy.Describe(inputFCcopy).OIDFieldName

for ifeature in range (count):
    arcpy.AddMessage("Generating cirque "+str(ifeature + 1)+" of "+str(count))
    query = FcID +" = "+str(ifeature+1)
    arcpy.Select_analysis(inputFCcs, FCselected, query)
    
    #arcpy.CopyFeatures_management(Streams, "c:\\testdata\\Streams.shp")
    #arcpy.CopyFeatures_management(FCselected, "c:\\testdata\\FCselected.shp")

    arcpy.Intersect_analysis([Streams, FCselected], Singlepoint, "#", "#", "POINT")

    pntcountResult = arcpy.GetCount_management(Singlepoint)
    pntcount = int(pntcountResult.getOutput(0))
    #arcpy.AddMessage("the number of point is:" + str(pntcount))
    if (pntcount == 0): ##if no intersect points, use the buffer to get the intersection points for another time
        tmpbuf = "in_memory\\tmpbuf"
        arcpy.Buffer_analysis(FCselected, tmpbuf, "5 Meters")
        arcpy.Intersect_analysis([Streams, tmpbuf], Singlepoint, "#", "#", "POINT")
        pntcountResult = arcpy.GetCount_management(Singlepoint)
        pntcount = int(pntcountResult.getOutput(0))
        #arcpy.AddMessage("the number of point is:" + str(pntcount))
        arcpy.Delete_management (tmpbuf)
        
    if (pntcount > 0):
        #Calculate Watershed
        outPour = SnapPourPoint(Singlepoint, facc, 0)
        outWs1 = Watershed(fdir, outPour)
        outWs1.save("c:\\testdata\\outWs1.tif")
        outWs = Con(outWs1 >= 0, 1)  ##Determine the highest flowaccumuation part and set others as Nodata
        arcpy.RasterToPolygon_conversion(outWs, SingleWs, "NO_SIMPLIFY", "VALUE")
        
        arcpy.Clip_analysis(Streams, SingleWs, Singlestream)
        
        singleDEM = ExtractByMask(fillDEM, outWs)
        ext_slope = ExtractByMask(outslope, outWs)

        singleBND, status = BoundaryRefineBySlope(singleDEM, ext_slope, Singlepoint, int(Max_backwall_slope), int(Min_backwall_slope), 120)

        polycountResult = arcpy.GetCount_management(singleBND)
        polycount = int(polycountResult.getOutput(0))

        if (polycount > 0 and status): ##meet the cirque standards
            arcpy.Append_management(singleBND, TotalBND, "NO_TEST")        

##Smooth polygons
#arcpy.CopyFeatures_(TotalBND, OutCirques, "PAEK", smoothtolerance)
#arcpy.CopyFeatures_management(TotalBND, OutCirques)
polycountResult = arcpy.GetCount_management(TotalBND)
polycount = int(polycountResult.getOutput(0))
if polycount > 0:
    arcpy.cartography.SmoothPolygon(TotalBND, OutCirques, "PAEK", smoothtolerance)
    arcpy.AddField_management(OutCirques, "CirIndex", "DOUBLE")
    with arcpy.da.UpdateCursor(OutCirques, ['SHAPE@LENGTH', 'SHAPE@AREA', 'CirIndex']) as cursor:
        for row in cursor:
            area = float(row[1])
            #if area > minCirqueArea: ##Do not use the min area for this tool
            length = float(row[0])
            row[2] = 4.0 * 3.14 * area / (length * length)
            cursor.updateRow(row)
            #else:
            #    cursor.deleteRow()
    del cursor, row
else:
    arcpy.AddMessage("Fail to delineate the polygons!")
##Delete intermidiate data
arcpy.Delete_management("in_memory") ### Empty the in_memory
'''
##It seems that arcpy.Delete_management ("in_memory") already remove all in_memory datasets
##Delete the created in-memory dataset
arcpy.Delete_management (FCselected)
arcpy.Delete_management (Singlepoint)
arcpy.Delete_management (Streams)
arcpy.Delete_management (Singlestream)
arcpy.Delete_management (SingleWs)
arcpy.Delete_management (singleBND)
arcpy.Delete_management (inputFCcopy)
arcpy.Delete_management (inputFCcs)
#Delete Raster datasets
arcpy.Delete_management (outslope)
arcpy.Delete_management (fillDEM)
arcpy.Delete_management (fdir)
arcpy.Delete_management (facc)
arcpy.Delete_management (outGreaterThan)
arcpy.Delete_management (outStreamLink)
'''
