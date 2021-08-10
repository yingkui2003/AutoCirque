#-------------------------------------------------------------------------------
# Tool Name: CirqueDelineationUsingAssignedFeatures.py
# Purpose: This tool delineates cirques based on input outlet features. The outlet
# feature can be points (outlet points) or polyline features, such as cross sections
# the cirque outlet or curved boundaries (or moraines) for the cirque.
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
        poly_layer = arcpy.MakeFeatureLayer_management(thinpolygon, "in_memory\\poly_layer")
        stream_layer = arcpy.MakeFeatureLayer_management(instreams, "in_memory\\stream_layer")
        arcpy.SelectLayerByLocation_management(poly_layer,"INTERSECT", stream_layer,"1 METERS","NEW_SELECTION","")

        selcount = int(arcpy.GetCount_management(poly_layer)[0])
        if selcount > 0:
            arcpy.Dissolve_management(poly_layer, outpoly)
        else:
            arcpy.CopyFeatures_management(wspoly, outpoly)

    except arcpy.ExecuteError:
        arcpy.CopyFeatures_management(wspoly, outpoly)

    return outpoly

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
        b = zip(*arr)
        points = np.array(list(zip(b[0],b[1])))

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
            print turn_count
        ##Make the new feature class
        numpy_array_to_features(new_line, arr, ['SHAPE@X', 'SHAPE@Y'], field)

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

def BoundaryRefineBySlope(dem, slopeRaster, stream, max_slope, min_slope, max_angle):
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
    #arcpy.CopyFeatures_management(polyline, "c:\\test\\merge_polyline.shp")

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
    
    return newpoly, True

    #except:
    #    return bndcleanpoly, False
        
    
##Main program
# Script arguments
InputDEM = arcpy.GetParameterAsText(0)
InputFC  = arcpy.GetParameterAsText(1) ##Input turning points or cross sections around the outlet points
StreamThresholdKM2 = arcpy.GetParameter(2)
Max_backwall_slope =  arcpy.GetParameter(3)
Min_backwall_slope =  arcpy.GetParameter(4)

OutCirques = arcpy.GetParameterAsText(5)

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

StreamThreshold = int(float(StreamThresholdKM2) * 1e6 / (cellsize_int * cellsize_int))
#minCirqueArea = StreamThresholdKM2 * 1e6

###temporay files
smoothtolerance = cellsize_int * 10
FCselected = "in_memory\\FCselected"  ##Set a in_memory file for each moraine feature
Singlepoint = "in_memory\\Singlepoint"
Singlestream = "in_memory\\Singlestream"
SingleWs = "in_memory\\SingleWs"
singleBND = "in_memory\\singleBND"
#TotalBND = "in_memory\\TotalBND"
tmpbuf = "in_memory\\tmpbuf"

###Get the count of cirque features
countResult = arcpy.GetCount_management(InputFC)
count = int(countResult.getOutput(0))

if count < 1:
    arcpy.AddMessage("There is no features to identify cirques")
    sys.exit()

FcID = arcpy.Describe(InputFC).OIDFieldName

###Step 1: Stream network
arcpy.AddMessage("Flow direction and accumulation analysis...")

#Hydro analysis
outslope = Slope(InputDEM)
fillDEM =Fill(InputDEM)  ##Fill the sink first
fdir = FlowDirection(fillDEM,"NORMAL") ##Flow direction
facc = FlowAccumulation(fdir) ##Flow accmulation

TotalBND = arcpy.CreateFeatureclass_management("in_memory", "TotalBND","POLYGON","","","",InputFC)

for ifeature in range (count):
    arcpy.AddMessage("Generating cirque "+str(ifeature + 1)+" of "+str(count))
    if FcID == "FID":
        query = FcID +" = "+str(ifeature)
    else:
        query = FcID +" = "+str(ifeature+1)
    arcpy.Select_analysis(InputFC, FCselected, query)

    if bFCLine > 0: ##if the inputFc is polyline
        ## use the small buffer of the line to make sure it can always get the maximum facc
        ##make a small buffer of the cross section to make sure the cross section get the highest fcc
        arcpy.Buffer_analysis(FCselected, tmpbuf, (str(cellsize_int/2)+ " Meter"))
        bufID = arcpy.Describe(tmpbuf).OIDFieldName
        
        outZonalStatistics = ZonalStatistics(tmpbuf, bufID, facc, "MAXIMUM") #Find the maximum flowaccumulation point on the moriane feature
        outPour = Con(facc == outZonalStatistics,facc)  ##Determine the highest flowaccumuation part
        
        #outZonalStatistics = ZonalStatistics(FCselected, FcID, facc, "MAXIMUM") #Find the maximum flowaccumulation point on the moriane feature
        #outPour = Con(facc == outZonalStatistics,facc)  ##Determine the highest flowaccumuation part
        arcpy.RasterToPoint_conversion(outPour, Singlepoint, "VALUE")
    else:
        outPour = SnapPourPoint(FCselected, facc, cellsize_int) ## Just create a pourpoint raster with the same extent of the input DEM
        arcpy.RasterToPoint_conversion(outPour, Singlepoint, "VALUE")          

    #Calculate Watershed
    outPour = SnapPourPoint(Singlepoint, facc, cellsize_int)
    outWs = Watershed(fdir, outPour)
    #outWs = Watershed(fdir, Singlepoint)
    
    arcpy.RasterToPolygon_conversion(outWs, SingleWs, "NO_SIMPLIFY", "VALUE")
    #arcpy.CopyFeatures_management(SingleWs, "c:\\test\\testsinglews.shp")
    
    # Process: Extract by Mask
    ExtraFcc = ExtractByMask(facc,outWs)
    # Process: Greater Than
    outGreaterThan = Con(ExtraFcc > StreamThreshold, 1,0)  ##Determine the highest flowaccumuation part

    # Process: Stream Link
    outStreamLink = StreamLink(outGreaterThan, fdir)
    
    # Process: Stream to Feature
    StreamToFeature(outStreamLink, fdir, Singlestream, "SIMPLIFY")

    #arcpy.CopyFeatures_management(Singlestream, "c:\\test\\Singlestream.shp")
    
    singleDEM = ExtractByMask(fillDEM, outWs)
    ext_slope = ExtractByMask(outslope, outWs)

    #singleBND = SkylineBoundaryExtract(singleDEM, Singlestream, SingleWs, 10, singleBND)
    singleBND, status = BoundaryRefineBySlope(singleDEM, ext_slope, Singlestream, int(Max_backwall_slope), int(Min_backwall_slope), 120)

    polycountResult = arcpy.GetCount_management(singleBND)
    polycount = int(polycountResult.getOutput(0))

    #arcpy.AddMessage("polgon count: " + str(polycount))
    #arcpy.AddMessage("status: " + str(status))
    
    if (polycount > 0 and status): ##meet the cirque standards
        #if ifeature < 1: ##The first loop
        #    arcpy.CopyFeatures_management(singleBND, TotalBND)
        #else:
        arcpy.Append_management(singleBND, TotalBND, "NO_TEST")        

##Smooth polygons
#arcpy.CopyFeatures_(TotalBND, OutCirques, "PAEK", smoothtolerance)
arcpy.CopyFeatures_management(TotalBND, OutCirques)

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

##Delete intermidiate data
arcpy.Delete_management("in_memory") ### Empty the in_memory
