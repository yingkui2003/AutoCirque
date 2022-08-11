# Introduction
The AutoCirque toolbox includes three GIS tools: 1) Cirque Auto-Delineation; 2) Cirque Delineation by Threshold Features; and 3) Cirque Potential Thresholds. All these tools are written in python and imported into ArcGIS as python script tools. These tools have been tested successfully for different ArcGIS versions after 10.5 and ArcGIS Pro 2.8. All python scripts have been imported to their corresponding tools; therefore, the user can directly download the *.tbx file and run it in ArcGIS or ArcGIS pro. The python codes are also included in the github site. The user can check the codes and comtinue imporving this toolbox. 

Cirque Auto-Delineation is a fully automated tool to delineate cirque outlines. The input parameters include DEM, a minimum elevation for cirques, the minimum source area to start a stream, the minimum total contributing area of a stream before joining another stream, the minimum circularity index, the minimum percentages of steep (>31°) and gentle (<20°) slope areas, and the minimum headwall slope, as well as its allowable flexibility, for the cirques. The DEM is recommended to have a UTM projection to ensure the correct calculations of DEM-related metrics without the specification of a Z-factor. The output parameters include the 3D profiles for each first order stream related to the potential cirques, the potential threshold points and cross sections along the profiles, and the delineated potential cirque outlines. 

![image](https://user-images.githubusercontent.com/24683137/140794960-6bcc276d-616e-414b-841c-d5d354537c7d.png)

Cirque Delineation By Threshold Features is the semi-automated tool to delineate the cirques based on DEM and user specified cirque threshold features (points or cross sections). The input parameters of this tool include DEM, the threshold points or cross sections, the minimum source area to start a stream, and the minimum headwall slope, as well as its allowable flexibility, for the cirques. The output of this tool is the delineated cirque outlines.

![image](https://user-images.githubusercontent.com/24683137/140795153-174b2e8f-c374-4877-88a0-4e08a5b3034e.png)

The cirque threshold points can be prepared by digitizing the points in ArcGIS or Google Earth (then transfer the KML file to ArcGIS). In terms of the cross sections, the user can digitize the linear features in ArcGIS or Google Earth to define the low margins of the cirques. The user can also use the interpolate 3D lines in the 3D Analyst Toolbar in ArcGIS to draw the potential cross sections and then convert the drawn graphic objects to feature classes using the Convert Graphics To Features tool in ArcGIS. 

A Cirque Potential Thresholds tool is also developed to identify potential cirque threshold points and cross sections along the stream profile. This tool generates a set of potential threshold points and cross sections along the first order stream, so that the user can select the appropriate threshold features to delineate cirques using the Cirque Delineation by Threshold Features tool. The input parameters of this tool include DEM, a minimum elevation for cirques, the minimum source area to start a stream, and the minimum total contributing area of a stream before joining another stream. The output parameters include the 3D profiles, the potential threshold points, and the potential threshold cross sections along the profiles.

![image](https://user-images.githubusercontent.com/24683137/140795019-bc30437d-7b92-4070-beaa-dc1b8103d471.png)

# How to download and use this toolbox in ArcGIS or ArcGIS Pro
The github site includes a AutoCirque.tbx file and a python folder, including three python files associated with the three tools. The user can click "Code" (green color) on the right side of the github page and choose Download Zip.

![image](https://user-images.githubusercontent.com/24683137/140193749-35671a0e-2664-4bca-b487-a1e2f071fc36.png)

A zip file of the while github folder will be downloaded to the local computer. Unzip this file will create a AutoCirque-main folder with both the tbx file and the python folder and three code files. The user can use this toolbox, check the codes, and comtinue imporving this toolbox.

The 'Cirque Auto-Delineation' tool may have some errors (unable to allocate memory) when processing a reltively large DEM (when processing >600 potential cirque thresholds). If this happens, try to run this tool in ArcGIS Pro or you can run the 'Cirque Potential Thresholds' tool to generate the threshold points/cross sections, select the suitable threshold points/cross sections, and use the 'Cirque Delineation By Threshold Features' tool to delineate the cirque outlines. 

Please report any errors or questions to Yingkui Li (yli32@utk.edu).

# Cite this work
Li Y. and Zhao, Z., 2021. AutoCirque: An automated method to delineate glacial cirque outlines from digital elevation models. Geomorphology, https://doi.org/10.1016/j.geomorph.2021.108059.

# Most recent updates
08/11/2022: Revise the codes to speed up the cirque delineation by turning off the plan_clos calculation. Plan_clos can be derived using ACME tool later after determining the cirque outlines.


# Contact info
Yingkui Li

Department of Geography & Sustainability

University of Tennessee

Knoxville, TN 37996

Email: yli32@utk.edu

Website: https://geography.utk.edu/about-us/faculty/dr-yingkui-li/

Google Scholar: https://scholar.google.com/citations?user=JoNuyCMAAAAJ&hl=en&oi=ao

