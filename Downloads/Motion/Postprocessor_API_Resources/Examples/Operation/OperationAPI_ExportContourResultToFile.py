# OperationAPI_ExportContourResultToFile.py
# -*- coding: utf-8 -*-
import sys

# OPERATION_API_MODULE_PATH : path to the folder containing OperationAPI.py
sys.path.append(OPERATION_API_MODULE_PATH)

# Import necessary modules
from OperationAPI import *
import os

# Start the headless application interface
application_handler = ApplicationHandler()

# Set array about result file - There may be multiple result files.
filepaths = List[str]()
# RESULT_FILE_PATH : .dfr result file path
filepaths.Add(RESULT_FILE_PATH)

# Open about result files
# When the result is first opened, a Page is created and an Animation View is created on that Page.
application_handler.AddDocument(filepaths)

# Get Active Page
# This retrieves the currently active page in the application.
page = application_handler.GetActivePage()

# Get Active View
# This retrieves the currently active view in the page.
activeView = page.GetActiveView()

# Loading Animation - Use this only when the current Active View is an Animation View.
# This method is used to load the animation in the active view.
application_handler.LoadingAnimation()

# Get instance of entity
# The GetViewModelByName method retrieves the target by its name.
febody = activeView.GetViewModelByName("FEBody_01")

# Check if the FEBody is found
if febody is None:
    sys.exit(1)

# OUTPUT_DIR : path to the folder where exported files are written
output_path = os.path.join(OUTPUT_DIR, "Displacement.txt")

# A boolean indicating whether to export the full frame or not.
# If True, the full frame will be exported; if False, only the current frame will be exported.
is_full_frame = True

# Set the febody fullname as the target for export
targetNames = List[str]()
targetNames.Add(febody.FullName)

# Set Contour Mapping Type
# ContourMappingType can be one of the following:
# - ContourMappingType.None - None
# - ContourMappingType.FENode - Node (Averaged across body)
# - ContourMappingType.FEElement - Element (Unaveraged)
# - ContourMappingType.FEElementNode - Node (Unaveraged)
# - ContourMappingType.FEMaterial - Node (Averaged within material)
# - ContourMappingType.BeamGroup - Beam Group
# - ContourMappingType.Contact - Contact
# - ContourMappingType.ChainedSystem - Chained System
# - ContourMappingType.Usersubroutine - Usersubroutine
contour_mapping_type = ContourMappingType.FENode

# Set Characteristic Path
characteristic_path = "Displacement"

# Set Component Path
component_path = "X"

# Combine Characteristic and Component Path
combined_characteristic_component_path = characteristic_path + "/" + component_path

# Export Contour Result to File - When the Animation View is the Active View
# RESULT_FILE_PATH - The path to the result file.
# output_path - The path to the output file where the contour results will be saved.
# targetNames - The list of target names for which the contour results will be exported.
# is_full_frame - A boolean indicating whether to export the full frame or not.
# contour_mapping_type - The type of contour mapping to be used.
# characteristic_component_path - The path to the characteristic and component for which the contour results will be exported.
# AnalysisType - The default value is AnalysisResultType.Dynamics. If you want to change it, please set the AnalysisType.
# FileFormatType - The format of the output file. The default is FileFormatType.Text. The format type can be set to TXT or BINARY.
Export.ExportContourResultToFile(RESULT_FILE_PATH, output_path, targetNames, is_full_frame, contour_mapping_type, combined_characteristic_component_path)

# Close the Pages
# Get all created pages.
pages = application_handler.GetPages()
for page in pages :
    page.Close()

# Close Document
application_handler.CloseDocument(RESULT_FILE_PATH)
