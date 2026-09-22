# OperationAPI_ImageCapture.py
import sys

# OPERATION_API_MODULE_PATH : path to the folder containing OperationAPI.py
sys.path.append(OPERATION_API_MODULE_PATH)

# Import necessary modules
from OperationAPI import *
import os

# Start the headless application interface
application_handler = ApplicationHandler()

# Set array about result file
filepaths = List[str]()
# RESULT_FILE_PATH : .dfr result file path
filepaths.Add(RESULT_FILE_PATH)

# Open about result files
# This will open the result file in the application.
# When the result is first opened, a Page is created and an Animation View is created on that Page.
application_handler.AddDocument(filepaths)

# Get all created pages.
# This retrieves all pages created in the application.
pages = application_handler.GetPages()

findViews = list()
for page in pages :
    views = page.GetViews()
    animationViews = [view for view in views if view.ViewType == ViewType.Animation]
    for animationView in animationViews :
        if animationView.DocumentFilePath == RESULT_FILE_PATH and animationView.AnalysisResultType == AnalysisResultType.Dynamics :
            findViews.append(animationView)

viewCount = len(findViews)
if viewCount > 0 :
    # RESULT_FILE_PATH - Get the document from the result file path.
    document = application_handler.GetDocument(RESULT_FILE_PATH)

    # This retrieves the analysis result from the document.
    # Types of Analysis Results
    # - Dynamics
    # - Eigenvalue
    analysis = document.GetAnalysisResultViewModel(AnalysisResultType.Dynamics)

    animationview = findViews[0]
    frame_count = 10
    analysis.SetAnimationFrame(frame_count)
    frame_index = 5
    analysis.MoveToAnimationFrame(frame_index) # MoveToAnimationFrame

    # OUTPUT_DIR : path to the folder where exported files are written
    export_filepath = os.path.join(OUTPUT_DIR, r'Image.png')

    # Image Caputre
    animationview.ExportImage(export_filepath, ImageFormat.Png, 1920, 1080)

# Get all created pages.
# This retrieves all pages created in the application.
pages = application_handler.GetPages()
for page in pages :
    page.Close()

# Close the Document
application_handler.CloseDocument(RESULT_FILE_PATH)
