# OperationAPI_CreatingAnimationViewBasedonType.py
import sys

# OPERATION_API_MODULE_PATH : path to the folder containing OperationAPI.py
sys.path.append(OPERATION_API_MODULE_PATH)

# Import necessary modules
from OperationAPI import *

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


# RESULT_FILE_PATH - Get the document from the result file path.
document = application_handler.GetDocument(RESULT_FILE_PATH)

# This retrieves the analysis result from the document.
    # Types of Analysis Results
    # - Dynamics
    # - Eigenvalue
dynamic_analysis = document.GetAnalysisResultViewModel(AnalysisResultType.Dynamics)
frame_count = 50
dynamic_analysis.SetAnimationFrame(frame_count)

viewCount = len(findViews)
if viewCount > 0 :
    page1 = application_handler.GetPage(findViews[0].GroupID)

    # This retrieves the analysis result from the document.
    # Types of Analysis Results
    # - Dynamics
    # - Eigenvalue
    eigenval_analysis = document.GetAnalysisResultViewModel(AnalysisResultType.Eigenvalue)

    # Create an Animation View on the active page
    # This will create an animation view based on the eigenvalue analysis.
    animation_view_name = "EigenvalueAnimation"
    eigenvalue_animation = page1.CreateAnimation(eigenval_analysis, animation_view_name)
    eigenval_analysis.Frame = 100
    
    # Get Sampling Times
    # This retrieves the sampling times from the eigenvalue analysis.
    times = eigenval_analysis.GetSamplingTimes()
    convertedtimes = list(times)

    # Set the target sampling time.
    eigenval_analysis.TargetSamplingTime = convertedtimes[0]

    # Get the frequency instance for the specified sampling time index and enable its mode-shape animation.
    frequency = eigenval_analysis.GetFrequency(0)
    frequency.Enable = True

# Close the Pages
# Get all created pages.
# This retrieves all pages created in the application.
pages = application_handler.GetPages()
for page in pages :
    page.Close()

# Close the Document
application_handler.CloseDocument(RESULT_FILE_PATH)
