# OperationAPI_CreatingAExpression.py
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

viewCount = len(findViews)
if viewCount > 0 :
    # RESULT_FILE_PATH - Get the document from the result file path.
    document = application_handler.GetDocument(RESULT_FILE_PATH)

    # This retrieves the analysis result from the document.
    # Types of Analysis Results
    # - Dynamics
    # - Eigenvalue
    analysis = document.GetAnalysisResultViewModel(AnalysisResultType.Dynamics)

    # Creating a Expression
    expression = analysis.CreateExpression("expression")
    expression.Expression = 'DM("Crank/CM")'

    # Get the expression and add its value as curve data.
    expression = analysis.GetExpression(expression.FullName)
    chart_name = "ExpressionChart"
    chart = page.CreateChart(chart_name)
    curve_paths = List[str]()
    # characteristic_path : characteristic/component path of the curve to add, in "Characteristic/Component" format
    characteristic_path = "Expression/Value"
    curve_paths.Add(characteristic_path)
    parameters = PlotParameters()
    parameters.Target = expression.FullName
    parameters.Paths = curve_paths
    curves = chart.AddCurves(RESULT_FILE_PATH, parameters)

    # Remove the expression from the dynamic analysis result.
    analysis.RemoveExpression(expression.FullName)

# Close the Page
page.Close()

# Close the Document
application_handler.CloseDocument(RESULT_FILE_PATH)
