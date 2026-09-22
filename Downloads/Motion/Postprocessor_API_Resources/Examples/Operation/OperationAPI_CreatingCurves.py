# OperationAPI_CreatingCurves.py
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

# Get Active Page
# This retrieves the currently active page in the application.
page = application_handler.GetActivePage()

# Creating a Chart
# Create a new Chart View on the page
# This will create a new chart view with the specified name.
chart_name = "Chart"
chartView = page.CreateChart(chart_name)

# Documentation Example
# Set array about combination of characteristic and component
curvePaths = List[str]()
# characteristic_path_magnitude : characteristic/component path of the curve to add, in "Characteristic/Component" format
characteristic_path_magnitude = r'Displacement/Magnitude'
curvePaths.Add(characteristic_path_magnitude)
# characteristic_path_z : characteristic/component path of the curve to add, in "Characteristic/Component" format
characteristic_path_z = r'Displacement/Z'
curvePaths.Add(characteristic_path_z)

# Create a PlotParameters object to specify the parameters for the plot.
# Set the Entity to Plot.
# The Target is the name of the target for which you want to retrieve the curves.
# Set the paths for the curves you want to retrieve.
# This is where you specify the characteristics and components you want to plot.
parameters = PlotParameters()
parameters.Paths = curvePaths
parameters.Target = "Crank"

# Add Curves (FilePath, Curve Parameter)
# FilePath - The path of the result to access.
# parameters - The class used as a parameter of the AddCurve function.
# The instance of the curve.
curves = chartView.AddCurves(RESULT_FILE_PATH, parameters)
count = len(curves)
if count > 0 :
    curves[0].SeriesName = r'Crank_DM'
    curves[0].Color = Colors.Blue
    curves[0].LineStyle = LineStyleType.Dot

# Close the Pages
page.Close()

# Close the Document
application_handler.CloseDocument(RESULT_FILE_PATH)
