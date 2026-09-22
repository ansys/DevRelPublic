# -*- coding: utf-8 -*-
# OperationAPI_CreateACoordinateSystem.py
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
    animationview = findViews[0]
    # Create a Coordinate System
    # Name - Set the name of instance.
    # ParentInfo - Specifies The path of an parent entity.
    csys_name = "Crank_CSYS"
    parent_name = "Crank"
    animationview.CreateCoordinateSystem(csys_name, parent_name)

    # Get instance of entity
    # The GetViewModelByName method retrieves the target by its name.
    csys = animationview.GetViewModelByName("Crank_CSYS")

    # Set the transformation offset position and angle for the coordinate system.
    csys.TransformationOffsetPosition = Vector(0,10,0)
    csys.TransformationOffsetAngle = Vector(1,0,0)
    csys.TransformationOffsetRotationType = RotationTypes.FixedAngle
    csys.TransformationOffsetRotationAxis = RotationAxes.XYX
    csys.MarkerSize = 10

    # The CurrentCoordinateSystemType property indicates Type attribute as displayed in the user interface.
    # Accepted input values are defined by the GeneralMarkerType enumeration.
    # The GeneralMarkerType ? property must be one of the following: CARTESIAN, CYLINDRICAL, SPHERICAL
    csys.CurrentCoordinateSystemType = GeneralMarkerType.CYLINDRICAL
    if csys.CurrentCoordinateSystemType == GeneralMarkerType.SPHERICAL:
        # The SphericalAxis1 property indicates Axis ?. Accepted input values are defined by the CoordinateType enumeration.
        # The Axis ? property must be one of the following: X, Y, Z
        csys.SphericalAxis1 = CoordinateType.X
        # The SphericalAxis2 property indicates Axis ?? Accepted input values are defined by the CoordinateType enumeration.
        # The Axis ??property must be one of the following: X, Y, Z
        csys.SphericalAxis2 = CoordinateType.Y
    elif csys.CurrentCoordinateSystemType == GeneralMarkerType.CYLINDRICAL:
        # The CylindricalAxisR property indicates Axis R. Accepted input values are defined by the CoordinateType enumeration.
        # The Axis R property must be one of the following: X, Y, Z
        csys.CylindricalAxisR = CoordinateType.Z
        # The CylindricalAxisZ property indicates Axis Z. Accepted input values are defined by the CoordinateType enumeration.
        # The Axis Z property must be one of the following: X, Y, Z    
        csys.CylindricalAxisZ = CoordinateType.X

# Close the Pages
# Get all created pages.
# This retrieves all pages created in the application.
pages = application_handler.GetPages()
for page in pages :
    page.Close()

# Close the Document
application_handler.CloseDocument(RESULT_FILE_PATH)
