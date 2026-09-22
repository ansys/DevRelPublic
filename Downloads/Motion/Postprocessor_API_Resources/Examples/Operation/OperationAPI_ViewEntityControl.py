# -*- coding: utf-8 -*-
# OperationAPI_ViewEntityControl.py
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

# Get Active View
# This retrieves the currently active view in the page.
animationView = page.GetActiveView()

# Get instance of entity
# The GetViewModelByName method retrieves the target by its name.
febody = animationView.GetViewModelByName("FEBody_01")

# Hides all entities except the specified ones. 
# In this case, we hide all entities except the FEBody_01.
# Input: febody.FullName 
# or Input: febody.ID
animationView.HideOthers(febody.ID)

# Fit the view to the current selection.
animationView.Fit()

# Show all entities in the view.
animationView.ShowAll()

# Get all created pages.
# This retrieves all pages created in the application.
pages = application_handler.GetPages()
for page in pages :
    page.Close()

# Close Document
application_handler.CloseDocument(RESULT_FILE_PATH)
