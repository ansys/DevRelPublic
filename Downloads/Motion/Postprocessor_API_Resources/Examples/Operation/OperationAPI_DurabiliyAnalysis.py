# OperationAPI_DurabiliyAnalysis.py
import sys

# OPERATION_API_MODULE_PATH : path to the folder containing OperationAPI.py
sys.path.append(OPERATION_API_MODULE_PATH)

# Import necessary modules
from OperationAPI import *

# Start the headless application interface
application_handler = ApplicationHandler()

# Prepare result file list
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
active_view = page.GetActiveView()

# Load animation
application_handler.LoadingAnimation()

# Get instance of entity
# The GetViewModelByName method retrieves the target by its name.
fe_body = active_view.GetViewModelByName("FEBody_01")
material = active_view.GetViewModelByName("2014_HV_O")

# Set FE body properties
fe_body.Material = material
fe_body.AnalysisType = FatigueAnalysisType.SN
fe_body.Scale = 2.0
fe_body.StressStrainCombination = StressStrainCombinationType.VonMises
fe_body.MeanStressCorrection = MeanStressCorrection.Neglect

# Set up fatigue analysis parameters
durability_param = DurabilityAnalysisParameter()
durability_param.DocumentFilePath = RESULT_FILE_PATH
# FATIGUE_RESULT_NAME : name of the fatigue result file in the Fatigue folder, without the .dffrf extension
durability_param.ResultName = FATIGUE_RESULT_NAME
durability_param.NoOfRepeatedLoad = 1
durability_param.Start = 1
durability_param.End = 10
durability_param.AddTarget(fe_body.ID)

# Run fatigue analysis
DurabilityAnalysis.RunFatigueAnalysis(durability_param)

# Create fatigue contours
# fatigue_characteristic_name : fixed characteristic name required by the fatigue contour API
fatigue_characteristic_name = "Fatigue"
# fatigue_life_cycle_component_name : fixed component name required for the fatigue life-cycle contour
fatigue_life_cycle_component_name = "Life Cycle"
# fatigue_damage_component_name : fixed component name required for the fatigue damage contour
fatigue_damage_component_name = "Damage"
document = application_handler.GetDocument(RESULT_FILE_PATH)
analysis_result = document.GetAnalysisResultViewModel(AnalysisResultType.Dynamics)
target_entities = List[str]()
target_entities.Add(fe_body.FullName)
analysis_result.CreateContour(target_entities, ContourMappingType.FEElement, fatigue_characteristic_name, fatigue_life_cycle_component_name)
analysis_result.CreateContour(target_entities, ContourMappingType.FEElement, fatigue_characteristic_name, fatigue_damage_component_name)

# Close all pages and document
# Get all created pages.
# This retrieves all pages created in the application.
for p in application_handler.GetPages():
    p.Close()
application_handler.CloseDocument(RESULT_FILE_PATH)
