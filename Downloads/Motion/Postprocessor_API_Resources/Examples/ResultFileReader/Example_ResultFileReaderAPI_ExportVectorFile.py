import sys

# BINARY_FOLDER_PATH : path to the folder containing the Postprocessor dll files
# RESULT_FILE_READER_API_MODULE_PATH : path to the folder containing ResultFileReaderAPI.py
sys.path.append(BINARY_FOLDER_PATH)
sys.path.append(RESULT_FILE_READER_API_MODULE_PATH)

# Import necessary modules
from ResultFileReaderAPI import *
import os

# RESULT_FILE_PATH : .dfr result file path
output_reader = OutputReader(RESULT_FILE_PATH)

state_ids = output_reader.GetStateIDArray()

# target - Specifies the name of vector displayable entity
# path - Specifies characteristc on vector display
target_name = "TJ_01"
# vector_path : fixed vector characteristic path for the Base Force result
vector_path = "Base Force"

targets = List[IVectorDisplay]()
vector_name = "vector"
vector = output_reader.CreateVector(vector_name, target_name, vector_path)
targets.Add(vector)
print ("===ExportVectorDisplayToFile===")

# OUTPUT_DIR : path to the folder where exported files are written
export_vector_file_path = os.path.join(OUTPUT_DIR, RFR_VECTOR_TARGETS_OUTPUT_FILENAME)
output_reader.ExportVectorDisplayToFile(export_vector_file_path, state_ids, targets, True, True, True, AnalysisResultType.Dynamics)

export_vector_file_path = os.path.join(OUTPUT_DIR, RFR_VECTOR_TARGET_OUTPUT_FILENAME)
output_reader.ExportVectorDisplayToFile(export_vector_file_path, state_ids, target_name, vector_path, True, True, True, AnalysisResultType.Dynamics)

# Close
output_reader.Close()
