import sys

# BINARY_FOLDER_PATH : path to the folder containing the Postprocessor dll files
# RESULT_FILE_READER_API_MODULE_PATH : path to the folder containing ResultFileReaderAPI.py
sys.path.append(BINARY_FOLDER_PATH)
sys.path.append(RESULT_FILE_READER_API_MODULE_PATH)

# Import necessary modules
from ResultFileReaderAPI import *

# RESULT_FILE_PATH : .dfr result file path
output_reader = OutputReader(RESULT_FILE_PATH)

print ("===GetVector===")
# target - Specifies the name of vector displayable entity
# path - Specifies characteristc on vector display
target_name = "TJ_01"
# vector_path : fixed vector characteristic path for the Action Force result
vector_path = "Action Force"
vectors = output_reader.GetVector(target_name, vector_path)
for vector in vectors:
    print(f"Vector : {vector.Key}")
    animation_data = vector.Value

    positions = len(list(animation_data.Positions))
    for i in range(positions):
        first_positions = len(list(animation_data.Positions[i]))
        for j in range(first_positions):
            second_positions = list(animation_data.Positions[i][j])
            print("Positions :", *second_positions, sep=',')

    vectors = len(list(animation_data.Vectors))
    for i in range(vectors):
        first_vectors = len(list(animation_data.Vectors[i]))
        for j in range(first_vectors):
            second_vectors = list(animation_data.Vectors[i][j])
            print("Vectors :", *second_vectors, sep=',')
# Close
output_reader.Close()
