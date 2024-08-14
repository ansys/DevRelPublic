import os, time, tempfile
import pandas as pd
from parsimonious.grammar import Grammar

directory = r'D:\\'

required_data = {'Version': ['BUILD=', '=', 'str'],
                 'Nodes': ['Total number of nodes', '=', 'float'],
                 'Nodes < v20': ['Number of total nodes', '=', 'float'],
                 'Elements': ['Total number of elements', '=', 'float'],
                 'Elements < v20': ['Number of total elements', '=', 'float'],
                 'DOF': ['Number of DOF', '=', 'float'],
                 'Solver': ['Equation solver used', ':', 'str'],
                 'Cores': ['Total number of cores requested', ':', 'float'],
                 'Steps': ['SOLVE FOR LS', 'OF', 'float'],
                 'memory_available': ['physical memory available', ':',
                                      'float'],
                 'memory_used_old': ['Sum of memory used on all processes',
                                     ':', 'float'],
                 'memory_used': ['Maximum total memory used', ':', 'float'],
                 'Time': ['Total CPU time summed for all threads', ':',
                          'float']}


def finder(all_text, search_term):
    s1 = all_text.index(search_term)
    return_line = all_text[s1:all_text.index('\n', s1)]
    return return_line


def find_lines_with_text(file_path, search_strings):
    return_data = []
    temp_df = {}
    with open(file_path, 'r') as file:
        text = file.read()

    grammar = Grammar(r"""
        file = line+
        line = ~".*?\n"
    """)

    tree = grammar.parse(text)
    all_text = tree.text
    for key in search_strings:
        search_term = search_strings[key][0]  # .lower()
        try:
            line = finder(all_text, search_term)
            datapoint = line.split(search_strings[key][1])[1].split()[0]
            temp_df[key] = [datapoint]
        except:
            datapoint = 0

        data_frame = pd.DataFrame(temp_df)
        return_data.append(datapoint)
    return return_data, data_frame


### Main loop through files
data = []
d2 = []
df = pd.DataFrame()
t1 = time.time()
for root, dirs, files in os.walk(directory):
    for file in files:
        if file == ("solve.out"):
            path2file = os.path.join(root, file)
            temp_data, data_frame = find_lines_with_text(path2file,
                                                         required_data)
            data.append(temp_data)
df = pd.DataFrame(data, columns=required_data.keys())
t2 = time.time()

print('Time to get all data for %s files = %s seconds' % (
str(len(data)), str(round(t2 - t1))))

# clean up dataframe
solver_dict = {'Sparse': 0, 'PCG': 1,
               0: -1}  # Simplified numeric data for solver (0,1,-1)
a = df['Elements'].str.isnumeric() != False
df = df[a]
for key in required_data.keys():
    data_type = required_data[key][-1]
    if data_type == 'float':
        df[key] = df[key].astype(float)
# get elements and nodes to one column
df['Elements'] = df[['Elements', 'Elements < v20']].max(axis=1)
df['Nodes'] = df[['Nodes', 'Nodes < v20']].max(axis=1)
df['memory_used'] = df[['memory_used', 'memory_used_old']].max(axis=1)
solver_data = [solver_dict[i] for i in df['Solver']]

df['Solver'] = solver_data
df['Solver'] = df['Solver'].astype(int)
df = df.drop(['Elements < v20', 'Nodes < v20', 'memory_used'], axis=1)

user = os.getlogin()
out_file_name = 'solves_out_' + user + '.json'
path2result = os.path.join(tempfile.gettempdir(), out_file_name)
df.to_json(path2result)
print(path2result)
print('Data saved to: ' + path2result)
