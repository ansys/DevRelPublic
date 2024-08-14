# imports
from sklearn.model_selection import train_test_split
from sklearn import ensemble
from sklearn.metrics import mean_absolute_error
from sklearn.model_selection import GridSearchCV
import pandas as pd
import time
import joblib

path2result = "D:\PYTHON\ResultFileGetData\jsons\solves_out_mcapella.json"
df = pd.read_json(path2result)
print(df.columns)
df_features = df.copy()
del df_features['Version']
del df_features['memory_available']
del df_features['memory_used_old']

# determine which column has the most data!
sort_results = []
columns_to_filter = ['Nodes','Elements','DOF','Time']
for column in columns_to_filter:
    count = (df_features[column]!=0).astype(int).sum()
    print(count)
    sort_results.append([column,count])
print(sort_results)

# Either nodes or elements have the max number of data points, delete DOF
del df_features['DOF']

# sort for nodes/elements/time ==0 and remove by index
index_count = []
for index,row in df_features.iterrows():
    if row['Nodes'] == 0 or row['Elements'] == 0 or row['Time'] == 0:
        index_count.append(index)

df_features.drop(axis = 0, index = index_count,inplace = True)

# Quick plot to see data
import matplotlib.pyplot as plt
x_feature = 'Elements'
y_feature = 'Time'
df_features.plot(kind = 'scatter', x = x_feature, y = y_feature, grid = True)
plt.show()

# scikit learn data set up
y = df_features['Time'].values
del df_features['Time']
X = df_features.values



test_size = 0.30
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=0)

# Fit regression model
model = ensemble.GradientBoostingRegressor(
    n_estimators = 1000,
    max_depth = 6,
    min_samples_leaf = 9,
    learning_rate = 0.1,
    max_features = 0.1,
    loss = 'huber')

model.fit(X_train,y_train)
## Find the error rate on the training set
mse = mean_absolute_error(y_train,model.predict(X_train))
print("Training Set Mean Absolute Error: %.4f" % mse)

# Find the error rate on the test set
mse = mean_absolute_error(y_test,model.predict(X_test))
print("Test Set Mean Absolute Error: %.4f" % mse)

n_cpus = 6
param_grid = {
    'n_estimators': [500, 1000, 2000 ],
    'max_depth': [3, 6, 12],
    'min_samples_leaf': [4, 9, 18],
    'learning_rate': [ 0.05, 0.1, 0.2],
    'max_features': [0.05, 0.1, 0.2],
    'loss': [ 'huber']
}

t1=time.time()
# Define the grid search we want to run. Run it with the maximum cpus in parallel 
gs_cv = GridSearchCV(model, param_grid, n_jobs = n_cpus)

# Run the grid search - on only the training data!
gs_cv.fit(X_train, y_train)

t2 = time.time()

# Print the parameters that gave us the best fit
print(gs_cv.best_params_)
print('time : ' + str(round(t2-t1)/60) + ' minutes' )

# Find the error rate on the training set using the best parameters
mse = mean_absolute_error(y_train, gs_cv.predict(X_train))
print("Training Set Mean Absolute Error: %.4f" % mse)
# Find the error rate on the test set using the best parameters
mse = mean_absolute_error(y_test, gs_cv.predict(X_test))
print("Test Set Mean Absolute Error: %.4f" % mse)

# Save the trained model to a file so we can use it in other instances
joblib.dump(model, 'ansys_data_best_parameters.pkl')

# Create a numpy array based on the model's feature importances
importance = model.feature_importances_

# Sort the feature labels based on the feature importance rankings from the model
feature_indexes_by_importance = importance.argsort()

# These are the feature labels from our data set
feature_labels = df_features.columns
print(feature_labels)
# Print each feature label, from most important to least important (reverse order)
for index in feature_indexes_by_importance:
    print("{} - {:.2f}%".format(feature_labels[index], (importance[index] * 100.0)))

# Load the model we trained previously
model = joblib.load('ansys_data_best_parameters.pkl')

### predict
nodes = 4e6
elements = 1e6
solver = 1
cores = 8
steps = 1

analysis_inputs = [
    # exact same order as training pandas column headers/features!
    nodes, # nodes
    elements, # Elements
    solver, # solver
    cores, # cores
    steps # Steps
    ]

analyses_inputs = [analysis_inputs]
predicted_time = model.predict(analyses_inputs)
predicted_value = predicted_time[0]

print("This analysis of {:,d} nodes, {:,d} elements,{:,d} cores has an estimated time of {:,.2f} seconds"
      .format(int(nodes),int(elements),cores,predicted_value))