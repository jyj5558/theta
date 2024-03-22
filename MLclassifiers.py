"""=================================================================================================
    Title: ThetaProject_MachingLearningClassifiers

    This python script will parse the input file and fit several machine learning classifier models
    to classify a Red List category based on GD metrics.

    Jong Yoon Jeon     January 08 2023

================================================================================================="""
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

#data load and check
df = pd.read_csv("C:/Users/jyj55/OneDrive - purdue.edu/DeWoody_Lab/Dissertation/Theta project/Statistics/rf_input.csv")
df.head()
df.tail() #too many NaN rows and columns -> clean NaN rows and columns as below
real_cols = [x for x in df.columns if not x.startswith("Unnamed: ")]
df = df[real_cols]
df = df.dropna(subset=['NCBI_name'])
df.head()
df.tail()
df.info()
print(df.isnull().sum())
df.describe().T

g = sns.pairplot(df, hue='IUCN_category')
g.fig.suptitle("Scatterplot and histogram of pairs of variables color coded by IUCN category",
               fontsize = 14, # defining the size of the title
               y=1.05) # y = definig title y position (height)

# Function to classify categories
def classify_iucn(row):
    if row['IUCN_category'] in ['CR', 'EN', 'VU']:
        return 'Threatened'
    elif row['IUCN_category'] in ['NT', 'LC']:
        return 'Non-Threatened'
    else:
        return 'Unknown'

# Apply the function to create a new column
df['IUCN.b'] = df.apply(classify_iucn, axis=1)
df.head()
df.tail()

df['IUCN.b'].unique()
df['IUCN.b'] = df['IUCN.b'].replace('Non-Threatened', 0).replace('Threatened', 1)

seed = 106

##K-Nearest Neighbors
#preprocessing for model
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split, GridSearchCV

#Xh = df.iloc[:,5] #define predictor (Het)
Xt = df.iloc[:,3] #define predictor (Theta)
Xf = df.iloc[:,5] #define predictor (F100)
y = df['IUCN.b'] #define response


Xt_train, Xt_test, y_train, y_test = train_test_split(Xt, y, stratify=y, test_size=0.3, random_state=seed) #random_state = seed
Xt_train.shape, Xt_test.shape # check the shape of X_train and X_test
Xf_train, Xf_test, y_train, y_test = train_test_split(Xf, y, stratify=y, test_size=0.3, random_state=seed) #random_state = seed
Xf_train.shape, Xf_test.shape # check the shape of X_train and X_test

Xt_train_2d = np.array(Xt_train).reshape((len(Xt_train),1)) # make a 2D numpy array
Xt_test_2d = np.array(Xt_test).reshape((len(Xt_test),1)) # make a 2D numpy array
Xf_train_2d = np.array(Xf_train).reshape((len(Xf_train),1)) # make a 2D numpy array
Xf_test_2d = np.array(Xf_test).reshape((len(Xf_test),1)) # make a 2D numpy array

# scaling predictors
scaler = StandardScaler()
scaler.fit(Xt_train_2d) #fit scaling only to training data to prevent data leakage (influnence of test data)
Xt_train_2d = scaler.transform(Xt_train_2d)
Xt_test_2d = scaler.transform(Xt_test_2d)
scaler.fit(Xf_train_2d) #fit scaling only to training data to prevent data leakage (influnence of test data)
Xf_train_2d = scaler.transform(Xf_train_2d)
Xf_test_2d = scaler.transform(Xf_test_2d)

#tune K-Nearest Neighbors model
from sklearn.neighbors import KNeighborsClassifier
#from sklearn.metrics import f1_score

param_grid = {
    'n_neighbors' : [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20],
    'weights' : ['uniform','distance'],
    'metric' : ['minkowski','euclidean','manhattan','hamming']}
# Base model
knn_grid = KNeighborsClassifier()
# Instantiate the grid search model
grid_knnt_search = GridSearchCV(estimator = knn_grid, param_grid = param_grid,
                          cv = 5, n_jobs = 8, verbose = 2, scoring='accuracy')
grid_knnt_search.fit(Xt_train_2d, y_train)
grid_knnt_search.best_score_
grid_knnt_search.best_estimator_
print(grid_knnt_search.score(Xt_train_2d, y_train))
print(grid_knnt_search.score(Xt_test_2d, y_test))

grid_knnf_search = GridSearchCV(estimator = knn_grid, param_grid = param_grid,
                          cv = 5, n_jobs = 8, verbose = 2, scoring='accuracy')
grid_knnf_search.fit(Xf_train_2d, y_train)
grid_knnf_search.best_score_
grid_knnf_search.best_estimator_
print(grid_knnf_search.score(Xf_train_2d, y_train))
print(grid_knnf_search.score(Xf_test_2d, y_test))

#K-Nearest Neighbors model fitting
KNNclf_t = KNeighborsClassifier(n_neighbors=17)
KNNclf_t.fit(Xt_train_2d, y_train)
yt_pred = KNNclf_t.predict(Xt_test_2d)

KNNclf_f = KNeighborsClassifier(n_neighbors=6)
KNNclf_f.fit(Xf_train_2d, y_train)
yf_pred = KNNclf_f.predict(Xf_test_2d)

#model evaluation
from sklearn.metrics import classification_report, confusion_matrix

print(KNNclf_t.score(Xt_train_2d, y_train))
print(KNNclf_t.score(Xt_test_2d, y_test))
#print("Accuracy:", KNNclf_h.accuracy_score(y_test, yh_pred))
print(KNNclf_f.score(Xf_train_2d, y_train))
print(KNNclf_f.score(Xf_test_2d, y_test))
#print("Accuracy:", KNNclf_f.accuracy_score(y_test, yf_pred))

print(classification_report(y_test,yt_pred))
print(classification_report(y_test,yf_pred))

cm_t = confusion_matrix(y_test, yt_pred)
sns.heatmap(cm_t, annot=True, fmt='d').set_title('IUCN category confusion matrix (0 = Non-Threatened, 1 = Threatened) by Het')
cm_f = confusion_matrix(y_test, yf_pred)
sns.heatmap(cm_f, annot=True, fmt='d').set_title('IUCN category confusion matrix (0 = Non-Threatened, 1 = Threatened) by F100')



##Support Vector Machine
#preprocessing for model
from sklearn.model_selection import train_test_split, GridSearchCV
X = df.iloc[:,[3,5]] #define predictors
y = df['IUCN.b'] #define response
X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y, test_size=0.3, random_state=seed) #random_state = seed
X_train.shape, X_test.shape # check the shape of X_train and X_test

# scaling predictors
scaler = StandardScaler()
scaler.fit(X_train) #fit scaling only to training data to prevent data leakage (influence of test data)
X_train = scaler.transform(X_train)
X_test = scaler.transform(X_test)

#tune linear Support Vector Machine model
from sklearn.svm import SVC
param_grid = {
    'C': [0.1,1,10,100,1000],
    'gamma': ['scale','auto',1,0.1,0.01,0.001,0.001],
    'kernel': ['linear'],
    'degree':[1,2,3,4,5,6]}

# Base model
svc_grid = SVC()
# Instantiate the grid search model
grid_svc_search = GridSearchCV(estimator = svc_grid, param_grid = param_grid,
                          cv = 5, n_jobs = 8, verbose = 2, scoring='accuracy')
grid_svc_search.fit(X_train, y_train)
grid_svc_search.best_score_
grid_svc_search.best_estimator_
print(grid_svc_search.score(X_train, y_train))
print(grid_svc_search.score(X_test, y_test))

#Support Vector Machine model fitting - linear
SVCclf = SVC(C=0.1, degree=1, kernel='linear')
SVCclf.fit(X_train, y_train)

# only available in case of linear kernel, provides the weight assigned to the features
SVCclf.coef_
# get the value of other attributes
#SVCclf.predict([[x,y]])
# the number of support vectors for each class
SVCclf.n_support_
# returns the support vectors
SVCclf.support_vectors_
# represents the independent term (constant) in decision function
SVCclf.intercept_
# the output would be 0 if it is correctly fitted. The output would be 1 if it is incorrectly fitted.
SVCclf.fit_status_

# make predictions
y_pred = SVCclf.predict(X_test)

# model evaluation
from sklearn.metrics import classification_report, confusion_matrix

print(classification_report(y_test,y_pred))
cm = confusion_matrix(y_test, y_pred)
sns.heatmap(cm, annot=True, fmt='d').set_title('IUCN category confusion matrix (0 = Non-Threatened, 1 = Threatened)')

print(SVCclf.score(X_train, y_train))
print(SVCclf.score(X_test, y_test))
#print("Accuracy:", SVCclf.accuracy_score(y_test, y_pred))



#tune Support Vector Machine model individually for Theta and F100
Xt = df.iloc[:,3] #define predictor (Theta)
Xf = df.iloc[:,5] #define predictor (F100)
y = df['IUCN.b'] #define response
Xt_train, Xt_test, y_train, y_test = train_test_split(Xt, y, stratify=y, test_size=0.3, random_state=seed) #random_state = seed
Xt_train.shape, Xt_test.shape # check the shape of X_train and X_test
Xf_train, Xf_test, y_train, y_test = train_test_split(Xf, y, stratify=y, test_size=0.3, random_state=seed) #random_state = seed
Xf_train.shape, Xf_test.shape # check the shape of X_train and X_test

Xt_train_2d = np.array(Xt_train).reshape((len(Xt_train),1)) # make a 2D numpy array
Xt_test_2d = np.array(Xt_test).reshape((len(Xt_test),1)) # make a 2D numpy array
Xf_train_2d = np.array(Xf_train).reshape((len(Xf_train),1)) # make a 2D numpy array
Xf_test_2d = np.array(Xf_test).reshape((len(Xf_test),1)) # make a 2D numpy array

# scaling predictors
scaler = StandardScaler()
scaler.fit(Xt_train_2d) #fit scaling only to training data to prevent data leakage (influnence of test data)
Xt_train_2d = scaler.transform(Xt_train_2d)
Xt_test_2d = scaler.transform(Xt_test_2d)
scaler.fit(Xf_train_2d) #fit scaling only to training data to prevent data leakage (influnence of test data)
Xf_train_2d = scaler.transform(Xf_train_2d)
Xf_test_2d = scaler.transform(Xf_test_2d)

param_grid = {
    'C': [0.1,1,10,100,1000],
    'gamma': ['scale','auto',1,0.1,0.01,0.001,0.001],
    'kernel': ['poly','sigmoid','rbf','linear'],
    'degree':[1,2,3,4,5,6]}

# Base model
svc_grid = SVC()
# Instantiate the grid search model
grid_svct_search = GridSearchCV(estimator = svc_grid, param_grid = param_grid,
                          cv = 5, n_jobs = 8, verbose = 2, scoring='accuracy')
grid_svct_search.fit(Xt_train_2d, y_train)
grid_svct_search.best_score_
grid_svct_search.best_estimator_
print(grid_svct_search.score(Xt_train_2d, y_train))
print(grid_svct_search.score(Xt_test_2d, y_test))

grid_svcf_search = GridSearchCV(estimator = svc_grid, param_grid = param_grid,
                          cv = 5, n_jobs = 8, verbose = 2, scoring='accuracy')

grid_svcf_search.fit(Xf_train_2d, y_train)
grid_svcf_search.best_score_
grid_svcf_search.best_estimator_
print(grid_svcf_search.score(Xf_train_2d, y_train))
print(grid_svcf_search.score(Xf_test_2d, y_test))

#Individual Support Vector Machine model fitting
SVCclf_t = SVC(C=1, degree=1, gamma='auto', kernel='sigmoid')
SVCclf_t.fit(Xt_train_2d, y_train)
yt_pred = SVCclf_t.predict(Xt_test_2d)

SVCclf_f = SVC(C=100, degree=1)
SVCclf_f.fit(Xf_train_2d, y_train)
yf_pred = SVCclf_f.predict(Xf_test_2d)

#model evaluation
print(SVCclf_t.score(Xt_train_2d, y_train))
print(SVCclf_t.score(Xt_test_2d, y_test))
#print("Accuracy:", KNNclf_h.accuracy_score(y_test, yh_pred))
print(SVCclf_f.score(Xf_train_2d, y_train))
print(SVCclf_f.score(Xf_test_2d, y_test))
#print("Accuracy:", KNNclf_f.accuracy_score(y_test, yf_pred))

print(classification_report(y_test,yt_pred))
print(classification_report(y_test,yf_pred))

cm_t = confusion_matrix(y_test, yt_pred)
sns.heatmap(cm_t, annot=True, fmt='d').set_title('IUCN category confusion matrix (0 = Non-Threatened, 1 = Threatened) by Het')
cm_f = confusion_matrix(y_test, yf_pred)
sns.heatmap(cm_f, annot=True, fmt='d').set_title('IUCN category confusion matrix (0 = Non-Threatened, 1 = Threatened) by F100')



##Random Forest
#preprocessing for model
from sklearn.model_selection import train_test_split, RandomizedSearchCV, GridSearchCV
from sklearn.ensemble import RandomForestClassifier

X = df.iloc[:,[3,5]] #define predictors
y = df['IUCN.b'] #define response
X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y, test_size=0.3, random_state=seed) #random_state = seed
X_train.shape, X_test.shape # check the shape of X_train and X_test

#tune Random Forest model
# Number of trees in random forest
n_estimators = np.linspace(10, 10000, int((10000-10)/10) + 1, dtype=int)
# Number of features to consider at every split
max_features = ['auto', 'sqrt']
# Maximum number of levels in tree
max_depth = [1, 2, 5, 10, 20, 30, 50, 75, 100]
# Minimum number of samples required to split a node
# min_samples_split = [int(x) for x in np.linespace(start = 2, stop = 10, num = 9)]
min_samples_split = [1, 2, 5, 10, 15, 20, 30, 50]
# Minimum number of samples required at each leaf node
min_samples_leaf = [1, 2, 3, 4, 5, 10]
# Method of selecting samples for training each tree
bootstrap = [True, False]
# Criterion
criterion=['gini', 'entropy']

random_grid = {'n_estimators': n_estimators,
               'max_features': max_features,
               'max_depth': max_depth,
               'min_samples_split': min_samples_split,
               'min_samples_leaf': min_samples_leaf,
               'bootstrap': bootstrap,
               'criterion': criterion}

rf_base = RandomForestClassifier()
rf_random = RandomizedSearchCV(estimator = rf_base,
                               param_distributions = random_grid,
                               n_iter = 30, cv = 5,
                               verbose=2,
                               random_state=seed, n_jobs = 4)
rf_random.fit(X_train, y_train)
rf_random.best_params_ #{'n_estimators': 920, 'min_samples_split': 15, 'min_samples_leaf': 4, 'max_features': 'sqrt', 'max_depth': 1, 'criterion': 'gini', 'bootstrap': True}
                       #{'n_estimators': 6410, 'min_samples_split': 1, 'min_samples_leaf': 4, 'max_features': 'sqrt', 'max_depth': 1, 'criterion': 'entropy', 'bootstrap': True}
                       #{'n_estimators': 290, 'min_samples_split': 15, 'min_samples_leaf': 4, 'max_features': 'sqrt', 'max_depth': 50, 'criterion': 'gini', 'bootstrap': False}
print(rf_random.score(X_train, y_train))
print(rf_random.score(X_test, y_test))

#param_grid = {
#    'n_estimators': np.linspace(500, 1000, 5, dtype = int),
#    'max_depth': [1, 2, 5, 10, 20, 30, 50],
#    'min_samples_split': [1, 2, 5, 10, 15, 20],
#    'min_samples_leaf': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
#    'bootstrap': [True, False],
#    'criterion': ['gini', 'entropy']}
#param_grid = {
#    'n_estimators': np.linspace(100, 5000, 50, dtype = int),
#    'max_depth': [1],
#    'min_samples_split': [1, 2, 5, 10, 15, 20],
#    'min_samples_leaf': [1, 2, 4],
#    'bootstrap': [True],
#    'criterion': ['gini', 'entropy']}
param_grid = {
    'n_estimators': np.linspace(100, 5000, 50, dtype = int),
    'max_depth': [1, 2, 5, 10, 20, 30, 50, 75, 100],
    'min_samples_split': [1, 2, 5, 10, 15, 20],
    'min_samples_leaf': [1, 2, 4],
    'bootstrap': [True, False],
    'criterion': ['gini', 'entropy']}
# Base model
rf_grid = RandomForestClassifier()
# Instantiate the grid search model
grid_rf_search = GridSearchCV(estimator = rf_grid, param_grid = param_grid,
                          cv = 5, n_jobs = 8, verbose = 2)
grid_rf_search.fit(X_train, y_train)

grid_rf_search.best_estimator_
print(grid_rf_search.score(X_train, y_train))
print(grid_rf_search.score(X_test, y_test))

#Random Forest model fitting
#RFclf = RandomForestClassifier(max_depth=1, min_samples_leaf=2, min_samples_split=5, n_estimators=750, random_state=seed) #random_state = seed
#RFclf = RandomForestClassifier(max_depth=1, min_samples_leaf=2, min_samples_split=5, random_state=seed) #random_state = seed
RFclf = RandomForestClassifier(bootstrap=False, max_depth=10, min_samples_split=1, n_estimators=200, random_state=seed) #random_state = seed
RFclf.fit(X_train, y_train)
y_pred =RFclf.predict(X_test)

#from sklearn import tree
#features = X.columns.values # The name of each column
#classes = ['1', '2'] # The name of each class

#for estimator in clf.estimators_[0:10]:
#    print(estimator)
#    plt.figure(figsize=(12,6))
#    tree.plot_tree(estimator,
#                   feature_names=features,
#                   class_names=classes,
#                   fontsize=8,
#                   filled=True,
#                   rounded=True)
#    plt.show()

#model evaluation
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score
cm = confusion_matrix(y_test, y_pred)
sns.heatmap(cm, annot=True, fmt='d').set_title('IUCN category confusion matrix (0 = Non-Threatened, 1 = Threatened)')
print(classification_report(y_test,y_pred))

print(RFclf.score(X_train, y_train))
print(RFclf.score(X_test, y_test))
print("Accuracy:", RFclf.accuracy_score(y_test, y_pred))

pd.DataFrame(RFclf.feature_importances_, index=X_train.columns).sort_values(by=0, ascending=False)
features_df = pd.DataFrame({'features': RFclf.feature_names_in_, 'importances': RFclf.feature_importances_ })
features_df_sorted = features_df.sort_values(by='importances', ascending=False)
print(features_df_sorted)
g = sns.barplot(data=features_df_sorted, x='importances', y ='features', palette="rocket")
sns.despine(bottom = True, left = True)
g.set_title('Feature importances')
g.set(xlabel=None)
g.set(ylabel=None)
g.set(xticks=[])
for value in g.containers:
    g.bar_label(value, padding=2)

#Classification report information
#accuracy = (true positives + true negatives) / (true positives + true negatives + false pasitives + false negatives)
#precision = true positives / (true positives + false positives)
#recall (sensitivity) = true positives / (true positives + false negatives)
#f1-score = 2 ∗ (precision ∗ recall) / (precision + recall)

#save the result
import joblib
joblib.dump(RFclf, "./random_forest_RedList.joblib")
