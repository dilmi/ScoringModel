__author__ = 'dilmiperera'
import numpy as np
import pandas as pd
from sklearn.ensemble import AdaBoostClassifier
from sklearn import cross_validation
from sklearn.metrics import confusion_matrix, roc_auc_score
from sklearn.grid_search import GridSearchCV


num_of_crossval = 10

mutations_df = pd.read_csv('Data/Data.csv', delimiter=',')
data_labels = ['GeneMappedUsingFANTOM5', 'DistanceToTSS', 'DHS', 'H3K4me1', 'H3K4me3', 'H3K27ac',
               'Conservation_at_mutation', 'BackgroundConservation', 'Fantom5_Promoter', 'Fantom5_Enhancer',
               'MotifScore']
X_all = mutations_df[data_labels].as_matrix()
y_all = mutations_df.Pathogenicity.values

# Split testing and training sets
X_train, X_test, y_train, y_test = cross_validation.train_test_split(X_all, y_all, test_size=0.33, random_state=42)

#To store output
classArray_train = np.zeros(len(y_train))
classArray_test = np.zeros((len(y_test), num_of_crossval))

kf = cross_validation.StratifiedKFold(y_train, n_folds=num_of_crossval, shuffle=False)

#Number of trees for grid search(drawing the graph to find the optimum parameter)
tune_params = [{'n_estimators': [50, 100, 300, 500]}]
#tune_params = [{'C': [10.0, 1.0, 0.1, 0.01, 0.001]}]
final_Results_mat = []
fold = 0
for train_index, test_index in kf:
    print('Fold : {0}'.format(fold))
    X_train_kf, X_test_kf, y_train_kf, y_test_kf = X_train[train_index], X_train[test_index], y_train[train_index], \
                                                   y_train[test_index]
    print('Training Starting ... ')
    gsCV = GridSearchCV(AdaBoostClassifier(n_estimators=300),tune_params,cv=len(tune_params[0]['n_estimators']), n_jobs=-1)
    print('Create model ... ')
    gsCV.fit(X_train_kf, y_train_kf)
    print('Training complete ... ')
    print(gsCV.best_params_)

    for params, mean_score, scores in gsCV.grid_scores_:
        print("%0.3f (+/-%0.03f) for %r" % (mean_score, scores.std() * 2, params))

    predictedClass = gsCV.predict(X_test_kf)
    cm = confusion_matrix(y_test_kf, predictedClass)
    print(cm)

    predictedClassTEST = gsCV.predict(X_test)
    cmTEST = confusion_matrix(y_test, predictedClassTEST)
    print(cmTEST)

    classArray_train[test_index] = predictedClass
    classArray_test[:, fold] = predictedClassTEST

    predictedProb = gsCV.predict_proba(X_test_kf)
    predClassesProb_test = gsCV.predict_proba(X_test)

    auc_TRAIN = roc_auc_score(y_test_kf, predictedProb[:, 1])
    auc_TEST = roc_auc_score(y_test, predClassesProb_test[:, 1])

    print(gsCV.best_estimator_.feature_importances_)

    fold+=1

print("Validation Final : ")
print(confusion_matrix(y_train, classArray_train))

classArray_test_final = np.mean(classArray_test, axis=1) > 0.5
print("Test Final : ")
print(confusion_matrix(y_test, classArray_test_final))

