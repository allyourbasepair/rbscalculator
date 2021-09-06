from stats import do_do_do_compute_stats
import numpy as np
import stats
from itertools import permutations, combinations
import pandas as pd
from tabulate import tabulate
import matplotlib.pyplot as plt
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import RFECV
from sklearn.linear_model import (
    LinearRegression, TheilSenRegressor, RANSACRegressor, HuberRegressor)
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline
from sklearn.metrics import mean_squared_error

def regression(csv_file, prot_mean_col_idx, prot_std_col_idx, energy_term_col_idxs):
  prot_mean, prot_std = np.genfromtxt(csv_file, delimiter=',',
      skip_header=1, usecols=(prot_mean_col_idx, prot_std_col_idx), unpack=True)
  prot_mean = np.log(prot_mean)

  energy_terms = np.genfromtxt(csv_file, delimiter=',',
      skip_header=1, usecols=energy_term_col_idxs)
  energy_terms_labels = np.genfromtxt(csv_file, delimiter=',', max_rows=1, usecols=energy_term_col_idxs, dtype=str)
  print(energy_terms_labels)

  X_train, X_test, y_train, y_test = train_test_split(energy_terms, prot_mean, test_size=0.33, random_state=42)

  #Create and fit the model
  # model = Pipeline([
  #   ('scaler', StandardScaler()),
  #   ('linear regression', LinearRegression())
  # ])
  # model = LinearRegression()
  estimators = [('LinearRegression', LinearRegression()),
              ('TheilSenRegressor', TheilSenRegressor(random_state=7)),
              ('RANSACRegressor', RANSACRegressor(random_state=7)),
              ('HuberRegressor', HuberRegressor())]
  for name, estimator in estimators:
    # model = make_pipeline(PolynomialFeatures(3), estimator)
    model = estimator
    # model = RANSACRegressor(random_state=0)
    print(name)
    #Fit the model using the training data
    model.fit(X_train,y_train)

    #Predict unseen data
    scores = model.score(X_test, y_test)

    print("R^2", scores)
    mse = mean_squared_error(model.predict(X_test), y_test)

    # print(metrics.mean_squared_error(y_test, y_predicted))
    print(mse)

  # # Create the RFE object and compute a cross-validated score.
  # svc = SVC(kernel="linear")
  # # The "accuracy" scoring is proportional to the number of correct
  # # classifications

  # min_features_to_select = 2  # Minimum number of features to consider
  # rfecv = RFECV(estimator=svc, step=1, cv=StratifiedKFold(2),
  #               scoring='accuracy',
  #               min_features_to_select=min_features_to_select)
  # rfecv.fit(energy_terms, prot_mean)

  # print("Optimal number of features : %d" % rfecv.n_features_)

  # # Plot number of features VS. cross-validation scores
  # plt.figure()
  # plt.xlabel("Number of features selected")
  # plt.ylabel("Cross validation score (nb of correct classifications)")
  # plt.plot(range(min_features_to_select,
  #               len(rfecv.grid_scores_) + min_features_to_select),
  #         rfecv.grid_scores_)
  # plt.show()

  # # col_idxs = range(0,len(energy_terms_labels))
  # # all_stats_dict = []
  # # for i in range(2, len(col_idxs)+1):
  # #   # if i == 3:
  # #   #   exit()
  # #   for comb in combinations(col_idxs, i):
  # #     energy_labels = [energy_terms_labels[comb[0]]]
  # #     energy_values = energy_terms[comb[0]]
  # #     for idx in comb[1:]:
  # #       energy_labels.append(energy_terms_labels[idx])
  # #       energy_values += energy_terms[idx]
  # #     reference = " + ".join(energy_labels)
  # #     print('processing {}'.format(reference))
  # #     stats_for_this_param_set = do_do_do_compute_stats(energy_values, prot_mean, prot_std, reference, linear_fit=True)
  # #     all_stats_dict.append(stats_for_this_param_set)
  # # pd.options.display.max_columns = None
  # # df = pd.DataFrame(all_stats_dict)
  # # print(tabulate(df,headers='keys',tablefmt='psql'))







if __name__ == "__main__":
  dataset_with_properties = "./salis_lab_v2_1/dataset_with_properties/train_with_salis_sd_algo.csv"
    # new
  prot_mean_idx, prot_std_idx = 4, 5
  energy_term_idxs = (12, 18, 22, 24, 25, 27, 30, 41)
  # energy_term_idxs = (12, 15, 19, 22, 24, 25, 27, 31, 32, 33, 34, 35, 36, 37, 41)


  regression(dataset_with_properties, prot_mean_idx, prot_std_idx, energy_term_idxs)

  # regression((1,2,3,45))
  # stats.linear_complete(x_valid,y_valid,std_valid,xScale,yScale,self.models[m].a1)
