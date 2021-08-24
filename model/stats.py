import numpy as np
import scipy
import scipy.stats
import pandas as pd
from tabulate import tabulate

'''
This file computes statistics for the models the RBS calculator.

Full credit for the functions herein go to Alex Reis (@reisalex) and The Salis
Lab. All functions in this file have been taken from the `SynBioMTS`
GitHub repo (https://github.com/reisalex/SynBioMTS), except for the function
`linear_complete` which has minor modifications.

Please read the documentation of the `linear_complete` function for details
on how to compute statistics of your RBS model.
'''




def area_under_ROC_curve(predictions, observations, cutoff=None, n=50):
    '''Calculate ROC curve and the area under the curve (auc)

    Returns:
        auroc (float) = area under the ROC curve
        fpr (array) = false positive rates
        tpr (array) =  true positive rates
        thresholds (array) = thresholds from ROC curve'''

    # define outcomes as TRUE/FALSE defined by cutoff
    if cutoff == None:
        cutoff = np.median(outcomes)
    outcomes = (observations > cutoff)

    # preallocate arrays and define n-many threshodls
    fpr = np.zeros(n, dtype=np.float64)
    tpr = np.zeros(n, dtype=np.float64)
    thresholds = np.linspace(0.9*min(observations),
                             1.1*max(observations), n, dtype=np.float64)

    for i in xrange(0, n):

        threshold = thresholds[i]
        scores = predictions > threshold

        # calculate confusion matrix
        a = sum(scores & outcomes)
        b = sum(scores & ~outcomes)
        c = sum(~scores & outcomes)
        d = sum(~scores & ~outcomes)

        # calculate fpr and tpr for each threshold
        fpr[i] = float(b)/(b+d)
        tpr[i] = float(a)/(a+c)

    # Remove nans
    nans = np.isnan(fpr) + np.isnan(tpr)
    x = np.array(fpr, dtype=np.float64)
    y = np.array(tpr, dtype=np.float64)
    x = x[~nans]
    y = y[~nans]

    # Sort by increasing x
    indx = np.argsort(x)
    x = x[indx]
    y = y[indx]

    # Remove repeated elements
    x, indx = np.unique(x, return_index=True)
    y = y[indx]

    # Add (0,0) or (1,1) if not present!
    if x[0] != 0.0:
        x = np.append([0.0], x)
        y = np.append([0.0], y)

    if x[-1] != 1.0:
        x = np.append(x, 1.0)
        y = np.append(y, 1.0)

    # Calculate area under the ROC curve
    fpr = x
    tpr = y
    auroc = np.trapz(tpr, fpr)

    return auroc, fpr, tpr, thresholds


def pdist(data, b, nbins=100, make_outliers_rand=True):
    '''Calculate the discrete probability distribution function
    of a set of data defined over the range (-b,b) with nbins
    Inputs:
        ignore_outliers (bool) = if True ignore outliers,
                                 if False add random information to pdf
    Returns:
        pdf (array) = discrete pdf with probability of event i as pdf[i]
        x (array) = midpoints of the bins of the histogram that defines pdf'''

    # remove outliers outside of range (-b,b)
    data2 = data[~(np.abs(data) > b)]

    # calculate Pset (the pdf of the set of data points within the range)
    edges = np.linspace(-b, b, nbins+1, dtype=np.float64)
    Pset, edges = np.histogram(data2, edges, density=True)
    fset = float(len(data2))/len(data)
    x = (edges[0:-1]+edges[1:])/2

    if make_outliers_rand:
        # calculate Poutliers (pdf corresponding to random guesses)
        Pout = np.ones(len(x), dtype=np.float64)/len(x)
        fout = 1 - fset
        pdf = fset*Pset + fout*Pout
    else:
        pdf = Pset

    return (pdf, x)


def entropy(pk, qk=None, base=None):
    '''Calculate the entropy of a distribution for given probability values

    Inputs:
        pk (array)   = probability distribution of a discrete distribution
        qk (array)   = pdf of a second sequence with which pk is compared against
        base (float) = the logarithmic base to use (default=e)
    Returns:
        S (float) = Shannon entropy if qk is None else Kullback-Leibler divergence'''
    return scipy.stats.entropy(pk, qk, base)


def normKLdiv(data, b):
    '''Calculate normalized Kullback-Leibler divergence

    normKLdiv uses data and calculates the normKLdiv relative to a
    random model (Q) over bin range (-b,b).
    Returns:
        NKLdiv (float)   = normalized Kullback-Leibler divergence
        KLdiv (float)    = KL-divergence
        KLdivmax (float) = KL-divergence for a perfect model (dirac delta fxn)'''
    (P, x) = pdist(data, b)
    Pmax = np.zeros(len(x), dtype=np.float64)
    Pmax[0] = 1.0
    Q = np.ones(len(x), dtype=np.float64)/len(x)
    KLdiv = entropy(P, Q)
    KLdivmax = entropy(Pmax, Q)
    NKLdiv = KLdiv/KLdivmax
    return (NKLdiv, KLdiv, KLdivmax)


def correlation(x, y, name="Pearson"):
    # Linear or rank correlation
    # Returns:
    #   r,rho, or tau (float) = the test statistic
    #   pvalue (float) = the p-value for the test

    assert len(x) == len(y), "Arrays x & y must be equal length."

    if name == "Pearson":
        (r, pvalue) = scipy.stats.pearsonr(x, y)
        return r, pvalue

    elif name == "Spearman":
        (rho, pvalue) = scipy.stats.spearmanr(x, y)
        return rho, pvalue

    elif name == "Kendall":
        (tau, pvalue) = scipy.stats.kendalltau(x, y)
        return tau, pvalue

    else:
        error = ("The {} correlation is not available."
                 "Please use 'Pearson', 'Spearman', or 'Kendall.'")
        error.format(name)
        raise ValueError(error)


def mad(x):
    '''Median absolute deviation
    Returns:
        MAD (float) = the median absolute deviation of the values in X
        diff (list) = a list of deviations of the values in X from the median of X'''
    median = np.median(x)
    diff = np.absolute(x - median)
    MAD = np.median(diff)
    return (MAD, diff)


''' Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
    Handle Outliers", The ASQC Basic References in Quality Control:
    Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.'''


def find_outliers(x, threshold=2):
    '''Outlier detection method using Iglewicz & Hoaglin's modified Z-score
    Returns a list of bools, where True is an outlier'''
    (MAD, diff) = mad(x)
    M = 0.6745 * diff / MAD
    outliers = M > threshold
    return outliers


def empirical_cdf(data, bins, rangemin=None, rangemax=None):
    '''Calculate the empirical cumulative distribution function for data
    Inputs:
        data (array) = data to calculate ecdf of
        bins (int/array) = same behavior as numpy.histogram (int defines
                           number of bins, array specifies bin edges)
        rangemin/rangemax (float) = must be specified if bins is of type=int,
                                    and rangemin < rangemax
    Returns:
        ecdf (array) = empirical cumulative distribution function
        edges (array) = edges used during binning of histogram'''

    if isinstance(bins, np.ndarray):
        edges = bins
    else:
        message = "bins must either be an int or a list corresponding to the edges"
        assert any(isinstance(bins, t) for t in [int, list]), message

        if isinstance(bins, int):
            assert bins > 1, "number of bins must be greater than 1."
            assert rangemin < rangemax, "rangemin must be less than rangemax"
            assert rangemin is not None and rangemax is not None, ("rangemin and "
                                                                   "rangemax must be specified if bins is an integer.")

            edges = np.linspace(rangemin, rangemax, bins+1)

        elif isinstance(bins, list):
            edges = np.array(bins)

    # print("data", data)
    # print("edges", edges)
    h = np.histogram(data, edges)
    # print(h[0])
    ecdf = np.cumsum(h[0], dtype=np.float64)/len(data)

    return ecdf, edges


def fit_linear_model(x, y, slope=None):
    '''Linear least squares (LSQ) linear regression (polynomial fit of degree=1)
    Returns:
        m (float) = slope of linear regression line of form (y = m*x + b)
        b (float) = intercept of linear regression line'''

    assert len(x) == len(y), ("Arrays x & Y must be equal length to fit "
                              "linear regression model.")
    if slope == None:
        (m, b) = scipy.polyfit(x, y, deg=1)
    else:
        def LSQ(b): return np.sum((y-(slope*x+b))**2.0)
        res = scipy.optimize.minimize(LSQ, x0=1, bounds=None)
        (m, b) = (slope, res.x[0])
    return (m, b)


def linear_complete(x_values, prot_mean, prot_std, linear_fit, x_scale='linear', y_scale='ln', slope=None, intercept=None):
    '''
    linear_complete returns all the statistics needed to analyse a RBS model.

    The usage of this function differs if your model is a free energy model or
    not.

    If your RBS model is a free-energy model (i.e. you compute a final total
    free energy value and then use that value to figure out the translation
    initiation rate), then you should set the `linear_fit` param is set to
    `True`.

    When the `linear_fit` param is set to `True`, this function expects that
    the `x_values` are the total free energy values. In this case, a line
    is fitted to `x_values` and `prot_mean` and the equation of the line is
    used to calculate the `y_predicted` (translation initiation rate) values.

    You should then update your `TranslationInitiationRate` function to use
    the outputted slope and intercept values to calculate translation
    initiation rate from your computed free energy values.

    Optionally, when `linear_fit` is set to `True`, you can pass the slope
    and/or the intercept values that will be used to calculate the translation
    initiation rate values (`y_predicted`) from your total free energy values
    (`x_values`).


    If your RBS model is not a free-energy model then the `x_values` should be
    the translation initiation rates. In this case, you should set `linear_fit`
    param is set to `False`. In this case, no line is fitted and only the
    statistics for your model are calculated.
    '''
    # Useful lambda functions
    def calc_x(a0, a1, y): return (y-a0)/a1
    def calc_y(a0, a1, x): return a1*x+a0

    if y_scale == 'log10':
        prot_mean = np.log10(prot_mean)
    elif y_scale == 'ln':
        prot_mean = np.log(prot_mean)
    elif y_scale == 'linear':
        pass
    else:
        raise ValueError(
            "Invalid input in ModelTest._linear_model_stats for yScale: {}".format(y_scale))

    if x_scale == 'log10':
        x_values = np.log10(x_values)
    elif x_scale == 'ln':
        x_values = np.log(x_values)
    elif x_scale == 'linear':
        pass
    else:
        raise ValueError(
            "Invalid input in ModelTest._linear_model_stats for xScale: {}".format(x_scale))

    outliers = 0
    if linear_fit:
        a1, a0 = 0, 0
        if slope != None and intercept != None:
            a1 = slope
            a0 = intercept
        else:
            # determine outliers with initial fit
            (a1, a0) = fit_linear_model(x_values, prot_mean, slope=slope)
            app_xVals = calc_x(a0, a1, prot_mean)
            abs_delta_x = np.absolute(x_values - app_xVals)
            outliers = find_outliers(abs_delta_x)
            keepers = np.invert(outliers)

            # reapply fit with trimmed dataset & calculate error
            xVals1 = x_values[keepers]
            yVals1 = prot_mean[keepers]
            (a1, a0) = fit_linear_model(xVals1, yVals1, slope=slope)

        y_predicted = calc_y(a0, a1, x_values)
    else:
        y_predicted = x_values

    residuals = prot_mean - y_predicted

    if y_scale == 'log10':
        yError = 10**(residuals)
    elif y_scale == 'ln':
        yError = np.exp(residuals)
    else:
        yError = prot_mean/y_predicted

    # Pearson/Spearman correlation coefficients
    (R, Pearson_p) = correlation(y_predicted, prot_mean)
    (rho, Spearman_p) = correlation(y_predicted, prot_mean)

    # Root Mean Square Error (RMSE)
    # n = len(yVals)
    # SSE = np.sum(res**2.0 for res in residuals)
    # RMSE = np.sqrt(SSE/n)
    RMSE = np.sqrt(1-R**2.0)*np.std(prot_mean)

    # One-sided model error cdfs
    yError1 = np.array([1/val if val < 1 else val for val in yError])
    bins = np.concatenate((np.linspace(1, 10, 10), np.linspace(
        20, 100, 9), np.linspace(200, 1000, 9)))
    # print(yError1)
    # print(bins)
    onesided_cdf, _ = empirical_cdf(yError1, bins)
    # print(onesided_cdf)

    # Kullback-Leibler divergence
    (NKLdiv, KLdiv, KLdivmax) = normKLdiv(yError, b=4)

    # AUC ROC
    threshold = (max(prot_mean) + min(prot_mean))/2.0
    auroc, fpr, tpr, thresholds = area_under_ROC_curve(
        y_predicted, prot_mean, cutoff=threshold)

    # Relative information gain (RIG) over uniform model
    # remove nans from ystd
    # if len(ystd)==0, then skip information theory analysis
    prot_std = np.array(prot_std, dtype=np.float64)
    nans = np.isnan(prot_std)
    prot_std = prot_std[~nans]
    prot_mean = prot_mean[~nans]

    if len(prot_std) == 0:
        CV = 0.0
        N = 0.0
        RIG = 0.0
    else:
        if y_scale == 'log10':
            prot_std = np.log10(prot_std)
        elif y_scale == 'ln':
            prot_std = np.log(prot_std)
        else:
            pass

        nonzero = np.nonzero(prot_std*prot_mean)
        prot_std = prot_std[nonzero]
        prot_mean = prot_mean[nonzero]

        negative = (prot_mean < 0.0) + (prot_std < 0.0)
        CV = np.mean(prot_std[~negative]/prot_mean[~negative])
        N = int(np.floor((max(prot_mean) - min(prot_mean))/CV))
        edges = np.linspace(-4, 4, N+1)
        # filter out residuals that don't fall in the bins
        rmv = (residuals < edges[0]) + (residuals > edges[-1])
        residuals = residuals[~rmv]

        [dist, _] = np.histogram(residuals, edges, density=True)
        dist_model = dist[dist != 0.0]
        Hmodel = entropy(dist_model)
        dist = np.ones(N)/N
        Hrandom = entropy(dist)
        RIG = 1 - Hmodel/Hrandom
    # WARNING yVals and y_predicted have been modified after MC calculations

    results = {
        "Pearson R-squared": R**2.0,
        "Pearson p": Pearson_p,
        "Spearman R-squared": rho**2.0,
        "Spearman p": Spearman_p,
        "N-outliers": np.sum(outliers),
        "slope": a1,
        "intercept": a0,
        "2-fold Error": onesided_cdf[1],
        "5-fold Error": onesided_cdf[4],
        "10-fold Error": onesided_cdf[9],
        "Normalized KL-divergence": NKLdiv,
        "KL-divergence": KLdiv,
        "Max KL-divergence": KLdivmax,
        "RMSE": RMSE,
        "RIG": RIG,
        "N-states": N,
        "CV": CV,
        "AUROC": auroc
    }

    return results, yError


def compute_stats(csv_file, x_values_col_idx, prot_mean_col_idx, prot_std_col_idx, linear_fit, slope=None, intercept=None, use_filter=False, filter_col_idx=None, filter_values=None):
    '''
    compute_stats computes the statistics needed to compare RBS models.

    Please see the documentation of `linear_complete` for more details on the
    `x_values_col_idx`, `linear_fit`, `slope` and `intercept` arguments.

    If `use_filter` is set to `True`, the function expects that the
    `filter_col_idx` and `filter_values` arguments are passed in.
    The `filter_values` argument should be a list of string of values available
    in `filter_col_idx` in the dataset.

    '''

    # read columns of csv file into the relevant variables
    if use_filter:
        if filter_col_idx == None or filter_values == None:
            raise ValueError("need `filter_col_idx` and `filter_values` arugments if `use_filter` is `True")

        filter_labels = np.genfromtxt(csv_file, delimiter=',',
                                    usecols=(filter_col_idx), dtype=str,
                                    skip_header=1)

    # print(np.unique(filter_labels))
    x_values, prot_mean, prot_std = np.genfromtxt(csv_file, delimiter=',',
        skip_header=1,
        usecols=(x_values_col_idx, prot_mean_col_idx, prot_std_col_idx),
        unpack=True)

    # Variable to keep track of the rows used from all datasets
    if use_filter:
        all_idxs = np.zeros(len(x_values), dtype=int)
    else:
        all_idxs = np.ones(len(x_values), dtype=int)

    all_stats_dict = []
    if use_filter:
        # calculate the
        for filter_value in filter_values:
            curr_filter_idxs = np.array(filter_labels == filter_value)
            all_idxs += curr_filter_idxs
            stats_dict = do_compute_stats(x_values, prot_mean, prot_std,
                curr_filter_idxs, filter_value, linear_fit, slope, intercept)
            all_stats_dict.append(stats_dict)

    # convert all_idxs from a integer array to boolean array
    all_idxs = np.array(all_idxs == 1)
    if use_filter:
        filter_value = " + ".join(filter_values)
    else:
        filter_value = "all"

    stats_dict = do_compute_stats(x_values, prot_mean, prot_std, all_idxs,
        filter_value, linear_fit, slope, intercept)
    all_stats_dict.append(stats_dict)

    pd.options.display.max_columns = None
    df = pd.DataFrame(all_stats_dict)
    print(tabulate(df,headers='keys',tablefmt='psql'))

    # only return the all stats dict
    return stats_dict


def do_compute_stats(total_free_energy, prot_mean, prot_std, idxs, filter_value, linear_fit, slope, intercept):
    print('processing {}'.format(filter_value))
    curr_dataset_total_free_energy = total_free_energy[idxs]
    curr_dataset_prot_mean = prot_mean[idxs]
    curr_dataset_prot_std = prot_std[idxs]

    # filter out no-prediction sequences (nan and inf values)
    invalid = np.isnan(curr_dataset_total_free_energy) + \
        np.isinf(curr_dataset_total_free_energy) + \
        (curr_dataset_total_free_energy == 0.0)
    valid_total_free_energy = curr_dataset_total_free_energy[~invalid]
    valid_prot_mean = curr_dataset_prot_mean[~invalid]
    valid_prot_std = curr_dataset_prot_std[~invalid]
    count_invalid = np.sum(invalid)

    dataset_stats, _yError = linear_complete(valid_total_free_energy,
                                             valid_prot_mean, valid_prot_std, linear_fit, slope=slope, intercept=intercept)
    dataset_stats["filter"] = filter_value
    dataset_stats["size"] = len(valid_total_free_energy)
    dataset_stats["invalid count"] = count_invalid
    slope, intercept = dataset_stats["slope"], dataset_stats["intercept"]
    dataset_stats["slope / intercept string"] = 'slope, intercept := {}, {}\n'.format(slope, intercept)

    return dataset_stats



if __name__ == "__main__":
    dataset_with_properties = "./salis_lab_v2_1/dataset_with_properties/train_1629043294.csv"

    dataset_idx = 7
    dataset_filter_values = ["Kosuri_PNAS_2013", "Goodman_Science_2013"]
    # organism_col_idx = 11
    # organism_filter_values = [
    #     "Escherichia coli BL21(DE3)", "Escherichia coli str. K-12 substr. DH10B", "Escherichia coli str. K-12 substr. MG1655"]
    # protein_col_idx = 5
    # protein_filter_values = ["sfGFP", "RFP", "RFP-GFP"]
    total_free_energy_idx, prot_mean_idx, prot_std_idx = 44, 4, 5

    stats_dict = compute_stats(dataset_with_properties, total_free_energy_idx,
        prot_mean_idx, prot_std_idx, linear_fit=True, use_filter=True,
        filter_col_idx=dataset_idx, filter_values=dataset_filter_values)
    # slope, intercept = stats_dict["slope"], stats_dict["intercept"]
    slope = -0.45
    # intercept = 7.11720550316
    stats_dict = compute_stats(dataset_with_properties, total_free_energy_idx,
        prot_mean_idx, prot_std_idx, linear_fit=True, slope=slope,
        use_filter=False, filter_col_idx=dataset_idx,
        filter_values=dataset_filter_values)
