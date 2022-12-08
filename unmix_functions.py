from panel_classes import *
from sklearn import linear_model
from sklearn.metrics import mean_squared_error


def unmix(measurement, spectra, method, **kwargs):
    if method == "LassoNN":
        reg = linear_model.Lasso( positive=True,\
                                 alpha = kwargs["alpha"], max_iter=kwargs["max_iter"])
        reg.fit(spectra.values, measurement)
        coef = reg.coef_
    elif method == "Lasso":
        reg = linear_model.Lasso( alpha = kwargs["alpha"]\
                                , max_iter=kwargs["max_iter"])
        reg.fit(spectra.values, measurement)
        coef = reg.coef_
    elif method == "OLS":
        reg = linear_model.LinearRegression()
        reg.fit(spectra.values, measurement)
        coef = reg.coef_
    elif method == "OLSNN":
        coef = nnls(spectra, measurement)[0]
    elif method == "Ridge":
        reg = linear_model.Ridge( alpha = kwargs["alpha"]\
                                , max_iter=kwargs["max_iter"])
        reg.fit(spectra.values, measurement)
        coef = reg.coef_
    color_dict ={}
    for i, fluor in enumerate(spectra.columns):
        color_dict[fluor] = coef[i]
    return color_dict
    
def recall(df_result, fluor, sample):
    positive_population = set(df_result.sort_values(fluor, ascending = False)["id"]\
                              [:int(SAMPLE_SIZE*POSITIVE_FRACTION)].index.to_list())
    match_size = len(set(sample.components[fluor]).intersection(set(positive_population)))
    recall = match_size/(SAMPLE_SIZE*POSITIVE_FRACTION)
    return recall