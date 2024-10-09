# Random Forest
"Two separate random forest regressors were trained using the
combinations of the plant properties (functional types, LAI,
LCC) and environmental conditions (meteorological vari-
ables, soil types, soil water content) with the optimized m and
V 25cmax, respectively. The plant functional types and soil types
were encoded with the one-hot encoder. The meteorological
variables and soil water content were collected according to
the time and location of each m and V 25cmax retrieval. Overall,
80 % of the data were randomly selected to train the random
forest regressors with a 5-fold cross-validation to determine
the hyperparameters, and 20 % of the data were used to eval-
uate the performance of the trained regressors in each round
of calibration." (from the [...]hourly[...] paper) 
If I understand correctly what they do in this papers: they have sparse data in space (on specific sites) and in time (weekly/monthly) for the value of VCMax (maximum rubisco capacity at 25°C), and they want to have estimations at the grid level. For the spatial grid, they train a random forest regressor to predict it.