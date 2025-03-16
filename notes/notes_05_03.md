# Idea to obtain a rs model compatible with T&C
* train rs model using usual pipeline
* plug-it in T&C
* collect new data from T&C using this model
* run the usual pipeline with this new data
* repeat this until convergence (model stays the same? loss stops decreasing?)

# clean way to do it
ideally, the best way to do it would be to translate all the modules between pb.m and the computation of Q_LE. if it's not too much / too complex, I could consider translating it to python. Or, find a way to have torch-like features inside of matlab.

# Investigate if Julia is easy to 
low priority side-project.