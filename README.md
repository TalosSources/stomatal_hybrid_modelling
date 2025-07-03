This is the code for the semester project titled *"Hybrid modeling of stomatal resistance using a complex ecohydrological model and eddy covariance measurements"*.

# Repository structure
## configs
This contains .yaml files describing a training run to perform. Config ```default.yaml``` can be used as a baseline to create new configurations. Config ```best_model.yaml``` trains the best performing model architecture with all training sites.

## results
Contains the output of running all the training configs discussed in the report. A result folder contains running time, losses as a numpy array, and the model output weights that can be used for inference.

## python
This contains the source code used in the entire project, responsible for loading the data, specifing the pipeline and the models, training and evaluating them, and plotting the results.

# How to run it
Prerequisites must be installed, preferably in a python environment like conda:
* pytorch
* numpy
* wandb
* scikit-learn
* scipy
* matplotlib
* omegaconf
* torchviz

To run a configuration, e.g. ```best_model.yaml```, from the root of the repository, use command:
```python python/main.py best_model.yaml```

# TODO:
* create conda environment using requirements
* where the T&C files are (need to download?)
* say that they need to run once T&C with the base model ??
* explain how to choose the sites
* explain how to run training (+ some infos about perf?)
* explain again the .yaml (and some other stuff from python)
