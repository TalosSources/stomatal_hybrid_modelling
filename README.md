yThis repository contains code for learning a hybrid model of stomatal resistance (r_s) within (Tethys-Chloris) T&C, using a separate python pipeline. Additionnaly, due to a mismatch in the python pipeline and T&C's equations, a calibration procedure is offered, which iteratively run T&C predictions, then learns an r_s model using those predictions, and then predict again using the new r_s model, and so on until convergence of the predictions.

# Repository structure
## configs
This contains .yaml files describing a training run to perform in the python pipeline. Config ```iterative_training.yaml``` should be used to run the iterative calibration procedure.

## results
Contains the output of a python training run. A result folder contains running time, losses as a numpy array, and the model output weights that can be used for inference.

## python
Contains the source code used in the python training pipeline, responsible for loading the data, specifing the pipeline and the models, training and evaluating them, and optionnally plotting the results.

## tandc
Contains the source code, inputs and outputs related to T&C.
### tandc-calibrated-parameters
Contains the parameter files necessary to run T&C on many observation sites.
### tandc-forcing
Should contain the observation data in ```.mat``` format, including "Data_..." and "Res_..." files.
### tandc-physics
Contains the source code, running scripts, and files specifying sites to run.
### tandc-physics-output
Contains the output files from T&C, used by the python pipeline.

## traced-models
Contains r_s models in the ```onnx``` format, readable by matlab.

# How to run it
Prerequisites must be installed, preferably in a python environment like conda.
They are present in the ```requirements.txt``` file.
The environment can be created with these commands:

```conda create --name stomatal```

```conda activate stomatal```

```conda install pip```

```pip install -r requirements.txt```


## Site choice
To choose the sites to use in the calibration, the file ```tandc/tandc-physics/site_tandc_fluxnet2015_ameriflux_final_v2_iterative_training_selection.csv``` should be modified. 

Each line corresponds to a site to run, and the lines can be copied from the file ```tandc/tandc-physics/site_tandc_fluxnet2015_ameriflux_final_v2.csv```.

## .yaml config
The configuration ```iterative_training.yaml``` defines (hyper)-parameters for the python pipeline. Deep model settings, learning rate, epochs, batch size etc. can be set there.

## Running calibration
Then, to run the iterative training script, the script ```iterative_training.sh``` should be run.
It will run either until convergence, or until ```max_iter``` iterations are performed, which is 5 by default.
Since at each iteration, the T&C model is run on potentially many observation sites on >100'000 timesteps, the script may take many hours to run.

After running, the final predictions can be found in ```tandc/tandc-physics-output/tandc_outputs_for_python``` and the r_s weights are saved under ```traced_models/iterative_training.onnx```, and will be used by default by our modified version of T&C if they're present at this path.

Then, T&C can be run as usual using the script ```tandc/tandc-physics/run_tandc_physics_fluxnet2015_ameriflux.sh```, replacing the ```file_site_list``` path with the ```.csv``` file containing the desired sites to run inference on. The predictions will be made with the learned r_s model.


# Modified T&C
In this project, we use a slightly modified version of the T&C sourcecode. Some small modifications are quality of life changes, but the main one is the use of the deep r_s model. In the file ```tandc/tandc-physics/tandc-model/src/photosynthesis_biochemical.m```, this is done starting at line 298, where the architecture and weights are loaded, if the weights path exists. Then, if that's the case, the model is used to compute r_s instead of the empirical one. This loading and inference code can be re-used to perform analyses and plots. The r_s model takes 7 predictors, which can be seen at line 315, and need to be rescaled first with specific scaling factors. Its output also need to be passed through the ```softplus``` function, and rescaled too. 

# Provided run
The repo comes with a site selection aiming to represent all possible combinations of cold/temperate/hot, arid/semi-arid/humid, and grassland/forest conditions. Since temperate/hot and arid forests don't really exist, 16 sites are used in total. The resulting traced models and outputs are available, with the large outputs being available on git LFS at path ```tandc/tandc-physics-output/tandc_outputs_for_python```.
