To correct for a bias in the training

# Plan
* train the model as usual using data from vanilla T&C
* repeat until convergence:
    * plug the trained model in T&C
    * generate new data by running T&C on all training sites
    * train the model using this new data

Can this be automated?
With smart handling of model and data paths, a bash script could easily iteratively start training and prediction rounds.
The only tricky part would be checking convergence. This would require saving the new data/model elsewhere than the old, so that both can co-exist and a comparison can be made. Or the comparison can be made before saving in memory. In any case, it would make sense to check only a small subset of the entire data for convergence, which would be equivalent to checking everything if the subset is representative and would save a lot of time/space. 
And to perform the actual check, it would be more convenient to do it in python, and output a flag that's readable by bash.

## Convergence
Can be defined in terms of model weights convergence, or predictions convergence. It seems the former always has significant stochasticity even when training goes well, so latter might make more sense.

# Roadmap
* worry somewhat about some of the magic numbers (including convergence threshold)

## Then
* hyperparams again?

## Goal
* select sites
    * try, for temperature and humidity, to get all combinations of low-moderate-high for each. Yields 9 sites. Then if possible, try to get for each combination, a grassland and a forest. might not be possible (no forest in hot arid places :( ). At maximum, we have 18 sites to find.
* clean the code and repo for Son DONE? need feedback from Akash, Sara and Son
    * an inference script? for T&C, that's quite automatic
    * refactor -> change/fix the t_and_c path, make sure everything is a parameter
    * write README, check requirements.txt is complete
    * readme should be usable by Sara and Akash as well, so that they can make modifications, change sites etc., make plots and analysis.
* check the installation and script works for a fresh system 
* run the iterative training calibration with the chosen sites

# Perf
## Knobs
* python: data size, epochs, batch size, network size
* matlab: iters, perhaps also data size
* sites count
## modified
* removed the python part
* changed iters and data size in python
* changed iters in matlab

# Problems
* no convergence -> consider using r_s predictions instead of An, or something else

# General notes
the script assumes that T&C was run before, with the modified final parameters. Perhaps, that should be mentioned somewhere, or the script should (offer the option to) run it before running python

# sites
* grass-cold-arid: CN-Dan (semiarid?), GL-ZaH
* grass-cold-semiarid: CN-Du2, CN-Cng, CN-Sw2
* grass-cold-humid: CH-Fru, IT-Mbo, IT-Tor, RU-Ha1 (semiarid?)
* grass-temperate-arid: AU-Asm (warm, and perhaps forest)
* grass-temperate-semiarid: IT-Noe (forest? humid?), US-Twt
* grass-temperate-humid: CH-Cha/CH-Oe1/CH-Oe2 (cold?)
* grass-hot-arid: SD-Dem (semiarid? forest?)
* grass-hot-semiarid: CG-Tch (forest? humid?), SN-Dhr (arid?) AU-Dry/Stp/TTE?? (wet&dry)
* grass-hot-humid: AU-DaP (tropical savanna, permanent wetland?), AU-DaS?
* forest-cold-arid: US-Prr?
* forest-cold-semiarid: FI-Sod? US-Ivo? US-Me2?
* forest-cold-humid: BE-Vie, CA-Qfo, CH-Dav, CH-Lae (not best), other DE/DK, IT-Col, IT-Lav, IT-Ren
* forest-temperate-arid: 
* forest-temperate-semiarid: ES-Amo?? (grass? arid?)
* forest-temperate-humid: AU-Wac, BE-Bra (semiarid?), CN-Din/CN-Qia (warm), IT-Cpz, IT-Ro2, IT-Sro, other DE
* forest-hot-arid: probably impossible
* forest-hot-semiarid: perhaps AU-Ade?
* forest-hot-humid: GF-Guy, GH-Ank, PA-SPn, AU-DaS? AU-How? (those 2 might be semi-arid/savannas), ZM-Mon

* wetlands: GL stuff, US-Ivo, US-Tw1 (good), SJ-Adv

* NO INFO: ES-Amo, RU-Vrk, SJ-Adv, AR-Vir, AU-Ade, AU-Asm, AU-Dry, AU-Stp, AU-TTE, CN-Cng, CN-Ha2, CN-HaM, CN-Sw2, JP-MBF

* needed and don't have:


All humid combinations are covered, I have issues for arid/semiarid. Some sites are not clearly assigned (shrubs? e.g. ES-Amo), also arid/semiarid depends on temperature? I don't really view a 300mm in Canada the same as in Sudan.