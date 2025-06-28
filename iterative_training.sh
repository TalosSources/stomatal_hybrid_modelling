#! /bin/bash

# This file contains a bash script that iteratively trains a stomatal resistance model, 
# and use it to run the T&C model and generate new training data until convergence.

# TODO: Make paths general?

git_dir="/home/talos/git_epfl/stomatal_hybrid_modelling/"
t_and_c_dir=$git_dir"t_and_c_alternative/"
t_and_c_physics_dir=$t_and_c_dir"tandc-physics/"
results_dir=$t_and_c_dir"/tandc-physics-output/fatichi_parameters_ismail/"
comparison_path=$git_dir"comparison_temp_dir/"
mkdir $comparison_path -p

comparison_site="CH-Dav" # TODO: choose several sites? choose a few time-steps for each site?

# 2 options for the loop: train then predict, or vice-versa. start simple
# for i in $(seq 1 10); # inclusive range
echo "Starting iterative training!"

max_iter=5 # stop iterative training after some time if no convergence 
converged="false"
i=0
while [ $i -lt $max_iter -a "$converged" = "false" ];
do

    i=$((i+1))
    echo "iteration: "$i

    # train an rs model
    cd $git_dir
    python python/main.py
    [ $? -eq  0 ] || exit 1 # Exit if the python training failed

    # copy subset of old predictions elsewhere for future comparison
    site_file_name="Results_"$comparison_site".mat"
    fresh_output_site_path=$results_dir$site_file_name
    cp $fresh_output_site_path $comparison_path # TODO: rearrange file structure 

    # run T&C to generate predictions
    cd $t_and_c_physics_dir
    ./run_tandc_physics_fluxnet2015_ameriflux_debug.sh
    [ $? -eq  0 ] || exit 1 # Exit if the matlab predictions failed

    # check stability condition
    cd $git_dir
    converged=$(python check_convergence.py $fresh_output_site_path $comparison_path$site_file_name)
    echo "converged="$converged
done

if [ "$converged" = "true" ]; then
    echo "calibration converged! exiting..."
else
    echo "calibration did not converge after "$i" iterations. exiting..."
fi

# cleanup files
cd $git_dir
rm -r $comparison_path

