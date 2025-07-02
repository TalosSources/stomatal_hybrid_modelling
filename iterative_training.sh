#! /bin/bash

# This file contains a bash script that iteratively trains a stomatal resistance model, 
# and use it to run the T&C model and generate new training data until convergence.

# TODO: Make paths general?

# git_dir="/home/talos/git_epfl/stomatal_hybrid_modelling/"
git_dir=$PWD"/"
t_and_c_dir=$git_dir"tandc/"
t_and_c_physics_dir=$t_and_c_dir"tandc-physics/"
results_dir=$t_and_c_dir"/tandc-physics-output/tandc_outputs_for_python/"
comparison_path=$git_dir"comparison_temp_dir/"
mkdir $comparison_path -p

comparison_sites_path=$t_and_c_physics_dir"/site_tandc_fluxnet2015_ameriflux_final_v2_iterative_training_selection.csv"

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
    python python/main.py "configs/iterative_training.yaml"
    [ $? -eq  0 ] || exit 1 # Exit if the python training failed

    # copy subset of old predictions elsewhere for future comparison
    while IFS=, read -r site_id sitename mod_param_name prova_name lat lon elev igbp
    do
        site_file_name="Results_"${site_id}".mat"
        fresh_output_site_path=$results_dir$site_file_name
        cp $fresh_output_site_path $comparison_path # TODO: rearrange file structure 
    done < <(tail -n +2 ${comparison_sites_path})

    # run T&C to generate predictions
    # cd $t_and_c_physics_dir
    # ./run_tandc_physics_fluxnet2015_ameriflux_debug.sh
    $t_and_c_physics_dir"run_tandc_physics_fluxnet2015_ameriflux_debug.sh" $git_dir
    [ $? -eq  0 ] || exit 1 # Exit if the matlab predictions failed

    # check stability condition
    cd $git_dir
    converged=$(python check_convergence.py $results_dir $comparison_path $comparison_sites_path)
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

