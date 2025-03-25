#!/bin/bash
# Description: Script to run the T&C model for FLUXNET2015 and Ameriflux
#              stations.
#              NOTE: This is the serial version. The computational time
#                    is not very high and this suffices for now. 
# Author: Akash Koppa
# Date: 29.05.2024
#
#
#set -e

# user defined configuration
root_directory="/home/talos/git_epfl/stomatal_hybrid_modelling/t_and_c_alternative" # user-specific path: should be specified outside?
physics_directory="$root_directory/tandc-physics"
directory_with_model="$physics_directory/tandc-model"
directory_with_src="$directory_with_model/src"
param_directory="$root_directory/tandc-calibrated-parameters/fatichi_final"
directory_with_mod_param="$param_directory/mod_param_files"
directory_with_prova="$param_directory/prova_files"
directory_with_input_mat="$root_directory/tandc-forcing"
directory_with_output_mat="$root_directory/tandc-physics-output/fatichi_parameters_ismail"
file_ca="$root_directory/CO2_Data/Ca_Data.mat"
file_site_list="$physics_directory/site_tandc_fluxnet2015_ameriflux_final_v2.csv"

while IFS=, read -r siteid sitename mod_param_name prova_name lat lon elev igbp
do

    cd ${physics_directory}

    echo ">>>> Station Under Process: ${siteid} <<<<"

    # construct the name of the input file
    file_input_mat="Data_${siteid}_run.mat"
    path_input_mat="${directory_with_input_mat}/${file_input_mat}"

    # construct the name of the output file
    file_output_mat="Results_${siteid}.mat"
    path_output_mat="${directory_with_output_mat}/${file_output_mat}"

    # construct the name of the mod_param file
    file_mod_param="MOD_PARAM_${mod_param_name}.m"
    path_mod_param="${directory_with_mod_param}/${file_mod_param}"
    echo "${path_mod_param}"
    
    # construct the name of the prova file
    file_prova="prova_${mod_param_name}.m"
    path_prova="${directory_with_prova}/${file_prova}"

    # copy the files to the main directory
    cp ${path_mod_param} "${directory_with_model}/mod_param.m"
    cp ${path_prova} "${directory_with_model}/prova.m"

    # modify the prova file to incorporate all the relevant file paths
    # comment out all lines starting with 'load' so that they can be 
    # replaced by actual paths
    echo "---- commenting out load statements ----"
    #LC_ALL=C sed -i '' "/[[:<:]]load[[:>:]]/s/^/%/g" "${directory_with_model}/prova.m"
    LC_ALL=C sed -i "/\bload\b/s/^/%/" "${directory_with_model}/prova.m"

    # proposition from GPT: to make it cross-platform
    # if [[ "$OSTYPE" == "darwin"* ]]; then
    #     # macOS (BSD sed)
    #     LC_ALL=C sed -i '' "/\<load\>/s/^/%/" "${directory_with_model}/prova.m"
    # else
    #     # Linux (GNU sed)
    #     LC_ALL=C sed -i "/\<load\>/s/^/%/" "${directory_with_model}/prova.m"
    # fi

    # input data
    echo "----- replacing input path ----"
    echo -e "load('${path_input_mat}')\n$(cat ${directory_with_model}/prova.m)" > ${directory_with_model}/prova.m

    # CO2 data
    echo "---- replacing CO2 data path ----"
    echo -e "load('${file_ca}')\n$(cat ${directory_with_model}/prova.m)" > ${directory_with_model}/prova.m
    
    # mod_param file
    echo "---- replacing the MOD_PARAM file path ----"
    LC_ALL=C sed -i.bak -e "/PARAM_IC/ s/^/%/" -e "/PARAM_IC/ a\\
PARAM_IC='${directory_with_model//\//\\/}\/mod_param.m'" ""${directory_with_model}/prova.m""

    # source code directory
    echo "---- replacing the source directory path ----"
    LC_ALL=C sed -i.bak -e "/Directory=/ s/^/%/" -e "/Directory=/ a\\
Directory='${directory_with_src//\//\\/}'" ""${directory_with_model}/prova.m""
    LC_ALL=C sed -i.bak -e "/Directory =/ s/^/%/" -e "/Directory =/ a\\
Directory ='${directory_with_src//\//\\/}'" ""${directory_with_model}/prova.m""

    # output file
    echo "---- replacing the output file path ----"
    LC_ALL=C sed -i.bak -e "/save(/ s/^/%/" -e "/save(/ a\\
save('${path_output_mat//\//\\/}')" ""${directory_with_model}/prova.m""

    # remove the end of line characters
    #LC_ALL=C sed -i.bak 's/\r$//' "${directory_with_model}/prova.m"
    LC_ALL=C sed -i.bak '/^%/d' "${directory_with_model}/prova.m"

    # run the tandc model (using the prova file)
    echo "---- running the T&C model ----"
    cd ${directory_with_model}
    matlab -batch "prova"

    cd ${physics_directory}

    # remove temporary mod_param and prova files
    rm ${directory_with_model}/mod_param.m
    rm ${directory_with_model}/prova.m

done < <(tail -n +2 ${file_site_list})

# remove temporary files
echo "---- removing temporary files ----"
rm -rf "${directory_with_model}/."* 2> /dev/null
