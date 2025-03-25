%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function to fill up the prova template
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = tmp2prova(station, prova_template, input_data_dir, param_file, ...
                        ca_file, model_src_dir, output_data_dir, output_prova_dir)
    %% Read in the prova template file
    prova = fileread(prova_template);

    %% Replace the required strings in the template
    % input data file
    input_data_file = join([input_data_dir, '/', 'Data_', station, '_run.mat'],"");
    prova = strrep(prova, '_inputdata_', input_data_file);

    % carbon dioxide file
    prova = strrep(prova, '_cadata_', ca_file);

    % model parameter file
    prova = strrep(prova, '_modparam_', param_file);

    % source code directory (which contains MAIN_FRAME)
    prova = strrep(prova, '_srccode_', model_src_dir);

    % output data file
    output_data_file = join([output_data_dir, '/', 'Results_', station, '.mat'], "");
    prova = strrep(prova, '_outputfile_', output_data_file);

    %% Write out the prova file
    output_prova_file = join([output_prova_dir, '/', 'prova.m'], "");
    fid = fopen(output_prova_file, 'wt');
    fwrite(fid, prova);
    fclose(fid);

end