%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function to fill up the mod_param template
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = tmp2param(station, param_file, param_template, output_dir)
    %% read in the parameter file
    param_data = readtable(param_file);

    %% extract the parameter values for the required station
    %row_index = find(param_data.SITE_ID == station);
    row_index = find(ismember(param_data.SITE_ID, station));

    %% subset the parameter file,
    param_data = param_data(row_index,:);

    %% remove the columns which are not needed (currently hard-coded)
    param_data = removevars(param_data, {'SITE_ID', ...
                                         'SITE_NAME' ...
                                         'LOCATION_LAT', ...
                                         'LOCATION_LONG', ...
                                         'LOCATION_ELEV', ...
                                         'IGBP', ...
                                         'DeltaGMT', ...
                                         'HorL'});
    
    %% read in the model template file
    mod_param = fileread(param_template);
    
    %% 
    param_names = param_data.Properties.VariableNames;
    param_names = param_names';

    %%
    for i = 1:length(param_names)
        % extract data 
        param_data_temp = string(table2array(param_data(1, [param_names{i}])));
        param_name_temp = param_names{i};
        param_name_temp = param_name_temp(2:end);
        %disp(join([param_name_temp, param_data_temp], " "))
        mod_param = strrep(mod_param, param_name_temp, param_data_temp);
    end

    %% write out the file
    output_file = join([output_dir,'/mod_param','.m'], "");
    fid = fopen(output_file, 'wt');
    fwrite(fid, mod_param);
    fclose(fid);

end
