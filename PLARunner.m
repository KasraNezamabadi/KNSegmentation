clear;
clc;
close all;
warning off;

%% INCLUDES
addpath('LoaderFunctions');
addpath('Segmentation');
addpath('Segmentation/QRS');
addpath('Segmentation/T');
addpath('Segmentation/P');
addpath('Segmentation/Utility');

addpath('Data');
addpath('Data/ECG');

%% Main
meta_path = "Data/ECGMeta.xlsx";

meta_table = readtable(meta_path);

for i = 1:size(meta_table, 1)

    ecg_id = table2array(meta_table(i, 1))
    frequency = table2array(meta_table(i, 4));

    ecg_file_name = strcat(num2str(ecg_id), '.csv');
    ecg = getECG(ecg_file_name);
    ecg = ecg';
    P_anns = [];
    QRS_anns = [];
    T_anns = [];
    for lead = 1:12
        %% Automatic Segmentation
        signal = ecg(lead, :);
        [ P, QRS, T ] = PLASegmentation( signal, frequency);
        P_anns = [P_anns, P];
        QRS_anns = [QRS_anns, QRS];
        T_anns = [T_anns, T];
    end
    save(strcat('Data/PLAAnnotation/', strcat(strcat(strcat(num2str(ecg_id), '_'), 'P'), '.mat')), 'P_anns')
    save(strcat('Data/PLAAnnotation/', strcat(strcat(strcat(num2str(ecg_id), '_'), 'QRS'), '.mat')), 'QRS_anns')
    save(strcat('Data/PLAAnnotation/', strcat(strcat(strcat(num2str(ecg_id), '_'), 'T'), '.mat')), 'T_anns')
end


function [ecg] = getECG(target_file_name)

path_to_folder = 'Data/ECG/';

ecg_dir = dir(fullfile(path_to_folder,'*.csv'));

for i = 1:numel(ecg_dir)
    file = ecg_dir(i);
    file_name = file.name;
    
    if strcmp(file_name, target_file_name)
        path_to_ecg = strcat(strcat(path_to_folder, '/'), file_name);
        ecg = csvread(path_to_ecg, 1);
        break
    end 
end

end





