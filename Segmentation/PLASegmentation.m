function [ P, QRS, T ] = PLASegmentation( ecg_signal, sample_rate )
%%% Accepts one lead ECG and sample rate (in Hz) and outputs three structs:
%%% P contains P.onset, P.offset, and P.peaks as lists
%%% QRS contains QRS.onset, QRS.offset, QRS.Q, QRS.R, QRS.R', and QRS.S as
%%% lists
%%% T contains T.onset, T.offset, T.peaks, and T.biphasic as lists

addpath('QRS');
addpath('Utility');

%% Denoising and baseline calculation
ecg_filtered = Denoise(ecg_signal, sample_rate);
baseline = repelem(mean(ecg_filtered),length(ecg_filtered));

%% Calculate threshold for Piecewise Linear Approximation
threshold = max(abs(ecg_filtered - baseline)) / 10;
if threshold < 25
    threshold = 25;
end

%% Step 1: Extract significant points using the Piecewise Linear Approximation algorithm
extracted_points = Extract_Significant_Points(ecg_filtered, threshold);

%% Step 2: Detect candidate QRS peaks using the significant point representation of ECG
[qrs_extremums, qrs_second_extremums, fragmented_qrs] = Detect_QRS_Peaks(...
    ecg_filtered,...
    baseline,...
    extracted_points,...
    sample_rate);

%% Step 3: Detect QRS onset, offset, and individual waves
EE_interval = Mean_EE_Interval(qrs_extremums);
[qrs_left_sets, qrs_right_sets, Q_list, R_list, R2_list, S_list] = Detect_QRSComplex( ecg_filtered, qrs_extremums, fragmented_qrs, EE_interval);

%% Step 4: Detect T onset, offset, and individual waves
[ T_peaks, T_left_sets, T_right_sets, T_biphasic] = delineate_T_waves( ecg_filtered, qrs_left_sets, qrs_right_sets, EE_interval, sample_rate);
%% Step 5: Detect P onset, offset, and peaks
[ P_peaks, P_left_sets, P_right_sets] = delineate_P_waves( ecg_filtered, T_right_sets, qrs_left_sets, EE_interval, sample_rate);

%% Prepare the output
P.peak = P_peaks;
P.onset = P_left_sets;
P.offset = P_right_sets;

QRS.peak = qrs_extremums;
QRS.fragmented = fragmented_qrs;
QRS.onset = qrs_left_sets;
QRS.offset = qrs_right_sets;
QRS.Q = Q_list;
QRS.R = R_list;
QRS.R2 = R2_list;
QRS.S = S_list;


T.peak = T_peaks;
T.biphasic = T_biphasic;
T.onset = T_left_sets;
T.offset = T_right_sets;

end

