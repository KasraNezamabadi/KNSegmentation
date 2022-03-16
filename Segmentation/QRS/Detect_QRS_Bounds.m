function [left_refractory_peak, onset, offset, right_refractory_peak, Q, R, S] = Detect_QRS_Bounds( signal, peak_index, avg_EE_interval)

left_refractory_peak = [];
onset = [];
offset = [];
right_refractory_peak = [];

%% Constants and Thresholds
threshold_baseline = 20; % Threshold 10 mV to be part of the baseline.

blank_period_refractory = round(0.12 * avg_EE_interval);
blank_period_delineation = round(0.12 * avg_EE_interval);
gap = 1;
smooth_window = 0.6;

%% Delineate QRS peak
[left_refractory_peak_candidate, ~, ~, right_refractory_peak_candidate, ~, ~] = delineate_wave( signal, peak_index, 1, avg_EE_interval, [], 0);

% plot(signal); hold on;
% plot(peak_index, signal(peak_index), 'o'); hold on;
% plot(left_refractory_peak_candidate, signal(left_refractory_peak_candidate), 'o'); hold on;
% plot(right_refractory_peak_candidate, signal(right_refractory_peak_candidate), 'o'); hold on;
% close all;

%% Establish smoothed Left Segment
lower_bound_left = peak_index - gap - blank_period_refractory;
if lower_bound_left < 1
    lower_bound_left = 1;
end
upper_bound = peak_index - gap;
left_segment_signal = signal(lower_bound_left : upper_bound);
left_segment_signal_smoothed = smooth(1 : length(left_segment_signal), left_segment_signal, smooth_window, 'loess');

left_refractory_peak_local = left_refractory_peak_candidate - (lower_bound_left - 1);

%% Establish smoothed Right Segment
lower_bound_right = peak_index + gap;
upper_bound = peak_index + gap + blank_period_refractory;
if upper_bound > length(signal)
    upper_bound = length(signal);
end

right_segment_signal = signal(lower_bound_right : upper_bound);
right_segment_signal_smoothed = smooth(1 : length(right_segment_signal), right_segment_signal, smooth_window, 'loess');

right_refractory_peak_local = right_refractory_peak_candidate - (lower_bound_right - 1);

% baseline_l = repelem(mean(signal),length(left_segment_signal));
% baseline_r = repelem(mean(signal),length(right_segment_signal));
% base_left = repelem(mean(left_segment_signal),length(left_segment_signal));
% base_right = repelem(mean(right_segment_signal),length(right_segment_signal));
% 
% close all;
% figure;
% plot(left_segment_signal); hold on;
% plot(left_refractory_peak_local, left_segment_signal(left_refractory_peak_local), 'o')
% plot(baseline_l); hold on;
% plot(base_left);
% 
% figure;
% plot(right_segment_signal); hold on;
% plot(right_refractory_peak_local, right_segment_signal(right_refractory_peak_local), 'o')
% plot(baseline_r); hold on;
% plot(base_right);
% 
% 
% temp = [left_segment_signal, right_segment_signal];
% baseline_temp = repelem(mean(temp),length(temp));
% baseline_temp_actual = repelem(mean(signal),length(temp));
% 
% figure;
% plot(temp); hold on;
% plot(baseline_temp); hold on;
% plot(baseline_temp_actual);

v = 9;

%% Calculate baseline and Left and righ refractory period amps and bases.

%baseline = repelem(mean([left_segment_signal, right_segment_signal]),length(signal));
[baseline, threshold] = calculate_segment_baseline([left_segment_signal, right_segment_signal]);
baseline = repelem(baseline, length(signal));
threshold_baseline = threshold;
left_amp = left_segment_signal(left_refractory_peak_local);
right_amp = right_segment_signal(right_refractory_peak_local);
left_base = baseline(left_refractory_peak_local);
right_base = baseline(right_refractory_peak_local);

%% Compare left and right with the baseline

left_condition_base = left_amp > left_base - threshold_baseline && left_amp < left_base + threshold_baseline;
left_condition_above = left_amp > left_base + threshold_baseline;
left_condition_below = left_amp < left_base - threshold_baseline;

right_condition_base = right_amp > right_base - threshold_baseline && right_amp < right_base + threshold_baseline;
right_condition_above = right_amp > right_base + threshold_baseline;
right_condition_below = right_amp < right_base - threshold_baseline;


signal_to_plot = [left_segment_signal, right_segment_signal];


% figure;
% plot(signal_to_plot); hold on;
% plot(baseline(1:length(signal_to_plot))); hold on;
% plot(left_refractory_peak_local, left_segment_signal(left_refractory_peak_local), 'o'); hold on;
% plot(right_refractory_peak_local + length(left_segment_signal)-1, right_segment_signal(right_refractory_peak_local), 'o'); hold on;

v = 9;



% figure;
% plot(signal); hold on;
% plot(peak_index, signal(peak_index), 'o', 'Color', [1 0 0]); hold on;
% plot(left_refractory_peak_candidate, signal(left_refractory_peak_candidate), 'o', 'Color', [0 0 1]); hold on;
% plot(right_refractory_peak_candidate, signal(right_refractory_peak_candidate), 'o', 'Color', [0 0 1]); hold on;
% 
% v = 9;



if left_condition_base % Left refractory is onset
    left_refractory_peak = left_refractory_peak_candidate;
    onset = left_refractory_peak_candidate;
else % onset is the onset of Left refractory
    [left_refractory_peak, ~, ~, ~, ~, ~] = delineate_wave( signal, left_refractory_peak_candidate, 0, avg_EE_interval, [],1);
    onset = left_refractory_peak;
end

if right_condition_base % Right refractory is offset
    right_refractory_peak = right_refractory_peak_candidate;
    offset = right_refractory_peak_candidate;
else % onset is the onset of Left refractory
    [~, ~, ~, right_refractory_peak, ~, ~] = delineate_wave( signal, right_refractory_peak_candidate, 1, avg_EE_interval, [], 0);
    offset = right_refractory_peak;
end

peak_amp = signal(peak_index);
peak_based = baseline(peak_index);

peak_condition_above = peak_amp > peak_based + threshold_baseline;
peak_condition_below = peak_amp < peak_based - threshold_baseline;


Q = [];
R = [];
S = [];

if peak_condition_above
    R = peak_index;
elseif left_condition_above
    R = left_refractory_peak_candidate;
elseif right_condition_above
    R = right_refractory_peak_candidate; 
end

if left_condition_below
    Q = left_refractory_peak_candidate;
elseif left_condition_base && peak_condition_below
    Q = peak_index;
end

if peak_condition_above && right_condition_below
    S = right_refractory_peak_candidate;
elseif left_condition_above && peak_condition_below
    S = peak_index;
end




% figure;
% plot(signal); hold on;
% plot(peak_index, signal(peak_index), 'o', 'Color', [1 0 0]); hold on;
% plot(left_refractory_peak, signal(left_refractory_peak), 'o', 'Color', [0 0 1]); hold on;
% plot(right_refractory_peak, signal(right_refractory_peak), 'o', 'Color', [0 0 1]); hold on;

v = 9;

% 
% 
% if diff_left < threshold_diff
%     onset = onset_candidate;
%     left_refractory_peak = left_refractory_peak_candidate;
% else
%     [left_refractory_peak, onset, ~, ~, ~, ~] = delineate_wave( signal, left_refractory_peak_candidate, 1, avg_EE_interval, []);
% end
% 
% if diff_right < threshold_diff
%     offset = offset_candidate;
%     right_refractory_peak = right_refractory_peak_candidate;
% else
%     [~, ~, offset, right_refractory_peak, ~, ~] = delineate_wave( signal, right_refractory_peak_candidate, 1, avg_EE_interval, []);
% end

end

