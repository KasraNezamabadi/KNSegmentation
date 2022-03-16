function [ P_peaks, P_left_sets, P_right_sets] = delineate_P_waves( ecg_filtered, T_right_sets, qrs_left_sets, avg_EE_interval, sample_rate)


P_peaks = [];
P_left_sets = {}; % [onset, refractory_peak]
P_right_sets = {}; % [offset, refractory_peak]
LineFittingThresholdFactor = 8;

for i = 2 : length(qrs_left_sets)
    
    end_index = cell2mat(qrs_left_sets(i));
    if i - 1 > length(T_right_sets)
        break;
    end
    start_index = cell2mat(T_right_sets(i-1));
    end_index = end_index(2);
    start_index = start_index(1);
    
    gap = 4;
    
    if end_index > start_index + gap
         
        segment = ecg_filtered(start_index + gap : end_index);
        segment_baseline = repelem(mean(segment),length(segment));
        
        segment_smoothed = smooth(1:length(segment), segment, round(0.3*length(segment)), 'loess');
        segment_smoothed = segment_smoothed';
        segment_baseline_smoothed = repelem(mean(segment_smoothed), length(segment_smoothed));
        
        threshold = max(abs(segment_smoothed - segment_baseline_smoothed)) / LineFittingThresholdFactor;
        extracted_points = lineFittingAlgorithmYDistance(segment_smoothed, threshold);
        
        corrected_extracted_points = [];
        [ A ] = get_all_peaks_valleys(segment_smoothed);
        if ~isempty(A)
            window_size = round(0.02 * sample_rate);
            for j = 1 : length(extracted_points)
                extracted_point = extracted_points(j);
                distance_to_peaks = abs(A - extracted_point);
                [min_distance, min_index] = min(distance_to_peaks);
                if min_distance <= window_size
                    corrected_extracted_points = [corrected_extracted_points, A(min_index)];
                end
                
            end
        end
        extracted_points = corrected_extracted_points;
        extracted_points = unique(extracted_points,'first');
        
        
        close all;
        plot(segment_smoothed); hold on;
        plot(segment_baseline_smoothed); hold on;
        plot(extracted_points, segment_smoothed(extracted_points), 'o'); hold on;
        
        if i == 12
            v = 9;
        end
        
        [p_peaks_result] = p_peak_detection_recursively(segment_smoothed, extracted_points, 10, []);
        p_peaks_result = p_peaks_result + start_index + gap - 1;
        
        for jj = 1 : length(p_peaks_result)
            [ left_refractory_peak, onset, offset, right_refractory_peak, ~, ~] = delineate_wave( ecg_filtered, p_peaks_result(jj), 0, avg_EE_interval, end_index, 0);
            P_left_sets{end+1} = [onset, left_refractory_peak];
            P_right_sets{end+1} = [offset, right_refractory_peak];
            P_peaks = [P_peaks, p_peaks_result(jj)];
        end
        
%         close all;
%         plot(segment_smoothed); hold on;
%         plot(segment_baseline_smoothed); hold on;
%         plot(p_peaks, segment_smoothed(p_peaks), 'o'); hold on;
        
    end
  
end

last_right = cell2mat(T_right_sets(end));
last_right = last_right(1);

end_index = length(ecg_filtered);
start_index = last_right;



gap = 4;
segment = ecg_filtered(start_index + gap : end_index - gap);
segment_baseline = repelem(mean(segment),length(segment));

threshold = max(abs(segment - segment_baseline)) / 4;
if threshold < 10
    threshold = 10;
end

extracted_points = Extract_Significant_Points(segment, threshold);

if ~isempty(extracted_points)
    [~, max_index] = max(abs(segment(extracted_points) - segment_baseline(extracted_points)));
    P_peak_index_local = extracted_points(max_index);
    P_peak_index = P_peak_index_local + start_index + gap - 1;
    
    [ left_refractory_peak, onset, offset, right_refractory_peak, ~, ~] = delineate_wave( ecg_filtered, P_peak_index, 0, avg_EE_interval, end_index, 0);
    
    P_left_sets{end+1} = [onset, left_refractory_peak];
    P_right_sets{end+1} = [offset, right_refractory_peak];
    P_peaks = [P_peaks, P_peak_index];
end

end