function [ T_peaks, T_left_sets, T_right_sets, T_biphasic] = delineate_T_waves( ecg_filtered, left_sets, right_sets, avg_EE_interval, sample_rate)


T_peaks = [];
T_biphasic = {};
T_left_sets = {}; % [onset, refractory_peak]
T_right_sets = {}; % [offset, refractory_peak]


LineFittingThresholdFactor = 8;
TSegmentTrimFactor = 8;

for i = 2 : length(left_sets)

    end_index = cell2mat(left_sets(i));
    start_index = cell2mat(right_sets(i-1));
    end_index = end_index(1);
    start_index = start_index(1);
    
    end_index = end_index - round((end_index - start_index)/TSegmentTrimFactor);
    
    if end_index - start_index > avg_EE_interval
        end_index = start_index + round(0.8*avg_EE_interval);
    end
    
    
    if end_index > start_index
        
        gap = 4;
        segment = ecg_filtered(start_index + gap : end_index - gap);
        segment_baseline = repelem(mean(segment),length(segment));
        
%         close all;
%         figure;
%         plot(segment); hold on;
%         plot(segment_baseline);
        
        if i == 13
            v = 9;
        end
        
        x = [1:length(segment)];
        segment_smoothed = smooth(x, segment, round(0.2*length(segment)), 'loess');
        segment_smoothed = segment_smoothed';
        segment_baseline_smoothed = repelem(mean(segment_smoothed), length(segment_smoothed));
                
        %threshold = max(abs(segment_smoothed - segment_baseline_smoothed)) / 4;
        
        [~,locsP] = findpeaks(segment_smoothed);
        [~,locsN] = findpeaks(-segment_smoothed);
        if ~isempty(locsP)
            if locsP(1) == length(segment_smoothed) - 1
                locsP(1) = [];
            end
        end
        if ~isempty(locsN)
            if locsN(1) == length(segment_smoothed) - 1
                locsN(1) = [];
            end
        end
        A = [locsP, locsN];
        A = sort(A);
        
%         if ~isempty(A)
%             
%             threshold = max(abs(segment_smoothed(A) - segment_baseline_smoothed(A))) / 2.5;
%             
%         else
%             threshold = max(abs(segment_smoothed - segment_baseline_smoothed)) / 4;
%         end
        
        threshold = max(abs(segment_smoothed - segment_baseline_smoothed)) / LineFittingThresholdFactor;
        
        %extracted_points = lineFittingAlgorithm(segment, threshold);
        extracted_points = lineFittingAlgorithmYDistance(segment_smoothed, threshold);
        
%         close all;
%         figure;
%         plot(segment_smoothed); hold on;
%         plot(segment_baseline); hold on;
%         plot(extracted_points, segment_smoothed(extracted_points), 'o'); hold on;
%         plot(A, segment_smoothed(A), 'o'); hold on;
        
        corrected_extracted_points = [];
        
        if ~isempty(A)
            % Correct extracted points
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
        
%         close all;
%         figure;
%         plot(segment_smoothed); hold on;
%         plot(segment_baseline); hold on;
%         plot(extracted_points, segment_smoothed(extracted_points), 'o');
        
        no_t_exists = 0;
        
        if isempty(extracted_points)
            % No T.
            no_t_exists = 1;
        else
            if length(extracted_points) == 1
                T_peak_index_local = extracted_points(1);
            else
                divergence_array = [];
                for j = 2 : length(extracted_points)
                    divergence_array(j-1) = segment_smoothed(extracted_points(j)) - segment_smoothed(extracted_points(j-1));
                end
                
                [~, max_ind] = max(abs(divergence_array));
                
                first_point = extracted_points(max_ind);
                second_point = extracted_points(max_ind+1);
                
                if max_ind+2 <= length(extracted_points) % Right exists
                    
                    third_point = extracted_points(max_ind+2);
                    max_divergence = divergence_array(max_ind);
                    right_divergence = segment_smoothed(third_point) - segment_smoothed(second_point);
                    
                    
                    if max_ind - 1 > 0 % Left exists too. Check for the biphasic T wave.
                        
                        zero_point = extracted_points(max_ind - 1);
                        left_divergence = segment_smoothed(first_point) - segment_smoothed(zero_point);
                        
                        d_base_first = segment_smoothed(first_point) - segment_baseline_smoothed(first_point);
                        d_base_second = segment_smoothed(second_point) - segment_baseline_smoothed(second_point);
                        
                        if d_base_first * d_base_second < 0 % Condition 1 for biphasic
                            
                            condition2 = (abs(left_divergence) < 0.8 * abs(max_divergence)) && (abs(right_divergence) < 0.8 * abs(max_divergence));
                            
                            if condition2 % Condition 2 for biphasic
                                
                                T_biphasic{end+1} = [first_point + start_index + gap - 1, second_point + start_index + gap - 1];
                                
                            end
                            
                        end
                        
                    end
                    
                    if max_divergence * right_divergence < 0
                        if (abs(right_divergence) > 0.7 * abs(max_divergence))
                            T_peak_index_local = second_point;
                        else
                            T_peak_index_local = first_point;
                        end
                        
                    else % right is T
                        
                        d_base_first = segment_smoothed(first_point) - segment_baseline_smoothed(first_point);
                        d_base_thrid = segment_smoothed(third_point) - segment_baseline_smoothed(third_point);
                        
                        if abs(d_base_first) > abs(d_base_thrid)
                            T_peak_index_local = first_point;
                        else
                            T_peak_index_local = third_point;
                        end
                        
                    end
                    
                else % No right. T is the one with max divergence from the baseline.
                    
                    d_base_first = segment_smoothed(first_point) - segment_baseline_smoothed(first_point);
                    d_base_second = segment_smoothed(second_point) - segment_baseline_smoothed(second_point);
                    
                    if abs(d_base_first) > abs(d_base_second)
                        T_peak_index_local = first_point;
                    else
                        T_peak_index_local = second_point;
                    end

                end
                
            end
        end
        
        if ~no_t_exists
            
            T_peak_index = T_peak_index_local + start_index + gap - 1;
            [ left_refractory_peak, onset, offset, right_refractory_peak, ~, ~] = delineate_wave( ecg_filtered, T_peak_index, 0, avg_EE_interval, [], 0);
            
            T_left_sets{end+1} = [onset, left_refractory_peak];
            T_right_sets{end+1} = [offset, right_refractory_peak];
            T_peaks = [T_peaks, T_peak_index];
            
            
            %         close all;
            %         figure;
            %         plot(segment_smoothed); hold on;
            %         plot(segment_baseline); hold on;
            %         plot(T_peak_index_local, segment_smoothed(T_peak_index_local), 'o');
        end
        

    end
    
end

last_right = cell2mat(right_sets(end));

if ~isempty(last_right)
    
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
        T_peak_index_local = extracted_points(max_index);
        T_peak_index = T_peak_index_local + start_index + gap - 1;
        
        if abs(length(ecg_filtered) - T_peak_index) < round(0.2 * avg_EE_interval)
            % Do nothing
        else
            
            [ left_refractory_peak, onset, offset, right_refractory_peak, ~, ~] = delineate_wave( ecg_filtered, T_peak_index, 0, avg_EE_interval, [],0 );
            
            T_left_sets{end+1} = [onset, left_refractory_peak];
            T_right_sets{end+1} = [offset, right_refractory_peak];
            T_peaks = [T_peaks, T_peak_index];
        end
    end
end

end