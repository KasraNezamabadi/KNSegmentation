function [p_peaks] = p_peak_detection_recursively(segment, extracted_points, previous_divergence, previous_sign)

p_peaks = [];
if length(extracted_points) < 2
    % No P peak with less than 2 extracted points. Stop Condition!
else
    segment_baseline = repelem(mean(segment),length(segment));
    if length(extracted_points) == 2
        
        first_pont = extracted_points(1);
        second_point = extracted_points(2);
        divergence = segment(second_point) - segment(first_pont);
        
        if abs(divergence) > 0.9 * abs(previous_divergence) % Must be in the range of previous peak to be a P peak. 
            
            first_divergence_from_baseline = segment(first_pont) - segment_baseline(first_pont);
            second_divergence_from_baseline = segment(second_point) - segment_baseline(second_point);
            
            if ~isempty(previous_sign)
                
                if abs(first_divergence_from_baseline) > abs(second_divergence_from_baseline) % First is the peak IF has the same sign as the previous peak
                    
                    if first_divergence_from_baseline * previous_sign > 0
                        p_peaks = first_pont;
                    end
                    
                else % Second is the peak IF has the same sign as the previous peak
                    if second_divergence_from_baseline * previous_sign > 0
                        p_peaks = second_point;
                    end
                end
            else
                if abs(first_divergence_from_baseline) > abs(second_divergence_from_baseline) % First is the peak
                    
                    p_peaks = first_pont;
                    
                else % Second is the peak
                    p_peaks = second_point;
                end
            end
            
        end
        
    else % Number of extracted points are at least 3
        
        divergence_array = [];
        for i = 2 : length(extracted_points)
            divergence = segment(extracted_points(i)) - segment(extracted_points(i-1));
            divergence_array(i-1) = divergence;
        end
        
        [~, max_index] = max(abs(divergence_array));
        
        first_point = extracted_points(max_index);
        second_point = extracted_points(max_index+1);
        max_divergence = segment(second_point) - segment(first_point);
        
        if abs(max_divergence) > 0.9 * abs(previous_divergence) % Must be in the range of previous peak to be a P peak.
            
            left_of_first_divergence = 0;
            right_of_second_divergence = 0;
            
            if max_index - 1 > 0
                left_of_first_divergence = segment(first_point) - segment(extracted_points(max_index - 1));
            end
            
            if max_index+2 <= length(extracted_points)
                right_of_second_divergence = segment(extracted_points(max_index+2)) - segment(second_point);
            end
            
            temp_result = [];
            temp_sign = [];
            
            
            if abs(left_of_first_divergence) > abs(right_of_second_divergence) % first_point is P peak IF sign matches.
                
                if abs(left_of_first_divergence) > 0.8 * abs(max_divergence)
                    temp_result = first_point;
                    temp_sign = left_of_first_divergence > 0;
                else
                    temp_result = second_point;
                    temp_sign = max_divergence > 0;
                end
                
                
                
            elseif abs(left_of_first_divergence) < abs(right_of_second_divergence) % second_point is P peak IF sign matches.
                
                if abs(right_of_second_divergence) > 0.8 * abs(max_divergence)
                    temp_result = second_point;
                    temp_sign = right_of_second_divergence < 0;
                else
                    temp_result = first_point;
                    temp_sign = max_divergence > 0;
                end
            end
            
            if temp_sign == 0
                temp_sign = -1;
            end
            
            if ~isempty(previous_sign)
                
                if temp_sign * previous_sign > 0 % signs matche. temp_result is p peak
                    
                    if abs(max_divergence) > abs(previous_divergence)
                        call_back_max_divergence = abs(max_divergence);
                    else
                        call_back_max_divergence = abs(previous_divergence);
                    end
                    
                    left_extracted_points = extracted_points(find(extracted_points<temp_result));
                    right_extracted_points = extracted_points(find(extracted_points>temp_result));
                    
                    [p_peaks_left] = p_peak_detection_recursively(segment, left_extracted_points, call_back_max_divergence, previous_sign);
                    [p_peaks_right] = p_peak_detection_recursively(segment, right_extracted_points, call_back_max_divergence, previous_sign);
                    
                    p_peaks = [p_peaks, p_peaks_left];
                    p_peaks = [p_peaks, temp_result];
                    p_peaks = [p_peaks, p_peaks_right];
                    
                end
                
            else % No previous sign exists. Current is P peak.
                
                if abs(max_divergence) > abs(previous_divergence)
                    call_back_max_divergence = abs(max_divergence);
                else
                    call_back_max_divergence = abs(previous_divergence);
                end
                
                left_extracted_points = extracted_points(find(extracted_points<temp_result));
                right_extracted_points = extracted_points(find(extracted_points>temp_result));
                
                [p_peaks_left] = p_peak_detection_recursively(segment, left_extracted_points, call_back_max_divergence, temp_sign);
                [p_peaks_right] = p_peak_detection_recursively(segment, right_extracted_points, call_back_max_divergence, temp_sign);
                
                p_peaks = [p_peaks, p_peaks_left];
                p_peaks = [p_peaks, temp_result];
                p_peaks = [p_peaks, p_peaks_right];
                
            end
            
        end
        
        
        
        
    end
end




end