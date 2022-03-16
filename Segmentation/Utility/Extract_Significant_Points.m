function [ inflection_points] = Extract_Significant_Points( signal, threshold )

x = 1 : length(signal);
inflection_points = [];

if length(signal) > 1
    
    xFit = [x(1), x(end)];
    yFit = [signal(1), signal(end)];
    
    coefficients = polyfit([xFit(1), xFit(2)], [yFit(1), yFit(2)], 1);    
    line = polyval(coefficients, xFit(1) : xFit(2));
    v1 = [xFit(1), yFit(1)];
    v2 = [xFit(2), yFit(2)];
    
    distance_array = zeros(1, length(signal));
    
    
    for index = 1 : length(signal)
        
        point = [x(index), signal(index)];
        distance = distance_point_to_line(point, v1, v2);
%         distance = abs(signal(index) - line(index));
        distance_array(index) = distance;
    end
    
  
    [max_values, max_indecies] = max(distance_array);
    max_value = max_values(1);
    max_index = max_indecies(1);
    
%     plot(x, signal); hold on;
%     plot(x, line); hold on;
%     scatter(max_index, signal(max_index)); hold off;

    diff_signal = diff(signal);
    
    diff_prev_prev = 100;
    diff_prev = 100;
    diff_next = 100;
    diff_next_next = 100;

    
    if max_index - 2 >= 1
        diff_prev_prev = diff_signal(max_index - 2);
    end
    
    if max_index - 1 >= 1
        diff_prev = diff_signal(max_index - 1);
    end
    
    if max_index + 1 <= length(diff_signal)
        diff_next = diff_signal(max_index + 1);
    end
    
    if max_index + 2 <= length(diff_signal)
        diff_next = diff_signal(max_index + 2);
    end
    
%     diff_current = diff_signal(max_index);
    
    
    eligible = 1;
    
    
%     eligible = 0;
%     
%     
%     if (diff_prev < 0 && diff_next > 0) || (diff_prev > 0 && diff_next < 0)
%         eligible = 1;
%     elseif abs(diff_current) < 2 || abs(diff_prev) < 2 || abs(diff_next) < 2
%         eligible = 1;
%     end

    
    if max_value > threshold && eligible == 1
        
        left_signal = signal(1 : max_index - 1);
        right_signal = signal(max_index + 1 : end);
        
        
        threshold_left = max(abs(left_signal)) / 10;
        if threshold_left < 25
            threshold_left = 25;
        end
        threshold_right = max(abs(right_signal)) / 10;
        if threshold_right < 25
            threshold_right = 25;
        end
        
        result_left = Extract_Significant_Points(left_signal, threshold_left);
        result_right = Extract_Significant_Points(right_signal, threshold_right);
        result_right = result_right + max_index;
        

        inflection_points = [inflection_points, result_left];
        inflection_points = [inflection_points, max_index];
        inflection_points = [inflection_points, result_right];
        
    end
    
end


end

