function [baseline, threshold] = calculate_segment_baseline(signal)

signal_avg = mean(signal);
range_threshold = 5;
points_above_mean = signal(find(signal > signal_avg));
points_below_mean = signal(find(signal < signal_avg));

baseline = repelem(signal_avg, length(signal));
dist_from_baseline = abs(signal - baseline);
threshold = round(1.5 * mean(dist_from_baseline));


if (length(points_above_mean) > length(points_below_mean) - range_threshold) && (length(points_above_mean) < length(points_below_mean) + range_threshold) % same range
    
    baseline = signal_avg;
    dist_from_baseline = abs(signal - baseline);
    threshold = round(1.5 * mean(dist_from_baseline));
    
elseif length(points_above_mean) > length(points_below_mean) + range_threshold % above are more points
    
    baseline = mean(points_above_mean);
    dist_from_baseline = abs(points_above_mean - baseline);
    threshold = round(1.5 * mean(dist_from_baseline));
    
elseif length(points_below_mean) > length(points_above_mean) + range_threshold % below are more people
    baseline = mean(points_below_mean);
    dist_from_baseline = abs(points_below_mean - baseline);
    threshold = round(1.5 * mean(dist_from_baseline));
end

end

