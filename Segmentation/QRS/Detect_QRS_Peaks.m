function [ qrs_peaks, significant_peaks, fragmented_qrs ] = Detect_QRS_Peaks( signal, baseLine, extracted_points, fs )

x = 1 : length(signal);
window_size = calculate_neighbouring_size( signal, extracted_points, fs );

left_neighbouring_points = [];
right_neighbouring_points = [];

max_peak_value = max(abs(signal(extracted_points)));


pivote_index = 14;

% [left_neighbouring_points] = get_left_neighbouring_points(extracted_points, pivote_index, window_size);
% [right_neighbouring_points] = get_right_neighbouring_points(extracted_points, pivote_index, window_size);
% 
% plot(x, signal); hold on; grid on;
% segment_signal = signal(extracted_points(pivote_index) - window_size : extracted_points(pivote_index) + window_size);
% segment_time = x(extracted_points(pivote_index) - window_size : extracted_points(pivote_index) + window_size);
% plot(segment_time, segment_signal, 'r-'); hold on;grid on
% scatter(x(extracted_points(pivote_index)), signal(extracted_points(pivote_index))); hold on;
% scatter(x(left_neighbouring_points), signal(left_neighbouring_points)); hold on;
% scatter(x(right_neighbouring_points), signal(right_neighbouring_points)); hold on;


qrs_peaks = [];
fragmented_qrs = {};
significant_peaks = [];

max_p_signal = max(signal(find(signal > 0)));
max_n_signal = max(abs(signal(find(signal < 0))));
max_abs_signal = max(abs(signal));

for index = 1 : length(extracted_points)
    
    
    left_neighbouring_points = get_left_neighbouring_points(extracted_points, index, window_size);
    right_neighbouring_points = fliplr(get_right_neighbouring_points(extracted_points, index, window_size));
    
    if isempty(left_neighbouring_points) || isempty(right_neighbouring_points) 
        %%  --- BOTH LEFT AND RIGHT NEIGHBORS HAVE TO HAVE AT LEAST 1 POINT. AT LEAST ONE OF THEM DO NOT HAVE POINT.
        if index ~= 1
            
            neighbor_point = 0;
            if ~isempty(left_neighbouring_points)
                neighbor_point = left_neighbouring_points(end);
            elseif ~isempty(right_neighbouring_points)
                neighbor_point = right_neighbouring_points(1);
            end
            
            if neighbor_point > 0 
                
                d_base_candidate = abs(signal(extracted_points(index)) - baseLine(extracted_points(index)));
                d_base_neighbor = abs(signal(neighbor_point) - baseLine(neighbor_point));
                
                if abs(signal(extracted_points(index))) > 0.5 * max_abs_signal && d_base_candidate > d_base_neighbor
                    qrs_peaks = check_and_add_to_results(signal, baseLine, qrs_peaks, extracted_points(index));
                end
            else
                if abs(signal(extracted_points(index))) > 0.5 * max_abs_signal
                    qrs_peaks = check_and_add_to_results(signal, baseLine, qrs_peaks, extracted_points(index));
                end
            end
        else
            if abs(signal(extracted_points(index))) > 0.7 * max_peak_value
                qrs_peaks = check_and_add_to_results(signal, baseLine, qrs_peaks, extracted_points(index));
            end
        end
        
    else
        %% --- BOTH LEFT AND RIGHT HAVE AT LEAST 1 POINT
        
        left_neighbours_contain_peak = check_left_contains_peak(left_neighbouring_points, qrs_peaks);
        
        if ~left_neighbours_contain_peak
            %% --- NO QRS EXTREMUM FOUND IN LEFT NEIGHBOR. GOOD TO GO.
            
            % CALCULATE D_LEFT, D_RIGHT, D_MAX_LEFT, D_MAX_RIGHT
            left_divergence = abs(signal(extracted_points(index)) - signal(left_neighbouring_points(end)));
            left_neighbour_max_divergence = calculate_maximum_divergence(signal, left_neighbouring_points);
            
            right_divergence = abs(signal(right_neighbouring_points(1)) - signal(extracted_points(index)));
            right_neighbour_max_divergence = calculate_maximum_divergence(signal, right_neighbouring_points);
            
            if left_neighbour_max_divergence == 0
                %% LEFT NEIGHBOR HAS ONLY ONE POINT. NO D_MAX_LEFT.
                
                if right_divergence >=  right_neighbour_max_divergence % This is Condition 3.
                    
                    if left_divergence > 0.5 * right_divergence % This is Condition 2.
                        qrs_peaks = check_and_add_to_results(signal, baseLine, qrs_peaks, extracted_points(index));
                    end
                    
                else
                    [max_peak_positive, max_peak_negetive] = calculate_maximum_peak_in_neighbor(signal, right_neighbouring_points);
                    if signal(extracted_points(index)) > 0
                        if signal(extracted_points(index)) > 0.4 * max_p_signal
                            significant_peaks = [significant_peaks, extracted_points(index)];
                        elseif signal(extracted_points(index)) > 0.5 * max_peak_positive
                            significant_peaks = [significant_peaks, extracted_points(index)];
                        end
                    else
                        if abs(signal(extracted_points(index))) > 0.4 * max_n_signal
                            significant_peaks = [significant_peaks, extracted_points(index)];
                        elseif abs(signal(extracted_points(index))) > 0.5 * abs(max_peak_negetive)
                            significant_peaks = [significant_peaks, extracted_points(index)];
                        end
                    end
                end
                
            else
                %% LEFT NEIGHBOR HAS MORE THAN ONE POINT. (ELABORATIVE CASE!)
                
                % CALCULATE LARGENESS FACTOR. IF IN LEFT NEIGHBOR A POINT HAD
                % EVERYTHING EXCEPT R>R_MAX OR R>R_R (IN ITS OWN ITERATION), LARGENESS
                % FACTOR IS 1. OTHERWISE IT IS 1.9.
                largeness_factor = 1.9;
                left_neighbours_contain_significant = check_left_contains_significant(left_neighbouring_points, significant_peaks);
                if left_neighbours_contain_significant
                    largeness_factor = 1;
                end
                
                % BUILDING UP CONDITIONS FOR PRIMERY CANDIDATE DETECTION
                condition1 = left_divergence >= largeness_factor * left_neighbour_max_divergence; % Necassary condition for a QRS extremum.
                condition2 = left_divergence > 0.5 * right_divergence; % Not to select first and small QRS peak.
                
                if condition1 && condition2
                    %% LIKELY TO HAVE QRS EXTREMUM. MUST LOOK AT THE RIGHT NEIGHBOR
                    
                    if length(right_neighbouring_points) == 1 % RIGHT NEIGHBOR HAS ONLY ONE POINT. NO D_MAX_RIGHT. CANDIDATE IS DEFENEITLY EXTREMUM.
                        qrs_peaks = check_and_add_to_results(signal, baseLine, qrs_peaks, extracted_points(index));
                        
                    else % RIGHT NEIGHBOR HAS MORE THAN ONE POINT (ALSO, LEFT NEIGHBOR HAS MORE THAN ONE POINT).
                        condition3 = right_divergence >= right_neighbour_max_divergence;
                        
                        if condition3
                            %% RIGHT POINT AND RIGHT_RIGHT POINT ARE IMPORTANT.
                            signed_right_divergence = signal(right_neighbouring_points(1)) - signal(extracted_points(index));
                            signed_right_right_divergence = signal(right_neighbouring_points(2)) - signal(right_neighbouring_points(1));
                            
                            if signed_right_divergence * signed_right_right_divergence > 0 % No ^ or ^ reverse.
                                qrs_peaks = check_and_add_to_results(signal, baseLine, qrs_peaks, extracted_points(index));
                            else % To the right we have ^ or ^ reveres
                                % We first check for fragmented QRS
                                d_time_right_right = abs(extracted_points(index) - right_neighbouring_points(2));
                                d_time_left_left = abs(extracted_points(index) - left_neighbouring_points(end-1));
                                d_time_left = abs(extracted_points(index) - left_neighbouring_points(end));
                                
                                d_time_left_right = abs(left_neighbouring_points(end) - right_neighbouring_points(1));
                                
                                signed_left_divergence = signal(extracted_points(index)) - signal(left_neighbouring_points(end));
                                signed_left_left_divergence = signal(left_neighbouring_points(end)) - signal(left_neighbouring_points(end-1));
                                
                                %% Fragmented QRS Detection
                                fragmented_qrs_found = 0;
                                
                                if signed_left_divergence > 0
                                    if length(left_neighbouring_points) > 2 % Type 1 Fragmented QRS has to have at least 3 points in its left
                                        signed_left_left_left_divergence = signal(left_neighbouring_points(end-1)) - signal(left_neighbouring_points(end-2));
                                        if signed_left_left_left_divergence > 0 && d_time_left_left < round(0.045 * fs)
                                            fragmented_qrs_found = 1;
                                            qrs_peaks = check_and_add_to_results(signal, baseLine, qrs_peaks, extracted_points(index));
                                            %fragmented_qrs{end+1} = [left_neighbouring_points(end-1), left_neighbouring_points(end), extracted_points(index), right_neighbouring_points(1)];
                                            fragmented_qrs{end+1} = [left_neighbouring_points(end-1), extracted_points(index)];
                                        end
                                    end
                                    
                                    if ~fragmented_qrs_found
                                        if length(right_neighbouring_points) > 2 % Type 2 Fragmented QRS has to have at least 3 points in its right
                                            signed_right_right_right_divergence = signal(right_neighbouring_points(3)) - signal(right_neighbouring_points(2));
                                            
                                            if signed_right_divergence < 0 && signed_right_right_divergence > 0 && signed_right_right_right_divergence < 0 && d_time_right_right < round(0.045 * fs)
                                                fragmented_qrs_found = 1;
                                                qrs_peaks = check_and_add_to_results(signal, baseLine, qrs_peaks, extracted_points(index));
                                                %fragmented_qrs{end+1} = [left_neighbouring_points(end), extracted_points(index), right_neighbouring_points(1), right_neighbouring_points(2)];
                                                fragmented_qrs{end+1} = [extracted_points(index), right_neighbouring_points(2)];
                                            end
                                        end
                                    end
                                    
                                else
                                    if signed_left_left_divergence > 0 && signed_right_right_divergence < 0 && d_time_left_right <= round(0.1 * fs) % Fragmented QRS is detected Type 3
                                        fragmented_qrs_found = 1;
                                        qrs_peaks = check_and_add_to_results(signal, baseLine, qrs_peaks, extracted_points(index));
                                        fragmented_qrs{end+1} = [left_neighbouring_points(end), extracted_points(index), right_neighbouring_points(1)];
                                    end
                                end
                                %% End Fragmented QRS Detection
                                
                                if ~fragmented_qrs_found % It is NOT fragmented QRS
                                    condition4 = abs(signed_right_right_divergence) <= 0.8 * abs(signed_right_divergence);
                                    if condition4
                                        qrs_peaks = check_and_add_to_results(signal, baseLine, qrs_peaks, extracted_points(index));
                                    else % 0.8*R < R_R
                                        significant_peaks = [significant_peaks, extracted_points(index)];
                                    end
                                end
                                
                                
%                                 if  d_time_left_left <= round(0.0625 * fs) && signed_left_divergence > 0 && signed_left_left_divergence < 0
%                                     
%                                     qrs_peaks = check_and_add_to_results(signal, baseLine, qrs_peaks, extracted_points(index));
%                                     fragmented_qrs{end+1} = [left_neighbouring_points(end), extracted_points(index), right_neighbouring_points(1), right_neighbouring_points(2)];
%                                     
%                                 elseif d_time_left <= round(0.03125 * fs) && signed_left_divergence < 0 && signed_left_left_divergence > 0 % Fragmented QRS is detected Type 2
%                                     
%                                     qrs_peaks = check_and_add_to_results(signal, baseLine, qrs_peaks, extracted_points(index));
%                                     fragmented_qrs{end+1} = [left_neighbouring_points(end-1), left_neighbouring_points(end), extracted_points(index), right_neighbouring_points(1)];
%                                     
%                                 elseif d_time_right_right <= round(0.0625 * fs) % Fragmented QRS is detected Type 1
%                                     qrs_peaks = check_and_add_to_results(signal, baseLine, qrs_peaks, extracted_points(index));
%                                     fragmented_qrs{end+1} = [left_neighbouring_points(end-1), left_neighbouring_points(end), extracted_points(index), right_neighbouring_points(1)];
%                                     
%                                 else % It is NOT fragmented QRS
%                                     
%                                     condition4 = abs(signed_right_right_divergence) <= 0.8 * abs(signed_right_divergence);
%                                     if condition4
%                                         qrs_peaks = check_and_add_to_results(signal, baseLine, qrs_peaks, extracted_points(index));
%                                     end
%                                 end
                            end
                            
                        else % R < R_MAX
                            significant_peaks = [significant_peaks, extracted_points(index)];   
                        end         
                    end
                end    
            end  
        else
            %% --- QRS EXTREMUM FOUND IN LEFT NEIGHBOR. CHECK FOR FRAGMENTED QRS.
            
%             if length(left_neighbouring_points) > 1
%                 if qrs_peaks(end) == left_neighbouring_points(end)
%                     d_time = abs(extracted_points(index) - left_neighbouring_points(end-1));
%                     if d_time <= round(0.0625 * fs)
%                         left_divergence = abs(signal(extracted_points(index)) - signal(left_neighbouring_points(end)));
%                         left_left_divergence = abs(signal(left_neighbouring_points(end)) - signal(left_neighbouring_points(end-1)));
%                         if left_divergence >= left_left_divergence
%                             fragmented_qrs{end+1} = [extracted_points(index), left_neighbouring_points(end), left_neighbouring_points(end-1)];
%                         end  
%                     end 
%                 end   
%             end   
        end    
    end
end

end

function qrs_peaks = check_and_add_to_results(signal, baseline, qrs_peaks, candidate_point_index)

if ~isempty(qrs_peaks)
    
    divergence_from_baseline = [];
    for index = 1 : length(qrs_peaks)
        divergence_from_baseline(index) = abs(signal(qrs_peaks(index)) - baseline(qrs_peaks(index)));
    end
    
    last_peak_index = qrs_peaks(end);
    last_peak_divergence_from_baseline = signal(last_peak_index) - baseline(last_peak_index);
    
    avg_peak_divergence_from_baseline = mean(divergence_from_baseline);
    
    candidate_peak_divergence_from_baseline = signal(candidate_point_index) - baseline(candidate_point_index);
    
    %if last_peak_divergence_from_baseline * candidate_peak_divergence_from_baseline > 0
    
    if abs(candidate_peak_divergence_from_baseline) > 0.5 * abs(avg_peak_divergence_from_baseline) && abs(signal(candidate_point_index)) > 0.2 * abs(signal(last_peak_index))
        qrs_peaks = [qrs_peaks, candidate_point_index];
    end
    %end
else
    qrs_peaks = [qrs_peaks, candidate_point_index];
end

end

function [max_divergence] = calculate_maximum_divergence(signal, neighbouring_points)

max_divergence = 0;
divergence_array = [];

if length(neighbouring_points) > 1
    for i = 2 : length(neighbouring_points)
        divergene = abs(signal(neighbouring_points(i)) - signal(neighbouring_points(i-1)));
        divergence_array = [divergence_array, divergene];
    end
end

if ~isempty(divergence_array)
    [value, ~] = max(divergence_array);
    max_divergence = value;
end

end

function [max_peak_positive, max_peak_negetive] = calculate_maximum_peak_in_neighbor(signal, neighbouring_points)

max_peak_positive = 0;
max_peak_negetive = 0;

if ~isempty(neighbouring_points)
    for i = 1 : length(neighbouring_points)
        peak = signal(neighbouring_points(i));
        if peak > 0
            if peak > max_peak_positive
                max_peak_positive = peak;
            end
        else
            if peak < max_peak_negetive
                max_peak_negetive = peak;
            end
        end
        
    end
end

end

function [found] = check_left_contains_peak(left_neighbouring_points, qrs_peaks)

found = 0;
for i = 1 : length(left_neighbouring_points)
    
    point = left_neighbouring_points(i);
    idx = find(qrs_peaks == point);
    if ~isempty(idx)
        found = 1;
        break
    end
end

end

function [found] = check_left_contains_significant(left_neighbouring_points, significant_peaks)

found = 0;
for i = 1 : length(left_neighbouring_points)
    
    point = left_neighbouring_points(i);
    idx = find(significant_peaks == point);
    if ~isempty(idx)
        found = 1;
        break
    end
end

end

function [left_neighbouring_points] = get_left_neighbouring_points(extracted_points, pivot_point_index, window_size)

left_neighbouring_points = [];

for index = 1 : length(extracted_points)
    
    if index < pivot_point_index
        
        if extracted_points(pivot_point_index) - extracted_points(index) < window_size
            
            left_neighbouring_points = [left_neighbouring_points, extracted_points(index)];
        end
        
    end
    
end

end

function [right_neighbouring_points] = get_right_neighbouring_points(extracted_points, pivot_point_index, window_size)

right_neighbouring_points = [];

for index = length(extracted_points) : -1 : 1
    
    if index > pivot_point_index
        
        if extracted_points(index) - extracted_points(pivot_point_index) < window_size
            
            right_neighbouring_points = [right_neighbouring_points, extracted_points(index)];
        end
        
    end
    
end
end

function [ neighbour_size ] = calculate_neighbouring_size( signal, extracted_points, fs )

peaks(:, 2) = signal(extracted_points);
peaks(:, 1) = extracted_points;

[~, idx] = max(abs(peaks(:, 2)));

peak_max(1) = peaks(idx, 1);
peak_max(2) = peaks(idx, 2);

left = get_left_inflection_point(signal, peak_max, fs);
right = get_right_inflection_point(signal, peak_max, fs);

neighbour_size = (right(1) - left(1)) * 4;
neighbour_size = round(0.35 * fs);

end


function left_inflection_point = get_left_inflection_point(signal, mid_point, fs)

blank_period = round(0.200*fs);
gap = 1;

peak_index = mid_point(1);

lower_bound = peak_index - gap - blank_period;
if lower_bound < 1
    lower_bound = 1;
end
upper_bound = peak_index - gap;

left_segment_signal = signal(lower_bound : upper_bound);
left_inflection_point = locate_left_inflection_point(left_segment_signal, fs);
index = left_inflection_point(2);
left_inflection_point(2) = left_inflection_point(1);
left_inflection_point(1) = index + lower_bound - 1;

end

function right_inflection_point = get_right_inflection_point(signal, mid_point, fs)

blank_period = round(0.200*fs);
gap = 1;

peak_index = mid_point(1);
upper_bound = peak_index + gap + blank_period;
if upper_bound > length(signal)
    upper_bound = length(signal);
end
lower_bound = peak_index + gap;

right_segment_signal = signal(lower_bound : upper_bound);
right_inflection_point = locate_right_inflection_point(right_segment_signal, [], []);
index = right_inflection_point(2);
right_inflection_point(2) = right_inflection_point(1);
right_inflection_point(1) = index + lower_bound - 1;

end


function q = locate_left_inflection_point(sig, fs)

q = [];
[~,locs] = findpeaks(sig);
[~,locs1] = findpeaks(-sig);

max_value = max(abs(sig));

if ~isempty(locs)
    if locs(end) == length(sig) - 1
        locs(end) = [];
    end
end

if ~isempty(locs1)
    if locs1(end) == length(sig) - 1
        locs1(end) = [];
    end
end

A = [locs, locs1];
A = sort(A);
B = A;

number_of_removed = 0;

if ~isempty(A)
    if ~isempty(fs)
        period = round(0.04 * fs);
        for i = 1 : length(A) - 1
            if A(i+1) - A(i) <= period && abs(sig(A(i+1)) - sig(A(i))) < 0.1 * max_value
                
                if i-1 > 0
                    
                    if abs(sig(A(i)) - sig(A(i-1))) > 0.1 * max_value
                        jj = i;
                        jj = jj - number_of_removed;
                        B(jj) = [];
                        B(jj) = [];
                        number_of_removed = number_of_removed + 2;
                    end
%                     if ~(A(i) - A(i-1) <= period && abs(sig(A(i)) - sig(A(i-1))) < 0.1 * max_value)
%                         jj = i;
%                         jj = jj - number_of_removed;
%                         B(jj) = [];
%                         B(jj) = [];
%                         number_of_removed = number_of_removed + 2;
%                     end
                end
            end
        end
    end
end

A = B;

if ~isempty(A)
    q(2) = max(A);
else
    
    diff_2_sig = diff(diff(sig));
    [~, diff_2_sig_zeros] = find(abs(diff_2_sig) <= 0.5);
    
    if ~isempty(diff_2_sig_zeros)
        
        diff_1_sig = diff(sig);
        diff_1_sig_candidates = diff_1_sig(diff_2_sig_zeros);
        [~, diff_min_loc] = min(abs(diff_1_sig_candidates));
        q(2) = diff_2_sig_zeros(diff_min_loc(1));
        
    else
        [~,locs]=findpeaks(diff(sig));
        if isempty(locs)
            [~,locs] = max(diff(sig,2));
        else
            locs = locs(1);
        end
        q(2) = locs;
    end
    
end
q(1) = sig(q(2));
end

function s = locate_right_inflection_point(sig, R_Notched, offset)

s = [];
[~,locs] = findpeaks(sig);
[~,locs1] = findpeaks(-sig);

if ~isempty(locs)
    if locs(1) == length(sig) - 1
        locs(1) = [];
    end
end

if ~isempty(locs1)
    if locs1(1) == length(sig) - 1
        locs1(1) = [];
    end
end


notched_found_at_local = [];
% if ~isempty(locs)
%     if ~isempty(find(R_Notched(:,1) == locs(1) + offset - 1))
%         notched_found_at_local = locs(1);
%         locs(1) = []; 
%     end
% end
% 
% if ~isempty(locs1)
%     if ~isempty(find(R_Notched(:,1) == locs1(1) + offset - 1))
%         notched_found_at_local = locs1(1);
%         locs1(1) = [];
%     end
% end

if ~isempty(notched_found_at_local)
    if ~isempty(locs)
        elements_to_be_deleted = [];
        for i = 1 : length(locs)
            if locs(i) <= notched_found_at_local
                elements_to_be_deleted = [elements_to_be_deleted, i];
            end
        end
        if ~isempty(elements_to_be_deleted)
            locs(elements_to_be_deleted) = [];
        end
    end
    
    if ~isempty(locs1)
        elements_to_be_deleted = [];
        for i = 1 : length(locs1)
            if locs1(i) <= notched_found_at_local
                elements_to_be_deleted = [elements_to_be_deleted, i];
            end
        end
        if ~isempty(elements_to_be_deleted)
            locs1(elements_to_be_deleted) = [];
        end
    end
end



if ~isempty(locs)
    locs = locs(1);
end
if ~isempty(locs1)
    locs1 = locs1(1);
end
A = [locs;locs1];

if ~isempty(A)
    s(2) = min(A);
     s(1) = sig(s(2));
    
else
    %     sig = -1 * sig;
    %     subplot(4,1,1)
    %     plot(sig)
    %     subplot(4,1,2)
    %     plot(diff(sig))
    %     subplot(4,1,3)
    %     plot(diff(diff(sig)))
    %     subplot(4,1,4)
    %     plot(diff(diff(diff(sig))))
    
    diff_2_sig = diff(diff(sig));
    [~, diff_2_sig_zeros] = find(abs(diff_2_sig) <= 0.5);
    
    if ~isempty(diff_2_sig_zeros)
        
        diff_1_sig = diff(sig);
        diff_1_sig_candidates = diff_1_sig(diff_2_sig_zeros);
        [~, diff_min_loc] = min(abs(diff_1_sig_candidates));
        s(2) = diff_2_sig_zeros(diff_min_loc(1));
        
    else
        [~,locs]=findpeaks(diff(sig));
        if isempty(locs)
            [~,locs] = max(diff(sig,2));
        else
            locs = locs(1);
        end
        s(2) = locs;
    end
    
    s(1) = sig(s(2));
end

end



