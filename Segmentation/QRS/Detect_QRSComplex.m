function [ left_sets, right_sets, Q_list, R_list, R2_list, S_list] = Detect_QRSComplex( signal, qrs_extremums, fragmented_qrs, EE_interval )


left_sets = {}; % [onset, refractory_peak]
right_sets = {}; % [offset, refractory_peak]

Q_list = [];
R_list = [];
R2_list = [];
S_list = [];



for i = 1 : length(qrs_extremums)
    extremum_index = qrs_extremums(i);
    
    [is_fragmented, qrs_fragments] = is_qrs_fragmented(extremum_index, fragmented_qrs);
    
    if is_fragmented
        first_peak = qrs_fragments(1);
        last_peak = qrs_fragments(end);
        
        [~, ~, offset, right_refractory_peak,~,~,~] = Detect_QRS_Bounds( signal, last_peak, EE_interval);
        right_sets{end+1} = [offset, right_refractory_peak];
        
        [left_refractory_peak, onset, ~, ~, ~, ~, ~] = Detect_QRS_Bounds( signal, first_peak, EE_interval);
        left_sets{end+1} = [onset, left_refractory_peak];
        
        R_list = [R_list, first_peak];
        R2_list = [R2_list, last_peak];
        
    else
        
        [left_refractory_peak, onset, offset, right_refractory_peak, Q, R, S] = Detect_QRS_Bounds( signal, extremum_index, EE_interval);
        
        left_sets{end+1} = [onset, left_refractory_peak];
        right_sets{end+1} = [offset, right_refractory_peak];
        
        Q_list = [Q_list, Q];
        R_list = [R_list, R];
        S_list = [S_list, S];
        
    end
    
end

Q_list = correct_QRS(signal, Q_list);
R_list = correct_QRS(signal, R_list);
S_list = correct_QRS(signal, S_list);


end

function [Q_list] = correct_QRS(signal, Q_list)

for i = 1 : length(Q_list)
    
    if i == 6
        v = 9;
    end
    
    candidate_peak = Q_list(i);
    period = 6;
    lower_bound = candidate_peak - period;
    upper_bound = candidate_peak + period;
    if lower_bound < 1
        lower_bound = 1;
    end
    if upper_bound > length(signal)
        upper_bound = length(signal);
    end
    
    if upper_bound > lower_bound + 2
        
        segment = signal(lower_bound : upper_bound);
        diff_signal = diff(signal);
        left_slop = diff_signal(candidate_peak - 1);
        right_slop = diff_signal(candidate_peak + 1);
        
%         figure
%         plot(segment);
%         
%         figure;
%         plot(signal); hold on;
%         plot(candidate_peak, signal(candidate_peak), 'o');
        
        
        candidate_peak = candidate_peak - (lower_bound-1);
        
        % You may want to delete this condition!! Some cases only 1 point
        % is to the top
        %if left_slop * right_slop > 0 % candidate_peak is not peak
        
        [~,locsP] = findpeaks(segment);
        [~,locsN] = findpeaks(-segment);
        if ~isempty(locsP)
            if locsP(1) == length(segment) - 1
                locsP(1) = [];
            end
        end
        if ~isempty(locsN)
            if locsN(1) == length(segment) - 1
                locsN(1) = [];
            end
        end
        A = [locsP, locsN];
        A = sort(A);
        
        if ~isempty(A)
            
            dists = abs(A - candidate_peak);
            [~, index] = min(dists);
            peak_corrected_local = A(index);
            peak_corrected = peak_corrected_local + lower_bound - 1;
            Q_list(i) = peak_corrected;
            
%             figure
%             plot(segment);
%             
%             figure;
%             plot(signal); hold on;
%             plot(peak_corrected, signal(peak_corrected), 'o');
            
        end
        
        % end
        
        
        
    end
    v = 9;
    %close all;
    
end
end

function [found, qrs_peaks] = is_qrs_fragmented(extremum_index, fragmented_qrs)

found = 0;
qrs_peaks = [];

for i = 1 : length(fragmented_qrs)
    qrs_fragments = cell2mat(fragmented_qrs(i));
    if sum(ismember(qrs_fragments, extremum_index)) > 0
        found = 1;
        qrs_peaks = qrs_fragments;
        break;
    end
end

end