function [ left_refractory_peak, onset, offset, right_refractory_peak, diff_left, diff_right] = delineate_wave( signal, peak_index, is_complex, avg_EE_interval, right_bound, is_QorS)

close all;

left_refractory_peak = [];
onset = [];
offset = [];
right_refractory_peak = [];
diff_left = [];
diff_right = [];

%% Thresholds and boundaries

if is_complex
    blank_period_refractory = round(0.12 * avg_EE_interval);
    blank_period_delineation = round(0.12 * avg_EE_interval);
    if blank_period_refractory > 50
        blank_period_refractory = 50;
    end
    if blank_period_delineation > 50
        blank_period_delineation = 50;
    end
%     blank_period_refractory = round(0.085 * avg_EE_interval);
%     blank_period_delineation = round(0.085 * avg_EE_interval);
    gap = 1;
    smooth_window = 0.6;
elseif is_QorS
    blank_period_refractory = round(0.08 * avg_EE_interval);
    blank_period_delineation = round(0.08 * avg_EE_interval);
    if blank_period_refractory > 30
        blank_period_refractory = 30;
    end
    if blank_period_delineation > 30
        blank_period_delineation = 30;
    end
    gap = 1;
    smooth_window = 0.6;
else
    blank_period_refractory = round(0.2 * avg_EE_interval);
    blank_period_delineation = round(0.2 * avg_EE_interval);
    gap = 5;
    smooth_window = 0.4;
end

%blank_period = round(0.200*fs);


%% Smoothing out the raw signal
% x = [1 : length(signal)];
% signalRaw = signal;
% signal = smooth(x, signal, 0.004, 'loess');

%% Onset. Left segment


lower_bound = peak_index - gap - blank_period_refractory;
if lower_bound < 1
    lower_bound = 1;
end
upper_bound = peak_index - gap;
left_segment_signal = signal(lower_bound : upper_bound);

% plot(left_segment_signal);
% v = 9;
% close all;

if length(left_segment_signal) > 2
    
    x = [1 : length(left_segment_signal)];
    signalRaw = left_segment_signal;
    left_segment_signal = smooth(x, left_segment_signal, smooth_window, 'loess');
    left_segment_signal = left_segment_signal';
    
    if is_complex % no smoothing
        %left_segment_signal = signalRaw;
    end
    
    %refractory_period_peak_local = lineFittingPrependecularAlgorithm(left_segment_signal);
    refractory_period_peak_local = lineFittingYDistAlgorithm(left_segment_signal);
    
    [~,locsP] = findpeaks(left_segment_signal);
    [~,locsN] = findpeaks(-left_segment_signal);
    if ~isempty(locsP)
        if locsP(1) == length(left_segment_signal) - 1
            locsP(1) = [];
        end
    end
    if ~isempty(locsN)
        if locsN(1) == length(left_segment_signal) - 1
            locsN(1) = [];
        end
    end
    A = [locsP, locsN];
    A = sort(A);
    
    min_dist_from_min = round(0.3 * length(left_segment_signal));
    index_of_closest_peak = 0;
    for i = 1 : length(A)
        dist = abs(A(i) - refractory_period_peak_local) ;
        if dist < min_dist_from_min && A(i) <= refractory_period_peak_local
            %refractory_period_peak_local = A(i);
            index_of_closest_peak = A(i);
            min_dist_from_min = dist;
        end
    end
    
    if index_of_closest_peak > 0
        refractory_period_peak_local = index_of_closest_peak;
    end
    
    left_refractory_peak_amp = left_segment_signal(refractory_period_peak_local);
    
    refractory_period_peak = refractory_period_peak_local + lower_bound - 1;
    left_refractory_peak = refractory_period_peak;
    
%     close all;
%     figure;
%     plot(signalRaw)
%     
%     figure;
%     x = [1 : length(left_segment_signal)];
%     plot(left_segment_signal); hold on;
%     plot(x(refractory_period_peak_local), left_segment_signal(refractory_period_peak_local), 'o');
    
    lower_bound = refractory_period_peak - gap - blank_period_delineation;
    if lower_bound < 1
        lower_bound = 1;
    end
    upper_bound = refractory_period_peak - gap;
    left_segment_signal = signal(lower_bound : upper_bound);
    
    if length(left_segment_signal) > 2
        
        x = [1 : length(left_segment_signal)];
        signalRaw = left_segment_signal;
        left_segment_signal = smooth(x, left_segment_signal, smooth_window, 'loess');
        left_segment_signal = left_segment_signal';
        if is_complex % no smoothing
            %left_segment_signal = signalRaw;
        end
        
        % figure;
        % plot(left_segment_signal)
        %
        % figure;
        % plot(signalRaw)
        
        if length(left_segment_signal) > 3
            
            [~,locsP] = findpeaks(left_segment_signal);
            [~,locsN] = findpeaks(-left_segment_signal);
            
            if ~isempty(locsP)
                if locsP(end) == length(left_segment_signal) - 1
                    locsP(end) = [];
                end
            end
            
            if ~isempty(locsN)
                if locsN(end) == length(left_segment_signal) - 1
                    locsN(end) = [];
                end
            end
            
            
            A = [locsP, locsN];
            A = sort(A);
            
            if ~isempty(A)
                onset_amp = left_segment_signal(A(end));
                onset = A(end) + lower_bound - 1;
            else
                diff_2_sig = diff(diff(left_segment_signal));
                [~, diff_2_sig_zeros] = find(abs(diff_2_sig) <= 0.5);
                
                if ~isempty(diff_2_sig_zeros)
                    
                    diff_1_sig = diff(left_segment_signal);
                    diff_1_sig_candidates = diff_1_sig(diff_2_sig_zeros);
                    [~, diff_min_loc] = min(abs(diff_1_sig_candidates));
                    onset_amp = left_segment_signal(diff_2_sig_zeros(diff_min_loc(1)));
                    onset = diff_2_sig_zeros(diff_min_loc(1)) + lower_bound - 1;
                else
                    [~,locs]=findpeaks(diff(left_segment_signal));
                    if isempty(locs)
                        [~,locs] = max(diff(left_segment_signal,2));
                    else
                        locs = locs(1);
                    end
                    onset_amp = left_segment_signal(locs);
                    onset = locs + lower_bound - 1;
                end
            end
            diff_left = abs(left_refractory_peak_amp - onset_amp);
        end
    end
    
end


%% Offset. Right segment

lower_bound = peak_index + gap;
if ~isempty(right_bound)
    upper_bound = right_bound - gap;
    if upper_bound - lower_bound < 3
        upper_bound = peak_index + gap + blank_period_refractory;
    end
else
    upper_bound = peak_index + gap + blank_period_refractory;
end
if upper_bound > length(signal)
    upper_bound = length(signal);
end

right_segment_signal = signal(lower_bound : upper_bound);

if length(right_segment_signal) > 2
    
    x = [1 : length(right_segment_signal)];
    signalRaw = right_segment_signal;
    %length(right_segment_signal)
    right_segment_signal = smooth(x, right_segment_signal, smooth_window, 'loess');
    right_segment_signal = right_segment_signal';
    if is_complex % no smoothing
        %right_segment_signal = signalRaw;
    end
    
    %refractory_period_peak_local = lineFittingPrependecularAlgorithm(right_segment_signal);
    refractory_period_peak_local = lineFittingYDistAlgorithm(right_segment_signal);
    
    [~,locsP] = findpeaks(right_segment_signal);
    [~,locsN] = findpeaks(-right_segment_signal);
    if ~isempty(locsP)
        if locsP(1) == length(right_segment_signal) - 1
            locsP(1) = [];
        end
    end
    if ~isempty(locsN)
        if locsN(1) == length(right_segment_signal) - 1
            locsN(1) = [];
        end
    end
    A = [locsP, locsN];
    A = sort(A);
    
    min_dist_from_min = round(0.3 * length(right_segment_signal));
    index_of_closest_peak = 0;
    for i = 1 : length(A)
        dist = abs(A(i) - refractory_period_peak_local) ;
        if dist < min_dist_from_min && A(i) >= refractory_period_peak_local
            %refractory_period_peak_local = A(i);
            index_of_closest_peak = A(i);
            min_dist_from_min = dist;
        end
    end
    
    if index_of_closest_peak > 0
        refractory_period_peak_local = index_of_closest_peak;
    end
    
    right_refractory_peak_amp = right_segment_signal(refractory_period_peak_local);
    
    refractory_period_peak = refractory_period_peak_local + lower_bound - 1;
    
    right_refractory_peak = refractory_period_peak;
    
    % figure;
    % plot(signalRaw)
    %
    % figure;
    % x = [1 : length(right_segment_signal)];
    % plot(right_segment_signal); hold on;
    % plot(x(refractory_period_peak_local), right_segment_signal(refractory_period_peak_local), 'o');
    
    upper_bound = refractory_period_peak + gap + blank_period_delineation;
    if upper_bound > length(signal)
        upper_bound = length(signal);
    end
    lower_bound = refractory_period_peak + gap;
    right_segment_signal = signal(lower_bound : upper_bound);
    
    if length(right_segment_signal) > 3
        
        x = [1 : length(right_segment_signal)];
        signalRaw = right_segment_signal;
        right_segment_signal = smooth(x, right_segment_signal, smooth_window, 'loess');
        right_segment_signal = right_segment_signal';
        if is_complex % no smoothing
            %right_segment_signal = signalRaw;
        end
        
        % figure;
        % plot(right_segment_signal)
        %
        % figure;
        % plot(signalRaw)
        
        [~,locsP] = findpeaks(right_segment_signal);
        [~,locsN] = findpeaks(-right_segment_signal);
        
        if ~isempty(locsP)
            if locsP(1) == length(right_segment_signal) - 1
                locsP(1) = [];
            end
        end
        
        if ~isempty(locsN)
            if locsN(1) == length(right_segment_signal) - 1
                locsN(1) = [];
            end
        end
        
        A = [locsP, locsN];
        A = sort(A);
        
        if ~isempty(A)
            offset_amp = right_segment_signal(A(1));
            offset = A(1) + lower_bound - 1;
        else
            diff_2_sig = diff(diff(right_segment_signal));
            [~, diff_2_sig_zeros] = find(abs(diff_2_sig) <= 0.5);
            
            if ~isempty(diff_2_sig_zeros)
                
                diff_1_sig = diff(right_segment_signal);
                diff_1_sig_candidates = diff_1_sig(diff_2_sig_zeros);
                [~, diff_min_loc] = min(abs(diff_1_sig_candidates));
                offset_amp = right_segment_signal(diff_2_sig_zeros(diff_min_loc(1)));
                offset = diff_2_sig_zeros(diff_min_loc(1)) + lower_bound - 1;
                
            else
                if length(diff(right_segment_signal)) > 3
                    [~,locs]=findpeaks(diff(right_segment_signal));
                    if isempty(locs)
                        [~,locs] = max(diff(right_segment_signal,2));
                    else
                        locs = locs(1);
                    end
                    offset_amp = right_segment_signal(locs);
                    offset = locs + lower_bound - 1;
                else
                    [~,locs] = max(diff(right_segment_signal,2));
                    offset_amp = right_segment_signal(locs);
                    offset = locs + lower_bound - 1;
                end
                
            end
            
        end
        
        
        diff_right = abs(right_refractory_peak_amp - offset_amp);
    else
    end
end

end























