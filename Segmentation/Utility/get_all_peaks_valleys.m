function [ A ] = get_all_peaks_valleys( segment )

A = [];

if length(segment) > 3
    
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
end


end

