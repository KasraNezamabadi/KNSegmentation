function [EE_interval] = Mean_EE_Interval(qrs_extremums)
sum = 0;
for i = 1 : length(qrs_extremums)-1 
    sum = sum + (qrs_extremums(i+1) - qrs_extremums(i));
end
EE_interval = round(sum / (length(qrs_extremums)-1));
end