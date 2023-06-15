function [] = shaded_patch_significant_timepoints(time_axis, is_significant, color)
%% This function takes a figure containing time series data and creates shaded patches
% over time windows that are deemed to be significant
%
% Inputs:
%
% time_axis: 1-d array denoting the the time_axis of the time series
%
% is_significant: 1-d boolean array that is ON at time points that are
% signifiacnt and OFF otherwise

ybounds = get(gca, 'YLim');
is_significant = reshape(is_significant,1,length(is_significant)); % set vector to be 1xn
time_axis = reshape(time_axis,1,length(time_axis));
if length(time_axis) ~= length(is_significant)
    error('Time_axis and Significance vectors must have same length.')
end
%% define significant time start/stop:
sig_onset = find(diff(is_significant) == 1)+1; % points where significant times start
sig_offset = find(diff(is_significant) == -1); % points where significant times stop
% fix end cases:
if is_significant(1)
    sig_onset = [1 sig_onset];
end
if is_significant(end)
    sig_offset = [sig_offset length(is_significant)];
end
%% Plot
hold on;
for k = 1:length(sig_onset)
    patch([time_axis(sig_onset(k)) time_axis(sig_offset(k)) time_axis(sig_offset(k)) time_axis(sig_onset(k))], [ybounds(1) ybounds(1) ybounds(2) ybounds(2)], color, 'FaceAlpha', 0.1, 'EdgeAlpha',0.1);
end
set(gca, 'YLim', ybounds);
end