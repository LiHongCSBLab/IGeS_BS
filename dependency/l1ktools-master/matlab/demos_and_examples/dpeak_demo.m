% DPEAK_DEMO Demo of peak detection.

lxbfile = fullfile(mortarpath, '../data/A10.lxb');
outfolder = fullfile(mortarpath, '../data/');
sample_analytes = [15, 25, 100, 200, 300];
% 
% % plot peak intensity distributions for a few analytes
fprintf ('Plotting distributions...\n');
l1kt_plot_peaks(lxbfile, sample_analytes);
drawnow
% close all

% detect peaks and save stats to a file
fprintf('Detecting peaks...\n');
[pkstats, raw] = l1kt_dpeak(lxbfile, 'out', outfolder, 'showfig', false);
fprintf('Done\n')
