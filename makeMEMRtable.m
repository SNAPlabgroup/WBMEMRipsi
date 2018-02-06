
resultsdir = 'D:\DATA\WBMEMRipsi\Results\';

csvfile = 'MEMR_BB_results.csv';
fid = fopen(csvfile, 'a');
%% Left Ear
files = dir(strcat(resultsdir, 'S*_L_BB.mat'));
for k = 1:numel(files)
    load(strcat(resultsdir, files(k).name));
    fprintf(fid, '%s, L, %f\n', files(k).name(1:4), MEMresults.cumulative(end));
end

%% Right Ear
files = dir(strcat(resultsdir, 'S*_R_BB.mat'));
for k = 1:numel(files)
    load(strcat(resultsdir, files(k).name));
    fprintf(fid, '%s, R, %f\n', files(k).name(1:4), MEMresults.cumulative(end));
end

fclose(fid);