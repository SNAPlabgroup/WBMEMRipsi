clear;
close all;
clc;

dirpath = uigetdir('.','Please pick MEM data directory');

subjects = dir(strcat(dirpath,'/S*'));

for ear = 1:2
    if ear == 1
        eartag = '_L';
        earname = 'left';
    else
        eartag = '_R';
        earname = 'right';
    end
    for k = 1:numel(subjects)
        s = subjects(k);
        fprintf(1, 'This is subject %s\n', s.name);
        
        MEMfiles = dir(strcat(dirpath, '/', s.name,...
            '/MEMR_', s.name, eartag, '*mat'));
        
        
        if(numel(MEMfiles) < 2)
            fprintf(2,...
                '%s has less than 2 files for %s ear! Skipping!\n',...
                s.name, earname);
            continue;
        end
        clear times;
        for m = 1:numel(MEMfiles)
            times(m) = MEMfiles(m).datenum; %#ok<SAGROW>
        end
        [~, ind] = sort(times, 'descend');
        HPfile = ind(1);
        BBfile = ind(2);
        fprintf(1, 'High-pass MEMR file is: %s\n', MEMfiles(HPfile).name);
        MEMfile = strcat(dirpath, '/', s.name, '/', MEMfiles(HPfile).name);
        savestem = strcat(s.name, eartag, '_HP');
        runSaveMEManalysis(MEMfile, savestem);
        fprintf(1, 'Broadband MEMR file is: %s\n', MEMfiles(BBfile).name);
        MEMfile = strcat(dirpath, '/', s.name, '/', MEMfiles(BBfile).name);
        savestem = strcat(s.name, eartag, '_BB');
        runSaveMEManalysis(MEMfile, savestem);
    end
end
