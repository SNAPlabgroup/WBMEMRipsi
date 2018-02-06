
function runSaveMEManalysis(MEMfile, savestem)


resultsdir = 'D:\DATA\WBMEMRipsi\Results\';
load(MEMfile);

endsamps = ceil(stim.clickwin*stim.Fs*1e-3);

freq = 10.^linspace(log10(200), log10(8000), 1024);
MEMband = [500, 2000];
ind = (freq >= MEMband(1)) & (freq <= MEMband(2));

for k = 1:stim.nLevels
    temp = reshape(squeeze(stim.resp(k, :, 2:end, 1:endsamps)),...
        (stim.nreps-1)*stim.Averages, endsamps);
    resp(k, :) = trimmean(temp, 20, 1); %#ok<*AGROW>
    resp_freq(k, :) = pmtm(resp(k, :), 4, freq, stim.Fs);
    
    %     if k < 3
    %         blevs = 1:k;
    %     else
    %         blevs = (k-2):k;
    %     end
    blevs = k;
    temp2 = squeeze(stim.resp(blevs, :, 1, 1:endsamps));
    
    if(numel(blevs) > 1)
        temp2 = reshape(temp2, size(temp2, 2)*numel(blevs), endsamps);
    end
    
    bline(k, :) = trimmean(temp2, 20, 1);
    bline_freq(k, :) = pmtm(bline(k, :), 4, freq, stim.Fs);
end



if(min(stim.noiseatt) == 6)
    elicitor = 94 - (stim.noiseatt - 6);
else
    elicitor = 94 - stim.noiseatt;
end

MEM = pow2db(resp_freq ./ bline_freq);
semilogx(freq / 1e3, MEM, 'linew', 2);
xlim([0.2, 3]);
ticks = [0.25, 0.5, 1, 2];
set(gca, 'XTick', ticks, 'XTickLabel', num2str(ticks'), 'FontSize', 16);
legend(num2str(elicitor'));
xlabel('Frequency (kHz)', 'FontSize', 16);
ylabel('Ear canal pressure (dB re: Baseline)', 'FontSize', 16);
print(gcf, strcat(resultsdir, savestem, '_spectrum.pdf'), '-dpdf');
close all;

figure;
growth = mean(abs(MEM(:, ind)), 2);
plot(elicitor,  growth, 'o-', 'linew', 2);
hold on;
xlabel('Elicitor Level (dB SPL)', 'FontSize', 16);
ylabel('Cumulative Absorbed Power Change (dB)', 'FontSize', 16);
set(gca,'FontSize', 16);
print(gcf, strcat(resultsdir, savestem, '_integrated.pdf'), '-dpdf');
close all;

MEMresults.spsectrum = MEM;
MEMresults.growth = growth;
MEMresults.cumulative = cumsum(growth); %#ok<*STRNU>
save(strcat(resultsdir, savestem, '.mat'), 'MEMresults');

