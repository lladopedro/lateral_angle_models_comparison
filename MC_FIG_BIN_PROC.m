% AUTHOR: Pedro Llado
% SCRIPT TO GENERATE THE PLOT SHOWING THE BINAURAL PROCESSING OUTPUTS

%% PARAMETERISATION
% Presentation level
lvl_dB = 70; % dB SPL
dboffset = 100; % ref = 10uPa --> set as in Breebaart, for the adaptation loop to work

% Middle ear
middleEarFc = [];%[500 2000] in Hz, empty if no BPF applied

% GTFB
fLow = 80; % in Hz
fHigh = 8000; % in Hz
spacingERB = 1; % in ERB

%%
randn('seed',13121);

%% Load SOFA file
%obj = SOFAload("aux_data/HRTFs/HUTUBS/pp55_HRIRs_simulated.sofa");
obj = SOFAload("aux_data/HRTFs/SADIE/D1_HRIR_SOFA/D1_48K_24bit_256tap_FIR_SOFA.sofa");
fs = obj.Data.SamplingRate;

%% Stimulus
noise = randn(fs/4,1);

w = cos(2*pi*25*[1/fs:1/fs:1/100]').^2;
noise(1:length(w)) = noise(1:length(w)) .* w(end:-1:1);
noise(end-length(w)+1:end) = noise(end-length(w)+1:end) .* w;
stim = noise;

stim = scaletodbspl(stim,lvl_dB);
bin_stim = func_binauralise_lsp_signals(stim', 30, 2, obj,fs);


%% LINDEMANN1986
hold on;
subplot(1,5,1)
% run Lindemann model on signal
c_s = 0.3; % stationary inhibition
w_f = 0.035; % monaural sensitivity
M_f = 6; % decrease of monaural sensitivity
T_int = inf; % integration time
N_1 = 1764; % sample at which first cross-correlation is calculated
[crosscorr,t,ild,cfreqs] = lindemann1986(bin_stim,fs,c_s,w_f,M_f,T_int,N_1);

crosscorr = squeeze(crosscorr);

MC_plot_allfreqch(crosscorr-min(crosscorr),cfreqs,'analysis_freqs',[125 250 500 1000 2000 4000],'color',[102,194,165]/256);
yline(1:2:11);
ylim([0.5 12.5])
yticks([1:2:11])
yticklabels([[125 250 500 1000 2000 4000]/1000]);
xlim([0 size(crosscorr,1)+1])
xticks([1:(size(crosscorr,1)-1)/8:size(crosscorr,1)])
xticklabels([-1:1/4:1])
xlabel("Lag (ms)")
ax = gca;
ax.FontSize = 20;
title("lindemann1986")


%% BREEBAART2001
subplot(1,5,2)
[ei_map,cfreq,outsigl,outsigr] = breebaart2001(bin_stim,fs,0,0,'flow',fLow,'fhigh',fHigh);

maxLag = 0.001;
iaccFuncts = zeros(round(2*maxLag*fs+1),length(cfreq));
for freqInd=1:length(cfreq)
    ccg(:,freqInd) = xcorr(outsigl(:,freqInd),...
        outsigr(:,freqInd),round(maxLag*fs),'coeff');
end

MC_plot_allfreqch(ccg - min(ccg),cfreq,'analysis_freqs',[125 250 500 1000 2000 4000],'color',[141,160,203]/256);
yline(1:2:11);
ylim([0.5 12.5])
yticks([1:2:11])
yticklabels([[125 250 500 1000 2000 4000]/1000]);
xlim([0 size(ccg,1)+1])
xticks([1:(size(ccg,1)-1)/8:size(ccg,1)])
xticklabels([-1:1/4:1])
xlabel("Lag (ms)")
ax = gca;
ax.FontSize = 20;
title("breebaart2001")

%% FALLER2004
subplot(1,5,3)
peripheral_out = MC_peripheral_gtfb(bin_stim,fs,fLow,fHigh,spacingERB,...
    'gtfb_order',4,'gtfb_type',"complex",'gtfb_may2011',...
    "false",'gtfb_compression_power',0.23);

peripheral_out = MC_peripheral_neuraltransduction(peripheral_out,'lpf_fc',...
    425,'lpf_order',4,'nt_compression_power',2,...
    'apply_envelope',"false");

%gaussian = randn(length(bin_stim),2);
% 9.4dB@2kHz and the rest of bands scaled according to ISO385
% [spl,freq] = iso226(11.9); %11.9 is hardcoded to guarantee the 9.4dB@2kHz
% spl_cfs = interp1(freq,spl,peripheral_out.cfs);
%gaussian = scaletodbspl(gaussian,11.9);
%gaussian_noise = MC_peripheral_gtfb(gaussian,fs,fLow,fHigh,spacingERB,...
%    'gtfb_order',4,'gtfb_type',"complex",'gtfb_may2011',...
%    "false",'gtfb_compression_power',0.23);

%gaussian_noise = MC_peripheral_neuraltransduction(gaussian_noise,'lpf_fc',...
%    425,'lpf_order',4,'nt_compression_power',2,...
%    'apply_envelope',"false");
%peripheral_out.gaussian_noise = gaussian_noise.ntout;

monaural_out = MC_monauralProcessing(peripheral_out,'mon_method','none');


maxlag_d = 48;% Size of IACC function in number of taps
frame_d = size(bin_stim,1)/4;
frameCount = 1;

ic_threshold = 0.95; % IC THRESHOLD ( 0 <= THETA_X <= 1)
alpha = 10; % EXP WIN TIME CONSTANT ( alpha_f >= 0 ) %in ms
ccg = prec_fallermerimaa(squeeze(monaural_out.ntout(:,1,:))',squeeze(monaural_out.ntout(:,2,:))',[],[],fs,ic_threshold,alpha,maxlag_d,frame_d,frameCount,[],[],[]);  

subplot(1,5,3)

MC_plot_allfreqch(ccg - min(ccg),monaural_out.cfs,'analysis_freqs',[125 250 500 1000 2000 4000],'color',[231,138,195]/256);
yline(1:2:11);
ylim([0.5 12.5])
yticks([1:2:11])
yticklabels([[125 250 500 1000 2000 4000]/1000]);
xlim([0 size(ccg,1)+1])
xticks([1:(size(ccg,1)-1)/8:size(ccg,1)])
xticklabels([-1:1/4:1])
xlabel("Lag (ms)")
ax = gca;
ax.FontSize = 20;
title("faller2004")

%% MAY2011
subplot(1,5,4)
out = may2011pl(bin_stim,fs);% the window size has been manipulated

MC_plot_allfreqch(out.xcorrel-min(out.xcorrel),out.cf,'analysis_freqs',[125 250 500 1000 2000 4000],'color',[166,216,84]/256);
yline(1:2:11);
ylim([0.5 12.5])
yticks([1:2:11])
yticklabels([[125 250 500 1000 2000 4000]/1000]);
xlim([0 size(out.xcorrel,1)+1])
xticks([1:(size(out.xcorrel,1)-1)/8:size(out.xcorrel,1)])
xticklabels([-1:1/4:1])
xlabel("Lag (ms)")
ax = gca;
ax.FontSize = 20;
title("may2011")

%% DESENA2020
subplot(1,5,5)
target = desena2020_estimateinterauralcues(bin_stim,fs);
MC_plot_allfreqch(target.iaccFuncts-min(target.iaccFuncts),target.cfs,'analysis_freqs',[125 250 500 1000 2000 4000],'color',[179,179,179]/256);
yline(1:2:11);
ylim([0.5 12.5])
yticks([1:2:11])
yticklabels([[125 250 500 1000 2000 4000]/1000]);
xlim([0 size(target.iaccFuncts,1)+1])
xticks([1:(size(target.iaccFuncts,1)-1)/8:size(target.iaccFuncts,1)])
xticklabels([-1:1/4:1])
xlabel("Lag (ms)")
ax = gca;
ax.FontSize = 20;
title("desena2020")

%% DIETZ2011
[fine,fc,ild,env,ihc_out,gtfb_out] = dietz2011pl(bin_stim,fs,'mod_center_frequency_hz',500);


figure;
subplot(1,2,1)
MC_plot_allfreqch([abs(fine.itf) abs(env.itf)],fc,'analysis_freqs',[250 500 1000 2000 4000],'color',[255,217,47]/256);

yline(1:2:9);
ylim([-0.5 10.5])
yticks([1:2:11])
yticklabels([[250 500 1000 2000 4000]/1000]);
xlim([0 size(fine.itf,1)+1])
xticks([1:2400:size(fine.itf,1)])
xticklabels(round([1:2400:size(fine.itf,1)]/fs,2)*1000)
xlabel("Time (ms)")
ax = gca;
ax.FontSize = 20;
title("amplitude")

subplot(1,2,2)
MC_plot_allfreqch([angle(fine.itf) angle(env.itf)],fc,'analysis_freqs',[250 500 1000 2000 4000],'normalise_lims',[-pi pi],'color',[255,217,47]/256);
yline(1.5:2:9.5);
ylim([0 11])
yticks([1.5:2:11.5])
yticklabels([[250 500 1000 2000 4000]/1000]);
xlim([0 size(fine.itf,1)+1])
xticks([1:2400:size(fine.itf,1)])
xticklabels(round([1:2400:size(fine.itf,1)]/fs,2)*1000)
xlabel("Time (ms)")
ax = gca;
ax.FontSize = 20;
title("phase")