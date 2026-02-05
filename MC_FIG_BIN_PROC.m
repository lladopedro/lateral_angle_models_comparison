% SCRIPT TO GENERATE THE PLOT SHOWING THE BINAURAL PROCESSING OUTPUTS

%CONVENTION:
%Left is positive angle, negative itd, negative ild.
%Listener translation to the left is positive angle

%% PARAMETERISATION

% Presentation level
lvl_dB = 70; % dB SPL
dboffset = 100; % ref = 10uPa --> set as in Breebaart, for the adaptation loop to work

% Middle ear
middleEarFc = [];%[500 2000];%[] % in Hz, empty if no BPF applied

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
%noise = [zeros(round(fs/1000),1);randn(fs/5,1);zeros(round(fs/1000),1)];
%noise = randn(fs/10,1);
noise = randn(fs/4,1);
%noise = randn(fs/50,1);

w = cos(2*pi*25*[1/fs:1/fs:1/100]').^2;
noise(1:length(w)) = noise(1:length(w)) .* w(end:-1:1);
noise(end-length(w)+1:end) = noise(end-length(w)+1:end) .* w;
stim = noise;

%stim = zeros(length(stim),1);
%stim(100) = 1;
% stim(200) = 1;

stim = scaletodbspl(stim,lvl_dB);%,dboffset);

bin_stim = func_binauralise_lsp_signals(stim', 30, 2, obj,fs);


%% Time dependency vs stationary
% Stationary for the moment

%% Middle ear filtering

%bin_stim = midEarFilt(bin_stim,fs,middleEarFc);

%% Peripheral processing I (frequency selectivity)
% Does the frequency range affect the results?
% Does the number of channels change the results substantially?
% Does compression affect? (Dietz 2011)
% 4th order phase compensated? (May 2011)
% DRNL?


% col_matrix = [%102,194,165; %lindemann
% 252,141,98; %takanen
% 141,160,203; %breebaart
% 231,138,195; %faller
% 166,216,84; %may
% 255,217,47; %dietz
% 229,196,148; %macpherson
% 179,179,179; %desena
% 0, 0, 0]/256; %saddler

col_matrix = [255,217,47; %dietz
    179,179,179; %desena
    166,216,84]/256; %may


%bin_stim = bin_stim(400:1400,:);
%figure; hold on;

%% Precomputed results for the same signal - my Matlab-Python interface doesn't work
load("FIGs/verhulst2012_y1_all.mat")
% verhulst2012's output is at 96000 sample rate
verhulst2012_y1_all = resample(verhulst2012_y1_all,48000,96000);
%MC_plot_allfreqch(verhulst2012_y1_all,[519 1054 2239 4458],'analysis_freqs',[500 1000 2000 4000],'color',[252,141,98]/256);

% %% LOAD precomputed output of Saddler's GTFB
% saddler2024_GTFB_out = load("saddler_GTFB_out.csv");
% MC_plot_allfreqch(saddler2024_GTFB_out',[519 1054 2239 4458],'analysis_freqs',[500 1000 2000 4000],'color',[0 0 0]);
%% LINDEMANN1986
%figure;
hold on;
subplot(1,5,1)
% run Lindemann model on signal
c_s = 0.3;%1;%0.5;%0.3; % stationary inhibition
w_f = 0.035;%0; % monaural sensitivity
M_f = 6; % decrease of monaural sensitivity
T_int = inf;%2; % integration time
N_1 = 1764; % sample at which first cross-correlation is calculated
[crosscorr,t,ild,cfreqs] = lindemann1986(bin_stim,fs,c_s,w_f,M_f,T_int,N_1);

crosscorr = squeeze(crosscorr);

MC_plot_allfreqch(crosscorr-min(crosscorr),cfreqs,'analysis_freqs',[125 250 500 1000 2000 4000],'color',[102,194,165]/256);
yline(1:2:11);%3)
ylim([0.5 12.5])%14.5])
yticks([1:2:11])%13])
yticklabels([[125 250 500 1000 2000 4000]/1000]);% 8000]/1000])
xlim([0 size(crosscorr,1)+1])
xticks([1:(size(crosscorr,1)-1)/8:size(crosscorr,1)])
xticklabels([-1:1/4:1])
xlabel("Lag (ms)")
ax = gca;
ax.FontSize = 20;
title("lindemann1986")

% %MACPHERSON1991
% [inoutsig,cfreq] = auditoryfilterbank(bin_stim,fs);
% inoutsig = ihcenvelope(inoutsig,fs,'ihc_lindemann1986');
% MC_plot_allfreqch(squeeze(inoutsig(:,:,1)),cfreq,'analysis_freqs',[500 1000 2000 4000],'color',[102,194,165]/256);
% ylim([0.5 8.5])
% xlim([0 fs/20])
% xticks([0:fs/100:fs/20])
% xticklabels([0:fs/100:fs/20]*1000/fs)
% xlabel("Time (ms)")
% ax = gca;
% ax.FontSize = 20;

%% BREEBAART2001
subplot(1,5,2)
[ei_map,cfreq,outsigl,outsigr] = breebaart2001(bin_stim,fs,0,0,'flow',fLow,'fhigh',fHigh);

maxLag = 0.001;
iaccFuncts = zeros(round(2*maxLag*fs+1),length(cfreq));
for freqInd=1:length(cfreq)
    ccg(:,freqInd) = xcorr(outsigl(:,freqInd),...
        outsigr(:,freqInd),round(maxLag*fs),'coeff'); %COEFF: 
                                %Normalizes the sequence so that the
                                %autocorrelations at zero lag equal 1
end

MC_plot_allfreqch(ccg - min(ccg),cfreq,'analysis_freqs',[125 250 500 1000 2000 4000],'color',[141,160,203]/256);
yline(1:2:11);%3)
ylim([0.5 12.5])%14.5])
yticks([1:2:11])%13])
yticklabels([[125 250 500 1000 2000 4000]/1000]);% 8000]/1000])
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

gaussian = randn(length(bin_stim),2);
% 9.4dB@2kHz and the rest of bands scaled according to ISO385
% I don't really understand, so I applied what I think should
% be very similar
% [spl,freq] = iso226(11.9); %11.9 is hardcoded to guarantee the 9.4dB@2kHz
% spl_cfs = interp1(freq,spl,peripheral_out.cfs);
gaussian = scaletodbspl(gaussian,11.9);
gaussian_noise = MC_peripheral_gtfb(gaussian,fs,fLow,fHigh,spacingERB,...
    'gtfb_order',4,'gtfb_type',"complex",'gtfb_may2011',...
    "false",'gtfb_compression_power',0.23);

gaussian_noise = MC_peripheral_neuraltransduction(gaussian_noise,'lpf_fc',...
    425,'lpf_order',4,'nt_compression_power',2,...
    'apply_envelope',"false");
peripheral_out.gaussian_noise = gaussian_noise.ntout;

%monaural_out = MC_monauralProcessing(peripheral_out,'mon_method','gaussianNoise');
monaural_out = MC_monauralProcessing(peripheral_out,'mon_method','none');


maxlag_d = 48;%48; % Size of IACC function in number of taps
frame_d = size(bin_stim,1)/4; % 
frameCount = 1;

ic_threshold = 0.95; % IC THRESHOLD ( 0 <= THETA_X <= 1)
%alpha = 0.01; % EXP WIN TIME CONSTANT ( alpha_f >= 0 )
alpha = 10; % EXP WIN TIME CONSTANT ( alpha_f >= 0 ) %in ms
ccg = prec_fallermerimaa(squeeze(monaural_out.ntout(:,1,:))',squeeze(monaural_out.ntout(:,2,:))',[],[],fs,ic_threshold,alpha,maxlag_d,frame_d,frameCount,[],[],[]);  

subplot(1,5,3)

MC_plot_allfreqch(ccg - min(ccg),monaural_out.cfs,'analysis_freqs',[125 250 500 1000 2000 4000],'color',[231,138,195]/256);
yline(1:2:11);%3)
ylim([0.5 12.5])%14.5])
yticks([1:2:11])%13])
yticklabels([[125 250 500 1000 2000 4000]/1000]);% 8000]/1000])
xlim([0 size(ccg,1)+1])
xticks([1:(size(ccg,1)-1)/8:size(ccg,1)])
xticklabels([-1:1/4:1])
xlabel("Lag (ms)")
ax = gca;
ax.FontSize = 20;
title("faller2004")

%% MAY2011
subplot(1,5,4)
out = may2011pl(bin_stim,fs);
% the window size has been manipulated to be 1 window only.
% xcorr = mean(out.xcorr,3)'; average over windows

MC_plot_allfreqch(out.xcorrel-min(out.xcorrel),out.cf,'analysis_freqs',[125 250 500 1000 2000 4000],'color',[166,216,84]/256);
yline(1:2:11);%3)
ylim([0.5 12.5])%14.5])
yticks([1:2:11])%13])
yticklabels([[125 250 500 1000 2000 4000]/1000]);% 8000]/1000])
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
yline(1:2:11);%3)
ylim([0.5 12.5])%14.5])
yticks([1:2:11])%13])
yticklabels([[125 250 500 1000 2000 4000]/1000]);% 8000]/1000])
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
%MC_plot_allfreqch([abs(fine.itf)],fc(1:12),'analysis_freqs',[250 500 1000],'normalise_lims',[0 5e-2],'color',[255,217,47]/256);
%MC_plot_allfreqch([abs(fine.itf) abs(env.itf)],fc,'analysis_freqs',[250 500 1000 2000 4000],'normalise_lims',[-1 1],'color',[255,217,47]/256);
%MC_plot_allfreqch([abs(fine.itf) abs(env.itf)],fc,'analysis_freqs',[250 500 1000 2000 4000],'color',[255,217,47]/256);

yline(1:2:9);%3)
ylim([-0.5 10.5])%14.5])
yticks([1:2:11])%13])
yticklabels([[250 500 1000 2000 4000]/1000]);% 8000]/1000])
xlim([0 size(fine.itf,1)+1])
xticks([1:2400:size(fine.itf,1)])
xticklabels(round([1:2400:size(fine.itf,1)]/fs,2)*1000)
xlabel("Time (ms)")
ax = gca;
ax.FontSize = 20;
title("amplitude")

subplot(1,2,2)
MC_plot_allfreqch([angle(fine.itf) angle(env.itf)],fc,'analysis_freqs',[250 500 1000 2000 4000],'normalise_lims',[-pi pi],'color',[255,217,47]/256);
yline(1.5:2:9.5);%3)
ylim([0 11])%14.5])
yticks([1.5:2:11.5])%13])
yticklabels([[250 500 1000 2000 4000]/1000]);% 8000]/1000])
xlim([0 size(fine.itf,1)+1])
xticks([1:2400:size(fine.itf,1)])
xticklabels(round([1:2400:size(fine.itf,1)]/fs,2)*1000)
xlabel("Time (ms)")
ax = gca;
ax.FontSize = 20;
title("phase")

% figure;
% subplot(1,2,1)
% MC_plot_allfreqch(abs(fine.itf),fc(1:size(fine.itf,2)),'analysis_freqs',[200 400 600 800 1000 1200],'color',[255,217,47]/256);
% yline(1:2:11);%3)
% ylim([0.5 12.5])%14.5])
% yticks([1:2:11])%13])
% yticklabels([[200 400 600 800 1000 1200]/1000]);% 8000]/1000])
% xlim([0 size(fine.itf,1)+1])
% xticks([1:1000:size(fine.itf,1)])
% xticklabels(round([1:1000:size(fine.itf,1)]/fs,2)*1000)
% xlabel("Time (ms)")
% ax = gca;
% ax.FontSize = 20;
% title("amplitude")
% 
% subplot(1,2,2)
% MC_plot_allfreqch(angle(fine.itf),fc(1:size(fine.itf,2)),'analysis_freqs',[200 400 600 800 1000 1200],'color',[255,217,47]/256);
% yline(1:2:11);%3)
% ylim([0.5 12.5])%14.5])
% yticks([1:2:11])%13])
% yticklabels([[200 400 600 800 1000 1200]/1000]);% 8000]/1000])
% xlim([0 size(fine.itf,1)+1])
% xticks([1:1000:size(fine.itf,1)])
% xticklabels(round([1:1000:size(fine.itf,1)]/fs,2)*1000)
% xlabel("Time (ms)")
% ax = gca;
% ax.FontSize = 20;
% title("phase")
% 
% sgtitle("dietz2011 (finestructure)")

%% TAKANEN2013
output = takanen2013pl(bin_stim,fs,1); %ComputationType should be 2, but in the code it seems it needs to be 0.

%output = takanen2013(bin_stim,fs,2);

% %%
% 
% gtfb_order = 4; % default = 4;
% gtfb_type = "complex"; % default = "complex"; "allpole";
% gtfb_may2011 = "true"; % It applies the all-poles 4th order gammatone filterbank from may2011
% %compression power is not used if gtfb_may2011 == true
% gtfb_compression_power = 1; % 0.23@FallerMerimaa; 0.4@Dietz; 1@noCompression
% 
% gtfb_order = 4; % default = 4;
% gtfb_type = "allpole"; % default = "complex"; "allpole";
% gtfb_may2011 = "false"; % It applies the all-poles 4th order gammatone filterbank from may2011
% %compression power is not used if gtfb_may2011 == true
% gtfb_compression_power = 0.4; % 0.23@FallerMerimaa; 0.4@Dietz; 1@noCompression
% 
% gtfb_order = 4; % default = 4;
% gtfb_type = "complex"; % default = "complex"; "allpole";
% gtfb_may2011 = "false"; % It applies the all-poles 4th order gammatone filterbank from may2011
% %compression power is not used if gtfb_may2011 == true
% gtfb_compression_power = 1; % 0.23@FallerMerimaa; 0.4@Dietz; 1@noCompression
% 
% 
% peripheral_out = MC_peripheral_gtfb(bin_stim,fs,fLow,fHigh,spacingERB,...
% 'gtfb_order',gtfb_order,'gtfb_type',gtfb_type,'gtfb_may2011',...
% gtfb_may2011,'gtfb_compression_power',gtfb_compression_power);
% 
% MC_plot_allfreqch(squeeze(peripheral_out.gtfb(:,1,:)),peripheral_out.cfs,'analysis_freqs',[500 1000 2000 4000],'color',col_matrix(i,:));
% xlim([0 fs/20])
% xticks([0:fs/100:fs/20])
% xticklabels([0:fs/100:fs/20]*1000/fs)
% xlabel("Time (ms)")
% ax = gca;
% ax.FontSize = 20;
% legend("Transmission line","","","","4th-ord phase-compensated allpole","","","","4th-ord allpole compression","","","","4th-ord complex",'Location','northoutside')
% 

%%


%%
% [V,Y,OAE,CF]=verhulst2012(bin_stim(:,1),fs,4458,lvl_dB);
% %[V,Y,OAE,CF]=verhulst2012(bin_stim(:,1),fs,[519 1054 2239 4458],lvl_dB);
% %
% %Run amtoolbox-1.6.0/environments/verhulst2012/run_cochlear_model.py
% %%
% y1 = load('y1.csv');
% plot(y1)
% 
% %y1_519 = y1;
% %y1_1054 = y1;
% %y1_2239 = y1;
% %y1_4458 = y1;
% 
% %verhulst2012_y1_all = [y1_519 y1_1054 y1_2239 y1_4458];
