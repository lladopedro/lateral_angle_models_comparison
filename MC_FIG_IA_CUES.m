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

target_angles = [-90:30:90];

figure;
newcolors = (turbo(length(target_angles)));
colororder(newcolors)


%% NORMALISE LEVEL BASED ON XX dB SPL on each ear for a frontal source
sig_uncorr = func_binauralise_lsp_signals(stim', 0, 2, obj,fs); %uncorrected
sig_uncorr_rms = rms(sig_uncorr);
sig_corr(:,1) = scaletodbspl(sig_uncorr(:,1),lvl_dB);
sig_corr(:,2) = scaletodbspl(sig_uncorr(:,2),lvl_dB);
sig_corr_rms = rms(sig_corr);
norm_factor = sig_corr_rms./sig_uncorr_rms;


for idir = 1:length(target_angles)

bin_stim = func_binauralise_lsp_signals(stim', target_angles(idir), 2, obj,fs);

bin_stim(:,1) = bin_stim(:,1) * norm_factor(1);
bin_stim(:,2) = bin_stim(:,2) * norm_factor(2);

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
%load("FIGs/verhulst2012_y1_all.mat")
% verhulst2012's output is at 96000 sample rate
%verhulst2012_y1_all = resample(verhulst2012_y1_all,48000,96000);
%MC_plot_allfreqch(verhulst2012_y1_all,[519 1054 2239 4458],'analysis_freqs',[500 1000 2000 4000],'color',[252,141,98]/256);

% %% LOAD precomputed output of Saddler's GTFB
% saddler2024_GTFB_out = load("saddler_GTFB_out.csv");
% MC_plot_allfreqch(saddler2024_GTFB_out',[519 1054 2239 4458],'analysis_freqs',[500 1000 2000 4000],'color',[0 0 0]);
%% LINDEMANN1986
%figure;
%hold on;
%subplot(1,5,1)
% run Lindemann model on signal
c_s = 0.3;%1;%0.5;%0.3; % stationary inhibition
w_f = 0.035;%0; % monaural sensitivity
M_f = 6; % decrease of monaural sensitivity
T_int = inf;%2; % integration time
N_1 = 1764; % sample at which first cross-correlation is calculated
[crosscorr,t,ild,cfreqs] = lindemann1986(bin_stim,fs,c_s,w_f,M_f,T_int,N_1);

%as per wierstorf2013_estimateazimuth.m
crosscorr = squeeze(crosscorr)';
% Calculate tau (delay line time) axes
tau = linspace(-1,1,size(crosscorr,2));
% find max in cc
itd = zeros(size(crosscorr,1),size(crosscorr,3));
itd_centroid = zeros(size(crosscorr,1),size(crosscorr,3));
for ii=1:size(crosscorr,1)
    tau_centroid(ii) = lindemann1986_centroid(crosscorr(ii,:));
    itd_centroid(ii) = tau_centroid(ii)/1000;
    for jj=1:size(crosscorr,3)
        [~,idx] = max(crosscorr(ii,:,jj));
        itd(ii,jj) = tau(idx)/1000;
    end
end

subplot(6,2,1); hold on;
%plot(itd*1000);%,'FaceColor',[1 1 1],'EdgeColor',[102,194,165]/256);
%hold on;
plot(cfreqs,itd_centroid*1000,'LineWidth',2);%,'FaceColor',[102,194,165]/256);
set(gca, 'XScale', 'log'); ylim([-1.05 1.05])
ylabel('ITD (ms)'); xlabel('Frequency (Hz)')
title("ITD - lindemann1986")
ax = gca;
ax.FontSize = 10;

subplot(6,2,2); hold on;
plot(cfreqs,ild,'LineWidth',2);%,'FaceColor',[102,194,165]/256)
set(gca, 'XScale', 'log'); ylim([-25 25])
ylabel('ILD (dB)'); xlabel('Frequency (Hz)')
title("ILD - lindemann1986")
ax = gca;
ax.FontSize = 10;





%% BREEBAART2001
[ei_map,cfreq,outsigl,outsigr] = breebaart2001(bin_stim,fs,0,0,'flow',fLow,'fhigh',fHigh);

clear itd ild crosscorr
maxLag = 0.001;
iaccFuncts = zeros(round(2*maxLag*fs+1),length(cfreq));
for freqInd=1:length(cfreq)
    crosscorr(:,freqInd) = xcorr(outsigl(:,freqInd),...
        outsigr(:,freqInd),round(maxLag*fs),'coeff'); %COEFF: 
                                %Normalizes the sequence so that the
                                %autocorrelations at zero lag equal 1
end

crosscorr = crosscorr';

% Calculate tau (delay line time) axes
tau = linspace(-1,1,size(crosscorr,2));
% find max in cc
itd = zeros(size(crosscorr,1),size(crosscorr,3));
itd_centroid = zeros(size(crosscorr,1),size(crosscorr,3));
for ii=1:size(crosscorr,1)
    %tau_centroid(ii) = lindemann1986_centroid(crosscorr(ii,:));
    %itd_centroid(ii) = tau_centroid(ii)/1000;
    for jj=1:size(crosscorr,3)
        [~,idx] = max(crosscorr(ii,:,jj));
        itd(ii,jj) = tau(idx)/1000;
    end
end

ild = dbspl(outsigr) - dbspl(outsigl);

subplot(6,2,3); hold on;
plot(cfreq,itd*1000,'LineWidth',2);%,'FaceColor',[141,160,203]/256);
set(gca, 'XScale', 'log'); ylim([-1.05 1.05])
ylabel('ITD (ms)'); xlabel('Frequency channel')
title("ITD - breebaart2001")
ax = gca;
ax.FontSize = 10;

subplot(6,2,4); hold on;
plot(cfreq,ild,'LineWidth',2);%,'FaceColor',[141,160,203]/256)
set(gca, 'XScale', 'log'); ylim([-25 25])
ylabel('ILD (dB)'); xlabel('Frequency channel')
title("ILD - breebaart2001")
ax = gca;
ax.FontSize = 10;


%% FALLER2004 %continue here

clear itd ild ccg crosscorr
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

crosscorr = ccg';

% Calculate tau (delay line time) axes
tau = linspace(-1,1,size(crosscorr,2));
% find max in cc
itd = zeros(size(crosscorr,1),size(crosscorr,3));
itd_centroid = zeros(size(crosscorr,1),size(crosscorr,3));
for ii=1:size(crosscorr,1)
    %tau_centroid(ii) = lindemann1986_centroid(crosscorr(ii,:));
    %itd_centroid(ii) = tau_centroid(ii)/1000;
    for jj=1:size(crosscorr,3)
        [~,idx] = max(crosscorr(ii,:,jj));
        itd(ii,jj) = tau(idx)/1000;
    end
end

%ild = dbspl(squeeze(monaural_out.bin_input(:,1,:))) - dbspl(squeeze(monaural_out.bin_input(:,2,:)));
ild = 10*log10(rms(squeeze(monaural_out.bin_input(:,2,:)),1)) - 10*log10(rms(squeeze(monaural_out.bin_input(:,1,:)),1));


subplot(6,2,5); hold on;
plot(peripheral_out.cfs,itd*1000,'LineWidth',2);%,'FaceColor',[231,138,195]/256);
set(gca, 'XScale', 'log'); ylim([-1.05 1.05])
ylabel('ITD (ms)'); xlabel('Frequency channel')
title("ITD - faller2004")
ax = gca;
ax.FontSize = 10;

subplot(6,2,6); hold on;
plot(peripheral_out.cfs,ild,'LineWidth',2);%,'FaceColor',[231,138,195]/256)
set(gca, 'XScale', 'log'); ylim([-25 25])
ylabel('ILD (dB)'); xlabel('Frequency channel')
title("ILD - faller2004")
ax = gca;
ax.FontSize = 10;


%% MAY2011
clear itd ild ccg crosscorr
cfs = 1000*[0.0800,0.1095,0.1418,0.1773,0.2161,0.2586,0.3052,0.3562,0.4121,0.4733,0.5404,0.6139,0.6945,0.7827,0.8794,0.9852,1.1013,1.2284,1.3676,1.5202,1.6873,1.8704,2.0710,2.2907,2.5315,2.7953,3.0842,3.4008,3.7477,4.1276,4.5439,5.0000];

out = may2011(bin_stim,fs);
% the window size has been manipulated to be 1 window only.
% xcorr = mean(out.xcorr,3)'; average over windows

subplot(6,2,7); hold on;
plot(cfs,mean(out.itd,2)*1000,'LineWidth',2);%,'FaceColor',[166,216,84]/256);
set(gca, 'XScale', 'log'); ylim([-1.05 1.05])
ylabel('ITD (ms)'); xlabel('Frequency channel')
title("ITD - may2011")
ax = gca;
ax.FontSize = 10;

subplot(6,2,8); hold on;
plot(cfs,mean(out.ild,2),'LineWidth',2);%,'FaceColor',[166,216,84]/256)
set(gca, 'XScale', 'log'); ylim([-25 25])
ylabel('ILD (dB)'); xlabel('Frequency channel')
title("ILD - may2011")
ax = gca;
ax.FontSize = 10;

%% DIETZ2011
clear itd ild ccg crosscorr
[fine,fc,ild,env] = dietz2011(bin_stim,fs);

subplot(6,2,9); hold on;
%fine.itd(fine.itd > 0.001) = fine.itd(fine.itd > 0.001)*NaN;
%fine.itd(fine.itd < -0.001) = fine.itd(fine.itd < -0.001)*NaN;
plot(fc,[nanmedian(fine.itd,1) nanmedian(env.itd,1)]*1000,'LineWidth',2);%,'FaceColor',[255,217,47]/256);
%plot(fc,[median(fine.itd,1) median(env.itd,1)]*1000,'LineWidth',2);%,'FaceColor',[255,217,47]/256);
%plot(fc,[median(fine.itd(1:round(end:2),:),1) median(env.itd(1:round(end:2),:),1)]*1000,'LineWidth',2);%,'FaceColor',[255,217,47]/256);

set(gca, 'XScale', 'log'); ylim([-1.05 1.05])
ylabel('ITD (ms)'); xlabel('Frequency channel')
title("ITD - dietz2011")
ax = gca;
ax.FontSize = 10;

subplot(6,2,10); hold on;
plot(fc,mean(ild,1),'LineWidth',2);%,'FaceColor',[255,217,47]/256)
set(gca, 'XScale', 'log'); ylim([-25 25])
ylabel('ILD (dB)'); xlabel('Frequency channel')
title("ILD - dietz2011")
ax = gca;
ax.FontSize = 10;
%% DESENA2020

target = desena2020_estimateinterauralcues(bin_stim,fs);

subplot(6,2,11); hold on;
plot(target.cfs,target.itd*1000,'LineWidth',2);%,'FaceColor',[179,179,179]/256);
set(gca, 'XScale', 'log'); ylim([-1.05 1.05])
ylabel('ITD (ms)'); xlabel('Frequency channel')
title("ITD - desena2020")
ax = gca;
ax.FontSize = 10;

subplot(6,2,12); hold on;
semilogx(target.cfs,target.ild,'LineWidth',2);%,'FaceColor',[179,179,179]/256)
set(gca, 'XScale', 'log'); ylim([-25 25])
ylabel('ILD (dB)'); xlabel('Frequency channel')
title("ILD - desena2020")
ax = gca;
ax.FontSize = 10;

end

%% TAKE CARE OF THE AXES
for isubplot = 1:12
    subplot(6,2,isubplot)
    xlim([100 20000])
    
    if mod(isubplot,2) == 0
        ylim([-35 35])
    else
        ylim([-1.2 1.2])
    end
    grid on;
end

%%% TAKANEN2013
%output = takanen2013pl(bin_stim,fs,1); %ComputationType should be 2, but in the code it seems it needs to be 0.

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
