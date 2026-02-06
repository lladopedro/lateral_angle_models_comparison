% AUTHOR: Pedro Llado
% SCRIPT TO GENERATE THE PLOT SHOWING THE FREQUENCY-DEPENDENT LOCALISATION 

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
obj = SOFAload("aux_data/HRTFs/SADIE/D1_HRIR_SOFA/D1_48K_24bit_256tap_FIR_SOFA.sofa");
fs = obj.Data.SamplingRate;

%% Stimulus
noise = randn(fs/4,1);

w = cos(2*pi*25*[1/fs:1/fs:1/100]').^2;
noise(1:length(w)) = noise(1:length(w)) .* w(end:-1:1);
noise(end-length(w)+1:end) = noise(end-length(w)+1:end) .* w;
stim = noise;

stim = scaletodbspl(stim,lvl_dB);%,dboffset);
figure;
target_dir = [-90 -60 -30 0 30 60 90];
newcolors = flipud(turbo(length(target_dir)));
colororder(newcolors)

%% Templates
load('MC_TEMPLATES.mat')
% Compute templates once (uncomment below) and store

% %%
% template_lindemann1986 = itd2angle_lookuptable_pl(obj,fs,'lindemann1986');
% % %%
% template_breebaart2001 = itd2angle_lookuptable_pl(obj,fs,'breebaart2001');
% % %%
% template_faller2004 = itd2angle_lookuptable_pl(obj,fs,'faller2004');
% % %%
% template_dietz2011 = itd2angle_lookuptable_pl(obj,fs,'dietz2011');
% % % %%
% template_desena2020 = desena2020_buildtemplate(obj);


%% NORMALISE LEVEL BASED ON XX dB SPL on each ear for a frontal source
sig_uncorr = func_binauralise_lsp_signals(stim', 0, 2, obj,fs); %uncorrected
sig_uncorr_rms = rms(sig_uncorr);
sig_corr(:,1) = scaletodbspl(sig_uncorr(:,1),lvl_dB);
sig_corr(:,2) = scaletodbspl(sig_uncorr(:,2),lvl_dB);
sig_corr_rms = rms(sig_corr);
norm_factor = sig_corr_rms./sig_uncorr_rms;


for idir = 1:length(target_dir)
    bin_stim = func_binauralise_lsp_signals(stim', target_dir(idir), 2, obj,fs);
    bin_stim(:,1) = bin_stim(:,1) * norm_factor(1);
    bin_stim(:,2) = bin_stim(:,2) * norm_factor(2);
    
    %% LINDEMANN1986
    % Calculate binaural parameters
    c_s = 0.3; % stationary inhibition
    w_f = 0.035; % monaural sensitivity
    M_f = 6; % decrease of monaural sensitivity
    T_int = 10; % integration time
    N_1 = 2400; % sample at which first cross-correlation is calculated

    [cc_tmp,dummy,ild,cfreqs] = lindemann1986(bin_stim,fs,c_s,w_f,M_f,T_int,N_1);
    clear dummy;

    for itw = 1:size(cc_tmp,1)
        cc = squeeze(cc_tmp(itw,:,:));
        for jj = 1:size(cc,2)
            itd(jj) = lindemann1986_centroid(cc(:,jj));
        end
    
        est_angle_time_freq(:,itw) = itd2angle(itd',template_lindemann1986);
    end

est_angle = nanmedian(est_angle_time_freq,2);

subplot(4,2,1)
plot(cfreqs,est_angle,'LineWidth',2);
set(gca, 'XScale', 'log');
hold on;
ylabel('Est. angle (°)'); xlabel('Frequency (Hz)')
title("lindemann1986")
ax = gca;
ax.FontSize = 14;

%% BREEBAART2001
[ei_map,cfreq,outsigl,outsigr] = breebaart2001(bin_stim,fs,0,0,'flow',fLow,'fhigh',fHigh);
    

window_size = round(fs*0.05); %50ms
hop_size_in_ms = 10;
hop_size = round(fs*hop_size_in_ms/1000); %5ms
Nwindows = ceil((size(bin_stim,1)-window_size)/hop_size);


clear itd ild crosscorr
maxLag = 0.001;
for itw = 1:Nwindows
    sample_start = 1 + (itw-1)*hop_size;
    sample_end = window_size + (itw-1)*hop_size;
    current_window_L = outsigl(sample_start:sample_end,:);
    current_window_R = outsigr(sample_start:sample_end,:);

    iaccFuncts = zeros(round(2*maxLag*fs+1),length(cfreq));
    clear crosscorr
    for freqInd=1:length(cfreq)
        crosscorr(:,freqInd) = xcorr(current_window_L(:,freqInd),...
            current_window_R(:,freqInd),round(maxLag*fs),'coeff');
    end
    
    crosscorr = crosscorr';
    
    % Calculate tau (delay line time) axes
    tau = linspace(-1,1,size(crosscorr,2));
    % find max in cc
    itd = zeros(size(crosscorr,1),size(crosscorr,3));
    itd_centroid = zeros(size(crosscorr,1),size(crosscorr,3));
    for ii=1:size(crosscorr,1)
        for jj=1:size(crosscorr,3)
            [~,idx] = max(crosscorr(ii,:,jj));
            itd(ii,jj) = tau(idx)/1000;
        end
    end
    ild = dbspl(current_window_R) - dbspl(current_window_L);
    est_angle_f_t(:,idir,itw) = itd2angle(itd,template_breebaart2001);
end

est_angle = squeeze(nanmedian(est_angle_f_t,3));

subplot(4,2,2)
plot(cfreq,est_angle,'LineWidth',2);
set(gca, 'XScale', 'log');
hold on;
ylabel('Est. angle (°)'); xlabel('Frequency (Hz)')
title("breebaart2001")
ax = gca;
ax.FontSize = 14;


%% FALLER2004

clear itd ild ccg crosscorr
peripheral_out = MC_peripheral_gtfb(bin_stim,fs,fLow,fHigh,spacingERB,...
    'gtfb_order',4,'gtfb_type',"complex",'gtfb_may2011',...
    "false",'gtfb_compression_power',0.23);

peripheral_out = MC_peripheral_neuraltransduction(peripheral_out,'lpf_fc',...
    425,'lpf_order',4,'nt_compression_power',2,...
    'apply_envelope',"false");

% gaussian = randn(length(bin_stim),2);
% 9.4dB@2kHz and the rest of bands scaled according to ISO385
% [spl,freq] = iso226(11.9); %11.9 is hardcoded to guarantee the 9.4dB@2kHz
% spl_cfs = interp1(freq,spl,peripheral_out.cfs);
% gaussian = scaletodbspl(gaussian,11.9);
% gaussian_noise = MC_peripheral_gtfb(gaussian,fs,fLow,fHigh,spacingERB,...
%     'gtfb_order',4,'gtfb_type',"complex",'gtfb_may2011',...
%     "false",'gtfb_compression_power',0.23);
% 
% gaussian_noise = MC_peripheral_neuraltransduction(gaussian_noise,'lpf_fc',...
%     425,'lpf_order',4,'nt_compression_power',2,...
%     'apply_envelope',"false");
% peripheral_out.gaussian_noise = gaussian_noise.ntout;


monaural_out = MC_monauralProcessing(peripheral_out,'mon_method','none');

maxlag_d = 48; % Size of IACC function in number of taps
frame_d = 50; % in ms

hopsize = fs*0.001*hop_size_in_ms;
frameCount = ceil((length(bin_stim)-(frame_d/1000)*fs)/(hopsize));

ic_threshold = 0.95; % IC THRESHOLD ( 0 <= THETA_X <= 1)
alpha = 10; % EXP WIN TIME CONSTANT ( alpha_f >= 0 ) %in ms

ccg = prec_fallermerimaa(squeeze(monaural_out.ntout(:,1,:))',squeeze(monaural_out.ntout(:,2,:))',[],[],fs,ic_threshold,alpha,maxlag_d,frame_d,frameCount,[],[],[]);  

for itw = 1:size(ccg,3)
    crosscorr = ccg(:,:,itw)';
    
    % Calculate tau (delay line time) axes
    tau = linspace(-1,1,size(crosscorr,2));
    % find max in cc
    itd = zeros(size(crosscorr,1),size(crosscorr,3));
    itd_centroid = zeros(size(crosscorr,1),size(crosscorr,3));
    for ii=1:size(crosscorr,1)
        for jj=1:size(crosscorr,3)
            [~,idx] = max(crosscorr(ii,:,jj));
            itd(ii,jj) = tau(idx)/1000;
        end
    end
    
    ild = 10*log10(rms(squeeze(monaural_out.bin_input(:,2,:)),1)) - 10*log10(rms(squeeze(monaural_out.bin_input(:,1,:)),1));


    est_angle_itd_t(:,itw) = itd2angle(itd,template_faller2004);
    est_angle_ild_t(:,itw) = -ild2angle(ild,template_faller2004)';
end
est_angle_itd = nanmedian(est_angle_itd_t,2);
est_angle_ild = nanmedian(est_angle_ild_t,2);

subplot(4,2,3)
plot(peripheral_out.cfs,est_angle_itd,'LineWidth',2);
set(gca, 'XScale', 'log');
hold on;
ylabel('Est. angle (°)'); xlabel('Frequency (Hz)')
title("faller2004 (ITD)")
ax = gca;
ax.FontSize = 14;

subplot(4,2,4)
plot(peripheral_out.cfs,est_angle_ild,'LineWidth',2);
set(gca, 'XScale', 'log');
hold on;
ylabel('Est. angle (°)'); xlabel('Frequency (Hz)')
title("faller2004 (ILD)")
ax = gca;
ax.FontSize = 14;

%% MAY2011
cfs = 1000*[0.0800,0.1095,0.1418,0.1773,0.2161,0.2586,0.3052,0.3562,0.4121,0.4733,0.5404,0.6139,0.6945,0.7827,0.8794,0.9852,1.1013,1.2284,1.3676,1.5202,1.6873,1.8704,2.0710,2.2907,2.5315,2.7953,3.0842,3.4008,3.7477,4.1276,4.5439,5.0000];
out = may2011(bin_stim,fs);

subplot(4,2,5)
plot(cfs,-median(out.azimuth,2),'LineWidth',2);
set(gca, 'XScale', 'log');

hold on;
ylabel('Est. angle (°)'); xlabel('Frequency (Hz)')
title("may2011")
ax = gca;
ax.FontSize = 14;

%% DIETZ2011
[fine,fc,ild,env] = dietz2011(bin_stim,fs);

est_angle = itd2angle([fine.itd env.itd],template_dietz2011);

subplot(4,2,6)
plot(fc,nanmedian(est_angle,1),'LineWidth',2);
set(gca, 'XScale', 'log');
hold on;
ylabel('Est. angle (°)'); xlabel('Frequency (Hz)')
title("dietz2011")
ax = gca;
ax.FontSize = 14;


%% TAKANEN2013
% 
% output = takanen2013pl(bin_stim,fs,1,0,0);
% 
% est_angle_takanen_left = mean(real(output.whereLeft),1);
% est_angle_takanen_right = mean(real(output.whereRight),1);
% est_angle_takanen = est_angle_takanen_right - est_angle_takanen_left;
% 
% 
% subplot(4,2,7)
% hold on;
% plot(output.fc,est_angle_takanen','LineWidth',2,'LineStyle','-','Color',newcolors(idir,:));
% 
% set(gca, 'XScale', 'log');
% ylabel('Est. angle (°)'); xlabel('Frequency (Hz)')
% title("takanen2013")
% ax = gca;
% ax.FontSize = 14;



%% DESENA2020
% Compute interaural cues for the binaural stimulus
target = desena2020_estimateinterauralcues(bin_stim,fs);
% Compute distance between the target and each direction in the template
distance = desena2020_distance(target,template_desena2020);
% Weight frequency band by loudness
distance = desena2020_loudnessweighting(distance,target);
% Estimate localisation
[~,argmax] = max(distance.likelihood_f);
estimated_localisation_f = template_desena2020.azimuths(argmax);


subplot(4,2,8)
plot(target.cfs,estimated_localisation_f,'LineWidth',2);
set(gca, 'XScale', 'log');

hold on;
ylabel('Est. angle (°)'); xlabel('Frequency (Hz)')
title("desena2020")
ax = gca;
ax.FontSize = 14;

end

%% TAKE CARE OF THE AXES
for isubplot = 1:8
    subplot(4,2,isubplot)
    xlim([80 20000])
    ylim([-100 100])
    yticks([-90:30:90])
    grid on;
end

