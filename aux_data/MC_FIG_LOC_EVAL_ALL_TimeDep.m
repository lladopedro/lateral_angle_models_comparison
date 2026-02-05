% SCRIPT TO GENERATE THE PLOT SHOWING ESTIMATED ANGLES PER FREQUENCY
% CHANNEL

%TODO
% NORMALISE THE LEVEL FOR THE GAUSSIAN NOISE FROM FALLER
% COMPUTE THE ITD-BASED LOCALISATION ONLY DURING/AFTER A MOMENT OF VERY
% HIGH COHERENCE?


% Include speech.
% Include low-passed and high-passed
% Time-dependent estimates (50ms window, 1ms stepsize)?
% ITD/ILD/FREQ DEPENDENT weights - include Ahrens2020?
% Check what's wrong with Faller's ild2angle.
% How to derive direction for takanen2013 properly
% 
% Does the frequency range affect the results?
% Does the number of channels change the results substantially?
% Does compression affect? (Dietz 2011)
% 4th order phase compensated? (May 2011)
% DRNL?

%CONVENTION:
%Left is positive angle, negative itd, negative ild.
%Listener translation to the left is positive angle

%% WHAT TO PLOT
%do_fig = "plot_FF_whiteNoise";
%do_fig = "plot_VBAP_whiteNoise";
%do_fig = "plot_LeadLag_whiteNoise";

do_fig = "plot_FF_speech";
%do_fig = "plot_VBAP_speech";
%do_fig = "plot_LeadLag_speech";

figure;

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
clear bin_stim_all est_angle est_angle_f est_angle_itd est_angle_ild X_VALUES

if do_fig == "plot_FF_whiteNoise"
    %noise = randn(fs/50,1);
    noise = randn(fs/4,1);

    w = cos(2*pi*25*[1/fs:1/fs:1/100]').^2;
    noise(1:length(w)) = noise(1:length(w)) .* w(end:-1:1);
    noise(end-length(w)+1:end) = noise(end-length(w)+1:end) .* w;
    stim = noise;
    %stim = scaletodbspl(stim,lvl_dB);%,dboffset);
    
    X_VALUES = -90:30:90;
    for idir = 1:length(X_VALUES)
        bin_stim_all(idir,:,:) = func_binauralise_lsp_signals(stim', X_VALUES(idir), 2, obj,fs);
    end
    
    newcolors = flipud(turbo(length(X_VALUES)));
    colororder(newcolors)

elseif do_fig == "plot_VBAP_whiteNoise"
    %noise = randn(fs/50,1);
    noise = randn(fs/4,1);
    w = cos(2*pi*25*[1/fs:1/fs:1/100]').^2;
    noise(1:length(w)) = noise(1:length(w)) .* w(end:-1:1);
    noise(end-length(w)+1:end) = noise(end-length(w)+1:end) .* w;
    stim = noise;
    %stim = scaletodbspl(stim,lvl_dB);%,dboffset);
    
    lsp_angles = [-30 30];
    X_VALUES = -20:5:20; %dB in VBAP
    for idir = 1:length(X_VALUES)
        bin_stim_all(idir,:,:) = func_binauralise_lsp_signals([stim'; stim'*10.^(X_VALUES(idir)/20)], lsp_angles, 2, obj,fs);
    end

    
    newcolors = flipud(cool(length(X_VALUES)));
    colororder(newcolors)

elseif do_fig == "plot_LeadLag_whiteNoise"
    % noise = randn(fs/50,1);
    noise = randn(fs/4,1);
    w = cos(2*pi*25*[1/fs:1/fs:1/100]').^2;
    noise(1:length(w)) = noise(1:length(w)) .* w(end:-1:1);
    noise(end-length(w)+1:end) = noise(end-length(w)+1:end) .* w;
    stim = noise;
    %stim = scaletodbspl(stim,lvl_dB);%,dboffset);
    
    lsp_angles = [30 -30];
    X_VALUES = [0 0.5 1:3:21]; %[0:0.25:2 5:5:20];
    extra_taps = ceil(max(X_VALUES*fs/1000));
    for idir = 1:length(X_VALUES)
        if X_VALUES(idir) == 0
            bin_stim_all(idir,:,:) = func_binauralise_lsp_signals([stim' zeros(1,extra_taps); stim' zeros(1,extra_taps)], lsp_angles, 2, obj,fs);
        else
            bin_stim_all(idir,:,:) = func_binauralise_lsp_signals([stim' zeros(1,extra_taps); zeros(1,round(X_VALUES(idir)*fs/1000)) stim' zeros(1,extra_taps - round(X_VALUES(idir)*fs/1000))], lsp_angles, 2, obj,fs);
        end
    end
    
    newcolors = hot(5+length(X_VALUES));
    colororder(newcolors)

elseif do_fig == "plot_FF_speech"
    [speech,stim_fs] = audioread('aux_data/Stimuli/EBU_FemaleSpeech.wav');
    assert(stim_fs == fs);
    speech = speech(30950:30950+fs-1);
    w = cos(2*pi*25*[1/fs:1/fs:1/100]').^2;
    speech(1:length(w)) = speech(1:length(w)) .* w(end:-1:1);
    speech(end-length(w)+1:end) = speech(end-length(w)+1:end) .* w;
    stim = speech;
    %stim = scaletodbspl(stim,lvl_dB);%,dboffset);
    
    X_VALUES = -90:30:90;
    for idir = 1:length(X_VALUES)
        bin_stim_all(idir,:,:) = func_binauralise_lsp_signals(stim', X_VALUES(idir), 2, obj,fs);
    end

    newcolors = flipud(turbo(length(X_VALUES)));
    colororder(newcolors)
elseif do_fig == "plot_VBAP_speech"
    [speech,stim_fs] = audioread('aux_data/Stimuli/EBU_FemaleSpeech.wav');
    assert(stim_fs == fs);
    speech = speech(30950:30950+fs-1);
    w = cos(2*pi*25*[1/fs:1/fs:1/100]').^2;
    speech(1:length(w)) = speech(1:length(w)) .* w(end:-1:1);
    speech(end-length(w)+1:end) = speech(end-length(w)+1:end) .* w;
    stim = speech;
    %stim = scaletodbspl(stim,lvl_dB);%,dboffset);
    
    lsp_angles = [-30 30];
    X_VALUES = -20:5:20; %dB in VBAP
    for idir = 1:length(X_VALUES)
        bin_stim_all(idir,:,:) = func_binauralise_lsp_signals([stim'; stim'*10.^(X_VALUES(idir)/20)], lsp_angles, 2, obj,fs);
    end

    newcolors = flipud(cool(length(X_VALUES)));
    colororder(newcolors)
    
elseif do_fig == "plot_LeadLag_speech"
    [speech,stim_fs] = audioread('aux_data/Stimuli/EBU_FemaleSpeech.wav');
    assert(stim_fs == fs);
    speech = speech(30950:30950+fs-1);
    w = cos(2*pi*25*[1/fs:1/fs:1/100]').^2;
    speech(1:length(w)) = speech(1:length(w)) .* w(end:-1:1);
    speech(end-length(w)+1:end) = speech(end-length(w)+1:end) .* w;
    stim = speech;
    %stim = scaletodbspl(stim,lvl_dB);%,dboffset);
    
    lsp_angles = [30 -30];
    X_VALUES = [0 0.5 1:3:21]; %[0:0.25:2 5:5:20];
    extra_taps = ceil(max(X_VALUES*fs/1000));
    for idir = 1:length(X_VALUES)
        if X_VALUES(idir) == 0
            bin_stim_all(idir,:,:) = func_binauralise_lsp_signals([stim' zeros(1,extra_taps); stim' zeros(1,extra_taps)], lsp_angles, 2, obj,fs);
        else
            bin_stim_all(idir,:,:) = func_binauralise_lsp_signals([stim' zeros(1,extra_taps); zeros(1,round(X_VALUES(idir)*fs/1000)) stim' zeros(1,extra_taps - round(X_VALUES(idir)*fs/1000))], lsp_angles, 2, obj,fs);
        end
    end
    newcolors = hot(5+length(X_VALUES));
    colororder(newcolors)
end



nStim = size(bin_stim_all,1);

%% Templates
load('MC_TEMPLATES.mat')
% template_lindemann1986 = itd2angle_lookuptable_pl(obj,fs,'lindemann1986');
% %%
% template_breebaart2001 = itd2angle_lookuptable_pl(obj,fs,'breebaart2001');
% %%
% template_faller2004 = itd2angle_lookuptable_pl(obj,fs,'faller2004');
% %%
% template_dietz2011 = itd2angle_lookuptable_pl(obj,fs,'dietz2011');
% %%
% template_desena2020 = desena2020_buildtemplate(obj);

% col_matrix = [%102,194,165; %lindemann
% 252,141,98; %takanen
% 141,160,203; %breebaart
% 231,138,195; %faller
% 166,216,84; %may
% 255,217,47; %dietz
% 229,196,148; %macpherson
% 179,179,179; %desena
% 0, 0, 0]/256; %saddler
%legend("lindemann1986","breebaart2001","faller2004","may2011","dietz2011","takanen2013","desena2020")

%bin_stim = bin_stim(400:1400,:);
%figure; hold on;

%%

window_size = round(fs*0.05); %50ms
hop_size_in_ms = 10;
hop_size = round(fs*hop_size_in_ms/1000); %5ms
Nwindows = ceil((size(bin_stim_all,2)-window_size)/hop_size);


%% NORMALISE LEVEL BASED ON XX dB SPL on each ear for a frontal source
sig_uncorr = func_binauralise_lsp_signals(stim', 0, 2, obj,fs); %uncorrected
sig_uncorr_rms = rms(sig_uncorr);
sig_corr(:,1) = scaletodbspl(sig_uncorr(:,1),lvl_dB);
sig_corr(:,2) = scaletodbspl(sig_uncorr(:,2),lvl_dB);
sig_corr_rms = rms(sig_corr);
norm_factor = sig_corr_rms./sig_uncorr_rms;
%%
%figure;
hold on;
%% LINDEMANN1986
% MODIFIED in wierstorf2013pl to force 5ms and start at 50ms, to have the
% same analysis as in the rest of the models


clear est_angle_time
for idir = 1:nStim
    bin_stim = squeeze(bin_stim_all(idir,:,:));
    
    bin_stim(:,1) = bin_stim(:,1) * norm_factor(1);
    bin_stim(:,2) = bin_stim(:,2) * norm_factor(2);

    % %Remove outliers don't provide meaningful results - only one center frequency left?
    % [~,~,itd,ild,cfreqs,azimuth_timefreq] = wierstorf2013pl_estimateazimuth(bin_stim,template_lindemann1986,'lindemann1986','include_outlier');%,'remove_outlier');
    % 
    % %est_angle_time(:,idir) = nanmedian(azimuth_timefreq,2);
    % est_angle_time(:,idir) = nanmean(azimuth_timefreq,2);

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
            % [v,idx] = max(cc_tmp(:,jj));
            % itd(ii,jj) = tau(idx)/1000;
            itd(jj) = lindemann1986_centroid(cc(:,jj));
        end
    
        est_angle_time(itw,idir) = nanmedian(itd2angle(itd',template_lindemann1986));
    end
    % subplot(2,3,idir);
    % image('CData',azimuth_timefreq','CDataMapping','scaled');
    % colorbar;
    % clim([-70 70])
    % xlabel("Time window")
    % ylabel("Freq. channel")
    % yticks([1:3:15])
    % yticklabels(cfreqs(1:3:15))
    % if do_fig == "plot_LeadLag_whiteNoise";
    %     title("Delay = " + X_VALUES(idir) + "ms")
    % end
end

subplot(1,7,1);
%figure;
%plot([1:6:48*6],est_angle_time,'LineWidth',2);%,'color',[102,194,165]/256
plot([1:hop_size_in_ms:hop_size_in_ms*length(est_angle_time)],est_angle_time,'LineWidth',2);%,'color',[102,194,165]/256
title("lindemann1986")
ax = gca;
ax.FontSize = 14;

% subplot(1,7,1);
% plot(X_VALUES,est_angle,'color',[102,194,165]/256,'LineWidth',3);
% title("lindemann1986")
% ax = gca;
% ax.FontSize = 14;

%% BREEBAART2001
clear est_angle_time

for idir = 1:nStim
    bin_stim = squeeze(bin_stim_all(idir,:,:));
    bin_stim(:,1) = bin_stim(:,1) * norm_factor(1);
    bin_stim(:,2) = bin_stim(:,2) * norm_factor(2);

    [ei_map,cfreq,outsigl,outsigr] = breebaart2001(bin_stim,fs,0,0,'flow',fLow,'fhigh',fHigh);
    
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
                current_window_R(:,freqInd),round(maxLag*fs),'coeff'); %COEFF: 
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
        ild = dbspl(current_window_R) - dbspl(current_window_L);
        est_angle_f = itd2angle(itd,template_breebaart2001);
        est_angle(idir) = nanmedian(est_angle_f);
        %est_angle(idir) = nanmean(est_angle_f);
        est_angle_time(itw,idir) = est_angle(idir);
    end
end

subplot(1,7,2);
%figure;
plot([1:hop_size_in_ms:hop_size_in_ms*length(est_angle_time)],est_angle_time,'LineWidth',2);%,'color',[102,194,165]/256
title("breebaart2001")
ax = gca;
ax.FontSize = 14;

% subplot(1,7,2);
% plot(X_VALUES,est_angle,'color',[141,160,203]/256,'LineWidth',3);
% title("breebaart2001")
% ax = gca;
% ax.FontSize = 14;


%% FALLER2004
% Faller's IC threshold is computed per sample. If there is at least one of
% the samples in a blocksize if IC above the threshold, the IC is the
% average of those samples with IC above the threshold.
clear est_angle_time

for idir = 1:nStim
    bin_stim = squeeze(bin_stim_all(idir,:,:));

    bin_stim(:,1) = bin_stim(:,1) * norm_factor(1);
    bin_stim(:,2) = bin_stim(:,2) * norm_factor(2);

    clear itd ild ccg crosscorr est_angle_itd est_angle_ild
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
    %gaussian(:,1) = gaussian(:,1) * norm_factor(1);
    %gaussian(:,2) = gaussian(:,2) * norm_factor(2);
    gaussian_noise = MC_peripheral_gtfb(gaussian,fs,fLow,fHigh,spacingERB,...
        'gtfb_order',4,'gtfb_type',"complex",'gtfb_may2011',...
        "false",'gtfb_compression_power',0.23);
    
    gaussian_noise = MC_peripheral_neuraltransduction(gaussian_noise,'lpf_fc',...
        425,'lpf_order',4,'nt_compression_power',2,...
        'apply_envelope',"false");
    peripheral_out.gaussian_noise = gaussian_noise.ntout;
    
    %monaural_out = MC_monauralProcessing(peripheral_out,'mon_method','gaussianNoise');
    monaural_out = MC_monauralProcessing(peripheral_out,'mon_method','none');
    
    
%     maxlag_d = 48;%48; % Size of IACC function in number of taps
%     frame_d = size(bin_stim,1)/4; % 
%     frameCount = 1;
% 
%     ic_threshold = 0.95; % IC THRESHOLD ( 0 <= THETA_X <= 1)
%     %alpha = 0.01; % EXP WIN TIME CONSTANT ( alpha_f >= 0 )
%     alpha = 10; % EXP WIN TIME CONSTANT ( alpha_f >= 0 ) %in ms
% 
%     for itw = 1:Nwindows
%         sample_start = 1 + (itw-1)*hop_size;
%         sample_end = window_size + (itw-1)*hop_size;
%         current_window_L = squeeze(monaural_out.ntout(sample_start:sample_end,1,:))';
%         current_window_R = squeeze(monaural_out.ntout(sample_start:sample_end,2,:))';
% 
% 
%         ccg = prec_fallermerimaa(current_window_L,current_window_R,[],[],fs,ic_threshold,alpha,maxlag_d,frame_d,frameCount,[],[],[]);  
% 
%         crosscorr = ccg';
% 
%         % Calculate tau (delay line time) axes
%         tau = linspace(-1,1,size(crosscorr,2));
%         % find max in cc
%         itd = zeros(size(crosscorr,1),size(crosscorr,3));
%         itd_centroid = zeros(size(crosscorr,1),size(crosscorr,3));
%         for ii=1:size(crosscorr,1)
%             %tau_centroid(ii) = lindemann1986_centroid(crosscorr(ii,:));
%             %itd_centroid(ii) = tau_centroid(ii)/1000;
%             for jj=1:size(crosscorr,3)
%                 [~,idx] = max(crosscorr(ii,:,jj));
%                 itd(ii,jj) = tau(idx)/1000;
%             end
%         end
% 
%         %ild = 10*log10(rms(squeeze(monaural_out.bin_input(:,1,:)),1)) - 10*log10(rms(squeeze(monaural_out.bin_input(:,2,:)),1));
%         ild = 10*log10(rms(squeeze(monaural_out.bin_input(:,2,:)),1)) - 10*log10(rms(squeeze(monaural_out.bin_input(:,1,:)),1));
% 
% 
%         est_angle_itd = itd2angle(itd,template_faller2004);
%         est_angle_ild = -ild2angle(ild,template_faller2004)';
% 
% %        est_angle(idir) = nanmean([est_angle_itd; est_angle_ild]);
%         est_angle(idir) = nanmedian([est_angle_itd; est_angle_ild]);
%         est_angle_time(itw,idir) = est_angle(idir);
%     end
% end

    maxlag_d = 48; %48; % Size of IACC function in number of taps
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
            %tau_centroid(ii) = lindemann1986_centroid(crosscorr(ii,:));
            %itd_centroid(ii) = tau_centroid(ii)/1000;
            for jj=1:size(crosscorr,3)
                [~,idx] = max(crosscorr(ii,:,jj));
                itd(ii,jj) = tau(idx)/1000;
            end
        end
        
        %ild = 10*log10(rms(squeeze(monaural_out.bin_input(:,1,:)),1)) - 10*log10(rms(squeeze(monaural_out.bin_input(:,2,:)),1));
        ild = 10*log10(rms(squeeze(monaural_out.bin_input(:,2,:)),1)) - 10*log10(rms(squeeze(monaural_out.bin_input(:,1,:)),1));
    
    
        est_angle_itd = itd2angle(itd,template_faller2004);
        est_angle_ild = -ild2angle(ild,template_faller2004)';
    
    %        est_angle(idir) = nanmean([est_angle_itd; est_angle_ild]);
        est_angle(idir) = nanmedian([est_angle_itd; est_angle_ild]);
        est_angle_time(itw,idir) = -est_angle(idir);
    end
end

subplot(1,7,3);
%figure;
plot([1:hop_size_in_ms:hop_size_in_ms*length(est_angle_time)],est_angle_time,'LineWidth',2);%,'color',[102,194,165]/256
title("faller2004")
ax = gca;
ax.FontSize = 14;

% subplot(1,7,3);
% plot(X_VALUES,est_angle,'color',[231,138,195]/256,'LineWidth',3);
% title("faller2004")
% ax = gca;
% ax.FontSize = 14;

%% MAY2011
clear est_angle_time
nSources = 1;
for idir = 1:nStim
    bin_stim = squeeze(bin_stim_all(idir,:,:));
    bin_stim(:,1) = bin_stim(:,1) * norm_factor(1);
    bin_stim(:,2) = bin_stim(:,2) * norm_factor(2);

    out = may2011pl(bin_stim,fs);

    for itw = 1:size(out.loglik,2)
        out_time.rangeAZ = out.rangeAZ;
        out_time.loglik = squeeze(out.loglik(:,itw,:));
        
        est_angle_time(itw,idir)=-may2011_estazimuthgmm(out_time,'HIST',nSources,0);
    end

end

subplot(1,7,4);
plot([1:hop_size_in_ms:hop_size_in_ms*length(est_angle_time)],est_angle_time,'LineWidth',2);%,'color',[102,194,165]/256
title("may2011")
ax = gca;
ax.FontSize = 14;

% subplot(1,7,4);
% plot(X_VALUES,est_angle,'color',[166,216,84]/256,'LineWidth',3);
% title("may2011")
% ax = gca;
% ax.FontSize = 14;


%% DIETZ2011

clear est_angle_time
for idir = 1:nStim
    bin_stim = squeeze(bin_stim_all(idir,:,:));
    bin_stim(:,1) = bin_stim(:,1) * norm_factor(1);
    bin_stim(:,2) = bin_stim(:,2) * norm_factor(2);
    
    [fine, cfreqs, ild_tmp, env] = dietz2011(bin_stim,fs);

    for itw = 1:Nwindows
        sample_start = 1 + (itw-1)*hop_size;
        sample_end = window_size + (itw-1)*hop_size;
        current_window_L = bin_stim(sample_start:sample_end,1,:);
        current_window_R = bin_stim(sample_start:sample_end,2,:);
        
        %est_angle_time(itw,idir) = wierstorf2013_estimateazimuth([current_window_L current_window_R],template_dietz2011,'dietz2011','remove_outlier');
        %est_angle_time(itw,idir) = wierstorf2013pl_estimateazimuth([current_window_L current_window_R],template_dietz2011,'dietz2011','remove_outlier');

        %[fine, cfreqs, ild_tmp, env] = dietz2011([current_window_L current_window_R],fs);
        % unwrap ITD
        % itd_tmp = dietz2011_unwrapitd(fine.itd,ild_tmp(:,1:12),fine.f_inst,2.5);
        % env_itd_tmp = dietz2011_unwrapitd(env.itd,ild_tmp(:,13:23),env.f_inst,2.5);
        itd_tmp = dietz2011_unwrapitd(fine.itd(sample_start:sample_end,:),ild_tmp(sample_start:sample_end,1:12),fine.f_inst(sample_start:sample_end,:),2.5);
        env_itd_tmp = dietz2011_unwrapitd(env.itd(sample_start:sample_end,:),ild_tmp(sample_start:sample_end,13:23),env.f_inst(sample_start:sample_end,:),2.5);

        
        % calculate the mean about time of the binaural parameters and store
        % them
        itd(1:12) = median(itd_tmp,1);
        itd(13:23) = median(env_itd_tmp,1);
        %ild(ii,:) = median(ild_tmp,1);
        
        est_angle_tmp = itd2angle(itd,template_dietz2011);

        est_angle_time(itw,idir) = nanmedian(itd2angle(itd,template_dietz2011));

    end
end

subplot(1,7,5);
plot([1:hop_size_in_ms:hop_size_in_ms*length(est_angle_time)],est_angle_time,'LineWidth',2);%,'color',[102,194,165]/256
title("dietz2011")
ax = gca;
ax.FontSize = 14;

% subplot(1,7,5);
% plot(X_VALUES,est_angle,'color',[255,217,47]/256,'LineWidth',3);
% title("dietz2011")
% ax = gca;
% ax.FontSize = 14;

%% TAKANEN2013

for idir = 1:nStim
    bin_stim = squeeze(bin_stim_all(idir,:,:));
    output(idir) = takanen2013pl(bin_stim,fs,1,0,0);
end

% 
% if do_fig == 'plot_FF_whiteNoise'
%     load('takanen_LOC_whiteNoise_FF.mat')
% elseif do_fig == 'plot_VBAP_whiteNoise'
%     load('takanen_LOC_whiteNoise_VBAP.mat')
% elseif do_fig == 'plot_LeadLag_whiteNoise'
%     load('takanen_LOC_whiteNoise_LeadLag.mat')
% end

%%
for idir = 1:nStim

    %est_angle_takanen_left(:,idir) = mean(real(output(idir).whereLeft),2);
    %est_angle_takanen_right(:,idir) = -mean(real(output(idir).whereRight),2);


    for itw = 1:Nwindows
        sample_start = 1 + (itw-1)*hop_size;
        sample_end = window_size + (itw-1)*hop_size;

        est_angle_takanen_left(itw,idir) = mean(real(output(idir).whereLeft(sample_start:sample_end,:)),"all");
        est_angle_takanen_right(itw,idir) = mean(real(output(idir).whereRight(sample_start:sample_end,:)),"all");
        est_angle_takanen(itw,idir) = est_angle_takanen_right(itw,idir) - est_angle_takanen_left(itw,idir);

    end
end
% %
% 
% %%
subplot(1,7,6);
%plot([1:hop_size_in_ms:hop_size_in_ms*length(est_angle_time)],est_angle_takanen_left,'LineWidth',2,'LineStyle','--');
%hold on;
%plot([1:hop_size_in_ms:hop_size_in_ms*length(est_angle_time)],est_angle_takanen_right,'LineWidth',2);
plot([1:hop_size_in_ms:hop_size_in_ms*length(est_angle_time)],est_angle_takanen,'LineWidth',2);


% plot(1000/fs:1000/fs:length(est_angle_takanen)*1000/fs,est_angle_takanen','LineWidth',3);
title("takanen2013")
% ax = gca;
% ax.FontSize = 14;

%% DESENA2020
clear est_angle_time
for idir = 1:nStim
    bin_stim = squeeze(bin_stim_all(idir,:,:));
    bin_stim(:,1) = bin_stim(:,1) * norm_factor(1);
    bin_stim(:,2) = bin_stim(:,2) * norm_factor(2);

    for itw = 1:Nwindows
        sample_start = 1 + (itw-1)*hop_size;
        sample_end = window_size + (itw-1)*hop_size;
        current_window_L = bin_stim(sample_start:sample_end,1,:);
        current_window_R = bin_stim(sample_start:sample_end,2,:);

        out = desena2020([current_window_L current_window_R],template_desena2020,fs);
        est_angle_time(itw,idir) = out.estimated_localisation;
    end
end

subplot(1,7,7);
plot([1:hop_size_in_ms:hop_size_in_ms*length(est_angle_time)],est_angle_time,'LineWidth',2);%,'color',[102,194,165]/256
title("desena2020")
ax = gca;
ax.FontSize = 14;

% subplot(1,7,7);
% plot(X_VALUES,est_angle,'color',[179,179,179]/256,'LineWidth',3);
% ylim([-50 50])
% title("desena2020")
% grid on;
% ax = gca;
% ax.FontSize = 14;


%% TAKE CARE OF THE AXES
if do_fig == "plot_FF_whiteNoise"
    for isubplot = 1:7
        subplot(1,7,isubplot)
        ylabel('Est. angle (°)');
        xlabel('Time (ms)')
        xlim([0 300])
        ylim([-100 100])
        xticks([0:100:200])
        yticks([-90:45:90])
        grid on;
        ax = gca;
        ax.FontSize = 22;
        %legend("-90°","-60°","-30°","0°","30°","60°","90°",'Location','eastoutside')
    end



elseif do_fig == "plot_VBAP_whiteNoise"
    for isubplot = 1:7
        subplot(1,7,isubplot)
        xlim([0 300])
        xticks([0:100:2000])
        ylim([-50 50])
        ylabel('Est. angle (°)');
        xlabel('Time (ms)')
        yticks([-45:15:45])
        grid on;
        ax = gca;
        ax.FontSize = 22;
        %legend("-20 dB","-15 dB","-10 dB","-5 dB","0 dB","5 dB","10 dB","15 dB","20 dB",'Location','eastoutside')

    end

elseif do_fig == "plot_LeadLag_whiteNoise"
    for isubplot = 1:7
        subplot(1,7,isubplot)
        xlim([0 300])
        xticks([0:100:200])
        ylabel('Est. angle (°)');
        xlabel('Time (ms)')
        ylim([-50 50])
        yticks([-45:15:45])
        grid on;
        ax = gca;
        ax.FontSize = 22;
    end
    %legend("0 ms","0.5 ms","1 ms","4 ms","7 ms","10 ms","13 ms","16 ms","19 ms",'Location','eastoutside')

elseif do_fig == "plot_FF_speech"
    for isubplot = 1:7
        subplot(1,7,isubplot)
        ylabel('Est. angle (°)');
        xlabel('Time (ms)')
        xlim([0 1100])
        ylim([-100 100])
        xticks([0:500:1000])
        yticks([-90:45:90])
        grid on;
        ax = gca;
        ax.FontSize = 22;
    end

elseif do_fig == "plot_VBAP_speech"
    for isubplot = 1:7
        subplot(1,7,isubplot)
        xlim([0 1100])
        xticks([0:500:1000])
        ylim([-50 50])
        ylabel('Est. angle (°)');
        xlabel('Time (ms)')
        yticks([-45:15:45])
        grid on;
        ax = gca;
        ax.FontSize = 22;
    end

elseif do_fig == "plot_LeadLag_speech"
    for isubplot = 1:7
        subplot(1,7,isubplot)
        xlim([0 1100])
        xticks([0:500:1000])
        ylabel('Est. angle (°)');
        xlabel('Time (ms)')
        ylim([-50 50])
        yticks([-45:15:45])
        grid on;
        ax = gca;
        ax.FontSize = 22;
    end
end





%%
%legend("lindemann1986","breebaart2001","faller2004","may2011","dietz2011","takanen2013","desena2020")
% title("")
% xticks([0:5:20])
% xlim([-1 21])
% xlabel("Lag (ms)")