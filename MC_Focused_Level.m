% AUTHOR: Pedro Llado
% SCRIPT TO GENERATE THE PLOT COMPARING DIFFERENT PRESENTATION LEVELS USING
% BREEBAART2001


%% WHAT TO PLOT
do_fig = "plot_FF_whiteNoise";
%do_fig = "plot_VBAP_whiteNoise";
%do_fig = "plot_LeadLag_whiteNoise";

%do_fig = "plot_FF_speech";
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
obj = SOFAload("aux_data/HRTFs/SADIE/D1_HRIR_SOFA/D1_48K_24bit_256tap_FIR_SOFA.sofa");
fs = obj.Data.SamplingRate;

%% Stimulus
clear bin_stim_all est_angle est_angle_f est_angle_itd est_angle_ild X_VALUES

if do_fig == "plot_FF_whiteNoise"
    noise = randn(fs/4,1);

    w = cos(2*pi*25*[1/fs:1/fs:1/100]').^2;
    noise(1:length(w)) = noise(1:length(w)) .* w(end:-1:1);
    noise(end-length(w)+1:end) = noise(end-length(w)+1:end) .* w;
    stim = noise;
    
    X_VALUES = 60;%-90:30:90;
    for idir = 1:length(X_VALUES)
        bin_stim_all(idir,:,:) = func_binauralise_lsp_signals(stim', X_VALUES(idir), 2, obj,fs);
    end

elseif do_fig == "plot_VBAP_whiteNoise"
    noise = randn(fs/4,1);
    w = cos(2*pi*25*[1/fs:1/fs:1/100]').^2;
    noise(1:length(w)) = noise(1:length(w)) .* w(end:-1:1);
    noise(end-length(w)+1:end) = noise(end-length(w)+1:end) .* w;
    stim = noise;
    
    lsp_angles = [-30 30];
    X_VALUES = -20:5:20; %dB in VBAP
    for idir = 1:length(X_VALUES)
        bin_stim_all(idir,:,:) = func_binauralise_lsp_signals([stim'; stim'*10.^(X_VALUES(idir)/20)], lsp_angles, 2, obj,fs);
    end

elseif do_fig == "plot_LeadLag_whiteNoise"
    noise = randn(fs/4,1);
    w = cos(2*pi*25*[1/fs:1/fs:1/100]').^2;
    noise(1:length(w)) = noise(1:length(w)) .* w(end:-1:1);
    noise(end-length(w)+1:end) = noise(end-length(w)+1:end) .* w;
    stim = noise;
    
    lsp_angles = [30 -30];
    X_VALUES = [0 0.5 1:3:21]; 
    extra_taps = ceil(max(X_VALUES*fs/1000));
    for idir = 1:length(X_VALUES)
        if X_VALUES(idir) == 0
            bin_stim_all(idir,:,:) = func_binauralise_lsp_signals([stim' zeros(1,extra_taps); stim' zeros(1,extra_taps)], lsp_angles, 2, obj,fs);
        else
            bin_stim_all(idir,:,:) = func_binauralise_lsp_signals([stim' zeros(1,extra_taps); zeros(1,round(X_VALUES(idir)*fs/1000)) stim' zeros(1,extra_taps - round(X_VALUES(idir)*fs/1000))], lsp_angles, 2, obj,fs);
        end
    end

elseif do_fig == "plot_FF_speech"
    [speech,stim_fs] = audioread('aux_data/Stimuli/EBU_FemaleSpeech.wav');
    assert(stim_fs == fs);
    speech = speech(30950:30950+fs-1);
    w = cos(2*pi*25*[1/fs:1/fs:1/100]').^2;
    speech(1:length(w)) = speech(1:length(w)) .* w(end:-1:1);
    speech(end-length(w)+1:end) = speech(end-length(w)+1:end) .* w;
    stim = speech;
    
    X_VALUES = -90:30:90;
    for idir = 1:length(X_VALUES)
        bin_stim_all(idir,:,:) = func_binauralise_lsp_signals(stim', X_VALUES(idir), 2, obj,fs);
    end

elseif do_fig == "plot_VBAP_speech"
    [speech,stim_fs] = audioread('aux_data/Stimuli/EBU_FemaleSpeech.wav');
    assert(stim_fs == fs);
    speech = speech(30950:30950+fs-1);
    w = cos(2*pi*25*[1/fs:1/fs:1/100]').^2;
    speech(1:length(w)) = speech(1:length(w)) .* w(end:-1:1);
    speech(end-length(w)+1:end) = speech(end-length(w)+1:end) .* w;
    stim = speech;
    
    lsp_angles = [-30 30];
    X_VALUES = -20:5:20; %dB in VBAP
    for idir = 1:length(X_VALUES)
        bin_stim_all(idir,:,:) = func_binauralise_lsp_signals([stim'; stim'*10.^(X_VALUES(idir)/20)], lsp_angles, 2, obj,fs);
    end
elseif do_fig == "plot_LeadLag_speech"
    [speech,stim_fs] = audioread('aux_data/Stimuli/EBU_FemaleSpeech.wav');
    assert(stim_fs == fs);
    speech = speech(30950:30950+fs-1);
    w = cos(2*pi*25*[1/fs:1/fs:1/100]').^2;
    speech(1:length(w)) = speech(1:length(w)) .* w(end:-1:1);
    speech(end-length(w)+1:end) = speech(end-length(w)+1:end) .* w;
    stim = speech;
    
    lsp_angles = [30 -30];
    X_VALUES = [0 0.5 1:3:21]; 
    extra_taps = ceil(max(X_VALUES*fs/1000));
    for idir = 1:length(X_VALUES)
        if X_VALUES(idir) == 0
            bin_stim_all(idir,:,:) = func_binauralise_lsp_signals([stim' zeros(1,extra_taps); stim' zeros(1,extra_taps)], lsp_angles, 2, obj,fs);
        else
            bin_stim_all(idir,:,:) = func_binauralise_lsp_signals([stim' zeros(1,extra_taps); zeros(1,round(X_VALUES(idir)*fs/1000)) stim' zeros(1,extra_taps - round(X_VALUES(idir)*fs/1000))], lsp_angles, 2, obj,fs);
        end
    end
end



newcolors = flipud(turbo(length(X_VALUES)));
colororder(newcolors)
nStim = size(bin_stim_all,1);

%% Templates
load('MC_TEMPLATES.mat')

% template_breebaart2001 = itd2angle_lookuptable_pl(obj,fs,'breebaart2001');
%%

window_size = round(fs*0.05); %50ms
hop_size_in_ms = 2;
hop_size = round(fs*hop_size_in_ms/1000); %5ms
Nwindows = ceil((size(bin_stim_all,2)-window_size)/hop_size);



%%
linsty = [":","-.","--","-"];
iplot = 1;

color_levels = [194,230,153;
                120,198,121;
                49,163,84;
                0,104,55]/256;
icolorlevel = 1;
icolorlevel2 = 1;

for lvl_dB = [30:20:90] % dB SPL

    %% NORMALISE LEVEL BASED ON XX dB SPL on each ear for a frontal source
    sig_uncorr = func_binauralise_lsp_signals(stim', 0, 2, obj,fs); %uncorrected
    sig_uncorr_rms = rms(sig_uncorr);
    sig_corr(:,1) = scaletodbspl(sig_uncorr(:,1),lvl_dB);
    sig_corr(:,2) = scaletodbspl(sig_uncorr(:,2),lvl_dB);
    sig_corr_rms = rms(sig_corr);
    norm_factor = sig_corr_rms./sig_uncorr_rms;

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
            est_angle_f = itd2angle(itd,template_breebaart2001);
            est_angle(idir) = nanmedian(est_angle_f);
            est_angle_time(itw,idir) = est_angle(idir);
    
            if itw == 1 && idir ==1
                figure(302); hold on;
                plot(1000*[1/fs:1/fs:length(outsigl)/fs],outsigl(:,14),'LineWidth',2,'Color',color_levels(icolorlevel,:)); %1 kHz
                icolorlevel = icolorlevel+1;
            end
        end
    end
    
    figure(301); hold on;
    plot([1:hop_size_in_ms:hop_size_in_ms*length(est_angle_time)],est_angle_time,'LineWidth',2,'Color',color_levels(icolorlevel2,:));
    icolorlevel2 = icolorlevel2+1;
    ax = gca;
    ax.FontSize = 14;
    legend("30 dB SPL","50 dB SPL","70 dB SPL","90 dB SPL")
    
    iplot = iplot+1;
end

%% TAKE CARE OF THE AXES
if do_fig == "plot_FF_whiteNoise"
        figure(301)
        ylabel('Est. angle (°)');
        xlabel('Time (ms)')
        xlim([0 300])
        ylim([0 75])
        xticks([0:100:200])
        yticks([0:30:60])
        grid on;
        ax = gca;
        ax.FontSize = 22;

        figure(302)
        ylabel('Amplitude');
        xlabel('Time (ms)')
        xlim([0 300])
        xticks([0:100:200])
        grid on;
        ax = gca;
        ax.FontSize = 22;
        legend("30 dB SPL","50 dB SPL","70 dB SPL","90 dB SPL");
 
end
