% AUTHOR: Pedro Llado
% SCRIPT TO GENERATE THE PLOT COMPARING DIFFERENT IC THRESHOLDS USING
% FALLER2004

%% WHAT TO PLOT
%do_fig = "plot_FF_whiteNoise"; 
%do_fig = "plot_VBAP_whiteNoise";
%do_fig = "plot_LeadLag_whiteNoise";

do_fig = "plot_FF_speech";
%do_fig = "plot_VBAP_speech";
%do_fig = "plot_LeadLag_speech";


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

    
    X_VALUES = -90:30:90;
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
    % noise = randn(fs/50,1);
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
    
    X_VALUES = 60;%-90:30:90;
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
% template_faller2004 = itd2angle_lookuptable_pl(obj,fs,'faller2004');


%%

window_size = round(fs*0.05); %50ms
hop_size_in_ms = 5;
hop_size = round(fs*hop_size_in_ms/1000); %5ms
Nwindows = ceil((size(bin_stim_all,2)-window_size)/hop_size);



%%
%figure;
hold on;

linsty = [":","--","-","-."];
iplot = 1;

%% NORMALISE LEVEL BASED ON XX dB SPL on each ear for a frontal source
sig_uncorr = func_binauralise_lsp_signals(stim', 0, 2, obj,fs); %uncorrected
sig_uncorr_rms = rms(sig_uncorr);
sig_corr(:,1) = scaletodbspl(sig_uncorr(:,1),lvl_dB);
sig_corr(:,2) = scaletodbspl(sig_uncorr(:,2),lvl_dB);
sig_corr_rms = rms(sig_corr);
norm_factor = sig_corr_rms./sig_uncorr_rms;

%% FALLER2004
clear est_angle_time

figure; hold on;

%newcolors = flipud(turbo(length(X_VALUES)));
newcolors = (turbo(length(X_VALUES)));

colororder(newcolors)

isubplot = 1;
ic_threshold_values = [0.95 0.99 0.999];
for ic_threshold = ic_threshold_values
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
    
        hopsize = fs*0.005;
        frameCount = ceil((length(bin_stim)-(frame_d/1000)*fs)/(hopsize));
    
        
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
                       [ic(ii,jj),idx] = max(crosscorr(ii,:,jj));
                        itd(ii,jj) = tau(idx)/1000;
                    end
                end
                
                ild = 10*log10(rms(squeeze(monaural_out.bin_input(:,2,:)),1)) - 10*log10(rms(squeeze(monaural_out.bin_input(:,1,:)),1));
                
                ic_t(:,itw) = ic(:,jj);
                itd_t(:,itw) = itd(:,jj);
                ild_t(:,itw) = ild';
                    
        end

        est_angle_itd_corr = itd2angle(itd_t',template_faller2004)';
        est_angle_ild_corr = -ild2angle(ild_t',template_faller2004)';

    end
    

    figure(255)
    subplot(1,3,isubplot);
    h = image('CData',est_angle_itd_corr,'CDataMapping','scaled');
    colormap(flipud(turbo))

    set(h, 'AlphaData', ~isnan(est_angle_itd_corr))
    clim([-90 90])
    est_angle_time(idir,:) = nanmedian([est_angle_itd_corr; est_angle_ild_corr]);
    title("IC_{threshold} = "+ ic_threshold_values(isubplot))
    xlim([0 200])
    ylim([0 32])
    xlabel('Time block')
    ylabel('Auditory filter')
    ax = gca;
    ax.FontSize = 16;
    isubplot = isubplot+1;
end