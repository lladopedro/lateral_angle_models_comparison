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

elseif do_fig == "plot_FF_speech"
    [speech,stim_fs] = audioread('aux_data/Stimuli/EBU_FemaleSpeech.wav');
    assert(stim_fs == fs);
    speech = speech(30950:30950+fs-1);
    w = cos(2*pi*25*[1/fs:1/fs:1/100]').^2;
    speech(1:length(w)) = speech(1:length(w)) .* w(end:-1:1);
    speech(end-length(w)+1:end) = speech(end-length(w)+1:end) .* w;
    stim = speech;
    %stim = scaletodbspl(stim,lvl_dB);%,dboffset);
    
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
    %stim = scaletodbspl(stim,lvl_dB);%,dboffset);
    
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
end



newcolors = flipud(turbo(length(X_VALUES)));
colororder(newcolors)
nStim = size(bin_stim_all,1);

%% Templates
load('MC_TEMPLATES.mat')
% template_lindemann1986 = itd2angle_lookuptable_pl(obj,fs,'lindemann1986');
% %%
% template_breebaart2001 = itd2angle_lookuptable_pl(obj,fs,'breebaart2001');
%
%template_breebaart2001_dau1997 = itd2angle_lookuptable_pl(obj,fs,'breebaart2001dau1997');
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
        
        maxlag_d = 48; %48; % Size of IACC function in number of taps
        frame_d = 50; % in ms
    
        hopsize = fs*0.005;
        frameCount = ceil((length(bin_stim)-(frame_d/1000)*fs)/(hopsize));
    
        
        %ic_threshold = 0.5; % IC THRESHOLD ( 0 <= THETA_X <= 1)
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
                       [ic(ii,jj),idx] = max(crosscorr(ii,:,jj));
                        itd(ii,jj) = tau(idx)/1000;
                    end
                end
                
                %ild = 10*log10(rms(squeeze(monaural_out.bin_input(:,1,:)),1)) - 10*log10(rms(squeeze(monaural_out.bin_input(:,2,:)),1));
                ild = 10*log10(rms(squeeze(monaural_out.bin_input(:,2,:)),1)) - 10*log10(rms(squeeze(monaural_out.bin_input(:,1,:)),1));
                
                ic_t(:,itw) = ic(:,jj);
                itd_t(:,itw) = itd(:,jj);
                ild_t(:,itw) = ild';
                    
            
            %     est_angle_itd = itd2angle(itd,template_faller2004);
            %     est_angle_ild = -ild2angle(ild,template_faller2004)';
            % 
            % %        est_angle(idir) = nanmean([est_angle_itd; est_angle_ild]);
            %     est_angle(idir) = nanmedian([est_angle_itd; est_angle_ild]);
            %     est_angle_time(itw,idir) = est_angle(idir);
        end

        est_angle_itd_corr = itd2angle(itd_t',template_faller2004)';
        est_angle_ild_corr = -ild2angle(ild_t',template_faller2004)';
            
            % It was already corrected in the function prec_fallermerimaa
            % est_angle_itd = itd2angle(itd_t',template_faller2004)';
            % est_angle_ild = -ild2angle(ild_t',template_faller2004)';
            % 
            % 
            % est_angle_itd_corr = NaN(size(est_angle_itd,1),size(est_angle_itd,2));
            % est_angle_ild_corr = NaN(size(est_angle_ild,1),size(est_angle_ild,2));
            % 
            % clear idx
            % if ic_threshold > 0
            %     idx = ic_t > ic_threshold;
            % else
            %     idx = ones(size(ic_t));
            % end
            % est_angle_itd_corr(idx == 1) = est_angle_itd(idx == 1);
            % est_angle_ild_corr(idx == 1) = est_angle_ild(idx == 1);

    end
    
    % figure;
    % plot([1:hop_size_in_ms:hop_size_in_ms*length(est_angle_time)],est_angle_time,'LineWidth',2,'LineStyle',linsty(iplot));%,'color',[102,194,165]/256
    % title("faller2004")
    % ax = gca;
    % ax.FontSize = 14;
    % iplot = iplot + 1;

    figure(255)
    subplot(1,3,isubplot);
    h = image('CData',est_angle_itd_corr,'CDataMapping','scaled');
    colormap(flipud(turbo))

    set(h, 'AlphaData', ~isnan(est_angle_itd_corr))
    clim([-90 90])
    %est_angle_time(idir,:) = nanmean([est_angle_itd_corr; est_angle_ild_corr]);
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