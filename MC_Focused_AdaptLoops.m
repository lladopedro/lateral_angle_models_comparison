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
% load('MC_TEMPLATES.mat')
% % template_lindemann1986 = itd2angle_lookuptable_pl(obj,fs,'lindemann1986');
% % %%
% % template_breebaart2001 = itd2angle_lookuptable_pl(obj,fs,'breebaart2001');
% %
% template_breebaart2001_dau1997 = itd2angle_lookuptable_pl(obj,fs,'breebaart2001dau1997');
% % %%
% % template_faller2004 = itd2angle_lookuptable_pl(obj,fs,'faller2004');
% % %%
% % template_dietz2011 = itd2angle_lookuptable_pl(obj,fs,'dietz2011');
% % %%
% % template_desena2020 = desena2020_buildtemplate(obj);

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
hop_size_in_ms = 1;%5;
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

        if itw == 1 && idir ==1
            figure(301);
            %subplot(2,1,1)
            hold on;
            plot(1000*[1/fs:1/fs:length(outsigl)/fs],outsigl(:,14),'LineWidth',2,'LineStyle',"-",'Color',[94,60,153]/256); %1kHz
        end
    end
end

%subplot(1,7,2);
figure(302); 
%subplot(2,1,2);
hold on;
newcolors = flipud(turbo(length(X_VALUES)));
colororder(newcolors)
plot([1:hop_size_in_ms:hop_size_in_ms*length(est_angle_time)],est_angle_time,'LineWidth',2,'LineStyle',"-",'color',[94,60,153]/256);
%title("breebaart2001")
ax = gca;
ax.FontSize = 14;

% subplot(1,7,2);
% plot(X_VALUES,est_angle,'color',[141,160,203]/256,'LineWidth',3);
% %title("breebaart2001")
% ax = gca;
% ax.FontSize = 14;

%% BREEBAART2001 ALT - dau1997
clear est_angle_time

for idir = 1:nStim
    bin_stim = squeeze(bin_stim_all(idir,:,:));
    bin_stim(:,1) = bin_stim(:,1) * norm_factor(1);
    bin_stim(:,2) = bin_stim(:,2) * norm_factor(2);

    [ei_map,cfreq,outsigl,outsigr] = breebaart2001(bin_stim,fs,0,0,'flow',fLow,'fhigh',fHigh,'adt_dau1997');
    
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
        est_angle_f = itd2angle(itd,template_breebaart2001_dau1997);
        est_angle(idir) = nanmedian(est_angle_f);
        %est_angle(idir) = nanmean(est_angle_f);
        est_angle_time(itw,idir) = est_angle(idir);
        if itw == 1 && idir ==1
            figure(301); hold on;
            %subplot(2,1,1)
            %subplot(1,3,2); hold on;
            %MC_plot_allfreqch(outsigl,cfreq,'analysis_freqs',[125 250 500 1000 2000 4000],'color',[0 0 0]/256);
            plot(1000*[1/fs:1/fs:length(outsigl)/fs],outsigl(:,14),'LineWidth',2,'LineStyle',"-",'Color',[230,97,1]/256); %1kHz

        end
    end
end

%subplot(1,7,2);
%figure;
figure(302);
%subplot(2,1,2)
plot([1:hop_size_in_ms:hop_size_in_ms*length(est_angle_time)],est_angle_time,'LineWidth',2,'LineStyle',"-",'Color',[230,97,1]/256);%,'color',[102,194,165]/256
%title("breebaart2001")
ax = gca;
ax.FontSize = 14;


% %% BREEBAART2001 ALT - relanoiborra2019
% clear est_angle_time
% 
% for idir = 1:nStim
%     bin_stim = squeeze(bin_stim_all(idir,:,:));
%     bin_stim(:,1) = bin_stim(:,1) * norm_factor(1);
%     bin_stim(:,2) = bin_stim(:,2) * norm_factor(2);
% 
%     [ei_map,cfreq,outsigl,outsigr] = breebaart2001(bin_stim,fs,0,0,'flow',fLow,'fhigh',fHigh,'adt_relanoiborra2019');%'adt_dau1997');
% 
%     clear itd ild crosscorr
%     maxLag = 0.001;
%     for itw = 1:Nwindows
%         sample_start = 1 + (itw-1)*hop_size;
%         sample_end = window_size + (itw-1)*hop_size;
%         current_window_L = outsigl(sample_start:sample_end,:);
%         current_window_R = outsigr(sample_start:sample_end,:);
% 
%         iaccFuncts = zeros(round(2*maxLag*fs+1),length(cfreq));
%         clear crosscorr
%         for freqInd=1:length(cfreq)
%             crosscorr(:,freqInd) = xcorr(current_window_L(:,freqInd),...
%                 current_window_R(:,freqInd),round(maxLag*fs),'coeff'); %COEFF: 
%                                         %Normalizes the sequence so that the
%                                         %autocorrelations at zero lag equal 1
%         end
% 
%         crosscorr = crosscorr';
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
%         ild = dbspl(current_window_R) - dbspl(current_window_L);
%         est_angle_f = itd2angle(itd,template_breebaart2001);
%         est_angle(idir) = nanmedian(est_angle_f);
%         %est_angle(idir) = nanmean(est_angle_f);
%         est_angle_time(itw,idir) = est_angle(idir);
%         if itw == 1 && idir ==5
%             figure(302); hold on;
%             %subplot(1,3,3); hold on;
%             %MC_plot_allfreqch(outsigl,cfreq,'analysis_freqs',[125 250 500 1000 2000 4000],'color',[0 0 0]/256);
%             plot(1000*[1/fs:1/fs:length(outsigl)/fs],outsigl(:,14),'LineWidth',2,'LineStyle',"-",'Color','k'); %1kHz
% 
%         end
%     end
% end
% subplot(1,7,2);
% figure(301);
% subplot(2,1,2)
% plot([1:hop_size_in_ms:hop_size_in_ms*length(est_angle_time)],est_angle_time,'LineWidth',2);%,'color',[102,194,165]/256
% %title("breebaart2001")
% ax = gca;
% ax.FontSize = 14;

%% TAKE CARE OF THE AXES

%subplot(2,1,1)
figure(301)
ylabel('Amplitude');
xlabel('Time (ms)')
xlim([0 200])
%ylim([-100 100])
xticks([0:50:250])
%yticks([-90:45:90])
grid on;
ax = gca;
ax.FontSize = 30;
legend("breebaart2001","dau1997")


%subplot(2,1,2)
figure(302)
ylabel('Est. angle (°)');
xlabel('Time (ms)')
xlim([0 200])
ylim([0 80])
xticks([0:50:250])
yticks([-60:30:60])
grid on;
ax = gca;
ax.FontSize = 30;
legend("breebaart2001","dau1997")
