% AUTHOR: Pedro Llado
% SCRIPT TO GENERATE THE PLOT COMPARING USING THE HILBERT TRANSFORM VS LPF
% FOR DESENA2020

%% WHAT TO PLOT
do_fig = "plot_FF_whiteNoise";
%do_fig = "plot_VBAP_whiteNoise";
%do_fig = "plot_LeadLag_whiteNoise";

figure(201);

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
end

newcolors = flipud(turbo(length(X_VALUES)));
colororder(newcolors)
nStim = size(bin_stim_all,1);

%% Templates
load('MC_TEMPLATES.mat')
% template_desena2020 = desena2020_buildtemplate(obj);
% %%
template_desena2020_altGTFB = desena2020_altGTFB_LPF_buildtemplate(obj);


%%

window_size = round(fs*0.05); %50ms
hop_size_in_ms = 5;
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
figure;
hold on;

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
        if itw == 1 && idir ==1
            figure(204);
            subplot(2,2,3); hold on;
            MC_plot_allfreqch(out.target.nt_output.left,out.target.cfs,'analysis_freqs',[500 1000 2000 4000],'color',[94,60,153]/256);

        end
        
        if itw == 1
            figure(204);hold on;

            newcolors = flipud(turbo(length(X_VALUES)));
            colororder(newcolors)
            subplot(2,2,1); hold on;

            plot(out.target.itd*1000,'LineWidth',2,'LineStyle',':');
            subplot(2,2,2); hold on;
            plot(out.target.ild,'LineWidth',2,'LineStyle',':');
        end
    end
end

figure(204); hold on;
subplot(2,2,4); hold on;

plot([1:hop_size_in_ms:hop_size_in_ms*length(est_angle_time)],est_angle_time,'LineWidth',2,'LineStyle',':');%,'color',[102,194,165]/256
ax = gca;
ax.FontSize = 14;

%% DESENA2020 ALT

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

        out = desena2020_altGTFB_LPF([current_window_L current_window_R],template_desena2020_altGTFB,fs);
        est_angle_time(itw,idir) = out.estimated_localisation;
        if itw == 1 && idir ==1

            figure(204); 
            subplot(2,2,3); hold on; 

            MC_plot_allfreqch(out.target.nt_output.left,out.target.cfs,'analysis_freqs',[500 1000 2000 4000],'color',[230,97,1]/256);
            xticks([0:480:length(out.target.nt_output.left)])
            xticklabels(([0:480:length(out.target.nt_output.left)])/fs*1000)
            xlabel("Time (ms)")
            yticklabels([500 1000 2000 4000]/1000)
            grid on;
            ax = gca;
            ax.FontSize = 22;
            legend("Hilbert Transform","","","","LPF (750 Hz)",'Location','northoutside')

        end
        if itw == 1
            figure(204);hold on;

            newcolors = flipud(turbo(length(X_VALUES)));
            colororder(newcolors)
            subplot(2,2,1); hold on;
            plot(out.target.itd*1000,'LineWidth',2);
            ax = gca;
            ax.FontSize = 22;
            ylabel("ITD (ms)")
            xlabel("Auditory band")
            legend("Hilbert Transform","","","","","","","LPF (750 Hz)",'Location','northoutside')            
            grid on;

            subplot(2,2,2); hold on;

            plot(out.target.ild,'LineWidth',2);
            ylabel("ILD (dB)")
            xlabel("Auditory band")
            ax = gca;
            ax.FontSize = 22;
            grid on;
            ylim([-40 40])
            yticks([-40:20:40])
            legend("Hilbert Transform","","","","","","","LPF (750 Hz)",'Location','northoutside')

        end
    end
end

figure(204);
subplot(2,2,4); hold on; 
plot([1:hop_size_in_ms:hop_size_in_ms*length(est_angle_time)],est_angle_time,'LineWidth',2);
ylim([-100 100])
yticks([-90:30:90])
legend("Hilbert Transform","","","","","","","LPF (750 Hz)",'Location','northoutside')
ylabel("Angle (°)")
xlabel("Time (ms)")
ax = gca;
ax.FontSize = 22;
grid on;


