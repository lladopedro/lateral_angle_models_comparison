% AUTHOR: Pedro Llado
% SCRIPT TO GENERATE THE PLOT SHOWING THE GTFB AND IHC OUTPUTS 

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
noise = randn(fs/50,1);

w = cos(2*pi*25*[1/fs:1/fs:1/100]').^2;
noise(1:length(w)) = noise(1:length(w)) .* w(end:-1:1);
noise(end-length(w)+1:end) = noise(end-length(w)+1:end) .* w;
stim = noise;

stim = scaletodbspl(stim,lvl_dB);%,dboffset);

bin_stim = func_binauralise_lsp_signals(stim', 30, 2, obj,fs);

%% Middle ear filtering

bin_stim = midEarFilt(bin_stim,fs,middleEarFc);

%% LINDEMANN1986
[inoutsig,cfreq] = auditoryfilterbank(bin_stim,fs);

figure(36)
subplot(7,2,1)
MC_plot_allfreqch(squeeze(inoutsig(:,:,1)),cfreq,'analysis_freqs',[500 1000 2000 4000],'color',[102,194,165]/256);
yticklabels("");ylim([0 8.5])

xlim([0 0.06*fs])
xticks([0:0.02*fs:0.04*fs])
xticklabels([0:20:40])
yticklabels([0.5 1 2 4])
xlabel("");("");
ax = gca;
ax.FontSize = 16;

inoutsig = ihcenvelope(inoutsig,fs,'ihc_lindemann1986');

figure(36);
subplot(7,2,2)
MC_plot_allfreqch(squeeze(inoutsig(:,:,1)),cfreq,'analysis_freqs',[500 1000 2000 4000],'color',[102,194,165]/256);
yticklabels("");ylim([0.5 8.5])
xlim([0 0.06*fs])
xticks([0:0.02*fs:0.04*fs])
xticklabels([0:20:40])
yticklabels([0.5 1 2 4])
xlabel("");("");
ax = gca;
ax.FontSize = 16;

%% BREEBAART2001
[~,fc,outsigl,outsigr,outsig_gtfb] = breebaart2001pl(bin_stim,fs,0,0);

figure(36)
subplot(7,2,3)
MC_plot_allfreqch(squeeze(outsig_gtfb(:,:,1)),fc,'analysis_freqs',[500 1000 2000 4000],'color',[141,160,203]/256);
yticklabels("");ylim([0 8.5])
xlim([0 0.06*fs])
xticks([0:0.02*fs:0.04*fs])
xticklabels([0:20:40])
yticklabels([0.5 1 2 4])
xlabel("");("");
ax = gca;
ax.FontSize = 16;

figure(36)
subplot(7,2,4)
MC_plot_allfreqch(outsigl,fc,'analysis_freqs',[500 1000 2000 4000],'color',[141,160,203]/256);
yticklabels("");ylim([0.5 8.5])
xlim([0 0.06*fs])
xticks([0:0.02*fs:0.04*fs])
xticklabels([0:20:40])
yticklabels([0.5 1 2 4])
xlabel("");("");
ax = gca;
ax.FontSize = 16;


%% FALLER2004

peripheral_out_gtfb = MC_peripheral_gtfb(bin_stim,fs,fLow,fHigh,spacingERB,...
    'gtfb_order',4,'gtfb_type',"complex",'gtfb_may2011',...
    "false",'gtfb_compression_power',0.23);

peripheral_out = MC_peripheral_neuraltransduction(peripheral_out_gtfb,'lpf_fc',...
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

monaural_out = MC_monauralProcessing(peripheral_out,'mon_method','none');%'gaussianNoise');

figure(36)
subplot(7,2,5)
MC_plot_allfreqch(squeeze(peripheral_out_gtfb.gtfb(:,1,:)),monaural_out.cfs,'analysis_freqs',[500 1000 2000 4000],'color',[231,138,195]/256);
yticklabels("");ylim([0 8.5])
xlim([0 0.06*fs])
xticks([0:0.02*fs:0.04*fs])
xticklabels([0:20:40])
yticklabels([0.5 1 2 4])
xlabel("");("");
ax = gca;
ax.FontSize = 16;

figure(36)
subplot(7,2,6)
MC_plot_allfreqch(squeeze(monaural_out.bin_input(:,1,:)),monaural_out.cfs,'analysis_freqs',[500 1000 2000 4000],'color',[231,138,195]/256);
yticklabels("");ylim([0.5 8.5])
xlim([0 0.06*fs])
xticks([0:0.02*fs:0.04*fs])
xticklabels([0:20:40])
yticklabels([0.5 1 2 4])
xlabel("");("");
ax = gca;
ax.FontSize = 16;


%% MAY2011
subplot(7,2,7)
[out, ihc_out] = may2011pl(bin_stim,fs);

figure(36)
subplot(7,2,7)
MC_plot_allfreqch(squeeze(out.bm(:,:,1)),peripheral_out.cfs,'analysis_freqs',[500 1000 2000 4000],'color',[166,216,84]/256);

yticklabels("");ylim([0 8.5])
xlim([0 0.06*fs])
xticks([0:0.02*fs:0.04*fs])
xticklabels([0:20:40])
yticklabels([0.5 1 2 4])
xlabel("");("");
ax = gca;
ax.FontSize = 16;

figure(36)
subplot(7,2,8)
MC_plot_allfreqch(squeeze(ihc_out(:,1,:)),peripheral_out.cfs,'analysis_freqs',[500 1000 2000 4000],'color',[166,216,84]/256);

yticklabels("");ylim([0.5 8.5])
xlim([0 0.06*fs])
xticks([0:0.02*fs:0.04*fs])
xticklabels([0:20:40])
yticklabels([0.5 1 2 4])
xlabel("");("");
ax = gca;
ax.FontSize = 16;

%% DIETZ2011
[~,fcs,~,~,ihc_out_dietz,gtfb_out_dietz] = dietz2011pl(bin_stim,fs);

figure(36)
subplot(7,2,9)
MC_plot_allfreqch(squeeze(gtfb_out_dietz(:,:,1)),fcs,'analysis_freqs',[500 1000 2000 4000],'color',[255,217,47]/256);
yticklabels("");ylim([0 8.5])
xlim([0 0.06*fs])
xticks([0:0.02*fs:0.04*fs])
xticklabels([0:20:40])
yticklabels([0.5 1 2 4])
xlabel("");("");
ax = gca;
ax.FontSize = 16;

figure(36)
subplot(7,2,10)
MC_plot_allfreqch(squeeze(ihc_out_dietz(:,:,1)),fcs,'analysis_freqs',[500 1000 2000 4000],'color',[255,217,47]/256);
yticklabels("");ylim([0.5 8.5])
xlim([0 0.06*fs])
xticks([0:0.02*fs:0.04*fs])
xticklabels([0:20:40])
yticklabels([0.5 1 2 4])
xlabel("");("");
ax = gca;
ax.FontSize = 16;


%% TAKANEN2013

output = takanen2013pl_periphery(bin_stim,fs);

figure(36)
subplot(7,2,11)
MC_plot_allfreqch(output.cochlear.velocityLeft,output.fc,'analysis_freqs',[500 1000 2000 4000],'color',[252,141,98]/256);

yticklabels("");ylim([0 8.5])
xlim([0 0.06*fs])
xticks([0:0.02*fs:0.04*fs])
xticklabels([0:20:40])
yticklabels([0.5 1 2 4])
xlabel("");("");
ax = gca;
ax.FontSize = 16;

figure(36)
subplot(7,2,12)
MC_plot_allfreqch(output.left,output.fc,'analysis_freqs',[500 1000 2000 4000],'color',[252,141,98]/256);

yticklabels("");ylim([0.5 8.5])
xlim([0 0.06*fs])
xticks([0:0.02*fs:0.04*fs])
xticklabels([0:20:40])
yticklabels([0.5 1 2 4])
xlabel("");("");
ax = gca;
ax.FontSize = 16;

%% DESENA2020

target = desena2020_estimateinterauralcues(bin_stim,fs);

figure(36)
subplot(7,2,13)
MC_plot_allfreqch(target.gtfb_output.left,target.cfs,'analysis_freqs',[500 1000 2000 4000],'color',[179,179,179]/256);
yticklabels("");ylim([0 8.5])
xlim([0 0.06*fs])
xticks([0:0.02*fs:0.04*fs])
xticklabels([0:20:40])
yticklabels([0.5 1 2 4])
xlabel("");("");
ax = gca;
ax.FontSize = 16;


figure(36)
subplot(7,2,14)
MC_plot_allfreqch(target.nt_output.left,target.cfs,'analysis_freqs',[500 1000 2000 4000],'color',[179,179,179]/256);
yticklabels("");ylim([0.5 8.5])
xlim([0 0.06*fs])
xticks([0:0.02*fs:0.04*fs])
xticklabels([0:20:40])
yticklabels([0.5 1 2 4])
xlabel("");("");
ax = gca;
ax.FontSize = 16;


%%
subplot(7,2,7);
ylabel('Characteristic frequency (kHz)')

subplot(7,2,13);
xlabel('Time (ms)')


subplot(7,2,14);
xlabel('Time (ms)')

for ii = 1:14
    subplot(7,2,ii);
    ax = gca;
    ax.FontSize = 16;
    
    if ii <13
        xticks([])
    end
    if mod(ii,2) == 0
        yticks([])
    end
end
