% SCRIPT TO GENERATE THE PLOT SHOWING THE IHC OUTPUTS

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
noise = randn(fs/50,1);

w = cos(2*pi*25*[1/fs:1/fs:1/100]').^2;
noise(1:length(w)) = noise(1:length(w)) .* w(end:-1:1);
noise(end-length(w)+1:end) = noise(end-length(w)+1:end) .* w;
stim = noise;
% stim = zeros(length(stim),1);
% stim(100) = 1;

stim = scaletodbspl(stim,lvl_dB);%,dboffset);

bin_stim = func_binauralise_lsp_signals(stim', 30, 2, obj,fs);


%% Time dependency vs stationary
% Stationary for the moment

%% Middle ear filtering

bin_stim = midEarFilt(bin_stim,fs,middleEarFc);

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
% load("FIGs/verhulst2012_y1_all.mat")
% % verhulst2012's output is at 96000 sample rate
% verhulst2012_y1_all = resample(verhulst2012_y1_all,48000,96000);
% %MC_plot_allfreqch(verhulst2012_y1_all,[519 1054 2239 4458],'analysis_freqs',[500 1000 2000 4000],'color',[252,141,98]/256);

% %% LOAD precomputed output of Saddler's GTFB
% saddler2024_GTFB_out = load("saddler_GTFB_out.csv");
% MC_plot_allfreqch(saddler2024_GTFB_out',[519 1054 2239 4458],'analysis_freqs',[500 1000 2000 4000],'color',[0 0 0]);
%% LINDEMANN1986
[inoutsig,cfreq] = auditoryfilterbank(bin_stim,fs);

figure(36)
subplot(7,2,1)
MC_plot_allfreqch(squeeze(inoutsig(:,:,1)),cfreq,'analysis_freqs',[500 1000 2000 4000],'color',[102,194,165]/256);
yticklabels("");ylim([0 8.5])
% xlim([0 fs/20])
% xticks([0:fs/100:fs/20])
% xticklabels([0:fs/100:fs/20]*1000/fs)
xlim([0 0.06*fs])
xticks([0:0.02*fs:0.04*fs])
xticklabels([0:20:40])
yticklabels([0.5 1 2 4])
xlabel("");("");%xlabel("Time (ms)")
ax = gca;
ax.FontSize = 16;
%title("lindemann1986")
%text(10,10,"lindemann1986","FontSize",20)

inoutsig = ihcenvelope(inoutsig,fs,'ihc_lindemann1986');

figure(36);
subplot(7,2,2)
MC_plot_allfreqch(squeeze(inoutsig(:,:,1)),cfreq,'analysis_freqs',[500 1000 2000 4000],'color',[102,194,165]/256);
yticklabels("");ylim([0.5 8.5])
% xlim([0 fs/20])
% xticks([0:fs/100:fs/20])
% xticklabels([0:fs/100:fs/20]*1000/fs)
xlim([0 0.06*fs])
xticks([0:0.02*fs:0.04*fs])
xticklabels([0:20:40])
yticklabels([0.5 1 2 4])
xlabel("");("");%xlabel("Time (ms)")
ax = gca;
ax.FontSize = 16;
%title("lindemann1986")

% %MACPHERSON1991
% [inoutsig,cfreq] = auditoryfilterbank(bin_stim,fs);
% inoutsig = ihcenvelope(inoutsig,fs,'ihc_lindemann1986');
% MC_plot_allfreqch(squeeze(inoutsig(:,:,1)),cfreq,'analysis_freqs',[500 1000 2000 4000],'color',[102,194,165]/256);
% yticklabels("");ylim([0.5 8.5])
% xlim([0 fs/20])
% xticks([0:fs/100:fs/20])
% xticklabels([0:fs/100:fs/20]*1000/fs)
% xlabel("");("");%xlabel("Time (ms)")
% ax = gca;
% ax.FontSize = 16;

%% BREEBAART2001
% [inoutsig,cfreq] = auditoryfilterbank(bin_stim,fs);
% inoutsig = ihcenvelope(inoutsig,fs,'ihc_breebaart2001');
[~,fc,outsigl,outsigr,outsig_gtfb] = breebaart2001pl(bin_stim,fs,0,0);

figure(36)
subplot(7,2,3)
MC_plot_allfreqch(squeeze(outsig_gtfb(:,:,1)),fc,'analysis_freqs',[500 1000 2000 4000],'color',[141,160,203]/256);
yticklabels("");ylim([0 8.5])
xlim([0 0.06*fs])
xticks([0:0.02*fs:0.04*fs])
xticklabels([0:20:40])
yticklabels([0.5 1 2 4])
xlabel("");("");%xlabel("Time (ms)")
ax = gca;
ax.FontSize = 16;
%title("breebaart2001")

figure(36)
subplot(7,2,4)
MC_plot_allfreqch(outsigl,fc,'analysis_freqs',[500 1000 2000 4000],'color',[141,160,203]/256);
yticklabels("");ylim([0.5 8.5])
xlim([0 0.06*fs])
xticks([0:0.02*fs:0.04*fs])
xticklabels([0:20:40])
yticklabels([0.5 1 2 4])
xlabel("");("");%xlabel("Time (ms)")
ax = gca;
ax.FontSize = 16;
%title("breebaart2001")

%% FALLER2004

peripheral_out_gtfb = MC_peripheral_gtfb(bin_stim,fs,fLow,fHigh,spacingERB,...
    'gtfb_order',4,'gtfb_type',"complex",'gtfb_may2011',...
    "false",'gtfb_compression_power',0.23);

peripheral_out = MC_peripheral_neuraltransduction(peripheral_out_gtfb,'lpf_fc',...
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

monaural_out = MC_monauralProcessing(peripheral_out,'mon_method','gaussianNoise');

figure(36)
subplot(7,2,5)
MC_plot_allfreqch(squeeze(peripheral_out_gtfb.gtfb(:,1,:)),monaural_out.cfs,'analysis_freqs',[500 1000 2000 4000],'color',[231,138,195]/256);
yticklabels("");ylim([0 8.5])
xlim([0 0.06*fs])
xticks([0:0.02*fs:0.04*fs])
xticklabels([0:20:40])
yticklabels([0.5 1 2 4])
xlabel("");("");%xlabel("Time (ms)")
ax = gca;
ax.FontSize = 16;

%title("faller2004")

figure(36)
subplot(7,2,6)
MC_plot_allfreqch(squeeze(monaural_out.bin_input(:,1,:)),monaural_out.cfs,'analysis_freqs',[500 1000 2000 4000],'color',[231,138,195]/256);
yticklabels("");ylim([0.5 8.5])
xlim([0 0.06*fs])
xticks([0:0.02*fs:0.04*fs])
xticklabels([0:20:40])
yticklabels([0.5 1 2 4])
xlabel("");("");%xlabel("Time (ms)")
ax = gca;
ax.FontSize = 16;

%title("faller2004")

%% MAY2011
subplot(7,2,7)
% peripheral_out = MC_peripheral_gtfb(bin_stim,fs,fLow,fHigh,spacingERB,...
%     'gtfb_order',[],'gtfb_type',[],'gtfb_may2011',...
%     "true",'gtfb_compression_power',1);
% 
% peripheral_out = MC_peripheral_neuraltransduction(peripheral_out,'lpf_fc',...
%     0,'lpf_order',5,'nt_compression_power',2,...
%     'apply_envelope',"false");
% MC_plot_allfreqch(squeeze(peripheral_out.ntout(:,1,:)),peripheral_out.cfs,'analysis_freqs',[500 1000 2000 4000],'color',[166,216,84]/256);


[out, ihc_out] = may2011pl(bin_stim,fs);

figure(36)
subplot(7,2,7)
MC_plot_allfreqch(squeeze(out.bm(:,:,1)),peripheral_out.cfs,'analysis_freqs',[500 1000 2000 4000],'color',[166,216,84]/256);

yticklabels("");ylim([0 8.5])
xlim([0 0.06*fs])
xticks([0:0.02*fs:0.04*fs])
xticklabels([0:20:40])
yticklabels([0.5 1 2 4])
xlabel("");("");%xlabel("Time (ms)")
ax = gca;
ax.FontSize = 16;
%title("may2011")

figure(36)
subplot(7,2,8)
MC_plot_allfreqch(squeeze(ihc_out(:,1,:)),peripheral_out.cfs,'analysis_freqs',[500 1000 2000 4000],'color',[166,216,84]/256);

yticklabels("");ylim([0.5 8.5])
xlim([0 0.06*fs])
xticks([0:0.02*fs:0.04*fs])
xticklabels([0:20:40])
yticklabels([0.5 1 2 4])
xlabel("");("");%xlabel("Time (ms)")
ax = gca;
ax.FontSize = 16;
%title("may2011")

%% DIETZ2011
%out = dietz2011(bin_stim,fs);
[~,fcs,~,~,ihc_out_dietz,gtfb_out_dietz] = dietz2011pl(bin_stim,fs);

figure(36)
subplot(7,2,9)
MC_plot_allfreqch(squeeze(gtfb_out_dietz(:,:,1)),fcs,'analysis_freqs',[500 1000 2000 4000],'color',[255,217,47]/256);
yticklabels("");ylim([0 8.5])
xlim([0 0.06*fs])
xticks([0:0.02*fs:0.04*fs])
xticklabels([0:20:40])
yticklabels([0.5 1 2 4])
xlabel("");("");%xlabel("Time (ms)")
ax = gca;
ax.FontSize = 16;
%title("dietz2011")


figure(36)
subplot(7,2,10)
MC_plot_allfreqch(squeeze(ihc_out_dietz(:,:,1)),fcs,'analysis_freqs',[500 1000 2000 4000],'color',[255,217,47]/256);
yticklabels("");ylim([0.5 8.5])
xlim([0 0.06*fs])
xticks([0:0.02*fs:0.04*fs])
xticklabels([0:20:40])
yticklabels([0.5 1 2 4])
xlabel("");("");%xlabel("Time (ms)")
ax = gca;
ax.FontSize = 16;
%title("dietz2011")




%% TAKANEN2013

output = takanen2013pl_periphery(bin_stim,fs);

% x=amt_load('takanen2013','cochleardelays.mat');
% cochlear.delaysV = x.velocitydelays;
% 
% load("FIGs/verhulst2012_v1_all.mat")
% processed.left = verhulst2012_v1_all;
% fc = [519 1054 2239 4458];
% nVals = size(verhulst2012_v1_all,1);
% nBands = size(verhulst2012_v1_all,2);
% % ------ Half-wave rectification -----------------------------------------
% processed.left= processed.left.*(processed.left >0);
% %periphOutput.ventralLeft = processed.right;
% 
% % ------ Impulse generation ----------------------------------------------
% %search of the local maximas of each half-wave separately for left and 
% %right ear signals
% for chanInd=1:nBands
%     left = processed.left(:,chanInd);
%     %right = processed.right(:,chanInd);
%     processed.left(:,chanInd) = zeros(nVals,1);
%     %processed.right(:,chanInd) = zeros(nVals,1);
% 
%     %start end end points for each half-wave
%     startingPoints = strfind((left'>0),[0 1]);
%     endpoints = [strfind((left'>0),[1 0])+1 length(left)];
%     if (isempty(startingPoints) && ~isempty(endpoints))
%         startingPoints = 1;
%     else
%         if startingPoints(1)>=endpoints(1)
%             startingPoints = [1 startingPoints];
%         end
%     end
% 
%     rmsVals=zeros(size(startingPoints));
%     locations = rmsVals;
%     nsamples = endpoints(1:length(startingPoints))-startingPoints+1;
%     %compute the rms values of each half-wave and position it at the local
%     %maximum
%     for locInd=1:length(startingPoints)
%         rmsVals(locInd)= norm(left(startingPoints(locInd):endpoints(locInd)))/sqrt(nsamples(locInd));
%         [unnecessary, locations(locInd)] = max(left(startingPoints(locInd):endpoints(locInd)));
%     end
%     processed.left(locations+startingPoints-1,chanInd) = rmsVals';
% 
%     % %processing of the right channel
%     % 
%     % %start end end points for each half-wave
%     % startingPoints = strfind((right'>0),[0 1]);
%     % endpoints = [strfind((right'>0),[1 0])+1 length(right)];
%     % if (isempty(startingPoints) && ~isempty(endpoints))
%     %     startingPoints = 1;
%     % else
%     %     if startingPoints(1)>=endpoints(1)
%     %         startingPoints = [1 startingPoints];
%     %     end
%     % end
%     % %compute the rms values of each half-wave and position it at the local
%     % %maximum
%     % rmsVals= zeros(size(startingPoints));locations = rmsVals;
%     % nsamples = endpoints(1:length(startingPoints))-startingPoints+1;
%     % for locInd=1:length(startingPoints)
%     %     rmsVals(locInd)= norm(right(startingPoints(locInd):endpoints(locInd)))/sqrt(nsamples(locInd));
%     %     [unnecessary, locations(locInd)] = max(right(startingPoints(locInd):endpoints(locInd)));
%     % end
%     % processed.right(locations+startingPoints-1,chanInd) = rmsVals';
% 
% end
% % ------ Convolution with Gaussian window-functions ----------------------
% % the width of the Gaussian window in the peripheral hearing model depends
% % on the center frequency
% N = zeros(size(fc));
% N(fc<800) = round((2./fc(fc<800))*fs);
% indMid = find(((800<=fc).*(fc<=2800))==1);N(indMid) = round(0.0024*(0.6+0.4*fc(indMid)./800)*fs);
% N(fc>2800) = round(0.0048*fs);
% % the constant alpha is set to 20
% alpha=20;
% 
% periphOutput.left = zeros(nVals,nBands);periphOutput.right = periphOutput.left;
% for chanInd=1:nBands
%     left = processed.left(:,chanInd);
%     %right = processed.right(:,chanInd);
%     n = -N(chanInd)/2:1:N(chanInd)/2;
%     winFunction = (exp(-0.5*(alpha*n/(N(chanInd)/2)).^2))';
%     %convolution with the gaussian window function
%     left = (1/sum(winFunction))*conv(left,winFunction,'same');
%     %right = (1/sum(winFunction))*conv(right,winFunction,'same');
% 
%     %% compensation for the cochlear model delays
%     periphOutput.left(:,chanInd) = [left(cochlear.delaysV(chanInd)+1:end);zeros(cochlear.delaysV(chanInd),1)];
%     %periphOutput.right(:,chanInd) = [right(cochlear.delaysV(chanInd)+1:end);zeros(cochlear.delaysV(chanInd),1)];
%     %periphOutput.ventralLeft(:,chanInd) = [periphOutput.ventralLeft(cochlear.delaysV(chanInd)+1:end,chanInd);zeros(cochlear.delaysV(chanInd),1)];
%     %periphOutput.ventralRight(:,chanInd) = [periphOutput.ventralRight(cochlear.delaysV(chanInd)+1:end,chanInd);zeros(cochlear.delaysV(chanInd),1)];
% end

figure(36)
subplot(7,2,11)
%MC_plot_allfreqch(verhulst2012_v1_all,[519 1054 2239 4458],'analysis_freqs',[500 1000 2000 4000],'color',[252,141,98]/256);
MC_plot_allfreqch(output.cochlear.velocityLeft,output.fc,'analysis_freqs',[500 1000 2000 4000],'color',[252,141,98]/256);

yticklabels("");ylim([0 8.5])
%xlim([0 fs/20])
xlim([0 0.06*fs])
xticks([0:0.02*fs:0.04*fs])
xticklabels([0:20:40])
yticklabels([0.5 1 2 4])
xlabel("");("");%xlabel("Time (ms)")
ax = gca;
ax.FontSize = 16;
%title("takanen2013")

figure(36)
subplot(7,2,12)
%MC_plot_allfreqch(periphOutput.left,[519 1054 2239 4458],'analysis_freqs',[500 1000 2000 4000],'color',[252,141,98]/256);
MC_plot_allfreqch(output.left,output.fc,'analysis_freqs',[500 1000 2000 4000],'color',[252,141,98]/256);

yticklabels("");ylim([0.5 8.5])
%xlim([0 fs/20])
xlim([0 0.06*fs])
xticks([0:0.02*fs:0.04*fs])
xticklabels([0:20:40])
yticklabels([0.5 1 2 4])
xlabel("");("");%xlabel("Time (ms)")
ax = gca;
ax.FontSize = 16;
%title("takanen2013")

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
xlabel("");("");%xlabel("Time (ms)")
ax = gca;
ax.FontSize = 16;
%title("desena2020")

figure(36)
subplot(7,2,14)
MC_plot_allfreqch(target.nt_output.left,target.cfs,'analysis_freqs',[500 1000 2000 4000],'color',[179,179,179]/256);
yticklabels("");ylim([0.5 8.5])
xlim([0 0.06*fs])
xticks([0:0.02*fs:0.04*fs])
xticklabels([0:20:40])
yticklabels([0.5 1 2 4])
xlabel("");("");%xlabel("Time (ms)")
ax = gca;
ax.FontSize = 16;
%title("desena2020")

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


%% SADDLER2024


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
% xlabel("");("");%xlabel("Time (ms)")
% ax = gca;
% ax.FontSize = 16;
% legend("Transmission line","","","","4th-ord phase-compensated allpole","","","","4th-ord allpole compression","","","","4th-ord complex",'Location','northoutside')
% 
% 
% %%
% 
% 
% %%
% %[V,Y,OAE,CF]=verhulst2012(bin_stim(:,1),fs,'all',lvl_dB);
% [V,Y,OAE,CF]=verhulst2012(bin_stim(:,1),fs,1054,lvl_dB);
% %[V,Y,OAE,CF]=verhulst2012(bin_stim(:,1),fs,[519 1054 2239 4458]',lvl_dB);
% % %
% % %Run amtoolbox-1.6.0/environments/verhulst2012/run_cochlear_model.py
% % %%
% % v1 = load('v1.csv');
% % plot(v1)
% % 
% % %v1_519 = v1;
% % %v1_1054 = v1;
% % %v1_2239 = v1;
% % %v1_4458 = v1;
% % 
% % %verhulst2012_v1_all = [v1_519 v1_1054 v1_2239 v1_4458];
