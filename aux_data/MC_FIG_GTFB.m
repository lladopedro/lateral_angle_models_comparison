% CHECK THE MC_FIG_IHC DIRECTLY. BOTH ARE PLOTTED.


%% OLD 
% % SCRIPT TO GENERATE THE PLOT SHOWING GTFB DIFFERENCES
% 
% %CONVENTION:
% %Left is positive angle, negative itd, negative ild.
% %Listener translation to the left is positive angle
% 
% % TODO: make takanen/verhulst work
% % TODO: Parameters for the models!
% % CHECK: 5th order butterworth filter - breebaart's approach is always
% % positive, but not sure how to compute the cutoff frequency.
% % CHECK MONAURAL: ADAPTLOOP -- DOUBLECHECK HOW TO CALCULATE MINLVLS --> for
% % the moment, milvls calculated as in the gaussian noise
% % TODO: Include the itd2lookup + ild2lookup in the decision
% % TBC: (should be ok with the correct SPL)Compression works different for values >1 and <1. How should they be?
% 
% 
% %% PARAMETERISATION
% 
% % Presentation level
% lvl_dB = 70; % dB SPL
% dboffset = 100; % ref = 10uPa --> set as in Breebaart, for the adaptation loop to work
% 
% % Middle ear
% middleEarFc = [];%[500 2000];%[] % in Hz, empty if no BPF applied
% 
% % GTFB
% fLow = 300; % in Hz
% fHigh = 8000; % in Hz
% spacingERB = 1; % in ERB
% 
% %%
% randn('seed',13121);
% 
% %% Load SOFA file
% %obj = SOFAload("aux_data/HRTFs/HUTUBS/pp55_HRIRs_simulated.sofa");
% obj = SOFAload("aux_data/HRTFs/SADIE/D1_HRIR_SOFA/D1_48K_24bit_256tap_FIR_SOFA.sofa");
% fs = obj.Data.SamplingRate;
% 
% %% Stimulus
% %noise = [zeros(round(fs/1000),1);randn(fs/5,1);zeros(round(fs/1000),1)];
% %noise = randn(fs/10,1);
% noise = randn(fs/50,1);
% 
% w = cos(2*pi*25*[1/fs:1/fs:1/100]').^2;
% noise(1:length(w)) = noise(1:length(w)) .* w(end:-1:1);
% noise(end-length(w)+1:end) = noise(end-length(w)+1:end) .* w;
% stim = noise;
% % stim = zeros(length(stim),1);
% % stim(100) = 1;
% 
% stim = scaletodbspl(stim,lvl_dB);%,dboffset);
% 
% bin_stim = func_binauralise_lsp_signals(stim', 30, 2, obj,fs);
% 
% 
% %% Time dependency vs stationary
% % Stationary for the moment
% 
% %% Middle ear filtering
% 
% bin_stim = midEarFilt(bin_stim,fs,middleEarFc);
% 
% %% Peripheral processing I (frequency selectivity)
% % Does the frequency range affect the results?
% % Does the number of channels change the results substantially?
% % Does compression affect? (Dietz 2011)
% % 4th order phase compensated? (May 2011)
% % DRNL?
% 
% 
% % col_matrix = [%102,194,165; %lindemann
% % 252,141,98; %takanen
% % 141,160,203; %breebaart
% % 231,138,195; %faller
% % 166,216,84; %may
% % 255,217,47; %dietz
% % 229,196,148; %macpherson
% % 179,179,179; %desena
% % 0, 0, 0]/256; %saddler
% 
% col_matrix = [255,217,47; %dietz
%     179,179,179; %desena
%     166,216,84]/256; %may
% 
% 
% %bin_stim = bin_stim(400:1400,:);
% figure; hold on;
% 
% %% Precomputed results for the same signal - my Matlab-Python interface doesn't work
% load("FIGs/verhulst2012_v1_all.mat")
% % verhulst2012's output is at 96000 sample rate
% verhulst2012_v1_all = resample(verhulst2012_v1_all,48000,96000);
% MC_plot_allfreqch(verhulst2012_v1_all,[519 1054 2239 4458],'analysis_freqs',[500 1000 2000 4000],'color',[252,141,98]/256);
% 
% % %% LOAD precomputed output of Saddler's GTFB
% % saddler2024_GTFB_out = load("saddler_GTFB_out.csv");
% % MC_plot_allfreqch(saddler2024_GTFB_out',[519 1054 2239 4458],'analysis_freqs',[500 1000 2000 4000],'color',[0 0 0]);
% %%
% for i = 1:3
%     %figure;
%     % if i == 1 % 
%     %     model_name = "desena2020";
%     %     gtfb_order = 3; % default = 4;
%     %     gtfb_type = "complex"; % default = "complex"; "allpole";
%     %     gtfb_may2011 = "false"; % It applies the all-poles 4th order gammatone filterbank from may2011
%     %     %compression power is not used if gtfb_may2011 == true
%     %     gtfb_compression_power = 1; % 0.23@FallerMerimaa; 0.4@Dietz; 1@noCompression
%     if i == 1
%         gtfb_order = 4; % default = 4;
%         gtfb_type = "complex"; % default = "complex"; "allpole";
%         gtfb_may2011 = "true"; % It applies the all-poles 4th order gammatone filterbank from may2011
%         %compression power is not used if gtfb_may2011 == true
%         gtfb_compression_power = 1; % 0.23@FallerMerimaa; 0.4@Dietz; 1@noCompression
% 
%     elseif i == 2
%         gtfb_order = 4; % default = 4;
%         gtfb_type = "allpole"; % default = "complex"; "allpole";
%         gtfb_may2011 = "false"; % It applies the all-poles 4th order gammatone filterbank from may2011
%         %compression power is not used if gtfb_may2011 == true
%         gtfb_compression_power = 0.4; % 0.23@FallerMerimaa; 0.4@Dietz; 1@noCompression
%     elseif i == 3
%         gtfb_order = 4; % default = 4;
%         gtfb_type = "complex"; % default = "complex"; "allpole";
%         gtfb_may2011 = "false"; % It applies the all-poles 4th order gammatone filterbank from may2011
%         %compression power is not used if gtfb_may2011 == true
%         gtfb_compression_power = 1; % 0.23@FallerMerimaa; 0.4@Dietz; 1@noCompression
%     end
% 
%     peripheral_out = MC_peripheral_gtfb(bin_stim,fs,fLow,fHigh,spacingERB,...
%         'gtfb_order',gtfb_order,'gtfb_type',gtfb_type,'gtfb_may2011',...
%         gtfb_may2011,'gtfb_compression_power',gtfb_compression_power);
% 
%     MC_plot_allfreqch(squeeze(peripheral_out.gtfb(:,1,:)),peripheral_out.cfs,'analysis_freqs',[500 1000 2000 4000],'color',col_matrix(i,:));
%     xlim([0 fs/20])
%     xticks([0:fs/100:fs/20])
%     xticklabels([0:fs/100:fs/20]*1000/fs)
%     xlabel("Time (ms)")
%     ax = gca;
%     ax.FontSize = 20;
%     legend("Transmission line","","","","4th-ord phase-compensated allpole","","","","4th-ord allpole compression","","","","4th-ord complex",'Location','northoutside')
% end
% 
% %%
% 
% 
% % %%
% % [V,Y,OAE,CF]=verhulst2012(bin_stim(:,1),fs,519,lvl_dB);
% % %[V,Y,OAE,CF]=verhulst2012(bin_stim(:,1),fs,[519 1054 2239 4458],lvl_dB);
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
