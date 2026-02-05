function [Hmin] = desena2020_minimumuncertainty(sofa_obj)

obj = sofa_obj;
%% 1) Load SOFA file and build template
obj = SOFAload("aux_data/HRTFs/HUTUBS/pp55_HRIRs_simulated.sofa");
template = desena2020_buildtemplate(obj);
%plot_desena2020_template(template);

%% 2) Generate test noise
randn('seed',013120);
fs = obj.Data.SamplingRate;
noise = randn(fs/4,1); % 250-ms noise
w = cos(2*pi*25*[1/fs:1/fs:1/100]').^2; % Generate a 10-ms cosine window for ramp on/off
noise(1:length(w)) = noise(1:length(w)) .* w(end:-1:1);
noise(end-length(w)+1:end) = noise(end-length(w)+1:end) .* w;
stim = noise;

%% 3) Binauralise test signal
test_direction = [240,60]; % [Azimuth, elevation]
testID = SOFAfind(obj,test_direction(1),test_direction(2));

test_bin_stim = zeros(length(stim)+length(obj.Data.IR(1,1,:))-1,2);
test_bin_stim(:,1) = conv(stim',squeeze(obj.Data.IR(testID,1,:)));
test_bin_stim(:,2) = conv(stim',squeeze(obj.Data.IR(testID,2,:)));

%% 4) Estimate localisation using desena2020
out = desena2020(test_bin_stim,template,fs);
plot_desena2020_template(template,'target',out.target);


end