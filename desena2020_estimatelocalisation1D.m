function [estimated_localisation] = desena2020_estimatelocalisation1D(distance,template)

%plot(distance.weightedLikelihood);

[~,argmax_itd] = max(distance.likelihood_f_itd);
estimated_localisation_f_itd = template.azimuths(argmax_itd);

[~,argmax_ild] = max(distance.likelihood_f_ild);
estimated_localisation_f_ild = template.azimuths(argmax_ild);

p = 0.7; %norm

% Freq. dependent (as explained by Enzo for Localisation, but not in the paper)
%threshold = 0.6;
threshold = 1;

%f_threshold = 1500;
f_threshold = 15000;

%threshold_id = find(template.cfs>f_threshold,1,"first");
weighting_function_itd_ild = [ interp1([1 length(template.cfs)],[0 threshold], 1:length(template.cfs)) ]; % 0 -> only ITD; 1 -> only ILD

estimated_localisation = mean(estimated_localisation_f_itd.*(1 - weighting_function_itd_ild') + estimated_localisation_f_ild .* weighting_function_itd_ild');%.^(1/p);


end
