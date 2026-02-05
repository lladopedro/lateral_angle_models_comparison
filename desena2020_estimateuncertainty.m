function [Hnorm] = desena2020_estimateuncertainty(distance)

n_angles = length(distance.likelihood_f);

alpha = sum(distance.weightedLikelihood.*cos(2.*[1:n_angles].*2.*pi./n_angles)');
beta = sum(distance.weightedLikelihood.*sin(2.*[1:n_angles].*2.*pi./n_angles)');

H = 1 - sqrt(alpha.^2 + beta.^2); % circular variance

Hnorm = (H - distance.Hmin) / (1 - distance.Hmin); % High Hnorm indicates high uncertainty
end
