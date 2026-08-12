% This subroutine is used to generate and plot the sub-spectra obtained by CCAGGMM,
% and calculate the relative error D, skewness, and kurtosis.
% Inputs:
    % CCAGGMM_para_new - Sub-spectrum parameters obtained by the CCAGGMM algorithm.
% Outputs:
    % D_value - Relative error defined by Eq.(43).
    % CCAGGMM_reconstructed - Reconstructed spectrum by CCAGGMM.
    % CCAGGMM_subspectrum - Sub-spectra obtained by CCAGGMM.
    % skew - Skewness of CCAGGMM sub-spectra.
    % kurt - Kurtosis of CCAGGMM sub-spectra.
% The global variables are defined in the 'spectrum_decomposition_main.m' function.
function [D_value, CCAGGMM_reconstructed, CCAGGMM_subspectrum, skew, kurt] = GP_CCAGGMM(CCAGGMM_para_new)
global dx xline yline length_yline spectrum_number
skew = zeros(1,spectrum_number); % Skewness.
kurt = zeros(1,spectrum_number); % Kurtosis.
%% Calculate sub-spectrum by Eq. (1).
fm = zeros(length_yline, spectrum_number); % Asymmetric generalized Gaussian distributions (AGGD) of sub-spectra.
for i = 1:spectrum_number
    u = CCAGGMM_para_new(2,i); % Expectation.
    al = CCAGGMM_para_new(3,i); % Left branch scale factor.
    ar = CCAGGMM_para_new(4,i); % Right branch scale factor.
    beta = CCAGGMM_para_new(5,i); % Shape factor.
    weight = CCAGGMM_para_new(1,i); % Weight.
    gamma_val = gamma(1/CCAGGMM_para_new(5,i)); % Gamma function.   
    for j = 1:length_yline
        x = xline(j);
        if x < u - 1e-10 % Calculate the left branch of AGGD by Eq. (1).
            fm(j,i) = weight * (beta / ((al + ar) * gamma_val)) * exp(-((u - x) / al)^beta);
        elseif x > u + 1e-10 % Calculate the right branch of AGGD by Eq. (1).
            fm(j,i) = weight * (beta / ((al + ar) * gamma_val)) * exp(-((x - u) / ar)^beta);
        else % Calculate the peak of AGGD.
            fm(j,i) = weight * (beta / ((al + ar) * gamma_val));
        end
    end
end
%% Calculate output spectra.
CCAGGMM_subspectrum = fm*sum(yline(:))*dx; % Sub-spectra obtained by CCAGGMM.
CCAGGMM_reconstructed = sum(fm, 2)*sum(yline(:))*dx; % Reconstructed spectrum obtained by CCAGGMM.
%% Calculate the skewness and kurtosis of each sub-spectrum.
for i = 1:spectrum_number
    [skew(1,i), kurt(1,i)]= Skewkurt(fm(:,i)/CCAGGMM_para_new(1, i), xline, dx);
end
%% Calculate D value by Eq. (43).
D_value = sum(abs(CCAGGMM_reconstructed - yline)) / sum(yline);
fprintf('CCAGGMM Fitting | D value = %.4f\n', D_value);
%% Plotting code.
figure();
plot(xline, yline, '-k', xline, CCAGGMM_reconstructed, '-r');
hold on;
for i = 1:spectrum_number
    plot(xline, CCAGGMM_subspectrum(:,i), '--b');
end
legend('Input Spectrum','CCAGGMM Spectrum','Sub-spectra','Location','best');
xlabel('Sample space', 'FontSize', 12, 'FontName', 'Arial');
ylabel('Proportion', 'FontSize', 12, 'FontName', 'Arial');
set(gca, 'FontSize', 12, 'FontName', 'Arial');
box on;
hold off;
xlim([xline(1), xline(end)]);
ylim([min(ylim), max(ylim)]);
end
% This subroutine calculates the skewness and kurtosis of a discrete probability distribution.
% Input:
    % yline - Proportion of samples.
    % xline - Sample space of sub-spectrum.
% Outputs:
    % skew - Skewness.
    % kurt - Kurtosis.
function [skew, kurt] = Skewkurt(yline, xline, dx)
EX = dx * sum(xline .* yline); % Mean.
xc = xline - EX; % Centralization.
m2 = dx * sum(xc.^2 .* yline); % Second central moment.
m3 = dx * sum(xc.^3 .* yline); % Third central moment.
m4 = dx * sum(xc.^4 .* yline); % Fourth central moment.
sigma3 = m2 ^ 1.5; % σ³.
sigma4 = m2 ^ 2; % σ⁴.
skew = abs(m3 / sigma3); % Skewness.
kurt = abs(m4 / sigma4); % Kurtosis.
end