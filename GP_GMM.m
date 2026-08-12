% This subroutine is used to generate and plot the sub-spectra obtained by GMM,
% and calculate the relative error D.
% Inputs:
    % GMM_para_new - Sub-spectrum parameters obtained by the GMM algorithm.
% Outputs:
    % D_value - Relative error defined by Eq.(43).
    % GMM_reconstructed - Reconstructed spectrum by GMM.
    % GMM_subspectrum - Sub-spectra obtained by GMM.
    % skew - Skewness of GMM sub-spectra.
    % kurt - Kurtosis of GMM sub-spectra.
function [D_value, GMM_reconstructed, GMM_subspectrum] = GP_GMM(GMM_para_new)
global length_yline spectrum_number xline yline dx
skew = zeros(1,spectrum_number); % Skewness.
kurt = zeros(1,spectrum_number); % Kurtosis.
%% Calculate sub-spectra by Gaussian distribution function.
fm = zeros(length_yline, spectrum_number); % Gaussian distributions (GD) of sub-spectra.
for M = 1:length_yline
    for spectrum = 1:spectrum_number
        fm(M, spectrum) = GMM_para_new(1, spectrum) * fg(xline(M), GMM_para_new(2, spectrum), GMM_para_new(3, spectrum)); % Gaussian distribution (GD) subroutine.
    end
end
%% Calculate output spectra.
GMM_subspectrum = fm*sum(yline(:))*dx; % Sub-spectra obtained by GMM.
GMM_reconstructed = sum(fm, 2)*sum(yline(:))*dx; % Reconstructed spectrum obtained by GMM.
%% Calculate D value by Eq. (43).
D_value = sum(abs(GMM_reconstructed - yline)) / sum(yline);
fprintf('GMM Fitting | D value = %.4f\n', D_value);
%% Plotting code.
figure();
plot(xline, yline, '-k', xline, GMM_reconstructed, '-r');
hold on;
for i = 1:spectrum_number
    plot(xline, GMM_subspectrum(:,i), '--b');
end
legend('Input Spectrum','GMM Spectrum','Sub-spectra','Location','best');
xlabel('Sample space', 'FontSize', 12, 'FontName', 'Arial');
ylabel('Proportion', 'FontSize', 12, 'FontName', 'Arial');
set(gca, 'FontSize', 12, 'FontName', 'Arial');
box on;
hold off;
end
% This subroutine generates the Gaussian distribution (GD).
% Input:
    % x - Sample space.
    % miu - Expectation.
    % theta - Standard deviation.
% Outputs:
    % fg - Gaussian distribution.
function fg = fg(x, miu, theta)
fg = 1 / sqrt(2*pi*theta) * exp(-(x - miu).^2 / (2*theta));
end