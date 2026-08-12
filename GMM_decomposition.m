% This subroutine implements the Gaussian mixtrue model (GMM) spectral decomposition method.
% Inputs:
    % yline - Proportion of input spectrum.
    % xline - Sample space of input spectrum.
    % spectrum_number - Number of sub-spectra.
    % GMM_para - Initial values of the sub spectrum parameters for GMM.
% Outputs:
    % GMM_para_new - Sub-spectrum parameters obtained by the GMM algorithm.
function GMM_para_new = GMM_decomposition(yline, xline, spectrum_number, GMM_para)
tic; % Start the clock.
% Spactrum sampling based on probability of sample space.
ns = 5000; % Number of drawn samples.
amp = yline(:,1); % Proportion of input spectrum.
amp(amp < 0) = 0; % Negative values are set to zero.
amp = amp / sum(amp); % Normalization.
ndata_x = randsample(xline, ns, true, amp); % Generate samples space.
spnmr = ns; % Size of sample space.
%% EM algorithm basing on GMM.
gama = zeros(spnmr, spectrum_number); % Sub-spctra matrix.
for idi = 1:100
    % E-step: Calculate posterior probability.
    for M = 1:spnmr
        x_real = ndata_x(M); % Samples.
        for N = 1:spectrum_number
            p(M,N) = fg(x_real, GMM_para(2,N), GMM_para(3,N)); % Gaussian function.
            gama(M,N) = p(M,N) * GMM_para(1,N); % GMM
        end
        lii(M) = p(M,:) * GMM_para(1,:)';
        gama(M,:) = gama(M,:) / lii(M); % Upddate posterior probability.
        gama(M, isnan(gama(M,:))) = 0;
    end
    for N = 1:spectrum_number
        nn(N) = sum(gama(:,N));
    end
    omega = nn / spnmr; % Updated sub-spectrum weight.
    % M-step: Update spectra parameters.
    miu = gama' * ndata_x ./ nn'; % Updated expectation.
    minus = zeros(spnmr, spectrum_number);
    for N = 1:spectrum_number
        minus(:,N) = gama(:,N) .* (ndata_x - miu(N)).^2;
    end
    epsilon = sum(minus) ./ nn; % Updated standard deviation.
    GMM_para = [omega; miu'; epsilon]; % Updated parameter matrix.
    GMM_para_new = GMM_para;
end
times = toc; % Stop the clock.
fprintf('GMM computation time: %.4f s\n', times);
end
% Gaussian probability density function.
% Input:
    % x - Samples space.
    % miu - Expectation.
    % theta - Standard deviation.
% Outputs:
    % fg - Gaussian probability density function value.
function fg = fg(x, miu, theta)
fg = 1 / sqrt(2*pi*theta) * exp(-(x - miu).^2 / (2*theta));
end