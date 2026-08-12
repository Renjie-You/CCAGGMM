% Copyright Renjie You, Wei Li, Xiangjun Liu, Jian Xiong, Yi Ding,  2026/08/05.
% Please report any bug to yourenjie@163.com or liwei2014@email.cugb.edu.cn
% Referred manuscript title: "Confidence-Constrained Asymmetric Generalized Gaussian Mixture Model for Spectral Decomposition in Well Logging".
% The CCAGGMM algorithm is short for the "Confidence-Constrained Asymmetric Generalized Gaussian Mixture Model" algorithm.
% The main function of the CCAGGMM algorithm is complex spectral decomposition. The program has 2 advantages. 
% (1) it employs an asymmetric generalized Gaussian distribution as the sub-spectral characterization function,
% which can more accurately describe the symmetry and kurtosis characteristics of the sub-spectra.
% (2) it introduces confidence constraints, which enhance the noise robustness of the spectral decomposition.
clc; clear; close all;
global dx yline length_yline key spectrum_number CCAGGMM_circle xline hi
%% Input spectrum.
data = importdata('C:\Example 3 (Figure.6(a))-input spectrum.txt'); % Load raw data from TXT file.
%% Initialize parameters
spectrum_number = 2; % Number of sub-spectras.
% Method selection parameter
key = 1; % Method selection parameter for spectral decomposition: =1 for CCAGGMM; % =2 for GMM; % =3 for executing both CCAGGMM and GMM separately.
light_switch = 1; % Numerical solver selection parameter: =1 Newton-Raphson method; % =2 Brute-force method.
% inputs for CCAGGMM
if key == 1 || key == 3
    CCAGGMM_para = zeros(5, spectrum_number); % Initial values of the sub spectrum parameters for CCAGGMM.
    confidence = 0.999; % Confidence level: = 1-pc.
    CCAGGMM_para(1,:) = [0.5 0.5]; % Weight.
    CCAGGMM_para(2,:) = [30 70]; % Expectation.
    CCAGGMM_para(3,:) = [5 10]; % Left branch scale factor.
    CCAGGMM_para(4,:) = [5 10]; % Right branch scale factor.
    CCAGGMM_para(5,:) = [2 5]; % Shape factor.
end
% inputs for GMM
if key == 2 || key == 3
    GMM_para = zeros(3, spectrum_number); % Initial values of the sub spectrum parameters for GMM.
    GMM_para(1,:) = [0.5 0.5]; % Weight.
    GMM_para(2,:) = [30 70]; % Expectation.
    GMM_para(3,:) = [5 10]; % Standard deviation.
end
% Iteration termination parameters
CCAGGMM_circle = 200; % Maximum number of iterations for the EM algorithm.
NR_circle = 100; % Maximum number of iterations for the Newton Raphson method.
NR_tol = 1e-10; % Tolerance of Newton Raphson method.
%% Input spectrum normalization.
xline = data(:,1); % Sample space of input spectrum.
yline = data(:,2); % Proportion of input spectrum.
dx = xline(2) - xline(1); % Sampling interval.
length_yline = length(yline(:,1)); % Length of proportion column vector.
% Ignore extremely small values of proportion column vector. 
for i = 1:length_yline
    if yline(i,1) < 1e-10
        yline(i,1) = 0;
    end
end
% Normalize the proportion column vector.
sum_yline = sum(yline(:,1));
hi = yline;
hi(:,1) = yline(:,1)/(sum_yline * dx); % Proportion in Eq. (4).
%% Spectral decomposition.
% skew = zeros(1,spectrum_number); % Skewness of the sub-spectra.
% kurt = zeros(1,spectrum_number); % Kurtosis of the sub-spectra.
D_value = zeros(1,2); % Relative area error D definded by Eq. (49).
% CCAGGMM spectral decomposition.
if key == 1 || key == 3
    CCAGGMM_para_new = CCAGGMM_decomposition(CCAGGMM_para, confidence, NR_circle, NR_tol, light_switch); % CCAGGMM decomposition subroutine.
    [D_value(1,1), CCAGGMM_hi_new, CCAGGMM_subspectrum, skew, kurt] = GP_CCAGGMM(CCAGGMM_para_new); % Plot the results.
end
% GMM spectral decomposition.
if key == 2 || key == 3
    GMM_para_new = GMM_decomposition(yline, xline, spectrum_number, GMM_para); % GMM decomposition subroutine.
    [D_value(1,2), GMM_hi_new, GMM_subspectrum] = GP_GMM(GMM_para_new); % Plot the results.
end