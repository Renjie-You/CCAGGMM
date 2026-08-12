% This subroutine implements the CCAGGMM spectral decomposition method.
% Inputs:
    % CCAGGMM_para - Initial values of the sub spectrum parameters for CCAGGMM.
    % confidence - Confidence level.
    % NR_circle - Maximum number of iterations for the Newton-Raphson method.
    % NR_tol - Tolerance of the Newton-Raphson method.
    % light_switch - Numerical solver selection parameter: =1 Newton-Raphson method; =2 Brute-force method.
% Outputs:
    % CCAGGMM_para_new - Sub-spectrum parameters obtained by the CCAGGMM algorithm.
% The global variables are defined in the 'spectrum_decomposition_main.m' function.
function [CCAGGMM_para_new] = CCAGGMM_decomposition(CCAGGMM_para, confidence, NR_circle, NR_tol, light_switch)
global dx length_yline hi spectrum_number CCAGGMM_circle xline
tic; % Start the clock.
CCAGGMM_para_new = CCAGGMM_para; % Sub-spectrum parameters obtained by the CCAGGMM algorithm.
x_lower = zeros(1, spectrum_number); % Lower quantiles of confidence intervals for sub-spectra.
x_upper = zeros(1, spectrum_number); % Upper quantiles of confidence intervals for sub-spectra.
%% CCAGGMM spectral decomposition
for circle = 1:CCAGGMM_circle
    %% Calculate the conditional probability by Eqs. (1) and (16).
    fprintf('\nIteration = %d\n', circle);
    fm = zeros(length_yline, spectrum_number); % Asymmetric generalized Gaussian distribution (AGGD).
    % Calculate the AGGD by Eq. (1).
    for i = 1:spectrum_number
        u = CCAGGMM_para(2,i); % Expectation.
        al = CCAGGMM_para(3,i); % Left branch scale factor.
        ar = CCAGGMM_para(4,i); % Right branch scale factor.
        beta = CCAGGMM_para(5,i); % Shape factor.
        gama = gamma(1/CCAGGMM_para(5,i)); % Gamma function.
        for j = 1:length_yline
            x = xline(j);
            if x < u - 1e-10 % Calculate the left branch of AGGD by Eq. (1).
                fm(j,i) = (beta / ((al + ar) * gama)) * exp(-((-x + u) / al)^beta);
            elseif x > u + 1e-10 % Calculate the right branch of AGGD by Eq. (1).
                fm(j,i) = (beta / ((al + ar) * gama)) * exp(-((x - u) / ar)^beta);
            else % Calculate the peak of AGGD.
                fm(j,i) = beta / ((al + ar) * gama);
            end
        end
    end
    % Calculate the conditional probability under the current parameters by Eq. (16).
    F_sub_spec = zeros(length_yline, spectrum_number); % Asymmetric generalized Gaussian mixture model (AGGMM).
    zim = zeros(length_yline, spectrum_number); % Conditional probability.
    for i = 1:spectrum_number
        for j = 1:length_yline
            F_sub_spec(j,i) = CCAGGMM_para(1,i) * fm(j,i); % Calculate the function of AGGMM by Eq. (3).
        end
    end
    sum_F = sum(F_sub_spec, 2); % Calculate the denominator of conditional probability in Eq. (16).
    for i = 1: length_yline
        % Avoid null values.
        if sum_F(i,1) < 10^(-10)
            zim(i,:) = 0;
        end
        if sum_F(i,1) > 10^(-10)
            zim(i,:) = F_sub_spec(i,:) / sum_F(i);
        end
    end
    %% Update paramenters by Eqs. (23)-(27).
    % The first parameter iteration uses asynchronous update method,
    % and the second and subsequent parameter iterations use synchronous update method.
    % The reason is explained at the end of Section 3.1 in the manuscript.
    if circle == 1
        for i = 1:spectrum_number
            con_hi_pre = hi; % Confidence-constrainted hi (Eq. (4)), and the confidence level is 1 now.
            pretend_zgm = zim(:, i) .* con_hi_pre(:,1); % Intermediate variables.
            % Newton-Raphson method.
            if light_switch == 1
                CCAGGMM_para_new(1,i) = sum(pretend_zgm(:,1)*dx); % Update weight by Eq. (23).
                CCAGGMM_para_new(2,i) = Update_miu_Newton_Raphson(CCAGGMM_para(2,i), CCAGGMM_para(3,i), CCAGGMM_para(4,i), CCAGGMM_para(5,i), con_hi_pre, NR_circle, NR_tol, zim(:,i)); % Update expectation by Eqs. (28)-(30).
                CCAGGMM_para_new(3,i) = Update_alphaL(CCAGGMM_para_new(2,i), CCAGGMM_para(5,i), con_hi_pre, zim(:,i)); % Update left branch scale factor by Eq. (25).
                CCAGGMM_para_new(4,i) = Update_alphaR(CCAGGMM_para_new(2,i), CCAGGMM_para(5,i), con_hi_pre, zim(:,i)); % Update right branch scale factor by Eq. (26).
                CCAGGMM_para_new(5,i) = Update_beta_Newton_Raphson(CCAGGMM_para_new(2,i), CCAGGMM_para_new(3,i), CCAGGMM_para_new(4,i), CCAGGMM_para(5,i), con_hi_pre, NR_circle, NR_tol, zim(:,i)); % Update shape factor by Eqs. (31)-(33).
            end
            % Brute-force method.
            if light_switch == 2
                CCAGGMM_para_new(1,i) = sum(pretend_zgm(:,1)*dx); % Update weight by Eq. (23).
                CCAGGMM_para_new(2,i) = Update_miu_brute_force(CCAGGMM_para(3,i), CCAGGMM_para(4,i), CCAGGMM_para(5,i), con_hi_pre, zim(:,i)); % Update expectation by Eq. (24).
                CCAGGMM_para_new(3,i) = Update_alphaL(CCAGGMM_para_new(2,i), CCAGGMM_para(5,i), con_hi_pre, zim(:,i)); % Update left branch scale factor by Eq. (25).
                CCAGGMM_para_new(4,i) = Update_alphaR(CCAGGMM_para_new(2,i), CCAGGMM_para(5,i), con_hi_pre, zim(:,i)); % Update right branch scale factor by Eq. (26).
                CCAGGMM_para_new(5,i) = Update_beta_brute_force(CCAGGMM_para_new(2,i), CCAGGMM_para_new(3,i), CCAGGMM_para_new(4,i), con_hi_pre, zim(:,i)); % Update shape factor by Eq. (27).
            end
            % Update upper and lower quantiles by Eqs. (41) and (42).
            [x_lower(1, i), x_upper(1, i)] = Confidence_interval(CCAGGMM_para_new(2,i), CCAGGMM_para_new(3,i), CCAGGMM_para_new(4,i), CCAGGMM_para_new(5,i), confidence);
        end
        % Constrain hi using the obtained upper and lower quantiles.
        mask = false(size(xline));
        for j = 1:numel(x_lower)
            mask = mask | (xline >= x_lower(j) & xline <= x_upper(j));
        end
        con_hi = zeros(size(hi));
        con_hi(mask) = hi(mask);
        CCAGGMM_para = CCAGGMM_para_new;
        continue;
    end
    for i = 1:spectrum_number
        con_hi_pre = con_hi; % Confidence-constrainted hi (Eq. (4)). The confidence level is defined by variable "confidence" now.
        pretend_zgm = zim(:, i) .* con_hi_pre(:,1); % Intermediate variables.
        % Newton-Raphson method.
        if light_switch == 1
            CCAGGMM_para_new(1,i) = sum(pretend_zgm(:,1)*dx); % Update weight by Eq. (23).
            CCAGGMM_para_new(2,i) = Update_miu_Newton_Raphson(CCAGGMM_para(2,i), CCAGGMM_para(3,i), CCAGGMM_para(4,i), CCAGGMM_para(5,i), con_hi_pre, NR_circle, NR_tol, zim(:,i)); % Update expectation by Eqs. (28)-(30).
            CCAGGMM_para_new(3,i) = Update_alphaL(CCAGGMM_para(2,i), CCAGGMM_para(5,i), con_hi_pre, zim(:,i)); % Update left branch scale factor by Eq. (25).
            CCAGGMM_para_new(4,i) = Update_alphaR(CCAGGMM_para(2,i), CCAGGMM_para(5,i), con_hi_pre, zim(:,i)); % Update right branch scale factor by Eq. (26).
            CCAGGMM_para_new(5,i) = Update_beta_Newton_Raphson(CCAGGMM_para(2,i), CCAGGMM_para(3,i), CCAGGMM_para(4,i), CCAGGMM_para(5,i), con_hi_pre, NR_circle, NR_tol, zim(:,i)); % Update shape factor by Eqs. (31)-(33).
        end
        % Brute-force method.
        if light_switch == 2
            CCAGGMM_para_new(1,i) = sum(pretend_zgm(:,1)*dx); % Update weight by Eq. (23).
            CCAGGMM_para_new(2,i) = Update_miu_brute_force(CCAGGMM_para(3,i), CCAGGMM_para(4,i), CCAGGMM_para(5,i), con_hi, zim(:,i)); % Update expectation by Eq. (24).
            CCAGGMM_para_new(3,i) = Update_alphaL(CCAGGMM_para(2,i), CCAGGMM_para(5,i), con_hi_pre, zim(:,i)); % Update left branch scale factor by Eq. (25).
            CCAGGMM_para_new(4,i) = Update_alphaR(CCAGGMM_para(2,i), CCAGGMM_para(5,i), con_hi_pre, zim(:,i)); % Update right branch scale factor by Eq. (26).
            CCAGGMM_para_new(5,i) = Update_beta_brute_force(CCAGGMM_para(2,i), CCAGGMM_para(3,i), CCAGGMM_para(4,i), con_hi_pre, zim(:,i)); % Update shape factor by Eq. (27).
        end
        % Update upper and lower quantiles by Eqs. (41) and (42).
        [x_lower(1, i), x_upper(1, i)] = Confidence_interval(CCAGGMM_para_new(2,i), CCAGGMM_para_new(3,i), CCAGGMM_para_new(4,i), CCAGGMM_para_new(5,i), confidence);
    end
    % Constrain hi using the obtained upper and lower quantiles.
    mask = false(size(xline));
    for j = 1:numel(x_lower)
        mask = mask | (xline >= x_lower(j) & xline <= x_upper(j));
    end
    con_hi = zeros(size(hi));
    con_hi(mask) = hi(mask);
    CCAGGMM_para = CCAGGMM_para_new;
end
runtime = toc; % Stop the clock.
fprintf('\nTotal runtime of AGGMM algorithm: %.4f seconds\n', runtime);
end