% This subroutine uses brute-force method to update the expectation by Eq. (24).
% Inputs:
    % al - Left branch scale factor.
    % ar - Right branch scale factor.
    % beta - Shape parameter.
    % con_hi - Confidence-constrainted hi (Eq. (4)).
    % zim - Conditional probability.
% Outputs:
    % u_new - Updated expectation after iteration.
% The global variables are defined in the 'spectrum_decomposition_main.m' function.
function [u_new] = Update_miu_brute_force(al, ar, beta, con_hi, zim)
global length_yline xline
F = zeros(1,length_yline); % Value of left side of Eq. (24).
for i = 1:length_yline
    u = xline(i); % Search range of expectation.
    sum_f1_L = 0; % Sum of first term in Eq. (24).
    sum_f1_R = 0; % Sum of second term in Eq. (24).
    for idx = 1:length_yline
        x = xline(idx);
        z_h = zim(idx, 1) * con_hi(idx, 1); % Intermediate variable.
        % Calculate each term in Eq. (24).
        if x < u - 1e-10
            sum_f1_L = sum_f1_L + z_h * (-x + u)^(beta-1); % First term. 
        end
        if x > u + 1e-10
            sum_f1_R = sum_f1_R + z_h * (x - u)^(beta-1); % Second term.
        end
    end
    % When x = u, the values is zero and can be ignored.
    F(i) = sum_f1_L/(al^beta) - sum_f1_R/(ar^beta); % Eq. (24).
end
%% Finding the optimal value of expectation.
[~, best_idx] = min(abs(F));
u_new = xline(best_idx);
end
