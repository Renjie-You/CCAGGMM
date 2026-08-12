% This subroutine uses Newton-Raphson method to update the expectation by Eqs. (28)-(30).
% Inputs:
    % u - Expectation.
    % al - Left branch scale factor.
    % ar - Right branch scale factor.
    % beta - Shape factor.
    % con_hi - Confidence-constrainted hi (Eq. (4)).
    % NR_circle - Maximum number of iterations for the Newton-Raphson method.
    % NR_tol - Tolerance of the Newton-Raphson method.
    % zim - Conditional probability.
% Outputs:
    % u_new - Updated expectation after iteration.
% The global variables are defined in the 'spectrum_decomposition_main.m' function.
function [u_new] = Update_miu_Newton_Raphson(u, al, ar, beta, con_hi, NR_circle, NR_tol, zim)
global length_yline xline;
u_new = u; % Expectation after iteration.
for k = 1:NR_circle
    u = u_new;
    sum_f1_L = 0; % Sum of the x < μ part of Eq. (28).
    sum_f1_R = 0; % Sum of the x >= μ part of Eq. (28).
    sum_f2_L = 0; % Sum of the x < μ part of Eq. (29).
    sum_f2_R = 0; % Sum of the x >= μ part of Eq. (29).
    for idx = 1:length_yline
        x = xline(idx);
        z_h = zim(idx, 1) * con_hi(idx,1); % Intermediate variable.
        % Calculate each term in Eq. (28).
        if x - u < -1e-10
            sum_f1_L = sum_f1_L + z_h * (-x + u)^(beta-1); % First term.
        elseif x - u > 1e-10
            sum_f1_R = sum_f1_R + z_h * (x - u)^(beta-1); % Second term.
        end
        % Calculate each term in Eq. (29).
        if x - u < -1e-10
            sum_f2_L = sum_f2_L + z_h * ((beta-1) * (-x + u)^(beta-2)); % First term.
        elseif x - u > 1e-10
            sum_f2_R = sum_f2_R + z_h * ((beta-1) * (x - u)^(beta-2)); % Second term.
        end
        % When x = u, the values of Eqs. (28) and (29) are zero and can be ignored.
    end
    %% Expectation iterition function by Eq. (30).
    f_1_sum = sum_f1_L / (al^beta) - sum_f1_R / (ar^beta); % Eq. (28).
    f_2_sum = sum_f2_L / (al^beta) + sum_f2_R / (ar^beta); % Eq. (29).
    u_new = u - f_1_sum / f_2_sum; % Eq. (30).
    % Termination condition for Newton-Raphson method.
    if abs(u_new - u) < NR_tol
        break;
    end
end
end