% This subroutine uses Newton-Raphson method to update the shape factor by Eqs. (31)-(33).
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
    % beta_new - Updated shape factor after iteration.
% The global variables are defined in the 'spectrum_decomposition_main.m' function.
function beta_new = Update_beta_Newton_Raphson(u, al, ar, beta, con_hi, NT_circle, NT_tol, zim)
global length_yline xline;
beta_new = beta;
for k = 1:NT_circle
    beta = beta_new;
    % Calculate Eq. (31).
    F_part1 = sum(zim .* con_hi) * (1/beta + psi(1/beta)/beta^2); % First term of Eq. (31).
    F_part2 = 0; % Initialize second term of Eq. (31).
    F_part3 = 0; % Initialize third term of Eq. (31).
    for idx = 1:length_yline
        x = xline(idx);
        z_h = zim(idx) * con_hi(idx,1); % Intermediate variable.
        if x < u - 1e-10
            F_part2 = F_part2 + z_h * ((-x + u) / al)^beta * (log((-x + u) / al)); % Second term in Eq.(31).
        elseif x > u + 1e-10
            F_part3 = F_part3 + z_h * ((x - u) / ar)^beta * (log((x - u) / ar)); % Third term in Eq.(31).
        end
        % When x = u, the values of Eq. (31) is zero and can be ignored.
    end
    F = F_part1 - F_part2 - F_part3; % Eq. (31).
    % Calculate Eq. (32).
    dF_part1 = sum(zim .* con_hi) * (-1/beta^2 + (-psi(1, 1/beta) - 2*beta*psi(1/beta)) / beta^4); % First term of Eq. (32).
    dF_part2 = 0; % Initialize second term of Eq. (32).
    dF_part3 = 0; % Initialize third term of Eq. (32).
    for idx = 1:length_yline
        x = xline(idx);
        z_h = zim(idx) * con_hi(idx,1); % Intermediate variable.
        if x < u - 1e-10
            dF_part2 = dF_part2 + z_h * (((-x + u) / al)^beta) * (log((-x + u) / al))^2; % Second term in Eq.(32).
        elseif x > u + 1e-10
            dF_part3 = dF_part3 + z_h * (((x - u) / ar)^beta) * (log((x - u) / ar))^2; % Third term in Eq.(32).
        end
        % When x = u, the values of Eq. (32) is zero and can be ignored.
    end
    dF = dF_part1 - dF_part2 - dF_part3; % Eq. (32).
    %% Shape factor iterition function by Eq. (33).
    beta_new = beta - F / dF;
    % Termination condition for Newton-Raphson method.
    if abs(beta_new - beta) < NT_tol
        beta_new = beta;
        break;
    end
end
end