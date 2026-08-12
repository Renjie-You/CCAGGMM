% This subroutine uses brute-force method to update the shape factor by Eq. (27).
% Inputs:
    % u - Expectation.
    % al - Left branch scale factor.
    % ar - Right branch scale factor.
    % con_hi - Confidence-constrainted hi (Eq. (4)).
    % zim - Conditional probability.
% Outputs:
    % beta_new - Updated shape factor after iteration.
% The global variables are defined in the 'spectrum_decomposition_main.m' function.
function beta_new = Update_beta_brute_force(u, al, ar, con_hi, zim)
global xline length_yline;
beta_range = 0.01:0.01:10; % Search range of shape factor.
length_bm = length(beta_range); % Length of the search area.
F_beta = zeros(1, length_bm);% Value of left side of Eq. (27).
for i = 1:length_bm
    bm = beta_range(i);
    F_part1 = sum(zim .* con_hi) * (1/bm + psi(1/bm)/bm^2); % First term in Eq. (27).
    F_part2 = 0; % Initialize second term of Eq. (27).
    F_part3 = 0; % Initialize third term of Eq. (27).
    for j = 1:length_yline
        x = xline(j);
        z_h = zim(j) * con_hi(j); % Intermediate variable.
        if x < u - 1e-10
            F_part2 = F_part2 + z_h * ((-x + u) / al)^bm * log((-x + u) / al); % Second term in Eq.(27).
        end
        if x > u + 1e-10
            F_part3 = F_part3 + z_h * ((x - u) / ar)^bm * log((x - u) / ar); % Third term in Eq.(27).
        end
    end
    % When x = u, the values of Eq. (27) is zero and can be ignored.
    F_beta(i) = F_part1 - F_part2 - F_part3; % Eq. (27).
end
%% Finding the optimal value of shape factor.
[~, best_idx] = min(abs(F_beta));
beta_new = beta_range(best_idx);
end
