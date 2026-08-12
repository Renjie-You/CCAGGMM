% This subroutine updates the right branch scale factor by Eq. (26).
% Inputs:
    % u - Expectation.
    % beta - Shape parameter.
    % con_hi - Confidence-constrainted hi (Eq. (4)).
    % zim - Conditional probability.
% Outputs:
    % ar_new - Updated right branch scale factor after iteration.
% The global variables are defined in the 'spectrum_decomposition_main.m' function.
function ar_new = Update_alphaR(u, beta, con_hi, zim)
global length_yline xline
f_ar_sum = 0; % Molecules of Eq. (26).
f_zh_sum = 0; % Denominator of Eq. (26).
% Calculate molecules and denominator of Eq. (26).
for idx = 1:length_yline
    x = xline(idx);
    if x > u + 1e-10
        f_zh_sum = f_zh_sum + zim(idx,1) * con_hi(idx,1); % Denominator.
        f_ar_sum = f_ar_sum + beta * zim(idx,1) * con_hi(idx,1) * (x - u)^beta; % Molecules.
    end
end
ar_new = (f_ar_sum / f_zh_sum)^(1/beta); % Eq. (26).
end