% This subroutine updates the left branch scale factor by Eq. (25).
% Inputs:
    % u - Expectation.
    % beta - Shape factor.
    % con_hi - Confidence-constrainted hi (Eq. (4)).
    % zim - Conditional probability.
% Outputs:
    % al_new - Updated left branch scale factor after iteration.
% The global variables are defined in the 'spectrum_decomposition_main.m' function.
function al_new = Update_alphaL(u, beta, con_hi, zim)
global length_yline xline
f_al_sum = 0; % Molecules of Eq. (25).
f_zh_sum = 0; % Denominator of Eq. (25).
% Calculate molecules and denominator of Eq. (25).
for idx = 1:length_yline
    x = xline(idx);
    if x < u - 1e-10
        f_zh_sum = f_zh_sum + zim(idx,1) * con_hi(idx,1); % Denominator.
        f_al_sum = f_al_sum + beta * zim(idx,1) * con_hi(idx,1) * (-x + u)^beta; % Molecules.
    end
end
al_new = (f_al_sum / f_zh_sum)^(1/beta); % Eq. (25).
end