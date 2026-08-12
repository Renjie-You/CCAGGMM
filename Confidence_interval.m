% This subfunction is used to calculate confidence intervals by Eqs. (41) and (42).
% Inputs:
    % u - Expectation.
    % al - Left branch scale factor.
    % ar - Right branch scale factor.
    % bm - Shape factor.
    % confidence - Confidence level.
% Outputs:
    % x_lower - Lower quantile of confidence interval.
    % x_upper - Upper quantile of confidence interval.
function [x_lower, x_upper] = Confidence_interval(u, al, ar, bm, confidence)
% Compute the lower quantile of confidence by Eq. (41).
x_lower = u - al * gammaincinv(confidence, 1/bm)^(1/bm);
% Compute the upper quantile of confidence by Eq. (42)
x_upper = u + ar * gammaincinv(confidence, 1/bm)^(1/bm);
end