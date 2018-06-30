function [y_original, z_original, Z_original, info1] = recoveryDual(res, info)
% It recovers the original dual solution from reduced dual solution.
%
% Inputs:
%   res:   Given by calling mosekopt
%           >> [~, res] = mosekopt('minimize info', probr);
%   info0: Output from SieveSDP:
%           >> [probr, info0] = SieveSDP(prob);
%
% Outputs:
%   y_original, z_original, Z_original: recovered dual solution
%   info.success = 1, if recovery succeeds (finds an optimal dual solution)
%                = 0, if recovery fails (reduced dual solution may be dual infeasible)
%   info.msg
%   info.iter
%   info.time

time = tic;

% check if no need for recovery
if ~isfield(info, 'DR')
    y_original    = res.sol.itr.y;
    z_original    = res.sol.itr.slc;
    Z_original    = res.sol.itr.bars;
    info1.success = 0;
    info1.msg     = 'Please set option.DR = 1 when preprocessing';
    info1.iter    = 0;
    info1.time    = toc(time);
    return;
end

if info.infeasible == 1
    y_original    = [];
    z_original    = [];
    Z_original    = [];
    info1.success = 0;
    info1.msg     = 'Primal is infeasible. Please solve original problem for dual solution';
    info1.iter    = 0;
    info1.time    = toc(time);
    return;
end

if info.reduction == 0
    y_original    = res.sol.itr.y;
    z_original    = res.sol.itr.slc;
    Z_original    = res.sol.itr.bars;
    info1.success = 1;
    info1.msg     = 'No reduction by Sieve. No need for dual recovery';
    info1.iter    = 0;
    info1.time    = toc(time);
    return;
end

if info.n_pre.l > 0
    [y_original, z_original, Z_original, info1] = recoveryDual_withPOS(res, info);
else
    [y_original, z_original, Z_original, info1] = recoveryDual_onlySDP(res, info);
end

info1.time = toc(time);

end
