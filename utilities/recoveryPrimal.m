function [x_original, X_original] = recoveryPrimal(res, info)
% It recovers the original primal solution from reduced primal solution.
% It always gives the correct recovery
%
% Inputs:
%   res:  Given by calling mosekopt
%           >> [~, res] = mosekopt('minimize', prob);
%   info: Output 'info' from SieveSDP:
%           >> [probr, info] = SieveSDP(prob);
%
% Output:
%   x_original: Original primal solution corresponding to linear variables
%   X_original: Original primal solution correpsonding to PSD variables

if info.infeasible == 1
    fprintf('Infeasible\n');
    x_original = [];
    X_original = [];
    return;
end

% Linear variables
x_reduced  = res.sol.itr.xx;
x_original = zeros(info.n_pre.l, 1);
if isfield(info, 'nonzero_pos')
    x_original(info.nonzero_pos) = x_reduced((info.n_post.f + 1):end);
end
x_original = [x_reduced(1:info.n_post.f), x_original];

% PSD variables
barx      = res.sol.itr.barx;
bardim    = info.n_post.s;
n_sdp     = length(bardim);
n_sdp_sum = sum(bardim);
subk      = cell(n_sdp, 1);
subl      = cell(n_sdp, 1);
val       = cell(n_sdp, 1);

Aind  = zeros(n_sdp + 1, 1);
Ainds = Aind;
% convert Mosek solution to n*n matrix
for j = 1:n_sdp
    Aind(j + 1)              = Aind(j) + bardim(j);
    Ainds(j + 1)             = Ainds(j) + bardim(j)*(bardim(j) + 1)/2;
    X                        = sparse(bardim(j), bardim(j));
    X(tril(true(bardim(j)))) = barx((Ainds(j) + 1):Ainds(j + 1));
    [row, col, vv]           = find(X);
    I_diag                   = row == col;
    vv(I_diag)               = vv(I_diag)/2;
    row                      = Aind(j) + row;
    col                      = Aind(j) + col;
    subk{j}                  = [row; col];
    subl{j}                  = [col; row];
    val{j}                   = [vv; vv];
end

row = vertcat(subk{1:n_sdp});
col = vertcat(subl{1:n_sdp});
val = vertcat(val{1:n_sdp});
if isempty(val)
    X_reduced = sparse(n_sdp_sum, n_sdp_sum);
else
    X_reduced = sparse(row, col, val, n_sdp_sum, n_sdp_sum);
end

n_sdp_sum                              = sum(info.n_pre.s);
X_original                             = sparse(n_sdp_sum, n_sdp_sum);
X_original(info.nonzero, info.nonzero) = X_reduced;

end
