function info = mosekCall(prob)

%%%%%%%%%%%
% Solving %
%%%%%%%%%%%

t = tic;
[~, info] = mosekopt('minimize info', prob);
if ~isfield(info, 'sol')
    info = []; fprintf('Mosek is failed to solve this problem.\n'); return;
end
info.time = toc(t);
if isfield(info.sol.itr, 'pobjval'), info.obj1 = info.sol.itr.pobjval; end
if isfield(info.sol.itr, 'dobjval'), info.obj2 = info.sol.itr.dobjval; end
if isnan(info.obj1) || strcmp(info.sol.itr.prosta, 'PRIMAL_INFEASIBLE')
    info.infeasible = 1;
else
    if isnan(info.obj2) || strcmp(info.sol.itr.prosta, 'DUAL_INFEASIBLE')
        info.infeasible = -1;
    else
        info.infeasible = 0;
    end
end

%%%%%%%%%%
% Dimacs %
%%%%%%%%%%

if isinf(info.obj1) || isinf(info.obj1) || isnan(info.obj1) || ...
        isinf(info.obj2) || isinf(info.obj2) || isnan(info.obj2)
    info.dimacsfull = [];
    info.DIMACS = NaN;
    return;
end

[A, b, c, K] = convert_mosek2sedumi(prob);
if size(A, 1) ~= length(b)
    A = A';
end
if size(b, 2) ~= 1
    b = b';
end
if size(c, 2) ~= 1
    c = c';
end
n_fre = 0;
if isfield(K, 'f') && isempty(K.f) && (K.f > 0)
    n_fre = K.f;
end
n_pos = 0;
if isfield(K, 'l') && isempty(K.l) && (K.l > 0)
    n_pos = K.l;
end

nj = prob.bardim;
njs = nj.^2;
n_sdp = length(nj);

xx = info.sol.itr.xx;
barx = info.sol.itr.barx;
if isempty(xx)
    x = [];
else
    x = xx;
end
s = info.sol.itr.slx - info.sol.itr.sux + info.sol.itr.snx;
bars = info.sol.itr.bars;
y = info.sol.itr.y;
xres = inf;
sres = inf;
if n_pos > 0
    xres = min(xx((1 + n_fre):(n_fre + n_pos)));
    sres = min(s((1 + n_fre):(n_fre + n_pos)));
end
count = 0;
for j = 1:n_sdp
   nn = nj(j)*(nj(j) + 1)/2;
   X = tril(ones(nj(j)));
   X(X == 1) = barx((count + 1):(count + nn));
   X = X + X' - diag(diag(X));
   x = [x; reshape(X, [njs(j), 1])];
   xres = min(xres, min(eig(X)));
   S = tril(ones(nj(j)));
   S(S == 1) = bars((count + 1):(count + nn));
   S = S + S' - diag(diag(S));
   s = [s; reshape(S, [njs(j), 1])];
   sres = min(sres, min(eig(S)));
   count = count + nn;
end

err1 = norm(b - A*x)/(1 + norm(b, inf));
err2 = max(0, -xres)/(1 + norm(b, inf));
err3 = norm(c - A'*y - s, 'fro')/(1 + norm(c, inf));
err4 = max(0, -sres)/(1 + max(abs(c)));
err5 = (info.obj1 - info.obj2)/(1 + abs(info.obj1) + abs(info.obj2));
err6 = x'*s/(1 + abs(info.obj1) + abs(info.obj2));
info.dimacsfull = [err1, err2, err3, err4, err5, err6];
info.DIMACS = max(abs(info.dimacsfull));

end
