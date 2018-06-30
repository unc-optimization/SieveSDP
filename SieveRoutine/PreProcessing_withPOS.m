function [probr, info] = PreProcessing_withPOS(prob, option)

% Convert from Mosek to Sieve format
convert_mosek2sieve;

time_preprocessing = tic;

% information for dual recovery
if option.DR == 1
    info.DR.A_convert   = A_convert;	% record A_convert
    info.DR.A_pos       = A_pos;
    info.DR.A_fre       = prob.a(:, 1:n_fre)';
    info.DR.c_convert   = c_convert;
    info.DR.c_pos       = c_pos;
    info.DR.c_fre       = prob.c(1:n_fre);
    info.DR.indices     = sparse(false(n_sdp_sum, m));  % record faces
    info.DR.indices_pos = sparse(false(n_pos, m));
    info.DR.pd          = ones(m, 1);   % 1 if deleted constraint is pd, -1 if nd
    info.DR.constr      = zeros(m, 1);  % the order when constraints are deleted
    j                   = 0;    % count iteration for facial reduction
end

% make b negative
neg    = find(b > 0);
b(neg) = -b(neg);
len    = length(neg);
for i = 1:len
    A_convert{neg(i)} = -A_convert{neg(i)};
end
A_pos(:, neg) = -A_pos(:, neg);

% initial nonzero indices of each constraint matrix
Ipos = true(n_pos, m + 1);
for i = 1:m
    Ipos(:, i) = any(A_pos(:, i), 2);
end
Ipos = sparse(Ipos);

I = true(n_sdp_sum, m + 1);
for i = 1:m
    I(:, i) = any(A_convert{i}, 2);
end
I = sparse(I);

constr_fre      = any(prob.a(:, 1:n_fre), 2);
undeleted       = ~constr_fre;    % Keep track of which constraints are deleted
undone          = 1;   % undone = 1 means that it still needs (re-)preprocessing
info.infeasible = 0;   % 1 means that we have found infeasibility, 0 o/w
info.reduction  = 0;   % 1 means that we have found reduction, 0 o/w
constr_indices  = find(undeleted);
constr_num      = length(constr_indices);
iter            = 0;
cholEPS         = option.cholEPS;
bn              = -option.sqrtEPS*max(1, norm(b, inf)); % b < 0 if b < -sqrt(eps)*max{1, ||b||}
bz              = bn*option.sqrtEPS;    % b = 0 if -eps*max{1, ||b||} < b <= 0

% Begin preprocessing:
while undone
    
    undone = 0;
    iter   = iter + 1;

    for ii = 1:constr_num
        i    = constr_indices(ii);
        Apos = A_pos(Ipos(:, i), i);    % get the nonzero constraint coefficients APOSi for pos vars
        
        if isempty(find(Apos, 1))
            At   = A_convert{i}(I(:, i), I(:, i));    % get the nonzero submatrix when APOSi = 0
            Iaux = any(At, 2);
            if find(Iaux == false, 1)
                I(I(:, i), i) = Iaux;
                At            = At(Iaux, Iaux);
            end
            
            if isempty(At)  % if APOSi = 0, Ai = 0, and bi < 0, then infeasible
                if b(i) < bn
                    info.infeasible         = 1;
                    info.iter               = iter;
                    info.time_preprocessing = toc(time_preprocessing);
                    return;
                end
                if b(i) > bz
                    undeleted(i) = 0;
                    continue;
                end
            end
            
            if b(i) < bn
                if cholEPS > 0
                    [~, pd_check] = chol(At - cholEPS*eye(size(At, 1)));
                else
                    [~, pd_check] = chol(At);
                end
                if pd_check == 0    % if APOSi = 0, Ai pd, and bi < 0, then infeasible
                    info.infeasible         = 1;
                    info.iter               = iter;
                    info.time_preprocessing = toc(time1);
                    return;
                end
            else
                if b(i) > bz
                    if cholEPS > 0
                        [~, pd_check] = chol(At - cholEPS*eye(size(At, 1)));
                    else
                        [~, pd_check] = chol(At);
                    end
                    if pd_check == 0    % if APOSi = 0, Ai pd, and bi = 0, then reduce
                        if option.DR == 1
                            j                     = j + 1;
                            info.DR.constr(j)     = i;
                            info.DR.indices(:, j) = I(:, i);
                        end
                        I(I(:, i), :) = false;
                        undeleted(i)  = 0;
                        undone        = 1;
                    else
                        if cholEPS > 0
                            [~, nd_check] = chol(-(At - cholEPS*eye(size(At, 1))));
                        else
                            [~, nd_check] = chol(-At);
                        end
                        if nd_check == 0    % if APOSi = 0, Ai nd, and bi = 0, then reduce
                            if option.DR == 1
                                j                     = j + 1;
                                info.DR.constr(j)     = i;
                                info.DR.indices(:, j) = I(:, i);
                                info.DR.pd(j)         = -1; % this constraint was nd, now it is pd
                            end
                            I(I(:, i), :) = false;
                            undeleted(i)  = 0;
                            undone        = 1;
                        end
                    end
                end
            end
        else
            if b(i) < bn
                if Apos > cholEPS
                    At = A_convert{i}(I(:, i), I(:, i));    % get nonzero submatrix when APOSi > 0, b < 0
                    Iaux = any(At, 2);
                    if find(Iaux == false, 1)
                        I(I(:, i), i) = Iaux;
                        At            = At(Iaux, Iaux);
                    end
                    if isempty(At)  % if APOSi > 0, Ai = 0, bi < 0, then infeasible
                        info.infeasible         = 1;
                        info.iter               = iter;
                        info.time_preprocessing = toc(time_preprocessing);
                        return;
                    else
                        if cholEPS > 0
                            [~, pd_check] = chol(At - epsilon*eye(size(At, 1)));
                        else
                            [~, pd_check] = chol(At);
                        end
                        if pd_check == 0    % if APOSi > 0, Ai pd, bi < 0, then reduce
                            info.infeasible         = 1;
                            info.iter               = iter;
                            info.time_preprocessing = toc(time_preprocessing);
                            return;
                        end
                    end
                end
            else
                if b(i) > bz
                    if Apos > cholEPS
                        At   = A_convert{i}(I(:, i), I(:, i));    % get nonzero submatrix when APOSi > 0, b = 0
                        Iaux = any(At, 2);
                        if find(Iaux == false, 1)
                            I(I(:, i), i) = Iaux;
                            At            = At(Iaux, Iaux);
                        end
                        if isempty(At)  % if APOSi > 0, Ai = 0, b = 0, reduce
                            if option.DR == 1
                                j                         = j + 1;
                                info.DR.constr(j)         = i;
                                info.DR.indices_pos(:, j) = Ipos(:, i);
                            end
                            Ipos(Ipos(:, i), :) = false;
                            undeleted(i)        = 0;
                            undone              = 1;
                        else
                            if cholEPS > 0
                                [~, pd_check] = chol(At - epsilon*eye(size(At, 1)));
                            else
                                [~, pd_check] = chol(At);
                            end
                            if pd_check == 0    % if APOSi > 0, Ai pd, b = 0, reduce
                                if option.DR == 1
                                    j                         = j + 1;
                                    info.DR.constr(j)         = i;
                                    info.DR.indices(:, j)     = I(:, i);
                                    info.DR.indices_pos(:, j) = Ipos(:, i);
                                end
                                Ipos(Ipos(:, i), :) = false;
                                I(I(:, i), :)       = false;
                                undeleted(i)        = 0;
                                undone              = 1;
                            end
                        end
                    else
                        if Apos < -cholEPS
                            At = A_convert{i}(I(:, i), I(:, i));    % get nonzero submatrix when APOSi < 0, b = 0
                            Iaux = any(At, 2);
                            if find(Iaux == false, 1)
                                I(I(:, i), i) = Iaux;
                                At            = At(Iaux, Iaux);
                            end
                            if isempty(At)  % if APOSi < 0, Ai = 0, b = 0, reduce
                                if option.DR == 1
                                    j                         = j + 1;
                                    info.DR.constr(j)         = i;
                                    info.DR.indices_pos(:, j) = Ipos(:, i);
                                    info.DR.pd(j)             = -1;
                                end
                                Ipos(Ipos(:, i), :) = false;
                                undeleted(i)        = 0;
                                undone              = 1;
                            else
                                if cholEPS > 0
                                    [~, nd_check] = chol(-(At - epsilon*eye(size(At, 1))));
                                else
                                    [~, nd_check] = chol(-At);
                                end
                                if nd_check == 0    % if APOSi < 0, Ai nd, b = 0, reduce
                                    if option.DR == 1
                                        j                         = j + 1;
                                        info.DR.constr(j)         = i;
                                        info.DR.indices(:, j)     = I(:, i);
                                        info.DR.indices_pos(:, j) = Ipos(:, i);
                                        info.DR.pd(j)             = -1;
                                    end
                                    Ipos(Ipos(:, i), :) = false;
                                    I(I(:, i), :)       = false;
                                    undeleted(i)        = 0;
                                    undone              = 1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    constr_indices = find(undeleted);
    constr_num     = length(constr_indices);
    
    % delete zero constraints after reaching maxiter
    if iter >= option.maxiter
        for ii = 1:constr_num
            i    = constr_indices(ii);
            Apos = A_pos(Ipos(:, i), i);
            if isempty(find(Apos, 1))
                At   = A_convert{i}(I(:, i), I(:, i));
                Iaux = any(At, 2);
                if find(Iaux == false, 1)
                    I(I(:, i), i) = Iaux;
                    At            = At(Iaux, Iaux);
                end
                if isempty(At)
                    if b(i) < bn
                        info.infeasible         = 1;
                        info.iter               = iter;
                        info.time_preprocessing = toc(time_preprocessing);
                        return;
                    end
                    if b(i) > bz
                        undeleted(i) = 0;
                    end
                end
            end
        end
        constr_indices = find(undeleted);
        constr_num     = length(constr_indices);
        undone         = 0;
    end
    
end
info.iter = iter;

% reverse signs
b(neg) = -b(neg);
% len = length(neg);
for i = 1:len
    A_convert{neg(i)} = -A_convert{neg(i)};
end
A_pos(:, neg) = -A_pos(:, neg);

% do reduction
Ipos_nonzero     = Ipos(:, m + 1);
I_nonzero        = I(:, m + 1);
info.nonzero     = I_nonzero;
info.nonzero_pos = Ipos_nonzero;
undeleted        = undeleted + constr_fre;
info.undeleted   = sparse(logical(undeleted));
undeleted        = find(undeleted);

info.n_post.f = n_fre;
info.n_post.l = nnz(Ipos_nonzero);
info.m_post   = length(undeleted);

if option.DR == 1
    info.DR.constr      = info.DR.constr(1:j);
    info.DR.indices     = info.DR.indices(:, 1:j);
    info.DR.indices_pos = info.DR.indices_pos(:, 1:j);
    info.DR.pd          = info.DR.pd(1:j);
end

% check reduction
if info.m_post < m
    info.reduction = 1;
    A_pos = A_pos(Ipos_nonzero, undeleted);
    c_pos = c_pos(Ipos_nonzero);
    for ii = 1:info.m_post
        i            = undeleted(ii);
        A_convert{i} = A_convert{i}(I_nonzero, I_nonzero);
    end
    c_convert = c_convert(I_nonzero, I_nonzero);
else
    probr                   = prob;
    info.n_post.s           = info.n_pre.s;
    info.time_preprocessing = toc(time_preprocessing);
    return;
end

info.time_preprocessing = toc(time_preprocessing);

% Convert from Sieve to Mosek format
convert_sieve2mosek;

end
