function [probr, info] = PreProcessing_withPOS(prob, option)

% Convert from Mosek to Sieve format
convert_mosek2sieve;

time_preprocessing = tic;

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

undeleted       = ones(m, 1);    % Keep track of which constraints are deleted
undone          = 1;   % undone = 1 means that it still needs (re-)preprocessing
info.infeasible = 0;   % 1 means that we have found infeasibility, 0 o/w
info.reduction  = 0;   % 1 means that we have found reduction, 0 o/w
constr_indices  = (1:m);
constr_num      = m;
iter            = 0;
cholEPS         = option.cholEPS;
bn              = -option.sqrtEPS*max(1, norm(b, inf)); % b < 0 if b < -sqrt(eps)*max{1, ||b||}
bz              = bn*option.sqrtEPS;    % b = 0 if -eps*max{1, ||b||} < b <= 0

% Begin preprocessing:
while undone
    
    undone = 0;
    iter = iter + 1;

    for ii = 1:constr_num
        i = constr_indices(ii);
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
                    info.infeasible = 1;
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
                        info.time_preprocessing = toc(time_preprocessing);
                        return;
                    end
                    if b(i) > bz
                        undeleted(i) = 1;
                    end
                end
            end
        end
        constr_indices = find(undeleted);
        constr_num     = length(constr_indices);
        undone         = 0;
    end
    
end

% do reduction
Ipos_nonzero   = Ipos(:, m + 1);
I_nonzero      = I(:, m + 1);
info.nonzero   = [Ipos_nonzero; I_nonzero];
info.undeleted = sparse(logical(undeleted));
undeleted      = constr_indices;
info.n_post    = nnz(Ipos_nonzero) + nnz(I_nonzero);
info.m_post    = constr_num;
if (info.n_post < n) || (info.m_post < m)
    info.reduction = 1;
    A_pos = A_pos(Ipos_nonzero, undeleted);
    c_pos = c_pos(Ipos_nonzero);
    for ii = 1:constr_num
        i            = undeleted(ii);
        A_convert{i} = A_convert{i}(I_nonzero, I_nonzero);
    end
    c_convert = c_convert(I_nonzero, I_nonzero);
else
    probr                   = prob;
    info.time_preprocessing = toc(time_preprocessing);
    return;
end

info.time_preprocessing = toc(time_preprocessing);

% Convert from Sieve to Mosek format
convert_sieve2mosek;

end
