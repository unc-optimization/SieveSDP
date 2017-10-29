function [probr, info] = PreProcessing_withFRE(prob, option)

% Convert from Mosek to Sieve format
convert_mosek2sieve;

time_preprocessing = tic;

% initial nonzero indices of each constraint matrix
I = true(n_sdp_sum, m + 1);
for i = 1:m
    I(:, i) = any(A_convert{i}, 2);
end
I = sparse(I);

constr_fre      = any(prob.a, 2);
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
    iter = iter + 1;

    for ii = 1:constr_num
        i  = constr_indices(ii);
        At = A_convert{i}(I(:, i), I(:, i));    % Get the nonzero submatrix
        
        Iaux = any(At, 2);
        if find(Iaux == false, 1)
            I(I(:, i), i) = Iaux;
            At            = At(Iaux, Iaux);
        end
        if isempty(At)  % if Ai = 0 and bi < 0, then infeasible
            if b(i) < bn
                info.infeasible         = 1;
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
            if pd_check == 0    % if Ai pd and bi < 0, then infeasible
                info.infeasible         = 1;
                info.time_preprocessing = toc(time_preprocessing);
                return;
            end
        else
            if b(i) > bz
                if cholEPS > 0
                    [~, pd_check] = chol(At - cholEPS*eye(size(At, 1)));
                else
                    [~, pd_check] = chol(At);
                end
                if pd_check == 0    % if Ai pd and bi = 0, then reduce
                    I(I(:, i), :) = false;
                    undeleted(i)  = 0;
                    undone        = 1;
                else
                    if cholEPS > 0
                        [~, nd_check] = chol(-(At - cholEPS*eye(size(At, 1))));
                    else
                        [~, nd_check] = chol(-At);
                    end
                    if nd_check == 0    % if Ai nd and bi = 0, then reuce
                        I(I(:, i), :) = false;
                        undeleted(i)  = 0;
                        undone        = 1;
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
            i  = constr_indices(ii);
            At = A_convert{i}(I(:, i), I(:, i));
            
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
        constr_indices = find(undeleted);
        constr_num     = length(constr_indices);
        undone         = 0;
    end
    
end

% do reduction
I_nonzero      = I(:, m + 1);
info.nonzero   = [ones(n_fre, 1); I_nonzero];
undeleted      = undeleted + constr_fre;
info.undeleted = sparse(logical(undeleted));
undeleted      = find(undeleted);
info.n_post    = n_fre + nnz(I_nonzero);
info.m_post    = length(undeleted);

if (info.n_post < n) || (info.m_post < m)
    info.reduction = 1;
    for ii = 1:info.m_post
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