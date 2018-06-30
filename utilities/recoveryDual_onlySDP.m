function [y_original, z_original, Z_original, info1] = recoveryDual_onlySDP(res, info)

y_reduced  = res.sol.itr.y;
z_original = res.sol.itr.slx;
alphaM     = max([100; 2*abs(y_reduced)]);   % largest step length, can be changed
epsilon    = 1e-4;   % determine Z >= 0 if Z + epsilon*I > 0

% initialization
y_original                 = zeros(info.m_pre, 1);
y_original(info.undeleted) = y_reduced;

% compute PSD slack
I_undeleted = find(info.undeleted);
len         = length(I_undeleted);
subk        = cell(len, 1);
subl        = cell(len, 1);
val         = cell(len, 1);
for i = 1:len
    [subk{i}, subl{i}, val{i}] = find(info.DR.A_convert{I_undeleted(i)});
    val{i}                     = y_reduced(i)*val{i};
end
subk       = vertcat(subk{1:len});
subl       = vertcat(subl{1:len});
val        = vertcat(val{1:len});
n_sdp_sum  = sum(info.n_pre.s);
Z_original = sparse(subk, subl, val, n_sdp_sum, n_sdp_sum);
Z_original = info.DR.c_convert - Z_original;

% check if the slack is in the largest dual cone
I3             = info.nonzero;
[~, psd_check] = chol(Z_original(I3, I3) + sparse(epsilon*eye(nnz(I3))));
if psd_check ~= 0
    info1.success = 1;    % dual solution is not even in the largest dual cone
    info1.msg     = 'Slack not in the largest dual cone';
    info1.iter    = 0;
    return;
end

% dual recovery
N = length(info.DR.constr);
for j = N:-1:1
    i = info.DR.constr(j);
    
    % see the paper for the following block division
    A     = info.DR.A_convert{i}*info.DR.pd(j);
    I2    = info.DR.indices(:, j);
    I23   = logical(I2 + I3);
    A22   = A(I2, I2);
    Z22   = Z_original(I2, I2);
    A2233 = A(I23, I23);
    Z2233 = Z_original(I23, I23);
    
    % look for alpha s.t. Z31 + alpha*A31 = 0
    I1       = ~I23;
    A31      = A(I3, I1);
    A31nzind = find(A31);
    if ~isempty(A31nzind)
        Z31    = Z_original(I3, I1);
        A31nz  = A31(A31nzind);
        Z31nz  = Z31(A31nzind);
        alphas = -Z31nz./A31nz;
        alphas = unique(alphas(alphas >= 0));
        if ~isempty(alphas)
            if length(alphas) == 1
                [~, psd_check] = chol(Z2233 + alphas*A2233 + sparse(epsilon*eye(nnz(I23))));
                if psd_check == 0
                    y_original(i) = -alphas*info.DR.pd(j);
                    Z_original    = Z_original + alphas*A;
                    continue;
                end
                alpha1 = alphas;
            else
                found = false;
                for k = 1:length(alphas)
                    [~, psd_check] = chol(Z2233 + alphas(k)*A2233 + sparse(epsilon*eye(nnz(I23))));
                    if psd_check == 0
                        if k == length(alphas)
                            found = true;
                            break;
                        end
                        [~, pd_check] = chol(Z22 + alphas(k)*A22);
                        if pd_check == 0
                            found = true;
                            break;
                        end
                    end
                end
                if found
                    y_original(i) = -alphas(k)*info.DR.pd(j);
                    Z_original    = Z_original + alphas(k)*A;
                    continue;
                end
                alpha1 = alphas(k);
            end
        end
    else
        [~, psd_check] = chol(Z2233 + sparse(epsilon*eye(nnz(I23))));
        if psd_check == 0
            [~, pd_check] = chol(Z22);
            if pd_check == 0
                continue;
            end
            y_original(i) = -info.DR.pd(j);
            Z_original    = Z_original + A;
            continue;
        end
        alpha1 = 0;
    end

    % find the smallest step size
    matrix_psd = Z2233 + alpha1*A2233 + sparse(epsilon*eye(nnz(I23)));
    matrix_pd  = Z22 + alpha1*A22;
    for alpha = 1:(alphaM - alpha1)
        matrix_psd = matrix_psd + A2233;
        matrix_pd  = matrix_pd + A22;
        [~, psd_check] = chol(matrix_psd);
        if psd_check == 0
            [~, pd_check] = chol(matrix_pd);
            if pd_check == 0
                y_original(i) = -(alpha1 + alpha)*info.DR.pd(j);
                Z_original    = Z_original + (alpha1 + alpha)*A;
                break;
            end
            y_original(i) = -(alpha1 + alpha + 1)*info.DR.pd(j);
            Z_original    = Z_original + (alpha1 + alpha + 1)*A;
            break;
        end
        
        % if alpha = 1 and 2 do not work, check if alphaM works
        if alpha == 2
            [~, psd_check] = chol(Z2233 + alphaM*A2233 + sparse(epsilon*eye(nnz(I23))));
            if psd_check ~= 0
                info1.success = 0;    % cannot find a dual feasible point
                info1.iter    = N - j + 1;
                info1.msg     = ['Dual recovery fails at iteration ', num2str(info1.iter), '/', num2str(N)];
                return;
            end
        end
    end
    I3 = I23;

end

% linear variables
z_original = info.DR.c_fre - info.DR.A_fre*y_original;

% success
info1.success = 1;
info1.msg     = 'Dual recovery succeeds';
info1.iter    = N;

end
