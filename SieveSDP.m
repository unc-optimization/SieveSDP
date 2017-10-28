function [probr, info] = SieveSDP(prob, option)
% Purpose: SieveSDP is a preprocessing routine for semidefinite programming
% ...
% Inputs: 
%      prob: This variable takes Mosek format. It is a structure and 
%            has 9 fields: bardim, c, barc, a, bara, blc, buc, blx, bux,
%            where blc = buc; blx has entries -inf or 0; bux = [].
%            For detail, see e.g. http://docs.mosek.com/8.1/toolbox/tutorial-sdo-shared.html
%   options: It is a structure containing three fiels: 
%            maxiter, epsilon, cholEPS.
%            The default values are: 
%                   option.maxiter = intmax;
%                   option.epsilon = eps;
%                   option.cholEPS = 0;
% Outputs:
%     probr: The reduced problem structure in Mosek format.
%      info: This is a structure containing the information about the
%            preproessing.
%     > info has fields:
%                   info.n_pre
%                   info.m_pre
%                   info.n_post (if not infeasible)
%                   info.m_post (if not infeasible)
%                   info.infeasible (= 0 or 1)
%                   info.reduction  (= 0 or 1)
%                   info.time_preprocessing
%                   info.time_total
%
% Information:
%    Created by Yuzixuan (Melody) Zhu, Department of Statistics and
%    Operations Research, UNC-Chapel Hill.
%    Joint work with Gabor Pataki and Quoc Tran Dinh, UNC-Chapel Hill.
%    Created date: August 31, 2017.
%    Last modified: September 5, 2017.
%    Contact: zyzx@live.unc.edu
%    More information: http://...

if nargin < 2
   option.maxiter = intmax;
   option.epsilon = eps;
   option.cholEPS  = 0;
end

time_total = tic;

if (length(fieldnames(prob)) >= 10) || ~isequal(prob.blc, prob.buc)
    fprintf('Please clean up the data!\n');
    info.infeasible = 0;
    info.reduction = 0;
    info.time_preprocessing = 0;
    info.time_total = toc(time_total);
end

n_fre = nnz(prob.blx == -inf);
n_pos = nnz(prob.blx == 0);
if (length(prob.blx) - n_fre - n_pos > 0) || (~isempty(prob.bux) && any(~isinf(prob.bux)))
    fprintf('Please clean up the data!\n');
    info.infeasible = 0;
    info.reduction = 0;
    info.time_preprocessing = 0;
    info.time_total = toc(time_total);
end

% Choose the best PreProcessing function
%addpath(genpath(pwd));
option.sqrtEPS = sqrt(option.epsilon);
if n_pos == 0
    if n_fre == 0
        [probr, info] = PreProcessing_onlySDP(prob, option);
    else
        [probr, info] = PreProcessing_withFRE(prob, option);
    end
else
    if n_fre > 0
        I_pos = find(prob.blx == 0);
        I_fre = find(prob.blx == -inf);
        prob.c = [prob.c(I_pos); prob.c(I_fre); -prob.c(I_fre)];
        prob.a = [prob.a(:, I_pos), prob.a(:, I_fre), -prob.a(:, I_fre)];
        prob.blx = zeros(1, n_pos + 2*n_fre);
    end
    [probr, info] = PreProcessing_withPOS(prob, option);
end

info.time_total = toc(time_total);

end

% SeiveSDP v.1.0 by Melody Zhu, Gabor Pataki and Quoc Tran-Dinh.
% Copyright 2017 Department of Statistics and Operations Research
%                UNC - Chapel Hill, USA.
% See the file LICENSE for full license information.
