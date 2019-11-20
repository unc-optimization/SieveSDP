function [probr, info] = SieveSDP(prob, option)
%
% Purpose: SieveSDP is a preprocessing routine for semidefinite programming
%   min. <C, X>
%   s.t. <Ai, X> == bi (i = 1, ..., m)
%              X >= 0
% where X >= 0 means that X is in positive semidefinite (PSD) cone. X may also
% include linear (free and/or nonnegative) variables:
%   min. <cf, xf> + <cl, xl> + <C, X>
%   s.t. <afi, xf> + <ali, xl> + <C, X> = bi (i = 1, ..., m)
%              xl >= 0, X >= 0
%
% Inputs: 
%      prob: This variable takes Mosek format. It is a structure with fields:
%            bardim:    sizes of PSD blocks
%            c:         obj coeficient for linear variables
%            barc:      obj coeficient for PSD variables
%            a:         constraint coeficient for linear variables
%            bara:      constraint coeficient for PSD variables
%            blc = buc: constraint rhs
%            blx:       -inf for free variables, 0 for nonnegative variables
%            bux:       [] for free and nonnegative linear variables
%            For detail, see e.g. Section 9.7 of
%            http://docs.mosek.com/7.0/toolbox/A_guided_tour.html
%   option: It is a structure containing fields: 
%            maxiter, epsilon, cholEPS and DR.
%            The default values are: 
%            maxiter = intmax:  maximum iteration number
%            epsilon = eps:     accuracy to check if bi == 0
%            cholEPS = 0:       accuracy to check psd-ness
%            DR = 0:            save information for dual recovery?
%
% Outputs:
%     probr: The reduced problem structure in Mosek format.
%      info: This is a structure containing the information about the
%            preproessing. info has fields:
%            n_pre:       order of X, including fields n_pre.f, .l and .s
%            n_post:      order of X_reduced, if not infeasible
%            m_pre:       number of constraints
%            m_post:      number of constraints in reduced problem, if not infeasible
%            nonzero:     indices of nonzero rows/columns of PSD variable X IN BINARY
%            nonzero_pos: indices of nonzero nonnegative linear variable x
%            undeleted:   indices of undeleted constraints IN BINARY, if not infeasible
%            infeasible:  0 or 1
%            reduction:   0 or 1
%            iter:        number of iterations
%            DR:          dual recovery information, if option.DR = 1
%            time_preprocessing
%            time_total
%
% Information:
%   Created by Yuzixuan Zhu
%   Joint work with Gabor Pataki and Quoc Tran-Dinh
%   Department of Statistics and Operations Research, UNC-Chapel Hill.
%   Created date: August 31, 2017.
%   Last modified: January 25, 2018.
%   Contact: zyzx@live.unc.edu
%   More information: https://arxiv.org/pdf/1710.08954.pdf

if nargin < 2
   option.maxiter = intmax;
   option.epsilon = eps;
   option.cholEPS = 0;
   option.DR      = 0;
else
    if ~isfield(option, 'maxiter'), option.maxiter = intmax; end
    if ~isfield(option, 'epsilon'), option.epsilon = eps;    end
    if ~isfield(option, 'cholEPS'), option.cholEPS = 0;      end
    if ~isfield(option, 'DR'),      option.DR      = 0;      end
end
option.sqrtEPS = sqrt(option.epsilon);

time_total = tic;

% Choose the best PreProcessing function
n_pos = nnz(prob.blx == 0);
if n_pos == 0
    [probr, info] = PreProcessing_onlySDP(prob, option);
else
    [probr, info] = PreProcessing_withPOS(prob, option);
end

if option.DR == 0 && isfield(probr, 'DR')
    probr = rmfield(probr, 'DR');
end

info.time_total = toc(time_total);

end

% SeiveSDP v.2.0 by Yuzixuan Zhu, Gabor Pataki and Quoc Tran-Dinh.
% Copyright 2017 Department of Statistics and Operations Research
%                UNC - Chapel Hill, USA.
% See the file LICENSE for full license information.
