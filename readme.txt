Sieve-SDP ([1]) is a preprocessing algorithm for semidefinite programming of the form

	min. <C, X>
	s.t. A(X) == b
	     X in K

where K is the direct product of R^p, R^q_+ and S^r_+ (Euclidean space, nonnegative orthant, and positive semidefinite cones). For detail, call

	>> help SieveSDP;



******* HOW TO USE SIEVE-SDP TO DO REDUCTION *******

To preprocess a problem using SieveSDP, call

	>> [probr, info] = SieveSDP(prob);

If info.infeasible == 1, then problem is infeasible. Otherwise, to solve a reduced problem by Mosek, call

	>> [rcode, res] = mosekopt(‘minimize’, prob);

The problem solution is saved in ‘res’, c.f. [2].



******* HOW TO RECOVER ORIGINAL SOLUTION ******

If solution of original (before-preprocessed) problem is desired, call

	>> [x_original, X_original] = recoveryPrimal(res, info);

where ‘info’ is given from the call to SieveSDP, ‘x_original’ corresponds to linear variables, and ‘X_original’ corresponds to PSD variables.

If solution of original dual problem

	max. <b, y>
	s.t. A^* (y) + Z = C
	     Z in K^*

is desired, set

	>> option.DR = 1;

and call SieveSDP by

	>> [probr, info] = SieveSDP(prob, option);

After solving the problem by Mosek, call

	[y_original, z_original, Z_original, info1] = recoveryDual(res, info);

Dual recovery may not always succeed ([5]). Its success status is saved in ‘info1’.



******* HOW TO FORMULATE A PROBLEM *******

To test a problem in our datasets, go to folder “test examples”. There are 20 datasets saved as .zip files. After unzipping a dataset, you will see SDP problems saved as .mat files. Sources of these problem sets are listed in [1]. They are in (Matlab-based) Mosek Matlab format ([2, Section 9.7]).

To load a problem ``Example1.mat”in Matlab, call

	>> prob = load(‘path/Example1.mat’);

To convert a problem to different formats supported by Mosek, e.g., Task format, in order to run it outside of Matlab ([3]), call

	>> mosekopt(['min write(Example1.task.gz)'], prob);

To convert a problem from other formats supported by Mosek, e.g. Task format, call

	>> [~, res] = mosekopt(['read(Example1.task.gz)']);
	>> prob = res.prob;

To convert a problem from Mosek format to SeDuMi format ([4]), call

	>> [A, b, c, K] = convert_mosek2sedumi(prob);

To convert a problem from SeDuMi format to Mosek format, call

	>> prob = convert_sedumi2mosek(A, b, c, K);



[1] Y. Zhu, G. Pataki, TD Quoc. Sieve-SDP: a simple facial reduction algorithm to preprocess semidefinite programs. https://arxiv.org/pdf/1710.08954.pdf
[2] http://docs.mosek.com/7.0/toolbox/A_guided_tour.html
[3] http://docs.mosek.com/8.1/matlabfusion/supported-file-formats.html
[4] http://sedumi.ie.lehigh.edu/sedumi/files/sedumi-downloads/SeDuMi_Guide_11.pdf
[5] G. Pataki. Bad semidefinite programs: they all look the same[J]. SIAM Journal on Optimization, 2017, 27(1): 146-172.
