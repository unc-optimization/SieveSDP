SieveSDP ([1]) is a preprocessing algorithm for semidefinite programming of the form
	min. <C, X>
	s.t. A(X) == b
	     X in K
where K is the direct product of R^p, R^q_+ and S^r_+ (Euclidean space, nonnegative orthant, and positive semidefinite cones). For detail of using the code, call
	>> help SieveSDP;
	
To preprocess a problem using SieveSDP, call
	>> [probr, info] = SieveSDP(prob);

To test a problem in our datasets, go to folder “test examples”. There are 20 datasets saved as .zip files. After unzipping a dataset, you may see one or more SDP problems consisted in this dataset and saved as .mat files. They are in (Matlab-based) Mosek input format ([2, Section 9.7]).

To load a problem ``Example1.mat”in Matlab, call
	>> prob = load(``path/Example1.mat”);

To solve this problem using Mosek in Matlab, call
	>> [rcode, res] = mosekopt(``minimize info”, prob);

To convert this problem to other formats supported by Mosek, e.g., Task format, in order to run it outside of Matlab (see [3]), call
	>> mosekopt(['min write(‘, path, ‘/Example1.task.gz)'], prob);

To convert a problem from other formats supported by Mosek, e.g. Task format, call
	>> [~, res] = mosekopt(['read(', path, ‘/Example1.task.gz)']);
	>> prob = res.prob;

To convert this problem to SeDuMi format ([4]) in Matlab, call
	>> [A, b, c, K] = convert_mosek2sedumi(prob);

To convert a problem from SeDuMi format, call
	>> prob = convert_sedumi2mosek(A, b, c, K);

Some notes about solution recovery after solving the reduced problem:
To recover original solution X_original from X_reduced, call
	>> X_original = zeros(n, n);
	>> X_original(info.nonzero, info.nonzero) = X_reduced;
where n is the order of X, and info.nonzero is an output of SieveSDP. This recovery is accurate.

The dual problem is
	max. <b, y>
	s.t. A^* (y) <= C
To recover original solution y_original from y_reduced, call
	>> y_original = zeros(m, 1);
	>> y_original(info.undeleted) = y_reduced;
where m is the length of y, and info.undeleted is an output of SieveSDP. This recovery may be infeasible for badly-behaved SDPs ([5]). We are currently developing more accurate dual recovery code.
	
[1] Y. Zhu, G. Pataki, TD Quoc. Sieve-SDP: a simple facial reduction algorithm to preprocess semidefinite programs. https://arxiv.org/pdf/1710.08954.pdf
[2] http://docs.mosek.com/7.0/toolbox/A_guided_tour.html
[3] http://docs.mosek.com/8.1/matlabfusion/supported-file-formats.html
[4] http://sedumi.ie.lehigh.edu/sedumi/files/sedumi-downloads/SeDuMi_Guide_11.pdf
[5] G. Pataki. Bad semidefinite programs: they all look the same[J]. SIAM Journal on Optimization, 2017, 27(1): 146-172.
