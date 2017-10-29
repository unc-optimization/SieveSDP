SieveSDP ([1]) is a preprocessing algorithm for semidefinite programming of the form
	min. <C, X>
	s.t. <A, X> == b
	     X in K
where K is the direct product of R^p, R^q_+ and S^r_+ (Euclidean space, nonnegative orthant, and positive semidefinite cone). For detail of using the code, call
	>> help SieveSDP

To test a problem in our datasets, go to folder “test examples”. There are 21 datasets saved as .zip files. After unzipping a dataset, you may see one or more SDP problems consisted in this dataset and saved as .mat files. They are in (Matlab-based) Mosek input format ([2, Section 9.7]).

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

[1] Sieve-SDP: a simple facial reduction algorithm to preprocess
semidefinite programs. https://arxiv.org/pdf/1710.08954.pdf
[2] http://docs.mosek.com/7.0/toolbox/A_guided_tour.html
[3] http://docs.mosek.com/8.1/matlabfusion/supported-file-formats.html
[4] http://sedumi.ie.lehigh.edu/sedumi/files/sedumi-downloads/SeDuMi_Guide_11.pdf
