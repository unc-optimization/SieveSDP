probr = [];
info  = [];
b     = prob.blc;
m     = length(b);  % Number of constraints

% Get the number of free variables
n_fre = nnz(prob.blx);

% Get the number of positive variables
n_pos = length(prob.blx) - n_fre;

% Get the orders of SDP variables
nj        = prob.bardim;
n_sdp     = length(nj);
n_sdp_sum = sum(nj);
n         = n_fre + n_pos + n_sdp_sum;
Aind      = zeros(1, n_sdp + 1);
for j = 1:n_sdp
    Aind(j + 1) = Aind(j) + nj(j);
end

info.n_pre = n;
info.m_pre = m;

% Convert c
c_pos = prob.c;
len = length(prob.barc.val);
subk  = zeros(1, len);
subl  = zeros(1, len);
val = prob.barc.val;
for j = 1:len
    addterm              = Aind(prob.barc.subj(j));
    subk(j)    = prob.barc.subk(j) + addterm;
    subl(j)    = prob.barc.subl(j) + addterm;
    if prob.barc.subk(j) == prob.barc.subl(j)
        val(j) = val(j)/2;
    end
end
subk1     = [subk, subl];
subl1     = [subl, subk];
val       = [val, val];
c_convert = sparse(subk1, subl1, val, n_sdp_sum, n_sdp_sum);

% Convert A
A_pos      = prob.a';
A_convert  = cell(m, 1);
[C, ia, ~] = unique(prob.bara.subi);
len        = length(C);
ia         = [ia; length(prob.bara.subi) + 1];
for i = 1:len
    len1  = ia(i + 1) - ia(i);
    subij = prob.bara.subj(ia(i):(ia(i + 1) - 1));
    subik = prob.bara.subk(ia(i):(ia(i + 1) - 1));
    subil = prob.bara.subl(ia(i):(ia(i + 1) - 1));
    val   = prob.bara.val(ia(i):(ia(i + 1) - 1));
    for j = 1:len1
       addterm  = Aind(subij(j));
       subik(j) = subik(j) + addterm;
       subil(j) = subil(j) + addterm;
       if subik(j) == subil(j)
           val(j) = val(j)/2;
       end
    end
    subik1 = [subik, subil];
    subil1 = [subil, subik];
    val    = [val, val];
    A_convert{C(i)} = sparse(subik1, subil1, val, n_sdp_sum, n_sdp_sum);
end
for i = 1:m
   if isempty(A_convert{i})
       A_convert{i} = sparse(n_sdp_sum, n_sdp_sum);
   end
end

% Convert b
I           = find(b > 0);
b(I)        = -b(I);
A_pos(:, I) = -A_pos(:, I);
len = length(I);
for i = 1:len
    A_convert{I(i)} = -A_convert{I(i)};
end
