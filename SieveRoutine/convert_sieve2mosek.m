% Get the orders of sdp variables
nj  = zeros(1, n_sdp);
for j = 1:n_sdp
    nj(j) = nnz(I_nonzero((Aind(j) + 1):Aind(j + 1)));
end
nj           = nj(nj > 0);
probr.bardim = nj;
n_sdp        = length(nj);
Aind         = zeros(1, n_sdp + 1);
for j = 1:n_sdp
    Aind(j + 1) = Aind(j) + nj(j);
end

% Convert c
if n_pos > 0
    probr.c = c_pos;
else
    probr.c = prob.c;
end
[row, col, v] = find(tril(c_convert));
count         = length(v);
subj          = zeros(1, count);
subk          = zeros(1, count);
subl          = zeros(1, count);
val           = v';

k = 1;
for vv = 1:count
    for j = k:n_sdp
        if col(vv) <= Aind(j + 1), break; end
    end
    subj(vv) = j;
    subk(vv) = row(vv) - Aind(j);
    subl(vv) = col(vv) - Aind(j);
    if j > k, k = j; end
end
probr.barc.subj = subj;
probr.barc.subk = subk;
probr.barc.subl = subl;
probr.barc.val  = val;

% Convert b
probr.blc = b(undeleted);
probr.buc = probr.blc;

% Convert A
if n_pos > 0
    probr.a = A_pos';
else
    probr.a = prob.a(undeleted, :);
end
subi    = cell(1, info.m_post);
subj    = cell(1, info.m_post);
subk    = cell(1, info.m_post);
subl    = cell(1, info.m_post);
val     = cell(1, info.m_post);
for ii = 1:info.m_post
    i = undeleted(ii);
    [row, col, v] = find(tril(A_convert{i}));
    count    = length(v);
    subi{ii} = ii*ones(1, count);
    subj{ii} = zeros(1, count);
    subk{ii} = zeros(1, count);
    subl{ii} = zeros(1, count);
    val{ii}  = v';
    k = 1;
    for vv = 1:count
        for j = k:n_sdp
            if col(vv) <= Aind(j + 1), break; end
        end
        subj{ii}(vv) = j;
        subk{ii}(vv) = row(vv) - Aind(j);
        subl{ii}(vv) = col(vv) - Aind(j);
        if j > k, k = j; end
    end
end
probr.bara.subi = horzcat(subi{1:info.m_post});
probr.bara.subj = horzcat(subj{1:info.m_post});
probr.bara.subk = horzcat(subk{1:info.m_post});
probr.bara.subl = horzcat(subl{1:info.m_post});
probr.bara.val  = horzcat(val{1:info.m_post});

% Convert x
if n_fre == 0
    if n_pos == 0
        probr.blx = [];
    else
        probr.blx = zeros(1, nnz(Ipos_nonzero));
    end
else
    probr.blx = -inf(1, n_fre);
end
probr.bux = [];
