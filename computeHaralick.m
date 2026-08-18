function F = computeHaralick(GLCM)
%COMPUTEHARALICK Haralick texture features from a grey-level co-occurrence matrix.
%
%   F = COMPUTEHARALICK(GLCM) returns 13 of the texture features defined by
%   Haralick, Shanmugam and Dinstein (1973) for the two-dimensional
%   co-occurrence matrix GLCM. GLCM is expected to be normalized so that its
%   entries sum to 1, i.e. to hold probabilities rather than raw counts.
%
%   Inputs
%     GLCM  - G-by-G normalized co-occurrence matrix, for example one slice
%             of the output of CALCULATE3DGLCM.
%
%   Outputs
%     F     - 1-by-13 row vector of features, in this order:
%               1  Energy (angular second moment)
%               2  Correlation
%               3  Inertia (contrast)
%               4  Entropy
%               5  Inverse difference moment (homogeneity)
%               6  Sum average
%               7  Sum variance
%               8  Sum entropy
%               9  Difference average
%              10  Difference variance
%              11  Difference entropy
%              12  Information measure of correlation 1
%              13  Information measure of correlation 2
%
%             This ordering is the contract between this function and the
%             CSV header written by GENERATEDATA_ALLFILTER. Do not reorder
%             the outputs without updating that header as well.
%
%   Example
%     vol  = normalizeAndScale(randn(64, 64, 64));
%     G    = calculate3DGLCM(vol, [1 0 0], 256);
%     F    = computeHaralick(G(:, :, 1));
%
%   Notes
%     - Only 13 features are returned. Haralick's original set has 14; the
%       maximal correlation coefficient, which requires an eigenvalue
%       decomposition, is not computed here.
%     - Entropy terms use SPFUN(@LOG2, ...) so that zero entries are skipped
%       rather than producing -Inf. Sparse-style indexing with FIND is used
%       throughout for the same reason, which is why the code operates on
%       the nonzero triplets [i, j, v] of the matrix.
%     - Grey level g is stored at index g+1, so the code subtracts 1 from
%       every index before using it as a grey-level value.
%     - If GLCM is not normalized to unit sum, the entropy-based and
%       correlation-based features are not on their intended scale. See the
%       normalization WARNING in CALCULATE3DGLCM.
%     - A constant or near-degenerate GLCM can make SIGMA_X or SIGMA_Y zero,
%       in which case Correlation is NaN.
%
%   Reference
%     R. M. Haralick, K. Shanmugam and I. Dinstein, "Textural Features for
%     Image Classification", IEEE Transactions on Systems, Man, and
%     Cybernetics, vol. SMC-3, no. 6, pp. 610-621, 1973.
%
%   See also CALCULATE3DGLCM, NORMALIZEANDSCALE, GRAYCOPROPS.

%   Authors:     Tingfei Men, Ahmed Ashraf, Sherif S. Sherif
%   Affiliation: Department of Electrical and Computer Engineering,
%                University of Manitoba, Winnipeg, Manitoba, Canada
%   Contact:     Tingfei Men <tingfeimen@gmail.com>
%
%   Citation:    If you use this code in a scientific publication, please cite
%                T. Men, A. Ashraf, and S. S. Sherif, "Generation of a large
%                three-dimensional texture image dataset for machine and deep
%                learning applications," accepted to the International Workshop
%                on New Approaches to Multivariate Statistical Analysis
%                (NAMSP 2026).
%   Project:     3D Synthetic Texture Dataset Generation
%   Created:     2024-12-09
%   Last update: 2026-08-17
%   License:     MIT. See the LICENSE file in the repository root.

% Calculating the texture features from the SGLD matrix
foo=GLCM;

%Entropy
entropy=sum(sum(-((full(spfun(@log2,foo))).*foo)));

% Energy:
energy=sum(sum(foo.*foo));

% Inertia:
[i,j,v]=find(foo);
inertia=sum((((i-1)-(j-1)).*((i-1)-(j-1))).*v);

% Inverse differnece moment:
inverse_diff=sum((1./(1+(((i-1)-(j-1)).*((i-1)-(j-1))))).*v);

% Correlation:
[m,n]=size(foo);

px=sum(foo,2);
[i,j,v]=find(px);
mu_x=sum((i-1).*v);
sigma_x=sum((((i-1)-mu_x).^2).*v);
h_x=sum(sum(-((full(spfun(@log2,px))).*px)));
temp1=repmat(px,[1 m]);

py=sum(foo,1);
[i,j,v]=find(py);
mu_y=sum((j-1).*v);
sigma_y=sum((((j-1)-mu_y).^2).*v);

h_y=sum(sum(-((full(spfun(@log2,py))).*py)));
temp2=repmat(py,[n 1]);


[i,j,v]=find(foo);
correlation=(sum(((i-1)-mu_x).*((j-1)-mu_y).*v))/sqrt(sigma_x*sigma_y);

% Information measures of correlation 1 and 2:
foo1=-(foo.*(((temp1.*temp2)==0)-1));
foo2=-((temp1.*temp2).*((foo1==0)-1));
[i1,j1,v1]=find(foo1);
[i2,j2,v2]=find(foo2);
h1=sum((sum(-(v1.*(log2(v2))))));
info_corr_1=(entropy-h1)/max(h_x,h_y);
[i,j,v]=find(temp1.*temp2);
h2=sum((sum(-(v.*(log2(v))))));
info_corr_2=sqrt((1-exp(-2*(h2-entropy))));

% Sum average, variance and entropy:
[i,j,v]=find(foo);
k=i+j-1;
pk_sum=zeros(max(k),1);
for l=min(k):max(k)
  pk_sum(l)=sum(v(find(k==l)));
end

[i,j,v]=find(pk_sum);
sum_avg=sum((i-1).*v);
sum_var=sum((((i-1)-sum_avg).^2).*v);
sum_entropy=sum(-((full(spfun(@log2,pk_sum))).*pk_sum));
 
% Difference average, variance and entropy:
[i,j,v]=find(foo);
k=abs(i-j);
pk_diff=zeros(max(k)+1,1);
for l=min(k):max(k)
   pk_diff(l+1)=sum(v(find(k==l)));
end

[i,j,v]=find(pk_diff);
diff_avg=sum((i-1).*v);
diff_var=sum((((i-1)-diff_avg).^2).*v);
diff_entropy=sum(-((full(spfun(@log2, pk_diff))).*pk_diff));   


F= [energy,correlation,inertia,entropy,inverse_diff,sum_avg,...
    sum_var,sum_entropy,diff_avg,diff_var,diff_entropy,info_corr_1,...
    info_corr_2];

% Alternative ordering used in earlier revisions, kept for reference:
% F = [energy,inertia,correlation,inverse_diff,sum_avg,sum_var,sum_entropy,...
%      entropy,diff_var,diff_entropy,info_corr_1,info_corr_2,diff_avg];

end
