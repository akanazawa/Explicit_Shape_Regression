function T = bestfit_nonreflective_similarity(uv, xy, fast)
% Adapted from Matlab's built-in cp2tform().
% Given Mx2 matrices uv and xy, with corresponding rows being corresponding points, returns 3x2
% matrix T such that [xy ones(M,1)] * T == uv, with T a nonreflective similarity
% matrix.  M must be at least 2.  If greater than 2, the best fit T is returned.
%
% Written by Thomas Berg
%%%%%%%%%%%%%%%%%%%%
if (~exist('fast', 'var')) || isempty(fast)
    fast = true;
end

M = size(xy,1);
x = xy(:,1);
y = xy(:,2);
X = [x   y  ones(M,1)   zeros(M,1);
     y  -x  zeros(M,1)  ones(M,1)  ];

u = uv(:,1);
v = uv(:,2);
U = [u; v];

% We know that X * r = U
% Skip this assert for performance.
if ~fast
    assert(rank(X) >= 4, ...
           'At least 2 unique points needed to infer nonreflective similarity transform.');
end
r = X \ U;

sc = r(1);
ss = r(2);
tx = r(3);
ty = r(4);

T = [sc -ss;
     ss  sc;
     tx  ty];
