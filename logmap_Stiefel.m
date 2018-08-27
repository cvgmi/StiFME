function [V] = logmap_Stiefel(X, Q)


%   logmap_sphere maps x on sphere to the tangent space T_{p}M.
%
%   v = logmap_sphere(p, x)
%
%   p and x are unit vectors (points on the unit sphere).
%   v is a tangent vector in T_{p}M.


if (norm(X-Q)) < 1e-10
    V = zeros(size(X,1),size(X,1));
    return;
end
[n,p]=size(X);
M = Q(1:p,1:p)'*X(1:p,1:p)+X(p+1:n,1:p)'*Q(p+1:n,1:p);
sk = 0.5*(M'-M);
A = 2*inv(X(1:p,1:p)'+Q(1:p,1:p)')*sk*inv(X(1:p,1:p)+Q(1:p,1:p));
B = (Q(p+1:n,1:p)-X(p+1:n,1:p))*inv(X(1:p,1:p)+Q(1:p,1:p));
V = [A,-B';B,zeros(n-p,n-p)];

end

