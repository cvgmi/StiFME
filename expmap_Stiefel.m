function [q] = expmap_Stiefel(X, V)


%   expmap_sphere maps the tangent vector v in T_{p}M to x on the unit
%   sphere.
%
%   x = expmap_sphere(p, v)
%
%   p and x are unit vectors (points on the unit sphere).
%   v is a tangent vector in T_{p}M.


if norm(V) < 1e-10
    q = X;
    return;
end

q = (eye(size(X,1))+V)/(eye(size(X,1))-V)*X;
end

