function [ dist ] = distance_stiefel( X, Y )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

v=logmap_Stiefel(X, Y);
dist=norm(reshape(v, [size(v,1)*size(v,2),1]),'fro');


end

