% This file is part of Presto-HF, a matlab toolbox to identify a circuit from
% its response.
%
% SPDX-License-Identifier: AGPL-3.0-or-later
%
% Copyright 2025 by
%   Centre Inria de l'Université Côte d'Azur
%   2004, route des Lucioles
%   06902 Sophia Antipolis Cedex
%
% and by
%   Mines Paris - Armines
%   60, boulevard Saint-Michel
%   75006 Paris
%
% Contributors: Fabien Seyfert, Jean-Paul Marmorat, Martine Olivi
%
% Presto-HF is free software: you can redistribute it and/or modify it under
% the terms of the GNU Affero General Public License as published by the Free
% Software Foundation, either version 3 of the License, or (at your option) any
% later version.
%
% Presto-HF is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
% A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
% details.
%
% You should have received a copy of the GNU Affero General Public License
% along with Presto-HF. If not, see <https://www.gnu.org/licenses/>.
%
function [l]=FirstNSmalestRadiusClusterOfNearestNeighbours(pv,M,N)

[si,sj]=size(M);
N=min(N,si);
Mf=[pv,M];
for k=1:N
    pind=SmalestRadiusClusterOfNearestNeighbours(Mf(:,1),Mf(:,2:end));
    l{k}=[Mf(sub2ind(size(Mf),pind,[1:sj+1]))];
    [siaux,sjaux]=size(Mf);
    Mf(sub2ind(size(Mf),pind,[1:sj+1]))=[];
    Mf=reshape(Mf,siaux-1,sjaux);
end
    

function [ppv_ind]=SmalestRadiusClusterOfNearestNeighbours(pv,M)
% Find the smallest raidus cluster of nearest neighbours of the points stored in the column pv and the points stored columnwise
% in M

l=NearestNeighbourClusters(pv,M);
r=[];
Mf=[pv,M];
for k=1:length(l)
    nump=length(l{k});
    r(k)=ClusterRadius(Mf(sub2ind(size(Mf),l{k},[1:nump])));
end
[aux,ind]=min(r);
ppv_ind=l{ind(1)};

function [l]=NearestNeighbourClusters(pv,M)
% Return the list of all the nearest neigbour clusters for all points running through pv
l={};
for k=1:length(pv)
    [ind]=NearestNeighbourCluster(pv(k),M);
    l{k}=[k,ind];
end

function [r]=ClusterRadius(pv)

r=0;
for k=1:length(pv)
    r=max([r,abs(pv-pv(k))]);
end
    

function [ind]=NearestNeighbourCluster(p,M)
% Returns de indices of the nearest neigbours of the set of points
% represented column by column by the matrix p
[aux,ind]=min(abs(M-p));