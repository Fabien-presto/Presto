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