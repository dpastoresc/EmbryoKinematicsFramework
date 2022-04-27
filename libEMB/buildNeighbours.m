%David Pastor Escuredo
%Mechanics Framework 2014
%(C) All rights reserved
%Code imported from analysis of Benoit Lombardot (CNRS, ISC-PIF, France)

%r is the spatial scale to consider the neighbouring

function [neigh2 dist2 foundNeighs] = buildNeighbours(posList, neighList, r, maxNeigh,spac)

%posList: centers of all data in this step
%neighList: valid ones plus yoursefl in the first row
%r is Xscale = 2*sigmaSpace
lPos = size(posList,2);

%Distance maps for neighs
%distance but only for the max of neighs we accept
dist2 = zeros(maxNeigh, lPos);
%matrix of neighs selected and all pos
neigh2 = 1*ones(maxNeigh+1, lPos);
neigh2(1,:)=0;%first is yoursefl

%one is embedded in the other... 
 List = [posList,neighList];
% size(List)
 xm = min(List(1,:));
 ym = min(List(2,:));
 zm = min(List(3,:));
 xM = max(List(1,:));
 yM = max(List(2,:));
 zM = max(List(3,:));

% xm = min(posList(1,:));
% ym = min(posList(2,:));
% zm = min(posList(3,:));
% xM = max(posList(1,:));
% yM = max(posList(2,:));
% zM = max(posList(3,:));

ncellmax = single( ceil(ceil(r/spac(1))*ceil(r/spac(2))*ceil(r/spac(3))/ceil(r^2)) );

%ncellmax
%size(ncellmax)
%r
%spac

sx = ceil((xM-xm)/r)+1+2;
sy = ceil((yM-ym)/r)+1+2;
sz = ceil((zM-zm)/r)+1+2;

indx3D = zeros([sx sy sz ncellmax],'uint16');

neighx = fix((neighList(1,:)-xm)/r)+1+1;
neighy = fix((neighList(2,:)-ym)/r)+1+1;
neighz = fix((neighList(3,:)-zm)/r)+1+1;

%indx3D is like a distribution is a scale of Xscale for all neighs
for i = 1:size(neighList,2)
    x = neighx(i);  y = neighy(i);  z = neighz(i);
    indx3D(x, y, z, 1 ) = indx3D(x,y,z,1)+1;
    indx3D(x, y, z, 1 + indx3D(x,y,z,1) ) = i;
end

% size(indx3D)
% A=indx3D(:,:,6,1);
% size(A)
% max(A(:))
% imagesc(A)
% % test = indx3D(:,:,:,1);
% % disp(num2str([ncellmax max(test(:))],' ; %d') );

posx = fix((posList(1,:)-xm)/r)+1+1;
posy = fix((posList(2,:)-ym)/r)+1+1;
posz = fix((posList(3,:)-zm)/r)+1+1;


%we only consider the neigh positions to take our samples but we can reuse
%the non valid vectors (might be not good idea though as they were not
%smoothed temporally)
% lpos=size(neighList,2)
% posx = fix((neighList(1,:)-xm)/r)+1+1;
% posy = fix((neighList(2,:)-ym)/r)+1+1;
% posz = fix((neighList(3,:)-zm)/r)+1+1;


for i = 1:lPos
    x=posx(i);  y=posy(i);  z=posz(i);
    %take the NR next to our ref.
    voisinage = indx3D(x-1:x+1,y-1:y+1,z-1:z+1,2:end);
    size(voisinage);
    neighloc = voisinage(voisinage(:)>0);
    size(neighloc);
    
    pos = posList(:,i);
    %pos=neightList(:,i);
    pos = pos*ones(1,length(neighloc));
    dist =  sum((pos-neighList(:,neighloc)).^2,1).^0.5;
    %find the neighs within r.
    %R is the Xave. Again we take neighs as far as Xave but Xave=2*sigma
    I = find(dist<r);
    
    if ~isempty(I)
        I = I';
        [dist permute ]= sort(dist(I));
        I = I(permute);
        foundNeighs=length(I);
        %maxNeigh
        %max number of neights is satured.
        lI = min(maxNeigh,length(I));
        I = I(1:lI);
        neighloc = neighloc(I);
        % [ lI ; neighloc ]
        neigh2(1:lI+1,i) = [ lI ; neighloc ];%the first value is the number of neighs. strange to keep this
        %lI
        %neigh2(1,i)
        %size(neigh2)
        %pause
        dist2(1:lI,i) = dist(1:lI);
    else
        neigh2(1,i) = 0 ;%no neighs->0
    end   
end
