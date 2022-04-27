function [traj, cellref1] = get_trajs(cellList, mother, child, tstart, tmin, tmax)
ncell = numel(cellList);
cellList = reshape(cellList,[1,ncell]);
traj = cellList;

%tmin
%tmax

motherList = cellList;
if tmin<tstart
    for t = tstart:-1:tmin+1
        motherList = mother{t}(motherList);
        traj = [ motherList ; traj ];
    end
end

childrenList1 = cellList;
cellref1 = 1:ncell;
%size(child{t})
%size(cellList)
%max(childrenList1)

if tmax>tstart
    for t2 = tstart:1:tmax-1
        t2
        size(child{t2})
        childrenList2 = child{t2}(2,childrenList1);
        childrenList1 = child{t2}(1,childrenList1);
        cellref2 = find(childrenList2>1);
        trajref = 1:size(traj,2);
        cellref1 = [ cellref1 , cellref1(cellref2) ];        
        trajref = [ trajref , cellref2 ];
        childrenList2 = childrenList2(cellref2);
        childrenList1 = [ childrenList1 , childrenList2 ];
        traj = [ traj(:,trajref) ; childrenList1 ];
    end
end

Straj=size(traj)

