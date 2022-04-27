%David Pastor Escuredo. 2013-2015 BIT-UPM
%(c) All rights reserved

%The closes neighbour
function cellIndex=get_closest_cell(my_pos, cand)

    if nargin>2
        my_pos=my_pos+my_vel*deltaT;
    end
    cells=size(cand,2);
    
    dx=cand(1,:)-my_pos(1);
    dy=cand(2,:)-my_pos(2);
    dz=cand(3,:)-my_pos(3);
    
    d=sqrt(dx.^2+dy.^2+dz.^2);
    
    md=min(d);
    cellIndex=find(d==md);
    cellIndex=cellIndex(1);
 