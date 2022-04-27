%David Pastor Escuredo. 2013 BIT-UPM
%(c) All rights reserved

%The closest neighbour to the predicted position
function childIndex=get_dynamic_child(my_pos, next_cand, my_vel, deltaT)

    if nargin>2
        my_pos=my_pos+my_vel*deltaT;
    end
    cells=size(next_cand,2);
    
    dx=next_cand(1,:)-my_pos(1);
    dy=next_cand(2,:)-my_pos(2);
    dz=next_cand(3,:)-my_pos(3);
    
    d=sqrt(dx.^2+dy.^2+dz.^2);
    
    md=min(d);
    childIndex=find(d==md);
    childIndex=childIndex(1);
 