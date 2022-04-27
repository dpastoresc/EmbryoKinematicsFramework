%David Pastor Escuredo. 2013 BIT-UPM
%(c) All rights reserved

%we get a fake child because no cell was found, so we take the closest...
%will be improving this.
function childIndex=get_dynamic_child_vel(my_pos, my_vel, next_cand, next_cand_vel, deltaT)

    if nargin==5
        %we work with the predicted position of the cell in the next step
          my_pos=my_pos+my_vel*deltaT;
    end
    cells=size(next_cand,2);
    my_vel=my_vel(1:3);
    my_vel=my_vel./norm(my_vel);
    
    dx=next_cand(1,:)-my_pos(1);
    dy=next_cand(2,:)-my_pos(2);
    dz=next_cand(3,:)-my_pos(3);
    
    d=sqrt(dx.^2+dy.^2+dz.^2);
    
    [dn id]=sort(d);
    cand=id(1:5);
    dn=dn(1:5);
    dnm=dn-min(dn);
    dn=dn(dnm<5);
    cand=cand(dnm<5);
    for i=1:length(cand)
        cvel=next_cand_vel(1:3,cand(i));      
        cvel=cvel./norm(cvel);
        da(i)=acosd(dot(my_vel(1:3)', cvel'));
    end
    
    dam=da-min(da);
    dn=dn(dam<3);
    cand=cand(dam<3);
    
    ii=find(dn==min(dn));
    childIndex=cand(ii(1));
    %da(childIndex)
    %d(childIndex)
  
    
    
    