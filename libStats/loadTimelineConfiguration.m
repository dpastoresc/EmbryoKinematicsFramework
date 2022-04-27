% David Pastor Escuredo. 2012/2015 BIT-UPM
% Tracking Kinematics Framework
% (C) All rights reserved

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Timeline %%%%%%%%%%%%%%%%%%%%%

% Some datasets may have been shifted to be at the same step in the
% temporal landmakr of alignment
shift_t0=0
% Shift dataset if necessary by changing t0 dinamically

%We establish limits by hpf times
if hini<0
    hini=6
end
if hfin<0
    hfin=13
end

if hpf_limits
    t_start=uint16((hini*3600-t0)/dt_sec)
    if t_start<tini
        t_start=tini
    end
    t_limit=uint16((hfin*3600-t0)/dt_sec)
    if t_limit>tfinal
        t_limit=tfinal
    end
end

if startInTref
    t_start=tinilag
    t_limit=t_start+140
end
