% David Pastor Escuredo. 2012/2015 BIT-UPM
% Tracking Kinematics Framework
% (C) All rights reserved

%Movit Selections Metadata for Zebrafish
function selCtag=selection2tissue(clusterSelec, extraTag)
tissues={'Hypoblast', 'Eye', 'Epiblast', 'All'}
selCtag=''
for sssc=1:length(clusterSelec)
    if sssc>1
        selCtag=[selCtag '-' tissues{clusterSelec(sssc)}]  
    else
        selCtag=[selCtag tissues{clusterSelec(sssc)}] 
    end
end
selCtag=[selCtag extraTag]