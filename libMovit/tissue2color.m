% David Pastor Escuredo. 2012/2015 BIT-UPM
% Tracking Kinematics Framework
% (C) All rights reserved

%Movit Selections Metadata for Zebrafish
function my_color=tissue2color(selectNumber)
tissues={'Hypoblast', 'Eye', 'Epiblast', 'All'}
tissues_colors=[255 0 0; 225 168 3; 0 0 255; 255 0 255]
%0 0 0.8; 0 1 0; 1 0.5 0; 1 0 1
my_color=tissues_colors(selectNumber,:)/255;
