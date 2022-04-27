% David Pastor (2012): 
% Mechanics Framework
% Taken from Benoit's framework.... redo this at some point


fsep = filesep;  
%datapath = load_centData('datapath')
embfilepath = [datapath 'embfiles' fsep];
embfile = [embfilepath embfilename];

% data dataset
dataset_info = struct(...
    't0',t0_hours, ... t en heure du début de l'acquisition
    'dt',dt_hours,... t en heure entre 2 acquisition successive
    'extent',sizePixel',... en pixel
    'spac',spacingPixel',... en micron
    'type',species,... embryon orientation (D: dorsal, AP: animal view, V: vegetal view) facultatif
    'nbstepValid',vsteps);%761);  % timestep pas de temps jusqu'où la reconstruction est correct, à défaut nombre de pas de temps 
 
dataset_info
dataset
datasetBioemergences
   
if indexZero
   % lancer la procédure pour formater un .emb file
   embDataFormat_index0(embfile,datasetBioemergences,dataset_info,datapath, numCol, zoff)
else
   embDataFormat(embfile,datasetBioemergences,dataset_info,datapath); 
end

