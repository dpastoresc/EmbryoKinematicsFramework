function prepareData(embfile, datasetFolder, t0, dt, extent1, extent2, extent3, sp1, sp2, sp3, validSteps, type)

addpath('../data/');
fsep = filesep;
datapath = load_centData('datapath');
% r�pertoire des donn�s en relatif par rapport au r�pertoire des scripts

% proc�dure � r�aliser pour un nouveau tracking
embfilepath = [datapath 'embfiles' fsep];
%embfilename = '070418a_tn10040t10040.emb';
%embfilename = '080528aF_tn2003t2003.emb';
%embfilename = '071226a_tn20046t20046.emb';
%embfilename= '110921aF_tn20184t20184.emb'
%embfilename= '111215aF_tn20186t20186.emb'
embfile = [embfilepath embfilename];

%dataset='070418a';
%dataset = '080528aF_tn2003'; % pr�fixe utilis� pour nommer les fichier et les appel�s ensuite 
%dataset = '111215aF_tn20186'; % pr�fixe utilis� pour nommer les fichier et les appel�s ensuite 

% data dataset
dataset_info = struct(...
    't0',6, ... t en heure du d�but de l'acquisition
    'dt',0.0386,... t en heure entre 2 acquisition successive
    'extent',[512;512;120],... en pixel
    'spac',[0.84; 0.84; 1.68],... en micron
    'type','L6-3',... embryon orientation (D: dorsal, AP: animal view, V: vegetal view) facultatif
    'nbstepValid',109);%761);  % timestep pas de temps jusqu'o� la reconstruction est correct, � d�faut nombre de pas de temps 
    

% lancer la proc�dure pour formater un .emb file
embDataFormat(embfile,datasetFolder,dataset_info,datapath);



