% Ce script cr�e les fichiers matlab associ�s avec de nouveaux jeu de donn�
% les fichier et donn�e cr�er pourront �tre appel� charg� via la fonction
% load_centData
%
% pour que le script fonctionne on doit se placer dans le r�pertoire ou se
% trouve les scripts

% placer le fichier emb dans le r�pertoire embfiles
% renseigner les variables ci-dessous (embfilename, dataset et dataset)
%

fsep = filesep;
datapath = load_centData('datapath');
% r�pertoire des donn�s en relatif par rapport au r�pertoire des scripts

% proc�dure � r�aliser pour un nouveau tracking
embfilepath = [datapath 'embfiles' fsep];
%embfilename = '070418a_tn10040t10040.emb';
embfilename = '080528aFr_tn20201t20201.emb';
%embfilename = '071226a_tn20046t20046.emb';
%embfilename= '110921aF_tn20184t20184.emb'
%embfilename= '111215aF_tn20197t20197.emb'
embfile = [embfilepath embfilename];

%dataset='070418a';
dataset = '080528aFr_tn20201'; % pr�fixe utilis� pour nommer les fichier et les appel�s ensuite 
%dataset = '111215aF_tn20197'; % pr�fixe utilis� pour nommer les fichier et les appel�s ensuite 

% data dataset
dataset_info = struct(...
    't0',1, ... t en heure du d�but de l'acquisition
    'dt',0.05,... t en heure entre 2 acquisition successive
    'extent',[512;512;144],... en pixel
    'spac',[0.75; 0.75; 1.5],... en micron
    'type','Phal',... embryon orientation (D: dorsal, AP: animal view, V: vegetal view) facultatif
    'nbstepValid',100);%761);  % timestep pas de temps jusqu'o� la reconstruction est correct, � d�faut nombre de pas de temps 
    

% lancer la proc�dure pour formater un .emb file
embDataFormat(embfile,dataset,dataset_info,datapath);



