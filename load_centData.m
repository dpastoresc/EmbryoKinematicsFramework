%Function to load tracking structure data
%*it was moved here to make sure it works properly

function varargout = load_centData(varargin)

% charge les donn�es associ�e au jeu dataset
dataset = varargin{1};
if iscell(dataset)
    prefix = dataset{1};
    dataset = dataset{2};
else
    prefix = '';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% d�finition de datapath reprise dans les autres programmes %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fsep = filesep;                        %
datapath = ['../../MechanicsData/'];            %
datapath = ['C://Users\Jose\Desktop\MechanicsData\'];  
%datapath = ['../MechanicsData/']; 
%datapath='/home/public/MechanicsData';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(dataset,'datapath');  varargout{1}=datapath;  return;  end

path2 = [datapath prefix dataset '_t/']
infofile = [path2 dataset '_info.mat']
if exist(infofile,'file')
    a = load(infofile);
    a= a.dataset_info;
    t0 = a.t0; % en heure
    dt = a.dt; % en heure
    extent = a.extent; % en pixel le long de chaque dimension
    spac = a.spac; % dimension d'un voxel en micron dans chaque dimension
    type = a.type; % nomenclature de Nadine pour l'orientation des donn�es
    if  isfield(a,'manualSph'), manualSph = a.manualSph; end % 4 valeurs, centre de la sph�re fittant les donn�es et sont rayon en deriere position (en microns)
    if  isfield(a,'manualRef'), manualRef = a.manualRef; end % matrice de rotation
    if  isfield(a,'axisRange'), axisRange = a.axisRange; end % 6 valeur xmin,xmax,... en microns
    if  isfield(a,'viewAngle'), viewAngle = a.viewAngle; end % azimuth elevation (voir 'view' function)
    if  isfield(a,'nbstepValid'), nbstepValid = a.nbstepValid; end % nombre de pas de temps, � partir de 1, estimer ok pour les mesures 
else
    error(['pas de donn�e correspondant � ' dataset]);
end
timeline = t0:dt:t0+(nbstepValid-1)*dt ; % temps en heure de chaque timestep

%% output preparation

for i = 2:length(varargin)
    if strcmp(varargin{i},'cent')
        varargin{i} = 'XYZ';
    end    
end

%path2
filelist = ls(path2);
%length(filelist)
filelist = [ filelist char(32*ones(size(filelist,1),1)) ]; % on rajoute un espace � la fin de chaque ligne (necessaire � la compatibilit� windows)
filelist = reshape(filelist',[1,numel(filelist)]);
filelist(strfind(filelist,char(10))) = ' '; % on remplace les retour chariot par des espace (compatibilit� mac)
filelist(strfind(filelist,char(9))) = ' '; % on remplace les retour chariot par des espace (compatibilit� linux)
filelist = strsplit(' ',filelist)

%length(filelist)
%filelist
%varargin
for j = 2:length(varargin)
    varargout{j-1} = -1;
    i=1;
    %[dataset '_' varargin{j} '.mat']

    while i <= length(filelist)
        name = filelist{i} 
        %name=name(1:24)
       % length(name)
        vv=[dataset '_' varargin{j} '.mat'];
       % length(vv)
        strcmp([dataset '_' varargin{j} '.mat'],name)   ;     
        
        
        if strcmp(vv,name)
            %'hola'
            [path2 name];
            a = load( [path2 name] );
            a = struct2cell(a);
            size(a);
           % pause
            varargout{j-1} = a{1};
            size(a{1});
        %else
        %   'no match'
        end
        i=i+1;
    end
    
    if exist(varargin{j},'var')==1
        eval(['varargout{j-1}=' varargin{j} ';']);
    end
    
end



