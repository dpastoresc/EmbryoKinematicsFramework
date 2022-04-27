function embDataFormat(embfile,dataset,dataset_info,savepath,numCol,zoff)

fsep = filesep;
savepath = [savepath dataset '_t' fsep];
mkdir(savepath) % est ce que ca marche sous mac?
copyfile(embfile,savepath,'f');
save([savepath dataset '_info.mat'],'dataset_info');

fid = fopen(embfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% il est possible que le format d'entrer change
% on utilise la cache de l'appli de thierry pour récupérer les données
% certaines versions de l'appli peuvent renvoyer des caches un peu différent
% à terme  ça serait mieux de se référer aux données telles qu'elles sont dans
% la base (ie. x,y,z en micron cellid_mother plutot que l'index des tableaux
% de donné de l'application de thierry.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cellid,cellnum,x,y,z,t,mothernumb,validCell
if numCol==8
centA = fscanf(fid,'%d;%d;%d;%d;%d;%d;%d;%d', [numCol inf]);
else
  %centA = fscanf(fid,'%d;%d;%d;%d;%d;%d;%d;%d;%d', [numCol inf]);  
  centA=dlmread(embfile, ';');
  centA=centA';
end
fclose(fid);
size(centA)
extent = dataset_info.extent
spac = dataset_info.spac

% pause
% max(centA(3:5,2:end)')
% min(centA(3:5,2:end)')
% pause

%nbstep=nbstep+1;
if zoff
centA(3:5,:) = centA(3:5,:) + fix(extent/2)*ones(1,size(centA,2));
centA(3:5,:) = centA(3:5,:) .* ( spac*ones(1,size(centA,2)) );
else
    extentxy=extent(1:2)
    spacxy=spac(1:2)
    centA(3:4,:) = centA(3:4,:) + fix(extentxy/2)*ones(1,size(centA,2));
    centA(3:4,:) = centA(3:4,:) .* ( spacxy*ones(1,size(centA,2) ));
    centA(5,:) = centA(5,:) .* spac(3);
end

%% get cent, mother, cellid, cellnum, cellvalid
size(centA)
nbstep = max(centA(6,:))
%unique(centA(6,:))
initstep = min(centA(6,:))


%nbstep=nbstep+1;
cent = cell(1,nbstep+1);
mother = cell(1,nbstep+1);
child = cell(1,nbstep+1);
cellid = cell(1,nbstep+1); 
cellnum = cell(1,nbstep+1);
cellvalid = cell(1,nbstep+1); 
 motherid = cell(1,nbstep+1);
for i  = initstep:nbstep
    dataaux = [];
    if ~isempty(centA)
        cellSelect = centA(6,:)==i; % il y a plus rapide à faire
        if ~isempty(cellSelect)
            dataaux = centA(:,cellSelect);
            dataaux(7,:) = dataaux(7,:)+1; % car la valeur référence les ligne d'un tableau en C (ref a 0 pour le 1er elt)
                                            % en matlab la ref au 1er elt se fait à 1.
            centA(:,cellSelect)= []; % il y a plus rapide à faire
        end
    end
    if numCol==8
    dataaux(:,1) = [-1 -1 0 0 0 i 1 -1];
    else
    dataaux(:,1) = [-1 -1 0 0 0 i 1 -1 0]; % now the mother global id is added   
    
    motherid{i+1}  = dataaux(9,:);
    savename = [dataset '_motherid.mat'];
    save([savepath savename],'motherid');
    end
    s=i+1;
    cent{s}      = dataaux(3:5,:);
    mother{s}    = dataaux(7,:);
    cellid{s}    = dataaux(1,:);
    cellnum{s}   = dataaux(2,:);
    cellvalid{s} = dataaux(8,:);
end
savename = [dataset '_XYZ.mat'];
save([savepath savename],'cent');
savename = [dataset '_mother.mat'];
save([savepath savename],'mother');
savename = [dataset '_cellid.mat'];
save([savepath savename],'cellid');
savename = [dataset '_cellnum.mat'];
save([savepath savename],'cellnum');  
savename = [dataset '_cellvalid.mat'];
save([savepath savename],'cellvalid');
clear('cellid','cellnum');


%% get child file
nbstep
for t = 1:nbstep
    if isscalar(cent{t}) || isscalar(mother{t+1})
        child{t} = [1;1];
        continue
    end
    
    childaux = ones(2,size(cent{t},2));
    t
    size(cent{t},2)
    %     first row is control
    for i = 2:size(mother{t+1},2);
        candidate = mother{t+1}(i);
        % 1 because it has been stored in a matlab basis
        if candidate>1 
            % it assigns a matlab index or two
            if childaux(1,candidate) == 1
                childaux(1,candidate) = i;
            else
                childaux(2,candidate) = i;
            end
        end
    end
    child{t} = childaux;

end
child{nbstep+1} = ones(2,size(cent{nbstep+1},2));

savename = [dataset '_child.mat'];
save([savepath savename],'child');



