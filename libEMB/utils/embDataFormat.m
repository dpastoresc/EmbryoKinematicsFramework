function embDataFormat(embfile,dataset,dataset_info,savepath)

fsep = filesep;
savepath = [savepath dataset '_t' fsep]
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
centA = fscanf(fid,'%d;%d;%d;%d;%d;%d;%d;%d', [8 inf]);
fclose(fid);

extent = dataset_info.extent;
spac = dataset_info.spac;

centA(3:5,:) = centA(3:5,:) + fix(extent/2)*ones(1,size(centA,2));
centA(3:5,:) = centA(3:5,:) .* ( spac*ones(1,size(centA,2)) );

%% get cent, mother, cellid, cellnum, cellvalid
nbstep = max(centA(6,:));
cent = cell(1,nbstep);
mother = cell(1,nbstep);
child = cell(1,nbstep);
cellid = cell(1,nbstep); 
cellnum = cell(1,nbstep);
cellvalid = cell(1,nbstep); 
for i  = 1:nbstep
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
    dataaux(:,1) = [-1 -1 0 0 0 i 1 -1];
    cent{i}      = dataaux(3:5,:);
    mother{i}    = dataaux(7,:);
    cellid{i}    = dataaux(1,:);
    cellnum{i}   = dataaux(2,:);
    cellvalid{i} = dataaux(8,:);
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
for t = 1:nbstep-1
    if isscalar(cent{t}) || isscalar(mother{t+1})
        child{t} = [1;1];
        continue
    end
    
    childaux = ones(2,size(cent{t},2));
    for i = 2:size(mother{t+1},2);
        candidate = mother{t+1}(i);
        if candidate>1
            if childaux(1,candidate) == 1
                childaux(1,candidate) = i;
            else
                childaux(2,candidate) = i;
            end
        end
    end
    child{t} = childaux;

end
child{nbstep} = ones(2,size(cent{nbstep},2));
savename = [dataset '_child.mat'];
save([savepath savename],'child');



