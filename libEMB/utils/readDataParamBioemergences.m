% David Pastor Escuredo BIT-UPM
% 2012. All rights reserved
% Script to get data from bioemergences server
function [sizePixel spacingPixel steps t0 dt species ]=readDataParamBioemergences(datasetName)

    % e.g. readDataParamBioemergences('070418a')
    %all in one line
    %name:071226a; id:6200; operator:Nadine Peyrieras; species:DA.RE; microscopetype:SP5; obective:20X/0.95NA Olympus; scanning:200; 
    %zoom:1x; averaging:; exposure:; pixelsx:512; pixelsxy:512; pixelsxz:120; micronsx:1.37; micronsy:1.37; micronsz:1.37; ratioz/x:1.00; 
    %totalmicronsx:700; totalmicronsy:700; totalmicronsz:164; timesteps:231; deltat:153939; animalage:19440;
    
   % http://bioemergences.iscpif.fr/workflow
    m_url=['http://bioemergences.iscpif.fr/workflow/album/galerie/reports/services/experimentDetails.php?name=' datasetName]
    m_url=['http://bioemergences.iscpif.fr/workflow/services/experimentDetails.php?name=' datasetName]
   % m_url=['https://bioemergences.iscpif.fr/workflow/services/experimentDetails.php?name=' datasetName]
    datainfo=urlread(m_url)
    
    %datainfo='name:071226a; id:6200; operator:Nadine Peyrieras; species:DA.RE; microscopetype:SP5; obective:SP5; scanning:200; zoom:1; averaging:0; exposure:; pixelsx:512; pixelsxy:512; pixelsxz:120; micronsx:1.37; micronsy:1.37; micronsz:1.37; ratioz/x:1.00; totalmicronsx:700; totalmicronsy:700; totalmicronsz:164; timesteps:231; deltat:153939; animalage:19440;'
    remain=datainfo  
    
    token='init';
    while strcmp(token,'')==0
     [token remain]= strtok(remain, ';');    
     [token2 remain2]= strtok(token, ':');
     if strcmp(token2,'animalage')
         t0=str2num(remain2(2:end));   
     elseif strcmp(token2,'timesteps')
         steps=str2num(remain2(2:end));
     elseif strcmp(token2,'deltat')
         dt=str2num(remain2(2:end));
     elseif strcmp(token2,'species')
         species=remain2(2:end);
     elseif strcmp(token2,'pixelsx')
         sx=str2num(remain2(2:end));
     elseif strcmp(token2,'pixelsxy')
         sy=str2num(remain2(2:end));
     elseif strcmp(token2,'pixelsxz')
         sz=str2num(remain2(2:end));
     elseif strcmp(token2,'micronsx')
         spx=str2num(remain2(2:end));
     elseif strcmp(token2,'micronsy')
         spy=str2num(remain2(2:end));
     elseif strcmp(token2,'micronsz')
         spz=str2num(remain2(2:end));
     end    
         
     remain=remain(3:end);
          
    end
   
    sizePixel=[sx sy sz];
    spacingPixel=[spx spy spz];