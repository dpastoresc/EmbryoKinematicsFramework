%%%Benoit Lombardot's (CREA, Polytechnique) calculation of Gradient of Velocity field abd other tensors. 2011
%%%Revised by David Pastor (BIT, UPMadrid) 2012

%Attention: this version includes radial correction and spatialAve flags

function ch2_call_MovStrainCalc_dav(dataset,tscale,Xscale,Rscale,maxNeigh,minNeigh,savepath, radialCorrection, spatialAve)

fsep = filesep; % separateur de fichier sur le system
%David: Here it reads the data from the structures. Don't touch this
[cent,mother,child,spac,sph,tl,tmax0] = load_centData(dataset,'cent','mother','child','spac','manualSph','timeline','nbstepValid');

%     params = struct('tscale',40,...     % minutes
%                     'Xscale',50,...     % microns
%                     'Rscale',15,...      % microns
%                     'XmaxNeigh',200,... % cell number
%                     'XminNeigh',5  ...  % cell number
%                     );

%%%David: This defines the window of analysis!
params = struct('tscale'   , tscale   ,...  % minutes, tscale =2*sigma_t
                'Xscale'   , Xscale   ,...  % microns, Xscale =2*sigma_X
                'Rscale'   , Rscale   ,...  % microns 1/2 largeur échelon R0+/-Rscale
                'XmaxNeigh', maxNeigh ,...  % cell number
                'XminNeigh', minNeigh  ...  % cell number
               );
           
%%%David: It calculates things in terms of hours. Time steps preparation
dt = tl(2)-tl(1); % temps entre 2 timesteps en heures
Dt = ceil( params.tscale/(60*dt) );
tmax = min( tmax0 , length(cent)-(Dt+1) );
tList = 1:tmax;
datas = {cent,mother,child,spac,sph,tl};
clear('cent','mother','child','spac','sph','tl');

%%%David: Initializing Mov / MeanMov / TenseurGrad
Mov=cell(1,tmax);
MeanMov = cell(1,tmax);
TenseurGrad = cell(1,tmax);
Neighall = cell(1,tmax);
for t=tList
    %tic
    %I added NeighAux to see the density changes in neighbours along time
    [MovAux MeanMovAux TenseurGradAux NeighAux] = ch2_MovStrainCalc2_dav(datas,t,params,radialCorrection, spatialAve);
    Mov{t} = MovAux;
    MeanMov{t} = MeanMovAux;
    %Here it calculates the gradient
    TenseurGrad{t} = TenseurGradAux;
    Neighall{t} = NeighAux;
    %disp([dataset ' , t=' int2str(t) ' , calc time = ' num2str(toc,'%f')]);
end
t1 = datas{6}(1);
dt1  = datas{6}(2)-datas{6}(1);

%%David: Here he does the estimation of error. I remove it because i dont
%%really get it and is according to 070418 GS. Also it gives erros
%meanMovErr = ch2_meanMovErr_calc(dt1,t1,MeanMov);
%[Eall Eerrall Rall Fall] = ch2_E_Eerr_R_calc(TenseurGrad,meanMovErr);
[Eall Rall Fall] = ch2_E_Eerr_R_calc_dav(TenseurGrad);%No error stuff

%David: save tensors as matlab cells
if spatialAve==1
    filename = [dataset '_Fave_T' int2str(params.tscale) '_X' int2str(params.Xscale) '.mat'];
else
    filename = [dataset '_F_T' int2str(params.tscale) '_X' int2str(params.Xscale) '.mat'];
end
save([savepath fsep dataset '_t' fsep filename],'Fall');
filename = [dataset '_Mov_T' int2str(params.tscale) '_X' int2str(params.Xscale) '.mat'];
save([savepath fsep dataset '_t' fsep filename],'Mov');
filename = [dataset '_MeanMov_T' int2str(params.tscale) '_X' int2str(params.Xscale) '.mat'];
save([savepath fsep dataset '_t' fsep filename],'MeanMov');
%filename = [dataset '_MeanMov_T' int2str(params.tscale) '_X' int2str(params.Xscale) '_Err.mat'];
%save([savepath fsep dataset '_t' fsep filename],'meanMovErr');
filename = [dataset '_TenseurGrad_T' int2str(params.tscale) '_X' int2str(params.Xscale) '.mat'];
save([savepath fsep dataset '_t' fsep filename],'TenseurGrad');
if spatialAve==1
    filename = [dataset '_Eave_T' int2str(params.tscale) '_X' int2str(params.Xscale) '.mat'];
else
    filename = [dataset '_E_T' int2str(params.tscale) '_X' int2str(params.Xscale) '.mat'];
end
save([savepath fsep dataset '_t' fsep filename],'Eall');
%filename = [dataset '_E_T' int2str(params.tscale) '_X' int2str(params.Xscale) '_Err.mat'];
%save([savepath fsep dataset '_t' fsep filename],'Eerrall');
if spatialAve==1
    filename = [dataset '_Rave_T' int2str(params.tscale) '_X' int2str(params.Xscale) '.mat'];
else
    filename = [dataset '_R_T' int2str(params.tscale) '_X' int2str(params.Xscale) '.mat'];
end
save([savepath fsep dataset '_t' fsep filename],'Rall');
filename = [dataset '_Neigh_T' int2str(params.tscale) '_X' int2str(params.Xscale) '.mat'];
save([savepath fsep dataset '_t' fsep filename],'Neighall');

