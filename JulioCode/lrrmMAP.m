%
% LINEAR REFERENCE REGON MODEL BY CARDENAS-RODRIGUEZ AND PAGEL
%
%         [MAPS,curves,ROIs]=lrrmMAP(DCE_MRI,time,mask) %%
%% INPUTS :
%       DCE_MRI:        3D Matrix of DCE-MRI data, with: 
%               dimension_1= rows
%               dimension_2=colums
%               dimension_3=concentration
%               dimension_4=time (minutes)
%
%       time:          column vector of time
%              The length of time should be equal to the lenght of
%              DCE_MRI(1,1,1,:)=dimension_4
%
%       mask:        Structure with the following fields: 
%               .RR= binary mask for the reference region
%               .TOI= binary mask for the TOI
%
%% OUTPUTS :
%       MAPS:    Structure with the following components
%      MAPS.RKtrans=    Parametric map of Relative Ktrans (RKtrans)
%      MAPS.kepTOI=     Parametric map of kep for the tissue of interest (TOI)
%      MAPS.kepRR=      Parametric map of kep for the region of interest (RR)
%
%       curves:    Structure with the following components:
%      curves.reference_region=     average signal/concentration for the RR
%      curves.tissue_of_interest=   average signal/concentration for the TOI
%      curves.RKtrans=              Relative Ktrans (RKtrans)
%      curves.kepTOI=               kep for the tissue of interest (TOI)
%      curves.kepRR=                kep for the region of interest (RR)
%      curves.RSQ=                 RSQ
%
%       mask:    Structure with the following components
%      mask.TOI=    Binary mask for the Tissue of Interest
%      mask.RR=     Binary mask for the Region of Interest;

%% Reference for this method:  
%   Cardenas-Rodriguez J, Howison CM, Pagel MD.
%   A linear algorithm of the reference region model for DCE-MRI is robust
%   and relaxes	requirements for temporal resolution.		
%   Magn Reson Imaging. 2013 May;31(4):497-507. 
%
%
% Julio Cardenas-Rodriguez,     Mark D. Pagel
% cardenaj@email.arizona.edu    mpagel@email.arizona.edu
% ver 2.0 July, 2013.
%%
function              [MAPS,curves,mask]=lrrmMAP(DCE_MRI,time,mask) 

%control of variable;
dce=DCE_MRI;
%gaussFilter = fspecial('gaussian',[7 7],1);
%dce= imfilter(dce,gaussFilter,'same');
% 

%% 1) Determines size of image data and stores values into row, cols and time points. 

[rows,cols,timepoints]=size(dce);


%% 2) Calculate Average concentrion for each ROI as function of time.

curves.reference_region=nan(timepoints,1);     %Preallocate
curves.tissue_of_interest=curves.reference_region;   %Preallocate

for j=1:timepoints;
    SLICE=dce(:,:,j);   
            curves.reference_region(j,1)=mean(SLICE(mask.RR)); 
                    curves.tissue_of_interest(j,1)=mean(SLICE(mask.TOI));
end

%% 3) Estimate parameters for the Curves using the LRRM
rr=smooth(curves.reference_region);
toi=smooth(curves.tissue_of_interest);

[pars,~]=lrrm_nonneg(time,rr,toi); 

curves.RKtrans=     pars(1);              
curves.kepTOI=      pars(2);                 
curves.kepRR=       pars(3);               
curves.RSQ=        pars(4);              

%% 4) Calculate Parametric Maps

% Preallocate Maps;
MAPS.RKtrans=    nan(rows,cols);
MAPS.kepTOI=     MAPS.RKtrans;
MAPS.kepRR=      MAPS.RKtrans;
MAPS.RSQ=        MAPS.RKtrans;

% Find indices for all voxes INSIDE that belong to the TOI
indices=find(mask.TOI>0);

% Reshape DCE Data Matrix to turn into a 2D Matrix
% Colums now have the same length than time= timepoints;
dce=reshape(dce,[],timepoints);

for q=1:length(indices)
    
    toi=squeeze(dce(indices(q),:))';
    
                                        [pars,~]=lrrm_nonneg(time,rr,toi);
    
            MAPS.RKtrans(indices(q))=     pars(1);
            MAPS.kepTOI(indices(q))=      pars(2);
            MAPS.kepRR(indices(q))=       pars(3);
            MAPS.RSQ(indices(q))=        pars(4);
    
end
%%
function [pars,prediction] = lrrm_lsq(time,Ct_RR,Ct_TOI)
%
%       [pars,Ct_TOI] = lrrm_lsq(time,Ct_RR,Ct_TOI)
%
%% INPUTS :
%       time:    column vector of time
%       Ct_RR:  concentration in the Reference Region as a function of time
%       Ct_TOI: concentration in the Tissue of Interest as a function of time
%
%% OUTPUTS :
%       R:    Structure with the following components
%           R(1)=  Relative Ktrans (RKtrans)
%           R(2)=  kep for the tissue of interest, kepTOI
%           R(3)=  kep for the reference region,   kepRR
%           R(4)=  RMSE for the fitting
%
%       Ct_TOI:    predicted curve for the TOI
% 
%

% Reference:
%   Cardenas-Rodriguez J, Howison CM, Pagel MD.
%   A linear algorithm of the reference region model for DCE-MRI is robust
%   and relaxes	requirements for temporal resolution.		
%   Magn Reson Imaging. 2013 May;31(4):497-507. 
%
%
% Julio Cardenas-Rodriguez,     Mark D. Pagel
% cardenaj@email.arizona.edu    mpagel@email.arizona.edu
% ver 2.0 July, 2013.
%
%% Define vectors and matrices for this method (see reference)

y=Ct_TOI;
    x1=Ct_RR;      
        x2=cumtrapz(time,Ct_RR); 
            x3=cumtrapz(time,Ct_TOI);
                X=[x1,x2,-x3,ones(size(y))];

%%  Estimate parameters using a LLSQ
%(Please not that this is a little different than using \ as described in the reference
%  but we haev discovered that is more reliabel for low SNE.
%B = robustfit(X,y);
 B=X\y;
 %B=lsqnonneg(X,y);
%%  Calculate each parameter of output pars (see reference) ##

pars(1,1)=B(1);             %RKtrans;
pars(2,1)=B(3);             %kepTOI;
pars(3,1)=(B(2))*B(1)^-1;   %kepRR;

%% Calculate predicted curve for the TOI

prediction=X*B; %#ok<*NBRAK>
%prediction=B*X; %#ok<*NBRAK>

% Calculate the the RMSE for the predicted curve and allocated to pars

pars(4,1)=rsquare(y,prediction);


%% Ithaca
% 
% When you set out for Ithaka
% ask that your way be long,
% full of adventure, full of instruction.
%
% The Laistrygonians and the Cyclops,
% angry Poseidon - do not fear them:
% such as these you will never find
% as long as your thought is lofty, as long as a rare
% emotion touch your spirit and your body.
% The Laistrygonians and the Cyclops,
% angry Poseidon - you will not meet them
% unless you carry them in your soul,
% unless your soul raise them up before you.  
%                                                   ....................
% Constantine P. Cavafy

end

%%

end