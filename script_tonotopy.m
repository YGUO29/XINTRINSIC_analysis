%% Tonotopy map
% -------------------
% Xindong's code below 
% -------------------
% R.PhPw_CycleAmp =      abs(        R.PhPwQt_FFTRaw(:,:,R.N_Qfc) )*2/R.N_Ft;
% R.PhPw_CycleAgl =      mod(angle(  R.PhPwQt_FFTRaw(:,:,R.N_Qfc) ), 2*pi);                      
% R.PtOne_CycleAmp =     reshape(R.PhPw_CycleAmp,   R.N_Pt, 1);
% R.PtOne_CycleAgl =     reshape(R.PhPw_CycleAgl, R.N_Pt, 1);  
% 
% % Compensate the pseudo-delay
% R.PtOne_Hue =          mod(    R.PtOne_CycleAgl - ...
%                                     T.PsuedoDelay/S.SesSoundDurTotal*2*pi, 2*pi)/(2*pi);
% % Reverse the hue for DOWN cycle
% if contains(lower(S.SesSoundFile), 'down')
%     R.PtOne_Hue =      1 - R.PtOne_Hue;
% end
% 
% % Phase cut    
%     R.PtOne_Hue =      (R.PtOne_Hue-S.TrlDurPreStim/S.TrlDurTotal) /...
%                                 (S.TrlDurStim/S.TrlDurTotal);
%                             % match the Stimulus ONSET / OFFSET to 0-1
%     R.PtOne_Hue =      min(max(R.PtOne_Hue, 0), 1);
%                             % limit the hue range as 0-HueLim
% % CONTINUOUS or DIScontinous hue    
% if contains(lower(S.SesSoundFile),  'sinusoidcycle')
%     R.HueLim =         1.0;
% else
%     R.HueLim =         0.8;
% end 
%     R.PtOne_Hue =      R.PtOne_Hue * R.HueLim;
% 
% % Construct the entire image       
% R.SaturationLim =  0.005;
% R.PtOne_Saturation =   min(R.PtOne_CycleAmp/R.SaturationLim,1);
% R.PtThree_TuneMap =    hsv2rgb([   R.PtOne_Hue,...
%                                         R.PtOne_Saturation,...
%                                         R.PtOne_Saturation]);
% R.PhPwThree_TuneMap =  reshape(R.PtThree_TuneMap, R.N_Ph, R.N_Pw,3);

%%