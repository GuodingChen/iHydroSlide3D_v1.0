function [cSM,channel,overland,interflow,infil,CRESTABLE3D_output,caET,cP,sr,cERO,cERI,cPOT, tot_infil] = ...
    iHydroSlide3D(rundata,cSM,Pin,PET,dt,step, tstep,cstep,compgrid,dx,...
    NEXTTIMEUNDER,TH,ALPHA,BETA,SLOPE_TH,tot_infil,LandslideData, DownscalingData)
%Coupled Routing Excess STorage (CREST) distributed hydrologic model
%One way overland coupling (from water balance to routing only)
%CRESTb is the version where both rainfall and potential evapotanspiration
%data are loaded within the model.
%
%Version 1.0 - March, 2014
%Written by Humberto Vergara - humber@ou.edu
%Hydrometeorology and Remote-Sensing (HyDROS) Laboratory
%Advanced Radar Research Center (ARRC)
%The University of Oklahoma
%
%Reference:
%Wang. J., Y. Hong, L. Li, J.J. Gourley, K. Yilmaz, S. I. Khan, F.S. 
%Policelli, R.F. Adler, S. Habib, D. Irwn, S.A. Limaye, T. Korme, 
%and L. Okello, 2011: The Coupled Routing and Excess STorage (CREST) 
%distributed hydrological model. Hydrol. Sciences Journal, 56, 84-98.
%
%Arguments
%rundata = Struct variable containing simulation run settings
%Pin = Precipitation grid
%PET = Potential Evapotranspiration grid
%dt = Time step in seconds
%tstep = Time step in hours
%compgrid = Computational grid (mask)
%dx = horizontal length (pixel size)

%CREST Parameters
PRain=rundata.Parameters.PRain;
PKE = rundata.Parameters.PKE; %Multiplier to convert PET
PIM = rundata.Parameters.PIM.*0.01; %The impervious area ratio, 0 - 1.
PWM = rundata.Parameters.PWM; %The maximum soil water capacity (depth integrated pore space) of the soil layer. (mm)
K = rundata.Parameters.PFC .* tstep; %The soil saturated hydraulic conductivity (Ksat ,mm/hr).
PB = rundata.Parameters.PB; %The exponent of the variable infiltration curve. 

%Routing
overland = rundata.States.overland;
interflow = rundata.States.interflow;
channel = rundata.States.channel;
rundata = rmfield(rundata, 'States');

%pre-allocate temporary variables
cI = zeros(size(compgrid));
cER = zeros(size(compgrid));
cERO = zeros(size(compgrid));
cERI = zeros(size(compgrid));
ctemX = zeros(size(compgrid));
cperc = zeros(size(compgrid));
W0 = zeros(size(compgrid));

%Water Balance
%if (P(cstep) <= aET(cstep))
%Convert precipitation from mm/hr into mm
cP = (Pin .* tstep .*PRain) + compgrid; 

%Convert PET from mm/hr into mm, and compute actual ET
cPOT = PKE.*PET;
caET = PKE.*(PET.*tstep) + compgrid; 

%This is the precip that makes it to the soil.
cPSoil = (cP-caET).* (1 - PIM);

%Deal with grids that don't get Precip (i.e. PSoil <= 0)
cI(cPSoil <= 0) = 0;
cER(cPSoil <= 0) = 0;
cERO(cPSoil <= 0) = 0;
cERI(cPSoil <= 0) = 0;
ctemX(cPSoil <= 0) = (caET(cPSoil <= 0) - cP(cPSoil <= 0)) .* ...
    cSM(cPSoil <= 0) ./ PWM(cPSoil <= 0);

W0((cPSoil <= 0) & ...
    (ctemX < cSM)) = cSM((cPSoil <= 0) & (ctemX < cSM)) - ctemX((cPSoil <= 0) & (ctemX < cSM));
W0(W0 > PWM) = PWM(W0 > PWM);

caET(cPSoil <= 0) = cSM(cPSoil <= 0) - W0(cPSoil <= 0);
cPSoil(cPSoil <= 0) = 0;

%else

%Deal with grids that do get Precip (i.e. PSoil > 0)
%A. Deal with "unsaturated grids"
vicPars.b = PB((cPSoil > 0) & (cSM < PWM));
vicPars.Wm = PWM((cPSoil > 0) & (cSM < PWM));
%cER - Excess Rainfall
[cI((cPSoil > 0) & (cSM < PWM)),cER((cPSoil > 0) & (cSM < PWM)),W0((cPSoil > 0) & ...
    (cSM < PWM)),cperc((cPSoil > 0) & (cSM < PWM))] = Derive_vic(cPSoil((cPSoil > 0) & ...
    (cSM < PWM)),vicPars,cSM((cPSoil > 0) & (cSM < PWM)), interflow((cPSoil > 0) & ...
    (cSM < PWM)));

%B. Deal with "saturated" grids
cER((cPSoil > 0) &(cSM >= PWM)) = cPSoil((cPSoil > 0) & (cSM >= PWM));
W0((cPSoil > 0) & (cSM >= PWM)) = PWM((cPSoil > 0) & (cSM >= PWM));
cI((cPSoil > 0) & (cSM >= PWM)) = 0;
cperc((cPSoil > 0) & (cSM >= PWM)) = 0;

%ER Partitioning into Overland Flow and Interflow
%EF5 Implementation
ctemX(cPSoil > 0) = ((cSM(cPSoil > 0) + W0(cPSoil > 0)) ./ ...
    (2.*PWM(cPSoil > 0))) .* K(cPSoil > 0); %Calculate how much water can infiltrate
cERI((cPSoil > 0) & (cER <= ctemX)) = cER((cPSoil > 0) & (cER <= ctemX));
cERI((cPSoil > 0) & (cER > ctemX)) = ctemX((cPSoil > 0) & (cER > ctemX));
cERI(cERI < 0) = 0; %check why getting negative values
cERO(cPSoil > 0) = cER(cPSoil > 0) - cERI(cPSoil > 0) + (cP(cPSoil > 0) - caET(cPSoil > 0)) .* ...
    PIM(cPSoil > 0);
cERO(cERO < 0) = 0; %check why getting negative values

%Values for States after Water Balance
cSM = W0;

%Perform distributed routing
%Invoke routing (Kinematic Wave + interflow linear reservoir)
[channel,overland,interflow] = Flow_routing(cERO.*(0.001/dt),cERI,channel,overland,NEXTTIMEUNDER,ALPHA,BETA,TH,SLOPE_TH,dt,dx,rundata);
infil = cI;
tot_infil = tot_infil + cI;
sr = cSM ./ PWM;
% fs=factorS_Macon(rundata.slopeDeg, rundata.Parameters.cohesion, ...
%     rundata.Parameters.friction, rundata.Parameters.porosity, ...
%     rundata.Parameters.hydrocond, cstep, tot_infil, sr);
%         Ncores = rundata.Ncores;
%         Parallel_environment = parpool(Ncores);


    if ismember(datestr(step, 'yyyymmddHH'), rundata.observed_time_serise) == 1
        SM_12_5 = Soil_downscaling(cSM, rundata, LandslideData, DownscalingData);
%         Ncores = rundata.Ncores;
%         Parallel_environment = parpool(Ncores);
        CRESTABLE3D_output = LandslideModel(rundata, cstep, tot_infil, SM_12_5, PWM, LandslideData);
    else
%         delete(Parallel_environment)
        CRESTABLE3D_output = [];
    end
    
end





