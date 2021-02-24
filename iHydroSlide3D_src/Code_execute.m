function [prediction] = Code_execute(rundata)
% use readgrids.m
    %The Main function to run the modeling system.
    %
    %
    %
    %Arguments
    %"rundata" is a struct variable containing all information about the
    %simulation run.

    switch lower(rundata.modelCore)
        case 'hydromodel'
            fprintf('-----------------------------------------------------\n');
            fprintf('      solo Hydromodel     \n');
            fprintf('  MATLAB R2018b platform - November 2020   \n');
            fprintf('        Version Beta 1.0   \n');
            fprintf('-----------------------------------------------------\n');
        case 'hydroslide'
            fprintf('-----------------------------------------------------\n');
            fprintf('  Coupled Hydromodel and SLIDE (HydroSlide)  \n');
            fprintf('     MATLAB R2018b platform - November 2020  \n');
            fprintf('       Version Beta 1.0    \n');
            fprintf('-----------------------------------------------------\n');
        case 'ihydroslide3d'
            fprintf('^.^ ^.^ Welcome to use iHydroSlide3D ^.^ ^.^ \n');
            fprintf('  Coupled Hydromodel and 3D Landslide model \n');
            fprintf('  MATLAB R2018b platform - November 2020 \n');
            fprintf('        Version Beta 1.0            \n');
            fprintf(['Parallel Ncores = ', num2str(rundata.Ncores)]);
            fprintf('-----------------------------------\n');
            
        otherwise
            error(['Fatal error: \n' , ...
                   'The model core specified by you doesn''t exist!\n', ...
                   'Your selection is ', rundata.modelCore '.']);
    end

    fprintf('Read and pre-process input data ...\n');

    %Read gridded data
    [rundata.dem,rundata.slope,rundata.facc,rundata.fdir,rundata.gridres ...
    ,rundata.Parameters,rundata.States,rundata.outlet, LandslideData, DownscalingData] = Hydro_BasicData_Read(rundata);
    
 
    %Were basic data previously generated?
    if (isfield(rundata, 'load_basicgrids') == 1)
        fprintf('Loading previously saved basic grids: ');
        load(rundata.load_basicgrids);

        %Store loaded data in rundata
        rundata.HorLen = HorLen;
        rundata.routto = routto;
        rundata.routfrom = routfrom;
        rundata.routfrom_n = routfrom_n;
        rundata.slope = slope;
        rundata.seqmatrix = seqmatrix;
        rundata.Parameters.ALPHA = ALPHA;
        rundata.Parameters.BETA = BETA;

        %Free some memory
        clear HorLen routto routfrom routfrom_n slope seqmatrix ALPHA BETA;

        fprintf('Done.\n');
    else
        %Generate computational grids
        [rundata.HorLen,rundata.routto,rundata.routfrom,rundata.routfrom_n ...
        ,rundata.slope,rundata.seqmatrix,rundata.Parameters.ALPHA ...
        ,rundata.Parameters.BETA] = HPro_GridDistance(rundata);
        
        if (isfield(rundata, 'save_basicgrids') == 1)
            fprintf('Saving basic grids: ');

            %Set temporary variables
    %         dem = rundata.dem;
    %         facc = rundata.facc;
    %         fdir = rundata.fdir;
    %         gridres = rundata.gridres;
    %         Parameters = rundata.Parameters; 
    %         States = rundata.States;
    %         outlet = rundata.outlet;

            HorLen = rundata.HorLen;
            routto = rundata.routto;
            routfrom = rundata.routfrom;
            routfrom_n = rundata.routfrom_n;
            slope = rundata.slope;
            seqmatrix= rundata.seqmatrix;
            ALPHA = rundata.Parameters.ALPHA;
            BETA = rundata.Parameters.BETA;

            %Write to disk then free some memory
            save(rundata.save_basicgrids, 'HorLen','routto','routfrom' ...
                ,'routfrom_n','slope','seqmatrix','ALPHA','BETA', '-v7.3');
            clear HorLen routto routfrom routfrom_n slope seqmatrix ALPHA BETA;

            fprintf('Done.\n');
        end
    end
    
    
    


    %Devired soil and land-cover related parameters for coupled hydro and
    %slide models only
    switch lower(rundata.modelCore)
        case {'hydromodel', 'hydroslide', 'ihydroslide3d'}   
            rundata.slopeDeg=atand(rundata.slope);
            %Soil land cover parameters
            [cohesion,friction,porosity,hydrocond]=readSoilLandPara(rundata);
            rundata.Parameters.cohesion=cohesion;
            rundata.Parameters.friction=friction;
            rundata.Parameters.porosity=porosity;
            rundata.Parameters.hydrocond=hydrocond;
    end


    %Simulation Settings
    tstep = rundata.tstep; %time step in hrs
    tperiod = rundata.tperiod; %time period for simulation
    dt = rundata.tstep * 3600; %rundata.tstep is in hours, dt in seconds
    %Stream length (grid side size)
    dx = rundata.HorLen; %in meters

    Parameters = rundata.Parameters;
    ALPHA=Parameters.ALPHA;
    BETA=Parameters.BETA;
    SLOPE_TH=Parameters.Slope_TH;
    TH=Parameters.TH;
    under=Parameters.UNDER;
    
    %Interflow crossing-time grid
    under = under.*(0.001/dt); %(conversion to meters per second)
    sq_slope = sqrt(rundata.slope);
    speedU = under.*sq_slope;
    nexTimeUnder = dx./speedU;
    
    [mr,mc,~] = size(rundata.routto);
    compgrid = reshape(rundata.routto(:,:,1),mr,mc); 
    compgrid(isnan(compgrid) == 0) = 0;

    States_0 = rundata.States;

    %Pre-allocate variables
    % prediction.Q = zeros([size(compgrid),length(tperiod)]);
    nextra = length(rundata.outlet.extra.y);
    prediction.R = zeros(1+nextra,length(tperiod));
    prediction.Rain = zeros(1+nextra,length(tperiod));
    prediction.EPot = zeros(1+nextra,length(tperiod));
    prediction.W = zeros(1+nextra,length(tperiod));
    prediction.SM = zeros(1+nextra,length(tperiod));
    prediction.RS = zeros(1+nextra,length(tperiod));
    prediction.RI = zeros(1+nextra,length(tperiod));
    prediction.EAct = zeros(1+nextra,length(tperiod));
    prediction.Infil = zeros(1+nextra,length(tperiod));
    prediction.FS = zeros(1+nextra,length(tperiod));
    prediction.Tot_infil = zeros(1+nextra,length(tperiod));
    %Initial condition
    inpars=rundata;
    inpars.sq_slope=sq_slope;
    
    switch lower(rundata.modelCore)
        case {'hydromodel', 'hydroslide', 'ihydroslide3d'}
            %cSM = States_0.IWU.*Parameters.PWM./100; 
            cSM = States_0.IWU;
            cSM(cSM>Parameters.PWM)=Parameters.PWM(cSM>Parameters.PWM);
            
            %Soil Moisture in mm
    end
    
    %Main Loop
    tot_infil=zeros(size(compgrid));
    
    cstep = 1;
    for step = tperiod(1):(tstep/24):tperiod(end)
        % Read Pin and PET
        if strcmpi(rundata.rainFormat,'asc')==1
            filename=[rundata.pp_path ...
                datestr(step,'yyyymmddHH') '.asc'];
            [Pin, ~] = read_ascii_gis(filename);
        elseif strcmpi(rundata.rainFormat,'txt')==1
            filename=[rundata.pp_path ...
                datestr(step,'yyyymmddHH') '.txt'];
            [Pin, ~] = read_ascii_gis(filename);
        elseif strcmpi(rundata.rainFormat,'tif')==1
            filename=[rundata.pp_path ...
                datestr(step,'yyyymmddHH') '.tif'];
            [Pin, ~] = geotiffread(filename); 
        end
        
        
        if strcmpi(rundata.petFormat,'asc')==1
            filename=[rundata.pet_path ...
                datestr(step,'yyyymmddHH') '.asc'];
            [PET, ~] = read_ascii_gis(filename); 
        elseif strcmpi(rundata.petFormat,'txt')==1
            filename=[rundata.pet_path ...
                datestr(step,'yyyymmddHH') '.txt'];
            [PET, ~] = read_ascii_gis(filename); 
        elseif strcmpi(rundata.petFormat,'tif')==1
            filename=[rundata.pet_path ...
                datestr(step,'yyyymmddHH') '.tif'];
            [PET, ~] = geotiffread(filename); 
        end
        
        
        switch lower(rundata.modelCore)
            case 'hydromodel'
                % Run solo hydrological model
                [cSM,channel,overland,interflow,infil,cP,sr,excS,excI,cPOT,caET] = Hydromodel(inpars,cSM,Pin,...
                    PET,dt,tstep,compgrid,dx,nexTimeUnder,TH,ALPHA, ...
                    BETA,SLOPE_TH);
                
                
                inpars.States.overland=overland;
                inpars.States.interflow=interflow;
                inpars.States.channel=channel;

                Qq = channel;
                Qq(inpars.facc < Parameters.TH) = overland(inpars.facc < Parameters.TH);  

                % Store other variables of interest
                prediction.R(1,cstep) = Qq(rundata.outlet.y,rundata.outlet.x);
                prediction.Rain(1,cstep) = mean(cP(rundata.mask>0));
                prediction.EPot(1,cstep) = mean(cPOT(rundata.mask>0));
                prediction.EAct(1,cstep) = mean(caET(rundata.mask>0));
%                 prediction.RS(1,cstep) = mean(overland(rundata.mask>0));
%                 prediction.RI(1,cstep) = mean(interflow(rundata.mask>0));
                prediction.W(1,cstep) = mean(cSM(rundata.mask>0));
                prediction.SM(1,cstep) = mean(sr(rundata.mask>0));
                prediction.Infil(1,cstep) = mean(infil(rundata.mask>0));
                
                for nout = 2:size(prediction.R,1)
                    prediction.R(nout,cstep) = Qq(rundata.outlet.extra.y(nout-1),rundata.outlet.extra.x(nout-1));
                    prediction.Rain(nout,cstep) = mean(cP(rundata.outlet.masks{nout-1}>0));
                    prediction.EPot(nout,cstep) = mean(cPOT(rundata.outlet.masks{nout-1}>0));
                    prediction.EAct(nout,cstep) = mean(caET(rundata.outlet.masks{nout-1}>0));
%                     prediction.RS(nout,cstep) = mean(overland(rundata.outlet.masks{nout-1}>0));
%                     prediction.RI(nout,cstep) = mean(interflow(rundata.outlet.masks{nout-1}>0));
                    prediction.W(nout,cstep) = mean(cSM(rundata.outlet.masks{nout-1}>0));
                    prediction.SM(nout,cstep) = mean(sr(rundata.outlet.masks{nout-1}>0));
                    prediction.Infil(nout,cstep) = mean(infil(rundata.outlet.masks{nout-1}>0));
                end
%                 prediction.channel = channel;
%                 prediction.overland = overland;
%                 prediction.interflow = interflow;   
            case 'hydroslide'
                % Run HydroSlide
                [cSM,channel,overland,interflow,infil,sr,fs,caET,cP,excS,excI,cPOT, tot_infil] = HydroSlide(...
                    inpars,cSM,Pin,PET,dt,tstep,cstep,compgrid,dx, ...
                    nexTimeUnder,TH,ALPHA,BETA,SLOPE_TH,tot_infil);

                if cstep==1
                    fsmin=fs;
                else
                    fsmin(fsmin>fs)=fs(fsmin>fs); 
                end
                
                inpars.States.overland=overland;
                inpars.States.interflow=interflow;
                inpars.States.channel=channel;

                Qq = channel;
                Qq(inpars.facc < Parameters.TH) = overland(inpars.facc < Parameters.TH);  
                % Store other variables of interest
                prediction.R(1,cstep) = Qq(rundata.outlet.y,rundata.outlet.x);
                prediction.Rain(1,cstep) = mean(cP(rundata.mask>0));
                prediction.EPot(1,cstep) = mean(cPOT(rundata.mask>0));
                prediction.EAct(1,cstep) = mean(caET(rundata.mask>0));
%                 prediction.RS(1,cstep) = mean(overland(rundata.mask>0));
%                 prediction.RI(1,cstep) = mean(interflow(rundata.mask>0));
                prediction.W(1,cstep) = mean(cSM(rundata.mask>0));
                prediction.SM(1,cstep) = mean(sr(rundata.mask>0));
                prediction.Infil(1,cstep) = mean(infil(rundata.mask>0));
                prediction.FS(1,cstep) = nanmean(fs(rundata.mask>0));
                
                prediction.Tot_infil(1,cstep) = nanmean(tot_infil(rundata.mask>0));
                
                for nout = 2:size(prediction.R,1)
                    prediction.R(nout,cstep) = Qq(rundata.outlet.extra.y(nout-1),rundata.outlet.extra.x(nout-1));
                    prediction.Rain(nout,cstep) = mean(cP(rundata.outlet.masks{nout-1}>0));
                    prediction.EPot(nout,cstep) = mean(cPOT(rundata.outlet.masks{nout-1}>0));
                    prediction.EAct(nout,cstep) = mean(caET(rundata.outlet.masks{nout-1}>0));
%                     prediction.RS(nout,cstep) = mean(overland(rundata.outlet.masks{nout-1}>0));
%                     prediction.RI(nout,cstep) = mean(interflow(rundata.outlet.masks{nout-1}>0));
                    prediction.W(nout,cstep) = mean(cSM(rundata.outlet.masks{nout-1}>0));
                    prediction.SM(nout,cstep) = mean(sr(rundata.outlet.masks{nout-1}>0));
                    prediction.Infil(nout,cstep) = mean(infil(rundata.outlet.masks{nout-1}>0));
                    prediction.FS(nout,cstep) = nanmean(fs(rundata.outlet.masks{nout-1}>0));
                    prediction.Tot_infil(nout,cstep) = nanmean(tot_infil(rundata.outlet.masks{nout-1}>0));
                end
                
                
            case 'ihydroslide3d'
               
                % Run iHydroSlide3D
                [cSM,channel,overland,interflow,infil,CRESTABLE3D_output,caET,cP,sr,excS,excI,cPOT, tot_infil] = ...
                    iHydroSlide3D(inpars,cSM,Pin,PET,dt,step, tstep,cstep,compgrid,dx, ...
                    nexTimeUnder,TH,ALPHA,BETA,SLOPE_TH,tot_infil,LandslideData, DownscalingData);


                inpars.States.overland=overland;
                inpars.States.interflow=interflow;
                inpars.States.channel=channel;

                Qq = channel;
                Qq(inpars.facc < Parameters.TH) = overland(inpars.facc < Parameters.TH);  
                % Store other variables of interest
                prediction.R(1,cstep) = Qq(rundata.outlet.y,rundata.outlet.x);
                prediction.Rain(1,cstep) = mean(cP(rundata.mask>0));
                prediction.EPot(1,cstep) = mean(cPOT(rundata.mask>0));
                prediction.EAct(1,cstep) = mean(caET(rundata.mask>0));
    %                 prediction.RS(1,cstep) = mean(overland(rundata.mask>0));
    %                 prediction.RI(1,cstep) = mean(interflow(rundata.mask>0));
                prediction.W(1,cstep) = mean(cSM(rundata.mask>0));
                prediction.SM(1,cstep) = mean(sr(rundata.mask>0));
                prediction.Infil(1,cstep) = mean(infil(rundata.mask>0));
                prediction.Tot_infil(1,cstep) = nanmean(tot_infil(rundata.mask>0));

                for nout = 2:size(prediction.R,1)
                    prediction.R(nout,cstep) = Qq(rundata.outlet.extra.y(nout-1),rundata.outlet.extra.x(nout-1));
                    prediction.Rain(nout,cstep) = mean(cP(rundata.outlet.masks{nout-1}>0));
                    prediction.EPot(nout,cstep) = mean(cPOT(rundata.outlet.masks{nout-1}>0));
                    prediction.EAct(nout,cstep) = mean(caET(rundata.outlet.masks{nout-1}>0));
    %                     prediction.RS(nout,cstep) = mean(overland(rundata.outlet.masks{nout-1}>0));
    %                     prediction.RI(nout,cstep) = mean(interflow(rundata.outlet.masks{nout-1}>0));
                    prediction.W(nout,cstep) = mean(cSM(rundata.outlet.masks{nout-1}>0));
                    prediction.SM(nout,cstep) = mean(sr(rundata.outlet.masks{nout-1}>0));
                    prediction.Infil(nout,cstep) = mean(infil(rundata.outlet.masks{nout-1}>0));
                    prediction.Tot_infil(nout,cstep) = nanmean(tot_infil(rundata.outlet.masks{nout-1}>0));
                end

    %                 prediction.channel = channel;
    %                 prediction.overland = overland;
    %                 prediction.interflow = interflow;  
    %                 prediction.aET=caET;

            %     prediction.channel(:,:,cstep) = flow_routed.channel;
            %     prediction.overland(:,:,cstep) = flow_routed.overland;
            %     prediction.overland_h(:,:,cstep) = flow_routed.overland_h;
            %     prediction.interflow(:,:,cstep) = flow_routed.interflow;
            %     prediction.QO(:,:,cstep) = QO;
            %     prediction.QI(:,:,cstep) = QI;
            %     prediction.P(:,:,cstep) = cP;
            %     prediction.PSoil(:,:,cstep) = cPSoil;
            %     prediction.aET(:,:,cstep) = caET;
            %     prediction.I(:,:,cstep) = cI;
            %     prediction.perc(:,:,cstep) = cperc;
            %     prediction.ER(:,:,cstep) = cER;
            %     prediction.ERI(:,:,cstep) = cERI;
            %     prediction.ERO(:,:,cstep) = cERO;
            %     prediction.temX(:,:,cstep) = ctemX;
    %                 prediction.infil(:,:,cstep)=infil;
    %                 prediction.sr(:,:,cstep)=sr;
    %                 prediction.fs(:,:,cstep)=fs;
    %                 prediction.caET(:,:,cstep)=caET;
    %                 prediction.P(:,:,cstep)=cP;
                
        end
        
        %Output results
        strdate = datestr(step,'yyyymmddHH');
        if rundata.output_runoff
            data=overland+interflow;
            data(isnan(data) | rundata.mask<0)=rundata.nodata;
            data(channel>0)=channel(channel>0);
            filename=[rundata.outPath 'GOVar_R_' strdate '.' lower(rundata.resultFormat)];
            output_results(filename,rundata.resultFormat,rundata.geoinfo,rundata.projinfo,rundata.info,data);
        end
        if rundata.output_W    
            data=cSM;
            data(isnan(data) | rundata.mask<0)=rundata.nodata;
            filename=[rundata.outPath 'GOVar_W_' strdate '.' lower(rundata.resultFormat)];
            output_results(filename,rundata.resultFormat,rundata.geoinfo,rundata.projinfo,rundata.info,data); 
        end
        if rundata.output_SM  
            data=sr;
            data(isnan(data) | rundata.mask<0)=rundata.nodata;
            filename=[rundata.outPath 'GOVar_SM_' strdate '.' lower(rundata.resultFormat)];
            output_results(filename,rundata.resultFormat,rundata.geoinfo,rundata.projinfo,rundata.info,data); 
        end
%         if rundata.output_ExcS 
%             data=excS;
%             data(isnan(data) | rundata.mask<0)=rundata.nodata;
%             filename=[rundata.outPath 'GOVar_ExcS_' strdate '.' lower(rundata.resultFormat)];
%             output_results(filename,rundata.resultFormat,rundata.geoinfo,rundata.projinfo,rundata.info,data); 
%         end
%         if rundata.output_ExcI  
%             data=excI;
%             data(isnan(data) | rundata.mask<0)=rundata.nodata;
%             filename=[rundata.outPath 'GOVar_ExcI_' strdate '.' lower(rundata.resultFormat)];
%             output_results(filename,rundata.resultFormat,rundata.geoinfo,rundata.projinfo,rundata.info,data);
%         end
        if rundata.output_rain  
            data=cP;
            data(isnan(data) | rundata.mask<0)=rundata.nodata;
            filename=[rundata.outPath 'GOVar_Rain_' strdate '.' lower(rundata.resultFormat)];
            output_results(filename,rundata.resultFormat,rundata.geoinfo,rundata.projinfo,rundata.info,data);
        end
        if rundata.output_epot  
            data=cPOT;
            data(isnan(data) | rundata.mask<0)=rundata.nodata;
            filename=[rundata.outPath 'GOVar_EPot_' strdate '.' lower(rundata.resultFormat)];
            output_results(filename,rundata.resultFormat,rundata.geoinfo,rundata.projinfo,rundata.info,data);
        end
        if rundata.output_eact  
            data=caET;
            data(isnan(data) | rundata.mask<0)=rundata.nodata;
            filename=[rundata.outPath 'GOVar_EAct_' strdate '.' lower(rundata.resultFormat)];
            output_results(filename,rundata.resultFormat,rundata.geoinfo,rundata.projinfo,rundata.info,data);
        end
        if rundata.output_rs    
            data=overland;
            data(isnan(data) | rundata.mask<0)=rundata.nodata;
            filename=[rundata.outPath 'GOVar_RS_' strdate '.' lower(rundata.resultFormat)];
            output_results(filename,rundata.resultFormat,rundata.geoinfo,rundata.projinfo,rundata.info,data);
        end
        if rundata.output_ri    
            data=interflow;
            data(isnan(data) | rundata.mask<0)=rundata.nodata;
            filename=[rundata.outPath 'GOVar_RI_' strdate '.' lower(rundata.resultFormat)];
            output_results(filename,rundata.resultFormat,rundata.geoinfo,rundata.projinfo,rundata.info,data);
        end
        if rundata.output_infil 
            data=infil;
            data(isnan(data) | rundata.mask<0)=rundata.nodata;
            filename=[rundata.outPath 'GOVar_Infil_' strdate '.' lower(rundata.resultFormat)];
            output_results(filename,rundata.resultFormat,rundata.geoinfo,rundata.projinfo,rundata.info,data);
        end
        
        if rundata.output_tot_infil 
            data=tot_infil;
            data(isnan(data) | rundata.mask<0)=rundata.nodata;
            filename=[rundata.outPath 'GOVar_tot_infil_' strdate '.' lower(rundata.resultFormat)];
            output_results(filename,rundata.resultFormat,rundata.geoinfo,rundata.projinfo,rundata.info,data);
        end
        
        if strcmpi(rundata.modelCore,'HydroSlide')==1
            if rundata.output_FS    
                data=fs;
                data(isnan(data) | rundata.mask<0)=rundata.nodata;
                filename=[rundata.outPath 'GOVar_FS_' strdate '.' lower(rundata.resultFormat)];
                output_results(filename,rundata.resultFormat,rundata.geoinfo,rundata.projinfo,rundata.info,data);
            end
        end      

        if (ismember(datestr(step, 'yyyymmddHH'), rundata.observed_time_serise) == 1) && ...
                    (strcmpi(rundata.modelCore,'iHydroSlide3D') == 1)
            if rundata.output_FS3D
                data = CRESTABLE3D_output.FS3D_whole_map;
                data(isnan(data) | LandslideData.mask<0)=rundata.nodata;
                filename=[rundata.outPath 'GOVar_FS3D_' strdate '.' lower(rundata.resultFormat)];
                output_results(filename,rundata.resultFormat,LandslideData.geoinfo,LandslideData.projinfo,rundata.info,data);
            end
            if rundata.output_PF
                data = CRESTABLE3D_output.PF_whole_map;
                data(isnan(data) | LandslideData.mask<0)=rundata.nodata;
                filename=[rundata.outPath 'GOVar_PF_' strdate '.' lower(rundata.resultFormat)];
                output_results(filename,rundata.resultFormat,LandslideData.geoinfo,LandslideData.projinfo,rundata.info,data);
            end
            if rundata.output_Volume
                data = CRESTABLE3D_output.Volume_whole_map;
                data(isnan(data) | LandslideData.mask<0)=rundata.nodata;
                filename=[rundata.outPath 'GOVar_Volume_' strdate '.' lower(rundata.resultFormat)];
                output_results(filename,rundata.resultFormat,LandslideData.geoinfo,LandslideData.projinfo,rundata.info,data);
            end
            if rundata.output_Area
                data = CRESTABLE3D_output.Area_whole_map;
                data(isnan(data) | LandslideData.mask<0)=rundata.nodata;
                filename=[rundata.outPath 'GOVar_Area_' strdate '.' lower(rundata.resultFormat)];
                output_results(filename,rundata.resultFormat,LandslideData.geoinfo,LandslideData.projinfo,rundata.info,data);
            end
        end
        
        
        
        disp(['- Simulating:  ' strdate]);
        cstep = cstep + 1;  
    end % step loop
    
    if strcmpi(rundata.modelCore,'HydroSlide') == 1
        fsmin(rundata.mask~=1)=NaN;
        prediction.fsmin=fsmin;
    end
    %save states
    
    if rundata.saveState
        data=cSM;
        data(isnan(data) | rundata.mask<0)=rundata.nodata;
        filename=[rundata.statePath 'State_' strdate '_W0.' lower(rundata.stateFormat)];
        output_results(filename,rundata.resultFormat,rundata.geoinfo,rundata.projinfo,rundata.info,data);
        
        
        data=overland;
        data(isnan(data) | rundata.mask<0)=rundata.nodata;
        filename=[rundata.statePath 'State_' strdate '_RS.' lower(rundata.stateFormat)];
        output_results(filename,rundata.resultFormat,rundata.geoinfo,rundata.projinfo,rundata.info,data);
        
        data=channel;
        data(isnan(data) | rundata.mask<0)=rundata.nodata;
        filename=[rundata.statePath 'State_' strdate '_RC.' lower(rundata.stateFormat)];
        output_results(filename,rundata.resultFormat,rundata.geoinfo,rundata.projinfo,rundata.info,data);
        
    end
    
end
