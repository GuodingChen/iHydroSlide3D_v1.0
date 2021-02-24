function [prediction] = Main_program(globalCtlFile)
% Developed by Ke Zhang kzhang@hhu.edu.cn %
% updated by Guoding Chen cgdwork@hhu.edu.cn

    % globalCtlFile is ControlFile, named by user and together with the
    % startwindow in your folder
    globalCtl = Variable_initialize(globalCtlFile);
    
    
    %rundata.States.overland_h = 0.01;
    
    rundata.outPath     = globalCtl.resultPath;
    rundata.statePath   = globalCtl.statePath;
    %rundata.nodata      = -9999.0;
    rundata.nodata      = nan;
    
    pth_rst=rundata.outPath;
    if ~exist(pth_rst,'dir')
        mkdir(pth_rst);
    end
    
    %CREST Configuration File
    rundata.modelCore  = globalCtl.modelCore;
    rundata.Hydro_basicFormat = globalCtl.Hydro_basicFormat;
    rundata.Landslide_basicFormat = globalCtl.Landslide_basicFormat;
    rundata.DownscalingBasicFormat = globalCtl.DownscalingBasicFormat;
    rundata.UseMaxLocalNcores = globalCtl.UseMaxLocalNcores;
    rundata.UserDefinedNcores = globalCtl.UserDefinedNcores;
    rundata.Landslide_density = globalCtl.Landslide_density;
    rundata.Divide_tile_number = globalCtl.Divide_tile_number;
    rundata.min_ae = globalCtl.min_ae;
    rundata.max_ae = globalCtl.max_ae;
    rundata.min_be = globalCtl.min_be;
    rundata.max_be = globalCtl.max_be;

    rundata.LandslideOutStart_data = globalCtl.LandslideOutStart_data;
    rundata.LandslideOutEnd_data = globalCtl.LandslideOutEnd_data;
    

    rundata.paramFormat = globalCtl.paramFormat;
    rundata.stateFormat = globalCtl.stateFormat;
    rundata.petFormat   = globalCtl.petFormat;
    rundata.icsFormat   = globalCtl.icsFormat;
    rundata.resultFormat= globalCtl.resultFormat;
    rundata.rainFormat  = globalCtl.rainFormat;
    
    rundata.loadState   = globalCtl.loadState;
    rundata.saveState   = globalCtl.saveState;
    
    
    rundata.output_runoff = globalCtl.output_runoff;
    rundata.output_W      = globalCtl.output_W;
    rundata.output_SM     = globalCtl.output_SM;
%     rundata.output_ExcS   = globalCtl.output_ExcS;
%     rundata.output_ExcI   = globalCtl.output_ExcI;
    rundata.output_FS     = globalCtl.output_FS;
    rundata.output_rain   = globalCtl.output_rain;
    rundata.output_tot_infil   = globalCtl.output_tot_infil;
%     rundata.output_pet    = globalCtl.output_pet;
    rundata.output_epot   = globalCtl.output_epot;
    rundata.output_eact   = globalCtl.output_eact;
    rundata.output_rs     = globalCtl.output_rs;
    rundata.output_ri     = globalCtl.output_ri;
    rundata.output_infil  = globalCtl.output_infil;
    rundata.output_FS3D  = globalCtl.output_FS3D;
    rundata.output_PF  = globalCtl.output_PF;
    rundata.output_Volume  = globalCtl.output_Volume;
    rundata.output_Area  = globalCtl.output_Area;
    %Sanity check
    good=strcmpi(rundata.modelCore,'Hydromodel')==1 || strcmpi(rundata.modelCore,'HydroSlide')...
        || strcmpi(rundata.modelCore,'iHydroSlide3D') == 1;
    if ~good
       error(['The ''model_core'' can only be crest or creslide! However,' ...
           ' you set it to ' rundata.model_core]);
    end
    % Paraller setting check
    if rundata.UseMaxLocalNcores
        rundata.Ncores = feature('numCores');
        if rundata.Ncores < 2
            warning('This is a single-core computer and does not have the conditions for parallel computing')
            rundata.Ncores = 1;
        end
    else
        rundata.Ncores = rundata.UserDefinedNcores;
        if rundata.Ncores > feature('numCores')
            error(['The number of cores set cannot exceed the maximum number of cores of the computer!' ...
           ' Please set the UserDefinedNcores <= ' feature('numCores')]);
        end
    end
    % Set the LandslideModel output time 
    start_date = datenum(rundata.LandslideOutStart_data);
    end_date = datenum(rundata.LandslideOutEnd_data);
    % compute date step: 1 h
    date_step = datenum('2020-06-28 01:00:00') - datenum('2020-06-28 00:00:00');
    dateNum_array = (start_date : date_step : end_date)';
    % prepare a dateStr array
    dateStr_array =  strings(length(dateNum_array),1);
    for i = 1 : length(dateNum_array)
        date_STR = datestr( dateNum_array(i), 'yyyy-mm-dd HH:MM:SS');
        new_date_STR = strcat( date_STR(1:4), date_STR(6:7), date_STR(9:10), date_STR(12:13));
        dateStr_array(i) = new_date_STR;
    end
    rundata.observed_time_serise = dateStr_array;
    
    
    
%     %basicFormat
%     good=strcmpi(rundata.basicFormat,'asc')==1 || ...
%          strcmpi(rundata.basicFormat,'txt')==1 || ...
%          strcmpi(rundata.basicFormat,'tif')==1;
%     if ~good
%         error(['The ''BasicFormat'' can only be asc, txt or tif! However,' ...
%            ' you set it to ' rundata.basicFormat]);
%         
%     end
%     
%     %paramFormat
%     good=strcmpi(rundata.paramFormat,'asc')==1 || ...
%          strcmpi(rundata.paramFormat,'txt')==1 || ...
%          strcmpi(rundata.paramFormat,'tif')==1;
%     if ~good
%         error(['The ''ParamFormat'' can only be asc, txt or tif! However,' ...
%            ' you set it to ' rundata.paramFormat]);
%         
%     end
%     
%     %stateFormat
%     good=strcmpi(rundata.stateFormat,'asc')==1 || ...
%          strcmpi(rundata.stateFormat,'txt')==1 || ...
%          strcmpi(rundata.stateFormat,'tif')==1;
%     if ~good
%         error(['The ''StateFormat'' can only be asc, txt or tif! However,' ...
%            ' you set it to ' rundata.stateFormat]);
%         
%     end
%     
%     %petFormat
%     good=strcmpi(rundata.petFormat,'asc')==1 || ...
%          strcmpi(rundata.petFormat,'txt')==1 || ...
%          strcmpi(rundata.petFormat,'tif')==1;
%     if ~good
%         error(['The ''PETFormat'' can only be asc, txt or tif! However,' ...
%            ' you set it to ' rundata.petFormat]);
%         
%     end
%     
%     %rainFormat
%     good=strcmpi(rundata.rainFormat,'asc')==1 || ...
%          strcmpi(rundata.rainFormat,'txt')==1 || ...
%          strcmpi(rundata.rainFormat,'tif')==1;
%     if ~good
%         error(['The ''RainFormat'' can only be asc, txt or tif! However,' ...
%            ' you set it to ' rundata.rainFormat]);
%         
%     end
%     
%     %resultFormat
%     good=strcmpi(rundata.resultFormat,'asc')==1 || ...
%          strcmpi(rundata.resultFormat,'txt')==1 || ...
%          strcmpi(rundata.resultFormat,'tif')==1;
%     if ~good
%         error(['The ''ResultFormat'' can only be asc, txt or tif! However,' ...
%            ' you set it to ' rundata.resultFormat]);
%         
%     end
%     
%     %icsFormat
%     good=strcmpi(rundata.icsFormat,'asc')==1 || ...
%          strcmpi(rundata.icsFormat,'txt')==1 || ...
%          strcmpi(rundata.icsFormat,'tif')==1;
%     if ~good
%         error(['The ''ICSFormat'' can only be asc, txt or tif! However,' ...
%            ' you set it to ' rundata.icsFormat]);
%         
%     end
    
    
    % Hydro_basicFormat
    good=strcmpi(rundata.Hydro_basicFormat,'asc')==1 || ...
         strcmpi(rundata.Hydro_basicFormat,'txt')==1 || ...
         strcmpi(rundata.Hydro_basicFormat,'tif')==1;
    if ~good
        error(['The ''Hydro_BasicFormat'' can only be asc, txt or tif! However,' ...
           ' you set it to ' rundata.Hydro_basicFormat]);
        
    end
    % Landslide basicFormat
    good=strcmpi(rundata.Landslide_basicFormat,'asc')==1 || ...
         strcmpi(rundata.Landslide_basicFormat,'txt')==1 || ...
         strcmpi(rundata.Landslide_basicFormat,'tif')==1;
    if ~good
        error(['The ''LandslideBasicFormat'' can only be asc, txt or tif! However,' ...
           ' you set it to ' rundata.Landslide_basicFormat]);
        
    end
    % soil downscalingFile format
    good=strcmpi(rundata.DownscalingBasicFormat,'asc')==1 || ...
         strcmpi(rundata.DownscalingBasicFormat,'txt')==1 || ...
         strcmpi(rundata.DownscalingBasicFormat,'tif')==1;
    if ~good
        error(['The ''DownscalingBasicFormat'' can only be asc, txt or tif! However,' ...
           ' you set it to ' rundata.DownscalingBasicFormat]);
        
    end
    % paramFormat
    good=strcmpi(rundata.paramFormat,'asc')==1 || ...
         strcmpi(rundata.paramFormat,'txt')==1 || ...
         strcmpi(rundata.paramFormat,'tif')==1;
    if ~good
        error(['The ''ParamFormat'' can only be asc, txt or tif! However,' ...
           ' you set it to ' rundata.paramFormat]);
        
    end
    
    % stateFormat
    good=strcmpi(rundata.stateFormat,'asc')==1 || ...
         strcmpi(rundata.stateFormat,'txt')==1 || ...
         strcmpi(rundata.stateFormat,'tif')==1;
    if ~good
        error(['The ''StateFormat'' can only be asc, txt or tif! However,' ...
           ' you set it to ' rundata.stateFormat]);
        
    end
    
    % petFormat
    good=strcmpi(rundata.petFormat,'asc')==1 || ...
         strcmpi(rundata.petFormat,'txt')==1 || ...
         strcmpi(rundata.petFormat,'tif')==1;
    if ~good
        error(['The ''PETFormat'' can only be asc, txt or tif! However,' ...
           ' you set it to ' rundata.petFormat]);
        
    end
    
    %rainFormat
    good=strcmpi(rundata.rainFormat,'asc')==1 || ...
         strcmpi(rundata.rainFormat,'txt')==1 || ...
         strcmpi(rundata.rainFormat,'tif')==1;
    if ~good
        error(['The ''RainFormat'' can only be asc, txt or tif! However,' ...
           ' you set it to ' rundata.rainFormat]);
        
    end
    
    %resultFormat
    good=strcmpi(rundata.resultFormat,'asc')==1 || ...
         strcmpi(rundata.resultFormat,'txt')==1 || ...
         strcmpi(rundata.resultFormat,'tif')==1;
    if ~good
        error(['The ''ResultFormat'' can only be asc, txt or tif! However,' ...
           ' you set it to ' rundata.resultFormat]);
        
    end
    
    %icsFormat
    good=strcmpi(rundata.icsFormat,'asc')==1 || ...
         strcmpi(rundata.icsFormat,'txt')==1 || ...
         strcmpi(rundata.icsFormat,'tif')==1;
    if ~good
        error(['The ''ICSFormat'' can only be asc, txt or tif! However,' ...
           ' you set it to ' rundata.icsFormat]);
        
    end
    
    %--------------------------------------------------------------------------
    % read the basic hydrological file
    mainfolder_HydroFile = globalCtl.Hydro_basicPath;
    switch lower(rundata.Hydro_basicFormat)
        case 'txt'
    
            rundata.fdir_path      = [mainfolder_HydroFile, 'FDR.txt'];

            rundata.facc_path      = [mainfolder_HydroFile, 'FAC.txt'];

            rundata.mask_path      = [mainfolder_HydroFile, 'MASK.txt'];

            rundata.dem_path       = [mainfolder_HydroFile, 'DEM.txt'];

            rundata.soil_path      = [mainfolder_HydroFile, 'SOIL.txt'];
 
            rundata.landcover_path = [mainfolder_HydroFile, 'LANDCOVER.txt']; 
        case 'asc'

            rundata.fdir_path      = [mainfolder_HydroFile, 'FDR.asc'];
  
            rundata.facc_path      = [mainfolder_HydroFile, 'FAC.asc'];

            rundata.mask_path      = [mainfolder_HydroFile, 'MASK.asc'];

            rundata.dem_path       = [mainfolder_HydroFile, 'DEM.asc'];

            rundata.soil_path      = [mainfolder_HydroFile, 'SOIL.asc'];

            rundata.landcover_path = [mainfolder_HydroFile, 'LANDCOVER.asc']; 
        case 'tif'

            rundata.fdir_path = [mainfolder_HydroFile, 'FDR.tif'];
  
            rundata.facc_path = [mainfolder_HydroFile, 'FAC.tif'];

            rundata.mask_path = [mainfolder_HydroFile, 'MASK.tif'];

            rundata.dem_path = [mainfolder_HydroFile, 'DEM.tif'];

            rundata.soil_path      = [mainfolder_HydroFile, 'SOIL.tif'];

            rundata.landcover_path = [mainfolder_HydroFile, 'LANDCOVER.tif'];
    end
    
    rundata.mask_path      = [mainfolder_HydroFile, 'MASK.tif'];
    
    % read the landslide basic data
    mainfolder_LandslideFile = globalCtl.Landslide_basicPath;

    switch lower(rundata.Landslide_basicFormat)
        case 'txt'
            
            rundata.mask_fine_path      = [mainfolder_LandslideFile, 'mask_fine.txt'];

            rundata.dem_fine_path       = [mainfolder_LandslideFile, 'DEM_fine.txt'];

            rundata.slope_fine_path      = [mainfolder_LandslideFile, 'slope_fine.txt'];
 
            rundata.aspect_fine_path = [mainfolder_LandslideFile, 'aspect_fine.txt']; 
        case 'asc'

            rundata.mask_fine_path      = [mainfolder_LandslideFile, 'mask_fine.asc'];

            rundata.dem_fine_path       = [mainfolder_LandslideFile, 'DEM_fine.asc'];

            rundata.slope_fine_path      = [mainfolder_LandslideFile, 'slope_fine.asc'];
 
            rundata.aspect_fine_path = [mainfolder_LandslideFile, 'aspect_fine.asc']; 
            
        case 'tif'


            rundata.mask_fine_path      = [mainfolder_LandslideFile, 'mask_fine.tif'];

            rundata.dem_fine_path       = [mainfolder_LandslideFile, 'DEM_fine.tif'];

            rundata.slope_fine_path      = [mainfolder_LandslideFile, 'slope_fine.tif'];
 
            rundata.aspect_fine_path = [mainfolder_LandslideFile, 'aspect_fine.tif']; 
    end
    
    rundata.mask_fine_path      = [mainfolder_LandslideFile, 'mask_fine.tif'];
    % read the soil downscaling file 
     mainfolder_downscaling = globalCtl.DownscalingBasicPath;

    switch lower(rundata.DownscalingBasicFormat)
        case 'txt'
            
            rundata.curvature_coarse_path      = [mainfolder_downscaling, 'curvature_coarse.txt'];

            rundata.curvature_fine_path       = [mainfolder_downscaling, 'curvature_fine.txt'];

            rundata.TWI_coarse_path      = [mainfolder_downscaling, 'TWI_coarse.txt'];
 
            rundata.TWI_fine_path = [mainfolder_downscaling, 'TWI_fine.txt']; 
        case 'asc'

            rundata.curvature_coarse_path      = [mainfolder_downscaling, 'curvature_coarse.asc'];

            rundata.curvature_fine_path       = [mainfolder_downscaling, 'curvature_fine.asc'];

            rundata.TWI_coarse_path      = [mainfolder_downscaling, 'TWI_coarse.asc'];
 
            rundata.TWI_fine_path = [mainfolder_downscaling, 'TWI_fine.asc']; 
            
        case 'tif'


            rundata.curvature_coarse_path      = [mainfolder_downscaling, 'curvature_coarse.tif'];

            rundata.curvature_fine_path       = [mainfolder_downscaling, 'curvature_fine.tif'];

            rundata.TWI_coarse_path      = [mainfolder_downscaling, 'TWI_coarse.tif'];
 
            rundata.TWI_fine_path = [mainfolder_downscaling, 'TWI_fine.tif']; 
            
            rundata.aspect_coarse_path = [mainfolder_downscaling, 'aspect_coarse.tif'];  
    end
    
    

    
    if ~exist(rundata.mask_path,'file')
        error(['The mask file you provided is not in the format of geotiff!' ...
            'You must provide a geotiff mask file!']);
    end
    rundata.projinfo  = geotiffinfo([mainfolder_HydroFile, 'MASK.tif']);
    
    [rundata.mask,rundata.geoinfo]=geotiffread([mainfolder_HydroFile, 'MASK.tif']);

    rundata.slope_path =[];
    
    rundata.info=[];

    switch lower(rundata.Hydro_basicFormat)
        case 'txt'
            [~,geoinfo]=read_ascii_gis([mainfolder_HydroFile, 'FDR.txt']);
            rundata.info=geoinfo;
            rundata.nodata=rundata.info.nodata;
        case 'asc'
            [~,geoinfo]=read_ascii_gis([mainfolder_HydroFile, 'FDR.asc']);
            rundata.info=geoinfo;
            rundata.nodata=rundata.info.nodata;
%         case 'tif'
%             rundata.projinfo = geotiffinfo([mainfolder, 'FDR.tif']);
%             [~,geoinfo]=geotiffread([mainfolder, 'FDR.tif']);
%             rundata.geoinfo=geoinfo;
    end
    
    %Input data (Forcing)
    %Period of rainfall and PET data has to match
    %Rainfall estimates
    %Define naming convention. Example: ST4.2006123100_02083500.tif;
    rundata.pp_path = globalCtl.rainPath;
    switch lower(rundata.rainFormat)
        case 'asc'
            rundata.pp_fname_suffix  = '.asc'; 
        case 'txt'
            rundata.pp_fname_suffix  = '.txt'; 
        case 'tif'
            rundata.pp_fname_suffix  = '.tif';
    end
    
    switch lower(rundata.petFormat)
        case 'asc'
            rundata.pet_fname_suffix = '.asc';
        case 'txt'
            rundata.pet_fname_suffix = '.txt';
        case 'tif'
            rundata.pet_fname_suffix = '.tif';
    end
    % set the time step intercal, in hours
    rundata.pp_tstep = globalCtl.timeStep/3600.0; %hours

    rundata.pet_path = globalCtl.petPath;
    rundata.pet_tstep = rundata.pp_tstep; %hours


    rundata.subfolder_search = 'No';


    rundata.tstep = rundata.pp_tstep; %Time Step in hours
    
    
    rundata.tperiod = datenum(globalCtl.startYear,globalCtl.startMonth, ...
        globalCtl.startDay,globalCtl.startHour,globalCtl.startMinute, ...
        globalCtl.startSecond):rundata.tstep/24: ...
        datenum(globalCtl.endYear,globalCtl.endMonth, ...
        globalCtl.endDay,globalCtl.endHour,globalCtl.endMinute, ...
        globalCtl.endSecond);

    


    if strcmpi(globalCtl.RainFactType,'uniform')
        rundata.Parameters.PRain = globalCtl.RainFact;
    else
        switch  lower(globalCtl.paramFormat)
            case 'asc'
                rundata.Parameters.PRain = [globalCtl.paramPath 'RainFact.asc'];
            case 'txt'
                rundata.Parameters.PRain = [globalCtl.paramPath 'RainFact.txt'];
            case 'tif'
                rundata.Parameters.PRain = [globalCtl.paramPath 'RainFact.tif'];
        end
    end
    

    if strcmpi(globalCtl.RainFactType,'uniform')
        rundata.Parameters.PKE = globalCtl.KE;
    else
        switch  lower(globalCtl.paramFormat)
            case 'asc'
                rundata.Parameters.PKE = [globalCtl.paramPath 'KE.asc'];
            case 'txt'
                rundata.Parameters.PKE = [globalCtl.paramPath 'KE.txt'];
            case 'tif'
                rundata.Parameters.PKE = [globalCtl.paramPath 'KE.tif'];
        end
    end
    
    

    if strcmpi(globalCtl.IMType,'uniform')==1
        rundata.Parameters.PIM = globalCtl.IM;
    else
        switch lower(globalCtl.paramFormat)
            case 'txt'
                rundata.Parameters.PIM = [globalCtl.paramPath, 'IM.txt'];
            case 'asc'
                rundata.Parameters.PIM = [globalCtl.paramPath, 'IM.asc'];
            case 'tif'
                rundata.Parameters.PIM = [globalCtl.paramPath, 'IM.tif'];
        end
    end
    rundata.Parameters.PIM_sc = 1; %Scaling factor
    
    
    

    if strcmpi(globalCtl.WMType,'uniform')==1
        rundata.Parameters.PWM = globalCtl.WM;
    else
        switch lower(globalCtl.paramFormat)
            case 'txt'
                rundata.Parameters.PWM = [globalCtl.paramPath, 'WM.txt'];
            case 'asc'
                rundata.Parameters.PWM = [globalCtl.paramPath, 'WM.asc'];
            case 'tif'
                rundata.Parameters.PWM = [globalCtl.paramPath, 'WM.tif'];
        end
    end
    rundata.Parameters.PWM_sc = 1; %Scaling factor
    

    if strcmpi(globalCtl.KsatType,'uniform')==1
        rundata.Parameters.PFC = globalCtl.Ksat;
    else
        switch lower(globalCtl.paramFormat)
            case 'txt'
                rundata.Parameters.PFC = [globalCtl.paramPath, 'FC.txt'];
            case 'asc'
                rundata.Parameters.PFC = [globalCtl.paramPath, 'FC.asc'];
            case 'tif'
                rundata.Parameters.PFC = [globalCtl.paramPath, 'FC.tif'];
        end
    end
    rundata.Parameters.PFC_sc = 1; %Scaling factor

    if strcmpi(globalCtl.BType,'uniform')==1
        rundata.Parameters.PB = globalCtl.B;
    else
        switch lower(globalCtl.paramFormat)
            case 'txt'
                rundata.Parameters.PB = [globalCtl.paramPath, 'B.txt'];
            case 'asc'
                rundata.Parameters.PB = [globalCtl.paramPath, 'B.asc'];
            case 'tif'
                rundata.Parameters.PB = [globalCtl.paramPath, 'B.tif'];
        end
    end
    rundata.Parameters.PB_sc = 1; %Scaling factor
    

    if strcmpi(globalCtl.coeMType,'uniform')==1
        rundata.Parameters.COEM = globalCtl.coeM;
    else
        switch lower(globalCtl.paramFormat)
            case 'txt'
                rundata.Parameters.COEM = [globalCtl.paramPath, 'COEM.txt'];
            case 'asc'
                rundata.Parameters.COEM = [globalCtl.paramPath, 'COEM.asc'];
            case 'tif'
                rundata.Parameters.COEM = [globalCtl.paramPath, 'COEM.tif'];
        end
    end
    rundata.Parameters.COEM_sc = 1; %Scaling factor
    

    rundata.Parameters.UNDER = rundata.Parameters.PFC;
    rundata.Parameters.UNDER_sc = 1.5; 

  
    rundata.Parameters.TH = 50; %In this case 5 = 113.41 km2 is the minimum for a grid to be considered as channel
    rundata.Parameters.Slope_TH = 0.012;


    %rundata.Parameters.lc=[mainfolder, ''];
    
    %Kinematic Wave rundata.Parameters
    %A = alpha*Q^beta
    %Define how to compute kinematic wave parameters alpha and beta above.
    %Options:
    %Input: The user independently estimates the values and provides gridded
    %data for alpha and beta.
    % rundata.Parameters.Kinematic_Mode = 'Input';
    % rundata.Parameters.BETA = [mainfolder, '02083500_BETA.tif'];
    % rundata.Parameters.BETA_sc = 1; %Scaling factor
    % rundata.Parameters.ALPHA = [mainfolder, '02083500_ALPHA.tif'];
    % rundata.Parameters.ALPHA_sc = 1; %Scaling factor

    %RatingCurve: Uses parameters found by power regression from data at the 
    %basin outlet. Propagates this values upstream using a simple linear 
    %function based on flow accumulation.
    rundata.Parameters.Kinematic_Mode = 'RatingCurve';
    rundata.Parameters.RIVER = 120; %100;
    rundata.Parameters.BETA = 0.92; %0.9; %good 0.9 %0.2; %default 0.2;%0.78;
    rundata.Parameters.ALPHA = 18.5; %good 18%10.0; %dafult 1.0;%1.54;

    % rundata.Parameters.Kinematic_Mode = 'Morphometry_2';
    % rundata.Parameters.RR = [mainfolder, 'Basics/Macon_reliefratio_250m.tif'];
    % rundata.Parameters.RR_sc = 1; %Scaling factor
    % rundata.Parameters.RIVER = 100;

    %Channel - Two sub-options: 
    % rundata.Parameters.Kinematic_Mode = 'ChannelO';
    %RIVER: Inverse of Roughness coefficient. For rivers literature suggests ~50
    % rundata.Parameters.RIVER = 100;
    %CHW: Channel Width. If a scalar value is provided
    % rundata.Parameters.CHW = 50; %Channel Width at the outlet

    % rundata.Parameters.Kinematic_Mode = 'ChannelB';
    %Provide path to channel width data over all channel grids
    % rundata.Parameters.CHW = '/Users/humbertovergara/Documents/Trabajo/HPro_Tool/Code/input_data/Tar_River/HRAP/02083500_RIVERWIDTH.tif';

    %Initial States
    %Can be given arbitrary values
    
    %CREST States
    %Initialize state
    if rundata.loadState
        str_date=datestr(rundata.tperiod(1)-rundata.tstep/24,'yyyymmddHH');
         
        switch lower(globalCtl.stateFormat)
            case 'txt'
                rundata.States.IWU = [globalCtl.statePath, 'State_', str_date, '_', 'W0.txt']; %depth
                %rundata.States.overland_h = [globalCtl.statePath, 'State_', str_date, '_', 'SS0.txt']; %depth
                rundata.States.overland = [globalCtl.statePath, 'State_', str_date, '_', 'RS.txt']; %cms
                rundata.States.interflow = [globalCtl.statePath, 'State_', str_date, '_', 'SI0.txt']; %depth
                rundata.States.channel = [globalCtl.statePath, 'State_', str_date, '_', 'RC.txt']; %cms
            case 'asc'
                rundata.States.IWU = [globalCtl.statePath, 'State_', str_date, '_', 'W0.asc'];
                %rundata.States.overland_h = [globalCtl.statePath, 'State_', str_date, '_', 'SS0.asc'];
                rundata.States.overland = [globalCtl.statePath, 'State_', str_date, '_', 'RS.asc'];
                rundata.States.interflow = [globalCtl.statePath, 'State_', str_date, '_', 'SI0.asc'];
                rundata.States.channel = [globalCtl.statePath, 'State_', str_date, '_', 'RC.asc'];
            case 'tif'
                rundata.States.IWU = [globalCtl.statePath, 'State_', str_date, '_', 'W0.tif'];
                %rundata.States.overland_h = [globalCtl.statePath, 'State_', str_date, '_', 'SS0.tif'];
                rundata.States.overland = [globalCtl.statePath, 'State_', str_date, '_', 'RS.tif'];
                rundata.States.interflow = [globalCtl.statePath, 'State_', str_date, '_', 'SI0.tif'];
                rundata.States.channel = [globalCtl.statePath, 'State_', str_date, '_', 'RC.tif'];
        end
        rundata.stateFormat=globalCtl.stateFormat;
    else
        if strcmpi(globalCtl.W0Type,'uniform')==1
            rundata.States.IWU = globalCtl.W0;
        else
            switch lower(globalCtl.icsFormat)
                case 'txt'
                    rundata.States.IWU = [globalCtl.icsPath, 'W0.txt'];
                case 'asc'
                    rundata.States.IWU = [globalCtl.icsPath, 'W0.asc'];
                case 'tif'
                    rundata.States.IWU = [globalCtl.icsPath, 'W0.tif'];
            end
        end

%         %Routing States
%         if strcmpi(globalCtl.SS0Type,'uniform')==1
%             rundata.States.overland_h = globalCtl.SS0;
%         else
%             switch lower(globalCtl.icsFormat)
%                 case 'txt'
%                     rundata.States.overland_h = [globalCtl.icsPath, 'SS0.txt'];
%                 case 'asc'
%                     rundata.States.overland_h = [globalCtl.icsPath, 'SS0.asc'];
%                 case 'tif'
%                     rundata.States.overland_h = [globalCtl.icsPath, 'SS0.tif'];
%             end
%         end

        if strcmpi(globalCtl.RSType,'uniform')==1
            rundata.States.overland = globalCtl.RS;
        else
            switch lower(globalCtl.icsFormat)
                case 'txt'
                    rundata.States.overland = [globalCtl.icsPath, 'RS.txt'];
                case 'asc'
                    rundata.States.overland = [globalCtl.icsPath, 'RS.asc'];
                case 'tif'
                    rundata.States.overland = [globalCtl.icsPath, 'RS.tif'];
            end
        end

        if strcmpi(globalCtl.SI0Type,'uniform')==1
            rundata.States.interflow = globalCtl.SI0;
        else
            switch lower(globalCtl.icsFormat)
                case 'txt'
                    rundata.States.interflow = [globalCtl.icsPath, 'SI0.txt'];
                case 'asc'
                    rundata.States.interflow = [globalCtl.icsPath, 'SI0.asc'];
                case 'tif'
                    rundata.States.interflow = [globalCtl.icsPath, 'SI0.tif'];
            end
        end


        if strcmpi(globalCtl.RType,'uniform')==1
            rundata.States.channel = globalCtl.R;
        else
            switch lower(globalCtl.icsFormat)
                case 'txt'
                    rundata.States.channel = [globalCtl.icsPath, 'R.txt'];
                case 'asc'
                    rundata.States.channel = [globalCtl.icsPath, 'R.asc'];
                case 'tif'
                    rundata.States.channel = [globalCtl.icsPath, 'R.tif'];
            end
        end
        
        rundata.stateFormat=globalCtl.icsFormat;
    end

    %Define lat/long coordinates of additional locations (optional)
    %"crestmodeling" always return the basin outlet as the point with highest
    %flow accumulation.
    
    if globalCtl.NOutPixs>0
        rundata.extra_locations_latlon=zeros(globalCtl.NOutPixs,2);
        for i=1:globalCtl.NOutPixs
            rundata.extra_locations_latlon(i,1)=globalCtl.OutPixLats(i);
            rundata.extra_locations_latlon(i,2)=globalCtl.OutPixLons(i);
        end
    else
        rundata.extra_locations_latlon=[];
    end
    
    rundata.extra_locations_name=globalCtl.OutPixNames;
    rundata.extra_locations_area=globalCtl.OutPixAreas;
    %rundata.extra_locations_latlon(1,:) = %[35.415,-83.568]; %02081500
    %rundata.extra_locations_latlon(2,:) = [35.3364,-83.5269]; %03503000


    %Drainage area can be provided for more accurate location
    %rundata.extra_locations_area(1) = 1722.5; %02081500
    %rundata.extra_locations_area(2) = 1129.2; %03503000

    %Get a mask from these extra locations?: Yes = 1, No = 0
    rundata.extra_locations_mask = 1;

    %Initiate model run
    rundata.run_mode = 'memory_efficient';
    % rundata.load_basicgrids = 'saved_basicgrids.mat';
    % rundata.save_basicgrids = 'saved_basicgrids.mat';
    %[prediction] = HPROLite(rundata);
    Ncores = rundata.Ncores;
    Parallel_environment = parpool(Ncores);
    [prediction] = Code_execute(rundata);

%    %Save results
%     for i=1:length(rundata.tperiod)
%         filename=[pth_rst '/FS_' datestr(rundata.tperiod(i),'yyyymmddHH') ...
%             '.tif'];
%         data_out=squeeze(prediction.fs(:,:,i));
%         geotiffwrite(filename,single(data_out),geoinfo, 'GeoKeyDirectoryTag'...
%             , rundata.projinfo.GeoTIFFTags.GeoKeyDirectoryTag...
%             , 'TiffTags', struct('Compression','Deflate'));
% 
%         filename=[pth_rst '/RI_' datestr(rundata.tperiod(i),'yyyymmddHH') ...
%             '.tif'];
%         data_out=squeeze(prediction.infil(:,:,i));
%         geotiffwrite(filename,single(data_out),geoinfo, 'GeoKeyDirectoryTag'...
%             , rundata.projinfo.GeoTIFFTags.GeoKeyDirectoryTag...
%             , 'TiffTags', struct('Compression','Deflate'));
%     end


%Save results
filename=[globalCtl.resultPath 'Outlet_Results.csv'];
fid=fopen(filename,'w');
tstep = rundata.tstep; %time step in hrs
tperiod = rundata.tperiod; %time period for simulation


cstep=1;
for step = tperiod(1):(tstep/24):tperiod(end)
    switch lower(rundata.modelCore)
        case 'crest'
            if cstep==1
                fprintf(fid,'DateTime,Rain,EPot,EAct,W,SM,Infil,R\n');
            end
            fprintf(fid,'%s,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f\n', ...
                datestr(step,'yyyy-mm-dd HH:MM:SS'),prediction.Rain(1,cstep), ...
                prediction.EPot(1,cstep),prediction.EAct(1,cstep), ...
                prediction.W(1,cstep),prediction.SM(1,cstep), ...
                prediction.Infil(1,cstep),prediction.R(1,cstep));
        case 'creslide'
            if cstep==1
                fprintf(fid,'DateTime,Rain,EPot,EAct,W,SM,Infil,FS,R, tot_infil\n');
            end
            fprintf(fid,'%s,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f\n', ...
                datestr(step,'yyyy-mm-dd HH:MM:SS'),prediction.Rain(1,cstep), ...
                prediction.EPot(1,cstep),prediction.EAct(1,cstep), ...
                prediction.W(1,cstep),prediction.SM(1,cstep), prediction.Infil(1,cstep),...
                prediction.FS(1,cstep),prediction.R(1,cstep), prediction.Tot_infil(1,cstep));
        case 'crestable3d'
            if cstep==1
                fprintf(fid,'DateTime,Rain,EPot,EAct,W,SM,Infil,R, tot_infil\n');
            end
            fprintf(fid,'%s,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f\n', ...
                datestr(step,'yyyy-mm-dd HH:MM:SS'),prediction.Rain(1,cstep), ...
                prediction.EPot(1,cstep),prediction.EAct(1,cstep), ...
                prediction.W(1,cstep),prediction.SM(1,cstep), prediction.Infil(1,cstep),...
                prediction.R(1,cstep), prediction.Tot_infil(1,cstep));
    end
    cstep=cstep+1;
end
fclose(fid);

prediction.Robs=prediction.R*NaN;
prediction.statistics.NSCE=size(prediction.Rain,1)*NaN;
prediction.statistics.CC=size(prediction.Rain,1)*NaN;
prediction.statistics.Bias=size(prediction.Rain,1)*NaN;
for id=1:size(prediction.Rain,1)-1
    %Read observations
    filename=[globalCtl.obsPath char(globalCtl.OutPixNames{id}) '_Obs.csv'];
    
    if exist(filename,'file')
        fid=fopen(filename,'r');
        data=textscan(fid,'%s %f','delimiter',',','headerlines',1);
        fclose(fid);
        data{1}=str2double(data{1});
        cstep=1;
        for step = tperiod(1):(tstep/24):tperiod(end)
            time=str2double(datestr(step,'yyyymmddHH'));
            ind = find(data{1} == time);
            if ind>0
                prediction.Robs(id+1,cstep)=data{2}(ind);
            end
            cstep=cstep+1;
        end
    end
    
    %Derive statistics
    ys=squeeze(prediction.Robs(id+1,:));
    xs=squeeze(prediction.R(id+1,:));
    
    
    xs(isnan(ys))=NaN;
    xs=xs(~isnan(xs));
   
    ys=ys(~isnan(ys));
%     xs = xs(418:554);
%     ys = ys(418:554);
    
    bias=mean(xs-ys)/mean(ys)*100;
    [r,~]=corrcoef(ys,xs);
    
    y_=mean(ys);
    ns=1-sum((ys-xs).*(ys-xs))/sum((ys-y_).*(ys-y_));
    
    filename=[globalCtl.resultPath 'Outpix_' char(globalCtl.OutPixNames{id}) '_Results_Statistics.csv'];
    fid=fopen(filename,'w');
    fprintf(fid,'%s,%f\n','NSCE',ns);
    fprintf(fid,'%s,%f\n','Bias(%)',bias);
    fprintf(fid,'%s,%f\n','CC',r(1,2));
    fclose(fid);
    
    prediction.statistics.NSCE(id)=ns;
    prediction.statistics.Bias(id)=bias;
    prediction.statistics.CC(id)=r(1,2);
    
    
    tmp=squeeze(prediction.Robs(id+1,:));
    tmp(isnan(tmp))=rundata.nodata;
    filename=[globalCtl.resultPath 'Outpix_' char(globalCtl.OutPixNames{id}) '_Results.csv'];
    fid=fopen(filename,'w');
    cstep=1;
    for step = tperiod(1):(tstep/24):tperiod(end)
        switch lower(rundata.modelCore)
            case 'crest'
                if cstep==1
                    fprintf(fid,'DateTime,Rain,EPot,EAct,W,SM,Infil,R,RObs\n');
                end
                fprintf(fid,'%s,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f\n', ...
                    datestr(step,'yyyy-mm-dd HH:MM:SS'),prediction.Rain(id+1,cstep), ...
                    prediction.EPot(id+1,cstep),prediction.EAct(id+1,cstep), ...
                    prediction.W(id+1,cstep),prediction.SM(id+1,cstep), prediction.Infil(id+1,cstep),...
                    prediction.R(id+1,cstep),tmp(cstep));
                cstep=cstep+1;
            case 'creslide'
                if cstep==1
                    fprintf(fid,'DateTime,Rain,EPot,EAct,W,SM,Infil,FS,R,RObs, tot_infil\n');
                end
                fprintf(fid,'%s,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f\n', ...
                    datestr(step,'yyyy-mm-dd HH:MM:SS'),prediction.Rain(id+1,cstep), ...
                    prediction.EPot(id+1,cstep),prediction.EAct(id+1,cstep), ...
                    prediction.W(id+1,cstep),prediction.SM(id+1,cstep),prediction.Infil(id+1,cstep), ...
                    prediction.FS(id+1,cstep),prediction.R(id+1,cstep),tmp(cstep), prediction.Tot_infil(id+1,cstep));
                cstep=cstep+1;
           case 'crestable3d'
                if cstep==1
                    fprintf(fid,'DateTime,Rain,EPot,EAct,W,SM,Infil,R,RObs, tot_infil\n');
                end
                fprintf(fid,'%s,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f\n', ...
                    datestr(step,'yyyy-mm-dd HH:MM:SS'),prediction.Rain(id+1,cstep), ...
                    prediction.EPot(id+1,cstep),prediction.EAct(id+1,cstep), ...
                    prediction.W(id+1,cstep),prediction.SM(id+1,cstep),prediction.Infil(id+1,cstep), ...
                    prediction.R(id+1,cstep),tmp(cstep), prediction.Tot_infil(id+1,cstep));
                cstep=cstep+1;
        end
    end
    fclose(fid);
end

    
         
                

    % figure(2);
    % clf;
    % extra.obsq = dlmread('/Users/goodganker/Zhang/OU/Codes/CRESLIDE_MATLAB/Example/Macon/OBS/03503000_Obs.csv', ',', 1,1);
    % extra.period = datenum('01-Jan-2003 00:00:00'):1/24:datenum('31-Dec-2007 23:00:00');
    % idxi_2 = find(extra.period == rundata.tperiod(1));
    % idxf_2 = find(extra.period == rundata.tperiod(end));
    % [bias,rmse,nsce] = nanhydrostat(extra.obsq(idxi_2:idxf_2),prediction.Q(2,:)');
    % rankcc = corr(extra.obsq(idxi_2:idxf_2),prediction.Q(2,:)', 'Type', 'Spearman');
    % cc = corrcoef(extra.obsq(idxi_2:idxf_2),prediction.Q(2,:));
    % plot(extra.period(idxi_2:idxf_2),extra.obsq(idxi_2:idxf_2), 'DisplayName', 'Observed', 'Color', 'k', 'LineWidth', 2);
    % hold all;
     %plot(rundata.tperiod,prediction.R(2,:), 'DisplayName', 'CREST', 'Color', 'r', 'LineWidth', 2);
    % text(rundata.tperiod(30),100, ['Bias (%) = ', num2str(bias, '%.2f')], 'FontSize', 14);
    % text(rundata.tperiod(30),90, ['RMSE (%) = ', num2str(rmse, '%.2f')], 'FontSize', 14);
    % text(rundata.tperiod(30),80, ['NSCE = ', num2str(nsce, '%.2f')], 'FontSize', 14);
    % text(rundata.tperiod(30),70, ['CC = ', num2str(cc(2), '%.2f')], 'FontSize', 14);
    % text(rundata.tperiod(30),60, ['Rank CC = ', num2str(rankcc(1), '%.2f')], 'FontSize', 14);
    % datetick('x','mmm');
    % set(gca, 'FontSize', 12); xlabel('Time (hours)', 'FontSize', 14); ylabel('Streamflow (cms)', 'FontSize', 14); legend('show', 'Location', 'NorthWest');
    % saveas(gcf, ['Test_Simulation_HPRO_', rundata.model_core, '_Inside.tif']);
    %close;
    % 
    % % %%Save Result
    % % fid=fopen('./Prediction.csv','w');
    % % fprintf(fid,'Time,Q(m^3/s)\n');
    % % for i=1:length(rundata.tperiod)
    % %    fprintf(fid,'%s, %8.2f\n',datestr(rundata.tperiod(i),'yyyy-mm-dd HH:MM:SS'), prediction.Q(2,i)); 
    % % end
    % % fclose(fid);
    delete(Parallel_environment)
    disp('------ model execution ends ------')
end
