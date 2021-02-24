function [dem,slope,facc,fdir,gridres,outpar,outstates,outlet,...
    LandslideData, DownscalingData] = Hydro_BasicData_Read(rundata)
%Function that reads basic data for Hydromodel
%Version 1.0 - February, 2014
%Written by Humberto Vergara - humber@ou.edu
%Hydrometeorology and Remote-Sensing (HyDROS) Laboratory - http://hydro.ou.edu
%Advanced Radar Research Center (ARRC) - http://arrc.ou.edu
%The University of Oklahoma - http://ou.edu
%
%Arguments
%rundata - Struct variable with all simulation info as defined in the
%configuration file

warning('off');

fprintf('Reading Grids: ');
% hydrological data need
Parameters = rundata.Parameters;
States     = rundata.States;
flowdir    = rundata.fdir_path;
flowacc    = rundata.facc_path;
mask       = rundata.mask_path;
dem_p      = rundata.dem_path;
mapinfo    = rundata.projinfo;

% landslideModel data need
mask_fine = rundata.mask_fine_path;
dem_fine = rundata.dem_fine_path;
aspect_fine = rundata.aspect_fine_path;
slope_fine = rundata.slope_fine_path;






% downscaling data need
curvature_coarse = rundata.curvature_coarse_path;
curvature_fine = rundata.curvature_fine_path;
TWI_coarse = rundata.TWI_coarse_path;
TWI_fine = rundata.TWI_fine_path;
aspect_coarse = rundata.aspect_coarse_path;



%% Basic files
%Read basic files and get geographic information
[mk_Z, ~] = geotiffread(mask); 

switch lower(rundata.Hydro_basicFormat)
    case {'txt','asc'}
        [fd_Z, ~] = read_ascii_gis(flowdir);
        [fa_Z, ~] = read_ascii_gis(flowacc);
        [sl_Z,~]  = read_ascii_gis(dem_p);
    case 'tif'
        [fd_Z, ~] = geotiffread(flowdir);
        [fa_Z, ~] = geotiffread(flowacc);
        [sl_Z,~]  = geotiffread(dem_p); 
end

% read the landslide need data
switch lower(rundata.Landslide_basicFormat)
    case {'txt','asc'}
        [LandslideData.mask, ~] = read_ascii_gis(mask_fine);
        [LandslideData.dem, ~] = read_ascii_gis(dem_fine);
        [LandslideData.aspect,~]  = read_ascii_gis(aspect_fine);
        [LandslideData.slope, ~] = read_ascii_gis(slope_fine);
    case 'tif'
        LandslideData.projinfo  = geotiffinfo(mask_fine);
        [LandslideData.mask, LandslideData.geoinfo] = geotiffread(mask_fine);
        [LandslideData.dem, ~] = geotiffread(dem_fine);
        [LandslideData.aspect, ~] = geotiffread(aspect_fine);
        [LandslideData.slope, ~] = geotiffread(slope_fine);
end

switch lower(rundata.DownscalingBasicFormat)
    case {'txt','asc'}
        [DownscalingData.curvature_coarse, ~]=read_ascii_gis(curvature_coarse);
        [DownscalingData.curvature_fine, ~] = read_ascii_gis(curvature_fine);
        [DownscalingData.TWI_coarse, ~] = read_ascii_gis(TWI_coarse);
        [DownscalingData.TWI_fine, ~] = read_ascii_gis(TWI_fine);
    case 'tif'
        [DownscalingData.curvature_coarse, DownscalingData.geoinfo]=geotiffread(curvature_coarse);
        [DownscalingData.curvature_fine, ~] = geotiffread(curvature_fine);
        [DownscalingData.TWI_coarse, ~] = geotiffread(TWI_coarse);
        [DownscalingData.TWI_fine, ~] = geotiffread(TWI_fine);
        [DownscalingData.aspect_90, ~] = geotiffread(aspect_coarse);
end
    LandslideData.Landslide_density = rundata.Landslide_density;
    LandslideData.Divide_tile_number = rundata.Divide_tile_number;
    LandslideData.min_ae = rundata.min_ae;
    LandslideData.max_ae = rundata.max_ae;
    LandslideData.min_be = rundata.min_be;
    LandslideData.max_be = rundata.max_be;






fd_Z = double(fd_Z); 
fd_Z(fd_Z<0) = NaN;
fa_Z = double(fa_Z); 
fa_Z(fa_Z<0) = NaN;
mk_Z = double(mk_Z); 
mk_Z(mk_Z<0) = NaN; 
mk_Z(isnan(mk_Z)==0) = 0;
dem = double(sl_Z); dem(dem<0) = NaN;
dem = dem+mk_Z;


if (isempty(rundata.slope_path) == 0)
    % in this way, the slope data has been input
    switch lower(rundata.Hydro_basicFormat)
        case {'txt','asc'}
            [sl_Z,~] = read_ascii_gis(rundata.slope_path); 
        case 'tif'
            [sl_Z,~] = geotiffread(rundata.slope_path); 
    end
    
    slope = double(sl_Z); 
    slope(slope<0) = 0.0001;
    slope = slope + mk_Z;
else
    % have no input slope map, the slope is computed by DEM
    slope = [];
end

%Mask to basin
fdir = fd_Z+mk_Z;
facc = fa_Z+mk_Z;


%% Read Parameters
% CREST Parameters
if (ischar(Parameters.PKE) == 1)
    switch lower(rundata.paramFormat)
        case {'txt','asc'}
            [cpar_Z,~] = read_ascii_gis(Parameters.PKE); 
        case 'tif'
            [cpar_Z,~] = geotiffread(Parameters.PKE); 
    end
    cpar_Z = double(cpar_Z); cpar_Z(cpar_Z<0) = NaN;
    outpar.PKE = (cpar_Z+mk_Z).*Parameters.PKE_sc;
else
    outpar.PKE = Parameters.PKE+mk_Z;
end

if (ischar(Parameters.PRain) == 1)
    switch lower(rundata.paramFormat)
        case {'txt','asc'}
            [cpar_Z,~] = read_ascii_gis(Parameters.PRain); 
        case 'tif'
            [cpar_Z,~] = geotiffread(Parameters.PRain); 
    end
    cpar_Z = double(cpar_Z); cpar_Z(cpar_Z<0) = NaN;
    outpar.PRain = (cpar_Z+mk_Z);
else
    outpar.PRain = Parameters.PRain+mk_Z;
end


if (ischar(Parameters.PWM) == 1)
    switch lower(rundata.paramFormat)
        case {'txt','asc'}
            [cpar_Z,~] = read_ascii_gis(Parameters.PWM); 
        case 'tif'
            [cpar_Z,~] = geotiffread(Parameters.PWM); 
    end   
    cpar_Z = double(cpar_Z); cpar_Z(cpar_Z<0) = NaN;
    outpar.PWM = (cpar_Z+mk_Z).*Parameters.PWM_sc;
else
    outpar.PWM = Parameters.PWM+mk_Z;
end


% if (ischar(Parameters.LEAKO) == 1)
%     [cpar_Z,~] = geotiffread(Parameters.LEAKO);
%     cpar_Z = double(cpar_Z); cpar_Z(cpar_Z<0) = NaN;
%     outpar.LEAKO = (cpar_Z+mk_Z).*Parameters.LEAKO_sc;
% else
%     outpar.LEAKO = Parameters.LEAKO+mk_Z;
% end

% if (ischar(Parameters.LEAKI) == 1)
%     [cpar_Z,~] = geotiffread(Parameters.LEAKI);
%     cpar_Z = double(cpar_Z); cpar_Z(cpar_Z<0) = NaN;
%     outpar.LEAKI = (cpar_Z+mk_Z).*Parameters.LEAKI_sc;
% else
%     outpar.LEAKI = Parameters.LEAKI+mk_Z;
% end

if (ischar(Parameters.PB) == 1)
    switch lower(rundata.paramFormat)
        case {'txt','asc'}
            [cpar_Z,~] = read_ascii_gis(Parameters.PB); 
        case 'tif'
            [cpar_Z,~] = geotiffread(Parameters.PB);
    end  
    
    cpar_Z = double(cpar_Z); cpar_Z(cpar_Z<0) = 0;
    outpar.PB = (cpar_Z+mk_Z).*Parameters.PB_sc;
else
    outpar.PB = Parameters.PB+mk_Z;
end

if (ischar(Parameters.PFC) == 1)
    switch lower(rundata.paramFormat)
        case {'txt','asc'}
            [cpar_Z,~] = read_ascii_gis(Parameters.PFC); 
        case 'tif'
            [cpar_Z,~] = geotiffread(Parameters.PFC);
    end  
    
    cpar_Z = double(cpar_Z); cpar_Z(cpar_Z<0) = NaN;
    outpar.PFC = (cpar_Z+mk_Z).*Parameters.PFC_sc;
else
    outpar.PFC = Parameters.PFC+mk_Z;
end

if (ischar(Parameters.PIM) == 1)
    switch lower(rundata.paramFormat)
        case {'txt','asc'}
            [cpar_Z,~] = read_ascii_gis(Parameters.PIM); 
        case 'tif'
            [cpar_Z,~] = geotiffread(Parameters.PIM);
    end
    
    cpar_Z = double(cpar_Z); cpar_Z(cpar_Z<0) = NaN;
    outpar.PIM = (cpar_Z+mk_Z).*Parameters.PIM_sc;
else
    outpar.PIM = Parameters.PIM+mk_Z;
end
        

        
if (ischar(Parameters.COEM) == 1)
    switch lower(rundata.paramFormat)
        case {'txt','asc'}
            [cpar_Z,~] = read_ascii_gis(Parameters.COEM); 
        case 'tif'
            [cpar_Z,~] = geotiffread(Parameters.COEM);
    end
    
    cpar_Z = double(cpar_Z); cpar_Z(cpar_Z<0) = NaN;
    outpar.COEM = (cpar_Z+mk_Z).*Parameters.COEM_sc;
else
    outpar.COEM = Parameters.COEM+mk_Z;
end

if (ischar(Parameters.UNDER) == 1)
    switch lower(rundata.paramFormat)
        case {'txt','asc'}
            [cpar_Z,~] = read_ascii_gis(Parameters.UNDER); 
        case 'tif'
            [cpar_Z,~] = geotiffread(Parameters.UNDER);
    end
    
    cpar_Z = double(cpar_Z); cpar_Z(cpar_Z<0) = NaN;
    outpar.UNDER = (cpar_Z+mk_Z).*Parameters.UNDER_sc;
else
    outpar.UNDER = Parameters.UNDER+mk_Z;
end


switch lower(rundata.modelCore)
    case 'hydroslide'
        switch lower(rundata.Hydro_basicFormat)
            case {'txt','asc'}
                [cpar_Z,~] = read_ascii_gis(rundata.soil_path); 
            case 'tif'
                [cpar_Z,~] = geotiffread(rundata.soil_path);
        end 
        
        cpar_Z = double(cpar_Z); cpar_Z(cpar_Z<0) = 0;

        outpar.soil=cpar_Z+mk_Z;

        switch lower(rundata.Hydro_basicFormat)
            case {'txt','asc'}
                [cpar_Z,~] = read_ascii_gis(rundata.landcover_path); 
            case 'tif'
                [cpar_Z,~] = geotiffread(rundata.landcover_path);
        end 
        cpar_Z = double(cpar_Z); cpar_Z(cpar_Z<0) = 0;
        outpar.landcover=cpar_Z+mk_Z;
    case 'ihydroslide3d'
    switch lower(rundata.Hydro_basicFormat)
        case {'txt','asc'}
            [cpar_Z,~] = read_ascii_gis(rundata.soil_path); 
        case 'tif'
            [cpar_Z,~] = geotiffread(rundata.soil_path);
    end 

    cpar_Z = double(cpar_Z); cpar_Z(cpar_Z<0) = 0;

    outpar.soil=cpar_Z+mk_Z;

    switch lower(rundata.Hydro_basicFormat)
        case {'txt','asc'}
            [cpar_Z,~] = read_ascii_gis(rundata.landcover_path); 
        case 'tif'
            [cpar_Z,~] = geotiffread(rundata.landcover_path);
    end 
    cpar_Z = double(cpar_Z); cpar_Z(cpar_Z<0) = 0;
    outpar.landcover=cpar_Z+mk_Z;
    otherwise %do nothing
end
            
            






%Lumped parameters
outpar.TH = Parameters.TH;
outpar.Slope_TH=Parameters.Slope_TH;

%% Kinematic Routing Parameters
outpar.Kinematic_Mode = Parameters.Kinematic_Mode;
switch Parameters.Kinematic_Mode
    case 'Input'
        %ALPHA parameter
        [cpar_Z,~] = geotiffread(Parameters.ALPHA);
        cpar_Z = double(cpar_Z); cpar_Z(cpar_Z<0) = NaN;
        outpar.ALPHA = (cpar_Z+mk_Z).*Parameters.ALPHA_sc;
        
        %BETA parameter
        [cpar_Z,~] = geotiffread(Parameters.BETA);
        cpar_Z = double(cpar_Z); cpar_Z(cpar_Z<0) = NaN;
        outpar.BETA = (cpar_Z+mk_Z).*Parameters.BETA_sc;
    case 'RatingCurve'
        outpar.ALPHA = Parameters.ALPHA; 
        outpar.BETA = Parameters.BETA;
    case 'ChannelO'
        outpar.CHW = Parameters.CHW;
        %Lumped parameters
        outpar.RIVER = Parameters.RIVER;
    case 'ChannelB'
        [cpar_Z,~] = geotiffread(Parameters.CHW);
        cpar_Z = double(cpar_Z); cpar_Z(cpar_Z<0) = NaN;
        outpar.CHW = (cpar_Z+mk_Z).*Parameters.CHW_sc;
        %Lumped parameters
        outpar.RIVER = Parameters.RIVER;
end

%% Read States
if (ischar(States.IWU) == 1)
    switch lower(rundata.stateFormat)
        case {'txt','asc'}
            [cpar_Z,~] = read_ascii_gis(States.IWU); 
        case 'tif'
            [cpar_Z,~] = geotiffread(States.IWU);
    end 
    
    cpar_Z = double(cpar_Z); cpar_Z(cpar_Z<0) = NaN;
    outstates.IWU = cpar_Z+mk_Z;
else
    outstates.IWU = States.IWU+mk_Z;
end

% if (ischar(States.ISO) == 1)
%     [cpar_Z,~] = geotiffread(States.ISO);
%     cpar_Z = double(cpar_Z); cpar_Z(cpar_Z<0) = NaN;
%     outstates.ISO = cpar_Z+mk_Z;
% else
%     outstates.ISO = States.ISO+mk_Z;
% end

% if (ischar(States.ISU) == 1)
%     [cpar_Z,~] = geotiffread(States.ISU);
%     cpar_Z = double(cpar_Z); cpar_Z(cpar_Z<0) = NaN;
%     outstates.ISU = cpar_Z+mk_Z;
% else
%     outstates.ISU = States.ISU+mk_Z;
% end    



if (ischar(States.overland) == 1)
    switch lower(rundata.stateFormat)
        case {'txt','asc'}
            [cpar_Z,~] = read_ascii_gis(States.overland); 
        case 'tif'
            [cpar_Z,~] = geotiffread(States.overland);
    end
    
    cpar_Z = double(cpar_Z); cpar_Z(cpar_Z<0) = NaN;
    outstates.overland = cpar_Z+mk_Z;
else
    outstates.overland = States.overland+mk_Z;
end

% if (ischar(States.overland_h) == 1)
%     switch lower(rundata.stateFormat)
%         case {'txt','asc'}
%             [cpar_Z,~] = read_ascii_gis(States.overland_h); 
%         case 'tif'
%             [cpar_Z,~] = geotiffread(States.overland_h);
%     end
%     
%     cpar_Z = double(cpar_Z); cpar_Z(cpar_Z<0) = NaN;
%     outstates.overland_h = cpar_Z+mk_Z;
% else
%     outstates.overland_h = States.overland_h+mk_Z;
% end

if (ischar(States.channel) == 1)
    switch lower(rundata.stateFormat)
        case {'txt','asc'}
            [cpar_Z,~] = read_ascii_gis(States.channel); 
        case 'tif'
            [cpar_Z,~] = geotiffread(States.channel);
    end
    
    cpar_Z = double(cpar_Z); cpar_Z(cpar_Z<0) = NaN;
    outstates.channel = cpar_Z+mk_Z;
else
    outstates.channel = States.channel+mk_Z;
end

if (ischar(States.interflow) == 1)
    switch lower(rundata.stateFormat)
        case {'txt','asc'}
            [cpar_Z,~] = read_ascii_gis(States.interflow); 
        case 'tif'
            [cpar_Z,~] = geotiffread(States.interflow);
    end
    
    cpar_Z = double(cpar_Z); cpar_Z(cpar_Z<0) = NaN;
    outstates.interflow = cpar_Z+mk_Z;
else
    outstates.interflow = States.interflow+mk_Z;
end
    
%% Return grid size in meters
switch mapinfo.SpatialRef.CoordinateSystemType
    case 'geographic'
        gridres.DeltaX = deg2km(mapinfo.SpatialRef.DeltaLon).*1000;
        gridres.DeltaY = deg2km(mapinfo.SpatialRef.DeltaLat).*1000;
    case 'planar'
        gridres.DeltaX = mapinfo.SpatialRef.DeltaX;
        gridres.DeltaY = mapinfo.SpatialRef.DeltaY;
end

%Identify basin outlet based on flow accumulation
[outlet.y,outlet.x] = find(facc == nanmax(facc(:)));

fprintf('Done.\n');

%Identify extra locations
if (isfield(rundata,'extra_locations_latlon') == 1)
    n_locs = size(rundata.extra_locations_latlon,1);
    fprintf('Estimating pixel coordinates of %g extra locations: ', n_locs);
    for locs = 1:n_locs
        given_outlet.latitude = rundata.extra_locations_latlon(locs,1);
        given_outlet.longitude = rundata.extra_locations_latlon(locs,2);
        tol = 20;
        max_r = 1000;
        if (isfield(rundata,'extra_locations_area') == 1)
            basin_area = rundata.extra_locations_area(locs);
            [outlet.masks{locs},outlet.extra.y(locs),outlet.extra.x(locs),~,~,~,~] = Extract_basin(facc,fdir,given_outlet,basin_area,mapinfo,rundata.extra_locations_mask,tol,max_r);
        else
            [outlet.masks{locs},outlet.extra.y(locs),outlet.extra.x(locs),~,~,~,~] = Extract_basin(facc,fdir,given_outlet,[],mapinfo,rundata.extra_locations_mask,tol,max_r);
        end
    end
    fprintf('Done.\n');
else
    outlet.extra.y = [];
    outlet.extra.x = [];
end