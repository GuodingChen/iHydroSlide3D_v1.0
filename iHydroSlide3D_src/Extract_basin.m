function [mask,outlet_row,outlet_col,user_outlet_row,user_outlet_col,rel_err,new_est_darea] = Extract_basin(facc,fdir,given_outlet,basin_area,mapinfo,mask_key,tolerance,max_radius)
%Function that identify drainage of given location (outlet) from a region
%Version 1.0 - February, 2014
%Written by Humberto Vergara - humber@ou.edu
%Hydrometeorology and Remote-Sensing (HyDROS) Laboratory - http://hydro.ou.edu
%Advanced Radar Research Center (ARRC) - http://arrc.ou.edu
%The University of Oklahoma - http://ou.edu

%% Locate correct outlet from approximate given location
lon = given_outlet.longitude; lat = given_outlet.latitude;

%% Return grid size in meters
switch mapinfo.SpatialRef.CoordinateSystemType
    case 'geographic'
        deltax=mapinfo.PixelScale(1);
        deltay=mapinfo.PixelScale(2);
        xmin=min(mapinfo.BoundingBox(:,1));
        xmax=max(mapinfo.BoundingBox(:,1));
        ymin=min(mapinfo.BoundingBox(:,2));
        ymax=max(mapinfo.BoundingBox(:,2));
        
        if strcmpi(mapinfo.SpatialRef.ColumnsStartFrom,'north')==1
            new_row=floor((ymax-lat)/deltay)+1;
            user_outlet_row =new_row;
        else
            new_row=floor((lat-ymin)/deltay)+1;
            user_outlet_row =new_row;
        end
        
        if strcmpi(mapinfo.SpatialRef.RowsStartFrom,'west')==1
            new_col=floor((lon-xmin)/deltax)+1;
            user_outlet_col =new_col;
        else
            new_col=floor((xmax-lon)/deltax)+1;
            user_outlet_col =new_col;
        end
        cellsize = deg2km(deltax)*1000;
    case 'planar'
        %Re-project outlet coordinates
        [x, y] = projfwd(mapinfo, lat, lon);

        %Estimate pixel location from coordinates
        [row,col] = map2pix(mapinfo.SpatialRef,x,y);
        new_row = round(row); new_col = round(col);
        user_outlet_row = round(row); user_outlet_col = round(col);

        %Estimate drainage area of user-provided outlet location
        %Use pixel size (assumed in meters for now)
        cellsize = mapinfo.PixelScale(1);
    otherwise
        error('The coordinate system type of the mask file is not supported!');
end
        
est_darea = facc(new_row,new_col)*((cellsize/1000)^2);

%Given Drainage Area in km^2
if (isempty(basin_area) == 0)
    darea = basin_area;
    %Compute the relative difference between given drainage area and
    %estimated drainage area to determined if outlet is located correctly
    rel_err = abs((darea - est_darea)/darea)*100;
else
    darea = [];
    new_est_darea = est_darea;
    rel_err = 0;
end

radius_distance = 2; %In pixels
while (rel_err > tolerance)
    %Look in vecinity
    search_pixels_rows = user_outlet_row-radius_distance:user_outlet_row+radius_distance;
    search_pixels_cols = user_outlet_col-radius_distance:user_outlet_col+radius_distance;
    %ignore out of bounds pixels
    search_pixels_rows = search_pixels_rows(search_pixels_rows>0);
    search_pixels_cols = search_pixels_cols(search_pixels_cols>0);
    max_row = max(search_pixels_rows);
    max_col = max(search_pixels_cols);
    %Check for mas allowed radius search
    if (radius_distance*(cellsize/1000) > max_radius || max_row > size(facc,1) || max_col > size(facc,2))
        %if the maximum radius search is exceeded, force choosing the
        %closest grid and issue a warning.
        fprintf('\n\nWARNING: Could not find an accurate location within allowed domain for outlet at:\nlat: %g lon: %g with known drainage area %g km^2',lat,lon,darea);
        fprintf('\nForcing location at row: %g col: %g with drainage area %g km^2 and %g percent error.\n\n', new_row, new_col, new_est_darea, rel_err);
        break;
    end
    
    radius_search = facc(search_pixels_rows,search_pixels_cols).*((cellsize/1000)^2);
    radius_search_diff = abs(radius_search - darea);
    [new_i,new_j] = find(radius_search_diff == min(radius_search_diff(:)));
    new_row = (new_i-radius_distance-1) + user_outlet_row; 
    new_col = (new_j-radius_distance-1) + user_outlet_col;
    %Get new drainage area in km2
    new_est_darea = facc(new_row,new_col)*((cellsize/1000)^2);

    %Re-compute relative difference
    rel_err = abs((darea - new_est_darea)/darea)*100;
    %Increase radius of search
    radius_distance = radius_distance + 2;
end

%Set the specified outlet coordinates
outlet_row = new_row;
outlet_col = new_col;
new_est_darea = facc(outlet_row,outlet_col)*((cellsize/1000)^2);
    
%% Track cells draining to outlet
%Which flow direction code convention?
if (max(fdir(isnan(fdir)==0)) == 128)
    option = 1;
else if (max(fdir(isnan(fdir)==0)) == 8)
        option = 2;
    else
        error('Flow direction coding not supported');
    end
end

%% MASK
if (mask_key == 1)
    %Pre-allocate basin mask
    mask = zeros(size(fdir));
    [ys,xs]=size(fdir);
    %Include catchment outlet in mask
    mask(outlet_row,outlet_col) = 1;

    %How many grids need to be accounted for?
    n_contributing_cells = facc(outlet_row,outlet_col);

    %Flow direction key convention (toward pixel of interest)
    %Order here strictly tied to variable "surroundings" defined below
    fdircodes(1,:) = [2, 4, 8, 16, 32, 64, 128, 1];
    fdircodes(2,:) = [2, 1, 8, 7, 6, 5, 4, 3];

    %Specify the indices of surrounding cells
    %The following order of pixels is strictly tied to the flow direction
    %convention variable "fdircodes" defined above
    surroundings = [new_row-1,new_col-1; %fdir key: 2(2)
                    new_row-1,new_col; %fdir key: 4(1)
                    new_row-1,new_col+1; %fdir key: 8(8)
                    new_row,new_col+1; %fdir key: 16(7)
                    new_row+1,new_col+1; %fdir key: 32(6)
                    new_row+1,new_col; %fdir key: 64(5)
                    new_row+1,new_col-1; %fdir key: 128(4)
                    new_row,new_col-1]; %fdir key: 1(3)

    %Ignore out of bounds grids
    [surr_r] = find(min(surroundings,[],2) > 0);
    surroundings = surroundings(surr_r,:);

    %Select only contributing cells from surroundings
    %Total of cells accounted for at this point
    n_accounted_cells = 0;
    cont = 0;
    for neig = 1:size(surroundings,1)
        if (fdir(surroundings(neig,1),surroundings(neig,2)) == fdircodes(option,neig))
            cont = cont + 1;
            n_accounted_cells = n_accounted_cells + 1;
            drain_group(cont,:) = surroundings(neig,:);
            mask(surroundings(neig,1),surroundings(neig,2)) = 1;
        end
    end

    %Loop through cells until all contributing grids are accounted for
    while (n_accounted_cells < n_contributing_cells)
        cont = 0;
        for cells = 1:size(drain_group,1)
            new_row = drain_group(cells,1);
            new_col = drain_group(cells,2);
            %Specify the indices of surrounding cells
            %The following order of pixels is strictly tied to the flow direction
            %convention variable "fdircodes" defined above
            surroundings = [new_row-1,new_col-1; %fdir key: 2(2)
                            new_row-1,new_col; %fdir key: 4(1)
                            new_row-1,new_col+1; %fdir key: 8(8)
                            new_row,new_col+1; %fdir key: 16(7)
                            new_row+1,new_col+1; %fdir key: 32(6)
                            new_row+1,new_col; %fdir key: 64(5)
                            new_row+1,new_col-1; %fdir key: 128(4)
                            new_row,new_col-1]; %fdir key: 1(3)

            %Ignore out of bounds grids
            [surr_r] = find(min(surroundings,[],2) > 0);
            surroundings = surroundings(surr_r,:);

            %Select only contributing cells from surroundings
            for neig = 1:size(surroundings,1)
                if (surroundings(neig,1)<=ys && surroundings(neig,2)<=xs && fdir(surroundings(neig,1),surroundings(neig,2)) == fdircodes(option,neig))
                    cont = cont + 1;
                    n_accounted_cells = n_accounted_cells + 1;
                    new_drain_group(cont,:) = surroundings(neig,:);
                    mask(surroundings(neig,1),surroundings(neig,2)) = 1;
                end
            end
        end
        %Update group of outlets
        drain_group = new_drain_group;
        new_drain_group = [];
    end
else
    mask = [];
end