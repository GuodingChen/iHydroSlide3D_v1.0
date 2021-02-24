

function [iHydroSlide3D_output] = LandslideModel(rundata, cstep, tot_infil, SM_12_5, PWM, LandslideData)

        % get the data
        % hydrological output data
        Hydro_geoinfo = rundata.geoinfo;
        % note cohesion used N/m^2 in readsoil.m, but here we need KN/m^2
        effective_cohesion_map = rundata.Parameters.cohesion ./ 1000;
        kst_map = rundata.Parameters.hydrocond; % m/s
        friction_map = rundata.Parameters.friction;
        porosity_map = rundata.Parameters.porosity;
        OB_time = cstep; % h
        RainIn_matrix = tot_infil;
        
        WM_matrix = PWM;
        
        % landslide data need
        Landslide_geoinfo = LandslideData.geoinfo;
        mask_fine = LandslideData.mask;
        z_matrix = LandslideData.dem;
        slope_matrix = LandslideData.slope;
        aspect_matrix = LandslideData.aspect;
        
        % landslide model basic setting

        ellipse_dendity = LandslideData.Landslide_density;
        total_tile_number = LandslideData.Divide_tile_number;
        min_ae = LandslideData.min_ae; 
        max_ae = LandslideData.max_ae;
        min_be = LandslideData.min_be; 
        max_be = LandslideData.max_be;
        
        
        % read data
        % note that DEM, slope, aspect is smaller than mask_90
        % basical Topographical Data



        % read the Time-series file
        % The model accept three ways to select the observed time
        % 1. Generate continuous time series with user-defined start and end time
        % 2. Directly define discontinuous time series matrix
        % 3. Define individual time you are interested 
        % Users can use one of the three methods by commenting the specific code

        % Way 1
        % -----------
%         start_date = datenum('2012-07-03 20:00:00');
%         end_date = datenum('2012-07-04 20:00:00');
%         % compute date step: 1 h
%         date_step = datenum('2020-06-28 01:00:00')-datenum('2020-06-28 00:00:00');
% 
%         dateNum_array = (start_date : date_step : end_date)';
%         % prepare a dateStr array
%         dateStr_array =  strings(length(dateNum_array),1);
%         for i = 1 : length(dateNum_array)
%             date_STR = datestr( dateNum_array(i), 'yyyy-mm-dd HH:MM:SS');
%             new_date_STR = strcat( date_STR(1:4), date_STR(6:7), date_STR(9:10), date_STR(12:13));
%             dateStr_array(i) = new_date_STR;
%         end
%         observed_time_serise = dateStr_array;

        % Way 2
        % observed_time_serise = ["2012070320";"2012070403";"2012070406";"2012070412"];

        % Way 3
        % observed_time_serise = ["2012070412"];

        

        % user define the density of ellipse
        
        
        % 
        

        % data adjustment----------------
        z_matrix( z_matrix<0 ) = 0;
        aspect_matrix( aspect_matrix < 0 ) = 1;
        slope_matrix(slope_matrix == 0) = 1;
        
        all_slope = (slope_matrix / 180) * pi;
        all_aspect = (aspect_matrix / 180) * pi;

        % coordinate information (90 m)
        % row_90 = R_90.RasterSize(1);
        % column_90 = R_90.RasterSize(2);
        x_min_90 = Hydro_geoinfo.LongitudeLimits(1);  
        x_max_90 = Hydro_geoinfo.LongitudeLimits(2);  
        y_min_90 = Hydro_geoinfo.LatitudeLimits(1);
        y_max_90 = Hydro_geoinfo.LatitudeLimits(2);
        cell_size_90 = Hydro_geoinfo.CellExtentInLatitude;

        % coordinate information (12.5 m)
        row_12_5 = Landslide_geoinfo.RasterSize(1);
        column_12_5 = Landslide_geoinfo.RasterSize(2);
        x_min_12_5 = Landslide_geoinfo.LongitudeLimits(1);  
        x_max_12_5 = Landslide_geoinfo.LongitudeLimits(2);  
        y_min_12_5 = Landslide_geoinfo.LatitudeLimits(1);
        y_max_12_5 = Landslide_geoinfo.LatitudeLimits(2);
        cell_size_12_5 = Landslide_geoinfo.CellExtentInLatitude;

        %-----------------generate origin coordinate (resolution = 90 m)
        x_ori = (x_min_90 + cell_size_90/2) : cell_size_90 : (x_max_90 - cell_size_90/2);
        y_ori = (y_max_90 - cell_size_90/2) : -cell_size_90 : (y_min_90 + cell_size_90/2); % 
        [x_ori_all, y_ori_all] = meshgrid(x_ori, y_ori);

        %-----------------generate big coordinate (resolution = 12.5 m)
        x_big = (x_min_12_5 + cell_size_12_5/2) : cell_size_12_5 : (x_max_12_5 - cell_size_12_5/2);
        y_big = (y_max_12_5 - cell_size_12_5/2) : -cell_size_12_5 : (y_min_12_5 + cell_size_12_5/2); % 
        [x_big_all, y_big_all] = meshgrid(x_big, y_big);




        %-----------------generate the user-defined coordinate (resolution = 12.5 m)
        % let the coordinate transformation occurs within a Cartesian coordinate system
        % it will be More convenient
        x_cell = 12.5;
        y_cell = 12.5;
        x_min_set = 0;
        x_max_set = x_cell * (column_12_5-1);
        y_min_set = 0; 
        y_max_set = y_cell * (row_12_5-1);
        x = x_min_set : x_cell : x_max_set;
        y = y_max_set : -y_cell : y_min_set; % y coordinate should be carefully defined
        [x_all, y_all] = meshgrid(x,y);

        % z_out = z_all;
        % x_out = x_all;
        % y_out = y_all;
        % out_slope_matrix = slope_matrix;
        % out_aspect_matrix = aspect_matrix;
        
        % define FS matrix
        %-----% split the study area into sevaral tiles------
        [i_mask, j_mask] = find(~isnan(mask_fine));
        mask_vector = [i_mask, j_mask];
        
        
        CellNumber_InEachTile = floor( length(mask_vector) / total_tile_number );
        A_s = CellNumber_InEachTile * 12.5 * 12.5; % extent of the study area
        
        % initial the whole FS_3D map
        % it is no matter what cell value in FS3D_whole_map due to that the it will
        % be replaced by each tile in final step.
        FS3D_whole_map = mask_fine.*2;
        Volume_whole_map = mask_fine;
        Area_whole_map = mask_fine;
        PF_UnstableCount_map = mask_fine - 1;
        PF_TotalCount_map = mask_fine - 1;
        Volume_whole_map = single(Volume_whole_map);
        Area_whole_map = single(Area_whole_map);
        FS3D_whole_map = single(FS3D_whole_map);
        PF_UnstableCount_map = single(PF_UnstableCount_map);
        PF_TotalCount_map = single(PF_TotalCount_map);
        


        % user-defined random ellipse geo restriction

        a_e_range = (min_ae : 5 : max_ae)';
        b_e_range = (min_be : 5 : max_be)';

        % initial condition
        % FS_compare = mask_fine;
        % FS_compare(~isnan(FS_compare)) = nan;
        
        % transfer double matrix to single
        mask_vector = single(mask_vector);
        kst_map = single(kst_map);
        mask_fine = single(mask_fine);
        x_all = single(x_all); 
        y_all = single(y_all);
        z_all = single(z_matrix);
        x_big_all = single(x_big_all); 
        y_big_all = single(y_big_all);
        x_ori_all = single(x_ori_all); 
        y_ori_all = single(y_ori_all);


        
        

        % generate the file name
%         observed_time_step = observed_time_serise(time_i);
%         SM_file_name = strcat('SM_12_5_', observed_time_step, '.tif');
%         TotInfile_file_name = strcat('GOVar_tot_infil_', observed_time_step, '.tif');
%         % read .tiff file
%         addpath('./data_need'); 



        SM_12_5 = single(SM_12_5);
        RainIn_matrix = single(RainIn_matrix);
        FS3D_cell = cell(total_tile_number, 1);
        Volume_cell = cell(total_tile_number, 1);
        Area_cell = cell(total_tile_number, 1);
        PF_count_cell = cell(total_tile_number, 1);
        PF_unstable_cell = cell(total_tile_number, 1);
        RandomIndex_tile = cell(total_tile_number, 1);



        parfor tile_class = 1 : total_tile_number
            Monitor_variable = ['tile = ', num2str(tile_class), '   &   step = ', num2str(cstep)];
            disp(Monitor_variable)
            if tile_class == 1
                random_matrix = mask_vector(1:CellNumber_InEachTile,:);
            
            else
                random_matrix = mask_vector( ((tile_class-1) * CellNumber_InEachTile) : (tile_class*CellNumber_InEachTile),: );
            end
            % make the small FS maxtrix
            % tile_result = zeros(length(random_matrix), 3);
            FS3D_tile_matrix = ones(length(random_matrix), 1) * 1000;
            Volume_tile_matrix = zeros(length(random_matrix), 1);
            Area_tile_matrix = zeros(length(random_matrix), 1);
            PF_unstable_matrix = zeros(length(random_matrix), 1);
            PF_count_matrix = zeros(length(random_matrix), 1);        
            FS3D_tile_final = FS3D_tile_matrix;
            Volume_tile_final = Volume_tile_matrix;
            Area_tile_final = Area_tile_matrix;
            random_matrix_MinColumn = min(random_matrix(:,2));
            random_matrix_MaxColumn = max(random_matrix(:,2));
            % each tile corresponding to a random_matrix_LineIndex, and there
            % is impossible to get the intersection between any two matrix
            random_matrix_LineIndex = sub2ind(size(mask_fine), random_matrix(:,1), random_matrix(:,2));
            %disp(strcat('tile = ', num2str(tile_class)))
            ellipsoid_number = ellipse_dendity * 16 * A_s / (pi * (min_ae + max_ae) * ...
                (min_be + max_be)); 
            ellipsoid_number = ceil(ellipsoid_number);
            % begin to generate the random ellipse surface
            for num_ellipsoid = 1 : ellipsoid_number

                    %disp(['tile : ', num2str(tile_class), 'elli_N = ', num2str(num_ellipsoid)]);
                    %disp(strcat('ellipsoid number = ', num2str(num_ellipsoid)))

                    center_index = randi(length(random_matrix));
                    % (i, j) is the index in whole study region 
                    i = random_matrix(center_index, 1);
                    j = random_matrix(center_index, 2);

                %---------------------------------------------------  
                
                    x_WGS = x_big_all(i,j);
                    y_WGS = y_big_all(i,j);

                    %-----------find the slope location: i,j, in 90m map
                    L_location_row = ceil( (y_max_90 - y_WGS) / cell_size_90 );
                    L_location_column = ceil( (x_WGS - x_min_90) / cell_size_90 );

                    % determine the parameter in 90m map
                    % unit: KN/m^2
                    effective_cohesion = effective_cohesion_map(L_location_row, L_location_column);
                    n = porosity_map(L_location_row, L_location_column);
                    friction = (friction_map(L_location_row, L_location_column) / 180) * pi;
                    RainIn = RainIn_matrix(L_location_row, L_location_column);
                    % compute sr = SM (mm) / WM (mm);

                    WM = WM_matrix(L_location_row, L_location_column);
        %             sr_90 = SM_90 / WM;
        %             sr_90(sr_90>1) = 1;
        %             sr_mean = sr_90(L_location_row, L_location_column);
                    k_sat = kst_map(L_location_row, L_location_column);

                    % determine the parameter in 12m map 
                    slope = all_slope(i, j);
                    landslide_aspect = all_aspect(i, j);


                    %---determine the ellipse center in 12.5m map
                    x_center = x_all(i, j);
                    y_center = y_all(i, j);
                    z_center = z_all(i, j);
                    
                    % make a smaller matrix which contain a single ellipse
                    x_single_all = x_all((i-30) : (i+30), (j-30) : (j+30));
                    y_single_all = y_all((i-30) : (i+30), (j-30) : (j+30));
                    z_single_all = z_all((i-30) : (i+30), (j-30) : (j+30));
                    slope_single_all = all_slope((i-30) : (i+30), (j-30) : (j+30));
                    aspect_single_all = all_aspect((i-30) : (i+30), (j-30) : (j+30));
                    SM_12_5_single = SM_12_5((i-30) : (i+30), (j-30) : (j+30));
                    all_x_transition = (x_single_all-x_center) * cos(landslide_aspect) ...
                                           + (y_single_all-y_center) * sin(landslide_aspect);
                    all_y_transition = (y_single_all-y_center) * cos(landslide_aspect) ...
                                           - (x_single_all-x_center) * sin(landslide_aspect); 

                %---------------------soil parameters setup-------------------------

                    %effective_cohesion = c_map ;  % unit: KN/m^2, we note the read soil C
                    %(c_map) value is too big
                    
                    %effective_cohesion = 10;
                    
                    theta_s = 0.5;
                    theta_r = 0.05;
                    theta_0 = 0.2; %initial soil water content
                    
                    gamma_d = 27; %  specific weight of dry regolith, Unit: KN/m^3
                    gamma_w = 9.81;  % specific weight of water, Unit: KN/m^3    
                    %----------------geometrical compute: a_e, b_e, MC method is used
                    %cstep = 24 * 17 + 21;
                    % random the length and width of ellipsoid
                    ae_random_index = randi([1,length(a_e_range)]);
                    be_random_index = randi([1,length(b_e_range)]);
                    a_e_pre = a_e_range(ae_random_index);
                    b_e = b_e_range(be_random_index);


                    %a_e = (max_ae-min_ae) *rand(1) + min_ae;
                    project_length = a_e_pre * cos(slope);
                    %b_e = (max_be-min_be) *rand(1) + min_be;
                    H_c = 1.5;          % Hc is capillary pressure, m
                    t = OB_time * 3600; % seconds
                    slip_depth = sqrt(2 * k_sat * H_c * t / (theta_s - theta_0));
                    W_t = RainIn / 1000;
                    c_e = slip_depth; 

        %             [FSR_3D, ellipse_matrix] = slope_stability_3d(z_out, x_out, y_out, out_slope_matrix, out_aspect_matrix,...
        %                 project_length, b_e, c_e , x_center, y_center, z_center, x_cell, y_cell, ...
        %                 friction, effective_cohesion, gamma_d, gamma_w, landslide_aspect, n, ...
        %                 SM_12_5, WM, W_t);
                    
                    %------find the involved ellipse---------
                    % jugge_location is the index in small single region which can 
                    % include a single ellipse
                    judge_location = (all_x_transition).^2/project_length^2+(all_y_transition).^2/b_e^2;
                    ellipse_location = find(judge_location<1);  % value<1 represent that the points are included in project ellipse
                    [ellipse_i, ellipse_j] = find(judge_location<1);
                    ellipse_matrix = [ellipse_i, ellipse_j];
                    ellipse_matrix = single(ellipse_matrix);
                  %---------------------------get Ellipsoidal region----------------
                    z_grid = z_single_all(ellipse_location);
                    x_transition = all_x_transition(ellipse_location);
                    y_transition = all_y_transition(ellipse_location);
                    slope_c = slope_single_all(ellipse_location);
                    aspect_c = aspect_single_all(ellipse_location);
                    % find soil moisture
                    SM = SM_12_5_single(ellipse_location);
                    sr = SM / WM; % SM and WM is mm
                    sr(sr>1) = 1;

                    main_slope = mean(slope_c,1);
                    main_aspect = mean(aspect_c,1);
                    a_e = project_length/cos(main_slope);
                    % ensure there are no nodata cell in SM matrix
                    if sum(isnan(SM)) == 0


                        [FSR_3D, Failure_volumn, Failure_area] = StabilityEQ_3D(a_e, b_e, c_e , z_center, x_cell, y_cell, ...
                            friction, effective_cohesion, gamma_d, gamma_w, landslide_aspect, n, ...
                          W_t, x_transition, y_transition, main_slope, main_aspect, ...
                            z_single_all, ellipse_location, ellipse_matrix, z_grid, judge_location, ...
                            ellipse_i, ellipse_j, sr);

                        % transform the ellipse_location index to (i,j) index map
                        [MovingWindow_center_i, ~] = find(y_single_all == y_center);
                        [~, MovingWindow_center_j] = find(x_single_all == x_center);
                        MovingWindow_center_i = unique(MovingWindow_center_i);
                        MovingWindow_center_j = unique(MovingWindow_center_j);

                        % update the ellipse_matrix into (i,j) index map (whole region)
                        ellipse_matrix(:,1) = ellipse_matrix(:,1) + (i - MovingWindow_center_i);
                        ellipse_matrix(:,2) = ellipse_matrix(:,2) + (j - MovingWindow_center_j);

                        % exclude the index which is out of the random_maxtrix,
                        % FS_ellipse_matrix represents the coordinate corresponding to
                        % the mask_fine (the total area)
                        %FS_ellipse_matrix = ellipse_matrix()
                        FS_column_index = find((ellipse_matrix(:,2) >= random_matrix_MinColumn)...
                            & (ellipse_matrix(:,2) <= random_matrix_MaxColumn));

                        FS_ellipse_matrix = ellipse_matrix(FS_column_index, :);
                        %FS_ellipse_matrix = intersect(random_matrix, ellipse_matrix, 'rows');

                        % FS_tile_index represents the normal index (begin with 1) in
                        % each tile of FS3D_tile_maxtrix
                        FS_ellipse_matrix_LineIndex = sub2ind(size(mask_fine), FS_ellipse_matrix(:,1), FS_ellipse_matrix(:,2));
                        %FS_tile_index = ismember(random_matrix, FS_ellipse_matrix, 'rows');

                        FS_tile_index = ismember(random_matrix_LineIndex, FS_ellipse_matrix_LineIndex);
                        FS3D_tile_matrix(FS_tile_index) = FSR_3D;



                        % count the cell considered into random test
                        PF_count_matrix(FS_tile_index) = PF_count_matrix(FS_tile_index) + 1;

                        % count the test which computed as unstable cell
                        if FSR_3D < 1
                            Volume_tile_matrix(FS_tile_index) = Failure_volumn;
                            Area_tile_matrix(FS_tile_index) = Failure_area;
                            PF_unstable_matrix(FS_tile_index) = PF_unstable_matrix(FS_tile_index) + 1;
                        end
                    end

                    % update the FS value

                    FS3D_tile_final = min(FS3D_tile_final, FS3D_tile_matrix);
                    Volume_tile_final = max(Volume_tile_final, Volume_tile_matrix);
                    Area_tile_final = max(Area_tile_final, Area_tile_matrix);

            end
            % 
            % adjust the nodata in raster

            FS3D_tile_final(FS3D_tile_final == 1000) = nan;
            Volume_tile_final(Volume_tile_final == 0) = nan;
            Area_tile_final(Area_tile_final == 0) = nan;
            % storage each tile in tile class
            FS3D_cell{tile_class, 1} = FS3D_tile_final;
            Volume_cell{tile_class, 1} = Volume_tile_final;
            Area_cell{tile_class, 1} = Area_tile_final;
            RandomIndex_tile{tile_class, 1} = random_matrix_LineIndex;
            PF_count_cell{tile_class, 1} = PF_count_matrix;
            PF_unstable_cell{tile_class, 1} = PF_unstable_matrix;

    %        FS3D_whole_map(random_matrix_LineIndex) = FS3D_tile_final;

            % total count, computed the percentage of the unstable area
    %         PF_TotalCount_map(random_matrix_LineIndex) = PF_count_matrix;
    %         PF_UnstableCount_map(random_matrix_LineIndex) = PF_unstable_matrix;
    %         
    %         PF_whole_map = PF_UnstableCount_map ./ PF_TotalCount_map;
    %         PF_whole_map(PF_TotalCount_map == 0) = 0;


        end  
        disp('# # # # # # End of current time calculation # # # # # # ')
        % reconstruct the tiles to the whole FS/PF map
        % In theory, the cell number of FS ~= nan will equal to PF ~= nan, but
        % this program initial the PF_unstable_matrix as zero matrix, and dose
        % not consider the cell which is not randomly selected. Therefor, PF ~= nan 
        % will have more cell number than that of FS ~= nan
        for i = 1 : total_tile_number
            EachTile_index = RandomIndex_tile{i};
            FS3D_whole_map(EachTile_index) = FS3D_cell{i};
            PF_TotalCount_map(EachTile_index) = PF_count_cell{i};
            PF_UnstableCount_map(EachTile_index) = PF_unstable_cell{i};
            Volume_whole_map(EachTile_index) = Volume_cell{i};
            Area_whole_map(EachTile_index) = Area_cell{i};

        end
        PF_whole_map = PF_UnstableCount_map ./ PF_TotalCount_map;
        PF_whole_map(PF_TotalCount_map == 0) = 0;
        iHydroSlide3D_output.FS3D_whole_map = FS3D_whole_map;
        iHydroSlide3D_output.PF_whole_map = PF_whole_map;
        iHydroSlide3D_output.Volume_whole_map = Volume_whole_map;
        iHydroSlide3D_output.Area_whole_map = Area_whole_map;
        % output the file in each observed time step
%         path_output = './FS_PF_output/';
%         FS_FileName = strcat(path_output, 'whole_FS_map', observed_time_step);
%         PF_FileName = strcat(path_output, 'whole_PF_map', observed_time_step);
%         Volume_Filename = strcat(path_output, 'Volume', observed_time_step);
%         Area_Filename = strcat(path_output, 'Area', observed_time_step);
%         geotiffwrite(FS_FileName, FS3D_whole_map, Landslide_geoinfo);  
%         geotiffwrite(PF_FileName, PF_whole_map, Landslide_geoinfo);
%         geotiffwrite(Volume_Filename, Volume_whole_map, Landslide_geoinfo);
%         geotiffwrite(Area_Filename, Area_whole_map, Landslide_geoinfo);
%         
%         disp('# # # # # # All time calculation ends # # # # # # ')
        
        
end



