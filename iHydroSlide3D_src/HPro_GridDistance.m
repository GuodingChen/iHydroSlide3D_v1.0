function [hlen,routto,routfrom,routfrom_n,slope,seqmatrix,alpha,beta] = HPro_GridDistance(rundata)
%Determine inter-grid distance for HPro
%Version 1.0 - February, 2014
%Written by Humberto Vergara - humber@ou.edu
%Hydrometeorology and Remote-Sensing (HyDROS) Laboratory - http://hydro.ou.edu
%Advanced Radar Research Center (ARRC) - http://arrc.ou.edu
%The University of Oklahoma - http://ou.edu
%
%Arguments
%rundata - simulation data

fprintf('Generating flow parameter grids: ');
flowdir = rundata.fdir;
facc = rundata.facc;
gridres = rundata.gridres;
dem_path = rundata.dem_path;
slope_path = rundata.slope_path;
outlet = rundata.outlet;

%For now, it is assumed projection is equidistance
hlen = zeros(size(flowdir))+(gridres.DeltaX*sqrt(2));
hlen(isnan(flowdir) == 1) = NaN;

if (isempty(slope_path) == 1)
    slope = zeros(size(flowdir)).*hlen;
else
    slope = rundata.slope;
end

if (isempty(dem_path) == 0)
    dem = rundata.dem;
end 

%Create grids of flow destination
routto = zeros([size(flowdir),2]);
routfrom = zeros([size(flowdir),7,2]);
routfrom_n = zeros(size(flowdir));

keycode = nanmax(flowdir(:));
switch keycode
    case 8 % 1-8 Direction keys    
        hlen(mod(flowdir,2) == 0) = gridres.DeltaX;
        for i = 1:size(flowdir,1)
            for j = 1:size(flowdir,2)
                if (isnan(flowdir(i,j)) == 0)
                    new_i = i;
                    new_j = j; 
                    switch flowdir(new_i,new_j)
                        case 1
                            new_i = new_i + 1;
                        case 2
                            new_i = new_i + 1;
                            new_j = new_j + 1;
                        case 3
                            new_j = new_j + 1;
                        case 4
                            new_i = new_i - 1;
                            new_j = new_j + 1;
                        case 5
                            new_i = new_i - 1;
                        case 6
                            new_i = new_i - 1;
                            new_j = new_j - 1;
                        case 7
                            new_j = new_j - 1;
                        case 8
                            new_i = new_i + 1;
                            new_j = new_j - 1;
                    end
                    
                    %Compute Slope if not provided
                    if (isempty(slope_path) == 1)
                        slope(i,j) = max(0.004,(dem(i,j)-dem(new_i,new_j))/hlen(i,j));
                    end

                    %Assign destination and departure grids corresponding row and column indices
                    routfrom_n(new_i,new_j) = routfrom_n(new_i,new_j) + 1;
                    routfrom(new_i,new_j,routfrom_n(new_i,new_j),1) = i;
                    routfrom(new_i,new_j,routfrom_n(new_i,new_j),2) = j;
                    if (isnan(flowdir(new_i,new_j))==1)
                        routto(i,j,1) = 0;
                        routto(i,j,2) = 0;
                    else
                        routto(i,j,1) = new_i;
                        routto(i,j,2) = new_j;
                    end
                else
                    routto(i,j,1) = NaN;
                    routto(i,j,2) = NaN;
                    routfrom(i,j,:,1) = NaN;
                    routfrom(i,j,:,2) = NaN;
                end
            end
        end
    case 128 %1-2-4-8-16-32-64-128 Direction keys
        hlen(flowdir == 1 | flowdir == 4 | flowdir == 16 | flowdir == 64) = gridres.DeltaX;
        for i = 1:size(flowdir,1)
            for j = 1:size(flowdir,2)
                if (isnan(flowdir(i,j)) == 0)
                    new_i = i;
                    new_j = j; 
                    switch flowdir(new_i,new_j)
                        case 4
                            new_i = new_i + 1;
                        case 2
                            new_i = new_i + 1;
                            new_j = new_j + 1;
                        case 1
                            new_j = new_j + 1;
                        case 128
                            new_i = new_i - 1;
                            new_j = new_j + 1;
                        case 64
                            new_i = new_i - 1;
                        case 32
                            new_i = new_i - 1;
                            new_j = new_j - 1;
                        case 16
                            new_j = new_j - 1;
                        case 8
                            new_i = new_i + 1;
                            new_j = new_j - 1;
                    end
                    
                    %Compute Slope if not provided
                    if (isempty(slope_path) == 1)
                        slope(i,j) = max(0.004,(dem(i,j)-dem(new_i,new_j))/hlen(i,j));
                    end

                    %Assign destination and departure grids corresponding row and column indices
                    routfrom_n(new_i,new_j) = routfrom_n(new_i,new_j) + 1;
                    routfrom(new_i,new_j,routfrom_n(new_i,new_j),1) = i;
                    routfrom(new_i,new_j,routfrom_n(new_i,new_j),2) = j;
                    if ((i == outlet.y) && (j == outlet.x))
                        routto(i,j,1) = 0;
                        routto(i,j,2) = 0;
                    else
                        routto(i,j,1) = new_i;
                        routto(i,j,2) = new_j;
                    end
                else
                    routto(i,j,1) = NaN;
                    routto(i,j,2) = NaN;
                    routfrom(i,j,:,1) = NaN;
                    routfrom(i,j,:,2) = NaN;

                end
            end
        end
    otherwise
        error('Flow Direction Coding not recognized');
end


%% Kinematic Wave Routing 
%Create up-to-downstream computational grid
% facc_classes = unique(facc(isnan(facc)==0));
% N_fc = length(facc_classes);
% template_i = zeros(N_fc,length(facc(isnan(facc)==0)));
% template_j = zeros(N_fc,length(facc(isnan(facc)==0)));
% prevL = 1;
% for fc = 1:N_fc
%     contours = facc.*NaN;
%     contours(facc == facc_classes(fc)) = 1;
%     c = 1;
%     is = [];
%     js = [];
%     for i = 1:size(facc,1)
%         for j = 1:size(facc,2)
%             if (contours(i,j) == 1);
%                 is(c) = i;
%                 js(c) = j;
%                 c = c + 1;
%             end
%         end
%     end
%     maxL = max(length(js),prevL);
%     prevL = maxL;
%     template_i(fc,1:length(is)) = is;
%     template_j(fc,1:length(js)) = js;
% end
% 
% seqmatrix(:,:,1) = template_i(:,1:maxL);
% seqmatrix(:,:,2) = template_j(:,1:maxL);
%------------------------
%changed by Sheng Wang
facc_classes = tabulate(facc(isnan(facc)==0));
N_fc = length(facc_classes(:,1));  
maxL = max(facc_classes(:,2));
seqmatrix = zeros(N_fc,maxL,2);

for fc = 1:N_fc
    contours = facc.*NaN;
    contours(facc == facc_classes(fc,1)) = 1;
    c = 1;
    is = [];
    js = [];
    for i = 1:size(facc,1)
        for j = 1:size(facc,2)
            if (contours(i,j) == 1)
                is(c) = i;
                js(c) = j;
                c = c + 1;
            end
        end
    end
    seqmatrix(fc,1:length(is),1) = is;
    seqmatrix(fc,1:length(js),2) = js;
end

%Compute parameters
Parameters = rundata.Parameters;
%Overland kinematic parameters are fixed - Model as a wide rectangular
%channel
%sq_hslope - hillslope, estimated from main flow direction slope: arbitrary value
hslope = slope;
hslope(facc < Parameters.TH) = slope(facc < Parameters.TH).*0.5;
sq_hslope = sqrt(hslope);
%Overland Conveyance Parameters
alpha_0 = 1./(Parameters.COEM.*sq_hslope);
beta_0 = facc.*0+3/5;

%Estimate channel kinematic wave parameters
switch Parameters.Kinematic_Mode
    case 'Input'
        %Don't do anything, it is done in "readgrids.m"
        alpha = Parameters.ALPHA;
        beta = Parameters.BETA;
    case 'RatingCurve'
        %Propagate parameter values upstream from the outlet
        beta = 0.6+(facc./nanmax(facc(:))).*(Parameters.BETA-0.6);
        alpha = (facc./nanmax(facc(:))).*Parameters.ALPHA;
    case 'ChannelO'
        %Estimate parameters from channel characteristics at the outlet
        beta = 0.6+(facc./nanmax(facc(:))).*0.2;
        %Channel Width at the outlet
        B = Parameters.CHW;
        %Estimate channel width for all channels in basin based on width at the outlet
        Bgrid = (facc./nanmax(facc(:))).*B;
        P = Bgrid; %Wetted Perimeter P = Bgrid+2.*y ??
        alpha = ((P.^2/3)./(Parameters.RIVER.*sqrt(slope))).^beta;
    case 'ChannelB'
        %Estimate parameters from channel characteristics across the basin
        %streams
        beta = 0.6+(facc./nanmax(facc(:))).*0.2;
        Bgrid = Parameters.CHW;
        P = Bgrid; %Wetted Perimeter P = Bgrid+2.*y ??
        alpha = ((P.^2/3)./(river.*sq_slope)).^beta;
end

%Assign overland grids the correct alpha and beta values
alpha(facc < Parameters.TH) = alpha_0(facc < Parameters.TH);
beta(facc < Parameters.TH) = beta_0(facc < Parameters.TH);

fprintf('Done\n');