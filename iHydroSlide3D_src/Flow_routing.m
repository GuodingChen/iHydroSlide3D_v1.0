function [cQ,q,interflow] = Flow_routing(Overland_H,Interflow_H,pQ,pq,nexTimeUnder,alpha,beta,th,slope_th,dt,dx,RD)
%Basic flow routing framework
%Channel (river) routing: Kinematic Wave (Implicit Linear Scheme)
%Overland routing: Kinematic Wave (Implicit Linear Scheme)
%Interflow (subsurface) routing: Simple linear storage
%
%Version 1.0 - March, 2014
%Written by Humberto Vergara - humber@ou.edu
%
%Arguments:
%
% Overland_H: Infiltration excess (i - f) grid in meters per second
% Interflow_H - interflow rate (subsurface flow) in millimeters
%
% Current State of stream network
% pQ: Current channel flow
% pq: Current overland flow per unit width
%
% PARAMETERS
% th: Drainage area threshold to separate overland from channel grids
% beta: Exponent of momentum equation A = alphaQ^beta
% alpha: Coefficient of momentum equation A = alphaQ^beta
% nexTimeUnder: Interflow crossing-time grid (seconds)
%
%Spatiotemporal discretization
% dt: Time step in seconds
% dx: Stream length (grid side size) in meters
%
% RD - Struct variable with all simulation info as defined in the
% configuration file
%
% Grids containing flow direction info
%  RD.routto: Grid indicating where grids flow to 
%  RD.routfrom: Grid indicating where grids receive flow from
%  RD.routfrom_n: How many grids flow into a given grid point
%
% Grids containing flow accumulation info
%  RD.facc: Flow accumulation grid
%  RD.seqmatrix: up-to-downstream computational grid

%% PRE-PROCESSING
%Pre-allocate output variables
%Grid containing resulting channel flow
cQ = zeros(size(pQ));
%Grid containing resulting overland flow
q = zeros(size(pQ));
%Grid containing resulting interflow
interflow = zeros(size(pQ));

slope = (RD.sq_slope).^2;

%% ROUTING COMPUTATIONS
%Loop through computational grid
for fc = 1:size(RD.seqmatrix,1) %Each row contains grids with equal flow accumulation
    col = 1; %Column counter for computational grid        
    while (RD.seqmatrix(fc,col,1) > 0 && col < size(RD.seqmatrix,2))
        %Coordinates (array indices) of current grid
        i = RD.seqmatrix(fc,col,1);
        j = RD.seqmatrix(fc,col,2); 
        
        %Determine whether the current grid is overland or channel
        if (RD.facc(i,j) < th && slope(i,j) < slope_th)
            %Overland grid

            %Start Overland Flow computations
            %Determine whether current grid is headwater or not
            if (RD.facc(i,j) == 0)
                %Headwater grid
                %No upstream flow
                up_q = 0;
            else
                %Interior grid
                %Sum of all upstream flow
                up_q = 0;
                for uc = 1:RD.routfrom_n(i,j)
                    up_q = up_q + q(RD.routfrom(i,j,uc,1),RD.routfrom(i,j,uc,2));
                end
            end

            %Finite difference computation of overland flow at current grid point
            avg_q = max((pq(i,j)+up_q)/2,1e-20); %this value need to be > 0, otherwise: 0^-a = Inf
            A = (dt/dx(i,j)) * up_q;
            B = alpha(i,j)*beta(i,j)*pq(i,j)*(avg_q)^(beta(i,j)-1);
            C = dt*Overland_H(i,j);
            D = dt/dx(i,j);
            E = alpha(i,j)*beta(i,j)*(avg_q)^(beta(i,j)-1);
            q(i,j) = (A + B + C)/(D + E); %Overland flow (m^3/m.s)
            %End overland flow computations
            %----------------------------------------------------------               

            %Start interflow computations
            currentsec = 0;
            previoussec = 0;
            previousxy = [];
            currentxy = [i,j];
            while (currentsec < dt)
                if (isnan(RD.routto(currentxy(1),currentxy(2),1)) == 0)
                    previoussec = currentsec;
                    previousxy = currentxy;
                    currentsec = currentsec + nexTimeUnder(currentxy(1),currentxy(2));
                    if (RD.routto(currentxy(1),currentxy(2),1) ~= 0)
                        currentxy = [RD.routto(currentxy(1),currentxy(2),1),RD.routto(currentxy(1),currentxy(2),2)];
                        %If reaching a channel cell
                        if (RD.facc(currentxy(1),currentxy(2)) >= th || slope(i,j) >= slope_th)
                            previoussec = 0;
                            currentsec = dt;
                            %Add interflow to overland infiltration excess at the channel grid
                            Overland_H(currentxy(1),currentxy(2)) = Overland_H(currentxy(1),currentxy(2)) + Interflow_H(i,j).*(0.001/dt);
                            Interflow_H(i,j) = 0;
                            break;
                        end
                    else
                        currentxy = [];
                        if (dt > currentsec)
                                previousxy = [];
                        end
                        break;
                    end
                else
                    error('Interflow water flowing out of basin from a point other that the basin outlet'); 
                end
            end
            cAmount = (dt - previoussec) / (currentsec - previoussec);
            if (isnan(cAmount) == 1)
                error('Interflow_H');
            end
            pAmount = 1 - cAmount;
            if (isempty(currentxy) == 0)
                interflow(currentxy(1),currentxy(2)) = interflow(currentxy(1),currentxy(2)) + Interflow_H(i,j)*cAmount;
            end

            if (isempty(previousxy) == 0)
                interflow(previousxy(1),previousxy(2)) = interflow(previousxy(1),previousxy(2)) + Interflow_H(i,j)*pAmount;
            end
            %End interflow computations
            %----------------------------------------------------------
        else
            %Channel grid

            %First do overland flow routing in channel grid
            %Sum of all upstream flow
            up_q_i = 0;
            up_facc = 0;
            up_slope = 0;
            for uc = 1:RD.routfrom_n(i,j)
                up_q_i(uc) = q(RD.routfrom(i,j,uc,1),RD.routfrom(i,j,uc,2));
                up_facc(uc) = RD.facc(RD.routfrom(i,j,uc,1),RD.routfrom(i,j,uc,2));
                up_slope(uc) = slope(RD.routfrom(i,j,uc,1),RD.routfrom(i,j,uc,2));
            end

            %Only get upstream overland flow from overland grids
            up_q = sum(up_q_i(up_facc < th & up_slope < slope_th));

            %Compute overland flow at current grid point
            avg_q = max((pq(i,j)+up_q)/2,1e-20); %this value need to be > 0, otherwise: 0^-a = Inf
            A = (dt/dx(i,j)) * up_q;
            B = alpha(i,j)*beta(i,j)*pq(i,j)*(avg_q)^(beta(i,j)-1);
            C = dt*Overland_H(i,j);
            D = dt/dx(i,j);
            E = alpha(i,j)*beta(i,j)*(avg_q)^(beta(i,j)-1);
            q(i,j) = (A + B + C)/(D + E); %Overland flow (m^3/m.s)
            %End overland routing computations in channel grid
            %----------------------------------------------------------

            %Channel Flow            
            %Get present Q value of upstream grid - Qi,j+1
            up_Q = 0;
            for uc = 1:RD.routfrom_n(i,j)
                %Sum of all upstream flow
                up_Q = up_Q + cQ(RD.routfrom(i,j,uc,1),RD.routfrom(i,j,uc,2));
            end

            %Compute Q at current grid point - i+1,j+1
            avgQ = max((pQ(i,j)+up_Q)/2,1e-20); %this value need to be > 0, otherwise: 0^-a = Inf
            A = (dt/dx(i,j)) * up_Q;
            B = alpha(i,j)*beta(i,j)*pQ(i,j)*(avgQ)^(beta(i,j)-1);
            C = dt*q(i,j);
            D = dt/dx(i,j);
            E = alpha(i,j)*beta(i,j)*(avgQ)^(beta(i,j)-1);
            cQ(i,j) = (A + B + C)/(D + E); %Channel flow (m^3/s)
            %End channel flow computations
            %----------------------------------------------------------
        end
        col = col + 1;
    end
end
