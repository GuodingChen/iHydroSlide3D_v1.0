function FoS_inf = StabilityEQ_1D(slope, effective_cohesion, friction, n, k_s, cstep, RainIn, SR)
% unit use in FS calculate (SLIDE model), as follow:
%t: cstep is in hr, t = cstep * 3600 ---->seconds
%A: KN/m^2, A is a parameter depending on the kind of soil 
%   and is linked to the peak shear stress at failure
%slope: (degree)бу
%friction: (degree)бу
%W_t: mm, should be convert to m
%z_t: slip depth, m
    t = cstep * 3600; 
    A = 60;          % ori = 20 KPa
    lambda = 25;        % ori = 0.4
    Sr = SR;  %the degree of saturation = W_t / W_m
    alpha = 5;      %3.4;
% W_t: Wt is total cell mean water of the three soil layers, 
%       an intermediate variable that is simulated from CREST
    W_t = RainIn / 1000; %rain convert mm to m
    H_c = 1.5;          % Hc is capillary pressure, m
%set the saturated soil water content and initial water content of the soil.
    theta_s = 0.5; %saturated soil water content
    theta_0 = 0.2; %initial soil water content
    gamma_s = 20;% natural weight, unit weight of soil, KN/m^3; dry or saturated? 
    gamma_w = 9.8;
    G_s = gamma_s / gamma_w;
    n_w = n .* (1 - Sr);
    Gamma = G_s .* (1 - n) + n .* Sr;
    z_t = sqrt(2 .* k_s .* H_c .* t / (theta_s - theta_0));
    m_t = W_t ./ (z_t .* n.* (1-Sr));
    OMG = 2 ./ (sind(2*slope) .* z_t .* gamma_w);
    m_t(m_t>1) = 1;                  
    m_t(Sr==1) = 1;
    % c2 is apparent cohesion
    c1 = effective_cohesion./1000; % note cohesion used N/m^2 in readsoil.m
    c2 = A * Sr .* (1-Sr).^lambda .* (1-m_t).^alpha;
    c_final = c1 + c2;
%----------FS------------------
    FoS_inf = ( cotd(slope) .* tand(friction) .* (Gamma + m_t .* (n_w - 1))...
                + c_final .* OMG ) ./ (Gamma + m_t .* n_w);
end
