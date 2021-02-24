function [I,R,W0,i] = Derive_vic(eP,Pars,W,RI)
%This function solves the variable infiltration curve (vic) model in CREST.
%I = Amount of water available for infiltration
%Excess Rainfall
%eP = Effective Precipitation (amount of precip that makes it to the soil)
%vic_Pars = Parameters for vic/CREST

%Parameters
Wm = Pars.Wm;
b = Pars.b;

%Equations
%Compute I
B = 1+b;
A = W./Wm;
im = Wm.*B;
i = im.*(1-(1-A).^(1./B));
I = zeros(size(W));
R = zeros(size(W));
W0 = zeros(size(W));

%if ((i+eP) >= im)
I((i+eP) >= im) = Wm((i+eP) >= im) - W((i+eP) >= im)+RI((i+eP) >= im);
R((i+eP) >= im) = eP((i+eP) >= im) - I((i+eP) >= im);
R(R < 0) = 0;
W0((i+eP) >= im) = Wm(((i+eP) >= im));

%else
I((i+eP) < im) = (Wm((i+eP) < im) - W((i+eP) < im)) - Wm((i+eP) < im).*(1-(i((i+eP) < im)+eP((i+eP) < im))./im((i+eP) < im)).^B((i+eP) < im)+RI((i+eP) < im);
R((i+eP) < im) = eP((i+eP) < im) - I((i+eP) < im);
R(R < 0) = 0;
W0((i+eP) < im) = W((i+eP) < im) + eP((i+eP) < im) - R((i+eP) < im);
W0(W0 > Wm) = Wm(W0 > Wm);

