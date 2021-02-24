function [Cohesion,FRI,N,K]=readSoilLandPara(RD)

%Read soil and land cover parameters
%Version 1.0 - April, 2014
%Written by Ke Zhang

%
%Reference:
%V.T. Chow, D. R. Maidment and L. W. Mays, Applied Hydrology, 1988, page 115. 
LC=RD.Parameters.landcover;
SOIL=RD.Parameters.soil;


Cohesion = zeros(size(LC));
FRI = zeros(size(LC));
N = zeros(size(LC));
K = zeros(size(LC));

for fc = 1:size(RD.seqmatrix,1) %Each row contains grids with equal flow accumulation
    col = 1; %Column counter for computational grid
    while (RD.seqmatrix(fc,col,1) > 0 && col < size(RD.seqmatrix,2))
        %Coordinates (array indices) of current grid
        i = RD.seqmatrix(fc,col,1);
        j = RD.seqmatrix(fc,col,2);
        lc = LC(i,j);
        soil = SOIL(i,j);
        
        switch lc
            case 0
                clc=0;
            case 1
                clc=10;
            case 2
                 clc=9;
            case 3
                 clc=8;
            case 4
                clc=7;
            case 5
                clc=6;
            case 6
                clc=5;
            case 7
                clc=4;
            case 8
                 clc=3;
            case 9
                clc=2;
            case 10
                clc=1;
            case 11
                clc=0;
            case 12
                clc=0;
            case 13
                clc=20;
            otherwise
                clc=0;
        end
        switch soil
            case 0
                c=NaN;fri=0;n=0;k=2.778*0*10e-6;
            case 1                                                          
                c=0;fri=40;n=0.43;k=2.446421e-5;%2.778*38.41*10e-6;
            case 2                                                         
                c=2000;fri=35;n=0.42;k=1.776770e-5;%2.778*10.87*10e-6;
            case 3                                                          % sandy loam, m/s
                c=4000;fri=30;n=0.453;k=1.022660e-5;%5.24*2.7778e-06;
            case 4
                c=6000;fri=29;n=0.46;k=2.501101e-6;%2.778*3.96*10e-6;
            case 5
                c=8000;fri=27;n=0.52;k=4.532431e-6;%2.778*8.59*10e-6;
            case 6                                                          % loam, m/s
                c=7000;fri=20;n=0.463;k=6.593731e-6;%1.97*2.7778e-06;
            case 7
                c=27000;fri=21;n=0.39;k=1.435262e-6;%2.778*2.4*10e-6;
            case 8
                c=28000;fri=21;n=0.48;k= 2.717260e-6;%2.778*4.57*10e-6;
            case 9
                c=29000;fri=19;n=0.46;k=4.314507e-6;%2.778*1.77*10e-6;
            case 10
                c=30000;fri=17;n=0.41;k=1.055191e-6;%2.778*1.19*10e-6;
            case 11
                c=35000;fri=17;n=0.49;k=1.307770e-6;%2.778*2.95*10e-6;
            case 12
                c=40000;fri=17;n=0.47;k=2.357930e-6;%2.778*3.18*10e-6;
                %case 13
                %    c=10000;fri=10;
                %case 14
                %    c=0;fri=1;
                %case 15
                %    c=100000;fri=40;
            otherwise
                c=25000;fri=30;n=0.4;k=2.778*2*10e-6;
        end
        
        Cohesion(i,j) = c+clc;
        FRI(i,j) = fri;
        N(i,j) = n;
        K(i,j) = k;
        
        col = col+1;
    end
end

end
