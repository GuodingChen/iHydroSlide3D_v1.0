clear;clear;
clc;
% First, you need specify the location of the iHydroSlide3D v1.0 executable
addpath('.\iHydroSlide3D_src'); 
% Second, you need provide a project file, which provides all necessary
% paths, files, initial conditions, parameter values and other information.
prediction = Main_program('ControlFile.ProjectWindows');



