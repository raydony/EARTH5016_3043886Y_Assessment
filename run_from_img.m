%***** RUN 2D MODEL FROM IMAGE ***********************************

% clear workspace
clear all; close all; %clc;

% load model setup from image, interpolate to target grid size
W       = 16e3;     % domain width (must correspond to width of image) [m]
Nx      = 200;      % target grid size z-direction
h       = W/Nx;     % grid spacing based on image width and target grid size
n_units = 9;        % number of rock units contained in image
[units,D,Nz] = ModelFromImage('section.tiff',n_units,W,Nx);


% material properties for each rock unit (update based on your calibration)
matprop = [
% unit  conductivity  density  heat capacity  heat production
   1        3.62        2679        1000        3.98
   2        4.02        2721        1000        4.27
   3        3.75        2692        1000        4.19
   4        3.41        2720        1000        4.20
   5        3.59        2676        1000        4.22
   6        2.97        2662        1000        5.76
   7        3.40        2726        1000        5.40
   8        2.90        2725        1000        5.48
   9        3.60        2701        1000        5.66];  % air/water

% get coefficient fields based on spatial distribution of rock units from image
% pay attention if any unit conversion is required!
rho    = reshape(matprop(units,3),Nz,Nx);
Cp     = reshape(matprop(units,4),Nz,Nx);
kT     = reshape(matprop(units,2),Nz,Nx);
Hr     = reshape(matprop(units,5),Nz,Nx);

% continue setting remaining model parameters, then call model routine