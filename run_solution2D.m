% 1. Load the section image and extract rock units
section_image = imread('section.tiff');  % Read the section image
rock_units = extract_rock_units(imresize(section_image, 0.2));  % Downsample and extract rock units

% 2. Define material properties for each rock unit (as per given calibration)
%He1,He2,Ms,Cm,Gr,Si,Sa,Bg,Air/water
matprop = [
    % unit  conductivity  density  heat capacity  heat production
    1        3.6788      2697.6     1000          4.172
    2        3.2197      2703.5     1000          5.575
    3        0.919       1905.90    1000          1
    4        0.924       2083.14    1000          1
    5        1.68        2648       1000          1
    6        0.249       1916       1000          1
    7        0.37        1942.3     1000          1
    8        1           2700       1000          1
    9        1e-6        1          1000          0];  % Air/water

% 3. Initialize the model grid and spatial dimensions (following run_transect2D.m style)
W = 16e3;    % Domain width [m]
Nx = 200;    % Number of grid cells in x-direction
h = W / Nx;  % Grid spacing based on domain width
n_units = 9; % number of rock units contained in image

[units, D, Nz] = ModelFromImage('section.tiff', 9, W, Nx); % Call the function

x_cells = h/2:h:W-h/2;  % Cell center positions in x-direction
z_cells = h/2:h:D-h/2;  % Cell center positions in z-direction

% Create 2D grids for cell centers
[Xc, Zc] = meshgrid(x_cells, z_cells);

% 4. Assign material properties to each rock unit (keeping similar to transect2D.m)
[m, n] = size(units);  % Dimensions of the model grid
rho = matprop(units(:), 3);      % Density [kg/m^3]
Cp = matprop(units(:), 4);       % Heat capacity [J/(kg*K)]
kT = matprop(units(:), 2);    % Thermal conductivity [W/(m*K)]
Hr = matprop(units(:), 5);       % Heat production [W/m^3]

% Reshape properties to match the grid dimensions
rho = reshape(rho, m, n);
Cp = reshape(Cp, m, n);
kT = reshape(kT, m, n);
Hr = reshape(Hr, m, n);

% Calculate diffusivity (k0) based on properties
k0 = kT ./ (rho .* Cp);

% remaining model parameters
dTdz_top = 0;                      % flux at top
dTdz_bot = 30/1000;                % flux at bottom
Ttop = 0;       % (C)              % air/water temperature
Tbot = dTdz_bot*5000;              % find temperature at model base


yr    = 3600*24*365;  % seconds per year [s]
tend = 1e7*yr;

CFL   = 1/5;         % Time step limiter
nop   = 50;          % output figure produced every 'nop' steps


T0 = 10; % surface air temperature

dTdz_boundaries = [0, 35/1000];

function rock_units = extract_rock_units(section_image)
    % Extract rock units from grayscale section image
    rock_units = round(double(section_image) / 255 * 9);  % Map grayscale to 1-9 units
end


% *****  RUN MODEL
run('solution2D.m');