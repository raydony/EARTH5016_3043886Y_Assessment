% set up index arrays for boundary conditions        
ix = [   1, 1:Nx, Nx  ];          % insulating x boundaries
iz = [   1, 1:Nz, Nz  ];          % isothermal start boundary

% set initial condition for temperature at cell centres

T = T0 + dTdz_boundaries(2)*Zc;  % initialise T array 

dt = CFL * (h/2)^2 / max(k0, [], "all"); % time step [s]

air = units == 9;% indices of air coordinates
k = 100; 
t = 0;    
tau = 0;  

figure;
while t <= t_end
    % canculate the temp
    R1 = diffusion(T, k0, h, [1, 1:Nx, Nx], [1, 1:Nz, Nz], dTdz_boundaries(2));
    R2 = diffusion(T + R1 * dt / 2, k0, h, [1, 1:Nx, Nx], [1, 1:Nz, Nz], dTdz_boundaries(2));
    R3 = diffusion(T + R2 * dt / 2, k0, h, [1, 1:Nx, Nx], [1, 1:Nz, Nz], dTdz_boundaries(2));
    R4 = diffusion(T + R3 * dt, k0, h, [1, 1:Nx, Nx], [1, 1:Nz, Nz], dTdz_boundaries(2));

    % Runge-Kutta
    T = T + dt * (R1 + 2 * R2 + 2 * R3 + R4) / 6 + Hr ./ (rho .* Cp);

    % time and step
    t = t + dt * k;  % add dt*k every time
    tau = tau + 1;

   % Visualize the temperature distribution every 'nop' steps
    if mod(tau, nop) == 0
        clf;%clear the fig

        imagesc(x_cells, z_cells, T); % Plot temperature distribution
        colorbar;
               
        hold on;
        contour(x_cells, z_cells, T, [150 150], 'LineWidth', 2, 'LineColor', 'r');% Addition of isotherm at 150°C
        hold off;
        t_years = t / yr;  % Calculate current time in years

        title(['Temperature Distribution at t = ' num2str(t_years) ' years']);% Update title with years
        xlabel('Lateral Distance [m]');
        ylabel('Depth [m]');   
        clim([min(T(:)), max(T(:))]);  % Set the color bar's min and max range
        drawnow;  % Refresh the figure window
    end
end

% 6. Visualize the temperature field
imagesc(x_cells, z_cells, T);  % Display temperature field
colorbar;
title('Temperature Distribution Simulation Result');
xlabel('Lateral Distance [m]');
ylabel('Depth [m]');

% Save the final temperature field
save('thermal_distribution.mat', 'T');

% Auxiliary functions
function dTdt = diffusion(T, k0, h, ix, iz, base_flux)
    % Calculates temperature diffusion rates
    kz = (k0(iz(1:end-1), :) + k0(iz(2:end), :)) / 2;
    kx = (k0(:, ix(1:end-1)) + k0(:, ix(2:end))) / 2;

    qz = -kz .* diff(T(iz, :), 1, 1) / h;  % Heat flux in z-direction
    qx = -kx .* diff(T(:, ix), 1, 2) / h;  % Heat flux in x-direction
    qz(end, :) = -k0(end, :) * base_flux;

    dTdt = -(diff(qz, 1, 1) + diff(qx, 1, 2)) / h;
end
