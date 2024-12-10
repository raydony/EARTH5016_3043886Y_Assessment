% set up index arrays for boundary conditions        
ix = [   1, 1:Nx, Nx  ];          % insulating x boundaries
iz = [   1, 1:Nz, Nz  ];          % isothermal start boundary

% set initial condition for temperature at cell centres

T = T0 + dTdz_boundaries(2)*Zc;  % initialise T array 

dt = CFL * (h/2)^2 / max(k0, [], "all"); % time step [s]

air = units == 9;% indices of air coordinates
k=100;
for t = 0:(dt*k):t_end
    % Calculate temperature changes using the diffusion function
    R1 = diffusion(T, k0, h, [1, 1:Nx, Nx], [1, 1:Nz, Nz], dTdz_boundaries(2));
    R2 = diffusion(T + R1 * dt / 2, k0, h, [1, 1:Nx, Nx], [1, 1:Nz, Nz], dTdz_boundaries(2));
    R3 = diffusion(T + R2 * dt / 2, k0, h, [1, 1:Nx, Nx], [1, 1:Nz, Nz], dTdz_boundaries(2));
    R4 = diffusion(T + R3 * dt, k0, h, [1, 1:Nx, Nx], [1, 1:Nz, Nz], dTdz_boundaries(2));
    
    % Update temperature with Runge-Kutta scheme
    T = T + dt * (R1 + 2 * R2 + 2 * R3 + R4) / 6 + Hr ./ (rho .* Cp);
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
