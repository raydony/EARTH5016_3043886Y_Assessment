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
% Function to calculate 2D advection rate
function dTdt = advection(f, u, w, dx, dz, ix, iz)
    % Split the velocities into positive and negative components
    upos = 0.5 * (u + abs(u));  % Positive velocity in x-direction
    uneg = 0.5 * (u - abs(u));  % Negative velocity in x-direction
    wpos = 0.5 * (w + abs(w));  % Positive velocity in z-direction
    wneg = 0.5 * (w - abs(w));  % Negative velocity in z-direction
    


    % Get values on stencil nodes for x-direction and z-direction
    fmmm = f(ix(1:end-6), iz);
    fmm  = f(ix(2:end-5), iz);
    fm   = f(ix(3:end-4), iz);
    fc   = f(ix(4:end-3), iz);
    fp   = f(ix(5:end-2), iz);
    fpp  = f(ix(6:end-1), iz);
    fppp = f(ix(7:end), iz);
   
    
    fmmmz = f(ix, iz(1:end-6));
    fmmz  = f(ix, iz(2:end-5));
    fmz   = f(ix, iz(3:end-4));
    fcz   = f(ix, iz(4:end-3));
    fpz   = f(ix, iz(5:end-2));
    fppz  = f(ix, iz(6:end-1));
    fpppz = f(ix, iz(7:end));

    % Calculate heat flux by advection in both directions (x and z)

    fppos_x = weno5poly(fmmm, fmm, fm, fc, fp); fpneg_x = weno5poly(fppp, fpp, fp, fc, fm);
    fmpos_x = weno5poly(fmmm, fmm, fm, fc, fp); fmneg_x = weno5poly(fpp, fp, fc, fm, fmm);

    fppos_z = weno5poly(fmmmz, fmmz, fmz, fcz, fpz); fpneg_z = weno5poly(fpppz, fppz, fpz, fcz, fmz);
    fmpos_z = weno5poly(fmmmz, fmmz, fmz, fcz, fpz); fmneg_z = weno5poly(fppz, fpz, fcz, fmz, fmmz);
    

    % Calculate flux balance for rate of change in both directions
    div_qpos_x = upos .* (fppos_x - fmpos_x) / dx;
    div_qneg_x = uneg .* (fpneg_x - fmneg_x) / dx;
    div_qpos_z = wpos .* (fppos_z - fmpos_z) / dz;
    div_qneg_z = wneg .* (fpneg_z - fmneg_z) / dz;

    div_q_x = div_qpos_x + div_qneg_x;
    div_q_z = div_qpos_z + div_qneg_z;

    dTdt = -(div_q_x + div_q_z);
end

% 5th order WENO polynomials function
function [fface] = weno5poly(fmm, fm, fc, fp, fpp)
    % 5th order polynomials
    p1 = (2*fmm - 7*fm + 11*fc) / 6;
    p2 = (-fm + 5*fc + 2*fp) / 6;
    p3 = (2*fc + 5*fp - fpp) / 6;

    % Smoothness measure
    b1 = 13/12 * (fmm - 2*fm + fc).^2 + 1/4 * (fmm - 4*fm + 3*fc).^2;
    b2 = 13/12 * (fm - 2*fc + fp).^2 + 1/4 * (fm - fp).^2;
    b3 = 13/12 * (fc - 2*fp + fpp).^2 + 1/4 * (3*fc - 4*fp + fpp).^2;

    % Weights
    g = [1/10, 6/10, 3/10];
    wp1 = g(1) ./ (b1.^2 + eps);
    wp2 = g(2) ./ (b2.^2 + eps);
    wp3 = g(3) ./ (b3.^2 + eps);

    % Normalize and compute flux
    fface = (wp1 .* p1 + wp2 .* p2 + wp3 .* p3) ./ (wp1 + wp2 + wp3);
end

