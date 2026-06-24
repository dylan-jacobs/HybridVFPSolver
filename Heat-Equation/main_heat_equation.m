%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ----- Heat Equation Temporal Accuracy Test ----- %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear variables; close all;

% initial rank
r0 = 10;

% time-stepping method: 1=B.Euler, 2=DIRK2, 3=DIRK3, 4=DIRK4
methods = ['1', '2', '3'];
tolerance = 1e-6;

lambdavals = (1:1:100)'; % for temporal accuracy test
lambdavals = 1;        % for all other tests
errors = zeros(numel(lambdavals), 4);

% mesh parameters
rmin = 0; rmax = 1;
zmin = 0; zmax = 1;
Nr = 40; Nz = 80;  % temporal accuracy test
Nr = 80; Nz = 160; % conservation, rank tests
[Rmat, Zmat, dr, dz] = GetRZ(rmin, rmax, zmin, zmax, Nr, Nz);
rvals = Rmat(:, 1);
zvals = Zmat(1, :)';

% final time
tf = 1; % temporal accuracy test
tf = 5; % conservation, rank tests

% thermal conductivity coefficient
alpha = 0.1;

Drr = gallery('tridiag', Nr, rvals(1:end-1)+(dr/2), -2*rvals, rvals(1:end-1)+(dr/2));
Drr(1, 1) = -(rvals(1) + (dr/2));
Drr(end, end-1) = ((1/3) * (rvals(end) + (dr/2))) + (rvals(end) - (dr/2));
Drr(end, end) =  ((-3) * (rvals(end) + (dr/2))) - (rvals(end) - (dr/2));
Drr = alpha*((1/(dr^2)) * (diag(1./rvals) * Drr));

Dzz = alpha*((1/dz^2)*gallery('tridiag', Nz, 1, -2, 1)); % centered nodes

for k = 1:numel(lambdavals)
    dt = lambdavals(k)/((1/dr) + (1/dz));
    tvals = (0:dt:tf)';
    if tvals(end) ~= tf
        tvals = [tvals; tf];
    end
    Nt = numel(tvals);

    % store rank, mass, momentum, energy, l1 decay, etc...
    mass = zeros(Nt, numel(methods));
    Jzvals = zeros(Nt, numel(methods));
    E = zeros(Nt, numel(methods));
    ranks = zeros(Nt, numel(methods));

    for methodIndex = 1:numel(methods)
        method = methods(methodIndex);
        disp(['Method: ', method, ', dt=', num2str(dt, 3), ', ', num2str(k), '/', num2str(numel(lambdavals))]);
            
        % % initial conditions
        j01 = 2.40482555769577; % first root of bessel function
        f0 = @(r, z) (besselj(0, (j01/(rmax-rmin))*r)) .* (sin((2*pi/(zmax-zmin))*z));
        f_exact = @(r, z, t) (exp(-(alpha^1)*t*(((j01/(rmax-rmin))^2) + (2*pi/(zmax-zmin))^2))*f0(r, z));
        f_exact = f_exact(Rmat, Zmat, tf);
        
        f = f0(Rmat, Zmat);

        rhoM = sum(sum(f.*Rmat))*2*pi*dr*dz;
        JzM  = sum(sum(f.*Rmat.*Zmat))*2*pi*dr*dz;
        kappaM   = sum(sum(f.*Rmat.*((Rmat.^2 + Zmat.^2)/2)))*2*pi*dr*dz;
        
        % init bases
        [Vr, S, Vz] = svd2(f, rvals);
        r0 = min(r0, size(Vr, 2));
        Vr = Vr(:, 1:r0); S = S(1:r0, 1:r0); Vz = Vz(:, 1:r0);
               
        mass(1, methodIndex) = 2*pi*dr*dz*sum(sum(Rmat .* f));
        Jzvals(1, methodIndex) = 2*pi*dr*dz*sum(sum(f .* Rmat .* Zmat));
        E(1, methodIndex) = 2*pi*dr*dz*sum(sum(f .* Rmat .* ((Rmat.^2 + Zmat.^2)/2)));
        ranks(1, methodIndex) = r0;
        
        % time-stepping loop
        for n = 2:Nt
            tval = tvals(n);
            dt = tval - tvals(n-1);
            switch(method)
                case '1'
                    [Vr, S, Vz, rank] = BackwardEulerTimestep(Vr, S, Vz, dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM);
                case '2'
                    [Vr, S, Vz, rank] = DIRK2Timestep(Vr, S, Vz, dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM);
                case '3'
                    [Vr, S, Vz, rank] = DIRK3Timestep(Vr, S, Vz, dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM);
            end
            
            f = Vr*S*Vz';
            mass(n, methodIndex) = 2*pi*dr*dz*sum(sum(Rmat .* f));    
            Jzvals(n, methodIndex) = 2*pi*dr*dz*sum(sum(f .* Rmat .* Zmat));
            E(n, methodIndex) = pi*dr*dz*sum(sum(f .* (Rmat.^2 + Zmat.^2) .* Rmat));
            ranks(n, methodIndex) = rank;
        end
        
        errors(k, methodIndex) = 2*pi*dr*dz*sum((sum(Rmat .* abs(f - f_exact)))); % L1 error
    end
end

%%
% 1. Final solution
figure(1); clf; surf(Rmat, Zmat, f);
colorbar; shading interp;
legend(sprintf('N_r = %s', num2str(Nr, 3)), 'Location','northwest');
xlabel('V_r'); ylabel('V_z'); zlabel('U'); 

% 2. Exact solution
figure(2); clf; surf(Rmat, Zmat, f_exact);
colorbar; shading interp;
xlabel('V_r'); ylabel('V_z'); zlabel('f(V_r, V_z, t)'); title([sprintf('f_{exact} at time t=%s', num2str(tf, 4))]);

% % 3. Temporal error plot
% c_blue   = [0.1216 0.4667 0.7059];
% c_orange = [1.0000 0.4980 0.0549];
% c_green  = [0.1725 0.6275 0.1725];
% c_red    = [0.8392 0.1529 0.1569];
% c_purple = [0.5804 0.4039 0.7412];
% c_brown  = [0.5490 0.3373 0.2941];
% c_pink   = [0.8902 0.4667 0.7608];
% c_gray   = [0.4980 0.4980 0.4980];
% c_black  = [0.0000 0.0000 0.0000];
% 
% dtvals = lambdavals./((1/dr) + (1/dz));
% figure(3); clf; 
% begin_cutoff_1 = ceil(0.1*numel(dtvals)); end_cutoff_1 = ceil(0.4*numel(dtvals));
% begin_cutoff_2 = ceil(0.25*numel(dtvals)); end_cutoff_2 = ceil(0.7*numel(dtvals));
% begin_cutoff_3 = ceil(0.45*numel(dtvals)); end_cutoff_3 = ceil(0.95*numel(dtvals));
% % begin_cutoff_4 = ceil(0.45*numel(dtvals)); end_cutoff_4 = ceil(1*numel(dtvals));
% 
% loglog(dtvals, errors(:, 1, 1), 'Color', c_blue,  'LineWidth', 1.5); hold on; % B. Euler
% loglog(dtvals, errors(:, 2, 1), 'Color', c_orange,  'LineWidth', 1.5); % DIRK2
% loglog(dtvals, errors(:, 3, 1), 'Color', c_green,  'LineWidth', 1.5); % DIRK3
% % loglog(dtvals, errors(:, 4, 1), 'Color', c_black,  'LineWidth', 1.5); % DIRK4
% 
% loglog(dtvals(begin_cutoff_1:end_cutoff_1), 0.1*dtvals(begin_cutoff_1:end_cutoff_1), '--', 'Color', c_blue, 'LineWidth', 1); % Order 1
% loglog(dtvals(begin_cutoff_2:end_cutoff_2), 0.045*dtvals(begin_cutoff_2:end_cutoff_2).^2, '--', 'Color', c_orange, 'LineWidth', 1); % Order 2
% loglog(dtvals(begin_cutoff_3:end_cutoff_3), 0.02*dtvals(begin_cutoff_3:end_cutoff_3).^3, '--', 'Color', c_green, 'LineWidth', 1); % Order 3
% % loglog(dtvals(begin_cutoff_4:end_cutoff_4), 0.035*dtvals(begin_cutoff_4:end_cutoff_4).^4, '--', 'Color', c_black, 'LineWidth', 1); % Order 4
% xlabel('\Deltat'); ylabel('L^1 Error');
% legend('Backward Euler', 'DIRK2', 'DIRK3', 'Order 1', 'Order 2', 'Order 3', 'location','northwest');
% ylim([9e-5, 7e-1]);
% fontsize(18,"points");
% set(gcf,'Units','pixels','Position',[100 100 800 500])
% saveas(gcf, './Plots/heat_eqn_temporal_error_single_spatial.fig');
% exportgraphics(gcf,'./Plots/heat_eqn_temporal_error.pdf','ContentType','vector');

% 4. Mass conservation
figure(4); clf; 
plot(tvals, mass, 'LineWidth', 1.5);
xlabel('t'); ylabel('Mass'); % title('Relative mass of numerical solution');
legend('Backward Euler', 'DIRK2', 'DIRK3');
fontsize(18,"points");
set(gcf,'Units','pixels','Position',[100 100 800 500])
saveas(gcf, './Plots/heat_eqn_mass.fig');
exportgraphics(gcf,'./Plots/heat_eqn_mass.pdf','ContentType','vector');

% 5. Rank plot
figure(5); clf;
plot(tvals, ranks, 'LineWidth', 1.5);
xlabel('t'); ylabel('Rank');
legend('Backward Euler', 'DIRK2', 'DIRK3');
fontsize(18,"points");
set(gcf,'Units','pixels','Position',[100 100 800 500])
saveas(gcf, './Plots/heat_eqn_rank_plot.fig');
exportgraphics(gcf,'./Plots/heat_eqn_rank_plot.pdf','ContentType','vector');















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ---- HELPER FUNCTIONS ---- %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Vr, S, Vz, rank] = BackwardEulerTimestep(Vr0, S0, Vz0, dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM)

    Nr = numel(rvals);
    Nz = numel(zvals);

    Vr0_star = Vr0;
    Vz0_star = Vz0;

    K0 = Vr0*S0;
    L0 = Vz0*S0';

    K1 = sylvester(eye(Nr) - (dt*Drr), -dt*(Dzz*Vz0_star)'*Vz0_star, K0);
    L1 = sylvester(eye(Nz) - (dt*Dzz), -dt*(Drr*Vr0_star)'*(rvals.*Vr0_star), L0);

    [Vr1_ddagger, ~] = qr2(K1, rvals);
    [Vz1_ddagger, ~] = qr(L1, 0);

    [Vr1_hat, Vz1_hat] = reduced_augmentation([Vr1_ddagger, Vr0], [Vz1_ddagger, Vz0], rvals);

    S1_hat = sylvester((speye(size(Vr1_hat, 2)) - (dt*((rvals .* Vr1_hat)')*(Drr*Vr1_hat))), -dt*(Dzz*Vz1_hat)'*Vz1_hat, ((rvals .* Vr1_hat)'*Vr0)*S0*((Vz0')*Vz1_hat));
    [Vr, S, Vz, rank] = LoMaC_mass_only(Vr1_hat, S1_hat, Vz1_hat, rvals, zvals, tolerance, rhoM);
end

function [Vr, S, Vz, rank] = DIRK2Timestep(Vr0, S0, Vz0, dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM)

    Nr = numel(rvals);
    Nz = numel(zvals);

    nu = 1-(sqrt(2)/2);

    % Stage 1: Backward Euler
    [Vr1, S1, Vz1, ~] = BackwardEulerTimestep(Vr0, S0, Vz0, nu*dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM);

    W0_Vr = [Vr0, Drr*Vr1, Vr1];
    W0_S  = blkdiag(S0, (1-nu)*dt*S1, (1-nu)*dt*S1);
    W0_Vz = [Vz0, Vz1, Dzz*Vz1];

    % Reduced Augmentation
    % Predict V_dagger using B. Euler for second stage
    [Vr1_dagger, ~, Vz1_dagger, ~] = BackwardEulerTimestep(Vr0, S0, Vz0, dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM);
    [Vr_star, Vz_star] = reduced_augmentation([Vr1_dagger, Vr1, Vr0], [Vz1_dagger, Vz1, Vz0], rvals);
   
    % K/L-Step
    W0_K1 = W0_Vr*W0_S*(W0_Vz' * Vz_star);
    W0_L1 = (((rvals .* Vr_star)' * W0_Vr)*W0_S*(W0_Vz'))';

    K1 = sylvester(eye(Nr) - (nu*dt*Drr), -nu*dt*(Dzz*Vz_star)'*Vz_star, W0_K1);
    L1 = sylvester(eye(Nz) - (nu*dt*Dzz), -nu*dt*(Drr*Vr_star)'*(rvals .* Vr_star), W0_L1);

    % Get bases
    [Vr_ddagger, ~] = qr2(K1, rvals); [Vz_ddagger, ~] = qr(L1, 0);

    % Reduced Augmentation
    [Vr1_hat, Vz1_hat] = reduced_augmentation([Vr_ddagger, Vr1, Vr0], [Vz_ddagger, Vz1, Vz0], rvals);

    % S-Step
    W0_S1 = ((rvals .* Vr1_hat)' * W0_Vr)*W0_S*(W0_Vz'*Vz1_hat);

    S1_hat = sylvester(eye(size(Vr1_hat, 2)) - (nu*dt*((rvals .* Vr1_hat)')*Drr*Vr1_hat), -nu*dt*(Dzz*Vz1_hat)'*Vz1_hat, W0_S1);
    [Vr, S, Vz, rank] = LoMaC_mass_only(Vr1_hat, S1_hat, Vz1_hat, rvals, zvals, tolerance, rhoM);
end

function [Vr3, S3, Vz3, rank] = DIRK3Timestep(Vr0, S0, Vz0, dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM)

    % RK butcher table values
    nu = 0.435866521508459;
    beta1 = -(3/2)*(nu^2) + (4*nu) - (1/4);
    beta2 = (3/2)*(nu^2) - (5*nu) + (5/4);
    
    % Stage 1: Backward Euler
    [Vr1, S1, Vz1, ~] = BackwardEulerTimestep(Vr0, S0, Vz0, nu*dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM);
    [Vr_dagger1, ~, Vz_dagger1, ~] = BackwardEulerTimestep(Vr0, S0, Vz0, ((1+nu)/2)*dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM);

    W1_Vr = [Vr0, Drr*Vr1, Vr1];
    W1_S  = blkdiag(S0, (dt*(1-nu)/2)*S1, (dt*(1-nu)/2)*S1);
    W1_Vz = [Vz0, Vz1, Dzz*Vz1];

    % Reduced Augmentation
    [Vr_star1, Vz_star1] = reduced_augmentation([Vr_dagger1, Vr1, Vr0], [Vz_dagger1, Vz1, Vz0], rvals);

    % Stage 2:  
    % K/L-Step
    W1_K2 = W1_Vr*W1_S*(W1_Vz' * Vz_star1);
    W1_L2 = (((rvals .* Vr_star1)' * W1_Vr)*W1_S*(W1_Vz'))';

    K2 = sylvester(eye(size(Drr)) - (nu*dt*Drr), -nu*dt*(Dzz*Vz_star1)'*Vz_star1, W1_K2);
    L2 = sylvester(eye(size(Dzz)) - (nu*dt*Dzz), -nu*dt*(Drr*(rvals .* Vr_star1))'*Vr_star1, W1_L2);

    % Get bases
    [Vr_ddagger2, ~] = qr2(K2, rvals); [Vz_ddagger2, ~] = qr(L2, 0);

    % Reduced Augmentation
    [Vr2_hat, Vz2_hat] = reduced_augmentation([Vr_ddagger2, Vr1, Vr0], [Vz_ddagger2, Vz1, Vz0], rvals);

    % S-Step
    W1_S2 = ((rvals .* Vr2_hat)' * W1_Vr)*W1_S*(W1_Vz'*Vz2_hat);

    S2_hat = sylvester(eye(size(Vr2_hat, 2)) - (nu*dt*(rvals .* Vr2_hat)'*Drr*Vr2_hat), -nu*dt*(Dzz*Vz2_hat)'*Vz2_hat, W1_S2);
    [Vr2, S2, Vz2, ~] = truncate_svd(Vr2_hat, S2_hat, Vz2_hat, tolerance);

    % Stage 3:
    % Predict V_dagger using B. Euler
    [Vr_dagger3, ~, Vz_dagger3, ~] = BackwardEulerTimestep(Vr0, S0, Vz0, dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM);

    W2_Vr = [Vr0, Drr*Vr1, Vr1, Drr*Vr2, Vr2];
    W2_S  = blkdiag(S0, beta1*dt*S1, beta1*dt*S1, beta2*dt*S2, beta2*dt*S2);
    W2_Vz = [Vz0, Vz1, Dzz*Vz1, Vz2, Dzz*Vz2];

    % Reduced augmentation
    [Vr_star3, Vz_star3] = reduced_augmentation([Vr_dagger3, Vr2, Vr1, Vr0], [Vz_dagger3, Vz2, Vz1, Vz0], rvals);

    % K/L-Step
    W2_K3 = W2_Vr*W2_S*(W2_Vz' * Vz_star3);
    W2_L3 = (((rvals .* Vr_star3)' * W2_Vr)*W2_S*(W2_Vz'))';
    
    K3 = sylvester(eye(size(Drr)) - (nu*dt*Drr), -nu*dt*(Dzz*Vz_star3)'*Vz_star3, W2_K3);
    L3 = sylvester(eye(size(Dzz)) - (nu*dt*Dzz), -nu*dt*(Drr*Vr_star3)'*(rvals .* Vr_star3), W2_L3);

    % Get bases
    [Vr_ddagger3, ~] = qr2(K3, rvals); [Vz_ddagger3, ~] = qr(L3, 0);

    % Reduced Augmentation
    [Vr3_hat, Vz3_hat] = reduced_augmentation([Vr_ddagger3, Vr2, Vr1, Vr0], [Vz_ddagger3, Vz2, Vz1, Vz0], rvals);

    % S-Step
    W2_S3 = ((rvals .* Vr3_hat)' * W2_Vr)*W2_S*(W2_Vz'*Vz3_hat);

    S3_hat = sylvester(eye(size(Vr3_hat, 2)) - (nu*dt*((rvals .* Vr3_hat)')*Drr*Vr3_hat), -nu*dt*(Dzz*Vz3_hat)'*Vz3_hat, W2_S3);
    [Vr3, S3, Vz3, rank] = LoMaC_mass_only(Vr3_hat, S3_hat, Vz3_hat, rvals, zvals, tolerance, rhoM);
end


function [Vr4, S4, Vz4, rank] = DIRK4Timestep(Vr0, S0, Vz0, dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM)

    % RK butcher table values
    a_kk = 0.5;
    
    % Stage 1: Backward Euler
    [Vr1, S1, Vz1, ~] = BackwardEulerTimestep(Vr0, S0, Vz0, a_kk*dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM);
    [Vr_dagger1, ~, Vz_dagger1, ~] = BackwardEulerTimestep(Vr0, S0, Vz0, 0.25*dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM);

    Y1 = (((Drr*Vr1*S1*(Vz1')) + (Vr1*S1*((Dzz*Vz1)'))));
    W1 = (Vr0*S0*(Vz0')) + (-(1/4)*dt*Y1);

    % Reduced Augmentation
    [Vr_star1, Vz_star1] = reduced_augmentation([Vr_dagger1, Vr1, Vr0], [Vz_dagger1, Vz1, Vz0], rvals);

    % Stage 2:  
    % K/L-Step
    K2 = sylvester(eye(size(Drr)) - (a_kk*dt*Drr), -a_kk*dt*(Dzz*Vz_star1)'*Vz_star1, W1*Vz_star1);
    L2 = sylvester(eye(size(Dzz)) - (a_kk*dt*Dzz), -a_kk*dt*(Drr*(rvals .* Vr_star1))'*Vr_star1, (W1')*(rvals .* Vr_star1));

    % Get bases
    [Vr_ddagger2, ~] = qr2(K2, rvals); [Vz_ddagger2, ~] = qr(L2, 0);

    % Reduced Augmentation
    [Vr2_hat, Vz2_hat] = reduced_augmentation([Vr_ddagger2, Vr1, Vr0], [Vz_ddagger2, Vz1, Vz0], rvals);

    % S-Step
    S2_hat = sylvester(eye(size(Vr2_hat, 2)) - (a_kk*dt*(rvals .* Vr2_hat)'*Drr*Vr2_hat), -a_kk*dt*(Dzz*Vz2_hat)'*Vz2_hat, ((rvals .* Vr2_hat)')*W1*Vz2_hat);
    [Vr2, S2, Vz2, ~] = truncate_svd(Vr2_hat, S2_hat, Vz2_hat, tolerance);

    % Stage 3:
    % Predict V_dagger using B. Euler
    [Vr_dagger3, ~, Vz_dagger3, ~] = BackwardEulerTimestep(Vr0, S0, Vz0, 1.5*dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM);
    Y2 = (((Drr*Vr2*S2*(Vz2')) + (Vr2*S2*((Dzz*Vz2)'))));
    W2 = (Vr0*S0*(Vz0')) + (-1*dt*Y1) + (2*dt*Y2);
      
    % Reduced augmentation
    [Vr_star3, Vz_star3] = reduced_augmentation([Vr_dagger3, Vr2, Vr1, Vr0], [Vz_dagger3, Vz2, Vz1, Vz0], rvals);

    % K/L-Step
    K3 = sylvester(eye(size(Drr)) - (a_kk*dt*Drr), -a_kk*dt*(Dzz*Vz_star3)'*Vz_star3, W2*Vz_star3);
    L3 = sylvester(eye(size(Dzz)) - (a_kk*dt*Dzz), -a_kk*dt*(Drr*Vr_star3)'*(rvals .* Vr_star3), (W2')*(rvals .* Vr_star3));

    % Get bases
    [Vr_ddagger3, ~] = qr2(K3, rvals); [Vz_ddagger3, ~] = qr(L3, 0);

    % Reduced Augmentation
    [Vr3_hat, Vz3_hat] = reduced_augmentation([Vr_ddagger3, Vr2, Vr1, Vr0], [Vz_ddagger3, Vz2, Vz1, Vz0], rvals);

    % S-Step
    S3_hat = sylvester(eye(size(Vr3_hat, 2)) - (a_kk*dt*((rvals .* Vr3_hat)')*Drr*Vr3_hat), -a_kk*dt*(Dzz*Vz3_hat)'*Vz3_hat, ((rvals .* Vr3_hat)')*W2*Vz3_hat);
    [Vr3, S3, Vz3, ~] = truncate_svd(Vr3_hat, S3_hat, Vz3_hat, tolerance);

    % Stage 4:
    % Predict V_dagger using B. Euler
    [Vr_dagger4, ~, Vz_dagger4, ~] = BackwardEulerTimestep(Vr0, S0, Vz0, dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM);
    Y3 = (((Drr*Vr3*S3*(Vz3')) + (Vr3*S3*((Dzz*Vz3)'))));
    W3 = (Vr0*S0*(Vz0')) + (-(1/12)*dt*Y1) + ((2/3)*dt*Y2) + (-(1/12)*dt*Y3);
      
    % Reduced augmentation
    [Vr_star4, Vz_star4] = reduced_augmentation([Vr_dagger4, Vr3, Vr2, Vr1, Vr0], [Vz_dagger4, Vz3, Vz2, Vz1, Vz0], rvals);

    % K/L-Step
    K4 = sylvester(eye(size(Drr)) - (a_kk*dt*Drr), -a_kk*dt*(Dzz*Vz_star4)'*Vz_star4, W3*Vz_star4);
    L4 = sylvester(eye(size(Dzz)) - (a_kk*dt*Dzz), -a_kk*dt*(Drr*Vr_star4)'*(rvals .* Vr_star4), (W3')*(rvals .* Vr_star4));

    % Get bases
    [Vr_ddagger4, ~] = qr2(K4, rvals); [Vz_ddagger4, ~] = qr(L4, 0);

    % Reduced Augmentation
    [Vr4_hat, Vz4_hat] = reduced_augmentation([Vr_ddagger4, Vr3, Vr2, Vr1, Vr0], [Vz_ddagger4, Vz3, Vz2, Vz1, Vz0], rvals);

    % S-Step
    S4_hat = sylvester(eye(size(Vr4_hat, 2)) - (a_kk*dt*((rvals .* Vr4_hat)')*Drr*Vr4_hat), -a_kk*dt*(Dzz*Vz4_hat)'*Vz4_hat, ((rvals .* Vr4_hat)')*W3*Vz4_hat);
    [Vr4, S4, Vz4, rank] = LoMaC_mass_only(Vr4_hat, S4_hat, Vz4_hat, rvals, zvals, tolerance, rhoM);
end



function [Vr, S, Vz, rank] = truncate_svd(Vr, S, Vz, tolerance)
    [U, Sigma, V] = svd(S, 0);
    rank = find(diag(Sigma) > tolerance, 1, 'last');
    if (sum(rank) == 0)
        rank = 1;
    end
    Vr = Vr*U(:, 1:rank);
    S = Sigma(1:rank, 1:rank);
    Vz = Vz*V(:, 1:rank);
end

function [Vr, Vz] = reduced_augmentation(Vr_aug, Vz_aug, rvals)
    tolerance = 1e-12;
    [Qr, Rr] = qr2(Vr_aug, rvals);
    [Qz, Rz] = qr(Vz_aug, 0);
    [Ur, Sr, ~] = svd(Rr, 0);
    [Uz, Sz, ~] = svd(Rz, 0);
    rr = find(diag(Sr) > tolerance, 1, 'last');
    rz = find(diag(Sz) > tolerance, 1, 'last');
    rank = max(rr, rz);
    rank = min(rank, min(size(Ur, 2), size(Uz, 2)));
    Vr = Qr*Ur(:, 1:rank);
    Vz = Qz*Uz(:, 1:rank);
end

function [Q, R] = qr2(X, rvals)
    [Q, R] = qr(sqrt(rvals) .* X, 0);
    Q = Q ./ sqrt(rvals);
end

function [U, S, V] = svd2(X, rvals)
    [U, S, V] = svd(sqrt(rvals) .* X, 0);
    U = U./sqrt(rvals);
end

function [Rmat, Zmat, dr, dz] = GetRZ(vmin, vmax, zmin, zmax, Nv, Nz)
    rvals = linspace(vmin, vmax, Nv+1)';
    zvals = linspace(zmin, zmax, Nz+1)';
    dr = rvals(2) - rvals(1);
    dz = zvals(2) - zvals(1);
    rmid = rvals(1:end-1) + (dr/2);
    zmid = zvals(1:end-1) + (dz/2);
    [Rmat, Zmat] = meshgrid(rmid, zmid);
    Rmat = Rmat';
    Zmat = Zmat';
end



% ------- LoMaC Truncation -------
function [Vr, S, Vz, rank] = LoMaC_mass_only(Vr, S, Vz, rvals, zvals, tolerance, rhoM)
    % LoMaC Truncates given maxwellian (assumed Low-Rank) to given tolerance while conserving
    % macroscopic quantities (in this case, only mass for the heat equation).

    Nr = numel(rvals); Nz = numel(zvals);
    dr = rvals(2) - rvals(1);
    dz = zvals(2) - zvals(1);

    % Step 1: Integrate to calculate macro quantities
    p = 2*pi*dr*dz*sum(Vr .* rvals)*S*(sum(Vz)');

    % Step 2: Scale by maxwellian to ensure inner product is well defined
    % (f -> 0 as v -> infinity)
    dropoff = 10;
    wr = exp(-dropoff*(rvals.^2));
    wz = exp(-dropoff*(zvals.^2));    

    % Step 3: Orthogonal projection
    w_norm_1_squared = 2*pi*dr*dz*sum(rvals .* wr)*sum(wz);
    
    f1_proj_S_mtx11 = (p / w_norm_1_squared);

    proj_basis_r = wr.*[ones(Nr, 1)];
    proj_basis_z = wz.*[ones(Nz, 1)];
    f1_proj_S_mtx   = f1_proj_S_mtx11;

    % f2 = f - f1 (do it via SVD)
    f2_U = [Vr, proj_basis_r];
    f2_S = blkdiag(S, -f1_proj_S_mtx);
    f2_V = [Vz, proj_basis_z];

    % QR factorize
    [f2_Vr, f2_S, f2_Vz, ~] = truncate(f2_U, f2_S, f2_V, rvals, tolerance);

    % compute Pn(Te(f)) to ensure moments are kept
    trun_f2_p = 2*pi*dr*dz*(sum(f2_Vr .* rvals)*f2_S*(sum(f2_Vz)'));
    trun_f2_proj_S_mtx11 = (trun_f2_p / w_norm_1_squared);
    trun_f2_proj_S_mtx = trun_f2_proj_S_mtx11;

    % compute fM
    fM_proj_S_mtx11 = (rhoM ./ w_norm_1_squared);
    fM_proj_S_mtx = fM_proj_S_mtx11;

    f_mass_S = fM_proj_S_mtx - trun_f2_proj_S_mtx;

    [Vr, S, Vz, rank] = truncate([proj_basis_r, f2_Vr], blkdiag(f_mass_S, f2_S), [proj_basis_z, f2_Vz], rvals, 1e-14);
end

function [Vr, S, Vz, rank] = truncate(Vr_aug, S_aug, Vz_aug, rvals, tolerance)
    [Qr, Rr] = qr2(Vr_aug, rvals); [Qz, Rz] = qr(Vz_aug, 0);
    [U, Sigma, V] = svd(Rr*S_aug*Rz', 0); 
    rank = find(diag(Sigma) > tolerance, 1, 'last');
    if numel(rank) == 0
        rank = 1;
    end
    Vr = Qr*U(:, 1:rank);
    S = Sigma(1:rank, 1:rank);
    Vz = Qz*V(:, 1:rank);
end



