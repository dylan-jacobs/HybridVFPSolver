% ------- QCM ------- %
function [f_inf, n_inf, uz_inf, T_inf] = QCM(rho0, Jz0, kappa0, R, rvals, zvals)

    dr = rvals(2) - rvals(1);
    dz = zvals(2) - zvals(1);

    tol = 1e-15;

    Mk = [pi^(3/2); 0; 3]; % init guess for moments at equilibrium
    n_k = Mk(1);
    u_para_k = Mk(2);
    T_k = Mk(3);

    exp_r = exp(-(rvals.^2)/(2*R*T_k));
    exp_z = exp(-((zvals - u_para_k).^2)/(2*R*T_k));

    Rk_n = rho0   - sum(exp_r .* rvals)*(n_k/((2*pi*R*T_k)^(3/2)))*sum(exp_z)*2*pi*dr*dz;
    Rk_u = Jz0    - sum(exp_r .* rvals)*(n_k/((2*pi*R*T_k)^(3/2)))*sum(exp_z .* zvals)*2*pi*dr*dz;
    Rk_T = kappa0 - (sum(exp_r .* (rvals.^2) .* rvals)*(n_k/((2*pi*R*T_k)^(3/2)))*sum(exp_z) + sum(exp_r .* rvals)*(n_k/((2*pi*R*T_k)^(3/2)))*sum(exp_z .* (zvals.^2)))*pi*dr*dz;

    Rk = [Rk_n; Rk_u; Rk_T];  

    iters = 0;

    while norm(Rk, 2) > tol

        exp_r = exp(-(rvals.^2)/(2*R*T_k));
        exp_z = exp(-((zvals - u_para_k).^2)/(2*R*T_k));

        Jk_nn = -sum(exp_r .* rvals)*(1/((2*pi*R*T_k)^(3/2)))*sum(exp_z)*2*pi*dr*dz;
        Jk_nu = -sum(exp_r .* rvals)*(n_k/((2*pi*R*T_k)^(3/2)))*sum(exp_z .* ((zvals - u_para_k)/(R*T_k)))*2*pi*dr*dz;
        Jk_nT = -(sum(exp_r .* rvals)*(n_k/((2*pi*R)^(3/2)))*(-1.5*T_k^(-5/2))*sum(exp_z) + sum(exp_r .* (rvals.^2) .* rvals)*(n_k/((2*R*T_k^(7/2))*(2*pi*R)^(3/2)))*sum(exp_z) + sum(exp_r .* rvals)*(n_k/((2*R*T_k^(7/2))*(2*pi*R)^(3/2)))*sum(exp_z .* ((zvals - u_para_k).^2)))*2*pi*dr*dz;

        Jk_un = -sum(exp_r .* rvals)*(1/((2*pi*R*T_k)^(3/2)))*sum(exp_z .* zvals)*2*pi*dr*dz;
        Jk_uu = -sum(exp_r .* rvals)*(n_k/((2*pi*R*T_k)^(3/2)))*sum(exp_z .* ((zvals - u_para_k)/(R*T_k) .* zvals))*2*pi*dr*dz;
        Jk_uT = -(sum(exp_r .* rvals)*(n_k/((2*pi*R)^(3/2)))*(-1.5*T_k^(-5/2))*sum(exp_z .* zvals) + sum(exp_r .* (rvals.^2) .* rvals)*(n_k/((2*R*T_k^(7/2))*(2*pi*R)^(3/2)))*sum(exp_z .* zvals) + sum(exp_r .* rvals)*(n_k/((2*R*T_k^(7/2))*(2*pi*R)^(3/2)))*sum(exp_z .* ((zvals - u_para_k).^2) .* zvals))*2*pi*dr*dz;

        Jk_Tn = -(sum(exp_r .* (rvals.^2) .* rvals)*(1/((2*pi*R*T_k)^(3/2)))*sum(exp_z) + sum(exp_r .* rvals)*(1/((2*pi*R*T_k)^(3/2)))*sum(exp_z .* (zvals.^2)))*pi*dr*dz;
        Jk_Tu = -(sum(exp_r .* (rvals.^2) .* rvals)*(n_k/((2*pi*R*T_k)^(3/2)))*sum(exp_z .* ((zvals - u_para_k)/(R*T_k))) + sum(exp_r .* rvals)*(n_k/((2*pi*R*T_k)^(3/2)))*sum(exp_z .* (zvals.^2) .* ((zvals - u_para_k)/(R*T_k))))*pi*dr*dz;
        Jk_TT = -(sum(exp_r .* (rvals.^2) .* rvals)*(n_k/((2*pi*R)^(3/2)))*(-1.5*T_k^(-5/2))*sum(exp_z) + sum(exp_r .* (rvals.^4) .* rvals)*(n_k/((2*R*T_k^(7/2))*(2*pi*R)^(3/2)))*sum(exp_z) + sum(exp_r .* (rvals.^2) .* rvals)*(n_k/((2*R*T_k^(7/2))*(2*pi*R)^(3/2)))*sum(exp_z .* ((zvals - u_para_k).^2)))*pi*dr*dz    +    -(sum(exp_r .* rvals)*(n_k/((2*pi*R)^(3/2)))*(-1.5*T_k^(-5/2))*sum(exp_z .* (zvals.^2)) + sum(exp_r .* (rvals.^2) .* rvals)*(n_k/((2*R*T_k^(7/2))*(2*pi*R)^(3/2)))*sum(exp_z .* (zvals.^2)) + sum(exp_r .* rvals)*(n_k/((2*R*T_k^(7/2))*(2*pi*R)^(3/2)))*sum(exp_z .* (zvals.^2) .* ((zvals - u_para_k).^2)))*pi*dr*dz;

        J = [Jk_nn, Jk_nu, Jk_nT;
            Jk_un, Jk_uu, Jk_uT;
            Jk_Tn, Jk_Tu, Jk_TT];

        dMk = -J\Rk;
        Mk = Mk + dMk;

        n_k = Mk(1);
        u_para_k = Mk(2); 
        T_k = Mk(3);

        exp_r = exp(-(rvals.^2)/(2*R*T_k));
        exp_z = exp(-((zvals - u_para_k).^2)/(2*R*T_k));
        Rk_n = rho0   - sum(exp_r .* rvals)*(n_k/((2*pi*R*T_k)^(3/2)))*sum(exp_z)*2*pi*dr*dz;
        Rk_u = Jz0    - sum(exp_r .* rvals)*(n_k/((2*pi*R*T_k)^(3/2)))*sum(exp_z .* zvals)*2*pi*dr*dz;
        Rk_T = kappa0 - (sum(exp_r .* (rvals.^2) .* rvals)*(n_k/((2*pi*R*T_k)^(3/2)))*sum(exp_z) + sum(exp_r .* rvals)*(n_k/((2*pi*R*T_k)^(3/2)))*sum(exp_z .* (zvals.^2)))*pi*dr*dz;

        Rk = [Rk_n; Rk_u; Rk_T];

        iters = iters + 1;
        if iters > 100
            break
        end
    end
    f_inf = (n_k/((2*pi*R*T_k)^(3/2)))*exp(-((zvals'-u_para_k).^2 + rvals.^2)/(2*R*T_k));
    n_inf = n_k;
    uz_inf = u_para_k;
    T_inf = T_k;

end