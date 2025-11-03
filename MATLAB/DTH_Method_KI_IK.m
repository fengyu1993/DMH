function [S, Gamma, delta_q, error] = DTH_Method_KI_IK(S_init, Gamma_init, delta_q_init, q_bar_list, T_d_list, epsilon, tau_max, sigma, flag)
%% Damped Traditional Halley Method for solving Kinematic Identification and Inverse Kinematics 
%% Input
% S_init: Initialized S represented in Eq.(17).
% Gamma_init: Initialized exponential coordinates of the end-effector pose relative to the frame {s} when the robot is at its zero position.
% delta_q_init: Initialized delta_q represented in Eq.(18).
% q_bar_list: Nominal values of the measured m set of joint variables q_bar represented in Eq.(18). q_bar_list = [q_bar_1, q_bar_2, ..., q_bar_m].
% T_d_list: Measured poses of m set of end effectors. T_d_list in 4x4xm.
% epsilon: Convergence tolerance.
% tau_max: Maximum iteration count.
% sigma: Damping factor.
% flag: Selection of identification parameters. 
    % flag = [flag_1 flag_2 flag_3] 
    % flag_1: Determine whether to identify the parameter S
    % flag_2: Determine whether to identify the parameter Gamma
    % flag_3: Determine whether to identify the parameter delta_q
    % e.g. flag = [1 1 0] -- Kinematic Identification: S Gamma
    %      flag = [0 0 1] -- Inverse Kinematics: delta_q
%% Output
% S, Gamma, delta_q: Result after identification. 
% error: Iterative process error, Xi^T Xi
%% License
% Copyright (c) 2025 Yuhan Chen (chenyuhan19930920@163.com)
% Licensed under the MIT License.
% See LICENSE file in the project root for full license information.
%% Initialization
[n, m] = size(q_bar_list);
Phi = SGqToPhi(S_init, Gamma_init, delta_q_init, flag); % Eq.(21)
l_Phi = length(Phi);
l_Xi = 6 * m;
J_Xi = zeros(l_Xi, l_Phi);
temp_H = zeros(l_Xi, l_Phi);
H_Xi = zeros(l_Xi, l_Phi, l_Phi);
I = eye(l_Phi);
error = zeros(tau_max, 1);
%% Repeat
for tau = 1 : tau_max
    % Eq.(21)(24)
    [S_tau , Gamma_tau, delta_q_tau] = PhiToSGq(Phi, n, flag, S_init, Gamma_init, delta_q_init);
    %% Repeat
    for i = 1 : m
        % Eq.(16)
        T_e_i = ForwardKinematic(S_tau, Gamma_tau, delta_q_tau, q_bar_list(:,i));
        % Eq.(23)
        T_xi_i = TransInv(T_e_i) * T_d_list(:, :, i);
        % Eq.(61)
        xi_i = Log6(T_xi_i);
        % Eq.(53)
        J_Phi = get_Jacobi_Phi(S_tau, Gamma_tau, delta_q_tau, q_bar_list(:,i), flag);
        J_xi_i = -get_dexp_inv(xi_i) * J_Phi;
        % Eq.(55)
        d_Phi_tilde_list = eye(l_Phi);
        %% Repeat
        for k = 1 : l_Phi
            % Eq.(55)
            d_Phi_k = d_Phi_tilde_list(:, k);
            % Eq.(21)(24)
            [d_S_tau , d_Gamma_tau, d_delta_q_tau] = PhiToSGq(d_Phi_k, n, flag, zeros(size(S_init)), zeros(size(Gamma_init)), zeros(size(delta_q_init)));
            % Eq.(44)
            d_xi_i = J_xi_i * d_Phi_k;
            % Eq.(60)
            d_J_Phi = get_d_J_Phi(S_tau, Gamma_tau, delta_q_tau, q_bar_list(:,i), d_S_tau , d_Gamma_tau, d_delta_q_tau, flag);
            d_J_xi_i = -get_dexp_inv_dt(xi_i, d_xi_i) * J_Phi - get_dexp_inv(xi_i) * d_J_Phi;
            % Eq.(58)(59)
            H_Xi(6*i-5 : 6*i,:,k) = d_J_xi_i;
        end
        % Eq.(22)(25)
        Xi(6*i-5 : 6*i, 1) = xi_i;
        % Eq.(43)
        J_Xi(6*i-5 : 6*i, :) = J_xi_i;
    end
    error(tau) = Xi' * Xi;
    if error(tau) < epsilon
        break
    end
    % Eq.(38)
    Delta_Phi_DNR = -inv(J_Xi' * J_Xi + sigma*I) * J_Xi' * Xi;
    % Eq.(39)
    temp_H = temp_H * 0;
    for i = 1 : l_Phi
        temp_H = temp_H + Delta_Phi_DNR(i) * H_Xi(:,:,i);
    end
    J_Xi_DNR = J_Xi + 1/2*temp_H;
    Delta_Phi_DTH = -inv(J_Xi_DNR' * J_Xi_DNR + sigma*I) * J_Xi_DNR' * Xi;
    % Eq.(42)
    Phi = Phi + Delta_Phi_DTH;
end
%% Output
% Eq.(21)(24)
[S , Gamma, delta_q] = PhiToSGq(Phi, n, flag, S_init, Gamma_init, delta_q_init);
error = error(1 : tau);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function invT = TransInv(T)
    [R, p] = TransToRp(T);
    invT = [R', -R' * p; 0, 0, 0, 1];
end
function  [R, p] = TransToRp(T)
    R = T(1:3,1:3);
    p = T(1:3,4);
end
function judge = NearZero(near)
    judge = norm(near) < 1e-7;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Eq.(2)
function so3mat = VecToso3(omg)
    so3mat = [0, -omg(3), omg(2); omg(3), 0, -omg(1); -omg(2), omg(1), 0];
end
function omg = so3ToVec(so3mat)
    omg = [so3mat(3,2); so3mat(1,3); so3mat(2,1)];
end
%% Eq.(3) 
function T = Exp6(lambda)
    lambda_w = lambda(1:3);
    lambda_v = lambda(4:6);
    norm_lambda_w = norm(lambda_w);
    if NearZero(norm_lambda_w)
        T = [eye(3), lambda_v; 0, 0, 0, 1];
    else
        s = sin(norm_lambda_w/2) / (norm_lambda_w/2);
        c = cos(norm_lambda_w/2);
        alpha = s * c;
        beta = s * s;
        lambda_w_hat = VecToso3(lambda_w);
        R = eye(3) + alpha*lambda_w_hat + beta/2*lambda_w_hat*lambda_w_hat;
        p = get_dexp_w(lambda_w) * lambda_v;
        T = [R, p; 0, 0, 0, 1];
    end
end
%% Eq.(4)
function dexp_w = get_dexp_w(lambda_w)
    norm_lambda_w = norm(lambda_w);
    norm_lambda_w_2 = norm_lambda_w * norm_lambda_w;
    s = sin(norm_lambda_w/2) / (norm_lambda_w/2);
    c = cos(norm_lambda_w/2);
    alpha = s * c;
    beta = s * s;
    lambda_w_hat = VecToso3(lambda_w);
    if NearZero(norm_lambda_w)
        dexp_w = eye(3);
    else
        dexp_w = eye(3) + beta/2*lambda_w_hat + (1-alpha)/norm_lambda_w_2 * lambda_w_hat * lambda_w_hat;
    end
end
function dexp_w_inv = get_dexp_w_inv(lambda_w)
    norm_lambda_w = norm(lambda_w);
    if NearZero(norm_lambda_w)
        dexp_w_inv = eye(3);
    else
        norm_lambda_w_2 = norm_lambda_w * norm_lambda_w;
        s = sin(norm_lambda_w/2) / (norm_lambda_w/2);
        c = cos(norm_lambda_w/2);
        if s < 1e-9     
            fprintf("dexp_w is invertible.");   
        end
        gamma = c / s; 
        lambda_w_hat = VecToso3(lambda_w);
        dexp_w_inv = eye(3) - 1/2*lambda_w_hat + (1-gamma)/norm_lambda_w_2 * lambda_w_hat * lambda_w_hat;
    end
end
%% Eq.(7)
function dexp_lambda = get_dexp(lambda)
    lambda_w = lambda(1:3);
    lambda_v = lambda(4:6);
    dexp_w = get_dexp_w(lambda_w);
    C_w = get_C_lambda_w(lambda_w, lambda_v);
    dexp_lambda = [dexp_w, zeros(3, 3);...
                    C_w, dexp_w];
end
function dexp_lambda_inv = get_dexp_inv(lambda)
    lambda_w = lambda(1:3);
    lambda_v = lambda(4:6);
    dexp_w_inv = get_dexp_w_inv(lambda_w);
    D_w = get_D_lambda_w(lambda_w, lambda_v);
    dexp_lambda_inv = [dexp_w_inv, zeros(3, 3);...
                       D_w, dexp_w_inv];
end
function C_w = get_C_lambda_w(lambda_w, y)
    norm_lambda_w = norm(lambda_w);
    norm_lambda_w_2 = norm_lambda_w * norm_lambda_w;
    s = sin(norm_lambda_w/2) / (norm_lambda_w/2);
    c = cos(norm_lambda_w/2);
    alpha = s * c;
    beta = s * s;
    lambda_w_hat = VecToso3(lambda_w);
    lambda_y_hat = VecToso3(y);
    if norm_lambda_w == 0
        C_w = 1/2*lambda_y_hat;
    else
        Lambda_1 = (1-alpha)/norm_lambda_w_2;
        Lambda_2 = (alpha - beta)/norm_lambda_w_2;
        Lambda_3 = - 1/norm_lambda_w_2*(3*(1-alpha)/norm_lambda_w_2 - beta/2);
        C_w = beta/2*lambda_y_hat + Lambda_1*(lambda_y_hat*lambda_w_hat + lambda_w_hat*lambda_y_hat) ...
                + Lambda_2*lambda_w'*y*lambda_w_hat ...
                + Lambda_3*lambda_w'*y*lambda_w_hat*lambda_w_hat;
    end
end
function D_w = get_D_lambda_w(lambda_w, y)
    norm_lambda_w = norm(lambda_w);
    norm_lambda_w_2 = norm_lambda_w * norm_lambda_w;
    s = sin(norm_lambda_w/2) / (norm_lambda_w/2);
    c = cos(norm_lambda_w/2);
    beta = s * s;
    if s < 1e-9     
            fprintf("dexp is invertible.");   
    end
    gamma = c / s;
    lambda_w_hat = VecToso3(lambda_w);
    lambda_y_hat = VecToso3(y);
    if norm_lambda_w == 0
        D_w = -1/2*lambda_y_hat;
    else
        Lambda_4 = (1-gamma)/norm_lambda_w_2;
        Lambda_5 = 1/norm_lambda_w_2*(1/beta+gamma-2)/norm_lambda_w_2;
        D_w = -1/2*lambda_y_hat + Lambda_4*(lambda_w_hat*lambda_y_hat + lambda_y_hat*lambda_w_hat) ...
                + Lambda_5 *lambda_w'*y*lambda_w_hat*lambda_w_hat;
    end
end
%% Eq.(11)
function AdT = Adjoint(T)
    [R, p] = TransToRp(T);
    AdT = [R, zeros(3); VecToso3(p) * R, R];
end
%% Eq.(12)
function AdT_Inv = AdjointInv(T)
    [R, p] = TransToRp(T);
    AdT_Inv = [R', zeros(3); -R' * VecToso3(p), R'];
end
function adV = ad(V)
    omgmat = VecToso3(V(1:3));
    adV = [omgmat, zeros(3); VecToso3(V(4:6)), omgmat];
end
%% Eq.(16)
function T = ForwardKinematic(S, Gamma, delta_q, q_bar)
    n = length(q_bar);
    Slist = reshape(S, 6, n);
    T = Exp6(Gamma);
    for i = n : -1 : 1
        T = Exp6(Slist(:,i) * (delta_q(i) + q_bar(i))) * T;
    end
end
%% Eq.(21)(24)
function [S, Gamma, delta_q] = PhiToSGq(Phi, n, flag, S_init, Gamma_init, delta_q_init)
    if flag(1) == 1
        S = Phi(1 : 6*n);
        Phi = Phi(6*n+1 : end);
    else
        S = S_init;
    end
    if flag(2) == 1
        Gamma = Phi(1 : 6);
        Phi = Phi(7 : end);
    else
        Gamma = Gamma_init;
    end
    if flag(3) == 1
        delta_q = Phi(1 : n);
    else
        delta_q = delta_q_init;
    end    
end
function Phi = SGqToPhi(S, Gamma, delta_q, flag)
    Phi = [];
    if flag(1) == 1
        Phi = [Phi; S];
    end
    if flag(2) == 1
        Phi = [Phi; Gamma];
    end
    if flag(3) == 1
        Phi = [Phi; delta_q];
    end
end
%% Eq.(22)
function Xi = get_Xi(S, Gamma, delta_q, q_bar_list, T_d_list)
    [~, m] = size(q_bar_list);
    Xi = zeros(6*m, 1);
    for i = 1 : m
        % Eq.(16)
        T_e_DTH_i = ForwardKinematic(S, Gamma, delta_q, q_bar_list(:,i));
        % Eq.(23)
        T_xi_DTH_i = TransInv(T_e_DTH_i) * T_d_list(:, :, i);
        % Eq.(60)
        xi_i = Log6(T_xi_DTH_i);
        % Eq.(22)
        Xi(6*i-5 : 6*i, 1) = xi_i;
    end
end
%% Eq.(53)
function J_Phi = get_Jacobi_Phi(S, Gamma, delta_q, q_bar, flag)
    n = length(q_bar);
    Slist = reshape(S, 6, n);
    J_temp = zeros(6, 6, n);
    T_temp = Exp6(Gamma);
    for i = n : -1 : 1
        J_temp(:, :, i) = Adjoint(TransInv(T_temp)) * get_dexp(-Slist(:, i) * (q_bar(i) + delta_q(i)));
        T_temp = Exp6(Slist(:, i) * (q_bar(i) + delta_q(i))) * T_temp;
    end
    J_Phi = [];
    if flag(1) == 1
        J_Phi_S = [];
        for i = 1 : n
            J_Phi_S = [J_Phi_S, J_temp(:, :, i) * (q_bar(i) + delta_q(i))];
        end
        J_Phi = [J_Phi, J_Phi_S];
    end
    if flag(2) == 1
        J_Phi_Gamma = get_dexp(-Gamma);
        J_Phi = [J_Phi, J_Phi_Gamma];
    end
    if flag(3) == 1
        J_Phi_delta_q = [];
        for i = 1 : n
            J_Phi_delta_q = [J_Phi_delta_q, J_temp(:, :, i) * Slist(:, i)];
        end
        J_Phi = [J_Phi, J_Phi_delta_q];
    end
end
%% Eq.(60)
function dJ_Phi_dt = get_d_J_Phi(S, Gamma, delta_q, q_bar, d_S, d_Gamma, d_delta_q, flag)
    n = length(q_bar);
    Slist = reshape(S, 6, n);
    dSlist = reshape(d_S, 6, n);
    %% T_list
    Tlist = zeros(4, 4, n+1);
    for i = 1 : n
        Tlist(:, :, i) = Exp6(Slist(:, i) * (q_bar(i) + delta_q(i)));
    end
    Tlist(:, :, n+1) = Exp6(Gamma);
    %% Ti_n1_list
    Ti_n1_list = zeros(4, 4, n+1);
    Ti_n1_list(:, :, n+1) = Tlist(:, :, n+1);
    for i = n : -1 : 2
        Ti_n1_list(:, :, i) = Tlist(:, :, i) * Ti_n1_list(:, :, i+1);
    end
    %% Vlist
    Vlist = zeros(6, n+1);
    for k = 1 : n
        Vlist(:,k) = get_dexp(-Slist(:, k) * (q_bar(k) + delta_q(k))) ...
            * (dSlist(:,k) * (q_bar(k) + delta_q(k)) + Slist(:, k) * d_delta_q(k));
    end
    Vlist(:,n+1) = get_dexp(-Gamma) * d_Gamma;
    %% dt_Ad_Tlist_inv
    dt_Ad_Tlist_inv = zeros(6, 6, n+1);
    for k = 2 : n+1
        dt_Ad_Tlist_inv(:, :, k) = -ad(Vlist(:,k)) * AdjointInv(Tlist(:, :, k));
    end
    %% dt_Ad_Ti_n1_inv
    dt_Ad_Ti_n1_inv = zeros(6, 6, n+1);
    dt_Ad_Ti_n1_inv(:, :, n+1) = dt_Ad_Tlist_inv(:, :, n+1);
    T_temp = Tlist(:, :, n+1);
    for k = n : -1 : 2
        dt_Ad_Ti_n1_inv(:, :, k) = dt_Ad_Ti_n1_inv(:, :, k+1) * AdjointInv(Tlist(:, :, k)) ...
                            + AdjointInv(T_temp) * dt_Ad_Tlist_inv(:, :, k);
        T_temp = Tlist(:, :, k) * T_temp;
    end
    %% dexp_list
    dexp_list = zeros(6, 6, n);
    for k = 1 : n
        dexp_list(:,:,k) = get_dexp(-Slist(:, k) * (q_bar(k) + delta_q(k)));
    end
    %% dt_dexp_list
    dt_dexp_list = zeros(6, 6, n);
    for k = 1 : n
        dt_dexp_list(:, :, k) = get_dexp_dt(-Slist(:, k) * (q_bar(k) + delta_q(k)), -dSlist(:, k) * (q_bar(k) + delta_q(k)) - Slist(:, k)*d_delta_q(k));
    end
    %% dt_Js
    dJ_Phi_dt = [];
    if flag(1) == 1
        dt_Js = [];
        for k = 1 : n
            dt_Jsk = dt_Ad_Ti_n1_inv(:, :, k+1) * dexp_list(:,:,k) * (q_bar(k) + delta_q(k))  ...
                + AdjointInv(Ti_n1_list(:, :, k+1)) * dt_dexp_list(:, :, k) * (q_bar(k) + delta_q(k)) ...
                + AdjointInv(Ti_n1_list(:, :, k+1)) * dexp_list(:,:,k) * d_delta_q(k);
            dt_Js = [dt_Js, dt_Jsk];
        end
        dJ_Phi_dt = [dJ_Phi_dt, dt_Js];
    end
    %% dt_Jgamma
    if flag(2) == 1
        dt_Jgamma = get_dexp_dt(-Gamma, -d_Gamma);
        dJ_Phi_dt = [dJ_Phi_dt, dt_Jgamma];
    end
    %% dt_Jdq
    if flag(3) == 1
        dt_Jdq = [];
        for k = 1 : n
            dt_Jdqk = dt_Ad_Ti_n1_inv(:, :, k+1) * Slist(:, k)  ...
                + AdjointInv(Ti_n1_list(:, :, k+1)) * dSlist(:, k);
            dt_Jdq = [dt_Jdq, dt_Jdqk];
        end        
        dJ_Phi_dt = [dJ_Phi_dt, dt_Jdq];
    end
end
%% Eq.(61)
function lambda = Log6(T)
    [R, p] = TransToRp(T);
    if NearZero(norm(R - eye(3)))
        lambda = [0; 0; 0; T(1:3,4)];
    else
        if NearZero(trace(R) + 1)
            if ~NearZero(1 + R(3,3))
                lambda_w = (pi / sqrt(2 * (1 + R(3,3)))) * [R(1,3); R(2,3); 1 + R(3,3)];
            elseif ~NearZero(1 + R(2,2))
                lambda_w = (pi / sqrt(2 * (1 + R(2,2)))) * [R(1,2); 1 + R(2,2); R(3,2)];
            else
                lambda_w = (pi / sqrt(2 * (1 + R(1,1)))) * [1 + R(1,1); R(2,1); R(3,1)];
            end
            lambda_v = get_dexp_w_inv(lambda_w) * p;
        else
            acosinput = (trace(R) - 1) / 2;
            if acosinput > 1
                acosinput = 1;
            elseif acosinput < -1
                acosinput = -1;
            end
            theta = acos(acosinput);
            if theta == 0
                lambda_w = zeros(3,1);
            else
                lambda_w = so3ToVec(theta/(2*sin(theta)) * (R - R'));
            end
            lambda_v = get_dexp_w_inv(lambda_w) * p;
        end
        lambda = [lambda_w; lambda_v];
    end
end
%% Eq.(78)
function dexp_dt = get_dexp_dt(lambda, d_lambda)
    lambda_w = lambda(1:3);
    d_lambda_w = d_lambda(1:3);
    dexp_dt = [get_C_lambda_w(lambda_w, d_lambda_w), zeros(3, 3);...
               get_dC_dt(lambda, d_lambda), get_C_lambda_w(lambda_w, d_lambda_w)];
end
function dC_dt = get_dC_dt(lambda, d_lambda)
    lambda_w = lambda(1:3);
    lambda_v = lambda(4:6);
    d_lambda_w = d_lambda(1:3);
    d_lambda_v = d_lambda(4:6);
    %%
    norm_lambda_w = norm(lambda_w);
    if norm_lambda_w == 0
        lambda_v_hat = VecToso3(lambda_v);
        d_lambda_w_hat = VecToso3(d_lambda_w);
        d_lambda_v_hat = VecToso3(d_lambda_v);
        dC_dt = 1/2*d_lambda_v_hat + 1/6*(d_lambda_w_hat*lambda_v_hat + lambda_v_hat*d_lambda_w_hat);
    else
        norm_lambda_w_2 = norm_lambda_w * norm_lambda_w;
        s = sin(norm_lambda_w/2) / (norm_lambda_w/2);
        c = cos(norm_lambda_w/2);
        alpha = s * c;
        beta = s * s;
        lambda_w_hat = VecToso3(lambda_w);
        lambda_v_hat = VecToso3(lambda_v);
        d_lambda_w_hat = VecToso3(d_lambda_w);
        d_lambda_v_hat = VecToso3(d_lambda_v);
        %%
        Lambda_1 = (1-alpha)/norm_lambda_w_2;
        Lambda_2 = (alpha - beta)/norm_lambda_w_2;
        Lambda_3 = - 1/norm_lambda_w_2*(3*(1-alpha)/norm_lambda_w_2 - beta/2);
        %%
        tau = lambda_w' * lambda_v * (lambda_w' * d_lambda_w) / norm_lambda_w_2;
        detla_0 = d_lambda_w'*lambda_v + lambda_w'*d_lambda_v;
        delta_1 = lambda_w'*lambda_v*d_lambda_w_hat + lambda_w'*d_lambda_w*lambda_v_hat ...
                   + (detla_0 - 4*tau)*lambda_w_hat;
        delta_2 = lambda_w'*lambda_v*(lambda_w_hat*d_lambda_w_hat + d_lambda_w_hat*lambda_w_hat) ...
                   + lambda_w'*d_lambda_w*(lambda_w_hat*lambda_v_hat + lambda_v_hat*lambda_w_hat) ...
                   + (detla_0 - 5*tau)*lambda_w_hat*lambda_w_hat;
        %%
        dC_dt = beta/2 * (d_lambda_v_hat - tau*lambda_w_hat) ...
              + Lambda_1 * (d_lambda_v_hat*lambda_w_hat + lambda_w_hat*d_lambda_v_hat + lambda_v_hat*d_lambda_w_hat + d_lambda_w_hat*lambda_v_hat + tau*lambda_w_hat) ...
               + Lambda_2 * (delta_1 + tau*lambda_w_hat*lambda_w_hat) ...
                + Lambda_3 * delta_2;
    end
end
%% Eq.(79)
function dexp_inv_dt = get_dexp_inv_dt(lambda, d_lambda)
    lambda_w = lambda(1:3);
    d_lambda_w = d_lambda(1:3);
    %%
    dexp_inv_dt = [get_D_lambda_w(lambda_w, d_lambda_w), zeros(3, 3);...
               get_dD_dt(lambda, d_lambda), get_D_lambda_w(lambda_w, d_lambda_w)];
end
function dD_dt = get_dD_dt(lambda, d_lambda)
    lambda_w = lambda(1:3);
    lambda_v = lambda(4:6);
    d_lambda_w = d_lambda(1:3);
    d_lambda_v = d_lambda(4:6);
    %%
    norm_lambda_w = norm(lambda_w);
    if norm_lambda_w == 0
        lambda_v_hat = VecToso3(lambda_v);
        d_lambda_w_hat = VecToso3(d_lambda_w);
        d_lambda_v_hat = VecToso3(d_lambda_v);
        dD_dt = -1/2*d_lambda_v_hat + 1/12*(d_lambda_w_hat*lambda_v_hat + lambda_v_hat*d_lambda_w_hat);
    else
        norm_lambda_w_2 = norm_lambda_w * norm_lambda_w;
        s = sin(norm_lambda_w/2) / (norm_lambda_w/2);
        c = cos(norm_lambda_w/2);
        beta = s * s;
        if s < 1e-9     
            fprintf("dexp is invertible.");   
        end
        gamma = c / s;
        lambda_w_hat = VecToso3(lambda_w);
        lambda_v_hat = VecToso3(lambda_v);
        d_lambda_w_hat = VecToso3(d_lambda_w);
        d_lambda_v_hat = VecToso3(d_lambda_v);
        %%
        Lambda_4 = (1-gamma)/norm_lambda_w_2;
        Lambda_5 = 1/norm_lambda_w_2*(1/beta+gamma-2)/norm_lambda_w_2;
        %%
        tau = lambda_w' * lambda_v * (lambda_w' * d_lambda_w) / norm_lambda_w_2;
        detla_0 = d_lambda_w'*lambda_v + lambda_w'*d_lambda_v;
        delta_3 = lambda_w'*lambda_v*(lambda_w_hat*d_lambda_w_hat + d_lambda_w_hat*lambda_w_hat) ...
                   + lambda_w'*d_lambda_w*(lambda_w_hat*lambda_v_hat + lambda_v_hat*lambda_w_hat) ...
                   + (detla_0 - 3*tau)*lambda_w_hat*lambda_w_hat;
        %%
        dD_dt = -1/2*d_lambda_v_hat + 2/norm_lambda_w_2 * (1 - gamma/beta)/norm_lambda_w_2*tau*lambda_w_hat*lambda_w_hat ...
                   + Lambda_4*(d_lambda_v_hat*lambda_w_hat + lambda_w_hat*d_lambda_v_hat + lambda_v_hat*d_lambda_w_hat + d_lambda_w_hat*lambda_v_hat) ...
                   + Lambda_5*delta_3;
    end
end
%% Eq.(82)
function d_AdT_inv_dt = get_d_AdT_inv_dt(T, V)
    d_AdT_inv_dt = -ad(V) * AdjointInv(T);
end



