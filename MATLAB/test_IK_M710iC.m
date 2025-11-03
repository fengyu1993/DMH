%% Inverse Kinematics for M710iC test
clc;
clear;
close all;
%% License
% Copyright (c) 2025 Yuhan Chen (chenyuhan19930920@163.com)
% Licensed under the MIT License.
% See LICENSE file in the project root for full license information.
%% Parameters
flag_robot = 1;
[S, Gamma] = get_M710iC_POE_Nominal();
epsilon = 1e-20;
n = length(S) / 6;
q_bar = zeros(n, 1); 
sigma = 1e-6;
tau_max = 30;
flag = [0 0 1];
noise_q = pi/4;
rng(1);
%% IK
% Desired pose and initial joints
q_d = (2*rand(n, 1) - 1) * pi;
T_d = ForwardKinematic(S, Gamma, q_d, q_bar);
delta_q_init = q_d + (2*rand(n, 1) - 1) * noise_q;
% IK
[~, ~, q_DMH, error_DMH] = DMH_Method_KI_IK(S, Gamma, delta_q_init, q_bar, T_d, epsilon, tau_max, sigma, flag);
[~, ~, q_DTH, error_DTH] = DTH_Method_KI_IK(S, Gamma, delta_q_init, q_bar, T_d, epsilon, tau_max, sigma, flag);
[~, ~, q_DNR, error_DNR] = DNR_Method_KI_IK(S, Gamma, delta_q_init, q_bar, T_d, epsilon, tau_max, sigma, flag);
%% Plot
figure;
semilogy(0:length(error_DMH)-1, error_DMH);
hold on;
semilogy(0:length(error_DTH)-1, error_DTH);
semilogy(0:length(error_DNR)-1, error_DNR);
legend('DMH', 'DTH', 'DNR');
xlabel("Number of iterations", 'Interpreter', 'latex');
ylabel("Pose error  $\Xi^T \Xi$", 'Interpreter', 'latex');





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S, Gamma] = get_M710iC_POE_Nominal()
    % Nominal value local POE
    p0_init = [0; 0; 0; 0.5; 0; 0.565];             T0_init = Exp6(p0_init);
    p1_init = [-0.5*pi; 0; 0; 0.15; 0; 0];          T1_init = Exp6(p1_init);
    p2_init = [0; 0; 0; 0; -0.87; 0];               T2_init = Exp6(p2_init);
    p3_init = [0; 0.5*pi; 0; 0.798; -0.17; 0.798];  T3_init = Exp6(p3_init);
    p4_init = [0; -0.5*pi; 0; 0; 0; 0];             T4_init = Exp6(p4_init);
    p5_init = [0; 0.5*pi; 0; 0; 0; 0];              T5_init = Exp6(p5_init);
    p6_init = [0; 0; 0; 0; 0; 0.3];                 T6_init = Exp6(p6_init);
    tau = [0; 0; 1; 0; 0; 0];   
    % Nominal value global POE
    S1_init = Adjoint(T0_init) * tau;
    S2_init = Adjoint(T0_init*T1_init) * tau;
    S3_init = Adjoint(T0_init*T1_init*T2_init) * tau;
    S4_init = Adjoint(T0_init*T1_init*T2_init*T3_init) * tau;
    S5_init = Adjoint(T0_init*T1_init*T2_init*T3_init*T4_init) * tau;
    S6_init = Adjoint(T0_init*T1_init*T2_init*T3_init*T4_init*T5_init) * tau;
    Gamma = Log6(T0_init*T1_init*T2_init*T3_init*T4_init*T5_init*T6_init);
    S = [S1_init; S2_init; S3_init; S4_init; S5_init; S6_init];
    
end
function T = ForwardKinematic(S, Gamma, delta_q, q_bar)
    n = length(q_bar);
    Slist = reshape(S, 6, n);
    T = Exp6(Gamma);
    for i = n : -1 : 1
        T = Exp6(Slist(:,i) * (delta_q(i) + q_bar(i))) * T;
    end
end
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
function judge = NearZero(near)
    judge = norm(near) < 1e-9;
end
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
function  [R, p] = TransToRp(T)
    R = T(1:3,1:3);
    p = T(1:3,4);
end
function dexp_w_inv = get_dexp_w_inv(lambda_w)
    norm_lambda_w = norm(lambda_w);
    if NearZero(norm_lambda_w)
        dexp_w_inv = eye(3);
    else
        norm_lambda_w_2 = norm_lambda_w * norm_lambda_w;
        s = sin(norm_lambda_w/2) / (norm_lambda_w/2);
        c = cos(norm_lambda_w/2);
        gamma = c / s;
        lambda_w_hat = VecToso3(lambda_w);
        dexp_w_inv = eye(3) - 1/2*lambda_w_hat + (1-gamma)/norm_lambda_w_2 * lambda_w_hat * lambda_w_hat;
    end
end
function invT = TransInv(T)
    [R, p] = TransToRp(T);
    invT = [R', -R' * p; 0, 0, 0, 1];
end
function AdT = Adjoint(T)
    [R, p] = TransToRp(T);
    AdT = [R, zeros(3); VecToso3(p) * R, R];
end
function so3mat = VecToso3(omg)
    so3mat = [0, -omg(3), omg(2); omg(3), 0, -omg(1); -omg(2), omg(1), 0];
end
function omg = so3ToVec(so3mat)
    omg = [so3mat(3,2); so3mat(1,3); so3mat(2,1)];
end
