%% Kinematic Identification for HP20D test without noises
%% License
% Copyright (c) 2025 Yuhan Chen (chenyuhan19930920@163.com)
% Licensed under the MIT License.
% See LICENSE file in the project root for full license information.
clc;
clear;
close all;
%% Nominal value
S1_init = [0; 0; 1; 0; 0; 0];
S2_init = [0; 1; 0; -505; 0; 150];
S3_init = [0; 1; 0; -1175; 0; 150];
S4_init = [1; 0; 0; 0; 1315; 0];
S5_init = [0; 1; 0; -1315; 0; 945];
S6_init = [1; 0; 0; 0; 1315; 0];
S_init = [S1_init; S2_init; S3_init; S4_init; S5_init; S6_init];
M_init = [0 -1 0 1070; 0 0 -1 0; 1 0 0 1415; 0 0 0 1];
Gamma_init = Log6(M_init);
%% Actual value
S1_d = [-0.00019; -0.00084; 0.999; 4.986; 2.460; 0.0030];
S2_d = [-0.00769; 0.999; 0.00084; -498.742; -3.964; 153.841];
S3_d = [-0.00769; 0.999; 0.00084; -1173.08; -9.161; 167.198];
S4_d = [0.999; 0.0077; -0.0368; -10.1473; 1322.304; 1.694];
S5_d = [-0.0065; 0.999; 0.0318; -1287.084; -38.923; 958.045];
S6_d = [0.999; 0.0058; 0.0232; -7.351; 1265.375; -0.135];
S_d = [S1_d; S2_d; S3_d; S4_d; S5_d; S6_d];
M_d = [-0.0230 -0.9997 0.0065 1087.27; -0.0320 -0.0058 -0.9995 13.013; 0.9992 -0.0232 -0.0318 1399.27; 0 0 0 1];
Gamma_d = Log6(M_d);
%% IK: DMH
% Preparation
rng(1);
flag = [1; 1; 0];
n = 6; m = 20;
delta_q = zeros(n,1);
q_bar_list = 2*pi*rand(n, m) - pi;
T_d_list = zeros(4, 4, m);
for i = 1 : m
    T_d_list(:,:,i) = ForwardKinematic(S_d, Gamma_d, delta_q, q_bar_list(:, i));
end
% DMH for KI
epsilon = 1e-15; sigma = 1e-6;
tau_max = 1; 
[S_DMH_1, Gamma_DMH_1, ~, error_DMH_1] = DMH_Method_KI_IK(S_init, Gamma_init, delta_q, q_bar_list, T_d_list, epsilon, tau_max, sigma, flag);
tau_max = 2; 
[S_DMH_2, Gamma_DMH_2, ~, error_DMH_2] = DMH_Method_KI_IK(S_init, Gamma_init, delta_q, q_bar_list, T_d_list, epsilon, tau_max, sigma, flag);
tau_max = 3; 
[S_DMH_3, Gamma_DMH_3, ~, error_DMH_3] = DMH_Method_KI_IK(S_init, Gamma_init, delta_q, q_bar_list, T_d_list, epsilon, tau_max, sigma, flag);
%% Veriﬁcation
num = 500;
T_d_list_verif = zeros(4, 4, num);
T_tau_0_verif = zeros(4, 4, num);
T_tau_1_verif = zeros(4, 4, num);
T_tau_2_verif = zeros(4, 4, num);
T_tau_3_verif = zeros(4, 4, num);
q_bar_list_verif = 2*pi*rand(n, num) - pi;
for i = 1 : num
    T_d_list_verif(:, :, i) = ForwardKinematic(S_d, Gamma_d, delta_q, q_bar_list_verif(:, i));
    T_tau_0_verif(:, :, i) = ForwardKinematic(S_init, Gamma_init, delta_q, q_bar_list_verif(:, i));
    T_tau_1_verif(:, :, i) = ForwardKinematic(S_DMH_1, Gamma_DMH_1, delta_q, q_bar_list_verif(:, i));
    T_tau_2_verif(:, :, i) = ForwardKinematic(S_DMH_2, Gamma_DMH_2, delta_q, q_bar_list_verif(:, i));
    T_tau_3_verif(:, :, i) = ForwardKinematic(S_DMH_3, Gamma_DMH_3, delta_q, q_bar_list_verif(:, i));
end
%% Data processing
error_0 = zeros(6, num);    error_1 = zeros(6, num);    error_2 = zeros(6, num);    error_3 = zeros(6, num);
error_Orient_0 = zeros(num, 1); error_Orient_1 = zeros(num, 1); error_Orient_2 = zeros(num, 1); error_Orient_3 = zeros(num, 1);
error_Posit_0 = zeros(num, 1); error_Posit_1 = zeros(num, 1); error_Posit_2 = zeros(num, 1); error_Posit_3 = zeros(num, 1);
for i = 1 : num
    error_0(:,i) = Log6(TransInv(T_tau_0_verif(:, :, i)) * T_d_list_verif(:, :, i));
    error_Orient_0(i) = norm(error_0(1:3,i));
    error_Posit_0(i) = norm(error_0(4:6,i));    
    error_1(:,i) = Log6(TransInv(T_tau_1_verif(:, :, i)) * T_d_list_verif(:, :, i));
    error_Orient_1(i) = norm(error_1(1:3,i));
    error_Posit_1(i) = norm(error_1(4:6,i));
    error_2(:,i) = Log6(TransInv(T_tau_2_verif(:, :, i)) * T_d_list_verif(:, :, i));
    error_Orient_2(i) = norm(error_2(1:3,i));
    error_Posit_2(i) = norm(error_2(4:6,i));   
    error_3(:,i) = Log6(TransInv(T_tau_3_verif(:, :, i)) * T_d_list_verif(:, :, i));
    error_Orient_3(i) = norm(error_3(1:3,i));
    error_Posit_3(i) = norm(error_3(4:6,i));     
end
error_Orient_max = [max(error_Orient_0); max(error_Orient_1); max(error_Orient_2); max(error_Orient_3)];
error_Posit_max = [max(error_Posit_0); max(error_Posit_1); max(error_Posit_2); max(error_Posit_3)];
error_Orient_mean = [mean(error_Orient_0); mean(error_Orient_1); mean(error_Orient_2); mean(error_Orient_3)];
error_Posit_mean = [mean(error_Posit_0); mean(error_Posit_1); mean(error_Posit_2); mean(error_Posit_3)];
%% Plot
figure; 
semilogy(0:length(error_Orient_max)-1, error_Orient_max);
hold on;
semilogy(0:length(error_Orient_mean)-1, error_Orient_mean);
xlabel("Number of iteration");
ylabel("Orientation error (rad)");
legend("Max DMH", "Mean DMH");
figure; 
semilogy(0:length(error_Posit_max)-1, error_Posit_max);
hold on;
semilogy(0:length(error_Posit_mean)-1, error_Posit_mean);
xlabel("Number of iteration");
ylabel("Position error (mm)");
legend("Max DMH", "Mean DMH");
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

function so3mat = VecToso3(omg)
    so3mat = [0, -omg(3), omg(2); omg(3), 0, -omg(1); -omg(2), omg(1), 0];
end

function omg = so3ToVec(so3mat)
    omg = [so3mat(3,2); so3mat(1,3); so3mat(2,1)];
end
