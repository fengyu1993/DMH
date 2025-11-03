%% Kinematic Identification for SCARA test without noises
%% License
% Copyright (c) 2025 Yuhan Chen (chenyuhan19930920@163.com)
% Licensed under the MIT License.
% See LICENSE file in the project root for full license information.
clc;
clear;
close all;
%% Nominal value
S1_init = [0; 0; 1; 0; 0; 0];
S2_init = [0; 0; 1; 0; -0.25; 0];
S3_init = [0; 0; 0; 0; 0; -1];
S4_init = [0; 0; -1; 0; 0.47; 0];
S_init = [S1_init; S2_init; S3_init; S4_init];
%% Actual value
S1_d = [0.01999; 0; 0.9998; 0; 0.01303; 0];
S2_d = [0; 0.0004; 1; -0.0003; -0.25399; 0.000102];
S3_d = [0; 0; 0; 0.02; 0.0196; -0.99961];
S4_d = [0.04077; 0.03917; -0.9984; -0.02668; 0.504558; 0.018706];
S_d = [S1_d; S2_d; S3_d; S4_d];
Gamma = [0; 0; 0; 0; 0; 0];
%% IK: DMH
% Preparation
rng(1);
flag = [1; 0; 0];
n = 4; m = 6;
delta_q = zeros(n,1);
q_bar_list = 2*pi*rand(n, m) - pi;
T_d_list = zeros(4, 4, m);
for i = 1 : m
    T_d_list(:,:,i) = ForwardKinematic(S_d, Gamma, delta_q, q_bar_list(:, i));
end
% DMH for KI
epsilon = 1e-10; tau_max = 1; sigma = 1e-6;
[S_DMH_1, ~, ~, error_1] = DMH_Method_KI_IK(S_init, Gamma, delta_q, q_bar_list, T_d_list, epsilon, tau_max, sigma, flag);
%% Plot
error_DMH_1 = norm(S_DMH_1 - S_d);
fprintf("DMH error of the first iteration：%d\n", error_DMH_1);
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
    judge = norm(near) < 1e-7;
end
function so3mat = VecToso3(omg)
    so3mat = [0, -omg(3), omg(2); omg(3), 0, -omg(1); -omg(2), omg(1), 0];
end
function omg = so3ToVec(so3mat)
    omg = [so3mat(3,2); so3mat(1,3); so3mat(2,1)];
end
