%% Kinematic Identification for FANUC M-710iC robot without noises
clc;
clear;
close all;
%% Nominal value local POE
p0_init = [0; 0; 0; 0.5; 0; 0.565];             T0_init = Exp6(p0_init);
p1_init = [-0.5*pi; 0; 0; 0.15; 0; 0];          T1_init = Exp6(p1_init);
p2_init = [0; 0; 0; 0; -0.87; 0];               T2_init = Exp6(p2_init);
p3_init = [0; 0.5*pi; 0; 0.798; -0.17; 0.798];  T3_init = Exp6(p3_init);
p4_init = [0; -0.5*pi; 0; 0; 0; 0];             T4_init = Exp6(p4_init);
p5_init = [0; 0.5*pi; 0; 0; 0; 0];              T5_init = Exp6(p5_init);
p6_init = [0; 0; 0; 0; 0; 0.3];                 T6_init = Exp6(p6_init);
tau = [0; 0; 1; 0; 0; 0];   
%% Nominal value global POE
S1_init = Adjoint(T0_init) * tau;
S2_init = Adjoint(T0_init*T1_init) * tau;
S3_init = Adjoint(T0_init*T1_init*T2_init) * tau;
S4_init = Adjoint(T0_init*T1_init*T2_init*T3_init) * tau;
S5_init = Adjoint(T0_init*T1_init*T2_init*T3_init*T4_init) * tau;
S6_init = Adjoint(T0_init*T1_init*T2_init*T3_init*T4_init*T5_init) * tau;
Gamma_init = Log6(T0_init*T1_init*T2_init*T3_init*T4_init*T5_init*T6_init);
S_init = [S1_init; S2_init; S3_init; S4_init; S5_init; S6_init];
%% Actual value local POE
p0_d = [0; 0; 0; 0.5; 0; 0.565];                                        T0_d = Exp6(p0_d);
p1_d = [-1.563699; -0.007441; 0.009085; 0.144251; 0.008838; -0.002457]; T1_d = Exp6(p1_d);
p2_d = [-0.009101; -0.006606; 0.002227; 0.009006; -0.865002; -0.005787];T2_d = Exp6(p2_d);
p3_d = [-0.008989; 1.578209; 0.009233; 0.788113; -0.163947; 0.801418];  T3_d = Exp6(p3_d);
p4_d = [0.00457; -1.577157; 0.009714; -0.007373; 0.009921; -0.003729];  T4_d = Exp6(p4_d);
p5_d = [-0.009715; 1.568830; -0.007609; 0.003373; -0.008376; 0.004125]; T5_d = Exp6(p5_d);
p6_d = [0.009377; 0.006014; 0.009199; -0.004681; -0.006531; 0.304876];  T6_d = Exp6(p6_d);
%% Actual value
S1_d = Adjoint(T0_d) * tau;
S2_d = Adjoint(T0_d*T1_d) * tau;
S3_d = Adjoint(T0_d*T1_d*T2_d) * tau;
S4_d = Adjoint(T0_d*T1_d*T2_d*T3_d) * tau;
S5_d = Adjoint(T0_d*T1_d*T2_d*T3_d*T4_d) * tau;
S6_d = Adjoint(T0_d*T1_d*T2_d*T3_d*T4_d*T5_d) * tau;
Gamma_d = Log6(T0_d*T1_d*T2_d*T3_d*T4_d*T5_d*T6_d);
S_d = [S1_d; S2_d; S3_d; S4_d; S5_d; S6_d]; 
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
epsilon = 1e-10; sigma = 1e-6;
tau_max = 1; 
[S_DMH_1, Gamma_DMH_1, ~, error_DMH_1] = DMH_Method_KI_IK(S_init, Gamma_init, delta_q, q_bar_list, T_d_list, epsilon, tau_max, sigma, flag);
%% Plot
error_S_DMH_0 = norm(S_init - S_d);
error_Gamma_DMH_0 = norm(Gamma_init - Gamma_d);
error_S_DMH_1 = norm(S_DMH_1 - S_d);
error_Gamma_DMH_1 = norm(Gamma_DMH_1 - Gamma_d);
fprintf("DMH error of the zeroth initation：%d, %d\n", error_S_DMH_0, error_Gamma_DMH_0);
fprintf("DMH error of the first iteration：%d, %d\n", error_S_DMH_1, error_Gamma_DMH_1);
figure;
plot(0:1, [error_S_DMH_0, error_S_DMH_1]);
plot(0:1, [error_Gamma_DMH_0, error_Gamma_DMH_1]);


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