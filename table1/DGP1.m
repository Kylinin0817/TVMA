function [ x_t ] = DGP1(T,h)

tau = @(t) [t/T];


epsilon_t = mvnrnd(zeros(2, 1), eye(2), T+h);


pi_T = epsilon_t(:, 1) + randn(T+h, 1);

omega = @(tau) [1.5 + 0.2 * exp(0.5 - tau), 0; 
                0.1 * exp(0.5 - tau), (1.5 + 0.5 * (tau - 0.5))^2];


a = @(tau,t) [0.5 * sin(2 * pi_T(t,1) * tau), 0.5 * cos(2 * pi_T(t,1) * tau)]';

A1 = @(tau,t) [0.8 * exp(-0.5 + tau), 0.8 * (tau - 0.5)^3;
             0.8 * (tau - 0.5)^3, 0.8 + 0.3 * sin(pi_T(t,1) * tau)];
         
A2 = @(tau,t) [-0.2 * exp(-0.5 + tau), 0.8 * (tau - 0.5)^2;
             0.8 * (tau - 0.5)^2, -0.4 + 0.3 * cos(pi_T(t,1) * tau)];


x_t = zeros(T+h, 2);
eta_t = zeros(T+h, 2);


for t = 1:(T+h)
    if t > 2
        eta_t(t, :) = (omega(tau(t)) * epsilon_t(t, :)')';
        x_t(t, :) = (a(tau(t),t) + A1(tau(t),t) * x_t(t-1, :)' + A2(tau(t),t) * x_t(t-2, :)')' + eta_t(t, :);
    elseif t > 1
        eta_t(t, :) = (omega(tau(t)) * epsilon_t(t, :)')';
        x_t(t, :) = (a(tau(t),t) + A1(tau(t),t) * x_t(t-1, :)')' + eta_t(t, :);
    else
        eta_t(t, :) = (omega(tau(t)) * epsilon_t(t, :)')';
        x_t(t, :) = a(tau(t),t)' + eta_t(t, :);
    end
end
end
