function [K] = K_h(tau_t,tau,h)
% Epanechnikov kernel
if abs((tau-tau_t)/h)<=1
    K=0.75*(1-power((tau-tau_t)/h,2));
else
    K=0;
end
end

