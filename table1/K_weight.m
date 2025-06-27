function [ K ] = K_weight(T,tau,h)
%return the T*T weighting matrix
tau_t=[1:T]/T;

K=zeros(T,T);             % kernel weighting matrix
for i=1:T 
    K(i,i)=K_h(tau_t(i),tau,h);
end
end

