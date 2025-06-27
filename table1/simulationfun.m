function [RMSE] = simulationfun(T,h,Pmax)
rng(1);
d = 2;
K = 1000;

pre_error = zeros(K, d, 9);
RMSE = zeros(d, 9);

parfor k = 1:K
    disp(k)
    X = DGP1(T, h);
    phat = NIC(X(1:T,:), Pmax);
    Y_pred = zeros(d, Pmax);
    Omega_hat = zeros(d, d, Pmax);
    eta_hat = zeros(T-Pmax, d, Pmax);
    AIC = zeros(1, Pmax);
    BIC = zeros(1, Pmax);
    HQ = zeros(1, Pmax);
    ww_sAIC = zeros(1, Pmax);
    ww_sBIC = zeros(1, Pmax);
    ww_sHQ = zeros(1, Pmax);

    for p = 1:Pmax
        X_lag = zeros((T-p+1), d*p+1);
        for t = 1:(T-p+1)
            temp = 1;
            for r = 1:p
                temp = [temp, X(p+t-r,:)];
            end
            X_lag(t,:) = temp;
        end
        X_p = X((p+1):T,:); % regresan
        X_lag_h = X_lag(1:T-p,:);
        h_opt = 2.34*sqrt(1/12)*(T-p)^(-0.2);
        Y_fit_p = zeros((T-p), d);
        eta_hat_p = zeros((T-p), d);
        A_hat_p = zeros(d, d*p+1, T-p);

        for t = 1:(T-p)
            tau_t = t/(T-p);
            K_1 = K_weight(T-p, tau_t, h_opt); % kernel weighting matrix
            Z = Z_data(X_lag_h, tau_t, h_opt);
            A_hat_p(:,:,t) = ([eye(d*p+1) zeros(d*p+1)]*((Z'*K_1*Z)\(Z'*K_1*X_p)))';
            Y_fit_p(t,:) = X_lag_h(t,:)*A_hat_p(:,:,t)';
            eta_hat_p(t,:) = X_p(t,:) - X_lag_h(t,:) * A_hat_p(:,:,t)'; % estimates of reduced form residuals
        end

        y_pred_tem = zeros(1, h*d);
        for i = 1:h
            if i <= p
                y_lagged = [X_lag(end,1), y_pred_tem(:,(end-d*(i-1)+1):end), X_lag(end,2:(1+d*p-d*i+2))];
            else
                y_lagged = [1, y_pred_tem(:, (end-d*(i-1)+1):(end-d*(i-p-1)))];
            end
            y_pred_tem(:,(end-d*i+1):(end-d*i+d)) = y_lagged * A_hat_p(:,:,T-p)';
        end
        
        Y_pred(:,p) = y_pred_tem(1,1:d)';
        tau_t_omega = 1;
        K_Ome = K_weight(T-p, tau_t_omega, h_opt);
        Omega_hat(:,:,p) = eta_hat_p(:,:)' * K_Ome * eta_hat_p(:,:) / (sum(diag(K_Ome)));
        eta_hat(:,:,p) = eta_hat_p((1+Pmax-p):end,:);
    end

    PVEC = 1:Pmax;
    h_opt_ma = 2.34*sqrt(1/12)*(T-Pmax)^(-0.2);
    QQ = zeros(Pmax, Pmax);
    K_ma = K_weight((T-Pmax), 1, h_opt_ma);

    for j = 1:(T-Pmax)
        QQ = QQ + K_ma(j, j) * squeeze(eta_hat(j,:,:))' * inv(squeeze(Omega_hat(:,:,Pmax))) * squeeze(eta_hat(j,:,:));
    end

    QQ = QQ + QQ';
    B = 2 * log((T-Pmax) * h_opt_ma) * d^2 * PVEC';
    ww = quadprog(QQ, B, zeros(1, Pmax), 0, ones(1, Pmax), 1, zeros(Pmax, 1), ones(Pmax, 1));
    ww = ww(ww > 0);
    ww = ww / sum(ww);

    for p = 1:Pmax
        Omega_residual = eta_hat(:,:,p)' * eta_hat(:,:,p) / (T-Pmax);
        AIC(1, p) = log(det(squeeze(Omega_residual))) + 2 * p * d^2 / (T-Pmax);
        BIC(1, p) = log(det(squeeze(Omega_residual))) + log(T-Pmax) * p * d^2 / (T-Pmax);
        HQ(1, p) = log(det(squeeze(Omega_residual))) + 2 * log(log(T-Pmax)) * p * d^2 / (T-Pmax);
    end

    [~, loc_AIC] = min(AIC);
    [~, loc_BIC] = min(BIC);
    [~, loc_HQ] = min(HQ);

    ww_sAIC(1,:) = exp(-AIC(1,:) / 2);
    ww_sBIC(1,:) = exp(-BIC(1,:) / 2);
    ww_sHQ(1,:) = exp(-HQ(1,:) / 2);
    ww_sAIC(1,:) = ww_sAIC(1,:) / sum(ww_sAIC(1,:));
    ww_sBIC(1,:) = ww_sBIC(1,:) / sum(ww_sBIC(1,:));
    ww_sHQ(1,:) = ww_sHQ(1,:) / sum(ww_sHQ(1,:));

    Y_true = X(end,:);
    Y_pred_ma = Y_pred * ww;
    Y_pred_NIC = Y_pred(:, phat);
    Y_pred_AIC = Y_pred(:, loc_AIC);
    Y_pred_BIC = Y_pred(:, loc_BIC);
    Y_pred_HQ = Y_pred(:, loc_HQ);
    Y_pred_sAIC = Y_pred * ww_sAIC';
    Y_pred_sBIC = Y_pred * ww_sBIC';
    Y_pred_sHQ = Y_pred * ww_sHQ';
    Y_pred_SA = Y_pred * ones(Pmax,1)/Pmax;

    pre_error_ma(k,:) = Y_true - Y_pred_ma';
    pre_error_NIC(k,:) = Y_true - Y_pred_NIC';
    pre_error_AIC(k,:) = Y_true - Y_pred_AIC';
    pre_error_BIC(k,:) = Y_true - Y_pred_BIC';
    pre_error_HQ(k,:) = Y_true - Y_pred_HQ';
    pre_error_sAIC(k,:) = Y_true - Y_pred_sAIC';
    pre_error_sBIC(k,:) = Y_true - Y_pred_sBIC';
    pre_error_sHQ(k,:) = Y_true - Y_pred_sHQ';
    pre_error_SA(k,:) = Y_true - Y_pred_SA';
end



pre_error(:,:,1) = pre_error_ma;
pre_error(:,:,2) = pre_error_NIC;
pre_error(:,:,3) = pre_error_AIC;
pre_error(:,:,4) = pre_error_BIC;
pre_error(:,:,5) = pre_error_HQ;
pre_error(:,:,6) = pre_error_sAIC;
pre_error(:,:,7) = pre_error_sBIC;
pre_error(:,:,8) = pre_error_sHQ;
pre_error(:,:,9) = pre_error_SA;

for i = 1:9
    for j = 1:d
        RMSE(j, i) = sqrt(pre_error(:, j, i)' * pre_error(:, j, i) / K);
    end
end
end