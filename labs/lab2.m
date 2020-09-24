function lab2()
    x = importdata('data1.txt');
    N = length(x);
    
    mu = find_expected_value(x);
    fprintf("mu = %.4f \n", mu);
    
    S2 = find_dispersion(x);
    fprintf("S2 = %.4f \n", S2);

    gamma = 0.9;  
    
    [lm, hm] = find_exp_v_interval(gamma, S2, mu, N);
    fprintf("MX interval for gamma=%.1f: (%.4f ; %.4f)\n", gamma, lm, hm);
    
    [lsigma, hsigma] = find_dispersion_interval(gamma, S2, N); 
    fprintf("DX interval for gamma=%.1f: (%.4f ; %.4f)\n", gamma, lsigma, hsigma);
    
    figure(1);
    hold on;
    expected_value_graph(x, N, gamma);
    
    figure(2);
    hold on;
    dispersion_graph(x, N, gamma); 
end

% оценка мат ожидания
function mu = find_expected_value(x)
    n = length(x);
    mu = sum(x)/n;
end
% оценка дисперсии
function S2 = find_dispersion(x)
    n = length(x);
    mu = find_expected_value(x);
    if n > 1
        S2 = sum((x - mu).^2) / (n-1);
    else
        S2 = 0;
    end
end
% вычисление границ гамма доверительного интервала для MX
function [lsigma, hsigma] = find_exp_v_interval(gamma, S2, mu, n)
    quantile = tinv((1 + gamma)/2, n-1);
    lsigma = mu - quantile * sqrt(S2) / sqrt(n);
    hsigma = mu + quantile * sqrt(S2) / sqrt(n);
end
% вычисление границ гамма доверительного интервала для DX
function [ld, hd] = find_dispersion_interval(gamma, S2, n)   
    quantilel = chi2inv((1 + gamma)/2, n-1);
    quantileh = chi2inv((1 - gamma)/2, n-1); 
    ld = S2*(n-1)/quantilel;
    hd = S2*(n-1)/quantileh;    
end
% построение прямой и графиков для мат ожидания
function expected_value_graph(x, N, gamma)
    mu = zeros(N,1);
    S2 = zeros(N,1);
    lmu = zeros(N,1);
    hmu = zeros(N,1);
    start = 1;
    for n = start:N
        newx = x(1:n);
        mu(n) = find_expected_value(newx);
        S2(n) = find_dispersion(newx);
        [lmu(n), hmu(n)] = find_exp_v_interval(gamma, S2(n), mu(n), n);
    end

    mu_line = zeros(N,1);
    mu_line(1:N) = mu(N);
    
    plot((start:N),mu_line(start:N), 'r');
    plot((start:N),lmu(start:N), 'b');
    plot((start:N),hmu(start:N), 'g');
    plot((start:N),mu(start:N), 'm');
    grid on;
    xlabel('n');
    ylabel('\mu');
    legend('\mu\^(x_N)','_{--}\mu^(x_n)', '^{--}\mu^(x_n)', '\mu\^(x_n)');
end
% построение прямой и графиков для дисперсии
function dispersion_graph(x, N, gamma)
    S2 = zeros(N,1);
    lsigma = zeros(N,1);
    hsigma = zeros(N,1);
    start = 1;
    for n = start:N
        newx = x(1:n);
        S2(n) = find_dispersion(newx);
        [lsigma(n), hsigma(n)] = find_dispersion_interval(gamma, S2(n), n);
    end

    S2_line = zeros(N,1);
    S2_line(start:N) = S2(N);
    
    plot((start:N), S2_line(start:N), 'r');
    plot((start:N), lsigma(start:N), 'b');
    plot((start:N), hsigma(start:N), 'g');
    plot((start:N), S2(start:N), 'm');
    grid on;
    xlabel('n');
    ylabel('\sigma');
    legend('S^2(x_N)','_{--}\sigma^2(x_n)', '^{--}\sigma^2(x_n)', 'S^2(x_n)');
end
