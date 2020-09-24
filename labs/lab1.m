function lab1()
    x = importdata('data1.txt');
    x = sort(x);
    
    Mmax = x(end);
    fprintf("Mmax = %.2f \n", Mmax);
    
    Mmin = x(1);
    fprintf("Mmin = %.2f \n", Mmin);
    
    R = Mmax - Mmin;
    fprintf("R = %.2f \n", R);
    
    mu = find_expected_value(x);
    fprintf("mu = %.4f \n", mu);
    
    S2 = find_dispersion(x);
    fprintf("S2 = %.4f \n", S2);
    
    m = find_n_intervals(x);
    fprintf("m = %d \n", m);
    
    intervals_table = group(x, m);
    for i=1:m-1
        fprintf("[%5.2f;%5.2f) ", intervals_table(1,i), intervals_table(1,i+1));
    end
    fprintf("[%5.2f;%5.2f]\n",intervals_table(1,m), Mmax);
    for i=1:m
        fprintf("%13d ", intervals_table(2,i));
    end
    fprintf("\n\n");

    grmin = -4;
    grmax = 4;
    figure(1);
    grid;
    hold on;
    hist_and_f(x,mu,S2,m,intervals_table,grmin,grmax)

    figure(2);
    grid;
    hold on;
    empir_and_F(x,mu,S2,grmin,grmax)
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
    S2 = sum((x - mu).^2) / (n - 1);
end
% кол-во интервалов
function m = find_n_intervals(x)
    n = length(x);
    m = floor(log2(n)) + 2;
end
% группировка значений выборки в m интервалов
function intervals_table = group(x, m)
    n = length(x);
    intervals_table = zeros(2,m);
    delta = (x(end) - x(1))/m; 
    for i = 0:m-1
        intervals_table(1, i+1) = x(1) + delta*i;
    end
    
    count = 0;
    for i = 1:n
        for j = 1:m-1
            if intervals_table(1,j) <= x(i) && x(i) < intervals_table(1, j+1)
                intervals_table(2, j) = intervals_table(2, j) + 1;
                count = count + 1;
                break;
            end    
        end    
    end
   intervals_table(2, m) = n - count;
end
% построение гистограммы и графика плотности
function hist_and_f(x,mu,S2,m,intervals_table,grmin,grmax)
    n = length(x);
    delta = (x(n)- x(1))/(n-1);

    graph = zeros(2,n+2);
    graph(1,1) = grmin;
    graph(2,1) = 0;
    for i = 1:n
        X = x(1) + delta*(i-1);
        graph(1,i+1) = X;
        graph(2,i+1) = f_normal(X,mu,S2);
    end
    graph(1,n+2) = grmax;
    graph(2,n+2) = 0;
   
    xx = zeros(m+4);
    yy = zeros(m+4);
    ddelta = (x(n)- x(1))/m;
    xx(1) = grmin;
    yy(1) = 0;
   
    for i = 1:m
        xx(i+1) = intervals_table(1,i);
        yy(i+1) = intervals_table(2,i) / (n * ddelta);
    end
    xx(m+2) = xx(m+1)+(xx(m+1)-xx(m));
    yy(m+2) = yy(m+1);
    xx(m+3) = xx(m+2);
    yy(m+3) = 0;
    xx(m+4)= grmax;
    yy(m+4)= 0;
    
    stairs(xx, yy, 'b'), grid;

    plot(graph(1,:), graph(2,:), 'r'),grid;  
end
% плотность нормального распределения
function y = f_normal(x,mx,dx) 
    y = exp(-((x-mx)^2)/2/dx)/sqrt(2*pi*dx);
end
% построение графиков эмпирической функции и функции распределения
function empir_and_F(x,mu,S2,grmin,grmax)
    n = length(x);
    delta = (grmax-grmin)/(n-1);

    graph = zeros(2,n);
    for i = 1:n
        X = grmin + delta*(i-1);
        graph(1,i) = X;        
        graph(2,i) = F_normal(X,mu,S2);
    end
    
    F_empir = zeros(n+2);    
    for i = 1:n
       F_empir(i+1) = empiric_F(x(i), x, n);
    end

    x =[grmin x grmax];
    F_empir(n+2) = 1;
    stairs(x, F_empir, 'b'),grid;
 
    plot(graph(1,:), graph(2,:), 'r'),grid;
end
% функция распределения нормальной случайной величины
function y = F_normal(x,mx,dx)
    y = 1/2 * (1 + erf((x - mx) /   sqrt(2*dx)));
end
% эмпирическая функция распределения
function Fi = empiric_F(X, x, n)
    count = 0;
    for i = 1:n
        if (x(i) <= X)
            count = count + 1;
        else
            break;
        end
    end
    Fi = count/n;
end