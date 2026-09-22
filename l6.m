clear; clc; close all;

alpha = 6;                  % параметр корреляционной функции
h = 0.01;                   % шаг дискретизации, с
Tf = 1/alpha;               % постоянная времени фильтра
kf = sqrt(2/(alpha*h));     % коэффициент усиления фильтра 
N_trans = ceil(3*Tf/h);     % число шагов переходного процесса (50)

% Моменты времени для проверки стационарности (за пределами переходного процесса)
t1 = 1.0;                   
t2 = 2.0;                   
% Индексы в массиве после исключения переходного процесса
idx1 = round((t1 - N_trans*h)/h) + 1;   
idx2 = round((t2 - N_trans*h)/h) + 1;   

% Количество реализаций для формирования выборок
M = 500;

%% Проверка стационарности
sample1 = zeros(1, M);      % значения в момент t1
sample2 = zeros(1, M);      % значения в момент t2

for k = 1:M
    % Генерация одной реализации
    N_total = N_trans + max(idx1, idx2);
    u = rand(1, N_total);
    u2 = rand(1, N_total);
    Z = sqrt(-2*log(u)) .* cos(2*pi*u2);
    X = zeros(1, N_total);
    for i = 1:N_total-1
        X(i+1) = X(i) + h*(-alpha*X(i) + kf*alpha*Z(i));
    end
    X_proc = X(N_trans+1:end);
    X_proc = (X_proc - mean(X_proc)) / std(X_proc);   % стандартизация
    U = normcdf(X_proc);
    Z_out = 1 + sqrt(U);
    
    sample1(k) = Z_out(idx1);
    sample2(k) = Z_out(idx2);
end

% Критерий Колмогорова–Смирнова для двух выборок
[delta_stat, lambda_stat] = ks2sample(sample1, sample2);
fprintf('Проверка стационарности:\n');
fprintf('  Δ = %.4f\n', delta_stat);
fprintf('  λ = %.4f\n', lambda_stat);
if lambda_stat <= 1.22
    fprintf('  λ ≤ 1.22 -> выборки однородны, процесс стационарен\n');
else
    fprintf('  λ > 1.22 -> выборки неоднородны, процесс нестационарен\n');
end

% Графики эмпирических функций распределения
figure('Name','Проверка стационарности');
plot_ecdf(sample1, sample2, 't_1 = 1 c', 't_2 = 2 c');
title('Проверка стационарности: функции распределения');

%% Проверка эргодичности
% Первая выборка: значения в момент t1 по множеству реализаций
% Вторая выборка: значения одной длинной реализации
N_long = N_trans + M;       % длина для M отсчётов после переходного процесса
u = rand(1, N_long);
u2 = rand(1, N_long);
Z = sqrt(-2*log(u)) .* cos(2*pi*u2);
X_long = zeros(1, N_long);
for i = 1:N_long-1
    X_long(i+1) = X_long(i) + h*(-alpha*X_long(i) + kf*alpha*Z(i));
end
X_proc_long = X_long(N_trans+1:end);
X_proc_long = (X_proc_long - mean(X_proc_long)) / std(X_proc_long);
U_long = normcdf(X_proc_long);
Z_out_long = 1 + sqrt(U_long);

sample_time = Z_out_long(1:M);   % M значений одной реализации

% Критерий Колмогорова–Смирнова
[delta_erg, lambda_erg] = ks2sample(sample1, sample_time);
fprintf('\nПроверка эргодичности:\n');
fprintf('  Δ = %.4f\n', delta_erg);
fprintf('  λ = %.4f\n', lambda_erg);
if lambda_erg <= 1.22
    fprintf('  λ ≤ 1.22 -> выборки однородны, процесс эргодичен\n');
else
    fprintf('  λ > 1.22 -> выборки неоднородны, процесс неэргодичен\n');
end

figure('Name','Проверка эргодичности');
plot_ecdf(sample1, sample_time, 'По множеству реализаций', 'По одной реализации');
title('Проверка эргодичности: функции распределения');

%% Вспомогательные функции
function [delta, lambda] = ks2sample(x, y)
    % Двухвыборочный критерий Колмогорова–Смирнова
    n = length(x);
    m = length(y);
    z = sort([x(:); y(:)]);
    Fx = zeros(size(z));
    Fy = zeros(size(z));
    for i = 1:length(z)
        Fx(i) = sum(x <= z(i)) / n;
        Fy(i) = sum(y <= z(i)) / m;
    end
    delta = max(abs(Fx - Fy));
    lambda = delta * sqrt(n*m/(n+m));
end

function plot_ecdf(x, y, label1, label2)
    % График эмпирических функций распределения
    z = sort([x(:); y(:)]);
    Fx = zeros(size(z));
    Fy = zeros(size(z));
    for i = 1:length(z)
        Fx(i) = sum(x <= z(i)) / length(x);
        Fy(i) = sum(y <= z(i)) / length(y);
    end
    plot(z, Fx, 'b-', 'LineWidth', 1.5); hold on;
    plot(z, Fy, 'r--', 'LineWidth', 1.5);
    xlabel('z'); ylabel('F(z)');
    legend(label1, label2, 'Location', 'best');
    grid on;
end