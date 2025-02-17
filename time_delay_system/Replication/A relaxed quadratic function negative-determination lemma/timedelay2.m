clc
clear
Tt = 0.3;
Tg = 0.1;
R = 0.05;
D = 1.0; 
bias = 21; 
M = 10;
Kp = 0.1;
Ki =0.05;


A = [-D/M,1/M,0,0;0,-1/Tt,1/Tt,0;-1/(R*Tg),0,-1/Tg,0;bias, 0,0 ,0;];
Ad = [0 0 0 0;0 0 0 0;-Kp*bias/Tg 0 0 -Ki/Tg;0 0 0 0];


h1 = 2;

n =4;

beta_values = 0:0.2:1; % Beta 的搜索范围
tolerance = 1e-4; % 二分法停止条件
%for beta = beta_values
beta = 0;
    h2_min = h1; % h2 必须大于 h1
        h2_max = 20; % 假设一个较大的上界
        feasible_h2 = -1; % 存储当前 beta 下的最大可行 h2
        
        while h2_max - h2_min > tolerance
            h2_test = (h2_max + h2_min) / 2; % 取中间值



delta1=0;
delta2 = beta;%beta自选
delta3= 0;
delta4=1-beta;

P = sdpvar(5*n,5*n,'symmetric');
Q1 = sdpvar(n,n,'symmetric');
Q2 = sdpvar(n,n,'symmetric');
R1 = sdpvar(n,n,'symmetric');
R2 = sdpvar(n,n,'symmetric');
L1 = sdpvar(2*n,n,'full');
L2 = sdpvar(2*n,n,'full');
N1 = sdpvar(12*n,3*n,'full');
N2 = sdpvar(12*n,3*n,'full');
h2 = sdpvar(1, 1); 

 h12 = h2_test - h1;

% 定义 e_i 向量的生成函数
e = @(i, n) [zeros(n, n * (i - 1)), eye(n), zeros(n, n * (12 - i))];

% 定义 E_i 矩阵的生成函数
E = @(i, n) vertcat( ...
    e(i, n) - e(i + 1, n), ...
    e(i, n) + e(i + 1, n) - 2 * e(i + 4, n), ...
    e(i, n) - e(i + 1, n) + 6 * e(i + 4, n) - 12 * e(i + 7, n) ...
);
E1 = E(1, n);
E2 = E(2, n);
E3 = E(3, n);

es= A*e(1,n)+Ad*e(3,n);
PI0 = [zeros(n, 12*n);zeros(n, 12*n);zeros(n, 12*n);zeros(n, 12*n);e(9,n)+e(10,n)];
PI1h1 = [e(1,n);h1*e(5,n);e(11,n)+e(12,n);(h1^2)*e(8,n);h12^2*e(10,n) + h12*e(11,n)];
PI1h2 = [e(1,n);h1*e(5,n);e(11,n)+e(12,n);(h1^2)*e(8,n);h12^2*e(9,n)];
Eb = h12*e(2,n)-e(11,n)-e(12,n);
PI2 = [es;e(1,n)-e(2,n);e(2,n)-e(4,n);h1*(e(1,n)-e(5,n));Eb];

y0 = PI0'*P*PI2+PI2'*P*PI0;
y1h1 = PI1h1'*P*PI2+(PI1h1'*P*PI2)';
y1h2 = PI1h2'*P*PI2+(PI1h2'*P*PI2)';
y2= e(1,n)'*Q1*e(1,n) -e(2,n)'*(Q1-Q2)*e(2,n)-e(4,n)'*Q2*e(4,n);

E = @(i, n) vertcat( ...
    e(i, n) - e(i + 1, n), ...
    e(i, n) + e(i + 1, n) - 2 * e(i + 4, n), ...
    e(i, n) - e(i + 1, n) + 6 * e(i + 4, n) - 12 * e(i + 7, n) ...
);
E1 = E(1, n);
E2 = E(2, n);
E3 = E(3, n);

R1_hat = [R1 zeros(n,n) zeros(n,n);zeros(n,n) 3*R1 zeros(n,n);zeros(n,n) zeros(n,n) 5*R1];
R2_hat = [R2 zeros(n,n) zeros(n,n);zeros(n,n) 3*R2 zeros(n,n);zeros(n,n) zeros(n,n) 5*R2];

y3 =es'*((h1^2)*R1+(h12^2)*R2)*es-E1'*R1_hat*E1;
y4h1=[E2;E3]'*[2*R2_hat zeros(3*n,3*n);zeros(3*n,3*n) R2_hat]*[E2;E3] ...
+ [E2;E3]'*[N1';zeros(12*n,3*n)'] + ([E2;E3]'*[N1';zeros(12*n,3*n)'])';

y4h2=[E2;E3]'*[R2_hat zeros(3*n,3*n);zeros(3*n,3*n) 2*R2_hat]*[E2;E3]...
+ [E2;E3]'*[zeros(12*n,3*n)';N2'] + ([E2;E3]'*[zeros(12*n,3*n)';N2'])';

y5h1 = [e(6,n)',e(11,n)']*L1*[-e(11,n)] +([e(6,n)',e(11,n)']*L1*[-e(11,n)])' ...
        +[e(7,n)',e(12,n)']*L2*[h12*e(7,n)-e(12,n)]+([e(7,n)',e(12,n)']*L2*[h12*e(7,n)-e(12,n)])';
y5h2 = [e(6,n)',e(11,n)']*L1*[h12*e(6,n)-e(11,n)] +([e(6,n)',e(11,n)']*L1*[h12*e(6,n)-e(11,n)])' ...
        +[e(7,n)',e(12,n)']*L2*[-e(12,n)]+([e(7,n)',e(12,n)']*L2*[-e(12,n)])';

yh1 = y1h1+y2+y3-y4h1+y5h1;
yh2 = y1h2+y2+y3-y4h2+y5h2;

LMI1 = [yh1 N2;N2' -R2_hat];
LMI2 = [yh1-(beta^2)*(h12^2)*y0 N2;N2' -R2_hat];
LMI3 = [yh2 N1;N1' -R2_hat];
LMI4 = [yh2-((1-beta)^2)*(h12^2)*y0,N1;N1' -R2_hat];



Fcond = [P>=0, Q1>=0, Q2>=0, R1>=0, R2>=0, LMI1<=0, LMI2<=0, LMI3<=0, LMI4<=0];
objective = -h2; % 最大化 h2
ops = sdpsettings('solver', 'sdpt3', 'verbose', 1);
    diagnostics = optimize(Fcond, [], ops);
[m,p] = checkset(Fcond);%%%验证可解性
tmin = min(m);
        if tmin > 0
            feasible_h2 = h2_test; % 更新当前最大可行 h2
            h2_min = h2_test; % 增大 h2 下界
        else
            h2_max = h2_test; % 减小 h2 上界
        end
        end
    % 记录结果
  if feasible_h2 > 0
        disp(['Beta = ', num2str(beta), ', Max h2 = ', num2str(feasible_h2)]);
       
    else
        disp(['Beta = ', num2str(beta), ' is infeasible']);
    end

%end
%%
clc
clear
Tt = 0.3;
Tg = 0.1;
R = 0.05;
D = 1.0; 
bias = 21; 
M = 10;
Kp = 0.1;
Ki =0.05;


A = [-D/M,1/M,0,0;0,-1/Tt,1/Tt,0;-1/(R*Tg),0,-1/Tg,0;bias, 0,0 ,0;];
Ad = [0 0 0 0;0 0 0 0;-Kp*bias/Tg 0 0 -Ki/Tg;0 0 0 0];


h1 = 2;

n =4;

tolerance = 1e-3; % 二分法停止条件

    h2_min = h1; % h2 必须大于 h1
        h2_max = 15; % 假设一个较大的上界
        feasible_h2 = -1; % 存储当前 beta 下的最大可行 h2
        
        while h2_max - h2_min > tolerance
            h2_test = (h2_max + h2_min) / 2; % 取中间值


P = sdpvar(5*n,5*n,'symmetric');
Q1 = sdpvar(n,n,'symmetric');
Q2 = sdpvar(n,n,'symmetric');
R1 = sdpvar(n,n,'symmetric');
R2 = sdpvar(n,n,'symmetric');
L1 = sdpvar(2*n,n,'full');
L2 = sdpvar(2*n,n,'full');
N1 = sdpvar(12*n,3*n,'full');
N2 = sdpvar(12*n,3*n,'full');
h2 = sdpvar(1, 1); 

 h12 = h2_test - h1;

% 定义 e_i 向量的生成函数
e = @(i, n) [zeros(n, n * (i - 1)), eye(n), zeros(n, n * (12 - i))];

% 定义 E_i 矩阵的生成函数
E = @(i, n) vertcat( ...
    e(i, n) - e(i + 1, n), ...
    e(i, n) + e(i + 1, n) - 2 * e(i + 4, n), ...
    e(i, n) - e(i + 1, n) + 6 * e(i + 4, n) - 12 * e(i + 7, n) ...
);
E1 = E(1, n);
E2 = E(2, n);
E3 = E(3, n);

es= A*e(1,n)+Ad*e(3,n);
PI0 = [zeros(n, 12*n);zeros(n, 12*n);zeros(n, 12*n);zeros(n, 12*n);e(9,n)+e(10,n)];
PI1h1 = [e(1,n);h1*e(5,n);e(11,n)+e(12,n);(h1^2)*e(8,n);h12^2*e(10,n) + h12*e(11,n)];
PI1h2 = [e(1,n);h1*e(5,n);e(11,n)+e(12,n);(h1^2)*e(8,n);h12^2*e(9,n)];
Eb = h12*e(2,n)-e(11,n)-e(12,n);
PI2 = [es;e(1,n)-e(2,n);e(2,n)-e(4,n);h1*(e(1,n)-e(5,n));Eb];

y0 = PI0'*P*PI2+PI2'*P*PI0;
y1h1 = PI1h1'*P*PI2+(PI1h1'*P*PI2)';
y1h2 = PI1h2'*P*PI2+(PI1h2'*P*PI2)';
y2= e(1,n)'*Q1*e(1,n) -e(2,n)'*(Q1-Q2)*e(2,n)-e(4,n)'*Q2*e(4,n);

E = @(i, n) vertcat( ...
    e(i, n) - e(i + 1, n), ...
    e(i, n) + e(i + 1, n) - 2 * e(i + 4, n), ...
    e(i, n) - e(i + 1, n) + 6 * e(i + 4, n) - 12 * e(i + 7, n) ...
);
E1 = E(1, n);
E2 = E(2, n);
E3 = E(3, n);

R1_hat = [R1 zeros(n,n) zeros(n,n);zeros(n,n) 3*R1 zeros(n,n);zeros(n,n) zeros(n,n) 5*R1];
R2_hat = [R2 zeros(n,n) zeros(n,n);zeros(n,n) 3*R2 zeros(n,n);zeros(n,n) zeros(n,n) 5*R2];

y3 =es'*((h1^2)*R1+(h12^2)*R2)*es-E1'*R1_hat*E1;
y4h1=[E2;E3]'*[2*R2_hat zeros(3*n,3*n);zeros(3*n,3*n) R2_hat]*[E2;E3] ...
+ [E2;E3]'*[N1';zeros(12*n,3*n)'] + ([E2;E3]'*[N1';zeros(12*n,3*n)'])';

y4h2=[E2;E3]'*[R2_hat zeros(3*n,3*n);zeros(3*n,3*n) 2*R2_hat]*[E2;E3]...
+ [E2;E3]'*[zeros(12*n,3*n)';N2'] + ([E2;E3]'*[zeros(12*n,3*n)';N2'])';

y5h1 = [e(6,n)',e(11,n)']*L1*[-e(11,n)] +([e(6,n)',e(11,n)']*L1*[-e(11,n)])' ...
        +[e(7,n)',e(12,n)']*L2*[h12*e(7,n)-e(12,n)]+([e(7,n)',e(12,n)']*L2*[h12*e(7,n)-e(12,n)])';
y5h2 = [e(6,n)',e(11,n)']*L1*[h12*e(6,n)-e(11,n)] +([e(6,n)',e(11,n)']*L1*[h12*e(6,n)-e(11,n)])' ...
        +[e(7,n)',e(12,n)']*L2*[-e(12,n)]+([e(7,n)',e(12,n)']*L2*[-e(12,n)])';

yh1 = y1h1+y2+y3-y4h1+y5h1;
yh2 = y1h2+y2+y3-y4h2+y5h2;

LMI1 = [yh1 N2;N2' -R2_hat];
LMI2 = [yh1-(h12^2)*y0 N2;N2' -R2_hat];
LMI3 = [yh2 N1;N1' -R2_hat];




Fcond = [P>=0, Q1>=0, Q2>=0, R1>=0, R2>=0, LMI1<=0, LMI2<=0, LMI3<=0];
objective = -h2; % 最大化 h2
ops = sdpsettings('solver', 'sdpt3', 'verbose', 1);
    diagnostics = optimize(Fcond, [], ops);
[m,p] = checkset(Fcond);%%%验证可解性
tmin = min(m);
        if tmin > 0
            feasible_h2 = h2_test; % 更新当前最大可行 h2
            h2_min = h2_test; % 增大 h2 下界
        else
            h2_max = h2_test; % 减小 h2 上界
        end
        end
    % 记录结果
    if feasible_h2 > 0
        disp([' Max h2 = ', num2str(feasible_h2)]);
      
    end
    

%%
clc;
clear;

% 定义参数
Tt = 0.3;
Tg = 0.1;
R = 0.05;
D = 1.0;
bias = 21;
M = 10;
Kp = 0.1;
Ki = 0.1;

A = [-D/M, 1/M, 0, 0;
      0, -1/Tt, 1/Tt, 0;
     -1/(R*Tg), 0, -1/Tg, 0;
      bias, 0, 0, 0];
Ad = [0, 0, 0, 0;
      0, 0, 0, 0;
      -Kp*bias/Tg, 0, 0, -Ki/Tg;
      0, 0, 0, 0];

% 定义延迟时间的区间
d_min = 2;
d_max = 11.091;

% 时间区间
t_span = [0, 60];  % 模拟时间范围

% 随机生成案例
num_cases = 200;
delays = d_min + (d_max - d_min) * rand(num_cases, 1); % 随机生成 d(t) ∈ [2, 11.091]
freq_deviations = -0.02 + (0.02 - (-0.02)) * rand(num_cases, 1); % 随机初始偏差

% 时间采样点
time_points = linspace(t_span(1), t_span(2), 100);

% 初始化存储状态响应的数组
x1_response = zeros(num_cases, length(time_points));
x2_response = zeros(num_cases, length(time_points));
x3_response = zeros(num_cases, length(time_points));
x4_response = zeros(num_cases, length(time_points));

% 遍历每个案例
for i = 1:num_cases
    % 当前延迟时间
    delay = delays(i);
    
    % 初始条件函数 phi(t)
    phi = @(t) [freq_deviations(i); 0; 0; 0];
    
    % 定义延迟微分方程
    dde = @(t, x, x_delayed) A * x + Ad * x_delayed;
    
    % 求解延迟微分方程
    sol = dde23(dde, delay, phi, t_span);
    
    % 插值到统一的时间点，确保维度匹配
    x_interpolated = deval(sol, time_points); % 返回大小为 [4 x length(time_points)]
    
    % 存储状态响应
    x1_response(i, :) = x_interpolated(1, :);
    x2_response(i, :) = x_interpolated(2, :);
    x3_response(i, :) = x_interpolated(3, :);
    x4_response(i, :) = x_interpolated(4, :);
end

% 定义不同颜色的色图
colors = jet(num_cases);

% 绘制图像
[X, Y] = meshgrid(time_points, 1:num_cases);

figure;

% x1 绘图
subplot(2, 2, 1);
hold on;
for i = 1:num_cases
    fill3([time_points, fliplr(time_points)], [i * ones(1, length(time_points)), fliplr(i * ones(1, length(time_points)))], ...
          [x1_response(i, :), zeros(1, length(time_points))], colors(i, :), 'EdgeColor', 'none');
end
view([-30 30]); % 调整视角
title('x_1(t)');
xlabel('Time [t]');
ylabel('Cases');
zlabel('x_1(t)');
grid on;
hold off;

% x2 绘图
subplot(2, 2, 2);
hold on;
for i = 1:num_cases
    fill3([time_points, fliplr(time_points)], [i * ones(1, length(time_points)), fliplr(i * ones(1, length(time_points)))], ...
          [x2_response(i, :), zeros(1, length(time_points))], colors(i, :), 'EdgeColor', 'none');
end
view([-30 30]);
title('x_2(t)');
xlabel('Time [t]');
ylabel('Cases');
zlabel('x_2(t)');
grid on;
hold off;

% x3 绘图
subplot(2, 2, 3);
hold on;
for i = 1:num_cases
    fill3([time_points, fliplr(time_points)], [i * ones(1, length(time_points)), fliplr(i * ones(1, length(time_points)))], ...
          [x3_response(i, :), zeros(1, length(time_points))], colors(i, :), 'EdgeColor', 'none');
end
view([-30 30]);
title('x_3(t)');
xlabel('Time [t]');
ylabel('Cases');
zlabel('x_3(t)');
grid on;
hold off;

% x4 绘图
subplot(2, 2, 4);
hold on;
for i = 1:num_cases
    fill3([time_points, fliplr(time_points)], [i * ones(1, length(time_points)), fliplr(i * ones(1, length(time_points)))], ...
          [x4_response(i, :), zeros(1, length(time_points))], colors(i, :), 'EdgeColor', 'none');
end
view([-30 30]);
title('x_4(t)');
xlabel('Time [t]');
ylabel('Cases');
zlabel('x_4(t)');
grid on;
hold off;

sgtitle('Responses of system states.');
