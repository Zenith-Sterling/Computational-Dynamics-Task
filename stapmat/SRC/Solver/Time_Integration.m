function [a,a_dot,a_dot2] = Time_Integration(N0,dt0,M,C,K,Q,a0,a0_dot)
%TIME_INTEGRATION 基于广义alpha法的时间积分方法

%% 参数设置
% 广义alpha值
alpha_f = 0.0;
alpha_m = 0.0;
beta = 0.25;
gama = 0.5;
% 最大迭代步数
N = N0;
% 时间步长
dt = dt0; %计算到N0*dt0时刻，时间维度上N+1
n = length(a0);
a = zeros(n,N+1);
a_dot = zeros(n,N+1);
a_dot2 = zeros(n,N+1);
% 赋初始值
a(:,1) = a0;
a_dot(:,1) = a0_dot;
a_dot2(:,1) = M\(Q-C*a0_dot-K*a0);


%% 中间参数
c_k = 1-alpha_f;
c_0 = (1-alpha_m)/(beta*dt*dt);
c_1 = c_k*gama/(beta*dt);
c_2 = dt*c_0;
c_3 = 0.5*c_2*dt-1;
c_4 = c_k*gama/beta-1;
c_5 = c_k*(gama/2/beta-1)*dt;


%% 循环迭代
for k = 1:N
    K_1 = c_k*K+c_0*M+c_1*C;
    Q_1 = Q-alpha_f*K*a(:,k)+M*(c_0*a(:,k)+c_2*a_dot(:,k)+c_3*a_dot2(:,k))+C*(c_1*a(:,k)+c_4*a_dot(:,k)+c_5*a_dot2(:,k));
    a(:,k+1) = K_1\Q_1;
    a_dot2(:,k+1) = (a(:,k+1)-a(:,k))/(beta*dt*dt)-a_dot(:,k)/(beta*dt)-(0.5/beta-1)*a_dot2(:,k);
    a_dot(:,k+1) = gama*(a(:,k+1)-a(:,k))/(beta*dt)+(1-gama/beta)*a_dot(:,k)+(1-0.5*gama/beta)*dt*a_dot2(:,k);
end


end

