function [a,a_dot1,a_dot2] = modified_alpha4(N0,dt0,M,C,K,Q,a0,a0_dot)
%TIME_INTEGRATION 基于广义alpha法的时间积分方法

%% 参数设置
% 广义alpha值
alpha_f = 0.6;
alpha_1 = 1.5;
alpha_2 = 1.1;

gama_1 = alpha_1 - 0.5;
gama_2 = 0.5 - alpha_f + alpha_2;

beta_1 = (1+4*gama_1+4*gama_1^2)^2/16;
beta_2 = (1+4*gama_2+4*gama_2^2)^2/16;
% 最大迭代步数
N = N0;
% 时间步长
dt = dt0; %计算到N0*dt0时刻，时间维度上N+1
% 储存空间预设
n = length(a0);
a = zeros(n,N+1);
a_dot1 = zeros(n,N+1);
a_dot2 = zeros(n,N+1);
% 赋初始值
a(:,1) = a0;
a_dot1(:,1) = a0_dot;
a_dot2(:,1) = M\(Q-C*a0_dot-K*a0);
A_dot1 = zeros(n,1);
A_dot2 = zeros(n,1);
A_dot3 = zeros(n,1);
% A_dot1 = M\(0-C*a_dot2(:,1)-K*a0_dot);
% A_dot2 = M\(0-C*A_dot1-K*a_dot2(:,1));
% A_dot3 = M\(0-C*A_dot2-K*A_dot1);



%% 中间参数
alpha_f0 = 0.0;
alpha_m0 = 0.0;
beta0 = 0.25;
gama0 = 0.5;
c_k = 1-alpha_f0;
c_0 = (1-alpha_m0)/(beta0*dt*dt);
c_1 = c_k*gama0/(beta0*dt);
c_2 = dt*c_0;
c_3 = 0.5*c_2*dt-1;
c_4 = c_k*gama0/beta0-1;
c_5 = c_k*(gama0/2/beta0-1)*dt;


for k = 1:2
    K_1 = c_k*K+c_0*M+c_1*C;
    Q_1 = Q-alpha_f*K*a(:,k)+M*(c_0*a(:,k)+c_2*a_dot1(:,k)+c_3*a_dot2(:,k))+C*(c_1*a(:,k)+c_4*a_dot1(:,k)+c_5*a_dot2(:,k));
    a(:,k+1) = K_1\Q_1;
    a_dot2(:,k+1) = (a(:,k+1)-a(:,k))/(beta*dt*dt)-a_dot1(:,k)/(beta*dt)-(0.5/beta-1)*a_dot2(:,k);
    a_dot1(:,k+1) = gama*(a(:,k+1)-a(:,k))/(beta*dt)+(1-gama/beta)*a_dot1(:,k)+(1-0.5*gama/beta)*dt*a_dot2(:,k);
end

A_dot1 = M\(0-C*a_dot2(:,2)-K*a_dot1(:,2));
A_dot2 = M\(0-C*A_dot1-K*a_dot2(:,2));
A_dot3 = M\(0-C*A_dot2-K*A_dot1);


%% 循环迭代
for k = 2:N
    M_1 = M + gama_1*dt*C + beta_1*dt^2*K;
    Q_1 = Q - M*(1-alpha_1)*(a_dot2(:,k) + dt*A_dot1 + 0.5*dt^2*A_dot2 + dt^3*A_dot3/6) - ...
        C*(a_dot1(:,k) + (1-gama_1)*dt*a_dot2(:,k) + (0.5-gama_1)*dt^2*A_dot1 + (1/6-gama_1/2)*dt^3*A_dot2 + (1/24-gama_1/6)*dt^4*A_dot3) - ...
        K*(a(:,k) + dt*a_dot1(:,k) + (0.5-beta_1)*dt^2*a_dot2(:,k) + (1/6-beta_1)*dt^3*A_dot1 + (1/24-beta_1/2)*dt^4*A_dot2 + (1/120-beta_1/6)*dt^5*A_dot3);
    a_dot2(:,k+1) = M_1\Q_1;
    P_n = a_dot2(:,k+1) - a_dot2(:,k) - dt*A_dot1 - 0.5*dt^2*A_dot2 - dt^3*A_dot3/6;
    a(:,k+1) = a(:,k) + dt*a_dot1(:,k) + 0.5*dt^2*a_dot2(:,k) + dt^3*A_dot1/6 + dt^4*A_dot2/24 + dt^5*A_dot3/120 + beta_1*dt^2*P_n;
    a_dot1(:,k+1) = a_dot1(:,k) + dt*a_dot2(:,k) + 0.5*dt^2*A_dot1 + dt^3*A_dot2/6 + dt^4*A_dot3/24 + gama_1*dt*P_n;
    M_2 = alpha_2*M + alpha_f*gama_2*dt*C + alpha_f*beta_2*dt^2*K;
    Q_2 = 0 - M*(1-alpha_2)*A_dot3 - C*(A_dot2+alpha_f*(1-gama_2)*dt*A_dot3) - ...
        K*(A_dot1 + alpha_f*dt*A_dot2 + alpha_f*dt^2*(0.5-beta_2)*A_dot3);
    A_dot3_new = M_2\Q_2;
    A_dot1_new = A_dot1 + dt*A_dot2 + 0.5*dt^2*A_dot3 + dt^2*beta_2*(A_dot3_new-A_dot3);
    A_dot2_new = A_dot2 + dt*A_dot3 + dt*gama_2*(A_dot3_new-A_dot3);
    A_dot1 = A_dot1_new;
    A_dot2 = A_dot2_new;
    A_dot3 = A_dot3_new;


end



end

