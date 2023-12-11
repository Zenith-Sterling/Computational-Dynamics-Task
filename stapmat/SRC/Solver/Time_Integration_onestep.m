function [a_new,a_new_dot,a_new_dot2] = Time_Integration_onestep(M,C,K,Q,a,a_dot,a_dot2)
%TIME_INTEGRATION 基于广义alpha法的时间积分方法(一步)

%% 参数设置
% 广义alpha值
alpha_f = 2.0;
alpha_m = 0.5;
beta = 0.25;
gama = 0.5;
% 时间步长
dt = 0.1;


%% 中间参数
c_k = 1-alpha_f;
c_0 = (1-alpha_m)/(beta*dt*dt);
c_1 = c_k*gama/(beta*dt);
c_2 = dt*c_0;
c_3 = 0.5*c_2*dt-1;
c_4 = c_k*gama/beta-1;
c_5 = c_k*(gama/2/beta-1)*dt;


%% 计算t+dt时刻
K_1 = c_k*K+c_0*M+c_1*C;
Q_1 = Q-alpha_f*K*a+M*(c_0*a+c_2*a_dot+c_3*a_dot2)+C*(c_1*a+c_4*a_dot+c_5*a_dot2);
a_new = K_1\Q_1;
a_new_dot2 = (a_new-a)/(beta*dt*dt)-a_dot/(beta*dt)-(0.5/beta-1)*a_dot2;
a_new_dot = gama*(a_new-a)/(beta*dt)+(1-gama/beta)*a_dot+(1-0.5*gama/beta)*dt*a_dot2;


end

