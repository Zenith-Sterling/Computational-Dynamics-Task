%clc,clear
%load DATA.mat
figure;
trisurf(ELE, NODE(:, 1), NODE(:, 2), NODE(:, 3), DIS(:,1));
title('Displacement Field');
xlabel('X');
ylabel('Y');
zlabel('Z');
colorbar; % 显示颜色条
view(3); % 调整视角
grid on
axis equal
% 调整颜色映射
colormap jet; % 可以使用不同的颜色映射，如 'jet', 'hot', 'cool', 等