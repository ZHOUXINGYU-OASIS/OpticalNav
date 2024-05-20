function plot_OD_estimated_error(Err, tk, relative_state)
%% Plot the estimated errors in the main code
tk = tk / 24 / 3600;
C2 = [102, 173, 194; 
      36,  59,  66;
      232, 69,  69;
      194, 148, 102;
      54,  43,  33] ./ 255;
C6 = [235, 75,  55;
      77,  186, 216;
      58,  84,  141;
      2,   162, 136;
      245, 155, 122] ./ 255;
%% Position estimated errors
figure(1)
set(gcf, 'unit', 'centimeters', 'Position', [5, 5, 12, 8]);
h(1) = plot(tk, Err(1,:), '-', 'color', 'b', 'linewidth', 2);
hold on;
h(2) = plot(tk, Err(2,:), '-', 'color', 'r', 'linewidth', 2);
h(3) = plot(tk, Err(3,:), '-', 'color', 'g', 'linewidth', 2);
axis([-inf, inf, -inf, inf])
Legend = legend(h, '$x$', '$y$', '$z$');
set(Legend, 'FontSize', 12, 'Location', 'best', 'FontName', ...
    'Times New Roman', 'Interpreter', 'latex', 'NumColumn', 1);
grid on;
box on;
set(gca,'FontName', 'Times New Roman', 'FontSize', 12, 'gridlinestyle', '--');
xlabel('Epoch$\,\rm(day)$', 'FontName', 'Times New Roman', ...
    'FontSize', 12, 'Interpreter', 'latex');
ylabel('Position estimated errors$\,\rm(km)$', 'FontName', 'Times New Roman', ...
    'FontSize', 12, 'Interpreter', 'latex');
%% Velocity estimated errors
figure(2)
set(gcf, 'unit', 'centimeters', 'Position', [5, 5, 12, 8]);
h(1) = plot(tk, Err(4,:) * 1e3, '-', 'color', 'b', 'linewidth', 2);
hold on;
h(2) = plot(tk, Err(5,:) * 1e3, '-', 'color', 'r', 'linewidth', 2);
h(3) = plot(tk, Err(6,:) * 1e3, '-', 'color', 'g', 'linewidth', 2);
axis([-inf, inf, -inf, inf])
Legend = legend(h, '$\dot{x}$', '$\dot{y}$', '$\dot{z}$');
set(Legend, 'FontSize', 12, 'Location', 'best', 'FontName', ...
    'Times New Roman', 'Interpreter', 'latex', 'NumColumn', 1);
grid on;
box on;
set(gca,'FontName', 'Times New Roman', 'FontSize', 12, 'gridlinestyle', '--');
xlabel('Epoch$\,\rm(day)$', 'FontName', 'Times New Roman', ...
    'FontSize', 12, 'Interpreter', 'latex');
ylabel('Velocity estimated errors$\,\rm(m/s)$', 'FontName', 'Times New Roman', ...
    'FontSize', 12, 'Interpreter', 'latex');
%% Relative coordinate system
len = length(tk);
for i = 1: len
    a = eye(3);
    a(:, 1) = relative_state(1:3, i);
    b = Schmidt_orthogonalization(a);
    Err(1:3, i) = b' * Err(1:3, i);
    Err(4:6, i) = b' * Err(4:6, i);
end
%% Position estimated errors
figure(3)
set(gcf, 'unit', 'centimeters', 'Position', [5, 5, 12, 8]);
h(1) = plot(tk, Err(1,:), '-', 'color', 'b', 'linewidth', 2);
hold on;
h(2) = plot(tk, Err(2,:), '-', 'color', 'r', 'linewidth', 2);
h(3) = plot(tk, Err(3,:), '-', 'color', 'g', 'linewidth', 2);
axis([-inf, inf, -inf, inf])
Legend = legend(h, '$\tilde{x}$', '$\tilde{y}$', '$\tilde{z}$');
set(Legend, 'FontSize', 12, 'Location', 'best', 'FontName', ...
    'Times New Roman', 'Interpreter', 'latex', 'NumColumn', 1);
grid on;
box on;
set(gca,'FontName', 'Times New Roman', 'FontSize', 12, 'gridlinestyle', '--');
xlabel('Epoch$\,\rm(day)$', 'FontName', 'Times New Roman', ...
    'FontSize', 12, 'Interpreter', 'latex');
ylabel('Position estimated errors$\,\rm(km)$', 'FontName', 'Times New Roman', ...
    'FontSize', 12, 'Interpreter', 'latex');
%% Velocity estimated errors
figure(4)
set(gcf, 'unit', 'centimeters', 'Position', [5, 5, 12, 8]);
h(1) = plot(tk, Err(4,:) * 1e3, '-', 'color', 'b', 'linewidth', 2);
hold on;
h(2) = plot(tk, Err(5,:) * 1e3, '-', 'color', 'r', 'linewidth', 2);
h(3) = plot(tk, Err(6,:) * 1e3, '-', 'color', 'g', 'linewidth', 2);
axis([-inf, inf, -inf, inf])
Legend = legend(h, '$\dot{\tilde{x}}$', '$\dot{\tilde{y}}$', '$\dot{\tilde{z}}$');
set(Legend, 'FontSize', 12, 'Location', 'best', 'FontName', ...
    'Times New Roman', 'Interpreter', 'latex', 'NumColumn', 1);
grid on;
box on;
set(gca,'FontName', 'Times New Roman', 'FontSize', 12, 'gridlinestyle', '--');
xlabel('Epoch$\,\rm(day)$', 'FontName', 'Times New Roman', ...
    'FontSize', 12, 'Interpreter', 'latex');
ylabel('Velocity estimated errors$\,\rm(m/s)$', 'FontName', 'Times New Roman', ...
    'FontSize', 12, 'Interpreter', 'latex');