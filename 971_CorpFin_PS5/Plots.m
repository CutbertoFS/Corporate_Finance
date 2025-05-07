% Plots

%--- Figure 1: Investor value function [b] ---
figure;
set(gca, 'FontSize', 14);
set(gcf, 'Position', [100, 100, 500, 500]);
plot(a_grid, b, 'LineWidth', 2, 'Color', 'k');
hold on;
plot(a_grid, b_FB, 'LineWidth', 2, 'Color', 'b');
xlabel('Agents payoff a', 'Interpreter', 'latex', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Investor''s payoff b', 'Interpreter', 'latex', 'FontSize', 18, 'FontWeight', 'bold');
xlim([10, 62]);
ylim([30, 100]);
legend({'Asymmetric Info', 'First Best'}, 'Location', 'northeast', 'FontSize', 14);
grid on;

%% 

figure;
set(gcf, 'Position', [100, 100, 1000, 900]);

%--- Subplot 1:  Subplot 1: a_H
subplot(3,2,1);
plot(a_grid, a_H, 'k', 'LineWidth', 2);
hold on;
plot(a_grid, a_H_FB, 'b', 'LineWidth', 2);
xlabel('$a$', 'Interpreter', 'latex');
ylabel('$a_H(a)$', 'Interpreter', 'latex');
title('Next Value $a_H(a)$', 'Interpreter', 'latex');
xlim([10, 110]);
ylim([0, 70]);
legend({'Asymmetric Info', 'First Best'}, 'FontSize', 9, 'Location', 'best');
grid on;

%--- Subplot 2: a_L
subplot(3,2,2);
plot(a_grid(2:end), a_L(2:end), 'k', 'LineWidth', 2);
hold on;
plot(a_grid(2:end), a_L_FB(2:end), 'b', 'LineWidth', 2);
xlabel('$a$', 'Interpreter', 'latex');
ylabel('$a_L(a)$', 'Interpreter', 'latex');
title('Next Value $a_L(a)$', 'Interpreter', 'latex');
xlim([10, 110]);
ylim([0, 70]);
legend({'Asymmetric Info', 'First Best'}, 'FontSize', 9, 'Location', 'best');
grid on;

%--- Subplot 3: p_H
subplot(3,2,3);
plot(a_grid, p_H, 'k', 'LineWidth', 2);
hold on;
plot(a_grid, p_H_FB, 'b', 'LineWidth', 2);
xlabel('$a$', 'Interpreter', 'latex');
ylabel('$p_H(a)$', 'Interpreter', 'latex');
title('Probability $p_H$', 'Interpreter', 'latex');
xlim([10, 110]);
ylim([0, 1]);
legend({'Asymmetric Info', 'First Best'}, 'FontSize', 9, 'Location', 'best');
grid on;

%--- Subplot 4: p_L
subplot(3,2,4);
plot(a_grid, p_L, 'k', 'LineWidth', 2);
hold on;
plot(a_grid, p_L_FB, 'b', 'LineWidth', 2);
xlabel('$a$', 'Interpreter', 'latex');
ylabel('$p_L(a)$', 'Interpreter', 'latex');
title('Probability $p_L$', 'Interpreter', 'latex');
xlim([10, 110]);
ylim([0, 1]);
legend({'Asymmetric Info', 'First Best'}, 'FontSize', 9, 'Location', 'best');
grid on;

%--- Subplot 5: d_H
subplot(3,2,5);
plot(a_grid, d_H, 'k', 'LineWidth', 2);
hold on;
plot(a_grid, d_H_FB, 'b', 'LineWidth', 2);
xlabel('$a$', 'Interpreter', 'latex');
ylabel('$d_H(a)$', 'Interpreter', 'latex');
title('Dividend $d_H$', 'Interpreter', 'latex');
xlim([10, 110]);
legend({'Asymmetric Info', 'First Best'}, 'FontSize', 9, 'Location', 'best');
grid on;

%--- Subplot 6: d_L
subplot(3,2,6);
plot(a_grid, d_L, 'k', 'LineWidth', 2);
hold on;
plot(a_grid, d_L_FB, 'b', 'LineWidth', 2);
xlabel('$a$', 'Interpreter', 'latex');
ylabel('$d_L(a)$', 'Interpreter', 'latex');
title('Dividend $d_L$', 'Interpreter', 'latex');
xlim([10, 110]);
legend({'Asymmetric Info', 'First Best'}, 'FontSize', 9, 'Location', 'best');
grid on;

set(gcf, 'PaperPositionMode', 'auto'); % Automatically adjust the position
set(gcf, 'PaperSize', [10 8]);         % Set paper size for export (in inches)

exportgraphics(gcf, 'Figure_02.pdf', 'ContentType', 'vector');V