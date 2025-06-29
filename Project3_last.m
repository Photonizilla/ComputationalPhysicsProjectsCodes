rng(0);
N_list = [16,32,64,128,256,512,1024];
M = 100000;

% p = 0.5
mean_x_p5 = zeros(size(N_list));
var_x_p5 = zeros(size(N_list));
figure(1); clf; colors = lines(length(N_list));

for i = 1:length(N_list)
    N = N_list(i);
    steps = 2*(rand(M, N) < 0.5) - 1;
    x = sum(steps, 2);
    
    mean_x_p5(i) = mean(x);
    var_x_p5(i) = mean(x.^2) - mean_x_p5(i)^2;
    
    [counts, edges] = histcounts(x, 'BinMethod', 'integers');
    x_vals = (edges(1:end-1) + edges(2:end)) / 2;
    P = counts / M;
    
    subplot(2,1,1); hold on;
    plot(x_vals, P, 'o-', 'Color', colors(i,:), 'MarkerSize', 3);
    
    subplot(2,1,2); hold on;
    sigma = sqrt(var_x_p5(i));
    gauss = exp(-(x_vals).^2/(2*sigma^2)) / (sqrt(2*pi)*sigma);
    plot(x_vals, gauss, '--', 'Color', colors(i,:));
end
subplot(2,1,1); title('P_N(x) (p=0.5)'); xlabel('x'); ylabel('P_N(x)'); 
legend(cellstr(num2str(N_list')), 'Location', 'best');
subplot(2,1,2); title('Gaussian Fit (p=0.5)'); xlabel('x'); ylabel('P(x)');

% p = 0.7
mean_x_p7 = zeros(size(N_list));
var_x_p7 = zeros(size(N_list));
for i = 1:length(N_list)
    N = N_list(i);
    steps = 2*(rand(M, N) < 0.7) - 1;
    x = sum(steps, 2);
    mean_x_p7(i) = mean(x);
    var_x_p7(i) = mean(x.^2) - mean_x_p7(i)^2;
end

% Log-log plot for scaling law
figure(2); clf;
loglog(N_list, var_x_p5, 'bo-', 'LineWidth', 1.5); hold on;
loglog(N_list, var_x_p7, 'rs-', 'LineWidth', 1.5);
xlabel('N'); ylabel('\langle \Delta x_N^2 \rangle');
legend('p=0.5', 'p=0.7', 'Location', 'best');
grid on;

% Fit exponent for p=0.5
p5_fit = polyfit(log(N_list), log(var_x_p5), 1);
exponent_p5 = p5_fit(1);
fprintf('Exponent for p=0.5: %.4f\n', exponent_p5);