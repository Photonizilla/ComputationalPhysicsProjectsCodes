function polymer_chain_simulation
    rng(0);
    N_list = [1:10, 20:10:100, 200:100:1000, 2000:1000:10000];
    M = min(10000, 1000000 ./ N_list);
    M(M < 100) = 100;
    
    [avg_R2_fj, avg_Rg2_fj, theory_fj] = freely_jointed(N_list, M);
    [avg_R2_fr, avg_Rg2_fr] = freely_rotating(N_list, M, 68*pi/180);
    
    figure;
    loglog(N_list, avg_R2_fj, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
    loglog(N_list, avg_Rg2_fj, 'rs-', 'LineWidth', 1.5, 'MarkerSize', 6);
    loglog(N_list, theory_fj, 'g--', 'LineWidth', 1.5);
    loglog(N_list, avg_R2_fr, 'c^-', 'LineWidth', 1.5, 'MarkerSize', 6);
    loglog(N_list, avg_Rg2_fr, 'mv-', 'LineWidth', 1.5, 'MarkerSize', 6);
    xlabel('N'); ylabel('Value'); legend('FJ: <R^2>', 'FJ: <R_g^2>', 'FJ: Theory', 'FR: <R^2>', 'FR: <R_g^2>');
    grid on; title('Random Flight vs Freely Rotating Chain');
    
    figure;
    plot(N_list, avg_R2_fj ./ theory_fj, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
    plot(N_list, avg_R2_fr ./ theory_fj, 'r^-', 'LineWidth', 1.5, 'MarkerSize', 6);
    xlabel('N'); ylabel('<R^2> / (b^2 N)'); legend('FJ Model', 'FR Model (\theta=68^\circ)');
    ylim([0.8, 1.2]); grid on; title('Scaling Behavior');
    
    N_hist = 100;
    M_hist = 100000;
    R_vals = generate_freely_jointed_chain_hist(N_hist, M_hist);
    R = linspace(0, max(R_vals), 200);
    P_theory = 4*pi*R.^2 .* (3/(2*pi*N_hist))^(3/2) .* exp(-3*R.^2/(2*N_hist));
    
    figure;
    histogram(R_vals, 100, 'Normalization', 'pdf', 'EdgeColor', 'none'); hold on;
    plot(R, P_theory, 'r-', 'LineWidth', 1.5);
    xlabel('R'); ylabel('P(R)'); legend('Simulation', 'Theory'); grid on;
    title(sprintf('End-to-End Distance Distribution (N=%d)', N_hist));
end

function [avg_R2, avg_Rg2, theory] = freely_jointed(N_list, M)
    avg_R2 = zeros(size(N_list));
    avg_Rg2 = zeros(size(N_list));
    for i = 1:length(N_list)
        N = N_list(i);
        R2_sum = 0;
        Rg2_sum = 0;
        for j = 1:M(i)
            vecs = randn(N, 3);
            vecs = vecs ./ vecnorm(vecs, 2, 2);
            coords = [0,0,0; cumsum(vecs, 1)];
            R = coords(end,:);
            R2_sum = R2_sum + sum(R.^2);
            cm = mean(coords, 1);
            Rg2_sum = Rg2_sum + mean(sum((coords - cm).^2, 2));
        end
        avg_R2(i) = R2_sum / M(i);
        avg_Rg2(i) = Rg2_sum / M(i);
    end
    theory = 6*(N_list+1)./(N_list+2) .* avg_Rg2;
end

function [avg_R2, avg_Rg2] = freely_rotating(N_list, M, theta)
    avg_R2 = zeros(size(N_list));
    avg_Rg2 = zeros(size(N_list));
    for i = 1:length(N_list)
        N = N_list(i);
        R2_sum = 0;
        Rg2_sum = 0;
        for j = 1:M(i)
            u1 = randn(1,3);
            u1 = u1 / norm(u1);
            vecs = zeros(N,3);
            vecs(1,:) = u1;
            for k = 2:N
                v_rand = randn(1,3);
                v_perp = v_rand - dot(v_rand, vecs(k-1,:)) * vecs(k-1,:);
                v_perp = v_perp / norm(v_perp);
                vecs(k,:) = cos(theta)*vecs(k-1,:) + sin(theta)*v_perp;
            end
            coords = [0,0,0; cumsum(vecs,1)];
            R = coords(end,:);
            R2_sum = R2_sum + sum(R.^2);
            cm = mean(coords,1);
            Rg2_sum = Rg2_sum + mean(sum((coords - cm).^2,2));
        end
        avg_R2(i) = R2_sum / M(i);
        avg_Rg2(i) = Rg2_sum / M(i);
    end
end

function R_vals = generate_freely_jointed_chain_hist(N, M)
    R_vals = zeros(M,1);
    for i = 1:M
        vecs = randn(N,3);
        vecs = vecs ./ vecnorm(vecs,2,2);
        R_vec = sum(vecs,1);
        R_vals(i) = norm(R_vec);
    end
end