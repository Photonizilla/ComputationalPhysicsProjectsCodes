function random_walk_3d
    rng(42);
    N_list_trad = [10,20,50,100,200,500,1000];
    num_chains_trad = 1000;
    [avg_re2_trad, avg_rg2_trad] = traditional_walk(N_list_trad, num_chains_trad);
    
    N_list_saw = [10,20,30,40,50];
    num_chains_saw = 200;
    max_attempts = 10000;
    [avg_re2_saw, avg_rg2_saw] = saw_walk(N_list_saw, num_chains_saw, max_attempts);
    
    analyze_results(N_list_trad, avg_re2_trad, avg_rg2_trad, 'Traditional');
    analyze_results(N_list_saw, avg_re2_saw, avg_rg2_saw, 'Self-Avoiding');
    
    check_transience(10000, 1000);
end

function [avg_re2, avg_rg2] = traditional_walk(N_list, num_chains)
    avg_re2 = zeros(size(N_list));
    avg_rg2 = zeros(size(N_list));
    directions = [1,0,0; -1,0,0; 0,1,0; 0,-1,0; 0,0,1; 0,0,-1];
    for i = 1:length(N_list)
        N = N_list(i);
        re2_sum = 0;
        rg2_sum = 0;
        for c = 1:num_chains
            pos = [0,0,0];
            positions = zeros(N+1,3);
            positions(1,:) = pos;
            for step = 1:N
                move = directions(randi(6),:);
                pos = pos + move;
                positions(step+1,:) = pos;
            end
            re2 = sum(pos.^2);
            cm = mean(positions,1);
            rg2 = sum(sum((positions - cm).^2,2)) / (N+1);
            re2_sum = re2_sum + re2;
            rg2_sum = rg2_sum + rg2;
        end
        avg_re2(i) = re2_sum / num_chains;
        avg_rg2(i) = rg2_sum / num_chains;
    end
end

function [avg_re2, avg_rg2] = saw_walk(N_list, num_chains, max_attempts)
    avg_re2 = zeros(size(N_list));
    avg_rg2 = zeros(size(N_list));
    directions = [1,0,0; -1,0,0; 0,1,0; 0,-1,0; 0,0,1; 0,0,-1];
    for i = 1:length(N_list)
        N = N_list(i);
        re2_sum = 0;
        rg2_sum = 0;
        chains_generated = 0;
        attempts = 0;
        while chains_generated < num_chains && attempts < max_attempts
            attempts = attempts + 1;
            pos = [0,0,0];
            positions = zeros(N+1,3);
            positions(1,:) = pos;
            visited = containers.Map(mat2str(pos), true);
            dead = false;
            for step = 1:N
                dirs_rand = directions(randperm(6),:);
                moved = false;
                for d = 1:6
                    new_pos = pos + dirs_rand(d,:);
                    key = mat2str(new_pos);
                    if ~isKey(visited, key)
                        pos = new_pos;
                        positions(step+1,:) = pos;
                        visited(key) = true;
                        moved = true;
                        break;
                    end
                end
                if ~moved
                    dead = true;
                    break;
                end
            end
            if ~dead
                re2 = sum(pos.^2);
                cm = mean(positions,1);
                rg2 = sum(sum((positions - cm).^2,2)) / (N+1);
                re2_sum = re2_sum + re2;
                rg2_sum = rg2_sum + rg2;
                chains_generated = chains_generated + 1;
            end
        end
        if chains_generated > 0
            avg_re2(i) = re2_sum / chains_generated;
            avg_rg2(i) = rg2_sum / chains_generated;
        else
            avg_re2(i) = NaN;
            avg_rg2(i) = NaN;
        end
    end
end

function analyze_results(N_list, avg_re2, avg_rg2, walk_type)
    valid_idx = ~isnan(avg_re2) & ~isnan(avg_rg2);
    N_list = N_list(valid_idx);
    avg_re2 = avg_re2(valid_idx);
    avg_rg2 = avg_rg2(valid_idx);
    logN = log(N_list);
    logRe2 = log(avg_re2);
    logRg2 = log(avg_rg2);
    fit_re = polyfit(logN, logRe2, 1);
    fit_rg = polyfit(logN, logRg2, 1);
    nu_re = fit_re(1)/2;
    nu_rg = fit_rg(1)/2;
    ratio = avg_rg2 ./ avg_re2;
    
    figure;
    subplot(2,1,1);
    loglog(N_list, avg_re2, 'o-', 'DisplayName', '\langle r_e^2 \rangle');
    hold on;
    loglog(N_list, exp(fit_re(2)) * N_list.^(2*nu_re), '--', 'DisplayName', sprintf('Fit: \\nu=%.3f', nu_re));
    loglog(N_list, avg_rg2, 's-', 'DisplayName', '\langle r_g^2 \rangle');
    loglog(N_list, exp(fit_rg(2)) * N_list.^(2*nu_rg), '--', 'DisplayName', sprintf('Fit: \\nu=%.3f', nu_rg));
    xlabel('N');
    ylabel('Value');
    title([walk_type ' Random Walk: Scaling']);
    legend;
    grid on;
    
    subplot(2,1,2);
    plot(N_list, ratio, 'd-');
    xlabel('N');
    ylabel('\langle r_g^2 \rangle / \langle r_e^2 \rangle');
    title('Ratio of Radius of Gyration to End-to-End Distance');
    grid on;
    
    fprintf('%s Walk:\n', walk_type);
    fprintf('  From <r_e^2>: ν = %.4f\n', nu_re);
    fprintf('  From <r_g^2>: ν = %.4f\n', nu_rg);
    fprintf('  Mean ratio <r_g^2>/<r_e^2> = %.4f\n', mean(ratio));
end

function check_transience(max_steps, num_chains)
    directions = [1,0,0; -1,0,0; 0,1,0; 0,-1,0; 0,0,1; 0,0,-1];
    returned = 0;
    for c = 1:num_chains
        pos = [0,0,0];
        visited = containers.Map(mat2str(pos), true);
        for step = 1:max_steps
            move = directions(randi(6),:);
            pos = pos + move;
            key = mat2str(pos);
            if isKey(visited, key) && all(pos == [0,0,0])
                returned = returned + 1;
                break;
            end
            visited(key) = true;
        end
    end
    p_return = returned / num_chains;
    fprintf('Traditional Walk Transience Check:\n');
    fprintf('  Returned to origin in %d steps: %.4f\n', max_steps, p_return);
    fprintf('  Transient since p_return < 1: %d\n', p_return < 1);
end
