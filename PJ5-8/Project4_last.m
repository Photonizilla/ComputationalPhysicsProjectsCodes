rng(0);
L_trad = [10,20,40,80,160,320];
L_saw = [4,8,16,32,64];
M_trad = 1000;
M_saw = 500;
max_trials = 10000;

results_cont = zeros(length(L_trad),3);
for i = 1:length(L_trad)
    L = L_trad(i);
    r_e2_sum = 0;
    r_g2_sum = 0;
    for j = 1:M_trad
        pos = [0,0];
        path = zeros(L+1,2);
        path(1,:) = pos;
        for k = 1:L
            theta = 2*pi*rand();
            pos = pos + [cos(theta), sin(theta)];
            path(k+1,:) = pos;
        end
        r_e2 = sum((path(end,:) - path(1,:)).^2);
        cm = mean(path,1);
        diff = path - cm;
        r_g2 = mean(sum(diff.^2,2));
        r_e2_sum = r_e2_sum + r_e2;
        r_g2_sum = r_g2_sum + r_g2;
    end
    results_cont(i,1) = r_e2_sum / M_trad;
    results_cont(i,2) = r_g2_sum / M_trad;
    results_cont(i,3) = results_cont(i,2) / results_cont(i,1);
end

results_latt = zeros(length(L_trad),3);
for i = 1:length(L_trad)
    L = L_trad(i);
    r_e2_sum = 0;
    r_g2_sum = 0;
    dirs = [1,0; -1,0; 0,1; 0,-1];
    for j = 1:M_trad
        pos = [0,0];
        path = zeros(L+1,2);
        path(1,:) = pos;
        for k = 1:L
            d = dirs(randi(4),:);
            pos = pos + d;
            path(k+1,:) = pos;
        end
        r_e2 = sum((path(end,:) - path(1,:)).^2);
        cm = mean(path,1);
        diff = path - cm;
        r_g2 = mean(sum(diff.^2,2));
        r_e2_sum = r_e2_sum + r_e2;
        r_g2_sum = r_g2_sum + r_g2;
    end
    results_latt(i,1) = r_e2_sum / M_trad;
    results_latt(i,2) = r_g2_sum / M_trad;
    results_latt(i,3) = results_latt(i,2) / results_latt(i,1);
end

results_saw = zeros(length(L_saw),3);
for i = 1:length(L_saw)
    L = L_saw(i);
    r_e2_sum = 0;
    r_g2_sum = 0;
    count = 0;
    trials = 0;
    dirs = [1,0; -1,0; 0,1; 0,-1];
    while count < M_saw && trials < max_trials
        pos = [0,0];
        path = zeros(L+1,2);
        path(1,:) = pos;
        visited = containers.Map();
        visited(num2str(pos)) = 1;
        success = true;
        for k = 1:L
            next_pos = pos + dirs;
            avail = [];
            for d = 1:4
                key = num2str(next_pos(d,:));
                if ~isKey(visited, key)
                    avail = [avail; next_pos(d,:)];
                end
            end
            if isempty(avail)
                success = false;
                break;
            end
            idx = randi(size(avail,1));
            pos = avail(idx,:);
            path(k+1,:) = pos;
            visited(num2str(pos)) = 1;
        end
        if success
            r_e2 = sum((path(end,:) - path(1,:)).^2);
            cm = mean(path,1);
            diff = path - cm;
            r_g2 = mean(sum(diff.^2,2));
            r_e2_sum = r_e2_sum + r_e2;
            r_g2_sum = r_g2_sum + r_g2;
            count = count + 1;
        end
        trials = trials + 1;
    end
    if count > 0
        results_saw(i,1) = r_e2_sum / count;
        results_saw(i,2) = r_g2_sum / count;
        results_saw(i,3) = results_saw(i,2) / results_saw(i,1);
    end
end

logL_cont = log(L_trad);
logRe_cont = log(results_cont(:,1));
logRg_cont = log(results_cont(:,2));
fit_e_cont = polyfit(logL_cont, logRe_cont, 1);
fit_g_cont = polyfit(logL_cont, logRg_cont, 1);
nu_e_cont = fit_e_cont(1)/2;
nu_g_cont = fit_g_cont(1)/2;

logL_latt = log(L_trad);
logRe_latt = log(results_latt(:,1));
logRg_latt = log(results_latt(:,2));
fit_e_latt = polyfit(logL_latt, logRe_latt, 1);
fit_g_latt = polyfit(logL_latt, logRg_latt, 1);
nu_e_latt = fit_e_latt(1)/2;
nu_g_latt = fit_g_latt(1)/2;

logL_saw = log(L_saw);
logRe_saw = log(results_saw(:,1));
logRg_saw = log(results_saw(:,2));
fit_e_saw = polyfit(logL_saw, logRe_saw, 1);
fit_g_saw = polyfit(logL_saw, logRg_saw, 1);
nu_e_saw = fit_e_saw(1)/2;
nu_g_saw = fit_g_saw(1)/2;

figure;
subplot(1,3,1);
loglog(L_trad, results_cont(:,1), 'o'); hold on;
loglog(L_trad, results_cont(:,2), 's');
loglog(L_trad, exp(fit_e_cont(2)) * L_trad.^fit_e_cont(1), '-');
loglog(L_trad, exp(fit_g_cont(2)) * L_trad.^fit_g_cont(1), '-');
title('Continuous Traditional RW');
xlabel('L'); ylabel('<r^2>');
legend('<r_e^2>', '<r_g^2>', 'Fit e', 'Fit g');

subplot(1,3,2);
loglog(L_trad, results_latt(:,1), 'o'); hold on;
loglog(L_trad, results_latt(:,2), 's');
loglog(L_trad, exp(fit_e_latt(2)) * L_trad.^fit_e_latt(1), '-');
loglog(L_trad, exp(fit_g_latt(2)) * L_trad.^fit_g_latt(1), '-');
title('Lattice Traditional RW');
xlabel('L'); ylabel('<r^2>');

subplot(1,3,3);
loglog(L_saw, results_saw(:,1), 'o'); hold on;
loglog(L_saw, results_saw(:,2), 's');
loglog(L_saw, exp(fit_e_saw(2)) * L_saw.^fit_e_saw(1), '-');
loglog(L_saw, exp(fit_g_saw(2)) * L_saw.^fit_g_saw(1), '-');
title('Lattice SAW');
xlabel('L'); ylabel('<r^2>');

disp('Continuous Traditional RW:');
disp(['nu_e = ', num2str(nu_e_cont), ', nu_g = ', num2str(nu_g_cont)]);
disp(['Mean ratio r_g^2/r_e^2 = ', num2str(mean(results_cont(:,3)))]);

disp('Lattice Traditional RW:');
disp(['nu_e = ', num2str(nu_e_latt), ', nu_g = ', num2str(nu_g_latt)]);
disp(['Mean ratio r_g^2/r_e^2 = ', num2str(mean(results_latt(:,3)))]);

disp('Lattice SAW:');
disp(['nu_e = ', num2str(nu_e_saw), ', nu_g = ', num2str(nu_g_saw)]);
disp(['Mean ratio r_g^2/r_e^2 = ', num2str(mean(results_saw(:,3)))]);
