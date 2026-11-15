clc;
clear;

%% ============================================================
%  1. Load simulation data
%% ============================================================

network_names = {'ER','WS','BA'};

file_names = { ...
    'polarization_ER.mat', ...
    'polarization_WS.mat', ...
    'polarization_BA.mat'};

alpha_selected = [2,4,6];

num_network = 3;
num_alpha = length(alpha_selected);

N = 100;

% Load the first network
S = load(file_names{1});

initial_polarization = S.initial_polarization;

ratio = [0,initial_polarization/N,1];

num_points = length(ratio);

data = zeros(num_network,num_alpha,num_points);

for ig = 1:num_network

    S = load(file_names{ig});

    [found,idx] = ismember(alpha_selected,S.alpha_values);

    assert(all(found), ...
        'Some alpha values are missing in %s.',file_names{ig});

    assert(isequal(S.initial_polarization, ...
        initial_polarization), ...
        'Initial polarization differs between networks.');

    data(ig,:,:) = S.polarization(idx,:);

end

fprintf('All three network datasets loaded.\n');


%% ============================================================
%  2. Search ranges
%% ============================================================

% Common game parameters
value = 0:0.1:1;

% Independent omega for each network and alpha
wList = 0:0.01:6;


%% ============================================================
%  3. Storage
%% ============================================================

bestLoss = inf;

bestParameter = zeros(1,4);

% Rows: ER, WS, BA
% Columns: alpha = 1,...,6
bestOmega = zeros(num_network,num_alpha);


%% ============================================================
%  4. Joint grid search
%% ============================================================

tic;

fprintf('\n========== Joint grid search ==========\n');

for a = value
    for b = value
        for c = value
            for d = value

                totalLoss = 0;

                omega_temp = zeros(num_network,num_alpha);

                %% Search omega independently for 18 datasets

                for ig = 1:num_network

                    for ia = 1:num_alpha

                        current_data = ...
                            squeeze(data(ig,ia,:))';

                        bestLoss_alpha = inf;

                        bestw = NaN;

                        for w = wList

                            % Neutral selection
                            if w == 0

                                phi = ratio;

                            else

                                phi = fixation_curve( ...
                                    a,b,c,d,w);

                            end

                            residual = phi(:)'-current_data;

                            loss = sum(residual.^2);

                            if loss < bestLoss_alpha

                                bestLoss_alpha = loss;

                                bestw = w;

                            end

                        end

                        omega_temp(ig,ia) = bestw;

                        totalLoss = totalLoss+bestLoss_alpha;

                    end

                end

                %% Update global optimum

                if totalLoss < bestLoss

                    bestLoss = totalLoss;

                    bestParameter = [a,b,c,d];

                    bestOmega = omega_temp;

                end

            end
        end
    end
end

toc;


%% ============================================================
%  5. Optimal common game parameters
%% ============================================================

besta = bestParameter(1);
bestb = bestParameter(2);
bestc = bestParameter(3);
bestd = bestParameter(4);

fprintf('\n========================================\n');
fprintf('       Optimal common game parameters\n');
fprintf('========================================\n');

fprintf('a = %.4f\n',besta);
fprintf('b = %.4f\n',bestb);
fprintf('c = %.4f\n',bestc);
fprintf('d = %.4f\n',bestd);

fprintf('\nTotal SSE = %.12g\n',bestLoss);


%% ============================================================
%  6. Optimal selection intensity
%% ============================================================

fprintf('\n========================================\n');
fprintf('       Optimal selection intensity\n');
fprintf('========================================\n');

for ig = 1:num_network

    fprintf('\n%s network:\n',network_names{ig});

    for ia = 1:num_alpha

        fprintf('alpha = %d, omega* = %.4f\n', ...
            alpha_selected(ia),bestOmega(ig,ia));

    end

end


%% ============================================================
%  7. Fitting curves and residuals
%% ============================================================

phi_fit = zeros(num_network,num_alpha,num_points);

residual_all = zeros(num_network,num_alpha,num_points);

SSE = zeros(num_network,num_alpha);

RMSE = zeros(num_network,num_alpha);

MAX = zeros(num_network,num_alpha);

for ig = 1:num_network

    for ia = 1:num_alpha

        current_data = squeeze(data(ig,ia,:))';

        w = bestOmega(ig,ia);

        if w == 0

            phi = ratio;

        else

            phi = fixation_curve( ...
                besta,bestb,bestc,bestd,w);

        end

        phi_fit(ig,ia,:) = phi(:)';

        residual_all(ig,ia,:) = ...
            phi_fit(ig,ia,:)-reshape(current_data,1,1,[]);

        current_residual = ...
            squeeze(residual_all(ig,ia,:));

        SSE(ig,ia) = sum(current_residual.^2);

        RMSE(ig,ia) = sqrt(mean(current_residual.^2));

        MAX(ig,ia) = max(abs(current_residual));

    end

end


%% ============================================================
%  8. Residual analysis
%% ============================================================

fprintf('\n========================================\n');
fprintf('             Residual analysis\n');
fprintf('========================================\n');

for ig = 1:num_network

    fprintf('\n========== %s network ==========\n', ...
        network_names{ig});

    for ia = 1:num_alpha

        fprintf('\nAlpha = %d:\n',alpha_selected(ia));

        fprintf('Omega* = %.4f\n',bestOmega(ig,ia));

        fprintf('SSE  = %.10g\n',SSE(ig,ia));

        fprintf('RMSE = %.10g\n',RMSE(ig,ia));

        fprintf('Max residual = %.10g\n',MAX(ig,ia));

    end

    fprintf('\nTotal SSE = %.12g\n',sum(SSE(ig,:)));

end

fprintf('\nAll-network SSE = %.12g\n',sum(SSE(:)));


%% ============================================================
%  9. Plot omega-alpha relationship
%% ============================================================

hold on;
box on;

colors = lines(num_network);

for ig = 1:num_network

    plot(alpha_selected,bestOmega(ig,:), ...
        '-o', ...
        'Color',colors(ig,:), ...
        'LineWidth',2, ...
        'MarkerFaceColor',colors(ig,:), ...
        'MarkerSize',7);

end

xlabel('$\alpha$', ...
    'Interpreter','latex','FontSize',18);

ylabel('$\omega^*$', ...
    'Interpreter','latex','FontSize',18);

legend(network_names,'Location','best');

xlim([1,6]);

set(gca,'FontSize',16);


%% ============================================================
%  10. Plot fitting curves for each network
%% ============================================================

alpha_colors = lines(num_alpha);

markers = {'o','s','^'};

for ig = 1:num_network

    subplot(2,3,ig);
    hold on;
    box on;

    % ========================================================
    % 1. Plot simulation results
    % ========================================================

    for ia = 1:num_alpha

        current_data = squeeze(data(ig,ia,:))';

        plot(ratio,current_data, ...
            markers{ia}, ...
            'Color',alpha_colors(ia,:), ...
            'LineWidth',2, ...
            'MarkerSize',6, ...
            'LineStyle','none', ...
            'DisplayName',sprintf( ...
            'Simulation: $\\alpha=%g$',alpha_selected(ia)));

    end

    % ========================================================
    % 2. Plot fitting results
    % ========================================================

    for ia = 1:num_alpha

        current_fit = squeeze(phi_fit(ig,ia,:))';

        plot(ratio,current_fit, ...
            '-', ...
            'Color',alpha_colors(ia,:), ...
            'LineWidth',1.5, ...
            'DisplayName',sprintf( ...
            'Fit: $\\alpha=%g$, $\\omega=%.2f$', ...
            alpha_selected(ia),bestOmega(ig,ia)));

    end

    % ========================================================
    % 3. Figure settings
    % ========================================================

    xlabel('Initial fraction of extreme opinions', ...
        'Interpreter','latex','FontSize',18);

    ylabel('Polarization/Fixation probability', ...
        'Interpreter','latex','FontSize',18);

    legend('show', ...
        'Interpreter','latex', ...
        'Location','best', ...
        'FontSize',13);

    title(sprintf('%s network',network_names{ig}));

    xlim([0,1]);
    ylim([0,1]);

    grid on;

    set(gca,'FontSize',16);

end


%% ============================================================
%  11. Plot residuals for each network
%% ============================================================

for ig = 1:num_network

    subplot(2,3,ig+3);
    hold on;
    box on;

    for ia = 1:num_alpha

        current_residual = ...
            squeeze(residual_all(ig,ia,:))';

        plot(ratio,current_residual, ...
           'LineStyle','-', ...
           'Marker',markers{ia}, ...
            'Color',alpha_colors(ia,:), ...
            'LineWidth',1.5, ...
            'MarkerSize',6);

    end

    yline(0,'k--');
    
    xlabel('Initial fraction of extreme opinions', ...
    'Interpreter','latex','FontSize',18);
    
    ylabel('Residual ${\rho }_{\alpha }^{fit}-P_{\alpha }^{sim}$', ...
    'Interpreter','latex','FontSize',18);

    legend(arrayfun(@(x) sprintf('\\alpha = %d',x), ...
        alpha_selected,'UniformOutput',false), ...
        'Location','best');

    title(sprintf('%s network',network_names{ig}));

    set(gca,'FontSize',14);

end


%% ============================================================
%  12. Save results
%% ============================================================

save('joint_network_fitting_results.mat', ...
    'bestParameter','bestOmega','bestLoss', ...
    'SSE','RMSE','MAX', ...
    'phi_fit','residual_all', ...
    'alpha_selected','ratio','data', ...
    'network_names');

% Export omega table

omega_table = table( ...
    alpha_selected(:), ...
    bestOmega(1,:)', ...
    bestOmega(2,:)', ...
    bestOmega(3,:)', ...
    'VariableNames', ...
    {'Alpha','Omega_ER','Omega_WS','Omega_BA'});

writetable(omega_table, ...
    'joint_network_omega_results.xlsx');

fprintf('\nAll results saved successfully.\n');
