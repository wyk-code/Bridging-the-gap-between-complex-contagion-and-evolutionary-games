clc;
clear;

%% ============================================================
%  1. Network settings
%% ============================================================

network_names = {'ER','WS','BA'};

file_names = { ...
    'polarization_ER.mat', ...
    'polarization_WS.mat', ...
    'polarization_BA.mat'};


num_network = length(network_names);

% Known alpha values
alpha_known = [2,4,6];

% Known omega values
% Rows: ER, WS, BA
% Columns: alpha = 2, 4, 6

bestOmega = [ ...
    2.74, 3.84, 4.17; ...
    2.56, 3.30, 3.55; ...
    2.00, 2.83, 2.94];

% Fixed payoff matrix
bestParameter = [0.2,0.3,0,0.6];

besta = bestParameter(1);
bestb = bestParameter(2);
bestc = bestParameter(3);
bestd = bestParameter(4);

N = 100;

wList = 0:0.01:10;

% Colors for ER, WS, BA
network_colors = [ ...
    0.0000, 0.4470, 0.7410; ...
    0.8500, 0.3250, 0.0980; ...
    0.4660, 0.6740, 0.1880];

% Markers for ER, WS, BA
network_markers = {'o','s','^'};


%% ============================================================
%  2. Load first network to initialize
%% ============================================================

S = load(file_names{1});

alpha_values = S.alpha_values(:);

num_alpha = length(alpha_values);

bestw_fitting = NaN(num_network,num_alpha);

RMSE_all = NaN(num_network,num_alpha);

Maxerror_all = NaN(num_network,num_alpha);

a_fit_all = NaN(num_network,1);
b_fit_all = NaN(num_network,1);

R2_all = NaN(num_network,1);


%% ============================================================
%  3. Loop over three networks
%% ============================================================

for ig = 1:num_network

    %% ========================================================
    %  Load simulation data
    %% ========================================================

    S = load(file_names{ig}, ...
        'polarization','alpha_values','initial_polarization');

    polarization = S.polarization;

    alpha_current = S.alpha_values(:);

    initial_polarization = S.initial_polarization;

    ratio = [0,initial_polarization(:)'/N,1];

    assert(isequal(alpha_current,alpha_values), ...
        'Alpha values differ between networks.');

    assert(size(polarization,1) == num_alpha && ...
           size(polarization,2) == length(ratio), ...
        'Polarization data dimensions do not match.');


    %% ========================================================
    %  Store known omega
    %% ========================================================

    known_idx = zeros(1,length(alpha_known));

    for k = 1:length(alpha_known)

        idx = find(alpha_values == alpha_known(k));

        assert(length(idx) == 1, ...
            'Known alpha value is missing or duplicated.');

        known_idx(k) = idx;

        bestw_fitting(ig,idx) = bestOmega(ig,k);

    end


    %% ========================================================
    %  Fit remaining omega
    %% ========================================================

    fit_idx = setdiff(1:num_alpha,known_idx);

    for j = fit_idx

        bestLoss1 = inf;
        bestw1 = NaN;

        current_data = polarization(j,:);

        for w1 = wList

            phi = fixation_curve( ...
                besta,bestb,bestc,bestd,w1);

            loss = sum((phi(:)-current_data(:)).^2);

            if loss < bestLoss1

                bestLoss1 = loss;
                bestw1 = w1;

            end

        end

        bestw_fitting(ig,j) = bestw1;

    end


    %% ========================================================
    %  Fit alpha-omega relationship
    %% ========================================================

    ft = fittype( ...
        'b*log(1+a*(x-1))', ...
        'independent','x', ...
        'coefficients',{'a','b'});

    [f,gof] = fit( ...
        alpha_values, ...
        bestw_fitting(ig,:)', ...
        ft, ...
        'StartPoint',[10,0.5], ...
        'Lower',[0,0]);

    a_fit_all(ig) = f.a;
    b_fit_all(ig) = f.b;

    R2_all(ig) = gof.rsquare;


    %% ========================================================
    %  Calculate polarization fitting errors
    %% ========================================================

    phi_fit = zeros(size(polarization));

    for j = 1:num_alpha

        phi = fixation_curve( ...
            besta,bestb,bestc,bestd, ...
            bestw_fitting(ig,j));

        phi_fit(j,:) = phi(:)';

    end

    residual = polarization - phi_fit;

    RMSE = sqrt(mean(residual.^2,2));

    Maxerror = max(abs(residual),[],2);

    RMSE_all(ig,:) = RMSE';
    Maxerror_all(ig,:) = Maxerror';


    %% ========================================================
    %  Display results
    %% ========================================================

    fprintf('\n====================================\n');

    fprintf('%s network\n',network_names{ig});

    fprintf('omega(alpha) = %.8f * log(1 + %.8f*(alpha-1))\n', ...
        b_fit_all(ig),a_fit_all(ig));

    fprintf('R-square = %.8f\n',gof.rsquare);

    fprintf('====================================\n');

    result_table = table( ...
        alpha_values, ...
        bestw_fitting(ig,:)', ...
        RMSE, ...
        Maxerror, ...
        'VariableNames', ...
        {'Alpha','Omega','RMSE','MaxError'});

    disp(result_table);

end


%% ============================================================
%  4. Plot three networks in one figure
%% ============================================================

% figure(1);
alpha_plot = linspace( ...
    min(alpha_values),max(alpha_values),500);


%% ============================================================
%  (a) Alpha-omega relationship
%% ============================================================

subplot(2,2,3);
hold on;
box on;

for ig = 1:num_network

    % Optimized omega
    plot(alpha_values,bestw_fitting(ig,:), ...
        'LineStyle','none', ...
        'Marker',network_markers{ig}, ...
        'Color',network_colors(ig,:), ...
        'MarkerSize',8, ...
        'LineWidth',2, ...
        'DisplayName',sprintf('%s: Optimized $\\omega$', ...
        network_names{ig}));

    % Fitted logarithmic curve
    omega_plot = b_fit_all(ig)*log( ...
        1+a_fit_all(ig)*(alpha_plot-1));

    plot(alpha_plot,omega_plot, ...
        '-', ...
        'Color',network_colors(ig,:), ...
        'LineWidth',2, ...
        'DisplayName',sprintf('%s: Fitting line', ...
        network_names{ig}));

end

xlabel('Nonlinear exponent $\alpha$', ...
    'Interpreter','latex','FontSize',18);

ylabel('Selection intensity $\omega$', ...
    'Interpreter','latex','FontSize',18);

legend('show', ...
    'Interpreter','latex', ...
    'Location','southeast', ...
    'FontSize',12);

xlim([min(alpha_values),max(alpha_values)]);
xticks(alpha_values);

set(gca,'FontSize',16);

grid on;


%% ============================================================
%  Display fitted equations
%% ============================================================

for ig = 1:num_network

    equationText = sprintf( ...
        '$\\omega_{\\mathrm{%s}}(\\alpha)=%.3f\\ln[1+%.3f(\\alpha-1)]$', ...
        network_names{ig}, ...
        b_fit_all(ig), ...
        a_fit_all(ig));

    text(0.03,0.95-(ig-1)*0.085, ...
        equationText, ...
        'Units','normalized', ...
        'Interpreter','latex', ...
        'FontSize',12, ...
        'BackgroundColor','w', ...
        'VerticalAlignment','top');

end


%% ============================================================
%  (b) Fitting errors
%% ============================================================

subplot(2,2,4);
hold on;
box on;

for ig = 1:num_network

    % RMSE: solid line
    plot(alpha_values,RMSE_all(ig,:), ...
        'LineStyle','-', ...
        'Marker',network_markers{ig}, ...
        'Color',network_colors(ig,:), ...
        'MarkerFaceColor',network_colors(ig,:), ...
        'LineWidth',2, ...
        'MarkerSize',7, ...
        'DisplayName',sprintf('%s: RMSE', ...
        network_names{ig}));

    % Maximum absolute error: dashed line
    plot(alpha_values,Maxerror_all(ig,:), ...
        'LineStyle','--', ...
        'Marker',network_markers{ig}, ...
        'Color',network_colors(ig,:), ...
        'LineWidth',2, ...
        'MarkerSize',7, ...
        'DisplayName',sprintf('%s: Max error', ...
        network_names{ig}));

end

xlabel('Nonlinear exponent $\alpha$', ...
    'Interpreter','latex','FontSize',18);

ylabel('Fitting error', ...
    'Interpreter','latex','FontSize',18);

legend('show', ...
    'Interpreter','latex', ...
    'Location','best', ...
    'FontSize',12);

xlim([min(alpha_values),max(alpha_values)]);
xticks(alpha_values);

set(gca,'FontSize',16);

grid on;


