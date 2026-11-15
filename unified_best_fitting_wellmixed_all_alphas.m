clc
clear

%% ============================================================
%  1. Load simulation data
%% ============================================================

load('polarization_wellmixed.mat', ...
    'polarization','alpha_values','initial_polarization');

bestParameter=[1,0.3,0.5,1,1.47,2.21,2.45];

N = 100;

ratio = [0,initial_polarization/N,1];

alpha_values = alpha_values(:);

num_alpha = length(alpha_values);


%% ============================================================
%  2. Fixed payoff matrix
%% ============================================================

besta = bestParameter(1);
bestb = bestParameter(2);
bestc = bestParameter(3);
bestd = bestParameter(4);

wList = 0:0.01:10;

% Known alpha and omega
alpha_known = [2,5,10];

omega_known = bestParameter(5:7);

% Store optimized omega
bestw_fitting = NaN(num_alpha,1);

% Store SSE
bestLoss = NaN(num_alpha,1);


%% ============================================================
%  3. Store known omega
%% ============================================================

known_idx = zeros(1,length(alpha_known));

for k = 1:length(alpha_known)

    idx = find(alpha_values == alpha_known(k));

    known_idx(k) = idx;

    bestw_fitting(idx) = omega_known(k);

end


%% ============================================================
%  4. Fit remaining omega
%% ============================================================

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

    bestw_fitting(j) = bestw1;
    bestLoss(j) = bestLoss1;

end


%% ============================================================
%  5. Display optimized omega
%% ============================================================

result_table = table( ...
    alpha_values, ...
    bestw_fitting, ...
    'VariableNames',{'Alpha','Omega'});

disp(result_table);


%% ============================================================
%  6. Fit alpha-omega relationship
%% ============================================================

ft = fittype( ...
    'b*log(1+a*(x-1))', ...
    'independent','x', ...
    'coefficients',{'a','b'});

[f,gof] = fit( ...
    alpha_values, ...
    bestw_fitting, ...
    ft, ...
    'StartPoint',[10,0.5], ...
    'Lower',[0,0]);

a_fit = f.a;
b_fit = f.b;

fprintf('\n====================================\n');

fprintf('Fitted alpha-omega relationship:\n');

fprintf('omega(alpha) = %.8f * log(1 + %.8f*(alpha-1))\n', ...
    b_fit,a_fit);

fprintf('R-square = %.8f\n',gof.rsquare);

fprintf('RMSE = %.8f\n',gof.rmse);

fprintf('====================================\n');


%% ============================================================
%  7. Plot alpha-omega relationship
%% ============================================================

alpha_plot = linspace(1,10,500);

omega_plot = b_fit*log(1+a_fit*(alpha_plot-1));

subplot(2,2,1)
hold on;
box on;

plot(alpha_values,bestw_fitting, ...
    'o', ...
    'MarkerSize',8, ...
    'LineWidth',2, ...
    'DisplayName','Well-mixed: Optimized $\omega$');

plot(alpha_plot,omega_plot, ...
    '-', ...
    'LineWidth',2, ...
    'DisplayName','Well-mixed: Fitting line');

xlabel('Nonlinear exponent $\alpha$', ...
    'Interpreter','latex','FontSize',18);

ylabel('Selection intensity $\omega$', ...
    'Interpreter','latex','FontSize',18);

legend('show', ...
    'Interpreter','latex', ...
    'Location','best', ...
    'FontSize',16);

% Generate fitted equation
equationText = sprintf( ...
    '$\\omega(\\alpha)=%.4f\\ln[1+%.4f(\\alpha-1)]$', ...
    b_fit,a_fit);

% Display equation on the figure
text(0.05,0.95,equationText, ...
    'Units','normalized', ...
    'Interpreter','latex', ...
    'FontSize',16, ...
    'VerticalAlignment','top');

xlim([1,10]);

set(gca,'FontSize',16);

grid on;


for j = 1:num_alpha

    phi = fixation_curve( ...
        besta,bestb,bestc,bestd,bestw_fitting(j));

    phi_fit(j,:) = phi(:)';

end

residual = polarization - phi_fit;
RMSE = sqrt(sum((residual.^2)')/51);
Maxerror = max(abs(residual'));
subplot(2,2,2)
hold on;
box on;

plot(alpha_values,RMSE,'-o', ...
    'LineWidth',2, ...
    'MarkerSize',8, ...
    'DisplayName','Well-mixed: RMSE');

plot(alpha_values,Maxerror,'-s', ...
    'LineWidth',2, ...
    'MarkerSize',8, ...
    'DisplayName','Well-mixed: Max error');

xlabel('Nonlinear exponent $\alpha$','Interpreter','latex','FontSize',18);

ylabel('Fitting error','Interpreter','latex','FontSize',18);

legend('show', ...
    'Interpreter','latex', ...
    'Location','best', ...
    'FontSize',12);

grid on;


set(gca,'FontSize',16);

xlim([1,10]);
xticks(1:10);