% clc;
clear;

%% ==================== Parameters ====================

% Select any three alpha values
alpha_selected = [2,5,10];

% Search ranges
value = 0:0.1:1;

wList = 0:0.01:10;

%% ==================== Load simulation data ====================

load('polarization_wellmixed.mat', ...
    'polarization','alpha_values','initial_polarization');

% Find the indices of selected alpha values
[found,idx] = ismember(alpha_selected,alpha_values);

% Extract the three datasets
data1 = polarization(idx(1),:);
data2 = polarization(idx(2),:);
data3 = polarization(idx(3),:);

% Initial polarization fraction
N = 100;

ratio = [0,initial_polarization/N,1];


%% ==================== Storage ====================

bestLoss = inf;

bestParameter = zeros(1,7);

%% ==================== Grid search ====================

tic;

for a = value
    for b = value
        for c = value
            for d = value

                %% ---------- First alpha ----------

                bestLoss1 = inf;
                bestw1 = NaN;

                for w1 = wList

                    phi = fixation_curve(a,b,c,d,w1);

                    loss = sum((phi(:)-data1(:)).^2);

                    if loss < bestLoss1

                        bestLoss1 = loss;
                        bestw1 = w1;

                    end

                end

                %% ---------- Second alpha ----------

                bestLoss2 = inf;
                bestw2 = NaN;

                for w2 = wList

                    phi = fixation_curve(a,b,c,d,w2);

                    loss = sum((phi(:)-data2(:)).^2);

                    if loss < bestLoss2

                        bestLoss2 = loss;
                        bestw2 = w2;

                    end

                end

                %% ---------- Third alpha ----------

                bestLoss3 = inf;
                bestw3 = NaN;

                for w3 = wList

                    phi = fixation_curve(a,b,c,d,w3);

                    loss = sum((phi(:)-data3(:)).^2);

                    if loss < bestLoss3

                        bestLoss3 = loss;
                        bestw3 = w3;

                    end

                end

                %% ---------- Total error ----------

                totalLoss = bestLoss1+bestLoss2+bestLoss3;

                if totalLoss < bestLoss

                    bestLoss = totalLoss;

                    bestParameter = ...
                        [a,b,c,d,bestw1,bestw2,bestw3];

                end

            end
        end
    end
end

toc;

%% ==================== Optimal parameters ====================

besta = bestParameter(1);
bestb = bestParameter(2);
bestc = bestParameter(3);
bestd = bestParameter(4);

bestw = bestParameter(5:7);

fprintf('\n========== Optimal parameters ==========\n');

fprintf('a = %.4f\n',besta);
fprintf('b = %.4f\n',bestb);
fprintf('c = %.4f\n',bestc);
fprintf('d = %.4f\n',bestd);

for k = 1:3

    fprintf('alpha = %g, omega = %.4f\n', ...
        alpha_selected(k),bestw(k));

end

fprintf('Total SSE = %.12g\n',bestLoss);

%% ==================== Fitted curves ====================

data = {data1,data2,data3};

phi_fit = cell(1,3);
residual = cell(1,3);

SSE = zeros(1,3);
RMSE = zeros(1,3);
MAX = zeros(1,3);

for k = 1:3

    phi_fit{k} = fixation_curve( ...
        besta,bestb,bestc,bestd,bestw(k));

    phi_fit{k} = phi_fit{k}(:)';

    residual{k} = phi_fit{k}-data{k};

    SSE(k) = sum(residual{k}.^2);

    RMSE(k) = sqrt(mean(residual{k}.^2));

    MAX(k) = max(abs(residual{k}));

end

%% ==================== Residual analysis ====================

fprintf('\n========== Residual analysis ==========\n');

for k = 1:3

    fprintf('\nAlpha = %g:\n',alpha_selected(k));

    fprintf('SSE  = %.10g\n',SSE(k));
    fprintf('RMSE = %.10g\n',RMSE(k));
    fprintf('Max absolute residual = %.10g\n',MAX(k));

end

fprintf('\nTotal SSE = %.10g\n',sum(SSE));

%% ==================== Plot fitting curves ====================

% figure;
subplot(1,2,1)
hold on;
box on;

colors = lines(3);

markers = {'o','s','^'};

legendText = cell(1,6);

for k = 1:3

    plot(ratio,data{k},markers{k}, ...
        'Color',colors(k,:), ...
        'MarkerSize',6, ...
        'LineStyle','none', ...
        'LineWidth',1.5);

    legendText{k} = sprintf( ...
        'Simulation: $\\alpha=%g$',alpha_selected(k));

end

for k = 1:3

    plot(ratio,phi_fit{k},'-', ...
        'Color',colors(k,:), ...
        'LineWidth',2);

    legendText{k+3} = sprintf( ...
        'Fit: $\\alpha=%g$, $\\omega=%.2f$', ...
        alpha_selected(k),bestw(k));

end

xlabel('Initial fraction of extreme opinions', ...
    'Interpreter','latex','FontSize',18);

ylabel('Polarization/Fixation probability', ...
    'Interpreter','latex','FontSize',18);

% legend(legendText, ...
%     'Interpreter','latex', ...
%     'FontSize',13, ...
%     'Location','best');

xlim([0,1]);
ylim([0,1]);

set(gca,'FontSize',16);

%% ==================== Plot residuals ====================

% figure;
subplot(1,2,2)
hold on;
box on;

legendResidual = cell(1,3);

for k = 1:3

    plot(ratio,residual{k}, ...
        'LineWidth',1.5, ...
        'Marker',markers{k}, ...
        'Color',colors(k,:));

    legendResidual{k} = sprintf( ...
        '$\\alpha=%g$',alpha_selected(k));

end

yline(0,'k--','LineWidth',1);

xlabel('Initial fraction of extreme opinions', ...
    'Interpreter','latex','FontSize',18);

ylabel('Residual ${\rho }_{\alpha }^{fit}-P_{\alpha }^{sim}$', ...
    'Interpreter','latex','FontSize',18);

legend(legendResidual, ...
    'Interpreter','latex', ...
    'Location','best');

set(gca,'FontSize',16);


% check codes
% value = 0:0.1:1;
% 
% wList = 0:0.01:3;
% 
% bestLoss = inf;
% 
% bestParameter = zeros(1,7);
% 
% for a = value
%     for b = value
%         for c = value
%             for d = value
% 
%                 % ===========================
%                 % 第一组数据
%                 % ===========================
% 
%                 bestLoss1 = inf;
% 
%                 bestw1 = NaN;
% 
%                 for w1 = wList
% 
%                     phi = fixation_curve(a,b,c,d,w1);
% 
%                     loss = sum((phi-polarization2).^2);
% 
%                     if loss<bestLoss1
% 
%                         bestLoss1 = loss;
% 
%                         bestw1 = w1;
% 
%                     end
% 
%                 end
% 
%                 % ===========================
%                 % 第二组
%                 % ===========================
% 
%                 bestLoss2 = inf;
% 
%                 bestw2 = NaN;
% 
%                 for w2 = wList
% 
%                     phi = fixation_curve(a,b,c,d,w2);
% 
%                     loss = sum((phi-polarization5).^2);
% 
%                     if loss<bestLoss2
% 
%                         bestLoss2 = loss;
% 
%                         bestw2 = w2;
% 
%                     end
% 
%                 end
% 
%                 % ===========================
%                 % 第三组
%                 % ===========================
% 
%                 bestLoss3 = inf;
% 
%                 bestw3 = NaN;
% 
%                 for w3 = wList
% 
%                     phi = fixation_curve(a,b,c,d,w3);
% 
%                     loss = sum((phi-polarization10).^2);
% 
%                     if loss<bestLoss3
% 
%                         bestLoss3 = loss;
% 
%                         bestw3 = w3;
% 
%                     end
% 
%                 end
% 
%                 % 总误差
% 
%                 totalLoss = bestLoss1 + bestLoss2 + bestLoss3;
% 
%                 if totalLoss<bestLoss
% 
%                     bestLoss = totalLoss;
% 
%                     bestParameter = [a b c d bestw1 bestw2 bestw3];
% 
%                 end
% 
%             end
%         end
%     end
% end
