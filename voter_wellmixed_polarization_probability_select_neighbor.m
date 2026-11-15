clc;
clear;

%% ==================== Parameters ====================

N = 100;

initial_polarization = 2:2:98;

alpha_values = 1:10;

num_runs = 5000;

num_alpha = length(alpha_values);
num_initial = length(initial_polarization);

%% ==================== Parallel pool ====================

delete(gcp('nocreate'));
parpool('local',10);

%% ==================== Storage ====================

x_polarization = zeros(num_alpha,num_runs,num_initial);

%% ==================== Monte Carlo ====================

tic;

for ia = 1:num_alpha

    alpha = alpha_values(ia);

    fprintf('Running alpha = %d\n',alpha);

    parfor j = 1:num_runs

        for j1 = 1:num_initial

            x_polarization(ia,j,j1) = ...
                voter_wellmixed_polarization_probability( ...
                alpha,initial_polarization(j1),N);

        end

    end

end

delete(gcp('nocreate'));

toc;

%% ==================== Mean polarization probability ====================

polarization = zeros(num_alpha,num_initial+2);

for ia = 1:num_alpha

    polarization(ia,:) = ...
        [0, squeeze(mean(x_polarization(ia,:,:),2))', 1];

end

save('polarization_wellmixed1214.mat','polarization','alpha_values','initial_polarization');

%% ==================== Plot ====================

figure;
hold on;
box on;

colors = lines(num_alpha);

for ia = 1:num_alpha

    plot([0,initial_polarization/N,1], ...
        polarization(ia,:), ...
        '-o', ...
        'Color',colors(ia,:), ...
        'LineWidth',1.5, ...
        'MarkerSize',5);

end

xlabel('Initial fraction of extreme opinions');
ylabel('Polarization probability');

legend(arrayfun(@(a) sprintf('\\alpha = %d',a), ...
    alpha_values,'UniformOutput',false), ...
    'Location','best');

xlim([0,1]);
ylim([0,1]);


%% ============================================================
%  Function: Conservative learning in well-mixed population
%% ============================================================

function x_polarization = ...
    voter_wellmixed_polarization_probability( ...
    alpha,initial_polarization,N)

    %% ==================== Initial opinions ====================

    size_A = randi(initial_polarization);
    size_B = initial_polarization-size_A;

    state = [ones(1,size_A), ...
             zeros(1,N-size_A-size_B), ...
             -ones(1,size_B)];

    state = state(randperm(N));

    %% ==================== Opinion counts ====================

    nPlus  = nnz(state==1);
    nMinus = nnz(state==-1);
    nZero  = nnz(state==0);

    %% ==================== Evolution ====================

    while true

        % Step 1: Randomly select a focal individual
        i = randi(N);

        % Step 2: Randomly select another individual
        % The focal individual cannot select itself
        j = randi(N-1);

        if j >= i
            j = j+1;
        end

        % Opinions of focal and selected individuals
        old_i = state(i);
        opinion_j = state(j);

        % Step 3: Calculate opinion fractions
        % Exclude the focal individual itself

        fA = (nPlus  - (old_i==1))/(N-1);
        fB = (nMinus - (old_i==-1))/(N-1);
        fC = (nZero  - (old_i==0))/(N-1);

        % Step 4: Conservative learning
        p = rand;

        % A -> C
        if old_i==1 && opinion_j==0

            if p <= fC^(alpha-1)

                state(i) = 0;

                nPlus = nPlus-1;
                nZero = nZero+1;

            end

        % B -> C
        elseif old_i==-1 && opinion_j==0

            if p <= fC^(alpha-1)

                state(i) = 0;

                nMinus = nMinus-1;
                nZero = nZero+1;

            end

        % C -> A
        elseif old_i==0 && opinion_j==1

            if p <= fA^(alpha-1)

                state(i) = 1;

                nZero = nZero-1;
                nPlus = nPlus+1;

            end

        % C -> B
        elseif old_i==0 && opinion_j==-1

            if p <= fB^(alpha-1)

                state(i) = -1;

                nZero = nZero-1;
                nMinus = nMinus+1;

            end

        end

        % A and B cannot directly learn from each other

        %% ==================== Absorbing states ====================

        if nZero==N || nZero==0

            % 1: Polarization edge
            % 0: Neutral consensus
            x_polarization = 1-nZero/N;

            break;

        end

    end

end