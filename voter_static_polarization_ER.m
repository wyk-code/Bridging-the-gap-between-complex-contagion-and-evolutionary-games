clc;
clear;

%% ==================== Parameters ====================

N = 100;

initial_polarization = 2:2:98;

alpha_values = 1:0.5:6;
% alpha_values=[2,4,6];
% alpha_values = [6.5,7,7.5,8];

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

    fprintf('Running alpha = %0.2f\n',alpha);

    parfor j = 1:num_runs

        for j1 = 1:num_initial

            x_polarization(ia,j,j1) = ...
                voter_static_polarization( ...
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

save('polarization_ER.mat','polarization','alpha_values','initial_polarization');

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

xlabel('Initial polarization');
ylabel('Polarization probability');

legend(arrayfun(@(a) sprintf('\\alpha = %d',a), ...
    alpha_values,'UniformOutput',false), ...
    'Location','best');

xlim([0,1]);
ylim([0,1]);

%% ============================================================
%  Function 1: Conservative learning on ER network
%% ============================================================

function x_plus = voter_static_polarization(alpha,polarization_initial,N)

    %% ---------- Generate network ----------

    kbar = 8;

    L = kbar*N/2;

    adj = generate_pairwise_ER(N,L);

    degree = sum(adj,2);

    neighbor = cell(N,1);

    for i = 1:N
        neighbor{i} = find(adj(i,:));
    end

    %% ---------- Initial opinions ----------

    size_A = randi(polarization_initial);
    size_B = polarization_initial-size_A;

    vec = [ones(1,size_A), ...
           zeros(1,N-size_A-size_B), ...
           -ones(1,size_B)];

    a_adj = vec(randperm(N));

    %% ---------- Opinion counts ----------

    nPlus  = nnz(a_adj==1);
    nMinus = nnz(a_adj==-1);
    nZero  = nnz(a_adj==0);

    %% ---------- Evolution ----------

    while true

        % Step 1: Randomly select a focal individual
        i = randi(N);

        % Step 2: Randomly select one of its neighbors
        nei = neighbor{i};

        deg = degree(i);

        if deg == 0
            continue;
        end

        j = nei(randi(deg));

        % Opinions of focal individual and selected neighbor
        old_i = a_adj(i);
        opinion_j = a_adj(j);

        % Step 3: Calculate neighborhood opinion fractions
        states = a_adj(nei);

        fA = nnz(states==1)/deg;
        fB = nnz(states==-1)/deg;
        fC = nnz(states==0)/deg;

        % Step 4: Conservative learning rule
        p = rand;

        % A -> C
        if old_i == 1 && opinion_j == 0

            if p <= fC^(alpha-1)

                a_adj(i) = 0;

                nPlus = nPlus-1;
                nZero = nZero+1;

            end

        % B -> C
        elseif old_i == -1 && opinion_j == 0

            if p <= fC^(alpha-1)

                a_adj(i) = 0;

                nMinus = nMinus-1;
                nZero = nZero+1;

            end

        % C -> A
        elseif old_i == 0 && opinion_j == 1

            if p <= fA^(alpha-1)

                a_adj(i) = 1;

                nZero = nZero-1;
                nPlus = nPlus+1;

            end

        % C -> B
        elseif old_i == 0 && opinion_j == -1

            if p <= fB^(alpha-1)

                a_adj(i) = -1;

                nZero = nZero-1;
                nMinus = nMinus+1;

            end

        end

        % A and B cannot directly learn from each other

        %% ---------- Absorbing states ----------

        if nZero == N || nZero == 0

            x_plus = 1-nZero/N;

            break;

        end

    end

end


%% ============================================================
%  Function 2: Generate ER network
%% ============================================================

function A = generate_pairwise_ER(N,L)

    % Step 1: Generate a random spanning tree
    perm = randperm(N);

    edges = zeros(N-1,2);

    for k = 2:N

        i = perm(k);
        j = perm(randi(k-1));

        edges(k-1,:) = [i,j];

    end

    % Step 2: Construct adjacency matrix
    A = zeros(N,N);

    for k = 1:N-1

        i = edges(k,1);
        j = edges(k,2);

        A(i,j) = 1;
        A(j,i) = 1;

    end

    edgeCount = N-1;

    % Step 3: Randomly add edges
    while edgeCount < L

        i = randi(N);
        j = randi(N);

        if i ~= j && A(i,j)==0

            A(i,j) = 1;
            A(j,i) = 1;

            edgeCount = edgeCount+1;

        end

    end

end
