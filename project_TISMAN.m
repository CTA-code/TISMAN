%% Script for Transcriptomics-Informed Stoichiometric Modelling and Network Analysis (TISMAN)
% Claudio Tomi-Andrino (2021)
% other files can be found in https://github.com/CTA-code/TISMAN

%% Load the COBRA-compatible model
% model: Human-GEM 1.3.0, as per 2019-12-19 https://github.com/SysBioChalmers/Human-GEM
% the model was made using RAVEN (10.1126/scisignal.aaz1482), so a couple
% of steps were necessary to make it COBRA-compatible
% https://sysbiochalmers.github.io/Human-GEM-guide/getting_started/#loading-human-gem-into-matlab
tic
initCobraToolbox(false)
load('HumanGEM.mat');
backup_model = ravenCobraWrapper(ihuman);
changeCobraSolver('ibm_cplex')
save backup_model.mat
% load('backup_model.mat')
backup_model.lb(7733) = -2 % upper limit for oxygen uptake rate, as done elsewhere (10.1038/srep41241 (2017))

%% rFASTCORMICS
% obtain a flux consistent model using 'fastcc'
% this function generates a vector 'A' comprising the indeces of flux
% consistent reactions. However, the generated model 'modelFlipped' does
% not really implement this information. Therefore, it is necessary to
% generate another model from it
model = backup_model
epsilon = 1e-4
printLevel = 0
[A, modelFlipped, V] = fastcc(model, epsilon, printLevel)
model_mod = model
clear modelFlipped
clear V

% generate an updated reaction list and S
num_rxns = size(A,1);
c_lb = zeros(num_rxns,1);
c_ub = zeros(num_rxns,1);
c_rev = zeros(num_rxns,1);
c_c = zeros(num_rxns,1);
c_rxnConfidenceScores = zeros(num_rxns,1);
c_rxnNames = num2cell(zeros(num_rxns,1));
c_subSystems = num2cell(zeros(num_rxns,1));
c_eccodes = num2cell(zeros(num_rxns,1));
c_rxnReferences = num2cell(zeros(num_rxns,1));
c_rxnFrom = num2cell(zeros(num_rxns,1));
c_rxns = num2cell(zeros(num_rxns,1));
c_S = zeros(size(model.mets,1),num_rxns);
c_rxnGeneMat = zeros(num_rxns,size(model.genes,1));

% convert sparse to double
S_ori = full(model_mod.S);
rxnGeneMat_ori = full(model_mod.rxnGeneMat);

for i = 1:size(A,1)
    tmp_indx = A(i,1)
    c_lb(i,1) = model_mod.lb(tmp_indx,1);
    c_ub(i,1) = model_mod.ub(tmp_indx,1);
    c_rev(i,1) = model_mod.rev(tmp_indx,1);
    c_c(i,1) = model_mod.c(tmp_indx,1);
    c_rxnConfidenceScores(i,1) = model_mod.rxnConfidenceScores(tmp_indx,1);
    c_rxnNames(i,1) = model_mod.rxnNames(tmp_indx,1);
    c_subSystems(i,1) = model_mod.subSystems(tmp_indx,1);
    c_rxnReferences(i,1) = model_mod.rxnReferences(tmp_indx,1);
    c_rxns(i,1) = model_mod.rxns(tmp_indx,1);
    c_S(:,i) = S_ori(:,tmp_indx);
    c_rxnGeneMat(i,:) = rxnGeneMat_ori(tmp_indx,:);
end
model_mod.lb = c_lb;
model_mod.ub = c_ub;
model_mod.rev = c_rev;
model_mod.c = c_c;
model_mod.rxnConfidenceScores = c_rxnConfidenceScores;
model_mod.rxnNames = c_rxnNames;
model_mod.subSystems = c_subSystems;
model_mod.rxnReferences = c_rxnReferences;
model_mod.rxns = c_rxns;
model_mod.S = sparse(c_S);
model_mod.rxnGeneMat = sparse(c_rxnGeneMat);
model_mod.rules = [];
positions = model_mod.rxns;
model_mod = createEmptyFields(model_mod, 'rules')
model_mod = creategrRulesField(model_mod, positions)
% save c_human_gem.mat
% load('c_human_gem.mat')
Cmodel = model_mod

%%% import medoids data, as generated in the R workflow
%%% Now we load the medoids previously calculated (check workflow in R)
% Medoid 1
%pre_medoid = readtable('human_gem_t_medoid_1.txt');
pre_medoid = readtable('human_gem_t_medoid_2.txt');
pre_medoid = table2cell(pre_medoid);
medoid_1(:,1) = pre_medoid(:,1);
for i = 1:size(medoid_1,1)
    tmp = pre_medoid(i,2);
    tmp_num = str2num(cell2mat(tmp));
    if tmp_num == 1
        medoid_1(i,2) = num2cell(1); % force active
    elseif tmp_num == 0
        medoid_1(i,2) = num2cell(0); % do not force
    else
        medoid_1(i,2) = num2cell(0); % do not force
    end
end
clear tmp
clear tmp_num

% get a binary list of reactions, flagging those to be forced to be active
% for medoid 1
bound_medoid = Cmodel.lb;
bound_medoid(:,2) = Cmodel.ub;
for i = 1:size(medoid_1,1)
    if cell2mat(medoid_1(i,2)) == 0
        tmp_indx = find(Cmodel.rxnGeneMat(:,i));
        for j = 1:size(tmp_indx,1) 
            tmp_rxn = tmp_indx(j,1);
            bound_medoid(tmp_rxn,1) = 0;
            bound_medoid(tmp_rxn,2) = 0;
        end  
    clear tmp_rxn
    end
end
rxns_medoid = Cmodel.rxns;
rxns_medoid(:,2) = num2cell(zeros(size(bound_medoid,1),1));
for i = 1:size(bound_medoid,1)
    tmp_lb = bound_medoid(i,1);
    tmp_ub = bound_medoid(i,2);
    if tmp_lb < 0 || tmp_ub > 0
        rxns_medoid(i,2) = num2cell(1);
    end
    clear tmp_lb
    clear tmp_ub
end
active_medoid = num2cell(zeros(nnz(cell2mat(rxns_medoid(:,2))),1));
j = 1;
for i = 1:size(rxns_medoid,1)
    tmp_act = cell2mat(rxns_medoid(i,2));
    if tmp_act == 1
        tmp_rxn = rxns_medoid(i,1);
        active_medoid(j,1) = tmp_rxn;
        j = j+1;
    end
    clear tmp_act
    clear tmp_rxn
end

% rFASTCORMICS for each ((c) Dr. Maria Pires Pacheco 2016_)
addpath(genpath(pwd)) % add all subfolders to path
% load colnames
%headers = readtable('headers_fpkm_fast_1_corr.txt');
headers = readtable('headers_fpkm_fast_2_corr.txt');
headers_table = table2cell(headers);
colnames = transpose(headers_table);
clear headers
clear headers_table
% load FPKM
%fpkm_data = readtable('re_fpkm_fast_cluster_1.txt','ReadVariableNames',false,'ReadRowNames',true);
fpkm_data = readtable('re_fpkm_fast_cluster_2.txt','ReadVariableNames',false,'ReadRowNames',true);
pre_fpkm = table2cell(fpkm_data);
fpkm = cell2mat(pre_fpkm);
clear pre_fpkm
% load rownames
rownames = fpkm_data.Properties.RowNames;
clear fpkm_data
figflag = 1; % set figure flag: 1 to output and save density plot
% Data discretization
% log2-transform the data
signal = log2(fpkm); % log2-transform the data
signal(isinf(signal)) = -10000; % transform all -Inf to -10000
for j = 1:size(colnames,2) % for each sample
    signal_sample = signal(:,j); % save only current sample
    data_keep = fpkm(:,j);
    signal_sample = signal_sample(signal_sample>-10000); % remove samples with low expression
% Density plot
    [probability_estimate,xi] = ksdensity((signal_sample)); % get the densities
% Create right curve
    peak_max_idx = (find(probability_estimate==max(probability_estimate))); % find the maximum of the density curve
    max_r_side = probability_estimate(peak_max_idx:end); % copy right-most side of density curve
    max_l_side = flip(max_r_side); % mirror the right side
    hybrid = [max_l_side(1:end-1), max_r_side];
    hybrid_curve = zeros(numel(probability_estimate),1); % create hybrid curve
    if numel(hybrid)> numel(probability_estimate)
        hybrid_curve(end+1-numel(hybrid(end-100+1:end)):end)=hybrid(end-100+1:end);
        x = zeros(numel(probability_estimate),1); % create new x-axis
    else
        hybrid_curve(end-numel(hybrid)+1:end) = hybrid;
        x = xi;
    end
    clear peak_max_idx max_l_side max_r_side hybrid
% Create left curve (main curve - rightmost curve)
    rest_curve = probability_estimate - hybrid_curve';
    rest_curve(find(rest_curve<0.0001))=0;
% fit the curves
    [fitresult_r, ~] = createFit_r(xi, hybrid_curve); % right curve
    [fitresult_l, ~] = createFit_l(xi, rest_curve); % left curve  
%  zFRMA transform the data and plot the density
    sigma1 = fitresult_r.c1/sqrt(2);
    mu1 = fitresult_r.b1; % x-value of right max
    zFRMA = (signal_sample-mu1)/sigma1;
    [yFRMA,xFRMA] = ksdensity(zFRMA);
    clear hybrid_curve probability_estimate rest_curve yFRMA xFRMA x xi ans
% Discretization
    zFRMA = (signal(:,j)-mu1)/sigma1; % calculate zFRMA 
    sigma2 = fitresult_l.c1/sqrt(2);
    mu2 = fitresult_l.b1; % x-value of left max
    mu2_t=(mu2-mu1)/sigma1;
    discretized = zFRMA;
    e = 0;
    ue = max(-3, mu2_t);
    zFRMA = reshape(zFRMA,size(data_keep,1),size(data_keep,2));
    exp_threshold = e;
    unexp_threshold = ue;
    exp=find(discretized>=exp_threshold);
    not_exp=find(discretized<=unexp_threshold);
    unknown = find(discretized<exp_threshold & discretized>unexp_threshold);
    discretized(unknown)    = 0;    % unknown expression is set to 0
    discretized(not_exp)    = -1;   % non-expression is set to -1
    discretized(exp)        = 1;    % expression is set to 1
    discretized = (reshape(discretized,size(data_keep,1),size(data_keep,2)));
    discretized_keep(j,:) = discretized'; % save for each sample
    clear fitresult_l fitresult_r data_keep e exp exp_threshold mu1 mu2 mu2_t not_exp sigma1 sigma2 signal2 signal_sample ue unexp_threshold unknown zFRMA
end
discretized=discretized_keep';
clear j discretized_keep signal
load medium_example % need to define medium for the cells used here
% overwrite dictionary
% recon.genes is based on Entrez IDs, whereas the example FPKM values use
% gene names; the FPKM values for GBM are identified by Ensemble IDs, so it
% is necessary to create a new dictionary
load dico_ML.mat % dictionary to map the rownname identifier to the genes in the model
dico_mod = [model_mod.genes model_mod.genes] 
dico = cell2table(dico_mod);
Cmodel.genes = regexprep(Cmodel.genes,'\.[0-9]+$',''); %removal of gene transcripts
consensus_proportion = 0.9; %gene has to be expressed in 90% of the cases in order to be included. Only relevant if you want to create one generic model from different samples
% Set optional settings such as:
unpenalizedSystems = {'Transport reactions'};
subsystems = Cmodel.subSystems
subsystems_2 = num2cell(zeros(size(subsystems,1),1))
for i = 1:size(subsystems,1)
    tmp_subs = subsystems{i,1};
    subsystems_2(i,1) = tmp_subs;
end
clear tmp_subs
Cmodel.subSystems = subsystems_2;
unpenalized = Cmodel.rxns(ismember(Cmodel.subSystems,unpenalizedSystems));
optional_settings.unpenalized = unpenalized
optional_settings.func = ['HMR_4358'; 'HMR_6552'; active_medoid] % forced reactions

% Create models
% single models:
for i = 1:numel(colnames) %for each sample
    [model_out{i}, A_keep{i}] = fastcormics_2018(Cmodel, discretized(:,i), ...
        rownames, dico , 0, consensus_proportion, epsilon, optional_settings);
end

% generic models:
[model_out_generic, A_keep_generic] = fastcormics_2018(Cmodel, discretized, ...
    rownames, dico , 0, consensus_proportion, epsilon, optional_settings);

backup_mod_model = model_out_generic % overwrite 
%save model_out_generic_forced_1.mat
save model_out_generic_forced_2.mat
rc = backup_mod_model.rxns
toc

% checking GPRs
% All genes are linked to at least 1 reaction, but only about half of the
% reactions are linked to at least 1 gene
gpr = backup_mod_model.rxnGeneMat;
% from the rxn point of view
gpr_list_rxn = num2cell(zeros(size(gpr,1),2));
gpr_list_rxn(:,1) = backup_mod_model.rxns;
for i = 1:size(gpr,1)
    tmp = num2cell(nonzeros(gpr(i,:)));
    tmp_size = size(tmp,1);
    
    if tmp_size > 0
        gpr_list_rxn(i,2) = num2cell(1);
        tmp_indx = find(gpr(i,:));
        gpr_list_rxn(i,3) = {num2cell(tmp_indx)};
        tmp_genes = [];
        for j = 1:size(tmp_indx,2)
            tmp_genes = [tmp_genes; backup_mod_model.genes(tmp_indx(1,j))];
        end
        gpr_list_rxn(i,4) = {tmp_genes};
    else
        gpr_list_rxn(i,2) = num2cell(0);
    end 
    clear tmp
    clear tmp_size
    clear tmp_indx
end
clear i
clear tmp_genes
rxn_with_gpr = nnz(cell2mat(gpr_list_rxn(:,2)));
rxn_with_gpr_percent = (rxn_with_gpr/size(backup_mod_model.rxns,1))*100;

% from the gene point of view
gpr_list_gene = num2cell(zeros(size(gpr,2),2));
gpr_list_gene(:,1) = backup_mod_model.genes;
for i = 1:size(gpr,2)
    tmp = num2cell(nonzeros(gpr(:,i)));
    tmp_size = size(tmp,1);
    
    if tmp_size > 0
        gpr_list_gene(i,2) = num2cell(1);
    else
        gpr_list_gene(i,2) = num2cell(0);
    end 
    clear tmp
    clear tmp_size
end
clear i
gene_with_gpr = nnz(cell2mat(gpr_list_gene(:,2)));
gene_with_gpr_percent = (gene_with_gpr/size(backup_mod_model.genes,1))*100;

% need to find the new indeces for rxns of interest when defining FBA
% ids of selected objective functions (rxns). ATP formation will be
% modelled indirectly. Production of pyruvate in the glycolysis will be
% maximised, and later on fluxes out of it considered (oxidative
% phosphorilation, anaerobic glycolysis or aerobic glycolysis)
list_rxns = model_out_generic.rxns;
biomass = 'biomass_human' % biomass_human
ATP_synth = 'HMR_4358' % HMR_4358; pyruvate kinase (https://metabolicatlas.org/explore/gem-browser/human1/reaction/HMR_4358)
migration = 'HMR_6552' % HMR_6552; synaptojanin 2, producing a cell membrane phospholipid regulating cell migration and invasion in anaplastic intracranial ependymomas (10.1016/j.cancergencyto.2007.12.008)
% obtain rxn index
for i = 1:size(list_rxns,1)
    tmp_rxn = list_rxns(i);
    if strcmp(tmp_rxn,biomass) == 1
        id_biomass = i
    elseif strcmp(tmp_rxn,ATP_synth) == 1
        id_ATP_synthesis = i 
    elseif strcmp(tmp_rxn,migration) == 1
        id_migration = i 
    end
end
clear tmp_rxn

%%% Simplex lattice design
% Generation of different objective functions, and selection of unique
% solutions
q = 3;  % number of components (i.e. objective functions)
n = 4;  % number of elements -> that's 5 values when also considering 0
N = factorial(n+q-1)/((factorial(q-1))*factorial(n));
matrix_mixtures = zeros(N,q);
% first coordinate + number rows starting with it = n + 1
k = n + 1;
% matlab starts counting positions on number 1 instead of 0, so some
% tweaking is necessary
k_p = k + 1;
j = 1;
p = 1;
for i = 1:k
    num_rows = k_p - i;
    while j <= num_rows
        matrix_mixtures(nnz((matrix_mixtures(:,1)))+1,1) = i;
        matrix_mixtures(nnz((matrix_mixtures(:,2)))+1,2) = p;
        j = j + 1;
        p = p + 1;
    end
    j = 1;
    p = 1;
end
% now we just need to shift all the numbers by 1, and calculate the last
% column
matrix_mixtures_shift(:,1) = matrix_mixtures(:,1) - 1;
matrix_mixtures_shift(:,2) = matrix_mixtures(:,2) - 1;
matrix_mixtures_shift(:,3) = n - matrix_mixtures_shift(:,1) - matrix_mixtures_shift(:,2);
% finally, we normalise the rows so that cells add up 1. Each column can
% only take the values 0, 0.25, 0.50, 0.75 and 1.
factor = 10/n;
matrix_mixtures_fin = matrix_mixtures_shift*factor/10;
clear q
clear n
clear i
clear j
clear p
clear k
clear k_p
% set a lb for the biomass reaction for all tests
FBA_pre = optimizeCbModel(backup_mod_model); 
min_obj = round(0.5*FBA_pre.f);
backup_mod_model.lb(id_biomass) = min_obj;

%%% to store the results
number_tests = size(matrix_mixtures_fin,1) + 1; % need to add min sum flux!
indx_vector_c = [1:number_tests];
testArray = sprintfc('%02d', indx_vector_c)';
for vector_c = 1:numel(testArray)
Varnames_model_specific{vector_c} = matlab.lang.makeValidName(strcat('model_specific_',testArray{vector_c}));
Varnames_FBA{vector_c} = matlab.lang.makeValidName(strcat('FBA_',testArray{vector_c}));
Varnames_essential_genes{vector_c} = matlab.lang.makeValidName(strcat('list_essential_genes',testArray{vector_c}));
Varnames_essential_reactions{vector_c} = matlab.lang.makeValidName(strcat('list_essential_rxns',testArray{vector_c}));
Varnames_pagerank{vector_c} = matlab.lang.makeValidName(strcat('pagerank_rxns_',testArray{vector_c}));
Varnames_nodenames{vector_c} = matlab.lang.makeValidName(strcat('node_names_',testArray{vector_c}));
Varnames_common_chokepoints_lethal_consumption_lethal_production{vector_c} = matlab.lang.makeValidName(strcat('common_chokepoints_lethal_consumption_lethal_production_',testArray{vector_c}));
Varnames_list_potential_targets{vector_c} = matlab.lang.makeValidName(strcat('list_potential_targets_',testArray{vector_c}));
end       
failed_track = [0]
%
for vector_c = 1:number_tests
    try 
model_specific = backup_mod_model; % load the pre-constrained model
model_specific.c = zeros(size(model_specific.c,1),1); % zeroing the weights for objective function(s)

if vector_c <= number_tests - 1
    model_specific.c(id_biomass) = matrix_mixtures_fin(vector_c,1);
    model_specific.c(id_ATP_synthesis) = matrix_mixtures_fin(vector_c,2);
    model_specific.c(id_migration) = matrix_mixtures_fin(vector_c,3);
else
    model_specific.c = ones(size(model_specific.c,1),1); % min sum fluxes
end
stored_model_specific.(Varnames_model_specific{vector_c}) = model_specific;

model = model_specific;
FBA = optimizeCbModel(model);
stored_FBA.(Varnames_FBA{vector_c}) = FBA;

%%% Essential reactions
% since GPR are not complete, we should have a look at reaction
% essentiality by itself. Focusing only on reactions with a known GPR would
% mean that other potential chokepoints (later on the code) would not be
% identified once the genome annotation improves. This way, we have a
% matrix that is later filtered, and focus here on the known ones
num_rxns = size(model_specific.rxns,1);
ko_w_feasible_solution_FBA = num2cell(zeros(num_rxns,3));
ko_w_feasible_solution_FBA(:,1) = num2cell(1:num_rxns)';
ko_w_feasible_solution_FBA(:,2) = model_specific.rxns;

for i = 1:num_rxns
    mymodel = model_specific;
    try
        mymodel.lb(i) = 0;
        mymodel.ub(i) = 0;
        solFBA = optimizeCbModel(mymodel);  % solFBA refers to these solutions, whereas FBA is the overall one
        ko_w_feasible_solution_FBA(i,3) = num2cell(solFBA.x(id_biomass));     
    catch
        ko_w_feasible_solution_FBA(i,3) = num2cell(-1);       % infeasible conditions will be flagged with -1
    end
    
    clear mymodel
end
clear i

% select a cut-off: values lower than 10e-12 are considered 0 (no growth,
% i.e. lethal). INSTEAD: values lower than 50% of FBA.f. 
% Surprisingly, sometimes shutting down certain reactions does improve a
% little bit the growth -> check pharma review
for i = 1:num_rxns
    if cell2mat(ko_w_feasible_solution_FBA(i,3)) <= 0.95*FBA.x(id_biomass) % instead of FBA.f, since more than one objective function is used, and the interesting part is ensuring growth
        ko_w_feasible_solution_FBA(i,4) = num2cell(1);      % flagging with '1' 
    else
        ko_w_feasible_solution_FBA(i,4) = num2cell(0); 
    end    
end
clear i

list_essential_rxns = [];
for i =1:size(ko_w_feasible_solution_FBA,1)
    if cell2mat(ko_w_feasible_solution_FBA(i,4)) == 1
        list_essential_rxns = [list_essential_rxns; ko_w_feasible_solution_FBA(i,:)];
    end    
end
clear i

% Add GPR data
for i = 1:size(list_essential_rxns,1)
    tmp_ind = cell2mat(list_essential_rxns(i,1));
    list_essential_rxns(i,5) = gpr_list_rxn(tmp_ind,4);   
end
clear i
stored_essential_reactions.(Varnames_essential_reactions{vector_c}) = list_essential_rxns;

% First, it is necessary to obtain a 'clean' binary S
% list of active rxns
model = model_specific;
rxns_num = size(model.rxns,1);

rxns_list = num2cell(zeros(rxns_num, 4));
rxns_list(:,1) = model.rxns;
rxns_list(:,2) = num2cell(FBA.x);
% we consider numbers e-14 to be numerical artifacts, so change it to 0
tmp = cell2mat(rxns_list(:,2));

for i = 1:rxns_num
    if tmp(i) > 0
        if tmp(i) < 1e-13
            tmp(i) = 0;
        end
    elseif tmp(i) < 0
        if tmp(i) > -1e-13
            tmp(i) = 0;
        end
    end
end
clear i
rxns_list(:,3) = num2cell(tmp);
clear tmp

% In an edge list for a directed graph, the edge goes from the first node
% to the second. Since some fluxes are negative, a correction factor will be
% necessary later when generating the binary S (saved independently)
tmp = cell2mat(rxns_list(:,3));
for i = 1:rxns_num
    if tmp(i) > 0
        tmp(i) = 1;
    elseif tmp(i) < 0
        tmp(i) = -1;
    else
        tmp(i) = 0;
    end
end
clear i
rxns_list(:,4) = num2cell(tmp);
rxns_list(:,5) = printRxnFormula(model_specific);
clear tmp

% now we start modifying the stoichiometric matrix
S = model.S;
S_corr_direc = S;

for i = 1:rxns_num
    if cell2mat(rxns_list(i,4)) == -1
        S_corr_direc(:,i) = -1*S_corr_direc(:,i);
    end
end
clear i

% we need to remove rxns that are not active
null_rxns = ones(rxns_num,1);

for i = 1:rxns_num
    if cell2mat(rxns_list(i,4)) == 0
        null_rxns(i,:) = 0;
    end
end
clear i

list_rxns(:,1) = model.rxns;
list_rxns(:,2) = num2cell(null_rxns);
list_rxns(:,3) = num2cell(1:size(list_rxns,1)); % keep track of the original numbering

% export active rxns
S_corr_direc_non_null_rxns = zeros(size(S_corr_direc,1),nnz(null_rxns));

j = 1;
for i = 1:rxns_num
    if null_rxns(i) == 1 %therefore active, either FW or RV
        S_corr_direc_non_null_rxns(:,j) = S_corr_direc(:,i);
        list_rxns_active(j,1) = list_rxns(i,1);
        list_rxns_active(j,2) = list_rxns(i,2); 
        list_rxns_active(j,3) = list_rxns(i,3); 
        j = j + 1;
    end
end
clear i
clear j

%%% identify metabolites that do not participate in any of the active
% reactions. If a met is inactive, then the sum of the row must be zero
num_met_S_corr_direc_non_null_rxns = size(S_corr_direc_non_null_rxns,1);
num_rxn_S_corr_direc_non_null_rxns = size(S_corr_direc_non_null_rxns,2);
active_met_S_corr_direc_non_null_rxns = zeros(num_met_S_corr_direc_non_null_rxns,1);

for i = 1:num_met_S_corr_direc_non_null_rxns
    active_met_S_corr_direc_non_null_rxns(i,1) = nnz(S_corr_direc_non_null_rxns(i,:));
end
clear i

% we need to remove mets that are not active
null_mets = ones(num_met_S_corr_direc_non_null_rxns,1);

for i = 1:num_met_S_corr_direc_non_null_rxns
    if active_met_S_corr_direc_non_null_rxns(i,1) == 0 % therefore inactive
        null_mets(i,:) = 0;
    end
end
clear i

list_mets(:,1) = model.mets;
list_mets(:,2) = num2cell(null_mets);

% export S with both active rxns and mets
S_clean = zeros(nnz(null_mets),nnz(null_rxns));

j = 1;
for i = 1:num_met_S_corr_direc_non_null_rxns
    if null_mets(i) == 1 %therefore active
        S_clean(j,:) = S_corr_direc_non_null_rxns(i,:);
        list_mets_active(j,1) = list_mets(i,1);
        list_mets_active(j,2) = list_mets(i,2);  
        j = j + 1;
    end
 end
clear i
clear j

list_mets_final = list_mets_active;

test_active_rxns = ones(1,size(S_clean,2));
for i = 1:num_rxn_S_corr_direc_non_null_rxns
    test_active_rxns(1,i) = nnz(S_clean(:,i));
end
nnz(test_active_rxns)
clear i

test_active_mets = ones(size(S_clean,1),1);
for i = 1:size(S_clean,1)
    test_active_mets(i,1) = nnz(S_clean(i,:));
end
nnz(test_active_mets)
clear i

%%% more cleaning
% for the sake of simplicity, we could remove transport reactions (but we
% are not; all fields will be '1', and there nothing will be changed for
% this particular example)
for i = 1:size(list_rxns_active, 1)
   list_rxns_active(i,4) = num2cell(1);
   list_rxns_active(i,5) = num2cell(1);
   list_rxns_active(i,6) = num2cell(1);
   list_rxns_active(i,7) = num2cell(1);
   list_rxns_active(i,8) = num2cell(1);
   list_rxns_active(i,9) = num2cell(1);
   list_rxns_active(i,10) = num2cell(1);
end
clear i

j = 1;
for i = 1:size(S_clean,2)
    if cell2mat(list_rxns_active(i,4)) == 1 && cell2mat(list_rxns_active(i,5)) == 1 && cell2mat(list_rxns_active(i,6)) == 1 && cell2mat(list_rxns_active(i,7)) == 1 && cell2mat(list_rxns_active(i,8)) == 1 && cell2mat(list_rxns_active(i,9)) == 1 && cell2mat(list_rxns_active(i,10)) == 1 
        S_clean_1b(:,j) = S_clean(:,i);
        list_rxns_active_1b(j,1) = list_rxns_active(i,1);
        list_rxns_active_1b(j,2) = list_rxns_active(i,2); 
        list_rxns_active_1b(j,3) = list_rxns_active(i,3);
        j = j + 1;
    end
end
clear i
clear j

S_clean_final = S_clean_1b;

try
   list_rxns_final = list_rxns_active_1b;   
end

% not removing side compounds in this script
list_mets_active(:,3) = num2cell(ones((size(list_mets_active, 1)),1));

try
j = 1;
for i = 1:size(S_clean_1b,1)
    if cell2mat(list_mets_active(i,3)) == 1 %therefore of interest
        S_clean_1c(j,:) = S_clean_1b(i,:);
        list_mets_active_1b(j,1) = list_mets_active(i,1);
        list_mets_active_1b(j,2) = list_mets_active(i,2); 
        j = j + 1;
    end
end
clear i
clear j
list_mets_final = list_mets_active_1b;

% now check for active rxns
% another cycle checking for active reactions is then necessary
num_rxn_S_clean_1c = size(S_clean_1c,2);
active_rxn_S_clean_1c = zeros(1,num_rxn_S_clean_1c);
for i = 1:num_rxn_S_clean_1c
    active_rxn_S_clean_1c(1,i) = nnz(S_clean_1c(:,i));
end
clear i

null_rxns_1c = ones(num_rxn_S_clean_1c,1);
for i = 1:num_rxn_S_clean_1c
    if active_rxn_S_clean_1c(1,i) == 0
        null_rxns_1c(i,:) = 0; % inactive
    end
end
clear i

if nnz(~null_rxns_1c) == 0 
    disp('no more cleaning steps are necessary')
    S_clean_final = S_clean_1c;
else
    disp('more cleaning steps ARE necessary')
    %return
end

%%% cleaning again
S_clean_2 = zeros(size(S_clean_1c,1),nnz(null_rxns_1c));

j = 1;
for i = 1:size(null_rxns_1c,1)
    if null_rxns_1c(i) == 1 %therefore active
        S_clean_2(:,j) = S_clean_1c(:,i);
        list_rxns_active_2(j,1) = list_rxns_active_1b(i,1);
        list_rxns_active_2(j,2) = list_rxns_active_1b(i,2); 
        list_rxns_active_2(j,3) = list_rxns_active_1b(i,3);
        j = j + 1;
    end
end
clear i
clear j

% if list_rxns_active_2 was generated, then overwrite the previous list
try
   list_rxns_final = list_rxns_active_2;   
end

test_active_rxns_2 = ones(1,size(S_clean_2,2));
for i = 1:size(S_clean_2,2)
    test_active_rxns_2(1,i) = nnz(S_clean_2(:,i));
end
nnz(test_active_rxns_2)
clear i

test_active_mets_2 = ones(size(S_clean_2,1),1);
for i = 1:size(S_clean_2,1)
    test_active_mets_2(i,1) = nnz(S_clean_2(i,:));
end
nnz(test_active_mets_2)
clear i

if nnz(test_active_rxns_2) ~= size(test_active_rxns_2,2) || nnz(test_active_mets_2) ~= size(test_active_mets_2,1)
    disp('more cleaning steps ARE necessary')
    S_clean_2b = zeros((nnz(test_active_mets_2)),size(S_clean_2,2));
    j = 1;
    for i = 1:size(test_active_mets_2,1)
        if test_active_mets_2(i) ~= 0 %therefore active
            S_clean_2b(j,:) = S_clean_2(i,:);
            list_mets_active_2(j,1) = list_mets_active_1b(i,1);
            list_mets_active_2(j,2) = list_mets_active_1b(i,2);   
            j = j + 1;
        end
    end
    clear i
    
    test_active_rxns_2b = ones(1,size(S_clean_2b,2));
    for i = 1:size(S_clean_2b,2)
        test_active_rxns_2b(1,i) = nnz(S_clean_2b(:,i));
    end
    nnz(test_active_rxns_2b)
    clear i

    test_active_mets_2b = ones(size(S_clean_2b,1),1);
    for i = 1:size(S_clean_2b,1)
        test_active_mets_2b(i,1) = nnz(S_clean_2b(i,:));
    end
    nnz(test_active_mets_2b)
    clear i  
else
    S_clean_final = S_clean_2;
end

try
if nnz(test_active_rxns_2b) ~= size(test_active_rxns_2b,2) || nnz(test_active_mets_2b) ~= size(test_active_mets_2b,1)
    disp('more cleaning steps ARE necessary')
    %return
else
    S_clean_final = S_clean_2b;
end
catch
end

% if list_mets_active_2 was generated, then overwrite the previous list
try
   list_mets_final = list_mets_active_2;   
end
end

% Chokepoints
% reactions in charge of exclusively produce or consume a certain
% metabolite. However, said reaction may catalyse other reactions too
% call final lists of rxns and mets!!
S = S_clean_final;
S_list_met = num2cell(zeros(size(S,1),8));
S_list_met(:,1) = num2cell(1:size(S,1));
S_list_met(:,2) = list_mets_final(:,1);

for i = 1:size(S,1)
    tmp = num2cell(nonzeros(S(i,:)));
    tmp_size = size(tmp,1);
    
    if tmp_size > 0
        if tmp_size == 1
            S_list_met(i,3) = num2cell(1);
        elseif tmp_size == 2
            S_list_met(i,3) = num2cell(2);
        else
            S_list_met(i,3) = num2cell(-1);     % meaning, any number > 2
        end      
    else
        S_list_met(i,3) = num2cell(0);          % not used
    end 
    clear tmp
    clear tmp_size
end
clear i

% now, for mets that are linked to 1 or 2 reactions, need to recall them
% (both index and name)
for i = 1:size(S_list_met,1)
    if cell2mat(S_list_met(i,3)) == 1
        tmp_indx_rxn = find(S(i,:));
        tmp_coeff = full(S(i,tmp_indx_rxn));
        if tmp_coeff < 0    % i.e. consumption
            S_list_met(i,4) = num2cell(tmp_indx_rxn);   % index in S_clean_final
            S_list_met(i,5) = list_rxns_final(tmp_indx_rxn,1);  % rxn name
            S_list_met(i,6) = list_rxns_final(tmp_indx_rxn,3);  % index in model_specific.S
        else                % i.e. production
            S_list_met(i,7) = num2cell(tmp_indx_rxn);   % index in S_clean_final
            S_list_met(i,8) = list_rxns_final(tmp_indx_rxn,1);  % rxn name
            S_list_met(i,9) = list_rxns_final(tmp_indx_rxn,3);  % index in model_specific.S
        end
          
    elseif cell2mat(S_list_met(i,3)) == 2
        tmp_indx_rxn = find(S(i,:));
        tmp_coeff = full(S(i,tmp_indx_rxn));
        if tmp_coeff(1,1) < 0 && tmp_coeff(1,2) > 0     % it cannot be 0!
            S_list_met(i,4) = num2cell(tmp_indx_rxn(1,1));
            S_list_met(i,5) = list_rxns_final(tmp_indx_rxn(1,1),1);
            S_list_met(i,6) = list_rxns_final(tmp_indx_rxn(1,1),3);
            S_list_met(i,7) = num2cell(tmp_indx_rxn(1,2));
            S_list_met(i,8) = list_rxns_final(tmp_indx_rxn(1,2),1);
            S_list_met(i,9) = list_rxns_final(tmp_indx_rxn(1,2),3);
        elseif tmp_coeff(1,1) > 0 && tmp_coeff(1,2) < 0
            S_list_met(i,4) = num2cell(tmp_indx_rxn(1,2));
            S_list_met(i,5) = list_rxns_final(tmp_indx_rxn(1,2),1);
            S_list_met(i,6) = list_rxns_final(tmp_indx_rxn(1,2),3);
            S_list_met(i,7) = num2cell(tmp_indx_rxn(1,1));
            S_list_met(i,8) = list_rxns_final(tmp_indx_rxn(1,1),1);
            S_list_met(i,9) = list_rxns_final(tmp_indx_rxn(1,1),3);
        end  
    end
end
clear i

list_chokepoints = [];
for i =1:size(S_list_met,1)
    if cell2mat(S_list_met(i,3)) == 2
        list_chokepoints = [list_chokepoints; S_list_met(i,:)];
    end    
end
clear i

%%% Intersection chokepoints and lethal reactions (i.e. essential)
final_lethal = cell2mat(list_essential_rxns(:,1));
final_chokepoints_consumption = cell2mat(list_chokepoints(:,6));
final_chokepoints_production = cell2mat(list_chokepoints(:,9));

common_chokepoints_lethal_consumption = intersect(final_lethal,final_chokepoints_consumption);
common_chokepoints_lethal_production = intersect(final_lethal,final_chokepoints_production);
common_chokepoints_lethal_consumption_lethal_production = num2cell(intersect(common_chokepoints_lethal_consumption,common_chokepoints_lethal_production));  % intersection rxns participating in chokepoints (at the same time consumption and production) with essential rxns

for i = 1:size(common_chokepoints_lethal_consumption_lethal_production,1)
    for j = 1:size(list_essential_rxns,1)
        if cell2mat(list_essential_rxns(j,1)) == cell2mat(common_chokepoints_lethal_consumption_lethal_production(i,1))
            common_chokepoints_lethal_consumption_lethal_production(i,2) = list_essential_rxns(j,2);    % rxn name
        end        
    end 
end
clear i

% the index in rxns_list corresponds to the original order, so we can
% easily obtain the reaction definition
for i = 1:size(common_chokepoints_lethal_consumption_lethal_production,1)
    tmp_indx = cell2mat(common_chokepoints_lethal_consumption_lethal_production(i,1));
    common_chokepoints_lethal_consumption_lethal_production(i,3) = rxns_list(tmp_indx,5); % rxn definition
end
clear tmp_indx
clear i

% get localisation of enzyme; there are 10 compartments in the model
local_enzymes = common_chokepoints_lethal_consumption_lethal_production;
local_enzymes(:,4) = num2cell(zeros); % s, extracellular
local_enzymes(:,5) = num2cell(zeros); % p, peroxisome
local_enzymes(:,6) = num2cell(zeros); % m, mitochondria
local_enzymes(:,7) = num2cell(zeros); % c, cytosol
local_enzymes(:,8) = num2cell(zeros); % l, lysosome
local_enzymes(:,9) = num2cell(zeros); % r, ER
local_enzymes(:,10) = num2cell(zeros); % g, Golgi
local_enzymes(:,11) = num2cell(zeros); % n, nucleus
local_enzymes(:,12) = num2cell(zeros); % x, boundary met?
local_enzymes(:,13) = num2cell(zeros); % i, inner mitochondria 

for i = 1:size(local_enzymes,1)
    tmp_string = char(local_enzymes(i,3));
    if contains(tmp_string,'[s]') == 1
        local_enzymes(i,4) = num2cell(1);
    end
    if contains(tmp_string,'[p]') == 1
        local_enzymes(i,5) = num2cell(1);
    end
    if contains(tmp_string,'[m]') == 1
        local_enzymes(i,6) = num2cell(1);
    end
    if contains(tmp_string,'[c]') == 1
        local_enzymes(i,7) = num2cell(1);
    end
    if contains(tmp_string,'[l]') == 1
        local_enzymes(i,8) = num2cell(1);
    end 
    if contains(tmp_string,'[r]') == 1
        local_enzymes(i,9) = num2cell(1);
    end
    if contains(tmp_string,'[g]') == 1
        local_enzymes(i,10) = num2cell(1);
    end
    if contains(tmp_string,'[n]') == 1
        local_enzymes(i,11) = num2cell(1);
    end
    if contains(tmp_string,'[x]') == 1
        local_enzymes(i,12) = num2cell(1);
    end
    if contains(tmp_string,'[i]') == 1
        local_enzymes(i,13) = num2cell(1);
    end
end
clear i

% add met consumption (by the reaction)
for i = 1:size(common_chokepoints_lethal_consumption_lethal_production,1)
    tmp_rxn_num = cell2mat(common_chokepoints_lethal_consumption_lethal_production(i,1));
    for j = 1:size(list_chokepoints,1)
        if cell2mat(common_chokepoints_lethal_consumption_lethal_production(i,1)) == cell2mat(list_chokepoints(j,6))
            common_chokepoints_lethal_consumption_lethal_production(i,4) = list_chokepoints(j,2);   % consumed met (considering directionalities)
        end        
    end
    clear tmp_rxn_num
end
clear i
clear j

% count how many compartments are related to each reaction. If > 1, then
% they are transport reactions
for i = 1:size(local_enzymes,1)
    local_enzymes(i,14) = num2cell(nnz(cell2mat(local_enzymes(i,4:13))));
end
common_chokepoints_lethal_consumption_lethal_production(:,5) = local_enzymes(:,14); % number of compartments used in the rxn
clear i

% In order to get a list of enzymes to further
% study, we need to check whether the reactions of interest are actually
% linked to known genes (and therefore, proteins). Gpr list follows the
% original order of rxns, so the index in
% common_chokepoints_lethal_consumption_lethal_production can be used
for i = 1:size(common_chokepoints_lethal_consumption_lethal_production,1)
    tmp_indx = cell2mat(common_chokepoints_lethal_consumption_lethal_production(i,1));
    common_chokepoints_lethal_consumption_lethal_production(i,6) = (gpr_list_rxn(tmp_indx,4));  % gene(s) associated to the rxn (via GPR)
end

stored_common_chokepoints_lethal_consumption_lethal_production.(Varnames_common_chokepoints_lethal_consumption_lethal_production{vector_c}) = common_chokepoints_lethal_consumption_lethal_production;

%%% important reactions via the adjoint matrix

% in Palsson's framework, biomass is not modelled as a sink (as opposed to
% Fell's), so biomass precursors are generated and 'going nowhere'. Thus,
% the biomass reaction will be modified and a pseudometabolite 'biomass'
% created. We need to add an extra row and a sij for it:
S_clean_final_biomass = [S_clean_final; zeros(1,size(S_clean_final,2))];
pre_id_biomass_rxn = strcmp({'biomass_human'}, list_rxns_final(:,1));
id_biomass_rxn = find(pre_id_biomass_rxn > 0);
S_clean_final_biomass(size(S_clean_final_biomass,1),id_biomass_rxn) = 1;
S_adjoint = S_clean_final_biomass';
S_adjoint = conj(S_adjoint).*(-1); % complex conjugate (and correction factor for the directions)
adjacency_matrix = zeros(size(S_adjoint,1),size(S_adjoint,1));
for i = 1:size(S_adjoint,2)
    tmp_column = S_adjoint(:,i);
    tmp_subs = find(tmp_column < 0);
    tmp_prods = find(tmp_column > 0);  
    for si = 1:size(tmp_subs,1)
        for pi = 1:size(tmp_prods,1)
            adjacency_matrix(tmp_subs(si,1),tmp_prods(pi,1)) = adjacency_matrix(tmp_subs(si,1),tmp_prods(pi,1)) + 1;
        end
    end
    clear tmp_column
    clear tmp_subs
    clear tmp_prods
    clear tmp_num_comb    
end
clear i
clear si

node_names = list_rxns_final(:,1);
G = digraph(adjacency_matrix,node_names,'omitselfloops');
pagerank_weight = centrality(G,'pagerank','Importance',G.Edges.Weight);
stored_node_pagerank.(Varnames_pagerank{vector_c}) = pagerank_weight
stored_node_names.(Varnames_nodenames{vector_c}) = node_names
clear pagerank_weight node_names adjacency_matrix G

clear tmp_indx
clear i
clear final_chokepoints_consumption
clear final_chokepoints_production
clear final_lethal
clear ko_w_feasible_solution_FBA
clear list_chokepoints
clear list_essential_genes
clear list_essential_rxns
clear list_mets
clear list_mets_active
clear list_mets_active_1b
clear list_mets_final
clear list_potential_targets
clear list_rxns
clear list_rxns_active
clear list_rxns_active_1b
clear list_rxns_active_2
clear list_rxns_final
clear model
clear model_specific
clear null_mets
clear null_rxns
clear null_rxns_1c
clear num_essential_genes
clear num_met_S_corr_direc_non_null_rxns
clear num_rows
clear num_rxn_S_clean_1c
clear num_rxn_S_corr_direc_non_null_rxns
clear num_rxns
clear S
clear S_clean
clear S_clean_1b
clear S_clean_1c
clear S_2
clear S_clean_final
clear S_corr_direc
clear S_corr_direc_non_null_rxns
clear S_list_met
clear tmp_coeff
clear tmp_ind
clear tmp_indx_rxn
clear tmp_string
clear tmp_target
clear common_chokepoints_lethal_consumption_lethal_production

    catch
clear final_chokepoints_consumption
clear final_chokepoints_production
clear final_lethal
clear ko_w_feasible_solution_FBA
clear list_chokepoints
clear list_essential_genes
clear list_essential_rxns
clear list_mets
clear list_mets_active
clear list_mets_active_1b
clear list_mets_final
clear list_potential_targets
clear list_rxns
clear list_rxns_active
clear list_rxns_active_1b
clear list_rxns_active_2
clear list_rxns_final
clear model
clear model_specific
clear null_mets
clear null_rxns
clear null_rxns_1c
clear num_essential_genes
clear num_met_S_corr_direc_non_null_rxns
clear num_rows
clear num_rxn_S_clean_1c
clear num_rxn_S_corr_direc_non_null_rxns
clear num_rxns
clear S
clear S_clean
clear S_clean_1b
clear S_clean_1c
clear S_2
clear S_clean_final
clear S_corr_direc
clear S_corr_direc_non_null_rxns
clear S_list_met
clear tmp_coeff
clear tmp_ind
clear tmp_indx_rxn
clear tmp_string
clear tmp_target
clear common_chokepoints_lethal_consumption_lethal_production

failed_track = [failed_track; vector_c]
    end
disp(vector_c)
end
toc

clear final_chokepoints_consumption
clear final_chokepoints_production
clear final_lethal
clear ko_w_feasible_solution_FBA
clear list_chokepoints
clear list_essential_genes
clear list_essential_rxns
clear list_mets
clear list_mets_active
clear list_mets_active_1b
clear list_mets_final
clear list_potential_targets
clear list_rxns
clear list_rxns_active
clear list_rxns_active_1b
clear list_rxns_active_2
clear list_rxns_final
clear model
clear model_specific
clear null_mets
clear null_rxns
clear null_rxns_1c
clear num_essential_genes
clear num_met_S_corr_direc_non_null_rxns
clear num_rows
clear num_rxn_S_clean_1c
clear num_rxn_S_corr_direc_non_null_rxns
clear num_rxns
clear S
clear S_clean
clear S_clean_1b
clear S_clean_1c
clear S_2
clear S_clean_final
clear S_corr_direc
clear S_corr_direc_non_null_rxns
clear S_list_met
clear tmp_coeff
clear tmp_ind
clear tmp_indx_rxn
clear tmp_string
clear tmp_target
clear common_chokepoints_lethal_consumption_lethal_production

%%% Select non-redundant solutions
% Two solutions are assumed to be the same when the flux values are the
% same. In order to reduce the effect of numerical artifacts, let's round the number to the 5th decimal
try
% need to extract them from the new variables
for i = 1:number_tests
    tmp = stored_FBA.(Varnames_FBA{i}).x;
    stored_FBA_all(:,i) = tmp;
    clear tmp
end

end
round_solFBA_all = round(stored_FBA_all,5,'significant');
% substitute numerical artifacts by 0
upper_b = 1e-4;
lower_b = -1e-7;

mod_round_solFBA_all = round_solFBA_all;
mod_round_solFBA_all(mod_round_solFBA_all(:,:) < upper_b) = 0;

% generate new matrices without rows that are always 0
mean_solFBA_all = zeros(size(mod_round_solFBA_all,1),1);
for i = 1:size(mean_solFBA_all,1)
    mean_solFBA_all(i,1) = mean(mod_round_solFBA_all(i,:));
end
flag_solFBA_all = ones(size(mod_round_solFBA_all,1),1);
for i = 1:size(flag_solFBA_all,1)
    if mean_solFBA_all(i,1) == 0
        flag_solFBA_all(i,1) = 0;
    end         
end
clean_solFBA_all = zeros(nnz(flag_solFBA_all),size(mod_round_solFBA_all,2));
j = 1;
for i=1:size(mod_round_solFBA_all,1)
    if flag_solFBA_all(i,1) == 1
        clean_solFBA_all(j,:) = mod_round_solFBA_all(i,:);
        j = j + 1;
    end
end
clear j

%%% Now we analyse the similarities
% number of succesful runs
succ_runs = size(clean_solFBA_all,2)
num_combinations = 0;

for p = 1:succ_runs
num_combinations = num_combinations + p;
end
clear p

list_comb = zeros(num_combinations,2); % combination of this list against itself

p = 1;
r = 1;

while p <= succ_runs
    for q = p:succ_runs
        list_comb(r,1) = p;
        list_comb(r,2) = q;
        r = r + 1;
    end
p = p + 1;
end

clear p
clear q
clear r

% now we calculate the Jaccard index for every combination
p = 1;
for p = 1:num_combinations
    X = [clean_solFBA_all(:,list_comb(p,1)) clean_solFBA_all(:,list_comb(p,2))]';
    JD = pdist(X,'jaccard');
    JI = 1 - JD;
    list_comb(p,3) = JI;
        
    clear X
    clear JD
    clear JI
end
clear p

% average and std dev of similarities
p = 1;
for p = 1:num_combinations
    list_comb(p,4) = mean(list_comb(p,3));
    list_comb(p,5) = std(list_comb(p,3));
end
clear p

% calculate pairs of unique average-std dev values
stats_list_comb(:,1) = unique(list_comb(:,4),'stable');
stats_list_comb(:,2) = unique(list_comb(:,5),'stable');

%%% I need a symmetric matrix to identify identical vectors (i.e. tests that
% produce the same flux/metabolic space solution) -> average of JIs (not
% actual flux values!)
ave_matrix = zeros(succ_runs,succ_runs);
for i = 1:size(list_comb,1)
    ave_matrix(list_comb(i,1),list_comb(i,2)) = list_comb(i,4);
    ave_matrix(list_comb(i,2),list_comb(i,1)) = list_comb(i,4);
end

space_solutions = unique(ave_matrix(:,:),'rows','stable');
identification_test = zeros(succ_runs,size(space_solutions,1)); % rows are the test number, column the case of space solution

for i = 1:size(space_solutions,1)
    for j = 1:succ_runs
        tmp = space_solutions(i,:) - ave_matrix(j,:);
        if tmp == 0
            identification_test(j,i) = 1;
        end
        clear tmp
    end
end

% Unify results and prioritise
unique_solutions_num = size(identification_test,2);
unique_solutions = zeros(unique_solutions_num,2);
unique_solutions(:,1) = [1:unique_solutions_num]';

% now we need to select one example of each solution from
% identification_test
for i = 1:unique_solutions_num
    for j = 1:number_tests
        if identification_test(j,i) == 1
            unique_solutions(i,2) = j;            
        end
    end
end

% pool potential targets from unique solutions
pooled_list_potential_targets = [];
for i = 1:unique_solutions_num
    j = unique_solutions(i,2);
    tmp_list =  stored_common_chokepoints_lethal_consumption_lethal_production.(Varnames_common_chokepoints_lethal_consumption_lethal_production{j});
    pooled_list_potential_targets = [pooled_list_potential_targets; tmp_list];
end
unique_rxns = unique(cell2mat(pooled_list_potential_targets(:,1)));
sorted_list = sort(cell2mat(pooled_list_potential_targets(:,1)));
table = tabulate(sorted_list);
for i = 1:size(unique_rxns,1)
    tmp_i = unique_rxns(i,1);
    for j = 1:size(table,1)
        tmp_j = table(j,1);
        if tmp_i == tmp_j
            unique_rxns(i,2) = table(j,2);
            unique_rxns(i,3) = table(j,3);
        end    
    end
end
clear tmp_i
clear tmp_j

unique_rxns_cell = num2cell(zeros(size(unique_rxns,1),size(pooled_list_potential_targets,2)+2));
unique_rxns_cell(:,1:3) = num2cell(unique_rxns(:,1:3));
for i = 1:size(unique_rxns_cell,1)
    tmp_i = cell2mat(unique_rxns_cell(i,1));
    for j = 1:size(pooled_list_potential_targets,1)
        tmp_j = cell2mat(pooled_list_potential_targets(j,1));
        if tmp_i == tmp_j
            unique_rxns_cell(i,4:end) = pooled_list_potential_targets(j,2:end);
        end    
    end    
end

%save medoid_1_analysis.mat
save medoid_2_analysis.mat

%% Analyse important proteins and integrate the results for both medoids (MATLAB R2017b)
% need to analyse PageRank (essential and upregulated)
clear all
clc
 
% for medoid_1
load medoid_1_analysis
% extract non-redundant data
for i = 1:unique_solutions_num
    j = unique_solutions(i,2);
    important_protein.node{i} =  stored_node_names.(Varnames_nodenames{j});
    important_protein.pr{i} =  stored_node_pagerank.(Varnames_pagerank{j});
    important_protein.essential{i} = stored_essential_reactions.(Varnames_essential_reactions{j});
end
% identify essential reactions, then keep top 10% of PR
% need to only consider reactions with an associated gene
for i = 1:unique_solutions_num
    tmp_joint_clean = [];
    tmp_essential = important_protein.essential{i};
    tmp_node = important_protein.node{i};
    tmp_pr = num2cell(important_protein.pr{i});
    tmp_joint = [tmp_node tmp_pr];
    tmp_joint(:,3) = num2cell(zeros(size(tmp_joint,1),1));
    for j = 1:size(tmp_joint,1)
        tmp_rxn = tmp_joint{j,1};
        for k = 1:size(tmp_essential,1)
            tmp_rxn_ess = tmp_essential{k,2};
            tmp_gene = tmp_essential{k,5};
            if strcmp(tmp_rxn, tmp_rxn_ess) == 1 && isempty(tmp_gene) == 0
                tmp_joint{j,3} = tmp_gene;
                tmp_joint_clean = [tmp_joint_clean; tmp_joint(j,:)];
            end
        end
    end
    important_protein.relevant{i} = tmp_joint_clean;
    clear tmp_essential tmp_node tmp_pr tmp_joint tmp_joint_clean
end

% load upregulated genes, as obtained in R (tumor v healthy, regardless of the cluster)
upgenes_raw_table = readtable('topGenes_interest_def.txt');
upgenes_raw = table2cell(upgenes_raw_table);
% we add this info to the previous table: still a reaction-centric view
for m = 1:unique_solutions_num    
    tmp_relevant = important_protein.relevant{m}
    tmp_relevant(:,7) = num2cell(zeros);
    for i = 1:size(tmp_relevant,1)
        tmp_genes = tmp_relevant{i,3};
        for j = 1:size(tmp_genes,1)
            tmp_genes_j = tmp_genes{j,1};
            for k = 1:size(upgenes_raw)
                tmp_genes_k = upgenes_raw{k,1};
                if strcmp(tmp_genes_j,tmp_genes_k) == 1
                     tmp_relevant{i,4} = [tmp_relevant{i,4}; upgenes_raw(k,1)]; % gene list
                     tmp_relevant{i,5} = [tmp_relevant{i,5}; upgenes_raw(k,2)]; % names
                     tmp_relevant{i,6} = [tmp_relevant{i,6}; upgenes_raw(k,4)]; % logFC
                end            
            end
        end
    end
    tmp_relevant_up = [];
    for n = 1:size(tmp_relevant,1)
        if cell2mat(tmp_relevant(n,2)) > 0.001 && isempty(tmp_relevant{n,4}) == 0
            tmp_relevant_up = [tmp_relevant_up; tmp_relevant(n,:)];
        end
    end
    important_protein.relevant_up{m} = tmp_relevant_up;
end

% pool the results, identify unique ones and add info on frequency
pool_important = [];
for i = 1:unique_solutions_num 
    pool_important = [pool_important; important_protein.relevant_up{i}];
end
pool_important(:,2) = num2cell(ones(size(pool_important,1),1));
% add up frequencies for repeated reactions
[~, uniqueIdx] = unique(pool_important(:,1));
unique_pooled = pool_important(uniqueIdx,:);
for i = 1:size(unique_pooled,1)
    tmp_rxn = unique_pooled{i,1};
    for j =1:size(pool_important,1)
        tmp_rxn_pre = pool_important{j,1};
        if strcmp(tmp_rxn, tmp_rxn_pre) == 1
            unique_pooled(i,2) = num2cell(cell2mat(unique_pooled(i,2))+cell2mat(pool_important(j,2))); 
        end
        clear tmp_rxn_pre
    end
    [genes_un, ~] = unique(unique_pooled{i,5});
    unique_pooled{i,7} = genes_un;
    clear genes_un
end

% make it gene centric
% gene centric view: since a gene could be participating in more than one
% reaction, we need to change the perspective. We are using a "gene
% participation" metrics, where the overall frequency needs to be accounted
% for
unique_pooled(:,8) = unique_pooled(:,7)
for i = 1:size(unique_pooled,1)
    tmp_freq = unique_pooled(i,2);
    for j = 1:size(unique_pooled{i,6},1)
        unique_pooled{i,8}(j) = (tmp_freq);
    end
end
clear tmp_freq

gene_list_a = []
gene_list_b = []
gene_list_c = []
gene_list_d = []
for i = 1:size(unique_pooled,1)
    tmp_genes = unique_pooled{i,4};
    tmp_names = unique_pooled{i,5};
    tmp_FC = unique_pooled{i,6};
    tmp_freq = unique_pooled{i,8};
    gene_list_a = [gene_list_a; tmp_genes];
    gene_list_b = [gene_list_b; tmp_names];
    gene_list_c = [gene_list_c; tmp_FC];
    gene_list_d = [gene_list_d; tmp_freq];
end
gene_list = [gene_list_a gene_list_b gene_list_c gene_list_d];
[genes_cod_un, ~] = unique(gene_list(:,1), 'stable');
[genes_names_un, ~] = unique(gene_list(:,2), 'stable');
[genes_FC_un, ~] = unique(gene_list(:,3), 'stable');
gene_list_un = [genes_names_un genes_cod_un genes_FC_un];
gene_list_un_sort = sortrows(gene_list_un);

% now we calculate the gene participation by considering the accumulated
% frequencies
gene_list_un_sort(:,4) = num2cell(zeros(size(gene_list_un_sort,1),1));
for i = 1:size(gene_list_un_sort,1)
    tmp_gene = gene_list_un_sort{i,1};
    for j = 1:size(gene_list,1)
        tmp_gene_all = gene_list{j,2};
        if strcmp(tmp_gene,tmp_gene_all) == 1
            gene_list_un_sort(i,4) = num2cell(cell2mat(gene_list_un_sort(i,4))+cell2mat(gene_list(j,4)));
        end
    end
    clear tmp_gene
end
gene_list_un_sort_pr_medoid_1 = gene_list_un_sort
clearvars -except gene_list_un_sort_pr_medoid_1

% for medoid_2
load medoid_2_analysis
% extract non-redundant data
for i = 1:unique_solutions_num
    j = unique_solutions(i,2);
    important_protein.node{i} =  stored_node_names.(Varnames_nodenames{j});
    important_protein.pr{i} =  stored_node_pagerank.(Varnames_pagerank{j});
    important_protein.essential{i} = stored_essential_reactions.(Varnames_essential_reactions{j});
end
% identify essential reactions, then keep top 10% of PR
% need to only consider reactions with an associated gene
for i = 1:unique_solutions_num
    tmp_joint_clean = [];
    tmp_essential = important_protein.essential{i};
    tmp_node = important_protein.node{i};
    tmp_pr = num2cell(important_protein.pr{i});
    tmp_joint = [tmp_node tmp_pr];
    tmp_joint(:,3) = num2cell(zeros(size(tmp_joint,1),1));
    for j = 1:size(tmp_joint,1)
        tmp_rxn = tmp_joint{j,1};
        for k = 1:size(tmp_essential,1)
            tmp_rxn_ess = tmp_essential{k,2};
            tmp_gene = tmp_essential{k,5};
            if strcmp(tmp_rxn, tmp_rxn_ess) == 1 && isempty(tmp_gene) == 0
                tmp_joint{j,3} = tmp_gene;
                tmp_joint_clean = [tmp_joint_clean; tmp_joint(j,:)];
            end
        end
    end
    important_protein.relevant{i} = tmp_joint_clean;
    clear tmp_essential tmp_node tmp_pr tmp_joint tmp_joint_clean
end

% load upregulated genes, as obtained in R (tumor v healthy, regardless of the cluster)
upgenes_raw_table = readtable('topGenes_interest_def.txt');
upgenes_raw = table2cell(upgenes_raw_table);
% we add this info to the previous table: still a reaction-centric view
for m = 1:unique_solutions_num    
    tmp_relevant = important_protein.relevant{m}
    tmp_relevant(:,7) = num2cell(zeros);
    for i = 1:size(tmp_relevant,1)
        tmp_genes = tmp_relevant{i,3};
        for j = 1:size(tmp_genes,1)
            tmp_genes_j = tmp_genes{j,1};
            for k = 1:size(upgenes_raw)
                tmp_genes_k = upgenes_raw{k,1};
                if strcmp(tmp_genes_j,tmp_genes_k) == 1
                     tmp_relevant{i,4} = [tmp_relevant{i,4}; upgenes_raw(k,1)]; % gene list
                     tmp_relevant{i,5} = [tmp_relevant{i,5}; upgenes_raw(k,2)]; % names
                     tmp_relevant{i,6} = [tmp_relevant{i,6}; upgenes_raw(k,4)]; % logFC
                end            
            end
        end
    end
    tmp_relevant_up = [];
    for n = 1:size(tmp_relevant,1)
        if cell2mat(tmp_relevant(n,2)) > 0.001 && isempty(tmp_relevant{n,4}) == 0
            tmp_relevant_up = [tmp_relevant_up; tmp_relevant(n,:)];
        end
    end
    important_protein.relevant_up{m} = tmp_relevant_up;
end

% pool the results, identify unique ones and add info on frequency
pool_important = [];
for i = 1:unique_solutions_num 
    pool_important = [pool_important; important_protein.relevant_up{i}];
end
pool_important(:,2) = num2cell(ones(size(pool_important,1),1));
% add up frequencies for repeated reactions
[~, uniqueIdx] = unique(pool_important(:,1));
unique_pooled = pool_important(uniqueIdx,:);
for i = 1:size(unique_pooled,1)
    tmp_rxn = unique_pooled{i,1};
    for j =1:size(pool_important,1)
        tmp_rxn_pre = pool_important{j,1};
        if strcmp(tmp_rxn, tmp_rxn_pre) == 1
            unique_pooled(i,2) = num2cell(cell2mat(unique_pooled(i,2))+cell2mat(pool_important(j,2))); 
        end
        clear tmp_rxn_pre
    end
    [genes_un, ~] = unique(unique_pooled{i,5});
    unique_pooled{i,7} = genes_un;
    clear genes_un
end

% make it gene centric
% gene centric view: since a gene could be participating in more than one
% reaction, we need to change the perspective. We are using a "gene
% participation" metrics, where the overall frequency needs to be accounted
% for
unique_pooled(:,8) = unique_pooled(:,7)
for i = 1:size(unique_pooled,1)
    tmp_freq = unique_pooled(i,2);
    for j = 1:size(unique_pooled{i,6},1)
        unique_pooled{i,8}(j) = (tmp_freq);
    end
end
clear tmp_freq

gene_list_a = []
gene_list_b = []
gene_list_c = []
gene_list_d = []
for i = 1:size(unique_pooled,1)
    tmp_genes = unique_pooled{i,4};
    tmp_names = unique_pooled{i,5};
    tmp_FC = unique_pooled{i,6};
    tmp_freq = unique_pooled{i,8};
    gene_list_a = [gene_list_a; tmp_genes];
    gene_list_b = [gene_list_b; tmp_names];
    gene_list_c = [gene_list_c; tmp_FC];
    gene_list_d = [gene_list_d; tmp_freq];
end
gene_list = [gene_list_a gene_list_b gene_list_c gene_list_d];
[genes_cod_un, ~] = unique(gene_list(:,1), 'stable');
[genes_names_un, ~] = unique(gene_list(:,2), 'stable');
[genes_FC_un, ~] = unique(gene_list(:,3), 'stable');
gene_list_un = [genes_names_un genes_cod_un genes_FC_un];
gene_list_un_sort = sortrows(gene_list_un);

% now we calculate the gene participation by considering the accumulated
% frequencies
gene_list_un_sort(:,4) = num2cell(zeros(size(gene_list_un_sort,1),1));
for i = 1:size(gene_list_un_sort,1)
    tmp_gene = gene_list_un_sort{i,1};
    for j = 1:size(gene_list,1)
        tmp_gene_all = gene_list{j,2};
        if strcmp(tmp_gene,tmp_gene_all) == 1
            gene_list_un_sort(i,4) = num2cell(cell2mat(gene_list_un_sort(i,4))+cell2mat(gene_list(j,4)));
        end
    end
    clear tmp_gene
end

gene_list_un_sort_pr_medoid_2 = gene_list_un_sort
clearvars -except gene_list_un_sort_pr_medoid_1 gene_list_un_sort_pr_medoid_2

% pool both lists
pooled_both_pr = [gene_list_un_sort_pr_medoid_1; gene_list_un_sort_pr_medoid_2];
% add up frequencies for repeated reactions
[~, uniqueIdx] = unique(pooled_both_pr(:,1));
unique_pooled = pooled_both_pr(uniqueIdx,:);
for i = 1:size(unique_pooled,1)
    tmp_rxn = unique_pooled{i,1};
    for j =1:size(pooled_both_pr,1)
        tmp_rxn_pre = pooled_both_pr{j,1};
        if strcmp(tmp_rxn, tmp_rxn_pre) == 1
            unique_pooled(i,4) = num2cell(cell2mat(unique_pooled(i,4))+cell2mat(pooled_both_pr(j,4))); 
        end
        clear tmp_rxn_pre
    end
 
end
important_genes = unique_pooled;
clearvars -except important_genes

%%% now we analyse the double chokepoints
load medoid_1_analysis
unique_medoid_1 = unique_rxns_cell
lower_thres_medoid_1 = floor(unique_solutions_num/2)
clearvars -except unique_medoid_1 important_genes lower_thres_medoid_1
load medoid_2_analysis
unique_medoid_2 = unique_rxns_cell
lower_thres_medoid_2 = floor(unique_solutions_num/2)
clearvars -except unique_medoid_1 unique_medoid_2 important_genes lower_thres_medoid_1 lower_thres_medoid_2
% consider the lower threshold
indx_med_1 = zeros(size(unique_medoid_1,1),1)
for i = 1:size(unique_medoid_1,1)
    tmp_f = cell2mat(unique_medoid_1(i,2));
    tmp_g = isempty(unique_medoid_1{i,8});
    if tmp_f >= lower_thres_medoid_1 & tmp_g == 0
        indx_med_1(i,1) = 1;
    end
    clear tmp_f tmp_g
end
rnx_interest_1 = nnz(indx_med_1)
indx_med_2 = zeros(size(unique_medoid_2,1),1)
for i = 1:size(unique_medoid_2,1)
    tmp_f = cell2mat(unique_medoid_2(i,2));
    tmp_g = isempty(unique_medoid_2{i,8});
    if tmp_f >= lower_thres_medoid_2 & tmp_g == 0
        indx_med_2(i,1) = 1;
    end
    clear tmp_f tmp_g
end
rnx_interest_2 = nnz(indx_med_2)

% get clean tables
unique_medoid_1_clean = num2cell(zeros(rnx_interest_1,size(unique_medoid_1,2)))
unique_medoid_2_clean = num2cell(zeros(rnx_interest_2,size(unique_medoid_2,2)))
j = 1
for i = 1:size(indx_med_1,1)
    if indx_med_1(i,1) == 1
        unique_medoid_1_clean(j,:) = unique_medoid_1(i,:);
        j = j+1;
    end
end
j = 1
for i = 1:size(indx_med_2,1)
    if indx_med_2(i,1) == 1
        unique_medoid_2_clean(j,:) = unique_medoid_2(i,:);
        j = j+1
    end
end

% pool tables: rxn index is not the same for both contextualised models.
% There are reactions appearing in both lists, so we are going to add up
% the frequency values. Since the rxnGeneMat were adapted earlier on, then
% the gene list associated to each reactions may be different for each list
% (e.g. lacking 1 gene)
pooled_rxns_pre_1 = unique_medoid_1_clean(:,2);
pooled_rxns_pre_1 = [pooled_rxns_pre_1 unique_medoid_1_clean(:,4)];
pooled_rxns_pre_1 = [pooled_rxns_pre_1 unique_medoid_1_clean(:,8)];
pooled_rxns_pre_2 = unique_medoid_2_clean(:,2);
pooled_rxns_pre_2 = [pooled_rxns_pre_2 unique_medoid_2_clean(:,4)];
pooled_rxns_pre_2 = [pooled_rxns_pre_2 unique_medoid_2_clean(:,8)];
pooled_rxns_pre = [pooled_rxns_pre_1; pooled_rxns_pre_2];
% when the same reaction pop, pick the one with the highest number of
% associated genes
for i = 1:size(pooled_rxns_pre,1)
    pooled_rxns_pre(i,4) = num2cell(size(pooled_rxns_pre{i,3},1));
end

% add up frequencies for repeated reactions
[~, uniqueIdx] = unique(pooled_rxns_pre(:,2));
unique_pooled = pooled_rxns_pre(uniqueIdx,:);
unique_pooled(:,5) = num2cell(zeros);
for i = 1:size(unique_pooled,1)
    tmp_rxn = unique_pooled{i,2};
    for j =1:size(pooled_rxns_pre,1)
        tmp_rxn_pre = pooled_rxns_pre{j,2};
        if strcmp(tmp_rxn, tmp_rxn_pre) == 1
            unique_pooled(i,5) = num2cell(cell2mat(unique_pooled(i,5))+cell2mat(pooled_rxns_pre(j,1)));
            unique_pooled{i,6} = [unique_pooled{i,3}; pooled_rxns_pre{j,3}];
        end
        clear tmp_rxn_pre
    end
    [genes_un, ~] = unique(unique_pooled{i,6});
    unique_pooled{i,7} = genes_un;
    clear genes_un
end

unique_pooled_clean_a = unique_pooled(:,2);
unique_pooled_clean_a(:,2) = unique_pooled(:,5);
unique_pooled_clean_a(:,3) = unique_pooled(:,7);
unique_pooled_clean_a(:,7) = num2cell(zeros);
% load upregulated genes, as obtained in R (tumor v healthy, regardless of the cluster)
upgenes_raw_table = readtable('topGenes_interest_def.txt');
upgenes_raw = table2cell(upgenes_raw_table);
% we add this info to the previous table: still a reaction-centric view
for i = 1:size(unique_pooled_clean_a,1)
    tmp_genes = unique_pooled_clean_a{i,3};
    for j = 1:size(tmp_genes,1)
        tmp_genes_j = tmp_genes{j,1};
        for k = 1:size(upgenes_raw)
            tmp_genes_k = upgenes_raw{k,1};
            if strcmp(tmp_genes_j,tmp_genes_k) == 1
                 unique_pooled_clean_a{i,4} = [unique_pooled_clean_a{i,4}; upgenes_raw(k,1)]; % gene list
                 unique_pooled_clean_a{i,5} = [unique_pooled_clean_a{i,5}; upgenes_raw(k,2)]; % names
                 unique_pooled_clean_a{i,6} = [unique_pooled_clean_a{i,6}; upgenes_raw(k,4)]; % logFC
            end            
        end
    end
end

% keep those which are upregulated, taking the subset of the gene list
% referring to upregulated genes (instead of the full list from the GSM)
unique_pooled_clean_b = [unique_pooled_clean_a(:,1) unique_pooled_clean_a(:,2) unique_pooled_clean_a(:,4) unique_pooled_clean_a(:,5) unique_pooled_clean_a(:,6)];
unique_pooled_clean_c = [];
for i = 1:size(unique_pooled_clean_b,1)
    if isempty(unique_pooled_clean_b{i,3}) == 0
        unique_pooled_clean_c = [unique_pooled_clean_c; unique_pooled_clean_b(i,:)];
    end
end
unique_pooled_clean_c(:,6) = unique_pooled_clean_c(:,3);
for i = 1:size(unique_pooled_clean_c,1)
    tmp_freq = unique_pooled_clean_c(i,2);
    for j = 1:size(unique_pooled_clean_c{i,6},1)
        unique_pooled_clean_c{i,6}(j) = (tmp_freq);
    end
end
clear tmp_freq

% gene centric view: since a gene could be participating in more than one
% reaction, we need to change the perspective. We are using a "gene
% participation" metrics, where the overall frequency needs to be accounted
% for
accum_size = 0;
for i = 1:size(unique_pooled_clean_c,1)
    tmp_size = size(unique_pooled_clean_c{i,3},1);
    accum_size = accum_size + tmp_size;
end
clear tmp_size accum_size

gene_list_a = []
gene_list_b = []
gene_list_c = []
gene_list_d = []
for i = 1:size(unique_pooled_clean_c,1)
    tmp_genes = unique_pooled_clean_c{i,3};
    tmp_names = unique_pooled_clean_c{i,4};
    tmp_FC = unique_pooled_clean_c{i,5};
    tmp_freq = unique_pooled_clean_c{i,6};
    gene_list_a = [gene_list_a; tmp_genes];
    gene_list_b = [gene_list_b; tmp_names];
    gene_list_c = [gene_list_c; tmp_FC];
    gene_list_d = [gene_list_d; tmp_freq];
end
gene_list = [gene_list_a gene_list_b gene_list_c gene_list_d];
[genes_cod_un, ~] = unique(gene_list(:,1), 'stable');
[genes_names_un, ~] = unique(gene_list(:,2), 'stable');
[genes_FC_un, ~] = unique(gene_list(:,3), 'stable');
gene_list_un = [genes_names_un genes_cod_un genes_FC_un];
gene_list_un_sort = sortrows(gene_list_un);

% now we calculate the gene participation by considering the accumulated
% frequencies
gene_list_un_sort(:,4) = num2cell(zeros(size(gene_list_un_sort,1),1));
for i = 1:size(gene_list_un_sort,1)
    tmp_gene = gene_list_un_sort{i,1}
    for j = 1:size(gene_list,1)
        tmp_gene_all = gene_list{j,2};
        if strcmp(tmp_gene,tmp_gene_all) == 1
            gene_list_un_sort(i,4) = num2cell(cell2mat(gene_list_un_sort(i,4))+cell2mat(gene_list(j,4)));
        end
    end
    clear tmp_gene
end

double_chokepoints = gene_list_un_sort
clearvars -except important_genes double_chokepoints

%%% combine both lists
full_genes = [double_chokepoints(:,1:3); important_genes(:,1:3)]
[full_genes_un, ~] = unique(full_genes(:,1));
for i = 1:size(full_genes_un,1)
    tmp_gene = full_genes_un{i,1};
    for j = 1:size(full_genes)
        tmp_gene_un = full_genes{j,1};
        if strcmp(tmp_gene,tmp_gene_un) == 1
            full_genes_un(i,2) = full_genes(j,2);
            full_genes_un(i,3) = full_genes(j,3);
        end
    end
end

%%% include gene participation values from each result (independently normalised to %)  
full_genes_un(:,4) = num2cell(zeros(size(full_genes_un,1),1)); % double choke points
full_genes_un(:,5) = num2cell(zeros(size(full_genes_un,1),1)); % important proteins
% gather values
for i = 1:size(full_genes_un,1)
    tmp_gene = full_genes_un{i,1};
    for j = 1:size(double_chokepoints,1)
        tmp_gene_un = double_chokepoints{j,1};
        if strcmp(tmp_gene,tmp_gene_un) == 1
            full_genes_un(i,4) = double_chokepoints(j,4);
        end
    end
end
for i = 1:size(full_genes_un,1)
    tmp_gene = full_genes_un{i,1};
    for j = 1:size(important_genes,1)
        tmp_gene_un = important_genes{j,1};
        if strcmp(tmp_gene,tmp_gene_un) == 1
            full_genes_un(i,5) = important_genes(j,4);
        end
    end
end
% normalise gene participation values
tmp_max_1 = max(cell2mat(full_genes_un(:,4)));
tmp_max_2 = max(cell2mat(full_genes_un(:,5)));
for i = 1:size(full_genes_un,1)
    full_genes_un(i,4) = num2cell(round(cell2mat(full_genes_un(i,4))*100/tmp_max_1),2);
    full_genes_un(i,5) = num2cell(round(cell2mat(full_genes_un(i,5))*100/tmp_max_2),2);
end
clear tmp_max_1 tmp_max_2

% genetic dependencies: there may be more than one entry for each gene,
% since info comes from different studies. 
% if there is one "Combined RNAi (Broad, Novartis, Marcotte)" pick it
full_genes_un(:,6) = num2cell(zeros(size(full_genes_un,1),1)); % 
full_genes_un(:,7) = num2cell(zeros(size(full_genes_un,1),1)); % 
full_genes_un(:,8) = num2cell(zeros(size(full_genes_un,1),1)); % 
depmap = importdata('depmap_GBM.xlsx');
data_name = depmap.textdata(2:end,3); % it includes headers!
data_db = depmap.textdata(2:end,4);
data_t_statistics = num2cell(depmap.data(:,1));
data_p_value = num2cell(depmap.data(:,2));
for i = 1:size(full_genes_un,1)
    tmp_gene = full_genes_un{i,1};
    tmp_db = []
    tmp_t_stat = []
    tmp_p_val = []
    for j = 1:size(data_name,1)
        tmp_gene_db = data_name{j,1};
        if strcmp(tmp_gene, tmp_gene_db) == 1
            tmp_db = [tmp_db; data_db(j)];
            tmp_t_stat = [tmp_t_stat; data_t_statistics(j)];
            tmp_p_val = [tmp_p_val; data_p_value(j)];
        end
    end
    full_genes_un{i,6} = tmp_db;
    full_genes_un{i,7} = tmp_t_stat;
    full_genes_un{i,8} = tmp_p_val;
end

% pick the 'combined' if available
for i = 1:size(full_genes_un,1)
    if isempty(full_genes_un{i,6}) == 0
        tmp_size = size(full_genes_un{i,6},1);
        if tmp_size == 1
            full_genes_un(i,9) = full_genes_un{i,6};
            full_genes_un(i,10) = full_genes_un{i,7};
            full_genes_un(i,11) = full_genes_un{i,8}; 
        else
            for j = 1:tmp_size
               tmp_data = full_genes_un{i,6}(j);
               if contains(tmp_data,'Combined') == 1
                    full_genes_un{i,9} = full_genes_un{i,6}(j);
                    full_genes_un(i,10) = full_genes_un{i,7}(j);
                    full_genes_un(i,11) = full_genes_un{i,8}(j); 
                    j = tmp_size;             
               end
            end
        end
    end
end

% rewrite the matrix
full_genes_un(:,6) = full_genes_un(:,9);
full_genes_un(:,7) = full_genes_un(:,10);
full_genes_un(:,8) = full_genes_un(:,11);

% drug data: interactions.tsv from https://dgidb.org/downloads
interactions_DGIdb = tdfread('interactions_DGIdb_Nov_2020.tsv','\t')

% get a gene-centric view, with one row per drug (therefore repeated genes
% will appear)
interactions_gene = []
interactions_source = []
interactions_type = []
interactions_drug_name = []
interactions_drug_concept_id = []
interactions_PMIDs = []

for i = 1:size(interactions_DGIdb.gene_name,1)
    interactions_gene = [interactions_gene; strtrim(convertCharsToStrings(interactions_DGIdb.gene_name(i,:)))];
    interactions_source = [interactions_source; strtrim(convertCharsToStrings(interactions_DGIdb.interaction_claim_source(i,:)))];
    interactions_type = [interactions_type; strtrim(convertCharsToStrings(interactions_DGIdb.interaction_types(i,:)))];
    interactions_drug_name = [interactions_drug_name; strtrim(convertCharsToStrings(interactions_DGIdb.drug_name(i,:)))];
    interactions_drug_concept_id = [interactions_drug_concept_id; strtrim(convertCharsToStrings(interactions_DGIdb.drug_concept_id(i,:)))];
    interactions_PMIDs = [interactions_PMIDs; strtrim(convertCharsToStrings(interactions_DGIdb.PMIDs(i,:)))];
end

genes_list_DGIdb_pre(:,1) = num2cell((1:size(interactions_DGIdb.gene_name,1))');
genes_list_DGIdb_pre(:,2) = cellstr(interactions_gene);
genes_list_DGIdb_pre(:,10) = cellstr(interactions_source);
genes_list_DGIdb_pre(:,11) = cellstr(interactions_type);
genes_list_DGIdb_pre(:,12) = cellstr(interactions_drug_name);
genes_list_DGIdb_pre(:,13) = cellstr(interactions_drug_concept_id);
genes_list_DGIdb_pre(:,14) = cellstr(interactions_PMIDs);
for i = 1:size(genes_list_DGIdb_pre,1)
    tmp_gene = genes_list_DGIdb_pre{i,2};
    for j = 1:size(full_genes_un,1)
        tmp_gene_un = full_genes_un{j,1};
        if strcmp(tmp_gene, tmp_gene_un) == 1
            genes_list_DGIdb_pre(i,3) = full_genes_un(j,2);
            genes_list_DGIdb_pre(i,4) = full_genes_un(j,3);
            genes_list_DGIdb_pre(i,5) = full_genes_un(j,4);
            genes_list_DGIdb_pre(i,6) = full_genes_un(j,5);
            genes_list_DGIdb_pre(i,7) = full_genes_un(j,6);
            genes_list_DGIdb_pre(i,8) = full_genes_un(j,7);
            genes_list_DGIdb_pre(i,9) = full_genes_un(j,8);
        end
    end
end

% keep only entries for genes of interest, which have a known interaction
% with a drug
genes_list_DGIdb_clean = []
for i = 1:size(genes_list_DGIdb_pre,1)
    if isempty(genes_list_DGIdb_pre{i,3}) == 0 && isempty(genes_list_DGIdb_pre{i,13}) == 0
        tmp_row = genes_list_DGIdb_pre(i,:);
        genes_list_DGIdb_clean = [genes_list_DGIdb_clean; tmp_row];
        clear tmp_row
    end
end
genes_list_DGIdb_clean = sortrows(genes_list_DGIdb_clean,2)

% obtain a list of unique drugs, with info of their participation in genes_list_DGIdb_clean
% since different DBs were use, there may duplicates (X gene, Y drug),so
% this needs to be considered too. 
genes_list_DGIdb_clean_nnrd = genes_list_DGIdb_clean(1,:)
for i = 1:size(genes_list_DGIdb_clean,1)
    repeated = 0;
    tmp_gene = genes_list_DGIdb_clean{i,2};
    tmp_drug = genes_list_DGIdb_clean{i,12};
    for j = 1:size(genes_list_DGIdb_clean_nnrd,1)
        tmp_gene_nnrd = genes_list_DGIdb_clean_nnrd{j,2};
        tmp_drug_nnrd = genes_list_DGIdb_clean_nnrd{j,12};
        if strcmp(tmp_gene,tmp_gene_nnrd) == 1 && strcmp(tmp_drug,tmp_drug_nnrd) == 1
            repeated = repeated + 1
        end
    end
    if repeated == 0
        genes_list_DGIdb_clean_nnrd = [genes_list_DGIdb_clean_nnrd; genes_list_DGIdb_clean(i,:)];
    end
    clear repeated
end

% get a drug-wise view, considering an accumulated gene participation (which was previously normalised to %)
[drug_wise_DGIdb, ~] = unique(genes_list_DGIdb_clean_nnrd(:,12));
drug_wise_DGIdb(:,2) = {[]};
drug_wise_DGIdb(:,3) = {'NB'};
drug_wise_DGIdb(:,4) = {[]};
drug_wise_DGIdb(:,5) = {[]};
drug_wise_DGIdb(:,6) = {[]};

for i = 1:size(drug_wise_DGIdb,1)
    tmp_gene = [];
    tmp_gene_partic_choke = 0;
    tmp_gene_partic_pagerank = 0;
    tmp_drug = drug_wise_DGIdb{i,1};
    for j = 1:size(genes_list_DGIdb_clean_nnrd)
        tmp_drug_j = genes_list_DGIdb_clean_nnrd{j,12};
        if strcmp(tmp_drug,tmp_drug_j) == 1
           tmp_gene = [tmp_gene; genes_list_DGIdb_clean_nnrd(j,2)];
           tmp_gene_partic_choke = tmp_gene_partic_choke + cell2mat(genes_list_DGIdb_clean_nnrd(j,5));
           tmp_gene_partic_pagerank = tmp_gene_partic_pagerank + cell2mat(genes_list_DGIdb_clean_nnrd(j,6));
        end
        drug_wise_DGIdb{i,2} = tmp_gene;

        drug_wise_DGIdb(i,4) = num2cell(tmp_gene_partic_choke);
        drug_wise_DGIdb(i,5) = num2cell(tmp_gene_partic_pagerank);
        drug_wise_DGIdb(i,6) = num2cell(tmp_gene_partic_choke+tmp_gene_partic_pagerank);
    end
end

% add the type of interaction
drug_wise_DGIdb(:,3) = drug_wise_DGIdb(:,2);
for i = 1:size(drug_wise_DGIdb,1)
    for j = 1:size(drug_wise_DGIdb{i,3},1)
        drug_wise_DGIdb{i,3}(j,:) = cellstr('NB');
    end
end
sto_inter = [];
for i = 1:size(drug_wise_DGIdb,1)
    tmp_drug = drug_wise_DGIdb{i,1};
    for j = 1:size(drug_wise_DGIdb{i,2},1)
       tmp_gene = drug_wise_DGIdb{i,2}(j);
       for k = 1:size(genes_list_DGIdb_clean_nnrd,1) 
           tmp_drug_nnr = genes_list_DGIdb_clean_nnrd{k,12};
           tmp_gene_nnr = genes_list_DGIdb_clean_nnrd{k,2};
           if strcmp(tmp_drug,tmp_drug_nnr) == 1 && strcmp(tmp_gene, tmp_gene_nnr) == 1
               drug_wise_DGIdb{i,3}(j,:) = cellstr(genes_list_DGIdb_clean_nnrd{k,11});
               sto_inter = [sto_inter; cellstr(genes_list_DGIdb_clean_nnrd{k,11})];
           end
       end
    end
end

% USE MATLABR2017b onwards!
% In column 14 there's the PMID from the DGIdb. Now we have a list of
% unique drugs (drug_wise), where the genes they interact with are in the
% 2nd column. We add a 3rd column for the automated results in pubmed of
% searching 'glioblastoma + DRUG', and a 4th column for 'cancer + DRUG'.
% Entries in #3 validate the capacity of identifying known/studied drug
% targets and drugs; those not appearing in #3 but appearing in #4 are of
% interest.
% only consider publications from the last 10 years
% if the combination yields nothing (because there are no papers for the
% drug) results will be just for GBM. Therefore, also get results from the
% drug by itself, and discard those for which no result was retrieved

% if the server tries to kick me out, pause the query for a minute and
% resume
 pause('on');
 n = 60;

drug_wise_DGIdb(:,7) = {[]};
drug_wise_DGIdb(:,8) = {[]};
drug_wise_DGIdb(:,9) = {[]};

secondTerm = ' '
for i = 1:size(drug_wise_DGIdb,1)
    try
    searchterm = drug_wise_DGIdb{i,1};
    [PMID] = getpubmed(searchterm, secondTerm);
    drug_wise_DGIdb{i,7} = PMID;  
    clear searchterm PMID
    catch
        pause(n)
        i = i-1
    end
end
clear i

secondTerm = 'glioblastoma'
for i = 1:size(drug_wise_DGIdb,1)
    try
    searchterm = drug_wise_DGIdb{i,1};
    [PMID] = getpubmed(searchterm, secondTerm);
    drug_wise_DGIdb{i,8} = PMID;  
    clear searchterm PMID
    catch
        pause(n)
        i = i-1
    end
end
clear i

secondTerm = 'cancer'
for i = 1:size(drug_wise_DGIdb,1)
    try
    searchterm = drug_wise_DGIdb{i,1};
    [PMID] = getpubmed(searchterm, secondTerm);
    drug_wise_DGIdb{i,9} = PMID;  
    clear searchterm PMID
    catch
        pause(n)
        i = i-1
    end
end

% identify drugs for which there were results for the drug itself and no PMIDs were recovered for 'GBM + drug', but
% there were results for 'cancer + drug', the point being it is more likely
% for a drug known to work for another cancer to work for GBM. In
% particular, only considercases where at least 8 papers where found (when
% it says 12, there may be more, since only takes the first page)
for i = 1:size(drug_wise_DGIdb,1)
    tmp_papers = size(drug_wise_DGIdb{i,9},2);
    if isempty(drug_wise_DGIdb{i,7}) == 0 && isempty(drug_wise_DGIdb{i,8}) == 1 && tmp_papers >= 0
        drug_wise_DGIdb(i,10) = num2cell(tmp_papers);
    end
end

drug_wise_interest_pre = []
for i = 1:size(drug_wise_DGIdb,1)
    if drug_wise_DGIdb{i,10} >= 9
        drug_wise_interest_pre = [drug_wise_interest_pre; drug_wise_DGIdb(i,:)];
    end
end

% create a cell per PMID so we can manage it manually later
entries = num2cell(zeros(size(drug_wise_interest_pre,1),10));
drug_wise_interest_pre = [drug_wise_interest_pre entries];
for i = 1:size(drug_wise_interest_pre,1)
    tmp_num = size(drug_wise_interest_pre{i,9},2);
    for j = 1:tmp_num
        drug_wise_interest_pre{i,10+j} = drug_wise_interest_pre{i,9}{j};    
    end
end

% flag whenever there is at least one hit for inhibitor/allosteric
% modulator
drug_wise_interest_pre(:,23) = {[]};
tmp_inhibitor = 'inhibitor';
tmp_allosteric = 'allosteric modulator';
for i = 1:size(drug_wise_interest_pre,1)
    for j = 1:size(drug_wise_interest_pre{i,3},1)
        tmp_inter = drug_wise_interest_pre{i,3}{j,:};
        if strcmp(tmp_inter,tmp_inhibitor) == 1 || strcmp(tmp_inter,tmp_allosteric) == 1
            drug_wise_interest_pre(i,23) = num2cell(1);
        end
    end
end

drug_wise_interest_DGI = num2cell(zeros(size(drug_wise_interest_pre,1),1));
drug_wise_interest_DGI(:,1) = ({'DGI'});
drug_wise_interest_DGI = [drug_wise_interest_DGI drug_wise_interest_pre];

%%% drug data: https://ctdbase.org/tools/batchQuery.go
% (already filtered, only human and for the unique genes from full_genes)
% the file was too big to handle automatically
CTD_data = table2cell(readtable('CTD_gene_chem_inter.txt'));
 
genes_list_CTDdb_pre(:,1) = num2cell((1:size(CTD_data,1))');
genes_list_CTDdb_pre(:,2) = (CTD_data(:,1));
genes_list_CTDdb_pre(:,10) = (CTD_data(:,3)); % drug
genes_list_CTDdb_pre(:,11) = (CTD_data(:,4)); % Chemical
genes_list_CTDdb_pre(:,12) = (CTD_data(:,5)); % CasRN
genes_list_CTDdb_pre(:,13) = (CTD_data(:,6)); % interaction
genes_list_CTDdb_pre(:,14) = (CTD_data(:,7)); % interaction action
genes_list_CTDdb_pre(:,15) = (CTD_data(:,8)); % PMID
for i = 1:size(genes_list_CTDdb_pre,1)
    tmp_gene = genes_list_CTDdb_pre{i,2};
    for j = 1:size(full_genes_un,1)
        tmp_gene_un = full_genes_un{j,1};
        if strcmp(tmp_gene, tmp_gene_un) == 1
            genes_list_CTDdb_pre(i,3) = full_genes_un(j,2);
            genes_list_CTDdb_pre(i,4) = full_genes_un(j,3);
            genes_list_CTDdb_pre(i,5) = full_genes_un(j,4);
            genes_list_CTDdb_pre(i,6) = full_genes_un(j,5);
            genes_list_CTDdb_pre(i,7) = full_genes_un(j,6);
            genes_list_CTDdb_pre(i,8) = full_genes_un(j,7);
            genes_list_CTDdb_pre(i,9) = full_genes_un(j,8);
        end
    end
end
% since different DBs were use, there may duplicates (X gene, Y drug),so
% this needs to be considered too. 
genes_list_CTD_clean_nnrd = genes_list_CTDdb_pre(1,:)
for i = 1:size(genes_list_CTDdb_pre,1)
    repeated = 0;
    tmp_gene = genes_list_CTDdb_pre{i,2};
    tmp_drug = genes_list_CTDdb_pre{i,11};
    for j = 1:size(genes_list_CTD_clean_nnrd,1)
        tmp_gene_nnrd = genes_list_CTD_clean_nnrd{j,2};
        tmp_drug_nnrd = genes_list_CTD_clean_nnrd{j,11};
        if strcmp(tmp_gene,tmp_gene_nnrd) == 1 && strcmp(tmp_drug,tmp_drug_nnrd) == 1
            repeated = repeated + 1
        end
    end
    if repeated == 0
        genes_list_CTD_clean_nnrd = [genes_list_CTD_clean_nnrd; genes_list_CTDdb_pre(i,:)];
    end
    clear repeated
end

% get a drug-wise view, considering an accumulated gene participation (which was previously normalised to %)
[drug_wise_CTDdb, ~] = unique(genes_list_CTD_clean_nnrd(:,10));
drug_wise_CTDdb(:,2) = {[]};
drug_wise_CTDdb(:,3) = {'NB'};
drug_wise_CTDdb(:,4) = {[]};
drug_wise_CTDdb(:,5) = {[]};
drug_wise_CTDdb(:,6) = {[]};

if isempty(drug_wise_CTDdb{1,1}) == 1
    drug_wise_CTDdb = drug_wise_CTDdb(2:end,:);
end
for i = 1:size(drug_wise_CTDdb,1)
    tmp_gene = [];
    tmp_gene_partic_choke = 0;
    tmp_gene_partic_pagerank = 0;
    tmp_drug = drug_wise_CTDdb{i,1};
    for j = 1:size(genes_list_CTD_clean_nnrd)
        tmp_drug_j = genes_list_CTD_clean_nnrd{j,10};
        if strcmp(tmp_drug,tmp_drug_j) == 1
           tmp_gene = [tmp_gene; genes_list_CTD_clean_nnrd(j,2)];
           tmp_gene_partic_choke = tmp_gene_partic_choke + cell2mat(genes_list_CTD_clean_nnrd(j,5));
           tmp_gene_partic_pagerank = tmp_gene_partic_pagerank + cell2mat(genes_list_CTD_clean_nnrd(j,6));
        end
        drug_wise_CTDdb{i,2} = tmp_gene;

        drug_wise_CTDdb(i,4) = num2cell(tmp_gene_partic_choke);
        drug_wise_CTDdb(i,5) = num2cell(tmp_gene_partic_pagerank);
        drug_wise_CTDdb(i,6) = num2cell(tmp_gene_partic_choke+tmp_gene_partic_pagerank);
    end
end

% add the type of interaction
drug_wise_CTDdb(:,3) = drug_wise_CTDdb(:,2);
for i = 1:size(drug_wise_CTDdb,1)
    for j = 1:size(drug_wise_CTDdb{i,3},1)
        drug_wise_CTDdb{i,3}(j,:) = cellstr('NB');
    end
end
sto_inter = [];
for i = 1:size(drug_wise_CTDdb,1)
    tmp_drug = drug_wise_CTDdb{i,1};
    for j = 1:size(drug_wise_CTDdb{i,2},1)
       tmp_gene = drug_wise_CTDdb{i,2}(j);
       for k = 1:size(genes_list_CTD_clean_nnrd,1) 
           tmp_drug_nnr = genes_list_CTD_clean_nnrd{k,10};
           tmp_gene_nnr = genes_list_CTD_clean_nnrd{k,2};
           if strcmp(tmp_drug,tmp_drug_nnr) == 1 && strcmp(tmp_gene, tmp_gene_nnr) == 1
               drug_wise_CTDdb{i,3}(j,:) = cellstr(genes_list_CTD_clean_nnrd{k,14});
               sto_inter = [sto_inter; cellstr(genes_list_CTD_clean_nnrd{k,11})];
           end
       end
    end
end

% USE MATLABR2017b onwards!
% In column 14 there's the PMID from the DGIdb. Now we have a list of
% unique drugs (drug_wise), where the genes they interact with are in the
% 2nd column. We add a 3rd column for the automated results in pubmed of
% searching 'glioblastoma + DRUG', and a 4th column for 'cancer + DRUG'.
% Entries in #3 validate the capacity of identifying known/studied drug
% targets and drugs; those not appearing in #3 but appearing in #4 are of
% interest.
% only consider publications from the last 10 years
% if the combination yields nothing (because there are no papers for the
% drug) results will be just for GBM. Therefore, also get results from the
% drug by itself, and discard those for which no result was retrieved

% if the server tries to kick me out, pause the query for a minute and
% resume
 pause('on');
 n = 60

drug_wise_CTDdb(:,7) = {[]};
drug_wise_CTDdb(:,8) = {[]};
drug_wise_CTDdb(:,9) = {[]};

secondTerm = ' '
for i = 1:size(drug_wise_CTDdb,1)
    try
    searchterm = drug_wise_CTDdb{i,1};
    [PMID] = getpubmed(searchterm, secondTerm);
    drug_wise_CTDdb{i,7} = PMID;  
    clear searchterm PMID
    catch
        pause(n)
        i = i-1
    end
end
clear i

secondTerm = 'glioblastoma'
for i = 1:size(drug_wise_CTDdb,1)
    try
    searchterm = drug_wise_CTDdb{i,1};
    [PMID] = getpubmed(searchterm, secondTerm);
    drug_wise_CTDdb{i,8} = PMID;  
    clear searchterm PMID
    catch
        pause(n)
        i = i-1
    end
end
clear i

secondTerm = 'cancer'
for i = 1:size(drug_wise_CTDdb,1)
    try
    searchterm = drug_wise_CTDdb{i,1};
    [PMID] = getpubmed(searchterm, secondTerm);
    drug_wise_CTDdb{i,9} = PMID;  
    clear searchterm PMID
    catch
        pause(n)
        i = i-1
    end
end

% since the CTDdb also includes pollutants, some hits will refer to
% chemical causing cancer, rather than therapeutics. Therefore, we add an
% extra search to filter out later on
secondTerm = 'therapy'
for i = 1:size(drug_wise_CTDdb,1)
    try
    searchterm = drug_wise_CTDdb{i,1};
    [PMID] = getpubmed(searchterm, secondTerm);
    drug_wise_CTDdb{i,10} = PMID;  
    clear searchterm PMID
    catch
        pause(n)
        i = i-1
    end
end

% identify drugs for which there were results for the drug itself and no PMIDs were recovered for 'GBM + drug', but
% there were results for 'cancer + drug', the point being it is more likely
% for a drug known to work for another cancer to work for GBM. In
% particular, only considercases where at least 9 papers where found (when
% it says 12, there may be more, since only takes the first page)
for i = 1:size(drug_wise_CTDdb,1)
    tmp_papers = size(drug_wise_CTDdb{i,9},2);
    if isempty(drug_wise_CTDdb{i,7}) == 0 && isempty(drug_wise_CTDdb{i,8}) == 1 && tmp_papers >= 0 && isempty(drug_wise_CTDdb{i,10}) == 0
        drug_wise_CTDdb(i,11) = num2cell(tmp_papers);
    end
end
clear tmp_papers
for i = 1:size(drug_wise_CTDdb,1)
    tmp_papers = size(drug_wise_CTDdb{i,10},2);
    if isempty(drug_wise_CTDdb{i,7}) == 0 && isempty(drug_wise_CTDdb{i,8}) == 1 && tmp_papers >= 0 && isempty(drug_wise_CTDdb{i,10}) == 0
        drug_wise_CTDdb(i,12) = num2cell(tmp_papers);
    end
end
clear tmp_papers
%
drug_wise_interest_pre_CTD = []
for i = 1:size(drug_wise_CTDdb,1)
    if drug_wise_CTDdb{i,11} < 9 % trick it so it goes through
    elseif drug_wise_CTDdb{i,12} >= 9 % a decent amount of hits for 'drug + cancer' and 'drug + therapy'
        drug_wise_interest_pre_CTD = [drug_wise_interest_pre_CTD; drug_wise_CTDdb(i,:)];
    end
end
%
% create a cell per PMID so we can manage it manually later
entries = num2cell(zeros(size(drug_wise_interest_pre_CTD,1),10));
drug_wise_interest_pre_CTD = [drug_wise_interest_pre_CTD entries];
for i = 1:size(drug_wise_interest_pre_CTD,1)
    tmp_num = size(drug_wise_interest_pre_CTD{i,9},2);
    for j = 1:tmp_num
        drug_wise_interest_pre_CTD{i,12+j} = drug_wise_interest_pre_CTD{i,9}{j};    
    end
end

% flag whenever something interesting appears
drug_wise_interest_pre_CTD(:,25) = {[]};
tmp_term_1 = 'decreases^reaction';
tmp_term_2 = 'decreases^activity';
tmp_term_3 = 'decreases^expression';
tmp_term_4 = 'affects^folding';
tmp_term_5 = 'affects^metabolic processing';

for i = 1:size(drug_wise_interest_pre_CTD,1)
    for j = 1:size(drug_wise_interest_pre_CTD{i,3},1)
        tmp_inter = drug_wise_interest_pre_CTD{i,3}{j,:};
        tmp_inter_term_1 = strfind(tmp_inter,tmp_term_1);
        tmp_inter_term_2 = strfind(tmp_inter,tmp_term_2);
        tmp_inter_term_3 = strfind(tmp_inter,tmp_term_3);
        tmp_inter_term_4 = strfind(tmp_inter,tmp_term_4);
        tmp_inter_term_5 = strfind(tmp_inter,tmp_term_5);
        
        if isempty(tmp_inter_term_1) == 0 || isempty(tmp_inter_term_2) == 0 || isempty(tmp_inter_term_3) == 0 || isempty(tmp_inter_term_4) == 0 || isempty(tmp_inter_term_5) == 0
            drug_wise_interest_pre_CTD(i,25) = num2cell(1);
        end
        clear tmp_inter_term_1 tmp_inter_term_2 tmp_inter_term_3 tmp_inter_term_4 tmp_inter_term_5
    end
end

drug_wise_interest_CTD = num2cell(zeros(size(drug_wise_interest_pre_CTD,1),1));
drug_wise_interest_CTD(:,1) = ({'CTD'});
drug_wise_interest_CTD = [drug_wise_interest_CTD drug_wise_interest_pre_CTD];

% drug data: Only data curated from articles by BindingDB
BD_inhibition = tabularTextDatastore('BindingDB_BindingDB_Inhibition_2020_11_01.txt');
BD_inhibition.SelectedVariableNames = {'BindingDBReactant_set_id','BindingDBLigandName','TargetNameAssignedByCuratorOrDataSource','PMID','ChEMBLIDOfLigand','DrugBankIDOfLigand','UniProt_SwissProt_EntryNameOfTargetChain','UniProt_SwissProt_PrimaryIDOfTargetChain'};
BD_inhibition.TextscanFormats(32) = {'%s'}; % need to change the format of CHEMbl and DrugBank
BD_inhibition.TextscanFormats(33) = {'%s'};
BD_inhibition_subset = table2cell(readall(BD_inhibition));
% only keep human genes
BD_inhibition_subset_human = []
tmp_human = 'HUMAN'
for i = 1:size(BD_inhibition_subset,1)
    tmp_gene = BD_inhibition_subset{i,7};
    j = strfind(tmp_gene,tmp_human);
    if isempty(j) == 0
        BD_inhibition_subset_human = [BD_inhibition_subset_human; BD_inhibition_subset(i,:)];
    end
    clear j
end

% import known Uniprot IDs (manually done, for full_genes_un). For some of
% the genes there are more than 1 entry
BD_inhibition_subset_human_uniprot = BD_inhibition_subset_human(:,8);
dictionary_gene_unique_uniprotID = table2cell(readtable('gene_uniprot_ID.txt'));
uniprot_ids_gene_unique = (dictionary_gene_unique_uniprotID(:,2));
[common_uniprot_ids,ia,ib] = intersect(BD_inhibition_subset_human_uniprot,uniprot_ids_gene_unique,'stable')
BD_inhibition_subset_human_interest = []
for j = 1:size(common_uniprot_ids,1)
    tmp_uni_j = common_uniprot_ids{j,1};
    for i=1:size(BD_inhibition_subset_human,1)
        tmp_uni_i = BD_inhibition_subset_human{i,8};
        if strcmp(tmp_uni_i,tmp_uni_j) == 1
            BD_inhibition_subset_human_interest = [BD_inhibition_subset_human_interest; BD_inhibition_subset_human(i,:)];
        end
    end
end
clear tmp_uni_j tmp_uni_i

% structure data similarly to previous cases
for i = 1:size(BD_inhibition_subset_human_interest,1)
    tmp_uni = BD_inhibition_subset_human_interest{i,8};
    for j = 1:size(dictionary_gene_unique_uniprotID,1)
        tmp_uni_dicc = dictionary_gene_unique_uniprotID{j,2};
        if strcmp(tmp_uni,tmp_uni_dicc) == 1
            BD_inhibition_subset_human_interest{i,9} = dictionary_gene_unique_uniprotID{j,1};            
        end
    end
end
genes_list_BD_inhibition_pre(:,1) = num2cell((1:size(BD_inhibition_subset_human_interest,1))');
genes_list_BD_inhibition_pre(:,2) = (BD_inhibition_subset_human_interest(:,9));
genes_list_BD_inhibition_pre(:,10) = (BD_inhibition_subset_human_interest(:,2)); % drug

% modify '::' byr ' OR ' so it the drug names can be used in a pubmed
% search
for i = 1:size(genes_list_BD_inhibition_pre,1)
    tmp_drug = genes_list_BD_inhibition_pre{i,10};
    tmp_drug_new = strrep(tmp_drug,'::',' OR ');
    genes_list_BD_inhibition_pre{i,10} = tmp_drug_new;
    clear tmp_drug tmp_drug_new
end

genes_list_BD_inhibition_pre(:,11) = (BD_inhibition_subset_human_interest(:,5)); % CHEMBL
genes_list_BD_inhibition_pre(:,12) = (BD_inhibition_subset_human_interest(:,6)); % Drugbank
genes_list_BD_inhibition_pre(:,13) = ({[]}); % 
genes_list_BD_inhibition_pre(:,14) = ({[]}); % 
genes_list_BD_inhibition_pre(:,15) = ({[]}); % 
for i = 1:size(genes_list_BD_inhibition_pre,1)
    tmp_gene = genes_list_BD_inhibition_pre{i,2};
    for j = 1:size(full_genes_un,1)
        tmp_gene_un = full_genes_un{j,1};
        if strcmp(tmp_gene, tmp_gene_un) == 1
            genes_list_BD_inhibition_pre(i,3) = full_genes_un(j,2);
            genes_list_BD_inhibition_pre(i,4) = full_genes_un(j,3);
            genes_list_BD_inhibition_pre(i,5) = full_genes_un(j,4);
            genes_list_BD_inhibition_pre(i,6) = full_genes_un(j,5);
            genes_list_BD_inhibition_pre(i,7) = full_genes_un(j,6);
            genes_list_BD_inhibition_pre(i,8) = full_genes_un(j,7);
            genes_list_BD_inhibition_pre(i,9) = full_genes_un(j,8);
        end
    end
end

% Focus just on those defined in the DrugBank
genes_list_BD_inhibition_clean_nnrd = [];
for i = 1:size(genes_list_BD_inhibition_pre,1)
    tmp_dbank = genes_list_BD_inhibition_pre{i,12};
    if isempty(tmp_dbank) == 0
        genes_list_BD_inhibition_clean_nnrd = [genes_list_BD_inhibition_clean_nnrd; genes_list_BD_inhibition_pre(i,:)];
    end
end

% get a drug-wise view, considering an accumulated gene participation (which was previously normalised to %)
[drug_wise_BD_inhibition, ~] = unique(genes_list_BD_inhibition_clean_nnrd(:,10));
drug_wise_BD_inhibition(:,2) = {[]};
drug_wise_BD_inhibition(:,3) = {'NB'};
drug_wise_BD_inhibition(:,4) = {[]};
drug_wise_BD_inhibition(:,5) = {[]};
drug_wise_BD_inhibition(:,6) = {[]};
if isempty(drug_wise_BD_inhibition{1,1}) == 1
    drug_wise_BD_inhibition = drug_wise_BD_inhibition(2:end,:);
end
for i = 1:size(drug_wise_BD_inhibition,1)
    tmp_gene = [];
    tmp_gene_partic_choke = 0;
    tmp_gene_partic_pagerank = 0;
    tmp_drug = drug_wise_BD_inhibition{i,1};
    for j = 1:size(genes_list_BD_inhibition_clean_nnrd)
        tmp_drug_j = genes_list_BD_inhibition_clean_nnrd{j,10};
        if strcmp(tmp_drug,tmp_drug_j) == 1
           tmp_gene = [tmp_gene; genes_list_BD_inhibition_clean_nnrd(j,2)];
           tmp_gene_partic_choke = tmp_gene_partic_choke + cell2mat(genes_list_BD_inhibition_clean_nnrd(j,5));
           tmp_gene_partic_pagerank = tmp_gene_partic_pagerank + cell2mat(genes_list_BD_inhibition_clean_nnrd(j,6));
        end
        drug_wise_BD_inhibition{i,2} = tmp_gene;

        drug_wise_BD_inhibition(i,4) = num2cell(tmp_gene_partic_choke);
        drug_wise_BD_inhibition(i,5) = num2cell(tmp_gene_partic_pagerank);
        drug_wise_BD_inhibition(i,6) = num2cell(tmp_gene_partic_choke+tmp_gene_partic_pagerank);
    end
end

% no interaction info was available for this DB (we can assume they have an
% inhibitory effect, as explained in the website)

% USE MATLABR2017b onwards!
% In column 14 there's the PMID from the DGIdb. Now we have a list of
% unique drugs (drug_wise), where the genes they interact with are in the
% 2nd column. We add a 3rd column for the automated results in pubmed of
% searching 'glioblastoma + DRUG', and a 4th column for 'cancer + DRUG'.
% Entries in #3 validate the capacity of identifying known/studied drug
% targets and drugs; those not appearing in #3 but appearing in #4 are of
% interest.
% only consider publications from the last 10 years
% if the combination yields nothing (because there are no papers for the
% drug) results will be just for GBM. Therefore, also get results from the
% drug by itself, and discard those for which no result was retrieved

% if the server tries to kick me out, pause the query for a minute and
% resume
 pause('on');
 n = 60

drug_wise_BD_inhibition(:,7) = {[]};
drug_wise_BD_inhibition(:,8) = {[]};
drug_wise_BD_inhibition(:,9) = {[]};

secondTerm = ' '
for i = 1:size(drug_wise_BD_inhibition,1)
    try
    searchterm = drug_wise_BD_inhibition{i,1};
    [PMID] = getpubmed(searchterm, secondTerm);
    drug_wise_BD_inhibition{i,7} = PMID;  
    clear searchterm PMID
    catch
        pause(n)
        i = i-1
    end
end
clear i

secondTerm = 'glioblastoma'
for i = 1:size(drug_wise_BD_inhibition,1)
    try
    searchterm = drug_wise_BD_inhibition{i,1};
    [PMID] = getpubmed(searchterm, secondTerm);
    drug_wise_BD_inhibition{i,8} = PMID;  
    clear searchterm PMID
    catch
        pause(n)
        i = i-1
    end
end
clear i

secondTerm = 'cancer'
for i = 1:size(drug_wise_BD_inhibition,1)
    try
    searchterm = drug_wise_BD_inhibition{i,1};
    [PMID] = getpubmed(searchterm, secondTerm);
    drug_wise_BD_inhibition{i,9} = PMID;  
    clear searchterm PMID
    catch
        pause(n)
        i = i-1
    end
end
 
% we add an extra search to filter out later on
secondTerm = 'therapy'
for i = 1:size(drug_wise_BD_inhibition,1)
    try
    searchterm = drug_wise_BD_inhibition{i,1};
    [PMID] = getpubmed(searchterm, secondTerm);
    drug_wise_BD_inhibition{i,10} = PMID;  
    clear searchterm PMID
    catch
        pause(n)
        i = i-1
    end
end
 % manually check results from DrugBank -> names too complex for Pubmed to
 % handle properly

%% PanDrugs
% since most drugs target just one gene (and if not, all genes are in the
% same cell), this time make it simpler by manually checking the relevance
% at the end
PAN_data = table2cell(readtable('GBM_pandrugs_2020_12_03.txt'));
 
% USE MATLABR2017b onwards!
% if the server tries to kick me out, pause the query for a minute and
% resume
 pause('on');
 n = 60

PAN_data(:,13) = {[]};
PAN_data(:,14) = {[]};
PAN_data(:,15) = {[]};
PAN_data(:,16) = {[]};

secondTerm = ' '
for i = 1:size(PAN_data,1)
    try
    searchterm = PAN_data{i,2};
    [PMID] = getpubmed(searchterm, secondTerm);
    PAN_data{i,13} = PMID;  
    clear searchterm PMID
    catch
        pause(n)
        i = i-1
    end
end
clear i

secondTerm = 'glioblastoma'
for i = 1:size(PAN_data,1)
    try
    searchterm = PAN_data{i,2};
    [PMID] = getpubmed(searchterm, secondTerm);
    PAN_data{i,14} = PMID;  
    clear searchterm PMID
    catch
        pause(n)
        i = i-1
    end
end
clear i

secondTerm = 'cancer'
for i = 1:size(PAN_data,1)
    try
    searchterm = PAN_data{i,2};
    [PMID] = getpubmed(searchterm, secondTerm);
    PAN_data{i,15} = PMID;  
    clear searchterm PMID
    catch
        pause(n)
        i = i-1
    end
end

% identify drugs for which there were results for the drug itself and no PMIDs were recovered for 'GBM + drug', but
% there were results for 'cancer + drug', the point being it is more likely
% for a drug known to work for another cancer to work for GBM. In
% particular, only considercases where at least 9 papers where found (when
% it says 12, there may be more, since only takes the first page)
for i = 1:size(PAN_data,1)
    tmp_papers = size(PAN_data{i,15},2);
    if isempty(PAN_data{i,13}) == 0 && isempty(PAN_data{i,14}) == 1 && tmp_papers >= 0 && isempty(PAN_data{i,15}) == 0
        PAN_data(i,16) = num2cell(tmp_papers);
    end
end
clear tmp_papers

% % create a cell per PMID so we can manage it manually later
PAN_data_interest = []
for i =1:size(PAN_data,1)
    tmp_papers = (PAN_data{i,16});
    if tmp_papers >= 9 
        PAN_data_interest = [PAN_data_interest; PAN_data(i,:)];
    end
end

for i = 1:size(PAN_data_interest,1)
    tmp_num = size(PAN_data_interest{i,15},2);
    for j = 1:tmp_num
        PAN_data_interest{i,16+j} = PAN_data_interest{i,15}{j};    
    end
end






































