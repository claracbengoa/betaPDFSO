function main

delete('*.txt')
delete('*.mat')
delete('*.eps')
delete('*.tif')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The beta-PDFSO problem is a new approach for fail-safe optimization of 
% structures, which takes into account two sources of uncertainty: 

%   - The inherent uncertainty affecting the structural response, 
%     characterized through random variables
%   - The available information on the probability of occurrence of 
%     different accidental scenarios. 

% The results lead to a less conservative and more appropriate design 
% compared to the traditional fail-safe RBDO optimization, as the actual 
% data of each accidental situation are included in the optimization 
% process.

% Input data must be modified in the function: defineInputParameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[d0,optPrm,modelPrm,limStatePrm,damConfPrm,rndPrm] = defineInputParameters;

auxiliaryParameters(modelPrm)

%--------------------------------------------------------------------------
options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter','Diagnostics','on',...
            'DiffMinChange',optPrm.DiffMinChange,...
            'DiffMaxChange',optPrm.DiffMaxChange,...
            'TolCon',optPrm.TolCon,'TolFun',optPrm.TolFun,...
            'OutputFcn',@outfun);

[d,fval,exitflag,output]=fmincon(@(d)Objfun(d),d0,[],[],[],[],optPrm.lb,optPrm.ub,...
                                 @(d)confun(d,modelPrm,optPrm,limStatePrm,damConfPrm,rndPrm),options)
%--------------------------------------------------------------------------
plotKey=1;
save('plotKey.mat','plotKey')

fval=Objfun(d);
[c,ceq] = confun(d,modelPrm,optPrm,limStatePrm,damConfPrm,rndPrm);
save('c.mat','c')

counter = 0;
nameFolder='opt_results_it';
writeOutputResults(nameFolder,exitflag,output)
writeActiveConstr(modelPrm,optPrm,nameFolder,counter,fval,d,c)
writeOvercomeLimitStateAndConfig(damConfPrm,nameFolder)

end

function [d0,optPrm,modelPrm,limStatePrm,damConfPrm,rndPrm] = defineInputParameters

% initial value of design variables:
% d0 = [1 1];
d0 = [2.5402 2.2277];

% number of design variables
nDV = length(d0);

% number of damaged configurations:
nDamages = 60;

% probability of occurrence of each damaged configuration
pDamages = [ones(1,20)'*0.025;ones(1,20)'*0.015;ones(1,20)'*0.010]; % must have nDamages components
% pDamages = ones(1,nDamages)'/nDamages; % all damaged configurations have the same probability

if abs(sum(pDamages)-1) > 1e-8 % two numbers equal within a tolerance.
   error('Error. The sum of probabilities must be equal to 1')
elseif length(pDamages) ~= nDamages
    error('Error. Length of pDamages must be equal to nDamages')
end

% target probability of failure for PDFSO approach:
pfPDFSO = 0.10;

%--------------------------------------------------------------------------
% random variables:
% x1
mu_x1 = 0.25;
sig_x1 = 0.1*mu_x1;
% x2
mu_x2 = 1.0;
sig_x2 = 0.1*mu_x2;

% target probability of failure for Reliability Analysis:
pfRA=1e-4;
beta_min=-norminv(pfRA,0,1);

% beta >= beta_min    G safety region: beta - beta_min >= 0     G: beta/beta_min - 1 >= 0
betaFun = @(betaValue,beta_min) betaValue/beta_min - 1;

%------------------------------------------

% number of Load Cases
nLC = 2;

% number of Limit-states:
nDcon = 2;

% lists of Load Cases and Limit-states
vLC =   [];
vDcon = [];
for h=1:nLC
    for k=1:nDcon
        vLC =   [vLC,h];
        vDcon = [vDcon,k];         
    end
end

% In this example, Coefficients for damaged configurations are randomly generated
rng(1);

%--------------------------------------------------------------------------
% Add here the limit-states and Load Cases to be considered:
% h = Load Case
% k = Limit-state
% response <= responseMax    G safety region: responseMax - response >= 0     G:1 - response/responseMax
h=1;
k=1;
resp{h,1}{k,1} = @(d,Cf,x1,x2) ((1+x1)/(Cf(1) + Cf(4)*d(1)*d(2) + Cf(5)*d(1)^2 + Cf(6)*d(2)^2 + x2));
respMax{h,1}{k,1} = 1/4;
GFun{h,1}{k,1} = @(resp,respMax) 1 - resp/respMax;
Cfd0{h,1}{k,1} = ones(1,6);
CfDamages{h,1}{k,1} = lhsdesign(nDamages,6);
dGFun_dresp{h,1}{k,1} = @(respMax)-1/respMax;
dresp_dx1{h,1}{k,1} = @(d,Cf,x1,x2) (1/(Cf(1) + Cf(4)*d(1)*d(2) + Cf(5)*d(1)^2 + Cf(6)*d(2)^2 + x2));
dresp_dx2{h,1}{k,1} = @(d,Cf,x1,x2) (-(1+x1)/(Cf(1) + Cf(4)*d(1)*d(2) + Cf(5)*d(1)^2 + Cf(6)*d(2)^2 + x2)^2);

h=1;
k=2;
resp{h,1}{k,1} = @(d,Cf,x1,x2) ((1+x1)/(Cf(1) + Cf(2)*d(1) + Cf(3)*d(2) + Cf(5)*d(1)^2 + Cf(6)*d(2)^2 + x2));
respMax{h,1}{k,1} = 1/4;
GFun{h,1}{k,1} = @(resp,respMax) 1 - resp/respMax;
Cfd0{h,1}{k,1} = ones(1,6);
CfDamages{h,1}{k,1} = lhsdesign(nDamages,6);
dGFun_dresp{h,1}{k,1} = @(respMax)-1/respMax;
dresp_dx1{h,1}{k,1} = @(d,Cf,x1,x2) (1/(Cf(1) + Cf(2)*d(1) + Cf(3)*d(2) + Cf(5)*d(1)^2 + Cf(6)*d(2)^2 + x2));
dresp_dx2{h,1}{k,1} = @(d,Cf,x1,x2) (-(1+x1)/(Cf(1) + Cf(2)*d(1) + Cf(3)*d(2) + Cf(5)*d(1)^2 + Cf(6)*d(2)^2 + x2)^2);

h=2;
k=1;
resp{h,1}{k,1} = @(d,Cf,x1,x2) ((1+x1)/(Cf(1) + Cf(2)*d(1) + Cf(3)*d(2) + Cf(4)*d(1)*d(2) + Cf(6)*d(2)^2 + x2));
respMax{h,1}{k,1} = 1/4;
GFun{h,1}{k,1} = @(resp,respMax) 1 - resp/respMax;
Cfd0{h,1}{k,1} = ones(1,6);
CfDamages{h,1}{k,1} = CfDamages{1,1}{1,1};
dGFun_dresp{h,1}{k,1} = @(respMax)-1/respMax;
dresp_dx1{h,1}{k,1} = @(d,Cf,x1,x2) (1/(Cf(1) + Cf(2)*d(1) + Cf(3)*d(2) + Cf(4)*d(1)*d(2) + Cf(6)*d(2)^2 + x2));
dresp_dx2{h,1}{k,1} = @(d,Cf,x1,x2) (-(1+x1)/(Cf(1) + Cf(2)*d(1) + Cf(3)*d(2) + Cf(4)*d(1)*d(2) + Cf(6)*d(2)^2 + x2)^2);

h=2;
k=2;
resp{h,1}{k,1} = @(d,Cf,x1,x2) ((1+x1)/(Cf(1) + Cf(2)*d(1) + Cf(3)*d(2) + Cf(4)*d(1)*d(2) + Cf(5)*d(1)^2 + x2));
respMax{h,1}{k,1} = 1/4;
GFun{h,1}{k,1} = @(resp,respMax) 1 - resp/respMax;
Cfd0{h,1}{k,1} = ones(1,6);
CfDamages{h,1}{k,1} = CfDamages{1,1}{2,1};
dGFun_dresp{h,1}{k,1} = @(respMax)-1/respMax;
dresp_dx1{h,1}{k,1} = @(d,Cf,x1,x2) (1/(Cf(1) + Cf(2)*d(1) + Cf(3)*d(2) + Cf(4)*d(1)*d(2) + Cf(5)*d(1)^2 + x2));
dresp_dx2{h,1}{k,1} = @(d,Cf,x1,x2) (-(1+x1)/(Cf(1) + Cf(2)*d(1) + Cf(3)*d(2) + Cf(4)*d(1)*d(2) + Cf(5)*d(1)^2 + x2)^2);

%--------------------------------------------------------------------------
lb=[1e-03,1e-03];   % lower bound values
ub=[]; % upper bound values

TolCon=1e-3; % constraint tolerance
TolFun=1e-3; % objective function tolerance
DiffMinChange=0.01; % minimum change in design variables for finite difference gradients
DiffMaxChange=0.02; % maximum change in design variables for finite difference gradients
%--------------------------------------------------------------------------

% admissible violations in probability constraints of damaged structure and beta to construct the CDF
admViolpf = 0.008; % admito 0.058
scaleFactorProbConstr = TolCon / ((pfPDFSO+admViolpf)/pfPDFSO - 1);

beta_min_viol_CDF = beta_min - 0.0190; % 3.70 instead of 3.7190
%--------------------------------------------------------------------------

optPrm = struct('lb',lb,'ub',ub,'TolCon',TolCon,'TolFun',TolFun,'DiffMinChange',DiffMinChange,'DiffMaxChange',DiffMaxChange);
modelPrm = struct('nDcon',nDcon,'nLC',nLC,'nDV',nDV,'vLC',vLC,'vDcon',vDcon);
limStatePrm = struct('resp',{resp},'respMax',{respMax},'GFun',{GFun},'Cfd0',{Cfd0},'CfDamages',{CfDamages},...
                     'dGFun_dresp',{dGFun_dresp},'dresp_dx1',{dresp_dx1},'dresp_dx2',{dresp_dx2});
damConfPrm = struct('nDamages',nDamages,'pDamages',pDamages,'pfPDFSO',pfPDFSO,...
    'scaleFactorProbConstr',scaleFactorProbConstr,'beta_min_viol_CDF',beta_min_viol_CDF);
rndPrm = struct('mu_x1',mu_x1,'sig_x1',sig_x1,...
                'mu_x2',mu_x2,'sig_x2',sig_x2,'betaFun',betaFun,'beta_min',beta_min);

save('modelPrm.mat','modelPrm')
end

function auxiliaryParameters(modelPrm)

% plot the CDFs at each iteration 0=no 1=yes
plotKey=0;
save('plotKey.mat','plotKey')

% finite difference gradients counter
n_CONT = 0;  
save('n_CONT.mat','n_CONT')

% hitory of design variables and objective function
history_d_opt = [];
history_fval_opt = [];
save('history_d_opt.mat','history_d_opt');
save('history_fval_opt.mat','history_fval_opt');
%--------------------------------------------------------------------------
% History of Results
counter = 'iter';
nameFolder= 'opt_results_it';
writeTitleObjFunAndDesVar(modelPrm.nDV,nameFolder,counter)

Results_history.designVariables = [];
Results_history.c = [];
Results_history.f = [];
Results_history.betaList = [];
Results_history.probList = [];

Results_history.RA_history.RA_evolution_d0 = {};
Results_history.RA_history.RA_evolution_danos = {};

Results_history.beta_converged_d0 = {};
Results_history.beta_converged_danos = {};

save('Results_history.mat','Results_history')

end

function fval=Objfun(d)

fval=d(1)+d(2);
save('fval.mat','fval')

end

function [c,ceq]=confun(d,modelPrm,optPrm,limStatePrm,damConfPrm,rndPrm)

load('plotKey.mat')
extension='.sh'; 

load('n_CONT.mat')
n_CONT = n_CONT + 1;  
save('n_CONT.mat','n_CONT')

ceq = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reliability_analysis(d,modelPrm,damConfPrm,limStatePrm,rndPrm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('RA_evolution_d0.mat')
load('RA_evolution_danos.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[beta_converged_d0,beta_converged_danos] = readBetaValues(modelPrm,damConfPrm,...
    RA_evolution_d0,RA_evolution_danos);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%----------------------------
% INTACT MODEL
%----------------------------
[constr_d0] = evaluateLimitStateIntactModel(modelPrm,rndPrm,beta_converged_d0);

[c] = evaluateConstraintsIntactModel(modelPrm,constr_d0);

%----------------------------
% DAMAGED CONFIGURATIONS
%----------------------------
[constrDamagedConf]=evaluateLimitStateDamagedConf(modelPrm,damConfPrm,rndPrm,...
    beta_converged_danos);

[x_CDFs,f_CDFs] = builtCDFs(modelPrm,damConfPrm,constrDamagedConf);
save('x_CDFs.mat','x_CDFs')
save('f_CDFs.mat','f_CDFs')

[p_CDFs] = obtainProbabilityFromCDFs(modelPrm,x_CDFs,f_CDFs,plotKey);

[c] = evaluateConstraintsDamagedConf(modelPrm,damConfPrm,c,p_CDFs);

%--------------------------------------------------------------------------
load('Results_history.mat')
load('beta_current_it.mat')
load('probabilities_current_it.mat')
load('fval.mat')
%---------------------------
Results_history.beta_converged_d0{n_CONT,1} = beta_converged_d0;
Results_history.beta_converged_danos{n_CONT,1} = beta_converged_danos;

Results_history.designVariables = [Results_history.designVariables;d];
Results_history.c = [Results_history.c;c];
Results_history.f = [Results_history.f;fval];
Results_history.betaList = [Results_history.betaList;beta_current_it];
Results_history.probList = [Results_history.probList;probabilities_current_it];
save('Results_history.mat','Results_history')

end

function stop = outfun(d,optimValues,state)

    stop = false;

     switch state
         case 'init'
%              hold on
         case 'iter'
            load('modelPrm.mat')
            load('history_d_opt.mat');
            load('history_fval_opt.mat');
            history_fval_opt = [history_fval_opt; optimValues.fval];
            history_d_opt = [history_d_opt; d];
            save('history_d_opt.mat','history_d_opt');
            save('history_fval_opt.mat','history_fval_opt');
            i=length(history_fval_opt);
            
            nameFolder='opt_results_it';
            writeObjFunAndDesVar(modelPrm.nDV,nameFolder,history_d_opt(i,:),history_fval_opt(i),i)
  
         case 'done'
%              hold off
         otherwise
     end
end

%--------------------------------------------------------------------------

function reliability_analysis(d,modelPrm,damConfPrm,limStatePrm,rndPrm)

[NumLimitStates] = setDataHistoryRA(modelPrm,damConfPrm,rndPrm);

%-----------------------------------
n_TOL=0; % numero de estados limite en los que beta converge 
TOL=5e-4;
save('TOL.mat','TOL') 

n_CONT_FORM = 0;
%-----------------------------------

% NumLimitStates --> se calcula en funcion de los datos de numero de danos y limit-states en cada dano
while n_TOL ~= NumLimitStates
    
    n_CONT_FORM = n_CONT_FORM + 1;
    save('n_CONT_FORM.mat','n_CONT_FORM') 
    
    %--------------------
    % MUY IMPORTANTE!!!!!
    %--------------------
    % Cargamos los ficheros .mat de la iteracion anterior dentro de cada
    % subrutina    
    %--------------------
   
    % Calculo x_new de cada estado limite   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    calculo_iteracion_FORM(d,modelPrm,damConfPrm,limStatePrm,rndPrm);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [n_TOL] = compruebo_convergencia_FORM(modelPrm,damConfPrm);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
end % end of while loop

end

function [NumLimitStates] = setDataHistoryRA(modelPrm,damConfPrm,rndPrm)

% INICIALMENTE NINGUN RA CONVERGE
Binary_converg_Matrix_d0 = {};
Binary_converg_Matrix_danos = {};

x1_itera_Matrix_d0 = {};
x2_itera_Matrix_d0 = {};
Estiff_itera_Matrix_d0 = {};

x1_itera_Matrix_danos = {};
x2_itera_Matrix_danos = {};
Estiff_itera_Matrix_danos = {};

RA_evolution_d0 = {};
RA_evolution_danos = {};

NumLimitStates = 0;
for g=0:damConfPrm.nDamages
    for iv=1:length(modelPrm.vLC)
        h=modelPrm.vLC(iv);
        k=modelPrm.vDcon(iv);
        hh=int2str(h);
        titulo2=strcat('LC',hh);
        kk=int2str(k);
        titulo3=strcat('LS',kk);
        if g == 0
            Binary_converg_Matrix_d0{1,1}.(titulo2).(titulo3) = 2;
            x1_itera_Matrix_d0{1,1}.(titulo2).(titulo3) = rndPrm.mu_x1;
            x2_itera_Matrix_d0{1,1}.(titulo2).(titulo3) = rndPrm.mu_x2;
            Gx_anterior_Matrix_d0{1,1}.(titulo2).(titulo3) = 1;
            RA_evolution_d0{1,1}.(titulo2).(titulo3) = struct('Gx',[],'beta',[],'x',[],'u',[],'c',[]);
        else
            Binary_converg_Matrix_danos{g,1}.(titulo2).(titulo3) = 2;
            x1_itera_Matrix_danos{g,1}.(titulo2).(titulo3) = rndPrm.mu_x1;
            x2_itera_Matrix_danos{g,1}.(titulo2).(titulo3) = rndPrm.mu_x2;
            Gx_anterior_Matrix_danos{g,1}.(titulo2).(titulo3) = 1;
            RA_evolution_danos{g,1}.(titulo2).(titulo3) = struct('Gx',[],'beta',[],'x',[],'u',[],'c',[]);
        end

        NumLimitStates = NumLimitStates + 1;
    end 
end

save('NumLimitStates.mat','NumLimitStates')

save('Binary_converg_Matrix_d0.mat','Binary_converg_Matrix_d0')
save('x1_itera_Matrix_d0.mat','x1_itera_Matrix_d0')
save('x2_itera_Matrix_d0.mat','x2_itera_Matrix_d0')
save('Gx_anterior_Matrix_d0.mat','Gx_anterior_Matrix_d0')
save('RA_evolution_d0.mat','RA_evolution_d0')

save('Binary_converg_Matrix_danos.mat','Binary_converg_Matrix_danos')
save('x1_itera_Matrix_danos.mat','x1_itera_Matrix_danos')
save('x2_itera_Matrix_danos.mat','x2_itera_Matrix_danos')
save('Gx_anterior_Matrix_danos.mat','Gx_anterior_Matrix_danos')
save('RA_evolution_danos.mat','RA_evolution_danos')

end

function calculo_iteracion_FORM(d,modelPrm,damConfPrm,limStatePrm,rndPrm)

load('Binary_converg_Matrix_d0.mat')
load('Binary_converg_Matrix_danos.mat')

load('x1_itera_Matrix_d0.mat')
load('x2_itera_Matrix_d0.mat')

load('x1_itera_Matrix_danos.mat')
load('x2_itera_Matrix_danos.mat')

load('RA_evolution_d0.mat')
load('RA_evolution_danos.mat')

Gx_nuevo_Matrix_d0 = {};
Gx_nuevo_Matrix_danos = {};

beta_Matrix_d0 = {};
beta_Matrix_danos = {};

mu_x = [rndPrm.mu_x1; rndPrm.mu_x2];
sig_x = [rndPrm.sig_x1; rndPrm.sig_x2];

%---------------------------------------------------------

for g=0:damConfPrm.nDamages
    for iv=1:length(modelPrm.vLC)
        h=modelPrm.vLC(iv);
        k=modelPrm.vDcon(iv);
        hh=int2str(h);
        titulo2=strcat('LC',hh);
        kk=int2str(k);
        titulo3=strcat('LS',kk);
            
        if g == 0
            BCM_g = Binary_converg_Matrix_d0{1,1}.(titulo2).(titulo3);
        else
            BCM_g = Binary_converg_Matrix_danos{g,1}.(titulo2).(titulo3);
        end

        if BCM_g == 1 % if the limit-state has converged
            ;
        else % an additional iteration is needed fot this limit-state

            dG_dresponse = limStatePrm.dGFun_dresp{h,1}{k,1}(limStatePrm.respMax{h,1}{k,1});
            
            if g == 0  
                x1_itera = x1_itera_Matrix_d0{1,1}.(titulo2).(titulo3);
                x2_itera = x2_itera_Matrix_d0{1,1}.(titulo2).(titulo3);
                
                resp = limStatePrm.resp{h,1}{k,1}(d,limStatePrm.Cfd0{h,1}{k,1},x1_itera,x2_itera);
                Gx = limStatePrm.GFun{h,1}{k,1}(resp,limStatePrm.respMax{h,1}{k,1});
                
                dresponse_dx1 = limStatePrm.dresp_dx1{h,1}{k,1}(d,limStatePrm.Cfd0{h,1}{k,1},x1_itera,x2_itera);
                dresponse_dx2 = limStatePrm.dresp_dx2{h,1}{k,1}(d,limStatePrm.Cfd0{h,1}{k,1},x1_itera,x2_itera);
                
            else
                x1_itera = x1_itera_Matrix_danos{g,1}.(titulo2).(titulo3);
                x2_itera = x2_itera_Matrix_danos{g,1}.(titulo2).(titulo3);
                
                resp = limStatePrm.resp{h,1}{k,1}(d,limStatePrm.CfDamages{h,1}{k,1}(g,:),x1_itera,x2_itera);
                Gx = limStatePrm.GFun{h,1}{k,1}(resp,limStatePrm.respMax{h,1}{k,1});
 
                dresponse_dx1 = limStatePrm.dresp_dx1{h,1}{k,1}(d,limStatePrm.CfDamages{h,1}{k,1},x1_itera,x2_itera);
                dresponse_dx2 = limStatePrm.dresp_dx2{h,1}{k,1}(d,limStatePrm.CfDamages{h,1}{k,1},x1_itera,x2_itera);
                
            end

            dG_dx1 = dG_dresponse * dresponse_dx1;
            dG_dx2 = dG_dresponse * dresponse_dx2;

            dG_dx = [dG_dx1; dG_dx2];
            dx_du = sig_x;
            dG_du = dG_dx .* dx_du;

            x_itera = [x1_itera; x2_itera]; 

            u = (x_itera - mu_x)./sig_x; 

            s1 = u' * dG_du;
            s2 = norm(dG_du);

            %Pf = func.distr.acum(-beta) -->  beta = - func.distr.acum.inv(Pf)
            beta = (Gx-s1)/s2;

            n1=(-dG_du(1))/s2;
            n2=(-dG_du(2))/s2;

            n=[n1;n2]; 

            %----------------------------------------------------------     
            u_anterior = u;

            u_nuevo = n*beta;

            c=3;
            u_new = u_anterior + (u_nuevo - u_anterior)/c;   
            %----------------------------------------------------------

            x_new = u_new .* sig_x + mu_x;

            load('n_CONT.mat')
            load('n_CONT_FORM.mat')
            
            if n_CONT_FORM==40
                if g == 0
                    Binary_converg_Matrix_d0{1,1}.(titulo2).(titulo3)=1;
                    save('Binary_converg_Matrix_d0.mat','Binary_converg_Matrix_d0')
                else
                    Binary_converg_Matrix_danos{g,1}.(titulo2).(titulo3)=1;
                    save('Binary_converg_Matrix_danos.mat','Binary_converg_Matrix_danos')
                end 
            end

            if g == 0
                Gx_nuevo_Matrix_d0{1,1}.(titulo2).(titulo3) = Gx;
                beta_Matrix_d0{1,1}.(titulo2).(titulo3) = beta;
                x1_itera_Matrix_d0{1,1}.(titulo2).(titulo3) = x_new(1);
                x2_itera_Matrix_d0{1,1}.(titulo2).(titulo3) = x_new(2);
                %----------------------
                RA_evolution_d0{1,1}.(titulo2).(titulo3).('Gx') = [RA_evolution_d0{1,1}.(titulo2).(titulo3).('Gx');Gx];
                RA_evolution_d0{1,1}.(titulo2).(titulo3).('beta') = [RA_evolution_d0{1,1}.(titulo2).(titulo3).('beta');beta];
                RA_evolution_d0{1,1}.(titulo2).(titulo3).('x') = [RA_evolution_d0{1,1}.(titulo2).(titulo3).('x');x_new'];
                RA_evolution_d0{1,1}.(titulo2).(titulo3).('u') = [RA_evolution_d0{1,1}.(titulo2).(titulo3).('u');u_new'];
                RA_evolution_d0{1,1}.(titulo2).(titulo3).('c') = [RA_evolution_d0{1,1}.(titulo2).(titulo3).('c');c];
            else
                Gx_nuevo_Matrix_danos{g,1}.(titulo2).(titulo3) = Gx;
                beta_Matrix_danos{g,1}.(titulo2).(titulo3) = beta;
                x1_itera_Matrix_danos{g,1}.(titulo2).(titulo3) = x_new(1);
                x2_itera_Matrix_danos{g,1}.(titulo2).(titulo3) = x_new(2);
                %----------------------
                RA_evolution_danos{g,1}.(titulo2).(titulo3).('Gx') = [RA_evolution_danos{g,1}.(titulo2).(titulo3).('Gx');Gx];
                RA_evolution_danos{g,1}.(titulo2).(titulo3).('beta') = [RA_evolution_danos{g,1}.(titulo2).(titulo3).('beta');beta];
                RA_evolution_danos{g,1}.(titulo2).(titulo3).('x') = [RA_evolution_danos{g,1}.(titulo2).(titulo3).('x');x_new'];
                RA_evolution_danos{g,1}.(titulo2).(titulo3).('u') = [RA_evolution_danos{g,1}.(titulo2).(titulo3).('u');u_new'];
                RA_evolution_danos{g,1}.(titulo2).(titulo3).('c') = [RA_evolution_danos{g,1}.(titulo2).(titulo3).('c');c];
            end 

        end            

    end
end

save('Gx_nuevo_Matrix_d0.mat','Gx_nuevo_Matrix_d0')
save('Gx_nuevo_Matrix_danos.mat','Gx_nuevo_Matrix_danos')

save('beta_Matrix_d0.mat','beta_Matrix_d0')
save('beta_Matrix_danos.mat','beta_Matrix_danos')

save('x1_itera_Matrix_d0.mat','x1_itera_Matrix_d0')
save('x2_itera_Matrix_d0.mat','x2_itera_Matrix_d0')

save('x1_itera_Matrix_danos.mat','x1_itera_Matrix_danos')
save('x2_itera_Matrix_danos.mat','x2_itera_Matrix_danos')

save('RA_evolution_d0.mat','RA_evolution_d0')
save('RA_evolution_danos.mat','RA_evolution_danos')

end

function [n_TOL] = compruebo_convergencia_FORM(modelPrm,damConfPrm)

load('n_CONT_FORM.mat') 
load('TOL.mat') 

load('Gx_anterior_Matrix_d0.mat')
load('Gx_anterior_Matrix_danos.mat')

load('Gx_nuevo_Matrix_d0.mat')
load('Gx_nuevo_Matrix_danos.mat')

load('Binary_converg_Matrix_d0.mat')
load('Binary_converg_Matrix_danos.mat')

n_TOL=0; % number of limit-states where Gx converges

for g=0:damConfPrm.nDamages
    gg=int2str(g);
    titulo1=strcat('d',gg);
    for iv=1:length(modelPrm.vLC)
        h=modelPrm.vLC(iv);
        k=modelPrm.vDcon(iv);
        hh=int2str(h);
        titulo2=strcat('LC',hh);
        kk=int2str(k);
        titulo3=strcat('LS',kk);
            
        if g == 0
            BCM_g = Binary_converg_Matrix_d0{1,1}.(titulo2).(titulo3);
        else
            BCM_g = Binary_converg_Matrix_danos{g,1}.(titulo2).(titulo3);
        end

        if BCM_g == 1 % if the limit-state has converged
            n_TOL=n_TOL+1;
        else % an additional iteration is needed fot this limit-state
            if g == 0
                Gx_nuevo    = Gx_nuevo_Matrix_d0{1,1}.(titulo2).(titulo3);
                Gx_anterior = Gx_anterior_Matrix_d0{1,1}.(titulo2).(titulo3);
            else
                Gx_nuevo    = Gx_nuevo_Matrix_danos{g,1}.(titulo2).(titulo3);
                Gx_anterior = Gx_anterior_Matrix_danos{g,1}.(titulo2).(titulo3);
            end 

            dif_abs = abs( Gx_nuevo - Gx_anterior );

            if (n_CONT_FORM>1) && (dif_abs <= TOL) && (abs(Gx_nuevo) < 1e-3)
                n_TOL = n_TOL+1;
                if g == 0
                    Binary_converg_Matrix_d0{1,1}.(titulo2).(titulo3) = 1; % indico en la matriz si el RA de ese estado limite converge
                else
                    Binary_converg_Matrix_danos{g,1}.(titulo2).(titulo3) = 1; % indico en la matriz si el RA de ese estado limite converge
                end
            end 


        end 
    end
end

save('Binary_converg_Matrix_d0.mat','Binary_converg_Matrix_d0')
save('Binary_converg_Matrix_danos.mat','Binary_converg_Matrix_danos')

Gx_anterior_Matrix_d0 = Gx_nuevo_Matrix_d0;
Gx_anterior_Matrix_danos = Gx_nuevo_Matrix_danos;
save('Gx_anterior_Matrix_d0.mat','Gx_anterior_Matrix_d0')
save('Gx_anterior_Matrix_danos.mat','Gx_anterior_Matrix_danos')

end

function [beta_converged_d0,beta_converged_danos] = readBetaValues(modelPrm,damConfPrm,...
    RA_evolution_d0,RA_evolution_danos)

beta_converged_d0 = {};
beta_converged_danos = {};

for g=0:damConfPrm.nDamages
    for iv=1:length(modelPrm.vLC)
        h=modelPrm.vLC(iv);
        k=modelPrm.vDcon(iv);
        hh=int2str(h);
        titulo2=strcat('LC',hh);
        kk=int2str(k);
        titulo3=strcat('LS',kk);
            
        if g == 0
            beta_converged_d0{1,1}.(titulo2).(titulo3) = RA_evolution_d0{1,1}.(titulo2).(titulo3).('beta')(end);
        else
            beta_converged_danos{g,1}.(titulo2).(titulo3) = RA_evolution_danos{g,1}.(titulo2).(titulo3).('beta')(end);
        end
    end
end

end

%--------------------------------------------------------------------------

function [constr_d0] = evaluateLimitStateIntactModel(modelPrm,rndPrm,beta_converged_d0)

beta_current_it = [];

for iv=1:length(modelPrm.vLC)
    h=modelPrm.vLC(iv);
    k=modelPrm.vDcon(iv);
    hh=int2str(h);
    titulo2=strcat('LC',hh);
    kk=int2str(k);
    titulo3=strcat('LS',kk);
    betaValue = beta_converged_d0{1,1}.(titulo2).(titulo3);
    beta_current_it = [beta_current_it,betaValue];
    constr_d0_h_k = rndPrm.betaFun(betaValue,rndPrm.beta_min);
    constr_d0{h,1}{k,1} = constr_d0_h_k;
end

save('beta_current_it.mat','beta_current_it')
save('constr_d0.mat','constr_d0')

end

function [c] = evaluateConstraintsIntactModel(modelPrm,constr_d0)

l=0;

for iv=1:length(modelPrm.vLC)
    h=modelPrm.vLC(iv);
    k=modelPrm.vDcon(iv);
    l=l+1;
    % constr>=0 --> c=-constr
    c(l) = -constr_d0{h,1}{k,1};
end

end

function [constrDamagedConf]=evaluateLimitStateDamagedConf(modelPrm,damConfPrm,rndPrm,beta_converged_danos)

for iv=1:length(modelPrm.vLC)
    h=modelPrm.vLC(iv);
    k=modelPrm.vDcon(iv);
    hh=int2str(h);
    titulo2=strcat('LC',hh);
    kk=int2str(k);
    titulo3=strcat('LS',kk);
    
    for g=1:damConfPrm.nDamages
        betaValue(g,:) = beta_converged_danos{g,1}.(titulo2).(titulo3); % leo en ese dano el beta asociado a ese LS
        if betaValue(g,:)<=rndPrm.beta_min && betaValue(g,:)>damConfPrm.beta_min_viol_CDF
            betaValue(g,:) = 1.01 * rndPrm.beta_min;
        end
    end
    
    constrDamagedConf_h_k = rndPrm.betaFun(betaValue,rndPrm.beta_min);
    constrDamagedConf{h,1}{k,1} = constrDamagedConf_h_k;
    betaDamagedConf{h,1}{k,1} = betaValue;

end

save('constrDamagedConf.mat','constrDamagedConf')
save('betaDamagedConf.mat','betaDamagedConf')

end

function [x_CDFs,f_CDFs] = builtCDFs(modelPrm,damConfPrm,constrDamagedConf)

for iv=1:length(modelPrm.vLC)
    h=modelPrm.vLC(iv);
    k=modelPrm.vDcon(iv);
    constrDamagedConf_h_k = constrDamagedConf{h,1}{k,1};
    [x_prov,index] = sort(constrDamagedConf_h_k,'ascend');

    f_prov = damConfPrm.pDamages(index);

    x_limitState_k = unique(x_prov); % sorted and non-repeated values of x_prov
    f_PDF_limitState_k = zeros(1,length(x_limitState_k))';
    for i=1:length(x_limitState_k)
         f_PDF_limitState_k(i) = sum(f_prov(x_limitState_k(i) == x_prov));
    end 

    f_CDF_limitState_k = cumsum(f_PDF_limitState_k);

    if abs(x_limitState_k(1))<=5e-4 % constr negtive but active. Modify the value to 1e-8
        x_limitState_k(1)=1e-8;
    end

    if length(x_limitState_k) == 1 % The CDF can not have only 1 value
        x_value=x_limitState_k;
        epsilon = 1e-8;
        x_limitState_k = [x_value-epsilon;x_value;x_value+epsilon];    
        f_CDF_limitState_k = [0;0.5;1];
    end     

    x_CDFs{h,1}{k,1} = x_limitState_k;
    f_CDFs{h,1}{k,1} = f_CDF_limitState_k;
end

end

function [p_CDFs] = obtainProbabilityFromCDFs(modelPrm,x_CDFs,f_CDFs,plotKey)

% Entro en CDF y calculo el valor de probabilidad
xq = 0;

for iv=1:length(modelPrm.vLC)
    h=modelPrm.vLC(iv);
    k=modelPrm.vDcon(iv);
    x_limitState_k = x_CDFs{h,1}{k,1};
    f_CDF_limitState_k = f_CDFs{h,1}{k,1};
    [fq] = obtainInterpolatedValue(x_limitState_k,f_CDF_limitState_k,xq);
    p_CDFs{h,1}{k,1} = fq;

    if plotKey ==1
        plot_CDF(k,h,xq,fq,x_limitState_k,f_CDF_limitState_k)
    end
end

end

function [fq] = obtainInterpolatedValue(x_limitState_k,f_CDF_limitState_k,xq)

fq = interp1(x_limitState_k,f_CDF_limitState_k,xq,'pchip'); % Previous neighbor interpolation
if xq < min(x_limitState_k)
    fq = 0 - min(x_limitState_k);
elseif xq >= max(x_limitState_k)
    fq = 1;
end

end

function [c] = evaluateConstraintsDamagedConf(modelPrm,damConfPrm,c,p_CDFs)

l=length(c);
probabilities_current_it = [];

for iv=1:length(modelPrm.vLC)
    h=modelPrm.vLC(iv);
    k=modelPrm.vDcon(iv);
    l=l+1;
    p = p_CDFs{h,1}{k,1};
    probabilities_current_it = [probabilities_current_it,p];
    % p<=pfPDFSO
    % scale the constraint to allow a violation of admViolpf using TolCon=1e-3 
    c(l) = (p/damConfPrm.pfPDFSO - 1)*damConfPrm.scaleFactorProbConstr;
    save('probabilities_current_it.mat','probabilities_current_it')
end

end

%--------------------------------------------------------------------------

function plot_CDF(LS,LC,xq,fq,x_limitState_k,f_CDF_limitState_k)

scrsz = get(0,'ScreenSize');   
h=figure('Position',[1,1,scrsz(3)*6/6,scrsz(4)*3/3]);
FIGTemp=gcf;
FigName=strcat('CDF_LC',num2str(LC),'_LS',num2str(LS),'_ex1');

p_piecewise = stairs(x_limitState_k,f_CDF_limitState_k,'LineStyle','-','color','b','LineWidth',0.5);
hold on

for j_points=1:length(x_limitState_k)
    p1=plot(x_limitState_k(j_points),f_CDF_limitState_k(j_points),'o','MarkerSize',2,'MarkerEdgeColor','r','MarkerFaceColor','r');
end

pinit = min(x_limitState_k);
pfin  = max(x_limitState_k);

pvector = linspace(pinit,pfin,1000);

[p_pchip] = plot_approximated_CDF(pvector,x_limitState_k,f_CDF_limitState_k);

% extend the axes limits
lb_factor = (1 - sign(pinit)*0.05);
ub_factor = (1 + sign(pfin)*0.05);

location_x = 0.33;
location_y = 0.70;

p2=plot([xq, xq],[0, 1],'LineStyle','--','color',[0.5,0.5,0.5],'LineWidth',0.3);
p3=plot([lb_factor*pinit, ub_factor*pfin],[fq, fq],'LineStyle','--','color',[0.5,0.5,0.5],'LineWidth',0.3);
p4=plot(xq,fq,'o','MarkerSize',2,'MarkerEdgeColor','g','MarkerFaceColor','g');

if fq<0 % truco para no poner probabilidad negativa en la leyenda
    fq=0;
end

[hleg1,hobj1]=legend([p_piecewise p1 p_pchip p4],'Piecewise CDF','Points from CDF','Interpolated CDF',strcat('$F_{\hat{H}_{',num2str(LC),',',num2str(LS),'}}(0)=P[\hat{H}_{',num2str(LC),',',num2str(LS),'} \le 0] =',num2str(fq),'$'),...
    'Location','east','Orientation','vertical', 'Interpreter', 'latex', 'fontsize', 10);

set(hleg1, 'Position', [location_x location_y 0.1 0.1]);

axis([lb_factor*pinit, ub_factor*pfin,0,1])
xlabel(strcat('$\hat{H}_{',num2str(LC),',',num2str(LS),'}$'),'FontSize',10,'Interpreter','latex')
ylabel(strcat('$F_{\hat{H}_{',num2str(LC),',',num2str(LS),'}}$'),'FontSize',10,'Interpreter','latex')
set(gca,'FontSize',10)
ax = gca;
ax.XAxis.MinorTick = 'off';
ax.YAxis.MinorTick = 'off';
ax.XAxis.MinorTickValues = lb_factor*pinit:(ub_factor*pfin-lb_factor*pinit)/8:ub_factor*pfin;
ax.YAxis.MinorTickValues = 0:0.1:1;

% % yaxis horizontal
% ylh = get(gca,'ylabel');
% ylp = get(ylh, 'Position');
% set(ylh, 'Rotation',0, 'Position',ylp-0.04, 'VerticalAlignment','middle', 'HorizontalAlignment','right')

% Paper position (to save files)
a1 = 20;
a2 = 12.2;
PaperPosition = [0 0 a1 a2];

set(FIGTemp,'PaperUnits','centimeters','PaperPosition',PaperPosition)

FileFormat1 = 'epsc';
FileFormat2 = 'tiff';
% FileFormat3 = 'png';
print(FIGTemp,strcat('-d',FileFormat1),'-r300', FigName);
print(FIGTemp,strcat('-d',FileFormat2),'-r300', FigName);
% print(FIGTemp,strcat('-d',FileFormat3),'-r300', FigName);
% saveas(FIGTemp,FigName,'fig');

close(FIGTemp);

end

function [p_pchip] = plot_approximated_CDF(pvector,x_limitState_k,f_CDF_limitState_k)

fqvector=[];
for j_point = pvector

    xq = j_point;
    fq = interp1(x_limitState_k,f_CDF_limitState_k,xq,'pchip'); % Previous neighbor interpolation
      
    fqvector=[fqvector,fq];
    
%     p1=plot(j_point,fq,'o','MarkerSize',3,'MarkerEdgeColor','b','MarkerFaceColor','none');

end

p_pchip = plot(pvector,fqvector,'LineStyle','-','color',[0.5,0.5,0.5],'LineWidth',0.5);

end

%--------------------------------------------------------------------------

function writeTitleObjFunAndDesVar(nDV,nameFolder,counter)

fid=fopen(strcat(nameFolder,'.txt'),'a'); % write design variables and fval for each iteration
    fprintf(fid,'%9s %12s', ...
                counter,'Obj. Fun');
    for k=1:nDV
        fprintf(fid,' %10s', ...
                     strcat('d',num2str(k)));
    end
    fprintf(fid,'\n');
fclose(fid);

end

function writeObjFunAndDesVar(nDV,nameFolder,d,fval,n_CONT) 

fid=fopen(strcat(nameFolder,'.txt'),'a');
    fprintf(fid,'%9d %12.4f', n_CONT, fval);
    for k=1:nDV
        fprintf(fid,' %10.4f', d(k));
    end 
    fprintf(fid,'\n');
fclose(fid);

end

function writeOutputResults(nameFolder,exitflag,output)

fid=fopen(strcat(nameFolder,'.txt'),'a');
    fprintf(fid,' \n')
fclose(fid);

fid=fopen(strcat(nameFolder,'.txt'),'a');
fprintf(fid,' \n %10s %16s \n', ...
                 'exitflag','constrviolation');
fprintf(fid,' %10i %16.4e \n', ...
             exitflag,output.constrviolation);             
fclose(fid);

end

function writeActiveConstr(modelPrm,optPrm,nameFolder,counter,fval,d,c)

load('n_CONT.mat')
load('Results_history.mat')

if counter ~= 0
    writeTitleObjFunAndDesVar(modelPrm.nDV,nameFolder,counter)
    writeObjFunAndDesVar(modelPrm.nDV,nameFolder,d,fval,n_CONT)
end
%--------------------------------------------------------------------------

matrixLimStateActViol = [];

% ACTIVE CONSTRAINTS
l=0;
fid=fopen(strcat(strcat(nameFolder,'.txt')),'a');
    fprintf(fid,'\n');
    
    fprintf(fid,'%30s \n','----- ACTIVE CONSTRAINTS -----');
    ll=0;
    g=0;
    for iv=1:length(modelPrm.vLC)
        h=modelPrm.vLC(iv);
        k=modelPrm.vDcon(iv);
        l=l+1;
        ll=ll+1;
        if abs(c(l)) <= optPrm.TolCon
            beta = Results_history.betaList(n_CONT,ll);
            fprintf(fid,'%24s', strcat('d0_LC',num2str(h),'_LS',num2str(k)));
            fprintf(fid,'%10s %12.4e %36s %12.4e', strcat('c(',num2str(l),')='),c(l),strcat('beta_d0_LC',num2str(h),'_LS',num2str(k),' (l=',num2str(l),')'),beta);
            fprintf(fid,'\n');
        end  
    end
    
    ll=0;
    % Damages
    for iv=1:length(modelPrm.vLC)
        h=modelPrm.vLC(iv);
        k=modelPrm.vDcon(iv);
        l=l+1;
        ll=ll+1;
        if abs(c(l)) <= optPrm.TolCon
            matrixLimStateActViol = [matrixLimStateActViol; [h,k]];
            prob = Results_history.probList(n_CONT,ll);
            fprintf(fid,'%24s', strcat('danos_LC',num2str(h),'_LS',num2str(k)));
            fprintf(fid,'%10s %12.4e %36s %12.4e', strcat('c(',num2str(l),')='),c(l),strcat('prob_failure_LC',num2str(h),'_LS',num2str(k),' (l=',num2str(l),')'),prob);
            fprintf(fid,'\n');
        end    
    end   
    
    fprintf(fid,'\n\n');
fclose(fid);

% VIOLATED CONSTRAINTS
l=0;
fid=fopen(strcat(strcat(nameFolder,'.txt')),'a');
    fprintf(fid,'\n');
    fprintf(fid,'%30s \n','----- VIOLATED CONSTRAINTS -----');
    ll=0;
    g=0;
    for iv=1:length(modelPrm.vLC)
        h=modelPrm.vLC(iv);
        k=modelPrm.vDcon(iv);
        l=l+1;
        ll=ll+1;
        if c(l) > optPrm.TolCon
            beta = Results_history.betaList(n_CONT,ll);
            fprintf(fid,'%24s', strcat('d0_LC',num2str(h),'_LS',num2str(k)));
            fprintf(fid,'%10s %12.4e %36s %12.4e', strcat('c(',num2str(l),')='),c(l),strcat('beta_d0_LC',num2str(h),'_LS',num2str(k),' (l=',num2str(l),')'),beta);
            fprintf(fid,'\n');
        end    
    end
    
    ll=0;
    % Damages
    for iv=1:length(modelPrm.vLC)
        h=modelPrm.vLC(iv);
        k=modelPrm.vDcon(iv);
        l=l+1;
        ll=ll+1;
        if c(l) > optPrm.TolCon
            matrixLimStateActViol = [matrixLimStateActViol; [h,k]];
            prob = Results_history.probList(n_CONT,ll);
            fprintf(fid,'%24s', strcat('danos_LC',num2str(h),'_LS',num2str(k)));
            fprintf(fid,'%10s %12.4e %36s %12.4e', strcat('c(',num2str(l),')='),c(l),strcat('prob_failure_LC',num2str(h),'_LS',num2str(k),' (l=',num2str(l),')'),prob);
            fprintf(fid,'\n');
        end    
    end  
    
    fprintf(fid,'\n\n');
fclose(fid);

save('matrixLimStateActViol.mat','matrixLimStateActViol')

end

function writeOvercomeLimitStateAndConfig(damConfPrm,nameFolder)

load('constrDamagedConf.mat')
load('betaDamagedConf.mat')
load('matrixLimStateActViol.mat')

[row,~]=size(matrixLimStateActViol);

fid=fopen(strcat(nameFolder,'.txt'),'a');
    fprintf(fid,' \n%69s \n', '----- DAMAGED CONFIGURATIONS WHERE THE LIMIT-STATE IS OVERCOME ----- '); 
fclose(fid);

for i=1:row
    h = matrixLimStateActViol(i,1);
    k = matrixLimStateActViol(i,2);
    constrDamagedConf_h_k = constrDamagedConf{h,1}{k,1};
    beta = betaDamagedConf{h,1}{k,1};
    for g=1:damConfPrm.nDamages
        if constrDamagedConf_h_k(g)<0
            constrDamagedConf_viol = constrDamagedConf_h_k(g);
            beta_viol = beta(g);
            ID_damageConfig_viol = g;
            p_viol = damConfPrm.pDamages(g);
            fid=fopen(strcat(nameFolder,'.txt'),'a');
                fprintf(fid,'\n %20s %5s %35s %25s %25s %25s ',strcat('LS=',num2str(k)),strcat('LC=',num2str(h)),strcat('constraintValue=',num2str(constrDamagedConf_viol)),strcat('DamagedConfiguration=',num2str(ID_damageConfig_viol)),strcat('BetaValue=',num2str(beta_viol)),strcat('ProbabilityValue=',num2str(p_viol))); 
            fclose(fid);
        end
    end
end

fid=fopen(strcat(nameFolder,'.txt'),'a');
    fprintf(fid,'\n'); 
fclose(fid);

end

