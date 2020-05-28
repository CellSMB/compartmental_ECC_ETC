%% Parameter sweep for NFAT model
% H Hunt 2020
% Explore IP3 concentration over a range of ip3r parameter sets

function results=NFAT_PS(isip3)
if isip3
    IP3 = 0:0.2:1.5; 
    k_f = 15;
else
    IP3 = 1; 
    k_f = 0:10:100;
end
K_c = [0.5 1 3 5 10 20 30 40 50];
K_h = 6;
K_t=10;
t_max=10;
total_runs=(size(IP3,2)*size(k_f,2)*size(K_c,2));
results=zeros(total_runs,7);
parfor ii=1:total_runs
    IP3it=ceil(ii/(size(k_f,2)*size(K_c,2)));
    kfit=ceil((mod(ii-1,size(k_f,2)*size(K_c,2))+1)/size(K_c,2));
    Kcit=mod(mod(ii-1,size(k_f,2)*size(K_c,2)),size(K_c,2))+1;
    model_params = [IP3(IP3it),K_t,K_c(Kcit),K_h,t_max,k_f(kfit),1,1,1,1,3e3];
    [~, STATES, ~, ~,~,~] = NFAT_Cooling2009(model_params);
    results(ii,:)=[IP3(IP3it),K_c(Kcit),k_f(kfit),STATES(end,:)];
end 
end