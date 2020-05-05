% Run parameter sweep and generate informative graphs from IncDelayModel.m
% Also requires getQuals.m
% HHunt 2018
namedesc='_1hzkf_recolrwp_thresh';
if 0
    numparRuns = 1;
% Choose parameters
% fluxvar     = 5*rand(numparRuns,1);
% t_max       = 1+1e9*rand(numparRuns,1);
% K_c         = 15*rand(numparRuns,1);
% K_h         = 15*rand(numparRuns,1);
% K_t         = 15*rand(numparRuns,1);

fluxvar     = 3*ones(numparRuns,1);
t_max       = 1e1*rand(numparRuns,1);
K_c         = 3.5*ones(numparRuns,1);
K_h         = 25*ones(numparRuns,1);
K_t         = 10*rand(numparRuns,1);
numRuns=numparRuns;

totalRuns = numRuns;
% Make room for results
resultsD=zeros(totalRuns,16);
resultsnD=zeros(totalRuns,16);
savepath='/home/hhunt1/MATLAB/Graphs/';
status = mkdir('/home/hhunt1/MATLAB/Graphs');
for it=1:numRuns
    params=[fluxvar(it),t_max(it),K_c(it),K_h(it),K_t(it)]
    [VOI, STATES, ALGEBRAIC, CONSTANTS] = ...
        IncDelayModel([fluxvar(it),t_max(it),K_c(it),K_h(it),K_t(it)]);
    period=CONSTANTS(:,7);
    lastTransient=find(VOI>=(VOI(end)-period));
	resultsD(it,:)=getQuals(VOI, STATES, ALGEBRAIC, CONSTANTS,params,lastTransient);
    [VOI2, STATES2, ALGEBRAIC2, CONSTANTS2] = ...
        IncDelayModel([fluxvar(it),1,K_c(it),K_h(it),1]);
    lastTransient2=find(VOI2>=(VOI2(end)-period));
	resultsnD(it,:)=getQuals(VOI2, STATES2, ALGEBRAIC2, CONSTANTS2,...
        [fluxvar(it),1,K_c(it),K_h(it),1],lastTransient2);
    
%     time=VOI>(VOI(end)-CONSTANTS(:,7)-500);
%     time2=VOI2>(VOI2(end)-CONSTANTS2(:,7)-500);
%     figure('pos',[600,0,600,600])
%     numa = 1;
%     numd =3;
%     subplot(numd,numa,1)
%     hold on
%     plot(VOI2(time2)-(VOI2(end)-CONSTANTS2(:,7)-500),STATES2(time2,1))
%     plot(VOI(time)-(VOI(end)-CONSTANTS(:,7)-500),STATES(time,1))
%     plot(VOI(1:(CONSTANTS(:,7)))+500,STATES(1:(CONSTANTS(:,7)),1))
%     legend({'Ca_i+IP_3','Ca_i+IP_3 w delay','Ca_i'})
%     title(strcat(num2str(fluxvar(it)),'-',num2str(t_max(it)),'-',...
%         num2str(K_c(it)),'-',num2str(K_h(it)),'-',num2str(K_t(it))));
%     hold on
%     subplot(numd,numa,2)
%     hold on
%     plot(VOI2(time2)-(VOI2(end)-CONSTANTS2(:,7)-500),ALGEBRAIC2(time2,66))
%     plot(VOI(time)-(VOI(end)-CONSTANTS(:,7)-500),ALGEBRAIC(time,66))
%     legend('IP_3R','IP_3R w delay')
%     subplot(numd,numa,3)
%     hold on
%     plot(VOI2(time2)-(VOI2(end)-CONSTANTS2(:,7)-500),ALGEBRAIC2(time2,65))
%     plot(VOI(time)-(VOI(end)-CONSTANTS(:,7)-500),ALGEBRAIC(time,65))
%     legend('P_0','P_0 w delay')
%     saveas(gcf,strcat(savepath,num2str(fluxvar(it)),...
%         '-',num2str(t_max(it)),'-',num2str(K_c(it)),'-',num2str(K_h(it)),...
%         '-',num2str(K_t(it)),'.png'))
end
% Calculate base qualities
transient = find(VOI<(CONSTANTS(:,7)));
baseresults=getQuals(VOI2, STATES2, ALGEBRAIC2, CONSTANTS2,[0,1,1,1,1],transient);
% Save results
results = [resultsnD;resultsD;baseresults];

end
% save('/home/hhunt1/MATLAB/paramsweep_delay','results');
% save('/Users/hhunt1/Documents/Compartmental Model/MATLAB/paramsweep_green','results');
%% Generate scatter plots comparing variables and dependent quantities

newfold = strcat('Graphs',namedesc);
savegin = strcat('/Users/hhunt1/Documents/Compartmental Model/MATLAB/Postsubmission/',newfold,'/');
st = mkdir(strcat('/Users/hhunt1/Documents/Compartmental Model/MATLAB/Postsubmission/',newfold));
ti = {'Amplitude','FDHM','FD at 90% of max','Time to peak','Baseline','Duty cycle',... %5+6=11
    'AmplitudeF','FDHMF','FD at 90% of maxF','Time to peak F','Duty cycle F',... %5+11=16
    'IP_3R flux','RyR flux','Diastolic [Ca^{2+}]_i','Diastolic [Ca^{2+}]_{SR}'};%20

yl={'(mM)','(ms)','(ms)','(ms)','(mM)','','(mM)','(ms)','(ms)','(ms)','','','','mM','mM'};
%switch point
numVars = 5;
numRuns=(size(results,1)-1);
% xl={'IP_3','t_{max}','K_c','K_h','K_t'};
xl={'k_f','t_{max}','K_c','K_h','K_t'};
% xl={'IP_3','t_{max}','K_c','K_h','K_t','max(I_{IP_3R})','N_{IP_3R}','N_{RyR}','N_{SERCA}','N_{NCX}'};
%%
if 0
for qual = [(numVars+1):11,12,16]
pttype = '.';
figure
numa=floor(sqrt(numVars));
numd=ceil(sqrt(numVars));
tvar=[1,2,3:numVars];
for plotit=1:(numVars)
subplot(numd,numa,plotit)
hold on
scatter(results(1:(end-1)/2,tvar(plotit)),results(1:(end-1)/2,qual),pttype,'MarkerEdgeAlpha',0.2)
scatter(results(((end-1)/2+1):(end-1),tvar(plotit)),results(((end-1)/2+1):(end-1),qual),pttype,'MarkerEdgeAlpha',0.2)
refline(0,results(end,qual))
xlabel(xl{tvar(plotit)})
ylabel(yl{qual-numVars})
end
subplot(numd,numa,1)
title(ti{qual-numVars})
subplot(numd,numa,2)
legend({'no delay','delay'})
saveas(gcf,strcat(savegin,...
    ti{qual-numVars},namedesc,'.png'))
end
end

%% I
% numVars=5;
numRuns=(size(results,1)-1)/2;
[a,b]=sort(results((numRuns+1):(2*numRuns),numVars+1));
resAmpD=results(b+numRuns,:);
[a,b]=sort(results(1:numRuns,numVars+1));
resAmpND=results(b,:);
[a,b]=sort(results((numRuns+1):(2*numRuns),numVars+5));
resbaseD=results(b+numRuns,:);
[a,b]=sort(results(1:numRuns,numVars+5));
resbaseND=results(b,:);
[a,b]=sort(results((numRuns+1):(2*numRuns),numVars+2));
resdur50D=results(b+numRuns,:);
[a,b]=sort(results(1:numRuns,numVars+2));
resdur50ND=results(b,:);
da=1e-4;db=2e-5;dd=20;
diffper=40;
dlowbound=1-diffper/100;
dhibound=1+diffper/100;
lowa = find(resAmpD(:,numVars+1)<(results(end,numVars+1)*dlowbound),1,'last');
lowna = find(resAmpND(:,numVars+1)<(results(end,numVars+1)*dlowbound),1,'last');
higha = find(resAmpD(:,numVars+1)>(results(end,numVars+1)*dhibound),1,'first');
highna = find(resAmpND(:,numVars+1)>(results(end,numVars+1)*dhibound),1,'first');
lowd=find(resdur50D(:,numVars+2)<results(end,numVars+2)*dlowbound,1,'last');
lownd=find(resdur50ND(:,numVars+2)<results(end,numVars+2)*dlowbound,1,'last');
highd=find(resdur50D(:,numVars+2)>results(end,numVars+2)*dhibound,1,'first');
highnd=find(resdur50ND(:,numVars+2)>results(end,numVars+2)*dhibound,1,'first');
lowb=find(resbaseD(:,numVars+5)<results(end,numVars+5)*dlowbound,1,'last');
lownb=find(resbaseND(:,numVars+5)<results(end,numVars+5)*dlowbound,1,'last');
highb = find(resbaseD(:,numVars+5)>results(end,numVars+5)*dhibound,1,'first');
highnb = find(resbaseND(:,numVars+5)>results(end,numVars+5)*dhibound,1,'first');
%%
% cm=[54 171 147]/255;cl=[22 70 217]/255;ch=[242 208 0]/255;
% cl=0.4*ones(3,1);cm=0.6*ones(3,1);ch=0.8*ones(3,1);
if var(resbaseND(:,3))~=0 & var(resbaseND(:,4))~=0
cl=[59 109 224]/255;cm=0.6*ones(3,1);ch=[0.9 0 0];
figure('position',[0,200,300,200])
hold on
pttype = 'o';sz=30;
h1=scatter(resbaseND(:,3),resbaseND(:,4),sz,pttype,'filled','MarkerEdgeColor',cm,'MarkerFaceColor',cm);
h2=scatter(resbaseND(1:lownb,3),resbaseND(1:lownb,4),sz,pttype,'filled','MarkerEdgeColor',cl,'MarkerFaceColor',cl);
h3=scatter(resbaseND(highnb:end,3),resbaseND(highnb:end,4),sz,pttype,'filled','MarkerEdgeColor',ch,'MarkerFaceColor',ch);
pttype = 'o';sz=40;lw=1;
h4=scatter(resdur50ND(:,3),resdur50ND(:,4),sz,pttype,'MarkerEdgeColor',cm,'LineWidth',lw);
h5=scatter(resdur50ND(highnd:end,3),resdur50ND(highnd:end,4),sz,pttype,'MarkerEdgeColor',ch,'LineWidth',lw);
pttype='.'; 
h6=scatter(resAmpND(:,3),resAmpND(:,4),pttype,'MarkerEdgeColor',cm);
h7=scatter(resAmpND(1:lowna,3),resAmpND(1:lowna,4),pttype,'MarkerEdgeColor',cl);
h8=scatter(resAmpND(highna:end,3),resAmpND(highna:end,4),pttype,'MarkerEdgeColor',ch);
legend([h6,h1,h4],{'amplitude','base [Ca^{2+}]','duration'},'Position',[0.7 0.8 0.2 0.1])
% legend boxoff
title('No delay')
ylabel('K_h')
xlabel('K_c')
saveas(gcf,strcat(savegin,...
    'nd',num2str(diffper),namedesc,'.png'))
end
%%
figure
scatter3(resAmpD(:,4),resAmpD(:,3),resAmpD(:,6))%,ones(500,1)*30,resAmpD(:,1))
view(90,90)
xlabel('K_h')
ylabel('K_c')
colorbar
%%
if var(resbaseND(:,3))~=0 & var(resbaseND(:,4))~=0
cm=[54 171 147]/255;cl=[22 70 217]/255;ch=[242 208 0]/255;
cl=0.4*ones(3,1);cm=0.6*ones(3,1);ch=0.8*ones(3,1);
cl=[59 109 224]/255;cm=0.6*ones(3,1);ch=[0.9 0 0];
figure('position',[0,200,300,200])
var1=3;var2=4;
hold on
pttype = 'o';sz=30;
h1=scatter(resbaseD(:,var1),resbaseD(:,var2),sz,pttype,'filled','MarkerEdgeColor',cm,'MarkerFaceColor',cm);
h2=scatter(resbaseD(1:lowb,var1),resbaseD(1:lowb,var2),sz,pttype,'filled','MarkerEdgeColor',cl,'MarkerFaceColor',cl);
h3=scatter(resbaseD(highb:end,var1),resbaseD(highb:end,var2),sz,pttype,'filled','MarkerEdgeColor',ch,'MarkerFaceColor',ch);
pttype = 'o';sz=40;lw=1;
h4=scatter(resdur50D(:,var1),resdur50D(:,var2),sz,pttype,'MarkerEdgeColor',cm,'LineWidth',lw);
h5=scatter(resdur50D(highd:end,var1),resdur50D(highd:end,var2),sz,pttype,'MarkerEdgeColor',ch,'LineWidth',lw);
pttype='.'; lw=1.5;
h6=scatter(resAmpD(:,var1),resAmpD(:,var2),pttype,'MarkerEdgeColor',cm);%,'LineWidth',lw);
h7=scatter(resAmpD(1:lowa,var1),resAmpD(1:lowa,var2),pttype,'MarkerEdgeColor',cl);%,'LineWidth',lw);
h8=scatter(resAmpD(higha:end,var1),resAmpD(higha:end,var2),pttype,'MarkerEdgeColor',ch);%,'LineWidth',lw);
legend([h6,h1,h4],{'amplitude','base [Ca^{2+}]','duration'},'Position',[0.7 0.8 0.2 0.1])
% legend boxoff
title('Delay')
xlabel('K_c')
ylabel('K_h')
saveas(gcf,strcat(savegin,...
    'dth',num2str(diffper),namedesc,'.png'))
end
if 0
%% Voronoi
% numVars=5;
pl=0;
[a,b]=sort(results((1:500)+pl,numVars+1));
resAmp=results(b+pl,:);
[a,b]=sort(results((1:500)+pl,numVars+5));
resbase=results(b+pl,:);
[a,b]=sort(results((1:500)+pl,numVars+2));
resdur50=results(b+pl,:);

lowa = find(resAmp(:,numVars+1)<(results(end,numVars+1)-2e-5),1,'last');
higha = find(resAmp(:,numVars+1)>(results(end,numVars+1)+2e-5),1,'first');
lowd=find(resdur50(:,numVars+2)<results(end,numVars+2)-400,1,'last');
highd=find(resdur50(:,numVars+2)>results(end,numVars+2)+400,1,'first');
lowb=find(resbase(:,numVars+5)<results(end,numVars+5)-0.1e-4,1,'last');
highb = find(resbase(:,numVars+5)>results(end,numVars+5)+0.1e-4,1,'first');
hmm3=intersect(resAmp(higha:end,3),intersect(resbase((lowb+1):(highb-1),3),resdur50(1:(highd-1),3),'stable'),'stable');
mmm3=intersect(resAmp((lowa+1):(higha-1),3),intersect(resbase((lowb+1):(highb-1),3),resdur50(1:(highd-1),3),'stable'),'stable');
mlm3=intersect(resAmp((lowa+1):(higha-1),3),intersect(resbase(1:lowb,3),resdur50(1:(highd-1),3),'stable'),'stable');
lhh3=intersect(resAmp(1:lowa,3),intersect(resbase(highb:end,3),resdur50(highd:end,3),'stable'),'stable');
lhm3=intersect(resAmp(1:lowa,3),intersect(resbase(highb:end,3),resdur50(1:(highd-1),3),'stable'),'stable');
lmh3=intersect(resAmp(1:lowa,3),intersect(resbase((lowb+1):(highb-1),3),resdur50(highd:end,3),'stable'),'stable');
lmm3=intersect(resAmp(1:lowa,3),intersect(resbase((lowb+1):(highb-1),3),resdur50(1:(highd-1),3),'stable'),'stable');
llh3=intersect(resAmp(1:lowa,3),intersect(resbase(1:lowb,3),resdur50(highd:end,3),'stable'),'stable');
llm3=intersect(resAmp(1:lowa,3),intersect(resbase(1:lowb,3),resdur50(1:(highd-1),3),'stable'),'stable');
hlm3=intersect(resAmp(higha:end,3),intersect(resbase(1:lowb,3),resdur50(1:(highd-1),3),'stable'),'stable');
K_c={hmm3,mmm3,mlm3,lhh3,lhm3,lmh3,lmm3,llh3,llm3,hlm3};
hmm4=intersect(resAmp(higha:end,4),intersect(resbase((lowb+1):(highb-1),4),resdur50(1:(highd-1),4),'stable'),'stable');
mmm4=intersect(resAmp((lowa+1):(higha-1),4),intersect(resbase((lowb+1):(highb-1),4),resdur50(1:(highd-1),4),'stable'),'stable');
mlm4=intersect(resAmp((lowa+1):(higha-1),4),intersect(resbase(1:lowb,4),resdur50(1:(highd-1),4),'stable'),'stable');
lhh4=intersect(resAmp(1:lowa,4),intersect(resbase(highb:end,4),resdur50(highd:end,4),'stable'),'stable');
lhm4=intersect(resAmp(1:lowa,4),intersect(resbase(highb:end,4),resdur50(1:(highd-1),4),'stable'),'stable');
lmh4=intersect(resAmp(1:lowa,4),intersect(resbase((lowb+1):(highb-1),4),resdur50(highd:end,4),'stable'),'stable');
lmm4=intersect(resAmp(1:lowa,4),intersect(resbase((lowb+1):(highb-1),4),resdur50(1:(highd-1),4),'stable'),'stable');
llh4=intersect(resAmp(1:lowa,4),intersect(resbase(1:lowb,4),resdur50(highd:end,4),'stable'),'stable');
llm4=intersect(resAmp(1:lowa,4),intersect(resbase(1:lowb,4),resdur50(1:(highd-1),4),'stable'),'stable');
hlm4=intersect(resAmp(higha:end,4),intersect(resbase(1:lowb,4),resdur50(1:(highd-1),4),'stable'),'stable');
K_h={hmm4,mmm4,mlm4,lhh4,lhm4,lmh4,lmm4,llh4,llm4,hlm4};
% K_c = [mean(resbase(1:lowb,3));mean(resbase(highb:end,3));mean(resbase((lowb+1):(highb-1),3));...
%     mean(resdur50(1:(highd-1),3));mean(resdur50(highd:end,3));...
%     mean(resAmp(1:lowb,3));mean(resAmp(highb:end,3));mean(resAmp((lowb+1):(highb-1),3))];
% K_h = [mean(resbase(1:lowb,4));mean(resbase(highb:end,4));mean(resbase((lowb+1):(highb-1),4));...
%     mean(resdur50(1:(highb-1),4));mean(resdur50(highb:end,4));...
%     mean(resAmp(1:lowb,4));mean(resAmp(highb:end,4));mean(resAmp((lowb+1):(highb-1),4))];

legtext={'hmm','mmm','mlm','lhh','lhm','lmh','lmm','llh','llm','hlm'};
figure
hold on
for ii=1:9
x=K_c{ii};
y=K_h{ii};
k=convhull(x,y);
plot(x(k),y(k))
% scatter(x,y,'b')
end
legend(legtext)
% K_h=15+(K_c-6.75)^2;
end
if 0%%
prms=[3,1000,4.5,25,1;...
    3,1000,3,26,1;...
    3,1000,4.5,13,1;...
    3,1000,2.8,17.5,1;...
    3,1000,1,0.75,1;...
    3,1000,1,18,1;...
    3,1000,3,2.7,1;...
    3,1000,3,10,1];
desc = {'hmh','hlh','mmh','mlh','lhh','lmh','lmm','llh'};
    simr=struct();
for ii=1:8
    [VOI, STATES, ALGEBRAIC, CONSTANTS]=IncDelayModel(prms(ii,:));
    simr(ii).ALGEBRAIC=ALGEBRAIC;
    simr(ii).CONSTANTS=CONSTANTS;
    simr(ii).STATES=STATES;
    simr(ii).VOI=VOI;
    simr(ii).desc=desc{ii};
end
end
%%
t_max=1e8;
K_t=1;
prms=[3,t_max,4.5,25,1;...
    3,t_max,3,26,1;...
    3,t_max,4.5,13,1;...
    3,t_max,2.8,17.5,1;...
    3,t_max,1,0.75,1;...
    3,t_max,1,18,1;...
    3,t_max,3,2.7,1;...
    3,t_max,3,10,1];
desc = {'hmh','hlh','mmh','mlh','lhh','lmh','lmm','llh'};
    simr=struct();
for ii=1:8
    [VOI, STATES, ALGEBRAIC, CONSTANTS]=IncDelayModel(prms(ii,:));
    simr(ii).ALGEBRAIC=ALGEBRAIC;
    simr(ii).CONSTANTS=CONSTANTS;
    simr(ii).STATES=STATES;
    simr(ii).VOI=VOI;
    simr(ii).desc=desc{ii};
end
%%
figure
    numa = 3;
    numd =3;
for ii=(1:8)
    VOI = simr(ii).VOI;
    CONSTANTS=simr(ii).CONSTANTS;
    ALGEBRAIC=simr(ii).ALGEBRAIC;
    STATES=simr(ii).STATES;

    time=VOI>(VOI(end)-CONSTANTS(:,7)-500);
    subplot(numd,numa,ii)
    hold on
    plot(VOI(1:(CONSTANTS(:,7)))+500,STATES(1:(CONSTANTS(:,7)),1),'-')
    plot(VOI(time)-(VOI(end)-CONSTANTS(:,7)-500),STATES(time,1),'-')
    ylabel('[Ca^{2+}] (mM)')
    xlabel('time (ms)')
    title(desc{ii})
    axis([0 3000 0 1.2e-3])
end

legend({'ECC [Ca^{2+}]','IP_3Rs activated'},'Position',[0.7 0.17 0.2 0.1])
legend boxoff

%%
%1
subset=[1,3,5,6];
subset=[2,8];
% subset=[4,6,7];
% %2
% subset=[1,2];
% %?
% subset=[1,3,5];

ii=subset(1);
VOI = simr(ii).VOI;
CONSTANTS=simr(ii).CONSTANTS;
ALGEBRAIC=simr(ii).ALGEBRAIC;
STATES=simr(ii).STATES;

time=VOI>(VOI(end)-CONSTANTS(:,7)-500);
    figure('pos',[600,0,600,600])
    numa = 1;
    numd =3;
    subplot(numd,numa,1)
    hold on
    plot(VOI(1:(CONSTANTS(:,7)))+500,STATES(1:(CONSTANTS(:,7)),1))
    plot(VOI(time)-(VOI(end)-CONSTANTS(:,7)-500),STATES(time,1))
%     legend({'Hinch model','Hinch+IP_3'})
    ylabel('[Ca^{2+}] (mM)')
%     title(var3);
    hold on
    subplot(numd,numa,2)
    hold on
    plot(VOI(1:(CONSTANTS(:,7)))+500,ALGEBRAIC(1:(CONSTANTS(:,7)),66))
    plot(VOI(time)-(VOI(end)-CONSTANTS(:,7)-500),ALGEBRAIC(time,66))
    ylabel('I_{IP_3R} (mM/ms)')
    subplot(numd,numa,3)
    hold on
    plot(VOI(1:(CONSTANTS(:,7)))+500,ALGEBRAIC(1:(CONSTANTS(:,7)),35))
    plot(VOI(time)-(VOI(end)-CONSTANTS(:,7)-500),ALGEBRAIC(time,35))
    ylabel('I_{RyR} (mM/ms)')
    xlabel('time (ms)')

for ii=subset(2:end)
    VOI = simr(ii).VOI;
    CONSTANTS=simr(ii).CONSTANTS;
    ALGEBRAIC=simr(ii).ALGEBRAIC;
    STATES=simr(ii).STATES;

    time=VOI>(VOI(end)-CONSTANTS(:,7)-500);
    subplot(numd,numa,1)
    plot(VOI(time)-(VOI(end)-CONSTANTS(:,7)-500),STATES(time,1))
    subplot(numd,numa,2)
    plot(VOI(time)-(VOI(end)-CONSTANTS(:,7)-500),ALGEBRAIC(time,66))
    subplot(numd,numa,3)
    plot(VOI(time)-(VOI(end)-CONSTANTS(:,7)-500),ALGEBRAIC(time,35))
end
legend([{'Hinch'},desc{subset}])

%% Setup graphs
numa = 4;
numd = 4;
% xl={'IP_3R flux','t_{max}','K_c','K_h','K_t'};
% numVars=5;
tvar=1:numVars;
varVals=cell(1,4);
dc=results(((end-1)/2+1):(end-1),11);
for ii=1:4
    varVals{ii}=results(((end-1)/2+1):(end-1),tvar(ii));
end

%% Circle graphs
for qual = [(numVars+1):11,17,18]
    dc=results(((end-1)/2+1):(end-1),qual)/max(results(((end-1)/2+1):(end-1),qual))/2;
 if prod(dc)~=0   
figure
for ii=1:4
    for jj=[1:ii-1,ii+1:4]
        subplot(numd,numa,numa*(ii-1)+jj)
%         [X,Y] = ndgrid(linspace(min(varVals{ii}),max(varVals{ii}),150),linspace(min(varVals{jj}),max(varVals{jj}),150));
%         Z = griddata(varVals{ii},varVals{jj},dc,X,Y,'cubic');
%         contourf(X,Y,Z)
%         hold on
        scatter(varVals{ii},varVals{jj},100*dc.^2,dc,'o')
        xlabel(xl(ii))
        ylabel(xl(jj))
    end
end
for ii=1:4
    subplot(numd, numa,(numa)*(ii-1)+ii)
    text(0.25,0.25,xl{ii},'FontSize',20); axis off
end
subplot(numd, numa,2)
title(ti{qual-numVars},'FontSize',25)
colormap copper
saveas(gcf,strcat(savegin,...
    'circle_varcomp',ti{qual-numVars},namedesc,'.png'))
 end
end
%% Envelope plots
for ii=1:4
    for jj=[1:ii-1,ii+1:4]
        subplot(numd,numa,numa*(ii-1)+jj)
        hold on
        scatter(varVals{ii},varVals{jj},'o')
        
        xlabel(xl(ii+1))
        ylabel(xl(jj+1))
    end
end
%% Generate scatter plots for just ip3 and K_c
ii=1;jj=3;
resrange=(numRuns+1):(2*numRuns);
% resrange=find(results(resrange,1)<5)+numRuns;
for qual = [(numVars+1):11,19:20]
    dc=results(resrange,qual)/results(end,qual)*100-100;
    figure
    scatter(results(resrange,(ii))*10,...
        results(resrange,(jj))*2e-1,1000,dc,'.')
    colorbar
    title(strcat('% change in',{' '},ti{qual-numVars}),'FontSize',25)
    xlabel(strcat(xl(ii),{' '},'(\muM)'),'FontSize',20)
ylabel(strcat(xl(jj),{' '},'(\muM)'),'FontSize',20)
%     caxis([-100 600])
    saveas(gcf,strcat(savegin,'scatter_ip3_kc_',...
    ti{qual-numVars},namedesc,'.eps'),'epsc')
end

%% Generate contour plots for just ip3/kf and K_c
ii=1;jj=3;
numRuns=400;
resrange=401:800;
% resrange=1:400;
isip3=0;
if isip3
resrange1=find(results(resrange,1)<=150)+numRuns;
resrange2=find(results(resrange,3)<=60)+numRuns;
resrange=intersect(resrange1,resrange2);
else
    resrange1=find(results(resrange,1)<=100)+numRuns;
    resrange2=find(results(resrange,3)<=100)+numRuns;
    resrange=intersect(resrange1,resrange2);
end
% un6=unique(results(:,6));
% un6=un6(2:end)';
% un3=unique(results(:,3));
% un3=un3(3:end)';
% for ll=un3
% for kk=un6
% resrange=1:400;

% resrange1=find(results(resrange,6)==kk);
% resrange2=find(results(resrange1,3)==ll);
% resrange=resrange(resrange1(resrange2));
quals=[numVars+1,numVars+2,numVars+5,numVars+6];
for qual = quals%[(numVars+1):(numVars+6),(numVars+14):(numVars+15)]
    dc=results(resrange,qual)/results(end,qual);
    dc=dc*100-100;
    dc(dc>160)=159.9;
    %switch pt
    if isip3
        f = fit([results(resrange,(ii))*10, 2e-1*results(resrange,(jj))],dc,'linearinterp');
    else
        f = fit([results(resrange,(ii))*3e-2, 2e-1*results(resrange,(jj))],dc,'linearinterp');
    end
    
%     f = fit([results(resrange,(ii)), results(resrange,(jj))],dc*100-100,'linearinterp');
    if ismember(qual,[numVars+14,numVars+15])
        figure('position',[100,200,560,450]*1.04)
    else
        figure('position',[0,200,560,450])
    end
    plot(f,'Style','contour')
    hold on
    if isip3
        xcross=10;
    else
        xcross=15*3e-2;
    end
    plot(xcross,40*2e-1,'x','MarkerSize',10,'LineWidth',4)
    plot(xcross,40*2e-1,'x','MarkerSize',10,'LineWidth',4)
    plot(xcross,32*2e-1,'x','MarkerSize',10,'LineWidth',4)
    if ~isip3
    plot(80*3e-2,48*2e-1,'x','MarkerSize',10,'LineWidth',4)
    else
        plot(xcross,1*2e-1,'x','MarkerSize',10,'LineWidth',4)
    end
    plot(xcross,1*2e-1,'x','MarkerSize',10,'LineWidth',4)
    
    plot(xcross,40*2e-1,'o','MarkerSize',15,'Color',[0 0 0],'LineWidth',2)
    plot(xcross,32*2e-1,'o','MarkerSize',15,'Color',[0 0 0],'LineWidth',2)
    plot(xcross,1*2e-1,'o','MarkerSize',15,'Color',[0 0 0],'LineWidth',2)
    if ~isip3
    plot(80*3e-2,48*2e-1,'o','MarkerSize',15,'Color',[0 0 0],'LineWidth',2)
    end
    set(gca,'FontSize',16)
    view(90,-90)
    axis('square')
    cb=colorbar;
    cb.Title.String = '+';
    cb.Title.Position(2) = 345;
    cb.Title.Position(1) = 50;
    title(strcat('Change in',{' '},ti{qual-numVars},{' '},'(%)'),'FontSize',25)
    ylabel(strcat(xl(jj),{' '},'(\muM)'),'FontSize',20)
    caxis([-60 160])
%     caxis([-100 600])
%     if qual==6
%         caxis([-60 300])
%     elseif qual==7
%         caxis([-60 300])
%     elseif qual==10
%         caxis([-60 300])
%     elseif qual==11
%         caxis([-60 300])
%     end
colormap(bluewhitered)
if isip3
    axis([0 14 0 10])
    xlabel(strcat(xl(ii),{' '},'(\muM)'),'FontSize',20)
      saveas(gcf,strcat(savegin,'contour_ip3_K_c_',...
    ti{qual-numVars},namedesc,'.eps'),'epsc')
saveas(gcf,strcat(savegin,'contour_ip3_K_c_',...
    ti{qual-numVars},namedesc,'.fig'),'epsc')
else
    axis([0 2.8 0 10])
    xlabel(strcat(xl(ii),{' '},'(\mum^3/ms)'),'FontSize',20)
    saveas(gcf,strcat(savegin,'contour_kf_K_c_',...
    ti{qual-numVars},namedesc,'.eps'),'epsc')
saveas(gcf,strcat(savegin,'contour_kf_K_c_',...
    ti{qual-numVars},namedesc,'.fig'),'epsc')
end
%     saveas(gcf,strcat(savegin,'contour_ser_ncx_kf',num2str(kk),'_Kc',num2str(round(ll)),...
%     ti{qual-numVars},namedesc,'.eps'),'epsc')
end
% end
% end
%% Generate contour plots for just ip3 and K_h
if 0
    ii=1;jj=4;
numRuns=400;
resrange=401:800;
resrange1=find(results(resrange,1)<=1.5)+numRuns;
resrange2=find(results(resrange,4)<=30)+numRuns;
resrange=intersect(resrange1,resrange2);
for qual = [(numVars+1):(numVars+6),(numVars+14):(numVars+15)]
    dc=results(resrange,qual)/results(end,qual);
    %switch pt
    f = fit([results(resrange,(ii))*10, 8e-5*results(resrange,(jj))],dc*100-100,'linearinterp');
    if ismember(qual,[numVars+14,numVars+15])
        figure('position',[100,200,560,450]*1.04)
    else
        figure('position',[0,200,560,450])
    end
    plot(f,'Style','contour')
    hold on
    plot(10,40*2e-1,'x','MarkerSize',10,'LineWidth',4)
    plot(10,40*2e-1,'x','MarkerSize',10,'LineWidth',4)
    plot(10,20*2e-1,'x','MarkerSize',10,'LineWidth',4)
    plot(10,1*2e-1,'x','MarkerSize',10,'LineWidth',4)
    plot(10,1*2e-1,'x','MarkerSize',10,'LineWidth',4)
    plot(10,40*2e-1,'o','MarkerSize',15,'Color',[0 0 0],'LineWidth',2)
    plot(10,20*2e-1,'o','MarkerSize',15,'Color',[0 0 0],'LineWidth',2)
    plot(10,1*2e-1,'o','MarkerSize',15,'Color',[0 0 0],'LineWidth',2)
    set(gca,'FontSize',16)
    view(90,-90)
    axis('square')
    colorbar
    title(strcat('Change in',{' '},ti{qual-numVars},{' '},'(%)'),'FontSize',25)
%     xlabel(strcat(xl(ii),{' '},'(\mum^3/ms)'),'FontSize',20)
    xlabel(strcat(xl(ii),{' '},'(\muM)'),'FontSize',20)
    ylabel(strcat(xl(jj),{' '},'(\muM)'),'FontSize',20)
%     caxis([-100 600])
%     saveas(gcf,strcat(savegin,'contour_ip3_K_h_',...
%     ti{qual-numVars},namedesc,'.eps'),'epsc')
%     saveas(gcf,strcat(savegin,'contour_ser_ncx_kf',num2str(kk),'_Kc',num2str(round(ll)),...
%     ti{qual-numVars},namedesc,'.eps'),'epsc')
end
 end
% end
%% Quantile regression
% For K_c
for qual = [(numVars+1):11,17,18]
numa=2;
numd=2;
figure
for ii=1:4
    if prod(varVals{ii}==varVals{ii}(1))~=1
    subplot(numd,numa,ii)
    [sv,I]=sort(varVals{ii});
    k=linspace(round(min(varVals{ii})),round(max(varVals{ii})),100);
    X = [ones(numRuns,1),varVals{ii}];
    Z = zeros(length(varVals{ii}),length(k));
    for j=1:length(k)
        Z(:,j)=max(X(:,2)-k(j),0);
    end
    lme = fitlmematrix(X,results(((end-1)/2+1):(end-1),qual),Z,[],'CovariancePattern','Isotropic');
    % X=[X Z];
    % lme_fixed = fitlmematrix(X,dc,[],[]);
    % compare(lme,lme_fixed,'NSim',500,'CheckNesting',true)
    R=response(lme);
    plot(sv,R(I),'o','MarkerFaceColor',[0.8,0.8,0.8],...
        'MarkerEdgeColor',[0.8,0.8,0.8],'MarkerSize',4);
    hold on
    F=fitted(lme);
    plot(sv,F(I))
    xlabel(xl{ii})
    ylabel(strcat(ti{qual-numVars},yl{qual-numVars}))
    end
end
subplot(numd,numa,2)
legend('Simulations','Regression')
saveas(gcf,strcat(savegin,...
    'regressions',ti{qual-numVars},namedesc,'.png'))
end
%% Box plot
col=@(x)reshape(x,numel(x),1);
boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});
numa=2;
numd=2;
figure
numbins=10;
for ii=1:4
    subplot(numd,numa,ii)
    [sv,I]=sort(varVals{ii});
    yval=dc(I);
    [N,edges]=histcounts(sv,numbins);
    N=[0 cumsum(N)];
    N=unique(N);
    numbins=size(N,2)-1;
    bins=cell(1,numbins);
    labels=bins;
    if ii==1
        sv=sv*1e-8;
    end
    for jj=1:numbins;
        bins{jj}=yval(N(jj)+1:N(jj+1));
%         labels{jj}=strcat(num2str(sv(N(jj)+1)),':',num2str(sv(N(jj+1))));
    labels{jj}=(round(mean(sv((N(jj)+1):N(jj+1))),2,'significant'));
    end
    boxplot2(bins,'Labels',labels)
    if ii>1
        xlabel(xl{ii+1})
    else
        xlabel(strcat(xl{ii+1},{' '},'(x10^8)'))
    end
    ylabel('Duty Cycle')
end