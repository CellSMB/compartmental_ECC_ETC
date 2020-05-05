%% Make figures for compartmental model paper
% Hilary Hunt 2019
savepath='/Users/hhunt1/Documents/Show/Apr17/';
mkdir(savepath(1:(end-1)));
if 0
%% Result 2: Sneyd IP3R model channel fluxes
chg=[1,3];
ik=2;
for kk=1%:0.02:0.15
    for jj=-2%:2
% Generate sims
if 1
% params=[40*kk*chg(ik),1,80,27,0.1];
params=[10/3*chg(ik),10,40,6,10];
[VOIb, STATESc, ALGEBRAICc, CONSTANTSc] = IncDelayModel([0,params(2:end)]);
[VOI, STATES, ALGEBRAIC, CONSTANTS] = IncDelayModel(params);
end
%% Plot fig 5
% needs 
titletext={'Cytosol (Model Resolution)','SR','I_{RyR}','I_{SERCA}','I_{IP_3R}','I_{NCX}'};
ydesc={'[Ca^{2+}] (\muuM)','[Ca^{2+}] (\muM)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)'};
xdesc='Time (ms)';
pretrans = 100;
tstart=find(VOIb>VOIb(end)-CONSTANTSc(:,7)-pretrans,1,'first');
timec=tstart:size(VOIb,1);
valuesc={STATESc(timec,1),STATESc(timec,2),ALGEBRAICc(timec,35),-ALGEBRAICc(timec,48),ALGEBRAICc(timec,66),ALGEBRAICc(timec,54),ALGEBRAICc(timec,46)};
tvalc=VOIb(timec)-VOIb(tstart);
tstart=find(VOI>VOI(end)-CONSTANTS(:,7)-pretrans,1,'first');
time=tstart:size(VOI,1);
values={STATES(time,1),STATES(time,2),ALGEBRAIC(time,35),-ALGEBRAIC(time,48),ALGEBRAIC(time,66),ALGEBRAIC(time,54),ALGEBRAIC(time,46)};
tval=VOI(time)-VOI(tstart);
% Convert to uM
valuesc=cellfun(@(x)1000*x,valuesc,'UniformOutput',0);
values=cellfun(@(x)1000*x,values,'UniformOutput',0);
% Reduce resolution
res = 65;
offset=46;
pts=ceil((round(max(tval))-offset)/res);
tvallr=zeros(pts,1);
tvalclr=tvallr;
count=1;
for ii=offset:res:max(tval)
   tvallr(count)=find(tval>=ii,1,'first');
   tvalclr(count)=find(tvalc>=ii,1,'first');
   count=count+1;
end
tvallr(tvallr==0)=[];
tvalclr(tvalclr==0)=[];

if 1
numa=2; numd=3;
figure
for ii=[1,2,3:numa*numd]
   subplot(numd,numa,ii)
   hold on
   plot(tvalc,valuesc{ii},'-','LineWidth',2)
   plot(tval,values{ii},'--','LineWidth',2)
   title(titletext{ii})
   if ii>numa*numd-numa
       xlabel(xdesc)
   end
   if mod(ii,2)==1
       ylabel(ydesc{ii})
   end
   axis([0 CONSTANTS(:,7) min([values{ii};valuesc{ii}]) max([values{ii};valuesc{ii}])])
end
% subplot(numd,numa,2)
% hold on
% plot(tvalc(tvalclr),valuesc{1}(tvalclr),'-x','LineWidth',0.5)
% plot(tval(tvallr),values{1}(tvallr),'-x','LineWidth',0.5)
% title('Cytosol (16Hz)')
% axis([0 CONSTANTS(:,7) min([values{1};valuesc{1}]) max([values{1};valuesc{1}])])

subplot(numd,numa,1)
legend({'ECC', 'IP_3Rs activated'})
strparams=strcat(num2str(params(1)),'-',num2str(params(2)),'-',num2str(params(3)),'-',num2str(params(4)),'-',num2str(params(5)));
if 1
saveas(gcf,strcat(savepath,'fluxhip',strparams,'.eps'))
saveas(gcf,strcat(savepath,'fluxhip',strparams,'.fig'))
end
end
        end
    end
    end
% if 0
% %% Plot fluxes ramos-franco2000
% params=[0.1,1,1];
% [VOIc, STATESc, ALGEBRAICc, CONSTANTSc] = Fbfluo([0,params(2:end)]);
% [VOI, STATES, ALGEBRAIC, CONSTANTS] = Fbfluo(params);
% %% Plot
% 
% titletext={'Cytosol','','SR','I_{RyR}','I_{SERCA}','I_{IP_3R}'};
% ydesc={'[Ca^{2+}] (\muuM)','','[Ca^{2+}] (\muM)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)'};
% xdesc='Time (ms)';
% pretrans = 100;
% tstart=find(VOIc>VOIc(end)-CONSTANTSc(:,7)-pretrans,1,'first');
% timec=tstart:size(VOIc,1);
% valuesc={STATESc(timec,1),0,STATESc(timec,2),ALGEBRAICc(timec,35),-ALGEBRAICc(timec,48),ALGEBRAICc(timec,60),ALGEBRAICc(timec,54)};
% tvalc=VOIc(timec)-VOIc(tstart);
% tstart=find(VOI>VOI(end)-CONSTANTS(:,7)-pretrans,1,'first');
% time=tstart:size(VOI,1);
% values={STATES(time,1),0,STATES(time,2),ALGEBRAIC(time,35),-ALGEBRAIC(time,48),ALGEBRAIC(time,60),ALGEBRAIC(time,54)};
% tval=VOI(time)-VOI(tstart);
% % Convert to uM
% valuesc=cellfun(@(x)1000*x,valuesc,'UniformOutput',0);
% values=cellfun(@(x)1000*x,values,'UniformOutput',0);
% 
% % Reduce resolution
% res = 65;
% offset=46;
% pts=ceil((round(max(tval))-offset)/res);
% tvallr=zeros(pts,1);
% tvalclr=tvallr;
% count=1;
% for ii=offset:res:max(tval)
%    tvallr(count)=find(tval>=ii,1,'first');
%    tvalclr(count)=find(tvalc>=ii,1,'first');
%    count=count+1;
% end
% tvallr(tvallr==0)=[];
% tvalclr(tvalclr==0)=[];
% 
% if 1
% numa=2; numd=3;
% figure
% for ii=[1,3:numa*numd]
%    subplot(numd,numa,ii)
%    hold on
%    plot(tvalc,valuesc{ii},'-','LineWidth',2)
%    plot(tval,values{ii},'--','LineWidth',2)
%    title(titletext{ii})
%    if ii>numa*numd-numa
%        xlabel(xdesc)
%    end
%    if mod(ii,2)==1
%        ylabel(ydesc{ii})
%    end
%    axis([0 CONSTANTS(:,7) min([values{ii};valuesc{ii}]) max([values{ii};valuesc{ii}])])
% end
% subplot(numd,numa,2)
% hold on
% plot(tvalc(tvalclr),valuesc{1}(tvalclr),'-x','LineWidth',0.5)
% plot(tval(tvallr),values{1}(tvallr),'-x','LineWidth',0.5)
% title(strcat(titletext{1},{' '},'(16Hz)'))
% axis([0 CONSTANTS(:,7) min([values{1};valuesc{1}]) max([values{1};valuesc{1}])])
% 
% subplot(numd,numa,1)
% legend({'ECC', 'IP_3Rs activated'})
% if 1
% strparams=strcat(num2str(params(1)),'-',num2str(params(2)),'-',num2str(params(3)));
% saveas(gcf,strcat(savepath,'fluxrf',strparams,'.eps'))
% saveas(gcf,strcat(savepath,'fluxrf',strparams,'.fig'))
% end
% end
% end
    
% if 0
% %% Plot low res
% res = 65;
% offset=1;
% pts=ceil((round(max(tval))-offset)/res);
% for offset=1:(res-1)
% tvallr=zeros(pts,1);
% tvalclr=tvallr;
% count=1;
% for ii=offset:res:max(tval)
%    tvallr(count)=find(tval>=ii,1,'first');
%    tvalclr(count)=find(tvalc>=ii,1,'first');
%    count=count+1;
% end
% tvallr(tvallr==0)=[];
% tvalclr(tvalclr==0)=[];
% ii=1;
% noise=5;
% outvc=awgn(valuesc{ii},noise,'measured');
% outv=awgn(values{ii},noise,'measured');
% numa=2; numd=1;
% figure
% subplot(numd,numa,ii)
% hold on
% plot(tvalc,valuesc{ii},'LineWidth',2)
% plot(tval,values{ii},'-','LineWidth',2)
% title(titletext{ii})
% if ii>numa*numd-numa
%    xlabel(xdesc)
% end
% if mod(ii,2)==1
%    ylabel(ydesc{ii})
% end
% axis([0 CONSTANTS(:,7) min([values{ii};valuesc{ii}]) max([values{ii};valuesc{ii}])])
% subplot(numd,numa,2)
% hold on
% plot(tvalc(tvalclr),valuesc{ii}(tvalclr),'x','LineWidth',2)
% plot(tval(tvallr),values{ii}(tvallr),'x','LineWidth',2)
% % plot(tval,outv,'-','LineWidth',2)
% % plot(tvalc,outvc,'LineWidth',2)
% title(titletext{ii})
% if ii>numa*numd-numa
%    xlabel(xdesc)
% end
% if mod(ii,2)==1
%    ylabel(ydesc{ii})
% end
% axis([0 CONSTANTS(:,7) min([values{ii};valuesc{ii}]) max([values{ii};valuesc{ii}])])
% legend({'ECC', 'IP_3Rs activated'})
% 
% end
% end

%% Plot envelope
if 0
load('/Users/hhunt1/Documents/MATLAB/datameans.mat');

% Calculate differences between model and data
cellNum = 1;
nanCut = 2199;
compareTo = mean(meanTP(1:nanCut,cellNum,:),3);
compareToC = mean(meanTC(1:nanCut,cellNum,:),3);
dataTime = (1:size(compareTo,1))'*1.04;
% compareTo2 = interp1q(dataTime,compareTo,(1:dataTime(end))');
% compareToC2 = interp1q(dataTime,compareToC,(1:dataTime(end))');
% maxComp = max(compareTo);
% maxCompC = max(compareToC);
% gradC = diff(compareToC(1:10:nanCut));

% Find envelope borders
stdMult = 2;
nanCut = 2115;
cellNum=2;
stdC = std(meanTC(1:nanCut,cellNum,:),[],3);
mTC=reshape(meanTC(1:nanCut,:,:),[nanCut,1,25]);
stdC = std(mTC,[],3);
upperC = compareToC(1:nanCut) + stdMult*stdC;
lowerC = compareToC(1:nanCut) - stdMult*stdC;

stdP = std(meanTP(1:nanCut,cellNum,:),[],3);
upperP = compareTo(1:nanCut) + stdMult*stdP;
lowerP = compareTo(1:nanCut) - stdMult*stdP;

trngDat = dataTime(1:nanCut)';
% Calculate fluorescence intensity
K_d_cyt = 1.1e-3; R_f_cyt= 7;
% caf = ((K_d_cyt+STATESc(timec,1)*R_f_cyt)-(K_d_cyt+min(STATESc(timec,1))))./...
%                 ((K_d_cyt+STATESc(timec,1))*(K_d_cyt+R_f_cyt*min(STATESc(timec,1))));
caf = (STATESc(timec,1)*R_f_cyt+K_d_cyt)*(K_d_cyt+min(STATESc(timec,1)))./...
    ((STATESc(timec,1)+K_d_cyt)*(K_d_cyt+R_f_cyt*min(STATESc(timec,1))));
fval=interp1q([0; tvalc],[min(STATESc(timec,1)); caf],(1:2201)');
if 1
    % Plot with data and envelope
    figure
    hold on
    % envelope
%     h1 = fill([trngDat fliplr(trngDat)],[upperP' fliplr(lowerP')],[255 222 173]/255,'EdgeColor','none');
    h2 = fill([trngDat fliplr(trngDat)],[upperC' fliplr(lowerC')],[135 206 250]/255,'EdgeColor','none');
    % data
%     h3 = plot(trngDat,compareToC(1:nanCut),'LineWidth',3,'Color',[0 0.4470 0.7410]);
%     h4 = plot(trngDat,compareTo(1:nanCut),'LineWidth',3,'Color',[0.8500 0.3250 0.0980]);
    % model
    colorHere = parula(5);
    h5 = plot([ones(150,1);fval(1:(end-150))],'-','LineWidth',2,'Color',colorHere(1,:));
    
%     h6 = plot(modelPts{bestComboLocInd2},'-.','LineWidth',2,'Color','r');
    % h5 = plot(runCtrl,'LineWidth',3);
    % colorHere = parula(5);
    % h6 = plot(runIP3,'LineWidth',3,'Color',colorHere(1,:));
    legend([h5,h2],{'control model','95% CI control data'})
    xlabel('time (ms)')
    ylabel('[Ca^{2+}] (F/F_0)')
    yl = ylim;
end
end

if 0
%% Result 3: Compare difference channel numbers
% Generate sims
if 0
    titletext={'Cytosol','','SR','I_{RyR}','I_{SERCA}','I_{IP_3R}'};
    ydesc={'[Ca^{2+}] (\muuM)','','[Ca^{2+}] (\muM)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)'};
    xdesc='Time (ms)';
    
params=[40,1,80,27,0.1,1e-2,3,0.7];
lowip3=5e-3;
numvars=9;
res = 65;
prm = [[0,params(2:end-2),1,1];[params(1:end-2),1,1];[params(1:end-1),1];params;...
    [0,params(2:end-3),lowip3,1,0.7];[params(1:end-3),lowip3,1,1];...
    [params(1:end-3),lowip3,3,1];[params(1:(end-3)),lowip3,params((end-1):end)];...
    [0,params(2:end-2),1,1]];
sims=cell(numvars,3);
simres=cell(numvars,3);
for kk=1:numvars
[VOI, STATES, ALGEBRAIC, CONSTANTS] = IncDelayModelNumbers(prm(kk,:));
sims{kk,1}=VOI;
sims{kk,2}=STATES;
sims{kk,3}=ALGEBRAIC;
if ismember(kk,[3,4,7,8,9])
    pretrans = 300;
    offset=51;
else
    pretrans = 100;
    offset=46;
end
tstart=find(VOI>VOI(end)-CONSTANTS(:,7)-pretrans,1,'first');
time=tstart:size(VOI,1);
% values={STATES(time,1),0,STATES(time,2),ALGEBRAIC(time,35),-ALGEBRAIC(time,48),ALGEBRAIC(time,66),ALGEBRAIC(time,54)};
tval=VOI(time)-VOI(tstart);
% Convert to uM
% values=cellfun(@(x)1000*x,values,'UniformOutput',0);
% Reduce resolution
pts=ceil((round(max(tval))-offset)/res);
tvallr=zeros(pts,1);
count=1;
for ii=offset:res:max(tval)
   tvallr(count)=find(tval>=ii,1,'first');
   count=count+1;
end
tvallr(tvallr==0)=[];
simres{kk,1}=tval;
simres{kk,2}=STATES(time,1)*1e6;
simres{kk,3}=tvallr;
end
end
if 0
%% Recreate Harzheim Fig
fpos=[400 250 560 140];
figure('pos',fpos)
subplot(1,2,1)
hold on
ik=6;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
ik=7;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
axis([0 CONSTANTS(:,7) 0 max(simres{3,2})])
xlabel('time (ms)')
title('5\muM IP_3')
ylabel('[Ca^{2+}] (nM)')
subplot(1,2,2)
hold on
ik=1;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
ik=9;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
axis([0 CONSTANTS(:,7) 0 max(simres{3,2})])
title('IP_3Rs inhibited')
xlabel('time (ms)')
legend({'ctrl','300% IP3Rs'})
ax=gca;
% ax.FontSize=12;

figure('pos',fpos)
subplot(1,2,1)
hold on
ik=6;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
ik=7;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
axis([0 CONSTANTS(:,7) 0 max(simres{3,2})])
xlabel('time (ms)')
title('5\muM IP_3')
ylabel('[Ca^{2+}] (nM)')
subplot(1,2,2)
hold on
ik=2;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
ik=3;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
title('10\muM IP_3')
axis([0 CONSTANTS(:,7) 0 max(simres{3,2})])
xlabel('time (ms)')
legend('ctrl','300% IP3Rs')

figure('pos',fpos)
subplot(1,2,1)
hold on
ik=6;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
ik=5;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
axis([0 CONSTANTS(:,7) 0 max(simres{3,2})])
title('5\muM IP_3')
ylabel('[Ca^{2+}] (nM)')
xlabel('time (ms)')
subplot(1,2,2)
hold on
ik=2;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
ik=4;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
title('10\muM IP_3')
axis([0 CONSTANTS(:,7) 0 max(simres{3,2})])
xlabel('time (ms)')
legend('ctrl','300% IP3Rs, 70% RyRs')
end
if 1
%% Recreate Harzheim Fig with dc
set(0,'defaultAxesFontSize',14)
auc=zeros(1,9);
maxauc=auc;
for ik=[1,2,3,4,5,6,7,9]
    auc(ik)=1e-5*(trapz(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3})));
    maxauc(ik)= 1e-5*max(simres{ik,1}(simres{ik,3}))*max(simres{1,2}(simres{1,3}));
end
legpos=[0.4 -0.015 0.253571428571429 0.0690476190476191];
    
figure
subplot(1,2,1)
hold on
ik=6;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
ik=7;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
axis([0 CONSTANTS(:,7) 0 max(simres{3,2})])
xlabel('time (ms)')
title('5\muM IP_3')
ylabel('[Ca^{2+}] (nM)')
legstr={strcat('\gamma =',32,num2str(round(auc(6)/maxauc(6),3)),32,...
    '\gamma_{inhibited} =',32,num2str(round(auc(1)/maxauc(1),3))),...
    strcat('\gamma =',32,num2str(round(auc(7)/maxauc(7),3)),32,...
    '\gamma_{inhibited} =',32,num2str(round(auc(9)/maxauc(9),3)))};
legend(legstr,'Units','normalized','Position',legpos,'Orientation','horizontal')
legend('boxoff')
subplot(1,2,2)
hold on
ik=1;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
ik=9;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
axis([0 CONSTANTS(:,7) 0 max(simres{3,2})])
title('IP_3Rs inhibited')
xlabel('time (ms)')
legend('ctrl','300% IP3Rs')

figure
subplot(1,2,1)
hold on
ik=6;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
ik=7;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
axis([0 CONSTANTS(:,7) 0 max(simres{3,2})])
xlabel('time (ms)')
title('5\muM IP_3')
ylabel('[Ca^{2+}] (nM)')
legstr={strcat('\gamma =',32,num2str(round(auc(6)/maxauc(6),3)),32,...
    '\gamma_{IP_3+} =',32,num2str(round(auc(2)/maxauc(2),3))),...
    strcat('\gamma =',32,num2str(round(auc(7)/maxauc(7),3)),32,...
    '\gamma_{IP_3+} =',32,num2str(round(auc(3)/maxauc(3),3)))};
legend(legstr,'Units','normalized','Position',legpos,'Orientation','horizontal')
legend('boxoff')
subplot(1,2,2)
hold on
ik=2;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
ik=3;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
title('10\muM IP_3')
axis([0 CONSTANTS(:,7) 0 max(simres{3,2})])
xlabel('time (ms)')
legend('ctrl','300% IP3Rs')

figure
subplot(1,2,1)
hold on
ik=6;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
ik=5;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
axis([0 CONSTANTS(:,7) 0 max(simres{3,2})])
title('5\muM IP_3')
ylabel('[Ca^{2+}] (nM)')
xlabel('time (ms)')
legstr={strcat('\gamma =',32,num2str(round(auc(6)/maxauc(6),3)),32,...
    '\gamma_{IP_3+} =',32,num2str(round(auc(2)/maxauc(2),3))),...
    strcat('\gamma =',32,num2str(round(auc(5)/maxauc(5),3)),32,...
    '\gamma_{IP_3+} =',32,num2str(round(auc(4)/maxauc(4),3)))};
legend(legstr,'Units','normalized','Position',legpos,'Orientation','horizontal')
legend('boxoff')
subplot(1,2,2)
hold on
ik=2;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
ik=4;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
title('10\muM IP_3')
axis([0 CONSTANTS(:,7) 0 max(simres{3,2})])
xlabel('time (ms)')
legend('ctrl','300% IP3Rs, 75% RyRs')
%% plot duty cycle/auc
figure
numa=3;numd=2;
auc=zeros(1,n);
maxauc=zeros(1,n);
legpos={[0.0891071428571429 0.475 0.253571428571429 0.0690476190476191],...
    [0.369464285714286 0.475 0.253571428571429 0.0690476190476191],...
    [0.654285714285714 0.475 0.253571428571429 0.0690476190476191],...
    [0.0891071428571429 0 0.253571428571429 0.0690476190476191],...
    [0.369464285714286 0 0.253571428571429 0.0690476190476191],...
    [0.654285714285714 0 0.253571428571429 0.0690476190476191]};
for ii=1:n
subplot(numd,numa,ii)
hold on
area(simres{ii,1},simres{ii,2},'FaceColor',[0.7 0.8 1])
plot(simres{ii,1},simres{ii,2},'LineWidth',2)
axis([0 CONSTANTS(:,7) 0 1000])
auc(ii)=1e-5*(trapz(simres{ii,1},simres{ii,2}));
maxauc(ii)= 1e-5*max(simres{ii,1})*max(simres{ii,2});
legstr={strcat('AUC =',32,num2str(auc(ii))...
    ,'x10^5nMms'),strcat('\gamma =',32,num2str(auc(ii)/maxauc(ii)))};
%     -trapz(simres{ii,1},min(simres{ii,2})*ones(size(simres{ii,1})))))
legend(legstr,'Units','normalized','Position',legpos{ii})
legend('boxoff')
end
end
end
if 0
%% Compare IP3 concs (sensitivity)
% [VOI, STATES, ALGEBRAIC, CONSTANTS] = IncDelayModelNumbers([40,1,80,27,0.1,3e-5,3,1]);
num=6;
ip=round([0, (logspace(-3.5,-1.5,num))],1,'significant');
% ip=[0, (linspace(3e-5,1e-1,num))];
tv=cell(num,1);
cv=cell(num,1);
% [0,1,80,27,0.1,1e-2,1,1]
%[40,1,60,8,0.1,ip(ii),3,1]
for ii = 1:num+1
    sh=ii*200;
    [VOI, STATES, ALGEBRAIC, CONSTANTS] = IncDelayModelNumbers([40,1,80,27,0.1,ip(ii),3,1]);
    tv{ii} = VOI(find(VOI>=VOI(end)-CONSTANTS(:,7)-sh,1,'first'):end)-VOI(end)+CONSTANTS(:,7)+sh;
    cv{ii} = 1e6*STATES(find(VOI>=VOI(end)-CONSTANTS(:,7)-sh,1,'first'):end,1);
end

figure('pos',[400 250 940 235]);
hold on
for ii=1:num+1
    plot(tv{ii},cv{ii},'LineWidth',2)
%     legend('no IP3Rs','[IP_3]=10\muM','[IP_3]=30nM')
    ylabel('[Ca^{2+}] (nM)')
    xlabel('time (ms)')
end
axis([0 CONSTANTS(:,7) 0 1600])
legend(arrayfun(@(x) strcat('[IP_3]=',num2str(1000*x),'\muM'),ip,'UniformOutput',0))
% ik=1;plot(simres{ik,1},simres{ik,2}(:,1),'LineWidth',2)
% ik=3;plot(simres{ik,1},simres{ik,2}(:,1),'LineWidth',2)
% sh=200;
% 
% plot(VOI(find(VOI>=VOI(end)-CONSTANTS(:,7)-sh,1,'first'):end)-VOI(end)+CONSTANTS(:,7)+sh,1e6*STATES(find(VOI>=VOI(end)-CONSTANTS(:,7)-sh,1,'first'):end,1),'LineWidth',2)
end
if 0
if 1
%% Result 5: Duty cycle
par=[40,1,76,10,1];
lowip=5e-3;
ams={[0,1,1],[lowip,1,1],[1e-2,1,1],[0,3,1],[lowip,3,1],[1e-2,3,1]};
n=size(ams,2);
sims=cell(n,3);
simres=cell(n,2);
pretrans=100;
offset=46;
res=65;
for ii=1:n
    [VOI, STATES, ALGEBRAIC, CONSTANTS] = IncDelayModelNumbers([par,ams{ii}]);
    sims{ii,1}=VOI;
    sims{ii,2}=STATES;
    sims{ii,3}=ALGEBRAIC;
    
    tstart=find(VOI>VOI(end)-CONSTANTS(:,7)-pretrans,1,'first');
    time=tstart:size(VOI,1);
    % values={STATES(time,1),0,STATES(time,2),ALGEBRAIC(time,35),-ALGEBRAIC(time,48),ALGEBRAIC(time,66),ALGEBRAIC(time,54)};
    tval=VOI(time)-VOI(tstart);
    % Convert to \muM
    % values=cellfun(@(x)1000*x,values,'UniformOutput',0);
    % Reduce resolution
    pts=ceil((round(max(tval))-offset)/res);
    tvallr=zeros(pts,1);
    simres{ii,1}=tval;
    simres{ii,2}=STATES(time,1)*1e6;
end
end
if 0
%% plot duty cycle/auc
set(0,'defaultAxesFontSize',14)
figure
numa=3;numd=2;
auc=zeros(1,n);
maxauc=zeros(1,n);
legpos={[0.0891071428571429 0.475 0.253571428571429 0.0690476190476191],...
    [0.369464285714286 0.475 0.253571428571429 0.0690476190476191],...
    [0.654285714285714 0.475 0.253571428571429 0.0690476190476191],...
    [0.0891071428571429 0 0.253571428571429 0.0690476190476191],...
    [0.369464285714286 0 0.253571428571429 0.0690476190476191],...
    [0.654285714285714 0 0.253571428571429 0.0690476190476191]};
for ii=1:n
subplot(numd,numa,ii)
hold on
area(simres{ii,1},simres{ii,2},'FaceColor',[0.7 0.8 1])
plot(simres{ii,1},simres{ii,2},'LineWidth',2)
axis([0 CONSTANTS(:,7) 0 1000])
auc(ii)=1e-5*(trapz(simres{ii,1},simres{ii,2}));
maxauc(ii)= 1e-5*max(simres{ii,1})*max(simres{4,2});
legstr={strcat('AUC =',32,num2str(auc(ii))...
    ,'x10^5nMms'),strcat('\gamma =',32,num2str(auc(ii)/maxauc(ii)))};
%     -trapz(simres{ii,1},min(simres{ii,2})*ones(size(simres{ii,1})))))
legend(legstr,'Units','normalized','Position',legpos{ii},'FontSize',12)
legend('boxoff')
end
% for ii=[1,4]
%     subplot(numd,numa,ii)
%     ylabel('[Ca^{2+}] (nM)')
% end
% for ii=4:6
%     subplot(numd,numa,ii)
%     xlabel('time (ms)')
% end
end
if 1
    %% plot duty cycle/auc
legpos={[0.0891071428571429 0.615 0.253571428571429 0.0690476190476191],...
    [0.369464285714286 0.615 0.253571428571429 0.0690476190476191],...
    [0.654285714285714 0.615 0.253571428571429 0.0690476190476191],...
    [0.0891071428571429 0 0.253571428571429 0.0690476190476191],...
    [0.369464285714286 0 0.253571428571429 0.0690476190476191],...
    [0.654285714285714 0 0.253571428571429 0.0690476190476191]
    };
    
set(0,'defaultAxesFontSize',14)
figure('pos',[400 250 560 630])
numa=3;numd=3;
dc1=[auc([1,4]);auc([2,5]);auc([3,6])]./[maxauc([1,4]);maxauc([2,5]);maxauc([3,6])];
for ii=1:n
subplot(numd,numa,ii+3*floor(ii/4))
hold on
area(simres{ii,1},simres{ii,2},'FaceColor',[0.7 0.8 1])
plot(simres{ii,1},simres{ii,2},'LineWidth',2)
axis([0 CONSTANTS(:,7) 0 1000])
legstr={strcat('AUC =',32,num2str(round(auc(ii),1))...
    ,'x10^5nMms'),strcat('\Delta\gamma =',32,'+',num2str(round(100*(dc1(ii)/dc1(1)-1))),'%')};
%     -trapz(simres{ii,1},min(simres{ii,2})*ones(size(simres{ii,1})))))
legend(legstr,'Units','normalized','Position',legpos{ii},'FontSize',12)
legend('boxoff')
% daspect manual
end
end
if 1
%% plot bar graph
% c={'(ctrl, no IP_3)','(ctrl, low IP_3)','(ctrl, high IP_3)',...
%     '(+123% IP_3Rs, no IP_3)','(+123% IP_3Rs, low IP_3)','(+123% IP_3Rs, high IP_3)'};
c={'(no IP_3)','(low IP_3)','(high IP_3)'};

auc1=[auc([1,4]);auc([2,5]);auc([3,6])];
figure
b=bar(auc1);
set(gca,'xticklabel',c)
ylabel('AUC (nMms)')
legend('ctrl','+123% IP_3Rs')

%% plot bar graph dc
% c={'(ctrl, no IP_3)','(ctrl, low IP_3)','(ctrl, high IP_3)',...
%     '(+123% IP_3Rs, no IP_3)','(+123% IP_3Rs, low IP_3)','(+123% IP_3Rs, high IP_3)'};
c={'(no IP_3)','(low IP_3)','(high IP_3)'};
set(0,'defaultAxesFontSize',14)
dc1=[auc([1,4]);auc([2,5]);auc([3,6])]./[maxauc([1,4]);maxauc([2,5]);maxauc([3,6])];
dc=100*(dc1/dc1(1)-1)
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple lines using matrix input to bar
bar1 = bar(dc);
set(bar1(1),'DisplayName','ctrl',...
    'FaceColor',[0 0.600000023841858 0.400000005960464]);
set(bar1(2),'DisplayName','300% IP_3Rs',...
    'FaceColor',[1 0.400000005960464 0]);

% Create ylabel
ylabel('Increase in duty cycle (%)','FontSize',14);

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'XTick',[1 2 3],'XTickLabel',...
    {'(no IP_3)','(low IP_3)','(high IP_3)'});
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.137090773809524 0.813468555453839 0.252604166666667 0.0861581920903954]);
end
if 0
%% plot
figure
numa=2;numd=1;
subplot(numd,numa,1)
for ii=1:3
hold on
plot(simres{ii,1},simres{ii,2},'LineWidth',2)
axis([0 CONSTANTS(:,7) 0 1000])
end
subplot(numd,numa,2)
for ii=4:6
hold on
plot(simres{ii,1},simres{ii,2},'LineWidth',2)
axis([0 CONSTANTS(:,7) 0 1000])
end
end
end
if 0
%% Plot IP3R gating as heat map
k=0.324;K_p=2e-4;
K_c=(1:5:100)*2e-4;
% K_h=0:0.1:60;
K_h=[1 20 27 60]*8e-5;

% ca=[0.1,0.5,1]*1e-3;
ca=0:1e-4:1.2e-3;
ip=[0.03,5,10]*1e-3;
alltog=size(K_h,2)*size(ip,2);
P_open=cell(1,alltog);
alpha=P_open;
beta=P_open;
for jj=1:size(ip,2)
    for kk=1:size(K_h,2)
        for ii=1:size(ca,2)
            B=ip(jj)^2./(K_p.^2+ip(jj)^2);
            m=ca(ii).^4./(K_c.^4+ca(ii).^4);
            h=K_h(kk).^4./(K_h(kk).^4+ca(ii).^4);
            beta{kk+3*jj-3}(:,ii)=B.*m;
            alpha{kk+3*jj-3}(:,ii)=(1-B).*(1-m*h);
            P_open{kk+4*jj-4}(:,ii)=beta{kk+3*jj-3}(:,ii)./(beta{kk+3*jj-3}(:,ii)+k*(beta{kk+3*jj-3}(:,ii)+alpha{kk+3*jj-3}(:,ii)));
        end
    end
end
%% plot
% set(0,'defaultAxesFontSize',14)
intvar=P_open;
figure
for ii=1:alltog
subplot(size(K_h,2),3,ii)
% imagesc(ca*1e3,K_c*1e3,intvar{ii})%,'edgecolor','none')
contourf(ca*1e3,K_c*1e3,intvar{ii})%,'edgecolor','none')
% mesh(ca*1e3,K_c,(P_open{ii}))
%view(-90,90)
%axis([0 1.2 0 100*2e-4 -60 0]) 
xlabel('Ca^{2+} (\muM)')
ylabel('K_c (\muM)')
zlabel('log(P_O)')
colorbar
caxis([0 0.5])
end
end
if 0
%% Result 3: Compare difference channel numbers
% Generate sims
if 1
    titletext={'cytosol','','SR','I_{RyR}','I_{SERCA}','I_{IP_3R}'};
    ydesc={'[Ca^{2+}] (uM)','','[Ca^{2+}] (uM)','flux (uM/ms)','flux (uM/ms)','flux (uM/ms)'};
    xdesc='time (ms)';
    
params=[40,1,80,27,0.1,1e-2,3,0.77];
lowip3=5e-3;
numvars=4;
res = 65;
prm = [[params(1:end-2),1,1];params;...
    [0,params(2:end-3),lowip3,1,0.77];[params(1:end-3),lowip3,1,1]];
sims=cell(numvars,3);
simres=cell(numvars,3);
for kk=1:numvars
[VOI, STATES, ALGEBRAIC, CONSTANTS] = IncDelayModelNumbers(prm(kk,:));
sims{kk,1}=VOI;
sims{kk,2}=STATES;
sims{kk,3}=ALGEBRAIC;
if ismember(kk,[2,3])
    pretrans = 300;
    offset=51;
else
    pretrans = 100;
    offset=46;
end
tstart=find(VOI>VOI(end)-CONSTANTS(:,7)-pretrans,1,'first');
time=tstart:size(VOI,1);
% values={STATES(time,1),0,STATES(time,2),ALGEBRAIC(time,35),-ALGEBRAIC(time,48),ALGEBRAIC(time,66),ALGEBRAIC(time,54)};
tval=VOI(time)-VOI(tstart);
% Convert to uM
% values=cellfun(@(x)1000*x,values,'UniformOutput',0);
% Reduce resolution
pts=ceil((round(max(tval))-offset)/res);
tvallr=zeros(pts,1);
count=1;
for ii=offset:res:max(tval)
   tvallr(count)=find(tval>=ii,1,'first');
   count=count+1;
end
tvallr(tvallr==0)=[];
simres{kk,1}=tval;
simres{kk,2}=STATES(time,1)*1e6;
simres{kk,3}=tvallr;
end
end
%% Recreate Harzheim Fig

figure
subplot(1,2,1)
hold on
ik=4;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
ik=3;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
axis([0 CONSTANTS(:,7) 0 max(simres{1,2})])
title('5uM IP_3')
ylabel('[Ca^{2+}] (nM)')
xlabel('time (ms)')
subplot(1,2,2)
hold on
ik=1;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
ik=2;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
title('10uM IP_3')
axis([0 CONSTANTS(:,7) 0 max(simres{1,2})])
xlabel('time (ms)')
legend('ctrl','+123% IP3Rs, -31% RyRs')

end

%% Result 3: Compare difference channel numbers
% Generate sims
if 1
    titletext={'Cytosol','','SR','I_{RyR}','I_{SERCA}','I_{IP_3R}'};
    ydesc={'[Ca^{2+}] (\muuM)','','[Ca^{2+}] (\muM)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)'};
    xdesc='Time (ms)';
lip3=1e-4;
paramsctl=[lip3,10,40,6,10,0,1,1,1,1];
paramshyp=[lip3,10,40,6,10,30,1,1,0.8,1.5];
paramshf=[lip3,10,40,6,10,10,1,0.7,0.63,1.4];

hip3=1e-2;
numvars=8;
res = 5;
% 1 ctrl no ip 2 ctrl hip 3 hyp hip 4 hf hip 5 hf lip 6 ctrl lip 7 hyp lip 9 hyp noip
prm = [[0,paramsctl(2:end)];[hip3,paramsctl(2:end)];[hip3,paramshyp(2:end)];...
    [hip3,paramshf(2:end)];...
    paramshf;paramsctl;...
    paramshyp;[0,paramshyp(2:end)]];
sims=cell(numvars,3);
simres=cell(numvars,3);
for kk=1:numvars
[VOI, STATES, ALGEBRAIC, CONSTANTS] = IncDelayModelNumbers(prm(kk,:));
sims{kk,1}=VOI;
sims{kk,2}=STATES;
sims{kk,3}=ALGEBRAIC;
if ismember(kk,[3,4,7,8,8])
    pretrans = 200;
    offset=51;
else
    pretrans = 100;
    offset=46;
end
tstart=find(VOI>VOI(end)-CONSTANTS(:,7)-pretrans,1,'first');
time=tstart:size(VOI,1);
% values={STATES(time,1),0,STATES(time,2),ALGEBRAIC(time,35),-ALGEBRAIC(time,48),ALGEBRAIC(time,66),ALGEBRAIC(time,54)};
tval=VOI(time)-VOI(tstart);
% Convert to uM
% values=cellfun(@(x)1000*x,values,'UniformOutput',0);
% Reduce resolution
pts=ceil((round(max(tval))-offset)/res);
tvallr=zeros(pts,1);
count=1;
for ii=offset:res:max(tval)
   tvallr(count)=find(tval>=ii,1,'first');
   count=count+1;
end
tvallr(tvallr==0)=[];
simres{kk,1}=tval;
simres{kk,2}=STATES(time,1)*1e6;
simres{kk,3}=tvallr;
end

%% Recreate Harzheim Fig
fpos=[400 250 560 140];
figure('pos',fpos)
subplot(1,2,1)
hold on
ik=6;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
ik=7;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
axis([0 CONSTANTS(:,7) 0 max(simres{3,2})])
xlabel('time (ms)')
title('5\muM IP_3')
ylabel('[Ca^{2+}] (nM)')
subplot(1,2,2)
hold on
ik=1;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
ik=8;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
axis([0 CONSTANTS(:,7) 0 max(simres{8,2})])
title('IP_3Rs inhibited')
xlabel('time (ms)')
legend({'ctrl','300% IP3Rs'})
ax=gca;
% ax.FontSize=12;
if 0
figure('pos',fpos)
subplot(1,2,1)
hold on
ik=6;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
ik=7;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
axis([0 CONSTANTS(:,7) 0 max(simres{3,2})])
xlabel('time (ms)')
title('5\muM IP_3')
ylabel('[Ca^{2+}] (nM)')
subplot(1,2,2)
hold on
ik=2;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
ik=3;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
title('10\muM IP_3')
axis([0 CONSTANTS(:,7) 0 max(simres{3,2})])
xlabel('time (ms)')
legend('ctrl','300% IP3Rs')

figure('pos',fpos)
subplot(1,2,1)
hold on
ik=6;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
ik=5;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
axis([0 CONSTANTS(:,7) 0 max(simres{3,2})])
title('5\muM IP_3')
ylabel('[Ca^{2+}] (nM)')
xlabel('time (ms)')
subplot(1,2,2)
hold on
ik=2;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
ik=4;plot(simres{ik,1}(simres{ik,3}),simres{ik,2}(simres{ik,3}),'LineWidth',2)
title('10\muM IP_3')
axis([0 CONSTANTS(:,7) 0 max(simres{3,2})])
xlabel('time (ms)')
legend('ctrl','300% IP3Rs, 70% RyRs')
end
end

%% Plot transients from all regions of Agne's map - ip3
% Generate sims
if 1
% strparams=strcat(num2str(params(1)),'-',num2str(params(2)),'-',num2str(params(3)),'-',num2str(params(4)),'-',num2str(params(5)));
strparams='varyingIP3';
paramsb=[0,1,1,1,1,0,1,1,1,1];
params1=[1,1,1,1,1,15,1,1,1,1];
params2=[1,10,32,6,10,15,1,1,1,1];
params3=[1,10,40,6,10,15,1,1,1,1];
params4=[];

% strparams='varyingkf';
% paramsb=[0,1,1,1,1,0,1,1,1,1];
% params1=[1,1,1,1,1,15,1,1,1,1];
% params2=[1,10,32,6,10,15,1,1,1,1];
% params3=[1,10,40,6,10,15,1,1,1,1];
% params4=[1,10,48,6,10,80,1,1,1,1];

maxRedo=1e3;
% Base case
[VOIb, STATESb, ALGEBRAICb, CONSTANTSb] = ...
    IncDelayModelNumbers_PS(paramsb,[]);
period=CONSTANTSb(:,7);
numOsc=round(VOIb(end)/period);
penUpts=unique([VOIb-(numOsc-2)*period,STATESb(:,1)*1e3],'rows','stable');     
dbleup=diff(penUpts(:,1));     
if sum(dbleup==0)>0         
    penUpts(dbleup==0,:)=[];     
end
penUlt=interp1(penUpts(:,1),penUpts(:,2),0:(period-1));
Upts=unique([VOIb-(numOsc-1)*period,STATESb(:,1)*1e3],'rows','stable');
dbleup=diff(Upts(:,1));
if sum(dbleup==0)>0
    Upts(dbleup==0,:)=[];
end
ult=interp1(Upts(:,1),Upts(:,2),0:(period-1));
oscDiff=sum((ult-penUlt).^2/range(ult).^2);
redoNum=0;
while oscDiff>1e-2 & redoNum<maxRedo
    [VOIb, STATESb, ALGEBRAICb, CONSTANTSb] = ...
    IncDelayModelNumbers_PS(paramsb,STATESb(end,:));
    redoNum=redoNum+1;
    penUpts=unique([VOIb-(numOsc-2)*period,STATESb(:,1)*1e3],'rows','stable');     
    dbleup=diff(penUpts(:,1));     
    if sum(dbleup==0)>0         
        penUpts(dbleup==0,:)=[];     
    end
    penUlt=interp1(penUpts(:,1),penUpts(:,2),0:(period-1));
    Upts=unique([VOIb-(numOsc-1)*period,STATESb(:,1)*1e3],'rows','stable');
    dbleup=diff(Upts(:,1));
    if sum(dbleup==0)>0
        Upts(dbleup==0,:)=[];
    end
    ult=interp1(Upts(:,1),Upts(:,2),0:(period-1));
    oscDiff=sum((ult-penUlt).^2/range(ult).^2);
end
if redoNum>=maxRedo
    warning(strcat('Not great clarity at base'));
    unstable(2*numRuns+1)=1;
end
lastTransientb=find(VOIb>=(VOIb(end)-1*period)&VOIb<=(VOIb(end)-0.8*period));
tvecb=VOIb(lastTransientb)-(VOIb(end)-period);
% 1st case
[VOI1, STATES1, ALGEBRAIC1, CONSTANTS1] = ...
    IncDelayModelNumbers_PS(params1,[]);
period=CONSTANTS1(:,7);
numOsc=round(VOI1(end)/period);
penUpts=unique([VOI1-(numOsc-2)*period,STATES1(:,1)*1e3],'rows','stable');     
dbleup=diff(penUpts(:,1));     
if sum(dbleup==0)>0         
    penUpts(dbleup==0,:)=[];     
end
penUlt=interp1(penUpts(:,1),penUpts(:,2),0:(period-1));
Upts=unique([VOI1-(numOsc-1)*period,STATES1(:,1)*1e3],'rows','stable');
dbleup=diff(Upts(:,1));
if sum(dbleup==0)>0
    Upts(dbleup==0,:)=[];
end
ult=interp1(Upts(:,1),Upts(:,2),0:(period-1));
oscDiff=sum((ult-penUlt).^2/range(ult).^2);
redoNum=0;
while oscDiff>1e-2 & redoNum<maxRedo
    [VOI1, STATES1, ALGEBRAIC1, CONSTANTS1] = ...
    IncDelayModelNumbers_PS(params1,STATES1(end,:));
    redoNum=redoNum+1;
    penUpts=unique([VOI1-(numOsc-2)*period,STATES1(:,1)*1e3],'rows','stable');     
    dbleup=diff(penUpts(:,1));     
    if sum(dbleup==0)>0         
        penUpts(dbleup==0,:)=[];     
    end
    penUlt=interp1(penUpts(:,1),penUpts(:,2),0:(period-1));
    Upts=unique([VOI1-(numOsc-1)*period,STATES1(:,1)*1e3],'rows','stable');
    dbleup=diff(Upts(:,1));
    if sum(dbleup==0)>0
        Upts(dbleup==0,:)=[];
    end
    ult=interp1(Upts(:,1),Upts(:,2),0:(period-1));
    oscDiff=sum((ult-penUlt).^2/range(ult).^2);
end
if redoNum>=maxRedo
    warning(strcat('Not great clarity at base'));
    unstable(2*numRuns+1)=1;
end
lastTransient1=find(VOI1>=(VOI1(end)-1*period)&VOI1<=(VOI1(end)-0.8*period));
tvec1=VOI1(lastTransient1)-(VOI1(end)-period);
% 2nd case
[VOI2, STATES2, ALGEBRAIC2, CONSTANTS2] = ...
    IncDelayModelNumbers_PS(params2,[]);
period=CONSTANTS2(:,7);
numOsc=round(VOI2(end)/period);
penUpts=unique([VOI2-(numOsc-2)*period,STATES2(:,1)*1e3],'rows','stable');     
dbleup=diff(penUpts(:,1));     
if sum(dbleup==0)>0         
    penUpts(dbleup==0,:)=[];     
end
penUlt=interp1(penUpts(:,1),penUpts(:,2),0:(period-1));
Upts=unique([VOI2-(numOsc-1)*period,STATES2(:,1)*1e3],'rows','stable');
dbleup=diff(Upts(:,1));
if sum(dbleup==0)>0
    Upts(dbleup==0,:)=[];
end
ult=interp1(Upts(:,1),Upts(:,2),0:(period-1));
oscDiff=sum((ult-penUlt).^2/range(ult).^2);
redoNum=0;
while oscDiff>1e-2 & redoNum<maxRedo
    [VOI2, STATES2, ALGEBRAIC2, CONSTANTS2] = ...
    IncDelayModelNumbers_PS(params2,STATES2(end,:));
    redoNum=redoNum+1;
    penUpts=unique([VOI2-(numOsc-2)*period,STATES2(:,1)*1e3],'rows','stable');     
    dbleup=diff(penUpts(:,1));     
    if sum(dbleup==0)>0         
        penUpts(dbleup==0,:)=[];     
    end
    penUlt=interp1(penUpts(:,1),penUpts(:,2),0:(period-1));
    Upts=unique([VOI2-(numOsc-1)*period,STATES2(:,1)*1e3],'rows','stable');
    dbleup=diff(Upts(:,1));
    if sum(dbleup==0)>0
        Upts(dbleup==0,:)=[];
    end
    ult=interp1(Upts(:,1),Upts(:,2),0:(period-1));
    oscDiff=sum((ult-penUlt).^2/range(ult).^2);
end
if redoNum>=maxRedo
    warning(strcat('Not great clarity at base'));
    unstable(2*numRuns+1)=1;
end
lastTransient2=find(VOI2>=(VOI2(end)-1*period)&VOI2<=(VOI2(end)-0.8*period));
tvec2=VOI2(lastTransient2)-(VOI2(end)-period);
% 3rd case
[VOI3, STATES3, ALGEBRAIC3, CONSTANTS3] = ...
    IncDelayModelNumbers_PS(params3,[]);
period=CONSTANTS3(:,7);
numOsc=round(VOI3(end)/period);
penUpts=unique([VOI3-(numOsc-2)*period,STATES3(:,1)*1e3],'rows','stable');     
dbleup=diff(penUpts(:,1));     
if sum(dbleup==0)>0         
    penUpts(dbleup==0,:)=[];     
end
penUlt=interp1(penUpts(:,1),penUpts(:,2),0:(period-1));
Upts=unique([VOI3-(numOsc-1)*period,STATES3(:,1)*1e3],'rows','stable');
dbleup=diff(Upts(:,1));
if sum(dbleup==0)>0
    Upts(dbleup==0,:)=[];
end
ult=interp1(Upts(:,1),Upts(:,2),0:(period-1));
oscDiff=sum((ult-penUlt).^2/range(ult).^2);
redoNum=0;
while oscDiff>1e-2 & redoNum<maxRedo
    [VOI3, STATES3, ALGEBRAIC3, CONSTANTS3] = ...
    IncDelayModelNumbers_PS(params3,STATES3(end,:));
    redoNum=redoNum+1;
    penUpts=unique([VOI3-(numOsc-2)*period,STATES3(:,1)*1e3],'rows','stable');     
    dbleup=diff(penUpts(:,1));     
    if sum(dbleup==0)>0         
        penUpts(dbleup==0,:)=[];     
    end
    penUlt=interp1(penUpts(:,1),penUpts(:,2),0:(period-1));
    Upts=unique([VOI3-(numOsc-1)*period,STATES3(:,1)*1e3],'rows','stable');
    dbleup=diff(Upts(:,1));
    if sum(dbleup==0)>0
        Upts(dbleup==0,:)=[];
    end
    ult=interp1(Upts(:,1),Upts(:,2),0:(period-1));
    oscDiff=sum((ult-penUlt).^2/range(ult).^2);
end
if redoNum>=maxRedo
    warning(strcat('Not great clarity at base'));
    unstable(2*numRuns+1)=1;
end
lastTransient3=find(VOI3>=(VOI3(end)-1*period)&VOI3<=(VOI3(end)-0.8*period));
tvec3=VOI3(lastTransient3)-(VOI3(end)-period);

% 4th case
if ~isempty(params4)
[VOI4, STATES4, ALGEBRAIC4, CONSTANTS4] = ...
    IncDelayModelNumbers_PS(params4,[]);
period=CONSTANTS4(:,7);
numOsc=round(VOI4(end)/period);
penUpts=unique([VOI4-(numOsc-2)*period,STATES4(:,1)*1e3],'rows','stable');     
dbleup=diff(penUpts(:,1));     
if sum(dbleup==0)>0         
    penUpts(dbleup==0,:)=[];     
end
penUlt=interp1(penUpts(:,1),penUpts(:,2),0:(period-1));
Upts=unique([VOI4-(numOsc-1)*period,STATES4(:,1)*1e3],'rows','stable');
dbleup=diff(Upts(:,1));
if sum(dbleup==0)>0
    Upts(dbleup==0,:)=[];
end
ult=interp1(Upts(:,1),Upts(:,2),0:(period-1));
oscDiff=sum((ult-penUlt).^2/range(ult).^2);
redoNum=0;
while oscDiff>1e-2 & redoNum<maxRedo
    [VOI4, STATES4, ALGEBRAIC4, CONSTANTS4] = ...
    IncDelayModelNumbers_PS(params4,STATES4(end,:));
    redoNum=redoNum+1;
    penUpts=unique([VOI4-(numOsc-2)*period,STATES4(:,1)*1e3],'rows','stable');     
    dbleup=diff(penUpts(:,1));     
    if sum(dbleup==0)>0         
        penUpts(dbleup==0,:)=[];     
    end
    penUlt=interp1(penUpts(:,1),penUpts(:,2),0:(period-1));
    Upts=unique([VOI4-(numOsc-1)*period,STATES4(:,1)*1e3],'rows','stable');
    dbleup=diff(Upts(:,1));
    if sum(dbleup==0)>0
        Upts(dbleup==0,:)=[];
    end
    ult=interp1(Upts(:,1),Upts(:,2),0:(period-1));
    oscDiff=sum((ult-penUlt).^2/range(ult).^2);
end
if redoNum>=maxRedo
    warning(strcat('Not great clarity at base'));
    unstable(2*numRuns+1)=1;
end
lastTransient4=find(VOI4>=(VOI4(end)-1*period)&VOI4<=(VOI4(end)-0.8*period));
tvec4=VOI4(lastTransient4)-(VOI4(end)-period);
end
end
%% Plot Fig 5

titletext={'Cytosol','SR','I_{RyR}','I_{SERCA}','I_{IP_3R}','I_{NCX}','I_{LCC}','I_{SR leak}','I_{CaB}','I_{pCa}'};
ydesc={'[Ca^{2+}] (\muM)','[Ca^{2+}] (\muM)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)'};
xdesc='Time (ms)';

valuesb={STATESb(lastTransientb,1),STATESb(lastTransientb,2),ALGEBRAICb(lastTransientb,35),-ALGEBRAICb(lastTransientb,48),ALGEBRAICb(lastTransientb,66),ALGEBRAICb(lastTransientb,46),...
    ALGEBRAICb(lastTransientb,44),ALGEBRAICb(lastTransientb,54),ALGEBRAICb(lastTransientb,53),ALGEBRAICb(lastTransientb,50)};
values1={STATES1(lastTransient1,1),STATES1(lastTransient1,2),ALGEBRAIC1(lastTransient1,35),-ALGEBRAIC1(lastTransient1,48),ALGEBRAIC1(lastTransient1,66),ALGEBRAIC1(lastTransient1,46),...
    ALGEBRAIC1(lastTransient1,44),ALGEBRAIC1(lastTransient1,54),ALGEBRAIC1(lastTransient1,53),ALGEBRAIC1(lastTransient1,50)};
values2={STATES2(lastTransient2,1),STATES2(lastTransient2,2),ALGEBRAIC2(lastTransient2,35),-ALGEBRAIC2(lastTransient2,48),ALGEBRAIC2(lastTransient2,66),ALGEBRAIC2(lastTransient2,46),...
    ALGEBRAIC2(lastTransient2,44),ALGEBRAIC2(lastTransient2,54),ALGEBRAIC2(lastTransient2,53),ALGEBRAIC2(lastTransient2,50)};
values3={STATES3(lastTransient3,1),STATES3(lastTransient3,2),ALGEBRAIC3(lastTransient3,35),-ALGEBRAIC3(lastTransient3,48),ALGEBRAIC3(lastTransient3,66),ALGEBRAIC3(lastTransient3,46),...
    ALGEBRAIC3(lastTransient3,44),ALGEBRAIC3(lastTransient3,54),ALGEBRAIC3(lastTransient3,53),ALGEBRAIC3(lastTransient3,50)};

% Convert to uM
valuesb=cellfun(@(x)1000*x,valuesb,'UniformOutput',0);
values1=cellfun(@(x)1000*x,values1,'UniformOutput',0);
values2=cellfun(@(x)1000*x,values2,'UniformOutput',0);
values3=cellfun(@(x)1000*x,values3,'UniformOutput',0);

if ~isempty(params4)
    values4={STATES4(lastTransient4,1),STATES4(lastTransient4,2),ALGEBRAIC4(lastTransient4,35),-ALGEBRAIC4(lastTransient4,48),ALGEBRAIC4(lastTransient4,66),ALGEBRAIC4(lastTransient4,46),...
        ALGEBRAIC4(lastTransient4,44),ALGEBRAIC4(lastTransient4,54),ALGEBRAIC4(lastTransient4,53),ALGEBRAIC4(lastTransient4,50)};
    values4=cellfun(@(x)1000*x,values4,'UniformOutput',0);
end

if 1
numa=2; numd=3;
figure
for ii=[1,2,3:6]
   subplot(numd,numa,ii)
   hold on
   plot(tvecb,valuesb{ii},'-','LineWidth',2)
   plot(tvec3,values3{ii},'-','LineWidth',2)
   plot(tvec2,values2{ii},'-','LineWidth',2)
   if ~isempty(params4)
       plot(tvec4,values4{ii},'-','LineWidth',2)
   else
       plot(tvec1,values1{ii},'-','LineWidth',2,'HandleVisibility','off')
   end
   plot(tvec1,values1{ii},'-','LineWidth',2)
%    xlim([0 40])
   title(titletext{ii})
   if ii>numa*numd-numa
       xlabel(xdesc)
   end
   if mod(ii,2)==1
       ylabel(ydesc{ii})
   end
   if ~isempty(params4)
   axis([0 200 min([0;values3{ii};valuesb{ii};values1{ii};values2{ii};values4{ii}]) max([values4{ii};valuesb{ii};values1{ii};values2{ii}])])
   end
end
% subplot(numd,numa,2)
% hold on
% plot(tvalc(tvalclr),valuesc{1}(tvalclr),'-x','LineWidth',0.5)
% plot(tval(tvallr),values{1}(tvallr),'-x','LineWidth',0.5)
% title('Cytosol (16Hz)')
% axis([0 CONSTANTS(:,7) min([values{1};valuesc{1}]) max([values{1};valuesc{1}])])

subplot(numd,numa,3)
if isempty(params4)
    legstr={'No IP3','K_c=6.4\muM','K_c=8\muM','IP3R type I parameters'};
else
    legstr={'No IP3','K_c=6.4\muM','K_c=8\muM','K_c=10\muM & k_f=2.4\mum^3/ms','IP3R type I parameters'};
end
% legend('No IP3','K_c=8\muM','K_c=6.4\muM','IP3R type I parameters')
% legend(legstr)
if 1
saveas(gcf,strcat(savepath,'fluxhip',strparams,'.eps'),'epsc')
saveas(gcf,strcat(savepath,'fluxhip',strparams,'.fig'))
end
end

%% Plot Fig 5 extended to all fluxes

titletext={'Cytosol','SR','I_{RyR}','I_{SERCA}','I_{IP_3R}','I_{NCX}','I_{LCC}','I_{SR}','I_{CaB}','I_{pCa}','I_{TRPN}','[CaTRPN]'};
ydesc={'[Ca^{2+}] (\muM)','[Ca^{2+}] (\muM)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)','(\muM)'};
xdesc='Time (ms)';

valuesb={STATESb(lastTransientb,1),STATESb(lastTransientb,2),ALGEBRAICb(lastTransientb,35),-ALGEBRAICb(lastTransientb,48),ALGEBRAICb(lastTransientb,66),ALGEBRAICb(lastTransientb,46),...
    ALGEBRAICb(lastTransientb,44),ALGEBRAICb(lastTransientb,54),ALGEBRAICb(lastTransientb,53),ALGEBRAICb(lastTransientb,50),ALGEBRAICb(lastTransientb,55),STATESb(lastTransientb,6)};
values1={STATES1(lastTransient1,1),STATES1(lastTransient1,2),ALGEBRAIC1(lastTransient1,35),-ALGEBRAIC1(lastTransient1,48),ALGEBRAIC1(lastTransient1,66),ALGEBRAIC1(lastTransient1,46),...
    ALGEBRAIC1(lastTransient1,44),ALGEBRAIC1(lastTransient1,54),ALGEBRAIC1(lastTransient1,53),ALGEBRAIC1(lastTransient1,50),ALGEBRAIC1(lastTransient1,55),STATES1(lastTransient1,6)};
values2={STATES2(lastTransient2,1),STATES2(lastTransient2,2),ALGEBRAIC2(lastTransient2,35),-ALGEBRAIC2(lastTransient2,48),ALGEBRAIC2(lastTransient2,66),ALGEBRAIC2(lastTransient2,46),...
    ALGEBRAIC2(lastTransient2,44),ALGEBRAIC2(lastTransient2,54),ALGEBRAIC2(lastTransient2,53),ALGEBRAIC2(lastTransient2,50),ALGEBRAIC2(lastTransient2,55),STATES2(lastTransient2,6)};
values3={STATES3(lastTransient3,1),STATES3(lastTransient3,2),ALGEBRAIC3(lastTransient3,35),-ALGEBRAIC3(lastTransient3,48),ALGEBRAIC3(lastTransient3,66),ALGEBRAIC3(lastTransient3,46),...
    ALGEBRAIC3(lastTransient3,44),ALGEBRAIC3(lastTransient3,54),ALGEBRAIC3(lastTransient3,53),ALGEBRAIC3(lastTransient3,50),ALGEBRAIC3(lastTransient3,55),STATES3(lastTransient3,6)};

% Convert to uM
valuesb=cellfun(@(x)1000*x,valuesb,'UniformOutput',0);
values1=cellfun(@(x)1000*x,values1,'UniformOutput',0);
values2=cellfun(@(x)1000*x,values2,'UniformOutput',0);
values3=cellfun(@(x)1000*x,values3,'UniformOutput',0);

if ~isempty(params4)
    values4={STATES4(lastTransient4,1),STATES4(lastTransient4,2),ALGEBRAIC4(lastTransient4,35),-ALGEBRAIC4(lastTransient4,48),ALGEBRAIC4(lastTransient4,66),ALGEBRAIC4(lastTransient4,46),...
        ALGEBRAIC4(lastTransient4,44),ALGEBRAIC4(lastTransient4,54),ALGEBRAIC4(lastTransient4,53),ALGEBRAIC4(lastTransient4,50)};
    values4=cellfun(@(x)1000*x,values4,'UniformOutput',0);
end

if 1
numa=2; numd=6;
figure('Position',[10,10,500,700])
for ii=1:size(titletext,2)
   subplot(numd,numa,ii)
   hold on
   plot(tvecb,valuesb{ii},'-','LineWidth',2)
   plot(tvec3,values3{ii},'-','LineWidth',2)
   plot(tvec2,values2{ii},'-','LineWidth',2)
   if ~isempty(params4)
       plot(tvec4,values4{ii},'-','LineWidth',2)
   else
       plot(tvec1,values1{ii},'-','LineWidth',2,'HandleVisibility','off')
   end
   plot(tvec1,values1{ii},'-','LineWidth',2)
%    xlim([0 40])
   title(titletext{ii})
   if ii>numa*numd-numa
       xlabel(xdesc)
   end
   if mod(ii,2)==1
       ylabel(ydesc{ii})
   end
   if ~isempty(params4)
   axis([0 200 min([0;values3{ii};valuesb{ii};values1{ii};values2{ii};values4{ii}]) max([values4{ii};valuesb{ii};values1{ii};values2{ii}])])
   end
end
ylabel(ydesc{ii})
% subplot(numd,numa,2)
% hold on
% plot(tvalc(tvalclr),valuesc{1}(tvalclr),'-x','LineWidth',0.5)
% plot(tval(tvallr),values{1}(tvallr),'-x','LineWidth',0.5)
% title('Cytosol (16Hz)')
% axis([0 CONSTANTS(:,7) min([values{1};valuesc{1}]) max([values{1};valuesc{1}])])

subplot(numd,numa,3)
if isempty(params4)
    legstr={'No IP3','K_c=6.4\muM','K_c=8\muM','IP3R type I parameters'};
else
    legstr={'No IP3','K_c=6.4\muM','K_c=8\muM','K_c=10\muM & k_f=2.4\mum^3/ms','IP3R type I parameters'};
end
% legend('No IP3','K_c=8\muM','K_c=6.4\muM','IP3R type I parameters')
% legend(legstr)
if 1
saveas(gcf,strcat(savepath,'fluxhip_all',strparams,'.eps'),'epsc')
saveas(gcf,strcat(savepath,'fluxhip_all',strparams,'.fig'))
end
end

%% New figure showing difference between Hinch and Hinch + Sneyd
% Calculate model
% 5th case
params5=[1,1,1,1,1,0.1,1,1,1,1];
[VOI5, STATES5, ALGEBRAIC5, CONSTANTS5] = ...
    IncDelayModelNumbers_PS(params5,[]);
period=CONSTANTS5(:,7);
numOsc=round(VOI5(end)/period);
penUpts=unique([VOI5-(numOsc-2)*period,STATES5(:,1)*1e3],'rows','stable');     
dbleup=diff(penUpts(:,1));     
if sum(dbleup==0)>0         
    penUpts(dbleup==0,:)=[];     
end
penUlt=interp1(penUpts(:,1),penUpts(:,2),0:(period-1));
Upts=unique([VOI5-(numOsc-1)*period,STATES5(:,1)*1e3],'rows','stable');
dbleup=diff(Upts(:,1));
if sum(dbleup==0)>0
    Upts(dbleup==0,:)=[];
end
ult=interp1(Upts(:,1),Upts(:,2),0:(period-1));
oscDiff=sum((ult-penUlt).^2/range(ult).^2);
redoNum=0;
while oscDiff>1e-2 & redoNum<maxRedo
    [VOI5, STATES5, ALGEBRAIC5, CONSTANTS5] = ...
    IncDelayModelNumbers_PS(params5,STATES5(end,:));
    redoNum=redoNum+1;
    penUpts=unique([VOI5-(numOsc-2)*period,STATES5(:,1)*1e3],'rows','stable');     
    dbleup=diff(penUpts(:,1));     
    if sum(dbleup==0)>0         
        penUpts(dbleup==0,:)=[];     
    end
    penUlt=interp1(penUpts(:,1),penUpts(:,2),0:(period-1));
    Upts=unique([VOI5-(numOsc-1)*period,STATES5(:,1)*1e3],'rows','stable');
    dbleup=diff(Upts(:,1));
    if sum(dbleup==0)>0
        Upts(dbleup==0,:)=[];
    end
    ult=interp1(Upts(:,1),Upts(:,2),0:(period-1));
    oscDiff=sum((ult-penUlt).^2/range(ult).^2);
end
if redoNum>=maxRedo
    warning(strcat('Not great clarity at base'));
    unstable(2*numRuns+1)=1;
end
lastTransient5=find(VOI5>=(VOI5(end)-1*period)&VOI5<=(VOI5(end)-0*period));
tvec5=VOI5(lastTransient5)-(VOI5(end)-period);
 
titletext={'Cytosol','SR','I_{RyR}','I_{SERCA}','I_{IP_3R}','I_{NCX}','I_{LCC}','I_{SR leak}','I_{CaB}','I_{pCa}','I_{TRPN}','[CaTRPN]'};
ydesc={'[Ca^{2+}] (\muM)','[Ca^{2+}] (\muM)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)','(\muM)'};
xdesc='Time (ms)';
% 
lastTransientb=find(VOIb>=(VOIb(end)-1*period)&VOIb<=(VOIb(end)-0*period));
tvecb=VOIb(lastTransientb)-(VOIb(end)-period);
valuesb={STATESb(lastTransientb,1),STATESb(lastTransientb,2),ALGEBRAICb(lastTransientb,35),-ALGEBRAICb(lastTransientb,48),ALGEBRAICb(lastTransientb,66),ALGEBRAICb(lastTransientb,46),...
    ALGEBRAICb(lastTransientb,44),ALGEBRAICb(lastTransientb,54),ALGEBRAICb(lastTransientb,53),ALGEBRAICb(lastTransientb,50),ALGEBRAICb(lastTransientb,55),STATESb(lastTransientb,6)};
values5={STATES5(lastTransient5,1),STATES5(lastTransient5,2),ALGEBRAIC5(lastTransient5,35),-ALGEBRAIC5(lastTransient5,48),ALGEBRAIC5(lastTransient5,66),ALGEBRAIC5(lastTransient5,46),...
    ALGEBRAIC5(lastTransient5,44),ALGEBRAIC5(lastTransient5,54),ALGEBRAIC5(lastTransient5,53),ALGEBRAIC5(lastTransient5,50),ALGEBRAIC5(lastTransient5,55),STATES5(lastTransient5,6)};
valuesb=cellfun(@(x)1000*x,valuesb,'UniformOutput',0);
values5=cellfun(@(x)1000*x,values5,'UniformOutput',0);
%%
figure('Position',[10,10,1200,200])
numa=4;
numd=1;
jj=[1,2,3,5];
for ii=1:4
    subplot(numd,numa,ii)
    
set(gca,'FontSize',16)
   hold on
   plot(tvecb,valuesb{jj(ii)},'-','LineWidth',2)
   plot(tvec5,values5{jj(ii)},'-','LineWidth',2)
   title(titletext{jj(ii)})
   if ii>numa*numd-numa
       xlabel(xdesc)
   end
   if mod(ii,2)==1
       ylabel(ydesc{ii})
   end
%    if ~isempty(params4)
%    axis([0 200 min([0;values3{ii};valuesb{ii};values1{ii};values2{ii};values4{ii}]) max([values4{ii};valuesb{ii};values1{ii};values2{ii}])])
%    end
if jj(ii)==1
    axis([0 1000 0 0.8])
elseif jj(ii)==2
axis([0 1000 0 800])

end
end
% lastTransientb=find(VOIb>=(VOIb(end)-1*period)&VOIb<=(VOIb(end)-0.8*period));
% tvecb=VOIb(lastTransientb)-(VOIb(end)-period);
if 1
saveas(gcf,strcat(savepath,'example_output_flux','.eps'),'epsc')
saveas(gcf,strcat(savepath,'example_output_flux','.fig'))
end