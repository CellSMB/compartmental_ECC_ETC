%% This is a script to generate figures for the paper Hunt, Tilunaite, 
% Bass, Roderick,Soeller,Rajagopal,Crampin 2020
% Written by H Hunt 2020
savepath='./';
%% Figure 2C
% Calculate model
% Base case
[VOIb, STATESb, ALGEBRAICb, CONSTANTSb] = ...
    ECC_ETC_PS(paramsb,[]);
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
    ECC_ETC_PS(paramsb,STATESb(end,:));
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
% Add in ETC with Sneyd et al. 2017 parameters
params5=[1,1,1,1,1,0.1,1,1,1,1];
[VOI5, STATES5, ALGEBRAIC5, CONSTANTS5] = ...
    ECC_ETC_PS(params5,[]);
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
    ECC_ETC_PS(params5,STATES5(end,:));
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

lastTransientb=find(VOIb>=(VOIb(end)-1*period)&VOIb<=(VOIb(end)-0*period));
tvecb=VOIb(lastTransientb)-(VOIb(end)-period); 

titletext={'Cytosol','SR','I_{RyR}','I_{SERCA}','I_{IP_3R}','I_{NCX}','I_{LCC}','I_{SR leak}','I_{CaB}','I_{pCa}','I_{TRPN}','[CaTRPN]'};
ydesc={'[Ca^{2+}] (\muM)','[Ca^{2+}] (\muM)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)','Flux (\muM/ms)','(\muM)'};
xdesc='Time (ms)';
% 
valuesb={STATESb(lastTransientb,1),STATESb(lastTransientb,2),ALGEBRAICb(lastTransientb,35),-ALGEBRAICb(lastTransientb,48),ALGEBRAICb(lastTransientb,66),ALGEBRAICb(lastTransientb,46),...
    ALGEBRAICb(lastTransientb,44),ALGEBRAICb(lastTransientb,54),ALGEBRAICb(lastTransientb,53),ALGEBRAICb(lastTransientb,50),ALGEBRAICb(lastTransientb,55),STATESb(lastTransientb,6)};
values5={STATES5(lastTransient5,1),STATES5(lastTransient5,2),ALGEBRAIC5(lastTransient5,35),-ALGEBRAIC5(lastTransient5,48),ALGEBRAIC5(lastTransient5,66),ALGEBRAIC5(lastTransient5,46),...
    ALGEBRAIC5(lastTransient5,44),ALGEBRAIC5(lastTransient5,54),ALGEBRAIC5(lastTransient5,53),ALGEBRAIC5(lastTransient5,50),ALGEBRAIC5(lastTransient5,55),STATES5(lastTransient5,6)};
valuesb=cellfun(@(x)1000*x,valuesb,'UniformOutput',0);
values5=cellfun(@(x)1000*x,values5,'UniformOutput',0);
% Plot results
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
if jj(ii)==1
    axis([0 1000 0 0.8])
elseif jj(ii)==2
axis([0 1000 0 800])

end
end
if 1
saveas(gcf,strcat(savepath,'example_output_flux','.eps'),'epsc')
saveas(gcf,strcat(savepath,'example_output_flux','.fig'))
end
%% Figure 3
% Plot approximate IP3R gating as contour map
k=0.324;K_p=2e-4;
K_c=(1:5:100)*2e-4;
K_h=[1 20 27 60]*8e-5;

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
% Plot figure
intvar=P_open;
figure
for ii=1:alltog
subplot(size(K_h,2),3,ii)
contourf(ca*1e3,K_c*1e3,intvar{ii})
xlabel('Ca^{2+} (\muM)')
ylabel('K_c (\muM)')
zlabel('log(P_O)')
colorbar
caxis([0 0.5])
end
%% Heatmaps in Figure 4, 6, 7, 8, S2
% Show effects of IP3/k_f vs K_c on transient qualities
% Change this to 0 to vary maximum IP3R flux (k_f) or 1 to vary IP3
% concentration.
isip3=1;
period=1000;
% Run the Calcium_dynamics_ECC_ETC_cytosol_PS to generate the table of
% simulations
results = Calcium_dynamics_ECC_ETC_cytosol_PS(isip3,period);

% Decide what to call outputs
ti = {'Amplitude','FDHM','FD at 90% of max','Time to peak','Baseline','Duty cycle',... %5+6=11
    'AmplitudeF','FDHMF','FD at 90% of maxF','Time to peak F','Duty cycle F',... %5+11=16
    'IP_3R flux','RyR flux','Diastolic [Ca^{2+}]_i','Diastolic [Ca^{2+}]_{SR}'};%20

yl={'(mM)','(ms)','(ms)','(ms)','(mM)','','(mM)','(ms)','(ms)','(ms)','','','','mM','mM'};
numVars = 5;
if isip3
    xl={'[IP_3]','t_{max}','K_c','K_h','K_t'};
else
    xl={'k_f','t_{max}','K_c','K_h','K_t'};
end

% Plot results as heat maps
ii=1;jj=3;
numRuns=400;
resrange=401:800;
isip3=1;
if isip3
resrange1=find(results(resrange,1)<=150)+numRuns;
resrange2=find(results(resrange,3)<=60)+numRuns;
resrange=intersect(resrange1,resrange2);
else
    resrange1=find(results(resrange,1)<=100)+numRuns;
    resrange2=find(results(resrange,3)<=100)+numRuns;
    resrange=intersect(resrange1,resrange2);
end
quals=[numVars+1,numVars+2,numVars+5,numVars+6];
for qual = quals
    dc=results(resrange,qual)/results(end,qual);
    dc=dc*100-100;
    % Threshold highest shown % change at ~160%. Higher values become
    % physiologically implausible
    dc(dc>160)=159.9; 
    %switch between IP3 and k_f on axis
    if isip3
        f = fit([results(resrange,(ii))*10, 2e-1*results(resrange,(jj))],dc,'linearinterp');
    else
        f = fit([results(resrange,(ii))*3e-2, 2e-1*results(resrange,(jj))],dc,'linearinterp');
    end
    
    if ismember(qual,[numVars+14,numVars+15])
        figure('position',[100,200,560,450]*1.04)
    else
        figure('position',[0,200,560,450])
    end
    plot(f','Style','contour')
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
    colormap(bluewhitered)
if isip3
    axis([0 14 0 10])
    xlabel(strcat(xl(ii),{' '},'(\muM)'),'FontSize',20)
else
    axis([0 2.8 0 10])
    xlabel(strcat(xl(ii),{' '},'(\mum^3/ms)'),'FontSize',20)
end
end

%% Plot of different parameter regions in Figure 4, 6, 7, 8, S2
% This code requires the same simulation results as the heatmaps above.
isip3=1;
period=1000;
% Run the Calcium_dynamics_ECC_ETC_cytosol_PS to generate the table of
% simulations
results = Calcium_dynamics_ECC_ETC_cytosol_PS(isip3,period);
base=results(401,:);

bg=401;
en=800;

if isip3
    ip3=reshape(results(bg:en,1),[20,20])*10;
else
    ip3=reshape(results(bg:en,1),[20,20])*3e-2;
end
Kc=reshape(results(bg:en,3),[20,20])*0.2;
Kh=reshape(results(bg:en,4),[20,20])*0.08;

amplitude=reshape((results(bg:en,6)-base(6))/base(6),[20,20]).*100;
fdhm=reshape((results(bg:en,7)-base(7))/base(7),[20,20]).*100;
diastolicca=reshape((results(bg:en,10)-base(10))/base(10),[20,20]).*100;
dc=reshape((results(bg:en,11)-base(11))/base(11),[20,20]).*100;
diastoliccasr=reshape((results(bg:en,20)-base(20))/base(20),[20,20]).*100;

figure;
% dc is always >=0, so we code remaining 3 parameters as rgb
rgb_increase_regions=zeros(20,20,3); %initialize
rgb_increase_regions(:,:,1)=1*(amplitude>0);
rgb_increase_regions(:,:,2)=1*(fdhm<0)-1*(fdhm>0);
rgb_increase_regions(:,:,3)=1*(diastolicca>0);
new_rgb=rgb_increase_regions;

c_p = [123 50 144]/255;
c_or = [217 83 25]/255;
c_blu = [23 114 189]/255;
c_y = [247 177 41]/255;
c_gr = [119 172 49]/255;
for ii=1:20
    for jj=1:20
        if permute(rgb_increase_regions(ii,jj,:),[1,3,2])==[1,1,0]
            new_rgb(ii,jj,:)=permute(c_p,[1,3,2]);
        elseif permute(rgb_increase_regions(ii,jj,:),[1,3,2])==[0,1,0]
            new_rgb(ii,jj,:)=permute(c_or,[1,3,2]);
        elseif permute(rgb_increase_regions(ii,jj,:),[1,3,2])==[0,1,1]
            new_rgb(ii,jj,:)=permute(c_blu,[1,3,2]);
        elseif permute(rgb_increase_regions(ii,jj,:),[1,3,2])==[0,0,0]
            new_rgb(ii,jj,:)=permute(c_blu,[1,3,2]);
        elseif permute(rgb_increase_regions(ii,jj,:),[1,3,2])==[0,-1,0]
            new_rgb(ii,jj,:)=permute(c_y,[1,3,2]);
        elseif permute(rgb_increase_regions(ii,jj,:),[1,3,2])==[0,-1,1]
            new_rgb(ii,jj,:)=permute(c_gr,[1,3,2]);
            elseif permute(rgb_increase_regions(ii,jj,:),[1,3,2])==[1,-1,0]
            new_rgb(ii,jj,:)=permute(c_or,[1,3,2]);
        else
            warning('t')
        end
    end
end
imagesc([Kc(1,1), Kc(end,end)],[ip3(1,1),ip3(end,end)],new_rgb)
    set(gca,'FontSize',16)
    axis('square')

view(0,-90)
xlabel('K_c (\mu M)')
if isip3
    ylabel('IP_3 (\mu M)','FontSize',20)
else
    ylabel('k_f (\mu m^3 /ms)','FontSize',20)
end
title('Change regions','FontSize',24)
%% Figure 5, S3
% Generate sims
% isip3=1 for Figure 5 and 0 for figure S3
isip3=1;
if isip3
    strparams='varyingIP3';
    kf=1.5;
    paramsb=[0,1,1,1,1,0,1,1,1,1,3e3];
    params1=[1,1,1,1,1,kf,1,1,1,1,3e3];
    params2=[1,10,28,6,10,kf,1,1,1,1,3e3];
    params3=[1,10,40,6,10,kf,1,1,1,1,3e3];
    params4=[];
else 
    strparams='varyingkf';
    paramsb=[0,1,1,1,1,0,1,1,1,1,3e3];
    params1=[1,1,1,1,1,15,1,1,1,1,3e3];
    params2=[1,10,32,6,10,15,1,1,1,1,3e3];
    params3=[1,10,40,6,10,15,1,1,1,1,3e3];
    params4=[1,10,48,6,10,80,1,1,1,1];
end
maxRedo=1e3;
% Base case
[VOIb, STATESb, ALGEBRAICb, CONSTANTSb] = ...
    ECC_ETC_PS(paramsb,[]);
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
    ECC_ETC_PS(paramsb,STATESb(end,:));
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
    ECC_ETC_PS(params1,[]);
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
    ECC_ETC_PS(params1,STATES1(end,:));
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
    ECC_ETC_PS(params2,[]);
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
    ECC_ETC_PS(params2,STATES2(end,:));
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
    ECC_ETC_PS(params3,[]);
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
    ECC_ETC_PS(params3,STATES3(end,:));
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
    ECC_ETC_PS(params4,[]);
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
    ECC_ETC_PS(params4,STATES4(end,:));
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

% Plot transients
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
   if ii==1
       axis([0 200 0.08 0.7])
   end
end

subplot(numd,numa,3)
if isempty(params4)
    legstr={'No IP3','K_c=6.4\muM','K_c=8\muM','IP3R type I parameters'};
else
    legstr={'No IP3','K_c=6.4\muM','K_c=8\muM','K_c=10\muM & k_f=2.4\mum^3/ms','IP3R type I parameters'};
end
if 1
saveas(gcf,strcat(savepath,'fluxhip',strparams,'.eps'),'epsc')
saveas(gcf,strcat(savepath,'fluxhip',strparams,'.fig'))
end
end

%% Figure 9
namedesc='NFAT2';
ti = {'Dephosphorylated NFAT_n','Phosphorylated NFAT_n','Dephosphorylated NFAT_c','Phosphorylated NFAT_c'};%20
numVars = 3;
numRuns=(size(results,1)-1);
xl={'IP_3','K_c','k_f'};
% isip3 should be 1 for Figure 9A and 0 for Figure 9B
isip3=1;
results = NFAT_PS(isip3);

ii=1;jj=2;
numRuns=size(results,1);
resrange=1:numRuns;

resrange1=find(results(resrange,1)<=2);
kfval=15;
resrange=intersect(find(results(resrange,3)==kfval),resrange1);
quals=[numVars+1];
for qual = quals
    dc=results(resrange,qual)/results(1,qual);
    dc=dc*100-100;
    if isip3
        f = fit([results(resrange,(ii))*10, 2e-1*results(resrange,(jj))],dc,'linearinterp');
    else
        f = fit([results(resrange,(ii))*3e-2, 2e-1*results(resrange,(jj))],dc,'linearinterp');
    end
    
    figure('position',[100,200,560,450]*1.04)
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
    cb.Title.Position(2) = 352;
    cb.Title.Position(1) = 50;
    title(strcat(ti{qual-numVars},{' '},'(%)'),'FontSize',25)
    ylabel(strcat(xl(jj),{' '},'(\muM)'),'FontSize',20)
    caxis([-60 160])
colormap(bluewhitered)
if isip3
    axis([0 14 1e-1 10])
    xlabel(strcat(xl(ii),{' '},'(\muM)'),'FontSize',20)
      saveas(gcf,strcat(savepath,'contour_ip3_K_c_pc',...
    ti{qual-numVars},namedesc,'.eps'),'epsc')
saveas(gcf,strcat(savepath,'contour_change_ip3_K_c_kf-',num2str(kfval),namedesc,'.fig'))
else
    axis([0 2.8 0 10])
    xlabel(strcat(xl(ii),{' '},'(\mum^3/ms)'),'FontSize',20)
    saveas(gcf,strcat(savepath,'contour_kf_K_c_',...
    ti{qual-numVars},namedesc,'.eps'),'epsc')
saveas(gcf,strcat(savepath,'contour_change_kf_K_c_',...
    ti{qual-numVars},namedesc,'.fig'),'epsc')
end
end

%% Figure S1
% Change the beat period in ECC_ETC_PS (CONSTANTS(:,7)) to
% change the amount of time simulated.
paramsb=[0,1,1,1,1,0,1,1,1,1];
[VOIb, STATESb, ALGEBRAICb, CONSTANTSb] = ...
    ECC_ETC_PS(paramsb,[]);
figure
plot(VOI,ALGEBRAIC(:,1),'LineWidth',2)
yaxis('V (mV)')
xaxis('time (ms)')
%% Figure S4
% Run the simulations for Figure 5
% Then generate this figure
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
subplot(numd,numa,3)
if isempty(params4)
    legstr={'No IP3','K_c=6.4\muM','K_c=8\muM','IP3R type I parameters'};
else
    legstr={'No IP3','K_c=6.4\muM','K_c=8\muM','K_c=10\muM & k_f=2.4\mum^3/ms','IP3R type I parameters'};
end
if 1
saveas(gcf,strcat(savepath,'fluxhip_all',strparams,'.eps'),'epsc')
saveas(gcf,strcat(savepath,'fluxhip_all',strparams,'.fig'))
end
end
