function results = Calcium_dynamics_ECC_ETC_cytosol_PS(isip3,period)
numparRuns = 400;
% Choose parameters
if isip3
    fluxvar     =1.5*repmat(0:(1/19):1,[1,20]);
else
    fluxvar     =100*repmat(0:(1/19):1,[1,20]);
end
t_max       = 10*ones(numparRuns,1);
K_c         = reshape(repmat(0:(50/19):50,[20,1]),400,1);%50*rand(numparRuns,1);
K_h         = 6*ones(numparRuns,1);
K_t         = 10*ones(numparRuns,1);
numRuns=numparRuns;
numStates=7;
totResults=18+numStates;
totalRuns = numRuns;
% Make room for results
resultsD=zeros(totalRuns,totResults);
resultsnD=zeros(totalRuns,totResults);
status = mkdir('/home/hhunt1/MATLAB/Graphs');
unstable=zeros(1,2*numRuns+1);
for it=1:numRuns
    params=[fluxvar(it),t_max(it),K_c(it),K_h(it),K_t(it)];
    if isip3
        [VOI, STATES, ALGEBRAIC, CONSTANTS] = ...
            ECC_ETC_ip3(params,period,[]);
    else
        [VOI, STATES, ALGEBRAIC, CONSTANTS] = ...
            ECC_ETC_kf(params,period,[]);
    end
    period=CONSTANTS(:,7);
    % stability
    numOsc=round(VOI(end)/period);
    penUlt=interp1(VOI-(numOsc-2)*period,STATES(:,1)*1e3,0:(period-1));
    ult=interp1(VOI-(numOsc-1)*period,STATES(:,1)*1e3,0:(period-1));
    oscDiff=sum((ult-penUlt).^2/range(ult).^2);
    redoNum=0;
    maxRedo=50;
    while oscDiff>1e-3 && redoNum<maxRedo
        if isip3
        [VOI, STATES, ALGEBRAIC, CONSTANTS] = ...
            ECC_ETC_ip3(params,period,STATES(end,:));
    else
        [VOI, STATES, ALGEBRAIC, CONSTANTS] = ...
            ECC_ETC_kf(params,period,STATES(end,:);
    end
        redoNum=redoNum+1;
        penUlt=interp1(VOI-(numOsc-2)*period,STATES(:,1)*1e3,0:(period-1));
        ult=interp1(VOI-(numOsc-1)*period,STATES(:,1)*1e3,0:(period-1));
        oscDiff=sum((ult-penUlt).^2/range(ult).^2);
    end
    if redoNum>=maxRedo
        warning(strcat('Not great clarity at it=',num2str(it)));
        unstable(it)=1;
    end
    %end stability
    lastTransient=find(VOI>=(VOI(end)-period));
	resultsD(it,:)=getQuals(VOI, STATES, ALGEBRAIC, CONSTANTS,params,lastTransient);
    params=[fluxvar(it),1,K_c(it),K_h(it),1];
    if isip3
        [VOI2, STATES2, ALGEBRAIC2, CONSTANTS2] = ...
            ECC_ETC_ip3(params,period,[]);
    else
        [VOI2, STATES2, ALGEBRAIC2, CONSTANTS2] = ...
            ECC_ETC_kf(params,period,[]);
    end
    % stability
    numOsc=round(VOI2(end)/period);
    penUlt=interp1(VOI2-(numOsc-2)*period,STATES2(:,1)*1e3,0:(period-1));
    ult=interp1(VOI2-(numOsc-1)*period,STATES2(:,1)*1e3,0:(period-1));
    oscDiff=sum((ult-penUlt).^2/range(ult).^2);
    redoNum=0;
    maxRedo=50;
    while oscDiff>1e-3 && redoNum<maxRedo
        if isip3
            [VOI2, STATES2, ALGEBRAIC2, CONSTANTS2] = ...
                ECC_ETC_ip3(params,period,STATES2(end,:));
        else
            [VOI2, STATES2, ALGEBRAIC2, CONSTANTS2] = ...
                ECC_ETC_kf(params,period,STATES2(end,:));
        end
        redoNum=redoNum+1;
        penUlt=interp1(VOI2-(numOsc-2)*period,STATES2(:,1)*1e3,0:(period-1));
        ult=interp1(VOI2-(numOsc-1)*period,STATES2(:,1)*1e3,0:(period-1));
        oscDiff=sum((ult-penUlt).^2/range(ult).^2);
    end
    if redoNum>=maxRedo
        warning(strcat('Not great clarity at it=',num2str(it)));
        unstable(it)=1;
    end
    %end stability
    lastTransient2=find(VOI2>=(VOI2(end)-period));
	resultsnD(it,:)=getQuals(VOI2, STATES2, ALGEBRAIC2, CONSTANTS2,...
        params,lastTransient2);
    
    time=VOI>(VOI(end)-CONSTANTS(:,7)-500);
    time2=VOI2>(VOI2(end)-CONSTANTS2(:,7)-500);
    if 0
    figure('pos',[600,0,600,600])
    numa = 1;
    numd =3;
    subplot(numd,numa,1)
    hold on
    plot(VOI2(time2)-(VOI2(end)-CONSTANTS2(:,7)-500),STATES2(time2,1))
    plot(VOI(time)-(VOI(end)-CONSTANTS(:,7)-500),STATES(time,1))
    plot(VOI(1:(CONSTANTS(:,7)))+500,STATES(1:(CONSTANTS(:,7)),1))
    legend({'Ca_i+IP_3','Ca_i+IP_3 w delay','Ca_i'})
    title(strcat(num2str(fluxvar(it)),'-',num2str(t_max(it)),'-',...
        num2str(K_c(it)),'-',num2str(K_h(it)),'-',num2str(K_t(it))));
    hold on
    subplot(numd,numa,2)
    hold on
    plot(VOI2(time2)-(VOI2(end)-CONSTANTS2(:,7)-500),ALGEBRAIC2(time2,66))
    plot(VOI(time)-(VOI(end)-CONSTANTS(:,7)-500),ALGEBRAIC(time,66))
    legend('IP_3R','IP_3R w delay')
    subplot(numd,numa,3)
    hold on
    plot(VOI2(time2)-(VOI2(end)-CONSTANTS2(:,7)-500),ALGEBRAIC2(time2,65))
    plot(VOI(time)-(VOI(end)-CONSTANTS(:,7)-500),ALGEBRAIC(time,65))
    legend('P_0','P_0 w delay')
    saveas(gcf,strcat('/home/hhunt1/MATLAB/Graphs/ip3',num2str(fluxvar(it)),...
        '-',num2str(t_max(it)),'-',num2str(K_c(it)),'-',num2str(K_h(it)),...
        '-',num2str(K_t(it)),'.png'))
    end
end
%% Calculate base qualities
params=[0,1,1,1,1];
if isip3
    [VOI3, STATES3, ALGEBRAIC3, CONSTANTS3] = ...
        ECC_ETC_ip3(params,period,[]);
else
    [VOI3, STATES3, ALGEBRAIC3, CONSTANTS3] = ...
        ECC_ETC_kf(params,period,[]);
end
period=CONSTANTS3(:,7);
    numOsc=round(VOI3(end)/period);
    penUlt=interp1(VOI3-(numOsc-2)*period,STATES3(:,1)*1e3,0:(period-1));
    ult=interp1(VOI3-(numOsc-1)*period,STATES3(:,1)*1e3,0:(period-1));
    oscDiff=sum((ult-penUlt).^2/range(ult).^2);
    redoNum=0;
    maxRedo=50;
    while oscDiff>1e-3 && redoNum<maxRedo
        if isip3
            [VOI3, STATES3, ALGEBRAIC3, CONSTANTS3] = ...
                ECC_ETC_ip3(params,period,STATES3(end,:));
        else
            [VOI3, STATES3, ALGEBRAIC3, CONSTANTS3] = ...
                ECC_ETC_kf(params,period,STATES3(end,:));
        end
        redoNum=redoNum+1;
        penUlt=interp1(VOI-(numOsc-2)*period,STATES(:,1)*1e3,0:(period-1));
        ult=interp1(VOI-(numOsc-1)*period,STATES(:,1)*1e3,0:(period-1));
        oscDiff=sum((ult-penUlt).^2/range(ult).^2);
    end
    if redoNum>=maxRedo
        warning(strcat('Not great clarity at base'));
        unstable(2*numRuns+1)=1;
    end
    lastTransient=find(VOI3>=(VOI3(end)-period));
    baseresults=getQuals(VOI3, STATES3, ALGEBRAIC3, CONSTANTS3,params,lastTransient);
%% Save results
results = [resultsnD;resultsD;baseresults];
end
