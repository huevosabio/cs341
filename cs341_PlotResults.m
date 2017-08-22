clear all; close all; clc;

fileslist={'MPC_Perf_5000v_rev3','MPC-LSTM5000v_rev3','REACTIVE5000v_rev3','MPC-EMPTY5000v_rev3'}%,'ReactiveSim_7500v_rev2','ReactiveSim_9800v_rev2','mpcperfect'};

for fid = 1:length(fileslist)
    load(fileslist{fid})
    tplot = linspace(0,24,t+1);
    figure()
    hold all
    plot(tplot,numVehiclesRebalancing,'g')
    plot(tplot,numVehiclesBusy-numVehiclesRebalancing,'b')
    plot(tplot,sum(custWaiting'),'r','LineWidth',3)
    axis([0,24,0,4500])
    legend('Rebalancing','Carrying passengers','Waiting customers','FontSize',16)
    xlabel('Time')
    ylabel('Count')
    title([num2str(v),' vehicles'],'FontSize',20)
    allWaitTimes = zeros(numDelivered, 1);
    ccc = 1;
    for i = 1:cc-1
        if customer(i).delivered
            allWaitTimes(ccc) = customer(i).waitTime;
            ccc = ccc+1;
        end
    end
    fprintf('%d vehicles\n',v)
    fprintf('Mean wait time:   %f\nMedian wait time: %f\n',mean(allWaitTimes),median(allWaitTimes))
    
    
end
