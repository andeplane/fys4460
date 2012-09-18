function moments
    close all;
    kmax = 10;
    
    numWalkers = 10000;
    
    moments = zeros(kmax,8);
    values = zeros(kmax,numWalkers);
    
    for k=1:kmax
        n = 2^k;
        
        for walker=1:numWalkers
            y = floor(3*rand(n,1)) - 1;
            values(k,walker) = sum(y);
        end
        
        for mom=1:8
           moments(k,mom) = moment(values(k,:),mom);
        end
    end
    x = log(2.^(1:kmax));
    plotMoment(moments,1,x,'b');
    plotMoment(moments,2,x,'r');
    plotMoment(moments,3,x,'g');
    plotMoment(moments,4,x,'k');
    plotMoment(moments,8,x,'c');
    legend('1st','2nd','3rd','4th','8th')
    
end

function plotMoment(moments,r,x,color)
    hold on;
    sprintf('Plotting the %dth moment',r)
    y = log(moments(:,r));
    
    A = [ones(size(y)) x'];
    B = regress(y,A);
    
    
    plot(x,y,color);
    axis('equal')
end