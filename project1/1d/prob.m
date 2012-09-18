function moments
    close all;
    kmax = 10;
    numWalkers = 100000;
    maxY = 0;
    
    values = zeros(kmax,numWalkers);
    P = zeros(kmax,2^(kmax+1)+1);
    
    x = linspace(-2^kmax,2^kmax,2^(kmax+1)+1);
    zeroIndex = find(x==0);
    
    for k=1:kmax
        n = 2^k;
        
        for walker=1:numWalkers
            y = floor(3*rand(n,1)) - 1;
            values(n,walker) = sum(y);
            index = zeroIndex + values(n,walker);
            
            P(k,index) = P(k,index) + 1;
        end
        P(k,:) = P(k,:)/numWalkers;
        
        maxY = max(maxY,max(P(k,:)));
    end
    
    colors = ['r ','k ','b ','c ','m ','r-', 'k-','b-','c-','m-'];
    for i=1:kmax
        plotProbability(P,i,x,colors(i),maxY)
    end
    
    figure;
    
    for k=1:kmax 
       n = 2^k;
       scaledX = x / sqrt(n);
       P(k,:) = P(k,:)*sqrt(n);
       maxY = max(maxY,max(P(k,:)));
       
       plotProbability(P,k,scaledX,colors(k),maxY)
    end
    
end

function plotProbability(P,r,x,color,maxY)
    hold on;
    
    plot(x,P(r,:),color);
    axis([-80 80 0 maxY]);
    
end