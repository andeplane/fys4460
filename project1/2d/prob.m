function prob
    close all;
    kmax = 10;
    numWalkers = 5000;
    maxY = 0;
    
    values = zeros(kmax,numWalkers);
    points = 2^kmax+1;
    x = linspace(0,2^kmax,points);
    P = zeros(kmax,points);
    
    for k=1:kmax
        n = 2^k;
        
        for walker=1:numWalkers
            walk = floor(4*rand(n,1));
            xp = length(find(walk==0));
            xn = length(find(walk==1));
            yp = length(find(walk==2));
            yn = length(find(walk==3));
            walkX = xp - xn;
            walkY = yp - yn;
            r = sqrt(sum(walkX)^2 + sum(walkY)^2);
            values(k,walker) = r;
            
            index = find(x>=values(k,walker),1);
            
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
       P(k,:) = P(k,:)*n;
       maxY = max(maxY,max(P(k,:)));
       
       plotProbability(P,k,scaledX,colors(k),maxY)
    end
    
end

function plotProbability(P,r,x,color,maxY)
    hold on;
    
    plot(x,P(r,:),color);
    axis([0 50 0 maxY]);
    
end