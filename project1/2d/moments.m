function moments
    close all;
    kmax = 10;
    
    numWalkers = 1000;
    
    moments = zeros(kmax,1);
    values = zeros(kmax,numWalkers);
    
    for k=1:kmax
        n = 2^k;
        
        for walker=1:numWalkers
            walkX = fix(floor((4*rand(n,1)-1))-0.5); %equal prob [-1 0 0 1]
            walkY = fix(floor((4*rand(n,1)-1))-0.5); %equal prob [-1 0 0 1]
            
            values(k,walker) = sqrt(sum(walkX)^2 + sum(walkY)^2);
        end
        
        moments(k) = moment(values(k,:),2);
    end
    
    x = log(2.^(1:kmax));
    y = log(moments);
    
    A = [ones(size(y)) x'];
    B = regress(y,A);
    
    
    plot(x,y);
    axis('equal')
end