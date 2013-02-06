function plotVelocityProfiles
    atomIndex = @(i,N) mod(i-1,N)+1;
    sampleIndex = @(i,N) floor((i-1)/N)+1;
    
    data = dlmread('velocities.dat');
    dataSize = size(data);
    
    N = data(1,1);
    
    numberOfSamples = (dataSize(1)-1)/N;
    
    v = zeros(numberOfSamples,N);
    
    for i=1:dataSize(1)-1
       v(sampleIndex(i,N),atomIndex(i,N)) = data(i+1,1);
       sprintf('sampleindex = %d, atomindex = %d',sampleIndex(i,N),atomIndex(i,N))
    end
    size(v)
    hist(v(20,:))
end