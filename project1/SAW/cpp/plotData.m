function plotData
    close all;

    A = load('data.dat');
    s = size(A);
    N_w = s(2);
    j_max = s(1);
    
    figure(1);

    hold('on');
    x = 0:1:40;

    for j=1:j_max
        if(sum(A(j,:)) > 0) 
            plot(x, histc(A(j,:), x)/N_w, 'color', rndclr());
        end
    end
    
    sprintf('A=1: %f',sum(A(j,:)==1)/length(A(j,:)))
    sprintf('A=2: %f',sum(A(j,:)==2)/length(A(j,:)))
    sprintf('A=sqrt2: %f',sum(A(find(A>1.4))<1.5)/length(A(j,:)))

    xlabel('r');
    ylabel('P(r)');
    legend('2','4','8','16','32','64','128', '256', '512', '1024')
end