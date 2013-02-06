function plotting
clear all;
close all;

energy = dlmread('build-debug/energy.dat');

temperature = dlmread('build-debug/temperature.dat');
pressure = dlmread('build-debug/pressure.dat');
diffusion = dlmread('build-debug/diffusion.dat');

t = energy(:,1);
Ek = energy(:,2);
Ep = energy(:,3);
E = energy(:,4);

figure(1)
subplot(4,1,1);
plot(t,Ek,'r');
hold on
plot(t,Ep,'g');
plot(t,E,'b');
legend('Kinetic energy','Potential energy','Total energy');
xlabel('Time')
ylabel('Energy')

subplot(4,1,2);
plot(t,temperature(:,2));
legend('Temperature');
xlabel('Time')
ylabel('Temperature')

subplot(4,1,3);
plot(t,pressure(:,2));
legend('Pressure');
xlabel('Time')
ylabel('Pressure')

subplot(4,1,4);
plot(diffusion(:,1),diffusion(:,2));
legend('Diffusion constant');
xlabel('Time')
ylabel('Diffusion constant')
    
end