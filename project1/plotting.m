function plotting
clear all;
close all;
format long;

energy = dlmread('release/energy.txt');
pressure = dlmread('release/pressure.txt');

time = energy(:,1);
kinetic_energy = energy(:,2);
potential_energy = energy(:,3);
total_energy = energy(:,4);
temperature = energy(:,5);

figure(1)
subplot(4,1,1);
plot(time,kinetic_energy,'r');
hold on
plot(time,potential_energy,'g');
plot(time,total_energy,'b');
legend('Kinetic energy','Potential energy','Total energy');
xlabel('Time')
ylabel('Energy')

subplot(4,1,2);
plot(time,temperature);
legend('Temperature');
xlabel('Time')
ylabel('Temperature')

subplot(4,1,3);
plot(time,pressure(:,2));
legend('Pressure');
xlabel('Time')
ylabel('Pressure')
end