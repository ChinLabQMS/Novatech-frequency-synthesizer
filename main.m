%% Initialization
clc;clear;close all;

saveDir = datetime("now",format='yyyy-MM-dd');
wavelength = 1.064e-6;
defaultFreq0 = 80.0;
defaultFreq1 = 80.0;

saveFile = false;
% saveFile = true;

%% Parameters
traj.para = [0.84, 1]; % Degree of "flatness" for flattop trajectory
traj.shape = @Cos_flat;
traj.dist = 279.809; % Total distance in mm (absolute value)
traj.step = 5; % Step size in the integer multiple of 100us
traj.time = 453.4; % Total time in ms

%% Generate curve - frequency sweep in MHz
amp = traj.dist/(wavelength*integral(@(x)abs(traj.shape(x,traj.para)),0,1)*traj.time);
curve = @(x) amp*traj.shape(x/traj.time,traj.para)*1e-6;

%% Generate and save plot/table/parameters
t = 0:0.1:traj.time;
numStep = length(t);

freq0 = defaultFreq0*ones(1,numStep);
freq1 = defaultFreq1*ones(1,numStep)+curve(fix(10*t/traj.step)*traj.step*0.1);

plot(t,freq0,t,freq1);
legend('freq0','freq1');
xlabel('time(ms)');
ylabel('frequency(MHz)');
title(traj_name(traj,' '));
a = gca;
a.Position(2) = 0.3;
a.Position(4) = 0.6;
max_freq_diff = max(abs(diff(freq1)));
txt = sprintf('Max frequency jump is %.2fkHz\nMax acceleration is %.2fm/s^2\nTotal distance is %.2fmm',...
    max_freq_diff*1e3,max_freq_diff*1e10*wavelength/traj.step,traj.dist);
annotation('textbox',[0.15,0.05,0.8,0.15],'String',txt,'EdgeColor','None');

if saveFile
    mkdir(saveDir);
    saveas(gcf,sprintf('%s/%s_%s.png',saveDir,datestr(now,'HH-MM'),traj_name(traj,'_')));
    fileName = sprintf('%s/%s_%s.m',saveDir,datestr(now,'HH-MM'),traj_name(traj,'_'));
    save(fileName,'traj');
    fileName = sprintf('%s/%s_%s.txt',saveDir,datestr(now,'HH-MM'),traj_name(traj,'_'));
    fid = fopen(fileName,'w');
    formatSpec = '0 %011.7f 00000 0696 0255\n1 %011.7f 00000 0820 0255\n';
    fprintf(fid,formatSpec,[freq0(1),freq1(1)]);
    formatSpec = '0 %011.7f 00000 0696 0001\n1 %011.7f 00000 0820 0001\n';
    fprintf(fid,formatSpec,[freq0(2:end-1);freq1(2:end-1)]);
    formatSpec = '0 %011.7f 00000 0696 0000\n1 %011.7f 00000 0820 0000\n';
    fprintf(fid,formatSpec,[freq0(end);freq1(end)]);
    fclose(fid);
end

%% 
function str = traj_name(traj,s)
    shape_name = strrep(func2str(traj.shape),'_',s);
    if contains(shape_name,'flat')
        shape_name = sprintf('%s%sa=%.1f',shape_name,s,traj.para(1));
    end
    if contains(shape_name,'asym')
        shape_name = sprintf('%s%sb=%.1f',shape_name,s,traj.para(2));
    end
    str = sprintf('%s%s%.1fmm%s%.1fms%s%dus',...
        shape_name,s,traj.dist,s,traj.time,s,100*traj.step);
end

%% Normalized trajectories shape
function y = Step(~,~)
    y = 1;
end

function y = Sin(x,~)
    y = sin(pi*x);
end

function y = Sin_round(x,~)
    y = sin(2*pi*x);
end

function y = Cos(x,~)
    y = (1-cos(2*pi*x))/2;
end

function y = Cos_flat(x,a)
    a = a(1);
    y = (x<(1-a)/2).*Cos(x/(1-a),a) +...
        (x>(1+a)/2).*Cos((1-x)/(1-a),a)+...
        ((x>=(1-a)/2) & (x<=(1+a)/2));
end

function y = Cos_round(x,~)
    y = (x<1/2).*Cos(2*x,0) - (x>=1/2).*Cos(2*x-1,0);
end

function y = Cos_flat_round(x,a) 
    y = (x<1/2).*Cos_flat(2*x,a) - (x>=1/2).*Cos_flat(2*x-1,a);
end

function y = Cos_flat_asym(x,a)
    b = a(2);
    a = a(1);
    y = (x<(1-a)/(1+b)).*Cos(x*(1+b)/(2*(1-a))) +...
        (x>(1+a*b)/(1+b)).*Cos((1-x)*(1+b)/(2*b*(1-a))) +...
        ((x>=(1-a)/(1+b)) & (x<=(1+a*b)/(1+b)));
        
end

function y = Cos_flat_asym_round(x,a)
    y = (x<1/2).*Cos_flat_asym(2*x,a) - (x>=1/2).*Cos_flat(2*x-1,a);
end