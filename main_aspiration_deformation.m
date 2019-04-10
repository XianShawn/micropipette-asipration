clc
clear
close all

readerobj=VideoReader('Mar19WT3.avi');

filename = 'cut_1.mat';
imgname = 'cut_1';


numFrames = get(readerobj, 'NumberOfFrames');
%%
for i = 1:10:numFrames
    vidFrames = read(readerobj,i);
    imshow(vidFrames(:,:,:,1));
    set(gcf,'outerposition',get(0,'screensize'));
    title('use left button on mouse to double click on target point.');
    [x,y] = getpts;
    position (i,1)=x;
    position (i,2)=y;
end

%%
% define position x and position y
position_x = position(:,1);
position_y = position(:,2);

% get rid of the zeros in x and y
position_x(position_x==0) = [];
position_y(position_y==0) = [];

% define new matrix for saving position information
position2(:,1) = position_x;
position2(:,2) = position_y;

save('Mar19WT3.mat','position2');

%% plot the distance/deformation over time
load('Mar19WT1.mat');

distance  = zeros(size(position2));
distance(:,2) = [];
distance(1,1) = 0;

for i = 2:numel(distance)
    distance(i,1)=((position2(i,1)-position2(1,1))^2+(position2(i,2)-position2(1,2))^2)^0.5;
    distance(i,1)=5*distance(i,1)/32;
    %covert pixel value to physical value 32 pixel for 5 micron
end

save('Mar19WT1.mat','position2','distance');

%% extract viscosity and stiffness from the deformation curve
time = 0:0.33:0.33*(numel(distance)-1); %0.33 is the time scale factor, 3 frame analyzed per sec
figure(1), plot(time,distance);
%curve fitting to calculate the values of E and N

%define geometry parameters of the micropipette y = c+a*exp(-b*x)
% model Viscoelastic behaviour of human mesenchymal stem cells, July 2008BMC Cell Biology 9(1):40
pipette.pressure = 100; %unit in pascal
pipette.radius = 0.0000025; %unit in meter

curvefit = aspirationFit(time, distance);

a = curvefit.a;
b = curvefit.b;
c = curvefit.c;

K1 = 2*pipette.pressure*pipette.radius/3.14/c
K2 = K1*a/(a+c);
N = (-1)*K1*a/c/b   % apparent viscosity

E_equ = 1.5*K1;      % equalibrium elastic modulus
E_ins = 1.5*(K1+K2); % instantous elastic modulus

