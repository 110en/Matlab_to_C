% Code for data loading: Ali Nilforoushan
% Code for plotting: Aditya Muppala (University of Michigan, Ann Arbor, ECE Department)

%% Clear All of the Params

clc
clear
clear path
close all 

%% Reading Arrays from Files

fhat3_dim = load('fhat3_dim.txt');
padded_kx = fhat3_dim(1);
padded_ky = fhat3_dim(2);
N = fhat3_dim(3);

fid = fopen('fft3_fhat.bin', 'rb');
raw = fread(fid, [2, padded_kx * padded_ky * N], 'double');
fclose(fid);

% Combine real and imaginary parts
fhat3 = complex(raw(1,:), raw(2,:));
fhat3 = reshape(fhat3, [padded_kx, padded_ky, N]);

%fhat3 = reshape(comp, [padded_kx, padded_ky, N]);

fid = fopen('xImg.txt','r');
xIm = fread(fid,'double');
fclose(fid);

fid = fopen('yImg.txt','r');
yIm = fread(fid,'double');
fclose(fid);

fid = fopen('distZ.txt','r');
distZ = fread(fid,'double');
fclose(fid);

%% Plotting MIP

for temp = 20:35
    f = abs(fhat3(:,:,temp));
    [X,Y] = ndgrid(xIm,yIm);
    figure(2); clf
    surf(X,Y,abs(f)/max(f(:)));
    shading interp;
    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('f(x,y) - Target Scene');
    title('z = ', num2str(temp));
    hold on;
    colormap gray;
    view([0 90])
    set(gca, 'clim', [0.1 0.8]);
    colorbar
    w = waitforbuttonpress;
end

%% Plotting Slices

f3D = abs(fhat3(:,:,19:22));
f = circshift(f,-20,1);
f = circshift(f,-5,2);
zslice = distZ(19:22);
[Z,Y,X] = meshgrid(zslice,xIm,yIm);
xslice = [];   
yslice = [];
f3D_2 = permute(f3D,[1 3 2]);
figure(2); clf
slice(Z.*100,Y.*100,X.*100,f3D_2./max(f3D(:)),zslice.*100,yslice,xslice)
shading interp;
xlabel('Z (mm)');
ylabel('X (mm)');
zlabel('Y (mm)');
title('Range cuts');
hold on;
pbaspect([2.8 1 1])
colormap gray;
view([47 -10])
set(gca, 'clim', [0 1]);
set(gcf,'position',[50,50,1700,920])
set(gca,'FontSize',20)
set(gca,'GridAlpha',1)
set(gca,'GridLineStyle','--')
set(gca,'fontname','times')
w = waitforbuttonpress;