% Holographic Illumination in paraxial approximation
% WALTHER, Mar 2014

clear all;

currentFolder = cd;
if ~exist(fullfile(currentFolder, 'Hologram'), 'dir'), error('Directory not found'); end
path(path, fullfile(currentFolder, 'Hologram'));
if ~exist(fullfile(currentFolder, 'Functions'), 'dir'), error('Functions folder missing in the active directory'); end
path(path, fullfile(currentFolder, 'Functions'));


% ******************** GENERAL PARAMETERS *******************************
% lengts in units of [mm]
lambda = 940e-6;  % wavelength in [mm]
nH2O = 1.33;

rObj = 15.12/2;  % [mm] radius of objective BFP
fObj = 8; %[mm] objective focal length
% fObj = 180 / 25; tube focal length OLYMPUS = 180 mm!

NA = 0.7;

f1 = 300;
f2 = 200;

Mlat = f2/f1;
Mang = f1/f2;

rBeam = 13/2;  % in [mm] AOD coordinates
rAOD = 15/2;

xtot = 20; % in [mm]
nx = 2^nextpow2(2000);

ytot = 20; % in [mm]
ny = 2^nextpow2(2000);

nz = 20;
zDefocus = linspace(-30,30,nz)*1e-3;   % in [mm]

% *********************************************************************
rainbow = [ zeros(51,1),zeros(51,1),linspace(0,1,51)';...
            zeros(51,1),linspace(0,1,51)',ones(51,1);...
            zeros(51,1),ones(51,1),linspace(1,0,51)';...
            linspace(0,1,51)',ones(51,1),zeros(51,1);...
            ones(51,1),linspace(1,0,51)',zeros(51,1)];
        
blueTOyellowTOred = vertcat( [linspace(0,1,2^7)', linspace(0,1,2^7)', linspace(1,0.8, 2^7)'],...
                    [ones(2^7,1),linspace(1,0,2^7)', linspace(0.8,0,2^7)'] );         
       
% ***** AOD-shaped light field in back focal plane (holographic plane) ********
% hy(:,1) = y = spatial ccordinate in y-direction
% hy(:,2) = 1D amplitude function (AODY power modulation)
% hy(:,3) = 1D phase function (AODY phase modulation)
% hy(:,4) = 1D frequency function (MHz)

% hx(:,1) = x =spatial coordinate in x-direction
% hx(:,2) = 1D amplitude function (AODX power modulation)
% hx(:,3) = 1D phase function (AODX phase modulation)
% hx(:,4) = 1D frequency function (MHz)

% 1D spatial coordinates in the ao plane
xao = xtot/Mlat * linspace(-0.5,0.5,nx)';  % in [mm]
yao = ytot/Mlat * linspace(-0.5,0.5,ny)';

% 1D spatial coordinates in the holographic plane
xh = xtot * linspace(-0.5,0.5,nx)';  % in [mm]
yh = ytot * linspace(-0.5,0.5,ny)';

% initialize input data
hy = cat(2, yh, ones(length(yh),1), zeros(length(yh),2));
hx = cat(2, xh, ones(length(xh),1), zeros(length(xh),2));

% *** Define input amplitude and phase ***

fileNameX = fullfile(currentFolder, 'Hologram', 'GS holo plane X_data.mat');
holoInput = importdata(fileNameX);
hx(:,2) = interp1(Mlat*holoInput(:,1),holoInput(:,2),hx(:,1), 'nearest',0);
hx(:,3) = interp1(Mlat*holoInput(:,1),holoInput(:,3),hx(:,1), 'nearest',1);
hx(:,4) = interp1(Mlat*holoInput(:,1),holoInput(:,6),hx(:,1), 'nearest',1);

fileNameY = fullfile(currentFolder, 'Hologram', 'GS holo plane Y_data.mat');
holoInput = importdata(fileNameY);
hy(:,2) = interp1(Mlat*holoInput(:,1),holoInput(:,2),hy(:,1), 'nearest',0);
hy(:,3) = interp1(Mlat*holoInput(:,1),holoInput(:,3),hy(:,1), 'nearest',1);
hy(:,4) = interp1(Mlat*holoInput(:,1),holoInput(:,6),hy(:,1), 'nearest',1);

% ***********************************************************************************
% ***** Trim 2-D Field in back focal plane (holographic plane) = input field ********
% hyx(:,:,1) = yy = spatial coordinates in y-direction
% hyx(:,:,2) = xx = spatial coordinates in x-direction
% hyx(:,:,3) = 2D amplitude function 
% hyx(:,:,4) = 2D phase function

[xx,yy] = meshgrid(xh,yh);
hyx = cat(3,yy,xx,ones(size(yy)),zeros(size(yy)));
hyx(:,:,3) = repmat(hy(:,2),1,size(hyx,2)) .* repmat(hx(:,2)',size(hyx,1),1);
hyx(:,:,4) = repmat(hy(:,3),1,size(hyx,2)) + repmat(hx(:,3)',size(hyx,1),1);

insideAOD = hyx(:,:,1).^2+hyx(:,:,2).^2 <= (Mlat*rAOD)^2;

% *** Apply objective aperture
% Hard Aperture
apert = ones(size(hyx(:,:,1)));
apert(~insideAOD) = 0;
hyx(:,:,3) = hyx(:,:,3) .* apert;

% *** Apply Beam Shape
% Gaussian Beam
ampli = exp(-(xx.^2+yy.^2)/(Mlat*rBeam)^2);
%hyx(:,:,3) = hyx(:,:,3) .* ampli;
hyx(:,:,3) = hyx(:,:,3) / max(hyx(cat(3,false([size(insideAOD),2]),insideAOD,false(size(insideAOD)))));

%% **********************************************************************
hWait = waitbar(0,'Calculating light field in object (Fourier) space');

% output calculation
% Spatial frequencies in holographic plane
X = nx/xtot * linspace(-0.5,0.5,nx);  % in [1/mm]
Y = ny/ytot * linspace(-0.5,0.5,ny);
[XX,YY] = meshgrid(X,Y);

% Fourier transform of h(y,x) in holographic plane
field = hyx(:,:,3) .* exp(1i * 2 * pi * hyx(:,:,4));
spec = fftshift(fft2(ifftshift(field),ny,nx));
HYX = cat(3,YY,XX,abs(spec),angle(spec)/(2*pi));
waitbar(1/(2+nz));
clear XX YY field spec;

% *************************************************************************
% ************* Field in the space domain ********************************
% fyx(:,:,1) = yy = spatial coordinates in y-direction
% fyx(:,:,2) = xx = spatial coordinates in x-direction
% fyx(:,:,3) = amplitude function
% fyx(:,:,4) = phase function

% spatial coordinates in front focal plane
x = nx*fObj*lambda/xtot * linspace(-0.5,0.5,nx);  % in [mm]
y = ny*fObj*lambda/ytot * linspace(-0.5,0.5,ny);
[xx,yy] = meshgrid(x,y); 

% Field f(y,x) in front focal plane
fyx = cat(3, yy, xx, HYX(:,:,3), HYX(:,:,4)); 

% spatial frequencies in front focal plane
X = xtot/(fObj*lambda) * linspace(-0.5,0.5,nx);  % in [1/mm]
Y = ytot/(fObj*lambda) * linspace(-0.5,0.5,ny);
[XX,YY] = meshgrid(X,Y);

% Fourier transform of f(y,x)
field = fyx(:,:,3) .* exp(1i * 2*pi * fyx(:,:,4));
spec = fftshift(fft2(ifftshift(field),ny,nx));
FYX = cat(3,YY,XX, abs(spec), angle(spec)/(2*pi));
waitbar(2/(2+nz));
clear x y xx yy XX YY field spec;

% *************************************************************************
% ************* Field propagation to defocus plane ******* ****************
% dyxz = amplitude function in defocused z-planes

dyxz = zeros([size(fyx(:,:,1)),nz],'single');

% claculation of input field
trueWave = (FYX(:,:,1).^2+FYX(:,:,2).^2) <= 1/lambda^2;
H = exp(-1i*2*pi*zDefocus(1)*sqrt((1/lambda^2-FYX(:,:,1).^2-FYX(:,:,2).^2).*trueWave));
spec = H.*(FYX(:,:,3) .* exp(1i * 2 * pi * FYX(:,:,4)));
field = fftshift(ifft2(ifftshift(spec), ny, nx));
dyxz(:,:,1) = sqrt(single(field .* conj(field)));   

for i = 2:nz
   dz = zDefocus(i)-zDefocus(i-1);
   H = exp(-1i*2*pi*dz*sqrt((1/lambda^2-FYX(:,:,1).^2-FYX(:,:,2).^2).*trueWave));
   spec = H .* fftshift(fft2(ifftshift(field),ny,nx));
   field = fftshift(ifft2(ifftshift(spec), ny, nx));
    
   dyxz(:,:,i) = sqrt(single(field .* conj(field)));
   waitbar(i/(nz-1));
    
end
close(hWait);
clear spec field;

%% **** FIGURE 1 ****
% New Coordinate System = AOD coordinates
yAOD = hyx(:,1,1)/Mlat;
xAOD = hyx(1,:,2)/Mlat;
yxAOD = hyx(:,:,[1,2])/Mlat;

xrange = find(abs(xAOD) <= rAOD);
yrange = find(abs(yAOD) <= rAOD);
yxrange = yxAOD(:,:,1).^2+yxAOD(:,:,2).^2 <= rAOD^2;

scrsize = get(groot, 'ScreenSize');
aspect = scrsize(4)/scrsize(3);
 
fig1 = figure;
set(fig1, 'units','normalized','outerposition',[0 0.25 0.5*aspect/2, 0.57]);
set(fig1, 'Name', 'Holographic Plane');

s1=subplot(2,1,1);
ampl = hyx(:,:,3);
ampl(~yxrange) = 0;
imMin = 0;
imMax = 1;
im = imresize(ampl(yrange,xrange), [512,512]);   
imagesc(xAOD(xrange),yAOD(yrange),im, [imMin,imMax]);  
colormap(s1, gray);
set(gca, 'FontSize', 12, 'TickDir', 'out', 'TickLength', [0.02,0.02]);
set(gca, 'XTick', [-6,-3,0,3,6], 'YTick', [-6,-3,0,3,6]);
xlabel('x-AOD (mm)', 'FontSize', 12);
ylabel('y-AOD (mm)', 'FontSize', 12);
title('Amplitude', 'FontSize', 12);
axis tight;

s2=subplot(2,1,2);
phase = hyx(:,:,4);
phase(~yxrange) = 0;
phase = phase - (min(phase(:)) + (max(phase(:))-min(phase(:)))/2);
imMin = -0.5;
imMax = 0.5;
im = phase;
im = mat2gray(im, double([imMin,imMax]));
[im,~] = gray2ind(im,255);
im = ind2rgb(im,blueTOyellowTOred);
im = uint8(255*im);
im(repmat(~yxrange,[1,1,3])) = 255;
im = imresize(im(yrange,xrange,:), [512,512]);
imagesc(xAOD(xrange),yAOD(yrange),im);  
set(gca, 'FontSize', 12, 'TickDir', 'out', 'TickLength', [0.02,0.02]);
set(gca, 'XTick', [-6,-3,0,3,6], 'YTick', [-6,-3,0,3,6]);
axis tight;
ylabel('y-AOD (mm)', 'FontSize', 12); 
xlabel('x-AOD (mm)','FontSize', 12);
title('Phase', 'FontSize', 12);

%% **** FIGURE 2 ****
fig2 = figure;
set(fig2, 'units','normalized','outerposition',[0.25 0.25 0.5*aspect 0.5]);
set(fig2, 'Name', 'Intensity Plots');
colorMAP = load('ocean_r.mat', 'cMap');

FOVyx = [15,15] * 1e-3;
yrange = find(abs(fyx(:,1,1)) <= FOVyx(1));
xrange = find(abs(fyx(1,:,2)) <= FOVyx(2));
zrange = find(abs(zDefocus) <= 30e-3);

subplot(2,2,1)
intensity = fyx(yrange,xrange,3).^2;
im = imresize(intensity, [512,512]);
imMin = 0;
%imMax = prctile(intensity(:),99.9);
imMax = 0.95 * max(intensity(:));
imagesc(1e3*fyx(1,xrange,2),1e3*fyx(yrange,1,1),im, [imMin,imMax]);
colormap(rainbow);
colormap(colorMAP.cMap);
set(gca, 'FontSize', 12, 'TickDir', 'out', 'TickLength', [0.02,0.02]);
%set(gca, 'XTick', [-20,-10,0,10,20], 'YTick', [-20,-10,0,10,20]);
set(gca, 'XTick', [-15,-10,-5,0,5,10,15], 'YTick', [-15,-10,-5,0,5,10,15]);
set(gca, 'FontSize', 12, 'TickDir', 'out');
ylabel('y (\mum)', 'FontSize', 12); 
xlabel('x (\mum)','FontSize', 12);
title('XY-Focal Plane', 'FontSize', 12, 'FontWeight', 'bold');

subplot(2,2,2);
intensity = max(dyxz(yrange,xrange,3:end),[],3).^2;
im = imresize(intensity, [512,512]);
imMin = 0;
imMax = prctile(intensity(:),99.9);
imMax = 0.95 * max(intensity(:));
imagesc(1e3*fyx(1,xrange,2),1e3*fyx(yrange,1,1), im,[imMin,imMax]);
colormap(rainbow);
colormap(colorMAP.cMap);
set(gca, 'FontSize', 12, 'TickDir', 'out', 'TickLength', [0.02,0.02] );
%set(gca, 'XTick', [-20,-10,0,10,20], 'YTick', [-20,-10,0,10,20]);
set(gca, 'XTick', [-15,-10,-5,0,5,10,15], 'YTick', [-15,-10,-5,0,5,10,15]);
xlabel('x (\mum)','FontSize', 12 );
ylabel('y (\mum)', 'FontSize', 12 ); 
title('XY-Plane & Max Z', 'FontSize', 12, 'FontWeight', 'bold' );

subplot(2,2,3);
xzero = find(fyx(1,:,2) >= 0, 1, 'first');
dzero = find(zDefocus >= 0,1,'first');
intensity = (squeeze(max(dyxz(:,xrange,zrange),[],1))).^2;
im = imresize(intensity, [512,512]);
imMin = 0;
imMax = prctile(intensity(:),99.9);
imMax = 0.95 * max(intensity(:));
imagesc(1e3*fyx(1,xrange,2), 1e3*zDefocus(zrange), im', [imMin,imMax]);
line([1e3*(min(squeeze(fyx(yrange,xzero,1)))),...
    1e3*(max(squeeze(fyx(yrange,xzero,1))))],zDefocus(dzero)*[1,1], 'LineStyle','--','Color', 'k', 'LineWidth', 1); 
colormap(rainbow);
colormap(colorMAP.cMap);
set(gca, 'FontSize', 12, 'TickDir', 'out', 'TickLength', [0.02,0.02]);
set(gca, 'XTick', [-20,-10,0,10,20], 'YTick', [-20,-10,0,10,20]);
set(gca, 'XTick', [-15,-10,-5,0,5,10,15], 'YTick', [-20,-10,0,10,20]);
xlabel('x (\mum)', 'FontSize', 12);
ylabel('z (\mum)', 'FontSize', 12); 
title('XZ-Plane & Max Y','FontSize', 12, 'FontWeight', 'bold' );

subplot(2,2,4);
yzero = find(fyx(:,1,1) >= 0, 1, 'first');
dzero = find(zDefocus>=0,1,'first');
intensity = squeeze(max(dyxz(yrange,:,zrange),[],2)).^2;
im = imresize(intensity, [512,512]);
imagesc(1e3*fyx(1,xrange,2), 1e3*zDefocus(zrange), im', [imMin,imMax]);
line([1e3*(min(squeeze(fyx(yrange,xzero,1)))),...
    1e3*(max(squeeze(fyx(yrange,xzero,1))))],zDefocus(dzero)*[1,1], 'LineStyle','--','Color', 'k', 'LineWidth', 1); 
colormap(fig2, colorMAP.cMap);
set(gca, 'FontSize', 12, 'TickDir', 'out', 'TickLength', [0.02,0.02]);
set(gca, 'XTick', [-20,-10,0,10,20], 'YTick', [-20,-10,0,10,20]);
set(gca, 'XTick', [-15,-10,-5,0,5,10,15], 'YTick', [-20,-10,0,10,20]);
xlabel('y (\mum)', 'FontSize', 12);
ylabel('z (\mum)', 'FontSize', 12); 
title('YZ-Plane & Max X','FontSize', 12, 'FontWeight', 'bold'); 


%% **** FIGURE 3 ****
fig3 = figure;
set(fig3, 'units','normalized','outerposition',[0.25 0.25 0.5*aspect 0.5]);
set(fig3, 'Name', 'Intensity**2 Plots');
colorMAP = load('gnuplot2.mat', 'cMap');

FOVyx = [15,15] * 1e-3;
yrange = find(abs(fyx(:,1,1)) <= FOVyx(1));
xrange = find(abs(fyx(1,:,2)) <= FOVyx(2));
zrange = find(abs(zDefocus) <= 30e-3);

subplot(2,2,1)
intensity = fyx(yrange,xrange,3).^4;
im = imresize(intensity, [512,512]);
imMin = 0;
%imMax = prctile(intensity(:),99.9);
imMax = 0.9 * max(intensity(:));
imagesc(1e3*fyx(1,xrange,2),1e3*fyx(yrange,1,1),im, [imMin,imMax]);
colormap(rainbow);
%colormap(fig3, cubehelix_r.cMap);
set(gca, 'FontSize', 12, 'TickDir', 'out', 'TickLength', [0.02,0.02]);
%set(gca, 'XTick', [-20,-10,0,10,20], 'YTick', [-20,-10,0,10,20]);
set(gca, 'XTick', [-15,-10,-5,0,5,10,15], 'YTick', [-15,-10,-5,0,5,10,15]);
set(gca, 'FontSize', 14, 'TickDir', 'out');
ylabel('y (\mum)', 'FontSize', 12); 
xlabel('x (\mum)','FontSize', 12);
title('XY-Focal Plane', 'FontSize', 12, 'FontWeight', 'bold' );

subplot(2,2,2);
intensity = max(dyxz(yrange,xrange,3:end),[],3).^4;
im = imresize(intensity, [512,512]);
imMin = 0;
%imMax = prctile(intensity(:),99.9);
imMax = 0.9 * max(intensity(:));
imagesc(1e3*fyx(1,xrange,2),1e3*fyx(yrange,1,1), im,[imMin,imMax]);
colormap(rainbow);
%colormap(fig3, cubehelix_r.cMap);
set(gca, 'FontSize', 12, 'TickDir', 'out', 'TickLength', [0.02,0.02] );
%set(gca, 'XTick', [-20,-10,0,10,20], 'YTick', [-20,-10,0,10,20]);
set(gca, 'XTick', [-15,-10,-5,0,5,10,15], 'YTick', [-15,-10,-5,0,5,10,15]);
xlabel('x (\mum)','FontSize', 12 );
ylabel('y (\mum)', 'FontSize', 12 ); 
title('XY-Plane & Max Z', 'FontSize', 12, 'FontWeight', 'bold', 'FontWeight', 'bold' );

subplot(2,2,3);
xzero = find(fyx(1,:,2) >= 0, 1, 'first');
dzero = find(zDefocus >= 0,1,'first');
intensity = (squeeze(max(dyxz(:,xrange,zrange),[],1))).^4;
im = imresize(intensity, [512,512]);
imMin = 0;
%imMax = prctile(intensity(:),99.9);
imMax = 0.9 * max(intensity(:));
imagesc(1e3*fyx(1,xrange,2), 1e3*zDefocus(zrange), im', [imMin,imMax]);
line([1e3*(min(squeeze(dyxz(yrange,xzero,1)))),...
    1e3*(max(squeeze(dyxz(yrange,xzero,1))))],zDefocus(dzero)*[1,1], 'LineStyle','--','Color', 'w', 'LineWidth', 1.2); 
colormap(rainbow);
%colormap(fig3, cubehelix_r.cMap);
set(gca, 'FontSize', 12, 'TickDir', 'out', 'TickLength', [0.02,0.02]);
set(gca, 'XTick', [-15,-10,-5,0,5,10,15], 'YTick', [-20,-10,0,10,20]);
xlabel('x (\mum)', 'FontSize', 12);
ylabel('z (\mum)', 'FontSize', 12); 
title('XZ-Plane & Max Y','FontSize', 12, 'FontWeight', 'bold' );

subplot(2,2,4);
yzero = find(fyx(:,1,1) >= 0, 1, 'first');
dzero = find(zDefocus>=0,1,'first');
intensity = squeeze(max(dyxz(yrange,:,zrange),[],2)).^4;
im = imresize(intensity, [512,512]);
imMin = 0;
imMax = prctile(intensity(:),99.9);
imMax = 0.9 * max(intensity(:));
imagesc(1e3*fyx(1,xrange,2), 1e3*zDefocus(zrange), im', [imMin,imMax]);
line([1e3*(min(squeeze(dyxz(yrange,xzero,1)))),...
    1e3*(max(squeeze(dyxz(yrange,xzero,1))))],zDefocus(dzero)*[1,1], 'LineStyle','--','Color', 'w', 'LineWidth', 1.2); 
colormap(fig3, colorMAP.cMap);
colormap(rainbow);
set(gca, 'FontSize', 12, 'TickDir', 'out', 'TickLength', [0.02,0.02]);
set(gca, 'XTick', [-15,-10,-5,0,5,10,15], 'YTick', [-20,-10,0,10,20]);
xlabel('y (\mum)', 'FontSize', 12);
ylabel('z (\mum)', 'FontSize', 12); 
title('YZ-Plane & Max X','FontSize', 12,'FontWeight', 'bold' );


%% **** FIGURE 4 ****
pause(0.001);
fig5 = figure;
set(fig5, 'units','normalized','outerposition',[0.6 0.5 0.5*aspect 0.5]);
set(fig5, 'Name', '3D Intensity Plot');
FOVyxz = [30, 30, 30] * 1e-3;   % yx field of view in fig 1 in [mm];
FOVyInd = abs(fyx(:,1,1)) <= FOVyxz(1)/2;
FOVxInd = abs(fyx(1,:,2)) <= FOVyxz(2)/2;
FOVzInd = abs(zDefocus) <= FOVyxz(3);
dThreshold = 1/exp(1);
[yyy,xxx,zzz] = meshgrid(1e3*fyx(FOVyInd,1,1),1e3*fyx(1,FOVxInd,2),1e3*zDefocus(FOVzInd));
dyxzSM = smooth3( dyxz(FOVyInd,FOVxInd,FOVzInd).^2,'gauss',[15,15,5]);
p = patch(isosurface(yyy,xxx,zzz,dyxzSM,dThreshold * max(dyxzSM(:))));
isonormals(yyy,xxx,zzz,dyxzSM,p);
set(p,'Facecolor', 'red', 'Edgecolor','none');
grid on;
daspect([1,1,1]);
set(gca,'ydir','reverse');
axis(1e3 * [-FOVyxz(2) FOVyxz(2) -FOVyxz(1) FOVyxz(2) -FOVyxz(3) FOVyxz(3)]);
set(gca, 'Zdir', 'reverse');
set(gca, 'FontSize', 12);
xlabel('X [\mum]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y [\mum]', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('Z [\mum]', 'FontSize', 12, 'FontWeight', 'bold');
set(gca,'YTick',[-10,-5,0,5,10]);
set(gca,'XTick',[-10,-5,0,5,10]);
set(gca,'ZTick',[-10,-5,0,5,10]);

view(60,20);

light('Position',[1,0,0],'Style','infinite');
light('Position',[0,0,-1],'Style','infinite');
lighting flat;

