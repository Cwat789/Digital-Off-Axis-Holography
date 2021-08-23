clear all; close all;
%%In V4 all images are displayed, and most prompts are suppressed

%Hard coded file name
IHO=imread('spray_1560bar_0001.tif','tif');

% figure(1); imshow(IHO);
% title('input image');
% axis off

IHO1 = IHO;
%initial parameters
pixelSize=0.00065; % pixel pitch 0.00065 cm (6.5 um)
lambda=400*10^-9; % wavelength 633 nm            Lambda: IMPORTANT PARAMETER
z=185; % z=M*deltax^2/w; % propagation distance   IMPORTANT PARAMETER


%gets size of the image, this will be useful later
[M,N] = size(IHO);  %in our setup,  image is 1372x1040 pixles

%converts image data matrix to a double
IHO=double(IHO);

%% Spectrum Interpretation

%creates spectrum image of fringe pattern
SP=fftshift(ifft2(fftshift(IHO)));

%displays spectrum to select box size
% figure(2); imshow(50.*mat2gray(abs(SP)));
% title('Hologram spectrum')
% axis off

%asks uset to identify the image in the spectral plane (x, y, and box size)
%Considering adding object identification in the future to circumvent this
% x = input('X val: ');
% y = input('Y val: ');
% boxsize = input('Box size: ');
% %while we are inputing values, might as well ask for z
% z = input('Propigation depth: ');

%these values can be hard coded for convienience in testing
x = 1080;
y = 495;
boxsize = 400;


%Create new spectral image with noise and shadow removed, and spray
%centered, and the rest is zero padded
newSP = zeros([M,N]);
newSP(M/2-boxsize/2+1:M/2+boxsize/2, N/2-boxsize/2+1:N/2+boxsize/2) = SP(y-boxsize/2+1:y+boxsize/2, x-boxsize/2+1: x+boxsize/2);

%displays new hologram spectrum
% figure(3); imshow(50.*mat2gray(abs(newSP)));
% title('New Hologram spectrum')
% axis off

%% Creating and Zero Padding New Fringe Pattern
%creates fringe pattern using new hologram spectrum image 
IHO = fftshift(fft2(fftshift(newSP)));

%displays new hologram
% figure(4); imshow(1.*mat2gray(abs(IHO)));
% title('new fringe, noise removed')
% axis off

%preformes zero padding
paddedimage=zeros(2*N);
paddedimage((N-M/2)+1:(N+M/2),1+N/2:N+N/2)= IHO;



%% Fresnal Defraction Method
% FDM: reconstruction (Fresnel diffraction)
r3=1:2*N;  %2*N
c3=1:2*N;
[C3, R3]=meshgrid(c3, r3);
THOR=((R3-N-1).^2+(C3-N-1).^2).^0.5;
A=THOR.*pixelSize/4;
QP=exp(1i*pi/lambda/z.*(A.^2));
%[x,y] = size(QP)
FTS=fftshift(fft2(fftshift(paddedimage.*QP)));
I2=FTS.*conj(FTS);

%Displays reconstructed image
%weight affects saturation of the image as far as I can tell
weight = 1;
% figure(5); imshow(weight.*mat2gray(I2));
% title('Reconstructed image')
% axis off


%% Final image cropping and saving
%crop zero padding out of reconstructed image
croppedimage = imcrop(I2, [N/2,N/2,N-1,M-1]); %double check that this is right

%self crop image
%croppedimage = imcrop(weight.*mat2gray(I2));

intensity = sum(sum(croppedimage));



% figure(6); imshow(weight.*mat2gray(croppedimage));
% title('Cropped Reconstructed image')
% axis off

figure(7);
subplot(2,3,1)
imshow(IHO1);
title('input image - Diesel Spray Hologram');
axis off


subplot(2,3,2)
imshow(50.*mat2gray(abs(SP)));
title('Hologram spectrum')
axis off


subplot(2,3,3)
imshow(50.*mat2gray(abs(newSP)));
title('New Hologram spectrum')
axis off

subplot(2,3,4)
imshow(1.*mat2gray(abs(IHO)));
title('new fringe, noise removed')
axis off

subplot(2,3,5)
imshow(weight.*mat2gray(I2));
title('Reconstructed image')
axis off

subplot(2,3,6)
imshow(weight.*mat2gray(croppedimage));
title('Cropped Reconstructed image')
axis off


%% Option to save turned off for this example


% saveyn = questdlg('Would you like to save this file?', ...
% 	'Options', ...
% 	'Yes','No ','No');
% if (saveyn == 'Yes')
%     filename = input('input file name: ','s');
%     imwrite(croppedimage,filename)
% end