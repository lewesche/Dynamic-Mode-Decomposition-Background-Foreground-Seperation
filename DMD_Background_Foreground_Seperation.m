% Leif Wesche
% Dynamic Mode Decompoisition Background and Foreground Seperation

%% Load video 1

close all
clear all
clc

test=[' Test 1 '];

res=[1080, 1920]/10;

vid=read(VideoReader('01040026Trim.MP4'));
vid=vid(:,:,:,1:2:end);
frame=size(vid);
frame=frame(4);
dt=9/frame;

for i=[1:frame]
video(:,:,i)=imresize(flipud(double(rgb2gray(vid(:,:,:,i)))), res);
end

%% Load video 2

close all
clear all
clc

test=[' Test 2 '];

res=[1080, 1920]/10;

vid=read(VideoReader('01040017Trim.MP4'));
vid=vid(:,:,:,1:2:end);
frame=size(vid);
frame=frame(4);
dt=12/frame;

for i=[1:frame]
video(:,:,i)=imresize(flipud(double(rgb2gray(vid(:,:,:,i)))), res);
end


%% Play video
close all
clc

figure
for i=[1:frame]
pcolor(video(:,:,i)), colormap(gray), shading interp
pause(dt)
end

%% Construct X
close all
clc

X=[];
for i=[1:frame-1]
dum=video(:,:,i);
X=[X, dum(:)];
end

X_tild=[];
for i=[2:frame]
dum=video(:,:,i);
X_tild=[X_tild, dum(:)];
end

%% DMD
close all
clc

[U, S, V]=svd(X, 'econ');

% figure
% bar((diag(S)/sum(sum(S))))

% U=U(:, 1:40);
% V=V(:, 1:40);
% S=S(1:40, 1:40);

A_tild=(U')*X_tild*V*inv(S);
%Eigen values of A and A_tilde the same

[W,D]=eig(A_tild);

phi=U*W;
%% Reconstruct Video
close all 
clc 

bar(abs(diag(D)));


dum=diag(D);
Omega=zeros(size(D));
for i=[1:frame-1]
Omega(i,i)=log(dum(i));
Omega(i,i)=Omega(i,i)/(dt);
end

b=phi\X(:,1);
t=(1:frame);
X_recon_whole=[];

for i=[1:frame-1]
X_recon_whole=[X_recon_whole, phi*D.^i*b];
end

%% Play Reconstructed Video 
close all
clc

figure
for i=[1:frame-1]
dum=abs(reshape(X_recon_whole(:,i), res));
pcolor(dum), colormap(gray), shading interp
pause(dt)
end

%% Select Eigenvalues corresponding to low Omega values
close all
clc

figure
subplot(2,2,1)
bar(abs(diag(D)));
title(['A)', test, 'Original A Eigen Values'])
ylabel('Eigen Values'); xlabel('Eigen Values #'); 

subplot(2,2,2)
bar(abs(diag(Omega)))
title(['B)', test, 'Mode Frequencies'])
ylabel('\omega'); xlabel('Mode #');  

low_cutoff=1;

D_low=D;
D_low(abs(Omega) > low_cutoff)=0;
subplot(2,2,3);
bar(abs(diag(D_low)));
title(['C)', test, 'Eigen Values Above Cutoff'])
ylabel('Eigen values'); xlabel('Eigen values #'); 

high_cutoff=1;

D_high=D;
D_high(abs(Omega) <= high_cutoff)=0;
subplot(2,2,4);
bar(abs(diag(D_high)));
title(['D)', test, 'Eigen Values Below Cutoff'])
ylabel('Eigen values'); xlabel('Eigen values #'); 

%% Reconstruct Video
X_recon_back=[];
for i=[1:frame-1]
X_recon_back=[X_recon_back, phi*D_low.^i*b];
end

X_recon_fore=[];
for i=[1:frame-1]
X_recon_fore=[X_recon_fore, phi*D_high.^i*b];
end
%% View Original/Background/Foreground
close all
clc

figure
for i=[1:frame-1]

dum1=abs(reshape(X_recon_whole(:,i), res));
subplot(2,2,1)
pcolor(dum1), colormap(gray), shading interp
caxis([0, 255]);
title(['A)', test, 'Original Video'])

dum2=abs(reshape(X_recon_back(:,i), res));
subplot(2,2,2)
pcolor(dum2), colormap(gray), shading interp
caxis([0, 255]);
title(['B)', test, 'Background'])

dum3=abs(reshape(X_recon_fore(:,i), res));
subplot(2,2,3)
pcolor(dum3), colormap(gray), shading interp
caxis([0, 255]);
title(['C)', test, 'Foreground'])

pause(0.01)

end
