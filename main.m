
clear all
close all
clc

addpath('Mex/')
mex -O -g CFLAGS="\$CFLAGS -Wall -std=c99" Mex/PM_2D.c  -outdir Mex/
mex -O -g CFLAGS="\$CFLAGS -Wall -std=c99" Mex/rebuilt.c  -outdir Mex/


pw = 2;  %patch_size
iter = 5;  %patchmatch iteration number

a = double(imread('./data/img_1.jpg'));
b = double(imread('./data/img_2.jpg'));

%a to b
tic;
[nnf, nnfd] = PM_2D(single(a), single(b), pw, iter);
toc;
disp('PM ended');

%rebuilt
[a_r] = rebuilt(single(b), int32(nnf(:,:,:,end)), pw);
a_r = double(a_r);

imwrite(uint8(a_r),'img_1_r.png');


figure, 
subplot 221
imagesc(uint8(a))
title('A')
subplot 222
imagesc(uint8(b))
title('B')
subplot 223
imagesc(uint8(a_r))
title('A rebuilt from B')
subplot 224
imagesc(sum((a-a_r).^2,3))
title(sprintf('Error (psnr=%1.3f)',psnr(uint8(a),uint8(a_r))));
