function [Z] = mainsolver(A,y,K,C2,C1,mu,maxiter,S)
%% initialization
[nr,nc,L] = size(S);
B = reshape(y,nr,nc);
sizeD = [nr,nc,L];
N = nr*nc*L;

Zt = zeros(nr,nc,L);
ZE = hyperConvert2D(Zt); 

V1=ZE;
V2=ZE;
V3=ZE;
G1=zeros(size(V1));
G2=zeros(size(V2));
G3=zeros(size(V3));

alpha=[0.5,0.5,0.5]';
x1   = zeros(N,1);
g11  = x1;
x2   = zeros(N,1);
g22  = x2;
x3   = zeros(N,1);
g33  = x3;

patchsize1=8; 
patchsize2=8;
overlap=4;
num1                =   (nr-patchsize1)/(patchsize1-overlap)+1;
num2                =   (nc-patchsize2)/(patchsize2-overlap)+1;
bparams.block_sz    =   [patchsize1,patchsize2];
bparams.overlap_sz  =   [overlap,overlap];
bparams.block_num   =   [num1,num2]; 

fkmeans_omaxpt.careful = 1;
predenoised_blocks = ExtractBlocks(B, bparams);
Y2   = Unfold(predenoised_blocks,size(predenoised_blocks),4);
[aa] = fkmeans(Y2,K,fkmeans_omaxpt);

%% iterations
for i=1:maxiter
%% updating Z
v1 = hyperConvert3D(V1,nr,nc);
v1 = v1(:);
v2 = hyperConvert3D(V2,nr,nc);
v2 = v2(:);
v3 = hyperConvert3D(V3,nr,nc);
v3 = v3(:);

g1 = hyperConvert3D(G1,nr,nc);
g1 = g1(:);
g2 = hyperConvert3D(G2,nr,nc);
g2 = g2(:);
g3 = hyperConvert3D(G3,nr,nc);
g3 = g3(:);

z  = Zt(:);
Zt0 = Zt;

z  = myPCGSolver(A,y,z,mu,v1,v2,v3,g1,g2,g3,x1,x2,x3,g11,g22,g33);

Zt = reshape(z,sizeD);
ZE = hyperConvert2D(Zt);

rmse(i) = getrmse(double(im2uint8(S)),double(im2uint8(Zt))) 
% relchg(i) = norm(Zt(:)-Zt0(:))/norm(Zt0(:));
% relerr(i) = norm(Zt(:)-S(:))/norm(S(:));
% psnr1(i)  = psnr(S,Zt);
% ssim1(i)  = ssim(S,Zt);
% save('changepsnr.mat','psnr1');
% save('changessim.mat','ssim1');
% save('relchg.mat','relchg');
% save('relerr.mat','relerr');

%% updating TV term
mode = 1;
[x1] = aistroTVThres(z,g11,sizeD,C1,mu,alpha(1),mode);
mode = 2;
[x2] = aistroTVThres(z,g22,sizeD,C1,mu,alpha(2),mode);
mode = 3;
[x3] = aistroTVThres(z,g33,sizeD,C1,mu,alpha(3),mode);

%% updating U, V, W (V1, V2, V3)
B1 = ZE-G1/(2*mu);
B1 = hyperConvert3D(B1,nr, nc );
predenoised_blocks1 = ExtractBlocks(B1, bparams);
Z_block1 = zeros( bparams.block_sz(1), bparams.block_sz(2),L, bparams.block_num(1)* bparams.block_num(2));
  
B2 = ZE-G2/(2*mu);
B2 = hyperConvert3D(B2,nr, nc );
predenoised_blocks2 = ExtractBlocks(B2, bparams);
Z_block2 = zeros( bparams.block_sz(1), bparams.block_sz(2),L, bparams.block_num(1)* bparams.block_num(2));
  
B3 = ZE-G3/(2*mu);
B3 = hyperConvert3D(B3,nr, nc );
predenoised_blocks3 = ExtractBlocks(B3, bparams);
Z_block3 = zeros( bparams.block_sz(1), bparams.block_sz(2),L, bparams.block_num(1)* bparams.block_num(2));
  
 for mn=1:max(aa)
    gg = find(aa==mn);
    XES = predenoised_blocks1(:,:,:,gg);
    [a,b,c,d] = size(XES);

    a1 = min(a*b,c*d);
    a2 = min(a*b*c,d);
    a3 = a;
    c1 = sqrt(a1)/(sqrt(a1)+sqrt(a2)+sqrt(a3));
    c2 = sqrt(a2)/(sqrt(a1)+sqrt(a2)+sqrt(a3));
    c3 = sqrt(a3)/(sqrt(a1)+sqrt(a2)+sqrt(a3));

    D1 = C2*c1;
    D2 = C2*c2;
    D3 = C2*c3;
    nsig2 = 1;
    XES = reshape(XES,[a*b c*d]);
    V11 = WNNM(XES, D1/2/mu, nsig2);
    V11 = reshape(V11,[a b c d]);
    Z_block1(:,:,:,gg) = V11;

    XES = predenoised_blocks2(:,:,:,gg);
    [a,b,c,d] = size(XES);
    XES = reshape(XES,[a*b*c d]);
    V22 = WNNM(XES, D2/2/mu, nsig2);
    V22 =reshape(V22,[a b c d]);
    Z_block2(:,:,:,gg) = V22;
    
    XES = predenoised_blocks3(:,:,:,gg);
    [a,b,c,d] = size(XES);
    XES = reshape(XES,[a b*c*d]);
    V33 = WNNM(XES,D3/2/mu,nsig2 );
    V33 = reshape(V33,[a b c d]);
    Z_block3(:,:,:,gg) = V33;
 end
    
V1 = JointBlocks(Z_block1,bparams);
V1 = hyperConvert2D(V1);
V2 = JointBlocks(Z_block2,bparams);
V2 = hyperConvert2D(V2);
V3 = JointBlocks(Z_block3,bparams);
V3 = hyperConvert2D(V3);

%% updating multipliers
G1 = G1+2*mu*(V1-ZE);
G2 = G2+2*mu*(V2-ZE);
G3 = G3+2*mu*(V3-ZE);
g11 = g11 - 2*mu*(z - x1);
g22 = g22 - 2*mu*(z - x2);
g33 = g33 - 2*mu*(z - x3);
end

%% output
Z=Zt;