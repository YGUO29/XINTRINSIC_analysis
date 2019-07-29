function [mask_sig, mask_consis, mask_corr] = PixSelect(para,DataMat,I,p_th,r_th,c_th,plot_on)

% DataMat = [para.nStim, para.nStim, para.height, para.width, para.nStim]
% I: surface image
% p_th: threshold for significance test, typical value: 0.0005;
% r_th: threshold for consistency index, typical value: 0.2
% c_th: threshold for correlation measurement, typical value: 0.6
% plot: 0 = no plots; 1 = plot masks; 2 = plot all figures

%% pixel selection: 1. select significantly responded pixels

mask_sig = false(para.height, para.width);

X = DataMat(:,:,:,:,1:floor(para.fr*para.preStim));
X = squeeze(mean(X,1));
X = permute(X,[2 3 4 1]);
X = reshape(X,[size(X,1),size(X,2),size(X,3)*size(X,4)]);
X = reshape(X,[size(X,1)*size(X,2),size(X,3)]);

Y = DataMat(:,:,:,:,floor(para.fr*para.preStim)+1:floor(para.fr*para.preStim+ + para.fr*para.durStim));
Y = squeeze(mean(Y,1));
Y = permute(Y,[2 3 4 1]);
Y = reshape(Y,[size(Y,1),size(Y,2),size(Y,3)*size(Y,4)]);
Y = reshape(Y,[size(Y,1)*size(Y,2),size(Y,3)]);

[~,p] = ttest2(X',Y');
mask1 = zeros(para.height,para.width);
mask1(p<p_th) = 1;

 
% mask_sig generated.

%% Pixel selection: 2. select reliably responded pixels

% construct a single-trial response matrix 
V = zeros(2,para.nStim,para.width*para.height);
ind = [1,2]; % select two reps to compare
for i = 1:2
    mov = squeeze(DataMat(ind(i),:,:,:,:));
    % mov_base & mov_amp: [para.nStim,para.height,para.width];
    mov_base = squeeze(mean(mov(:,:,:,1:floor(para.fr*para.preStim)),4)); % first second: baseline, for one rep
    mov_amp = squeeze(mean(mov(:,:,:,floor(para.fr*para.preStim)+1:floor(para.fr*para.preStim+ + para.fr*para.durStim)),4));
    % figure,imshow(img_base,[])
    % figure,imshow(img_amp,[])
    mov_relamp = (mov_amp - mov_base)./mov_base; 
    V(i,:,:) = reshape(mov_relamp,para.nStim,para.width*para.height);
end
% for two repetitions:
V1 = squeeze(V(1,:,:)); 
V2 = squeeze(V(2,:,:)); 

%% ==== calculate consistency index r ====
r = zeros(para.width*para.height,1);
for i = 1 : para.width*para.height
    v1 = V1(:,i); v2 = V2(:,i);
    proj = v2.*(v2'*v1)/norm(v2);
    r(i) = 1-norm(v1-proj)/norm(v1);
end


mask_consis = zeros(para.height,para.width);
mask_consis((r>r_th)) = 1;

% ==== calculate correlation index c ====
c = zeros(para.width*para.height,1);
for i = 1:para.width*para.height
    v1 = V1(:,i); v2 = V2(:,i);
    c(i) = corr(v1,v2);
end

mask_corr = zeros(para.height,para.width);
mask_corr(c > c_th) = 1;


% plot masks on top of surface image
if plot_on == 1
I_norm = (I - min(min(I)))./(max(max(I)) - min(min(I)));

img = repmat(I_norm,1,1,3); % three layers, representing R,G,B 
img(:,:,1) = img(:,:,1) + mask_sig;
figure,imagesc(img),axis image, axis off, title('Significance mask')

img = repmat(I_norm,1,1,3); % three layers, representing R,G,B 
img(:,:,1) = img(:,:,1) + mask_consis;
figure,imagesc(img),axis image, axis off, title('Consistency mask')

img = repmat(I_norm,1,1,3); % three layers, representing R,G,B 
img(:,:,1) = img(:,:,1) + mask_corr;
figure,imagesc(img),axis image, axis off, title('Correlation mask')
figure,imagesc(reshape(c,para.height,para.width)), axis image, axis off, colorbar, title('Correlation value map')

end

end