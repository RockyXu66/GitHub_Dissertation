function comp_ratio=image_cell_PCA(directory,test_directory,dist_dir,file_pre,dataset,threshold,image_width,block_width)
max_num_bands=200;
uniform=1;
clear comp_data rcnst_data; H=0; P=[]; n=1; b=1; m_l=0; m_eig=0;





first_file='frame_1.mhd'; %[V,info]=ReadData3D([directory '\' first_file]); 
N=length(dir(directory))-2;
tic
k=1; band=[]; 
M=30; %number of eigenvectors
V = double(reshape(imread([directory '\' file_pre num2str(2) '.bmp']),1,[]));
lengthV=length(V);
sample_mean=zeros(1,lengthV);
numbers=1:N;
for i=numbers
info.DataFile=[file_pre num2str(i) '.raw']; k=k+1; clear V;
V = double(reshape(imread([directory '\' file_pre num2str(i) '.bmp']),1,[])); sample_mean=sample_mean+V;
end;
sample_mean=sample_mean/N;
sample_var=zeros(1,lengthV);
for i=numbers
info.DataFile=[file_pre num2str(i) '.raw']; k=k+1; clear V;
V = double(reshape(imread([directory '\' file_pre num2str(i) '.bmp']),1,[])); sample_var=sample_var+(V-sample_mean).^2;
end;
num_bands=(image_width/block_width)^2;
dist_dir=[dist_dir '\' dataset '\image_band_' num2str(num_bands) '_' num2str(N)];

mkdir(dist_dir);
sample_var=sample_var/N;
save([dist_dir '\image_sample_mean_' num2str(N) '.mat'],'sample_mean');
save([dist_dir '\image_sample_var_' num2str(N) '.mat'],'sample_var');
min_m=min(sample_mean);
max_m=max(sample_mean);
min_v=min(sample_var);
max_v=max(sample_var);
variab=[];
for i=1:image_width/block_width
for j=1:image_width/block_width
B=zeros(image_width,image_width,3);
B((1:block_width)+(j-1)*block_width,(1:block_width)+(i-1)*block_width,:)=ones(block_width,block_width,3);
band=find(reshape(B,1,[])==1);
variab(j+(i-1)*image_width/block_width)=((sum(sample_var(band))+sum((sample_mean(band)-mean(sample_mean(band))).^2))/length(band));
save([dist_dir '\image_band_' num2str(num_bands) '_' num2str(N) '_' num2str(j+(i-1)*image_width/block_width) '.mat'],'band');
end; 
end;
variab=variab/sum(variab);
toc
for j=1:num_bands
clear tv_volume C daul_eigenvectors dual_eigenvalues explained eigenvectors volumes;
load([dist_dir '\image_band_' num2str(num_bands) '_' num2str(N) '_' num2str(j) '.mat']);
if (variab(j)/(length(band)/lengthV)>0.1)
status=['processing band: ' num2str(j) ' percentage ' num2str(100*length(band)/lengthV) '%']
status='computing covariance matrix'
tic
for i=1:N
clear V; V = double(reshape(imread([directory '\' file_pre num2str(i) '.bmp']),1,[])); volumes(i,:)=V(band)-sample_mean(band); 
end;
C=volumes*volumes';
toc
status='computing eigenspace'
tic
[dual_eigenvectors,dual_eigenvalues,explained] = pcacov(C/(N-1)); eigenvectors=volumes'*dual_eigenvectors;
toc
status='normalizing eigenvectors'
tic
for i=1:N-1
eigenvectors(:,i)=eigenvectors(:,i)./((N-1)*var(eigenvectors(:,i)'*volumes'))^(1/4); 
end;

save([dist_dir '\image_band_' num2str(num_bands) '_' num2str(N) '_' num2str(j) '.mat'],'band','eigenvectors','explained');
status=['band ' num2str(j) ' processed and saved successfully']
%save(['C:\ARTIVVIS Gitlab\Compact Volume Renderer\Compact Volume Renderer\time_varying\supernova\band_based_PCA\CovarMatOfBand_' num2str(j) '.mat'],'band','C');
end;
end;




z=1;
%test_directory=directory;
Ntest=length(dir(test_directory))-2;
load([dist_dir '\image_sample_mean_' num2str(N)]);
status='saving reconstructed test samples'
tic
for i=1:Ntest
V = double(reshape(imread([test_directory '\' file_pre num2str(i) '.bmp']),1,[]))-sample_mean;    
rcnst_data=zeros(1,lengthV);   
for b=1:num_bands
b
num_eig(b)=0;
l(b)=0;
load([dist_dir '\image_band_' num2str(num_bands) '_' num2str(N) '_' num2str(b) '.mat']);
if (variab(b)/(length(band)/lengthV)>0.1)
m_l=m_l+length(band); k=1; while (sum(explained(1:k))<threshold) k=k+1; end; num_eig(b)=k; l(b)=length(band); n=n+1;
comp_data=eigenvectors(:,1:k)'*V(band)'; 
rcnst_data(band)=comp_data'*eigenvectors(:,1:k)'; 
end;
b=b+1;
end;  rcnst_data=rcnst_data+sample_mean; 

imwrite(reshape(round(rcnst_data)/255,size(imread([test_directory '\' file_pre num2str(i) '.bmp']))),[dist_dir '\bandPCAfig_' num2str(num_bands) '_' num2str(N) '_' num2str(threshold) '_' num2str(i) '.bmp']);
end;
toc
% status='saving reconstructed volumes'
% tic
% reconstructed_volume=zeros(1,numel(V)); for i=1:N
% reconstructed_volume=rcnst_data(i,:);
% fid=fopen(['C:\ARTIVVIS Gitlab\Compact Volume Renderer\Compact Volume Renderer\time_varying\supernova\reconstructed\' num2str(i) '.raw'],'w'); %reconstructed_volume=single(reconstructed_volume);
% cnt=fwrite(fid,reshape((reconstructed_volume),size(mha_read_volume(info))),'float');
% fclose(fid); end;
% for i=1:N
% fid=fopen(['C:\ARTIVVIS Gitlab\Compact Volume Renderer\Compact Volume Renderer\time_varying\supernova\reconstructed\frame_' num2str(i) '.mhd'],'w');
% fprintf(fid,['NDims = 3 \nDimSize = 128 128 128 \nElementType = MET_FLOAT \nElementSpacing = 1.0 1.0 1.0 \nElementByteOrderMSB = False \nElementDataFile = ' num2str(i) '.raw']);
% fclose(fid);
% end; toc
% m_l=m_l/num_bands; m_eig=m_eig/num_bands; comp_ratio(max_th-90+1,num_bands-1)=N*length(V)/(num_bands*(m_l+N*m_eig));
% mean_quality(num_bands-1)=mean_quality(num_bands-1)/length(V);
% %load('C:\ARTIVVIS Gitlab\Compact Volume Renderer\Compact Volume Renderer\time_varying\supernova\band_based_PCA\sample_var');
% %p=var(sample_var(band))/var(sample_var);
% salah_image_mean_quality(max_th-90+1,num_bands-1)=salah_image_mean_quality(max_th-90+1,num_bands-1); P=P/H;
comp_ratio=sum(l)/sum(num_eig);