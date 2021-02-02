%% regression with sound features
% plot coefficients maps 
% need: R (response matrix of each component), K (number of components)
% load('D:\=code=\Sound_analysis\F_test.mat') 
% load('D:\SynologyDrive\=data=\F_halfcosine_marm_4reps.mat') 
load('D:\SynologyDrive\=data=\F_halfcosine_marm.mat') 
% load('D:\SynologyDrive\=data=\F_halfcosine_marm_LZVoc_Allmix.mat') 

%%
K = 10;
R = Rs{K};

%
nFeat       = size(F.F_mat,1);
result_p = zeros(nFeat,K);
result_r = zeros(nFeat,K);
% ============= Regression ===============
for iComp = 1:K
    yy = R(:,iComp);
    for iFeat = 1:nFeat
        xx = F.F_mat(iFeat,:)';
%         [p,rsq,~] = RSquared(xx,yy);
%         result_p(iFeat,iComp) = p(1);
%         result_r(iFeat,iComp) = sqrt(rsq);
        rr = corrcoef(xx,yy);
        result_r(iFeat,iComp) = rr(1,2);
    end
end

%% ===========================================
% plot feature correlations
scr_size = get(0,'ScreenSize');
scr_size = scr_size(3:4); %width and height
% f1 = figure,set(gcf,'position',[1,scr_size(2)./2,scr_size(1),scr_size(2)./3]);
% f2 = figure,set(gcf,'position',[1,1,scr_size(1),scr_size(2)./3]);
f = figurex([1,1,0.8.*scr_size]); 

% ind = [1 2 5 6 4 3];
% ind = [1 2 3 4 5 6]; % plot order
% ind = [3 2 4 6 5 1]; % 80Z
ind = 1:K;
for i = 1:K
    iComp = ind(i);
% for iComp = 1:K
    % correlation with frequency power
%     figure(f1)    
    subplot(2,K,i),
    
    plot(F.FreqBounds(1:end-1), result_r(1:F.nFreq,ind(i)), 'linewidth',4,'Marker','x')
    hold on, plot([F.FreqBounds(1:2), F.FreqBounds(end-1)], [0 0 0],'linestyle','--','color','k')    
    set(gca,'xscale','log');
%     ymax = max(max(  abs( result_r(1:F.nFreq,1:K) )  )); ylim([-0.77,0.77]);
    ymax = max(max(  abs( result_r(1:F.nFreq,1:K) )  )); ylim([-0.77,1]);

    set(gca,'xtick',F.FreqBounds(1:2:end))
    set(gca,'ytick',-0.6:0.2:1)
    set(gca,'xticklabels',arrayfun(@num2str,F.FreqBounds(1:2:end-1),'UniformOutput',false))
    set(gca,'fontsize',24);
    axis square
%     xlabel('Frequency','fontsize',14),
%     title(['Correlation coefficient, component #',num2str(ind(iComp))],'fontsize',14)
    xtickangle(45)
%     figure(f2)

    subplot(2,K,i + K), 
    % use averaged modulation maps
%     spectemp_r = reshape(result_r(F.nFreq+F.nTemp+F.nSpec+1:F.nFreq++F.nTemp+F.nSpec+F.nSpectemp,ind(i)), ...
%                         size(F.spectemp_mod,1), size(F.spectemp_mod,2));

    % use weighted average modulation maps
    spectemp_r = reshape(result_r(F.nFreq+F.nTemp+F.nSpec+F.nSpectemp+1:F.nFreq++F.nTemp+F.nSpec+2*F.nSpectemp,ind(i)), ...
        size(F.spectemp_mod,1), size(F.spectemp_mod,2));
    cmax = max(max(  abs( result_r(F.nFreq+1 : F.nFreq+F.nSpectemp,1:K) )  )); 
    imagesc(flipud(spectemp_r),[-cmax, cmax]), colormap('jet')
    colorbar
    set(gca,'ytick',1:2:length(F.spec_mod_rates))
    set(gca,'yticklabels',arrayfun(@num2str,fliplr(F.spec_mod_rates(1:2:end)),'UniformOutput',false))
    set(gca,'xtick',1:2:length(F.temp_mod_rates)-1)
    set(gca,'xticklabels',arrayfun(@num2str,F.temp_mod_rates(2:2:end),'UniformOutput',false))
    axis square
    set(gca,'fontsize',24);

%     xlabel('Temporal modulation rate (cycles/s)','fontsize',14),
%     ylabel('Spectral modulation rate (cycles/octave)','fontsize',14),
%     title(['Correlation coefficient, component #',num2str(ind(iComp))],'fontsize',14)

%   correlation with the full spectrotemporal modulation power map
%     subplot(2,K,iComp + K), 
%     spectemp_r = reshape(result_r(F.nFreq+1:F.nFreq+F.nSpectemp_full,ind(iComp)), ...
%                         size(F.spectemp_mod_full,1), size(F.spectemp_mod_full,2));
%     cmax = max(max(  abs( result_r(F.nFreq+1 : F.nFreq+F.nSpectemp_full,1:K) )  )); 
%     imagesc(flipud(spectemp_r),[-cmax, cmax]),colorbar, colormap('jet')
%     set(gca,'yticklabels',arrayfun(@num2str,fliplr(F.spec_mod_rates),'UniformOutput',false), 'fontsize',12)
%     set(gca,'xticklabels',arrayfun(@num2str,F.temp_mod_rates,'UniformOutput',false), 'fontsize',12)
%     xlabel('Temporal modulation rate (cycles/s)','fontsize',14),
%     ylabel('Spectral modulation rate (cycles/octave)','fontsize',14),
%     title(['Correlation coefficient, component #',num2str(ind(iComp))],'fontsize',14)

%     saveas(f,['80Z_session1_component_',num2str(iComp),'.png'])

end
% saveas(f1,'132D_session1_FreqPower.png')
% saveas(f2,'132D_session1_SpecTemp.png')
% saveas(f,'80Z_session1_reg.png')
