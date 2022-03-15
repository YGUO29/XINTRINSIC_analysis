%% regression with sound features
% plot coefficients maps 
% need: R (response matrix of each component), K (number of components)
% load('D:\=code=\Sound_analysis\F_test.mat') 
load('D:\SynologyDrive\=data=\F_halfcosine_marm_NatJM.mat') 
% load('D:\SynologyDrive\=data=\F_halfcosine_marm_NatPNAS-168.mat') 
% load('D:\SynologyDrive\=data=\F_halfcosine_marm_NatCustomize.mat') 

Tbl = readtable('D:\SynologyDrive\=code=\_python\sound_analysis_python\eGeMAPSv02_Nat165_functionals.csv');
%% test: correlation with openSMILE features
close all
iK = 1;
K = size(Rs{iK},2);
R = Rs{iK};
% for model-matched
% R = R(1:2:end,:);

nFeat       = size(Tbl,2);
result_p    = zeros(nFeat,K);
result_r    = zeros(nFeat,K);
% ============= Regression ===============
for iComp = 1:K
    yy = R(:,iComp);
%     yy = yy(16:end);

    for iFeat = 3:nFeat
        xx = table2array(Tbl(:,iFeat));
%         xx = xx(16:end);
        rr = corrcoef(xx,yy);
        result_r(iFeat,iComp) = rr(1,2);
        
%         [p,rsq,~] = RSquared(xx,yy);
%         result_p(iFeat,iComp) = p(1);
%         result_r(iFeat,iComp) = sqrt(rsq);
    end
end
%% correlation with features from auditory models
% frequency + spectrotemporal modulation
% F = F_0214;
Decomp = Decomp_102D;
close all
iK = 1;
K = size(Decomp.Rs{iK},2);
R = Decomp.Rs{iK};
% for model-matched
% R = R(1:2:end,:);

nFeat       = size(F.F_mat,1);
result_p    = zeros(nFeat,K);
result_r    = zeros(nFeat,K);
% ============= Regression ===============
for iComp = 1:K
    yy = R(:,iComp);
%     yy = yy(16:end);

    for iFeat = 1:nFeat
        xx = F.F_mat(iFeat,:)';
%         xx = xx(16:end);
        rr = corrcoef(xx,yy);
        result_r(iFeat,iComp) = rr(1,2);
        
%         [p,rsq,~] = RSquared(xx,yy);
%         result_p(iFeat,iComp) = p(1);
%         result_r(iFeat,iComp) = sqrt(rsq);
    end
end

% ===========================================
% plot feature correlations
scr_size = get(0,'ScreenSize'); % left, bottom, width and height
panel_width = scr_size(3)/10;
f = figurex([1, 1, K*panel_width, scr_size(4)*0.8]);
fig_pos = get(f,'position');

% ind = [1 2 5 6 4 3];
% ind = [1 2 3 4 5 6]; % plot order
% ind = [3 2 4 6 5 1]; % 80Z
ind = 1:K;
spectemp_r = cell(1,K);
for i = 1:K
    iComp = ind(i);

    % 1st row: plot components 
%     cutoff = 0.15;
    cutoff      = mean(Decomp.Ws{iK}(i,:)) + 7*std(Decomp.Ws{iK}(i,:)); % variable cutoff values for each components
    comp_temp   = zeros(para.height, para.width);
    comp_temp(para.ind_save) = Decomp.Ws{iK}(i,:);
%     comp_temp = comp(:,:,iComp);
    subplot(3, K, i);
    imagesc(comp_temp,cutoff.*[-1 1]); axis image 
%     cmap = cbrewer('div', 'RdYlGn', 256);
    colormap(gca, 'jet'), 
    colorbar
    if para.mirror
        set(gca,'XDir','reverse')
    else 
    end
    axis off
    
    % 2nd row: correlation with frequency power
    subplot(3, K, i+K);
    plot(F.FreqBounds(1:end-1), result_r(1:F.nFreq, ind(i)), 'linewidth',4,'Marker','x')
    hold on, plot([F.FreqBounds(1:2), F.FreqBounds(end-1)], [0 0 0],'linestyle','--','color','k')    
    set(gca,'xscale','log');
%     ymax = max(max(  abs( result_r(1:F.nFreq,1:K) )  )); ylim([-0.77,0.77]);
    ymax = max(max(  abs( result_r(1:F.nFreq,1:K) )  )); ylim([-1,1]);
    set(gca,'xtick',F.FreqBounds(1:2:end))
    set(gca,'ytick',-0.9:0.3:0.9)
    set(gca,'xticklabels',arrayfun(@num2str,F.FreqBounds(1:2:end-1),'UniformOutput',false))
    set(gca,'fontsize',18);
    axis square
    if i == 1
        xlabel('Frequency','fontsize',20),
        ylabel('Correlation Coeff.','fontsize',20),
    end
    xtickangle(45)

    % 3rd row: correlation with spectrotemporal modulation
    subplot(3, K, i + 2*K), 
    % use averaged modulation maps
    spectemp_r{i} = reshape(result_r(F.nFreq+F.nTemp+F.nSpec+1:F.nFreq++F.nTemp+F.nSpec+F.nSpectemp,ind(i)), ...
                        size(F.spectemp_mod,1), size(F.spectemp_mod,2));
    
    cmax = max(max(  abs( result_r(F.nFreq+F.nTemp+F.nSpec+1 : F.nFreq+F.nTemp+F.nSpec+F.nSpectemp,1:K) )  )); 
    imagesc(flipud(spectemp_r{i}),[-cmax, cmax]); 
    colormap(gca, 'jet')
%     colorbar
    set(gca,'ytick',1:2:length(F.spec_mod_rates))
    set(gca,'yticklabels',arrayfun(@num2str,fliplr(F.spec_mod_rates(1:2:end)),'UniformOutput',false))
    set(gca,'xtick',1:2:length(F.temp_mod_rates)-1)
    set(gca,'xticklabels',arrayfun(@num2str,F.temp_mod_rates(2:2:end),'UniformOutput',false))
    axis square
    set(gca,'fontsize',18);
    if i == 1
        xlabel('Temp. Mod. (Hz)','fontsize',20),
        ylabel('Spec. Mod. (cyc/oct)','fontsize',20),
    end
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

getColorbar(cmax,'number')
% saveas(f1,'132D_session1_FreqPower.png')
% saveas(f2,'132D_session1_SpecTemp.png')
% saveas(f,'80Z_session1_reg.png')
%% test: try using weighted average modulation maps
figure,
ind = 1:K;
for i = 1:K
    iComp = ind(i);

    % 1st row: temporal rates
    ax{i+K} = subplot(3, K, i);
    plot(F.temp_mod_rates(2:end), result_r(F.nFreq+1 : F.nFreq+F.nTemp, ind(i)), 'linewidth',4,'Marker','x')
    hold on, plot([F.temp_mod_rates(2), F.temp_mod_rates(end)], [0 0],'linestyle','--','color','k')    
    set(gca,'xscale','log');
%     ymax = max(max(  abs( result_r(1:F.nFreq,1:K) )  )); ylim([-0.77,0.77]);
    ymax = max(max(  abs( result_r(1:F.nFreq,1:K) )  )); ylim([-1,1]);
    set(gca,'xtick',F.FreqBounds(1:2:end))
    set(gca,'ytick',-0.9:0.3:0.9)
    set(gca,'xticklabels',arrayfun(@num2str,F.FreqBounds(1:2:end-1),'UniformOutput',false))
    set(gca,'fontsize',24);
    axis square
    if i == 1
    xlabel('Temporal modulation','fontsize',24),
    ylabel('Correlation Coeff.','fontsize',24),
    end
    xtickangle(45)
    
    % 2nd row: correlation with frequency power
    ax{i+K} = subplot(3, K, i+K);
    plot(F.spec_mod_rates, result_r(F.nFreq+F.nTemp+1 : F.nFreq+F.nTemp+F.nSpec, ind(i)), 'linewidth',4,'Marker','x')
    hold on, plot([F.spec_mod_rates(1), F.spec_mod_rates(end)], [0 0],'linestyle','--','color','k')    
    set(gca,'xscale','log');
%     ymax = max(max(  abs( result_r(1:F.nFreq,1:K) )  )); ylim([-0.77,0.77]);
    ymax = max(max(  abs( result_r(1:F.nFreq,1:K) )  )); ylim([-1,1]);
    set(gca,'xtick',F.FreqBounds(1:2:end))
    set(gca,'ytick',-0.9:0.3:0.9)
    set(gca,'xticklabels',arrayfun(@num2str,F.FreqBounds(1:2:end-1),'UniformOutput',false))
    set(gca,'fontsize',24);
    axis square
    if i == 1
    xlabel('Temporal modulation','fontsize',24),
    ylabel('Correlation Coeff.','fontsize',24),
    end
    xtickangle(45)

    % 3rd row: use weighted average modulation maps
    ax{i+K} = subplot(3, K, i+2*K);
    spectemp_r = reshape(result_r(F.nFreq+F.nTemp+F.nSpec+F.nSpectemp+1:F.nFreq++F.nTemp+F.nSpec+2*F.nSpectemp,ind(i)), ...
        size(F.spectemp_mod,1), size(F.spectemp_mod,2));
    
    cmax = max(max(  abs( result_r(F.nFreq+1 : F.nFreq+F.nSpectemp,1:K) )  )); 
    imagesc(flipud(spectemp_r),[-cmax, cmax]), colormap('jet')
%     colorbar
    set(gca,'ytick',1:2:length(F.spec_mod_rates))
    set(gca,'yticklabels',arrayfun(@num2str,fliplr(F.spec_mod_rates(1:2:end)),'UniformOutput',false))
    set(gca,'xtick',1:2:length(F.temp_mod_rates)-1)
    
    set(gca,'xticklabels',arrayfun(@num2str,F.temp_mod_rates(2:2:end),'UniformOutput',false))
    axis square
    set(gca,'fontsize',24);

end