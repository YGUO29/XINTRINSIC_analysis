function plotTonotopyWithComponents(iK, Decomp, T, para)

% input: a comp*pixel matrix (weight matrix)
close all
% iK = 1;
nComp = size(Decomp.Ws{iK},1);
ct = getContour_compmponent(Decomp.Ws{iK}, Decomp.comp{iK});

% verification: plot components contour with components
figurex;
for iComp = 1:nComp
    ax = subplot(1, nComp, iComp);
    max_amp = max(max( abs(Decomp.comp{iK}(:,:,iComp)) ));
    imagesc(Decomp.comp{iK}(:,:,iComp), [-max_amp, max_amp]), axis image, axis off, colormap(jet)
    if para.mirror
        set(gca,'XDir','reverse')
    else 
    end
    plotContour(ct(iComp),'RdBu');
    
end

% align: plot components contour with tonotopy
figurex;
for iComp = 1:nComp
    subplot(1, nComp, iComp)
    imagesc(T.tonotopy{3}); axis image, axis off
    if para.mirror
        set(gca,'XDir','reverse')
    else 
    end
    plotContour(ct(iComp), 'RdBu');
end