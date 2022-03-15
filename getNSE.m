function [nse_map, noise_map_x, noise_map_y] = getNSE(X, Y, noise_correct)
% based on Norman-Haignere 2018, measure similarity of two vectors 
% (e.g. for one pixel, the actural response vs. the predicted response to all sounds)
if ~noise_correct 
    % without noise correction
    % x, y: n x 1 vectors
    nPix = size(X,2);
    nse_map = zeros(1, nPix);

    for i = 1:nPix
        x = X(:,i); y = Y(:,i);
        
        a = mean(x.^2) + mean(y.^2) - 2*mean(x.*y);
        b = mean(x.^2) + mean(y.^2) - 2*mean(x).*mean(y);
        nse_map(i) = a./b;
        % 1: x and y are completely independent
        % 0: x and y are the same
    end
    noise_map_x = [];
    noise_map_y = [];
    
else % with noise correction
    assert(length(size(X)) == 3, 'Make sure input matrix X has at least 2 repetitions')
    assert(length(size(Y)) == 3, 'Make sure input matrix Y has at least 2 repetitions')
    
    nPix = size(X,2);
    nse_map = zeros(1, nPix);
    noise_map_x = nse_map;
    noise_map_y = nse_map;
    for i = 1:nPix
        x1 = X(:,i,1); x2 = X(:,i,2); 
        y1 = Y(:,i,1); y2 = Y(:,i,2);
        mu_sx2 = mean(x1.^2)/2 + mean(x2.^2)/2 - mean((x1-x2).^2)/2;
        mu_sy2 = mean(y1.^2)/2 + mean(y2.^2)/2 - mean((y1-y2).^2)/2;
        mu_sxsy = (mean(x1.*y1) + mean(x1.*y2) + mean(x2.*y1) + mean(x2.*y2))/4;
        mu_sx = mean(x1)/2 + mean(x2)/2;
        mu_sy = mean(y1)/2 + mean(y2)/2;
        a = mu_sx2 + mu_sy2 - 2 * mu_sxsy;
        b = mu_sx2 + mu_sy2 - 2 * mu_sx .* mu_sy;
        nse_map(i) = a./b;
        
        noise_map_x(i) = mean((x1-x2).^2)/2;
        noise_map_y(i) = mean((y1-y2).^2)/2;
        % 1: x and y are completely independent
        % 0: x and y are the same
    end

end


end