function output = myRegression(z,Y,type)
switch type
    case 'regular'
        [b,~,~,~,stats] = regress(z, Y); 
        output.coef = b'; % row vector

        % stats: R-squared (uncorrected), F-stat, p-value, error variance 
        output.Rsquared = stats(1);
        
    case 'lasso'
        % watch out here!!! do not add column of 'ones' in predictors for
        % lass regression; Add back the intercept later
        [b, FitInfo] = lasso(Y(:,2:end), z, 'CV',5);
        % default lambda search region:
        % If you do not supply Lambda, then lasso calculates the largest value 
        % of Lambda that gives a nonnull model. default lambda ratio = 1e-4
        b0 = FitInfo.Intercept(FitInfo.Index1SE); % IndexMinMSE or Index1SE
        output.coef = [b0; b(:,FitInfo.Index1SE)]';
        
        SSresid = sum((Y*output.coef' - z).^2);
        SStotal = (length(z) - 1) * var(z);
        output.Rsquared = 1 - SSresid/SStotal; % uncorrected r-squared
     
%         lassoPlot(b,FitInfo,'PlotType','CV');
%         legend('show') % Show legend
    otherwise 
end
assert(size(output.coef,1) == 1)
end

%%
% [B, FitInfo] = lasso(Y, z,'CV', 3);
% b = B(:, FitInfo.IndexMinMSE); % IndexMinMSE or Index1SE
% b0 = FitInfo.Intercept(FitInfo.IndexMinMSE);
% zpredict = Y*b + b0;
% rr = corrcoef(zpredict, z);
% Corr(k,i) = rr(1,2);
% Coef{i}(k,:) = b';