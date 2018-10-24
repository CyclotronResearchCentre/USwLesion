function h = crc_MSchalResults(val,lab_title)
% Small function to automatize the display of results

figure, h = gcf;

SZ = size(val);
if numel(SZ)<3
    % Single plot, for 1D/2D array
    m_val = mean(val);
    s_val = std(val);
    plot_error_bar(m_val,s_val,lab_title)
elseif numel(SZ)==3
    % mulitple plots
    for ii=1:SZ(3)
        m_val = mean(squeeze(val(:,:,ii)));
        s_val = std(squeeze(val(:,:,ii)));
        subplot(1,SZ(3),ii)
        plot_error_bar(m_val,s_val,lab_title(ii,:))
    end
else
    error('Wrong input format, should be 2D or 3D array')
end

end

%% SUB-FUNCTION
function plot_error_bar(mv,sv,lab)
hold on
bar(1:3,mv)
errorbar(1:3,mv,sv,'.')
title(lab)
end
