txtsz=24;
set(gca,'FontSize',txtsz)
screen_size = get(0, 'ScreenSize');
h_plotszz = get(gca, 'title');
kpp_plotszz = get(gca, 'xlabel');
l_plotszz = get(gca, 'ylabel');
m_plotszz = get(gca, 'zlabel');
set(h_plotszz, 'FontName', 'Helvetica', 'FontSize', txtsz*1.5)
set(kpp_plotszz, 'FontName', 'Helvetica', 'FontSize', txtsz)
set(l_plotszz, 'FontName', 'Helvetica', 'FontSize', txtsz)
set(m_plotszz, 'FontName', 'Helvetica', 'FontSize', txtsz)
set(gcf, 'PaperOrientation', 'portrait');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperType', 'A4');
% set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
% set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','normalized');
% set(gcf,'PaperPosition', [0 0 1 1]);
%  set(gcf, 'PaperPosition', [1 1 27 19]);
 set(gcf, 'Position',  0.95*[screen_size(1) screen_size(2) screen_size(3) screen_size(4) ] );
 
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + 3*ti(1);
% bottom = outerpos(2) + 1.9*ti(2);
% ax_width = outerpos(3) - 1*ti(1) - 1*ti(3);
% ax_height = outerpos(4) - 1*ti(2) - 1*ti(4);
% ax.Position = [left bottom ax_width ax_height];