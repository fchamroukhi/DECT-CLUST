
function save_pdf(fig_hdl, fn)

set(fig_hdl,'Units','Inches');
pos = get(fig_hdl,'Position');
set(fig_hdl,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig_hdl,fn,'-dpdf','-r0')
