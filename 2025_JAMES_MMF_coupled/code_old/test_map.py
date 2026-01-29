import ngl
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs

fig_file = 'figs/test_map'

wks = ngl.open_wks('png',fig_file)
plot = [None]*2
res = ngl.Resources()
res.nglDraw                      = False
res.nglFrame                     = False
res.mpProjection                 = 'Robinson'

res.mpPerimOn = True
res.mpPerimLineThicknessF = 10
res.mpPerimLineColor = 'red'
res.mpPerimDrawOrder = 'PostDraw'
plot[0] = ngl.map(wks, res)

res.mpPerimOn = False
# res.lbPerimOn = False
# res.lgPerimOn = False
# res.txPerimOn = False
# res.mpPerimLineColor = 'white'
# res.mpOutlineOn = False
# res.mpPerimLineThicknessF = 0.001
# res.mpPerimDrawOrder = 'PreDraw'
# res.tmXTOn = False
# res.tmYROn = False

plot[1] = ngl.map(wks, res)

# ngl.draw(plot[1])
# ngl.frame(wks)

ngl.panel(wks,plot,[len(plot),1],False)
ngl.end()

print(f'\n  {fig_file}.png\n')
# hc.trim_png(fig_file)
