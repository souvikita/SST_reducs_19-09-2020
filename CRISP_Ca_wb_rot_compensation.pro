dpath_wb = '/mn/stornext/d18/lapalma/reduc/2020/2020-09-19/CRISP/cubes_wb/'
data_ca_wb = readfits(dpath_wb+'wb_8542_2020-09-19T08:23:27_scans=0-149_corrected_im.fits')
nx_ca = 1232
ny_ca = 1219
nw_ca = 9
nt_ca = 150
rot_wb_ca=fltarr(nx_ca,ny_ca,nt_ca)

for t=0, nt_ca-1 do begin
	rot_wb_ca[*,*,t] = rot(reform(data_ca_wb[*,*,0,0,t]),-0.4,/interp)
endfor

lp_write, rot_wb_ca, dpath_wb+'wb_8542_2020-09-19T08:23:27_scans=0-149_rot_corrected_im.fcube'
end

