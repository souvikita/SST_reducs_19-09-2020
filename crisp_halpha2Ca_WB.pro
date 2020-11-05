;;%%%%%%%%%%% CRISP H-alpha to Ca-8542 WB alignment %%%%%%%%%%%%%%%%%%%''

dpath_wb = '/mn/stornext/d18/lapalma/reduc/2020/2020-09-19/CRISP/cubes_wb/'
wb_halpha = readfits(dpath_wb+'wb_6563_2020-09-19T08:23:27_scans=0-149_corrected_im.fits', h1)
wb_ca = readfits(dpath_wb+'wb_8542_2020-09-19T08:23:27_scans=0-149_corrected_im.fits',h2)
shifts = fltarr(2,150);for storing the shifts
tiles = [10,20,30,40]
clips = [12,6,3,1]
destretched_halpha = fltarr(1232,1219,150);destretched cube H-alpha stored 
median_WB = median(wb_halpha)
for t=0,149 do begin
	tmp = wb_halpha[*,*,0,0,t]
	i_nan = where(~finite(tmp), /null)
	tmp[i_nan] = median_WB
	shifts[*,t]=trk(reform(tmp[300:700,400:700]),reform(wb_ca[300:700,400:700,0,0,t])); the X and Y regions indicate the region for performing X-corr
	shifted_wb_halpha=red_rotation(tmp,0,shifts[0,t], shifts[1,t]); Compensating for the shifts obtaing after X-corr
	grid = red_dsgridnest(reform(wb_ca[*,*,0,0,t]),shifted_wb_halpha[0:1231,0:1218], tiles, clips);computing the grid + the new X and Y show the crop 
	destretched_halpha[*,*,t] = rot(red_stretch(shifted_wb_halpha[0:1231,0:1218], grid),-0.4,/interp); directly compensating for the rotation compensation of 0.4 deg
	print,string(13b)+'Compensating the X Y offsets and destretching. Takes time you can have coffee!:  % finished: ',float(t)*100./(149),format='(a,f4.0,$)'
endfor

save, shifts, filename=dpath_wb+'H_alpha_XY_offsets.sav'
lp_write, destretched_halpha, dpath_wb+'wb_6563_2020-09-19T08:23:27_scans=0-149_rot_corrected_im.fcube'
end

;cropped_halpha = reform(wb_halpha[73:1231+73,10:1218+10,0,0,10]) ;just for a given time step index 10
;ca_ref = reform(wb_ca[*,*,0,0,10]) ;; For the reference Ca WB
;tiles = [15,30,45,60]
;clips = [12,6,3,1]
;grid = red_dsgridnest(ca_ref, cropped_halpha, tiles, clips)
;destretched_halpha = red_stretch(cropped_halpha, grid)

;window,0

;plot_image, destretched_halpha, title='6563'

;window,1

;plot_image, ca_ref, title='8542'

;blink, [0,1]

;end
