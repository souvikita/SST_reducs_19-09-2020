dpath_nb='/mn/stornext/d18/lapalma/reduc/2020/2020-09-19/CRISP/cubes_nb/'
dpath_wb='/mn/stornext/d18/lapalma/reduc/2020/2020-09-19/CRISP/cubes_wb/'

nx_ca = 1232
ny_ca = 1219
nw_ca = 9
nt_ca = 150
dat_nb_ca = readfits(dpath_nb+'nb_8542_2020-09-19T08:23:27_scans=0-149_corrected_im.fits',header)
dat_nb_H = readfits(dpath_nb+'nb_6563_2020-09-19T08:23:27_scans=0-149_corrected_im.fits',header)

restore,dpath_wb+'H_alpha_XY_offsets.sav';shifts Array[2,150]
restore,dpath_wb+'Grids_after_destretch_Halpha.sav'; grid Array[2, 40, 40, 150]
tiles = [10,20,30,40] ; for the destretch
clips = [12,6,3,1]
destretched_halpha_NB =fltarr(nx_ca,ny_ca,31,nt_ca);Because we want the NB h-alpha cubes to be aligned with Ca 8542

for scans=0,nt_ca-1 do begin
	for wav=0,30 do begin
		tmp = dat_nb_H[*,*,wav,0,scans]
		i_nan = where(~finite(tmp), /null)
		tmp[i_nan] = median(tmp,/double);2.3272820e-08;;Median value of the NB H-alpha cube. Chosen this because the WB borders has this value. 
		shifted_NB_halpha=red_rotation(tmp,0,shifts[0,scans], shifts[1,scans])
		destretched_halpha_NB[*,*,wav,scans] = rot(red_stretch(shifted_NB_halpha[0:1231,0:1218], grid[*,*,*,scans]),-0.4,/interp)
	endfor
	print,string(13b)+'Compensating the X Y offsets and destretching. Takes time you can have coffee!:  % finished: ',float(scans)*100./(149),format='(a,f4.0,$)'
endfor
;create the header
head = red_unpol_lpheader(nx_ca,ny_ca,nt_ca*31,/float)
; open the file and write the header
ofile = dpath_nb+'nb_6563_2020-09-19T08:23:27_scans=0-149_rot_corrected_im.fcube'
openw,lun,ofile,/get_lun
writeu,lun,head

;Write the data using assoc

a = assoc(lun,fltarr(nx_ca,ny_ca,31,/nozero),512)
for tt=0, nt_ca-1 do a[tt] = destretched_halpha_NB[*,*,*,tt]
;;close file
free_lun, lun
;;create the spectral .sp cube
red_flipthecube_unpol,ofile,nt=nt_ca,nw=31
end
