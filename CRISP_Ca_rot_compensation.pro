nx_ca = 1232
ny_ca = 1219
nw_ca = 9
nt_ca = 150
dpath_nb = '/mn/stornext/d18/lapalma/reduc/2020/2020-09-19/CRISP/cubes_nb/'
dat_nb_ca = readfits(dpath_nb+'nb_8542_2020-09-19T08:23:27_scans=0-149_corrected_im.fits',header)
rot_nb_ca = fltarr(nx_ca,ny_ca,nw_ca,nt_ca)

for t=0, nt_ca-1 do begin
	for wav=0, nw_ca-1 do begin
		rot_nb_ca[*,*,wav,t] = rot(reform(dat_nb_ca[*,*,wav,0,t]),-0.4,/interp)
	endfor
endfor
;create the header
head = red_unpol_lpheader(nx_ca,ny_ca,nt_ca*nw_ca,/float)
; open the file and write the header
ofile = dpath_nb+'nb_8542_2020-09-19T08:23:27_scans=0-149_rot_corrected_im.fcube'
openw,lun,ofile,/get_lun
writeu,lun,head

;Write the data using assoc

a = assoc(lun,fltarr(nx_ca,ny_ca,nw_ca,/nozero),512)
for tt=0, nt_ca-1 do a[tt] = rot_nb_ca[*,*,*,tt]
;;close file
free_lun, lun

;;create the spectral (.sp) cube

;ofile= dpath_nb+'nb_8542_2020-09-19T08:23:27_scans=0-149_rot_corrected_sp.fcube'

red_flipthecube_unpol,ofile,nt=nt_ca,nw=nw_ca

end

