;######################## IRIS time to seconds conversion #################3
function time_iris_secs,time_str
  time_secs=fltarr(n_elements(time_str))
  for i=0,n_elements(time_str)-1 do begin

   time = strmid(time_str[i],11,12)
   str1=strsplit(time,':',/extract)
   time_secs[i]=float(str1[0])*60*60.+float(str1[1])*60+float(str1[2])
  endfor
return, time_secs
end

;######################## New SST cubes already have times in second #########
dpath_wb = '/mn/stornext/d18/lapalma/reduc/2020/2020-09-19/CRISP/cubes_wb/'
red_fitscube_getwcs, dpath_wb+'wb_8542_2020-09-19T08:23:27_scans=0-149_corrected_im.fits', coordinates = coordinates
time_secs_ca = coordinates.time[0,0,0,*] ; This is already in seconds from the headers.

nx_ca = 1232
ny_ca = 1219
nw_ca = 9
nt_ca = 150
dpath_nb = '/mn/stornext/d18/lapalma/reduc/2020/2020-09-19/CRISP/cubes_nb/'
;openr, lun, dpath_nb+'nb_8542_2020-09-19T08:23:27_scans=0-149_corrected_im.fits', /get_lun

;dat_nb_ca = assoc(lun, fltarr(nx_ca,ny_ca,nw_ca,/nozer),512)
dat_nb_ca = readfits(dpath_nb+'nb_8542_2020-09-19T08:23:27_scans=0-149_corrected_im.fits',header)
dpath_l2_iris ='/mn/stornext/d10/HDC2/iris/data/level2_decompressed/2020/09/19/20200919_065936_3633109417/'

f1= file_search(dpath_l2_iris+'*_SJI_*.fits')
read_iris_l2, f1[1], index, data ;Reading IRIS Mg SJI
time_iris = index[*].date_obs ;Entire time series of the dataset
data_Mg = iris_dustbuster(index,data); removing the dust particles
time_secs_Mg = time_iris_secs(time_iris)
trk2 =fltarr(2,nt_ca)
blown_iris = fltarr(2084,2179,nt_ca)
cropped_iris = fltarr(nx_ca,ny_ca,nt_ca)
time_steps_Mg = strarr(nt_ca)

;########### Looping through #########

for j =0, nt_ca-1 do begin
	target_time_sst = time_secs_ca[j]
	minimizer = abs(time_secs_Mg-target_time_sst)
	t1 = data_Mg[*,*,where(minimizer eq min(minimizer))]
	time_steps_Mg = time_iris[where(minimizer eq min(minimizer))]
	;print, time_steps_Mg
	;stop
	Ca_IR = reform(dat_nb_ca[*,*,*,0,j])
	Ca_IR_ref = mean(Ca_IR,dimension=3)
	rot_Ca_IR = rot(Ca_IR_ref,-0.4,cubic=-0.5)
	blown_iris[*,*,j] = reform(congrid(t1,fix(371/0.178),fix(388/0.178),cubic=-0.5))
	blown_iris1 = reform(blown_iris[*,*,j])
	w1 = where(blown_iris1 le 0.)
	blown_iris1[w1] = 0.
	roi_iris = blown_iris1[764:964,843:1343]
	roi_Ca_IR = rot_Ca_IR[250:450,400:900]
	thres = mean(roi_iris[*,0:300])+6.*stddev(roi_iris[*,0:300])
	w_pixels = where(roi_iris[*,0:300] ge thres)
	;print, n_elements(w_pixels)
	;window,0
	;plot_image, roi_iris
	;window,1
	;plot_image, roi_Ca_IR
	;stop
	if (n_elements(w_pixels) ge 10.) then begin
		trk2[*,j] = !VALUES.F_NAN
		;print, 'Bad Frames @:' +''+ string(j)
	endif else begin
		trk2[*,j] = trk(roi_Ca_IR[*,0:300], roi_iris[*,0:300])
		;print, 'Good Frames @:'+''+string(j)
	endelse
	print,string(13b)+'Calculating shifts:  % finished: ',float(j)*100./(nt_ca-1),format='(a,f4.0,$)'
endfor
;print, 'Shifts computed graciously!!'

bad_x = Where(finite(trk2[0,*]) eq 0, nbadx, COMPLEMENT=goodx, NCOMPLEMENT=ngoodx)
vector_x = reform(trk2[0,*])
if nbadx gt 0 && ngoodx gt 1 then vector_x[bad_x] = interpol(vector_x[goodx], goodx, bad_x)
bad_y = Where(finite(trk2[1,*]) eq 0, nbady, COMPLEMENT=goody, NCOMPLEMENT=ngoody)
vector_y = reform(trk2[1,*])
if nbady gt 0 && ngoody gt 1 then vector_y[bad_y] = interpol(vector_y[goody], goody, bad_y)

print,'Done with interpolation for the badframes with SAA'
window,2
cgplot, vector_x, xtitle='Frames', ytitle='Shifts'
window,3
cgplot, vector_y, xtitle='Frames', ytitle='Shifts'

for j=0, nt_ca-1 do begin
	t2 = blown_iris[*,*,j]
	w2 = where(t2 le 0.)
	t2[w2]=0.
	shifted_iris = red_rotation(t2,0,-vector_x[j],-vector_y[j])
	cropped_iris[*,*,j] = shifted_iris[514:1745,443:1661]; this is to get the IRIS aligned to SST. It possible that the pointing of SST is off by 0.4 degrees
	print,string(13b)+'Compensating the X Y offsets:  % finished: ',float(j)*100./(nt_ca-1),format='(a,f4.0,$)'
endfor

window,0
plot_image, rot(reform(dat_nb_ca[*,*,8,0,10]),-0.4)
window,1
plot_image, reform(cropped_iris[*,*,10])
save, vector_x, vector_y, filename=dpath_nb+'Mg_SJI_XY_offsets.sav'
;lp_write, cropped_iris, dpath_nb+'IRIS_Mg_SJI_2020-09-19T08:23:27_scans=0-149_aligned_8542.fcube' ;saving in native lapalma file format
;blink,[0,1]
end

