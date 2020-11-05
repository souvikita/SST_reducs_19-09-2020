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
;######################## Reading the precomputed shifts ##########
print, 'Reading the X Y offsets (vector_x, vector_y) from the Mg SJI alignment: Please wait'
restore, '/mn/stornext/d18/lapalma/reduc/2020/2020-09-19/CRISP/cubes_nb/Mg_SJI_XY_offsets.sav'
;######################## New SST cubes already have times in second #########
dpath_wb = '/mn/stornext/d18/lapalma/reduc/2020/2020-09-19/CRISP/cubes_wb/'
red_fitscube_getwcs, dpath_wb+'wb_8542_2020-09-19T08:23:27_scans=0-149_corrected_im.fits', coordinates = coordinates
time_secs_ca = coordinates.time[0,0,0,*] ; This is already in seconds from the headers.
nx_ca = 1232
ny_ca = 1219
nw_ca = 9
nt_ca = 150
dpath_nb = '/mn/stornext/d18/lapalma/reduc/2020/2020-09-19/CRISP/cubes_nb/'
dpath_l2_iris ='/mn/stornext/d10/HDC2/iris/data/level2_decompressed/2020/09/19/20200919_065936_3633109417/'

f1= file_search(dpath_l2_iris+'*_SJI_*.fits')
read_iris_l2, f1[0], index, data ;Reading IRIS Si SJI
time_iris = index[*].date_obs ;Entire time series of the dataset
data_Si = iris_dustbuster(index,data); removing the dust particles
time_secs_Si = time_iris_secs(time_iris)
blown_iris = fltarr(2084,2179,nt_ca)
cropped_iris = fltarr(nx_ca,ny_ca,nt_ca)
;########### Looping through #########

for j =0, nt_ca-1 do begin
	target_time_sst = time_secs_ca[j]
	minimizer = abs(time_secs_Si-target_time_sst)
	t1 = data_Si[*,*,where(minimizer eq min(minimizer))]
	blown_iris[*,*,j] = reform(congrid(t1,fix(371/0.178),fix(388/0.178),cubic=-0.5))
	t2 = blown_iris[*,*,j]
	w2 = where(t2 le 0.)
	t2[w2]=0.
	shifted_iris = red_rotation(t2,0,-(vector_x[j]+2.69),-(vector_y[j]+1.33))
	cropped_iris[*,*,j] = shifted_iris[514:1745,443:1661]
	print,string(13b)+'Compensating the X Y offsets:  % finished: ',$
	   float(j)*100./(nt_ca-1),format='(a,f4.0,$)'
endfor
;********The additional Offsets of (2.69,1.33) is because of misalignment between the Si $
; and Mg Slit jaws. 
;window,0
;plot_image, rot(reform(dat_nb_ca[*,*,8,0,10]),-0.4)
;window,1
;plot_image, reform(cropped_iris[*,*,10])
;blink, [0,1]
lp_write, cropped_iris, dpath_nb+'IRIS_Si_SJI_2020-09-19T08:23:27_scans=0-149_aligned_8542.fcube'
end
