PRO jcmtreduction

idnames = ['L108054','L109102','L109038', 'L110038', 'L109066', 'L108045', 'L109010', 'L101037', 'L109028']
idnums = [108054,109102,109038, 110038, 109066, 108045, 109010, 101037, 109028]
readcol, 'LOWMASS_catalog.dat', v1, v2, v3, v4, v5, v6, v7, v8, v9, v10
files = FILE_SEARCH('*.fits')
filesnew = files[where(STRPOS(files, '850') GT 0)]

!P.MULTI = [0, 1, 1]
DEVICE, RETAIN=3, DECOMPOSED=0
entry_device = !D.NAME & help, entry_device 
set_plot, 'ps' 
device, XSIZE=5, YSIZE=5, /INCHES, file="jcmt_proposal_upper_limits.eps"

x = [8.63, 9.08, 9.47, 9.99, 10.51 ]
y = [-2.07, -2.17, -2.3, -2.77, -3.01]
plot, x, y, yrange=[-3.4, -2.0], xtitle='Stellar Mass (Log M)', ytitle=TeXtoIDL("Log(M_{dust}/M_{Stellar})"), LINESTYLE=0.0, thick=4.0,  CHARTHICK=3.0, xthick=5, ythick=5, charsize=1.0
image_size = 21

for n=0, n_elements(idnames)-1 do begin
	imagenew=fltarr(image_size,image_size)
	print, idnames[n]
	
	for m=0, n_elements(filesnew)-1 do begin
		image = mrdfits(filesnew[m], 0, hdr, /silent)
		naxisi = sxpar(hdr,'OBJECT')
		exposure_time = sxpar(hdr,'EXP_TIME')
		units = sxpar(hdr,'BUNIT')
			
		if (strtrim(naxisi) EQ idnames[n]) then begin
		
			ra = v2[where(v1 EQ idnums[n])]
			dec = v3[where(v1 EQ idnums[n])]
			print, ra, dec
			d=size(image,/dimension)

			nx=d[0]

			ny=d[1]
			print, filesnew[m],"  ", idnames[n],"  ", exposure_time,"s","  ", units, "  ", d
			ADXY, hdr, ra, dec, x_center, y_center  ; coordinates of center in pixel coordinates
			print, x_center, y_center, 'centers of the image'

			imagenew =  image[x_center-((image_size-1.0)/2.0):x_center+((image_size-1.0)/2.0),x_center-((image_size-1.0)/2.0):x_center+((image_size-1.0)/2.0)] 
		
	ra = v2[where(v1 EQ idnums[n])]
	dec = v3[where(v1 EQ idnums[n])]
	redshift = v4[where(v1 EQ idnums[n])]
	
;########################Now to do the actual photometry##################################################
	d=size(imagenew,/dimension)

	nx=d[0]

	ny=d[1]
	galaxysize = 30.0          ;this is roughly in arcseconds the diameter of our galaxies
	print, nx, ny
	radius_in_pixels = 2.0
	x_center = ((image_size+1.0)/2.0)
	y_center = ((image_size+1.0)/2.0)
	
	radin=radius_in_pixels

	radout=3.0*radius_in_pixels
	areacircle = 3.14159*radin*radin
	areashell = (3.14159*radout*radout) - areacircle


	totalflux=0.0

	background=0.0

	for i=0, nx-1 do begin
		for j=0, ny-1 do begin
			if (sqrt((i-x_center)^2.0 + (j-y_center)^2.0) LE radout and sqrt((i-x_center)^2.0 + (j-y_center)^2.0) GT radin) then begin
				background = background + (imagenew[i, j])
			endif else begin 
				totalflux = totalflux
				background = background
			endelse
		endfor
	endfor

	skylevel = ((background)/(areashell))

	
	error = 0.0
	for i=0, nx-1 do begin
		for j=0, ny-1 do begin
			if (sqrt((i-x_center)^2.0 + (j-y_center)^2.0) LE radout and sqrt((i-x_center)^2.0 + (j-y_center)^2.0) GT radin) then begin
				error = error + (imagenew[i, j] - skylevel)^2.0
			endif else begin 
				error = error 
			endelse
		endfor
	endfor

	flux_error = sqrt(error/areashell)

	for i=0, nx-1 do begin
		for j=0, ny-1 do begin
			if (sqrt((i-x_center)^2.0 + (j-y_center)^2.0) LE radin) then begin
				totalflux = totalflux + ((imagenew[i, j]-skylevel))
				;print, imagenew[i,j], skylevel
			endif else begin 
				totalflux = totalflux
				background = background	
			endelse
		endfor
	endfor

	print, 'Sky level   ', skylevel
	print, 'Detection  ',  flux_error*3.0
	print, 'calculated rms equals  ', flux_error
	print, 'signal to noise', (totalflux/(flux_error*4.0))
	print, 'upper flux limit', flux_error*6.0
	
	wavelength = 850E-06
	c = 2.99E+08					;ms-1
	plank_constant = 6.626E-34      ;Js
	frequency = (c/wavelength)		;s-1
	boltzman_constant = 1.380E-23   ;JK-1
	dust_temperature = 20.0			;K
	absoprtion_coefficient = 0.07   ;m2/kg
	flux = flux_error*3.0*1E-29   ;Jm-2
	distance = lumdist(redshift)*(1E+06)*(3.0856E+16) ;m
	steradian = 1.0;(3.14159*((galaxysize/2.0)^2.0))/(4.25E+10)                ;(3.14159*galaxysize*galaxysize*4.848*4.848E-12)/4.0
	intermediate_freq_plank = (1E+22)*frequency*plank_constant
	intermdiate_dust_temp_boltzmann = dust_temperature*boltzman_constant*1E+22
	
	dust_upper_limit = (distance*distance*flux)/((2.0*absoprtion_coefficient*frequency*frequency*frequency*plank_constant/(c*c))*(1/(exp((intermediate_freq_plank/intermdiate_dust_temp_boltzmann))-1.0)))*(1/steradian)
	dust_upper_limit_solar_units = (dust_upper_limit)/(1.9891E+30)
	print, 'dust mass', alog10(dust_upper_limit_solar_units)
	print, 'stellar mass', v5[where(v1 EQ idnums[n])]
	print, 'dust_to_stars_ratio', alog10(dust_upper_limit_solar_units)-v5[where(v1 EQ idnums[n])]
	PLOTSYM, 1 ,3.0, /FILL, thick=2.0
	oplot, v5[where(v1 EQ idnums[n])], alog10(dust_upper_limit_solar_units)-v5[where(v1 EQ idnums[n])], psym=8.0, symsize=1.0
			
;#########################################################################################################
		endif 
	endfor
endfor

oplot, range(8.63, 9.75, 10), range(-2.9, -2.9, 10), LINESTYLE=2.0, thick=4.0
oplot, range(9.70, 10.51, 10), range(-3.3, -3.3, 10), LINESTYLE=2.0, thick=4.0

al_legend,[TeXtoIDL("Cortese et al. (2012) Scaling Relation"), TeXtoIDL("Our Lower Bounds")], linestyle=[0, 2], /bottom, /left, charsize=0.8, symsize= 0.7, charthick=3

device, /close
set_plot, entry_device

spawn, 'epstopdf jcmt_proposal_upper_limits.eps'

end




