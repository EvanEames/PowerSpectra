PROGRAM PowerSpectra

!This is a program to find the Power Spectrum of a lightcone at various redshifts. It takes a section of the lightcone at the specified redshift and with a sepcified width. It then adjusts the section so all pixels are equal comoving distance. This cube is Fourier Transformed using the MKL 3DFFT. Lastly spherical shells are used to calculate the power spectrum, which is output. A FITS file can also be created of various stages of the cube (Section 5).

!BEFORE USING SEE SECTION 7 FOR USAGE NOTES

!====================================================================================================
!SECTION 0 - INDEX
!====================================================================================================
! 1 - Preamble and Formalities
!	1.1 - Set Parameters
!	1.2 - Initiate Other Variables
!	1.3 - Setting the File Names and the Input Directory
!	1.4 - Define Constants
!	1.5 - Read in the Files and Initiate a Loop Over Them
!	1.6 - Read in Temperature Lightcone
!	1.7 - Clip Errant Tb Values and Find Any Empty Cells
!	1.8 - Initiate a loop over the desired redshifts
!
! 2 - Preparing the Cube
!	2.1 - Calculate the Comoving Distance to Each Slice
!	2.2 - Find slices that fall within the specified distance (L_ps) of the target redshift (z_ps)
!	2.3 - Reshape the array so that comoving distances are constant
!	2.4 - Calculate the variance of the reshaped cube
!
! 3 - Fourier Transform
!	3.1 - Fourier Transform the Cube
!	3.2 - Organize the Fourier Transform Output
!
! 4 - Creating the Power Spectrum
!	4.1 - Set up radii at which spherical shells will be placed to calculate the power spectrum
!	4.2 - Figure out which radii bin each pixel falls into and add its value to that bin's sum
!	4.3 - Write out the Power Spectrum
!
! 5 - Preparing the FITS file
!	5.1 - Choose what to output
!
! 6 - Writing Out the Multi-Redshift Power Spectrum for Each File
!
! 7 - Usage Notes

Use MKL_DFTI
implicit none

!====================================================================================================
!SECTION 1 - PREAMBLE AND FORMALITIES
!====================================================================================================

!§1.1 - SET PARAMETERS
integer, parameter :: grid_size = 1024, nz = grid_size*8 !Size of the lightcone (physical and redshift, which is assumed to be 8x physical)
real, parameter :: z_i=15,z_f=6 !Inital and final redshifts of the cube
real, parameter :: cut_off = 50 !If clip = true, then this is the mK at which values will be clipped (above this they are considered errant).
real, parameter :: dr = grid_size/2 !Number of Spherical shells used in construcing the power spectrum. Default = grid_size/2.
real, parameter :: Mpc = 200 !Physical slice of the size in Mpc*h^(-1)
real(KIND=4), dimension(21), parameter :: z_list = (/7., 7.25, 7.5, 7.75, 8., 8.25, 8.5, 8.75, 9., 9.25, 9.5, 9.75, 10., 10.25, 10.5, 10.75, 11., 11.25, 11.5, 11.75, 12./) !Use this to calculate the PS at multiple z if many_z = true (also change dimension)
real, parameter :: L_ps = 200 !Thickness of the cube in which the power spectrum will be calculated (Mpc/h)
logical, parameter :: clip = .true. !Clip any errant Tb values (as specified by cut-off) and look for empty cells
logical, parameter :: variance = .false. !Calculate the variance of the cube before performing the FT (as a test of the PS)
logical, parameter :: make_fits = .false. !Output a fits file in 2D or 3D of the cube before or after being Fourier Transformed
logical, parameter :: chatty = .true. !False suppresses all message output
logical, parameter :: print_bins = .false. !True prints information about the spherical bins used to created the power spectrum
logical, parameter :: many_z = .true. !If true then the power spectrum will be calculated for all z specified in z_list
logical, parameter :: many_files = .false. !If true then runs for every file specified in the file PS_files.dat
logical, parameter :: equal_volume = .false. !If true the spherical bins used when calculating the power spectrum will be equal volume.
logical, parameter :: linear_k = .true. !If true the spherical bins used when calculating the power spectrum will have equal width in k_space. This one seems to work the best fyi.
logical, parameter :: log_bins = .false. !If true the spherical bins used when calculating the power spectrum will be logarithmic

!§1.2 - INITIATE OTHER VARIABLES
real, parameter :: a_i=1./(1.+z_i),a_f=1./(1.+z_f)
real, parameter :: da=(a_f-a_i)/(nz-1)
integer(KIND=4) :: control_block
integer, parameter :: WP = selected_real_kind(6,37)
real(KIND=4), dimension(grid_size,grid_size,nz) :: tk
real(KIND=4), dimension(2,dr,size(z_list)) :: PS_data
real(WP), dimension(grid_size,grid_size,grid_size) :: ps_cube, tmp_cube, reshaped_cube
real(KIND=4), dimension(:,:), allocatable :: map
real(KIND=4), dimension(nz) :: coDist
real(KIND=8), dimension(dr) :: ps_bins, pixels_in_bins, ps_sums
complex(WP), dimension(grid_size/2+1,grid_size,grid_size) :: ps_complex
real(WP), dimension(grid_size/2,grid_size,grid_size) :: real_cube
integer(KIND=4), dimension(3) :: L
integer(KIND=4), dimension(4) :: rstrides, cstrides
integer :: i, j, k, f, m, ir, main_loop, file_loop, fits_status, unit, blocksize, bitpix, group, fpixel, naxes(3), nelements, slice_ps, a, b, xi, xj, xk, xf, NaN_pixels, a_file, b_file
real(kind=4) :: exp_fact, rad, D_c, dx, dz, pi, c, H_0, O_M, O_lambda, lambda_21, z, z_tmp, fi, ff, dist, avg, stdev, var, z_ps, k_min, k_max
integer :: Status, unitary_cells, middle_pixels, fitschoice
character filename*80, filename2*80, main_file*160, in_dir*80, out_file*80
logical simple, extend, map3D
Character (len=3), dimension(5), parameter :: fx = (/'0.1','0.3','1','3','10'/)
Character (len=3), dimension(3), parameter :: HXR = (/'0','0.5','1'/)
Character (len=3), dimension(3), parameter :: xae = (/'0.5','1','2'/)
Character (len=3), dimension(3), parameter :: direction = (/'x','y','z'/)
Character (len=80), dimension(size(fx)*size(HXR)*size(xae)*size(direction)) :: PS_files
type(DFTI_DESCRIPTOR), Pointer :: plan_forwards, plan_backwards
L = [grid_size, grid_size, grid_size]

!§1.3 - SETTING THE FILE NAMES AND THE INPUT DIRECTORY
main_file = 'fx=1_HXR=0.5_xae=1_y_lightcone_dtb_fullres_onebigquasar.dat' !The file whose Power Spectrum you want (if many_file = .false.)
in_dir = '../../fbolgar/lightcone_suite/for_power_spectrum/' !The directory where the file (or files) you wish to  are stored
out_file = 'fx=1_HXR=0.5_xae=1_y_onebigquasar_PS.dat' !The filename in which you'd like to store the file (if many_file = .false.)
main_file = TRIM(in_dir) // TRIM(main_file)


!§1.4 - DEFINE CONSTANTS
 pi = 3.14
 c = 3E8 !km/s
 H_0 = 67.80 !(km/s)/Mpc
 O_M = 0.316
 O_lambda = 0.684
 lambda_21 = 0.21 !m
 dx = L_ps / grid_size !Thickness of each pixel for the power spectrum cube (Mpc/h)
 z_ps = 9 !Redshift at which the power spectrum will be calculated (should be between z_i and z_f) if many_z = False

!§1.5 - READ IN ALL THE FILES AND INITIATE A LOOP OVER THEM
if (many_files) then
	m = 1
	do i = 1, size(fx)
	do j = 1, size(HXR)
	do k = 1, size(xae)
	do f = 1, size(direction)
		PS_files(m) = 'fx=' // TRIM(fx(i)) // '_HXR=' // TRIM(HXR(j)) // '_xae=' // TRIM(xae(k)) // '_' // TRIM(direction(f)) // '_lightcone_dtb_fullres.dat'
		m = m + 1
	enddo
	enddo
	enddo
	enddo
endif
do file_loop = 1, size(PS_files)
	if (many_files) main_file = TRIM(PS_files(file_loop))
	if (chatty .and. many_files) print*, 'Treating file ', file_loop, ' of ', size(PS_files), ': ',main_file


!§1.6 - READ IN TEMPERATURE LIGHTCONE
if (chatty) print*,'Reading temperature lightcone'
if (many_files) open(23,file=TRIM(in_dir) // TRIM(main_file),status='old',form='unformatted')
if (.not. many_files) open(23,file=main_file,status='old',form='unformatted')
	read(23) tk
 close(23)


!§1.7 - CLIP ERRANT Tb VALUES AND FIND ANY EMPTY CELLS
if (clip) then
	if (chatty) print*,"Clipping errant Tb values and fixing empty cells in the lightcone..."
	NaN_pixels = 0
	do i = 1, nz
		do j = 1, grid_size
			do k = 1, grid_size
				if (tk(k,j,i) > cut_off) then
					tk(k,j,i) = cut_off
				endif
				if (isnan(tk(k,j,i))) then
					!print*,'NaN at',i,j,k
					tk(k,j,i) = 0
					NaN_pixels = NaN_pixels + 1
				endif
			enddo
		enddo
	enddo
	print*,"There were ",NaN_pixels," NaN pixels."
endif

!§1.8 - INITIATE A LOOP OVER THE DESIRED REDSHIFTS
do main_loop = 1, size(z_list)
if (many_z) z_ps = z_list(main_loop)
print*,"Calculating Power Spectra at redshift ",z_ps," (",main_loop,"/",size(z_list),")"

!====================================================================================================
!SECTION 2 - PREPARING THE CUBE
!====================================================================================================
!§2.1 - CALCULATE THE COMOVING DISTANCE OF EACH SLICE
if (chatty) print*,"Calculating the comoving distance to each slice..."
dz = z_f/1000.
D_c = 0
do i = 0,1000
	z = dz*i
	D_c = D_c + dz*(1/sqrt(O_M*(1+(z+dz))**3+O_lambda) + 1/sqrt(O_M*(1+z)**3+O_lambda))/2
enddo
D_c = D_c*(3E5/H_0)
coDist(nz) = D_c
z_tmp = z_f

slice_ps = 1
do k = 1,nz-1
	z = 1/(a_f - da*k) - 1
	if (z > z_ps) then
		slice_ps = slice_ps + 1
	endif
	D_c = D_c + (3E5/H_0)*(z-z_tmp)*(1/sqrt(O_M*(1+z)**3+O_lambda) + 1/sqrt(O_M*(1+z_tmp)**3+O_lambda))/2
	coDist(nz-k) = D_c
	z_tmp = z
enddo

!§2.2 - FIND THE SLICES THAT FALL WITHIN L_ps/2 OF z_ps
a = nz
b = 1
i = 1
do while (coDist(slice_ps) - coDist(a) > L_ps/2) 
	a = a - 1
enddo
do while (coDist(b) - coDist(slice_ps) > L_ps/2)
	b = b + 1
enddo
if (chatty) print*,'The power spectrum will be calculated between pixels ',a,' and ',b,' which corresponds to redshift ',z_ps,' and width ',L_ps,' MHz/h.'

!§2.3 - RESHAPE THE ARRAY SO COMOVING DISTANCES ARE CONSTANT (AND SIZE IS GRID_SIZE**3)
if (chatty) print*,"Rescaling array so comoving distances are constant..."
xi = a
xf = a
do i = 1, grid_size
	!Find the tk pixel in which the beginning of the ith cell in the ps_cube will fall into
	do while ((i-1)*dx > coDist(xi)-coDist(a))
		xi = xi - 1
	enddo
	!Find the tk pixel in which the end of the ith pixel in the ps_cube will fall into
	do while (i*dx > coDist(xf)-coDist(a))
		xf = xf - 1
	enddo
	if (xi == xf) then
		!If the beginning and end fall into the same tk pixel, then simply find the ratio of comoving distances, and scale the Tb accordingly
		ps_cube(:,:,grid_size - i + 1) = tk(:,:,xi)
	else
		!Otherwise, find the fraction of the first and last pixels in tk that will contribute to the ith ps_cube pixel
		fi = (coDist(xi)-coDist(a) - dx*(i-1))/dx
		ff = (dx*i-(coDist(xf+1)-coDist(a)))/dx
		ps_cube(:,:,grid_size - i + 1) = tk(:,:,xi)*fi + tk(:,:,xf)*ff
		middle_pixels = 1
		!Also add any intermediary pixels (assuming the ith ps_cube pixel spans more than 2 pixels in tk)
		do while (xi + middle_pixels < xf)
			ps_cube(:,:,grid_size - i + 1) = ps_cube(:,:,grid_size - i) + tk(:,:,xi+middle_pixels)
			middle_pixels = middle_pixels + 1
		enddo
	endif
enddo
reshaped_cube = ps_cube

!§2.4 - CALCULATE THE VARIANCE OF THE RESHAPED CUBE (IF DESIRED)
if (variance) then
	avg = sum(ps_cube)/grid_size**3
	var = 0
	do i = 1, grid_size
		do j = 1, grid_size
			do k = 1, grid_size
				stdev = (ps_cube(k,j,i) - avg)**2
				var = var + stdev
			enddo
		enddo
	enddo
	var = var/grid_size**3
	print*,'variance = ', var
endif

!====================================================================================================
!SECTION 3 - FOURIER TRANSFORM
!====================================================================================================
!§3.1 - CALCULTE THE FOURIER TRANSFORM OF THE CUBE
print*,"Performing 3D Fourier Transform..."
status = DftiCreateDescriptor(plan_forwards, DFTI_SINGLE, DFTI_REAL, 3, L)
status = DftiSetValue(plan_forwards, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
status = DftiSetValue(plan_forwards, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
	!MKL doesn't set up the required data spacing for 3D FT, so it has to be done manually
	cstrides = [0, 1, INT(L(1)/2.0)+1, L(2)*(INT(L(1)/2.0)+1)]
	rstrides = [0, 1, L(1), L(2)*L(1)]
status = DftiSetValue(plan_forwards, DFTI_INPUT_STRIDES, rstrides)
status = DftiSetValue(plan_forwards, DFTI_OUTPUT_STRIDES, cstrides)
status = DftiCommitDescriptor(plan_forwards)
status = DftiComputeForward(plan_forwards, ps_cube(:,1,1), ps_complex(:,1,1))
status = DftiFreeDescriptor(plan_forwards)

!§3.2 - ORGANIZE THE FOURIER TRANSFORM OUTPUT
print*,"Rearranging Fourier Parameter Space..."
do i = 1, grid_size
	do j = 1, grid_size
		do k = 1, grid_size/2
			!Take the modulus of each value and normalize (divide by cube volume)
			real_cube(k, j, i) = sqrt(real(ps_complex(k,j,i))**2 + aimag(ps_complex(k,j,i))**2)/sqrt(real(grid_size**3))*Mpc**3
		enddo
	enddo
enddo

!====================================================================================================
!SECTION 4 - CREATING THE POWER SPECTRUM
!====================================================================================================
!§4.1 - SET UP THE RADII AT WHICH SPHERICAL SHELLS WILL BE TAKEN FOR THE POWER SPECTRUM
print*,"Creating Power Spectrum..."
k_min = 2*pi/Mpc
k_max = k_min*(grid_size/2)
do i = 1, dr
	if (equal_volume) then
		ps_bins(i) = real(grid_size/2)*(i/dr)**(1./3.) !Choose the bin limits (equal volume spacing)
		ps_bins(i) = ps_bins(i)*k_min !This switches from physical pixels to k-space
	else if (linear_k) then
		ps_bins(i) = (real(i)+0.5)*k_min !Choose bin limits. Here we assume even spacing in k-space
	else if (log_bins) then
		ps_bins(i) = 2.*(real(grid_size/2)/2.)**(i/dr) !Choose the bin limits (logarithmic spacing)
		ps_bins(i) = ps_bins(i)*k_min !This switches from physical pixels to k-space
	else
		print*,"Please make one of the following true: equal_volume, linear_k, or log_bins"
	end if
	ps_sums(i) = 0
	pixels_in_bins(i) = 0
enddo
!§4.2 - FIGURE OUT WHICH RADII BIN EACH PIXEL FALLS INTO AND ADD ITS VALUE TO THE SUM
do i = 1, grid_size
	do j = 1, grid_size
		do k = 1, grid_size
			!Because the FFT only outputs a half cube, with the most bright points at the four k=1 corners,
			!we have to take advantage of some symmetries to properly decide how far a pixel is from the centre
			if (i <= grid_size/2) xi = i-1
			if (i > grid_size/2) xi = grid_size - i + 1
			if (j <= grid_size/2) xj = j-1
			if (j > grid_size/2) xj = grid_size - j + 1
			if (k <= grid_size/2) xk = k-1
			if (k > grid_size/2) xk = grid_size - k + 1
			!Calculate the distance to the pixel and decide which bin it falls into
			dist = sqrt(real(xi**2 + xj**2 + xk**2))
			dist = dist*(2*pi/Mpc)
			if (dist < real(grid_size/2)*k_min) then
				ir = int((dist-0.5*k_min)/k_min)+1
				ps_sums(ir) = ps_sums(ir) + real_cube(xk + 1,j,i)**2
				pixels_in_bins(ir) = pixels_in_bins(ir) + 1
			endif
		enddo
	enddo
enddo
!The central pixel tends to be errant, so we remove it
ps_sums(1) = ps_sums(1) - real_cube(1,1,1)**2
pixels_in_bins(1) = pixels_in_bins(1) - 1
if (print_bins) then
	print*,'The bin limits (in k_space) are:'
	print*,ps_bins(1:10)
	print*,'The number of pixels that fall into each bin are as follows:'
	print*,pixels_in_bins(1:10)
	print*,'The sum of Tb for all pixels in each bin are as follows:'
	print*,ps_sums(1:10)
endif
do i = 1, dr
	!Average all the radii bins
	ps_sums(i) = ps_sums(i)/pixels_in_bins(i)
enddo
if (print_bins) then 
	print*,'The average Tb in each bin are as follows:'
	print*,ps_sums(1:10)
endif

!§4.3 - WRITE OUT THE POWER SPECTRA
if (chatty) print*,"Writing Power Spectrum..."
if (.not. many_z) then
	if (z_ps<10) then
		write (filename2, "(A15,F3.1,A4)") "PowerSpectra_z=", z_ps, '.dat'
	else
		write (filename2, "(A15,F4.1,A4)") "PowerSpectra_z=", z_ps, '.dat'
	endif
	open(20,file=filename2,status='replace',form='formatted')
	do i = 1, dr
		if (i == 1) then
			write(20,*)	(ps_bins(i)*(2./3.)), &
	                     	 	ps_sums(i)*(k_min)**3/(2*pi**2)/Mpc**3
		else
			write(20,*)	((ps_bins(i) - ps_bins(i-1))/2 + ps_bins(i-1)), &
					ps_sums(i)*((ps_bins(i) - ps_bins(i-1))/2 + ps_bins(i-1))**3/(2*pi**2)/Mpc**3
		endif
	enddo
	close(20)
else
	do i = 1,dr
		if (i==1) ps_data(:,i,main_loop) = (/ps_bins(i)*(2./3.), ps_sums(i)*(ps_bins(i)/2)**3/(2*pi**2)/)
		if (i/=1) ps_data(:,i,main_loop) = (/((ps_bins(i) - ps_bins(i-1))/2 + ps_bins(i-1)), ps_sums(i)*((ps_bins(i) - ps_bins(i-1))/2 + ps_bins(i-1))**3/(2*pi**2)/)
	enddo
endif


!====================================================================================================
!SECTION 5 - PREPARING THE FITS FILE
!====================================================================================================
if (make_fits) then
	if (chatty) print*,'Generating FITS file...'

!§5.1 - CHOOSE WHAT TO OUTPUT IN THE FITS FILE
	fitschoice = 3 !1 = initial section of Tb map, 2 = resized Tb map, 3 = Tb map FT
	map3D = .true. !Create a 3D fits file if true. Otherwise fits will be a 2D slice

	!Take a slice of what is requested
	if (.not. map3D) then
		if (fitschoice == 1) then
			allocate(map (a-b, grid_size))
			do i = 1, a-b
				do j = 1, grid_size
					map(i,j) = tk(j,grid_size/2,b+i-1)
				enddo
			enddo
		else
			allocate(map (grid_size,grid_size))
			if (fitschoice == 2) then
				do i = 1, grid_size
					do j = 1, grid_size
						map(i,j) = reshaped_cube(j,grid_size/2,i)
					enddo
				enddo
			else
				do i = 1, grid_size/2
					do j = 1, grid_size
						map(i,j) = 0
						map(i+grid_size/2,j) = real_cube(i,grid_size/2,j)
					enddo
				enddo
			endif
		endif
	endif

	!Name of the FITS file to be created:
	filename = 'ReshapedCube.fits'

	!Get an unused Logical Unit Number to create the FITS file
	fits_status = 0
	call ftgiou(unit,fits_status)

	!create the new empty FITS file
	blocksize = 1
	call ftinit(unit,filename,blocksize,fits_status)

	!initialize parameters about the FITS image 
	simple = .true.
	bitpix = -32
	naxes(1) = grid_size
	naxes(2) = grid_size
	naxes(3) = grid_size
	if (fitschoice == 1) then
		if (map3D) then
			naxes(3) = a-b
		else
			naxes(1) = a-b
		endif
	else if (fitschoice == 3) then
		naxes(1) = grid_size/2
	endif
	extend = .true.

	!write the required header keywords
	if (map3D) then
		call ftphpr(unit,simple,bitpix,3,naxes,0,1,extend,fits_status)
	else
		call ftphpr(unit,simple,bitpix,2,naxes(1:2),0,1,extend,fits_status)
	endif

	!write the array to the FITS file
	group = 1
	fpixel = 1
	if (map3D) then
		nelements = naxes(1)*naxes(2)*naxes(3)
		if (fitschoice == 1) call ftppre(unit,group,fpixel,nelements,tk(:,:,b:a),fits_status)
		if (fitschoice == 2) call ftppre(unit,group,fpixel,nelements,reshaped_cube,fits_status)
		if (fitschoice == 3) call ftppre(unit,group,fpixel,nelements,real_cube,fits_status)
	else
		nelements = naxes(1)*naxes(2)
		call ftppre(unit,group,fpixel,nelements,map,fits_status)
	endif

	!close the file and free the unit number
	call ftclos(unit, fits_status)
	call ftfiou(unit, fits_status)
endif

if (.not. many_z) stop
enddo !This closes the loop over the different files (see §1.7)

if (chatty) print*,'FINISHED'

!====================================================================================================
!SECTION 6 - WRITING OUT THE MULTI-REDSHIFT POWER SPECTRUM FOR EACH FILE
!====================================================================================================
b_file = SCAN(main_file,'l', .False.)
if (many_files) filename2 = 'lightcone_suite_PS/' // main_file(1:b_file-1) // 'PS.dat'
if (.not. many_files) filename2 = out_file
if (chatty) print*,'Writing final power spectrum file to ',filename2
open(20,file=filename2,status='replace',form='formatted')
do i = 1, size(z_list)
	do j = 1, dr
		write(20,*) z_list(i), ps_data(1,j,i), ps_data(2,j,i)
	enddo
	write(20,*)
enddo
 close(20)

if (.not. many_files) stop
enddo !This closes the loop over the different files (see §1.4)

!====================================================================================================
!SECTION 7 - USE NOTES
!====================================================================================================

!The compile command looks like this:
!ifort PowerSpectra.f90 -o PowerSpectra.out -mcmodel=large -mkl -L/usr/lib64 -L/scratch/semelin/cfitsio -lcfitsio -lnsl
!Change the directories as needed to find lib64 and cfitsio
!Be sure that you have the following in the same directory as this file to assure the MKL FT works:
!	mkl_dfti.f90
!	mkl_dfti.h
!	mkl_dfti.mod
!	mkl_dfti.o
!	mkl_dft_type.mod
!
!To start, choose the input lightcone in §1.3, and then adjust the parameters as desired in §1.1. If you want to calculate the PS at just one redshift, make sure the many_z variable in §1.1 is False, and set the redshift you choose by adjusting the z_ps variable in §1.3. Otherwise, you can set many_z to True, and give a list of redshifts in the z_list parameter (also §1.1).
!If many_z = .false.
!	The program will output a file called 'PowerSpectra_z=#.dat', and will automatically replace # with the redshift at which the FT has been performed.
!	With Gnuplot, you can then plot the PS as follows:
!		set logscale
!		set xrange [0.2:6]
!		plot 'PowerSpectra_z=10.0.dat' u 1:2 with lines lw 2 title "z = 10"
!
!If many_z = .true.
!	The program will output a file whose name will be whatever you set in §1.3
!	With gnuplot, you can ghen plot the PS as follows:
!		splot 'YOUR_FILE.dat' u 1:2:3
!		
!If many_z = .true. AND many_files = .true.
!	The program will output a file called 'fx=#_HXR=#_xae=#_x_PS.dat' (where x could also be y or z). You can plot each individual file with gnuplot using the same splot command
!
!If the FT isn't working it is possible that there are NaN values in the lightcone. This can be fixed by turning on the 'clip' option, which will fill any NaN values with 0.


END program
