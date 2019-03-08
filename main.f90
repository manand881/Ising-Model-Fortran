	program ising  
	implicit  none
	! The input parameters for this program are in "ising.in"
	! Variable declarations:
	integer :: i, j, m, n, m2, n2,i1	! dummy integers

	integer, allocatable :: A(:,:)  ! matrix containing spins
	integer :: nrows, ncols	! number  of rows  and  cols of A
	real :: temp, beta	! temperature, inverse temperature
	integer :: ConfigType	! starting configuration type
	integer  :: npass	! number of passes for MC algorithm
	integer :: ipass	! the current pass number
	integer  :: nequil	! number of equilibration steps
	integer :: trial_spin	! values of changed spin
	real :: high_temp	! starting temp for scan
	real :: low_temp	! final temp for scan
	real :: temp_interval	! interval between scan points
	integer :: nscans	! number of scans (each at diff T)
	integer :: iscan,iscan1	! current scan number
	real :: deltaU	! change in energy between 2 configs
	real :: deltaU1, deltaU2	! energy changes for lattice gas
	real :: log_eta	! log of random number to compare to
	real :: magnetization	! magnetization of all spins in lattice
	real :: magnetization_ave	! cumulative average magnetization
	real :: magnetization2_ave	! cumulative average of mag. squared
	real :: energy	! energy of all spins in lattice
	real :: energy_ave	! cumulative average of energy
	real :: energy2_ave	! cumulative average of energy squared
	integer :: output_count	! # times things have been added to averages

	real*8 :: ran0, rand_uniform

	print*,  "	MONTE CARLO 2D ISING MODEL	"
	print*, "Monte Carlo Statistics for 2D Ising Model with"
	print*, "  periodic boundary conditions."
	print*, "The critical temperature is approximately 2.3, as seen on Chandler p. 123."

	! Read in input parameters from file "ising.in"
	open(unit=11,file='ising.in',status='old', action='read')
	read(11,*);
	read(11,*) nrows
	read(11,*);
	read(11,*) ncols 
	read(11,*);
	read(11,*) npass
	read(11,*);
	read(11,*) nequil
	read(11,*);
	read(11,*) high_temp
	read(11,*);
	read(11,*) low_temp
	read(11,*);
	read(11,*) temp_interval
	read(11,*);
	read(11,*) ConfigType
	
	close(unit=11)
	 
	! Set the dimensions of the matrix of spin arrays. This program uses
	! periodic boundary conditions, so the first two rows and columns are
	! the same as the last two.
	allocate(A(nrows+2,ncols+2))
	
	! Open output files:
	open(unit=32,file='spin-array',status='replace',action='write')
	write(32,*) nrows
	write(32,*) ncols
	nscans  =  int((high_temp  -  low_temp)/temp_interval)  +  1
	
	open(unit=33,file='magnetization',status='replace',action='write')
	write(33,*) "temp	ave_magnetization	ave_magnetization^2    susceptibility"
	open(unit=34,file='energy',status='replace',action='write')
	write(34,*)   "temp	ave_energy	ave_energy^2	C_v"
	scan_loop: do iscan = 1, nscans
	temp = high_temp - temp_interval*(iscan-1)
	print*, "Running program for T =", temp
		 
	! Initialize variables
	beta  =  1.0/temp
	output_count   =   0
	energy_ave  =  0.0
	energy2_ave  =  0.0
	magnetization_ave  =  0.0
	magnetization2_ave  =  0.0
	
	! Set up the initial spin configuration.
	select case(ConfigType)
	case(1)  ! checkerboard setup
	A(1,1)  =  1
	do i = 1, nrows+1
	A(i+1,1)  =  -A(i,1)
	enddo
	do j = 1, ncols+1
	A(:,j+1)  =  -A(:,j)
	enddo
	case default
	print*, "Error! Check ConfigType parameter in ising.in"
	stop
	end select
	

	
	! Main loop containing Monte Carlo algorithm:
	MC_passes: do ipass = 0, npass
	
	
	! If ipass is greater than nequil (the number of equilibration steps),
	! calculate the magnetization and energy:
	if (ipass > nequil) then
	output_count  =  output_count  +  1
	magnetization     =     sum(A(2:nrows+1,2:nrows+1))/(ncols*nrows*1.0)
	magnetization_ave  =  magnetization_ave  +  magnetization
	magnetization2_ave   =   magnetization2_ave   +   magnetization**2
	energy  =  0.0
	do i = 2, nrows + 1
	do j = 2, ncols + 1
	energy = energy - A(m,n)*(A(m-1,n)+A(m+1,n)+A(m,n-1)+A(m,n+1))
	enddo
	enddo
	
	! Divide the energy by the total number of spins to get the ave
	! energy per spin, and divide by 2 to account for double counting.
	energy  =  energy/(ncols*nrows*2.0)
	energy_ave  =  energy_ave  +  energy
	energy2_ave  =  energy2_ave  +  energy**2
	endif
	
	! Randomly choose a spin to change:
	ran0=rand_uniform()
	m  = nint((nrows-1)*ran0 + 2)   ! choose a random row
	ran0=rand_uniform()
	n = nint((ncols-1)*ran0 + 2)   ! choose a random column
	trial_spin = -A(m,n)	! trial spin value
	! Find change in energy (deltaU) due to trial move.
	! If exp(-beta*deltaU) > eta, where eta is random, accept move:
	deltaU    =    -trial_spin*(A(m-1,n)+A(m+1,n)+A(m,n-1)+A(m,n+1))*2
	ran0=rand_uniform()
	log_eta = dlog(ran0 + 1.0d-10)   ! random number 0-1 (+ tiny offset)
	if (-beta*deltaU  >  log_eta)  then
	A(m,n)  =  trial_spin
	if (m == 2) A(nrows+2,n) = trial_spin
	if (m == nrows+1) A(1,n) = trial_spin
	if (n == 2) A(m,ncols+2) = trial_spin
	if (n == ncols+1) A(m,1) = trial_spin
	endif
	enddo MC_passes
	
	
	write(33,*)   temp,   abs(magnetization_ave/output_count),   & 
	magnetization2_ave/output_count,         &
	beta*(magnetization2_ave/output_count - (magnetization_ave/output_count)**2)
	write(34,*)    temp,    energy_ave/output_count,    energy2_ave/output_count,    &
	(beta**2)*(energy2_ave/output_count - (energy_ave/output_count)**2)
	enddo  scan_loop
	close(32)
	close(33)
	close(34)
	open(40,file='temperature.dat',status= 'replace') 
	write(40,*)    "temperature", "	i" ,"	 j","	spin A(i,j)"
	scan_loop1:do iscan1 = 1, nscans
	do i = 2, nrows + 1
	do j = 2, ncols + 1
        temp = high_temp - temp_interval*(iscan1-1)
	print*,  temp, i , j, A(i,j)
	write(40,*)   temp, i , j, A(i,j)
	enddo 
	enddo
	enddo scan_loop1
	close(40) 
	print*, "Program ising.f90 complete!"

	end program ising


!!########### UNIFORM RANDOM NUMBER GENERATOR ############

	function rand_uniform()

	  IMPLICIT REAL*8(A-H,O-Z)

! ----- variables for portable seed setting -----
	  INTEGER, DIMENSION(:), ALLOCATABLE :: ia_seed
	  INTEGER, DIMENSION(1:8) :: idt_seed
! ----- end of variables for seed setting -----

! ----- Set up random seed portably -----
	  CALL RANDOM_SEED(size=i_seed)
	  ALLOCATE(ia_seed(1:i_seed))
	  CALL RANDOM_SEED(get=ia_seed)
	  CALL DATE_AND_TIME(values=idt_seed)
	  ia_seed(i_seed)=idt_seed(8) 
	  ia_seed(1)=idt_seed(8)*idt_seed(7)*idt_seed(6)
	  CALL RANDOM_SEED(put=ia_seed)
	  DEALLOCATE(ia_seed)
! ----- Done setting up random seed -----

	  CALL RANDOM_NUMBER(r)
	  rand_uniform = r
	  
	  return
 	  
	end function rand_uniform




