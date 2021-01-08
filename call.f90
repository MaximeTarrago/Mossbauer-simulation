module procedures
!implicit none
IMPLICIT REAL*8 (A-H,O-Z)
contains
subroutine  vector_builder(orb_pop,nstate,vector_matrix,ndummy,nconfig,multiplicity)
! CALCULATES THE SPIN-EIGENSTATES FROM THE GIVEN ELECTRONIC CONFIGURATIONS DEFINED BY THE USER, LATER TO BE USED AS BASIS TO CONSTRUCT THE HAMILTONIAN MATRICES
! INPUT PARAMETERS:
!     ORB_POP (INPUT,INTEGER, 2-DIMENSIONAL [1:NSTATE,1:5]) : THE ELECTRONIC CONFIGURATIONS. IT IS A TABLE CONTAINING ONE SET OF INTEGERS FOR EACH ELECTRONIC CONFIGURATION.  
! THOSE INTEGERS RANGE BETWEEN 0 AND 2 FOR EACH OF THE 5 D-ORBITALS (DXY, DXZ, DYZ, DZ2, DX2-Y2 IN ORDER).
! 	NSTATE (INPUT, INTEGER) : THE NUMBER OF ELECTRONIC CONFIGURATIONS DEFINED BY THE USER. IT DEFINES THE LENGTH OF ONE DIMENSION OF VECTOR_MATRIX
!  	VECTOR_MATRIX (OUTPUT, INTEGER, 4-DIMENSION [1:NSTATE,1;MULTIPLICITY,1:NCONFIG,1:10]) : THE SPIN-EIGENSTATES. IT TAKES THE FORM OF SEVERAL SET OF 10 INTEGER OF VALUE 0 OR 1, EACH CORRESPONDING 
! TO A SPINORBITAL POPULATION. EACH SET OF 10 SPINORBITAL POPULATIONS AMOUNTS TO A SPINORBITAL CONFIGURATION.  THERE ARE SEVERAL SUCH SPINORBITAL POPULATIONS FOR A GIVEN MS, 
! AND SEVERAL MS FOR A GIVEN MULTIPLICITY/ELECTRONIC CONFIGURATION, ENDING UP IN A 4-DIMENSIONAL TABLE.
! THE FIRST DIMENSION CORRESPONDS TO THE CHOSEN EL. CONFIGURATIONS, THE SECOND CORRESPONDS TO THE MS IN A GIVEN CONFIGURATION, THE THIRD TO THE SPINORBITAL CONFIGURATIONS IN A GIVEN MS, AND THE FOURTH 
! THE SPINORBITAL POPULATIONS IN EACH SPINORBITAL CONFIGURATION.
!	NDUMMY (OUTPUT, INTEGER, 1-DIMENSIONAL [1:NSTATE]) : STORES THE MULTIPLICITY OF MS FOR EACH ELECTRONIC CONFIGURATION (TAKING THE MAXIMUM S VALUE). IS USED IN LATER SUBROUTINES    
!	NCONFIG (OUTPUT,INTEGER) : STORES THE MAXIMUM NUMBER OF SPINORBITAL CONFIGURATIONS ENTERING IN THE COMPOSITION OF A SPIN-EIGENSTATES. IT DEFINES THE LENGTH OF ONE DIMENSION OF VECTOR_MATRIX.
! IS USED IN LATER SUBROUTINES.
! MULTIPLICITY (OUTPUT, INTEGER) : STORES THE MAXIMUM NUMBER OF MS EIGENSTATES ENTERING IN A MULTIPLICITY. IT DEFINES THE LENGTH OF ONE DIMENSION OF VECTOR_MATRIX.

implicit none
integer :: nstate,i,j,multiplicity,nconfig,k,l,aa,b,c,d,e,X,Ms
integer, dimension (1:nstate) :: ndummy
integer, dimension (1:nstate,1:5), intent(in) :: orb_pop
integer, dimension (1:6) :: ai, af,a,Y
integer, allocatable, dimension(:,:,:,:) :: vector_matrix
integer, allocatable, dimension (:,:) :: config_sample
!do i=1,nstate
!         print*, orb_pop(i,1), orb_pop(i,2), orb_pop(i,3), orb_pop(i,4), orb_pop(i,5)  	!TEST
!enddo

multiplicity=0
do i=1,nstate			! the following loop sets up the variables ndummy and multiplicity (multiplicities and maximum multiplicity)
	ndummy(i)=1
	do j=1,5
		if (orb_pop(i,j)==1) then
			ndummy(i)=ndummy(i)+1
		endif
	enddo	
	if (multiplicity < ndummy(i)) then
		multiplicity=ndummy(i)
	endif
enddo


if (multiplicity == 1) then
	print*, "Careful! All the configurations are closed-shell. The system is not paramagnetic!"
endif

nconfig=0
if (multiplicity==1.or.multiplicity==2) then	! The following test sets up the maximum number of spinorbital configurations associated with a Ms 
	nconfig=1
elseif (multiplicity==3) then
	nconfig=2
elseif (multiplicity==4) then
	nconfig=3
elseif (multiplicity==5) then
	nconfig=6
elseif (multiplicity==6) then
	nconfig=10
endif

allocate (vector_matrix(1:nstate,1:multiplicity,1:nconfig,1:10))

do i=1,nstate
	do j=1,multiplicity
		do k=1,nconfig
			do l=1,10
				vector_matrix(i,j,k,l)=0
			enddo
		enddo
	enddo
enddo


!!!!! The following procedure is fairly complicated but ends up in the construction of the vector_matrix. For each configuration, the procedure builds all possible 
!!!!! spinorbital configurations involved and stores them in a 2-dimensional array config_sample(i,k). Then, the procedure writes these configurations stored in V in appropriate places in vector_matrix,
!!!!! in function of the associated Ms of these configurations. Then repeat for the next electronic configuration


!!!!! Let us build here a matrix with entries corresponding to electronic configurations x entries corresponding to the orbital populations.
!!!! First, we define a set of integer ai(i), af(i) with i=1,5. For doubly populated orbitals ai,bi...=1 and af,bf..=1. For singly-populated orbitals ai,bi..=0, af,bf..=1
!!!! Finally, for empty orbitals, ai,bi...=af,bf...=2. Those define (in non-straightforward manner) the different values that the associated spinorbital populations may take in config_sample.


do j=1,nstate
	do i=1,5
		if (orb_pop(j,i)==2) then
			ai(i)=1
			af(i)=1
		elseif (orb_pop(j,i)==1) then
			ai(i)=0
			af(i)=1
		elseif (orb_pop(j,i)==0) then
			ai(i)=2
			af(i)=2
		endif
	enddo

!!!!!!! After this is done one may define a 2-dimensional array config_sample(i,k')  corresponding to the spinorbital configurations (total dimension 2**(ndummy-1)) and the spinorbitals pop. (1:10)
	allocate (config_sample(1:2**(ndummy(j)-1),1:10))
!!!!!! Now, we can fill the entries of this vector. For each i, we fill the k' entries in multiple manners depending on ai and af. 
!!!!!! That should give an overall combination of all spinorbital configurations possible. 

	i=0
	do while (i<2**(ndummy(j)-1))
		do aa=ai(1),af(1)
			do b=ai(2),af(2)
				do c=ai(3),af(3)
					do d=ai(4),af(4)
						do e=ai(5),af(5)
							a(1)=aa
							a(2)=b
							a(3)=c
							a(4)=d
							a(5)=e
							i=i+1
							do k=1,5
								if (ai(k)==1) then
									config_sample(i,2*k-1)=a(k)
									config_sample(i,2*k)=a(k)
								elseif (ai(k)==0) then
									config_sample(i,2*k-1)=1-a(k)
               	                                	         	config_sample(i,2*k)=a(k)
								elseif (ai(k)==2) then
									config_sample(i,2*k-1)=2-a(k)
									config_sample(i,2*k)=2-a(k)
								endif
							enddo 
						enddo
					enddo
				enddo
			enddo
		enddo
	enddo


!!!!! The vectors are now created. They must now be discriminated by Ms and integrated to the final matrix vector_matrix(j,X,Y,k) where X and Y are indexing the Ms and the configuration, respectively.
!!!!! The jist of it is calculating for each vector the number of unpaired electrons, then associating this number to the variable X (~Ms) 
!!!!! Then, the vector is written in the position corresponding to Y(X)=1(Y is the spinorbital configuration index)  in vector_matrix(j,X,Y(X),k), and Y is incremented by 1, 
!!!!! thus ensuring each configuration is written in a different slot. Y is a 1-dimensional array of rank X (max. 6 for S=5/2)   
!!!!! First, Y is initiated
	do k=1,6
		Y(k)=1
	enddo
!!!!! Now the number 2Ms (variable name Ms) is extracted from the configuration in the vector config_sample(i,k)
	do k=1,i
		Ms=0
		do l=1,5
			Ms=Ms+config_sample(k,2*l-1)-config_sample(k,2*l)
		enddo 
!!!!!! The following formula links the Ms of the considered vector (through the variable Ms) to the rank of X of the configuration matrix vector_matrix
		X=(ndummy(j)-Ms+2)/2
		do l=1,10
			vector_matrix(j,X,Y(X),l)=config_sample(k,l)
		!	print*, vector_matrix(j,X,Y(X),l),j,X,Y(X),k,l
		enddo
!!!!!!!!!!!!!!Now that the vector is written, we must increment the rank of Y(X) in order to write the next configuration of same Ms in a different array
		Y(X)=Y(X)+1
	enddo
	deallocate (config_sample)
enddo

!!!!!! Ooof... We went through it. The worst is behind us..
end subroutine




subroutine matrix_builder(vector_matrix,multiplicity,ndummy,nstate,nconfig,maxdim,Lx2,Ly2,Lz2,Sx2,Sy2,Sz2,SOC2,SDx2,SDy2,SDz2)
! CALCULATES THE SPIN ORBIT COUPLING, X,Y,Z COMPONENT OF THE SPIN, ORBITAL MOMENTUM AND SPIN-DIPOLE (SEE README FILE) MATRICES USED IN LATER SUBROUTINES FROM THE 
! SPINORBITAL CONFIGURATIONS STORED IN VECTOR_MATRIX.
!	VECTOR_MATRIX (INPUT, INTEGER, 4-DIMENSION [1:NSTATE,1;MULTIPLICITY,1:NCONFIG,1:10]) : THE SPIN-EIGENSTATES. IT TAKES THE FORM OF SEVERAL SET OF 10 INTEGER OF VALUE 0 OR 1, EACH CORRESPONDING
! TO A SPINORBITAL POPULATION. EACH SET OF 10 SPINORBITAL POPULATIONS AMOUNTS TO A SPINORBITAL CONFIGURATION.  THERE ARE SEVERAL SUCH SPINORBITAL POPULATIONS FOR A GIVEN MS,
! AND SEVERAL MS FOR A GIVEN MULTIPLICITY/ELECTRONIC CONFIGURATION, ENDING UP IN A 4-DIMENSIONAL TABLE.
! THE FIRST DIMENSION CORRESPONDS TO THE CHOSEN EL. CONFIGURATIONS, THE SECOND CORRESPONDS TO THE MS IN A GIVEN CONFIGURATION, THE THIRD TO THE SPINORBITAL CONFIGURATIONS IN A GIVEN MS, AND THE FOURTH
! THE SPINORBITAL POPULATIONS IN EACH SPINORBITAL CONFIGURATION. CALCULATED FROM THE VECTOR_BUILDE SUBROUTINE.
!	MULTIPLICITY (INPUT,INTEGER): MAX. NUMBER OF MS IN A GIVEN CONFIGURATION. CALCULATED FROM THE VECTOR_BUILDER SUBROUTINE
!	NDUMMY (INPUT, INTEGER, 1-DIMENSIONAL [1:NSTATE]) : INDIVIDUAL MULTIPLICITY FOR EACH CONFIGURATION. CALCULATED IN THE VECTOR_BUILDER SUBROUTINE
!	NSTATE (INPUT, INTEGER) : NUMBER OF ELECTRONIC CONFIGURATIONS
!	NCONFIG (INPUT, INTEGER): MAX. NUMBER OF SPINORBITAL CONFIGURATIONS IN ONE MS. CALCULATED IN THE VECTOR_BUILDER SUBROUTINE
!	MAXDIM (OUTPUT, INTEGER): DIMENSION OF THE MATRICES. IS USED IN LATER SUBROUTINES
!	LX2,LY2,LZ2,SX2,SY2,SZ2,SOC2,SDX2,SDY2,SDZ2 (OUTPUT, COMPLEX, 2-DIMENSIONAL [1:MAXDIM,1:MAXDIM]): THE X,Y,Z COMPONENTS OF THE L MOMENTUM, X,Y,Z COMPONENT OF THE S MOMENTUM, 
!SPIN-ORBIT COUPLING AND X,Y,Z COMPONENT OF THE SPIN-DIPOLE VECTOR. CALCULATED IN THE BASIS OF THE {S2,MS} EIGENSTATES CORRESPONDING TO THE USER-DEFINED ELECTRONIC CONFIGURATION.


implicit none
integer, dimension (1:nstate) :: ndummy
integer, allocatable, dimension(:,:,:,:) :: vector_matrix
integer :: i,j,k,ll,multiplicity,nstate,m,n,index1,index2,flag0,flag00,p,q,count1,count2,flag2,&
&flag1,maxdim,o,spin1,spin2,nconfig,position1,position2,stopcount1,stopcount2
complex*16 :: soc,Lx1,Ly1,Lz1,Sx1,Sy1,Sz1,SDx,SDy,SDz
complex*16, allocatable, dimension (:,:) :: SOC2,Lx2,Ly2,Lz2,Sx2,Sy2,Sz2,SDx2,SDy2,SDz2
character(len=1) :: xyz
real*8 :: flipfactor

p=0
q=0
maxdim=0
do i=1,nstate				
	maxdim=maxdim+ndummy(i)			!! Builds the total dimension of the matrix, i.e.sum of all multiplicities 
enddo
allocate (SOC2(1:maxdim,1:maxdim))
allocate (Lx2(1:maxdim,1:maxdim))
allocate (Ly2(1:maxdim,1:maxdim))
allocate (Lz2(1:maxdim,1:maxdim))
allocate (Sx2(1:maxdim,1:maxdim))
allocate (Sy2(1:maxdim,1:maxdim))
allocate (Sz2(1:maxdim,1:maxdim))
allocate (SDx2(1:maxdim,1:maxdim),SDy2(1:maxdim,1:maxdim),SDz2(1:maxdim,1:maxdim))
p=0

!!!!!!!!!!!!!!! The following procedure (1) develops each interacting pair of {S2,MS} eigenstate into a sum of interacting spinorbital configuration, 
!!!!!!!!!!!!!!! (2) calculates the coupling of each pair of spinorbital configuration via reducing associated determinant coupling into a sum of pair of orbital couplings (using Slater-Condon rules), 
!!!!!!!!!!!!!!  (3) calculates the orbital couplings, (4) sums all configuration-wise couplings to get the total coupling  (matrix elements) of the interacting {S2,MS} eigenstates 
do i=1,nstate
	do j=1,multiplicity
		p=p+1	! p and q, the index of the final matrices, correspond to the interacting spin eigenstates (space x spin function)
		q=0
		do k=1,nstate
			do ll=1,multiplicity
				m=1
				n=1
				count1=0
				count2=0
				soc=0.0
				q=q+1
				if (p<=q) then   ! All considered operators are hermitian so this reduces the operation to the upper triangle
				!	print*, 'p=', p, 'q=',q		!TEST
					soc=0				! All these correspond to the coupling of UNNORMALIZED matrix elements 
					SDx=0.0
					SDy=0.0
					SDz=0.0
					Lx1=0
					Ly1=0
					Lz1=0
					Sx1=0
					Sy1=0
					Sz1=0
					do  m=1,nconfig			!  
						do n=1,nconfig		! Developps the coupling spin-eigenstates into sum of coupling configurations
							flag0=0		! This variable signals if the bra is an empty configuration
							flag00=0	! signals if the ket is an empty configuration
							flag1=0		! tracks the number of different orbitals 
							index1=0   ! marks the position of the orbitals involved in excitation from the bra to the ket
							index2=0    ! idem
							position1=0    ! counts the position of the aformentionned orbitals AMONG OCCUPIED ORBITALS for the bra
							stopcount1=0	! signals that the counting for position1 must be stopped (i.e. the orbital of involved in exc. has been reached)
							position2=0	! same for the ket
							stopcount2=0	! same for the ket
							do o=1,10	! This loops analyzes the differences of pop. in the two configurations
								if (vector_matrix(i,j,m,o)==1) then
									flag0=1  !!! Signals that the configuration is non-zero. Useful for not adding empty configurations(all 0) into the SOC
										 !!!  matrix element
									if (stopcount1==0) then
										position1=position1+1   ! tracks the position of the excited spinorb. in terms of occ. spinorb in the bra
									endif
								endif
								if (vector_matrix(k,ll,n,o)==1) then
       		                       		                          flag00=1  !!! Signals that the configuration is non-zero. Useful for not adding empty configurations into the SOC matrix element
									if (stopcount2==0) then
										position2=position2+1	! tracks the position of the excited orb. in terms of occ. orb in the ket
									endif
 	      	                               		          endif
							!	print*, vector_matrix(i,j,m,o), vector_matrix(k,ll,n,o)
								if (vector_matrix(i,j,m,o)/=vector_matrix(k,ll,n,o)) then
									flag1=flag1+1			! Counts the number of excitations between bra and ket
									if (vector_matrix(i,j,m,o)==1) then
										index1=o		! index 1 is written with the position (among the 10 spinorbitals) of the spinorbital
													! involved in the excitation (in the bra)
										stopcount1=1		! position1 will stop incrementing, and its final value is the position (among occ. orb) of the excited orb.
									elseif (vector_matrix(k,ll,n,o)==1) then
										index2=o		! Same but for the ket
										stopcount2=1		
									endif
								endif    
							enddo  

							if (flag0==0.or.flag00==0) then		! If one of the two configuration is empty, then it must not be counted in the sum
								flag0=0
								flag00=0
							else !!!! Cases where both configurations are non-empty
							!!!!! The next step is to discriminate cases in fonction of the number of different spinorbitals
								if (flag1==2) then  !! In this case there is one excitation
									if (modulo(index1,2)==0) then   !!! In the following tests we separate the spin from the space part in all spinorbs. 
										spin1=1
										index1=index1/2
									elseif (modulo(index1,2)==1) then
										spin1=0
										index1=index1+1
										index1=index1/2
									endif   !!!! extract the spin and space part of the spinorbital
                               	       		 		        if (modulo(index2,2)==0) then
                               	 	       		           	        spin2=1
                               	        	       		    	        index2=index2/2
                                        	       			 elseif (modulo(index2,2)==1) then
                                        	               		 	spin2=0
                                        	                		index2=index2+1
                                               	        			index2=index2/2
                                               				 endif  
									 flipfactor=(-1.0)**(real(position2-position1))			! Factor relative to the exchange of orb. positions in the Slater determinant in order to get two determinant differing by one orbital at the same position (necessary to apply Slater-condon rules)
									soc=soc+(l(index1,index2,'x')*s(spin1,spin2,'x')+l(index1,index2,'y')*s(spin1,spin2,'y')&
&+l(index1,index2,'z')*s(spin1,spin2,'z'))*flipfactor									! Slater-Condon reduction of 2 determinants interacting via one-electron operator
									SDx=SDx+(vv(index1,index2,'x','x')*s(spin1,spin2,'x')+vv(index1,index2,'x','y')*s(spin1,spin2,'y')&
&+vv(index1,index2,'x','z')*s(spin1,spin2,'z'))*flipfactor
									SDy=SDy+(vv(index1,index2,'y','x')*s(spin1,spin2,'x')+vv(index1,index2,'y','y')*s(spin1,spin2,'y')&
&+vv(index1,index2,'y','z')*s(spin1,spin2,'z'))*flipfactor
									SDz=SDz+(vv(index1,index2,'z','x')*s(spin1,spin2,'x')+vv(index1,index2,'z','y')*s(spin1,spin2,'y')&
&+vv(index1,index2,'z','z')*s(spin1,spin2,'z'))*flipfactor
									if (spin1==spin2) then				! interaction via L is non-zero only if same ms
										Lx1=Lx1+l(index1,index2,'x')*flipfactor
										Ly1=Ly1+l(index1,index2,'y')*flipfactor
										Lz1=Lz1+l(index1,index2,'z')*flipfactor
									endif
									if (index1==index2) then			! interaction with S non 0 only if same ml
										Sx1=Sx1+s(spin1,spin2,'x')*flipfactor
										Sy1=Sy1+s(spin1,spin2,'y')*flipfactor
										Sz1=Sz1+s(spin1,spin2,'z')*flipfactor
									endif
									
								elseif (flag1==0) then				! In this case there are no differences between bra and ket
								!	print*, 'diagonal elements!'		! TEST
									do o=1,10				!sum over populated orbitals (Slater-Condon rule) 
										if (vector_matrix(i,j,m,o)==1) then ! only populated spin-orbs participate to the sum
											if (modulo(o,2)==0) then   !!! In the following we separate the spin from the space part in all spinorbs.
                                                                        	       		spin1=1
                                                                                		index1=o/2
                                                                        		elseif (modulo(o,2)==1) then
                                                                                		spin1=0
                                                                                		index1=o+1
                                                                                		index1=index1/2
                                                                        		endif
											soc=soc+l(index1,index1,'x')*s(spin1,spin1,'x')+l(index1,index1,'y')*s(spin1,spin1,'y')+&
&l(index1,index1,'z')*s(spin1,spin1,'z')						!Slater condon rule for two identical determinant interacting via a one-electron operator						
											SDx=SDx+vv(index1,index1,'x','x')*s(spin1,spin1,'x')+vv(index1,index1,'x','y')*s(spin1,spin1,'y')&
&+vv(index1,index1,'x','z')*s(spin1,spin1,'z')
											SDy=SDy+vv(index1,index1,'y','x')*s(spin1,spin1,'x')+vv(index1,index1,'y','y')*s(spin1,spin1,'y')&
&+vv(index1,index1,'y','z')*s(spin1,spin1,'z')
											SDz=SDz+vv(index1,index1,'z','x')*s(spin1,spin1,'x')+vv(index1,index1,'z','y')*s(spin1,spin1,'y')&
&+vv(index1,index1,'z','z')*s(spin1,spin1,'z')


											Lx1=Lx1+l(index1,index1,'x')
											Ly1=Ly1+l(index1,index1,'y')
											Lz1=Lz1+l(index1,index1,'z')
											Sx1=Sx1+s(spin1,spin1,'x')
											Sy1=Sy1+s(spin1,spin1,'y')
											Sz1=Sz1+s(spin1,spin1,'z')
										endif
										
									enddo
								else
								!	print*, 'more than two orbitals differing'	!TEST
								endif	
							endif
						enddo
					enddo
				SOC2(p,q)=soc/(NORM(j,ndummy(i))*NORM(ll,ndummy(i))) !! All these correspond to the coupling of NORMALIZED spin-eigenstates
				Lx2(p,q)=Lx1/(NORM(j,ndummy(i))*NORM(ll,ndummy(i)))
				Ly2(p,q)=Ly1/(NORM(j,ndummy(i))*NORM(ll,ndummy(i)))
				Lz2(p,q)=Lz1/(NORM(j,ndummy(i))*NORM(ll,ndummy(i)))
				Sx2(p,q)=Sx1/(NORM(j,ndummy(i))*NORM(ll,ndummy(i)))
				Sy2(p,q)=Sy1/(NORM(j,ndummy(i))*NORM(ll,ndummy(i)))
				Sz2(p,q)=Sz1/(NORM(j,ndummy(i))*NORM(ll,ndummy(i)))
				SDx2(p,q)=SDx/(NORM(j,ndummy(i))*NORM(ll,ndummy(i)))
				SDy2(p,q)=SDy/(NORM(j,ndummy(i))*NORM(ll,ndummy(i)))
				SDz2(p,q)=SDz/(NORM(j,ndummy(i))*NORM(ll,ndummy(i)))

				endif
			enddo
		enddo
	enddo
enddo


do i=1,maxdim			!! Setting the lower triangle of the matrix
	do j=1,maxdim
		if (i<j) then
			SOC2(j,i)=conjg(SOC2(i,j))
			Lx2(j,i)=conjg(Lx2(i,j))
			Ly2(j,i)=conjg(Ly2(i,j))
			Lz2(j,i)=conjg(Lz2(i,j))
			Sx2(j,i)=conjg(Sx2(i,j))
			Sy2(j,i)=conjg(Sy2(i,j))
			Sz2(j,i)=conjg(Sz2(i,j))
			SDx2(j,i)=conjg(SDx2(i,j))
			SDy2(j,i)=conjg(SDy2(i,j))
			SDz2(j,i)=conjg(SDz2(i,j))
		endif
	enddo
enddo



end subroutine


subroutine susc_calc(Ntemp,T,gridlevel,job_name,DE,SOC2,Lx2,Ly2,Lz2,Sx2,Sy2,Sz2,maxdim,ndummy,lambda,B)		

!CALCULATES THE MAGNETIZATION/SUSCEPTIBILITY/CHI*T FOR EACH SPECIFIED TEMPERATURE POINT, USING THE MATRICES DETERINED IN MATRIX_BUILDER SUBROUTINE.
! IN DETAILS FOR EACH DIRECTION OF THE MAGNETIC FIELD, THE SET OF HAMILTONIAN EIGENVALUES IS STORED IN A MATRIX.
! THEN FOR EACH TEMPERATURE, THE SUBROUTINE CALCULATES THE DERIVATE OF THE PARTITION FUNCTION AND INFERS FROM IT THE DIRECTION-DEPENDENT MAGNETIZATION ALONG THE APPLIED FIELD DIRECTION.
! THE AVERAGE MAGNETIZATION IS OBTAINED BY NUMERICAL INTEGRATION OF THE DIRECTION-DEPENDENT MAGNETIZATION.
!	NTEMP (INPUT, INTEGER): NUMBER OF TEMPERATURE POINTS FOR WHICH TO CALCULATE THE MAGNETIZATION
!	T (INPUT, REAL, 1-DIMENSIONAL [1:NTEMP]): TEMPERATURE POINTS.
!	GRIDELEVEL (INPUT, CHARACTER): DEFINES THE TYPE OF GRID FOR THE INTEGRATION. 'LEB' USES LEBEDEV QUADRATURE WITH 38 DIRECTION-DEPENDENT MAGNETIZATIONS CALCULATED FOR EACH TEMPERATURE.
!'NORMAL' INTEGRATES NUMERICALLY WITH A POLAR AND AZIMUTHAL ANGLE INCREMENT OF 0.01*PI, WHICH ENDS UP IN A CALCULATION OF 20 000 DIRECTION-DEPENDENT MAGNETIZATIONS. 'FINE' INTEGRATES WITH AN ANGLE 
! INCREMENT OF 0.005*PI, WHICH AMOUNTS TO 80 000 DIRECTION-DEPENDENT MAGNETIZATIONS. IT SEEMS  THE LEBEDEV QUADRATURE YIELDS A REASONABLE ACCURACY. 
!	JOB_NAME (INPUT, CHARACTER): SPECIFIES WHETHER TO PRINT OUT THE MAGNETIC SUSCEPTIBILITY ('X'), CHI*T ('XT'), OR THE MAGNETIZATION ('MAG')
!	DE (INPUT, REAL, 1-DIMENSION [1:NSTATE]): THE VERTICAL EXCITATION ENERGIES ASSOCIATED WITH EACH ELECTRONIC CONFIGURATION, RELATIVE TO THE FIRST DEFINED ELECTRONIC CONFIGURATION.  
!	SOC2,LX2,LY2,LZ2,SX2,SY2,SZ2 (INPUT, COMPLEX, 2-DIMENSIONAL [1:MAXDIM,1:MAXDIM]): SPIN-ORBIT COUPLING, X,Y,Z COMPONENTS OF THE ORBITAL AND SPIN MOMENTA MATRICES. 
! TOGETHER WITH DE, THEY ARE USED TO CONSTRUCT THE TOTAL HAMILTONIAN.
! 	MAXDIM (INPUT, INTEGER) : THE TOTAL DIMENSION OF THE HAMILTONIAN
!	NDUMMY (INPUT, INTEGER, 1-DIMENSION [1:NSTATE]): MULTIPLICITY OF EACH ELECTRONIC CONFIGURATION
!	LAMBDA (INPUT, REAL) : SPIN-ORBIT COUPLING CONSTANT OF THE TRANSITION METAL
!	B (INPUT,REAL): NORM OF THE MAGNETIC FIELD

implicit none													
integer :: Ntemp,i,maxdim,u,v,j,k,ll,m
real*16, allocatable, dimension (:) :: T
real*16 :: kB=0.695039,mag,dmag,ex,ey,ez,lambda,B,pi=3.1415927,increment
real*16, dimension (0:37) :: Coeff
character (len=3) :: job_name
character (len=100) :: gridlevel
real*16, dimension (-1:1) :: Z
real*16, dimension (0:37,1:3) :: emat
real*16, allocatable, dimension (:) :: DE
integer, allocatable, dimension (:) :: ndummy
real*8, allocatable, dimension (:) :: eigenval
complex*16, allocatable, dimension(:,:) :: SOC2,Lx2,Ly2,Lz2,Sx2,Sy2,Sz2,result_micro, HAM
complex*16, allocatable, dimension (:,:,:,:) :: magnetization_matrix,store_matrix

gridlevel=trim(gridlevel)

if (gridlevel=='normal') then	! sets up the number and positions of the nodes in the quadrature for the powder-averaging
	increment=0.01*pi
	u=100
	v=200
elseif (gridlevel=='fine') then
	increment=0.005*pi
	u=200
	v=400
elseif (gridlevel=='leb') then
	u=37
	v=0
	call lebedev(Coeff,emat)	! Sets up the directional vectors and coefficients associated with each node
endif
allocate (store_matrix(0:u,0:v,1:maxdim,-1:1))
do i=0,u
	do j=0,v
		ez=cos(real(i)*increment)                               
                ey=sin(real(i)*increment)*sin(real(j)*increment)
                ex=sin(real(i)*increment)*cos(real(j)*increment)
		if (gridlevel=='leb') then
			ex=emat(i,1)
			ey=emat(i,2)
			ez=emat(i,3)
		endif
		call diag_Ham(DE,SOC2,Lx2,Ly2,Lz2,Sx2,Sy2,Sz2,ex,ey,ez,maxdim,ndummy,lambda,B,HAM,eigenval)	! diagonalizes the electronic Hamiltonian and returns the eigenvalues
		do k=1,maxdim
			store_matrix(i,j,k,0)=eigenval(k)							! stores the eigenvalues corresponding to the applied field direction
		enddo
		deallocate (HAM,eigenval)
		B=B-0.05
		call diag_Ham(DE,SOC2,Lx2,Ly2,Lz2,Sx2,Sy2,Sz2,ex,ey,ez,maxdim,ndummy,lambda,B,HAM,eigenval)	! repeat for B-dB (for later derivation of log(Z))
                do k=1,maxdim
                        store_matrix(i,j,k,-1)=eigenval(k)
                enddo
		B=B+0.1
		deallocate (HAM,eigenval)
		call diag_Ham(DE,SOC2,Lx2,Ly2,Lz2,Sx2,Sy2,Sz2,ex,ey,ez,maxdim,ndummy,lambda,B,HAM,eigenval)	! repeat for B+dB (for later derivation of log(Z)
                do k=1,maxdim
                        store_matrix(i,j,k,1)=eigenval(k)
                enddo
		B=B-0.05
		deallocate (HAM,eigenval)
	enddo
enddo

!! Note: we calculate three set of eigenfunction: 1 for B, one for B+dB, one for B-dB. We leave it like that because we may want to introduce later an alternative definition of the magnetic susceptibility

!! as X=d2E/dB2 instead of M/B. For now, the calculation for B is useless, but especially with Lebedev quadrature, the time wasted is negligible.

!!!!!!!! At this point all direction-dependent energies have been stored. Now for each temperature the direction-dependent magnetization is calculated 
!!!!!!!! from the partition function (calculated from the eigenvalues) and averaged over the sphere. Note. This step is decoupled from the last one to avoid 
!!!!!!!! rediagonalizing the Hamiltonian for each temperature, which is completely unnecessary and excruciatingly slow. 

do i=1,NTemp
	mag=0
	do j=0,u
		do k=0,v
			Z(-1)=0.0
			Z(0)=0.0
			Z(1)=0.0
			do ll=1,maxdim			! setting up the partition functions
				do m=-1,1
					Z(m)=Z(m)+exp(-(store_matrix(j,k,ll,m)-store_matrix(j,k,1,1))/(kB*T(i)))
				enddo
			enddo
			if (gridlevel=='leb') then	! two different formula depending on the type of quadrature
				dmag=T(i)*1.38*6.02*(10.0**7.0)*((log(Z(1))-log(Z(-1)))/1000.0)		!dmag is the direction-dependent magnetization
				mag=mag+(dmag*Coeff(j))		!mag is the final, powder-averaged function after incrementation of all direction-dependent function
			else
				dmag=T(i)*1.38*6.02*(10.0**7.0)*((log(Z(1))-log(Z(-1)))/1000.0)*(increment**2.0)*sin(real(j)*increment)		!dmag is the direction-dependent magnetization
				mag=mag+dmag
			endif
		enddo
	enddo
	if (gridlevel=='leb') then
	!	mag=mag
	else
		mag=mag/(4.0*pi)		! a factor 1/4pi must be added to get the averaged function over the sphere direction. this factor is already included in the Lebedev quadrature
	endif
	if (job_name=='mag') then		! different functions of mag must be printed, depending on whether the user specified magnetization, chi or chi*T
		print*, mag
	elseif (job_name=='X') then
		mag=mag/(B*10000.0)
		print*, mag
	elseif (job_name=='XT') then
		mag=mag*T(i)/(B*10000.0)
		print*, mag
	endif
enddo

end subroutine

subroutine lebedev(Coeff,emat)
! RETRUNS THE DIRECTIONAL VECTORS AND COEFFCIENTS ASSOCIATED WITH THE 38 NODES IN LEBEDEV QUADRATURE
!	COEFF (OUTPUT, REAL, 1-DIMENSIONAL [0:37]): RETURNS THE COEFFICIENT ASSOCIATED WITH THE NODE
!	EMAT (OUTPUT, REAL, 2-DIMENSIONAL [0:37,1:3]): EMAT(I,J) RETURNS THE JTH COMPONENT (X,Y, OR Z IN THAT ORDER) OF THE DIRECTIONAL VECTOR ASSOCIATED WITH THE NODE I
implicit none
real*16, dimension (0:37) :: Coeff
real*16, dimension (0:37,1:3) :: emat    ! x=1, y=2, z=3
real*16 :: p,q
integer :: i
p=0.8880738
q=0.4597008
emat(0,1)=0.0
emat(0,2)=0.0
emat(0,3)=1.0
emat(1,1)=0.0
emat(1,2)=0.0
emat(1,3)=-1.0
emat(2,1)=0.0
emat(2,2)=1.0
emat(2,3)=0.0
emat(3,1)=0.0
emat(3,2)=-1.0
emat(3,3)=0.0
emat(4,1)=1.0
emat(4,2)=0.0
emat(4,3)=0.0
emat(5,1)=-1.0
emat(5,2)=0.0
emat(5,3)=0.0
emat(6,1)=1.0/(3.0**0.5)
emat(6,2)=1.0/(3.0**0.5)
emat(6,3)=1.0/(3.0**0.5)
emat(7,1)=-1.0/(3.0**0.5)
emat(7,2)=1.0/(3.0**0.5)
emat(7,3)=1.0/(3.0**0.5)
emat(8,1)=1.0/(3.0**0.5)
emat(8,2)=-1.0/(3.0**0.5)
emat(8,3)=1.0/(3.0**0.5)
emat(9,1)=1.0/(3.0**0.5)
emat(9,2)=1.0/(3.0**0.5)
emat(9,3)=-1.0/(3.0**0.5)
emat(10,1)=-1.0/(3.0**0.5)
emat(10,2)=-1.0/(3.0**0.5)
emat(10,3)=1.0/(3.0**0.5)
emat(11,1)=-1.0/(3.0**0.5)
emat(11,2)=1.0/(3.0**0.5)
emat(11,3)=-1.0/(3.0**0.5)
emat(12,1)=1.0/(3.0**0.5)
emat(12,2)=-1.0/(3.0**0.5)
emat(12,3)=-1.0/(3.0**0.5)
emat(13,1)=-1.0/(3.0**0.5)
emat(13,2)=-1.0/(3.0**0.5)
emat(13,3)=-1.0/(3.0**0.5)
emat(14,1)=p
emat(14,2)=q
emat(14,3)=0.0
emat(15,1)=-p
emat(15,2)=q
emat(15,3)=0.0
emat(16,1)=p
emat(16,2)=-q
emat(16,3)=0.0
emat(17,1)=-p
emat(17,2)=-q
emat(17,3)=0.0
emat(18,1)=p
emat(18,2)=0.0
emat(18,3)=q
emat(19,1)=-p
emat(19,2)=0.0
emat(19,3)=q
emat(20,1)=p
emat(20,2)=0.0
emat(20,3)=-q
emat(21,1)=-p
emat(21,2)=0.0
emat(21,3)=-q
emat(22,1)=0.0
emat(22,2)=p
emat(22,3)=q
emat(23,1)=0.0
emat(23,2)=-p
emat(23,3)=q
emat(24,1)=0.0
emat(24,2)=p
emat(24,3)=-q
emat(25,1)=0.0
emat(25,2)=-p
emat(25,3)=-q
emat(26,1)=q
emat(26,2)=p
emat(26,3)=0.0
emat(27,1)=-q
emat(27,2)=p
emat(27,3)=0.0
emat(28,1)=q
emat(28,2)=-p
emat(28,3)=0.0
emat(29,1)=-q
emat(29,2)=-p
emat(29,3)=0.0
emat(30,1)=q
emat(30,2)=0.0
emat(30,3)=p
emat(31,1)=-q
emat(31,2)=0.0
emat(31,3)=p
emat(32,1)=q
emat(32,2)=0.0
emat(32,3)=-p
emat(33,1)=-q
emat(33,2)=0.0
emat(33,3)=-p
emat(34,1)=0.0
emat(34,2)=q
emat(34,3)=p
emat(35,1)=0.0
emat(35,2)=-q
emat(35,3)=p
emat(36,1)=0.0
emat(36,2)=q
emat(36,3)=-p
emat(37,1)=0.0
emat(37,2)=-q
emat(37,3)=-p

do i=0,37
	if (i<=5) then
		Coeff(i)=0.00952
	elseif (i>5.and.i<=13) then
		Coeff(i)=0.03214
	elseif (i>13.and.i<=37) then
		Coeff(i)=0.02857
	endif
enddo


end subroutine

subroutine mosleb(Coeff,emat)
! RETRUNS THE DIRECTIONAL VECTORS AND COEFFCIENTS ASSOCIATED WITH THE 86 NODES IN LEBEDEV QUADRATURE
!       COEFF (OUTPUT, REAL, 1-DIMENSIONAL [0:37]): RETURNS THE COEFFICIENT ASSOCIATED WITH THE NODE
!       EMAT (OUTPUT, REAL, 2-DIMENSIONAL [0:37,1:3]): EMAT(I,J) RETURNS THE JTH COMPONENT (X,Y, OR Z IN THAT ORDER) OF THE DIRECTIONAL VECTOR ASSOCIATED WITH THE NODE I
implicit none
real*16, dimension (0:85) :: Coeff
real*16, dimension (0:85,1:3) :: emat    ! x=1, y=2, z=3
real*16 :: p,q,m1,l1,m2,l2
integer :: i

p=0.9169
q=0.3992
m1=0.87852
m2=0.3643
l1=0.3378
l2=0.6585

emat(0,1)=0.0
emat(0,2)=0.0
emat(0,3)=1.0
emat(1,1)=0.0
emat(1,2)=0.0
emat(1,3)=-1.0
emat(2,1)=0.0
emat(2,2)=1.0
emat(2,3)=0.0
emat(3,1)=0.0
emat(3,2)=-1.0
emat(3,3)=0.0
emat(4,1)=1.0
emat(4,2)=0.0
emat(4,3)=0.0
emat(5,1)=-1.0
emat(5,2)=0.0
emat(5,3)=0.0
emat(6,1)=1.0/(3.0**0.5)
emat(6,2)=1.0/(3.0**0.5)
emat(6,3)=1.0/(3.0**0.5)
emat(7,1)=-1.0/(3.0**0.5)
emat(7,2)=1.0/(3.0**0.5)
emat(7,3)=1.0/(3.0**0.5)
emat(8,1)=1.0/(3.0**0.5)
emat(8,2)=-1.0/(3.0**0.5)
emat(8,3)=1.0/(3.0**0.5)
emat(9,1)=1.0/(3.0**0.5)
emat(9,2)=1.0/(3.0**0.5)
emat(9,3)=-1.0/(3.0**0.5)
emat(10,1)=-1.0/(3.0**0.5)
emat(10,2)=-1.0/(3.0**0.5)
emat(10,3)=1.0/(3.0**0.5)
emat(11,1)=-1.0/(3.0**0.5)
emat(11,2)=1.0/(3.0**0.5)
emat(11,3)=-1.0/(3.0**0.5)
emat(12,1)=1.0/(3.0**0.5)
emat(12,2)=-1.0/(3.0**0.5)
emat(12,3)=-1.0/(3.0**0.5)
emat(13,1)=-1.0/(3.0**0.5)
emat(13,2)=-1.0/(3.0**0.5)
emat(13,3)=-1.0/(3.0**0.5)
emat(14,1)=p
emat(14,2)=q
emat(14,3)=0.0
emat(15,1)=-p
emat(15,2)=q
emat(15,3)=0.0
emat(16,1)=p
emat(16,2)=-q
emat(16,3)=0.0
emat(17,1)=-p
emat(17,2)=-q
emat(17,3)=0.0
emat(18,1)=p
emat(18,2)=0.0
emat(18,3)=q
emat(19,1)=-p
emat(19,2)=0.0
emat(19,3)=q
emat(20,1)=p
emat(20,2)=0.0
emat(20,3)=-q
emat(21,1)=-p
emat(21,2)=0.0
emat(21,3)=-q
emat(22,1)=0.0
emat(22,2)=p
emat(22,3)=q
emat(23,1)=0.0
emat(23,2)=-p
emat(23,3)=q
emat(24,1)=0.0
emat(24,2)=p
emat(24,3)=-q
emat(25,1)=0.0
emat(25,2)=-p
emat(25,3)=-q
emat(26,1)=q
emat(26,2)=p
emat(26,3)=0.0
emat(27,1)=-q
emat(27,2)=p
emat(27,3)=0.0
emat(28,1)=q
emat(28,2)=-p
emat(28,3)=0.0
emat(29,1)=-q
emat(29,2)=-p
emat(29,3)=0.0
emat(30,1)=q
emat(30,2)=0.0
emat(30,3)=p
emat(31,1)=-q
emat(31,2)=0.0
emat(31,3)=p
emat(32,1)=q
emat(32,2)=0.0
emat(32,3)=-p
emat(33,1)=-q
emat(33,2)=0.0
emat(33,3)=-p
emat(34,1)=0.0
emat(34,2)=q
emat(34,3)=p
emat(35,1)=0.0
emat(35,2)=-q
emat(35,3)=p
emat(36,1)=0.0
emat(36,2)=q
emat(36,3)=-p
emat(37,1)=0.0
emat(37,2)=-q
emat(37,3)=-p

emat(38,1)=l1
emat(38,2)=l1
emat(38,3)=m1
emat(39,1)=-l1
emat(39,2)=l1
emat(39,3)=m1
emat(40,1)=l1
emat(40,2)=-l1
emat(40,3)=m1
emat(41,1)=l1
emat(41,2)=l1
emat(41,3)=-m1
emat(42,1)=-l1
emat(42,2)=-l1
emat(42,3)=m1
emat(43,1)=-l1
emat(43,2)=l1
emat(43,3)=-m1
emat(44,1)=l1
emat(44,2)=-l1
emat(44,3)=-m1
emat(45,1)=-l1
emat(45,2)=-l1
emat(45,3)=-m1
emat(46,1)=l1
emat(46,2)=m1
emat(46,3)=l1
emat(47,1)=-l1
emat(47,2)=m1
emat(47,3)=l1
emat(48,1)=l1
emat(48,2)=-m1
emat(48,3)=l1
emat(49,1)=l1
emat(49,2)=m1
emat(49,3)=-l1
emat(50,1)=-l1
emat(50,2)=-m1
emat(50,3)=l1
emat(51,1)=-l1
emat(51,2)=m1
emat(51,3)=-l1
emat(52,1)=l1
emat(52,2)=-m1
emat(52,3)=-l1
emat(53,1)=-l1
emat(53,2)=-m1
emat(53,3)=-l1
emat(54,1)=m1
emat(54,2)=l1
emat(54,3)=l1
emat(55,1)=-m1
emat(55,2)=l1
emat(55,3)=l1
emat(56,1)=m1
emat(56,2)=-l1
emat(56,3)=l1
emat(57,1)=m1
emat(57,2)=l1
emat(57,3)=-l1
emat(58,1)=-m1
emat(58,2)=-l1
emat(58,3)=l1
emat(59,1)=-m1
emat(59,2)=l1
emat(59,3)=-l1
emat(60,1)=m1
emat(60,2)=-l1
emat(60,3)=-l1
emat(61,1)=-m1
emat(61,2)=-l1
emat(61,3)=-l1


emat(62,1)=l2
emat(62,2)=l2
emat(62,3)=m2
emat(63,1)=-l2
emat(63,2)=l2
emat(63,3)=m2
emat(64,1)=l2
emat(64,2)=-l2
emat(64,3)=m2
emat(65,1)=l2
emat(65,2)=l2
emat(65,3)=-m2
emat(66,1)=-l2
emat(66,2)=-l2
emat(66,3)=m2
emat(67,1)=-l2
emat(67,2)=l2
emat(67,3)=-m2
emat(68,1)=l2
emat(68,2)=-l2
emat(68,3)=-m2
emat(69,1)=-l2
emat(69,2)=-l2
emat(69,3)=-m2
emat(70,1)=l2
emat(70,2)=m2
emat(70,3)=l2
emat(71,1)=-l2
emat(71,2)=m2
emat(71,3)=l2
emat(72,1)=l2
emat(72,2)=-m2
emat(72,3)=l2
emat(73,1)=l2
emat(73,2)=m2
emat(73,3)=-l2
emat(74,1)=-l2
emat(74,2)=-m2
emat(74,3)=l2
emat(75,1)=-l2
emat(75,2)=m2
emat(75,3)=-l2
emat(76,1)=l2
emat(76,2)=-m2
emat(76,3)=-l2
emat(77,1)=-l2
emat(77,2)=-m2
emat(77,3)=-l2
emat(78,1)=m2
emat(78,2)=l2
emat(78,3)=l2
emat(79,1)=-m2
emat(79,2)=l2
emat(79,3)=l2
emat(80,1)=m2
emat(80,2)=-l2
emat(80,3)=l2
emat(81,1)=m2
emat(81,2)=l2
emat(81,3)=-l2
emat(82,1)=-m2
emat(82,2)=-l2
emat(82,3)=l2
emat(83,1)=-m2
emat(83,2)=l2
emat(83,3)=-l2
emat(84,1)=m2
emat(84,2)=-l2
emat(84,3)=-l2
emat(85,1)=-m2
emat(85,2)=-l2
emat(85,3)=-l2



do i=0,85
        if (i<=5) then
                Coeff(i)=0.011544
        elseif (i>5.and.i<=13) then
                Coeff(i)=0.011944
        elseif (i>13.and.i<=37) then
                Coeff(i)=0.011812
	elseif (i>37.and.i<=61) then
		Coeff(i)=0.011111
	elseif (i>61.and.i<=85) then
		Coeff(i)=0.011812
        endif

enddo

end subroutine

subroutine mag_micro(HAM,maxdim,ndummy,Lx2,Ly2,Lz2,Sx2,Sy2,Sz2,SDx2,SDy2,SDz2,result_micro)   
! CALCULATES THE EXPECTATION VALUES OF THE SPIN, ORBITAL MOMENTA AND SPIN-DIPOLE VECTOR OPERATORS FOR EACH HAMILTONIAN EIGENSTATES. 
! IT CALCULATES THEM FROM THE MATRIX IN THE BASIS OF THE {S2,MS} EIGENSTATES AND THE
! EIGENVECTORS STATE COMPOSITION IN THE BASIS OF THE {S2,MS} EIGENSTATES. IN DETAILS, THE BASIS IN WHICH THE RELEVANT OPERATORS IS SHIFTED FROM THE {S2,MS} EIGENSTATES
! BASIS TO THE EIGENVECTOR BASIS USING O'=C(TRANS)OC WHERE C IS THE EIGENVECTOR BASIS, C(TRANS) THE TRANSPOSED MATRIX OF C AND O IS THE OPERATOR IN THE {S2,MS} EIGENSTATES
! BASIS. THEN, THE DIAGONAL VALUES OF THE REULTNG O' MATRICES (THE OPERATOR EXPECTATION VALUE OF THE EIGENSTATE) IS STORED IN A SMALLER MATRIX RESULT_MICRO
!	HAM (INPUT, COMPLEX, 2-DIMENSIONAL [1:MAXDIM,1:MAXDIM]): MATRIX CONTAINING THE EIGENVECTORS (HAM(I,J) IS THE COEFFICIENT OF THE ITH {S2,MS} EIGENSTATES IN THE JTH HAM. EIGENSTATE
!	MAXDIM (INPUT, INTEGER): DIMENSION OF THE HAMILTONIAN
!	Lx2,Ly2,Lz2,Sx2,Sy2,Sz2,SDx2,SDy2,SDz2 (INPUT, COMPLEX, 2-DIMENSIONAL [1:MAXDIM,1:MAXDIM]): X,Y,Z COMPONENTS OF THE ORBITAL SPIN, MOMENTA AND SPIN-DIPOLE MATRICES.
!	RESULT_MICRO (OUTPUT, COMPLEX, 2-DIMENSIONAL[1:MAXDIM,1:9]): RESULT_MICRO(I,J) CONTAINS THE EXPECTATION VALUES OF LX(J=1),LY(J=2),LZ(J=3),SX(J=4),SY(J=5),
!	SZ(J=6),SDX(J=7),SDY(J=8),SDZ(J=9) FOR THE ITH HAMILTONIAN EIGENVECTOR 

	   
implicit none
real*16, allocatable, dimension (:) :: DE
integer, allocatable, dimension (:) :: ndummy
real*16 :: ex,ey,ez,lambda,B
complex*16, allocatable, dimension (:,:) :: SOC2,Lx2,Ly2,Lz2,Sx2,Sy2,Sz2,TLx2,SDx2,SDy2,SDz2,TLy2,TLz2,TSx2,TSy2,TSz2,HAM,darkHAM&
&,ELx2,ELy2,ELz2,ESx2,ESy2,ESz2,Tsdx,Tsdy,Tsdz,Esdx,Esdy,Esdz,result_micro
real*8, allocatable, dimension (:) :: eigenval
integer :: i,j,maxdim

! Allocates the required matrices
allocate (darkHAM(1:maxdim,1:maxdim))
allocate (TLx2(1:maxdim,1:maxdim),TLy2(1:maxdim,1:maxdim),TLz2(1:maxdim,1:maxdim),TSx2(1:maxdim,1:maxdim),&
&Tsy2(1:maxdim,1:maxdim),Tsz2(1:maxdim,1:maxdim),Tsdx(1:maxdim,1:maxdim),Tsdy(1:maxdim,1:maxdim),Tsdz(1:maxdim,1:maxdim))
allocate (ELx2(1:maxdim,1:maxdim),ELy2(1:maxdim,1:maxdim),ELz2(1:maxdim,1:maxdim),ESx2(1:maxdim,1:maxdim),&
&Esy2(1:maxdim,1:maxdim),Esz2(1:maxdim,1:maxdim),Esdx(1:maxdim,1:maxdim),Esdy(1:maxdim,1:maxdim),Esdz(1:maxdim,1:maxdim))

do i=1,maxdim
	do j=1,maxdim
		darkHAM(i,j)=conjg(HAM(j,i)) ! The transpose of the Hamiltonian is necessary for the shift of basis.
	enddo
enddo

!!!!! The following procedure operates the shift of basis for all the relevant operators, the matrices in the new basis being ELx, ELy,ELz,ESx, etc.. 

call multmat(Lx2,HAM,TLx2,maxdim)	! Straightforward matrix multiplication subroutine
call multmat(darkHAM,TLx2,ELx2,maxdim)

call multmat(Ly2,HAM,TLy2,maxdim)
call multmat(darkHAM,TLy2,ELy2,maxdim)

call multmat(Lz2,HAM,TLz2,maxdim)
call multmat(darkHAM,TLz2,ELz2,maxdim)

call multmat(Sx2,HAM,TSx2,maxdim)
call multmat(darkHAM,TSx2,ESx2,maxdim)

call multmat(Sy2,HAM,TSy2,maxdim)
call multmat(darkHAM,TSy2,ESy2,maxdim)

call multmat(Sz2,HAM,TSz2,maxdim)
call multmat(darkHAM,TSz2,ESz2,maxdim)

call multmat(SDx2,HAM,Tsdx,maxdim)
call multmat(darkHAM,Tsdx,Esdx,maxdim)

call multmat(SDy2,HAM,Tsdy,maxdim)
call multmat(darkHAM,Tsdy,Esdy,maxdim)

call multmat(SDz2,HAM,Tsdz,maxdim)
call multmat(darkHAM,Tsdz,Esdz,maxdim)

!!!!!! Now all the relevant physical constants have been written. The values are tabulated in result_micro(i,j) (dim maxdim x 9) . For j, from 1 to 9, expec. of Lx,Ly,Lz,Sx,Sy,Sz,SDx,SDy,SDz

do i=1,maxdim
	result_micro(i,1)=ELx2(i,i)
	result_micro(i,2)=ELy2(i,i)
	result_micro(i,3)=ELz2(i,i)
	result_micro(i,4)=ESx2(i,i)
	result_micro(i,5)=ESy2(i,i)
	result_micro(i,6)=ESz2(i,i)
	result_micro(i,7)=Esdx(i,i)
	result_micro(i,8)=Esdy(i,i)
	result_micro(i,9)=Esdz(i,i)
enddo
end subroutine

subroutine diag_Ham(DE,SOC2,Lx2,Ly2,Lz2,Sx2,Sy2,Sz2,ex,ey,ez,maxdim,ndummy,lambda,B,HAM,eigenval)

! CALCULATES AND DIAGONALIZES THE TOTAL HAMILTONIAN FROM THE VARIOUS MATRICES DETERMINED IN MATRIX_BUILDER SUBROUTINE AND THE CHOSEN APPLED FIELD INTENSITY AND DIRECTIONAL VECTORS
!	DE (INPUT, REAL, 1-DIMENSION [1:NSTATE]): THE VERTICAL EXCITATION ENERGIES ASSOCIATED WITH EACH ELECTRONIC CONFIGURATION, RELATIVE TO THE FIRST DEFINED ELECTRONIC CONFIGURATION.
!       SOC2,LX2,LY2,LZ2,SX2,SY2,SZ2 (INPUT, COMPLEX, 2-DIMENSIONAL [1:MAXDIM,1:MAXDIM]): SPIN-ORBIT COUPLING, X,Y,Z COMPONENTS OF THE ORBITAL AND SPIN MOMENTA MATRICES.
!	EX,EY,EZ (INPUT,REAL): DIRECTIONAL VECTORS OF THE APPLIED FIELD
! 	MAXDIM (INPUT, INTEGER): DIMENSION OF THE HAMILTONIAN
!	NDUMMY (INPUT, INTEGER, 1-DIMENSIONAL[1:NSTATE]): MULTIPLICITY ASSOCIATED WITH EACH CONFIGURATION
!	LAMBDA (INPUT, REAL): SPIN-ORBIT COUPLING CONSTANT OF THE TRANSITION METAL.
!	HAM (OUTPUT, COMPLEX, 2-DIMENSIONAL [1:MAXDIM,1:MAXDIM]): BEFORE DIAGONALIZATION, TOTAL HAMILTONIAN MATRIX; AFTER DIAGONALIZATION, CONTAINS THE EIGENVECTORS
!	EIGENVAL (OUTPUT, REAL, 1-DIMENSIONAL [1:MAXDIM]): EIGENVALUES OF THE TOTAL HAMILTONIAN

implicit none
real*16, allocatable, dimension (:) :: DE
real*8, allocatable, dimension (:) :: eigenval
integer, allocatable, dimension (:) :: ndummy
real*16 :: ex,ey,ez,lambda,beta=0.46686437,B
complex*16, allocatable, dimension (:,:) :: SOC2,Lx2,Ly2,Lz2,Sx2,Sy2,Sz2,HAM
complex*16, dimension (1:200) :: WORK
integer :: i,j,maxdim,diagcount,a,LWORK,INFO
real*8, dimension(1:250) :: RWORK
character (len=1) :: JOBZ,UPLO
external zheev
JOBZ='V'
UPLO='U'
LWORK=1000
diagcount=0
a=1
allocate (HAM(1:maxdim,1:maxdim))
allocate (eigenval(1:maxdim))
do i=1,maxdim
	eigenval(i)=0.0
	do j=1,maxdim
		HAM(i,j)=lambda*SOC2(i,j)+beta*B*(ex*(Lx2(i,j)+2.00*Sx2(i,j))+ey*(Ly2(i,j)+2.0*Sy2(i,j))+ez*(Lz2(i,j)+2.0*Sz2(i,j)))	! Setting up all parts of the H except the orbital excitation energies (DE)
		if (i==j) then								! DE only appears on diagonal elements
			HAM(i,i)=HAM(i,i)+cmplx(DE(a),0.0)
			diagcount=diagcount+1
			if (diagcount<ndummy(a)) then
				a=a
			elseif (diagcount==ndummy(a)) then
				a=a+1
				diagcount=0
			endif
		endif
	enddo
enddo
!!!!!!!!!!!! NOW THE HAMILTONIAN HAS BEEN SET UP. CALLING ZHEEV (EXTERNAL FUNCTION)

call zheev (JOBZ,UPLO,maxdim,HAM,maxdim,eigenval,WORK,LWORK,RWORK,INFO)	! EXTERNAL FUNCTION: RETURNS EIGENVECTORS AND EIGENVALUES (LAPACK LIBRARY)


end subroutine

subroutine multmat(mat1,mat2,mat3,maxdim)
!MATRIX MULTIPLICATION SUBROUTINE MAT1 X MAT2 = MAT3
!	MAT1,MAT2 (INPUT, COMPLEX, 2-DIMENSIONAL [1:MAXDIM]): INPUT MATRICES TO BE MULTIPLICATED MAT1 X MAT2 = MAT3
!	MAT3 (OUTPUT, COMPLEX, 2-DIMENSIONAL [1:MAXDIM]): RESULT MATRIX OF THE OPERATION MAT1 X MAT2 = MAT3
!	MAXDIM (INPUT, INTEGER) : DIMENSION OF THE MATRICES

implicit none
integer :: i,j,k,l,maxdim
complex*16, dimension (1:maxdim,1:maxdim) :: mat1, mat2, mat3


do i=1,maxdim
        do j=1,maxdim
                mat3(i,j)=0.0
        enddo
enddo



do i=1,maxdim
        do j=1,maxdim
                do k=1,maxdim
                                mat3(i,j)=mat3(i,j)+mat1(i,k)*mat2(k,j)


                enddo
        enddo
enddo


end subroutine

subroutine mos_calc(B,DEq,eta,d,T,limit,DE,SOC2,Lx2,Ly2,Lz2,Sx2,Sy2,Sz2,SDx2,SDy2,SDz2,maxdim,ndummy,lambda,gridlevel,Imax,lorentz&
&,rd,AFC,extend,resolution)
! CALCULATES THE POWDER-AVERAGED MOSSBAUER SPECTRUM FROM THE VARIOUS PARAMETERS. IN DETAILS, FOR EACH DIRECTION OF THE APPLIED FIELD, 
! THE HYPERFINE FIELDS ARE CALCULATED FROM THE ELECTRONIC EIGENSTATES AND EIGENVALUES AFTER DIAGONALIZING THE ELECTRONIC HAMILTONIAN (SEE README FILE). 
! THEN, THE NUCLEAR HAMILTONIAN IS DIAGONALIZED. THE APPLIED-FIELD DIRECTION- AND GAMMA-WAVE PROPAGATION-VECTOR- DEPENDENT ABSORPTION FUNCTION IS CALCULATED. THEN, THE ABSORPTION FUNCTION
! IS AVEARGED OVER ALL POSSIBLE DIRECTIONS OF THE PROPAGATION VECTOR AND APPLIED FIELD. FINALLY, THE RESULTING FUNCTION IS SUBSTRACTED TO 1.0, TO GET THE USUAL TRANSMISSION MODE
!	B (INPUT, REAL): NORM OF THE APPLIED FIELD 
!	DEQ (INPUT,REAL): QUADRUPOLE SPLITTING VALUE
!	ETA (INPUT,REAL): ISOMER SHIFT
!	DE (INPUT, REAL, 1-DIMENSION [1:NSTATE]): THE VERTICAL EXCITATION ENERGIES ASSOCIATED WITH EACH ELECTRONIC CONFIGURATION, RELATIVE TO THE FIRST DEFINED ELECTRONIC CONFIGURATION.
!       SOC2,LX2,LY2,LZ2,SX2,SY2,SZ2,SDx2,SDy2,SDz2 (INPUT, COMPLEX, 2-DIMENSIONAL [1:MAXDIM,1:MAXDIM]): SPIN-ORBIT COUPLING, X,Y,Z COMPONENTS OF THE ORBITAL SPIN MOMENTA AND SPIN-DIPOLE MATRICES.
!	MAXDIM (INPUT,INTEGER): DIMENSION OF THE HAMILTONIAN
!	NDUMMY (INPUT,INTEGER, 1-DIMENSIONAL[1:NSTATE]): MULTIPLICITY OF EACH CONFIGURATION
!	LAMBDA (INPUT,REAL): SPIN-ORBIT COUPLING CONSTANT OF THE IRON
!	GRIDLEVEL (INPUT, CHARACTER): TYPE AND ACCURACY OF THE QUADRATURE.'LEB' FOR LEBEDEV QUADRATURE, WITH 86 NODES FOR THE DIRECTION OF THE APPLIED FIELD, PLUS 50 NODES FOR THE GAMMA-WAVE
!PROPAGATION VECTOR, WHICH ACCOUNTS FOR A TOTAL OF 4300 SPECTRA TO CALCULATE. 'NORMAL' CORRESPONDS TO A NUMERICAL INTEGRATION WITH AN INCREMENT OF POLAR ANGLE BY 0.02*PI, INCREMENT
!OF AZIMUTHAL ANGLE OF 0.05*PI AND AN INCREMENT OF 0.02*PI FOR THE ANGLE DESCRIBING THE POSITION OF THE GAMMA-WAVE PROPAGATION VECTOR. THIS AMOUNTS TO A TOTAL OF 100 000 SPECTRA TO CALCULATE.
! 'FINE' CORRESPONDS TO AN INCREMENT OF THE POLAR ANGLE BY 0.01*PI, AN INCREMENT OF THE AZIMUTHAL ANGLE BY 0.02*PI AND OF THE GAMMA-WAVE PROPAGATION VECTOR ANGLE OF 0.01*PI. THIS AMOUNTS TO A TOTAL OF 
!10**6 SPECTRA TO CALCULATE! UNFORTUNATELY, AS OF NOW, THE LEBEDEV GRID LEADS TO THE APPARITION OF INTRUDER PEAKS, AND A FINE GRID LEVEL IS REQUIRED FOR A FINAL ACCURATE CALCULATION OF THE SPECTRUM
!	IMAX (INPUT,REAL): THE MAXIMUM INTENSITY OF THE SPECTRUM (ARBITRARY UNIT)
!	LORENTZ (INPUT, REAL, 1-DIMENSIONAL[1:6]): SPECTRUM WIDTH. AS IS CODED, SAME WIDTH FOR ALL PEAKS, BUT THIS CAN BE MODIFIED BY SPECIFYING INDIVIDUAL VALUES IN THE CODE
!	RD (INPUT,REAL): AVERAGE RD-3 VALUE,USED TO CALCULATE THE HYPERFINE FIELDS
!	AFC (INPUT, REAL): PROPORTIONALITY FACTOR BETWEEN THE FERMI-CONTACT FIELD AND SPIN MOMENTUM VECTOR
!	EXTEND (INPUT, INTEGER): NUMBER OF POINT CONSITUTING HALF WIDTH OF THE SPECTRUM
!	RESOLUTION (INPUT, REAL): DISTANCE BETWEEN TWO CALCULATED POINTS
!	



implicit none
real*16 :: B,DEq,eta,d,T,lambda,pi=3.1415927,kB=0.695039,increment,incrementprime,increment2,ex,ey,ez,rd,AFC,Z
real*8 :: exg,eyg,ezg, phi,theta,resolution,Imax
real :: start,finish
character (len=100) :: limit,gridlevel
complex*16, allocatable, dimension (:,:) :: SOC2,Lx2,Ly2,Lz2,Sx2,Sy2,Sz2,SDx2,SDy2,SDz2,HAM,result_micro
complex*16, dimension (1:6,1:6) :: Hnuc
complex*16 :: Bintx, Binty, Bintz
real*16, allocatable, dimension (:) :: DE
real*16, dimension (0:85) :: Coeff
real*16, dimension (0:85,1:3) :: emat
real*8, allocatable, dimension (:) :: eigenval
real*8, dimension (1:8, 1:2) :: intens
real*8, dimension (1:8) :: lorentz 	! A UNIQUE VALUE IS USED BY DEFAULT, BUT THEY CAN BE MODULATED HERE 
real*8, dimension (1:6) :: Nuc_eigenval
complex*8, dimension (-1000:1000) :: dummyspectrum, spectrumb, powderspec
integer :: maxdim,u,v,i,j,k,w,extend,ll
integer, allocatable, dimension (:) :: ndummy

allocate (result_micro(1:maxdim,1:9))
if (gridlevel=='normal') then		! setting up the number and type of nodes for the powder-integration
	u=50
	v=40
	w=50
	increment=0.02*pi
	incrementprime=0.05*pi
	increment2=0.02*pi
elseif (gridlevel=='fine') then
	u=100
        v=100
        w=100
        increment=0.01*pi
        incrementprime=0.02*pi
        increment2=0.01*pi

elseif (gridlevel=='leb') then		 
	u=85
	v=0
	w=50
	call mosleb(Coeff,emat)		! sets up a tight set of node for lebedev quadrature
	increment2=0.02*pi		
endif

do i=1,8
	do j=1,2
		intens(i,j)=0.0
	enddo
enddo

do i=-extend,extend
	powderspec(i)=0.0
enddo

do i=0,u
!	print*, i, '/',u
	do j=0,v
		ez=cos(real(i)*increment)
		ex=sin(real(i)*increment)*cos(real(j)*incrementprime)
		ey=sin(real(i)*increment)*sin(real(j)*incrementprime)
		if (gridlevel=='leb') then
			ex=emat(i,1)
			ey=emat(i,2)
			ez=emat(i,3)
		endif
		call diag_Ham(DE,SOC2,Lx2,Ly2,Lz2,Sx2,Sy2,Sz2,ex,ey,ez,maxdim,ndummy,lambda,B,HAM,eigenval)	! Calculates and diagonalizes the electronic Hamiltonian
		call mag_micro(HAM,maxdim,ndummy,Lx2,Ly2,Lz2,Sx2,Sy2,Sz2,SDx2,SDy2,SDz2,result_micro)	! calculates the expectation values of the S,L and spin-dipole operators 
		Bintx=0.0
		Binty=0.0
		Bintz=0.0
		Z=0.0
		do k=1,maxdim		! for fast electron relaxation (only implemented for now). Calculates the Boltzmann averaged hyperfine field acting on the nucleus
			Bintx=Bintx+(12.5*rd*(result_micro(k,7)-result_micro(k,1))+AFC*result_micro(k,4))*exp(-(eigenval(k)-eigenval(1))/(kB*T))
			Binty=Binty+(12.5*rd*(result_micro(k,8)-result_micro(k,2))+AFC*result_micro(k,5))*exp(-(eigenval(k)-eigenval(1))/(kB*T))
			Bintz=Bintz+(12.5*rd*(result_micro(k,9)-result_micro(k,3))+AFC*result_micro(k,6))*exp(-(eigenval(k)-eigenval(1))/(kB*T))
			Z=Z+exp(-(eigenval(k)-eigenval(1))/(kB*T))
		enddo
		Bintx=Bintx/Z
		Binty=Binty/Z
		Bintz=Bintz/Z

		call diag_nuc(DEq,eta,d,Hnuc,Nuc_eigenval,B,ex,ey,ez,Bintx,Binty,Bintz)		! diagonalizes the nuclear Hamiltonian
		
		theta=real(i)*increment
		phi=real(j)*incrementprime
		if (gridlevel=='leb') then
			call findangles(ex,ey,ez,theta,phi)	! for lebedev grid only, the angle is not directly accessible and must be calculated from ex,ey and ez
		endif			
			
		do k=-extend,extend
			spectrumb(k)=0.0
			dummyspectrum(k)=0.0
		enddo
		do k=0,w
			exg=cos(real(k)*increment2)			!! Calculates the directional vector in the lab frame (B applied along z)
			eyg=sin(real(k)*increment2)
			ezg=0.0
			call rotate(exg,eyg,ezg,phi,theta) 		!!! Rotates the lab frame to the molecular frame
			call intensitybuilder(Hnuc,exg,eyg,ezg,Nuc_eigenval,intens)	! Calculates the intensities from the eigenstates of the nuclear Hamiltonian
			call spectrumbuilder(intens,lorentz,resolution,extend,Imax,dummyspectrum)	! Calculates the absorption function from transition energies and intensities 
			do ll=-extend,extend		! averages the absorption function over all possible directions of the gamma-wave propagation vector for a given applied field direction
							! in the molecular framework
				spectrumb(ll)=spectrumb(ll)+dummyspectrum(ll)*increment2/(2.0*pi)
			enddo
		enddo
		do k=-extend,extend		! Averages the function over all possible directions of the applied field in the molecular framework.
			if (gridlevel=='leb') then
				powderspec(k)=powderspec(k)+spectrumb(k)*Coeff(i)*10.0 !the factor 10 doesn't matter since the overall absorption unit is arbitrary (introduced so Lebedev and non-lebedev
											! grids are easily comparable)
			else
				powderspec(k)=powderspec(k)+spectrumb(k)*increment*incrementprime*sin(theta)
			endif
		enddo
		deallocate (HAM)
		deallocate (eigenval)
	enddo
enddo

do i=-extend,extend
	print*, 1.00-real(powderspec(i))	! switches to transmission mode
enddo

end subroutine

subroutine findangles(ex,ey,ez,theta,phi)		!! Numerical instabilities relative to the arccos and arcsin made it difficult to perform this task intelligently
! CALCULATES THE POLAR AND AZIMUTHAL ANGLES ASSOCIATED WITH EX,EY AND EZ
!	ex,ey,ez (INPUT, REAL): DIRECTIONAL VECTORS OF THE APPLIED FIELDS IN HTE MOLECULAR FRAMEWORK
!	THETA,PHI (OUTPUT,REAL): POLAR AND AZIMUTHAL ANGLE IN THE MOLECULAR FRAMEWORK 

implicit none					 
real*16 :: ex,ey,ez,pi=3.1415927
real*8 :: theta, phi,phi1,phi1prime,phi1second,phi2,phi2prime
integer :: i

theta=acos(ez)
phi=0.0
i=0

if (abs(theta)<0.0001.or.abs(theta-pi)<0.0001) then
	phi=0.0
else
	phi1=ex/sin(theta)
	if (phi1>1.0) then		! Was necessary because in extreme cases ex/sin(theta) were slightly above 1 because of numerical instability. Acos then returns NaN
		phi1=1.0
	elseif (phi1<-1.0) then
		phi1=-1.0
	endif	
	phi1=acos(phi1)
	phi1second=-phi1
	if (abs(sin(theta)*sin(phi1)-ey)<0.01) then
		phi=phi1
	elseif (abs(sin(theta)*sin(phi1second)-ey)<0.01) then
		phi=phi1second
	else
		print*, 'error!'
		print*, ex,ey,ez, phi1
		print*, sin(theta)*sin(phi1)
		stop
	endif
endif

end subroutine

subroutine mosleb2(Coeff,emat)
!! TIGHTER SET OF NODES FOR LEBEDEV QUADRATURE. NOT USED IN THE PRESENT VERSION



implicit none
real*16, dimension (0:109) :: Coeff
real*16, dimension (0:109,1:3) :: emat    ! x=1, y=2, z=3
real*16 :: p,q,m1,l1,m2,l2,l3,m3
integer :: i

p=0.878159
q=0.478369
m1=0.965124
m2=0.840256
m3=0.21596
l1=0.185116
l2=0.383386
l3=0.690421



emat(0,1)=0.0
emat(0,2)=0.0
emat(0,3)=1.0
emat(1,1)=0.0
emat(1,2)=0.0
emat(1,3)=-1.0
emat(2,1)=0.0
emat(2,2)=1.0
emat(2,3)=0.0
emat(3,1)=0.0
emat(3,2)=-1.0
emat(3,3)=0.0
emat(4,1)=1.0
emat(4,2)=0.0
emat(4,3)=0.0
emat(5,1)=-1.0
emat(5,2)=0.0
emat(5,3)=0.0
emat(6,1)=1.0/(3.0**0.5)
emat(6,2)=1.0/(3.0**0.5)
emat(6,3)=1.0/(3.0**0.5)
emat(7,1)=-1.0/(3.0**0.5)
emat(7,2)=1.0/(3.0**0.5)
emat(7,3)=1.0/(3.0**0.5)
emat(8,1)=1.0/(3.0**0.5)
emat(8,2)=-1.0/(3.0**0.5)
emat(8,3)=1.0/(3.0**0.5)
emat(9,1)=1.0/(3.0**0.5)
emat(9,2)=1.0/(3.0**0.5)
emat(9,3)=-1.0/(3.0**0.5)
emat(10,1)=-1.0/(3.0**0.5)
emat(10,2)=-1.0/(3.0**0.5)
emat(10,3)=1.0/(3.0**0.5)
emat(11,1)=-1.0/(3.0**0.5)
emat(11,2)=1.0/(3.0**0.5)
emat(11,3)=-1.0/(3.0**0.5)
emat(12,1)=1.0/(3.0**0.5)
emat(12,2)=-1.0/(3.0**0.5)
emat(12,3)=-1.0/(3.0**0.5)
emat(13,1)=-1.0/(3.0**0.5)
emat(13,2)=-1.0/(3.0**0.5)
emat(13,3)=-1.0/(3.0**0.5)


emat(14,1)=l3
emat(14,2)=l3
emat(14,3)=m3
emat(15,1)=l3
emat(15,2)=l3
emat(15,3)=-m3
emat(16,1)=l3
emat(16,2)=-l3
emat(16,3)=m3
emat(17,1)=-l3
emat(17,2)=l3
emat(17,3)=m3
emat(18,1)=l3
emat(18,2)=-l3
emat(18,3)=-m3
emat(19,1)=-l3
emat(19,2)=l3
emat(19,3)=-m3
emat(20,1)=-l3
emat(20,2)=-l3
emat(20,3)=m3
emat(21,1)=-l3
emat(21,2)=-l3
emat(21,3)=-m3
emat(22,1)=l3
emat(22,2)=m3
emat(22,3)=l3
emat(23,1)=l3
emat(23,2)=m3
emat(23,3)=-l3
emat(24,1)=l3
emat(24,2)=-m3
emat(24,3)=l3
emat(25,1)=-l3
emat(25,2)=m3
emat(25,3)=l3
emat(26,1)=l3
emat(26,2)=-m3
emat(26,3)=-l3
emat(27,1)=-l3
emat(27,2)=m3
emat(27,3)=-l3
emat(28,1)=-l3
emat(28,2)=-m3
emat(28,3)=l3
emat(29,1)=-l3
emat(29,2)=-m3
emat(29,3)=-l3
emat(30,1)=m3
emat(30,2)=l3
emat(30,3)=l3
emat(31,1)=m3
emat(31,2)=l3
emat(31,3)=-l3
emat(32,1)=m3
emat(32,2)=-l3
emat(32,3)=l3
emat(33,1)=-m3
emat(33,2)=l3
emat(33,3)=l3
emat(34,1)=m3
emat(34,2)=-l3
emat(34,3)=-l3
emat(35,1)=-m3
emat(35,2)=l3
emat(35,3)=-l3
emat(36,1)=-m3
emat(36,2)=-l3
emat(36,3)=l3
emat(37,1)=-m3
emat(37,2)=-l3
emat(37,3)=-l3

emat(38,1)=l1
emat(38,2)=l1
emat(38,3)=m1
emat(39,1)=-l1
emat(39,2)=l1
emat(39,3)=m1
emat(40,1)=l1
emat(40,2)=-l1
emat(40,3)=m1
emat(41,1)=l1
emat(41,2)=l1
emat(41,3)=-m1
emat(42,1)=-l1
emat(42,2)=-l1
emat(42,3)=m1
emat(43,1)=-l1
emat(43,2)=l1
emat(43,3)=-m1
emat(44,1)=l1
emat(44,2)=-l1
emat(44,3)=-m1
emat(45,1)=-l1
emat(45,2)=-l1
emat(45,3)=-m1
emat(46,1)=l1
emat(46,2)=m1
emat(46,3)=l1
emat(47,1)=-l1
emat(47,2)=m1
emat(47,3)=l1
emat(48,1)=l1
emat(48,2)=-m1
emat(48,3)=l1
emat(49,1)=l1
emat(49,2)=m1
emat(49,3)=-l1
emat(50,1)=-l1
emat(50,2)=-m1
emat(50,3)=l1
emat(51,1)=-l1
emat(51,2)=m1
emat(51,3)=-l1
emat(52,1)=l1
emat(52,2)=-m1
emat(52,3)=-l1
emat(53,1)=-l1
emat(53,2)=-m1
emat(53,3)=-l1
emat(54,1)=m1
emat(54,2)=l1
emat(54,3)=l1
emat(55,1)=-m1
emat(55,2)=l1
emat(55,3)=l1
emat(56,1)=m1
emat(56,2)=-l1
emat(56,3)=l1
emat(57,1)=m1
emat(57,2)=l1
emat(57,3)=-l1
emat(58,1)=-m1
emat(58,2)=-l1
emat(58,3)=l1
emat(59,1)=-m1
emat(59,2)=l1
emat(59,3)=-l1
emat(60,1)=m1
emat(60,2)=-l1
emat(60,3)=-l1
emat(61,1)=-m1
emat(61,2)=-l1
emat(61,3)=-l1


emat(62,1)=l2
emat(62,2)=l2
emat(62,3)=m2
emat(63,1)=-l2
emat(63,2)=l2
emat(63,3)=m2
emat(64,1)=l2
emat(64,2)=-l2
emat(64,3)=m2
emat(65,1)=l2
emat(65,2)=l2
emat(65,3)=-m2
emat(66,1)=-l2
emat(66,2)=-l2
emat(66,3)=m2
emat(67,1)=-l2
emat(67,2)=l2
emat(67,3)=-m2
emat(68,1)=l2
emat(68,2)=-l2
emat(68,3)=-m2
emat(69,1)=-l2
emat(69,2)=-l2
emat(69,3)=-m2
emat(70,1)=l2
emat(70,2)=m2
emat(70,3)=l2
emat(71,1)=-l2
emat(71,2)=m2
emat(71,3)=l2
emat(72,1)=l2
emat(72,2)=-m2
emat(72,3)=l2
emat(73,1)=l2
emat(73,2)=m2
emat(73,3)=-l2
emat(74,1)=-l2
emat(74,2)=-m2
emat(74,3)=l2
emat(75,1)=-l2
emat(75,2)=m2
emat(75,3)=-l2
emat(76,1)=l2
emat(76,2)=-m2
emat(76,3)=-l2
emat(77,1)=-l2
emat(77,2)=-m2
emat(77,3)=-l2
emat(78,1)=m2
emat(78,2)=l2
emat(78,3)=l2
emat(79,1)=-m2
emat(79,2)=l2
emat(79,3)=l2
emat(80,1)=m2
emat(80,2)=-l2
emat(80,3)=l2
emat(81,1)=m2
emat(81,2)=l2
emat(81,3)=-l2
emat(82,1)=-m2
emat(82,2)=-l2
emat(82,3)=l2
emat(83,1)=-m2
emat(83,2)=l2
emat(83,3)=-l2
emat(84,1)=m2
emat(84,2)=-l2
emat(84,3)=-l2
emat(85,1)=-m2
emat(85,2)=-l2
emat(85,3)=-l2

emat(86,1)=p
emat(86,2)=q
emat(86,3)=0.0
emat(87,1)=p
emat(87,2)=-q
emat(87,3)=0.0
emat(88,1)=-p
emat(88,2)=q
emat(88,3)=0.0
emat(89,1)=-p
emat(89,2)=-q
emat(89,3)=0.0
emat(90,1)=p
emat(90,2)=0.0
emat(90,3)=q
emat(91,1)=p
emat(91,2)=0.0
emat(91,3)=-q
emat(92,1)=-p
emat(92,2)=0.0
emat(92,3)=q
emat(93,1)=-p
emat(93,2)=0.0
emat(93,3)=-q
emat(94,1)=0.0
emat(94,2)=p
emat(94,3)=q
emat(95,1)=0.0
emat(95,2)=p
emat(95,3)=-q
emat(96,1)=0.0
emat(96,2)=-p
emat(96,3)=q
emat(97,1)=0.0
emat(97,2)=-p
emat(97,3)=-q

emat(98,1)=p
emat(98,2)=q
emat(98,3)=0.0
emat(99,1)=p
emat(99,2)=-q
emat(99,3)=0.0
emat(100,1)=-p
emat(100,2)=q
emat(100,3)=0.0
emat(101,1)=-p
emat(101,2)=-q
emat(101,3)=0.0
emat(102,1)=p
emat(102,2)=0.0
emat(102,3)=q
emat(103,1)=p
emat(103,2)=0.0
emat(103,3)=-q
emat(104,1)=-p
emat(104,2)=0.0
emat(104,3)=q
emat(105,1)=-p
emat(105,2)=0.0
emat(105,3)=-q
emat(106,1)=0.0
emat(106,2)=p
emat(106,3)=q
emat(107,1)=0.0
emat(107,2)=p
emat(107,3)=-q
emat(108,1)=0.0
emat(108,2)=-p
emat(108,3)=q
emat(109,1)=0.0
emat(109,2)=-p
emat(109,3)=-q












do i=0,109
        if (i<=5) then
                Coeff(i)=0.003828
        elseif (i>5.and.i<=13) then
                Coeff(i)=0.009886
        elseif (i>13.and.i<=37) then
                Coeff(i)=0.009943
        elseif (i>37.and.i<=61) then
                Coeff(i)=0.008441
        elseif (i>61.and.i<=85) then
                Coeff(i)=0.009595
	elseif (i>85.and.i<=109) then
		Coeff(i)=0.009695
        endif

enddo






end subroutine

subroutine intensitybuilder(mos,exg,eyg,ezg,eigenval,intens)
! CALCULATES THE 8 TRANSITION INTENSITIES AND ENERGIES FROM THE EIGENVECTORS OF THE NUCLEAR HAMILTONIAN AND THE DIRECTIONAL VECTORS OF THE INCIDENT GAMMA-WAVE, STORED IN INTENS
! IN DETAILS, PROVIDES THE RIGHT INPUT FORM FOR THE SUBROUTINE INT57, WHICH ACTUALLY CALCULATES THE INTENSITIES. 
! THE SUBROUTINE WAS INITIALLY WRITTEN IN FORTRAN77 BY A PREVIOUS AUTHOR. AT THE TIME THERE WAS NO COMPLEX NUMBERS AND IMAGINARY/REAL PARTS
! HAD TO BE STORED IN SEPARATE VARIABLES
!	MOS (INPUT, COMPLEX, 2-DIMENSION [1:6,1:6]): EIGENSTATES OF THE NUCLEAR HAMILTONIAN
!	EXG,EYG,EZG (INPUT,REAL): DIRECTIONAL VECTORS OF THE GAMMA-WAVE PROPAGATION VECTOR
!	EIGENVAL (INPUT,REAL, 1-DIMENSION [1:6]): EIGENVALUES OF THE NUCLEAR HAMILTONIAN
!	INTENS (OUTPUT, REAL, 2-DIMENSION [1:8,1:2]): INTENS(I,J) RETURNS THE TRANSITION INTENSITY (J=1) AND ENERGY (J=2) OF THE ITH TRANSITION 
implicit none
complex*16, dimension (1:6, 1:6) :: mos
real*8, dimension (1:6) :: eigenval
real*8 :: C1, C2R, C2I,a,b,exgp,eygp,ezgp
real*8 :: exg,eyg,ezg
real*8, dimension (4,4) :: UR, UI
integer :: k,i,flag,IPOL,j
real*8, dimension(8) :: XNT
real*8, dimension (1:8,1:2) :: intens
flag=0


!! The following procedure separates imaginary and real parts of the eigenvectors
if (abs(real(mos(6,1)))<0.0001) then
        a=0.0
        b=1.0
        flag=1
        do i=1,6
              mos(i,1)=mos(i,1)*cmplx(a,b)
        enddo
else
        do while (abs(aimag(mos(6,1)))>0.0001)
                a=(real(mos(6,1))**2.0)/((aimag(mos(6,1))**2.0)+(real(mos(6,1))**2.0))
                b=(1.0-a**2.0)**0.5
                flag=1
                do i=1,6
                        mos(i,1)=mos(i,1)*cmplx(a,b)
                enddo
        enddo
endif

if (abs(real(mos(5,2)))<0.0001) then
        a=0.0
        b=1.0
        flag=1
        do i=1,6
                mos(i,2)=mos(i,2)*cmplx(a,b)
        enddo
else
        do while (abs(aimag(mos(5,2)))>0.0001)
                 a=(real(mos(5,2))**2.0)/((aimag(mos(5,2))**2.0)+(real(mos(5,2))**2.0))
                 b=(1.0-a**2.0)**0.5
                 flag=1
                do i=1,6
                   mos(i,2)=mos(i,2)*cmplx(a,b)
                enddo
        enddo
endif



C1=real(mos(6,1))
C2R=real(mos(5,1))
C2I=aimag(mos(5,1))
!print*, mos(6,1), mos(5,2)
!print*, mos(5,1), mos(6,2)
do i=1,4
        do j=1,4
                UR(j,i)=real(mos(5-j,i+2))
                UI(j,i)=aimag(mos(5-j,i+2))
        enddo
enddo
IPOL=0
call INT57(UR,UI,C1,C2R,C2I,exg,eyg,ezg,IPOL,XNT)	! calculates the intensities
do i=1,4			! stores the intensities in intens
        intens(i,1)=XNT(i)
enddo
!!!! the energies are now stored in intens. The -100 removes an artificial separation of the I=1/2 and I=3/2 states introduced in diag_nuc
intens(1,2)=eigenval(3)-eigenval(1)-100    ! 1->3
intens(2,2)=eigenval(4)-eigenval(1)-100    ! 1->4
intens(3,2)=eigenval(5)-eigenval(1)-100    ! 1->5
intens(4,2)=eigenval(6)-eigenval(1)-100    ! 1->6
intens(5,2)=eigenval(3)-eigenval(2)-100    !2->3
intens(6,2)=eigenval(4)-eigenval(2)-100    ! 2->4
intens(7,2)=eigenval(5)-eigenval(2)-100    ! 2->5
intens(8,2)=eigenval(6)-eigenval(2)-100    ! 2->6



do i=1,8
        intens(i,1)=XNT(i)
!        print*, XNT(i), intens(i,2)
enddo

end subroutine

!*****************************************************************
      SUBROUTINE INT57(UR,UI,C1,C2R,C2I,EX,EY,EZ,IPOL,XNT)
!*****************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)

!    ***
!    ***  THE SUBROUTINE INT57 CALCULATES THE INTENSITIES FROM THE
!    ***  57 FE NUCLEAR EIGENSTATES AND THE DIRECTION COSINES OF THE
!    ***  gamm-RAY.
!    ***
!    ***  INPUT PARAMETERS
      DIMENSION UR(4,4),UI(4,4)
!    ***  UR(I,K) AND UI(I,K) ARE THE REAL PART AND THE IMAGINARY PART
!    ***   OF THE I-TH COMPONENT OF THE EIGENVECTOR K FOR THE 3/2 STATE
!    ***  C1,C2R,C2I=COEFFICIENTS OF THE EIGENSTATES OF THE 57FE NUCLEAR
!    ***   GROUND STATE, AS DEFINED IN SUBROUTINE N57G I.E.
!    ***         1)=C1 *  1/2) + (C2R+I*C2I) *  -1/2)
!    ***         2)=-(C2R-I*C2I) *  1/2) + C1 *  -1/2)
!    ***        I=DSQRT(-1) IN THE ABOVE FORMULA
!    ***  EX,EY,EZ=VECTOR CHARACTERIZING THE ELECTROMAGNETIC WAVE. TO
!    ***   BE USED AS DEFINED BY IPOL
!    ***  IPOL  CHARACTERIZES THE TYPE OF INTENSITY CALCULATION
!    ***   IPOL=-3   INTENSITY TENSOR AND INTENSITY VECTOR ARE TO BE
!    ***             CALCULATED. NO INPUT FOR EX,EY,EZ NECESSARY
!    ***   IPOL=-2   POWDER MEAN VALUE; NO INPUT FOR EX,EY,EZ NECESSARY
!    ***   IPOL=-1   PERPENDICULAR CONFIGURATION (SPLIT COIL)
!    ***             EX,EY,EZ IS TO BE PARALLEL TO THE EXTERNAL FIELD
!    ***   IPOL=0    UNPOLARIZED gamm-RAY PARALLEL TO EX,EY,EZ
!    ***   IPOL=1    CIRCULAR POLARISATION,WAVE VECTOR PARALLEL EX,EY,EZ
!    ***   IPOL=2    LINEAR POLARISATION; IN 57FE (EX,EY,EZ) IS THE VECTOR
!    ***             PERPENDICULAR TO THE WAVE VECTOR AND PERPENDICULAR TO
!    ***             THE VECTOR POTENTIAL A
!    ***
!    ***  OUTPUT PARAMETERS
      DIMENSION XNT(1)
!    ***  THE INTENSITIES BELONGING TO THE  1) TRANSITION ARE STORED
!    ***  IN XNT(1 4), THOSE BELONGING TO THE  2) TRANSITION
!    ***  IN XNT(5 8)
!    ***  FOR IPOL=-3   INTENSITY TENSOR IS STORED IN THE FOLLOWING WAY
!    ***  IXX IN 1 8, IYY IN 9 16, IZZ IN 17 24,IXY IN 25 32,IXZ IN 33 40,
!    ***  IYZ IN 41 48,THE INTENSITY VECTOR IS STORED IN THE FOLLOWING WAY
!    ***  IX IN 49 56, IY IN 57 64, IZ IN 65 72
!    ***
      SQRT3=1.732050808
      IF (IPOL.LT.-1) GOTO 5
      HXX=EX*EX
      HYY=EY*EY
      HZZ=EZ*EZ
      HXY=EX*EY
      HXZ=EX*EZ
      HYZ=EY*EZ
      ESQ=HXX+HYY+HZZ
      IF (IPOL.EQ.1) ESQR=SQRT(ESQ)
   5  DO 10 K=1,4
      UR1=UR(1,K)
      UR2=UR(2,K)
      UR3=UR(3,K)
      UR4=UR(4,K)
      UI1=UI(1,K)
      UI2=UI(2,K)
      UI3=UI(3,K)
      UI4=UI(4,K)
      VR2=C2R*UR2+C2I*UI2
      VR3=C2R*UR3+C2I*UI3
      VR4=C2R*UR4+C2I*UI4
      VI2=C2R*UI2-C2I*UR2
      VI3=C2R*UI3-C2I*UR3
      VI4=C2R*UI4-C2I*UR4
      R1=C1*UR1*SQRT3+VR2
      R2=C1*UR2+VR3
      R3=C1*UR3+VR4*SQRT3
      U1=C1*UI1*SQRT3+VI2
      U2=C1*UI2+VI3
      U3=C1*UI3+VI4*SQRT3
      I=K
  20  CONTINUE
      X22=R2*R2+U2*U2
      X33=R3*R3+U3*U3
      XR23=R2*R3+U2*U3
      XI23=U2*R3-R2*U3
      X11=R1*R1+U1*U1
      XR12=R1*R2+U1*U2
      XR13=R1*R3+U1*U3
      XI12=U1*R2-R1*U2
      XI13=U1*R3-R1*U3
      H0=(X11+X33)/4.0
      H1=H0+X22
      H2=XR13/2.0
      RXX=H1+H2
      RYY=H1-H2
      RZZ=2.0*H0
      IF (IPOL.NE.-2) GOTO 60
      XNTEN=(RXX+RYY+RZZ)/3.
      GOTO 50
  60  CONTINUE
      RXY=-XI13
      RXZ=XR12-XR23
      RYZ=-XI12+XI23
      IF (IPOL.NE.-3) GOTO 48
      XNT(I)=RXX
      II=I+8
      XNT(II)=RYY
      II=I+16
      XNT(II)=RZZ
      II=I+24
      XNT(II)=RXY/2.
      II=I+32
      XNT(II)=RXZ/2.
      II=I+40
      XNT(II)=RYZ/2.
      II=I+48
      XNT(II)=XR12+XR23
      II=I+56
      XNT(II)=-XI12-XI23
      II=I+64
      XNT(II)=(X11-X33)/2.0
      GOTO 53
  48  XNTEN=(RXX*HXX+RYY*HYY+RZZ*HZZ+RXY*HXY+RXZ*HXZ+RYZ*HYZ)/ESQ
      IF (IPOL) 49,50,51
  49  XNTEN=0.5*(RXX+RYY+RZZ-XNTEN)
      GOTO 50
  51  IF (IPOL.EQ.2) GOTO 52
      RX=XR12+XR23
      RY=-XI12-XI23
      RZ=(X11-X33)/2.0
      XNTEN=XNTEN+(RX*EX+RY*EY+RZ*EZ)/ESQR
      GOTO 50
  52  XNTEN=RXX+RYY+RZZ-2.*XNTEN
  50  XNT(I)=XNTEN
  53  IF (I.GT.4) GO TO 30
      I=K+4
      VR1=-C2R*UR1+C2I*UI1
      VR2=-C2R*UR2+C2I*UI2
      VR3=-C2R*UR3+C2I*UI3
      VI1=-C2R*UI1-C2I*UR1
      VI2=-C2R*UI2-C2I*UR2
      VI3=-C2R*UI3-C2I*UR3
      R1=VR1*SQRT3+C1*UR2
      R2=VR2+C1*UR3
      R3=VR3+C1*UR4*SQRT3
      U1=VI1*SQRT3+C1*UI2
      U2=VI2+C1*UI3
      U3=VI3+C1*UI4*SQRT3
      IF (I.LT.9) GO TO 20
  30  CONTINUE
  10  CONTINUE
      II=8
      IF (IPOL.EQ.-3) II=72
      DO 100 I=1,II
  100  XNT(I)=XNT(I)/4.
      RETURN
      END

subroutine spectrumbuilder(intens,lorentz,resolution,extend,Imax,spectrum)
! CALCULATES THE SPECTRAL SHAPE OF THE ABSORPTION FUNCTION FROM THE TRANSITION INTENSITIES AND ENERGIES. THE OVERALL FUNCTION IS A SUM OF 8 LORENTZIAN (ONE FOR EACH TRANSITION)
!	INTENS (INPUT, REAL, 2-DIMENSION [1:8,1:2]): INTENS(I,J) RETURNS THE TRANSITION INTENSITY (J=1) AND ENERGY (J=2) OF THE ITH TRANSITION
!	LORENTZ (INPUT,REAL, 1-DIMENSION [1:8]): THE SPECTRAL WIDTH OF EACH TRANSITION INDIVIDUAL LORENTZIAN 
!	RESOLUTION (INPUT, REAL): THE DISTANCE BETWEEN TWO POINTS IN THE ABSORPTION FUNCTION
!	EXTEND (INPUT, INTEGER): THE NUMBER OF POINTS CONSITUTING HALF OF THE SPECTRUM
!	IMAX (INPUT, REAL): THE OVERALL INTENSITY OF THE ABSORPTION FUNCTION (ARBITRARY UNIT)
!	SPECTRUM (OUTPUT, COMPLEX, 1-DIMENSION [-1000:1000]): THE ABSORPTION FUNCTION


implicit none
real*8, dimension (1:8,1:2) :: intens
real*8 :: resolution, Imax
real*8, dimension (1:8) :: lorentz
integer :: extend,i,j,ii
complex*8, dimension (-1000:1000) :: spectrum

do i=-extend,extend
        spectrum(i)=0
enddo


do ii=-extend,extend
        do j=1,8
                spectrum(ii)=spectrum(ii)+Imax*intens(j,1)/(((intens(j,2)-real(ii)*resolution)**2.0)+((0.5*lorentz(j))**2.0))
        enddo
enddo


end subroutine


subroutine rotate(exg,eyg,ezg,phi,teta)
! ROTATES THE DIRECTIONAL VECTORS OF THE INCIDENT GAMMA-WAVE FROM THE LAB FRAME TO THE MOLECULAR FRAME
!	EXG,EYG,EZG (INPUT,OUTPUT): DIRECTIONAL VECTORS IN THE LAB FRAME; RETURNS DIRECTIONAL VECTORS IN THE MOLECULAR FRAME (OUTPUT)
!	TETA, PHI (INPUT) POLAR AND AZIMUTHAL ANGLES OF THE MAGNETIC FIELD DIRECTIONAL VECTOR

implicit none
real*8 :: exg, eyg, ezg,phi,teta
real*16, dimension (1:3,1:3) :: Ry, Rz
real*16, dimension (1:3) :: workvec, eg
integer :: i,j

Ry(1,1)=cos(teta)
Ry(1,2)=0.0
Ry(1,3)=sin(teta)
Ry(2,1)=0.0
Ry(2,2)=1.0
Ry(2,3)=0.0
Ry(3,1)=-sin(teta)
Ry(3,2)=0.0
Ry(3,3)=cos(teta)

Rz(1,1)=cos(phi)
Rz(1,2)=-sin(phi)
Rz(1,3)=0.0
Rz(2,1)=sin(phi)
Rz(2,2)=cos(phi)
Rz(2,3)=0.0
Rz(3,1)=0.0
Rz(3,2)=0.0
Rz(3,3)=1.0

eg(1)=exg
eg(2)=eyg
eg(3)=ezg

do i=1,3
        workvec(i)=0.0
        do j=1,3
                workvec(i)=workvec(i)+Ry(i,j)*eg(j)
        enddo
enddo

do i=1,3
        eg(i)=0.0
        do j=1,3
                eg(i)=eg(i)+Rz(i,j)*workvec(j)
        enddo
enddo

exg=eg(1)
eyg=eg(2)
ezg=eg(3)

end subroutine

subroutine diag_nuc(DEq,eta,d,Hnuc,Nuc_eigenval,B,ex,ey,ez,Bintx,Binty,Bintz)
!CALCULATES AND DIAGONALIZES THE NUCLEAR HAMILTONIAN COMPOSED OF THE ELECTRIC, MAGNETIC HYPERFINE EFFECTS PLUS THE NUCLEAR ZEEMAN EFFECT
!	DEQ (INPUT,REAL) : QUADRUPOLE SPLITTING
!	ETA (INPUT,REAL) : ASYMMETRY PARAMETER
!	D (INPUT,REAL) : ISOMER SHIFT
!	HNUC(OUTPUT,COMPLEX, 2-DIMENSIONAL [1:6,1:6]): BEFORE DIAGONALIZATION, NUCLEAR HAMILTONIAN; AFTER DIAGONALIZATION, CONTAINS THE NUCLEAR EIGENVECTORS
! 	NUC_EIGENVAL (OUTPUT,REAL,1-DIMENSIONAL [1:6]): RETURNS THE NUCLEAR EIGENVALUES 
! 	B (INPUT, REAL): APPLIED FIELD NORM
!	EX,EY,EZ (INPUT,REAL): APPLIED FIELD DIRECTIONAL VECTORS
!	BINTX,BINTY,BINTZ (INPUT,COMPLEX):INTERNAL FIELD VECTORS

implicit none
real*16 :: DEq,eta,d,B,ex,ey,ez,Ae,Ag,Bx,By,Bz
complex*16 :: Bintx,Binty,Bintz
complex*16, dimension (1:6,1:6) :: Hnuc
complex*16, dimension (1:11) :: WORK
real*8, dimension (1:6) :: Nuc_eigenval
real*8, dimension (1:16) :: RWORK
character (len=1) :: JOBZ,UPLO
integer :: i,j,LWORK,INFO
external zheev


LWORK=11
JOBZ='V'
UPLO='U'
!! Filling the Hamiltonian here. Starting from 3/2, -3/2; 3/2, -1/2... and 1/2,-1/2; 1/2,+1/2
Ae=0.068
Ag=-0.119
Bx=B*ex+real(Bintx)
By=B*ey+real(Binty)
Bz=B*ez+real(Bintz)

Hnuc(5,5)=-Ag*0.5*Bz
Hnuc(5,6)=cmplx(0.5*Ag*Bx,0.5*Ag*By)
Hnuc(6,6)=Ag*0.5*Bz
Hnuc(1,1)=100.0+d-Ae*1.5*Bz+DEq/(2.0*((1.0+((eta**2.0)/3.0))**0.5))
Hnuc(1,2)=cmplx(0.866*Ae*Bx,0.866*Ae*By)
Hnuc(1,3)=((3.0**0.5)/6.0)*eta*DEq/((1.0+((eta**2.0)/3.0))**0.5)
Hnuc(1,4)=0.0
Hnuc(2,2)=100.0+d-Ae*0.5*Bz-DEq/(2.0*((1.0+((eta**2.0)/3.0))**0.5))
Hnuc(2,3)=cmplx(Ae*Bx,Ae*By)
Hnuc(2,4)=((3.0**0.5)/6.0)*eta*DEq/((1.0+((eta**2.0)/3.0))**0.5)
Hnuc(3,3)=100.0+d+Ae*0.5*Bz-DEq/(2.0*((1.0+((eta**2.0)/3.0))**0.5))
Hnuc(3,4)=cmplx(0.866*Ae*Bx,0.866*Ae*By)
Hnuc(4,4)=100.0+d+Ae*1.5*Bz+DEq/(2.0*((1.0+((eta**2.0)/3.0))**0.5))

! Note: the number 100 is introduced to separate the I=1/2 and I=3/2 eigenvalues for convenient mathematical treatment.
!  It is artificial and will be removed in the subroutine intensitybuilder. 
do i=1,6		!LOWER TRIANGLE
	do j=1,6
		if (j>i) then
			Hnuc(j,i)=conjg(Hnuc(i,j))
		endif
	enddo
enddo


call zheev (JOBZ,UPLO,6,Hnuc,6,Nuc_eigenval,WORK,LWORK,RWORK,INFO) 	! DIAGONALIZES THE HAMILTONIAN AND RETURNS EIGENVECTORS ND EIGENVALUES

end subroutine





complex*8 function NORM(a,b)
integer :: a,b
!!! RETURNS THE NORMALIZATION FACTOR ASSOCIATED WITH EACH {S2,MS} EIGENSTATE, WHERE EACH SUCH EIGENSTATE IS A SUM OF NORMALIZED SPINORBITAL CONFIGURATIONS, NORMALIZED BY 1/N
!!!! a is the Ms. b is the multiplicity, i.e. 1,2, 3,4,5 or 6. 
NORM=1
if (b==2) then
	if (a==1) then
		NORM=1	
	elseif (a==2) then
		NORM=1
	endif
elseif (b==3) then
	if (a==1) then
		NORM=1
	elseif (a==2) then
		NORM=cmplx(2.0**0.5,0.0)
	elseif (a==3) then
		NORM=1
	endif
elseif (b==4) then
	if (a==1) then
		NORM=1.0
	elseif (a==2) then
		NORM=cmplx(3.0**0.5,0.0)
	elseif (a==3) then
		NORM=cmplx(3.0**0.5,0.0)
	elseif (a==4) then
		NORM=1.0
	endif
elseif (b==5) then
	if (a==1) then
		NORM=1.0
	elseif (a==2) then
		NORM=2.0
	elseif (a==3) then
		NORM=cmplx(6.0**0.5,0.0)
	elseif (a==4) then
		NORM=2.0
	elseif (a==5) then
		NORM=1.0
	endif
elseif (b==6) then
	if (a==1) then
		NORM=1.0
	elseif (a==2) then
		NORM=cmplx(5.0**0.5,0.0)
	elseif (a==3) then
		NORM=cmplx(10.0**0.5,0.0)
	elseif (a==4) then
		NORM=cmplx(10.0**0.5,0.0)
	elseif (a==5) then
		NORM=cmplx(5.0**0.5,0.0)
	elseif (a==6) then
		NORM=1.0
	endif

endif


end function



complex*8 function l(i,j,a)
integer :: i,j,k,ll
character(len=1) :: a
complex*8, dimension (1:5,1:5) :: lx,ly,lz
! RETURNS THE VALUE OF THE L MATRIX ELEMENTS IN THE BASIS OF THE REAL D-ORBITALS
!	I (INPUT, INTEGER) : ORBITAL IN THE BRA (DXY,DXZ,DYZ,DZ2,DX2-Y2 IN THAT ORDER)
!	J (INPUT, INTEGER) : ORBITAL IN THE KET (DXY,DXZ,DYZ,DZ2,DX2-Y2 IN THAT ORDER)
!	A (INPUT, CHARACTER): X,Y OR Z COMPONENT
!!!!!!!!!!Entering here matrix elements !!!!! 
lx(1,1)=0.0
lx(1,2)=cmplx(0.0,-1.0)
lx(1,3)=0.0
lx(1,4)=0.0
lx(1,5)=0.0
lx(2,2)=0.0
lx(2,3)=0.0
lx(2,4)=0.0
lx(2,5)=0.0
lx(3,3)=0.0
lx(3,4)=cmplx(0.0,-(3.0**0.5))
lx(3,5)=cmplx(0.0,-1.0)
lx(4,4)=0.0
lx(4,5)=0.0
lx(5,5)=0.0

ly(1,1)=0.0
ly(1,2)=0.0
ly(1,3)=cmplx(0.0,1.0)
ly(1,4)=0.0
ly(1,5)=0.0
ly(2,2)=0.0
ly(2,3)=0.0
ly(2,4)=cmplx(0.0,3.0**0.5)
ly(2,5)=cmplx(0.0,-1.0)
ly(3,3)=0.0
ly(3,4)=0.0
ly(3,5)=0.0
ly(4,4)=0.0
ly(4,5)=0.0
ly(5,5)=0.0

lz(1,1)=0.0
lz(1,2)=0.0
lz(1,3)=0.0
lz(1,4)=0.0
lz(1,5)=cmplx(0.0,2.0)
lz(2,2)=0.0
lz(2,3)=cmplx(0.0,-1.0)
lz(2,4)=0.0
lz(2,5)=0.0
lz(3,3)=0.0
lz(3,4)=0.0
lz(3,5)=0.0
lz(4,4)=0.0
lz(4,5)=0.0
lz(5,5)=0.0
do k=1,5
	do ll=1,5
		if (ll<k) then
			lx(k,ll)=conjg(lx(ll,k))
			ly(k,ll)=conjg(ly(ll,k))
			lz(k,ll)=conjg(lz(ll,k))
		endif
	enddo
enddo
if (a=='x') then
	l=lx(i,j)
elseif (a=='y') then
	l=ly(i,j)
elseif (a=='z') then
	l=lz(i,j)
endif


end function

complex*8 function vv(i,j,a,b)

! RETURNS THE VALUE OF THE ELECTRIC FIELD GRADIENT TENSOR COMPONENT IN THE BASIS OF THE REAL D-ORB.
!  I (INPUT, INTEGER) : ORBITAL IN THE BRA (DXY,DXZ,DYZ,DZ2,DX2-Y2 IN THAT ORDER)
!       J (INPUT, INTEGER) : ORBITAL IN THE KET (DXY,DXZ,DYZ,DZ2,DX2-Y2 IN THAT ORDER)
!	A,B (CHARACTER, INPUT): COMPONENT OF THE VAB (X,Y OR Z FOR A AND B) 

integer :: i,j,k,ll
character (len=1) :: a,b
complex*8, dimension (1:5,1:5) :: vxx,vxy,vxz,vyy,vyz,vzz

do k=1,5
	do ll=1,5
		vxx(k,ll)=0.0
		vxy(k,ll)=0.0
		vxz(k,ll)=0.0
		vyy(k,ll)=0.0
		vyz(k,ll)=0.0
		vzz(k,ll)=0.0
	enddo
enddo
vxx(1,1)=-2.0/7.0
vxx(2,2)=-2.0/7.0
vxx(3,3)=4.0/7.0
vxx(4,4)=2.0/7.0
vxx(5,5)=-2.0/7.0
vxx(4,5)=(2.0/7.0)*(3.0**0.5)
vxx(5,4)=(2.0/7.0)*(3.0**0.5)

vxy(1,4)=(2.0/7.0)*(3**0.5)
vxy(4,1)=(2.0/7.0)*(3**0.5)
vxy(2,3)=-3.0/7.0
vxy(3,2)=-3.0/7.0

vxz(1,3)=-3.0/7.0
vxz(3,1)=-3.0/7.0
vxz(2,5)=-3.0/7.0
vxz(5,2)=-3.0/7.0
vxz(2,4)=-(3.0**0.5)/7.0
vxz(4,2)=-(3.0**0.5)/7.0

vyy(1,1)=-2.0/7.0
vyy(2,2)=4.0/7.0
vyy(3,3)=-2.0/7.0
vyy(4,4)=2.0/7.0
vyy(5,5)=-2.0/7.0
vyy(4,5)=-(2.0/7.0)*(3.0**0.5)
vyy(5,4)=-(2.0/7.0)*(3.0**0.5)

vyz(1,2)=-3.0/7.0
vyz(2,1)=-3.0/7.0
vyz(3,4)=-(3.0**0.5)/7.0
vyz(4,3)=-(3.0**0.5)/7.0
vyz(3,5)=3.0/7.0
vyz(5,3)=3.0/7.0

vzz(1,1)=4.0/7.0
vzz(2,2)=-2.0/7.0
vzz(3,3)=-2.0/7.0
vzz(4,4)=-4.0/7.0
vzz(5,5)=4.0/7.0

if (a=='x') then
	if (b=='x') then
		vv=vxx(i,j)
	elseif (b=='y') then
		vv=vxy(i,j)
	elseif (b=='z') then
		vv=vxz(i,j)
	endif
elseif (a=='y') then
	if (b=='x') then
		vv=vxy(i,j)
	elseif (b=='y') then
		vv=vyy(i,j)
	elseif (b=='z') then
		vv=vyz(i,j)
	endif
elseif (a=='z') then
	if (b=='x') then
		vv=vxz(i,j)
	elseif (b=='y') then
		vv=vyz(i,j)
	elseif (b=='z') then
		vv=vzz(i,j)
	endif
endif



end function


complex*8 function s(i,j,a)
integer :: i,j,k,l
character(len=1) :: a
complex*8, dimension (0:1,0:1) :: sx,sy,sz

! RETURNS THE VALUE OF THE SPIN VECTOR COMPONENTS IN THE BASIS OF THE 1/2,+1/2; 1/2,-1/2 SPIN EIGENSTATES
! I (INPUT, INTEGER): SPIN-UP (0) OR SPIN-DOWN (1) IN THE BRA
! J (INPUT, INTEGER): SPIN-UP (0) OR SPIN-DOWN (1) IN THE KET
! A (INPUT, CHARACTER): COMPONENT OF THE SPIN VECTOR (X,Y OR Z)

sx(0,0)=0.0
sx(0,1)=0.5
sx(1,1)=0.0
sy(0,0)=0.0
sy(0,1)=cmplx(0.0,-0.5)
sy(1,1)=0.0
sz(0,0)=0.5
sz(0,1)=0.0
sz(1,1)=-0.5

do k=0,1
        do l=0,1
                if (l<k) then
                        sx(k,l)=conjg(sx(l,k))
                        sy(k,l)=conjg(sy(l,k))
                        sz(k,l)=conjg(sz(l,k))
                endif
        enddo
enddo

if (a=='x') then
        s=sx(i,j)
elseif (a=='y') then
        s=sy(i,j)
elseif (a=='z') then
        s=sz(i,j)
endif

end function

end module



program call			
! READS THE INPUT FILE AND LAUNCHES THE APPROPRIATE PROCEDURES. IN ORDER, BUILDS THE {MS,S2} EIGENSTATES FROM THE USER-DEFINED ELECTRONIC CONFIGURATION (S2 IS ALWAYS MAXIMUMFOR A GIVEN CONFIGURATION)
! THEN USES the eigenstates to BUILD THE APPROPRIATE SOC, ORBITAL MOMENTA, SPIN MOMENTA AND SPIN DIPOLE MATRICES NECESSARY FOR BUILDING THE FIELD-DEPENDENT HAMILTONIAN IN THE BASIS OF THESE EIGENSTATES.
! THEN, EITHER CALCULATES THE MAGNETIZATION/SUSCEPTIBILITY/CHI*T OR THE FE MOSSBAUER SPECTRUM ACCORDING TO A SET OF USER DEFINED CONSTRAINTS.  
! THE END OF THE PROGRAM GIVES INFORMATIONS ON THE ELECTRONIC STRUCTURE SUCH AS PRINTING THE SOC MATRIX, THE EIGENVALUES, THE EIGENSTATES IN THE BASIS OF THE {MS,S2} SPIN EIGENSTATES OF THE CHOSEN BASIS.
use procedures
implicit none
integer :: nstate,i,j,k,multiplicity,nconfig,maxdim,Ntemp,extend
integer, allocatable, dimension (:,:) :: orb_pop
integer, allocatable, dimension(:,:,:,:) :: vector_matrix
integer, allocatable, dimension (:) :: ndummy
complex*16, allocatable, dimension (:,:) :: HAM,Lx2,Ly2,Lz2,Sx2,Sy2,Sz2,SOC2,SDx2,SDy2,SDz2,result_micro
real*8, dimension (1:8) :: lorentz
real*8, allocatable, dimension (:) :: eigenval
real*16, allocatable, dimension (:) :: DE,T
real*16 :: lambda,B,Bdummy, DEq,eta,d,lorentz2,AFC,rd,ex,ey,ez
real*8 :: Imax,resolution,bound,extendreal
character (len=3) :: job_name
character (len=100) :: temp_input,gridlevel, limit,info


read(*,*) info
info=trim(info)
read(*,*) nstate
allocate (orb_pop(1:nstate,1:5))
allocate (ndummy(1:nstate))

do i=1,nstate
	read(*,*) orb_pop(i,1), orb_pop(i,2), orb_pop(i,3), orb_pop(i,4), orb_pop(i,5)
enddo


call vector_builder(orb_pop,nstate,vector_matrix,ndummy,nconfig,multiplicity)	!! Calculates the MS,S2 eigenstates from the configurations stored in orb_pop
call matrix_builder(vector_matrix,multiplicity,ndummy,nstate,nconfig,maxdim,Lx2,Ly2,Lz2,Sx2,Sy2,Sz2,SOC2,SDx2,SDy2,SDz2)	!! Uses the Ms, S2 eigenstates to calculate the appropriate matrices

!!!! Now storing the energy differences from the input: CAREFUL! FIRST STATE ENERGY FIXED TO 0 
allocate (DE(1:nstate))

do i=2,nstate
	read(*,*) DE(i)
enddo
DE(1)=0.0

!!! FINALLY getting the SOC constant for the metal center
read (*,*) lambda

!!!!!! DONE. Now, let us recognize the type of job. 'Mos' for Mossbauer and "mag" for magnetic susceptibility measurement
read(*,*) job_name
job_name=trim(job_name)
if (job_name=='mag'.or.job_name=='X'.or.job_name=='XT') then    	!!! Magnetization, magnetic susceptibility or chi*T
	read(*,*) B
	read(*,*) temp_input						
	temp_input=trim(temp_input)
	open(unit=1,file=temp_input)
	read(1,*) Ntemp
	allocate (T(1:Ntemp))
	do i=1,Ntemp
		read(1,*) T(i)
	enddo
	close(unit=1)
	read(*,*) gridlevel
	gridlevel=trim(gridlevel)
	call susc_calc(Ntemp,T,gridlevel,job_name,DE,SOC2,Lx2,Ly2,Lz2,Sx2,Sy2,Sz2,maxdim,ndummy,lambda,B)	! calculates and prints the magnetization or associated functions for each temp.
elseif (job_name=='mos') then		! Mossbauer section. Input B, DEq, d, T(1), limit 
	allocate (T(1:1))
	read(*,*) B			! (quadrupole splitting,asymmetry parameter, isomer shift, temperature, fast or slow magnetic relaxation limit)
	read(*,*) DEq,eta, d
	read(*,*) T(1)
	read(*,*) limit 		! limit=fast or limit=slow
	limit=trim(limit)
	read(*,*) gridlevel
	gridlevel=trim(gridlevel)
	read(*,*) Imax
	read(*,*) lorentz2
	do i=1,8
		lorentz(i)=lorentz2
	enddo
	read(*,*) rd
	read(*,*) AFC
	read(*,*) bound
	read(*,*) resolution
	extendreal=bound/resolution
	extend=int(extendreal)
	if (real(extend)<extendreal) then
		extend=extend+1
	endif
	call mos_calc(B,DEq,eta,d,T(1),limit,DE,SOC2,Lx2,Ly2,Lz2,Sx2,Sy2,Sz2,SDx2,SDy2,SDz2,maxdim,ndummy,lambda,gridlevel,&
&Imax,lorentz,rd,AFC,extend,resolution)			! calculates and prints out the mossbauer spectrum associated with the user-specified set of parameters 


else
	print*, 'job type unavailable'
endif
if (info=='yes') then						! prints a bunch of relevant values to understand the magnetism at a microscopic level
	print*, '**** SPIN ORBIT COUPLING MATRIX*****'
	do i=1,maxdim
		do j=1,maxdim
			print*, 'SOC(',i,',',j,')',SOC2(i,j)
		enddo
	enddo
	Bdummy=0
	allocate(result_micro(1:maxdim,1:9))
	call diag_Ham(DE,SOC2,Lx2,Ly2,Lz2,Sx2,Sy2,Sz2,ex,ey,ez,maxdim,ndummy,lambda,Bdummy,HAM,eigenval)
       ! call mag_micro(HAM,maxdim,ndummy,Lx2,Ly2,Lz2,Sx2,Sy2,Sz2,SDx2,SDy2,SDz2,result_micro)	
	print*, '****HAMILTONIAN  EIGENVALUES *******'
	do i=1,maxdim
		print*,i, eigenval(i)
	enddo	
	print*, '****HAMILTONIAN EIGENSTATES  *******'
	print*, 'now printing the composition of the eigenstates in the basis of the |S,Ms> eigenstates of each chosen &
&configuration (descending order of Ms)'
	do i=1,maxdim
		print*, 'state', i 
		do j=1,maxdim
			print*, HAM(j,i)
		enddo
	enddo
	print*, '**** MAGNETIZATION ALONG Z  *******'	
	print*, 'spin and orbital momenta for each states under a field of',B,'T along z'
	deallocate (HAM,eigenval)
	ez=1.0
	ex=0.0
	ey=0.0
	call diag_Ham(DE,SOC2,Lx2,Ly2,Lz2,Sx2,Sy2,Sz2,ex,ey,ez,maxdim,ndummy,lambda,B,HAM,eigenval)
	call mag_micro(HAM,maxdim,ndummy,Lx2,Ly2,Lz2,Sx2,Sy2,Sz2,SDx2,SDy2,SDz2,result_micro)
	do i=1,maxdim
		print*, 'state',i
		print*, '<Lz>=',real(result_micro(i,3)), '<Sz>=',real(result_micro(i,6))
	enddo 	
	print*, '**** MAGNETIZATION ALONG x  *******'
        print*, 'spin and orbital momenta for each states under a field of',B,'T along x'
        deallocate (HAM,eigenval)
        ez=0.0
        ex=1.0
        ey=0.0
        call diag_Ham(DE,SOC2,Lx2,Ly2,Lz2,Sx2,Sy2,Sz2,ex,ey,ez,maxdim,ndummy,lambda,B,HAM,eigenval)
        call mag_micro(HAM,maxdim,ndummy,Lx2,Ly2,Lz2,Sx2,Sy2,Sz2,SDx2,SDy2,SDz2,result_micro)
        do i=1,maxdim
                print*, 'state',i
                print*, '<Lx>=',real(result_micro(i,1)), '<Sx>=',real(result_micro(i,4))
        enddo
	print*, '**** MAGNETIZATION ALONG y  *******'
        print*, 'spin and orbital momenta for each states under a field of',B,'T along y'
        deallocate (HAM,eigenval)
        ez=0.0
        ex=0.0
        ey=1.0
        call diag_Ham(DE,SOC2,Lx2,Ly2,Lz2,Sx2,Sy2,Sz2,ex,ey,ez,maxdim,ndummy,lambda,B,HAM,eigenval)
        call mag_micro(HAM,maxdim,ndummy,Lx2,Ly2,Lz2,Sx2,Sy2,Sz2,SDx2,SDy2,SDz2,result_micro)
        do i=1,maxdim
                print*, 'state',i
                print*, '<Ly>=',real(result_micro(i,2)), '<Sy>=',real(result_micro(i,5))
        enddo

endif


end program

