!****************************************************************************!
!* Project           : Kinematic fitting
!*
!* Program name      : Fast UNiversal Kinematic FITting (FUNK_FIT.F90)
!*
!* Author            : Robin Smith
!*
!* Institution       : University of Birmingham
!*
!* Date created      : 22/05/2017
!*
!* Purpose           : Provide a general kinematic fitting platform that can
!*                     easily be adapted by code users for their purpose
!*
!* Revision History  :
!*
!* Date        Author      Ref    Revision (Date in DD/MM/YYYY format) 
!* N/A           N/A       N/A                  N/A
!*
!* Compile by typing : gfortran FUNK_FIT.F90
!*
!****************************************************************************!
!******************************** MODULES ***********************************!

MODULE CONSTANTS ! Some constants useful in nuclear physics
                 ! Some arrays used in the kinematic fitting code

 Double precision, parameter, public 			:: pi=3.14159274
 Double precision, parameter, public 			:: hbar=197.326 !MeV c
 Double precision, parameter, public 			:: me=0.511 !MeV / c^2
 Double precision, parameter, public 			:: mn=939.565 !MeV / c^2
 Double precision, parameter, public 			:: mp=938.272 !MeV / c^2
 Double precision, parameter, public 			:: Q12C3A=-7.275 !MeV
 Double precision, parameter, public 			:: ma=4.0 !MeV
 Double precision, parameter, public 			:: m12c=12.0 !MeV
 Double precision, parameter, public 			:: m8be=8.0 !MeV
 Double precision, dimension(:), allocatable, public 	:: parameters,errors,parametersnew,parametersvar
 Double precision, dimension(:), allocatable, public 	:: origparameters,deltaparameters
 integer, public					:: Nparam,Nconst

END MODULE CONSTANTS

!****************************************************************************!
!***************************** START PROGRAM ********************************!

PROGRAM BASICFITTING
 USE CONSTANTS
 IMPLICIT NONE

!------------------------------ FUNCTIONS -----------------------------------!
!-------------ANY NEW CONSTRAINT EQUATIONS SHOULD BE ADDED HERE--------------!
 Double precision, external 			:: SimpleConstraints

!---------------------------- SYSTEM CONSTANTS ------------------------------!
!---------------VARIABLES USED IN THE KINEMATIC FITTING CODE-----------------!
 integer :: i,j,k,n,Nevent,event,DeAllocateStatus,EventNumber,iter,maxiter,dum
 Double precision,  allocatable, dimension(:) 	:: d1,Lambda,Mult2
 Double precision,  allocatable, dimension(:,:) :: Va,D,DT,Mult,VDinv1,VDinv,VD,Vanew
 Double precision,  allocatable, dimension(:,:) :: C1,C2,C3,C4
 Double precision				:: deltaparam
 Double precision				:: soften

 Nevent=10000 ! Total number of events in the input file

 ! Open files for input data and output fitted data
 open(unit=11,file='/Users/rxs883/Documents/KinematicFitting/FUNK_FIT/TestData.txt')
 open(unit=12,file='/Users/rxs883/Documents/KinematicFitting/FUNK_FIT/TestDataFit.txt')

 ! Loop over all events in the file
 do event=1,Nevent

  Nparam = 3 ! Number of parameters in each event
  Nconst = 2 ! Number of constraint equations

  ! Allocate matrices for calculations in memory
  allocate ( origparameters(Nparam) )
  allocate ( deltaparameters(Nparam) )
  allocate ( parameters(Nparam) )
  allocate ( parametersnew(Nparam) )
  allocate ( parametersvar(Nparam) )
  allocate ( errors(Nparam) )
  allocate ( D(Nconst,Nparam) )
  allocate ( d1(Nconst) )
  allocate ( Va(Nparam,Nparam) )

  ! Zero the arrays
  do i=1,Nconst
   d1(i)=0.0
  end do

  do i=1,Nparam
   parameters(i)=0.0
   parametersnew(i)=0.0
   errors(i)=0.0
   do k=1,Nconst
    D(k,i)=0.0
   end do
   do j=1,Nparam
    Va(i,j)=0.0
   end do
  end do

  ! Set original parameters and errors from the input file
  do i=1,Nparam
   read(11,*)origparameters(i),errors(i)
   parameters(i)=origparameters(i)
   errors(i)=errors(i)
  end do

  ! Set up Va matrix
  do i=1,Nparam
   Va(i,i)=errors(i)**2
  end do

  maxiter=1000 ! Number of iterations required to achieve convergence
  do iter=1,maxiter

   ! Set up d matrix
   do i=1,Nconst
    d1(i)=SimpleConstraints(i,1,0)
   end do

   ! Set up D matrix numerically
   do i=1,Nconst
    do j=1,Nparam
     deltaparam=2.0*abs(errors(j)/1000.0)
     if (errors(j).ne.0) then
      D(i,j)=(SimpleConstraints(i,j,1) - SimpleConstraints(i,j,-1))/deltaparam
     else
      D(i,j)=(SimpleConstraints(i,j,1) - SimpleConstraints(i,j,-1))
     end if
    end do
   end do
 
   DT=transpose(D)

   ! Calculate VDinv matrix
   VDinv1=matmul(Va,DT)
   VDinv=matmul(D,VDinv1)
   ! Calculate VD matrix
   allocate ( VD(size(VDinv,DIM=1),size(VDinv,DIM=2)) )
   call inverse(VDinv,VD,size(VDinv,DIM=1))

   ! Calculate Lambda Lagrange multiplier matrix
   deltaparameters=parameters-origparameters
   Lambda=matmul(D,deltaparameters)+matmul(VD,d1)

   ! Now put all of these together to calculate the new parameters
   Mult=matmul(Va,DT)
   Mult2=matmul(Mult,Lambda)
   soften=0.05 ! variable to slow down the convergence
   ! soften = 1
   parameters=origparameters-(soften*Mult2)

   origparameters=parameters

   ! Calculate the new covariance matrix
   C1=matmul(D,Va)
   C2=matmul(DT,VD)
   C3=matmul(C2,C1)
   C4=matmul(Va,C3)
   Vanew = Va - C4

   deallocate (VD, STAT = DeAllocateStatus)

  end do
 
  ! Write the parameters after kinematic fitting to file, along with their reduced uncertainties
  write(12,*)parameters(1),Vanew(1,1)
  write(12,*)parameters(2),Vanew(2,2)
  write(12,*)parameters(3),Vanew(3,3)

  ! **REMEMBER TO DEALLOCATE ALL ARRAYS BEFORE READING IN A NEW EVENT**
  deallocate (origparameters, STAT = DeAllocateStatus)
  deallocate (deltaparameters, STAT = DeAllocateStatus)
  deallocate (parameters, STAT = DeAllocateStatus)
  deallocate (parametersnew, STAT = DeAllocateStatus)
  deallocate (parametersvar, STAT = DeAllocateStatus)
  deallocate (errors, STAT = DeAllocateStatus)
  deallocate (D, STAT = DeAllocateStatus)
  deallocate (d1, STAT = DeAllocateStatus)
  deallocate (Va, STAT = DeAllocateStatus)
  deallocate (Vanew, STAT = DeAllocateStatus)
  deallocate (VD, STAT = DeAllocateStatus)

 end do

 ! Close files after all events have been read
 close(unit=11)
 close(unit=12) 

END PROGRAM BASICFITTING


!-------------------------------------------------------------------------

 FUNCTION SimpleConstraints (ConstNum,ParamNum,var) RESULT (r)
 USE CONSTANTS

 ! Function evaluates the value of the constraint equations
 ! ConstNum denotes which constraint equation to evaluate
 ! ParamNum signifies which parameter to vary if 
 ! var says whether to vary the parameter

 integer ConstNum,ParamNum,var
 Double precision constant1,constant2,r
 Double precision VarParameters(Nparam)
 Double precision e(3)

 exhoyle=7.6542 ! MeV

 ! First we have to vary the parameters if necessary for numerical differentiation
 VarParameters=parameters ! Not varied
 ! If parameter var = +1 of -1 then 
 if (var.eq.1) then
  VarParameters(ParamNum)=VarParameters(ParamNum)+(abs(errors(ParamNum))/1000.0)
 end if
 if (var.eq.-1) then
  VarParameters(ParamNum)=VarParameters(ParamNum)-(abs(errors(ParamNum))/1000.0)
 end if

 ! Now link the parameters with their physical quantities
 do i=1,3
  e(i)=VarParameters( i )
 end do

 if ( ConstNum.eq.1 ) then ! Constraint Eqn 1
  r = e(1) + e(2) + e(3) - 29
 end if

 if ( ConstNum.eq.2 ) then ! Constraint Eqn 2
  r = e(1)**2 + e(2)**2 + e(3)**2 - 353
 end if

 if ( ConstNum.eq.3 ) then ! Constraint Eqn 3
  r = sqrt(e(1)) + sqrt(e(2)) + sqrt(e(3)) - 9
 end if

 END FUNCTION SimpleConstraints

!-------------------------------------------------------------------------
  subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
 
integer n
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k
! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0
! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do
! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do
! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse
!-------------------------------------------------------------------------
