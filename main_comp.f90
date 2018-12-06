
      module c_module
        INTERFACE
          subroutine INTERPOlATE() BIND (C, NAME='interpolate')
            USE ISO_C_BINDING
            implicit none
          end subroutine INTERPOLATE
        END INTERFACE
      end module c_module


!--------------------------------------------------------------------------------
!
!     Compressible channel code (x is direction normal to the wall)
!
!--------------------------------------------------------------------------------

  program inital
  use decomp_2d
  implicit none
  include 'param.txt'
  include 'mpif.h'

  real*8, dimension(iold,jold/p_row,kold/p_col)             :: rold,uold,vold,wold,eold
  real*8, dimension(inew,jnew/p_row,knew/p_col)             :: rnew,unew,vnew,wnew,enew
  real*8, dimension(0:iold+1,0:jold/p_row+1,0:kold/p_col+1) :: rbo,ubo,vbo,wbo,ebo,tmpo
  real*8, dimension(0:inew+1,0:jnew/p_row+1,0:knew/p_col+1) :: rbn,ubn,vbn,wbn,ebn,tmpn
  real*8, dimension(iold)                                   :: xold 
  real*8, dimension(kold)                                   :: zold 
  real*8, dimension(inew)                                   :: xnew 
  real*8, dimension(knew)                                   :: znew 
  
  integer ierr,sstep,i,j,k

  character*6 cha

  sstep = 2i275
 
  call grid_init(3.0,xold,zold,iold,jold,kold) 
  call grid_init(3.0,xnew,znew,inew,jnew,knew) 
  call mpi_init(ierr)
  call decomp_2d_init(iold,jold,kold,p_row,p_col)
  if(nrank.eq.0) write(*,*) "Reading old files" 
  call loadRestart(sstep,rold,uold,vold,wold,eold,iold,jold,kold,p_row,p_col)
  call output(nrank,rold,uold,vold,wold,eold,iold,jold,kold,xold,p_row,p_col,0)
  call output2D(nrank,rold,uold,vold,wold,eold,iold,jold,kold,xold,zold,p_row,p_col,0)
  rbo= 0
  ubo= 0
  vbo= 0
  wbo= 0
  ebo= 0
  rbn= 0
  ubn= 0
  vbn= 0
  wbn= 0
  ebn= 0
  do i=1,iold
   do j=1,jold/p_row
    do k=1,kold/p_col
     rbo(i,j,k) = rold(i,j,k)
     ubo(i,j,k) = uold(i,j,k)
     vbo(i,j,k) = vold(i,j,k)
     wbo(i,j,k) = wold(i,j,k)
     ebo(i,j,k) = eold(i,j,k)
    enddo
   enddo
  enddo
  write(*,*) "interpolating density";  call INTERPOLATE(rbo, rbn); !, iold, jold, kold, inew, jnew, knew);
  write(*,*) "interpolating u-vel";    call INTERPOLATE(ubo, ubn); !, iold, jold, kold, inew, jnew, knew);
  write(*,*) "interpolating v-vel";    call INTERPOLATE(vbo, vbn); !, iold, jold, kold, inew, jnew, knew);
  write(*,*) "interpolating w-vel";    call INTERPOLATE(wbo, wbn); !, iold, jold, kold, inew, jnew, knew);
  write(*,*) "interpolating energy";   call INTERPOLATE(ebo, ebn); !, iold, jold, kold, inew, jnew, knew);
  write(*,*) "Finished interpolating, sending values"
 
  rnew=rbn(1:inew,1:jnew/p_row,1:knew/p_col)
  unew=ubn(1:inew,1:jnew/p_row,1:knew/p_col)
  vnew=vbn(1:inew,1:jnew/p_row,1:knew/p_col)
  wnew=wbn(1:inew,1:jnew/p_row,1:knew/p_col)
  enew=ebn(1:inew,1:jnew/p_row,1:knew/p_col)
  call output(nrank,rnew,unew,vnew,wnew,enew,inew,jnew,knew,xnew,p_row,p_col,1)
  call output2D(nrank,rnew,unew,vnew,wnew,enew,inew,jnew,knew,xnew,znew,p_row,p_col,1)
  call decomp_2d_finalize

  
  call decomp_2d_init(inew,jnew,knew,p_row,p_col)
  if(nrank.eq.0) write(*,*) "Writing down new files" 
  call saveRestart(sstep,rnew,unew,vnew,wnew,enew,inew,jnew,knew,p_row,p_col)
  call decomp_2d_finalize
  call mpi_finalize(ierr)

end program

subroutine loadRestart(istep,rnew,unew,vnew,wnew,enew,iold,jold,kold,p_row,p_col)
  use decomp_2d
  use decomp_2d_io
  implicit none
  integer iold,jold,kold,p_row,p_col,i,j,k
  real*8 rnew(iold,jold/p_row,kold/p_col)
  real*8 enew(iold,jold/p_row,kold/p_col)
  real*8 unew(iold,jold/p_row,kold/p_col)
  real*8 vnew(iold,jold/p_row,kold/p_col)
  real*8 wnew(iold,jold/p_row,kold/p_col)
  real*8 a(iold,2),rho_0,U_0
  integer istep
  character*6 cha
  write(cha,'(i5.5)') istep
  call decomp_2d_read_one(1,unew,'data-old/u.'//cha//'.dns')
  call decomp_2d_read_one(1,vnew,'data-old/v.'//cha//'.dns')
  call decomp_2d_read_one(1,wnew,'data-old/w.'//cha//'.dns')
  call decomp_2d_read_one(1,rnew,'data-old/c.'//cha//'.dns')
  call decomp_2d_read_one(1,enew,'data-old/e.'//cha//'.dns')
end subroutine loadRestart


subroutine saveRestart(istep,rnew,unew,vnew,wnew,enew,iold,jold,kold,p_row,p_col)
  use decomp_2d
  use decomp_2d_io
  implicit none
  integer iold,jold,kold,p_row,p_col
  real*8 rnew(iold,jold/p_row,kold/p_col)
  real*8 enew(iold,jold/p_row,kold/p_col)
  real*8 unew(iold,jold/p_row,kold/p_col)
  real*8 vnew(iold,jold/p_row,kold/p_col)
  real*8 wnew(iold,jold/p_row,kold/p_col)
  integer istep
  character*6 cha
  write(cha,'(I5.5)') istep
  call decomp_2d_write_one(1,unew,'data-new/u.'//cha//'.dns')
  call decomp_2d_write_one(1,vnew,'data-new/v.'//cha//'.dns')
  call decomp_2d_write_one(1,wnew,'data-new/w.'//cha//'.dns')
  call decomp_2d_write_one(1,rnew,'data-new/c.'//cha//'.dns')
  call decomp_2d_write_one(1,enew,'data-new/e.'//cha//'.dns')
end subroutine saveRestart


subroutine output(rank,rnew,unew,vnew,wnew,enew,iold,jold,kold,x,p_row,p_col,flag)
  use decomp_2d
  implicit none
  include 'mpif.h'
  integer rank,ierr,i,j,k,p_row,p_col,iold,jold,kold,flag
  character*6 cha
  character*6 cha2

  real*8 Re

  real*8, dimension(iold,jold/p_row,kold/p_col) :: rnew,unew,vnew,wnew,enew

  real*8, dimension(iold) :: rm,um,vm,wm,tm,pm,em,    &
                             rr,ur,vr,wr,tr,pr,uw,    &
                             ufm,vfm,wfm,tfm,         &
                             ufr,vfr,wfr,tfr,         &
                             tmpiold,stress,x

  rm=0
  um=0
  vm=0
  wm=0
  tm=0
  pm=0
  do i=1,iold
    rm(i) = sum(rnew(i,:,:))
    um(i) = sum(unew(i,:,:))
    vm(i) = sum(vnew(i,:,:))
    wm(i) = sum(wnew(i,:,:))
    em(i) = sum(enew(i,:,:))
  enddo
  call mpi_allreduce(rm,tmpiold,iold,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr);   rm = tmpiold/(jold*kold)
  call mpi_allreduce(um,tmpiold,iold,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr);   um = tmpiold/(jold*kold)
  call mpi_allreduce(vm,tmpiold,iold,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr);   vm = tmpiold/(jold*kold)
  call mpi_allreduce(wm,tmpiold,iold,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr);   wm = tmpiold/(jold*kold)
  call mpi_allreduce(em,tmpiold,iold,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr);   em = tmpiold/(jold*kold)


  rr=0
  ur=0
  vr=0
  wr=0
  uw = 0
  do k=1,kold/p_col
     do j=1,jold/p_row
        do i=1,iold
           rr(i) = rr(i) + (rnew(i,j,k)-rm(i))**2
           ur(i) = ur(i) + (unew(i,j,k)-um(i))**2
           vr(i) = vr(i) + (vnew(i,j,k)-vm(i))**2
           wr(i) = wr(i) + (wnew(i,j,k)-wm(i))**2
           uw(i) = uw(i) + (wnew(i,j,k)-wm(i))*(unew(i,j,k)-um(i))
        enddo
     enddo
  enddo
  call mpi_allreduce(rr,tmpiold,iold,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr);  rr = tmpiold/(jold*kold)
  call mpi_allreduce(ur,tmpiold,iold,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr);  ur = tmpiold/(jold*kold)
  call mpi_allreduce(vr,tmpiold,iold,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr);  vr = tmpiold/(jold*kold)
  call mpi_allreduce(wr,tmpiold,iold,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr);  wr = tmpiold/(jold*kold)
  call mpi_allreduce(uw,tmpiold,iold,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr);  uw = tmpiold/(jold*kold)


  ufm=0
  vfm=0
  wfm=0
  tfm=0
  do k=1,kold/p_col
     do j=1,jold/p_row
        do i=1,iold
           ufm(i) = ufm(i) + rnew(i,j,k)*unew(i,j,k)/(jold*kold)
           vfm(i) = vfm(i) + rnew(i,j,k)*vnew(i,j,k)/(jold*kold)
           wfm(i) = wfm(i) + rnew(i,j,k)*wnew(i,j,k)/(jold*kold)
        enddo
     enddo
  enddo
  call mpi_allreduce(ufm,tmpiold,iold,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr);  ufm = tmpiold/rm/(jold*kold)
  call mpi_allreduce(vfm,tmpiold,iold,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr);  vfm = tmpiold/rm/(jold*kold)
  call mpi_allreduce(wfm,tmpiold,iold,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr);  wfm = tmpiold/rm/(jold*kold)


  ufr=0
  vfr=0
  wfr=0
  do k=1,kold/p_col
     do j=1,jold/p_row
        do i=1,iold
           ufr(i) = ufr(i) + (unew(i,j,k)-ufm(i))**2
           vfr(i) = vfr(i) + (vnew(i,j,k)-vfm(i))**2
           wfr(i) = wfr(i) + (wnew(i,j,k)-wfm(i))**2
        enddo
     enddo
  enddo
  call mpi_allreduce(ufr,tmpiold,iold,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr);   ufr = tmpiold/(jold*kold)
  call mpi_allreduce(vfr,tmpiold,iold,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr);   vfr = tmpiold/(jold*kold)
  call mpi_allreduce(wfr,tmpiold,iold,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr);   wfr = tmpiold/(jold*kold)


  if (rank.eq.0) then
     if(flag.eq.0)  open(17,file = 'old_prof')
     if(flag.eq.1)  open(17,file = 'new_prof')
     do i=1,iold
        write(17,'(18E16.8)') x(i),um(i),vm(i),wm(i),rm(i),em(i), &   ! 1,2,3,4,5,6,7,8
                              sqrt(ur(i)),    &  ! 11
                              sqrt(vr(i)),    &  ! 12
                              sqrt(wr(i)),    &  ! 13
                              uw(i),          &  ! 14
                              sqrt(rr(i))        ! 15
     enddo
     close(17)
  endif
end subroutine output


subroutine output2D(rank,rnew,unew,vnew,wnew,enew,iold,jold,kold,x,z,p_row,p_col,flag)
  use decomp_2d
  implicit none
  include 'mpif.h'
  integer rank,ierr,i,j,k,p_row,p_col,iold,jold,kold,flag
  character*6 cha
  character*6 cha2

  real*8 Re

  real*8, dimension(iold,jold/p_row,kold/p_col) :: rnew,unew,vnew,wnew,enew,dwdx

  real*8, dimension(iold,kold/p_col) :: rm,um,vm,wm,tm,pm,em,    &
                                        rr,ur,vr,wr,tr,pr,uw,    &
                                        ufm,vfm,wfm,tfm,         &
                                        ufr,vfr,wfr,tfr,         &
                                        tmpiold,stress
  real*8 x(iold),z(kold/p_col)

  rm=0
  um=0
  vm=0
  wm=0
  tm=0
  pm=0
  do i=1,iold
   do k=1,kold/p_col
    rm(i,k) = sum(rnew(i,:,k))
    um(i,k) = sum(unew(i,:,k))
    vm(i,k) = sum(vnew(i,:,k))
    wm(i,k) = sum(wnew(i,:,k))
    em(i,k) = sum(enew(i,:,k))
   enddo
  enddo
  call mpi_allreduce(rm,tmpiold,iold*kold/p_col,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr);   rm = tmpiold/(jold)
  call mpi_allreduce(um,tmpiold,iold*kold/p_col,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr);   um = tmpiold/(jold)
  call mpi_allreduce(vm,tmpiold,iold*kold/p_col,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr);   vm = tmpiold/(jold)
  call mpi_allreduce(wm,tmpiold,iold*kold/p_col,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr);   wm = tmpiold/(jold)
  call mpi_allreduce(em,tmpiold,iold*kold/p_col,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr);   em = tmpiold/(jold)

  if (rank.eq.0) then
     if(flag.eq.0)  open(17,file = 'old_prof2D')
     if(flag.eq.1)  open(17,file = 'new_prof2D')
     do i=1,iold
      do k=1,kold/p_col
        write(17,'(18E16.8)') x(i),16*atan(1.0)/kold*(k-0.5),um(i,k),vm(i,k),wm(i,k),rm(i,k),em(i,k)    ! 1,2,3,4,5,7
      enddo
     enddo
     close(17)
  endif
end subroutine output2D

subroutine grid_init(fact,xp,zp,imax,jmax,kmax)
  implicit none
  integer imax,jmax,kmax,i,rank,k
  real*8 x,fact,xu(0:imax),xp(imax),mr_c(0:imax),mr_s(imax),dmr_s(imax),zp(kmax)
  real*8 dx,dy,dz,Ly,Lz,pi

  pi = 4.*atan(1.)
  Ly = 2.*pi
  Lz = 4.*pi
  dx = 2./imax
  dy = Ly/jmax
  dz = Lz/kmax

  do i=0,imax
     x       = 1.*i/imax - 0.5
     xu(i)   = 1.0 + tanh(fact*x)/tanh(fact*0.5)
     mr_c(i) = 0.5*fact/tanh(fact/2.)/cosh(fact*x)**2.0
  enddo

  do i=1,imax
     x        =  1.*(i-0.5)/imax - 0.5
     xp(i)    =  1.0 + tanh(fact*x)/tanh(fact*0.5)
  enddo
 
  do k=1,kmax
     zp(k)    = (k-0.5)*dz
  enddo 

end subroutine grid_init

subroutine init_inflow(rank,rnew,unew,vnew,wnew,enew,imax,jmax,kmax,p_row,p_col)
  implicit none
  integer rank,imax,jmax,kmax,p_row,p_col,i,j,k
  real*8  temp
  real*8, dimension(imax,jmax/p_row,kmax/p_col) :: rnew,enew,unew,vnew,wnew


  rnew = 1.0
  unew = 0.0
  vnew = 0.0
  wnew = 0.0

  do i=1,imax
     wnew(i,:,:) = 1.5*(i)*(imax-(i))
     unew(i,:,:) = wnew(i,:,:)/(1.0*i)
     vnew(i,:,:) = wnew(i,:,:)/(1.0*i)+(i)
  enddo
  do i =1,imax
     temp = 1. 
     enew(i,:,:) = temp/0.4 + 0.5*(unew(i,:,:)**2 + vnew(i,:,:)**2 + wnew(i,:,:)**2)
  enddo

  if (rank.eq.0) then
    open(11,file = 'Init.txt')
    do i=1,imax
       write(11,'(I5,5E20.10)') i,rnew(i,1,1),unew(i,1,1),vnew(i,1,1),wnew(i,1,1),enew(i,1,1)
    enddo
    close(11)
  endif

end subroutine init_inflow

