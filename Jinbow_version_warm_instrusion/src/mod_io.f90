subroutine save3d(Nx,Ny,Nz,var,filename)
   integer :: Nx,Ny,Nz
   character(len=*) :: filename
   real*8 :: var(Nx,Ny,Nz)
   open(3333,file=trim(filename),form='unformatted',access='direct',recl=Nx*Ny*Nz,status='replace')
   write(3333,rec=1) real(var,4)
   close(3333)
end subroutine 

subroutine save2d(Nx,Ny,var,filename)
   !Nx,Ny are not necessary in zonal and meridional directions.
   integer :: Nx,Ny
   character(len=*) :: filename
   real*8 :: var(Nx,Ny)

   open(3333,file=trim(filename),form='unformatted',access='direct',recl=Nx*Ny,status='replace')
   write(3333,rec=1) real(var,4)
   close(3333)
end subroutine 

subroutine w_pickup(filename)
   use header, only : h,u,v,w,uf,vf,wf,T,S
   character(len=*) :: filename
   open(3333,file=trim(filename),form='unformatted',access='stream',status='replace')
   write(3333) h,u,v,w,uf,vf,wf,T,S
   close(3333)
   open(3333,file=trim(filename)//'.meta',status='replace')
   write(3333,100) 'h',size(h,1),size(h,2),0,0
   write(3333,100) 'u',size(u,1),size(u,2),size(u,3),size(u,4)
   write(3333,100) 'v',size(v,1),size(v,2),size(v,3),size(v,4)
   write(3333,100) 'w',size(w,1),size(w,2),size(w,3),size(w,4)
   write(3333,100) 'uf',size(uf,1),size(uf,2),size(uf,3),0
   write(3333,100) 'vf',size(vf,1),size(vf,2),size(vf,3),0
   write(3333,100) 'wf',size(wf,1),size(wf,2),size(wf,3),0
   write(3333,100) 'T',size(T,1),size(T,2),size(T,3),size(T,4)
   write(3333,100) 'S',size(S,1),size(S,2),size(S,3),size(S,4)
!   write(3333,100) 'pcorr',size(pcorr,1),size(pcorr,2),size(pcorr,3)
   close(3333)
100 format (A,4I)
end subroutine w_pickup

subroutine r_pickup(step)
   use header, only : h,u,v,w,uf,vf,wf,T,S
   integer :: step
   character(len=10) :: stepchar


   open(3333,file='op.pickup.'//stepchar(step)//'.bin',form='unformatted',access='stream',status='old')
   write(3333) h,u,v,w,uf,vf,wf,T,S
   close(3333)
   print*, '# pickup '//'op.pickup.'//stepchar(step)//'.bin'
end subroutine r_pickup 


subroutine load_dy(dyM)
   use header, only : NJ
   real*8, intent(out) :: dyM(0:NJ+1)

   call io_error('dyM.data')

   open(3333,file='dyM.data',form='unformatted',access='stream',status='old')
   read(3333) dyM
   close(3333)
   print*, dyM
   print*, '# load dyM.data'
end subroutine load_dy

subroutine io_error(filename)
   character(len=*) filename
   logical :: lexist
   
   inquire(file=trim(filename),exist=lexist)

   if (.not. lexist) then
      print*, "# error: file "//trim(filename)//" does not exsit."
      stop
   endif
end subroutine io_error
