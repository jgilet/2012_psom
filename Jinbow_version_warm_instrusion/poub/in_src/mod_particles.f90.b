MODULE particles
  USE header, ONLY : NI,NJ,NK,uf,vf,wf,Jifc,Jjfc,J2d,ux,vy,wz,NP,PI,dtf,vor,shear,rho,strain
  ! define the class for particles
  TYPE particle
     REAL*8 :: i,j,k,x,y,z,u,v,w,s,t,u0,v0,w0,id,vor,strain,shear,rho,time
  END TYPE particle
  ! NP is the particle number, which is assigned in the namelist.
  !  integer :: NP
  TYPE (particle), DIMENSION(:), ALLOCATABLE :: parti
  REAL*8 :: pcx,pcy,pcz,pcr ! the initial position of the particles

CONTAINS

  SUBROUTINE ini_particles(time)
    USE header, ONLY : rho
    IMPLICIT NONE
    INTEGER :: i,time
    REAL*8 :: rand,r,theta,xcenter,ycenter,isovalue,ywidth,xwidth
    INTEGER,PARAMETER :: seed = 86456  
    logical :: isolower(0:NI+1,0:NJ+1,0:NK+1)
    CALL RANDOM_SEED()
    Do i=1, NP
          parti(i)%id = i
          parti(i)%time = dble(time)
          parti(i)%u0=0d0
          parti(i)%v0=0d0
          parti(i)%w0=0d0
          parti(i)%u=0d0
          parti(i)%v=0d0
          parti(i)%w=0d0
          parti(i)%t=0d0
          parti(i)%s=0d0
          parti(i)%rho=0d0
          parti(i)%vor=0d0
          parti(i)%shear=0d0
          parti(i)%strain=0d0
    enddo 

    SELECT CASE (1)
    CASE (0)
       DO i = 1, NP
          parti(i)%j = REAL(i/250)*3 + real(NJ/2-1)
          parti(i)%i = Real(mod(i,250))
          PRINT*, "rand",theta
          parti(i)%k = REAL(NK)*pcz
       ENDDO
    CASE (1)
       DO i = 1, NP
          parti(i)%id = i
          CALL RANDOM_NUMBER(r)
          CALL RANDOM_NUMBER(theta)
          theta  = (theta-0.5d0)*2d0*PI
          parti(i)%i = REAL(NI/2) + pcx*REAL(NI)*r*SIN(theta)
          parti(i)%j = REAL(NJ/2) + pcy*REAL(NJ)*r*COS(theta)
          parti(i)%k = REAL(NK)*pcz
       ENDDO
    CASE (2)
       DO i=1,NP
          parti(i)%j = REAL(i/50) + REAL(NJ/3)
          parti(i)%i = REAL(MOD(i,50))
          parti(i)%k = 15
       ENDDO
    case (3)
       xcenter = real(NI/2)
       ycenter = real(NJ/2)
       xwidth  = 10d0
       ywidth  = 10d0
       isovalue = 1025.2
       !isolower = where(rho>isovalue)
       !print*, isolower
       Do i=1,NP
         !parti(i)%i = real(mod(i,ywidth))
         parti(i)%j = real(i/ywidth) - parti(i)%i
         !to be finished, locate z at any iso
       
       enddo
          

    END SELECT

  END SUBROUTINE ini_particles

  SUBROUTINE get_parti_vel(time)
    IMPLICIT NONE
    INTEGER :: i,j,k,ip,ic,jc,kc,ifc,jfc,kfc,time
    real*8 :: dic,djc,dkc,dif,djf,dkf
    REAL*8, DIMENSION(    0:NI,0:NJ+1        )        :: uxf
    REAL*8, DIMENSION(    0:NI+1,0:NJ        )        :: vyf
    REAL*8, DIMENSION(    0:NI+1,0:NJ+1,0:NK )        :: wzf
    REAL*8, DIMENSION(      NI,  0:NJ,     NK)        :: vfp
    REAL*8, DIMENSION(    0:NI+1,0:NJ,   0:NK+1)      :: vf_ex
    REAL*8, DIMENSION(    0:NI,    NJ,     NK)        :: ufp
    REAL*8, DIMENSION(    0:NI,  0:NJ+1, 0:NK+1)      :: uf_ex
    REAL*8, DIMENSION(      NI,    NJ,   0:NK)        :: wfp
    REAL*8, DIMENSION(    0:NI+1,0:NJ+1, 0:NK)        :: wf_ex

    !rearrange the ux and vy to face grids
    !uxf = 0.5d0*(ux(0:NI,:)+ux(1:NI+1,:))
    !vyf = 0.5d0*(vy(:,0:NJ)+vy(:,1:NJ+1))
    wzf = 0.5d0*(wz(:,:,0:NK) + wz(:,:,1:NK+1))

    !calculate the face velocity
    k=0
    wfp(:,:,k) = wf(:,:,k)/J2d(1:NI,1:NJ)*wzf(1:NI,1:NJ,k)

    DO k = 1, NK
       ufp(:,:,k) = uf(:,:,k)/Jifc(:,:,k)
       vfp(:,:,k) = vf(:,:,k)/Jjfc(:,:,k)
       wfp(:,:,k) = wf(:,:,k)/J2d(1:NI,1:NJ)*wzf(1:NI,1:NJ,k)
    ENDDO
    uf_ex=0d0
    uf_ex(:,1:NJ,1:NK) = ufp
    !=== vertical extrapolation
    uf_ex(:,:,NK+1) = 2*uf_ex(:,:,NK)-uf_ex(:,:,NK-1) ! extrapolation

    vf_ex=0d0
    vf_ex(1:NI,:,1:NK) = vfp
    !=== zonally periodic
    vf_ex(0,:,:) = vf_ex(NI,:,:)
    vf_ex(NI+1,:,:)=vf_ex(1,:,:)
    !=== vertical extrapolation
    vf_ex(:,:,NK+1) = 2*vf_ex(:,:,NK)-vf_ex(:,:,NK-1)

    wf_ex=0d0
    wf_ex(1:NI,1:NJ,:) = wfp
    !===zonal periodic condition
    wf_ex(0,:,:) = wf_ex(NI,:,:)
    wf_ex(NI+1,:,:) = wf_ex(1,:,:)
    !===
    DO ip = 1, NP
       parti(ip)%time=dble(time)
       IF (parti(ip)%j < NJ .AND. parti(ip)%j > 0 .AND. &
            parti(ip)%k < NK .AND. parti(ip)%k > 0) THEN
          !ic, jc, kc, is the integer index of the particle relative to 
          !the grids center. Use these values for variables wih the ghost points.
          !ifc, jfc, and kfc is the index relative to the coordinates of grid faces. 
          !Use these values for variables on faces.
          ic = INT(parti(ip)%i+0.5d0)
          jc = INT(parti(ip)%j+0.5d0)
          kc = INT(parti(ip)%k+0.5d0)

          ifc = INT(parti(ip)%i)
          jfc = INT(parti(ip)%j)
          kfc = INT(parti(ip)%k)

          dif = parti(ip)%i - ifc
          djf = parti(ip)%j - jfc
          dkf = parti(ip)%k - kfc

          dic = parti(ip)%i - ic + 0.5d0
          djc = parti(ip)%j - jc + 0.5d0
          dkc = parti(ip)%k - kc + 0.5d0
          !calcuate the zonal velocity 
          !i = INT(parti(ip)%i)
          !j = INT(parti(ip)%j+0.5d0)
          !k = INT(parti(ip)%k+0.5d0)
          !di = parti(ip)%i - ifc
          !dj = parti(ip)%j - jc + 0.5d0
          !dk = parti(ip)%k - kc + 0.5d0
          CALL interp_trilinear(dif,djc,dkc,uf_ex(ifc:ifc+1,jc:jc+1,kc:kc+1),parti(ip)%u)
          !calcuate the meridional velocity
          !i = INT(parti(ip)%i+0.5d0)
          !j = INT(parti(ip)%j)
          !k = INT(parti(ip)%k+0.5d0)
          !di = parti(ip)%i - ic + 0.5d0
          !dj = parti(ip)%j - jfc
          !dk = parti(ip)%k - kc + 0.5d0
          CALL interp_trilinear(dic,djf,dkc,vf_ex(ic:ic+1,jfc:jfc+1,kc:kc+1),parti(ip)%v)
          !calcuate the vertical velocity
          !i = INT(parti(ip)%i+0.5d0)
          !j = INT(parti(ip)%j+0.5d0)
          !k = INT(parti(ip)%k)
          !di = parti(ip)%i - ic + 0.5d0
          !dj = parti(ip)%j - jc + 0.5d0
          !dk = parti(ip)%k - kfc
          CALL interp_trilinear(dic,djc,dkf,wf_ex(ic:ic+1,jc:jc+1,kfc:kfc+1),parti(ip)%w)
          !diagnose other properties
          CALL interp_trilinear(dic,djc,dkc,vor(ic:ic+1,jc:jc+1,kc:kc+1),parti(ip)%vor)
          CALL interp_trilinear(dic,djc,dkc,rho(ic:ic+1,jc:jc+1,kc:kc+1),parti(ip)%rho)
          CALL interp_trilinear(dic,djc,dkc,shear(ic:ic+1,jc:jc+1,kc:kc+1),parti(ip)%shear)
          CALL interp_trilinear(dic,djc,dkc,strain(ic:ic+1,jc:jc+1,kc:kc+1),parti(ip)%strain)
       ELSE
          parti(ip)%u=0d0
          parti(ip)%v=0d0
          parti(ip)%w=0d0
       ENDIF
    ENDDO
    ! get the zonal velocity u

  END SUBROUTINE get_parti_vel

  SUBROUTINE parti_forward()
    IMPLICIT NONE
    INTEGER :: i
    DO i = 1, NP
       parti(i)%i = parti(i)%i + 0.5d0 * dtf * (3d0 * parti(i)%u - parti(i)%u0)
       IF (parti(i)%i >NI)   parti(i)%i = parti(i)%i - REAL(NI)
       IF (parti(i)%i <0d0 ) parti(i)%i = parti(i)%i + REAL(NI)

       if (parti(i)%j>NJ-1 .and. parti(i)%v>0) then
       parti(i)%j = parti(i)%j + parti(i)%v * dtf / (1d0 + (parti(i)%v * dtf)/(dble(NJ) - parti(i)%j) )
       else if (parti(i)%j<1 .and. parti(i)%v<0) then
       parti(i)%j = parti(i)%j + parti(i)%v * dtf / ( 1d0 - dtf/parti(i)%j )
       else
       parti(i)%j = parti(i)%j + 0.5d0 * dtf * (3d0 * parti(i)%v - parti(i)%v0)
       endif

       if (parti(i)%k>NK-1 .and. parti(i)%w>0) then
       parti(i)%k = parti(i)%k + parti(i)%w * dtf / (1d0 + (parti(i)%w * dtf)/(dble(NK) - parti(i)%k) )
       else if (parti(i)%k<1 .and. parti(i)%w<0) then
       parti(i)%k = parti(i)%k + parti(i)%w * dtf / ( 1d0 - dtf/parti(i)%k ) 
       else
       parti(i)%k = parti(i)%k + 0.5d0 * dtf * (3d0 * parti(i)%w - parti(i)%w0)
       endif
   
       !debug part
       if (parti(i)%j<0d0 .or. parti(i)%j>NJ .or. parti(i)%k>NK .or. parti(i)%k<0d0 ) then
          print*, "particles coordinates are wrong"
          stop
       endif

       parti(i)%u0 = parti(i)%u
       parti(i)%v0 = parti(i)%v
       parti(i)%w0 = parti(i)%w
    ENDDO

  END SUBROUTINE parti_forward


  SUBROUTINE interp_trilinear(di,dj,dk,var,velp)
    IMPLICIT NONE
    REAL*8, INTENT(in) :: di,dj,dk
    REAL*8, INTENT(in), DIMENSION(    2,    2  ,   2  )      :: var
    REAL*8, INTENT(out) :: velp
    REAL*8 :: i1,i2,i3,i4,j1,j2,ti,tj,tk

    ! calcuate the Trilinear interpolation
    i1 = (var(2,1,  1)   - var(1,1,  1))*di + var(1,1,  1)
    i2 = (var(2,1,  2) - var(1,1,2))*di + var(1,1,  2)
    i3 = (var(2,2,2) - var(1,2,2))*di +var(1,2,2)
    i4 = (var(2,2,1)   - var(1,2,1))*di + var(1,2,1)

    j1 = (i3 - i2)*dj + i2
    j2 = (i4 - i1)*dj + i1

    velp = (j1 - j2) * dk + j2
  END SUBROUTINE interp_trilinear

  SUBROUTINE  get_parti_vel_ana()
    INTEGER :: ip
    DO ip = 1, NP
       parti(ip)%u = -1*SIN(pi*parti(ip)%i/REAL(NI))*COS(pi*parti(ip)%j/REAL(NJ))
       parti(ip)%v = COS(pi*parti(ip)%i/REAL(NI))*SIN(pi*parti(ip)%j/REAL(NJ))
    ENDDO
  END SUBROUTINE get_parti_vel_ana

  SUBROUTINE save_parti()
    IMPLICIT NONE
    INTEGER :: i
    DO i = 1, NP
       WRITE(125) parti(i)%id, parti(i)%i,   parti(i)%j,     parti(i)%k, &
                               parti(i)%u,   parti(i)%v,     parti(i)%w, &
                               parti(i)%rho, parti(i)%s,     parti(i)%t, &
                               parti(i)%vor, parti(i)%shear, parti(i)%strain
    ENDDO
  END SUBROUTINE save_parti


END MODULE particles
