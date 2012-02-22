SUBROUTINE ini_rho()
  !     --------------------                                              
  USE header, ONLY: NI,NJ,NK,rho,s,T,pi,zc,DL,LEN,yc,xc,dx,dy,total_depth,T_ref
  USE grids
  IMPLICIT NONE

  ! Initializes s,T from Andrey's section gridded for model grid
  ! s_init.dat and T_init.dat  written in ascii - NK+2 x NJx2

  INTEGER n,i,j,k,Nu,k0
  REAL*8 :: A,B,C,D,G,potdens,x,y(0:NJ+1),z,tmp,yfront,mldepth,width,slfac
  real*8 :: tmp1, tmp2, tmp3, tmp4
  INTEGER,PARAMETER :: seed = 86456  
  print*, 'ini_grids'
  CALL RANDOM_SEED()

  ! assign the density or temperature and salinity by either analytic functions or 
  ! any particular hydrographic section.
  CALL ini_grids()

  ! t=A/B exp(A*z)erf(B*y) + G*(1+z)
  ! z=-1:0, y = -0.5,0.5, 

  ! ===========================
  ! cases 
  !-2: eddy
  ! 0: jet
  ! 1: idealized surface jet
  ! 2: krushio
  SELECT CASE (1)
  CASE(1)
     A = 4.0d0
     B = 80d0 !jet width in km as yc is in km
     C = 2d0
     G = 19d0
     D = 1d0  !linear temperature difference from south to north 
     !Nu =NK-9
     mldepth = 100d0
     s(:,:,:,0) = 34d0
     y = yc - (yc(NJ+1)-yc(0))/2d0
     DO k = NK, 0, -1
        DO j = 0, NJ+1
           DO i = 0, NI+1
              CALL RANDOM_NUMBER(tmp)
              !y = -(REAL(j)-REAL(NJ+1)/2d0)/REAL(NJ+1)
              !y = -(dble(j)-dble(NJ+1)/2d0)*dy
              ! y = - (yc(j) - (yc(NJ+1)-yc(0))/2d0)
              !IF (k>Nu) THEN
              IF (zc(1,1,k)*DL>-mldepth) THEN
                 z=0d0
                 Nu=k
              ELSE
                 z = dble(k-Nu)/dble(Nu)
              ENDIF
              T(i,j,k,0) = C * EXP(A*z) * erf(-y(j)/B) - D*y/maxval(y)*1d0  + G*(1+z) !+ 0.02*tmp ! + 1d-2*SIN(DBLE(i)/DBLE(NI+1)*2d0*pi*1d0)
              !T(i,j,k,0) = C * EXP(A*z) * erf(B*y) + G*(1+z) + 1d-3*SIN(DBLE(i)/DBLE(NI+1)*2d0*pi*3d0)
              rho(i,j,k) = potdens(s(0,j,k,0),T(0,j,k,0))
              !print*, "t=",t(1,j,k,0),"rho=",rho(1,j,k)
           ENDDO
        ENDDO
     ENDDO
     T(:,0,:,0)= T(:,1,:,0)
     T(:,NJ+1,:,0)= T(:,NJ,:,0)
     s(:,0,:,0)= s(:,1,:,0)
     s(:,NJ+1,:,0)= s(:,NJ,:,0)
     rho(:,0,:)= rho(:,1,:)
     rho(:,NJ+1,:)= rho(:,NJ,:)
     T(:,:,NK+1,:)=T(:,:,NK,:)
     s(:,:,NK+1,:)=s(:,:,NK,:)
     rho(:,:,NK+1)=rho(:,:,NK)
     DO k = 1, NK
        rho(:,:,k) = 0.25d0*(rho(:,:,k-1)+2d0*rho(:,:,k)+rho(:,:,k+1))
        s(:,:,k,0) = 0.25d0*(s(:,:,k-1,0)+2d0*s(:,:,k,0)+s(:,:,k+1,0))
        T(:,:,k,0) = 0.25d0*(T(:,:,k-1,0)+2d0*T(:,:,k,0)+T(:,:,k+1,0))
     ENDDO
  case(-4)
     s(:,:,:,:) = 34d0
     DO k = 0, NK+1
        z = zc(1,1,k)*DL
        T(:,:,k,:) = 28d0 + 5d0 * z/800d0
        DO i = 0, NI+1
           DO j = 0, NJ+1
              rho(i,j,k) = potdens(s(i,j,k,0),T(i,j,k,0))
           enddo
        enddo
      enddo

  case(-3)
     call stprofile
  CASE(-2)
     A = 1025.3
     mldepth = 50d0
     width=100d3
     B = 200d0
     C = 300d0
     G = 1025.3d0
     D = 0d0
     slfac = 0.01
     s(:,:,:,0) = 34d0
     yfront = 0.5d0*(yc(NJ/2)+yc(NJ/2+1))*dy
     tmp1 = 0.5d0*(xc(NI/2)+xc(NI/2+1))*dx
     PRINT*, "# yfront = ", yfront
     PRINT*, "zc(1,1,:)=",zc(1,1,:)*DL
y = yc - (yc(NJ+1)-yc(0))/2d0
     DO j = 0, NJ+1
        DO i = 0, NI+1
           DO k = NK+1,0,-1
              z = zc(1,1,k)*DL
              IF (z > -mldepth) THEN
                 s(i,j,k,0) = G  - 2e-5*z
              ELSEIF (z > -C) THEN
                 z = zc(1,1,k)*DL + mldepth
                 s(i,j,k,0) = G  - 1e-4*z&
                      + slfac * (EXP(z/B)) 
                 k0 = k
              ELSE             
                 z = zc(1,1,k)*DL + 300
                 tmp = (1d0-EXP(z/50d0))*EXP(z/100d0)
                 s(i,j,k,0) = s(i,j,k0,0)  - 1e-5 * z
              ENDIF
              CALL RANDOM_NUMBER(tmp)
           ENDDO
        ENDDO
     ENDDO
     DO j = 0, NJ+1
        y=yc(j)*dy
        DO i = 0, NI+1
           x=xc(i)*dx
           DO k = NK+1,0,-1
              z = zc(1,1,k)*DL
                 s(i,j,k,0) = s(i,j,k,0) + &
                 ( slfac*exp(z/B)*(1-exp(z/B))* &
                 exp(-((y-yfront)/width/2d0)**2- &
                 ((x-tmp1+z/total_depth*100d3)/width/2d0)**2))
           ENDDO
        ENDDO
     ENDDO
     DO k = 1, NK
        s(:,:,k,0) = 0.25d0*(s(:,:,k-1,0)+2d0*s(:,:,k,0)+s(:,:,k+1,0))
     ENDDO
     DO k = 1, NK
        s(:,:,k,0) = 0.25d0*(s(:,:,k-1,0)+2d0*s(:,:,k,0)+s(:,:,k+1,0))
     ENDDO
     s(:,0,:,:)= s(:,1,:,:)
     s(:,NJ+1,:,:)= s(:,NJ,:,:)
     s(:,:,NK+1,:)=s(:,:,NK,:)
     s(:,:,0,:) = s(:,:,1,:)
     rho = s(:,:,:,0)
     call save3d(NI+2,NJ+2,NK+2,rho,'op.rho.init.bin')
  CASE(0)
     s(:,:,:,0) = 34d0
     DO k = 0, NK
        DO j = 0, NJ+1
           CALL RANDOM_NUMBER(tmp)
           y = (REAL(j)-REAL(NJ+1)/2d0)/REAL(NJ+1)
           T(:,j,k,0) = 6d0 * erf(9d0*y) + 15d0 + 1d-2*tmp
           rho(:,j,k) = potdens(s(0,j,k,0),T(0,j,k,0))
        ENDDO
     ENDDO
     T(:,:,NK+1,:)=T(:,:,NK,:)
     s(:,:,NK+1,:)=s(:,:,NK,:)
     rho(:,:,NK+1)=rho(:,:,NK)
!     s(:,:,:,0) = rho
     !OPEN(212, file='tmp-rho-ini.bin',form='unformatted',access='stream')
!     call save3d(NI+2,NJ+2,NK+2,rho,'op.rho.init.bin')
     !WRITE(212) rho
     !CLOSE(212)

  CASE(-1)
     A = 1025.3
     mldepth = 200d0
     width=5d3
     B = 200d0
     C = 500d0
     G = 1025.3d0
     slfac = 0.05
     s(:,:,:,0) = 34d0
     yfront = 0.5d0*(yc(NJ/2)+yc(NJ/2+1))*dy
     PRINT*, "# yfront = ", yfront
     PRINT*, "zc(1,1,:)=",zc(1,1,:)*DL
     DO j = 0, NJ+1
        y=yc(j)*dy
        DO i = 0, NI+1
           DO k = NK+1,0,-1
              z = zc(1,1,k)*DL
              IF (z > -mldepth) THEN
                 s(i,j,k,0) = G +  slfac * erf((y-yfront)/width/2d0)
              ELSEIF (z > -C) THEN
                 z = zc(1,1,k)*DL + mldepth
                 s(i,j,k,0) = G  - 5e-4*z&
                      + slfac * (EXP(z/B)) * erf((y-yfront)/width/2d0)
                 k0 = k
              ELSE             
                 z = zc(1,1,k)*DL + 300
                 tmp = (1d0-EXP(z/50d0))*EXP(z/100d0)
                 s(i,j,k,0) = s(i,j,k0,0)  - 1e-5 * z
              ENDIF
              CALL RANDOM_NUMBER(tmp)
           ENDDO
        ENDDO
     ENDDO
     DO k = 1, NK
        s(:,:,k,0) = 0.25d0*(s(:,:,k-1,0)+2d0*s(:,:,k,0)+s(:,:,k+1,0))
     ENDDO
     s(:,0,:,:)= s(:,1,:,:)
     s(:,NJ+1,:,:)= s(:,NJ,:,:)
     s(:,:,NK+1,:)=s(:,:,NK,:)
     s(:,:,0,:) = s(:,:,1,:)
     rho = s(:,:,:,0)
     OPEN(212, file='tmp-rho-ini.bin',form='unformatted',access='stream')
     WRITE(212) REAL(rho,4)
     CLOSE(212)
  CASE(2)
     n=0
     OPEN(31,file='s_kuroshio65.dat')
     DO k=0,NK
        READ(31,*) (s(0,j,k,n),j=1,NJ)
        s(0,0,k,n)= s(0,1,k,n)
        s(0,NJ+1,k,n)= s(0,NJ,k,n)
     END DO
     PRINT*,"-----"
     CLOSE(31)
     OPEN(31,file='T_kuroshio65.dat')
     DO k=0,NK
        READ(31,*) (T(0,j,k,n),j=1,NJ)
        T(0,0,k,n)= T(0,1,k,n)
        T(0,NJ+1,k,n)= T(0,NJ,k,n)
     END DO
     CLOSE(31)

     i=0
     DO j=0,NJ+1
        s(i,j,NK+1,n)= s(i,j,NK,n) 
        T(i,j,NK+1,n)= T(i,j,NK,n) 
     END DO

     DO k=0,NK+1
        DO j=0,NJ+1
           DO i=1,NI+1
              s(i,j,k,n)= s(0,j,k,n)
              T(i,j,k,n)= T(0,j,k,n)
           END DO
        END DO
     END DO
     CALL evalrho(rho,n)
  CASE(3)

     n=0
     OPEN(unit=31,file='s_kuroshio65.dat')
     READ(31,*) ((s(0,j,k,n),j=1,NJ),k=0,NK)
     OPEN(unit=32,file='T_kuroshio65.dat')
     READ(32,*) ((T(0,j,k,n),j=1,NJ),k=0,NK)
     DO k=0,NK
        DO j=1,NJ
           CALL RANDOM_NUMBER(tmp)
           rho(0,j,k) = potdens(s(0,j,k,0),T(0,j,k,0)) + 1d-2*tmp
        ENDDO
     END DO
     CLOSE(31)
     CLOSE(32)

     DO k=0,NK
        DO j=1,NJ
           s(:,j,k,0) = s(0,j,k,0)
           T(:,j,k,0) = T(0,j,k,0)
           rho(:,j,k) = rho(0,j,k)
        ENDDO
     ENDDO

     T(:,0,:,0)= T(:,1,:,0)
     T(:,NJ+1,:,0)= T(:,NJ,:,0)
     s(:,0,:,0)= s(:,1,:,0)
     s(:,NJ+1,:,0)= s(:,NJ,:,0)
     rho(:,0,:)= rho(:,1,:)
     rho(:,NJ+1,:)= rho(:,NJ,:)

     T(:,:,NK+1,0)=T(:,:,NK,0)
     s(:,:,NK+1,0)=s(:,:,NK,0)
     rho(:,:,NK+1)=rho(:,:,NK)
  CASE(4)
     CALL stprofile
     call save3d(NI+2,NJ+2,NK+2,rho,'initial-rho.bin')
  CASE(5)
     CALL stprofile
  END SELECT

  OPEN(1,file='tmp-rho-t0.bin',form='unformatted',access='stream',status='replace')
  PRINT*, "write the initial density field"
  WRITE(1) rho
  CLOSE(1)
  T_ref = T(:,:,:,0)

END SUBROUTINE ini_rho
