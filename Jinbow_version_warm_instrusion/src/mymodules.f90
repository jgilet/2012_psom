MODULE header
#include "cppdefs.h"

  INCLUDE 'size.h'

  integer :: NPR,ini_particle_time,parti_file_num

  LOGICAL, PARAMETER :: rect = .TRUE., periodicew = .TRUE.
!*** S0 and T0 are the mean salinity and temp. sreal= S0+s,Treal=T0+T  
  REAL*8,  PARAMETER :: S0=35.7d0, T0=15.d0 ,R0=1027.d0
!*** EPS is the Rosby number
!*** OMEGA is the angular velocity of rotation of the earth
  REAL*8,  PARAMETER :: EPS= 0.1d0, AL=1.d7, FPAR=1.d-4,                     &
                        PI=3.14159265358979323846, OMEGA=7.272d-5,           &
                        gpr= 0.981,apr=0.6371,dztop=1.d-3,z0= 0.2d0,zm= 0.1d0
!*** delta is the ratio of horizontal to vertical length scale = L/D   
!*** delinv= 1/delta
  real*8 :: delta,delinv
!*** L, DL are the characteristic length and depth scales.  
!*** lexp is a constant in the equations denoting the length scale     
!*** LEN= 10**(4+lexp),  0<=lexp<=1        
  REAL*8, PARAMETER :: LEN= 1.d5, DL=1.d3
!*** dtime : time step
!*** dtf   : dimensionless time step
!*** fnhhy : 0 for hydrostatic, 1 for nonhydrostatic
  REAL*8  :: dtf,fnhhy,dx,dy
  integer :: pickup_step
!  REAL*8  :: dtf,fnhhy,dx,dy,dtime
!*** sigrelease is the isopycnal of tracer release
!*** dztop is the depth of the top layer non-dim by DL
!*** pfac is the grid stretching factor in z, used in findzall and sigma
  REAL*8  :: phi0deg
  REAL*8  :: yfront,dyfront,sclwidth,tightness,mldepth,total_depth
  REAL*8  :: stressmax,pfac,sigrelease(ntr)
  INTEGER :: conv(0:NI+1,0:NJ+1,0:NK+1)
  INTEGER :: con100(0:NI+1,0:NJ+1,0:NK+1)
  INTEGER :: ktest(3)
!*** qpr is the ratio of non-hyd to hyd pressure = Q/P
!*** qprinv = 1/qpr
!*** lambda is the ration of hor length scale to earth's rad = L/A 
!*** kappah is the implcitness parameter in the Crank-Nicolson scheme f
  REAL*8               :: sbackgrnd, beta,P1,qpr,lambda, &
        &                 UL,WL,TL,HL,HDL,kappah,kaphinv

  real*8 :: Kx_TS, Ky_TS_bdy,Ky_TS_int, Kx_m, Ky_m
  real*8 :: r_sponge(0:NJ+1)

  REAL*8, dimension(3) :: mldn2,mldn2init,zbtop,zbbot,frictop,diatop,     &
        &                 zbtopinit,zbbotinit,advecpv,friction,           &
        &                 diabatic,diabot,fricbot
 
  REAL*8, dimension(    0:NI+1,0:NJ+1, 0:NK+1,nconsume) :: consump 
  REAL*8, target, dimension(ntr,0:NI+1,0:NJ+1, 0:NK+1,0:1)  :: Tr
  REAL*8, pointer, dimension(:,:,:)  :: Tr_p
  REAL*8, dimension(ntr,0:NI+1,0:NJ+1, 0:NK+1    )  :: wtr
  REAL*8, dimension(    0:NI+1,0:NJ+1, 0:NK+1,0:1)  :: u,v,w,s,T
  REAL*8, dimension(    0:NI+1,0:NJ+1, 0:NK,  3  )  :: gqk
  REAL*8, dimension(    0:NI,    NJ,     NK,  2  )  :: gi,gqi
  REAL*8, dimension(      NI,  0:NJ,     NK,  2  )  :: gj,gqj
  REAL*8, dimension(    0:NI+1,0:NJ+1,-1:NK+1)      :: zf
  REAL*8, dimension(    0:NI+1,0:NJ+1, 0:NK+1)      :: zc,wx,wy,wz,p,strain,shear,Jac,rho,&
        &                                              vor,pv,freqN2,si,sj,sk,cx,cy,cz,rp,T_ref
  REAL*8, dimension(    0:NI+1,0:NJ+1, 0:NK  )      :: wt,wzk,skfc
  REAL*8, dimension(      NI,    NJ,   0:NK  )      :: czf,Kz,wf
  REAL*8, dimension(           0:NJ+1,   NK  )      :: ueast,uwest
  REAL*8, dimension(      NI,  0:NJ,     NK  )      :: hyn,sjfc,cyf,gj3,gqj3,Jjfc,vf
  REAL*8, dimension(    0:NI,    NJ  ,   NK  )      :: hxn,sifc,cxf,gi3,gqi3,Jifc,uf
  REAL*8, dimension(      NI,    NJ,     NK  )      :: uvis,vvis,wvis,fricu,fricv, &
        &                                              fricw,fricb,rhoadv,rhoprev
  REAL*8, dimension(             NJ,     NK  )      :: ufbce,ufbcw,trinit,divreyn, &
        &                                              divmean,dcdt,prod
  REAL*8, dimension(      NI,            NK  )      :: vfbcn,vfbcs,vnorth,vsouth,ssouth
  REAL*8, dimension(    0:NI+1,0:NJ+1,   2   )      :: gradhn
  REAL*8, dimension(    0:NI+1,0:NJ+1        )      :: ux,uy,vx,vy,ffc,bbc,oldh,h,hdt,D, &
        &                                              J2d,Ddx,Ddy,g11,g22,g12
  REAL*8, dimension(      NI  ,0:NJ          )      :: bbj,ffj
  REAL*8, dimension(    0:NI  ,  NJ          )      :: ffi,bbi
  REAL*8, dimension(      NI,    NJ          )      :: wfbcb
  REAL*8, dimension(           0:NJ+1        )      :: yc,latrad
  REAL*8, dimension(    0:NI+1               )      :: xc
  REAL*8, dimension(             NJ          )      :: stressx

!rpgrads
  REAL*8 :: drpx(NI,NJ,NK),drpy(NI,NJ,NK),                                  &
            grpifc(0:NI,NJ,NK),grpjfc(NI,0:NJ,NK)

END MODULE header

