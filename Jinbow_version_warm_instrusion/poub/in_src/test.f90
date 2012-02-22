MODULE mod_param
 INTEGER                                   :: IMT, JMT, KM
 INTEGER                                   :: JMAX, LBT, NTRACMAX
 INTEGER, PARAMETER                        :: MR=1001
 INTEGER                                   :: NEND
 INTEGER, PARAMETER                        :: NST=2,NNRJ=8,NTRJ=7
#ifdef streamts
 INTEGER, PARAMETER                        :: LOV=3
#else
 INTEGER, PARAMETER                        :: LOV=1
#endif
 INTEGER                                   :: ncoor,kriva,iter,ngcm
 REAL*8                                    :: tseas,tday,tyear,dtmin,voltr
 REAL*8                                    :: tstep,dstep,tss,partQuant
 REAL*8, PARAMETER                         :: UNDEF=1.d20 
 REAL*8, PARAMETER                         :: PI = 3.14159265358979323846d0
ENDMODULE mod_param


allocate ( uflux(imt,jmt,km,nst), vflux(imt,0:jmt,km,nst) )

