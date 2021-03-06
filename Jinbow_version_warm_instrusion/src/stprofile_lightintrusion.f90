subroutine stprofile 
  !     --------------------                                              
  USE header
  !     initializes s as pot.density with a single vertical profile       
  !     Across-front use the TANH profile with tightness as a measure of f
  !     spread                                                            
  !     Larger the factor, tighter the front and larger b_xx.             
  !     tightness=4, gives the needed b_xx to make b_xx H/ (2f^2) < 1.    
  !     tightness = 10 represents a very tight front, =1 loose front      
  !     mldepth is the desired ML depth                                   
  !========================================
  integer  i,j,k,it,iseed,npm1,iter,iter2 
  double precision rhovert(np),Tvert(np),dep(np),depoff(np),drho(np) 
  double precision bs1(np),cs1(np),ds1(np),bT1(np),cT1(np),dT1(np), &
       &     z,seval,zmm,sbkgrnd,z1,z2,zprev,slfacmin     
  double precision slfac0,slfac,dum,dscl,rdscl,yarg,ex2y,thy,ran3,  &
       perturb,slfacnew,dz,bfsqbkgrnd,wiggles,amplitude             
  !     Average rho profile for 40N from DATA/ts/Levitus_Atlas.nc         
  data dep / -5500.d0, -5000.d0, -4500.d0, -4000.d0, -3500.d0,      &
       -3000.d0, -2500.d0, -2000.d0, -1750.d0, -1500.d0, -1400.d0,  &
       -1300.d0, -1200.d0, -1100.d0, -1000.d0, -900.d0, -800.d0,    &
       -700.d0, -600.d0, -500.d0, -400.d0, -300.d0, -250.d0,        &
       -200.d0, -150.d0, -125.d0, -100.d0, -75.d0, -50.d0,          &
       -30.d0, -20.d0, -10.d0, -0.d0 /                             
  data rhovert /  1027.7694, 1027.7852, 1027.7894, 1027.7886,       &
       1027.7897, 1027.8050, 1027.8165, 1027.7848, 1027.7630,       &
       1027.7297, 1027.7115, 1027.6846, 1027.6727, 1027.6372,       &
       1027.5988, 1027.5540, 1027.5026,                             &
       1027.4337, 1027.3380, 1027.2524, 1027.1506, 1027.0222,       &
       1026.9653, 1026.8773, 1026.7458, 1026.6521,                  &
       1026.5293, 1026.3756, 1026.0901, 1025.7001, 1025.5029,       &
       1025.3558, 1025.3215 /                                       

!     double n2  (temporarily)                                          
!==   ---------------------                                             
!      do k=np-1,1,-1                                                   
!         drho(k)= rhovert(k+1)-rhovert(k)                              
!      end do                                                           
!      do k=np-1,1,-1                                                   
!         rhovert(k)= rhovert(k+1)-2.d0*drho(k)                         
!      end do                                                           
!     -----------------                                                 
!     Specify MLD                                                       
!=======================
!      mldepth= 100.d0
!      mldepth= 150.d0  
      mldepth= 200.d0 
!      mldepth= 250.d0  
!      mldepth= 300.d0                                                  

!       mldepth= 20.d0
                                                                        
!      tightness= 0.3d0   !used in inith_fixdh, for larger domain
      tightness= 0.03d0 
      bfsqbkgrnd = 1.d-6 
      write(6,*) 'ML background stratific N2 = ', bfsqbkgrnd 
      do k=1,np 
         depoff(k)= dep(k)- mldepth 
      end do 
                                                                        
      call spline (np,depoff,rhovert,bs1,cs1,ds1) 
                                                                        
!     sclwidth is the scale width of the channel (i.e. orig width = 48km
!     yfront is the distance of the front from the coast                
      sclwidth = 48.0 
      yfront = 0.5*(yc(NJ+1) +yc(0)) 
!     ADD a PERTURBATION to the width of the front                      
      iseed= 44294 
      dum = ran3(iseed) 
                                                                        
!     z1 is the depth of the ml (diff rho on both sides above this depth
!     z2 is the vertical extent of the density anamoly (it is gradually 
!        linearly anihillated with depth).                              
!     Orig vals. z1= 50. z2=250.                                        
      z1= mldepth - 10.d0 
      z2= mldepth + 10.d0 
                      
           
      slfac0=0.1
      slfac=slfac0
!      slfac=3.0*slfac0/4.0    
!      slfac=2.0*slfac0  
!      slfac=3.0*slfac0  
!      slfac=4.0*slfac0  
!       slfac=5.0*slfac0  
!      slfac=slfac0/5.0    

!      slfac=slfac0/4.0  
!      slfac=slfac0/2.0 

!      slfac= 0.1  !drho= 2*slfac    ! for 100-200 m deep ML,usual val by=1e-7
!      slfac=4.0*slfac0
!      slfacmin= 2.d0*slfac0   ! adding bz in ML
!      slfac= 0.5 ! slfac=0.02 for byby5, drho= 2*slfac for 100-200 m deep ML 
!      slfac=0.02             
      dfront= slfac        ! dfront will be saved in mymodules,for use in inith

      do j=0,NJ+1 
         do i=0,NI+1 
            do k=0,NK+1 
               z= DL*zc(i,j,k) 
                  if (z.ge.-mldepth) then 
                     sbkgrnd =                                          &
     &               seval(np,-1.*mldepth,depoff,rhovert,bs1,cs1,ds1)   &
     &                       - (z+mldepth)*bfsqbkgrnd*R0/(gpr*10.)      
                  else 
                   sbkgrnd =                                          &
                          seval(np,z,depoff,rhovert,bs1,cs1,ds1)        
                  end if
                  s(i,j,k,0)=  sbkgrnd 
!                  if (k.eq.NK) then                                    
!                     slfacnew= slfac*2.d0                              
!                     if ((-1.d0*slfacnew*thy).gt.0.d0)                 
!     &                    s(i,j,k,0)= -slfacnew*thy + s(i,j,k,0)       
!                  end if                                               
!                                                                       
!-               else                                                   
!-                  s(i,j,k,0)=                                         
!-     &                 seval(np,z,dep,svert,bs1,cs1,ds1) -S0          
!-               end if                                                 
            end do 
         end do 
      end do 
                                                                        
                                                                        
!      do iter2=1,40                                                    
!         call conadjust(0)                                             
!      end do                                                           
                                                                        
                                                                        
!     WIGGLE in FRONT                                                   
      wiggles=1.d0    ! number of wiggles
!      amplitude= 2.0d0  !THIS IS WHAT I USED IN MOST EXPTS with the WIGGLE
      amplitude = 0.d0   !FLAT FRONT
      do j=0,NJ+1 
         do i=0,NI+1 
            yfront= 0.5d0*(yc(NJ+1) +yc(0))                             &
     &           + amplitude*dsin(2.d0*PI*xc(i)/                        &
     &           (0.5*(xc(NI+1)+xc(NI))/wiggles)  )                     
            thy = tanh(tightness*(yc(j)-yfront)*PI) 
            !Same form is used in inith
                                                                        
            do k=1,NK 
               z= DL*zc(i,j,k) 
               if (z.gt.-mldepth) then 
                                                                        
                  if (z.ge.-z1) then 
                     slfacnew = slfac
!strat on light side          slfacnew = (slfac-slfacmin)*(1.d0 +z/mldepth)+slfacmin
                  else if (z.ge.-z2) then 
                     slfacnew = slfac*(z+z2)/(z2-z1) 
!strat onlight side                    slfacnew = slfacmin*(z+z2)/(z2-z1) 
                  else 
                     slfacnew = 0.d0 
                  end if 
                                                                        
                  s(i,j,k,0)= slfacnew*(thy-1.d0) + s(i,NJ,k,0) 
                                                                        
!     if (z.ge.-z1) then                                                
!                  sbkgrnd =                                            
!     &              seval(np,-1.*z1,depoff,rhovert,bs1,cs1,ds1)        
!     &              - (z+z1)*bfsqbkgrnd*R0/(gpr*10.)                   
!               s(i,j,k,0)= slfac*(1.d0+thy) + sbkgrnd                  
!     &               seval(np,-1.*mldepth,depoff,rhovert,bs1,cs1,ds1)  
!     &                    -(z+mldepth)*bfsqbkgrnd*R0/(gpr*10.)         
                                                                        
               end if 
            end do 
         end do 
      end do 
                                                                        
      do k=0,NK+1 
         do i=1,NI 
            s(i,0,k,0)= s(i,1,k,0) 
            s(i,NJ+1,k,0)= s(i,NJ,k,0) 
         end do 
!     periodicew                                                        
         do j=0,NJ+1 
            s(0,j,k,0)= s(NI,j,k,0) 
            s(NI+1,j,k,0)= s(1,j,k,0) 
         end do 
      end do 
                                                                        
                                                                        
      return 
      END                                           
