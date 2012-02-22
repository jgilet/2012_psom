
PSOM: Process Study Ocean Model
===============================

Description of this version: 2012/02/22 version, obtained by JBG, from JBW code. 

Modifications from JBW version:

  copied stprofile_X

  copied viscous

  changed namelist

  changed grid.h

  got rid of sponge

  got rid of load_dy

  got rid of load_r

  got rid of relaxation

  ****************************
  added:

    229 !     CALL ini_rho()

    230      call stprofile

    231      call evalrho(rho,0)

  ****************************
  shifted from variable to constant grid

    100             ydv(i,j)= dy/LEN !- used for constant dy

    101 !           ydv(i,j)= dyM(j)/LEN 

  ***************************
  added in mod_grids:

    yc(0) = -0.5*dy*1.d-3 

    do j=1,NJ+1 

      yc(j)= yc(j-1) + dy*1.d-3 

    end do  

    xc(0) = -0.5*dx*1.d-3 

    do i=1,NI+1 

      xc(i)= xc(i-1) + dx*1.d-3 

    end do  

    declared them properly

  ***************************
  added write*.f90 subroutines for visualization purposes.

  in these, I replaced USE dimensions by : USE header,ONLY : NI,NJ,NK,ntr,nconsume


  *************************
  pb in geostroph at t=0 

  tried to replace by old geostroph, solved by adding the complilation option -real-size 64.


  **************************
  Shifted from f to beta-plane

  Modified init_tr

    56 DO  j=0,NJ+1

    57      latrad(j)= phi0 + DBLE(j-jmid)*dphi 

    58      !latrad(j)= phi0    

    59   enddo

  ************************
  Set RR=0 in diffusion_wind  


====================================

Description of the model:
PSOM, pronounced "soam" (the nectar derived from the churning of the oceans in Indian mythology), stands for Process Study Ocean Model. It is a versatile, three-dimensional, non-hydrostatic, computational fluid dynamical model for oceanographic (as well as other) applications (Mahadevan et al., 1996a,b). The model uses the finite volume method on a structured grid with boundary fitted coordinates (topography conforming sigma grid in the vertical, and boundary conforming in the horizontal). The model has a free-surface. It can be used for large- and small-scale phenomena and can be run in hydrostatic or non-hydrostatic mode (Mahadevan, 2006). It uses a highly efficient solution procedure that exploits the skewness arising from the small geometrical aspect (depth to length scale) ratio of the ocean to speed up the solution of the non-hydrostatic pressure, which is solved by the multigrid method. The model has been used for a number of process studies, including investigation of the vertical transport of nutrients for phytoplankton production (Mahadevan and Archer, 2000) and the dynamics of submesoscale processes (Mahadevan and Tandon, 2006; Mahadevan, Tandon and Ferrari, 2010). Since the non-hydrostatic model is well-posed with open boundaries, it can be used as a nested high-resolution model with time-varying boundary conditions applied from a coarser resolution general circulation model (Mahadevan and Archer, 1998). The model is thus ideally suited for high-resolution, limited-region modeling studies. 

