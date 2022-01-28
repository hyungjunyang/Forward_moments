c-----------------------------------------------------------------------
c     Construction of Penality Function                                !
c     **********************************                               !
c                                                                      !
c     This subroutine construct the penality function for optimization !
c     The penality function is constructed from pressure difference and!
c     sensitivity coefficients                                         !
c===================================================================70=c
c2345678901234567890123456789012345678901234567890123456789012345678901!
c     *************************************                            !
c     Compute                                                          !
c----------------------------------------------------------------------!
c     PARAMETER LIST                                                   !
c                                                                      !
c     mxtpe    maximum number of parameters in optimization            !
c              (=mxtmed+mxmp)                                          !
c              mxtmed :  maximum number of measured ln(k) data         !
c              mxmp   :  maximum number of master points               !
c     nmpx     no. of columms to be generated for master points        !
c     nmpy     no. of rows    to be generated for master points        !
c     nmp=nmpx*nmpy                                                    !
c     ntmed    number of permeability measurements                     !
c     xtms     the x-coordinate of the master points                   !
c     ytms     the y-coordinate of the master points                   !
c     nxpm,    the number of nodes in the x-, y- directions            !
c     xm1,xm2  the coordinates of x and y directions of perm filed     !
c                                                                      !
c     PERMY(Nr,NXPM,NYPM) = ln K (Nr,NXPM,NYPM)                        !
c     PERMF(Nr,NX*NY,5) = K(i+0.5,j)/<K>                               !
c     htj   :  sensitivity coefficients                                !
c                                                                      !
c     input :  htj, tppe, wpms, pms, pcal, tms, vtms                   !
c     ouput :  cpq, ppq, apq1, bpq, xpq                                !
c----------------------------------------------------------------------!
c     GRID DISCRETIZATION                                              !
c     Pointed Distributed                                              !
c                                                                      !
c=====7================================================================c
      subroutine buildfd(mxwell,mxtpe,mxstep,ntmed,nmp,nwell,nstep,             
     +                   htj,htj1,wpms,sumdiv,                                  
     +                   pms,pcal,pstd,                                         
     +                   am_y,tms,vtms,                                         
     +                   cpq,ppq,apq1,bpq,tppe,xpq                              
     +                  )                                                       
      implicit real*8 (a-h,o-z)                                                 
      dimension htj(mxwell,mxtpe,mxstep)                                        
      dimension htj1(mxwell,mxtpe,mxstep)                                       
      dimension wpms(mxwell,mxstep)                                             
      dimension pms(mxwell,mxstep),pcal(mxwell,mxstep)                          
      dimension pstd(mxwell,mxstep)                                             
      dimension tms(mxtpe),vtms(mxtpe)                                          
      dimension cpq(mxtpe,mxtpe),                                               
     *          ppq(mxtpe), tppe(mxtpe), xpq(mxtpe),                            
     *          apq1(2 * mxtpe), bpq(2 * mxtpe)                                 
                                                                                
C1)   ===== calculating cpq or c matrix =====
      do i = 1, ntmed + nmp    
         do j = i, ntmed + nmp
            ac = 0.0
            do itt = 1, nstep
               do k = 1, nwell
                  ac = ac + (htj (k,i,itt) * htj (k,j,itt) * wpms(k,itt)
     *                      +htj1(k,i,itt) * htj1(k,j,itt) * wpms(k,itt)
     *                      )
     *                      /sumdiv/sumdiv
               enddo
            enddo
            cpq(i,j) = ac
            cpq(j,i) = ac
         enddo
      enddo
C1)   ===== End of calculating cpq matrix =====

C2)   ===== calculating ppq or d Vector =====
      do i = 1, ntmed + nmp
         ac = 0.0
         do itt = 1, nstep
            do k = 1, nwell
               ac = ac +( (pcal(k,itt) - pms(k,itt) ) * htj (k,i,itt)
     *                   + pstd(k,itt)                * htj1(k,i,itt)
     *                  ) * wpms(k,itt) / sumdiv /sumdiv
            enddo
         enddo
         ppq(i) = 2. * ac
      enddo                                                                     
C2)   ===== End of calculating ppq Vector =====

c                                     
C     ===== Calculate apq1, bpq, and xpq ===============================        
c     Constraints of permrability perturbation                                  
c     In the measured locations for the case of conditional field               
      do i = 1, ntmed                                                           
         ii = 2 * i - 1                                                         
         apq1(ii  ) = -1.0                                                      
         apq1(ii+1) =  1.0                                                      
         dinfe = tms( i ) - am_y * sqrt( vtms(i) )                              
         supe  = tms( i ) + am_y * sqrt( vtms(i) )                              

         bpq (ii  ) =  tppe(i) - dinfe                                          
         bpq (ii+1) = -tppe(i) + supe                                           
         xpq (i   ) = (-bpq(ii+1) + bpq(ii) )/2./apq1(ii)                       
      enddo                                                                     
C     ===== Calculate apq1, bpq, and xpq ===============================        
c     Constraints of permrability perturbation                                  
c     In the Master locations for the case of conditional field                 
c     A small presentage of master points whose values are out                  
c     of the constraint inverval                                                
c      tms(i) =  the ensemble Y mean     at location i                 !        
c     vtms(i) =  the ensemble Y variance at location i                 !        
c     tppe(i) = {Y^(0)}, supe = Y_max, and dinfe = Y_min                        
c     bpq(i  )= {Y^(0)} - Y_min                                                 
c     bpq(i+1)= Y_max   - {Y^(0)}                                               
c     xpq(i)  =-[{Y^(0)} -(Y_min + Y_max)/2 ]                                   
c     equation (25) and (26) of page 170                                        
      do i = ntmed + 1, ntmed + nmp                                             
         ii  = 2 * i - 1                                                        
         apq1(ii    ) = -1.0      ! - [ I ]                                     
         apq1(ii + 1) =  1.0      !   [ I ]                                     
         dinfe = tms( i ) - am_y * sqrt( vtms(i) )                              
         supe  = tms( i ) + am_y * sqrt( vtms(i) )                              

         if    ( tppe(i) .gt. supe  ) then                                      
            supe  = tppe(i)                                                     
         elseif( tppe(i) .lt. dinfe ) then                                      
            dinfe = tppe(i)                                                     
         endif                                                                  
                                                                                
         bpq(ii    ) =  tppe(i) - dinfe  ! P_min                                
         bpq(ii + 1) =  supe - tppe(i)   ! P_max                                
         xpq(i)    = ( -bpq(ii+1) + bpq(ii) )/2./apq1(ii)                       
      enddo                                                                     
                                                                                
      return                                                                    
      end                                                                       
