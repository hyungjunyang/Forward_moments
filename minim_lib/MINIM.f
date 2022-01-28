c============ Parameter List ===========================================
c     mxtpe    maximum number of parameters in optimization(mxtmed+mxmp)
c     mxtmed   maximum number of measured ln(k) data
c     mxmp     maximum number of master points
c     nkrig    nkrig = ntmed + nmp
c     it_min   min. number of iterations
c     ifobj    no. of times that the difference of two consecutive obj. 
c              < eps5
c              is a dummy one, can be replaced by nact1
c     idbg     debug level
c     ldbg     the IO parameter of output bebugging contents   
c     fobj0    
c     sumwall   
c     eps1     for activation of constraints of perturbations of K
c     eps3     for compare norm 1
c     eps5     the difference in obj. function in two consecutive itera
c              Parameters used in each iteration of minimization:    
c     cpq      ( or [C]) matrix of problem in optimization
c     ppq      ( or {d}) vector of problem in optimization
c     apq1     ( or [H]) coefficients of constraints to the perturbation of T
c     pbq      
c     xpq      ( or {p}) results from optimizations
C
C     u        stands for s in the paper
C     alf_k    stands for beta
C     alf_c    stands for beta_opt
C     alf_mi   stands for beta_min
C     alf_ma   stands for beta_max
C
C=======================================================================
      subroutine minim2(nkrig,it_min,ifobj,idbg,ldbg,
     +                  mxtpe,fobj0,sumwall,eps1,eps3,eps5,
     +                  cpq,ppq,apq1,bpq,xpq)
      implicit real*8 (a-h,o-z)  
c       
c     =======    Optimization related ===========
c       
      DIMENSION xpq(mxtpe),cpq(mxtpe,mxtpe),
     +          ppq(mxtpe),apq1(2*mxtpe),bpq(2*mxtpe)
      
c          
c     =======    Local Variables     ===========
      integer   ia_1(mxtpe)
      DIMENSION grad(mxtpe),u(mxtpe),hh_t1(mxtpe),uprim(mxtpe)
      real*8    minfob
      real*8    alf_c,alf_mi,alf_ma,alf_k,acum,acum2,acum3,
     +          norma1
c       
      data maxite/9999/,maxit2/1000/
      data minfob/0.005/
c              
      nvar  = nkrig
      fobj  = fobj0 / sumwall

c      write(*,*) 'Obj. before optim: ', fobj * sumwall
       
c     Initial solution (k = iteration counter)
      k = 0   
100   continue
c                                                                                
      if(k .gt. maxit2 .and. fobj .lt. minfob) goto 600                         
      if(k .eq. maxite) goto 600                                                
                                                                                
      ielre = 0                                                                 
                                                                                
C     ===== Step 2 =====================================================        
C     compute {delta g} = 2[C]*{p} + {d}                                        
      do i = 1, nvar                                                            
         acum = 0.0                                                             
         do j = 1, nvar                                                         
            acum = acum + cpq(i,j) * xpq(j)                                     
         enddo                                                                  
         grad(i) = acum * 2.0 + ppq(i)                                          
      enddo                                                                     
150   continue                                                                  
c                                                                               
c     control the active restriction and calculate vector u:                    
c     .. active restriction of first type                                       
c     i  = 1, 3, 5, 7, ...                                                      
c     i+1= 2, 4, 6, 8, ...                                                      
c     ii = 1, 2, 3, 4, ...                                                      
c     apq1 = [H]                                                                
c     apq1(i  ) = -1.                                                           
c     apq1(i+1) =  1.                                                           
c     bpq(i  )= {Y^(0)} - Y_min                                                 
c     bpq(i+1)= Y_max   - {Y^(0)}                                               
c     xpq(i)  =-[{Y^(0)} -(Y_min + Y_max)/2 ]                                   
c     comment: diesac is not used!!!
      j = 1                                                                     
      do i = 1, 2 * nkrig, 2                                                    
         ii  = (i + 1)/2                                                        
         df1 = -apq1(i  ) * xpq(ii) + bpq(i  )                                  
         df2 = -apq1(i+1) * xpq(ii) + bpq(i+1)                                  
         if(    df1 .le. 0.) then                                               
            ia_1(j) = i                                                         
            j = j + 1                                                           
            xpq(ii) = bpq(i  ) / apq1(i  )                                      
         elseif(df2 .le. 0.) then                                               
            ia_1(j) = i + 1                                                     
            j = j + 1                                                           
            xpq(ii) = bpq(i+1) / apq1(i+1)                                      
         endif                                                                  
      enddo                                                                     
      nact1  = j - 1                                                            
c     .. evaluate the vector u                                                  
200   continue                                                                  
C     ===== Step 3 =====================================================        
      do i = 1, nvar                                                            
         u(i) = - grad(i)                                                       
      enddo                                                                     
c     if no active restriction: goto 400                                        
      if(nact1 .eq. 0) goto 400                                                 

c     there are active restrictions                                             
c     initialize {s'} or {u'}                                                   
      do i = 1, nact1                                                           
         uprim(i) = 0.0                                                         
      enddo                                                                     
c     calculate {s'} (or {u'})                                                  
c     {s'} = ([H][H]^{T})^{-1} * [H] * {delta g}                                
      do i = 1, nact1                                                           
         hh_t1(i) = apq1( ia_1(i) )**2.0                                        
      enddo                                                                     
c                                                                               
      do i = 1, nact1                                                           
         i1 = ia_1(i)                                                           
         if(mod(i1,2) .ne. 1) then                                              
            i1 = i1/2                                                           
         else                                                                   
            i1 = (i1+1)/2                                                       
         endif                                                                  
         uprim(i) = ( apq1( ia_1(i) ) * grad(i1) ) / hh_t1(i)                   
      enddo                                                                     
c                                                                               
c     calculate {s} (or {u})                                                    
c     {s} = -f * {delta g} + [H]^{T} * {s'} where f seems to be 1?              
      do i = 1, nact1                                                           
         ii = ia_1(i)                                                           
         if(mod(ii,2) .ne. 1) then                                              
          i1 = ii/2                                                             
         else                                                                   
          i1 = (ii+1)/2                                                         
         endif                                                                  
         u(i1) = apq1(ii) * uprim(i) + u(i1)                                    
      enddo                                                                     
c                                                                               
400   continue                                                                  
C     ===== Step 4 =====================================================        
C     calculate beta_min (alf_mi) and beta_max (alf_ma)                         
      alf_mi = -1.0d+20                                                         
      alf_ma = +1.0d+20                                                         
      m = 1                                                                     
      do 401 i = 1, 2 * nkrig                                                   
         if(m .le. nact1 .and. ia_1(m) .eq. i) then                             
            m = m + 1                                                           
            goto 401                                                            
         endif                                                                  
c        restriction type 1 is not active                                       
         i1 = i                                                                 
         if(mod(i,2) .ne. 1) then                                               
            i1 = i1/2                                                           
         else                                                                   
            i1 = (i1 + 1)/2                                                     
         endif                                                                  
         acum  = apq1(i) * u(i1)                                                
         acum2 = bpq(i) - apq1(i) * xpq(i1) + eps1                              
         if(acum .gt. 0.) then
            alf_ma = dmin1(alf_ma, acum2/acum)
         elseif(acum .lt. 0.) then
            alf_mi = dmax1(alf_mi, acum2/acum)
         endif
401   continue
c     calculate beta_opt (alf_c)
c     acum3 stands for {s}^{T} * {Delta g}
c     acum2 stands for {s}^{T} * [C] * {s}
      acum2 = 0.0
      acum3 = 0.0
      do i = 1, nvar
         acum = 0.0
         do j1 = 1, nvar
            acum = acum + cpq(i,j1) * u(j1)
         enddo
         acum2 = acum2 + acum * u(i)
         acum3 = acum3 + u(i) * grad(i)
      enddo
c
      if( acum3 .eq. 0. .or. acum2 .eq. 0.) then
          alf_c = 0.0
      else
          alf_c = - ( acum3 / ( 2.0 * acum2) )
      endif
c     .. para que la funcion objetivo disminuya se debe cumplir
c     .. alfa <= 2*alf_c
      alf_ma = dmin1(alf_ma, 2. * alf_c)
c     .. comprueba m x. y min.
      if( alf_mi .gt. alf_ma) then
c        .. este problema debe venir originado por tenerse alguna restricci¢n
c        .. activa no considerada. Normalmente ser  al volver del proceso de
c        .. comprobacion de u'
         if( idbg .ge. 3 ) write(ldbg,*) ' * Problem with Alfa ',
     1                     alf_mi, alf_ma, alf_c
         ielre = ielre + 1
         if( ielre .gt. 10) goto 600
         goto 150
      endif

C     ===== select beta (or alf_k, in which alf_c stands for beta_opt)
      if(alf_c .le. alf_mi) then
         alf_k = alf_mi
      elseif(alf_c .ge. alf_ma) then
         alf_k = alf_ma
      else
         alf_k = alf_c
      endif
c                
c     ===== Update solution {p}(i+1) = {p}(i) + beta * {s}
      do i = 1, nvar
         xpq(i) = xpq(i) + alf_k * u(i)
      enddo
c     ===== Calculate the object function
      acum2 = 0.
      do i = 1, nvar
         acum = 0.0
         do j = 1, nvar
            acum = acum + cpq(i,j) * xpq(j)
         enddo
         acum2 = acum2 + xpq(i) * acum + xpq(i) * ppq(i)
      enddo
      fobja = fobj
      fobj = (acum2 + fobj0)/sumwall
c     
      if( abs(fobj - fobja) .lt. eps5) ifobj = ifobj - 1
c      
c     .. puede pasar a comprobacion de optimo
c        
c     calculate the norm of the vector u
      norma1 = 0.0
      do i = 1, nkrig
         norma1 = norma1 + u(i) * u(i)
      enddo
      norma1 = dsqrt(norma1) * dabs(alf_k)
      if( k .gt. it_min .and. norma1 .lt. 10. * eps3 ) goto 500
c         
c     ..pasa a siguiente iteracion
      k = k + 1
      goto 100
C     ===== END of Step 4 ==============================================
500   continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     comprobacion de optimo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     se busca el maximo de las componentes de u'
c       
      uprimx = 0.0
C     ===== Step 5 =====================================================
      if(nact1 .eq. 0) goto 600
      i_umx = 0
      do i = 1, nact1
         if(uprim(i) .gt. uprimx) then
            uprimx = uprim(i)
            i_umx = i
         endif
      enddo
      if(i_umx .eq. 0 .and. norma1 .lt. eps3) goto 600
c     elimina la fila correspondiente a i_umx de la matriz H
      if(i_umx .gt. 0 .and. i_umx .le. nact1) then
c        es del primer tipo. se altera ia_1() y nact1
         do i = i_umx, nact1-1
            ia_1(i) = ia_1(i+1)
         enddo
         nact1 = nact1 - 1
      elseif(i_umx .eq. 0 .and. norma1 .lt. 5.*eps3 .and.
     1       ifobj .le. 0) then
c        ..considera encontrado el optimo
         goto 600
      elseif(norma1 .lt. 5. * eps3 .and. ifobj .le. 0) then
c        ..considera encontrado el optimo
         goto 600 
      else 
c        ..pasa a siguiente iteracion
         k = k + 1 
         goto 100 
      endif 
c     ..vuelve a evaluar el vector u
      goto 200 
C     ====== END of Step 5 =============================================        
600   continue                                                                  
c     ====== Optimia found =============================================        
c     output options:                                                           

      if(idbg .ge. 3) then                                                      
         write(ldbg,*)                                                          
         write(ldbg,*) 'output of optimization:'                                
         write(ldbg,'(a13,i7,5x,a10,e13.5)') ' Iterations:',k,                  
     *   ' f. obj.: ', fobj * sumwall                                           
         write(ldbg,'(a12,i4)') ' Active Restr.: ',  nact1                   
      endif                                                                     
c     write(*,*) 'Obj. after optim: ', fobj * sumwall                           
c                                                                               
      return                                                                    
      end                                                                       
