c     driver for routine simplx
      Program xsimplx1
      implicit real*8 (a-h,o-z)

      integer n, m, np, mp, m1, m2, m3, nm1m2
c     m = m1 + m2 + m3
      parameter(n  = 4, m  = 2, 
     *          np = 5, mp = 4,
     *          m1 = 0, m2 = 0, m3 = 2, 
     *          nm1m2 = n + m1 + m2)

      integer i, icase, j, jj, jmax, izrov(N), iposv(m)
      dimension a(mp, np), anum(NP)
      character txt(nm1m2)*2, alpha(np)*2
      logical rite
      data txt/'x1','x2','x3','x4'/
      data a/ 0.0, 2.0, 8.0, 0.0,
     *        0.0,-1.0, 0.0, 0.0,
     *        2.0,-6.0, 3.0, 0.0,
     *       -4.0, 1.0,-4.0, 0.0,
     *        0.0, 0.0,-1.0, 0.0/
      do i = 1, mp
         write(*,101) (a(i,j), j = 1, np)
      enddo
101   format(1x,5(f10.2,2x))

      call simplx(a, M, N, Mp, Np, M1, M2, M3, icase, izrov, iposv)
      write(*,*) 
      write(*,*) 'after simplx'
      do i = 1, mp
         write(*,101) (a(i,j), j = 1, np)
      enddo
      write(*,*) 'izrov= ',(izrov(i), i = 1, N)
      write(*,*) 'iposv= ',(iposv(i), i = 1, m)
      
      if(     icase .eq.  1) then
         write(*,*) 'Unbounded objective function'
      else if(icase .eq. -1) then
         write(*,*) 'No solution satisfy constraints given'
      else
         jj = 1
         do i = 1, N
            if(izrov(i) .le. nm1m2) then
               alpha(jj) = txt( izrov(i) )
               jj = jj + 1
            endif
         enddo
         jmax = jj - 1
         write(*,'(/3x,5a10)') ' ',(alpha(jj), jj = 1, jmax)
         do i = 1, m + 1
            if(i .eq. 1) then
               alpha(1) = '  '
               rite = .true.
            elseif( iposv(i-1) .le. nm1m2 ) then
               alpha(1) = txt(iposv(i-1))
               rite = .true.
            else
               rite = .false.
            endif
            if(rite) then
               anum(1) = a(i,1)
               jj = 2
               do j = 2, n + 1
                  if( izrov(j - 1) .le. nm1m2 ) then
                      anum(jj) = a(i,j)
                      jj = jj + 1
                  endif
               enddo
               jmax = jj - 1
               write(*,'(1x,a3,(5f10.2))') alpha(1),
     *                                     (anum(jj), jj = 1, jmax)
            endif
         enddo
      endif
c
c
      write(*,*) 
      write(*,*) "the program ends normally!"
      stop
      end
      include "simplx.f"
      include "simp1.f"
      include "simp2.f"
      include "simp3.f"
