!     !!!!!!!!!!!!!!!!!!!!!!!!! CALDER Flux Average SUBROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine fluxavg(input,output)
  use gemx_com
  use equil

  implicit none 

  ! Input
  real(8), dimension(0:imx,0:jmx,0:kmx) :: input !3D array to be flux averaged
  ! Outputs
  real(8), dimension(0:imx, 0:jmx) :: output ! 2D interpolated array

  ! Local variables
  real(8), dimension(:), allocatable :: phiavg1d, psi1d, psi1d_private,phiavg1d_private
  integer :: gi, xix, yjy, miw, psi_zero, store, k, line_large, line_small
  real(8) :: weightinput,weightinput3D, phiavggi, psival, wmx0, wmx1, psi_private_min
  allocate(phiavg1d(0:101),psi1d(0:101),psi1d_private(0:101),phiavg1d_private(0:101))

  phiavg1d         = 0
  phiavg1d_private = 0
  psi1d            = 0
  psi1d_private    = 0
  
  line_large = 1
  line_small = 0
  psi_zero   = 1
  !input      = 1

  !Computation
  do line = 1, 101
      phiavggi = 0.0d0 
      store = 0
      do gi = 1, num_lines
        if (gindex(gi) == line-1) then
              weightinput = (weight00(gi)*input(iarray(gi), jarray(gi),0) + &
                          weight10(gi)*input(iarray(gi)+1, jarray(gi),0) + &
                          weight01(gi)*input(iarray(gi), jarray(gi)+1,0) + &
                          weight11(gi)*input(iarray(gi)+1, jarray(gi)+1,0))
  
              if (i3D == 0) then                  
                 phiavggi = phiavggi + (weightinput*jacobian(gi))/deno(gi)
        
              else
                 do k=1, kmx
                    weightinput3D = (weight00(gi)*input(iarray(gi), jarray(gi),k) + &
                                   weight10(gi)*input(iarray(gi)+1, jarray(gi),k) + &
                                   weight01(gi)*input(iarray(gi), jarray(gi)+1,k) + &
                                   weight11(gi)*input(iarray(gi)+1, jarray(gi)+1,k))
                    weightinput = weightinput3D + weightinput
                 enddo
                 phiavggi = phiavggi +(weightinput*jacobian(gi))/(deno(gi)*(kmx+1))
              end if

           if (priv(gi) == 0) then
              store = gi
           end if
        end if

        !Remove redundancy from closed loop integration process
        if (phiavggi /= 0) then
           if (gindex(gi) /= line-1) then
              if (i3D == 0) then
                 phiavggi = phiavggi - (weightinput*jacobian(gi-1))/deno(gi-1)
              else
                 phiavggi = phiavggi - (weightinput*jacobian(gi-1))/(deno(gi-1)*(kmx+1))
              end if
              exit
           end if
        end if
     end do
     
     if (store /= 0) then
        phiavg1d(psi_zero) = phiavggi
        psi1d(psi_zero)    = psitab(store)
        psi_zero = psi_zero + 1
     else 
        if (line_small == 0) then
           psi_private_min = psitab(line)
        end if
        phiavg1d_private(line_small) = phiavggi
        psi1d_private(line_small) = psitab(gi)
        line_small = line_small + 1
     end if
  end do
  phiavg1d(0) = phiavg1d(1)
  
  open(unit=11, file = 'testphiavg1d',status='unknown',action='write')
  write(11,*) phiavg1d
  close(11)



  
  ! phiavg1d(0) = input(268,254,0)
  
  !Save psi1d and timesteps of phiavg1d to understand convergence
  ! if (timestep == 10) then
  !    open(unit=11, file = 'psi1d',status='unknown',action='write')
  !                write(11,*) psi1d(:)
  !             close(11)
  ! endif

  ! open(unit=11, file = 'phiavg1d',status='unknown',position='append')                
  ! write(11,*) phiavg1d(:)
  ! close(11)

  !Initialize output to zero
  output(xix,yjy) = 0.0
  
  !INTERPOLATION
  do xix = 0, nx
      do yjy = 0, nz
          psival = psi_p(xix,yjy)
          if (mask(xix,yjy) < 0.99) then 
              output(xix,yjy) = 0
          else
              if (yjy < 75 .and. xix < 150 .and. psival > 0.29 .and. psival<0.31) then !Private region under X-point
                 miw  = int((psival-psi_private_min)/(psi1d(2)-psi1d(1)))
                 wmx0 = ((miw+1)*(psi1d_private(2)-psi1d_private(1))-psival)/(psi1d(2)-psi1d(1))
                 wmx1 = 1.-wmx0
                 output(xix,yjy) = wmx0*phiavg1d_private(miw) + wmx1*phiavg1d_private(miw+1)
              else 
                 miw  = int(psival/(psi1d(2)-psi1d(1)))
                 wmx0 = ((miw+1)*(psi1d(2)-psi1d(1))-psival)/(psi1d(2)-psi1d(1))
                 wmx1 = 1.-wmx0  
                 output(xix,yjy) = wmx0*phiavg1d(miw) + wmx1*phiavg1d(miw+1)
              end if
          end if
      enddo
   enddo
  !  output = output*32
end subroutine fluxavg