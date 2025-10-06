            i=int(xt/dx)
            j=int(zt/dz)
	    k=int(zeta/dzeta)


            wx0=float(i+1)-xt/dx 
            wx1=1.-wx0
            wy0=float(j+1)-zt/dz
            wy1=1.-wy0
	    wz0=float(k+1)-zeta/dzeta
            wz1=1.-wz0

	      k_plus_1=k+1
	      if(k==kmx) k_plus_1=0  
            exp1=exp1 + wx0*wy0*wz0*ex(i,j,k) + wx1*wy0*wz0*ex(i+1,j,k) &
            + wx0*wy1*wz0*ex(i,j+1,k) + wx1*wy1*wz0*ex(i+1,j+1,k) + &
            wx0*wy0*wz1*ex(i,j, k_plus_1) + wx1*wy0*wz1*ex(i+1,j, k_plus_1) + &
            wx0*wy1*wz1*ex(i,j+1, k_plus_1) + wx1*wy1*wz1*ex(i+1,j+1, k_plus_1)

            ezp=ezp + wx0*wy0*wz0*ez(i,j,k) + wx1*wy0*wz0*ez(i+1,j,k) &
            + wx0*wy1*wz0*ez(i,j+1,k) + wx1*wy1*wz0*ez(i+1,j+1,k) + &
            wx0*wy0*wz1*ez(i,j, k_plus_1) + wx1*wy0*wz1*ez(i+1,j, k_plus_1) + &
            wx0*wy1*wz1*ez(i,j+1, k_plus_1) + wx1*wy1*wz1*ez(i+1,j+1, k_plus_1)

            ezetap =ezetap + wx0*wy0*wz0*ezeta(i,j,k) + wx1*wy0*wz0*ezeta(i+1,j,k) &
            + wx0*wy1*wz0*ezeta(i,j+1,k) + wx1*wy1*wz0*ezeta(i+1,j+1,k) + &
            wx0*wy0*wz1*ezeta(i,j, k_plus_1) + wx1*wy0*wz1*ezeta(i+1,j, k_plus_1) + &
            wx0*wy1*wz1*ezeta(i,j+1, k_plus_1) + wx1*wy1*wz1*ezeta(i+1,j+1, k_plus_1)

            delbxp =delbxp + wx0*wy0*wz0*delbx(i,j,k)  &
            + wx1*wy0*wz0*delbx(i+1,j,k) &
            + wx0*wy1*wz0*delbx(i,j+1,k) &
            + wx1*wy1*wz0*delbx(i+1,j+1,k) & 
            + wx0*wy0*wz1*delbx(i,j, k_plus_1) &
            + wx1*wy0*wz1*delbx(i+1,j, k_plus_1) & 
            + wx0*wy1*wz1*delbx(i,j+1, k_plus_1) &
            + wx1*wy1*wz1*delbx(i+1,j+1, k_plus_1)

            delbzp =delbzp + wx0*wy0*wz0*delbz(i,j,k) &
            + wx1*wy0*wz0*delbz(i+1,j,k) &
            + wx0*wy1*wz0*delbz(i,j+1,k)  &
            + wx1*wy1*wz0*delbz(i+1,j+1,k)  &
            + wx0*wy0*wz1*delbz(i,j, k_plus_1)  &
            + wx1*wy0*wz1*delbz(i+1,j, k_plus_1)  &
            + wx0*wy1*wz1*delbz(i,j+1, k_plus_1)  &
            + wx1*wy1*wz1*delbz(i+1,j+1, k_plus_1)


!	     write(*,*) exp1, ezp,ezetap,delbxp,delbzp
			   
