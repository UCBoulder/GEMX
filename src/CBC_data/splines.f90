MODULE splines
!inr模块的功能是实现一维的三次分段多项式样条插值。
!实际运用时，直接调用的一维三次样条插值子程序或函数：inrcsd,inrcsp,inrcsnak,inrcsval,inrcsitg

!SUBROUTINE inrcsd(x,y,dldf,yp1,ypn,coef)
!边界条件取端点的一或二阶导数，求分段多项式的系数。
!输出结果coef(1:4,i)，表示区间(x(i),x(i+1))上插值多项式x^3,x^2,x^1,x^0项的系数。
!输入：实数数组x(n),y(n)；整数dldf(11或12或21或22)；实数yp1,ypn
!输出：实数数组coef(4,n-1)

!SUBROUTINE inrcsp(x,y,coef)
!取周期边界条件，求分段多项式的系数。
!输出结果coef(1:4,i)，表示区间(x(i),x(i+1))上插值多项式x^3,x^2,x^1,x^0项的系数。
!输入：实数数组x(n),y(n)
!输出：实数数组coef(4,n-1)

!SUBROUTINE inrcsnak(x,y,coef)
!边界条件取not-a-knot，即区间(x1,x2)和(x2,x3)上的多项式取为相同，区间(x(n-2),x(n-1))和(x(n-1),x(n))上的多项式相同，求分段多项式的系数。
!输出结果coef(1:4,i)，表示区间(x(i),x(i+1))上插值多项式x^3,x^2,x^1,x^0项的系数。
!输入：实数数组x(n),y(n)
!输出：实数数组coef(4,n-1)

!FUNCTION inrcsval(xa,coef,x,i)
!该函数求解一维样条插值的0，1，2阶导数。依赖inrcsd,inrcsp,inrcsnak求解出系数。
!输入：实数数组xa(n)；实数矩阵coef(4,n-1)；实数x；整数i(0或1或2)
!输出：实数inrcsval

!FUNCTION inrcsitg(x,coef,xa,xb)
!该函数求区间(xa,xb)上的积分值。依赖inrcsd,inrcsp,inrcsnak求解出系数。
!输入：实数数组x(n)；实数矩阵coef(4,n-1)；实数xa,xb
!输出：实数inrcsitg

!cdbbs模块的功能是实现一维和二维的样条基插值。
!实际运用时，直接调用的样条基插值子程序或函数：cdbbsnak,cdbbscoef,cdbbsval,cdbbscoef2,cdbbsval2d

!SUBROUTINE cdbbsnak(x,k,knot)
!功能：对于数组x(n)，采用not-a-knot方法生成k阶节点数组knot(n+k)。
!输入：实数数组x(n)；整数k
!输出：实数数组knot(n+k)

!SUBROUTINE cdbbscoef(x,y,knot,k,bcoef)
!功能：给定自变量数组x(n)，函数值数组y(n)，k阶的节点数组knot(n+k)，求解出样条基插值的系数bcoef(n)。
!节点数组需要事先生成，可用cdbbsnak生成。
!输入：实数数组x(n),y(n),knot(n+k)；整数k
!输出：实数数组bcoef(n)

!FUNCTION cdbbsval(knot,k,bcoef,x,jderiv)
!功能：给出节点序列knot(n+k)，阶数k，系数bcoef(n)，求出点x的jderiv阶导数。
!需要cdbbscoef事先求解出系数。
!输入：实数数组knot(n+k)；整数k；实数数组becoef(n)；实数x；整数jderiv
!输出：实数cdbbsval

!SUBROUTINE cdbbscoef2d(x,y,fxy,knotx,knoty,kx,ky,coef2d)
!功能：给出自变量数组x(nx)，y(ny)，函数值数组fxy(nx,ny)，x方向的节点数组knotx(nx+kx)，y方向的节点数组knoty(ny+ky)，阶数kx，ky，求出样条基的系数coef2d(nx,ny)
!节点数组需要事先生成，可用cdbbsnak生成。
!输入：实数数组x(nx),y(ny)；实数矩阵fxy(nx,ny)；实数数组knotx(nx+kx),knoty(ny+ky)；整数kx，ky
!输出：实数矩阵coef2d(nx,ny)

!FUNCTION cdbbsval2d(knotx,knoty,kx,ky,coef2d,x,y,derivx,derivy)
!功能：给出节点数组knotx(nx+kx)，knoty(ny+ky)，阶数kx，ky，系数coef2d(nx,ny)，求出点(x,y)上x方向derivx阶导y方向derivy阶导的值。
!需要cdbbscoef2d事先求解出系数。
!输入：实数数组knotx(nx+kx),knoty(ny+ky)；整数kx，ky；实数矩阵coef2d(nx,ny)；实数x,y；整数derivx,derivy
!输出：实数cdbbsval2d
	INTERFACE inrtridag
		MODULE PROCEDURE inrtridag_ser
	END INTERFACE
CONTAINS
	SUBROUTINE inrbandec(a,m1,m2,al,indx,d)
	USE inrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
	INTEGER(I4B), INTENT(IN) :: m1,m2
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: al
	INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
	REAL(SP), INTENT(OUT) :: d
	REAL(SP), PARAMETER :: TINY=1.0e-20_sp
	INTEGER(I4B) :: i,k,l,mdum,mm,n
	REAL(SP) :: dum
	REAL(SP) :: swaptemp(m1+m2+1)
	INTEGER(I4B) :: vec(m1),imax(1)
	if (size(a,1) == size(al,1) .and. size(al,1) == size(indx)) then
		n=size(a,1)
	else
		write (*,*) 'error: matrix dimension check error, inrbandec: n'
		STOP 'program terminated by inrbandec'
	end if
	if (size(a,2) == m1+m2+1 ) then
		mm=size(a,2)
	else
		write (*,*) 'error: matrix dimension check error, inrbandec: mm'
		STOP 'program terminated by inrbandec'
	end if
	if (size(al,2) == m1 ) then
		mdum=size(al,2)
	else
		write (*,*) 'error: matrix dimension check error, inrbandec: mdum'
		STOP 'program terminated by inrbandec'
	end if
	do i=1,m1
		vec(i) = m1+1-i
	end do
	a(1:m1,:)=eoshift(a(1:m1,:),dim=2,shift=vec)
	d=1.0
	do k=1,n
		l=min(m1+k,n)
		imax=maxloc(abs(a(k:l,1)))
		i=imax(1)+k-1
		dum=a(i,1)
		if (dum == 0.0) a(k,1)=TINY
		indx(k)=i
		if (i /= k) then
			d=-d
			swaptemp=a(k,1:mm)
			a(k,1:mm)=a(i,1:mm)
			a(i,1:mm)=swaptemp
		end if
		do i=k+1,l
			dum=a(i,1)/a(k,1)
			al(k,i-k)=dum
			a(i,1:mm-1)=a(i,2:mm)-dum*a(k,2:mm)
			a(i,mm)=0.0
		end do
	end do
	END SUBROUTINE inrbandec


	SUBROUTINE inrbanbks(a,m1,m2,al,indx,b)
	USE inrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: a,al
	INTEGER(I4B), INTENT(IN) :: m1,m2
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
	INTEGER(I4B) :: i,k,l,mdum,mm,n
	REAL(SP) :: swapdum
	if (size(a,1) == size(al,1) .and. size(al,1) == size(b) .and. size(b) == size(indx)) then
		n=size(a,1)
	else
		write (*,*) 'error: matrix dimension check error, inrbanbks: n'
		STOP 'program terminated by inrbanbks'
	end if
	if (size(a,2) == m1+m2+1 ) then
		mm=size(a,2)
	else
		write (*,*) 'error: matrix dimension check error, inrbanbks: mm'
		STOP 'program terminated by inrbanbks'
	end if
	if (size(al,2) == m1 ) then
		mdum=size(al,2)
	else
		write (*,*) 'error: matrix dimension check error, inrbanbks: mdum'
		STOP 'program terminated by inrbanbks'
	end if
	do k=1,n
		l=min(n,m1+k)
		i=indx(k)
		if (i /=k) then
			swapdum=b(i)
			b(i)=b(k)
			b(k)=swapdum
		end if
		b(k+1:l)=b(k+1:l)-al(k,1:l-k)*b(k)
	end do
	do i=n,1,-1
		l=min(mm,n-i+1)
		b(i)=(b(i)-dot_product(a(i,2:l),b(1+i:i+l-1)))/a(i,1)
	end do
	END SUBROUTINE inrbanbks


	SUBROUTINE inrtridag_ser(a,b,c,r,u)
	USE inrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,c,r
	REAL(SP), DIMENSION(:), INTENT(OUT) :: u
	REAL(SP), DIMENSION(size(b)) :: gam
	INTEGER(I4B) :: n,j
	REAL(SP) :: bet
        if (size(a)+1 == size(b) .and. size(b) == size(c)+1 .and.  size(c)+1==size(r) .and. size(r)==size(u)) then
		n=size(a)+1
	else
		write (*,*) 'error: dimension check error, inrtridag_ser'
		STOP 'program terminated by inrtridag_ser'
	end if
	bet=b(1)
	if (bet == 0.0) then
		write (*,*) 'error: inrtridag_ser: Error at code stage 1'
		STOP 'program terminated by inrtridag_ser'
	end if
	u(1)=r(1)/bet
	do j=2,n
		gam(j)=c(j-1)/bet
		bet=b(j)-a(j-1)*gam(j)
		if (bet == 0.0) then
			write (*,*) 'error: inrtridag_ser: Error at code stage 2'
			STOP 'program terminated by inrtridag_ser'
		end if
		u(j)=(r(j)-a(j-1)*u(j-1))/bet
	end do
	do j=n-1,1,-1
		u(j)=u(j)-gam(j+1)*u(j+1)
	end do
	END SUBROUTINE inrtridag_ser


	SUBROUTINE inrcyclic(a,b,c,alpha,beta,r,x)
	USE inrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN):: a,b,c,r
	REAL(SP), INTENT(IN) :: alpha,beta
	REAL(SP), DIMENSION(:), INTENT(OUT):: x
	INTEGER(I4B) :: n
	REAL(SP) :: fact,gamma
	REAL(SP), DIMENSION(size(x)) :: bb,u,z
        if (size(a)+1 == size(b) .and. size(b) == size(c)+1 .and. size(c)+1 == size(r) .and. size(r) == size(x)) then
		n=size(a)+1
	else
		write (*,*) 'error: matrix dimension check error, inrcyclic'
		STOP 'program terminated by inrcyclic'
	end if
	if (.not. n>2) then
		write (*,*) 'rerror: matrix dimension is too short, inrcyclic'
		STOP 'program terminated by inrcyclic'
	end if
	gamma=-b(1)
	bb(1)=b(1)-gamma
	bb(n)=b(n)-alpha*beta/gamma
	bb(2:n-1)=b(2:n-1)
	call inrtridag(a,bb,c,r,x)
	u(1)=gamma
	u(n)=alpha
	u(2:n-1)=0.0
	call inrtridag(a,bb,c,u,z)
	fact=(x(1)+beta*x(n)/gamma)/(1.0_sp+z(1)+beta*z(n)/gamma)
	x=x-fact*z
	END SUBROUTINE inrcyclic


	SUBROUTINE inrcsd(x,y,dldf,yp1,ypn,coef)
	USE inrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
	INTEGER(I4B), INTENT(IN) :: dldf
	REAL(SP), INTENT(IN) :: yp1,ypn
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: coef
	INTEGER(I4B) :: n
	REAL(SP), DIMENSION(size(x)) :: a,b,c,r
	REAL(SP), DIMENSION(size(x)) :: y2
	REAL(SP), DIMENSION(size(x)-1) :: tmp1,tmp2,tmp3,a1,a2,a3,a4,a5,a6
        if (size(x) == size(y) .and. size(y) == size(coef,2)+1 .and. size(coef,1) == 4) then
		n=size(x)
	else
		write (*,*) 'error: dimension check error, inrcsd'
		STOP 'program terminated by inrcsd'
	end if
	if (.not. (dldf==11 .or. dldf==12 .or. dldf==21 .or. dldf==22)) then
		write(*,*) 'error: dldf should be 11, or 12, or 21, or 22, inrcsd'
		STOP 'program terminated by inrcsd'
	end if
	c(1:n-1)=x(2:n)-x(1:n-1)
	r(1:n-1)=6.0_sp*((y(2:n)-y(1:n-1))/c(1:n-1))
	r(2:n-1)=r(2:n-1)-r(1:n-2)
	a(2:n-1)=c(1:n-2)
	b(2:n-1)=2.0_sp*(c(2:n-1)+a(2:n-1))
	b(1)=1.0
	b(n)=1.0
	if (dldf==21 .or. dldf==22) then
		r(1)=yp1
		c(1)=0.0
	else
		r(1)=(3.0_sp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
		c(1)=0.5
	end if
	if (dldf==12 .or. dldf==22) then
		r(n)=ypn
		a(n)=0.0
	else
		r(n)=(-3.0_sp/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
		a(n)=0.5
	end if
	call inrtridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
	a1=x(2:n)
	a2=x(1:n-1)
	a3=y2(2:n)
	a4=y2(1:n-1)
	a5=y(2:n)
	a6=y(1:n-1)
	tmp1=a2-a1
	tmp2=a1+a2
	tmp3=a1*a2
	coef(1,:)=(a4-a3)/tmp1/6.0_sp
	coef(2,:)=0.5_sp*(a2*a3-a1*a4)/tmp1
	coef(3,:)=(a6-a5)/tmp1-tmp2*(a4+a3)/3.0_sp+(a2*a4*(tmp2+a1)-a1*a3*(tmp2+a2))/tmp1/6.0_sp
	coef(4,:)=(a2*a5-a1*a6)/tmp1+a1*a2*(a4*(tmp1-a1)+a3*(tmp1+a2))/tmp1/6.0_sp
	END SUBROUTINE inrcsd


	SUBROUTINE inrcsp(x,y,coef)
	USE inrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: coef
	INTEGER(I4B) :: n
	REAL(SP), DIMENSION(size(x)) :: a,b,c,r
	REAL(SP) alpha,beta
	REAL(SP), DIMENSION(size(x)) :: y2
	REAL(SP), DIMENSION(size(x)-1) :: tmp1,tmp2,tmp3,a1,a2,a3,a4,a5,a6
        if (size(x) == size(y) .and. size(y) == size(coef,2)+1 .and. size(coef,1) == 4) then
		n=size(x)
	else
		write (*,*) 'error: dimension check error, inrcsp'
		STOP 'program terminated by inrcsp'
	end if
	c(1:n-1)=x(2:n)-x(1:n-1)
	r(1:n-1)=6.0_sp*((y(2:n)-y(1:n-1))/c(1:n-1))
	r(2:n-1)=r(2:n-1)-r(1:n-2)
	a(2:n-1)=c(1:n-2)
	b(2:n-1)=2.0_sp*(c(2:n-1)+a(2:n-1))
	alpha=c(n-1)
	beta=alpha
	b(1)=2.0_sp*(c(1)+c(n-1))
	r(1)=6.0_sp*((y(2)-y(1))/c(1)-(y(n)-y(n-1))/c(n-1))
	call inrcyclic(a(2:n-1),b(1:n-1),c(1:n-2),alpha,beta,r(1:n-1),y2(1:n-1))
	y2(n)=y2(1)
	a1=x(2:n)
	a2=x(1:n-1)
	a3=y2(2:n)
	a4=y2(1:n-1)
	a5=y(2:n)
	a6=y(1:n-1)
	tmp1=a2-a1
	tmp2=a1+a2
	tmp3=a1*a2
	coef(1,:)=(a4-a3)/tmp1/6.0_sp
	coef(2,:)=0.5_sp*(a2*a3-a1*a4)/tmp1
	coef(3,:)=(a6-a5)/tmp1-tmp2*(a4+a3)/3.0_sp+(a2*a4*(tmp2+a1)-a1*a3*(tmp2+a2))/tmp1/6.0_sp
	coef(4,:)=(a2*a5-a1*a6)/tmp1+a1*a2*(a4*(tmp1-a1)+a3*(tmp1+a2))/tmp1/6.0_sp
	END SUBROUTINE inrcsp


	SUBROUTINE inrcsnak(x,y,coef)
	USE inrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: coef
	INTEGER(I4B) :: n
	REAL(SP), DIMENSION(size(x)-2) :: xx,yy
	REAL(SP), DIMENSION(size(x)-2) :: a,b,c,r
	REAL(SP), DIMENSION(size(x)-2) :: y2
	REAL(SP), DIMENSION(size(x)-3) :: tmp1,tmp2,tmp3,a1,a2,a3,a4,a5,a6
	REAL(SP) :: x32,x31,x21
	REAL(SP), DIMENSION(4, size(x)-3) :: coef2
        if (size(x) == size(y) .and. size(y) == size(coef,2)+1 .and. size(coef,1) == 4) then
		n=size(x)
	else
		write (*,*) 'error: dimension check error, inrcsnak'
		STOP 'program terminated by inrcsnak'
	end if
	xx(1)=x(1)
	xx(2:n-3)=x(3:n-2)
	xx(n-2)=x(n)
	yy(1)=y(1)
	yy(2:n-3)=y(3:n-2)
	yy(n-2)=y(n)
	c(1:n-3)=xx(2:n-2)-xx(1:n-3)
	r(1:n-3)=6.0_sp*((yy(2:n-2)-yy(1:n-3))/c(1:n-3))
	r(2:n-3)=r(2:n-3)-r(1:n-4)
	a(2:n-3)=c(1:n-4)
	b(2:n-3)=2.0_sp*(c(2:n-3)+a(2:n-3))
	x32=x(3)-x(2)
	x21=x(2)-x(1)
	x31=x(3)-x(1)
	b(1)=x32*((x32-x31)*(x32+x31))/6.0_sp
	c(1)=x21*((x21+x31)*(x21-x31))/6.0_sp
	r(1)=y(2)*x31-y(1)*x32-y(3)*x21
	x32=x(n)-x(n-1)
	x21=x(n-1)-x(n-2)
	x31=x(n)-x(n-2)
	a(n-2)=x32*((x32-x31)*(x32+x31))/6.0_sp
	b(n-2)=x21*((x21+x31)*(x21-x31))/6.0_sp
	r(n-2)=y(n-1)*x31-y(n-2)*x32-y(n)*x21
	call inrtridag(a(2:n-2),b(1:n-2),c(1:n-3),r(1:n-2),y2(1:n-2))
	a1=xx(2:n-2)
	a2=xx(1:n-3)
	a3=y2(2:n-2)
	a4=y2(1:n-3)
	a5=yy(2:n-2)
	a6=yy(1:n-3)
	tmp1=a2-a1
	tmp2=a1+a2
	tmp3=a1*a2
	coef2(1,:)=(a4-a3)/tmp1/6.0_sp
	coef2(2,:)=0.5_sp*(a2*a3-a1*a4)/tmp1
	coef2(3,:)=(a6-a5)/tmp1-tmp2*(a4+a3)/3.0_sp+(a2*a4*(tmp2+a1)-a1*a3*(tmp2+a2))/tmp1/6.0_sp
	coef2(4,:)=(a2*a5-a1*a6)/tmp1+a1*a2*(a4*(tmp1-a1)+a3*(tmp1+a2))/tmp1/6.0_sp
	coef(:,1)=coef2(:,1)
	coef(:,2)=coef2(:,1)
	coef(:,3:n-3)=coef2(:,2:n-4)
	coef(:,n-2)=coef2(:,n-3)
	coef(:,n-1)=coef2(:,n-3)
	END SUBROUTINE inrcsnak


	FUNCTION inrlocate(xx,x)
	USE inrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: xx
	REAL(SP), INTENT(IN) :: x
	INTEGER(I4B) :: inrlocate
	INTEGER(I4B) :: n,jl,jm,ju
	LOGICAL :: ascnd
	n=size(xx)
	ascnd = (xx(n) >= xx(1))
	jl=0
	ju=n+1
	do
		if (ju-jl <= 1) exit
		jm=(ju+jl)/2
		if (ascnd .eqv. (x >= xx(jm))) then
			jl=jm
		else
			ju=jm
		end if
	end do
	if (x == xx(1)) then
		inrlocate=1
	else if (x == xx(n)) then
		inrlocate=n-1
	else
		inrlocate=jl
	end if
	END FUNCTION inrlocate


	FUNCTION inrcsval(xa,coef,x,i)
	USE inrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: xa
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: coef
	REAL(SP), INTENT(IN) :: x
	INTEGER(I4B), INTENT(IN) :: i
	REAL(SP) :: inrcsval
	INTEGER(I4B) :: khi,klo,n
	REAL(SP) :: h
        if (size(xa) == size(coef,2)+1 .and. size(coef,1) == 4) then
		n=size(xa)
	else
		write (*,*) 'error: dimension check error, inrcsval'
		STOP 'program terminated by inrcsval'
	end if
	klo=max(min(inrlocate(xa,x),n-1),1)
	khi=klo+1
	h=xa(khi)-xa(klo)
	if (h == 0.0) then
		write (*,*) 'error: bad xa input in inrcsval'
		STOP 'program terminated by inrcsval'
	end if
	if (i==0) then
		inrcsval=((coef(1,klo)*x+coef(2,klo))*x+coef(3,klo))*x+coef(4,klo)
	else if (i==1) then
		inrcsval=(3.0_sp*coef(1,klo)*x+2.0_sp*coef(2,klo))*x+coef(3,klo)
	else if (i==2) then
		inrcsval=6.0_sp*coef(1,klo)*x+2.0_sp*coef(2,klo)
	else
		write (*,*) 'error: i should be 0, or 1, or 2.'
		STOP 'program terminated by inrcsval'
	end if
	END FUNCTION inrcsval

	FUNCTION inrcsitg(x,coef,xa,xb)
	USE inrtype
	implicit none
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: coef
	REAL(SP), INTENT(IN) :: xa,xb
	REAL(SP) :: inrcsitg
	INTEGER(I4B) :: n,kloa,khia,klob,khib,i
        if (size(x) == size(coef,2)+1 .and. size(coef,1) == 4) then
		n=size(x)
	else
		write (*,*) 'error: dimension check error, inrcsitg'
		STOP 'program terminated by inrcsitg'
	end if
	kloa=max(min(inrlocate(x,xa),n-1),1)
	khia=kloa+1
	klob=max(min(inrlocate(x,xb),n-1),1)
	khib=klob+1
	if (kloa==klob) then
		inrcsitg = (((0.25_sp*coef(1,kloa)*xb+coef(2,kloa)/3.0_sp)*xb+0.5_sp*coef(3,kloa))*xb+coef(4,kloa))*xb - &
			   (((0.25_sp*coef(1,kloa)*xa+coef(2,kloa)/3.0_sp)*xa+0.5_sp*coef(3,kloa))*xa+coef(4,kloa))*xa
	else
		inrcsitg = (((0.25_sp*coef(1,kloa)*x(khia)+coef(2,kloa)/3.0_sp)*x(khia)+0.5_sp*coef(3,kloa))*x(khia)+coef(4,kloa))*x(khia) - &
			   (((0.25_sp*coef(1,kloa)*xa+coef(2,kloa)/3.0_sp)*xa+0.5_sp*coef(3,kloa))*xa+coef(4,kloa))*xa
		inrcsitg = inrcsitg + &
			   (((0.25_sp*coef(1,klob)*xb+coef(2,klob)/3.0_sp)*xb+0.5_sp*coef(3,klob))*xb+coef(4,klob))*xb - &
			   (((0.25_sp*coef(1,klob)*x(klob)+coef(2,klob)/3.0_sp)*x(klob)+0.5_sp*coef(3,klob))*x(klob)+coef(4,klob))*x(klob)
		do i=khia,klob-1
			inrcsitg = inrcsitg + &
				   (((0.25_sp*coef(1,i)*x(i+1)+coef(2,i)/3.0_sp)*x(i+1)+0.5_sp*coef(3,i))*x(i+1)+coef(4,i))*x(i+1) - &
				   (((0.25_sp*coef(1,i)*x(i)+coef(2,i)/3.0_sp)*x(i)+0.5_sp*coef(3,i))*x(i)+coef(4,i))*x(i)
		end do
	end if
	end FUNCTION inrcsitg


	FUNCTION cdblocate(xx,x)
	USE inrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: xx
	REAL(SP), INTENT(IN) :: x
	INTEGER(I4B) :: cdblocate
	INTEGER(I4B) :: n,jl,jm,ju,i
	LOGICAL :: ascnd
	n=size(xx)
	ascnd = (xx(n) >= xx(1))
	jl=0
	ju=n+1
	do
		if (ju-jl <= 1) exit
		jm=(ju+jl)/2
		if (ascnd .eqv. (x >= xx(jm))) then
			jl=jm
		else
			ju=jm
		end if
	end do
	if (x == xx(1)) then
		do i=2,n
			if(xx(i)==x) then
				cycle
			else
				cdblocate=i-1
				exit
			end if
		end do
	else if (x == xx(n)) then
		do i=n-1,1,-1
			if(xx(i)==x) then
				cycle
			else
				cdblocate=i
				exit
			end if
		end do
	else
		cdblocate=jl
	end if
	END FUNCTION cdblocate


	SUBROUTINE cdbbsnak(x,k,knot)
	USE inrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: X
	INTEGER(I4B), INTENT(IN) :: K
	REAL(SP), DIMENSION(:), INTENT(OUT) :: knot	
	INTEGER(I4B) :: i,n
	if (size(x)+k == size(knot)) then
		n=size(x)
	else
		write (*,*) 'error: dimension check error, cdbbsnak'
		STOP 'program terminated by cdbbsnak'
	end if
	!!!knot(1:k)=x(1)
	knot(1:k)=x(1)-0.01_sp*(x(2)-x(1))
	!!!knot(n+1:n+k)=x(n)
	knot(n+1:n+k)=x(n)+0.01_sp*(x(n)-x(n-1))
	if (mod(k,2)==0) then
		do i=k+1,n
			knot(i)=x(i-k/2)
		end do
	else
		do i=k+1,n
			knot(i)=0.5_sp*(x(i-(k-1)/2)+x(i-1-(k-1)/2))
		end do
	end if
	END SUBROUTINE cdbbsnak


	SUBROUTINE cdbbsvalvec(knot,k,x,left,val)
	USE inrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: knot
	INTEGER(I4B), INTENT(IN) :: k, left
	REAL(SP), INTENT(IN) :: x
	REAL(SP), DIMENSION(:), INTENT(OUT) :: val
	INTEGER(I4B) n1,n2,i,j
	REAL(SP) deltar(k),deltal(k),term,temp

	if (size(val) /= k*(k+1)/2) then
		write (*,*) 'error: dimension of val is wrong, cdbbsvalvec'
		STOP 'program terminated by cdbbsvalvec'
	end if
	if (knot(left) >= knot(left+1)) then
		write(*,*) 'error: knot(left) < knot(left+1) is not satisfied, cdbbsvalvec'
		STOP 'program terminated by cdbbsvalvec'
	end if
	val(1)=1.0_sp
	n1=0
	do j=1,k-1
		deltar(j) = knot(left+j) - x
		deltal(j) = x - knot(left+1-j)
		temp=0.0_sp
		n2=n1+j
		do i = 1, j
			term = val(n1+i) / ( deltar(i) + deltal(j+1-i) )
			val(n2+i) = temp + deltar(i) * term
			temp = deltal(j+1-i) * term
		end do
		val(n2+j+1) = temp
		n1=n2
	end do
	END SUBROUTINE cdbbsvalvec


	SUBROUTINE cdbbscoef(x,y,knot,k,bcoef)
	USE inrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
	REAL(SP), DIMENSION(:), INTENT(IN) :: knot
	INTEGER(I4B), INTENT(IN) :: k
	REAL(SP), DIMENSION(:), INTENT(OUT) :: bcoef	
	INTEGER(I4B) :: i,n,left,right,n0
	INTEGER(I4B), DIMENSION(size(x)) :: indx
	REAL(SP) :: val(k*(k+1)/2),b,a(size(x),2*k-1), al(size(x),k-1)	
	if (size(x)==size(y) .and. size(y)==size(bcoef) .and. size(bcoef)+k==size(knot) ) then
		n=size(x)
	else
		write (*,*) 'error: dimension check error, cdbbscoef'
		STOP 'program terminated by cdbbscoef'
	end if
	left=k
	n0=k*(k-1)/2
	a=0.0_sp
	do i=1,n
		left=max(left,i)
		right = min ( i + k, n + 1 )
		if(x(i) < knot(left)) then
			write(*,*) 'error, x(i) >= knot(i) is not satisfied, i=',i
			STOP 'program terminated by cdbbscoef'
		end if		
		do while ( knot(left+1) <= x(i) )
			left = left + 1
			if ( left < right ) then
				cycle
			end if
			if ( knot(left) < x(i) ) then
      				write(*,*) 'error, x(i) < knot(i+k) is not satisfied, i=',i
				STOP 'program terminated by cdbbscoef'
      			end if
      			left = left - 1
			exit
		end do
		call cdbbsvalvec(knot,k,x(i),left,val)
		a(i,left-i+1:left-i+k)=val(n0+1:n0+k)
	end do
	call inrbandec(a,k-1,k-1,al,indx,b)
	bcoef=y
	call inrbanbks(a,k-1,k-1,al,indx,bcoef)
	END SUBROUTINE cdbbscoef


	FUNCTION cdbbsval(knot,k,bcoef,x,jderiv)
	USE inrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: knot, bcoef
	INTEGER(I4B), INTENT(IN) :: k
	REAL(SP), INTENT(IN) :: x
	INTEGER(I4B), INTENT(IN) :: jderiv
	REAL(SP) :: cdbbsval
	INTEGER(I4B) :: n,i,ilo,j,jc,jcmax,jcmin,jj
	real(SP) :: aj(k), dl(k), dr(k)
	if (size(knot) == size(bcoef)+k) then
		n=size(bcoef)
	else
		write(*,*) 'error: size of knot and bcoef is wrong. cdbbsval'
		STOP 'program terminated by cdbbsval.'
	end if
	cdbbsval = 0.0_sp
	if ( k <= jderiv ) then
		return
	end if
	i=cdblocate(knot,x)
	if ( k <= 1 ) then
		cdbbsval = bcoef(i)
		return
	end if
	jcmin = 1
	if ( k <= i ) then
		do j = 1, k-1
			dl(j) = x - knot(i+1-j)
		end do
	else
		jcmin = 1 - ( i - k )
		do j = 1, i
			dl(j) = x - knot(i+1-j)
		end do
		do j = i, k-1
			aj(k-j) = 0.0_sp
			dl(j) = dl(i)
		end do
	end if
	jcmax = k
	if ( n < i ) then
		jcmax = k + n - i
		do j = 1, k + n - i
			dr(j) = knot(i+j) - x
		end do
		do j = k+n-i, k-1
			aj(j+1) = 0.0_sp
			dr(j) = dr(k+n-i)
		end do
	else
		do j = 1, k-1
			dr(j) = knot(i+j) - x
		end do
	end if
	do jc = jcmin, jcmax
		aj(jc) = bcoef(i-k+jc)
	end do
	do j = 1, jderiv
		ilo = k - j
		do jj = 1, k - j
			aj(jj) = ( ( aj(jj+1) - aj(jj) ) / ( dl(ilo) + dr(jj) ) ) * real ( k - j, kind = sp )
			ilo = ilo - 1
		end do

	end do
	do j = jderiv+1, k-1
		ilo = k-j
		do jj = 1, k-j
			aj(jj) = ( aj(jj+1) * dl(ilo) + aj(jj) * dr(jj) ) / ( dl(ilo) + dr(jj) )
			ilo = ilo - 1
		end do
	end do
	cdbbsval = aj(1)
	END FUNCTION cdbbsval


	SUBROUTINE cdbbscoef2d(x,y,fxy,knotx,knoty,kx,ky,coef2d)
	USE inrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,knotx,knoty
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: fxy
	INTEGER(I4B), INTENT(IN) :: kx,ky
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: coef2d
	INTEGER(I4B) :: nx,ny,i,j
	REAL(SP) :: coef(size(x)),coefy(size(y)),temp(size(x),size(y))
	INTEGER(I4B) :: left,n0,right,indxx(size(x)),indxy(size(y))
	REAL(SP) :: b,alx(size(x),kx-1),aly(size(y),ky-1), valx(kx*(kx+1)/2),valy(ky*(ky+1)/2), &
		    ax(size(x),2*kx-1),ay(size(y),2*ky-1)
	nx=size(x)
	ny=size(y)
	if(size(fxy,1)/=nx .or. size(fxy,2)/=ny .or. size(coef2d,1)/=nx .or. size(coef2d,2)/=ny &
	   .or. nx+kx/=size(knotx) .or. ny+ky/=size(knoty)) then
		write(*,*) 'error: check the size of input and output, cdbbscoef2d'
		STOP 'program terminated by cdbbscoef2d'
	end if
	left=kx
	n0=kx*(kx-1)/2
	ax=0.0_sp
	do i=1,nx
		left=max(left,i)
		right = min ( i + kx, nx + 1 )
		if(x(i) < knotx(left)) then
			write(*,*) 'error, x(i) >= knotx(i) is not satisfied, i=',i
			STOP 'program terminated by cdbbscoef2d'
		end if
		do while ( knotx(left+1) <= x(i) )
			left = left + 1
			if ( left < right ) then
				cycle
			end if
			if ( knotx(left) < x(i) ) then
      				write(*,*) 'error, x(i) < knotx(i+kx) is not satisfied, i=',i
				STOP 'program terminated by cdbbscoef2d'
      			end if
      			left = left - 1
			exit
		end do
		call cdbbsvalvec(knotx,kx,x(i),left,valx)
		ax(i,left-i+1:left-i+kx)=valx(n0+1:n0+kx)
	end do
	call inrbandec(ax,kx-1,kx-1,alx,indxx,b)
	do i=1,ny
		coef=fxy(1:nx,i)
		call inrbanbks(ax,kx-1,kx-1,alx,indxx,coef)
		temp(1:nx,i)=coef(1:nx)
	end do
	left=ky
	n0=ky*(ky-1)/2
	ay=0.0_sp
	do i=1,ny
		left=max(left,i)
		right = min ( i + ky, ny + 1 )
		if(y(i) < knoty(left)) then
			write(*,*) 'error, y(i) >= knoty(i) is not satisfied, i=',i
			STOP 'program terminated by cdbbscoef2d'
		end if
		do while ( knoty(left+1) <= y(i) )
			left = left + 1
			if ( left < right ) then
				cycle
			end if
			if ( knoty(left) < y(i) ) then
      				write(*,*) 'error, y(i) < knoty(i+ky) is not satisfied, i=',i
				STOP 'program terminated by cdbbscoef2d'
      			end if
      			left = left - 1
			exit
		end do
		call cdbbsvalvec(knoty,ky,y(i),left,valy)
		ay(i,left-i+1:left-i+ky)=valy(n0+1:n0+ky)
	end do
	call inrbandec(ay,ky-1,ky-1,aly,indxy,b)
	do i=1,nx
		coefy=temp(i,1:ny)
		call inrbanbks(ay,ky-1,ky-1,aly,indxy,coefy)
		coef2d(i,1:ny)=coefy(1:ny)
	end do
	END SUBROUTINE cdbbscoef2d		


	FUNCTION cdbbsval2d(knotx,knoty,kx,ky,coef2d,x,y,derivx,derivy)
	USE inrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: knotx,knoty
	INTEGER(I4B), INTENT(IN) :: kx,ky
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: coef2d
	REAL(SP), INTENT(IN) :: x,y
	INTEGER(I4B), INTENT(IN) :: derivx,derivy
	REAL(SP) :: cdbbsval2d
	INTEGER(I4B) :: nx,ny,i,lefty
	REAL(SP) :: temp(ky)
	if(size(knotx)==kx+size(coef2d,1) .and. size(knoty)==ky+size(coef2d,2)) then
		nx=size(coef2d,1)
		ny=size(coef2d,2)
	else
		write(*,*) 'error: check the size of input and output, cdbbsval2d.'
		STOP 'program terminated by cdbbsval2d'		
	end if
	lefty=cdblocate(knoty,y)
	do i=1,ky
		temp(i)=cdbbsval(knotx,kx,coef2d(1:nx,lefty-ky+i),x,derivx)
	end do
	cdbbsval2d = cdbbsval(knoty(lefty-ky+1:lefty+ky),ky,temp,y,derivy)
	END FUNCTION cdbbsval2d

	FUNCTION cdbppval(break, coef, k, x, jderiv)
	USE inrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: break
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: coef
	INTEGER(I4B), INTENT(IN) :: k, jderiv
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: cdbppval
	INTEGER(I4B) :: i,n,itmp,m
	REAL(SP) :: h
	if (size(break)==size(coef,2)+1 .and. size(coef,1)==k) then
		n=size(break)
	else
		write (*,*) 'error: matrix dimension check error, cdbppval'
		STOP 'program terminated by cdbppval'
	end if
	cdbppval=0.0_sp
	itmp=k-jderiv
	if (itmp <= 0) return
	i=max(min(inrlocate(break,x),n-1),1)
	h=x-break(i)
	do m=k,jderiv+1,-1
		cdbppval = (cdbppval/real(itmp,kind=sp))*h + coef(m,i)
		itmp = itmp-1
	end do
	END FUNCTION cdbppval

	SUBROUTINE cdbbs2pp(t, bcoef, k, break, coef, l)
	USE inrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: t, bcoef
	INTEGER(I4B), INTENT(IN) :: k
	REAL(SP), DIMENSION(:), INTENT(OUT) :: break
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: coef
	INTEGER(I4B), INTENT(OUT) :: l
	INTEGER(I4B) :: lsofar,left,n,i,j,jp1
	REAL(SP) :: scrtch(k,k),diff,val(k*(k+1)/2)
	n=size(bcoef)
	if (.not.(size(t)==n+k .and. size(break)==n-k+2 .and. size(coef,1)==k .and. size(coef,2)==n-k+1)) then
		write (*,*) 'error: matrix dimension check error, cdbbs2pp'
		STOP 'program terminated by cdbbs2pp'
	end if
	lsofar = 0
	break(1) = t(k)
	do left = k, n
		if ( t(left+1) == t(left) ) cycle
		lsofar = lsofar + 1
		break(lsofar+1) = t(left+1)
		if ( k <= 1 ) then
			coef(1,lsofar) = bcoef(left)
			cycle
		end if
		do i = 1, k
			scrtch(i,1) = bcoef(left-k+i)
		end do
		do jp1 = 2, k
			j = jp1 - 1
			do i = 1, k - j
				diff = t(left+i) - t(left+i-(k-j))
				if ( diff > 0.0_sp ) then
					scrtch(i,jp1) = ((scrtch(i+1,j)-scrtch(i,j))/diff)*real(k-j, kind=sp)
				end if
			end do
		end do
		call cdbbsvalvec(t,k,t(left),left,val)
		coef(k,lsofar) = scrtch(1,k)
		do jp1 = 2, k
			coef(k+1-jp1,lsofar) = dot_product(val(jp1*(jp1-1)/2+1 : jp1*(jp1+1)/2), scrtch(1:jp1,k+1-jp1) ) 
		end do
	end do
	l = lsofar
	if (l/=n-k+1) write(*,*) 'warning: knot sequnce t(k:n) is not monotonic, by cdbbs2pp.'
	END SUBROUTINE cdbbs2pp

	SUBROUTINE cdbbs2pp2d0(knotx,knoty,bscoef2d,kx,ky,breakx,breaky,ppcoef2d,lx,ly)
	USE inrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: knotx,knoty
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: bscoef2d
	INTEGER(I4B), INTENT(IN) :: kx,ky
	REAL(SP), DIMENSION(:), INTENT(OUT) :: breakx,breaky
	REAL(SP), DIMENSION(:,:,:,:), INTENT(OUT) :: ppcoef2d
	INTEGER(I4B), INTENT(OUT) :: lx,ly
	INTEGER(I4B) :: nx,ny,lxmax,lymax,i,j
	REAL(SP), ALLOCATABLE :: ppcoefx(:,:),ppcoefy(:,:),temp(:,:,:)
	nx=size(bscoef2d,1)
	ny=size(bscoef2d,2)
	lxmax=nx-kx+1
	lymax=ny-ky+1
	if (.not.(size(knotx)==nx+kx .and. size(knoty)==ny+ky .and. size(breakx)==lxmax+1 .and. size(breaky)==lymax+1 &
		& .and. size(ppcoef2d,1)==kx .and. size(ppcoef2d,2)==lxmax &
		& .and. size(ppcoef2d,3)==ky .and. size(ppcoef2d,4)==lymax)) then
		write (*,*) 'error: matrix dimension check error, cdbbs2pp2d0'
		STOP 'program terminated by cdbbs2pp2d0'
	end if
	allocate(ppcoefx(kx,lxmax),ppcoefy(ky,lymax),temp(kx,lxmax,ny))
	do i=1,ny
		call cdbbs2pp(knotx,bscoef2d(:,i),kx,breakx,ppcoefx,lx)
		temp(1:kx,1:lx,i)=ppcoefx(1:kx,1:lx)
	end do
	do i=1,kx
		do j=1,lx
			call cdbbs2pp(knoty,temp(i,j,:),ky,breaky,ppcoefy,ly)
			ppcoef2d(i,j,1:ky,1:ly)=ppcoefy(1:ky,1:ly)
		end do
	end do
	deallocate(ppcoefx,ppcoefy,temp)
	END SUBROUTINE cdbbs2pp2d0

	FUNCTION cdbppval2d(breakx,breaky,kx,ky,ppcoef2d,x,y,derivx,derivy)
	USE inrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: breakx,breaky
	INTEGER(I4B), INTENT(IN) :: kx,ky
	REAL(SP), DIMENSION(:,:,:,:), INTENT(IN) :: ppcoef2d
	REAL(SP), INTENT(IN) :: x,y
	INTEGER(I4B), INTENT(IN) :: derivx,derivy
	REAL(SP) :: cdbppval2d
	INTEGER(I4B) :: lx,ly,i,lefty
	REAL(SP) :: temp(ky,1)
	if(size(ppcoef2d,1)==kx .and. size(ppcoef2d,2)+1==size(breakx) .and. &
		& size(ppcoef2d,3)==ky .and. size(ppcoef2d,4)+1==size(breaky)) then
		lx=size(breakx)-1
		ly=size(breaky)-1
	else
		write(*,*) 'error: check the size of input and output, cdbppval2d.'
		STOP 'program terminated by cdbppval2d'
	end if
	lefty=max(min(inrlocate(breaky,y),ly),1)
	do i=1,ky
		temp(i,1)=cdbppval(breakx,ppcoef2d(:,:,i,lefty),kx,x,derivx)
	end do
	cdbppval2d = cdbppval(breaky(lefty:lefty+1),temp,ky,y,derivy)
	END FUNCTION cdbppval2d
END MODULE splines
