module inrtype
	INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
	INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
	INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
	INTEGER, PARAMETER :: SP = KIND(1.0D0)
	INTEGER, PARAMETER :: DP = KIND(1.0D0)
	INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
	INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
	INTEGER, PARAMETER :: LGT = KIND(.true.)
	REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
	REAL(SP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_sp
	REAL(SP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_sp
	REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sp
	REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_sp
	REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
	REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_dp
	REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp
	real(SP), parameter :: mp = 1.6726485e-27_sp , E_charge = 1.6021892e-19_sp, mu0 = 4._sp*pi*1.e-7_sp, &
				eps0 =8.854187817e-12_sp, me=9.10938215e-31_sp
end module inrtype

module opes
   use inrtype
   implicit none

	type :: vector
		real(sp) :: r
		real(sp) :: p
		real(sp) :: z
	end type vector

   interface operator(.dot.)
      module procedure dot
   end interface

   interface operator(.cross.)
      module procedure cross
   end interface

   contains
      function dot(vec1,vec2)
         implicit none
         type(vector), intent(in) :: vec1,vec2
         real(sp) dot
         dot = vec1%r*vec2%r + vec1%p*vec2%p + vec1%z*vec2%z
         return
      end function dot

      function cross(vec1,vec2)
         implicit none
         type(vector), intent(in) :: vec1,vec2
         type(vector) cross
         cross%r = (vec1%p*vec2%z - vec1%z*vec2%p)
         cross%p = (vec1%z*vec2%r - vec1%r*vec2%z)
         cross%z = (vec1%r*vec2%p - vec1%p*vec2%r)
         return
      end function cross
end module opes
