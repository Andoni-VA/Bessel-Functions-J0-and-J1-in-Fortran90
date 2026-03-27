module bessel

        use mcf_tipos

        public:: J0

        contains

        function J0(rho,alpha)

        real(kind=dp) , intent(in) :: rho , alpha
        real(kind=dp) :: J0
        real(kind=dp) , parameter :: pi=acos(-1.0_dp)

        J0 = cos( rho * cos(alpha))

        J0 = J0  / 2.0_dp / pi

        end function J0

end module bessel

module trap

        use mcf_tipos
        use bessel

        public :: trapez

        contains

                function trapez(rho)

                real(kind=dp) , intent(in) :: rho
                real(kind=dp) :: alpha , a , b , h
                real(kind=dp) , parameter :: pi=acos(-1.0_dp)
                integer, parameter :: N=100
                integer :: i
                real(kind=dp) :: trapez

                a=0.0_dp
                b=2*pi

                h = (b-a)/N

                trapez = J0(rho,a)
                trapez = trapez + J0(rho,b)

                do i=2,N

                trapez = 2.0_dp*J0(rho,a+((i-1)*h)) + trapez

                end do

                trapez = trapez * h / 2.0_dp

                end function trapez

end module trap

module strangerthings

        use mcf_tipos
        use bessel
        use trap

        public :: J1

        contains

                function J1(rho)

                real(kind=dp) , intent(in) :: rho
                real(kind=dp) :: J1

                J1 = rho * trapez(rho)

                end function J1

end module strangerthings

program once

        use mcf_tipos
        use bessel
        use trap
        use strangerthings
        use mcf_cuadratura

        real(kind=dp) :: x , y , h , Z , emaitza
        integer , parameter :: N=300
        integer :: i

        !GRAFICAR LAS FUNCIONES J0 y J1

        open(unit=11,file="Bessel0.dat",action="write",status="replace")

        h = 20.0_dp / N

        do i=1,N
        x = -10.0_dp + (h*(i-1))
        y = trapez(x)
        write(unit=11,fmt="(2f16.8)") x , y
        end do

        close(11)

        open(unit=18,file="Bessel1.dat",action="write",status="replace")

        do i=1,N+1

        Z = -10_dp + h*(i-1)

        call romberg(J1,0.0_dp,Z,emaitza,1e-6_dp)

        emaitza = emaitza /Z

        write(unit=18,fmt="(2es16.8)") Z , emaitza

        end do

        close(18)

end program once



