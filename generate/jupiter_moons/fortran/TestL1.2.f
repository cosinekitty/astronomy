       PROGRAM  TESTL1_2

      !********************************************************************
      !*     Semi-analytic theory of motion of Galilean satellites        *
      !*          fitted to observations between 1891 and 2003            *
      !*                                                                  *
      !*                        (Version 1.2)                             *
      !*                                                                  *
      !*------------------------------------------------------------------*
      !*                   L.Duriez (Luc.Duriez@imcce.fr)                 *
      !*                   V.Lainey (lainey@oma.be)                       *
      !*                   A.Vienne (vienne@imcce.fr)                     *
      !*                                                                  *
      !*    IMCCE - 77, Avenue Denfert-Rochereau 75014 Paris (France)     *
      !*__________________________________________________________________*
      !*                                                                  *
      !*         Example of monitor for computing ephemerides             *
      !*         using the subroutine  DL1_2(ET,nsat,iv,ELEM,XVE) :       *
      !*                                                                  *
      !*                                                                  *
      !*  ET : Julian Day  ie  JD in the Terrestrial Time scale (TT)      *
      !*       ET  must be limited to dates between about                 *
      !*       years 1140 and 2760  (these limits concern accurate        *
      !*       ephemerides, else see comments in the subroutine DL1_2)    *
      !*                                                                  *
      !* nsat: satellite number, from 1 (Io) to  4 (Callisto)             *
      !*                                                                  *
      !*  iv : choice of output :                                         *
      !*       0 = orbital elements referred to the center of Jupiter in  *
      !*           a fixed reference frame close to the jovian equatorial *
      !*           frame for J2000.0  epoch,                              *
      !*       1 = cartesian jovicentric coordinates in the fixed         *
      !*           celestial reference frame : mean equator and mean      *
      !*           equinox for J2000.0 epoch.                             *
      !*                                                                  *
      !*         To compute geocentric ephemerides (and other), add the   *
      !*         Earth-Jupiter vector given by  DE406                     *
      !*                                                                  *
      !* ELEM(6) : orbital elements output :                              *
      !*            a, L, Re(z), Im(z), Re(zeta), Im(zeta)                *
      !*                                                                  *
      !* XVE(6)  : cartesian coordinates output :                         *
      !*                       x, y, z, vx, vy, vz                        *
      !********************************************************************
      !*                          15/02/2006                              *
      !********************************************************************

       implicit real*8 (a-h,o-z)
       parameter(ua=149597870.691d0)     !! UA value, from DE406
       real*8 ELEM(6),XVE(6),ELEM1(6),XVE1(6)
       real*8 axkm(4),fc(6),A0(4),MU(4),AN(4)
       common /axes/ A0,MU

       pi=4*atan(1.D0)
       pi2=2*pi
*
*  reading the series  L1.2 from BisL1.2.dat :
*
       call INITbis
*
*  OR reading the series  L1.2 from GalileanL1.2.dat :
*
*       call INIT

       do i=1,4
         axkm(i)= A0(i)*ua
         AN(i)=DSQRT(MU(i)/A0(i)**3.D0)
       end do

       open(12,file='TestL1.2.res')
       
       write(12,*)'semi-major axes in km '
       write(12,*) axkm
       write(12,*)'mean motions in rad/day '
       write(12,*) AN
       write(12,*)'periods in days '
       write(12,*) (pi2/AN(i),i=1,4)
       write(12,*)


      !* ELEM(6) : orbital elements in order :
      !*            a, L, Re(z), Im(z), Re(zeta), Im(zeta)
      !*            a in AU   while  other elements in rad
      !* XVE(6)  : cartesian coordinates output, in order :
      !*                    x, y, z,  vx, vy, vz
      !*              x, y, z in AU   and   vx, vy, vz  in AU/day

       ET0=2451545.0d0
       pas=12.345d0

       write(12,*)'positions velocities (in AU and AU/day) :'
       iv=1
       do nsat=1,4
         write(12,*)'satellite',nsat
         do k=-2000,2000,400
           ET=ET0+k*pas
           CALL DL1_2(ET,nsat,iv,ELEM,XVE)
           write(12,'(f14.3,6d20.10)') ET,((XVE(i)),i=1,6)
         end do
       end do
       write(12,*)

       write(12,*)'differences on elements (in km) :'
       iv=0
       do nsat=1,4
         write(12,*)'satellite',nsat
         fc(1)=ua
         do i=2,6
          fc(i)=axkm(nsat)
         end do
         do k=-2000,2000,400
           ET=ET0+k*pas
           CALL DL1_2(ET,nsat,iv,ELEM,XVE)
           CALL DL1_2FRC(ET,nsat,0,ELEM1,XVE1)
           write(12,'(i8,6f20.6)') k,((ELEM1(i)-ELEM(i))*fc(i),i=1,6)
         end do
       end do
       close(12)
       stop
       end

       include  'L1.2.f'
