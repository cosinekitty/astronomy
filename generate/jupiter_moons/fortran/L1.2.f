       subroutine DL1_2(ET,ks,iv,ELEM,XVE)

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
      !*  ET : Julian Day  ie  JD in the Terrestrial Time scale (TT)      *
      !*       ET  must be limited to dates between about                 *
      !*             years 1140 and 2760                                  *
      !*  ks : satellite number, from 1 (Io) to  4 (Callisto)             *
      !*  iv : choice of output :                                         *
      !*       0 = orbital elements referred to the center of Jupiter in  *
      !*           a fixed reference frame close to the jovian equatorial *
      !*           frame for J2000.0  epoch,                              *
      !*       1 = cartesian jovicentric coordinates in the fixed         *
      !*           celestial reference frame : mean equator and mean      *
      !*           equinox for J2000.0 epoch.                             *
      !* ELEM(6) : orbital elements output, in order :                    *
      !*            a, L, Re(z), Im(z), Re(zeta), Im(zeta)                *
      !*            a in AU   while  other elements in rad                *
      !* XVE(6)  : cartesian coordinates output, in order :               *
      !*                    x, y, z,  vx, vy, vz                          *
      !*              x, y, z in AU   and   vx, vy, vz  in AU/day         *
      !********************************************************************
      !*                          15/02/2006                              *
      !********************************************************************

       implicit real*8 (a-h,o-z)
       real*8 ampl(100,4,4),phas(100,4,4),freq(100,4,4),al(2,4)
       real*8 elem(6),MU(4),XV(6),XVE(6),c(9,5,4),val(5),TN(0:8)
       integer nbterm(4,4)
       logical AccurateEphem
       common /serL1/nbterm,T0,T1,T2,al,ampl,phas,freq,c,MU,ome,ainc

       pi=4.d0*datan(1.d0)
       pi2=pi*2.d0
       AccurateEphem=.true.

       T=ET-T0  ! (T0=2433282.5d0 => JD from 01/01/1950 at 0H00 (TT)
       val=0.D0
*
*  computing corrections to elements from Tchebycheff polynomials :
*  These corrections represent very long period perturbations which
*  are not computed in form of trigonometric terms but only approximated
*  by Tchebycheff polynomials, hence limiting the interval of validity
*  of ephemerides to about [-810,810] years centered on J1950.
*
*  To compute approximate ephemerides at longer time, skip
*  these corrections by giving a .false. value to AccurateEphem
*
       if( AccurateEphem ) then
        a= T1    !  -819.727638594856D0
        b= T2    !   812.721806990360D0
        x=(T/365.25D0-0.5D0*(b+a))/(0.5D0*(b-a))
	if(abs(x).gt.1.d0) then
	  write(*,*)'date',T,' out of interval of validity'
	  return
	end if
        TN(0)=1.D0
        TN(1)=x
        do it=2,8
          TN(it)=2.D0*x*TN(it-1)-TN(it-2)
        end do
        do nv=1,5
          do it=0,8
          val(nv)=val(nv)+c(it+1,nv,ks)*TN(it)
          end do
          val(nv)=val(nv)-0.5D0*c(1,nv,ks)
        end do
       end if
*
*  computing elements from series
*
        kv=1
         s=0.d0
        do k=1,nbterm(ks,kv)
         arg=phas(k,ks,kv)+freq(k,ks,kv)*T
         s=s+ampl(k,ks,kv)*cos(arg)
	end do
	elem(1)=s

        kv=2
         s=al(1,ks) +al(2,ks)*T
        do k=1,nbterm(ks,kv)
         arg=phas(k,ks,kv)+freq(k,ks,kv)*T
         s=s+ampl(k,ks,kv)*sin(arg)
	end do
	s=mod(s+val(1),pi2)
	if(s.lt.0) s=s+pi2
	elem(2)=s

        kv=3
         s1=0.d0
         s2=0.d0
        do k=1,nbterm(ks,kv)
         arg=phas(k,ks,kv)+freq(k,ks,kv)*T
         s1=s1+ampl(k,ks,kv)*cos(arg)
         s2=s2+ampl(k,ks,kv)*sin(arg)
	end do
	elem(3)=s1+val(2)
	elem(4)=s2+val(3)

        kv=4
         s1=0.d0
         s2=0.d0
        do k=1,nbterm(ks,kv)
         arg=phas(k,ks,kv)+freq(k,ks,kv)*T
         s1=s1+ampl(k,ks,kv)*cos(arg)
         s2=s2+ampl(k,ks,kv)*sin(arg)
	end do
	elem(5)=s1+val(4)
	elem(6)=s2+val(5)
*
*  computing cartesian coordinates from elements
*
        if(iv.EQ.1)then
          call ELEM2PV(MU(ks),elem,XV)
*
*  then transforming to cartesian coordinates in the
*      J2000 mean equator and mean equinox reference frame
*
          XVE(1)= XV(1)*DCOS(ome)-XV(2)*DSIN(ome)*DCOS(ainc)
     &           +XV(3)*DSIN(ainc)*DSIN(ome)
          XVE(2)= XV(1)*DSIN(ome)+XV(2)*DCOS(ome)*DCOS(ainc)
     &           -XV(3)*DSIN(ainc)*DCOS(ome)
          XVE(3)= XV(2)*DSIN(ainc)+XV(3)*DCOS(ainc)
          XVE(4)= XV(4)*DCOS(ome)-XV(5)*DSIN(ome)*DCOS(ainc)
     &           +XV(6)*DSIN(ainc)*DSIN(ome)
          XVE(5)= XV(4)*DSIN(ome)+XV(5)*DCOS(ome)*DCOS(ainc)
     &           -XV(6)*DSIN(ainc)*DCOS(ome)
          XVE(6)= XV(5)*DSIN(ainc)+XV(6)*DCOS(ainc)
        end if
        return
       end
***************************************************
       subroutine INIT
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
      !*   Read all constants and series from the file  GalileanL1.2.dat  *
      !*                                                                  *
      !*        Caution : this file contains comments which are skipped   *
      !*                  by appropriate read instructions                *
      !*        In case of problems in skipping the comments, use         *
      !*        the subroutine INITbis below  with the file  BisL1.2.dat  *
      !*                                                                  *
      !********************************************************************
      !*                          15/02/2006                              *
      !********************************************************************
       implicit real*8 (a-h,o-z)
       character*188 ligne
       real*8 ampl(100,4,4),phas(100,4,4),freq(100,4,4),al(2,4)
       real*8 c(9,20),c3(9,5,4)
       equivalence(c(1,1),c3(1,1,1))
       real*8 dfr(100,4,4),dph(100,4,4)
       integer nbterm(4,4),icomb(17,100,4,4)
       real*8 phasc(100,4,4),freqc(100,4,4)
       real*8 frq0(17),pha0(17)
       real*8 MU(4),A0(4),masses(4)    !,amu(4)
       common /serL1/nbterm,T0,T1,T2,al,ampl,phas,freq,c,MU,ome,ainc
       common /serL12/ frq0,pha0,phasc,freqc
       common /axes/ A0 ,masses    !,psi,aI

       pi=4.d0*datan(1.d0)
*
*  read the file GalileanL1.2.dat , skipping the comments
*
       open(15,file='GalileanL1.2.dat')
       read(15,'(a)') ligne ! 'Series for ephemerides ...'
       read(15,*)
       read(15,*) T0   ! 'T0 and Fundamental arguments ...'
       do i=1,17
         read(15,*) ii,pha0(i),frq0(i)
       end do
       read(15,*)
       read(15,'(a)') ligne ! 'masses [G x (Jupiter + satellite) ...'
       read(15,*) (MU(i),i=1,4)
       read(15,*)
       read(15,'(a)') ligne ! 'rotations ... '
       read(15,*) ome, ainc
       read(15,*)
       read(15,'(a)') ligne ! 'amplitude (AU) ... '
       read(15,*)
       do ks=1,4
        do kv=1,4
         read(15,'(a)') ligne  ! 'sat  ...'
         read(15,*) nbterm(ks,kv)
	 if (kv.EQ.2) then
         read(15,'(a)') ligne
	 read(ligne(4:54),*)  al(1,ks),al(2,ks)
	 end if
	 do nb=1,nbterm(ks,kv)
	 read(15,*)nnb,amp,phas(nb,ks,kv),freq(nb,ks,kv),ampkm,
     &     (icomb(k,nb,ks,kv),k=1,17),dfr(nb,ks,kv),dph(nb,ks,kv)
         ampl(nb,ks,kv)=amp
500    format(i3,f20.16,2d22.14,17i4,2d12.3)   !,f13.3,17i4,2d12.3)
	 end do
	 read(15,*)
	end do
       end do

       read(15,*) T1,T2   !  'T1, T2 (years) Tchebycheff polynomials...'
       DO ks=1,4
        read(15,'(a)') ligne  !  'satellite',ks
        read(15,'(a)') ligne  ! ' L','Re(z)','Im(z)', ...
	do j=1,9
	  read(15,'(i3,5d25.15)') jj,(c3(j,kv,ks),kv=1,5)
	end do
	read(15,*)
       ENDDO
       close(15)

       do i=1,4
         A0(i)=ampl(1,i,1)
         masses(i)=MU(i)
       end do
*
* computing phases and frequencies from the integer combinations
*     when such a combination is available
*
       do ks=1,4
        do kv=1,4
         do nb=1,nbterm(ks,kv)
          if(icomb(k,nb,ks,kv).ne. 99) then
           sp=0
           sf=0
           do k=1,17
            sp=sp+icomb(k,nb,ks,kv)*pha0(k)
            sf=sf+icomb(k,nb,ks,kv)*frq0(k)
           end do
           phasc(nb,ks,kv)=sp
           freqc(nb,ks,kv)=sf
          end if
          if(dph(nb,ks,kv).eq. 999.d0) then
           phasc(nb,ks,kv)=phas(nb,ks,kv)
           freqc(nb,ks,kv)=freq(nb,ks,kv)
          end if
          if(dfr(nb,ks,kv).eq. 999.d0) then
           freqc(nb,ks,kv)=freq(nb,ks,kv)
          end if
         end do
        end do
       end do
       return
       end
***************************************************
       subroutine INITbis
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
      !*   Read all constants and series from the file  BisL1.2.dat       *
      !*                                                                  *
      !*                                                                  *
      !********************************************************************
      !*                          15/02/2006                              *
      !********************************************************************
       implicit real*8 (a-h,o-z)
       character*188 ligne
       real*8 ampl(100,4,4),phas(100,4,4),freq(100,4,4),al(2,4)
       real*8 c(9,20),c3(9,5,4)
       equivalence(c(1,1),c3(1,1,1))
       real*8 dfr(100,4,4),dph(100,4,4)
       integer nbterm(4,4),icomb(17,100,4,4)
       real*8 phasc(100,4,4),freqc(100,4,4)
       real*8 frq0(17),pha0(17)
       real*8 MU(4),A0(4),masses(4)    !,amu(4)
       common /serL1/nbterm,T0,T1,T2,al,ampl,phas,freq,c,MU,ome,ainc
       common /serL12/ frq0,pha0,phasc,freqc
       common /axes/ A0 ,masses

       pi=4.d0*datan(1.d0)

*
*  read the file BisL1.2.dat , same as GalileanL1.2.dat without comments
*
       open(16,file='BisL1.2.dat')
       read(16,'(D18.10)') T0
       do i=1,17
         read(16,'(i3,2D23.15)') ii,pha0(i),frq0(i)
       end do
       read(16,'(4D25.16)') (MU(i),i=1,4)
       read(16,'(2D25.16)') ome, ainc
       do ks=1,4
        do kv=1,4
         read(16,'(i3)') nbterm(ks,kv)
	 if (kv.EQ.2) then
	 read(16,'(2D25.16)')  al(1,ks),al(2,ks)
	 end if
	 do nb=1,nbterm(ks,kv)
         read(16,500)nnb,amp,phas(nb,ks,kv),freq(nb,ks,kv),
     &     (icomb(k,nb,ks,kv),k=1,17),dfr(nb,ks,kv),dph(nb,ks,kv)
         ampl(nb,ks,kv)=amp
500    format(i3,f20.16,2d22.14,17i4,2d12.3)
	 end do
	end do
       end do
       read(16,'(2D25.15)') T1,T2
       DO ks=1,4
	do j=1,9
	  read(16,'(i3,5d25.15)') jj,(c3(j,kv,ks),kv=1,5)
	end do
       ENDDO
       close(16)

       do i=1,4
         A0(i)=ampl(1,i,1)
         masses(i)=MU(i)
       end do
*
* computing phases and frequencies from the integer combinations
*     when such a combination is available
*
       do ks=1,4
        do kv=1,4
         do nb=1,nbterm(ks,kv)
          if(icomb(k,nb,ks,kv).ne. 99) then
           sp=0
           sf=0
           do k=1,17
            sp=sp+icomb(k,nb,ks,kv)*pha0(k)
            sf=sf+icomb(k,nb,ks,kv)*frq0(k)
           end do
           phasc(nb,ks,kv)=sp
           freqc(nb,ks,kv)=sf
          end if
          if(dph(nb,ks,kv).eq. 999.d0) then
           phasc(nb,ks,kv)=phas(nb,ks,kv)
           freqc(nb,ks,kv)=freq(nb,ks,kv)
          end if
          if(dfr(nb,ks,kv).eq. 999.d0) then
           freqc(nb,ks,kv)=freq(nb,ks,kv)
          end if
         end do
        end do
       end do
       return
       end
***************************************************
       subroutine DL1_2FRC(ET,ks,iv,ELEM,XVE)

      !********************************************************************
      !*     Semi-analytic theory of motion of Galilean satellites        *
      !*          fitted to observations between 1891 and 2003            *
      !*                                                                  *
      !*     Computes ephemerides from series where arguments are         *
      !*     recomputed from the integer combination of fundamental       *
      !*     arguments (when these combinations are available)            *
      !*                                                                  *
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
      !*  ET : Julian Day  ie  JD in the Terrestrial Time scale (TT)      *
      !*       ET  must be limited to dates between about                 *
      !*             years 1140 and 2760                                  *
      !*  ks : satellite number, from 1 (Io) to  4 (Callisto)             *
      !*  iv : choice of output :                                         *
      !*       0 = orbital elements referred to the center of Jupiter in  *
      !*           a fixed reference frame close to the jovian equatorial *
      !*           frame for J2000.0  epoch,                              *
      !*       1 = cartesian jovicentric coordinates in the fixed         *
      !*           celestial reference frame : mean equator and mean      *
      !*           equinox for J2000.0 epoch.                             *
      !* ELEM(6) : orbital elements output, in order :                    *
      !*            a, L, Re(z), Im(z), Re(zeta), Im(zeta)                *
      !*            a in AU   while  other elements in rad                *
      !* XVE(6)  : cartesian coordinates output, in order :               *
      !*                    x, y, z,  vx, vy, vz                          *
      !*              x, y, z in AU   and   vx, vy, vz  in AU/day         *
      !********************************************************************
      !*                          15/02/2006                              *
      !********************************************************************

       implicit real*8 (a-h,o-z)
       real*8 ampl(100,4,4),phas(100,4,4),freq(100,4,4),al(2,4)
       real*8 phasc(100,4,4),freqc(100,4,4)
       real*8 elem(6),MU(4),XV(6),XVE(6),c(9,5,4),val(5),TN(0:8)
       integer nbterm(4,4)
       logical AccurateEphem
       real*8 frq0(17),pha0(17)
       integer icomb(17,100,4,4)
       common /serL1/nbterm,T0,T1,T2,al,ampl,phas,freq,c,MU,ome,ainc
       common /serL12/ frq0,pha0,phasc,freqc

       pi = 4.d0*datan(1.d0)
       pi2 = pi*2.d0
       AccurateEphem = .true.

       T=ET-T0  ! (T0=2433282.5d0 => jours juliens a partir du 01/01/1950 a 0H00 (TT)
       val=0.D0
*
*  computing corrections to elements from Tchebycheff polynomials :
*  These corrections represent very long period perturbations which
*  are not computed in form of trigonometric terms but only approximated
*  by Tchebycheff polynomials, hence limiting the interval of validity
*  of ephemerides to about [-810,810] years centered on J1950.
*
*  To compute approximate ephemerides at longer time, skip
*  these corrections by giving a .false. value to AccurateEphem
*
       if( AccurateEphem ) then
        a= T1    !  -819.727638594856D0
        b= T2    !   812.721806990360D0
        x=(T/365.25D0-0.5D0*(b+a))/(0.5D0*(b-a))
	if(abs(x).gt.1.d0) then
	  write(*,*)'date',T,' out of interval of validity'
	  return
	end if
        TN(0)=1.D0
        TN(1)=x
        do it=2,8
          TN(it)=2.D0*x*TN(it-1)-TN(it-2)
        end do
        do nv=1,5
          do it=0,8
          val(nv)=val(nv)+c(it+1,nv,ks)*TN(it)
          end do
          val(nv)=val(nv)-0.5D0*c(1,nv,ks)
        end do
       end if
*
*  computing elements from series
*
        kv=1
        s=0.d0
        do k=1,nbterm(ks,kv)
         arg=phasc(k,ks,kv)+freqc(k,ks,kv)*T
         s=s+ampl(k,ks,kv)*cos(arg)
	end do
	elem(1)=s

        kv=2
         s=al(1,ks) +al(2,ks)*T
        do k=1,nbterm(ks,kv)
         arg=phasc(k,ks,kv)+freqc(k,ks,kv)*T
         s=s+ampl(k,ks,kv)*sin(arg)
	end do
	s=mod(s+val(1),pi2)
	if(s.lt.0) s=s+pi2
	elem(2)=s

        kv=3
         s1=0.d0
         s2=0.d0
        do k=1,nbterm(ks,kv)
         arg=phasc(k,ks,kv)+freqc(k,ks,kv)*T
         s1=s1+ampl(k,ks,kv)*cos(arg)
         s2=s2+ampl(k,ks,kv)*sin(arg)
	end do
	elem(3)=s1+val(2)
	elem(4)=s2+val(3)

        kv=4
         s1=0.d0
         s2=0.d0
        do k=1,nbterm(ks,kv)
         arg=phasc(k,ks,kv)+freqc(k,ks,kv)*T
         s1=s1+ampl(k,ks,kv)*cos(arg)
         s2=s2+ampl(k,ks,kv)*sin(arg)
	end do
	elem(5)=s1+val(4)
	elem(6)=s2+val(5)
*
*  computing cartesian coordinates from elements
*
        if(iv.EQ.1)then
          call ELEM2PV(MU(ks),elem,XV)
*
*  then transforming to cartesian coordinates in the
*      J2000 mean equator and mean equinox reference frame
*
          XVE(1)= XV(1)*DCOS(ome)-XV(2)*DSIN(ome)*DCOS(ainc)
     &           +XV(3)*DSIN(ainc)*DSIN(ome)
          XVE(2)= XV(1)*DSIN(ome)+XV(2)*DCOS(ome)*DCOS(ainc)
     &           -XV(3)*DSIN(ainc)*DCOS(ome)
          XVE(3)= XV(2)*DSIN(ainc)+XV(3)*DCOS(ainc)
          XVE(4)= XV(4)*DCOS(ome)-XV(5)*DSIN(ome)*DCOS(ainc)
     &           +XV(6)*DSIN(ainc)*DSIN(ome)
          XVE(5)= XV(4)*DSIN(ome)+XV(5)*DCOS(ome)*DCOS(ainc)
     &           -XV(6)*DSIN(ainc)*DCOS(ome)
          XVE(6)= XV(5)*DSIN(ainc)+XV(6)*DCOS(ainc)
        end if
        return
       end
************************************************************
      SUBROUTINE ELEM2PV(MU,ELEM,XV)

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
      !*       Convert orbital elements to cartesian coordinates          *
      !*       in the same frame                                          *
      !*                                                                  *
      !*  MU     : gravitationnal constant of the keplerian motion        *
      !*           Here, values are in UA^3 / day^2                       *
      !* ELEM(6) : orbital elements input, in order :                     *
      !*            a, L, Re(z), Im(z), Re(zeta), Im(zeta)                *
      !*            a in AU   while  other elements in rad                *
      !* XVE(6)  : cartesian coordinates output, in order :               *
      !*                    x, y, z,  vx, vy, vz                          *
      !*              x, y, z in AU   and   vx, vy, vz  in AU/day         *
      !********************************************************************
      !*                          15/02/2006                              *
      !********************************************************************
*
*  given at time  t  :
*   x   \           /   a    semi-major axe
*   y   |    mu     |   L    mean longitude (= OMEGA+omega+M)
*   z   |  /-----\  |   k    e cos(OMEGA+omega)
*  vx   |  \-----/  |   h    e sin(OMEGA+omega)
*  vy   |           |   q    sin(i/2) cos(OMEGA)
*  vz   /           \   p    sin(i/2) sin(OMEGA)
*
*  mu is the gravitational constant of keplerian elliptic motion (focus F)
*  a  is the semi-major axe
*  M  is the mean anomaly at time  t
*     (M=n(t-to) with  n  mean osculating motion such as  n^2a^3=mu
*                and  to  time of passage at pericentre)
*  e  is the eccentricity
*  i  is the inclination of orbital plane on Fxy plane
*  OMEGA  is the longitude of ascending node on Fxy plane
*  omega  is the argument of pericentre such as
*         OMEGA + omega represents the longitude of pericentre
*     and  L (= OMEGA+omega+M) représents the mean longitude
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(PI2=3.141592653589793D0*2.d0)
      DOUBLE PRECISION MU,K
      DIMENSION XV(6),ELEM(6),G(3),E(3)
      K=ELEM(3)
      H=ELEM(4)
      Q=ELEM(5)
      P=ELEM(6)
      A=ELEM(1)
      AL=ELEM(2)
      AN=DSQRT(MU/A**3.D0)
      EE=AL+K*DSIN(AL)-H*DCOS(AL)
20    CONTINUE
      CE=DCOS(EE)
      SE=DSIN(EE)
      DE=(AL-EE+K*SE-H*CE)/(1.D0-K*CE-H*SE)
      EE=EE+DE
      IF (ABS(DE).GE.1.D-12) GOTO 20
      CE=DCOS(EE)
      SE=DSIN(EE)
      DLE=H*CE-K*SE
      RSAM1=-K*CE-H*SE
      ASR=1.D0/(1.D0+RSAM1)
      PHI=DSQRT(1.D0-K*K-H*H)
      PSI=1.D0/(1.D0+PHI)
      X1=A*(CE-K-PSI*H*DLE)
      Y1=A*(SE-H+PSI*K*DLE)
      VX1=AN*ASR*A*(-SE-PSI*H*RSAM1)
      VY1=AN*ASR*A*( CE+PSI*K*RSAM1)
      F2=2.D0*DSQRT(1.d0-Q*Q-P*P)
      P2=1.D0-2.D0*P*P
      Q2=1.D0-2.D0*Q*Q
      PQ=2.D0*P*Q
      XV(1)=X1*P2+Y1*PQ
      XV(2)=X1*PQ+Y1*Q2
      XV(3)=(Q*Y1-X1*P)*F2
      XV(4)=VX1*P2+VY1*PQ
      XV(5)=VX1*PQ+VY1*Q2
      XV(6)=(Q*VY1-VX1*P)*F2
      RETURN
      END
