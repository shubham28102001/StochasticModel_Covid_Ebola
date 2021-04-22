      program EVD
      implicit none

      integer*8 l,j,m,nodes,nsteps,xi,np,ncounts,n
      parameter (nodes = 500000, nsteps = 2000000000, np = 16, !Define Population Size (nodes) and End Conditions
     $     ncounts = 1)
      real*8 mu,bi,bd,bh,k,sigma,gir,mue,tau,delta,ghr,muh,alpha
      real*8 Ot,Of,Os,Rz,Eu,Er
      real*8 dt,t,dk
      integer*8 S,E,I,R,D,H
      real*8 ran2,rand1,rand2,ran3,dummy(3)
      real*8  p(np),a0,atemp,tol
      integer seeds(ncounts),seed,ext_count
      character(len=20) :: filename

            !  parameters
               mu=0.00005d0
               Bi=0.5
               Bd=0.6d0
               Bh=0.00016d0
               k=0.000000005d0
               sigma=0.10d0
               gir=0.07d0
               mue=0.12d0
               tau=0.20d0
               delta=0.33d0
               ghr=0.10d0
               muh=0.12d0


               write (filename,"('Simulation_data')")
               open (111,file=filename,status = 'new')




c              seed generation
               print *,'Enter the seed'
               read(*,*) seed
               seed = -seed


                  ! This is the do that determines how many simulations are
                  ! done with a single parameter set.

                  do n = 1,ncounts  !A
                     print *,n
                     ext_count = 0


      ! Build Initial Conditions
      Ot=sigma/(gir+tau+mue+mu)
      Of=tau*sigma/((gir+tau+mue+mu)*(ghr+mue+mu))
      Os=(gir*Ot+ghr*Of)/mu
      Rz=Ot*(Bi+Bd*mue/(delta+mu)+Bh*tau/(ghr+mue+mu))/(mu+sigma)
      Eu=(Rz-1-k*alpha/mu)
      Er=(Eu+sqrt(Eu**2+4*Rz*k*alpha/mu))*mu*nodes/(2*Rz*(mu+sigma))
      E=int((Eu+sqrt(Eu**2+4*Rz*k*alpha/mu))*mu*nodes/(2*Rz*(mu+sigma)))
      S=int(nodes-(1+sigma/mu)*dble(Er))
      I=int(Ot*dble(Er))
      D=int(mue*Ot*dble(Er)/(delta+mu))
      H=int(tau*Ot*dble(Er)/(ghr+mue+mu))
      R=int(Os*dble(Er))




                  t = 0.d0
                  dk = 0.2d0    !This determines how often the disease state will be recorded (in days)




 100     format(f15.3,2x,7(I13,2x))


            ! This do loop is what actually runs the simulation
                  do l = 1,nsteps   !B



!                       This is extinction event end condition
                        if(ext_count.gt.300) then
                              goto 300
                        endif

                        !Generate Random numbers so that an event and a time step size can be chosen randomly
                        rand1 = ran3(seed)
                        rand2 = ran3(seed)

                        do while (rand1.eq.0.d0)
                              rand1 = ran3(seed)
                        enddo
                        do while (rand2.eq.0.d0)
                              rand2 = ran3(seed)
                        enddo


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Define Transition Rates (These will effectively act as relative probabilities)
c Each p() is a possible event.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                  p(1) = mu*dble(nodes)  ! Birth N -> S
                  p(2) = mu*dble(S)   ! death while susceptible
                  p(3) = (Bi*dble(I)+Bd*dble(D) + Bh*dble(h))*dble(s)
     $                     /dble(nodes)        ! S -> E
                  p(4) = k*dble(S) ! exposed from reservoir
                  p(5) = mu*dble(E) ! natural death while exposed
                  p(6) = sigma*dble(E) ! exposed -> infected
                  p(7) = mue*dble(I) ! I -> D
                  p(8) = tau*dble(I) ! I -> H
                  p(9) = gir*dble(I) ! I -> R
                  p(10) = delta*dble(D) ! D - >
                  p(11) = muh*dble(H) ! H ->
                  p(12) = mu*dble(R) ! R ->
                  p(13) = ghr*dble(H) ! H -> R
                  p(14) = mu*dble(H) ! H ->
                  p(15) = mu*dble(I) ! I ->
                  p(16) = mu*dble(D) ! D ->


                        ! sum up all of the different transition rates.
                        a0 = 0.d0
                        do j = 1,np
                              a0 = a0 + p(j)
                        enddo


!                       time step
                        dt = 1.d0/a0*log(1.d0/rand1) ! exponential time step
                        t = t + dt

!                       attuating the even
                        atemp = 0.d0
                        do j = 1,np

!                             conditions on the probabilities
                              if (p(j).eq.0.d0) then
                                    goto 50
                              endif

                              ! Determine which event actually happens
                              atemp = atemp + p(j)
                              if(rand2*a0.le.atemp) then
                                    xi = j
                                    goto 200
                              endif

 50                           continue
                        enddo

                        !Change the size of different state variables based on which event happened.
 200                    continue
                              if(xi.eq.1) then    ! 0 -> S
                                    S = S + 1
                              else if(xi.eq.2) then ! S ->
                                    s = s - 1
c                       infection
                              else if(xi.eq.3) then ! S -> E
                                    S = S - 1
                                    E = E + 1
c                       Exposed from reservoir
                              else if(xi.eq.4) then ! S -> E
                                    S = S - 1
                                    E = E + 1
                              else if(xi.eq.5) then ! E ->
                                    E = E - 1
                              else if(xi.eq.6) then ! E -> I
                                    E = E - 1
                                    I = I + 1
                              else if(xi.eq.7) then ! I -> D
                                    I = I - 1
                                    D = D + 1
                              else if(xi.eq.8) then ! I -> H
                                    I = I - 1
                                    H = H + 1
                              else if(xi.eq.9) then ! I -> R
                                    I = I - 1
                                    R = R + 1
                              else if(xi.eq.10) then ! D ->
                                    D = D - 1

                              else if(xi.eq.11) then ! H ->
                                    H = h -1

                              else if(xi.eq.12) then ! R ->
                                    R = R - 1
                              else if(xi.eq.13) then ! H -> R
                                    H = H  - 1
                                    R = R  + 1
                              else if(xi.eq.14) then !  H ->
                                    h = h - 1
                              else if(xi.eq.15) then ! I ->
                                    i = i - 1
                              else if(xi.eq.16) then !  D ->
                                    d = d - 1

                              else
                                    print *,'Error: xi =  ', xi
                                    stop
                              endif


c There are many different things that you might track during simulations. This is
c an example of a script that keeps track of the number of extinction events. The
c code checks every .2 days to find out what the size of the infected compartment is, and if
c the size is zero then it increases the extinction count.
c It will also track the system progress as well.


                        if (t.gt.dk) then
                              if (I.eq.0) then
                                    ext_count = ext_count + 1
                              endif
                              write(111,*) t, S, E, I, D, H, R
                              dk = dk + 0.2
!                             This is 'max time' end condition
!                              print *,t
                              if(t.gt.1000) then
                                  goto 300
                              endif
                        endif


                  enddo
 300              continue
                  open (unit=110,file='Extinction_Times',status = 'new')
                  write(110,*) ext_count
                  close(110)
                  enddo   !A
                  close(111)
      END

      FUNCTION ran3(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
C     REAL*8 MBIG,MSEED,MZ
      REAL*8 ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.d0/MBIG)
C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
C     REAL*8 mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
         iff=1
         mj=MSEED-iabs(idum)
         mj=mod(mj,MBIG)
         ma(55)=mj
         mk=1
         do 11 i=1,54
            ii=mod(21*i,55)
            ma(ii)=mk
            mk=mj-mk
            if(mk.lt.MZ)mk=mk+MBIG
            mj=ma(ii)
 11      continue
         do 13 k=1,4
            do 12 i=1,55
               ma(i)=ma(i)-ma(1+mod(i+30,55))
               if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
 12         continue
 13      continue
         inext=0
         inextp=31
         idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      return
      END



