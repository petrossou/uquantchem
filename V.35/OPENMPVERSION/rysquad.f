      SUBROUTINE RysQuad(NROOTS,XX,rts,wts)
      implicit real*8(a-h,o-z)
c
      parameter (mxquad=25)
      dimension rts(mxquad),wts(mxquad),wrk(mxquad),
     *          alpha(0:mxquad-1),beta(0:mxquad-1),
     *          xasymp(mxquad)
c
      parameter (mxaux=100)
      dimension rtsaux(mxaux),wtsaux(mxaux),
     *          alpaux(0:mxaux),betaux(0:mxaux),
     *          p0(mxaux),p1(mxaux),p2(mxaux),
     *          rgrid(mxaux),wgrid(mxaux),nauxs(mxaux)
c
      character*10 type,code,auxqd
      logical RPprt
c
c         Rys' RMS agreement between asymptotic code and Gautschi code
c         goes below 1d-13 for both roots and weights above these, as
c         it seems Ishida-san was a bit too conservative about this.
c
      data (xasymp(I),I=1,13)/
     *    29.0d+00,  37.0d+00,  43.0d+00,  49.0d+00,  55.0d+00,
     *    60.0d+00,  65.0d+00,  71.0d+00,  76.0d+00,  81.0d+00,
     *    86.0d+00,  91.0d+00,  96.0d+00/
c
c         Rys auxiliary grid size determined at the xasymp values,
c         eyeballing against using a naux=80 "converged" grid size.
c         It is likely that at lower X values, a smaller auxiliary
c         quadrature would be adequate.
c
      data (nauxs(I),I=1,13)/
     *    20,25,30,30,35,  40,40,40,45,50,  50,55,55/
c
c     Program to develop Gaussian quadrature roots and weights
c       assembled/written by Mike Schmidt in November 2005.
c
c     User can select the type of quadrature here, the choices are
c       Hermite, Laguerre, Legendre, or Rys
c
      type = 'Rys'
c
c     ================================================================
c     for the theory of Orthogonal Polynomials, see the following:
c
c             "Calculation of Gauss Quadrature Rules"
c                 Gene H. Golub, John H. Welsch
c           Mathematics of Computation, 23, 221-230(1969)
c
c     "Algorithm 726: ORTHPOL - A Package of Routines for Generating
c         Orthogonal Polynomials and Gauss-Type Quadrature Rules"
c                         Walter Gautschi
c       ACM Transactions on Mathematical Software  20, 21-62(1994)
c
c         Orthogonal Polynomials: Computation and Approximation
c            W.Gautschi, Oxford University Press, Oxford 2004
c              see Theorem 1.31 and Table 1.1 in this book.
c     
c     The orthogonal polynomial family is defined on an interval [a,b]
c     with respect to a weight function, w(t), 
c         integral p (t) * p (t) * w(t) = delta(i,j) == <p |p >
c           a->b    i       j                             i  j
c     The inner product notation <|> implies the use of the weight
c     function (metric) in the integration.  Any such family will
c     obey the same recursion formula, namely,
c        p   (t) = (t-a )p (t) - b p   (t)      k=0,1,2,...
c         k+1          k  k       k k-1
c     where a and b are called alpha,beta in this program, and are
c                <t*p |p >               <p |p >
c                    k  k                  k  k
c           a  = -----------      b  = -----------
c            k    <p |p >          k   <p   |p   >
c                   k  k                 k-1  k-1
c     We must define the first two members in each family as
c        p  (t) = 0    p (t) = 1
c         -1            0
c     and state that k >= 0 in the a-k formula, while k >=1 in the b-k
c     formula.  Then, stating that b-0 = <p-0|p-0>, e.g. the integral
c     of the weight function only, completes the definition of the
c     orthogonal polynomial.
c
c     Thus, forming alpha,beta recurrence coefficients is tantamount
c     to forming the orthogonal polynomial.  Furthermore, use of the
c     Golub-Walsh library routine 'gauss' will give the roots as the
c     eigenvalues of a tridiagonal matrix (alpha=diagonal, beta=next
c     to diag), and the weights are the first element of each vector.
c     ================================================================
c
c      if(type.eq.'Rys') then
c         write(6,900)
c         read(5,*) NROOTS,XX
c      else
c         write(6,910) type
c         read(5,*) NROOTS
c      end if
      if(NROOTS.gt.mxquad) stop 'fix maximum quadrature dimension'
c
      PI=3.14159265358979323846264338327950288419716939937510D+00
      zero= 0.0D+00
      one = 1.0D+00
      two = 2.0D+00
      four= 4.0D+00
c          epsilon is a relative error tolerance in the QL routine
      eps = 1.0D-14
c
      if(type.eq.'Rys') then
         do i=14,mxquad
            xasymp(i) = 1.0d+05
            nauxs(i) = mxaux
         enddo
      end if
c
c        Gauss-Hermite quadrature
c        interval (-inf,+inf)    weight exp(-t*t)
c
      if(type.eq.'Hermite') then
         do i=0,NROOTS-1
            alpha(i) = zero
         enddo
         beta(0)=sqrt(pi)
         do i=1,NROOTS-1
            beta(i) = i/two
         end do
         call gauss(NROOTS,alpha,beta,eps,rts,wts,ierr,wrk)
      end if
c
c        Gauss-Laguerre quadrature
c        interval [0,+inf)       weight exp(-t)
c
      if(type.eq.'Laguerre') then
         do i=0,NROOTS-1
            alpha(i) = two*i + one
         enddo
         beta(0)= one
         do i=1,NROOTS-1
            beta(i) = i*i
         end do
         call gauss(NROOTS,alpha,beta,eps,rts,wts,ierr,wrk)
      end if
c
c        Gauss-Legendre quadrature
c        interval [-1,+1]        weight unity
c
      if(type.eq.'Legendre') then
         do i=0,NROOTS-1
            alpha(i) = zero
         enddo
         beta(0)= two
         do i=1,NROOTS-1
            beta(i) = one/(four-(one/(i*i)))
         end do
         call gauss(NROOTS,alpha,beta,eps,rts,wts,ierr,wrk)
      end if
c
c        Gauss-Rys quadrature
c        interval [0,+1]         weight exp(-XX*t*t)
c        Introduced by John Rys, Michel Dupuis, Harry King in 1976,
c        to evaluate Cartesian Gaussian integrals in quantum chemistry:
c           M.Dupuis, J.Rys, H.F.King  J.Chem.Phys. 65, 111-116(1976)
c           J.Rys, M.Dupuis, H.F.King  J.Comput.Chem. 4, 154-157(1983)
c        The numerical parameter XX appearing in the weight function
c        depends on the distance between basis function centers, and
c        the exponents in the four Gaussians, 0 <= XX < infinite.
c
c        Rys polynomial computation was first discussed by
c          1. H.F.King, M.Dupuis  J.Comput.Phys. 21, 144-165(1976)
c        Other papers about the computation of Rys polynomials:
c          2. J.D.Augspurger, D.E.Bernholdt, C.E.Dykstra
c                J.Comput.Chem. 11, 972-977(1990)
c          3. K.Ishida  J.Chem.Phys. 95, 5198-5205(1991)
c                       J.Chem.Phys. 98, 2176-2181(1993)
c          4. R.B.Sagar, V.H.Smith Int.J.Quantum Chem. 42, 827-836(1992)
c          5. R.C.Y.Chin  J.Comput.Phys. 99, 321-336(1992)
c          6. P.Carsky, M.Polasek  J.Comput.Phys. 143, 266-277(1998)
c          7. B.I.Schneider, N.Nygaard J.Phys.Chem.A 106, 10773-6(2002)
c        of which, Sagar and Smith's is particularly illuminating.
c
c        It is known from paper 1 that in the large XX limit, the
c        Rys polynomial approaches a scaled Hermite polynomial.
c        See the 2nd paper by Ishida for more information on the
c        large XX limit.
c        The limit XX=0 is manifestly the shifted Legendre case.
c
c        The approach here is to provide a number of different
c        implementations, allowing convenient testing of their
c        speed and numerical accuracy.  
c
c        Sample numerical results can be found in references
c        1, 4, and 7, for NROOTS=5, 13, and 20.  Pictures in 1.
c
c        Using the Golub-Walsh procedure to obtain the roots and
c        weights, from recurrence coefficients -alpha- and -beta-
c        that are generated by the so-called "bootstrap method",
c        using a Stieltjes discretization of the Rys polynomial on
c        an auxiliary grid, is considered the "correct answer":
c                     code=Gautschi.
c
      if(type.eq.'Rys') then
c
c        Choose code from: Dupuis, Sierra, Schmidt, or Gautschi.
c            The last two can print the Rys polynomial, see RPprt.
c            The last two codes need an auxiliary quadrature on
c               our interval [0,1], select from: Fejer or ShiftLeg.
c            The first two codes are dimension limited to NROOTS=13.
c
         code='Gautschi'
         RPprt=.false.
c
         call rysset
         call setRaiz128
         call setPeso128
c
         if(XX.ge.xasymp(NROOTS)) then
            code='Asymptotic'
            call rysasy(NROOTS,XX,rts,wts)
         end if
c          
c            For the standard test case, namely NROOTS=13 and XX=10.0,
c            ShiftLeg is much more accurate than Fejer around naux=25,
c            and converges to 14 significant digits about naux=30.
c            Fejer converges to 14 figures only near naux=60 or so.
c            The auxiliary grid must be denser for larger XX, and the
c            values stored in -nauxs- were determined at the XX for
c            which we are able to shift to the asymptotic code.
c            These are therefore overly accurate for smaller XX!
c 
         if(code.eq.'Schmidt'  .or.  code.eq.'Gautschi') then
            auxqd = 'ShiftLeg'
            naux = nauxs(NROOTS)
         else
            auxqd = 'none'
            naux = 0
            if(NROOTS.gt.13) stop 'chosen Rys code stops at NROOTS=13'
         end if
         if(naux.gt.mxaux) stop 'fix dimension for auxiliary quad.'
c
         if(code.eq.'Schmidt'  .or.  code.eq.'Gautschi') then
c
c              Set up auxiliary "shifted Legendre quadrature",
c              which is weight function 1, on interval [0,1].
c              Note that choosing naux=20 exactly reproduces
c              the slightly wrong values in the first column
c              of Tables 3 and 4 given by Schneider/Nygaard,
c              and using naux=40 will give their Tables 1 and 2.
c
            if(auxqd.eq.'ShiftLeg') then
               do i=0,naux-1
                  alpaux(i) = one/two
               enddo
               betaux(0)= one
               quart = one/four
               do i=1,naux-1
                  betaux(i) = quart/(four-(one/(i*i)))
               end do
               call gauss(naux,alpaux,betaux,eps,rtsaux,wtsaux,ierr,p0)
            end if
c
c              Set up auxiliary "Fejer quadrature", which is a
c              generic choice for roots and weights.
c              naux=40 gives about 9 significant figure agreement to
c              Sagar/Smith's table IV, and naux=80 is almost perfect.
c
            if(auxqd.eq.'Fejer') then
               call fejer(naux,rtsaux,wtsaux)
            end if
         end if
c
c
c            This is the original HONDO routine, which loses
c            accuracy rapidly with NROOTS, even if run in the
c            quadruple precision.  In double precision, the
c            RMS error compared to the Gautschi code is largest
c            at small XX, and gradually decreases with XX:
c                         XX=0    XX=xasymp
c                    6    1d-11    1d-14
c                    7    1d-09    1d-12
c                    8    5d-08    1d-12
c                    9    2d-06    1d-10
c                   10    5d-05    1d-10
c                   11    3d-03    1d-7
c
         if(code.eq.'Dupuis') then
            call dupuis(NROOTS,XX,rts,wts)
         end if
c
c            Jose Sierra's arbitrary precision fitting method.
c            Note that m2dR=32 and m2dW=128 handle NROOTS=13, XX=10.0
c            but that serious numerical problems turn up at larger XX.
c
         if(code.eq.'Sierra') then
            m2dR= 128
            m2dW= 128
            do i=1,NROOTS
               if(m2dR.eq.128) call Raiz128(rts(i),i,NROOTS,XX,m2dR-2)
               if(m2dW.eq.128) call Peso128(wts(i),i,NROOTS,XX,m2dW-2)
            enddo
         end if
c
c           Generate alpha,beta recurrence coefficients, as well
c           as the Rys Polynomials, by Stieltjes "bootstrapping".
c           This code choice will also print the Rys polynomial.
c           The only thing "Schmidt" about this is the programming,
c           and it is accurate if executed in quadruple precision.
c
         if(code.eq.'Schmidt') then
            call RPboot(rtsaux,wtsaux,naux,NROOTS,XX,alpha,beta,RPprt)
            call gauss(NROOTS,alpha,beta,eps,rts,wts,ierr,wrk)
         end if
c
c           Gautschi's discretized Stieltjes "bootstrapping", which
c           bootstraps the numerical values of the Rys polynomials on
c           the grid, yielding numerically stable alpha,beta values.
c           Of course, the Rys polynomial can be generated using the
c           final alpha,beta values in the recurrence.  If you
c           want to be sure about having the polynomial exactly,
c           run in quadruple precision, at least in RysGen.
c
         if(code.eq.'Gautschi') then
            do i=1,naux
               t2 = rtsaux(i)*rtsaux(i)
               rgrid(i) = t2
               wgrid(i) = wtsaux(i)*exp(-XX*t2)
            enddo
            call sti(NROOTS,naux,rgrid,wgrid,alpha,beta,ierr,p0,p1,p2)
            if(ierr.ne.0) stop 'discretized Stieltjes failed'
            if(RPprt) call RysGen(alpha,beta,NROOTS)
            call gauss(NROOTS,alpha,beta,eps,rts,wts,ierr,wrk)
         end if
      end if
c
c     the following converts to everyone's convention for Tables.
c
      if(type.eq.'Rys') then
         do i=1,NROOTS
            rts(i) = sqrt(rts(i))
         enddo
      end if
c
c      if(type.eq.'Rys') then
c         write(6,920) XX
c      else
c         write(6,930) type
c      end if
c      do i=1,NROOTS
c         write(6,940) i,rts(i),wts(i)
c      enddo
c      if(type.eq.'Rys') then
c         if(code.eq.'Asymptotic') write(6,950) code
c         if(code.eq.'Dupuis')     write(6,950) code
c         if(code.eq.'Sierra')     write(6,960) code,m2dR,m2dW
c         if(code.eq.'Schmidt')    write(6,970) code,auxqd,naux
c         if(code.eq.'Gautschi')   write(6,970) code,auxqd,naux
c      end if
      return
c
c  900 format('Enter desired Gauss-Rys NROOTS and X: ',$)
c  910 format('Enter desired order of Gauss-',a,' quadrature: ',$)
c  920 format('Roots and Weights for Gauss-Rys quadrature, at X=',f8.3)
c  930 format('Roots and Weights for Gauss-',a,' quadrature.')
c  940 format(1x,'i=',i4,1P,' root=',e21.14,' weight=',e21.14)
c  950 format('using Code=',a)
c  960 format('using Code=',a,' m2dR=',i5,'   m2dW=',i5)
c  970 format('using Code=',a,' AuxQd=',a,'   Naux=',i4)
      end
c
      subroutine gauss(n,alpha,beta,eps,zero,weight,ierr,wrk)
      implicit real*8(a-h,o-z)
      dimension alpha(n),beta(n),zero(n),weight(n),wrk(n)
c
c         this routine is a supplement to the journal article,
c      Algorithm 726: ORTHPOL - A Package of Routines for Generating
c        Orthogonal Polynomials and Gauss-Type Quadrature Rules.
c                         Walter Gautschi, 
c      ACM Transactions on Mathematical Software, 20, 21-62(1994)
c
c            A footnote in Gautschi's article reveals that
c            this particular routine came from Gene Golub.
c
c Given  n  and a measure  dlambda, this routine generates the n-point
c Gaussian quadrature formula
c 
c     integral over supp(dlambda) of f(x)dlambda(x)
c
c        = sum from k=1 to k=n of w(k)f(x(k)) + R(n;f).
c 
c The nodes are returned as  zero(k)=x(k) and the weights as 
c weight(k)=w(k), k=1,2,...,n. The user has to supply the recursion 
c coefficients  alpha(k), beta(k), k=0,1,2,...,n-1, for the measure 
c dlambda. The routine computes the nodes as eigenvalues, and the 
c weights in term of the first component of the respective normalized 
c eigenvectors of the n-th order Jacobi matrix associated with  dlambda.
c It uses a translation and adaptation of the algol procedure  imtql2,
c Numer. Math. 12, 1968, 377-383, by Martin and Wilkinson, as modified 
c by Dubrulle, Numer. Math. 15, 1970, 450. See also Handbook for 
c Autom. Comput., vol. 2 - Linear Algebra, pp.241-248, and the eispack
c routine  imtql2.
c
c        Input:  n - - the number of points in the Gaussian quadrature	
c                      formula; type integer
c                alpha,beta - - arrays of dimension  n  to be filled 
c                      with the values of  alpha(k-1), beta(k-1), k=1,2,
c                      ...,n
c                eps - the relative accuracy desired in the nodes
c                      and weights
c
c        Output: zero- array of dimension  n  containing the Gaussian 
c                      nodes (in increasing order)  zero(k)=x(k), k=1,2,
c                      ...,n
c                weight - array of dimension  n  containing the 
c                      Gaussian weights  weight(k)=w(k), k=1,2,...,n
c                ierr- an error flag equal to  0  on normal return,
c                      equal to  i  if the QR algorithm does not
c                      converge within 30 iterations on evaluating the 
c                      i-th eigenvalue, equal to  -1  if  n  is not in
c                      range, and equal to  -2  if one of the beta's is 
c                      negative.
c
c The array  wrk  is needed for working space.
c
      if(n.lt.1) then
        ierr=-1
        return
      end if
c
      ierr=0
      zero(1)=alpha(1)
      if(beta(1).lt.0.0D+00) then
        ierr=-2
        return
      end if
      weight(1)=beta(1)
      if (n.eq.1) return
c
      weight(1)=1.0D+00
      wrk(n)=0.0D+00
      do 100 k=2,n
        zero(k)=alpha(k)
        if(beta(k).lt.0.0D+00) then
          ierr=-2
          return
        end if
        wrk(k-1)=sqrt(beta(k))
        weight(k)=0.0D+00
  100 continue
c
      do 240 l=1,n
        j=0
c
c Look for a small subdiagonal element.
c
  105   do 110 m=l,n
          if(m.eq.n) go to 120
          if(abs(wrk(m)).le.eps*(abs(zero(m))+abs(zero(m+1)))) 
     *      go to 120
  110   continue
  120   dp=zero(l)
        if(m.eq.l) go to 240
        if(j.eq.30) go to 400
        j=j+1
c
c Form shift.
c
        dg=(zero(l+1)-dp)/(2.0D+00*wrk(l))
        dr=sqrt(dg*dg+1.0D+00)
        dg=zero(m)-dp+wrk(l)/(dg+sign(dr,dg))
        ds=1.0D+00
        dc=1.0D+00
        dp=0.0D+00
        mml=m-l
c
c For i=m-1 step -1 until l do ...
c
        do 200 ii=1,mml
          i=m-ii
          df=ds*wrk(i)
          db=dc*wrk(i)
          if(abs(df).lt.abs(dg)) go to 150
          dc=dg/df
          dr=sqrt(dc*dc+1.0D+00)
          wrk(i+1)=df*dr
          ds=1.0D+00/dr
          dc=dc*ds
          go to 160
  150     ds=df/dg
          dr=sqrt(ds*ds+1.0D+00)
          wrk(i+1)=dg*dr
          dc=1.0D+00/dr
          ds=ds*dc
  160     dg=zero(i+1)-dp
          dr=(zero(i)-dg)*ds+2.0D+00*dc*db
          dp=ds*dr
          zero(i+1)=dg+dp
          dg=dc*dr-db
c
c Form first component of vector.
c
          df=weight(i+1)
          weight(i+1)=ds*weight(i)+dc*df
          weight(i)=dc*weight(i)-ds*df
  200   continue
        zero(l)=zero(l)-dp
        wrk(l)=dg
        wrk(m)=0.0D+00
        go to 105
  240 continue
c         
c Order eigenvalues and eigenvectors.
c
      do 300 ii=2,n
        i=ii-1
        k=i
        dp=zero(i)
        do 260 j=ii,n
          if(zero(j).ge.dp) go to 260
          k=j
          dp=zero(j)
  260   continue
        if(k.eq.i) go to 300
        zero(k)=zero(i)
        zero(i)=dp
        dp=weight(i)
        weight(i)=weight(k)
        weight(k)=dp
  300 continue
      do 310 k=1,n
        weight(k)=beta(1)*weight(k)*weight(k)
  310 continue
      return
c
c Set error - no convergence to an eigenvalue after 30 iterations.
c
  400 ierr=l
      return
      end
c
      subroutine sti(n,ncap,x,w,alpha,beta,ierr,p0,p1,p2)
      implicit real*8(a-h,o-z)
      dimension x(ncap),w(ncap),alpha(n),beta(n),
     *          p0(ncap),p1(ncap),p2(ncap)
c
c         this routine is a supplement to the journal article,
c      Algorithm 726: ORTHPOL - A Package of Routines for Generating
c        Orthogonal Polynomials and Gauss-Type Quadrature Rules.
c                         Walter Gautschi, 
c      ACM Transactions on Mathematical Software, 20, 21-62(1994)
c
c This routine applies Stieltjes's procedure (cf. Section 2.1 of
c W. Gautschi,On generating orthogonal polynomials'', SIAM J. Sci.
c Statist. Comput. 3, 1982, 289-317) to generate the recursion
c coefficients  alpha(k), beta(k) , k=0,1,...,n-1, for the discrete
c (monic) orthogonal polynomials associated with the inner product
c
c     (f,g)=sum over k from 1 to ncap of w(k)*f(x(k))*g(x(k)).
c
c The integer  n  must be between  1  and  ncap, inclusive; otherwise,
c there is an error exit with  ierr=1. The results are stored in the
c arrays  alpha, beta; the arrays  p0, p1, p2  are working arrays.
c
c If there is a threat of underflow or overflow in the calculation
c of the coefficients  alpha(k)  and  beta(k), the routine exits with
c the error flag  ierr  set equal to  -k  (in the case of underflow) 
c or  +k  (in the case of overflow), where  k  is the recursion index 
c for which the problem occurs. The former [latter] can often be avoided
c by multiplying all weights  w(k)  by a sufficiently large [small]
c scaling factor prior to entering the routine, and, upon exit, divide
c the coefficient  beta(0)  by the same factor.
c
c This routine should be used with caution if  n  is relatively close 
c to  ncap, since there is a distinct possibility of numerical 
c instability developing. (See W. Gautschi,Is the recurrence relation 
c for orthogonal polynomials always stable?'', BIT, 1993, to appear.) 
c In that case, the routine  lancz  should be used.
c
      tiny = 1.0d-40
      huge = 1.0d+40
c
      ierr=0
      if(n.le.0 .or. n.gt.ncap) then
        ierr=1
        return
      end if
      nm1=n-1
c
c Compute the first alpha- and beta-coefficient.
c
      sum0=0.0d+00
      sum1=0.0d+00
      do 10 m=1,ncap
        sum0=sum0+w(m)
        sum1=sum1+w(m)*x(m)
   10 continue
      alpha(1)=sum1/sum0
      beta(1)=sum0
      if(n.eq.1) return
c
c Compute the remaining alpha- and beta-coefficients.
c
      do 20 m=1,ncap
        p1(m)=0.0d+00
        p2(m)=1.0d+00
   20 continue
      do 40 k=1,nm1
        sum1=0.0d+00
        sum2=0.0d+00
        do 30 m=1,ncap
c
c The following statement is designed to avoid an overflow condition
c in the computation of  p2(m)  when the weights  w(m)  go to zero 
c faster (and underflow) than the  p2(m)  grow. 
c
          if(w(m).eq.0.0d+00) goto 30
          p0(m)=p1(m)
          p1(m)=p2(m)
          p2(m)=(x(m)-alpha(k))*p1(m)-beta(k)*p0(m)
c
c Check for impending overflow.
c
          if(abs(p2(m)).gt.huge .or. abs(sum2).gt.huge) then
            ierr=k
            return
          end if
          t=w(m)*p2(m)*p2(m)
          sum1=sum1+t
          sum2=sum2+t*x(m)
   30   continue
c
c Check for impending underflow.
c
        if(abs(sum1).lt.tiny) then
          ierr=-k
          return
        end if
        alpha(k+1)=sum2/sum1
        beta(k+1)=sum1/sum0
        sum0=sum1
   40 continue
c
      return
      end
C
      subroutine fejer(nh,x,w)
      implicit real*8(a-h,o-z)
      dimension x(nh),w(nh)
c
c This routine generates the n-point Fejer quadrature rule.
c Again, this came from Walter Gautschi's package, but was
c modified by MWS to store only the quadrature points on
c the positive half of the usual interval, which is [-1,+1].
c
c         input:   nh  - the number of quadrature nodes at positive x.
c         output:  x,w - arrays of dimension -nh- holding the quadrature
c                        nodes and weights, respectively; the nodes
c                        are ordered increasingly
c
      n=nh+nh
      pi=4.0d+00*atan(1.0d+00)
      np1h=(n+1)/2
      fn=n
      do 10 k=1,nh
        x(nh+1-k)=cos(0.5d+00*(2*k-1)*pi/fn)
   10 continue
c     if(2*nh.ne.n) x(np1h)=0.0d+00
      do 30 k=1,np1h
        c1=1.0d+00
        c0=2.0d+00*x(k)*x(k)-1.0d+00
        t=2.0d+00*c0
        sum=c0/3.0d+00
        do 20 m=2,nh
          c2=c1
          c1=c0
          c0=t*c1-c2
          sum=sum+c0/(4*m*m-1)
   20   continue
        w(k)=2.0d+00*(1.0d+00-2.0d+00*sum)/fn
   30 continue
      return
      end
c
      subroutine RPboot(rtsaux,wtsaux,naux,NROOTS,XX,alpha,beta,RPprt)
      implicit real*8(a-h,o-z)
      dimension rtsaux(naux),wtsaux(naux),
     *          alpha(0:NROOTS-1),beta(0:NROOTS-1)
      parameter (mxquad=25)
      dimension t2exp(0:2*mxquad),RP(0:mxquad,0:mxquad),rpint(0:mxquad)
      logical RPprt
      parameter (zero=0.0d+00, one=1.0d+00)
c
c        This routine written by Mike Schmidt, Department of Chemistry,
c        Iowa State University, November 2005.  It does not control
c        the errors in alpha,beta recurrence coefficients well enough.
c        Running in Q.P. is decent, however, for NROOTS=13, XX=10.0.
c
c        Integrate all even powers of t times the Rys weight function,
c        using our auxiliary quadrature.
c
      do k=0,2*NROOTS
         t2exp(k) = zero
      enddo
      do i=1,naux
         t   = rtsaux(i)
         t2  = t*t
         dum = wtsaux(i)*exp(-XX*t2)
	 fac = one
	 do k=0,2*NROOTS
            t2exp(k) = t2exp(k) + fac*dum
            fac = fac*t2
         enddo
      enddo
c
      do k=0,NROOTS
         do j=0,NROOTS
            RP(j,k) = zero
         enddo
      enddo
c
c        by definition, the zero-th Rys Polynomial is unity
c
      RP(0,0) = one
c
c        generate the first Rys polynomial as a special case.
c          note: beta(0) = sqrt(pi/4*x) * derf(sqrt(x))
c
      rpint(0)= t2exp(0)
      alpha(0)= t2exp(1)/rpint(0)
      beta(0) = rpint(0)
      RP(0,1) = -alpha(0)
      RP(1,1) = one
      if(NROOTS.eq.1) return
c
      do k=1,NROOTS-1
c            numerical integration gives the recursion coefficients
         suma = zero
         sumb = zero
c            the multiplications and additions of alternating
c            signs in the loop below must be the roundoff problem.
         do i=0,k
            do j=0,k
               suma = suma + RP(i,k) * RP(j,k) * t2exp(i+j+1)
               sumb = sumb + RP(i,k) * RP(j,k) * t2exp(i+j)
            enddo
         enddo
         rpint(k) = sumb
         alpha(k) = suma/rpint(k)
         beta(k)  = rpint(k)/rpint(k-1)
c            upward recurrence generates the next Rys Polynomial
         do j=0,k
            RP(j,k+1) = -alpha(k)*RP(j,k) - beta(k)*RP(j,k-1)
         enddo
         do j=0,k
            RP(j+1,k+1) = RP(j+1,k+1) + RP(j,k)
         enddo
      enddo
c
c        print the Rys Polynomial coefficients
c
      if(RPprt) write(6,900) RP(0,NROOTS),(RP(j,NROOTS),2*j,j=1,NROOTS)
  900 format(1P,e15.7,7x,2(e15.7,'*t**',i2,1x)/
     *               (1P,3(e15.7,'*t**',i2,1x)))
c
      return
      end 
c
      subroutine RysGen(alpha,beta,NROOTS)
      implicit real*8(a-h,o-z)
      dimension alpha(0:*),beta(0:*)
      parameter (mxquad=25)
      dimension RP(0:mxquad,0:mxquad)
      parameter (zero=0.0d+00, one=1.0d+00)
c
      do k=0,NROOTS
         do j=0,NROOTS
            RP(j,k) = zero
         enddo
      enddo
c
c        by definition, the zero-th Rys Polynomial is unity
c
      RP(0,0) = one
c
c        then set the first Rys Polynomial
c
      RP(0,1) = -alpha(0)
      RP(1,1) = one
      if(NROOTS.eq.1) go to 200
c
c        upward recurrence generates the next Rys Polynomial, at k+1
c
      do k=1,NROOTS-1
         do j=0,k
            RP(j,k+1) = -alpha(k)*RP(j,k) - beta(k)*RP(j,k-1)
         enddo
         do j=0,k
            RP(j+1,k+1) = RP(j+1,k+1) + RP(j,k)
         enddo
      enddo
c
  200 continue
      write(6,900) RP(0,NROOTS),(RP(j,NROOTS),2*j,j=1,NROOTS)
      return
  900 format(1P,e15.7,7x,2(e15.7,'*t**',i2,1x)/
     *               (1P,3(e15.7,'*t**',i2,1x)))
      end
c
      subroutine dupuis(NROOTS,XX,rts,wts)
      implicit real*8(a-h,o-z)
      dimension rts(NROOTS),wts(NROOTS)
      COMMON /ROOT  / x,uf(13),wf(13),n
c
c        compute Rys Polynomial roots and weights,
c        using the original method of Michel Dupuis
c
      n= NROOTS
      x= XX
      if(n.gt.13) stop 'exceeded dimension limit for Dupuis code'
c
      call root6
      do i=1,NROOTS
         rts(i) = uf(i)
         wts(i) = wf(i)
      enddo
      return
      end
C*MODULE RYSPOL  *DECK ROOT6
      SUBROUTINE ROOT6
C
      implicit real*8(a-h,o-z)
C
      DIMENSION C(14,14),S(14,14),A(14),RT(14),R(13,13),W(13,13),FF(27)
C
      COMMON /ROOT  / XX,UF(13),WF(13),NROOTS
C
      PARAMETER (PT5=0.5D+00, ZERO=0.0D+00, ONE=1.0D+00, FOUR=4.0D+00)
C
C     I-TH ROOT OF THE J-TH RYS POLYNOMIAL IS RETURNED IN R(I,J)
C     WITH THE CORRESPONDING WEIGHT FACTOR IN W(I,J).   J=1,2,...,N
C     THIS VERSION USES CHRISTOFFEL FORMULA FOR WEIGHTS, AND IT IS
C     COMPLETELY GENERAL UP TO NROOTS=9.  SEE THE REFERENCE
C        "NUMERICAL INTEGRATION USING RYS POLYNOMIALS"
C                   H.F.KING AND M.DUPUIS
C              J.COMPUT.PHYS. 21, 144-165(1976)
C
      N=NROOTS
      X=XX
C
C     X LARGE (ASYMPTOTIC SOLUTION)              NROOTS=6,7,8,9
c     MWS: November 2005, I no longer think the asymptotic
c                         routines are very accurate!
C
c---  IF(X.GT.64.0D+00) THEN
      IF(X.GT.1000.0D+00) THEN
         CALL ROOTSA
         RETURN
      ENDIF
C
      IF(N.LT.2) N=2
      N1=N+1
      NN=N+N
      CALL RYSFUN(X,NN,FF)
C
      DO 11 I=1,N1
         DO 10 J=1,N1
            S(I,J)=FF(I+J-1)
   10    CONTINUE
   11 CONTINUE
      CALL RYSSMT(C,S,N1)
C
      DO 21 I=1,N
         DO 20 J=1,I
            W(I,J)= ZERO
            R(I,J)= ZERO
   20    CONTINUE
   21 CONTINUE
C
      WSUM=FF(1)
      W(1,1)=WSUM
      R(1,1)=FF(2)/WSUM
      DUM= SQRT(C(2,3)**2-FOUR *C(1,3)*C(3,3))
      R(1,2)= PT5*(-C(2,3)-DUM)/C(3,3)
      R(2,2)= PT5*(-C(2,3)+DUM)/C(3,3)
      IF(N.EQ.2) GO TO 70
      DO 25 I=3,N1
   25 RT(I)= ONE
      RT(1)=R(1,2)
      RT(2)=R(2,2)
C
      DO 60 K=3,N
         K1=K+1
         DO 30 I=1,K1
            A(I)=C(I,K1)
   30    CONTINUE
         CALL RYSNOD(A,RT,K)
         DO 50 I=1,K
            R(I,K)=RT(I)
   50    CONTINUE
   60 CONTINUE
C
   70 CONTINUE
      DO 150 K=2,N
         JMAX=K-1
         DO 140 I=1,K
            ROOT=R(I,K)
            DUM = ONE/FF(1)
            DO 110 J=1,JMAX
               J1=J+1
               POLY=C(J1,J1)
               DO 100 M=1,J
                  POLY=POLY*ROOT+C(J1-M,J1)
  100          CONTINUE
               DUM=DUM+POLY*POLY
  110       CONTINUE
            W(I,K) = ONE/DUM
  140    CONTINUE
  150 CONTINUE
C
C       GAMESS normally wants the reciprocal form to do an
C       actual integral, but here we want the root directly.
C
      DO 160 K=1,NROOTS
C--      DUM=R(K,NROOTS)
C--      UF(K)=DUM/(ONE-DUM)
         uf(k)=r(k,nroots)
         WF(K)=W(K,NROOTS)
  160 CONTINUE
      RETURN
      END
C*MODULE RYSPOL  *DECK RYSFUN
      SUBROUTINE RYSFUN(X,N,FF)
      implicit real*8(a-h,o-z)
      DIMENSION FF(27)
      DIMENSION XMAX(13)
      PARAMETER (ZERO=0.0D+00, PT5=0.5D+00, ONE=1.0D+00, TWO=2.0D+00)
C          THE FIRST FIVE ROOTS ARE DONE IN SPECIALIZED ROUTINES.
      DATA XMAX/5*0.0D+00,
     *          9000000.0D+00, 750000.0D+00, 150000.0D+00, 50000.0D+00,
     *           4*40000.0D+00/
C
      TOL   = 1.0D-14
      XX    = X+X
      FACMIN= XX
C
C        NOTE ADDED BY MWS, SEPTEMBER 1990:
C        VALUES X.GT.3200 AND NROOTS.EQ.6 GENERATE NEGATIVE ARGUMENTS
C        TO THE SQRT FUNCTION IN RYSSMT WHEN ONE USES E=10**(-35).
C        THE PROBLEM GOES AWAY, AND THE ROOTS AND WEIGHTS STABILIZE
C        FOR X OUT TO 5000 WHEN 10**(-40) TO 10**(-60) IS USED.
C        SINCE VAX UNDERFLOWS AT AROUND 10**(-38) WE HAVE TO BE MORE
C        CAREFUL THAN JUST SETTING E A FEW POWERS OF TEN DOWNWARD.
C        MICHEL'S ORIGINAL CODE IS SHOWN FIRST.   THE VALUE OF 360
C        IS WHERE IBMS GET INTO UNDERFLOWS (ABOUT 10**(-78)).  THE
C        1990 CHANGES INCLUDE THE COMPUTATION OF E IN THE ASYMPTOTIC
C        CODE A FEW LINES DOWN.
C
C--   E     = 1.0D-35
C--   IF(FACMIN.LT.360.0D+00) E=EXP(-X)
C--   IF(FACMIN.GT.80.0D+00) GO TO 100
C
      E = ZERO
      IF(FACMIN.LT.360.0D+00) E=EXP(-X)
      IF(FACMIN.GT.80.0D+00) GO TO 100
C
      TERM= ONE
      SUM = ONE
      FAC = N
      FAC = FAC+PT5
C
   10 FAC = FAC+ONE
      TERM= TERM*X/FAC
      SUM = SUM+TERM
      IF(FAC.LE.FACMIN) GO TO 10
      T=TERM
      S=SUM
      IF(T.GT.S*TOL) GO TO 10
C
      FAC=N+N+1
      FF(N+1)=SUM*E/FAC
      M=N-1
      FAC=M+M+1
C
   20 IF(M.LT.0) RETURN
      FF(M+1)=(E+XX*FF(M+2))/FAC
      M=M-1
      FAC=FAC-TWO
      GO TO 20
C
C     NOTE BY MWS, NOVEMBER 1997:  THIS CODE HAD BEEN SETTING -E- TO
C     1.0D-35 WITH DOWNWARDS ADJUSTMENTS TO 1.0D-45 ON MOST MACHINES.
C     MICHEL USUALLY USES A VALUE OF 1.0D-78 HERE, WHERE THE LOWER
C     END OF THE IBM MAINFRAME'S DYNAMIC RANGE LIES.  MIKE'S HACK AT
C     TRYING TO ACCOMODATE BOTH THE VAX, AND THE NEED FOR A LOWER GUESS
C     AT E IS SHOWN COMMENTED OUT.  INSTEAD WE WILL NOW TEST FOR THE
C     MAXIMUM X FOR WHICH 1D-78 WILL RETURN 9 DIGIT ACCURACY ROOTS AND
C     WEIGHTS FOR.  ANYONE WHO FINDS THEMSELVES HERE MIGHT BE ABLE TO
C     INTEGRATE ACCURATELY BY CHOOSING A STILL SMALLER VALUE OF -E-,
C     BUT SHOULD BE CAREFUL TO TEST THE NUMERICS OF THE SITUATION
C     INSTEAD OF JUST BLINDLY LOWERING THE THRESHHOLD FURTHER.
C
C  TO CONVERGE 9 SIGNIFICANT DIGITS (IS THIS ENOUGH???) FOR NROOTS=6
C          X    ROOT CONV.  WEIGHT CONV
C        5,000      -43        -45       (NUMBER IS POWER OF TEN FOR E)
C       10,000      -47        -48
C       20,000      -49        -50
C       50,000      -54        -54
C      100,000      -57        -58
C      200,000      -60        -61
C    1,000,000      -67        -69
C    9,000,000      -76        -78
C   10,000,000      -79        -79
C     SIMILAR TESTS SHOW THAT HIGHER POLYNOMIALS GO TO SMALLER X VALS
C
C     NO. OF SIGNIFICANT DIGITS ACCURACY USING E=1.0D-45 FOR NROOTS=6
C          X       ROOTS    WEIGHTS
C        5,000        9         9
C       10,000        8         8
C       15,000        7         5
C       20,000        4         4
C       30,000        3         2
C       40,000        2         2
C       50,000        1         0
C    THAT IS, FOR X=40000, ROOT 1 IS 2.442D-6 BUT SHOULD BE 2.469D-6,
C    SO ONLY TWO DIGITS OF ACCURACY IS ATTAINED.  AT X=58000, THE
C    GUESS OF 1D-45 FOR -E- FAILS TO FIND THE ROOTS AT ALL.
C
C     --- USE ASYMPTOTIC EXPANSION FOR LARGE ARGUMENTS ---
C
  100 CONTINUE
C---      IF(E.EQ.ZERO) THEN
C---         ESAVE = 1.0D-35
C---         DO 5 I=1,10
C---            E = ESAVE*0.1D+00
C---            IF(E.EQ.ZERO) GO TO 6
C---            ESAVE = E
C---    5    CONTINUE
C---    6    CONTINUE
C---         E = ESAVE
C---      END IF
C
      IF(E.EQ.ZERO) THEN
         NROOTS = N/2
         IF(X.GT.XMAX(NROOTS)) THEN
            WRITE(6,9010) NROOTS,X,XMAX(NROOTS)
            STOP 'error termination in RYSFUN'
         END IF
         E = 1.0D-78
      END IF
C
      A= SQRT(0.7853981633974483096156608D+00/X)
      TMAX= A*TOL/E
      TERM= ONE/XX
      SUM = TERM
      FAC = ONE
C
  110 FAC = FAC-TWO
      TERM= FAC*TERM/XX
      SUM = TERM+SUM
      T   = TERM
      IF(ABS(T).GT.TMAX) GO TO 110
C
      FF(1)= A-E*SUM
      FAC  = -ONE
      M    = 0
C
  120 IF(M.EQ.N) RETURN
      M      = M+1
      FAC    = FAC+TWO
      FF(M+1)= (FAC*FF(M)-E)/XX
      GO TO 120
C
 9010 FORMAT(1X,'RYSFUN: UNABLE TO FIND ROOTS OF RYS POLYNOMIAL.'/
     *       1X,'NROOTS=',I4,' X=',F20.5,' XMAX=',F20.5)
      END
C*MODULE RYSPOL  *DECK RYSNOD
      SUBROUTINE RYSNOD(A,RT,K)
      implicit real*8(a-h,o-z)
      DIMENSION A(14),RT(14)
      PARAMETER (ZERO=0.0D+00)
C
C     RETURNS IN RT(I) THE ITH ROOT OF A POLYNOMIAL OF ORDER K
C     WHOSE MTH COEFFICIENT IS STORED IN A(M+1).  IT IS ASSUMED
C     THAT THE INITIAL VALUES IN RT BRACKET THE FINAL VALUES.
C
      TOL=1.0D-11
      K1=K+1
      R2= ZERO
      P2=A(1)
      DO 100 M=1,K
         R1=R2
         P1=P2
         R2=RT(M)
         P2=A(K1)
         DO 10 I=1,K
            P2=P2*R2+A(K1-I)
   10    CONTINUE
         PROD=P1*P2
         IF(PROD.GE.ZERO) THEN
            WRITE(6,15) M,K
            STOP 'error termination in RYSNOD'
         END IF
         R5=R1
         P5=P1
         R6=R2
         P6=P2
   30    CONTINUE
         R3=R5
         P3=P5
         R4=R6
         P4=P6
         R =(R3*P4-R4*P3)/(P4-P3)
         DR=R4-R3
         DELTA=DR
         IF(ABS(DELTA).LT.TOL) GO TO 90
         DR=0.0625D+00*DR
         R5=R-DR
         IF(R5.LT.R3) R5=R3
         R6=R+DR
         IF(R6.GT.R4) R6=R4
         P5=A(K1)
         P6=P5
         DO 40 I=1,K
            P5=P5*R5+A(K1-I)
            P6=P6*R6+A(K1-I)
   40    CONTINUE
   45    CONTINUE
         PROD=P5*P6
         IF(PROD.LT. ZERO) GO TO 30
         PROD=P3*P5
         IF(PROD.GT. ZERO) GO TO 60
         R5=0.25D+00*R3+0.75D+00*R5
         P5=A(K1)
         DO 50 I=1,K
            P5=P5*R5 + A(K1-I)
   50    CONTINUE
         GO TO 45
C
   60    CONTINUE
         R6=0.25D+00*R4 + 0.75D+00*R6
         P6=A(K1)
         DO 70 I=1,K
            P6=P6*R6 + A(K1-I)
   70    CONTINUE
         GO TO 45
C
   90    CONTINUE
         RT(M)=R
  100 CONTINUE
      RETURN
C
   15 FORMAT(//1X,'RYSNOD: ROOT NUMBER ',I3,
     *            ' WAS NOT FOUND FOR POLYNOMIAL OF ORDER ',I3,
     *            ' ON NODE',I4//)
      END
C*MODULE RYSPOL  *DECK RYSSMT
      SUBROUTINE RYSSMT(C,S,N)
      implicit real*8(a-h,o-z)
      DIMENSION C(14,14),S(14,14),V(14),Y(14)
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00)
C
C     ROUTINE RETURNS AN N BY N TRIANGULAR MATRIX C SUCH THAT
C     C(TRANSPOSE)SC=I,  WHERE I IS AN N BY N IDENTITY MATRIX.
C
      DO 10 I=1,N
      DO 10 J=1,I
   10 C(I,J)= ZERO
C
      DO 100 J=1,N
         KMAX=J-1
         FAC=S(J,J)
         IF(KMAX.EQ.0) GO TO 60
C
         DO 20 K=1,KMAX
            V(K)= ZERO
            Y(K)=S(K,J)
   20    CONTINUE
         DO 50 K=1,KMAX
            DOT= ZERO
            DO 30 I=1,K
               DOT=C(I,K)*Y(I)+DOT
   30       CONTINUE
            DO 40 I=1,K
               V(I)=V(I)-DOT*C(I,K)
   40       CONTINUE
            FAC=FAC-DOT*DOT
   50    CONTINUE
C
   60    CONTINUE
C
C               THE ORIGINAL CODE LEADS TO DIVIDE BY
C               ZERO ON VAX COMPUTERS FOR X>3100 OR SO.
C--      FAC=ONE/SQRT(FAC)
C
         IF(FAC.GT.ZERO) THEN
            FAC=  ONE/SQRT(FAC)
         ELSE
            FAC = 1.0D+25
         END IF
C
         C(J,J)=FAC
         IF(KMAX.EQ.0) GO TO 100
         DO 70 K=1,KMAX
            C(K,J)=FAC*V(K)
   70    CONTINUE
  100 CONTINUE
      RETURN
      END
C*MODULE RYSPOL  *DECK ROOTSA
      SUBROUTINE ROOTSA
C
      implicit real*8(a-h,o-z)
C
      COMMON /ROOT  / X,U(13),W(13),NROOTS
C
C     COMPUTE THE ROOTS AND WEIGHTS OF THE RYS POLYNOMIALS OF ORDER
C     NROOTS=1 UP TO NROOTS=13 FOR ASYMPTOTICALLY LARGE -X- VALUES.
C
C     AUTHOR: JOSE M. SIERRA, SYNSTAR CORPORATION, SEPTEMBER 2004
C
C     THE CODE IS CORRECT FOR ANY NROOTS, BUT MAY BE A LITTLE BIT SLOWER
C     THAN THE ORIGINAL SMALL NROOTS VALUE EXPANSIONS OF DUPUIS AND RYS.
C
      PARAMETER (MR=13, MX=(MR*(MR+1))/2)
      DIMENSION AR(MX),AW(MX)
C
      PARAMETER (ONE=1.0D+00, PI4=7.85398163397448D-01)
C
C    NEW VALUES (OBTAINED FROM 'ANALYTICAL' POLYNOMIALS)
C
C NROOTS=1
      DATA (AR(I),I= 1, 1)/  5.00000000000000D-01/
C NROOTS=2
      DATA (AR(I),I= 2, 3)/  2.75255128608411D-01, 2.72474487139159D+00/
C NROOTS=3
      DATA (AR(I),I= 4, 6)/  1.90163509193488D-01, 1.78449274854325D+00,
     3 5.52534374226326D+00/
C NROOTS=4
      DATA (AR(I),I= 7,10)/  1.45303521503317D-01, 1.33909728812636D+00,
     4 3.92696350135829D+00, 8.58863568901203D+00/
C NROOTS=5
      DATA (AR(I),I=11,15)/  1.17581320211778D-01, 1.07456201243690D+00,
     5 3.08593744371755D+00, 6.41472973366203D+00, 1.18071894899717D+01/
C NROOTS=6
      DATA (AR(I),I=16,21)/  9.87470140684812D-02, 8.98302834569617D-01,
     6 2.55258980266817D+00, 5.19615253005447D+00, 9.12424803753117D+00,
     6 1.51299597811081D+01/
C NROOTS=7
      DATA (AR(I),I=22,28)/  8.51154429975940D-02, 7.72137920042777D-01,
     7 2.18059188845046D+00, 4.38979288673077D+00, 7.55409132610392D+00,
     7 1.19899930398239D+01, 1.852827749585248D+01/
C NROOTS=8
      DATA (AR(I),I=29,36)/  7.47918825968183D-02, 6.77249087649289D-01,
     8 1.90511363503143D+00, 3.80947636148493D+00, 6.48314542862698D+00,
     8 1.00933236752216D+01, 1.49726270884262D+01, 2.19842728409627D+01/
C NROOTS=9
      DATA (AR(I),I=37,45)/  6.67022309581943D-02, 6.03236357081748D-01,
     9 1.69239507979318D+00, 3.36917627024335D+00, 5.69442334295712D+00,
     9 8.76975673027048D+00, 1.27718253548667D+01, 1.80465054677307D+01,
     9 2.54859791660990D+01/
C NROOTS=10
      DATA (AR(I),I=46,55)/  6.01920631495866D-02, 5.43867500294645D-01,
     A 1.52294410540445D+00, 3.02251337645400D+00, 5.08490775008944D+00,
     A 7.77743923154111D+00, 1.12081302043343D+01, 1.55611633321971D+01,
     A 2.11938920962985D+01, 2.90249503402365D+01/
C NROOTS=11
      DATA (AR(I),I=56,66)/  5.48398695788185D-02, 4.95174123350356D-01,
     B 1.38465574008460D+00, 2.74191994010634D+00, 4.59773770048869D+00,
     B 6.99939746952188D+00, 1.00189082759681D+01, 1.37693058660907D+01,
     B 1.84411196809850D+01, 2.44019612423853D+01, 3.25949800914407D+01/
C NROOTS=12
      DATA (AR(I),I=67,78)/  5.03618891172940D-02, 4.54506681563779D-01,
     C 1.26958994010396D+00, 2.50984809722837D+00, 4.19841564490072D+00,
     C 6.36997538792228D+00, 9.07543423132649D+00, 1.23904479630736D+01,
     C 1.64321950885250D+01, 2.13967559356126D+01, 2.76611087800458D+01,
     C 3.61913603605784D+01/
C NROOTS=13
      DATA (AR(I),I=79,91)/  4.65600832450248D-02, 4.20027406401213D-01,
     D 1.17231077327777D+00, 2.31454086434948D+00, 3.86458503822818D+00,
     D 5.84873481130518D+00, 8.30455348999401D+00, 1.12857509934872D+01,
     D 1.48709603775928D+01, 1.91809194854928D+01, 2.44166923332210D+01,
     D 3.09639382746126D+01, 3.98104260687924D+01/
C
C    VALUES (OBTAINED USING GAMESS FOR X=128)
C
C NROOTS=1
      DATA (AW(I),I= 1, 1)/  1.00000000000000D+00/
C NROOTS=2
      DATA (AW(I),I= 2, 3)/  1.00000000000000D+00, 9.17517095361369D-02/
C NROOTS=3
      DATA (AW(I),I= 4, 6)/  1.00000000000000D+00, 1.77231492083829D-01,
     3 5.11156880411248D-03/
C NROOTS=4
      DATA (AW(I),I= 7,10)/  1.00000000000000D+00, 2.34479815323517D-01,
     4 1.92704402415764D-02, 2.25229076750736D-04/
C NROOTS=5
      DATA (AW(I),I=11,15)/  1.00000000000000D+00, 2.70967405960535D-01,
     5 3.82231610015404D-02, 1.51614186862443D-03, 8.62130526143657D-06/
C NROOTS=6
      DATA (AW(I),I=16,21)/  1.00000000000000D+00, 2.93934096090658D-01,
     6 5.82333758247259D-02, 4.40676137506590D-03, 9.67436984518118D-05,
     6 2.99985433527357D-07/
C NROOTS=7
      DATA (AW(I),I=22,28)/  1.00000000000000D+00, 3.08166679684882D-01,
     7 7.73002176482838D-02, 8.85783821383684D-03, 4.00679107517203D-04,
     7 5.32198268804275D-06, 9.73632251525836D-09/
C NROOTS=8
      DATA (AW(I),I=29,36)/  1.00000000000000D+00, 3.16676745499791D-01,
     8 9.45695047037012D-02, 1.45338752009040D-02, 1.05196985297312D-03,
     8 3.06000643175543D-05, 2.61894643174802D-07, 2.99562944518145D-10/
C NROOTS=9
      DATA (AW(I),I=37,45)/  1.00000000000000D+00, 3.21370607772183D-01,
     9 1.09793264917933D-01, 2.10330354860907D-02, 2.13096958967959D-03,
     9 1.03597922686261D-04, 2.04310479048232D-06, 1.18109769258849D-08,
     9 8.83317751281791D-12/
C
C     X LARGE (ASYMPTOTIC SOLUTION)              NROOTS=1,..,13
C
      N0=(NROOTS*NROOTS-NROOTS)/2
      R = ONE/X
      FW= SQRT(PI4*R)
      SW=-FW
      DO 110 M=1,NROOTS
         TMS = AR(N0+M)*R
         U(M)= TMS/(ONE-TMS)
         W(M)= AW(N0+M)*FW
         SW= SW+W(M)
  110 CONTINUE
      W(1)= FW-SW
      RETURN
      END
c
      subroutine rysset
      implicit real*8(a-h,o-z)
      PARAMETER (MR=13, MR2=2*MR)
      dimension rts(MR2),wts(MR2),wrk(MR2),
     *          alpha(0:MR2-1),beta(0:MR2-1)
      common /ryspar/ rtsasy(MR,MR),wtsasy(MR,MR)
c
      PI=3.14159265358979323846264338327950288419716939937510D+00
      zero= 0.0D+00
      two = 2.0D+00
      eps = 1.0d-14
c
c     set up asymptotic root computation, by generating
c     the roots of Hermite polynomials of order 2N.
c
      do i=1,MR
         n=i+i
         do j=0,n-1
            alpha(j) = zero
         enddo
         beta(0)=sqrt(pi)
         do j=1,n-1
            beta(j) = j/two
         end do
         call gauss(n,alpha,beta,eps,rts,wts,ierr,wrk)
         do j=1,i
            rtsasy(j,i) = rts(i+j)*rts(i+j)
            wtsasy(j,i) = wts(i+j)
         enddo
      enddo
      return
      end
c
      subroutine rysasy(NROOTS,XX,rts,wts)
      implicit real*8(a-h,o-z)
      dimension rts(NROOTS),wts(NROOTS)
      PARAMETER (MR=13)
      common /ryspar/ rtsasy(MR,MR),wtsasy(MR,MR)
      parameter (one=1.0d+00)
C
C     The roots and weights for X = big are given by:
c       t*t = s*s/XX    w = v/sqrt(XX)
c     where s and v are roots and weights of the
c     Hermite polynomials of order 2*NROOTS.
c     See Ishida-san's second paper for details, although he appears
c     to have been a bit too conservative about when to go asymptotic.
C
      factr = one/XX
      factw = sqrt(factr)
      DO I=1,NROOTS
         rts(I)= factr * rtsasy(I,NROOTS)
         wts(I)= factw * wtsasy(I,NROOTS)
      enddo
      RETURN
      END
c
      subroutine Raiz128(S,I,N,X,M)
      implicit real*8(a-h,o-z)
      PARAMETER (MPO=13, M2D=128, M2L=M2D-1)
      common /jms128/ rtsexp(0:M2L,MPO,MPO),wtsexp(0:M2L,MPO,MPO)
c
C                     2       m             k
C Calculate  S (x) = t (x) = Sum [sin(k) * x ]
C             in      in     k=0
C
      if(M.gt.0) stop "Sierra code is truncated"
c
      L= MIN(M,M2D)-1
      S= rtsexp(L,I,N)
      DO K=L-1, 0,-1
         S= S*X+rtsexp(K,I,N)
      ENDDO
C
      RETURN
      END
c
      subroutine Peso128(W,I,N,X,M)
      implicit real*8(a-h,o-z)
      PARAMETER (MPO=13, M2D=128, M2L=M2D-1)
      common /jms128/ rtsexp(0:M2L,MPO,MPO),wtsexp(0:M2L,MPO,MPO)
C
C                     m             k
C Calculate  W (x) = Sum [win(k) * x ]
C             in     k=0
C
      L= MIN(M,M2D)-1
      W= wtsexp(L,I,N)
      DO K=L-1, 0,-1
         W= W*X+wtsexp(K,I,N)
      ENDDO
C
      RETURN
      END
      subroutine setRaiz128
      return
      end
      subroutine setPeso128
      return
      end
