      real function scphi(x)
      implicit none
      real x
      double precision myderf
      DOUBLE PRECISION FITPAD(29),FITFUN
      COMMON/HCFITD/FITPAD,FITFUN
      double precision p3g
      double precision y1, y2, y3, y4 ,y5, y6, y7, y8, y9, y10,
     & y11, y12, w1, w2, w3
      double precision s1, ss, sm, s2, ys, sig, ds, xd, 
sigs, sigg
      double precision p(4)
      double precision x1(4), x2(4), x3(4)
      double precision y(4)
      double precision sxi(12), dsi(3), ssi(3), sxim, 
dsxim, agt0
      integer ni, i, ind, ind1, ind2,j
      common/rims/ smin, smax
      double precision smin, smax
      double precision ssj(5), dsj(5), ssjt, dsjt, der
      vector xi(10)
      vector si(3)
      der=2.40385
      ni=0
      do i=1,3
        if (si(i).ne.-1000) then
          ni=ni+1
        endif
      enddo
      ds=1.5/120.0*180/3.1415927
      xd=x
      y1=abs(fitpad(1))
      y2=abs(fitpad(2))
      y3=abs(fitpad(3))
      y4=abs(fitpad(4))
      y5=abs(fitpad(5))
      y6=abs(fitpad(6))
      y7=abs(fitpad(7))
      y8=abs(fitpad(8))
      y9=abs(fitpad(9))
      y10=abs(fitpad(10))
      y11=abs(fitpad(11))
      y12=abs(fitpad(29))
      s1=abs(fitpad(12))
      ss=abs(fitpad(13))
      sm=abs(fitpad(14))
      s2=abs(fitpad(15))
      ys=abs(fitpad(16))
      sig=abs(fitpad(17))
      sigs=abs(fitpad(18))
      sigg=abs(fitpad(19)) c sigg=sigs
      dsi(1)=abs(fitpad(20)/2.0d0)
      dsi(2)=abs(fitpad(21)/2.0d0)
      dsi(3)=abs(fitpad(22)/2.0d0)
      ssi(1)=si(1)+fitpad(23)
      ssi(2)=si(2)+fitpad(24)
      ssi(3)=si(3)+fitpad(25)
      w1=abs(fitpad(26))
      w2=abs(fitpad(27))
      w3=abs(fitpad(28)) c
      smin=s1
      smax=s2
      ind=0
      if (((abs(smin-0).lt.10).or.
     & (abs(smin-120).lt.10).or.
     & (abs(smin-240).lt.10).or.
     & (abs(smin-360).lt.10))) then
        ind=ind+1
        dsj(ind)=0.001
        ssj(ind)=smin+der
      endif
      if (((abs(smax-0).lt.10).or.
     & (abs(smax-120).lt.10).or.
     & (abs(smax-240).lt.10).or.
     & (abs(smax-360).lt.10))) then
        ind=ind+1
        dsj(ind)=0.001
        ssj(ind)=smax-der
      endif
      ni=ni+ind
      do i=1,3
        ind=ind+1
        dsj(ind)=dsi(i)
        ssj(ind)=ssi(i)
      enddo
      do i=1,ni-1
        do j=i+1,ni
          if (ssj(i).gt.ssj(j)) then
            ssjt=ssj(i)
            dsjt=dsj(i)
            ssj(i)=ssj(j)
            dsj(i)=dsj(j)
            ssj(j)=ssjt
            dsj(j)=dsjt
          endif
        enddo
      enddo c
      FITFUN=0 c ds=abs(fitpad(19))/120.0*180/3.1415927 c 
if (xd.lt.sm) then
        x1(1)=s1
        x1(2)=xi(2)
        x1(3)=xi(3)
        x1(4)=ss
        y(1)=y1
        y(2)=y2
        y(3)=y3
        y(4)=y4
        call m4(x1,y,p) c 
FITFUN=p3g(xd,s1,sig,p)-p3g(xd,ss-ds,sigs,p)
        ind=1
        ind1=ind
        ind2=0
        sxi(ind1)=s1
        do i=1,ni
          if 
((ssj(i)-dsj(i).gt.s1).and.(ssj(i)-dsj(i).lt.ss-ds)) then
            ind=ind+1
            sxi(ind)=ssj(i)-dsj(i)
          endif
          if 
((ssj(i)+dsj(i).gt.s1).and.(ssj(i)+dsj(i).lt.ss-ds)) then
            ind=ind+1
            sxi(ind)=ssj(i)+dsj(i)
          endif
          if 
((s1.gt.ssj(i)-dsj(i)).and.(s1.lt.ssj(i)+dsj(i))) then
            ind1=ind
          endif
          if ((ss-ds.gt.ssj(i)-dsj(i)).and.
     & (ss-ds.lt.ssj(i)+dsj(i))) then
            ind2=ind
          endif
        enddo
        if (ind2.eq.0) then
          ind=ind+1
          ind2=ind
          sxi(ind2)=ss-ds-w1
        endif
        do i=ind1,ind2,2
          if (i.eq.ind1.and.i+1.eq.ind2) then
            
FITFUN=FITFUN+p3g(xd,sxi(i),sig,p)-p3g(xd,sxi(i+1),sigs,p)
          else if (i.eq.ind1.and.i+1.ne.ind2) then
            
FITFUN=FITFUN+p3g(xd,sxi(i),sig,p)-p3g(xd,sxi(i+1),sigg,p) 
c sxim=(sxi(i+2)+sxi(i+1))/2 c 
dsxim=(sxi(i+2)-sxi(i+1))/2/3 c 
FITFUN=FITFUN+(p3g(xd,sxim,sig,p)-p3g(xd,sxim,sigs,p))/2/2* 
c & (myderf(((xd-(sxim-dsxim))/sqrt(2.0)/sigg))- c & 
myderf(((xd-(sxim+dsxim))/sqrt(2.0)/sigg)))
          else if (i.ne.ind1.and.i+1.eq.ind2) then
            
FITFUN=FITFUN+p3g(xd,sxi(i),sigg,p)-p3g(xd,sxi(i+1),sigs,p) 
c sxim=(sxi(i)+sxi(i-1))/2 c dsxim=(sxi(i)-sxi(i-1))/2/3 c 
FITFUN=FITFUN+(p3g(xd,sxim,sig,p)-p3g(xd,sxim,sigs,p))/2/2* 
c & (myderf(((xd-(sxim-dsxim))/sqrt(2.0)/sigg))- c & 
myderf(((xd-(sxim+dsxim))/sqrt(2.0)/sigg)))
          else
            
FITFUN=FITFUN+p3g(xd,sxi(i),sigg,p)-p3g(xd,sxi(i+1),sigg,p) 
c sxim=(sxi(i)+sxi(i-1))/2 c dsxim=(sxi(i)-sxi(i-1))/2/3 c 
FITFUN=FITFUN+(p3g(xd,sxim,sig,p)-p3g(xd,sxim,sigs,p))/2/2* 
c & (myderf(((xd-(sxim-dsxim))/sqrt(2.0)/sigg))- c & 
myderf(((xd-(sxim+dsxim))/sqrt(2.0)/sigg))) c 
sxim=(sxi(i+2)+sxi(i+1))/2 c dsxim=(sxi(i+2)-sxi(i+1))/2/3 
c 
FITFUN=FITFUN+(p3g(xd,sxim,sig,p)-p3g(xd,sxim,sigs,p))/2/2* 
c & (myderf(((xd-(sxim-dsxim))/sqrt(2.0)/sigg))- c & 
myderf(((xd-(sxim+dsxim))/sqrt(2.0)/sigg)))
          endif
        enddo c
        x2(1)=ss
        x2(2)=xi(5)
        x2(3)=xi(6)
        x2(4)=sm
        y(1)=y5
        y(2)=y6
        y(3)=y7
        y(4)=y8
        call m4(x2,y,p) c 
FITFUN=FITFUN+p3g(xd,ss+ds,sigs,p)-p3g(xd,sm,sigs,p)
        ind=1
        ind1=ind
        ind2=0
        sxi(ind1)=ss+ds+w2
        do i=1,ni
          if 
((ssj(i)-dsj(i).gt.ss+ds).and.(ssj(i)-dsj(i).lt.sm)) then
            ind=ind+1
            sxi(ind)=ssj(i)-dsj(i)
          endif
          if 
((ssj(i)+dsj(i).gt.ss+ds).and.(ssj(i)+dsj(i).lt.sm)) then
            ind=ind+1
            sxi(ind)=ssj(i)+dsj(i)
          endif
          if ((ss+ds.gt.ssj(i)-dsj(i)).and.
     & (ss+ds.lt.ssj(i)+dsj(i))) then
            ind1=ind
          endif
          if 
((sm.gt.ssj(i)-dsj(i)).and.(sm.lt.ssj(i)+dsj(i))) then
            ind2=ind
          endif
        enddo
        if (ind2.eq.0) then
          ind=ind+1
          ind2=ind
          sxi(ind2)=sm-w3/2.0d0
        endif
        do i=ind1,ind2,2
          if (i.eq.ind1.and.i+1.eq.ind2) then
            
FITFUN=FITFUN+p3g(xd,sxi(i),sigs,p)-p3g(xd,sxi(i+1),sigs,p)
          else if (i.eq.ind1.and.i+1.ne.ind2) then
            
FITFUN=FITFUN+p3g(xd,sxi(i),sigs,p)-p3g(xd,sxi(i+1),sigg,p)
          else if (i.ne.ind1.and.i+1.eq.ind2) then
            
FITFUN=FITFUN+p3g(xd,sxi(i),sigg,p)-p3g(xd,sxi(i+1),sigs,p)
          else
            
FITFUN=FITFUN+p3g(xd,sxi(i),sigg,p)-p3g(xd,sxi(i+1),sigg,p)
          endif
        enddo c else
        x3(1)=sm
        x3(2)=xi(8)
        x3(3)=xi(9)
        x3(4)=s2
        y(1)=y12
        y(2)=y9
        y(3)=y10
        y(4)=y11
        call m4(x3,y,p) c 
FITFUN=FITFUN+p3g(xd,sm,sigs,p)-p3g(xd,s2,sig,p) cc 
FITFUN=p3g(xd,s2,-sig,p)
        ind=1
        ind1=ind
        ind2=0
        sxi(ind1)=sm+w3/2.0d0
        do i=1,ni
          if 
((ssj(i)-dsj(i).gt.sm).and.(ssj(i)-dsj(i).lt.s2)) then
            ind=ind+1
            sxi(ind)=ssj(i)-dsj(i)
          endif
          if 
((ssj(i)+dsj(i).gt.sm).and.(ssj(i)+dsj(i).lt.s2)) then
            ind=ind+1
            sxi(ind)=ssj(i)+dsj(i)
          endif
          if 
((sm.gt.ssj(i)-dsj(i)).and.(sm.lt.ssj(i)+dsj(i))) then
            ind1=ind
          endif
          if 
((s2.gt.ssj(i)-dsj(i)).and.(s2.lt.ssj(i)+dsj(i))) then
            ind2=ind
          endif
        enddo
        if (ind2.eq.0) then
          ind=ind+1
          ind2=ind
          sxi(ind2)=s2
        endif
        do i=ind1,ind2,2
          if (i.eq.ind1.and.i+1.eq.ind2) then
            
FITFUN=FITFUN+p3g(xd,sxi(i),sigs,p)-p3g(xd,sxi(i+1),sig,p)
          else if (i.eq.ind1.and.i+1.ne.ind2) then
            
FITFUN=FITFUN+p3g(xd,sxi(i),sigs,p)-p3g(xd,sxi(i+1),sigg,p)
          else if (i.ne.ind1.and.i+1.eq.ind2) then
            
FITFUN=FITFUN+p3g(xd,sxi(i),sigg,p)-p3g(xd,sxi(i+1),sig,p)
          else
            
FITFUN=FITFUN+p3g(xd,sxi(i),sigg,p)-p3g(xd,sxi(i+1),sigg,p)
          endif
        enddo c endif
        FITFUN=FITFUN+agt0(ys)/2.0d0*
     & (myderf(((xd-(ss-ds))/sqrt(2.0)/sigs))-
     & myderf(((xd-(ss+ds))/sqrt(2.0)/sigs)))
        FITFUN=FITFUN+agt0(0.0d0)-agt0(0.0d0)/2.0d0*
     & (myderf(((xd-s1)/sqrt(2.0)/sig))-
     & myderf(((xd-s2)/sqrt(2.0)/sig))) c if 
(FITFUN.gt.1.0d-5) then c FITFUN = 
FITFUN/(1.0d0-exp(-FITFUN)) c else c FITFUN = 
2.0d0/(2.0d0-FITFUN) c endif
      if (FITFUN.gt.1000) then
        FITFUN = 1.0d0
      endif
      scphi=FITFUN
      return
      end

      double precision function agt0(a)
      implicit none
      double precision a, w, a0
      vector pix(3) 
      vector wp(1)
      if (wp(1).eq.1) then
        w=pix(1)
        a0=pix(2)
        agt0 = (w*a0+a)/(1.0d0-(1.0d0-w)*exp(-a))
      else
        agt0 = a
      endif
      return
      end

      double precision function p3g(x,xr,sr,p)
      implicit none
      double precision x,xr,sr,p(4)
      double precision er, der
      double precision p3g0, p3gl
      common/rims/ smin, smax
      double precision smin, smax
      er=0.821691
      der=2.40385
      if (((abs(smin-0).lt.10).or.
     & (abs(smin-120).lt.10).or.
     & (abs(smin-240).lt.10).or.
     & (abs(smin-360).lt.10)).and.
     & ((xr-smin).lt.der)) then
        p3g=p3gl(x,xr,sr,er,smin,der,p)
        return
      endif
      if (((abs(smax-0).lt.10).or.
     & (abs(smax-120).lt.10).or.
     & (abs(smax-240).lt.10).or.
     & (abs(smax-360).lt.10)).and.
     & ((smax-xr).lt.der)) then
        p3g=p3gl(x,xr,sr,er,smax,-der,p)
        return
      endif
      p3g=p3g0(x,xr,sr,p)
      return
      end
      double precision function p3g0(x,xr,sr,p)
      implicit none
      double precision x,xr,sr,p(4),p0,p1,p2,p3
      double precision srr, arg, PI
      double precision myderf
      PI=3.1415926535
      p0=p(1)
      p1=p(2)
      p2=p(3)
      p3=p(4)
      srr=0.5+(sr-0.5)*(
     & myderf((x-(xr-5*sr))/sqrt(2.0)/sr)-
     & myderf((x-(xr+5*sr))/sqrt(2.0)/sr)
     & )/2.0
      arg=(x-xr)/sqrt(2.0)/srr
      p3g0=exp(-arg**2)*srr/sqrt(2.0*PI)*
     & (p1+p2*(x+xr)+p3*(x*x+x*xr+xr*xr+2*srr*srr))+
     & (1+myderf((arg)))/2*
     & (p0+x*(p1+x*(p2+p3*x))+(p2+3*p3*x)*srr*srr)
      return
      end
      double precision function p3gl(x,xr,sr,er,xer,der,p)
      implicit none
      double precision x,xr,sr,er,xer,der,
     & p(4),p0,p1,p2,p3
      double precision earg, arg, PI, srr
      double precision myderf
      double precision part1, part2, part3, part4, part5
      PI=3.1415926535
      p0=p(1)
      p1=p(2)
      p2=p(3)
      p3=p(4)
      srr=0.5+(sr-0.5)*(
     & myderf((x-(xr-5*sr))/sqrt(2.0)/sr)-
     & myderf((x-(xr+5*sr))/sqrt(2.0)/sr)
     & )/2.0 c srr=max(sr,0.001)
      arg=(x - xr)/(sqrt(2.0)*srr)
      earg=exp(-arg**2)
      p3gl=(2*(der*er*(sqrt(2*PI)*(p0 + p2*(srr**2 + x**2) 
+
     & x*(p1 + p3*(3*srr**2 + x**2))) +
     & srr*(p1 + p2*(x + xr) +
     & p3*(2*srr**2 + x**2 + x*xr + xr**2))*earg) -
     & (-1 + er)*(sqrt(2*PI)*
     & (p1*srr**2 + 3*p3*srr**4 + p0*x + 3*p2*srr**2*x + 
p1*x**2 +
     & 6*p3*srr**2*x**2 + p2*x**3 + p3*x**4 -
     & (p0 + p2*(srr**2 + x**2) +
     & x*(p1 + p3*(3*srr**2 + x**2)))*xer) +
     & srr*(p0 + p1*(x - xer + xr) +
     & p2*(2*srr**2 + x**2 - x*xer + x*xr - xer*xr + xr**2) 
+
     & p3*(x**2*(x - xer) + x*(x - xer)*xr + (x - 
xer)*xr**2 +
     & xr**3 + srr**2*(5*x - 2*xer + 3*xr)))*earg)) +
     & sqrt(2*PI)*(-(der*er*(p0 + p2*(srr**2 + x**2) +
     & x*(p1 + p3*(3*srr**2 + x**2)))) +
     & (-1 + er)*(p1*srr**2 + 3*p3*srr**4 + p0*x +
     & 3*p2*srr**2*x + p1*x**2 +
     & 6*p3*srr**2*x**2 + p2*x**3 + p3*x**4 -
     & (p0 + p2*(srr**2 + x**2) +
     & x*(p1 + p3*(3*srr**2 + x**2)))*xer))*
     & (1 - myderf(arg)))/(2.*der*sqrt(2*PI))
      return
      end
      subroutine m4(x,y,p)
      implicit none
      double precision x(4), y(4), p(4)
      double precision d
      double precision x1, x2, x3, x4
      double precision y1, y2, y3, y4
      double precision agt0
      x1=x(1)
      x2=x(2)
      x3=x(3)
      x4=x(4)
      y1=agt0(abs(y(1)))
      y2=agt0(abs(y(2)))
      y3=agt0(abs(y(3)))
      y4=agt0(abs(y(4)))
      d=(x1-x2)*(x1-x3)*(x2-x3)*(x1-x4)*(x2-x4)*(x3-x4)
      p(1)=(x1*(x1-x3)*x3*(x1-x4)*(x3-x4)*x4*y2+
     & x2*x4**2*(-x3**3*y1+x3**2*x4*y1+x1**2*(x1-x4)*y3)+
     & x1**2*x2*x3**2*(-x1+x3)*y4+
     & x2**3*(x4*(-x3**2*y1+x3*x4*y1+x1*(x1-x4)*y3)+
     & x1*x3*(-x1+x3)*y4)+
     & x2**2*(x1*x4*(-x1**2+x4**2)*y3+
     & x3**3*(x4*y1-x1*y4)+x3*(-x4**3*y1+x1**3*y4)))/d
      p(2)=(x1**2*(x1-x4)*x4**2*(y2-y3)+
     & x3**3*(x4**2*(y1-y2)+x1**2*(y2-y4))+
     & x2**2*(x4**3*(y1-y3)+x1**3*(y3-y4)+x3**3*(-y1+y4))+
     & x3**2*(x4**3*(-y1+y2)+x1**3*(-y2+y4))+
     & 
x2**3*(x4**2*(-y1+y3)+x3**2*(y1-y4)+x1**2*(-y3+y4)))/d
      p(3)=(-x1*(x1-x4)*x4*(x1+x4)*(y2-y3)+
     & x3*(x4**3*(y1-y2)+x1**3*(y2-y4))+
     & x3**3*(-x4*y1-x1*y2+x4*y2+x1*y4)+
     & x2**3*(x4*y1+x1*y3-x4*y3-x1*y4+x3*(-y1+y4))+
     & x2*(x4**3*(-y1+y3)+x3**3*(y1-y4)+x1**3*(-y3+y4)))/d
      p(4)=(x1*(x1-x4)*x4*(y2-y3)+
     & x3**2*(x4*y1+x1*y2-x4*y2-x1*y4)+
     & x2**2*(-x4*y1-x1*y3+x4*y3+x3*(y1-y4)+x1*y4)+
     & x2*(x4**2*(y1-y3)+x1**2*(y3-y4)+x3**2*(-y1+y4))+
     & x3*(x4**2*(-y1+y2)+x1**2*(-y2+y4)))/d
      end
