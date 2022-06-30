macro mylib mylibinit='imp'
if ([mylibinit].eq.'imp') then
  gl/imp mylibinit
  mylibinit='no'
endif
*mylibinit='yes'
if ([mylibinit].ne.'yes') then
  gl/cre mylibinit yes
  appli comis quit
  !file mydspin2.sl
  !file mydspcd2.sl
  !file mydspps2.sl
  !file mydspap2.sl
  !file mydspap1.sl
  !file mydspin1.sl
  !file mydspps1.sl
  !file mydspcd1.sl
  !file mydsimps.sl
  !file myderf.sl
*  !file accmappl.sl
  !file mhdiff.sl
  !file mydgamma.sl
  !file accmappl_blocks.sl
  !file sndvcopy.sl
  !file /personal/l3-1-205-1/konctbel/libs/mylib.sl
  quit
  exec /personal/l3-1-205-1/konctbel/kskl/k#mylib
else
  mess MyLib already is set
endif
x=$call('test.f(0.)')
exec hsigma test = &test( array(101,-1.#1.) )
x=$call('mydsimps.sl(test,-1.,1.,100)')
mess [x]
*x=$call('myderf.sl(0.0)')
ve/del test
mess $call('myderf.sl(0.0)')
mess $call('mydgamma.sl(1.0)') 
return

macro cuts fun=zf mylibinit='imp'
exec mapcal#mylib mylibinit=[mylibinit]

cut 0 -
cut $1 nc.ge.2.and.np.ge.nc
cut $2 energy(1)/beam<0.8.and.energy(2)/beam<0.6
cut $3 act.eq.0
dphi=$sigma(3/180*3.1415927)
cut $4 abs(abs(phi(1)-phi(2))-3.1415927)<[dphi]
dtheta=$sigma(5/180*3.1415927)
cut $5 abs(theta(1)+theta(2)-3.1415927)<[dtheta]
cut $6 abs(z0(1)-z0(2))<2.and.abs(z0(1)+z0(2))<14
cut $7 (e1/14.7<2.5).and.(e6/14.7<2.5)
cut $8 energy(1)/beam>0.6.and.energy(2)/beam>0.6
cut $9 abs(d0(1)+d0(2))<0.4

zmax=10
*mess $call('zf.f(1,1)')
*mess $call('zfo.f(1,1)')
*cut $10 abs([fun](1,1))<10
*cut $20 abs([fun](2,1))<10
dphi=2.5
ve/cre tmin(9) r 1500 1000 1100 1150 975 950 1025 800 925
ve/cre tmax(9) r 1700 1200 1300 1350 1175 1150 1225 1000 1125
do i=1,9
  ci=10+[i]
	dphi1=40*([i]-0.5)+[dphi]
	dphi2=40*([i]-0.5)-5+[dphi]
  cut $[ci] abs([fun](1,2)-[dphi1])<17.5.and.abs([fun](1,2)-[dphi2])>3
  ci=20+[i]
  cut $[ci] abs([fun](2,2)-[dphi1])<17.5.and.abs([fun](2,2)-[dphi2])>3
	t1=tmin([i])
	t2=tmax([i])
	ci=30+[i] 
	cut $[ci] [t1]<tch([i])<[t2]
enddo
cut $30 $31.or.$32.or.$33.or.$34.or.$35.or.$36.or.$37.or.$38.or.$39
cut $89 l1+l2+l3+l4+l5+l6>0.and.7<(e1+e2+e3+e4+e5+e6)/(l1+l2+l3+l4+l5+l6)<17
return

macro data
n=0
n=[n]+1; ebeam[n]=508.7 ; il[n]=0
n=[n]+1; ebeam[n]=509.8 ; il[n]=0
n=[n]+1; ebeam[n]=510.5 ; il[n]=0
n=[n]+1; ebeam[n]=525   ; il[n]=363.2
n=[n]+1; ebeam[n]=537.5 ; il[n]=546.1
n=[n]+1; ebeam[n]=550   ; il[n]=485.7
n=[n]+1; ebeam[n]=562.5 ; il[n]=517.3
n=[n]+1; ebeam[n]=575   ; il[n]=419.1
n=[n]+1; ebeam[n]=587.5 ; il[n]=534.4
n=[n]+1; ebeam[n]=600   ; il[n]=489.8999940
n=[n]+1; ebeam[n]=612.5 ; il[n]=546.7000120
n=[n]+1; ebeam[n]=625   ; il[n]=440.5000000
n=[n]+1; ebeam[n]=637.5 ; il[n]=496.2000120
n=[n]+1; ebeam[n]=650   ; il[n]=456.1000060
n=[n]+1; ebeam[n]=662.5 ; il[n]=526.2999880
n=[n]+1; ebeam[n]=675   ; il[n]=559.7999880
n=[n]+1; ebeam[n]=687.5 ; il[n]=576.2000120
n=[n]+1; ebeam[n]=700   ; il[n]=577.0999760
n=[n]+1; ebeam[n]=712.5 ; il[n]=593.2000120
n=[n]+1; ebeam[n]=725   ; il[n]=432.2999880
n=[n]+1; ebeam[n]=737.5 ; il[n]=600.0999760
n=[n]+1; ebeam[n]=750   ; il[n]=694.7000120
n=[n]+1; ebeam[n]=762.5 ; il[n]=487.2999880
n=[n]+1; ebeam[n]=775   ; il[n]=546.5000000
n=[n]+1; ebeam[n]=787.5 ; il[n]=502.3999940
n=[n]+1; ebeam[n]=800   ; il[n]=439.3999940
n=[n]+1; ebeam[n]=812.5 ; il[n]=509.8999940
n=[n]+1; ebeam[n]=825   ; il[n]=475.6000060
n=[n]+1; ebeam[n]=850   ; il[n]=462.2000120
n=[n]+1; ebeam[n]=862.5 ; il[n]=504.6000060
n=[n]+1; ebeam[n]=875   ; il[n]=506.3999940
n=[n]+1; ebeam[n]=887.5 ; il[n]=525.7000120
n=[n]+1; ebeam[n]=900   ; il[n]=384.1000060
n=[n]+1; ebeam[n]=912.5 ; il[n]=481.3999940
n=[n]+1; ebeam[n]=925   ; il[n]=401.2999880
n=[n]+1; ebeam[n]=935   ; il[n]=637.5000000
n=[n]+1; ebeam[n]=945   ; il[n]=577.0999760
n=[n]+1; ebeam[n]=950   ; il[n]=451.2999880
n=[n]+1; ebeam[n]=962.5 ; il[n]=561.9000240
n=[n]+1; ebeam[n]=987.5 ; il[n]=467.5000000
n=[n]+1; ebeam[n]=1000  ; il[n]=538.0999760
n2=[n]

cuts=$1.and.$3.and.$4.and.$5.and.$6.and.$8.and.$9

dir=/work/users/konctbel/SepPar

ebeam=[ebeam3]
exec [dir]/sp#hists [ebeam]
nt/pl //exp/1.eton [cuts].and.eton<1.2
exec [dir]/sp#hists [ebeam]
ver=0
mod=1
cnt=1
ct=0
gl/imp ntup
do q=1,[ntup]
  exec [dir]/sp#hists [ebeam] n=[q] mod=[mod]
  hi/del 2
  ct=[ct]+1
  exec [dir]/rejectm //expi/1 test_[ct].hbook 2 [cuts]
  wait ! 5
  shell ls -ltr test_[ct].hbook > ls1_[cnt].txt
  cnt=[cnt]+1
  chain -tmp
  chain tmp test_[ct].hbook
  hi/del 1
  exec [dir]/rejectm //tmp/2 exp[q]_[ver]_[ebeam].hbook 1 [cuts]
  wait ! 5
  shell ls -ltr exp[q]_[ver]_[ebeam].hbook > ls1_[cnt].txt
  cnt=[cnt]+1
  shell rm test_[ct].hbook
enddo
return


macro mhad2011plx
exp=mhad2011
ebeam=510.5
gl/cre nxg 11
gl/cre nyg 11
gl/cre mapmd 1
ncnt=1
exec mapcal#prep parpl tmp [ncnt] [ebeam] [exp] no
ve/inp xi(1) $sigma(xi(1)+1000)
shell gfortran -shared -o accmappl.sl accmappl.f `cernlib mathlib packlib` -llapack -lm
ve/inp yr(1) yi(1)
ve/inp xr(1) -9.79
ve/cre yci(500) r
ve/cre zci(500) r
do i=1,500
  ve/inp yr(1) $sigma(([i]-1)*0.1)
  ve/inp yci([i]) yr(1)
  ve/inp zci([i]) $call('accmappl.sl(xr,yr,3)')
enddo
ve/cre dzc(500) r
set pmci 2
** exec $PER/s#vpl zci dzc yci dzc sz=0.01 ll=-1
set pmci 4
** exec $PER/s#vpl zc dzc yc dzc sz=0.01 ll=-1 o=s
return

macro mhad2011pl
exp=mhad2011
ebeam=510.5
gl/cre nxg 11
gl/cre nyg 11
gl/cre mapmd 1
file=[exp]_[ebeam]_mapcal.txt
if ($fexist([file]).ne.0) then
  shell rm [file]
endif
for/file  20 [file] new
close 20
ve/cre xr(1) r
ve/cre yr(1) r
do ncnt=1,9
  if ($fexist('map.tmp').ne.0) then
    shell rm map.tmp
  endif
  exec mapcal#prep parpl tmp [ncnt] [ebeam] [exp] no
  ve/inp xi(1) $sigma(xi(1)+1000)
  shell gfortran -shared -o accmappl.sl accmappl.f `cernlib mathlib packlib` -llapack -lm
  fmess 'Counter' [file]
  fmess [ncnt] [file]
  mess $call('accmappl.sl(xr,yr,3)')
  shell cat [file] map.tmp > tmp.tmp
  shell cp tmp.tmp [file]
*  fmess 'xi' [file]
*  fmess [nxg] [file]
*  do i=1,[nxg]
*    xit=xi([i])
*    xit=$format([xit],f10.2)
*    ve/inp xi([i]) [xit]
*    fmess $sigma(xi([i])) [file]
*  enddo
*  fmess 'yi' [file]
*  fmess [nyg] [file]
*  do i=1,[nyg]
*    yit=yi([i])
*    yit=$format([yit],f10.2)
*    ve/inp yi([i]) [yit]
*    fmess $sigma(yi([i])) [file]
*  enddo
*  fmess 'map' [file]
*  do j=1,[nyg]
*    ve/inp yr(1) yi([j])
*    do i=1,[nxg]
*      ve/inp xr(1) xi([i])
*      fmess $call('accmappl.sl(xr,yr,3)') [file]
*    enddo
*  enddo
enddo
return


macro mhad2011cal
exp=mhad2011
ebeam=510.5
gl/cre nxg 11
gl/cre nyg 11
gl/cre mapmd 1
do ncnt=5,8
*  exec mapcal#prep parn tmp [ncnt] [ebeam] [exp] no
enddo
ncnt=8
exec mapcal#prep parn tmp [ncnt] [ebeam] [exp] no
return



macro prep p=parn fout=mapfit.out.old ncnt=1 ebeam=510.2 exp=mhad2011 edata=no ver=0

chain -exp
chain exp mhad2011_510.5-2_col_p1.hbook
chain exp mhad2011_510.5-2_col_p2.hbook
chain exp mhad2011_510.5-2_col_p3.hbook
chain exp mhad2011_510.5-2_col_p4.hbook
chain exp -P /work/users/konctbel/exp/MHAD2011-2/col/

*chain -exp
*chain exp exp1_0_510.5.hbook
*chain exp exp2_0_510.5.hbook
*chain exp exp3_0_510.5.hbook
*chain exp exp4_0_510.5.hbook

mess shell ./mkcal.sh [exp] [ebeam] [ncnt]
shell ./mkcal.sh [exp] [ebeam] [ncnt]
dir=[exp]_[ebeam]/counter_[ncnt]
datafile=[exp]_[ebeam]_nc[ncnt].dat
datafile=mapfit.txt

gl/cre sc $sigma(180/3.1415927)
r1=10.666
r2=13.834
r=([r1]+[r2])/2

*-----------------------------

*ncnt=1
mnz=80
mzmin=-20
mzmax=20
mnf=80
mfmax=$sigma((40*[ncnt]-20)+40)
mfmin=$sigma((40*[ncnt]-20)-40)
dz=$sigma(([mzmax]-[mzmin])/[mnz])
df=$sigma(([mfmax]-[mfmin])/[mnf])
ve/cre xvi([mnz]) r
ve/cre yvi([mnf]) r
l0=[mzmin]+[dz]/2.
r0=[mzmax]-[dz]/2.
sigma xvi = array([mnz],[l0]#[r0])
l0=[mfmin]+[df]/2.
r0=[mfmax]-[df]/2.
sigma yvi = array([mnf],[l0]#[r0])

*goto 1
if ([edata].eq.'yes') then

cuts=$1.and.$3.and.$4.and.$5.and.$6.and.$8.and.$9.and.eventtime.eq.0

profile 1000 ! [mnz] [mzmin] [mzmax] ! !
profile 2000 ! [mnz] [mzmin] [mzmax] ! !
profile 3000 ! [mnz] [mzmin] [mzmax] ! !

ve/cre  avi([mnz]) r
ve/cre davi([mnz]) r
ve/cre  av([mnz],[mnf]) r
ve/cre dav([mnz],[mnf]) r

do i=1,[mnf]

  phi1=$sigma([mfmin]+[df]*([i]-0.5))
  phi2=[phi1]+360
  phi0=[phi1]-360

  nt/pl //exp/1.ach([ncnt])*sin(theta(1))%(z0(1)+[r]/tan(theta(1))) _
  [cuts].and.(abs(phi(1)*[sc]-([phi0]))<0.5*[df].or.abs(phi(1)*[sc]-([phi1]))<0.5*[df].or.abs(phi(1)*[sc]-([phi2]))<0.5*[df]) idh=1000
  
  nt/pl //exp/1.ach([ncnt])*sin(theta(2))%(z0(2)+[r]/tan(theta(2))) _
  [cuts].and.(abs(phi(2)*[sc]-([phi0]))<0.5*[df].or.abs(phi(2)*[sc]-([phi1]))<0.5*[df].or.abs(phi(2)*[sc]-([phi2]))<0.5*[df]) idh=2000
  
  hi/op/add 1000 2000 3000 1 1 e
	
  hi/get/cont 3000  avi
  hi/get/err  3000 davi
	
  ve/copy  avi(1:[mnz])  av(1:[mnz],[i])
  ve/copy davi(1:[mnz]) dav(1:[mnz],[i])

enddo
ve/del avi,davi
hi/del 1000
hi/del 2000
hi/del 3000

*ve/write av,dav [dir]/mapfit.txt 2g15.6
ve/write av,dav [dir]/[datafile] 2g15.6

endif
1:
*
ve/del av,dav
ve/cre  av([mnz],[mnf]) r
ve/cre dav([mnz],[mnf]) r
*ve/read av,dav [dir]/mapfit.txt 2g15.6
ve/read av,dav [dir]/[datafile] 2g15.6

nsc=[mnz]*[mnf]
ve/del av1,dav1
ve/read av1,dav1 [dir]/[datafile] 2g15.6

2d 1000 ! [mnz] [mzmin] [mzmax] [mnf] [mfmin] [mfmax]

hi/put/cont 1000  av
hi/put/err  1000 dav

*------------------

r1=10.666
r2=13.834
r=([r1]+[r2])/2
d1=0.5
d2=1
s0=0
mode=1
md1=1
md2=1
pk=3
pn=180
pm=40
zk=3
zn=100
zm=20
kx=3
ky=3
nx=11
ny=11
nx0=[nx]-2
ny0=[ny]-2
nxy0=[nx0]*[ny0]
npar0=21
ind0=[npar0]+[nx0]
nx0st=1
nx0sp=[nx0]
npar=[npar0]+[nx0]+[nx0]*[ny0]

*------------------

ve/cre s1f(9) r 3.9 43.5 82.6 122.6 163.4 203.9 245.2 285.4 324.7
ve/cre s2f(9) r 42.6 81.4 120.2 162.6 202.1 242.5 284.1 323.6 362.1
ve/cre ssf(9) r 17.6 57.6 96.3 136.3 178.0 217.7 259.1 299.5 338.8

ds0=2
*exec mapcal#s1s2 $sigma([ds0]+4.2+40*([ncnt]-1)) $sigma([ds0]+42.2+40*([ncnt]-1)) $sigma([ds0]+17.5+40*([ncnt]-1))
exec mapcal#s1s2 $sigma(s1f([ncnt])) $sigma(s2f([ncnt])) $sigma(ssf([ncnt]))
gl/imp s1g
gl/imp s2g
gl/imp ssg
gl/imp sigs1g
gl/imp sigs2g
gl/imp sigssg
exec mapcal#z1z2 $sigma(ssf([ncnt]))
gl/imp z1g
gl/imp z2g
*z1g=-10
*z2g=10

*------------------

z0=0; dz0=0.1
z1=[z1g]; dz1=0.1
z2=[z2g];  dz2=0.1
s1=[s1g];  ds1=[sigs1g]/10
s2=[s2g];  ds2=[sigs2g]/10
sigz=0.6; dsigz=0.001
sigs=([sigs1g]+[sigs2g])/2; dsigs=([sigs1g]+[sigs2g])/10
sigs0=[sigssg]; dsigs0=0.01
pds=0;    dpds=0.001
rPMT=0.86; drPMT=0.001
sPMT=[ssg]; dsPMT=0.01
apmt=14; daPMT=0.1
zs1=[z1g]; dzs1=0.1
zs2=[z2g]; dzs2=0.1
dzds=0; ddzds=0.00001
ss=[ssg]; dss=0.01
alpha=0; dalpha=0.001
ps1=0.05; dps1=0.01
ps2=0.05; dps2=0.01
ts1=40; dts1=1
ts2=40; dts2=1
ve/cre xi([nx]) r
ve/cre yi([ny]) r
l0=[z1]
r0=[z2]
sigma xi = array([nx],[l0]#[r0])
l0=[s1]
r0=[s2]
sigma yi = array([ny],[l0]#[r0])

ve/cre wls([nx0]) r [nx0]*15
ve/cre dwls([nx0]) r [nx0]*0.3
ve/cre map([nx0],[ny0]) r [nxy0]*2
ve/cre dmap([nx0],[ny0]) r [nxy0]*0.1
exec mapcal#mapc [ssg]
sigma  map =  mape/([r2]-[r1])
sigma dmap = dmape/([r2]-[r1])

ve/cre parn([npar]) r [z1] [z2] [s1] [s2] [sigz] [sigs] [pds] [rpmt] [apmt] [zs1] [zs2] [dzds] [ss] [z0] [sPMT] [sigs0] [ps1] [ps2] [ts1] [ts2]
i1=[ind0]-[nx0]+1
i2=[ind0]
ve/copy wls(1:[nx0]) parn([i1]:[i2])
do j=1,[ny0]
  i1=[i1]+[nx0]
  i2=[i2]+[nx0]
  ve/copy map(1:[nx0],[j]) parn([i1]:[i2])
enddo


if ([p].eq.'parpl') then

  n=0
  n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out
  n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out
  n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out
  n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out
  n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out.cut
  n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out.cut
  n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out.cut
  n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out.cut
  n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out
  fout=[file[ncnt]]
  
*  shell /work/snd2000/users/konctbel/SNDSimACC/readlifstat_new mapfit.out
  shell /work/snd2000/users/konctbel/SNDSimACC/readlifstat [fout]
  ve/del pars0,dpars0
  ve/read pars0,dpars0 [fout].pars

  fout=[dir]/mapfit.out
  dir=.
  ve/del xi,yi
  shell fgrep -e "map0,xi:" [fout] >& tmp.txt
  shell $unquote('cat tmp.txt | sed "s/map0,xi:/ /g" > tmp1.txt')
  ve/read xi tmp1.txt
  shell fgrep -e "map0,yi:" [fout] >& tmp.txt
  shell $unquote('cat tmp.txt | sed "s/map0,yi:/ /g" > tmp1.txt')
  ve/read yi tmp1.txt

  ve/inp xi(1) $sigma(xi(1) - 1000)
  gl/cre nxg $vlen(xi)
  gl/cre nyg $vlen(yi)
  nx=[nxg]
  ny=[nyg]
  nx0=[nx]-2
  ny0=[ny]-2
  ind0=[npar0]+[nx0]
  nx0st=1
  nx0sp=[nx0]
  npar=[npar0]+[nx0]+[nx0]*[ny0]
  ve/cre wls([nx0]) r [nx0]*5
  ve/cre dwls([nx0]) r [nx0]*0.3
  nxy0=[nx0]*[ny0]
  ve/cre map([nx0],[ny0]) r [nxy0]*2
  ve/cre dmap([nx0],[ny0]) r [nxy0]*0.1

  ve/cre parn([npar]) r

  i1=2
  i2=[npar]+1
  ve/copy pars0([i1]:[i2]) parn(1:[npar])
  i=1
  i=[i]+1; z1=pars0([i]);    dz1=dpars0([i])
  i=[i]+1; z2=pars0([i]);    dz2=dpars0([i])
  i=[i]+1; s1=pars0([i]);    ds1=dpars0([i])
  i=[i]+1; s2=pars0([i]);    ds2=dpars0([i])
  i=[i]+1; sigz=pars0([i]);  dsigz=dpars0([i])
  i=[i]+1; sigs=pars0([i]);  dsigs=dpars0([i])
  i=[i]+1; pds=pars0([i]);   dpds=dpars0([i])
  i=[i]+1; rPMT=pars0([i]);  drPMT=dpars0([i])
  i=[i]+1; apmt=pars0([i]);  daPMT=dpars0([i])
  i=[i]+1; zs1=pars0([i]);   dzs1=dpars0([i])
  i=[i]+1; zs2=pars0([i]);   dzs2=dpars0([i])
  i=[i]+1; dzds=pars0([i]);  ddzds=dpars0([i])
  i=[i]+1; ss=pars0([i]);    dss=dpars0([i])
  i=[i]+1; alpha=pars0([i]); dalpha=dpars0([i])
  i=[i]+1; z0=pars0([i]);    dz0=dpars0([i])
  i=[i]+1; sPMT=pars0([i]);    dsPMT=dpars0([i])
  i=[i]+1; sigs0=pars0([i]);  dsigs0=dpars0([i])
  i=[i]+1; ps1=pars0([i]);  dps1=dpars0([i])
  i=[i]+1; ps2=pars0([i]);  dps2=dpars0([i])
  i=[i]+1; ts1=pars0([i]);  dts1=dpars0([i])
  i=[i]+1; ts2=pars0([i]);  dts2=dpars0([i])
  i1=[i]+1;
  i2=[i1]+[nx0]-1
  ve/copy  pars0([i1]:[i2])  wls(1:[nx0])
*  ve/copy dpars0([i1]:[i2]) dwls(1:[nx0])
  do i=1,[ny0]
    i1=[i2]+1
    i2=[i1]+[nx0]-1
    ve/copy  pars0([i1]:[i2])  map(1:[nx0],[i])
*    ve/copy dpars0([i1]:[i2]) dmap(1:[nx0],[i])
  enddo
endif

if ([p].eq.'pars0') then
*  shell /work/snd2000/users/konctbel/SNDSimACC/readlifstat_new mapfit.out
  shell /work/snd2000/users/konctbel/SNDSimACC/readlifstat [fout]
  ve/del pars0,dpars0
  ve/read pars0,dpars0 [fout].pars

  xx=4
  if ([xx].eq.4) then
  gl/imp nxg
  gl/imp nyg
  nx=[nxg]
  ny=[nyg]
  nx0=[nx]-2
  ny0=[ny]-2
  ind0=[npar0]+[nx0]
  nx0st=1
  nx0sp=[nx0]
  npar=[npar0]+[nx0]+[nx0]*[ny0]
  ve/cre wls([nx0]) r [nx0]*5
  ve/cre dwls([nx0]) r [nx0]*0.3
  nxy0=[nx0]*[ny0]
  ve/cre map([nx0],[ny0]) r [nxy0]*2
  ve/cre dmap([nx0],[ny0]) r [nxy0]*0.1
  ve/cre xi([nx]) r
  ve/cre yi([ny]) r
  sigma xi = array([nx],-12#12)
  l0=5+40*([ncnt]-1)
  r0=43+40*([ncnt]-1)
  sigma yi = array([ny],[l0]#[r0])
  ve/cre parn([npar]) r
  endif

  i1=2
  i2=[npar]+1
  ve/copy pars0([i1]:[i2]) parn(1:[npar])
  i=1
  i=[i]+1; z1=pars0([i]);    dz1=dpars0([i])
  i=[i]+1; z2=pars0([i]);    dz2=dpars0([i])
  i=[i]+1; s1=pars0([i]);    ds1=dpars0([i])
  i=[i]+1; s2=pars0([i]);    ds2=dpars0([i])
  i=[i]+1; sigz=pars0([i]);  dsigz=dpars0([i])
  i=[i]+1; sigs=pars0([i]);  dsigs=dpars0([i])
  i=[i]+1; pds=pars0([i]);   dpds=dpars0([i])
  i=[i]+1; rPMT=pars0([i]);  drPMT=dpars0([i])
  i=[i]+1; apmt=pars0([i]);  daPMT=dpars0([i])
  i=[i]+1; zs1=pars0([i]);   dzs1=dpars0([i])
  i=[i]+1; zs2=pars0([i]);   dzs2=dpars0([i])
  i=[i]+1; dzds=pars0([i]);  ddzds=dpars0([i])
  i=[i]+1; ss=pars0([i]);    dss=dpars0([i])
  i=[i]+1; alpha=pars0([i]); dalpha=dpars0([i])
  i=[i]+1; z0=pars0([i]);    dz0=dpars0([i])
  i=[i]+1; sPMT=pars0([i]);    dsPMT=dpars0([i])
  i=[i]+1; sigs0=pars0([i]);  dsigs0=dpars0([i])
  i=[i]+1; ps1=pars0([i]);  dps1=dpars0([i])
  i=[i]+1; ps2=pars0([i]);  dps2=dpars0([i])
  i=[i]+1; ts1=pars0([i]);  dts1=dpars0([i])
  i=[i]+1; ts2=pars0([i]);  dts2=dpars0([i])
  i1=[i]+1;
  i2=[i1]+[nx0]-1
  ve/copy  pars0([i1]:[i2])  wls(1:[nx0])
*  ve/copy dpars0([i1]:[i2]) dwls(1:[nx0])
  do i=1,[ny0]
    i1=[i2]+1
    i2=[i1]+[nx0]-1
    ve/copy  pars0([i1]:[i2])  map(1:[nx0],[i])
*    ve/copy dpars0([i1]:[i2]) dmap(1:[nx0],[i])
  enddo
*  exec mapcal#mapc [ssg]
*  sigma  map =  mape/([r2]-[r1])
*  sigma dmap = dmape/([r2]-[r1])
*  s1=[s1g];  ds1=[sigs1g]/10
*  s2=[s2g];  ds2=[sigs2g]/10
*  sigs=([sigs1g]+[sigs2g])/2; dsigs=([sigs1g]+[sigs2g])/10
*  alpha=0; dalpha=0.0001
*  dzds=0; ddzds=0.0001
*  z1=[z1g];
*  z2=[z2g];
*
  i2=[ind0]
  do i=1,-[ny0]
    i1=[i2]+1
    i2=[i1]+[nx0]-1
*    ve/copy  map(1:[nx0],[i]) parn([i1]:[i2])
  enddo
endif

if ([p].eq.'chmap') then
*  shell /work/snd2000/users/konctbel/SNDSimACC/readlifstat_new mapfit.out
  shell /work/snd2000/users/konctbel/SNDSimACC/readlifstat [fout]
  ve/del pars0,dpars0
  ve/read pars0,dpars0 [fout].pars

  i1=2
  i2=[npar]+1
  ve/copy pars0([i1]:[i2]) parn(1:[npar])
  i=1
  i=[i]+1; z1=pars0([i]);    dz1=dpars0([i])
  i=[i]+1; z2=pars0([i]);    dz2=dpars0([i])
  i=[i]+1; s1=pars0([i]);    ds1=dpars0([i])
  i=[i]+1; s2=pars0([i]);    ds2=dpars0([i])
  i=[i]+1; sigz=pars0([i]);  dsigz=dpars0([i])
  i=[i]+1; sigs=pars0([i]);  dsigs=dpars0([i])
  i=[i]+1; pds=pars0([i]);   dpds=dpars0([i])
  i=[i]+1; rPMT=pars0([i]);  drPMT=dpars0([i])
  i=[i]+1; apmt=pars0([i]);  daPMT=dpars0([i])
  i=[i]+1; zs1=pars0([i]);   dzs1=dpars0([i])
  i=[i]+1; zs2=pars0([i]);   dzs2=dpars0([i])
  i=[i]+1; dzds=pars0([i]);  ddzds=dpars0([i])
  i=[i]+1; ss=pars0([i]);    dss=dpars0([i])
  i=[i]+1; alpha=pars0([i]); dalpha=dpars0([i])
  i=[i]+1; z0=pars0([i]);    dz0=dpars0([i]); dz0=0.1
  i=[i]+1; sPMT=pars0([i]);    dsPMT=dpars0([i])
  i=[i]+1; sigs0=pars0([i]);  dsigs0=dpars0([i])
  i=[i]+1; ps1=pars0([i]);  dps1=dpars0([i])
  i=[i]+1; ps2=pars0([i]);  dps2=dpars0([i])
  i=[i]+1; ts1=pars0([i]);  dts1=dpars0([i])
  i=[i]+1; ts2=pars0([i]);  dts2=dpars0([i])
  i1=[i]+1;
  i2=[i1]+[nx0]-1
  ve/copy  pars0([i1]:[i2])  wls(1:[nx0])
  ve/copy dpars0([i1]:[i2]) dwls(1:[nx0])
  do i=1,[ny0]
    i1=[i2]+1
    i2=[i1]+[nx0]-1
    ve/copy  pars0([i1]:[i2])  map(1:[nx0],[i])
    ve/copy dpars0([i1]:[i2]) dmap(1:[nx0],[i])
  enddo
  gl/imp nxg
  gl/imp nyg
  nx=[nxg]
  ny=[nyg]
  nx0=[nx]-2
  ny0=[ny]-2
  ind0=[npar0]+[nx0]
  nx0st=1
  nx0sp=[nx0]
  npar=[npar0]+[nx0]+[nx0]*[ny0]
  ve/cre wls([nx0]) r [nx0]*5
  ve/cre dwls([nx0]) r [nx0]*0.3
  nxy0=[nx0]*[ny0]
  ve/cre map([nx0],[ny0]) r [nxy0]*2
  ve/cre dmap([nx0],[ny0]) r [nxy0]*0.1
  ve/cre xi([nx]) r
  ve/cre yi([ny]) r
  sigma xi = array([nx],-12#12)
  l0=5+40*([ncnt]-1)
  r0=43+40*([ncnt]-1)
  sigma yi = array([ny],[l0]#[r0])
  ve/cre xr(1) r
  ve/cre yr(1) r
  do i=1,[nx0]
    ve/inp xr(1) $sigma(xi([i]))
    ve/inp wls([i]) $call('accmappl.sl(xr,yr,2)')
    do j=1,[ny0]
      ve/inp yr(1) $sigma(yi([j]))
      ve/inp map([i],[j]) $call('accmappl.sl(xr,yr,1)')
    enddo
  enddo
  ve/cre parn([npar]) r [z1] [z2] [s1] [s2] [sigz] [sigs] [pds] [rpmt] [apmt] [zs1] [zs2] [dzds] [ss] [z0] [sPMT] [sigs0]
  i1=[ind0]-[nx0]+1
  i2=[ind0]
  ve/copy wls(1:[nx0]) parn([i1]:[i2])
  do j=1,[ny0]
    i1=[i1]+[nx0]
    i2=[i2]+[nx0]
    ve/copy map(1:[nx0],[j]) parn([i1]:[i2])
  enddo
endif

fname=mapfit

file=[dir]/[fname].inc
if ($fexist([file]).ne.0) then
*  shell rm [file]
  shell mv -v [file] [file].old
endif
for/file  20 [file] new
close 20

fmess '      double precision r1, r2, r, s1, s2, ss, z1, z2, d1, d2, z0, s0,' [file]
fmess '     &  sigz, sigs, mode, pds, dss, aPMT, zPMT, rPMT, alpha, sPMT,' [file]
fmess '     &  sigs0, rs1, rs2, zs1, zs2, dsdz, ts1, ts2, ps1, ps2' [file]
fmess '      common /accmap/ r1, r2, r, s1, s2, ss, z1, z2, d1, d2, z0, s0,' [file]
fmess '     &  sigz, sigs, mode, pds, dss, aPMT, zPMT, rPMT, alpha, sPMT,' [file]
fmess '     &  sigs0, rs1, rs2, zs1, zs2, dsdz, ts1, ts2, ps1, ps2' [file]
*fmess '      common/mapfitc/md1,md2,xi(20),yi(20),zi(20,20),' [file]
txt='      common/mapfitc/md1,md2,'
txt=$unquote([txt])xi([nx]),yi([ny]),zi([nx],[ny]),
fmess [txt] [file]
fmess '     &  nx0start,nx0stop,ind0,kx,ky,nx,ny,' [file]
fmess '     &  nx0,ny0,kag,kwls,kpds,kpmt' [file]
fmess '      integer kx, ky, nx, ny, nx0, ny0, nx0start, nx0stop, ind0' [file]
fmess '      double precision xi, yi, zi, md1, md2' [file]
*fmess '      common/wlsfit/yis(20)' [file]
txt='      common/wlsfit/'
txt=$unquote([txt])yis([nx])
fmess [txt] [file]
fmess '      double precision yis' [file]
*fmess '      common/ptfit/pxi(180),pyi(180),pk,pn,pm' [file]
txt='      common/ptfit/'
txt=$unquote([txt])pxi([pn]),pyi([pn]),pk,pn,pm
fmess [txt] [file]
fmess '      integer pk,pn,pm' [file]
fmess '      double precision pxi,pyi' [file]
*fmess '      common/pzfit/zxi(100),zyi(100),zk,zn,zm' [file]
txt='      common/pzfit/'
txt=$unquote([txt])zxi([zn]),zyi([zn]),zk,zn,zm
fmess [txt] [file]
fmess '      integer zk,zn,zm' [file]
fmess '      double precision zxi,zyi' [file]
*fmess '      data r1,r2/10.666,13.834/' [file]
txt='      data r1,r2'
txt=$unquote([txt])/[r1],[r2]/
fmess [txt] [file]
*fmess '      data d1,d2,s0,mode,md1,md2/0.5,1,0,1,1,1/' [file]
txt='      data d1,d2,s0,mode,md1,md2'
txt=$unquote([txt])/[d1],[d2],[s0],[mode],[md1],[md2]/
fmess [txt] [file]
*fmess '      data pk,pn,pm/3,180,40/' [file]
txt='      data pk,pn,pm'
txt=$unquote([txt])/[pk],[pn],[pm]/
fmess [txt] [file]
*fmess '      data zk,zn,zm/3,100,20/' [file]
txt='      data zk,zn,zm'
txt=$unquote([txt])/[zk],[zn],[zm]/
fmess [txt] [file]
fmess '      integer kag,kwls,kpds,kpmt,ktail' [file]
fmess '      data kag,kwls,kpds,kpmt,ktail/1,1,1,1,1/' [file]
*fmess '      data kx,ky,nx,ny,nx0,ny0,ind0/3,3,7,7,5,5,18/' [file]
txt='      data kx,ky,nx,ny,nx0,ny0,ind0'
txt=$unquote([txt])/[kx],[ky],[nx],[ny],[nx0],[ny0],[ind0]/
fmess [txt] [file]
*fmess '      data nx0start,nx0stop/1,5/' [file]
txt='      data nx0start,nx0stop'
txt=$unquote([txt])/[nx0st],[nx0sp]/
fmess [txt] [file]
fmess '      integer npar0' [file]
txt='      data npar0 '
txt=$unquote([txt])/[npar]/
fmess [txt] [file]
*fmess '      common/parnc/parn(34)' [file]
txt='      common/parnc/'
txt=$unquote([txt])parn([npar])
fmess [txt] [file]
fmess '      double precision parn' [file]
fmess '      data parn /' [file]
exec mapcal#vewrite [file] parn f10.4 [npar] 4
fmess '     &  /' [file]
fmess '' [file]
fmess '' [file]
fmess '' [file]
fmess '' [file]

fmess '      data xi /' [file]
exec mapcal#vewrite [file] xi f10.4 [nx] 4
fmess '     &  /' [file]

fmess '      data yi /' [file]
exec mapcal#vewrite [file] yi f10.4 [ny] 4
fmess '     &  /' [file]

*---------------

file=[dir]/[fname]_data.inc
if ($fexist([file]).eq.1) then
*  shell rm -v [file]
  shell mv -v [file] [file].old
endif
for/file  20 [file] new
close 20

fmess '      integer nvx, nvy' [file]
txt='      data nvx,nvy'
txt=$unquote([txt])/[mnz],[mnf]/
fmess [txt] [file]

txt='      common/mapc/'
txt=$unquote([txt])av([mnz],[mnf]),dav([mnz],[mnf]),xvi([mnz]),yvi([mnf])
fmess [txt] [file]
fmess '      double precision av,dav,xvi,yvi' [file]

d=$sigma(([mzmax]-[mzmin])/[mnz]/2)
ve/cre xvi([mnz]) r
l=[mzmin]+[d]
r=[mzmax]-[d]
sigma xvi = array([mnz],[l]#[r])
fmess '      data xvi /' [file]
exec mapcal#vewrite [file] xvi f10.2 [mnz] 4
fmess '     &  /' [file]

d=$sigma(([mfmax]-[mfmin])/[mnf]/2)
ve/cre yvi([mnf]) r
l=[mfmin]+[d]
r=[mfmax]-[d]
sigma yvi = array([mnf],[l]#[r])
fmess '      data yvi /' [file]
exec mapcal#vewrite [file] yvi f10.2 [mnf] 4
fmess '     &  /' [file]

fmess '      data av /' [file]
exec mapcal#vewrite [file] av1 g10.3 [nsc] 4
fmess '     &  /' [file]
fmess '      data dav /' [file]
exec mapcal#vewrite [file] dav1 g10.3 [nsc] 4
fmess '     &  /' [file]

*---------------

file=[dir]/[fname].dat
if ($fexist([file]).eq.1) then
  shell rm -v [file]
endif
for/file  20 [file] new
close 20

fmess 'SET TITLE' [file]
ni=0
fmess 'counter ununiformity fitting' [file]
fmess 'PARAMETERS' [file]
ls=[z1]-1
rs=[z1]+1
ni=[ni]+1; pname=z1; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5); tmp4=$format([ls],e13.5); tmp5=$format([rs],e13.5)
tmp=$unquote([tmp1])$quote([pname])      $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
fmess [tmp] [file]
ls=[z2]-1
rs=[z2]+1
ni=[ni]+1; pname=z2; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5); tmp4=$format([ls],e13.5); tmp5=$format([rs],e13.5)
tmp=$unquote([tmp1])$quote([pname])      $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
fmess [tmp] [file]
ls=[s1]-1
rs=[s1]+1
ni=[ni]+1; pname=s1; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5); tmp4=$format([ls],e13.5); tmp5=$format([rs],e13.5)
tmp=$unquote([tmp1])$quote([pname])      $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5]) 
fmess [tmp] [file]
ls=[s2]-1
rs=[s2]+1
ni=[ni]+1; pname=s2; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5); tmp4=$format([ls],e13.5); tmp5=$format([rs],e13.5)
tmp=$unquote([tmp1])$quote([pname])      $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
fmess [tmp] [file]
ni=[ni]+1; pname=sigz; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5)
tmp=$unquote([tmp1])$quote([pname])    $unquote([tmp2])$unquote([tmp3]) 0 1
fmess [tmp] [file]
ni=[ni]+1; pname=sigs; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5)
tmp=$unquote([tmp1])$quote([pname])    $unquote([tmp2])$unquote([tmp3]) 0 1
fmess [tmp] [file]
ni=[ni]+1; pname=pds; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5)
tmp=$unquote([tmp1])$quote([pname])     $unquote([tmp2])$unquote([tmp3]) -0.01 0.01
fmess [tmp] [file]
ni=[ni]+1; pname=rPMT; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5)
tmp=$unquote([tmp1])$quote([pname])    $unquote([tmp2])$unquote([tmp3]) 0 1
fmess [tmp] [file]
ni=[ni]+1; pname=aPMT; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5)
tmp=$unquote([tmp1])$quote([pname])    $unquote([tmp2])$unquote([tmp3]) 0 30
fmess [tmp] [file]
ls=[z1]-2
rs=[z1]+1
ni=[ni]+1; pname=zs1; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5); tmp4=$format([ls],e13.5); tmp5=$format([rs],e13.5)
tmp=$unquote([tmp1])$quote([pname])     $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
fmess [tmp] [file]
ls=[z2]-1
rs=[z2]+2
ni=[ni]+1; pname=zs2; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5); tmp4=$format([ls],e13.5); tmp5=$format([rs],e13.5)
tmp=$unquote([tmp1])$quote([pname])     $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
fmess [tmp] [file]
ni=[ni]+1; pname=dzds; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5)
tmp=$unquote([tmp1])$quote([pname])    $unquote([tmp2])$unquote([tmp3]) -0.1 0.1
fmess [tmp] [file]
ls=[ss]-1
rs=[ss]+1
ni=[ni]+1; pname=ss; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5); tmp4=$format([ls],e13.5); tmp5=$format([rs],e13.5)
tmp=$unquote([tmp1])$quote([pname])      $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
fmess [tmp] [file]
ni=[ni]+1; pname=alpha; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5)
tmp=$unquote([tmp1])$quote([pname])   $unquote([tmp2])$unquote([tmp3]) -20 20
fmess [tmp] [file]
ni=[ni]+1; pname=z0; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5)
tmp=$unquote([tmp1])$quote([pname])      $unquote([tmp2])$unquote([tmp3]) -1 1
fmess [tmp] [file]
ls=[sPMT]-1
rs=[sPMT]+1
ni=[ni]+1; pname=sPMT; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5); tmp4=$format([ls],e13.5); tmp5=$format([rs],e13.5)
tmp=$unquote([tmp1])$quote([pname])    $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
fmess [tmp] [file]
ni=[ni]+1; pname=sigs0; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5)
tmp=$unquote([tmp1])$quote([pname])   $unquote([tmp2])$unquote([tmp3]) 0 1
fmess [tmp] [file]
ni=[ni]+1; pname=ps1; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5)
tmp=$unquote([tmp1])$quote([pname])     $unquote([tmp2])$unquote([tmp3]) 0 10
fmess [tmp] [file]
ni=[ni]+1; pname=ps2; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5)
tmp=$unquote([tmp1])$quote([pname])     $unquote([tmp2])$unquote([tmp3]) 0 10
fmess [tmp] [file]
ni=[ni]+1; pname=ts1; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5)
tmp=$unquote([tmp1])$quote([pname])     $unquote([tmp2])$unquote([tmp3]) 5 1000
fmess [tmp] [file]
ni=[ni]+1; pname=ts2; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5)
tmp=$unquote([tmp1])$quote([pname])     $unquote([tmp2])$unquote([tmp3]) 5 1000
fmess [tmp] [file]
nif=[ni]
do i=1,[nx0]
  tmp0=wls_[i]
  wlsi=$sigma(wls([i]))
  wlsil=0.8*[wlsi]
  wlsiu=1.2*[wlsi]
  dwlsi=$sigma(dwls([i]))
  ni=[ni]+1; pname=wls_[i]; tmp1=$format([ni],i-4); tmp2=$format([wlsi],e13.5); tmp3=$format([dwlsi],e13.5); tmp4=tmp2=$format([wlsil],e13.5); tmp5=tmp2=$format([wlsiu],e13.5);
  tmp=$unquote([tmp1])$quote([pname])   $unquote([tmp2])$unquote([tmp3]) 0 30
  fmess [tmp] [file]
enddo
nig=[ni]
do j=1,[ny0]
  do i=1,[nx0]
    tmp0=n_x[i]y[j]
    mapij=map([i],[j])
    mapijl=0.5*[mapij]
    mapiju=1.5*[mapij]
    dmapij=dmap([i],[j])
    ni=[ni]+1; pname=n_x[i]y[j]; tmp1=$format([ni],i-4); tmp2=$format([mapij],e13.5); tmp3=$format([dmapij],e13.5); tmp4=$format([mapijl],e13.5); tmp5=$format([mapiju],e13.5);
    tmp=$unquote([tmp1])$quote([pname])  $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
*    tmp=$unquote([tmp1])$quote([pname])  $unquote([tmp2])$unquote([tmp3]) 0 10
    fmess [tmp] [file]
  enddo
enddo
fmess '' [file]
fmess '' [file]
*fmess 'set err 0.5' [file]
*fmess 'set eps 1e-9' [file]
*fmess 'set strategy 2' [file]
*fmess '' [file]

dnx0=3

do i=1,[npar0]
  txt=fix [i]
  fmess [txt] [file]
enddo
do i=$sigma([npar0]+1),[npar]
  txt=release [i]
  fmess [txt] [file]
enddo
fmess '' [file]
fmess 'mini' [file]
fmess '' [file]

ve/cre fix0([npar0]) r [npar0]*1
ve/cre fix0([npar0]) r 17*1 4*0
ve/cre fix1([npar0]) r 1 0 1 1 1 1 1 0 0 1 0 1 1 1 0 0 1
ve/cre fix2([npar0]) r 0 0 1 1 0 1 1 0 0 0 0 1 1 1 0 0 1
ve/cre fix3([npar0]) r 0 1 1 1 1 1 1 1 1 0 1 1 1 1 0 1 1


im=$sigma([nx0]-[dnx0])
do i=0,[im]
  if ([i].eq.0) then
    fixv=fix0
  endif
  if (([i].gt.0).and.([i].lt.[im])) then
    fixv=fix0
  endif
  if ([i].eq.[im]) then
    fixv=fix0
  endif
  do j=1,[npar0]
    if ($sigma([fixv]([j])).eq.1) then
      txt=release [j]
    else
      txt=fix [j]
    endif
    fmess [txt] [file]
  enddo
  i1=[i]+1
  i2=[i1]+[dnx0]-1
  ind=[npar0]
  do j=0,[ny0]
    do k=1,[nx0]
      ind=[ind]+1
      if (([k].ge.[i1]).and.([k].le.[i2])) then
        txt=release [ind]
      else
        txt=fix [ind]
      endif
      fmess [txt] [file]
    enddo
  enddo
  fmess '' [file]
  fmess 'mini' [file]
  fmess '' [file]
enddo

do i=1,[npar]
  if (([i].ge.18).and.([i].le.21)) then
    txt=release [i]
  else
    txt=fix [i]
  endif
  fmess [txt] [file]
enddo
fmess '' [file]
fmess 'mini' [file]
fmess '' [file]

do i=1,[npar]
  txt=release [i]
  fmess [txt] [file]
enddo
fmess '' [file]
fmess 'mini' [file]
fmess 'mini' [file]
fmess 'mini' [file]
fmess '' [file]

fmess 'ret' [file]

*---------------------------------
if ([p].ne.'parpl') then
  shell ./mk.sh [dir]
endif

return

macro fmtprep fmtn=f8.1 nsc=10 nn=3 lch=' '
n0=$sigma(int(([nsc])/[nn]))
dn0=$sigma(([nsc])-[nn]*[n0])
if ([dn0].eq.0) then
  n0=[n0]-1
  dn0=[nn]
endif
if ([dn0].ne.1) then
  dn0=[dn0]-1
  if ([n0].ne.0) then
    gl/cre fmt [n0](5x,"&",2x,[nn]([fmtn],",",1x,)/),5x,"&",2x,[dn0]([fmtn],",",1x,)[fmtn]$unquote([lch])
    gl/cre fmt0 [n0](5x,"&",2x,[nn]([fmtn],",",1x,)/)
    gl/cre fmt1 5x,"&",2x,[dn0]([fmtn],",",1x,)[fmtn]$unquote([lch])
    gl/cre nn0 $sigma([n0]*[nn])
    gl/cre nn1 $sigma([dn0])
  else
    gl/cre fmt 5x,"&",2x,[dn0]([fmtn],",",1x,)[fmtn]$unquote([lch])
  endif
else
  if ([n0].ne.0) then
    gl/cre fmt [n0](5x,"&",2x,[nn]([fmtn],",",1x,)/),5x,"&",2x,[fmtn]$unquote([lch])
    gl/cre fmt0 [n0](5x,"&",2x,[nn]([fmtn],",",1x,)/)
    gl/cre fmt1 5x,"&",2x,[fmtn]$unquote([lch])
    gl/cre nn0 $sigma([n0]*[nn])
    gl/cre nn1 $sigma([dn0])
  else
    gl/cre fmt 5x,"&",2x,[fmtn]$unquote([lch])	  
  endif
endif	
return

macro vewrite file=vewrite.txt vname=v fmtn=f8.1 nsc=10 nn=3 lch=' '
if ($vexist([vname]).ne.0) then
  if ($fexist(tmp.f).ne.0) then
    shell rm tmp.f
  endif
  if ($fexist(tmp0.dat).ne.0) then
    shell rm tmp0.dat
  endif
  if ($fexist(tmp1.dat).ne.0) then
    shell rm tmp1.dat
  endif
*  for/file  20 tmp.f new 
  exec mapcal#fmtprep fmtn=[fmtn] nsc=[nsc] nn=[nn] lch=[lch]
  gl/imp fmt
  gl/imp fmt0
  gl/imp fmt1
  gl/imp nn0
  gl/imp nn1
  ve/del v0,v1
  ve/copy [vname](1:[nn0]) v0
  i1=[nn0]+1
  ve/copy [vname]([i1]:[nsc]) v1
*  ve/write [vname] tmp.f [fmt] ' '
  ve/write v0 tmp0.dat [fmt0]
  ve/write v1 tmp1.dat [fmt1]
*  close 20
*  shell cat [file] tmp.f > tmp2.f
  shell cat [file] tmp0.dat tmp1.dat > tmp2.f
  shell cat tmp2.f > [file]
endif
return


macro mapplx x0=0.
n=$vlen(xi)
xmin=xi(1)
xmax=xi([n])
ix0=$sigma(int(0.5+([n]-1)*([x0]-[xmin])/([xmax]-[xmin])))
mess [ix0]
n=$vlen(yi)
ny0=[n]-2
ve/cre xm([ny0]) r
ve/cre ym([ny0]) r
ve/cre dm([ny0]) r
i1=2
i2=[ny0]+1
ve/copy yi([i1]:[i2]) xm(1:[ny0])
ve/copy map([ix0],1:[ny0]) ym(1:[ny0])
r1=10.666
r2=13.834
sigma ym = ym*([r2]-[r1])
*
n=$vlen(xvi)
xmin=xvi(1)
xmax=xvi([n])
ix=$sigma(int(1.5+([n]-1)*([x0]-[xmin])/([xmax]-[xmin])))
ny=$vlen(yvi)
ve/cre dyvi([ny]) r
ve/del ye,dye
ve/copy  av([ix],1:[ny])  ye
ve/copy dav([ix],1:[ny]) dye
n=200
ymin=$sigma(vmin(yvi))
ymax=$sigma(vmax(yvi))
ve/cre y([n]) r
ve/cre ys([n]) r
ve/cre ya([n]) r
sigma x = array([n],[ymin]#[ymax])
ve/cre xr(1) r
ve/cre yr(1) r
do i=1,[n]
  ve/inp yr(1) $sigma(x([i]))
  ve/inp xr(1) $sigma(xvi([ix]))
  ve/inp y([i]) $call('accmappl_blocks.sl(xr,yr,0)')
*  ve/inp y([i]) $call('accmappl.sl(xr,yr,0)')
*  ve/inp ys([i]) $call('accmappl.sl(xr,yr,2)')
*  ve/inp ya([i]) $call('accmappl.sl(xr,yr,1)')
enddo
*sigma y = y*([r2]-[r1])
u=$sigma(1.1*max(vmax(y),vmax(ye+dye)))
null [ymin] [ymax] 0 [u]
ve/cre d([n]) r
* exec $PER/s#vpl y d x d sz=0.01 o=s
graph [n] x y sl
* exec $PER/s#vpl ys d x d sz=0.01 o=s
set plci 2
set ltyp 14
set basl 0.005
graph [n] x ys sl
* exec $PER/s#vpl ya d x d sz=0.01 o=s
set ltyp 13
set basl 0.01
set plci 1
graph [n] x ya sl
set ltyp 1
set plci 1
* exec $PER/s#vpl ye dye yvi dyvi o=s iatt=20 sz=0.05
set pmci 2
** exec $PER/s#vpl ym dm xm dm o=s
set pmci 1
return

macro mapply y0=17.
n=$vlen(yvi)
ymin=yvi(1)
ymax=yvi([n])
iy=$sigma(int(1.5+([n]-1)*([y0]-[ymin])/([ymax]-[ymin])))
nx=$vlen(xvi)
ve/cre dxvi([nx]) r
ve/del ye,dye
ve/copy  av(1:[nx],[iy])  ye
ve/copy dav(1:[nx],[iy]) dye
n=300
xmin=$sigma(vmin(xvi))
xmax=$sigma(vmax(xvi))
ve/cre y([n]) r
sigma x = array([n],[xmin]#[xmax])
ve/cre xr(1) r
ve/cre yr(1) r
do i=1,[n]
  ve/inp xr(1) $sigma(x([i]))
  ve/inp yr(1) $sigma(yvi([iy]))
*  ve/inp y([i]) $call('accmappl.sl(xr,yr,0)')
  ve/inp y([i]) $call('accmappl_blocks.sl(xr,yr,0)')
enddo
u=$sigma(1.1*max(vmax(y),vmax(ye+dye)))
null [xmin] [xmax] 0 [u]
ve/cre d([n]) r
set pmci 2
* exec $PER/s#vpl y d x d sz=0.01 o=s
set plci 2
graph [n] x y sl
set pmci 1
* exec $PER/s#vpl ye dye xvi dxvi o=s iatt=20 sz=0.05
return

macro mapplxy
nx=80
xmin=-20
xmax=20
ny=80
ymin=70
ymax=130
goto 1
l=$sigma([xmin]+([xmax]-[xmin])/[nx]/2)
r=$sigma([xmax]-([xmax]-[xmin])/[nx]/2)
sigma x = array([nx],[l]#[r])
l=$sigma([ymin]+([ymax]-[ymin])/[ny]/2)
r=$sigma([ymax]-([ymax]-[ymin])/[ny]/2)
sigma y = array([ny],[l]#[r])
ve/cre  mapc([nx],[ny]) r
ve/cre dmapc([nx],[ny]) r
ve/cre xr(1) r
ve/cre yr(1) r
do i=1,[nx]
  ve/inp xr(1) $sigma(x([i]))
  do j=1,[ny]
    ve/inp yr(1) $sigma(y([j]))
*    ve/inp mapc([i],[j]) $call('accmappl.sl(xr,yr,0)')
    ve/inp mapc([i],[j]) $call('accmappl_blocks.sl(xr,yr,0)')
  enddo
enddo
1:
2d 2000 ! [nx] [xmin] [xmax] [ny] [ymin] [ymax]
ve/cre z([ny]) r
ve/cre d([ny]) r
do i=1,-[nx]
  ve/copy mapc([i],1:[ny]) z(1:[ny])
  * exec $PER/s#vpl z d y d sz=0.01
  graph [ny] y z sl
  read xy
enddo
*
do i=1,[nx]
  do j=1,[ny]
    tmp=mapc([i],[j])
    if ($sigma([tmp]).lt.0.001) then
      ve/inp mapc([i],[j]) 0
    endif
  enddo
enddo
hi/put/cont 2000  mapc
hi/put/err  2000 dmapc
return



macro s1s2 l=4.2 r=43.1 s=17.5 i1=1 i2=nx idh=0
ny=$vlen(yvi)
ymin=$sigma(vmin(yvi))
ymax=$sigma(vmax(yvi))
1d 100 ! [ny] [ymin] [ymax]
ve/cre  yes([ny]) r
ve/cre dyes([ny]) r
ve/cre  ye([ny]) r
ve/cre dye([ny]) r
ve/cre  yt([ny]) r
ve/cre dyt([ny]) r
nx=$vlen(xvi)
if ([i2].eq.'nx') then
  i2=[nx]
endif
do i=[i1],[i2]
  ve/copy  av([i],1:[ny])  yt
  ve/copy dav([i],1:[ny]) dyt
  exec hsigma ye = ye + yt/dyt**2
  exec hsigma dye = dye + 1/dyt**2
  exec hsigma yes = yes + yt
  exec hsigma dyes = dyes + dyt**2
enddo
sigma ye = ye / dye
sigma dye = 1/sqrt(dye)
di=[i2]-[i1]+1
exec hsigma dyes = sqrt(dyes)
sigma yes = yes/[di]
sigma dyes = dyes/[di]
ve/copy yes ye
ve/copy dyes dye
if ([idh].eq.0) then
  hi/put/cont 100  ye
  hi/put/err  100 dye
else
  hi/copy [idh] 100
  hi/get/cont 100 ye
  hi/get/err 100 dye
endif
ve/cre dyvi([ny]) r
** exec $PER/s#vpl ye dye yvi dyvi iatt=20 sz=0.1 ll=-1
exec vpl#pl0 ye dye yvi dyvi iatt=20 sz=0.1 ll=-1
ve/cre p5(5) r [l] $sigma(8) 0.0 0.0 0.5
ve/cre s5(5) r 0.1 0.1 0.01 0 0.01
lf=[l]-2.01
rf=[l]+5.01
hi/pl 100([lf]:[rf])
do i=1,5
  hi/fit 100([lf]:[rf]) fs.f s 5 p5 s5
  lf=$sigma(p5(1)-1.5*p5(5))
enddo
hi/pl 100([lf]:[rf])
*read x
gl/cre s1g $sigma(p5(1))
gl/cre sigs1g $sigma(abs(p5(5)))
ve/cre p5(5) r [r] $sigma(5) 0 0 0.5
lf=[r]-5.01
rf=[r]+2.01
hi/pl 100([lf]:[rf])
mess [r] [lf] [rf]
*read x
do i=1,5
  hi/fit 100([lf]:[rf]) fsx.f s 5 p5 s5
  rf=$sigma(p5(1)+1.5*p5(5))
enddo
hi/pl 100([lf]:[rf])
gl/cre s2g $sigma(p5(1))
gl/cre sigs2g $sigma(abs(p5(5)))
*read x
ve/cre p7(7) r [s] 0.5 $sigma(25) $sigma(8) 0 -0.1 0
lf=[s]-5.01
rf=[s]+5.01
hi/pl 100([lf]:[rf])
ve/cre s7(7) r 0.1 0 1 0 0 0
hi/fit 100([lf]:[rf]) fss.f sb 7 p7 s7
do i=1,5
  hi/fit 100([lf]:[rf]) fss.f s 7 p7
enddo
gl/cre ssg $sigma(p7(1))
gl/cre sigssg $sigma(abs(p7(2)))
*read x
return


macro z1z2 fs=17
nx=$vlen(xvi)
xmin=$sigma(vmin(xvi))
xmax=$sigma(vmax(xvi))
1d 100 ! [nx] [xmin] [xmax]
ve/cre  xe([nx]) r
ve/cre dxe([nx]) r
ve/cre  xt([nx]) r
ve/cre dxt([nx]) r
ny=$vlen(yvi)
yvl=[fs]-11
yvr=[fs]+11
do i=1,[ny]
  yv=yvi([i])
  if (([yv].lt.[yvl]).or.([yv].gt.[yvr])) then
    ve/copy  av(1:[nx],[i])  xt
    ve/copy dav(1:[nx],[i]) dxt
    exec hsigma xe = xe + xt
    exec hsigma dxe = dxe + dxt**2
  endif
enddo
exec hsigma dxe = sqrt(dxe)
hi/put/cont 100  xe
hi/put/err  100 dxe
ve/cre dxvi([nx]) r
* exec $PER/s#vpl xe dxe xvi dxvi iatt=20 sz=0.1 ll=-1
ve/cre p4(4) r 10 100 0 0
ve/cre s4(4) r 0.1 1 0.01 0
lf=4.01
rf=15.01
hi/pl 100([lf]:[rf])
hi/fit 100([lf]:[rf]) fz.f sb 4 p4 ! ! ! s4
hi/fit 100([lf]:[rf]) fz.f s 4 p4
hi/fit 100([lf]:[rf]) fz.f s 4 p4
gl/cre z2g $sigma(p4(1))
*read x
*
ve/cre p4(4) r -10 100 0 0
lf=-15.1
rf=-4.1
hi/pl 100([lf]:[rf])
hi/fit 100([lf]:[rf]) fz.f sb 4 p4 ! ! ! s4
hi/fit 100([lf]:[rf]) fz.f s 4 p4
hi/fit 100([lf]:[rf]) fz.f s 4 p4
gl/cre z1g $sigma(p4(1))
*read x
return


macro mapc y0=17.5
nx0=$vlen(xi)
ny0=$vlen(yi)
xmin=$sigma(vmin(xi))
xmax=$sigma(vmax(xi))
ymin=$sigma(vmin(yi))
ymax=$sigma(vmax(yi))
dx=$sigma(([xmax]-[xmin])/([nx0]-1)/2)
dy=$sigma(([ymax]-[ymin])/([ny0]-1)/2)

nx0=[nx0]-2
ny0=[ny0]-2

nx=$vlen(xvi)
ny=$vlen(yvi)
xvmin=$sigma(vmin(xvi))
xvmax=$sigma(vmax(xvi))
yvmin=$sigma(vmin(yvi))
yvmax=$sigma(vmax(yvi))

ve/cre  mape([nx0],[ny0]) r
ve/cre dmape([nx0],[ny0]) r
do i=1,[nx0]
  ii=[i]+1
  xl=xi([ii])-[dx]
  i1=$sigma(int(1+[nx]*([xl]-[xvmin])/([xvmax]-[xvmin])))
  xr=xi([ii])+[dx]
  i2=$sigma(int(1+[nx]*([xr]-[xvmin])/([xvmax]-[xvmin])))
  do j=1,[ny0]
    jj=[j]+1
    yl=yi([jj])-[dy]
    j1=$sigma(int(1+[ny]*([yl]-[yvmin])/([yvmax]-[yvmin])))
    yr=yi([jj])+[dy]
    j2=$sigma(int(1+[ny]*([yr]-[yvmin])/([yvmax]-[yvmin])))
    
    nxt=[i2]-[i1]+1
    nyt=[j2]-[j1]+1
    ve/del avt,davt
    ve/cre  avt([nxt],[nyt]) r
    ve/cre davt([nxt],[nyt]) r
    do jj=[j1],[j2]
      jjj=[jj]-[j1]+1
      ve/copy  av([i1]:[i2],[jj])  avt(1:[nxt],[jjj])
      ve/copy dav([i1]:[i2],[jj]) davt(1:[nxt],[jjj])
    enddo
    mean=$sigma(vsum(avt))
    sig=$sigma(sqrt(vsum(davt**2)))
    ve/inp  mape([i],[j]) $sigma([mean]/[nxt]/[nyt])
    ve/inp dmape([i],[j]) $sigma([sig]/[nxt]/[nyt])
  enddo
enddo

n=$vlen(yi)
ymin=yi(1)
ymax=yi([n])
iy=$sigma(int(1.5+([n]-1)*([y0]-[ymin])/([ymax]-[ymin])))
i0=[iy]-1
i1=[i0]-1
i2=[i0]+1
n=$vlen(yi)
ve/del mape0,mape1,mape2
ve/copy mape(1:[nx0],[i1]) mape1
ve/copy mape(1:[nx0],[i2]) mape2
sigma mape0 = (mape1+mape2)/2
ve/copy mape0 mape(1:[nx0],[i0])

return


macro glmapx
do i=1,9
  exec mapcal#glmap [i]
enddo
return


macro glmap ncnt=1

chain -exp
chain exp mhad2011_1000_col_p1.hbook           
chain exp mhad2011_1000_col_p2.hbook           
chain exp mhad2011_1000_col_p3.hbook           
chain exp mhad2011_508.7_col.hbook             
chain exp mhad2011_509.8_col.hbook             
chain exp mhad2011_510.5-2_col_p1.hbook        
chain exp mhad2011_510.5-2_col_p2.hbook        
chain exp mhad2011_510.5-2_col_p3.hbook        
chain exp mhad2011_510.5-2_col_p4.hbook        
chain exp mhad2011_525_col_p1.hbook            
chain exp mhad2011_525_col_p2.hbook            
chain exp mhad2011_525_col_p3.hbook            
chain exp mhad2011_525_col_p4.hbook            
chain exp mhad2011_537.5_col_p1.hbook          
chain exp mhad2011_537.5_col_p2.hbook          
chain exp mhad2011_537.5_col_p3.hbook          
chain exp mhad2011_550_col_p1.hbook            
chain exp mhad2011_550_col_p2.hbook            
chain exp mhad2011_550_col_p3.hbook            
chain exp mhad2011_562.5_col_p1.hbook          
chain exp mhad2011_562.5_col_p2.hbook          
chain exp mhad2011_575_col_p1.hbook            
chain exp mhad2011_575_col_p2.hbook            
chain exp mhad2011_587.5_col_p1.hbook          
chain exp mhad2011_587.5_col_p2.hbook          
chain exp mhad2011_600_col_p1.hbook            
chain exp mhad2011_600_col_p2.hbook            
chain exp mhad2011_612.5_col_p1.hbook          
chain exp mhad2011_612.5_col_p2.hbook          
chain exp mhad2011_625_col_p1.hbook            
chain exp mhad2011_625_col_p2.hbook            
chain exp mhad2011_637.5_col_p1.hbook          
chain exp mhad2011_637.5_col_p2.hbook          
chain exp mhad2011_650_col_p1.hbook            
chain exp mhad2011_650_col_p2.hbook            
chain exp mhad2011_650_col_p3.hbook            
chain exp mhad2011_662.5_col_p1.hbook          
chain exp mhad2011_662.5_col_p2.hbook          
chain exp mhad2011_675_col_p1.hbook            
chain exp mhad2011_675_col_p2.hbook            
chain exp mhad2011_687.5_col_p1.hbook          
chain exp mhad2011_687.5_col_p2.hbook          
chain exp mhad2011_700_col_p1.hbook            
chain exp mhad2011_700_col_p2.hbook            
chain exp mhad2011_712.5_col_p1.hbook          
chain exp mhad2011_712.5_col_p2.hbook          
chain exp mhad2011_725_col_p1.hbook            
chain exp mhad2011_725_col_p2.hbook            
chain exp mhad2011_737.5_col_p1.hbook          
chain exp mhad2011_737.5_col_p2.hbook          
chain exp mhad2011_750-1_col_p1.hbook          
chain exp mhad2011_750-1_col_p2.hbook          
chain exp mhad2011_750_col_p1.hbook            
chain exp mhad2011_750_col_p2.hbook            
chain exp mhad2011_762.5_col_p1.hbook          
chain exp mhad2011_762.5_col_p2.hbook          
chain exp mhad2011_775_col_p1.hbook            
chain exp mhad2011_775_col_p2.hbook            
chain exp mhad2011_787.5_col_p1.hbook          
chain exp mhad2011_787.5_col_p2.hbook          
chain exp mhad2011_800_col_p1.hbook            
chain exp mhad2011_800_col_p2.hbook            
chain exp mhad2011_812.5_col.hbook             
chain exp mhad2011_825_col_p1.hbook            
chain exp mhad2011_825_col_p2.hbook            
chain exp mhad2011_825_col_p3.hbook            
chain exp mhad2011_837.5.hbook                 
chain exp mhad2011_850_col_p1.hbook            
chain exp mhad2011_850_col_p2.hbook            
chain exp mhad2011_862.5_col.hbook             
chain exp mhad2011_875_col_p1.hbook            
chain exp mhad2011_875_col_p2.hbook            
chain exp mhad2011_887.5_col.hbook             
chain exp mhad2011_900_col.hbook               
chain exp mhad2011_912.5_col.hbook             
chain exp mhad2011_925_col_p1.hbook            
chain exp mhad2011_925_col_p2.hbook            
chain exp mhad2011_935_col_p1.hbook            
chain exp mhad2011_935_col_p2.hbook            
chain exp mhad2011_945_col_p1.hbook            
chain exp mhad2011_945_col_p2.hbook            
chain exp mhad2011_950_col_p1.hbook            
chain exp mhad2011_950_col_p2.hbook            
chain exp mhad2011_950_col_p3.hbook            
chain exp mhad2011_962.5_col_p1.hbook          
chain exp mhad2011_962.5_col_p2.hbook          
chain exp mhad2011_962.5_col_p3.hbook          
chain exp mhad2011_987.5_col_p1.hbook          
chain exp mhad2011_987.5_col_p2.hbook          
chain exp mhad2011_987.5_col_p3.hbook          
chain exp -P /work/users/konctbel/exp/MHAD2011-2/col/

dir=.
datafile=global_map_counter_[ncnt].dat

gl/cre sc $sigma(180/3.1415927)
r1=10.666
r2=13.834
r=([r1]+[r2])/2

mnz=80
mzmin=-20
mzmax=20
mnf=240
mfmax=$sigma((40*[ncnt]-20)+40)
mfmin=$sigma((40*[ncnt]-20)-40)
dz=$sigma(([mzmax]-[mzmin])/[mnz])
df=$sigma(([mfmax]-[mfmin])/[mnf])
ve/cre xvi([mnz]) r
ve/cre yvi([mnf]) r
l0=[mzmin]+[dz]/2.
r0=[mzmax]-[dz]/2.
sigma xvi = array([mnz],[l0]#[r0])
l0=[mfmin]+[df]/2.
r0=[mfmax]-[df]/2.
sigma yvi = array([mnf],[l0]#[r0])

cuts=$1.and.$3.and.$4.and.$5.and.$6.and.$8.and.$9.and.eventtime.eq.0

profile 1000 ! [mnz] [mzmin] [mzmax] ! !
profile 2000 ! [mnz] [mzmin] [mzmax] ! !
profile 3000 ! [mnz] [mzmin] [mzmax] ! !

ve/cre  avi([mnz]) r
ve/cre davi([mnz]) r
ve/cre  av([mnz],[mnf]) r
ve/cre dav([mnz],[mnf]) r

do i=1,[mnf]

  phi1=$sigma([mfmin]+[df]*([i]-0.5))
  phi2=[phi1]+360
  phi0=[phi1]-360

  nt/pl //exp/1.ach([ncnt])*sin(theta(1))%(z0(1)+[r]/tan(theta(1))) _
  [cuts].and.(abs(phi(1)*[sc]-([phi0]))<0.5*[df].or.abs(phi(1)*[sc]-([phi1]))<0.5*[df].or.abs(phi(1)*[sc]-([phi2]))<0.5*[df]) idh=1000
  
  nt/pl //exp/1.ach([ncnt])*sin(theta(2))%(z0(2)+[r]/tan(theta(2))) _
  [cuts].and.(abs(phi(2)*[sc]-([phi0]))<0.5*[df].or.abs(phi(2)*[sc]-([phi1]))<0.5*[df].or.abs(phi(2)*[sc]-([phi2]))<0.5*[df]) idh=2000
  
  hi/op/add 1000 2000 3000 1 1 e
	
  hi/get/cont 3000  avi
  hi/get/err  3000 davi
	
  ve/copy  avi(1:[mnz])  av(1:[mnz],[i])
  ve/copy davi(1:[mnz]) dav(1:[mnz],[i])

enddo
ve/del avi,davi
hi/del 1000
hi/del 2000
hi/del 3000

*ve/write av,dav [dir]/mapfit.txt 2g15.6
ve/write av,dav [dir]/[datafile] 2g15.6

*
ve/del av,dav
ve/cre  av([mnz],[mnf]) r
ve/cre dav([mnz],[mnf]) r
*ve/read av,dav [dir]/mapfit.txt 2g15.6
ve/read av,dav [dir]/[datafile] 2g15.6

nsc=[mnz]*[mnf]
ve/del av1,dav1
ve/read av1,dav1 [dir]/[datafile] 2g15.6

2d 1000 ! [mnz] [mzmin] [mzmax] [mnf] [mfmin] [mfmax]

hi/put/cont 1000  av
hi/put/err  1000 dav

*------------------

return


macro mapx fout=mhad2011_510.5/counter_1/mapfit.out
*  shell /work/snd2000/users/konctbel/SNDSimACC/readlifstat_new mapfit.out
  shell /work/snd2000/users/konctbel/SNDSimACC/readlifstat [fout]
  ve/del pars0,dpars0
  ve/read pars0,dpars0 [fout].pars

  i1=2
  i2=[npar]+1
  ve/copy pars0([i1]:[i2]) parn(1:[npar])
  i=1
  i=[i]+1; z1=pars0([i]);    dz1=dpars0([i])
  i=[i]+1; z2=pars0([i]);    dz2=dpars0([i])
  i=[i]+1; s1=pars0([i]);    ds1=dpars0([i])
  i=[i]+1; s2=pars0([i]);    ds2=dpars0([i])
  i=[i]+1; sigz=pars0([i]);  dsigz=dpars0([i])
  i=[i]+1; sigs=pars0([i]);  dsigs=dpars0([i])
  i=[i]+1; pds=pars0([i]);   dpds=dpars0([i])
  i=[i]+1; rPMT=pars0([i]);  drPMT=dpars0([i])
  i=[i]+1; apmt=pars0([i]);  daPMT=dpars0([i])
  i=[i]+1; zs1=pars0([i]);   dzs1=dpars0([i])
  i=[i]+1; zs2=pars0([i]);   dzs2=dpars0([i])
  i=[i]+1; dzds=pars0([i]);  ddzds=dpars0([i])
  i=[i]+1; ss=pars0([i]);    dss=dpars0([i])
  i=[i]+1; alpha=pars0([i]); dalpha=dpars0([i])
  i=[i]+1; z0=pars0([i]);    dz0=dpars0([i]); dz0=0.1
  i=[i]+1; sPMT=pars0([i]);    dsPMT=dpars0([i])
  i=[i]+1; sigs0=pars0([i]);  dsigs0=dpars0([i])
  i=[i]+1; ps1=pars0([i]);  dps1=dpars0([i])
  i=[i]+1; ps2=pars0([i]);  dps2=dpars0([i])
  i=[i]+1; ts1=pars0([i]);  dts1=dpars0([i])
  i=[i]+1; ts2=pars0([i]);  dts2=dpars0([i])
  i1=[i]+1;
  i2=[i1]+[nx0]-1
  ve/copy  pars0([i1]:[i2])  wls(1:[nx0])
  ve/copy dpars0([i1]:[i2]) dwls(1:[nx0])
  do i=1,[ny0]
    i1=[i2]+1
    i2=[i1]+[nx0]-1
    ve/copy  pars0([i1]:[i2])  map(1:[nx0],[i])
    ve/copy dpars0([i1]:[i2]) dmap(1:[nx0],[i])
  enddo
  gl/imp nxg
  gl/imp nyg
  nx=[nxg]
  ny=[nyg]
  nx0=[nx]-2
  ny0=[ny]-2
  ind0=[npar0]+[nx0]
  nx0st=1
  nx0sp=[nx0]
  npar=[npar0]+[nx0]+[nx0]*[ny0]
  ve/cre wls([nx0]) r [nx0]*5
  ve/cre dwls([nx0]) r [nx0]*0.3
  nxy0=[nx0]*[ny0]
  ve/cre map([nx0],[ny0]) r [nxy0]*2
  ve/cre dmap([nx0],[ny0]) r [nxy0]*0.1
  ve/cre xi([nx]) r
  ve/cre yi([ny]) r
  sigma xi = array([nx],-12#12)
  l0=5+40*([ncnt]-1)
  r0=43+40*([ncnt]-1)
  sigma yi = array([ny],[l0]#[r0])
  ve/cre xr(1) r
  ve/cre yr(1) r
  do i=1,[nx0]
    ve/inp xr(1) $sigma(xi([i]))
    ve/inp wls([i]) $call('accmappl.sl(xr,yr,2)')
    do j=1,[ny0]
      ve/inp yr(1) $sigma(yi([j]))
      ve/inp map([i],[j]) $call('accmappl.sl(xr,yr,1)')
    enddo
  enddo
  ve/cre parn([npar]) r [z1] [z2] [s1] [s2] [sigz] [sigs] [pds] [rpmt] [apmt] [zs1] [zs2] [dzds] [ss] [z0] [sPMT] [sigs0]
  i1=[ind0]-[nx0]+1
  i2=[ind0]
  ve/copy wls(1:[nx0]) parn([i1]:[i2])
  do j=1,[ny0]
    i1=[i1]+[nx0]
    i2=[i2]+[nx0]
    ve/copy map(1:[nx0],[j]) parn([i1]:[i2])
  enddo
return



macro cmp ncnt=1

*goto 1
chain -exp
chain exp mhad2011_510.5-2_col_p1.hbook
chain exp mhad2011_510.5-2_col_p2.hbook
chain exp mhad2011_510.5-2_col_p3.hbook
chain exp mhad2011_510.5-2_col_p4.hbook
chain exp -P /work/users/konctbel/exp/MHAD2011-2/col/

chain -sim
do i=2500,2519
*  chain sim ../snd2k/R005-999/ee_t10e10-510-10943-[i]-50000.hbook
enddo
do i=2500,2519
*  chain sim ../snd2k/R005-999/ee_t35e10-510-10943-[i]-100000_corr.hbook
  chain sim ../snd2k/R005-999/ee_t35e10-510-10943-[i]-100000.hbook
enddo

chain -smun
do i=101,110
  chain smun ../snd2k/R005-999/smun_nrc-105500-10943-[i]-100000.hbook
enddo
chain -smup
do i=111,120
  chain smup ../snd2k/R005-999/smup_nrc-105500-10943-[i]-100000.hbook
enddo

chain -mumu
do i=1271,1280
  chain mumu ../snd2k/R005-999/mumu_ps_nrc-510-10943-[i]-100000.hbook
enddo

sc=180/3.1415927
dphi=([ncnt]-1)*40

exec mapcal#cuts

cut $70 $1.and.$3.and.$4.and.$5.and.$6.and.$8.and.$9.and.eventtime.eq.0
cut $72 $1.and.$4.and.$5.and.$6.and.$9.and.eventtime.eq.0
cut $71 200<energy(1)<400

nz=120

profile 110 ! [nz] -15 15 ! !
profile 210 ! [nz] -15 15 ! !
profile 310 ! [nz] -15 15 ! !
profile 410 ! [nz] -15 15 ! !
profile 111 ! 80 $sigma(([ncnt]-1)*40-20) $sigma([ncnt]*40+20) ! !
profile 211 ! 80 $sigma(([ncnt]-1)*40-20) $sigma([ncnt]*40+20) ! !
profile 311 ! 80 $sigma(([ncnt]-1)*40-20) $sigma([ncnt]*40+20) ! !
profile 411 ! 80 $sigma(([ncnt]-1)*40-20) $sigma([ncnt]*40+20) ! !
profile 120 ! [nz] -15 15 ! !
profile 220 ! [nz] -15 15 ! !
profile 320 ! [nz] -15 15 ! !
profile 420 ! [nz] -15 15 ! !
profile 121 ! 80 $sigma(([ncnt]-1)*40-20) $sigma([ncnt]*40+20) ! !
profile 221 ! 80 $sigma(([ncnt]-1)*40-20) $sigma([ncnt]*40+20) ! !
profile 321 ! 80 $sigma(([ncnt]-1)*40-20) $sigma([ncnt]*40+20) ! !
profile 421 ! 80 $sigma(([ncnt]-1)*40-20) $sigma([ncnt]*40+20) ! !
profile 100 ! [nz] -15 15 ! !
profile 200 ! [nz] -15 15 ! !
profile 300 ! [nz] -15 15 ! !
profile 400 ! [nz] -15 15 ! !
profile 101 ! 80 $sigma(([ncnt]-1)*40-20) $sigma([ncnt]*40+20) ! !
profile 201 ! 80 $sigma(([ncnt]-1)*40-20) $sigma([ncnt]*40+20) ! !
profile 301 ! 80 $sigma(([ncnt]-1)*40-20) $sigma([ncnt]*40+20) ! !
profile 401 ! 80 $sigma(([ncnt]-1)*40-20) $sigma([ncnt]*40+20) ! !

nt/pl //exp/1.ach([ncnt])*sin(theta(1))%phi(1)*[sc] $70.and.abs(z0(1)+12/tan(theta(1)))<10 idh=111
nt/pl //sim/1.ach([ncnt])*sin(theta(1))%phi(1)*[sc] $70.and.abs(z0(1)+12/tan(theta(1)))<10 idh=211
nt/pl //mumu/1.ach([ncnt])*sin(theta(1))%phi(1)*[sc] $72.and.abs(z0(1)+12/tan(theta(1)))<10 idh=411
nt/pl //smun/1.ach([ncnt])*sin(theta(1))%phi(1)*[sc] $71.and.abs(z0(1)+12/tan(theta(1)))<10 idh=311
nt/pl //exp/1.ach([ncnt])*sin(theta(2))%phi(2)*[sc] $70.and.abs(z0(2)+12/tan(theta(2)))<10 idh=121
nt/pl //sim/1.ach([ncnt])*sin(theta(2))%phi(2)*[sc] $70.and.abs(z0(2)+12/tan(theta(2)))<10 idh=221
nt/pl //mumu/1.ach([ncnt])*sin(theta(2))%phi(2)*[sc] $72.and.abs(z0(2)+12/tan(theta(2)))<10 idh=421
nt/pl //smup/1.ach([ncnt])*sin(theta(1))%phi(1)*[sc] $71.and.abs(z0(1)+12/tan(theta(1)))<10 idh=321

nt/pl //exp/1.ach([ncnt])*sin(theta(1))%z0(1)+12/tan(theta(1)) $70.and.abs(phi(1)*[sc]-22.5-[dphi])<16.and.abs(phi(1)*[sc]-17.5-[dphi])>3 idh=110
nt/pl //sim/1.ach([ncnt])*sin(theta(1))%z0(1)+12/tan(theta(1)) $70.and.abs(phi(1)*[sc]-22.5-[dphi])<16.and.abs(phi(1)*[sc]-17.5-[dphi])>3 idh=210
nt/pl //mumu/1.ach([ncnt])*sin(theta(1))%z0(1)+12/tan(theta(1)) $72.and.abs(phi(1)*[sc]-22.5-[dphi])<16.and.abs(phi(1)*[sc]-17.5-[dphi])>3 idh=410
nt/pl //smun/1.ach([ncnt])*sin(theta(1))%z0(1)+12/tan(theta(1)) $71.and.abs(phi(1)*[sc]-22.5-[dphi])<16.and.abs(phi(1)*[sc]-17.5-[dphi])>3 idh=310
nt/pl //exp/1.ach([ncnt])*sin(theta(2))%z0(2)+12/tan(theta(2)) $70.and.abs(phi(2)*[sc]-22.5-[dphi])<16.and.abs(phi(2)*[sc]-17.5-[dphi])>3 idh=120
nt/pl //sim/1.ach([ncnt])*sin(theta(2))%z0(2)+12/tan(theta(2)) $70.and.abs(phi(2)*[sc]-22.5-[dphi])<16.and.abs(phi(2)*[sc]-17.5-[dphi])>3 idh=220
nt/pl //mumu/1.ach([ncnt])*sin(theta(2))%z0(2)+12/tan(theta(2)) $72.and.abs(phi(2)*[sc]-22.5-[dphi])<16.and.abs(phi(2)*[sc]-17.5-[dphi])>3 idh=420
nt/pl //smup/1.ach([ncnt])*sin(theta(1))%z0(1)+12/tan(theta(1)) $71.and.abs(phi(1)*[sc]-22.5-[dphi])<16.and.abs(phi(1)*[sc]-17.5-[dphi])>3 idh=320

*nt/pl //exp/1.ach([ncnt])*sin(theta(1))%z0(1)+12/tan(theta(1))  $70.and.(phi(1)*[sc]-22.5-[dphi])<0.and.abs(phi(1)*[sc]-22.5-[dphi])<16.and.abs(phi(1)*[sc]-17.5-[dphi])>3 idh=110
*nt/pl //sim/1.ach([ncnt])*sin(theta(1))%z0(1)+12/tan(theta(1))  $70.and.(phi(1)*[sc]-22.5-[dphi])<0.and.abs(phi(1)*[sc]-22.5-[dphi])<16.and.abs(phi(1)*[sc]-17.5-[dphi])>3 idh=210
*nt/pl //smun/1.ach([ncnt])*sin(theta(1))%z0(1)+12/tan(theta(1)) $71.and.(phi(1)*[sc]-22.5-[dphi])<0.and.abs(phi(1)*[sc]-22.5-[dphi])<16.and.abs(phi(1)*[sc]-17.5-[dphi])>3 idh=310
*nt/pl //exp/1.ach([ncnt])*sin(theta(2))%z0(2)+12/tan(theta(2))  $70.and.(phi(2)*[sc]-22.5-[dphi])<0.and.abs(phi(2)*[sc]-22.5-[dphi])<16.and.abs(phi(2)*[sc]-17.5-[dphi])>3 idh=120
*nt/pl //sim/1.ach([ncnt])*sin(theta(2))%z0(2)+12/tan(theta(2))  $70.and.(phi(2)*[sc]-22.5-[dphi])<0.and.abs(phi(2)*[sc]-22.5-[dphi])<16.and.abs(phi(2)*[sc]-17.5-[dphi])>3 idh=220
*nt/pl //smup/1.ach([ncnt])*sin(theta(1))%z0(1)+12/tan(theta(1)) $71.and.(phi(1)*[sc]-22.5-[dphi])<0.and.abs(phi(1)*[sc]-22.5-[dphi])<16.and.abs(phi(1)*[sc]-17.5-[dphi])>3 idh=320

hi/op/add 110 120 100
hi/op/add 210 220 200
hi/op/add 310 320 300
hi/op/add 410 420 400
hi/op/add 111 121 101
hi/op/add 211 221 201
hi/op/add 311 321 301
hi/op/add 411 421 401

1d 112 ! 200 -5 35
1d 212 ! 200 -5 35
1d 312 ! 200 -5 35
1d 412 ! 200 -5 35
1d 122 ! 200 -5 35
1d 222 ! 200 -5 35
1d 322 ! 200 -5 35
1d 422 ! 200 -5 35
1d 102 ! 200 -5 35
1d 202 ! 200 -5 35
1d 302 ! 200 -5 35
1d 402 ! 200 -5 35

nt/pl //exp/1.ach([ncnt])*sin(theta(1)) $70.and.abs(z0(1)+12/tan(theta(1)))<10.and.abs(phi(1)*[sc]-22.5-[dphi])<16.and.abs(phi(1)*[sc]-17.5-[dphi])>3 idh=112
nt/pl //sim/1.ach([ncnt])*sin(theta(1)) $70.and.abs(z0(1)+12/tan(theta(1)))<10.and.abs(phi(1)*[sc]-22.5-[dphi])<16.and.abs(phi(1)*[sc]-17.5-[dphi])>3 idh=212
nt/pl //mumu/1.ach([ncnt])*sin(theta(1)) $72.and.abs(z0(1)+12/tan(theta(1)))<10.and.abs(phi(1)*[sc]-22.5-[dphi])<16.and.abs(phi(1)*[sc]-17.5-[dphi])>3 idh=412
nt/pl //smun/1.ach([ncnt])*sin(theta(1)) $71.and.abs(z0(1)+12/tan(theta(1)))<10.and.abs(phi(1)*[sc]-22.5-[dphi])<16.and.abs(phi(1)*[sc]-17.5-[dphi])>3 idh=312
nt/pl //exp/1.ach([ncnt])*sin(theta(2)) $70.and.abs(z0(2)+12/tan(theta(2)))<10.and.abs(phi(2)*[sc]-22.5-[dphi])<16.and.abs(phi(2)*[sc]-17.5-[dphi])>3 idh=122
nt/pl //sim/1.ach([ncnt])*sin(theta(2)) $70.and.abs(z0(2)+12/tan(theta(2)))<10.and.abs(phi(2)*[sc]-22.5-[dphi])<16.and.abs(phi(2)*[sc]-17.5-[dphi])>3 idh=222
nt/pl //mumu/1.ach([ncnt])*sin(theta(2)) $72.and.abs(z0(2)+12/tan(theta(2)))<10.and.abs(phi(2)*[sc]-22.5-[dphi])<16.and.abs(phi(2)*[sc]-17.5-[dphi])>3 idh=422
nt/pl //smup/1.ach([ncnt])*sin(theta(1)) $71.and.abs(z0(1)+12/tan(theta(1)))<10.and.abs(phi(1)*[sc]-22.5-[dphi])<16.and.abs(phi(1)*[sc]-17.5-[dphi])>3 idh=322

hi/op/add 112 122 102
hi/op/add 212 222 202
hi/op/add 312 322 302
hi/op/add 412 422 402
exec hsigma @1202 = @202/vsum(@202)*vsum(@102)
exec hsigma %1202 = %202/vsum(@202)*vsum(@102)
exec hsigma @1302 = @302/vsum(@302)*vsum(@102)
exec hsigma %1302 = %302/vsum(@302)*vsum(@102)

zone 2 1
exec seteps 0
opt ngrid 
opt nstat
set mtyp 20

set pmci 4
hi/pl 101 
set pmci 2
hi/pl 201 s
set pmci 1
hi/pl 301 s

set pmci 4
hi/pl 100 
set pmci 2
hi/pl 200 s
set pmci 1
hi/pl 300 s

exec hsigma @1200 = @200/@100
exec hsigma %1200 = @1200*sqrt((%100/@100)**2+(%200/@200)**2)
ve/cre p0(1) r 1
ve/cre dp0(1) r 0
hi/fit 1200(-8.:8.) p0 s 1 p0 ! ! ! dp0
ve/inp p0i([ncnt]) p0(1)
ve/inp dp0i([ncnt]) dp0(1)

set pmci 1
1:
exec mapcal#readlif [ncnt]

ve/cre mod(2) r 1 1
s1=pars0(4)
s2=pars0(5)
ss=pars0(14)
sm=$sigma(([s2]+[ss])/2)
sig=pars0(7)
ve/cre xi(10) r [s1] $sigma([s1]+([ss]-[s1])/3) $sigma([s1]+([ss]-[s1])*2/3) _
                [ss] $sigma([ss]+([sm]-[ss])/3) $sigma([ss]+([sm]-[ss])*2/3) _
                [sm] $sigma([sm]+([s2]-[sm])/3) $sigma([sm]+([s2]-[sm])*2/3) [s2]
exec hsigma ym = vmax(@101)                
ve/cre p16(16) r 6 8 9 11 9 7 4.7 4.7 4.3 3 [s1] [ss] [sm] [s2] $sigma(ym(1)) [sig]
do i=1,5
  hi/fit 101 scphi.f s 16 p16
enddo
do i=1,5
  l=$sigma([s1]-0*[sig])
  r=$sigma([s2]+0*[sig])
  s1=p16(11)
  ss=p16(12)
  sm=p16(13)
  s2=p16(14)
  ve/cre xi(10) r [s1] $sigma([s1]+([ss]-[s1])/3) $sigma([s1]+([ss]-[s1])*2/3) _
                  [ss] $sigma([ss]+([sm]-[ss])/3) $sigma([ss]+([sm]-[ss])*2/3) _
                  [sm] $sigma([sm]+([s2]-[sm])/3) $sigma([sm]+([s2]-[sm])*2/3) [s2]
  hi/fit 101([l]:[r]) scphi.f s 16 p16
enddo
do i=1,10
  hi/fit 101([l]:[r]) scphi.f s 16 p16
enddo

return

macro cmpx
ve/cre  p0i(9) r
ve/cre dp0i(9) r
do i=1,9
  exec mapcal#cmp [i]
enddo
ve/cre ns(9) r  1 2 3 4 5 6 7 8 9
ve/cre dns(9) r
* exec $PER/s#vpl p0i dp0i ns dns sz=0.1
ve/fit ns p0i dp0i p0 s
return

macro rimsx

ve/cre  s1(9) r
ve/cre ds1(9) r
ve/cre  s2(9) r
ve/cre ds2(9) r
ve/cre  ss(9) r
ve/cre dss(9) r
n=0
n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out
n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out
n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out
n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out
n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out.cut
n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out.cut
n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out.cut
n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out.cut
n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out
do i=1,9
  fout=[file[i]]
  mess [fout]
  shell /work/snd2000/users/konctbel/SNDSimACC/readlifstat [fout]
  ve/del pars0,dpars0
  ve/read pars0,dpars0 [fout].pars
  ve/inp s1([i]) pars0(4)
  ve/inp s2([i]) pars0(5)
  ve/inp ss([i]) pars0(14)
  ve/inp ds1([i]) dpars0(4)
  ve/inp ds2([i]) dpars0(5)
  ve/inp dss([i]) dpars0(14)
enddo
ve/cre ns(9) r 1 2 3 4 5 6 7 8 9
ve/cre dns(9) r

return


macro readlif ncnt=1
n=0
n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out
n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out
n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out
n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out
n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out.cut
n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out.cut
n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out.cut
n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out.cut
n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out
fout=[file[ncnt]]
shell /work/snd2000/users/konctbel/SNDSimACC/readlifstat [fout]
ve/del pars0,dpars0
ve/read pars0,dpars0 [fout].pars
return



macro mappr n1s=1 n2s=41
n=0
n=[n]+1; ebeam[n]=508.7 ; il[n]=0
n1=[n]
n=[n]+1; ebeam[n]=509.8 ; il[n]=0
n=[n]+1; ebeam[n]=510.5 ; il[n]=0
n=[n]+1; ebeam[n]=525   ; il[n]=363.2
n=[n]+1; ebeam[n]=537.5 ; il[n]=546.1
n=[n]+1; ebeam[n]=550   ; il[n]=485.7
n=[n]+1; ebeam[n]=562.5 ; il[n]=517.3
n=[n]+1; ebeam[n]=575   ; il[n]=419.1
n=[n]+1; ebeam[n]=587.5 ; il[n]=534.4
n=[n]+1; ebeam[n]=600   ; il[n]=489.8999940
n=[n]+1; ebeam[n]=612.5 ; il[n]=546.7000120
n=[n]+1; ebeam[n]=625   ; il[n]=440.5000000
n=[n]+1; ebeam[n]=637.5 ; il[n]=496.2000120
n=[n]+1; ebeam[n]=650   ; il[n]=456.1000060
n=[n]+1; ebeam[n]=662.5 ; il[n]=526.2999880
n=[n]+1; ebeam[n]=675   ; il[n]=559.7999880
n=[n]+1; ebeam[n]=687.5 ; il[n]=576.2000120
n=[n]+1; ebeam[n]=700   ; il[n]=577.0999760
n=[n]+1; ebeam[n]=712.5 ; il[n]=593.2000120
n=[n]+1; ebeam[n]=725   ; il[n]=432.2999880
n=[n]+1; ebeam[n]=737.5 ; il[n]=600.0999760
n=[n]+1; ebeam[n]=750   ; il[n]=694.7000120
n=[n]+1; ebeam[n]=762.5 ; il[n]=487.2999880
n=[n]+1; ebeam[n]=775   ; il[n]=546.5000000
n=[n]+1; ebeam[n]=787.5 ; il[n]=502.3999940
n=[n]+1; ebeam[n]=800   ; il[n]=439.3999940
n=[n]+1; ebeam[n]=812.5 ; il[n]=509.8999940
n=[n]+1; ebeam[n]=825   ; il[n]=475.6000060
n=[n]+1; ebeam[n]=850   ; il[n]=462.2000120
n=[n]+1; ebeam[n]=862.5 ; il[n]=504.6000060
n=[n]+1; ebeam[n]=875   ; il[n]=506.3999940
n=[n]+1; ebeam[n]=887.5 ; il[n]=525.7000120
n=[n]+1; ebeam[n]=900   ; il[n]=384.1000060
n=[n]+1; ebeam[n]=912.5 ; il[n]=481.3999940
n=[n]+1; ebeam[n]=925   ; il[n]=401.2999880
n=[n]+1; ebeam[n]=935   ; il[n]=637.5000000
n=[n]+1; ebeam[n]=945   ; il[n]=577.0999760
n=[n]+1; ebeam[n]=950   ; il[n]=451.2999880
n=[n]+1; ebeam[n]=962.5 ; il[n]=561.9000240
n=[n]+1; ebeam[n]=987.5 ; il[n]=467.5000000
n=[n]+1; ebeam[n]=1000  ; il[n]=538.0999760
n2=[n]

cut $70 $1.and.$3.and.$4.and.$5.and.$6.and.$8.and.$9.and.eventtime.eq.0
cut $71 $70.and.abs(theta(1)-1.5708)<0.02
cut $72 $70.and.abs(theta(2)-1.5708)<0.02
cut $71 $70
cut $72 $70
*cut $72 $1.and.$4.and.$5.and.$6.and.$9.and.eventtime.eq.0
*cut $71 200<energy(1)<400
dir=/work/users/konctbel/SepPar
nz=300
nf=320

sc=180/3.1415927

phir1=phi(1)*[sc]
phir2=phi(2)*[sc]

phir1=phir.f(1)*[sc]
phir2=phir.f(2)*[sc]

rm=12.3

do i=[n1s],[n2s]
  ebeam=[ebeam[i]]
  exec [dir]/sp#hists [ebeam]
  nt/pl //exp/1.eton $1.and.eton<1.3
  exec [dir]/sp#hists [ebeam]
  hfile=hist_[ebeam]_prof_v2.his
  if ($fexist([hfile]).eq.1) then
    shell rm -v [hfile]
  endif
  hi/file 20 [hfile] ! N
  do ncnt=1,9
    ve/del p16
    ve/read p16 mhad2011_amplitude_vs_phi_counter[ncnt]_fit.par 'f15.10'
    s1=p16(12)
    ss=p16(13)
    sm=p16(14)
    s2=p16(15)
    sig=p16(17)
    sigs=p16(18)
    x1=$sigma([s1]+2*[sig])
    x2=$sigma([s2]-2*[sig])
    x1s=$sigma([ss]-1.5/120.0*180/3.1415927-3*[sigs])
    x2s=$sigma([ss]+1.5/120.0*180/3.1415927+3*[sigs])
*
    dphi=([ncnt]-1)*40
    nh0=100+[ncnt]
    profile [nh0] ! [nz] -15 15 ! !
    nh1=110+[ncnt]
    profile [nh1] ! [nz] -15 15 ! !
    nh2=120+[ncnt]
    profile [nh2] ! [nz] -15 15 ! !
    ng0=200+[ncnt]
    profile [ng0] ! [nf] $sigma(([ncnt]-1)*40-20) $sigma([ncnt]*40+20) ! !
    ng1=210+[ncnt]
    profile [ng1] ! [nf] $sigma(([ncnt]-1)*40-20) $sigma([ncnt]*40+20) ! !
    ng2=220+[ncnt]
    profile [ng2] ! [nf] $sigma(([ncnt]-1)*40-20) $sigma([ncnt]*40+20) ! !
    nt/pl //exp/1.ach([ncnt])*sin(theta(1))%[phir1] $71.and.abs(z0(1)+[rm]/tan(theta(1)))<10 idh=[ng1]
    nt/pl //exp/1.ach([ncnt])*sin(theta(2))%[phir2] $72.and.abs(z0(2)+[rm]/tan(theta(2)))<10 idh=[ng2]
    nt/pl //exp/1.ach([ncnt])*sin(theta(1))%z0(1)+[rm]/tan(theta(1)) $71.and.abs([phir1]-22.5-[dphi])<18.and.abs([phir1]-17.5-[dphi])>3 idh=[nh1]
    nt/pl //exp/1.ach([ncnt])*sin(theta(2))%z0(2)+[rm]/tan(theta(2)) $72.and.abs([phir2]-22.5-[dphi])<18.and.abs([phir2]-17.5-[dphi])>3 idh=[nh2]
    nh01=1100+[ncnt]
    profile [nh01] ! [nz] -15 15 ! !
    nh11=1110+[ncnt]
    profile [nh11] ! [nz] -15 15 ! !
    nh21=1120+[ncnt]
    profile [nh21] ! [nz] -15 15 ! !
    nt/pl //exp/1.ach([ncnt])*sin(theta(1))%z0(1)+[rm]/tan(theta(1)) $71.and.[x1]<[phir1]<[x1s] idh=[nh11]
    nt/pl //exp/1.ach([ncnt])*sin(theta(2))%z0(2)+[rm]/tan(theta(2)) $72.and.[x1]<[phir2]<[x1s] idh=[nh21]
    hi/op/add [nh11] [nh21] [nh01]    
    nh02=2100+[ncnt]
    profile [nh02] ! [nz] -15 15 ! !
    nh12=2110+[ncnt]
    profile [nh12] ! [nz] -15 15 ! !
    nh22=2120+[ncnt]
    profile [nh22] ! [nz] -15 15 ! !
    nt/pl //exp/1.ach([ncnt])*sin(theta(1))%z0(1)+[rm]/tan(theta(1)) $71.and.[x2s]<[phir1]<[sm] idh=[nh12]
    nt/pl //exp/1.ach([ncnt])*sin(theta(2))%z0(2)+[rm]/tan(theta(2)) $72.and.[x2s]<[phir2]<[sm] idh=[nh22]
    hi/op/add [nh12] [nh22] [nh02]    
    nh03=3100+[ncnt]
    profile [nh03] ! [nz] -15 15 ! !
    nh13=3110+[ncnt]
    profile [nh13] ! [nz] -15 15 ! !
    nh23=3120+[ncnt]
    profile [nh23] ! [nz] -15 15 ! !
    nt/pl //exp/1.ach([ncnt])*sin(theta(1))%z0(1)+[rm]/tan(theta(1)) $71.and.[sm]<[phir1]<[x2] idh=[nh13]
    nt/pl //exp/1.ach([ncnt])*sin(theta(2))%z0(2)+[rm]/tan(theta(2)) $72.and.[sm]<[phir2]<[x2] idh=[nh23]
    hi/op/add [nh13] [nh23] [nh03]    
    if ([ncnt].eq.1) then
      ng1x=1210+[ncnt]
      profile [ng1x] ! [nf] $sigma(([ncnt]-1)*40-20) $sigma([ncnt]*40+20) ! !
      ng2x=1220+[ncnt]
      profile [ng2x] ! [nf] $sigma(([ncnt]-1)*40-20) $sigma([ncnt]*40+20) ! !
      nt/pl //exp/1.ach([ncnt])*sin(theta(1))%([phir1]-360) $71.and.abs(z0(1)+[rm]/tan(theta(1)))<10 idh=[ng1x]
      nt/pl //exp/1.ach([ncnt])*sin(theta(2))%([phir2]-360) $72.and.abs(z0(2)+[rm]/tan(theta(2)))<10 idh=[ng2x]
      hi/op/add [ng1] [ng1x] [ng1]
      hi/op/add [ng2] [ng2x] [ng2]
    endif
    if ([ncnt].eq.9) then
      ng1x=1210+[ncnt]
      profile [ng1x] ! [nf] $sigma(([ncnt]-1)*40-20) $sigma([ncnt]*40+20) ! !
      ng2x=1220+[ncnt]
      profile [ng2x] ! [nf] $sigma(([ncnt]-1)*40-20) $sigma([ncnt]*40+20) ! !
      nt/pl //exp/1.ach([ncnt])*sin(theta(1))%([phir1]+360) $71.and.abs(z0(1)+[rm]/tan(theta(1)))<10 idh=[ng1x]
      nt/pl //exp/1.ach([ncnt])*sin(theta(2))%([phir2]+360) $72.and.abs(z0(2)+[rm]/tan(theta(2)))<10 idh=[ng2x]
      hi/op/add [ng1] [ng1x] [ng1]
      hi/op/add [ng2] [ng2x] [ng2]
    endif
    hi/op/add [nh1] [nh2] [nh0]
    hi/op/add [ng1] [ng2] [ng0]
    hrout [nh0]
    hrout [nh1]
    hrout [nh2]
    hrout [ng0]
    hrout [ng1]
    hrout [ng2]
    hrout [nh01]
    hrout [nh11]
    hrout [nh21]
    hrout [nh02]
    hrout [nh12]
    hrout [nh22]
    hrout [nh03]
    hrout [nh13]
    hrout [nh23]
  enddo
  close 20
enddo
return

macro mapprplzx n1=1 n2=9
do i=[n1],[n2]
  nh=1100+[i]
  exec ../MinuitTest/mapcal#mapprplz [nh]
enddo
do i=[n1],[n2]
  nh=2100+[i]
  exec ../MinuitTest/mapcal#mapprplz [nh]
enddo
do i=[n1],[n2]
  nh=3100+[i]
  exec ../MinuitTest/mapcal#mapprplz [nh]
enddo
return

macro mapprplz nh=[nh]
exec ../MinuitTest/mapcal#mapprpl [nh] ver=v0
hi/copy 10000 20000
exec ../MinuitTest/mapcal#mapprpl [nh] ver=v1
set ksiz 0.1
null -10 10 0 20
set mtyp 20
set pmci 1
hi/pl 10000 s
set mtyp 20
set pmci 2
hi/pl 20000 s
set mtyp 20
set pmci 1
atitle 'z?r!, cm' 'A, pe'
exec save mhad2011_amplitude_vs_zr_counter[nh]_gaps.eps f
return

macro mapprplx n1=1 n2=9
do i=[n1],[n2]
  nh=1100+[i]
  exec ../MinuitTest/mapcal#mapprpl [nh]
  exec ../MinuitTest/mapcal#mapprpl [nh]
enddo
do i=[n1],[n2]
  nh=2100+[i]
  exec ../MinuitTest/mapcal#mapprpl [nh]
  exec ../MinuitTest/mapcal#mapprpl [nh]
enddo
do i=[n1],[n2]
  nh=3100+[i]
  exec ../MinuitTest/mapcal#mapprpl [nh]
  exec ../MinuitTest/mapcal#mapprpl [nh]
enddo
return

macro mapprply n1=1 n2=41 i1=1
do j=[n1],[n2]
  do i=[i1],9
    nh=200+[i]
    exec ../MinuitTest/mapcal#mapprpl [nh] v2 [j] [j]
  enddo
  i1=1
enddo
return

macro mapprplya
do i=1,9
  nh=200+[i]
  exec ../MinuitTest/mapcal#mapprpl [nh] v0 1 41
enddo
return

macro mapprplypl ncnt=1 ng=1 ver=v0
np=41
ve/cre  s1([np]) r
ve/cre ds1([np]) r
ve/cre  ss([np]) r
ve/cre dss([np]) r
ve/cre  s2([np]) r
ve/cre ds2([np]) r
ve/cre  sig([np]) r
ve/cre dsig([np]) r
ve/cre  sigs([np]) r
ve/cre dsigs([np]) r
ve/cre chi2([np]) r
ve/cre ndf([np]) r
ve/cre  g1([np]) r
ve/cre dg1([np]) r
ve/cre  g2([np]) r
ve/cre dg2([np]) r
ve/cre  g3([np]) r
ve/cre dg3([np]) r
ve/cre ag1([np]) r 
ve/cre ag2([np]) r
ve/cre ag3([np]) r
ve/cre dag([np]) r
ve/cre ord([np]) r 39 40 41 2 38 3 37 4 36 5 35 6 34 7 33 8 32 9 29 10 30 1 31 11 28 12 27 13 14 26 15 25 16 24 17 23 22 18 19 20 21
do i=1,[np]
  fname=mhad2011_amplitude_vs_phi_counter[ncnt]_p[i]_[ver]
  ve/del p16,dp16
  ve/read p16,dp16 [fname]_fit.par '2f15.10'
  ve/inp s1([i]) p16(12)
  ve/inp ds1([i]) dp16(12)
  ve/inp ss([i]) p16(13)
  ve/inp dss([i]) dp16(13)
  ve/inp s2([i]) p16(15)
  ve/inp ds2([i]) dp16(15)
  ve/inp sig([i]) p16(17)
  ve/inp dsig([i]) dp16(17)
  ve/inp sigs([i]) p16(18)
  ve/inp dsigs([i]) dp16(18)
  ve/inp g1([i]) p16(20)
  ve/inp dg1([i]) dp16(20)
  ve/inp g2([i]) p16(21)
  ve/inp dg2([i]) dp16(21)
  ve/inp g3([i]) p16(22)
  ve/inp dg3([i]) dp16(22)
  ve/inp chi2([i]) p16(26)
  ve/inp ndf([i]) dp16(26)
*
  s1=p16(12)
  ss=p16(13)
  sm=p16(14)
  s2=p16(15)
  i1=$sigma(int([s1]/15))
  i2=$sigma(int([s2]/15))
  ni=[i2]-[i1]
  i1=[i1]+1
  f1=[i1]*15
  f2=[i2]*15
  ve/del si0,si,mod,phi,xi
  ve/cre si0([ni]) r
  sigma si0 = array([ni],[f1]#[f2])
  ve/cre si(3) r 3*-1000
  ve/copy si0(1:[ni]) si(1:[ni])
  ve/cre mod(4) r 4*1
  ve/copy p16 pars
  ve/del p16
  ve/copy pars(1:25) p16
  ve/cre xi(10) r [s1] $sigma([s1]+([ss]-[s1])/3) $sigma([s1]+([ss]-[s1])*2/3) _
                  [ss] $sigma([ss]+([sm]-[ss])/3) $sigma([ss]+([sm]-[ss])*2/3) _
                  [sm] $sigma([sm]+([s2]-[sm])/3) $sigma([sm]+([s2]-[sm])*2/3) [s2]
  do j=1,[ni]
    k=22+[j]
    ve/del phi
    ve/cre phi(1) r $sigma(si([j])+p16([k]))
    a0=$call('scphip.f(phi)')
    k=19+[j]
    ve/inp p16([k]) 0
    a1=$call('scphip.f(phi)')
    ve/inp ag[j]([i]) $sigma([a1]-[a0])
  enddo
enddo
*s1=p16(12)
*ss=p16(13)
*sm=p16(14)
*s2=p16(15)
*sig=p16(17)
*sigs=p16(18)
ve/cre  vnp([np]) r
ve/cre vdnp([np]) r
sigma vnp = array([np],1#[np])
ve/cre  qs1([np]) r
ve/cre  qss([np]) r
ve/cre  qs2([np]) r
sigma qs1 = s1-vsum(s1)/[np]-1
sigma qss = ss-vsum(ss)/[np]+0
sigma qs2 = s2-vsum(s2)/[np]+1
*sigma qs1 = order(qs1,ord)
*sigma qss = order(qss,ord)
*sigma qs2 = order(qs2,ord)
*sigma ds1 = order(ds1,ord)
*sigma dss = order(dss,ord)
*sigma ds2 = order(ds2,ord)
*sigma vnp = order(vnp,ord)
null -1 43 -3 3
sigma ds1 = ds1*sqrt(chi2)
sigma dss = dss*sqrt(chi2)
sigma ds2 = ds2*sqrt(chi2)
* exec $PER/s#vpl qs2 ds2 ord vdnp iatt=22 o=s
* exec $PER/s#vpl qss dss ord vdnp iatt=24 o=s
* exec $PER/s#vpl qs1 ds1 ord vdnp iatt=20 o=s
*read x
sigma dsig = dsig*sqrt(chi2)
sigma dsigs = dsigs*sqrt(chi2)
* exec $PER/s#vpl  sig  dsig ord vdnp iatt=20
* exec $PER/s#vpl sigs dsigs ord vdnp iatt=24 o=s
*read x
sigma g1 = abs(g1)
sigma dg1 = dg1*sqrt(chi2)
sigma g2 = abs(g2)
sigma dg2 = dg2*sqrt(chi2)
sigma g3 = abs(g3)
sigma dg3 = dg3*sqrt(chi2)
*fun/pl scphip.f -20 40 
*
if ($sigma(si([ng])).ne.-1000) then
  * exec $PER/s#vpl g[ng] dg[ng] ord vdnp iatt=20 ll=-1
  line -10 0 100 0
  atitle 'energy point number' 'gap, degree'
  gapn=$sigma(si([ng])/15)
  txt=mhad2011: gap [gapn]
  exec $PER/s#tf 0.6 0.9 [txt]
  exec save mhad2011_gap[gapn]_vs_time_[ver].eps f
*  read x
  * exec $PER/s#vpl ag[ng] dag ord vdnp iatt=20
  atitle 'energy point number' 'gap amplitude, degree'
  exec save mhad2011_gap[gapn]_amplitude_vs_time_[ver].eps f
endif
** exec $PER/s#vpl g2 dg2 ord vdnp iatt=24 o=s
** exec $PER/s#vpl g3 dg3 ord vdnp iatt=22 o=s
return


macro phigaps ver=v0
do i=1,9
  do j=1,3
    exec mapcal#mapprplypl [i] [j] ver=[ver]
  enddo
enddo
return



macro phiposic ver=v0
do i=1,-9
  exec mapcal#phipos [i] ver=[ver]
  gl/imp phii
  gl/imp dphii
  gl/imp phic
  if ([i].eq.1) then
    ve/cre phic(9) r
    ve/cre phir(9) r
    ve/cre dphir(9) r
  endif
  ve/inp phic([i])  [phic]
  ve/inp phir([i])  [phii]
  ve/inp dphir([i]) [dphii]
enddo
ve/cre ni(9) r
sigma ni = array(9,1#9)
ve/cre dni(9) r
null 0 10 0 3
set hcol 1
set plci 1
* exec $PER/s#vpl phir dphir ni dni o=s
* exec $PER/s#vpl phic dni ni dni o=s iatt=24
atitle 'walls' 'walls width, degree'
sphic=$sigma(vsum(phic))
sphir=$sigma(vsum(phir))
dsphir=$sigma(sqrt(vsum(dphir**2)))
txt=[Sf]?i! = ( $sigma(int([sphir]*100+0.5)/100) [\261] $sigma(int([dsphir]*100+0.5)/100) )[\260]
exec $PER/s#tf 0.1 0.9 [txt]
txt=[Sf]?i! =   $sigma(int([sphic]*100+0.5)/100)[\260]
exec $PER/s#tf 0.1 0.8 [txt]
exec save mhad2011_walls_hist_[ver].eps f
return


macro phiposi ver=v0
do i=1,9
  exec mapcal#phipos [i] ver=[ver]
  gl/imp ss0
  gl/imp dss0
  if ([i].eq.1) then
    ve/cre phi(9) r
    ve/cre dphi(9) r
  endif
  ve/inp phi([i]) [ss0]
  ve/inp dphi([i]) [dss0]
enddo
ve/cre ni(9) r
sigma ni = array(9,1#9)
ve/cre dni(9) r
null 0 10 -30 350
set hcol 1
set plci 1
* exec $PER/s#vpl phi dphi ni dni o=s
npar=2
ve/cre chi2(2) r
ve/cre paru([npar]) r
ve/cre dparu([npar]) r
ve/cre covu([npar],[npar]) r
ve/cre p2(2) r -40 40
ve/cre s2(2) r -1 0
ve/fit ni phi dphi p1 sb 2 p2 s2
call covm.f(1)
call covmpen(chi2,[npar],paru,dparu)
call covmcov([npar],covu)
atitle 'counter' '[f]?wall!, degree'
exec save mhad2011_phi_walls_[ver].eps f

b=p2(1)
k=p2(2)
sigma wphi = phi-(([b])+([k])*ni)
* exec $PER/s#vpl wphi dphi ni dni ll=-1
atitle 'counter' 'd[f]?wall!, degree'
exec save mhad2011_dphi_walls_[ver].eps f
return


macro phipos iw=1 ver=v0
do i=1,-9
  exec mapcal#mapprpl $sigma(200+[i]) ver=[ver] nfit=1
  hi/copy 10000 10000[i]
  exec mapcal#mapprplypl [i] 1 ver=[ver]
  sigma  s1 = order(s1,ord)
  sigma ds1 = order(ds1,ord)
  sigma  s2 = order(s2,ord)
  sigma ds2 = order(ds2,ord)
  sigma  ss = order(ss,ord)
  sigma dss = order(dss,ord)
  ve/copy  s1  sl[i]
  ve/copy ds1 dsl[i]
  ve/copy  s2  sr[i]
  ve/copy ds2 dsr[i]
  ve/copy  ss  sm[i]
  ve/copy dss dsm[i]
  ve/del par[i],dpar[i]
  ve/read par[i],dpar[i] mhad2011_amplitude_vs_phi_counter[i]__p1_p41_[ver]_fit.par '2f15.10'
enddo
*iw=1
ic1=[iw]-1
df=0
if ([ic1].eq.0) then
  ic1=9
  df=360
endif
sigma s1 = sr[ic1]-([df])
sigma ds1 = dsr[ic1]
ic2=[iw]
sigma s2 = sl[ic2]
sigma ds2 = dsl[ic2]
np=$vlen(s1)
lvl=$sigma(vsum((s1+s2))/[np]/2)

sigma ds12 = s2-s1
sigma dds12 = sqrt(ds2**2+ds1**2)
npar=1
ve/cre chi2(2) r
ve/cre paru([npar]) r
ve/cre dparu([npar]) r
ve/cre covu([npar],[npar]) r
ve/cre g

sigma ss12 = (s2+s1)/2
sigma dss12 = sqrt(ds2**2+ds1**2)/2
ve/fit vnp ss12 dss12 p0 ! 1 g
call covm.f(1)
call covmpen(chi2,[npar],paru,dparu)
call covmcov([npar],covu)
gl/cre ss0 $sigma(paru(1))
gl/cre dss0 $sigma(dparu(1)*sqrt(chi2(1)))

ve/fit vnp ds12 dds12 p0 ! 1 g
call covm.f(1)
call covmpen(chi2,[npar],paru,dparu)
call covmcov([npar],covu)
sigma dparu = dparu*sqrt(vsum(chi2))
phi=$sigma(int(paru(1)*100+0.5)/100)
dphi=$sigma(int(dparu(1)*100+0.5)/100)
text=[d f]?[iw]! = ( [phi] "A# [dphi] )^o!
*"

gl/cre  phii  [phi]
gl/cre dphii [dphi]

null 0 $sigma([np]+1) $sigma([lvl]-5) $sigma([lvl]+5)
set hcol 1
set pmci 1
* exec $PER/s#vpl s1 ds1 vnp vdnp iatt=20 o=s
* exec $PER/s#vpl s2 ds2 vnp vdnp iatt=24 o=s
atitle 'energy point' '[f], degree'
exec $PER/s#tf 0.1 0.8 [text]
exec save mhad2011_wall[iw]_vs_time_[ver].eps f
*read x
*
  ve/del p16
  l=[lvl]-20
  r=[lvl]+20
  set ksiz 0.1
  set mtyp 20
*
  null [l] [r] 0 20
*  null [l] [r] 0 20 sab
  pi=3.1415926535
  if ($sigma(mod([iw]-1,3)).eq.0) then
    dw=$sigma(1/123*180/[pi])
    dt=$sigma(0.66/123*180/[pi])
    s0=$sigma((par[ic2](12)+par[ic1](15)-[df])/2)
    set bord 1
    set fais 1
    set plci 1
    set faci 5
    bl=$sigma([s0])
    br=$sigma([s0]+[dw])
    box [bl] [br] 0 15
    bl=$sigma([s0]-[dw])
    br=$sigma([s0])
    box [bl] [br] 0 15
    set fais 0
    bl=$sigma([s0]+[dw])
    br=$sigma([s0]+[dw]+[dt])
    box [bl] [br] 0 15
    bl=$sigma([s0]-[dw]-[dt])
    br=$sigma([s0]-[dw])
    box [bl] [br] 0 15
    dwc=$sigma(int(200*([dw]+[dt])+0.5)/100)
  else
    dw=$sigma(0.5/123*180/[pi])
    dt=$sigma(0.66/123*180/[pi])
    s0=$sigma((par[ic2](12)+par[ic1](15)-[df])/2)
    set bord 1
    set fais 1
    set plci 1
    set faci 5
    bl=$sigma([s0]-[dw]/2)
    br=$sigma([s0]+[dw]/2)
    box [bl] [br] 0 15
    set fais 0
    bl=$sigma([s0]+[dw]/2)
    br=$sigma([s0]+[dw]/2+[dt])
    box [bl] [br] 0 15
    bl=$sigma([s0]-[dw]/2-[dt])
    br=$sigma([s0]-[dw]/2)
    box [bl] [br] 0 15
    dwc=$sigma(int(200*([dw]/2+[dt])+0.5)/100)
  endif
  gl/cre phic [dwc]
*
  ve/copy par[ic2] p16
  s1=p16(12)
  ss=p16(13)
  sm=p16(14)
  s2=p16(15)
  i1=$sigma(int([s1]/15))
  i2=$sigma(int([s2]/15))
  ni=[i2]-[i1]
  i1=[i1]+1
  f1=[i1]*15
  f2=[i2]*15
  mess [s1] [ss] [sm] [s2] [i1] [i2]
  ve/del si0,si,mod,xi
  ve/cre si0([ni]) r
  sigma si0 = array([ni],[f1]#[f2])
  ve/cre si(3) r 3*-1000
  ve/copy si0(1:[ni]) si(1:[ni])
  ve/cre mod(4) r 4*1
  ve/cre xi(10) r [s1] $sigma([s1]+([ss]-[s1])/3) $sigma([s1]+([ss]-[s1])*2/3) _
                  [ss] $sigma([ss]+([sm]-[ss])/3) $sigma([ss]+([sm]-[ss])*2/3) _
                  [sm] $sigma([sm]+([s2]-[sm])/3) $sigma([sm]+([s2]-[sm])*2/3) [s2]
  set hcol 1
  hi/pl 10000[ic2] s
  set hcol 2                
  fun/pl scphip.f $sigma(max([s1]-5,[l])) $sigma(min([s2]+5,[r])) s
*
  ve/copy par[ic1] p16
*  ve/inp p16(12) $sigma(p16(12)-360)
*  ve/inp p16(13) $sigma(p16(13)-360)
*  ve/inp p16(14) $sigma(p16(14)-360)
*  ve/inp p16(15) $sigma(p16(15)-360)
  s1=p16(12)
  ss=p16(13)
  sm=p16(14)
  s2=p16(15)
  i1=$sigma(int([s1]/15))
  i2=$sigma(int([s2]/15))
  ni=[i2]-[i1]
  i1=[i1]+1
  f1=[i1]*15
  f2=[i2]*15
  mess [s1] [ss] [sm] [s2] [i1] [i2]
  ve/del si0,si,mod,xi
  ve/cre si0([ni]) r
  sigma si0 = array([ni],[f1]#[f2])
  ve/cre si(3) r 3*-1000
  ve/copy si0(1:[ni]) si(1:[ni])
  ve/cre mod(4) r 4*1
  ve/cre xi(10) r [s1] $sigma([s1]+([ss]-[s1])/3) $sigma([s1]+([ss]-[s1])*2/3) _
                  [ss] $sigma([ss]+([sm]-[ss])/3) $sigma([ss]+([sm]-[ss])*2/3) _
                  [sm] $sigma([sm]+([s2]-[sm])/3) $sigma([sm]+([s2]-[sm])*2/3) [s2]
  null $sigma([l]+[df]) $sigma([r]+[df]) 0 20 sab
  set hcol 1
  hi/pl 10000[ic1] s
  set hcol 2
  fun/pl scphip.f $sigma(max([s1]-5,[l]+[df])) $sigma(min([s2]+5,[r]+[df])) s
*
  exec $PER/s#tf 0.1 0.9 [text]
  txt=[d f]?[iw]! = [dwc] ^o! (calc.)
  exec $PER/s#tf 0.1 0.8 [txt]
  atitle '[f], dergee' 'A, pe.'
  exec save mhad2011_wall[iw]_[ver].eps f
return


macro mapprpl nh=201 ver=v0 n10=! n20=! nfit=0
n=0
n=[n]+1;np[n]=39; ebeam[n]=508.7 ; ebeam[n]=508.7 ; il[n]=0
n1=[n]
n=[n]+1;np[n]=40; ebeam[n]=509.8 ; ebeam[n]=509.8 ; il[n]=0
n=[n]+1;np[n]=41; ebeam[n]=510.5 ; ebeam[n]=510.5 ; il[n]=0
n=[n]+1;np[n]=2;  ebeam[n]=525   ; ebeam[n]=525   ; il[n]=363.2
n=[n]+1;np[n]=38; ebeam[n]=537.5 ; ebeam[n]=537.5 ; il[n]=546.1
n=[n]+1;np[n]=3;  ebeam[n]=550   ; ebeam[n]=550   ; il[n]=485.7
n=[n]+1;np[n]=37; ebeam[n]=562.5 ; ebeam[n]=562.5 ; il[n]=517.3
n=[n]+1;np[n]=4;  ebeam[n]=575   ; ebeam[n]=575   ; il[n]=419.1
n=[n]+1;np[n]=36; ebeam[n]=587.5 ; ebeam[n]=587.5 ; il[n]=534.4
n=[n]+1;np[n]=5;  ebeam[n]=600   ; ebeam[n]=600   ; il[n]=489.8999940
n=[n]+1;np[n]=35; ebeam[n]=612.5 ; ebeam[n]=612.5 ; il[n]=546.7000120
n=[n]+1;np[n]=6;  ebeam[n]=625   ; ebeam[n]=625   ; il[n]=440.5000000
n=[n]+1;np[n]=34; ebeam[n]=637.5 ; ebeam[n]=637.5 ; il[n]=496.2000120
n=[n]+1;np[n]=7;  ebeam[n]=650   ; ebeam[n]=650   ; il[n]=456.1000060
n=[n]+1;np[n]=33; ebeam[n]=662.5 ; ebeam[n]=662.5 ; il[n]=526.2999880
n=[n]+1;np[n]=8;  ebeam[n]=675   ; ebeam[n]=675   ; il[n]=559.7999880
n=[n]+1;np[n]=32; ebeam[n]=687.5 ; ebeam[n]=687.5 ; il[n]=576.2000120
n=[n]+1;np[n]=9;  ebeam[n]=700   ; ebeam[n]=700   ; il[n]=577.0999760
n=[n]+1;np[n]=29; ebeam[n]=712.5 ; ebeam[n]=712.5 ; il[n]=593.2000120
n=[n]+1;np[n]=10; ebeam[n]=725   ; ebeam[n]=725   ; il[n]=432.2999880
n=[n]+1;np[n]=30; ebeam[n]=737.5 ; ebeam[n]=737.5 ; il[n]=600.0999760
n=[n]+1;np[n]=1;  ebeam[n]=750   ; ebeam[n]=750   ; il[n]=694.7000120
n=[n]+1;np[n]=31; ebeam[n]=762.5 ; ebeam[n]=762.5 ; il[n]=487.2999880
n=[n]+1;np[n]=11; ebeam[n]=775   ; ebeam[n]=775   ; il[n]=546.5000000
n=[n]+1;np[n]=28; ebeam[n]=787.5 ; ebeam[n]=787.5 ; il[n]=502.3999940
n=[n]+1;np[n]=12; ebeam[n]=800   ; ebeam[n]=800   ; il[n]=439.3999940
n=[n]+1;np[n]=27; ebeam[n]=812.5 ; ebeam[n]=812.5 ; il[n]=509.8999940
n=[n]+1;np[n]=13; ebeam[n]=825   ; ebeam[n]=825   ; il[n]=475.6000060
n=[n]+1;np[n]=14; ebeam[n]=850   ; ebeam[n]=850   ; il[n]=462.2000120
n=[n]+1;np[n]=26; ebeam[n]=862.5 ; ebeam[n]=862.5 ; il[n]=504.6000060
n=[n]+1;np[n]=15; ebeam[n]=875   ; ebeam[n]=875   ; il[n]=506.3999940
n=[n]+1;np[n]=25; ebeam[n]=887.5 ; ebeam[n]=887.5 ; il[n]=525.7000120
n=[n]+1;np[n]=16; ebeam[n]=900   ; ebeam[n]=900   ; il[n]=384.1000060
n=[n]+1;np[n]=24; ebeam[n]=912.5 ; ebeam[n]=912.5 ; il[n]=481.3999940
n=[n]+1;np[n]=17; ebeam[n]=925   ; ebeam[n]=925   ; il[n]=401.2999880
n=[n]+1;np[n]=23; ebeam[n]=935   ; ebeam[n]=935   ; il[n]=637.5000000
n=[n]+1;np[n]=22; ebeam[n]=945   ; ebeam[n]=945   ; il[n]=577.0999760
n=[n]+1;np[n]=18; ebeam[n]=950   ; ebeam[n]=950   ; il[n]=451.2999880
n=[n]+1;np[n]=19; ebeam[n]=962.5 ; ebeam[n]=962.5 ; il[n]=561.9000240
n=[n]+1;np[n]=20; ebeam[n]=987.5 ; ebeam[n]=987.5 ; il[n]=467.5000000
n=[n]+1;np[n]=21; ebeam[n]=1000  ; ebeam[n]=1000  ; il[n]=538.0999760
n2=[n]

if ([n10].ne.'!') then
  n1=[n10]
endif
if ([n20].ne.'!') then
  n2=[n20]
endif

dir=/work/users/konctbel/MinuitTest

ci=0
hi/del [nh]
do i=[n1],[n2]
  ebeam=[ebeam[i]]
  if ([ver].eq.'v0') then
    hfile=hist_[ebeam]_prof.his
  endif
  if ([ver].eq.'v1') then
    hfile=hist_[ebeam]_prof_v1.his
  endif
  if ([ver].eq.'v2') then
    hfile=hist_[ebeam]_prof_v2.his
  endif
  hi/file 20 [dir]/[hfile]
  ci=$sigma(mod([i],5)+1)
  set pmci [ci]
  if ([i].eq.[n1]) then
    hi/pl //lun20/[nh]
    hi/copy [nh] 10000
  else
    hi/pl //lun20/[nh] s
    hi/op/add [nh] 10000 10000
  endif
  close 20
enddo

set pmci 1
*read x
hi/pl 10000
if (([nh].lt.200).or.([nh].gt.1000)) then
*ve/cre p30(31) r -11.629 -10.100 -8.4205 -5.9000 -4.4165 -2.7000 -1.6125 0 1.4896 3.5677 4.3000 6.5673 7.5000 10.100 10.700 _
ve/cre p30(34) r -10.295 -9 -7.8826 -6.1744 -4.6141 -3.0539 -1.5129 0 1.5366 2.2765 4.2113 5.5521 7.8676 9.3502 10.422 _
5.5740 6.5977 7.5585 6.8670 7.9323 7.4458 7.7306 8.0075 8.2367 8.2778 9.0746 8.6534 8.5318 4.0953 0 1 1 11 0.3
ve/cre p30(34) r -11.681 -10.3075 -8.6210 -7.0100 -5.3990 -3.3885 -1.7770 -0.39974 1.0450 2.6560 4.2670 6.2445 7.8890 9.1000 10.321 _
7.5298 7.9029 9.2221 8.9689 9.3056 9.4394 9.2645 9.8690 10.270 10.567 11.350 11.154 10.873 9.9358 0.0000 0.66154 4.1908 10.816 0.63145
if ([ver].eq.'v0') then
  pfile=mhad2011_h[nh]_fit.pars
endif
if ([ver].eq.'v1') then
  pfile=mhad2011_h[nh]_fit_v1.pars
endif
if ($fexist([pfile])) then
  ve/del p30,s30,pmin,pmax
  ve/read p30,s30,pmin,pmax [pfile] '4f15.10'
endif
ve/inp p30(30) 0.0
ve/cre xi(15) r
sigma xi = array(14,-11.643#9.3000)
ve/inp xi(15) 10.590
*ve/copy xi(1:15) p30(1:15)
ve/cre s30(34) r 15*0.01 14*0.1 0 0.001 0.1 0.01 0.01
ve/cre pmin(34) r 30*0 -2 0 10.5 0.2
ve/cre pmax(34) r 30*20 5 10 12 1
ve/cre xi(15) r
ve/copy p30(1:15) xi(1:15)
ve/cre ximin(15) r
sigma ximin = xi-0.3
ve/cre ximax(15) r
sigma ximax = xi+0.3
ve/copy ximin(1:15) pmin(1:15)
ve/copy ximax(1:15) pmax(1:15)
do i=1,10
*  hi/fit 10000(-13.:12.) /work/users/konctbel/SepPar/uzfit.f b 34 p30 s30 pmin pmax
enddo
ve/inp p30(30) 2.0
*hi/fit 10000(-13.:13.) /work/users/konctbel/SepPar/uzfit.f b 34 p30 s30 pmin pmax
do i=1,15
  line $sigma(p30([i])) 0 $sigma(p30([i])) 20
enddo
*shell cp -v [pfile] [pfile].old
*ve/write p30,s30,pmin,pmax [pfile] '4f15.10'

*************************************************************

*************************************************************

else 

if ([nfit].eq.0) then

ncnt=[nh]-200
exec mapcal#readlif [ncnt]
point 1000

if ([ver].ne.'v2') then

hi/copy 10000 20000
if ([ncnt].eq.1) then
  exec hsigma @10001 = @10000*($10000 gt 14 and $10000 lt 16)
  exec hsigma @10001 = @10001/@10001
  exec hsigma @10001 = mod(@10001+1,2)
  hi/op/mult 20000 10001 20000
endif
if ([ncnt].eq.5) then
  exec hsigma @10001 = @10000*($10000 gt 179.4 and $10000 lt 180.6)
  exec hsigma @10002 = @10000*($10000 gt 194 and $10000 lt 196)
  exec hsigma @10001 = @10001/@10001+@10002/@10002
  exec hsigma @10001 = mod(@10001+1,2)
  hi/op/mult 20000 10001 20000
endif

ve/cre mod(2) r 1 1
s1=pars0(4)
s2=pars0(5)
ss=pars0(14)
sm=$sigma(([s2]+[ss])/2)
sig=pars0(7)
ve/cre xi(10) r [s1] $sigma([s1]+([ss]-[s1])/3) $sigma([s1]+([ss]-[s1])*2/3) _
                [ss] $sigma([ss]+([sm]-[ss])/3) $sigma([ss]+([sm]-[ss])*2/3) _
                [sm] $sigma([sm]+([s2]-[sm])/3) $sigma([sm]+([s2]-[sm])*2/3) [s2]
exec hsigma ym = vmax(@10000)                
*
i1=$sigma(int([s1]/15))
i2=$sigma(int([s2]/15))
ni=[i2]-[i1]
i1=[i1]+1
f1=[i1]*15
f2=[i2]*15
ve/cre si0([ni]) r
sigma si0 = array([ni],[f1]#[f2])
ve/cre si(3) r 3*-1000
ve/copy si0(1:[ni]) si(1:[ni])
dni=3-[ni]
hi/copy 10000 20000
do i=1,[ni]
  l=si0([i])-2.1
  r=si0([i])+2.1
  exec hsigma @1000[i] = @10000*($10000 gt [l] and $10000 lt [r])
enddo
hi/copy 10000 10004
hi/op/res 10004
exec hsigma @10004 = 0*@10000
do i=1,[ni]
  exec hsigma @10004 = @10004+@1000[i]/@1000[i]
enddo
exec hsigma @10005 = mod(@10004+1,2)
hi/op/mult 20000 10005 20000
*
npar=25
ve/cre chi2(2) r
ve/cre paru([npar]) r
ve/cre dparu([npar]) r
ve/cre covu([npar],[npar]) r
*
ve/cre p16([npar]) r 6 8 9 11 11 9 7 4.7 4.7 4.3 3 [s1] [ss] [sm] [s2] $sigma(ym(1)) 0.3 0.3 1.0 [ni]*0.5
ve/cre dp16([npar]) r
if ([dni].eq.0) then
  ve/cre s16([npar]) r 11*0.5 0.2 0 0 0.2 1*1 0.01 0.01 0.0 [ni]*0.00 [ni]*0
else
  ve/cre s16([npar]) r 11*0.5 0.2 0 0 0.2 1*1 0.01 0.01 0.0 [ni]*0.00 [dni]*0 [ni]*0 [dni]*0 
endif
ve/cre pmin([npar]) r 11*0 $sigma([s1]-2) $sigma([ss]-2) $sigma([sm]-2) $sigma([s2]-2) 0 3*0.2 3*0 3*-2
ve/cre pmax([npar]) r 11*50 $sigma([s1]+2) $sigma([ss]+2) $sigma([sm]+2) $sigma([s2]+2) 50 3*1.0 3*1 3*2
hi/pl 20000 e
do i=1,3
  hi/fit 20000 scphi.f sb [npar] p16 s16 pmin pmax dp16
enddo
ve/inp s16(23) 0.01
ve/inp s16(24) 0.01
ve/inp s16(25) 0.01
do i=1,3
  hi/fit 10000 scphi.f sb [npar] p16 s16 pmin pmax dp16
enddo
ve/inp s16(19) 0.01
ve/inp s16(20) 0.01
ve/inp s16(21) 0.01
ve/inp s16(22) 0.01
do i=1,3
  hi/fit 10000 scphi.f sb [npar] p16 s16 pmin pmax dp16
enddo
ve/inp s16(13) 0.02
ve/inp s16(14) 0.02
do i=1,3
  hi/fit 10000 scphi.f sb [npar] p16 s16 pmin pmax dp16
enddo
*read x
hi/pl 10000
i=0
if ([i].gt.1) then
sg=$sigma(abs(p16(19)))
dg=$sigma(abs(p16(20)))
ag=$sigma([dg]/sqrt(2*3.1415927)/[sg]*exp(-([dg]/[sg])**2/8))
if (([ag].gt.0.1).and.([dg].gt.0.05)) then
  ve/inp s16(23) 0.01
endif
dg=$sigma(abs(p16(21)))
ag=$sigma([dg]/sqrt(2*3.1415927)/[sg]*exp(-([dg]/[sg])**2/8))
if (([ag].gt.0.1).and.([dg].gt.0.05)) then
  ve/inp s16(24) 0.01
endif
dg=$sigma(abs(p16(22)))
ag=$sigma([dg]/sqrt(2*3.1415927)/[sg]*exp(-([dg]/[sg])**2/8))
if (([ag].gt.0.1).and.([dg].gt.0.05)) then
  ve/inp s16(25) 0.01
endif
ve/inp s16(17) 0.01
ve/inp s16(18) 0.01
ve/inp s16(19) 0.01
do i=1,5
  hi/fit 10000 scphi.f sb [npar] p16 s16 pmin pmax dp16
enddo
endif
do i=1,5
  l=$sigma([s1]-1*[sig])
  r=$sigma([s2]+1*[sig])
  s1=p16(12)
  ss=p16(13)
  sm=p16(14)
  s2=p16(15)
  ve/cre xi(10) r [s1] $sigma([s1]+([ss]-[s1])/3) $sigma([s1]+([ss]-[s1])*2/3) _
                  [ss] $sigma([ss]+([sm]-[ss])/3) $sigma([ss]+([sm]-[ss])*2/3) _
                  [sm] $sigma([sm]+([s2]-[sm])/3) $sigma([sm]+([s2]-[sm])*2/3) [s2]
  ve/inp s16(17) 0.0
  ve/inp s16(18) 0.0
  ve/inp s16(19) 0.0
  ve/inp s16(23) 0.0
  ve/inp s16(24) 0.0
  ve/inp s16(25) 0.0
  ve/inp s16(19) 0.0
  ve/inp s16(20) 0.0
  ve/inp s16(21) 0.0
  ve/inp s16(22) 0.0
  hi/fit 10000([l]:[r]) scphi.f sb [npar] p16 s16 pmin pmax dp16
  ve/inp s16(17) 0.01
  ve/inp s16(18) 0.01
  ve/inp s16(19) 0.01
  ve/inp s16(23) 0.01
  ve/inp s16(24) 0.01
  ve/inp s16(25) 0.01
  ve/inp s16(19) 0.01
  ve/inp s16(20) 0.01
  ve/inp s16(21) 0.01
  ve/inp s16(22) 0.01
  hi/fit 10000([l]:[r]) scphi.f sb [npar] p16 s16 pmin pmax dp16
enddo
do i=1,5
  hi/fit 10000([l]:[r]) scphi.f sb [npar] p16 s16 pmin pmax dp16
enddo

else

  npar=25
  fname=mhad2011_amplitude_vs_phi_counter[ncnt]_p[n1]
  fname=mhad2011_amplitude_vs_phi_counter[ncnt]__p1_p41_v2
  ve/read pars,dpars [fname]_fit.par '2f15.10'
  ve/del p16
  ve/copy pars(1:25) p16
  ve/cre dp16([npar]) r
  s1=p16(12)
  ss=p16(13)
  sm=p16(14)
  s2=p16(15)
  ve/cre xi(10) r [s1] $sigma([s1]+([ss]-[s1])/3) $sigma([s1]+([ss]-[s1])*2/3) _
                  [ss] $sigma([ss]+([sm]-[ss])/3) $sigma([ss]+([sm]-[ss])*2/3) _
                  [sm] $sigma([sm]+([s2]-[sm])/3) $sigma([sm]+([s2]-[sm])*2/3) [s2]
  i1=$sigma(int([s1]/15))
  i2=$sigma(int([s2]/15))
  ni=[i2]-[i1]
  i1=[i1]+1
  f1=[i1]*15
  f2=[i2]*15
  ve/cre si0([ni]) r
  sigma si0 = array([ni],[f1]#[f2])
  ve/cre si(3) r 3*-1000
  ve/copy si0(1:[ni]) si(1:[ni])
  dni=3-[ni]
  if ([dni].ne.0) then
    ve/cre s16([npar]) r 11*0.5 4*0.2 1*1 3*0.01 [ni]*0.001 [dni]*0 [ni]*0.01 [dni]*0
  else
    ve/cre s16([npar]) r 11*0.5 4*0.2 1*1 3*0.01 [ni]*0.001 [ni]*0.01
  endif
  ve/cre pmin([npar]) r 11*0 $sigma([s1]-2) $sigma([ss]-2) $sigma([sm]-2) $sigma([s2]-2) 0 3*0.2 3*0 3*-2
  ve/cre pmax([npar]) r 11*50 $sigma([s1]+2) $sigma([ss]+2) $sigma([sm]+2) $sigma([s2]+2) 50 3*1.0 3*1 3*2
  do i=1,3
    sig=p16(17)
    l=$sigma([s1]-1*[sig])
    r=$sigma([s2]+1*[sig])
    s1=p16(12)
    ss=p16(13)
    sm=p16(14)
    s2=p16(15)
    sig=p16(17)
    ve/cre xi(10) r [s1] $sigma([s1]+([ss]-[s1])/3) $sigma([s1]+([ss]-[s1])*2/3) _
                    [ss] $sigma([ss]+([sm]-[ss])/3) $sigma([ss]+([sm]-[ss])*2/3) _
                    [sm] $sigma([sm]+([s2]-[sm])/3) $sigma([sm]+([s2]-[sm])*2/3) [s2]
    ve/inp s16(17) 0.0
    ve/inp s16(18) 0.0
    ve/inp s16(19) 0.0
    ve/inp s16(23) 0.0
    ve/inp s16(24) 0.0
    ve/inp s16(25) 0.0
    ve/inp s16(19) 0.0
    ve/inp s16(20) 0.0
    ve/inp s16(21) 0.0
    ve/inp s16(22) 0.0
    hi/fit 10000([l]:[r]) scphi.f sb [npar] p16 s16 pmin pmax dp16
    ve/inp s16(17) 0.01
    ve/inp s16(18) 0.01
    ve/inp s16(19) 0.01
    ve/inp s16(23) 0.01
    ve/inp s16(24) 0.01
    ve/inp s16(25) 0.01
    ve/inp s16(19) 0.01
    ve/inp s16(20) 0.01
    ve/inp s16(21) 0.01
    ve/inp s16(22) 0.01
    hi/fit 10000([l]:[r]) scphi.f sb [npar] p16 s16 pmin pmax dp16
  enddo
  do i=1,3
    hi/fit 10000([l]:[r]) scphi.f sb [npar] p16 s16 pmin pmax dp16
  enddo
  ve/cre chi2(2) r
  ve/cre paru([npar]) r
  ve/cre dparu([npar]) r
  ve/cre covu([npar],[npar]) r

endif

exec hsigma @10001 = @10000
exec hsigma %10001 = %10000
hi/fit 10000([l]:[r]) scphi.f sb [npar] p16 s16 pmin pmax dp16
call covm.f(1)
call covmpen(chi2,[npar],paru,dparu)
call covmcov([npar],covu)
exec seteps 0
opt ngrid
opt nstat
*set hcol 2
*set pmci 2
set ksiz 0.1
set mtyp 20
set fcol 2
set hwid 1
hi/pl 10000
atitle '[f], degree' 'A, pe'
xbins=$hinfo(10000,'xbins')
xmin=$hinfo(10000,'xmin')
xmax=$hinfo(10000,'xmax')
ndf=$sigma(int([xbins]/([xmax]-[xmin])*([r]-[l]))-[npar])
chi=$sigma(int(chi2(1)*[ndf]*1000+0.5)/1000)
txt=[h]^2!/ndf = [chi] / [ndf]
ve/cre ndf(1) r [ndf]
exec $PER/s#tf 0.05 0.9 [txt]

i=0
if ([i].eq.1) then
  set hcol 1
  set pmci 2
  ve/cre mod(4) r 1 0 0 0
  fun/pl scphip.f $hinfo(10000,'xmin') $hinfo(10000,'xmax') s
  ve/cre mod(4) r 0 1 0 0
  fun/pl scphip.f $hinfo(10000,'xmin') $hinfo(10000,'xmax') s
  ve/cre mod(4) r 0 0 1 0
  fun/pl scphip.f $hinfo(10000,'xmin') $hinfo(10000,'xmax') s
  ve/cre mod(4) r 0 0 0 1
  fun/pl scphip.f $hinfo(10000,'xmin') $hinfo(10000,'xmax') s
*hi/pl 20000 s
  s1=p16(12)
  ss=p16(13)
  sm=p16(14)
  s2=p16(15)
  sig=p16(17)
  sigs=p16(18)
  x1=$sigma([s1]+2*[sig])
  x2=$sigma([s2]-2*[sig])
  line [x1] 0 [x1] 20
  line [x2] 0 [x2] 20
  x1s=$sigma([ss]-1.5/120.0*180/3.1415927-3*[sigs])
  x2s=$sigma([ss]+1.5/120.0*180/3.1415927+3*[sigs])
  line [x1s] 0 [x1s] 20
  line [x2s] 0 [x2s] 20
  mess [x1] [x1s] [x2s] [x2] $sigma(([x2]-[x1])/2) $sigma(([x2s]-[x1s])/2)
endif

if ([n1].ne.[n2]) then
  fname=mhad2011_amplitude_vs_phi_counter[ncnt]__p[n1]_p[n2]_[ver]
else
  fname=mhad2011_amplitude_vs_phi_counter[ncnt]_p[n1]_[ver]
endif
exec save [fname].eps f 
exec vappend p16 chi2
exec vappend dp16 ndf
ve/write p16,dp16 [fname]_fit.par '2f15.10'
endif
endif
return


macro mapprplsx n1=1 n2=9
do i=[n1],[n2]
  nh=1100+[i]
  exec ../MinuitTest/mapcal#mapprpls [nh]
  exec ../MinuitTest/mapcal#mapprpls [nh]
enddo
do i=[n1],[n2]
  nh=2100+[i]
  exec ../MinuitTest/mapcal#mapprpls [nh]
  exec ../MinuitTest/mapcal#mapprpls [nh]
enddo
do i=[n1],[n2]
  nh=3100+[i]
  exec ../MinuitTest/mapcal#mapprpls [nh]
  exec ../MinuitTest/mapcal#mapprpls [nh]
enddo
return

macro mapprpls nh=201
n=0
n=[n]+1; ebeam[n]=508.7 ; il[n]=0
n1=[n]
n=[n]+1; ebeam[n]=509.8 ; il[n]=0
n=[n]+1; ebeam[n]=510.5 ; il[n]=0
n=[n]+1; ebeam[n]=525   ; il[n]=363.2
n=[n]+1; ebeam[n]=537.5 ; il[n]=546.1
n=[n]+1; ebeam[n]=550   ; il[n]=485.7
n=[n]+1; ebeam[n]=562.5 ; il[n]=517.3
n=[n]+1; ebeam[n]=575   ; il[n]=419.1
n=[n]+1; ebeam[n]=587.5 ; il[n]=534.4
n=[n]+1; ebeam[n]=600   ; il[n]=489.8999940
n=[n]+1; ebeam[n]=612.5 ; il[n]=546.7000120
n=[n]+1; ebeam[n]=625   ; il[n]=440.5000000
n=[n]+1; ebeam[n]=637.5 ; il[n]=496.2000120
n=[n]+1; ebeam[n]=650   ; il[n]=456.1000060
n=[n]+1; ebeam[n]=662.5 ; il[n]=526.2999880
n=[n]+1; ebeam[n]=675   ; il[n]=559.7999880
n=[n]+1; ebeam[n]=687.5 ; il[n]=576.2000120
n=[n]+1; ebeam[n]=700   ; il[n]=577.0999760
n=[n]+1; ebeam[n]=712.5 ; il[n]=593.2000120
n=[n]+1; ebeam[n]=725   ; il[n]=432.2999880
n=[n]+1; ebeam[n]=737.5 ; il[n]=600.0999760
n=[n]+1; ebeam[n]=750   ; il[n]=694.7000120
n=[n]+1; ebeam[n]=762.5 ; il[n]=487.2999880
n=[n]+1; ebeam[n]=775   ; il[n]=546.5000000
n=[n]+1; ebeam[n]=787.5 ; il[n]=502.3999940
n=[n]+1; ebeam[n]=800   ; il[n]=439.3999940
n=[n]+1; ebeam[n]=812.5 ; il[n]=509.8999940
n=[n]+1; ebeam[n]=825   ; il[n]=475.6000060
n=[n]+1; ebeam[n]=850   ; il[n]=462.2000120
n=[n]+1; ebeam[n]=862.5 ; il[n]=504.6000060
n=[n]+1; ebeam[n]=875   ; il[n]=506.3999940
n=[n]+1; ebeam[n]=887.5 ; il[n]=525.7000120
n=[n]+1; ebeam[n]=900   ; il[n]=384.1000060
n=[n]+1; ebeam[n]=912.5 ; il[n]=481.3999940
n=[n]+1; ebeam[n]=925   ; il[n]=401.2999880
n=[n]+1; ebeam[n]=935   ; il[n]=637.5000000
n=[n]+1; ebeam[n]=945   ; il[n]=577.0999760
n=[n]+1; ebeam[n]=950   ; il[n]=451.2999880
n=[n]+1; ebeam[n]=962.5 ; il[n]=561.9000240
n=[n]+1; ebeam[n]=987.5 ; il[n]=467.5000000
n=[n]+1; ebeam[n]=1000  ; il[n]=538.0999760
n2=[n]

dir=/work/users/konctbel/MinuitTest

ci=0
hi/del [nh]
do i=[n1],[n2]
  ebeam=[ebeam[i]]
  hfile=hist_[ebeam]_prof.his
  hi/file 20 [dir]/[hfile]
  ci=$sigma(mod([i],5)+1)
  set pmci [ci]
  if ([i].eq.[n1]) then
    hi/pl //lun20/[nh]
    hi/copy [nh] 10000
  else
    hi/pl //lun20/[nh] s
    hi/op/add [nh] 10000 10000
  endif
  close 20
enddo

set pmci 1
set fcol 2
set fwid 3
*read x
hi/pl 10000
if (([nh].lt.200).or.([nh].gt.1000)) then
pfile=mhad2011_h[nh]_fit.pars
if ($fexist([pfile])) then
  ve/del p30,s30,pmin,pmax
  ve/read p30,s30,pmin,pmax [pfile] '4f15.10'
endif
ve/cre s30(34) r 34*0
hi/fit 10000(-13.:13.) /work/users/konctbel/SepPar/uzfit.f b 34 p30 s30 pmin pmax
do i=1,15
  line $sigma(p30([i])) 0 $sigma(p30([i])) 20
enddo
atitle 'z?r!, cm' 'A, pe'
*exec save mhad2011_amplitude_vs_zr_counter[nh].eps f
else
endif
return


macro mapprfitxx sm=u opt=f
exec mapcal#mapprfitx 1 1 9 [sm] [opt]
exec mapcal#mapprfitx 2 1 9 [sm] [opt]
exec mapcal#mapprfitx 3 1 9 [sm] [opt]
exec mapcal#mapprfitx 4 1 9 [sm] [opt]
return


macro mapprfitx map=1 n1=1 n2=9 sm=0 opt=f
*
if ([map].eq.1) then
  gl/cre mname mapcal1
  gl/cre nh1 11
  gl/cre nh2 11
  gl/cre nb 14
  npt=0
endif
if ([map].eq.2) then
  gl/cre mname mapcal2
  gl/cre nh1 12
  gl/cre nh2 12
  gl/cre nb 14
  npt=1
endif
if ([map].eq.3) then
  gl/cre mname mapcal3
  gl/cre nh1 13
  gl/cre nh2 13
  gl/cre nb 10
  npt=2
endif
if ([map].eq.4) then
  gl/cre mname mapcal4
  gl/cre nh1 14
  gl/cre nh2 14
  gl/cre nb 14
  npt=3
endif
*
if ([opt].eq.'f') then
*
  ks=0
*
  ve/del ptl
  ve/read ptl /work/users/konctbel/MinuitTest/phicuts_test.txt
  ve/cre ic(9) r
  ve/cre if1(9) r
  ve/cre if2(9) r
  ve/cre if3(9) r
  ve/cre if4(9) r
  ve/cre if5(9) r
*
  do i=1,9
    ind=$sigma(128*[npt]+14*([i]-1)+14)
    dxl=$sigma(ptl([ind]))
    do j=1,5
      ind=$sigma(128*[npt]+14*([i]-1)+8+[j])
      xl=$sigma(ptl([ind])-[dxl]*[ks])
      ve/inp if[j]([i]) [xl]
    enddo
  enddo
*
  cmd=.testrelease/.mainrelease/Offline/submit.sh -q clusters,180
*
  do layer=1,3

    do i=[n1],[n2]
      nh=$sigma([layer]100+[i])
	  mess im here
      exec ../MinuitTest/mapcal#mapprfit [nh] [sm]
    enddo
  
    do i=[n1],[n2]
      nh=$sigma([layer]100+[i])
      fname=[mname]_amplitude_vs_zr_counter[nh]_[sm]
      if ($fexist([fname])) then
        shell rm [mname]_amplitude_vs_zr_counter[nh]_[sm]
      endif
      shell rm [mname]/[fname].o*
*      shell nice -19 FIT [fname] >& [fname].log &
*      shell $unquote([cmd]) FIT [fname] >& [fname]
      shell nice -19 FIT [fname] norun
      shell $unquote([cmd]) -o [mname] [fname]
    enddo
  
  enddo
 
endif
*
if ([opt].eq.'c') then

  do layer=1,3

    do i=[n1],[n2]
    
      nh=$sigma([layer]100+[i])
      fname=[mname]_amplitude_vs_zr_counter[nh]_[sm]
      shell mv [fname]{.dat,_exp.his,.for,_thr.his} [mname]/
      shell mv [mname]/[fname].o* [mname]/[fname].log
      
      shell cp [mname]/[fname]_exp.his tmp.txt
      shell $unquote('cat tmp.txt | sed "s/.00 0.05-0./.00 0.05 0./g" > tmp1.txt')
      shell cp tmp1.txt [mname]/[fname]_exp.his      
      
    enddo
  
  enddo

endif
*
*
if ([opt].eq.'p') then

  fname=[mname]/[mname]_z_profiles_[sm].tex
  if ($fexist([fname]).eq.1) then
    shell rm [fname]
  endif
  for/file  20 [fname] new
  close 20
  
  do layer=1,3

    do i=[n1],[n2]
    
      nh=$sigma([layer]100+[i])
      flname=[mname]_amplitude_vs_zr_counter[nh]_[sm]
      epsfile=[flname].eps
      
      fkname=mapprfit[i].kumac
      if ($fexist([fkname]).eq.1) then
        shell rm [fkname]
      endif
      for/file 20 [fkname] new
      close 20
      
      txt=exec crsfit load u [mname]/[flname]
      fmess [txt] [fkname]
      txt=opt nstat
      fmess [txt] [fkname]
      txt=opt ngrid
      fmess [txt] [fkname]
      txt=set plci 1
      fmess [txt] [fkname]
      txt=set lwid 1 
      fmess [txt] [fkname]
      txt=set fcol 2
      fmess [txt] [fkname]
      txt=set fwid 3
      fmess [txt] [fkname]
      txt=exec crsfit u 1 c
      fmess [txt] [fkname]
      txt=atitle 'z?mid!, cm' 'amplitude, pe.'
      fmess [txt] [fkname]
      txt=exec save [mname]/[epsfile] f
      fmess [txt] [fkname]
      shell pawbigX11 -b [fkname]
  
      fmess '\begin{figure}[ht!b]' [fname]
      fmess '  \begin{minipage}{\textwidth}' [fname]
      fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
      txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
      fmess [txt] [fname]
      txt=$unquote('     ')\caption{[mname]: counter \No [i],  layer \No [layer]}
      fmess [txt] [fname]
      fmess '   \end{minipage}' [fname]
      fmess '\end{figure}' [fname]
    enddo
    fmess '\clearpage' [fname]
  
  enddo

  fname=[mname]_latex.sh
  if ($fexist([fname]).eq.1) then
    shell rm [fname]
  endif
  for/file  20 [fname] new
  close 20
  
  txt=cd [mname]
  fmess [txt] [fname]
  txt=latex [mname]_amplitude_vs_zr_[sm]
  fmess [txt] [fname]
  txt=latex [mname]_amplitude_vs_zr_[sm]
  fmess [txt] [fname]
  txt=dvips [mname]_amplitude_vs_zr_[sm]
  fmess [txt] [fname]
  
  shell chmod +x [fname]
  shell ./[fname]
endif
*
return



macro mapprfit nh=201 sm=0 prof=1
n=0
n=[n]+1; ebeam[n]=508.7 ; il[n]=0
n1=[n]
n=[n]+1; ebeam[n]=509.8 ; il[n]=0
n=[n]+1; ebeam[n]=510.5 ; il[n]=0
n=[n]+1; ebeam[n]=525   ; il[n]=363.2
n=[n]+1; ebeam[n]=537.5 ; il[n]=546.1
n=[n]+1; ebeam[n]=550   ; il[n]=485.7
n=[n]+1; ebeam[n]=562.5 ; il[n]=517.3
n=[n]+1; ebeam[n]=575   ; il[n]=419.1
n=[n]+1; ebeam[n]=587.5 ; il[n]=534.4
n=[n]+1; ebeam[n]=600   ; il[n]=489.8999940
n=[n]+1; ebeam[n]=612.5 ; il[n]=546.7000120
n=[n]+1; ebeam[n]=625   ; il[n]=440.5000000
n=[n]+1; ebeam[n]=637.5 ; il[n]=496.2000120
n=[n]+1; ebeam[n]=650   ; il[n]=456.1000060
n=[n]+1; ebeam[n]=662.5 ; il[n]=526.2999880
n=[n]+1; ebeam[n]=675   ; il[n]=559.7999880
n=[n]+1; ebeam[n]=687.5 ; il[n]=576.2000120
n=[n]+1; ebeam[n]=700   ; il[n]=577.0999760
n=[n]+1; ebeam[n]=712.5 ; il[n]=593.2000120
n=[n]+1; ebeam[n]=725   ; il[n]=432.2999880
n=[n]+1; ebeam[n]=737.5 ; il[n]=600.0999760
n=[n]+1; ebeam[n]=750   ; il[n]=694.7000120
n=[n]+1; ebeam[n]=762.5 ; il[n]=487.2999880
n=[n]+1; ebeam[n]=775   ; il[n]=546.5000000
n=[n]+1; ebeam[n]=787.5 ; il[n]=502.3999940
n=[n]+1; ebeam[n]=800   ; il[n]=439.3999940
n=[n]+1; ebeam[n]=812.5 ; il[n]=509.8999940
n=[n]+1; ebeam[n]=825   ; il[n]=475.6000060
n=[n]+1; ebeam[n]=850   ; il[n]=462.2000120
n=[n]+1; ebeam[n]=862.5 ; il[n]=504.6000060
n=[n]+1; ebeam[n]=875   ; il[n]=506.3999940
n=[n]+1; ebeam[n]=887.5 ; il[n]=525.7000120
n=[n]+1; ebeam[n]=900   ; il[n]=384.1000060
n=[n]+1; ebeam[n]=912.5 ; il[n]=481.3999940
n=[n]+1; ebeam[n]=925   ; il[n]=401.2999880
n=[n]+1; ebeam[n]=935   ; il[n]=637.5000000
n=[n]+1; ebeam[n]=945   ; il[n]=577.0999760
n=[n]+1; ebeam[n]=950   ; il[n]=451.2999880
n=[n]+1; ebeam[n]=962.5 ; il[n]=561.9000240
n=[n]+1; ebeam[n]=987.5 ; il[n]=467.5000000
n=[n]+1; ebeam[n]=1000  ; il[n]=538.0999760
n2=[n]

dir=/work/users/konctbel/MinuitTest

gl/imp mname
gl/imp nb

ncnt=$sigma([nh]-10*int([nh]/10))
mode=$sigma(int([nh]/1000))

if ([prof].eq.0) then

  if ([ncnt].ne.0) then

    ci=0
    hi/del [nh]
    do i=[n1],[n2]
      ebeam=[ebeam[i]]
      hfile=hist_[ebeam]_prof.his
      hi/file 20 [dir]/[hfile]
      ci=$sigma(mod([i],5)+1)
      set pmci [ci]
      if ([i].eq.[n1]) then
        hi/pl //lun20/[nh]
        hi/copy [nh] 10000
      else
        hi/pl //lun20/[nh] s
        hi/op/add [nh] 10000 10000
      endif
      close 20
    enddo

  else

    ci=0
    hi/del [nh]
    do i=[n1],[n2]
      ebeam=[ebeam[i]]
      hfile=hist_[ebeam]_prof.his
      hi/file 20 [dir]/[hfile]
      ci=$sigma(mod([i],5)+1)
      set pmci [ci]
      do j=1,9
        nhj=[nh]+[j]
        if (([i].eq.[n1]).and.([j].eq.1)) then
          hi/pl //lun20/[nhj]
          hi/copy [nhj] 10000
        else
          hi/pl //lun20/[nhj] s
          hi/op/add [nhj] 10000 10000
        endif
      enddo
      close 20
    enddo

  endif

else

*  ve/del ic,if1,if2,if3,if4,if5
*  ve/read ic,if1,if2,if3,if4,if5 phi_cuts.txt
  if ([mode].eq.1) then
    f1=if1([ncnt])
    f2=if2([ncnt])
  endif
  if ([mode].eq.2) then
    f1=if3([ncnt])
    f2=if4([ncnt])
  endif
  if ([mode].eq.3) then
    f1=if4([ncnt])
    f2=if5([ncnt])
  endif
  gl/imp nh1
  gl/imp nh2
  mess [nh1] [nh2] [f1] [f2]
  
  exec mapcal#histprepf [ncnt] [nh1] [nh2] [f1] [f2]
  idh=10000
  hi/copy 200 [idh]
  
*  mnz=300
*  mzmin=-15
*  mzmax=15
*  mnf=320
*  mfmax=[ncnt]*40+20
*  mfmin=[ncnt]*40-60
*  i1=$sigma([mnf]*([f1]-[mfmin])/([mfmax]-[mfmin]))
*  j=[i]+1
*  i2=$sigma([mnf]*([f2]-[mfmin])/([mfmax]-[mfmin]))
*  ic=$sigma(int(([i1]+[i2])/2+0.5))
*  di=$sigma(int(([i2]-[i1])/2+0.5))
*  mess [ic] [di] [i1] [i2] [f1] [f2]

*  hi/file 30 [mname]_profilex_exp_pds.his

*  idh=10000
*  if ($hexist([idh]).eq.1) then
*    hi/del [idh]
*  endif
*  i1=$sigma(int([i1]+0.5))
*  i2=$sigma(int([i2]+0.5))
*  do j=[i1],[i2]

*    idhj=[j]*100+10*0+[ncnt]
    
*    if ($hexist([idh]).eq.0) then
*      hi/copy [idhj] [idh]
*    else
*      hi/op/add [idhj] [idh] [idh]
*    endif
    
*  enddo

*  close 30

endif

set pmci 1
set fcol 2
set fwid 3
*read x
hi/pl 10000

if (([nh].lt.200).or.([nh].gt.1000)) then

*  pfile=/work/users/konctbel/SepPar/mhad2011_h[nh]_fit.pars
*  if ($fexist([pfile])) then
*    ve/del p30,s30,pmin,pmax
*    ve/read p30,s30,pmin,pmax [pfile] '4f15.10'
*  endif
*  ve/cre s30(34) r 34*0
*  hi/fit 10000(-13.:13.) /work/users/konctbel/SepPar/uzfit.f b 34 p30 s30 pmin pmax
*  do i=1,15
*    line $sigma(p30([i])) 0 $sigma(p30([i])) 20
*  enddo
*  atitle 'z?r!, cm' 'A, pe'
*  exec save mhad2011_amplitude_vs_zr_counter[nh].eps f
*  read x
*=====================================================

  exec ../SepPar/sp#brn 10000 2 100 ! 1
*  exec hsigma @200 = %10000**2
*  exec ../SepPar/sp#brn 200 2 100
*  exec hsigma @10200 = @10100/2
*  exec hsigma %10200 = sqrt(@300)/2

  hexp=10100

  if ([sm].eq.'u') then
    z0=0
    dz0=0
*    sigz0=p30(30)
    sigz0=2
    dsigz0=0.1
*    alpha=p30(31)
*    alpha=0.55
    alpha=2
    beta=100
    dbeta=0.0
*    ag=p30(32)
    ag=2
    dag=0.1
*    mg=p30(33)
    mg=11
    dmg=0.01
*    sg=p30(34)
    sg=0.46
    dsg=0.01
    pds=0
    dpds=0.01
*    nb=14
    nr=[nb]+1
    ng=[nb]-1
    zil=-11.6
    zir=10.5
    ve/cre  zi([nr]) r -11.6 [ng]*1000 10.5
    ve/cre dzi([nr]) r  0.01 [ng]*0     0.01
*    ve/copy p30(1:15) zi(1:15)
*    zil=zi(1)
*    zir=zi(15)
    lambda=-60 ; dlambda=1
    xr1=-15    ; dxr1=1
    sr1=5      ; dsr1=1
    xr2=13     ; dxr2=1
    sr2=3      ; dsr2=1
*    sigma zi = array([nr],[zil]#[zir])
*    i=2; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0 
*    i=3; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0 
*    i=4; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0 
*    i=5; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0
*    i=6; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0 
*    i=7; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0 
*    i=8; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0 
*    i=9; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0 
*    i=10; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0
*    i=11; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0 
*    i=12; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0 
*    i=13; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0 
*    i=14; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0 
    ve/cre ziw([nr]) r [nr]*2
    ve/cre  ai([nb]) r [nb]*7
    ve/cre dai([nb]) r [nb]*0.1 
    ve/cre  ki([nb]) r 
    ve/cre dki([nb]) r 3*0 8*0.01 3*0
    ve/cre kiw([nb]) r 3*0 8*1 3*0
*    ve/copy p30(16:29) ai(1:14)
    ve/cre  gzi([ng]) r 
    ve/cre dgzi([ng]) r 
    ve/cre  alpi([nb]) r 1*[alpha] [ng]*1000
    ve/cre dalpi([nb]) r 1*0.01 [ng]*0
*ve/cre  alpi([nb]) r [nb]*[alpha]
*ve/cre dalpi([nb]) r [nb]*0.01
*    exec hsigma  ampi = @[hexp]
*    exec hsigma dampi = %[hexp]
*    exec hsigma   zri = $[hexp]
    np=$hinfo([hexp],'xbins')
    xmin=$hinfo([hexp],'xmin')
    xmax=$hinfo([hexp],'xmax')
    ve/cre ampi([np]) r
    hi/get/cont [hexp] ampi
    ve/cre dampi([np]) r
    hi/get/err [hexp] dampi
    l=$sigma([xmin]+([xmax]-[xmin])/[np]/2)
    r=$sigma([xmax]-([xmax]-[xmin])/[np]/2)
    ve/cre zri([np]) r
    sigma zri=array([np],[l]#[r])
    np=$vlen(ampi)
    ve/cre dzri([np]) r [np]*0.05
    nsc=1
    ve/cre inpoie([nsc]) i [np]
    sigma ai = ai/3.4
    ve/cre  ai([nb]) r 2.5 [ng]*1000
    ve/cre dai([nb]) r 0.1 [ng]*0
    
  endif

  if ([sm].eq.0) then

    nh0=$sigma(10*int([nh]/10))

    fname0=[mname]/[mname]_amplitude_vs_zr_counter[nh]_u
    shell /work/snd2000/users/konctbel/SNDSimACC/readlifstat [fname0].log
    ve/del pars0,dpars0
    ve/read pars0,dpars0 [fname0].log.pars


    z0=0
    dz0=0
*    sigz0=p30(30)
    sigz0=2
    dsigz0=0.1
*    alpha=p30(31)
*    alpha=0.55
    alpha=2
    beta=100
    dbeta=0.0
*    ag=p30(32)
    ag=2
    dag=0.1
*    mg=p30(33)
    mg=11
    dmg=0.01
*    sg=p30(34)
    sg=0.46
    dsg=0.01
    pds=0
    dpds=0.01
*    nb=14
    nr=[nb]+1
    ng=[nb]-1
    ve/cre  zi([nr]) r -11.6 [ng]*1000 10.5
    ve/cre dzi([nr]) r  0.01 [ng]*0     0.01
*    ve/copy p30(1:15) zi(1:15)
*    zil=zi(1)
*    zir=zi(15)
    ip=5*[nb]+8
    ip=[ip]+1; lambda=pars0([ip]) ; dlambda=5
    ip=[ip]+1; xr1=pars0([ip]) ;    dxr1=1
    ip=[ip]+1; sr1=pars0([ip]) ;    dsr1=0.5
    ip=[ip]+1; xr2=pars0([ip]) ;    dxr2=1
    ip=[ip]+1; sr2=pars0([ip]) ;    dsr2=0.5
*    sigma zi = array([nr],[zil]#[zir])
*    i=2; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0 
*    i=3; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0 
*    i=4; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0 
*    i=5; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0
*    i=6; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0 
*    i=7; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0 
*    i=8; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0 
*    i=9; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0 
*    i=10; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0
*    i=11; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0 
*    i=12; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0 
*    i=13; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0 
*    i=14; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0 
    ve/cre ziw([nr]) r [nr]*2
    ve/cre  ai([nb]) r [nb]*7
    ve/cre dai([nb]) r [nb]*0.1
    dd=3
    dnb=[nb]-2*[dd]
    ve/cre  ki([nb]) r 
    ve/cre dki([nb]) r [dd]*0 [dnb]*0.01 [dd]*0
    ve/cre kiw([nb]) r [dd]*0 [dnb]*1    [dd]*0
*    ve/copy p30(16:29) ai(1:14)
    ve/cre  gzi([ng]) r 
    ve/cre dgzi([ng]) r 
    ve/cre  alpi([nb]) r 1*[alpha] [ng]*1000
    ve/cre dalpi([nb]) r 1*0.01 [ng]*0
*ve/cre  alpi([nb]) r [nb]*[alpha]
*ve/cre dalpi([nb]) r [nb]*0.01
*    exec hsigma  ampi = @[hexp]
*    exec hsigma dampi = %[hexp]
*    exec hsigma   zri = $[hexp]
    np=$hinfo([hexp],'xbins')
    xmin=$hinfo([hexp],'xmin')
    xmax=$hinfo([hexp],'xmax')
    ve/cre ampi([np]) r
    hi/get/cont [hexp] ampi
    ve/cre dampi([np]) r
    hi/get/err [hexp] dampi
    l=$sigma([xmin]+([xmax]-[xmin])/[np]/2)
    r=$sigma([xmax]-([xmax]-[xmin])/[np]/2)
    ve/cre zri([np]) r
    sigma zri=array([np],[l]#[r])
    np=$vlen(ampi)
    ve/cre dzri([np]) r [np]*0.05
    nsc=1
    ve/cre inpoie([nsc]) i [np]
    sigma ai = ai/3.4
    ip=[nr]+[ng]+2
    ve/cre  ai([nb]) r [nb]*$sigma(pars0([ip]))
    ve/cre dai([nb]) r [nb]*$sigma(dpars0([ip])*3)

  endif

  fname=[mname]_amplitude_vs_zr_counter[nh]_[sm]

  if ([sm].eq.1) then
    exec mapcal#fitpl [nh]

    fname0=[mname]/[mname]_amplitude_vs_zr_counter[nh]_0
    shell /work/snd2000/users/konctbel/SNDSimACC/readlifstat [fname0].log
    ve/del pars0,dpars0
    ve/read pars0,dpars0 [fname0].log.pars

    ip=5*[nb]+1
    ip=[ip]+1; z0=pars0([ip])
    dz0=dpars0([ip])
    ip=[ip]+1; sigz0=pars0([ip])
    dsigz0=dpars0([ip])
    ip=[ip]+1; ag=pars0([ip])
    dag=dpars0([ip])
    ip=[ip]+1; mg=pars0([ip])
    dmg=dpars0([ip])
    ip=[ip]+1; sg=pars0([ip])
    dsg=dpars0([ip])
    ip=[ip]+1; pds=pars0([ip])
    dpds=dpars0([ip])
    ip=[ip]+1; beta=pars0([ip])
    dbeta=10
    dbeta=0.1
    ip=[ip]+1; lambda=pars0([ip]) ; dlambda=5
    ip=[ip]+1; xr1=pars0([ip]) ;    dxr1=0
    ip=[ip]+1; sr1=pars0([ip]) ;    dsr1=0
    ip=[ip]+1; xr2=pars0([ip]) ;    dxr2=0
    ip=[ip]+1; sr2=pars0([ip]) ;    dsr2=0
    nr=[nb]+1
    ng=[nb]-1
    ve/cre  zi([nr]) r
    ve/cre dzi([nr]) r
    lp0=1    ; lp=[lp0]+1
    rp0=[nr] ; rp=[rp0]+1
    ve/copy  pars0([lp]:[rp])  zi([lp0]:[rp0])
    ve/copy dpars0([lp]:[rp]) dzi([lp0]:[rp0])
*i=2; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0 
*i=3; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0 
*i=4; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0 
    do i=1,[nr]
      if ($sigma(zi([i])).eq.1000) then
        ve/inp zi([i]) $sigma(zi(1)+(zi([nr])-zi(1))*([i]-1)/([nr]-1))
        ve/inp dzi([i]) 0.01
      endif
      ve/inp dzi([i]) 0.01
    enddo
*i=12; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0 
*i=13; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0 
*i=14; ve/inp zi([i]) 1000; ve/inp dzi([i]) 0 
    ve/cre ziw([nr]) r 1*1 2*0.2 9*0.5 2*0.2 1*1
    ve/cre ziw([nr]) r 15*0.2
    ve/cre  ai([nb]) r 
    ve/cre dai([nb]) r 
    ve/copy  pars0(30:43)  ai(1:14)
    ve/copy dpars0(30:43) dai(1:14)
    ve/cre  ki(14) r 
    ve/cre dki(14) r 
    ve/copy  pars0(44:57)  ki(1:14)
    ve/copy dpars0(44:57) dki(1:14)
    ve/inp ki(1) 0
    ve/inp ki(2) 0
    ve/inp ki(3) 0
    ve/inp ki(12) 0
    ve/inp ki(13) 0
    ve/inp ki(14) 0
    ve/inp dki(1) 0.1
    ve/inp dki(2) 0.1
    ve/inp dki(3) 0.1
    ve/inp dki(12) 0.1
    ve/inp dki(13) 0.1
    ve/inp dki(14) 0.1
    ve/cre kiw(14) r 3*0.5 8*1 3*0.5
    ve/cre  gzi(13) r 
    ve/cre dgzi(13) r
    ve/copy  pars0(17:29)  gzi(1:13)
    ve/copy dpars0(17:29) dgzi(1:13)
    ve/cre  alpi(14) r 
    ve/cre dalpi(14) r
    ve/copy  pars0(58:71)  alpi(1:14)
    ve/copy dpars0(58:71) dalpi(1:14)
*ve/cre  alpi(14) r 14*$sigma(pars0(58))
*ve/cre dalpi(14) r 14*$sigma(abs(dpars0(58)))
*exec hsigma alpi = sqrt(alpi)
*    exec hsigma  ampi = @[hexp]
*    exec hsigma dampi = %[hexp]
*    exec hsigma   zri = $[hexp]
    np=$hinfo([hexp],'xbins')
    xmin=$hinfo([hexp],'xmin')
    xmax=$hinfo([hexp],'xmax')
    ve/cre ampi([np]) r
    hi/get/cont [hexp] ampi
    ve/cre dampi([np]) r
    hi/get/err [hexp] dampi
    l=$sigma([xmin]+([xmax]-[xmin])/[np]/2)
    r=$sigma([xmax]-([xmax]-[xmin])/[np]/2)
    ve/cre zri([np]) r
    sigma zri=array([np],[l]#[r])
    np=$vlen(ampi)
    ve/cre dzri([np]) r [np]*0.05
    nsc=1
    ve/cre inpoie([nsc]) i [np]

  endif


*----- for-file ------

file=[fname].for
if ($fexist([file]).eq.1) then
  shell rm [file]
endif
for/file  20 [file] new
close 20

fmess '      implicit none' [file]
tmp=$unquote('      ')include 'fitdata.inc'
fmess [tmp] [file]
tmp=$unquote('      ')include 'fitdata1.inc'
fmess [tmp] [file]
fmess '' [file]

fmess '      integer inproc' [file]
tmp=$unquote('      ')parameter ( inproc =  [nsc] )
fmess [tmp] [file]
fmess '      integer iarrlen' [file]
tmp=$unquote('      ')parameter ( iarrlen =  [np] )
fmess [tmp] [file]
fmess '      integer istatmode' [file]
fmess '      parameter ( istatmode =  1 )' [file]
fmess '      character*(*) cminfile' [file]
tmp=$unquote('      ')parameter ( cminfile = $quote([fname]) )
fmess [tmp] [file]
fmess '      character*(*) chisfile' [file]
tmp=$unquote('      ')parameter ( chisfile = $quote([fname]) )
fmess [tmp] [file]
fmess '      integer initer / 100000000 /' [file]
fmess '' [file]

fmess '      real*8 demaxrg(inproc) /1.0/' [file]
tmp=$unquote('      ')integer inpoie(inproc) / [np] /
fmess [tmp] [file]

fmess '      integer ifitmode(inproc) / 2001 /' [file]

fmess '      common/pars0/shift,kx,nb' [file]
fmess '      integer nb' [file]
fmess '      real*8 shift,kx' [file]
fmess '      integer i' [file]

fmess '' [file]
fmess '      real*8 decm0(iarrlen)' [file]
fmess '      real*8 dsige0(iarrlen)' [file]
fmess '      real*8 dlum(iarrlen)' [file]
fmess '      real*8 derrlum(iarrlen)' [file]
fmess '      real*8 dlivtim(iarrlen)' [file]
fmess '      real*8 dnev(iarrlen)' [file]
fmess '      real*8 derrnevdn(iarrlen)' [file]
fmess '      real*8 derrnevup(iarrlen)' [file]
fmess '      real*8 deff(iarrlen)' [file]
fmess '      real*8 derreff(iarrlen)' [file]
fmess '      real*8 derre0(iarrlen)' [file]
fmess '' [file]

fmess '      data decm0 /' [file]
exec mapcal#vewrite [file] zri f8.3 $vdim(zri) 5
fmess '     &  /' [file]

fmess '      data dsige0 /' [file]
exec mapcal#vewrite [file] dzri f8.3 $vdim(dzri) 5
fmess '     &  /' [file]

fmess '      data dlum / iarrlen*1 /' [file]
fmess '      data derrlum / iarrlen*0 /' [file]
fmess '      data dlivtim / iarrlen*1 /' [file]

fmess '      data dnev /' [file]
exec mapcal#vewrite [file] ampi f8.3 $vdim(ampi) 5
fmess '     &  /' [file]

sigma dampi = dampi + (1-dampi/dampi)*1.0
fmess '      data derrnevup /' [file]
exec mapcal#vewrite [file] dampi f8.3 $vdim(dampi) 5
fmess '     &  /' [file]

fmess '      data derrnevdn /' [file]
exec mapcal#vewrite [file] dampi f8.3 $vdim(dampi) 5
fmess '     &  /' [file]

fmess '      data deff / iarrlen*1 /' [file]
fmess '      data derreff / iarrlen*0 /' [file]
fmess '      data derre0 / iarrlen*0.05 /' [file]
fmess '' [file]

tmp=$unquote('      ')nb =  $vdim(ai)
fmess [tmp] [file]
tmp=$unquote('      ')shift =  2000
fmess [tmp] [file]
tmp=$unquote('      ')kx =  100
fmess [tmp] [file]
fmess '' [file]

fmess '      do i=1,iarrlen' [file]
fmess '        decm0(i)=kx*decm0(i)+shift' [file]
fmess '        dsige0(i)=kx*dsige0(i)' [file]
*fmess '        dnev(i)=dnev(i)+10' [file]
fmess '      enddo' [file]
fmess '' [file]

fmess '      call inifit' [file]
fmess '' [file]
fmess '      call fitbegin(' [file]
fmess '     &  inproc,' [file]
fmess '     &  initer,' [file]
fmess '     &  ifitmode,' [file]
fmess '     &  inpoie,' [file]
fmess '     &  iarrlen,' [file]
fmess '     &  dnev,' [file]
fmess '     &  decm0,' [file]
fmess '     &  derre0,' [file]
fmess '     &  dsige0,' [file]
fmess '     &  dlum,' [file]
fmess '     &  derrlum,' [file]
fmess '     &  dlivtim,' [file]
fmess '     &  deff,' [file]
fmess '     &  derreff,' [file]
fmess '     &  derrnevup,' [file]
fmess '     &  derrnevdn,' [file]
fmess '     &  demaxrg,' [file]
fmess '     &  cminfile,' [file]
fmess '     &  chisfile )' [file]
fmess '' [file]
fmess '      end' [file]
fmess '' [file]
tmp=$unquote('      ')include 'uzfit_tail.for'
fmess [tmp] [file]


*--- dat-file ------

*shell /work/snd2000/users/konctbel/SNDSimACC/readlifstat dedx_[ver0].lif
*ve/del pars,dpars
*ve/read pars,dpars dedx_[ver0].lif.pars

file=[fname].dat
if ($fexist([file]).eq.1) then
  shell rm -v [file]
endif
for/file  20 [file] new
close 20

fmess 'SET TITLE' [file]
ni=0
fmess 'Z-profile fit of amplitude' [file]
fmess 'PARAMETERS' [file]
ni1z=[ni]+1
do i=1,$vdim(zi)
  tmp0=z_[i]
  zix=$sigma(zi([i]))
  if (([i].gt.1).and.([i].lt.$vdim(zi))) then
    zixl=[zix]-0.4
    zixu=[zix]+0.4
  else
    zixl=[zix]-2.
    zixu=[zix]+2.
  endif
  zixl=[zix]-ziw([i])
  zixu=[zix]+ziw([i])
  dzix=$sigma(dzi([i]))
  ni=[ni]+1; pname=z_[i]; tmp1=$format([ni],i-4); tmp2=$format([zix],e13.5); tmp3=$format([dzix],e13.5); tmp4=$format([zixl],e13.5); tmp5=$format([zixu],e13.5);
  tmp=$unquote([tmp1])$quote([pname])   $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
  fmess [tmp] [file]
  ni2z=[ni]
enddo
ni1g=[ni]+1
do i=1,$vdim(gzi)
  tmp0=gz_[i]
  gzix=$sigma(gzi([i]))
  gzixl=0.0
  gzixu=0.05
  dgzix=$sigma(dgzi([i]))
  ni=[ni]+1; pname=gz_[i]; tmp1=$format([ni],i-4); tmp2=$format([gzix],e13.5); tmp3=$format([dgzix],e13.5); tmp4=$format([gzixl],e13.5); tmp5=$format([gzixu],e13.5);
  tmp=$unquote([tmp1])$quote([pname])   $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
  fmess [tmp] [file]
  ni2g=[ni]
enddo
do i=1,$vdim(ai)
  tmp0=a_[i]
  aix=$sigma(ai([i]))
  aixl=0.0
  aixu=40.0
  daix=$sigma(dai([i]))
  ni=[ni]+1; pname=a_[i]; tmp1=$format([ni],i-4); tmp2=$format([aix],e13.5); tmp3=$format([daix],e13.5); tmp4=$format([aixl],e13.5); tmp5=$format([aixu],e13.5);
  tmp=$unquote([tmp1])$quote([pname])   $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
  fmess [tmp] [file]
enddo
ni1k=[ni]+1
do i=1,$vdim(ki)
  tmp0=k_[i]
  kix=$sigma(ki([i]))
  kixl=-kiw([i])
  kixu=+kiw([i])
  dkix=$sigma(dki([i]))
  ni=[ni]+1; pname=k_[i]; tmp1=$format([ni],i-4); tmp2=$format([kix],e13.5); tmp3=$format([dkix],e13.5); tmp4=$format([kixl],e13.5); tmp5=$format([kixu],e13.5);
  tmp=$unquote([tmp1])$quote([pname])   $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
  fmess [tmp] [file]
  ni2k=[ni]
enddo
ni1b=[ni]+1
do i=1,$vdim(alpi)
  tmp0=alpha_[i]
  alpix=$sigma(alpi([i]))
  alpixl=-10.0
  alpixu=10.0
  dalpix=$sigma(dalpi([i]))
  ni=[ni]+1; pname=alpha_[i]; tmp1=$format([ni],i-4); tmp2=$format([alpix],e13.5); tmp3=$format([dalpix],e13.5); tmp4=$format([alpixl],e13.5); tmp5=$format([alpixu],e13.5);
  tmp=$unquote([tmp1])$quote([pname])   $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
  fmess [tmp] [file]
  ni2b=[ni]
enddo
z0l=[z0]-1
z0u=[z0]+1
ni=[ni]+1; pname=z0; tmp1=$format([ni],i-4); tmp2=$format([z0],e13.5); tmp3=$format([dz0],e13.5); tmp4=$format([z0l],e13.5); tmp5=$format([z0u],e13.5);
tmp=$unquote([tmp1])$quote([pname])   $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
fmess [tmp] [file]
niz0=[ni]
sigz0l=0.
sigz0u=3.
ni=[ni]+1; pname=sigz0; tmp1=$format([ni],i-4); tmp2=$format([sigz0],e13.5); tmp3=$format([dsigz0],e13.5); tmp4=$format([sigz0l],e13.5); tmp5=$format([sigz0u],e13.5);
tmp=$unquote([tmp1])$quote([pname])   $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
fmess [tmp] [file]
nidz0=[ni]
agl=0
agu=30
ni=[ni]+1; pname=ag; tmp1=$format([ni],i-4); tmp2=$format([ag],e13.5); tmp3=$format([dag],e13.5); tmp4=$format([agl],e13.5); tmp5=$format([agu],e13.5);
tmp=$unquote([tmp1])$quote([pname])   $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
fmess [tmp] [file]
mgl=[mg]-1
mgu=[mg]+1
ni=[ni]+1; pname=mg; tmp1=$format([ni],i-4); tmp2=$format([mg],e13.5); tmp3=$format([dmg],e13.5); tmp4=$format([mgl],e13.5); tmp5=$format([mgu],e13.5);
tmp=$unquote([tmp1])$quote([pname])   $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
fmess [tmp] [file]
sgl=0
sgu=1
ni=[ni]+1; pname=sg; tmp1=$format([ni],i-4); tmp2=$format([sg],e13.5); tmp3=$format([dsg],e13.5); tmp4=$format([sgl],e13.5); tmp5=$format([sgu],e13.5);
tmp=$unquote([tmp1])$quote([pname])   $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
fmess [tmp] [file]
sgl=-0.1
sgu=0.1
ni=[ni]+1; pname=pds; tmp1=$format([ni],i-4); tmp2=$format([pds],e13.5); tmp3=$format([dpds],e13.5); tmp4=$format([sgl],e13.5); tmp5=$format([sgu],e13.5);
tmp=$unquote([tmp1])$quote([pname])   $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
fmess [tmp] [file]
sgl=0.001
sgu=10.0
ni=[ni]+1; pname=beta; tmp1=$format([ni],i-4); tmp2=$format([beta],e13.5); tmp3=$format([dbeta],e13.5); tmp4=$format([sgl],e13.5); tmp5=$format([sgu],e13.5);
tmp=$unquote([tmp1])$quote([pname])   $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
fmess [tmp] [file]
sgl=-100
sgu=100
ni=[ni]+1; pname=lambda; tmp1=$format([ni],i-4); tmp2=$format([lambda],e13.5); tmp3=$format([dlambda],e13.5); tmp4=$format([sgl],e13.5); tmp5=$format([sgu],e13.5);
tmp=$unquote([tmp1])$quote([pname])   $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
fmess [tmp] [file]
sgl=-30
sgu=-10
ni=[ni]+1; pname=xr1; tmp1=$format([ni],i-4); tmp2=$format([xr1],e13.5); tmp3=$format([dxr1],e13.5); tmp4=$format([sgl],e13.5); tmp5=$format([sgu],e13.5);
tmp=$unquote([tmp1])$quote([pname])   $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
fmess [tmp] [file]
sgl=1
sgu=10
ni=[ni]+1; pname=sr1; tmp1=$format([ni],i-4); tmp2=$format([sr1],e13.5); tmp3=$format([dsr1],e13.5); tmp4=$format([sgl],e13.5); tmp5=$format([sgu],e13.5);
tmp=$unquote([tmp1])$quote([pname])   $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
fmess [tmp] [file]
sgl=10
sgu=30
ni=[ni]+1; pname=xr2; tmp1=$format([ni],i-4); tmp2=$format([xr2],e13.5); tmp3=$format([dxr2],e13.5); tmp4=$format([sgl],e13.5); tmp5=$format([sgu],e13.5);
tmp=$unquote([tmp1])$quote([pname])   $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
fmess [tmp] [file]
sgl=1
sgu=10
ni=[ni]+1; pname=sr2; tmp1=$format([ni],i-4); tmp2=$format([sr2],e13.5); tmp3=$format([dsr2],e13.5); tmp4=$format([sgl],e13.5); tmp5=$format([sgu],e13.5);
tmp=$unquote([tmp1])$quote([pname])   $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
fmess [tmp] [file]

fmess '' [file]
fmess '' [file]
*fmess 'set err 0.5' [file]
*fmess 'set eps 1e-9' [file]
*fmess 'set strategy 2' [file]

tmp=fix [niz0]; fmess [tmp] [file]
tmp=fix [nidz0]; fmess [tmp] [file]

if ([sm].eq.'u') then
  do i=[ni1g],[ni2g]
    tmp=fix [i]
    fmess [tmp] [file]
  enddo
  do i=[ni1k],[ni2k]
    tmp=fix [i]
    fmess [tmp] [file]
  enddo
  do i=[ni1b],[ni2b]
    tmp=fix [i]
    fmess [tmp] [file]
  enddo
  fmess 'mini' [file]
  do i=[ni1b],[ni2b]
    tmp=release [i]
    fmess [tmp] [file]
  enddo
  fmess 'mini' [file]
  do i=[ni1z],[ni2z]
    tmp=fix [i]
    fmess [tmp] [file]
  enddo
  do i=[ni1k],[ni2k]
    tmp=release [i]
    fmess [tmp] [file]
  enddo
  fmess 'mini' [file]
endif

if ([sm].eq.0) then
  do i=[ni1g],[ni2g]
    tmp=fix [i]
    fmess [tmp] [file]
  enddo
  do i=[ni1k],[ni2k]
    tmp=fix [i]
    fmess [tmp] [file]
  enddo
  do i=[ni1b],[ni2b]
    tmp=fix [i]
    fmess [tmp] [file]
  enddo
  fmess 'mini' [file]
  do i=[ni1b],[ni2b]
    tmp=release [i]
    fmess [tmp] [file]
  enddo
  fmess 'mini' [file]
  do i=[ni1z],[ni2z]
    tmp=fix [i]
    fmess [tmp] [file]
  enddo
  do i=[ni1k],[ni2k]
    tmp=release [i]
    fmess [tmp] [file]
  enddo
  fmess 'mini' [file]
endif

if ([sm].eq.1) then
  do i=[ni1g],[ni2g]
    tmp=fix [i]
    fmess [tmp] [file]
  enddo
  do i=[ni1k],[ni2k]
    tmp=fix [i]
    fmess [tmp] [file]
  enddo
  do i=[ni1b],[ni2b]
    tmp=fix [i]
    fmess [tmp] [file]
  enddo
  fmess 'mini' [file]
  do i=[ni1g],[ni2g]
    tmp=release [i]
    fmess [tmp] [file]
  enddo
  do i=[ni1k],[ni2k]
    tmp=release [i]
    fmess [tmp] [file]
  enddo
  do i=[ni1b],[ni2b]
    tmp=release [i]
    fmess [tmp] [file]
  enddo
  fmess 'mini' [file]
endif

fmess 'ret' [file]

*=====================================================

else
endif
return


macro fitres
*exec crsfit load cs7 mhad2011_amplitude_vs_zr_counter2107
exec crsfit cs7 1 xD vx vex vy vdn vup
n1=100
n2=250
ve/del vxr,vexr,vyr,vdnr,vupr
ve/copy vx([n1]:[n2]) vxr
ve/copy vex([n1]:[n2]) vexr
ve/copy vy([n1]:[n2]) vyr
ve/copy vdn([n1]:[n2]) vdnr
ve/copy vup([n1]:[n2]) vupr
* exec $PER/s#vpl vyr vdnr vxr vexr ll=-1 iatt=20 sz=0.1
ve/fit vxr vyr vdnr p0 s
return

macro fitsum nhl=2 sm=0 corr=0

set mtyp 20
set ksiz 0.05
set hcol 1
set pmci 1
set hwid 1
set lwid 1

do i=1,9

  nh=[nhl]*1000+100+[i]
  fname=mhad2011_amplitude_vs_zr_counter[nh]_[sm]

  if ([corr].eq.1) then
    shell cp [fname]_exp.his tmp.txt
    shell $unquote('cat tmp.txt | sed "s/.00 0.05-0./.00 0.05 0./g" > tmp1.txt')
    shell cp tmp1.txt [fname]_exp.his
  endif

  ve/del vcs
  exec crsfit load vcs [fname]
  exec crsfit vcs 1 xc vx vex vy vdn vup
  if ([i].eq.1) then
    ve/del zr,dzr,az,daz,uaz
    ve/copy vx zr
    ve/copy vex dzr
    np=$vlen(zr)
    ve/cre az([np]) r
    ve/cre daz([np]) r
    ve/cre uaz([np]) r
  endif
  sigma az = az + vy
  sigma daz = daz + vdn**2
  sigma uaz = uaz + vup**2
enddo
sigma az = az/9
sigma daz = sqrt(daz)/9
sigma uaz = sqrt(uaz)/9

* exec $PER/s#vpl az daz zr dzr sz=0.1

return


macro fitplx n1=1 n2=9 sm=0
do i=[n1],[n2]
  nh=1100+[i]
  exec ../MinuitTest/mapcal#fitpl [nh] [sm]
enddo
do i=[n1],[n2]
  nh=2100+[i]
  exec ../MinuitTest/mapcal#fitpl [nh] [sm]
enddo
do i=[n1],[n2]
  nh=3100+[i]
  exec ../MinuitTest/mapcal#fitpl [nh] [sm]
enddo
return



macro fitpl nh=2107 sm=0 corr=0

fname=mhad2011-4_amplitude_vs_zr_counter[nh]_[sm]

if ([corr].eq.1) then
  shell cp [fname]_exp.his tmp.txt
  shell $unquote('cat tmp.txt | sed "s/.00 0.05-0./.00 0.05 0./g" > tmp1.txt')
  shell cp tmp1.txt [fname]_exp.his
endif

set mtyp 20
set ksiz 0.05
set hcol 1
set pmci 1
set hwid 1
set lwid 1

*file=i.am.busy
*while ($fexist([file]).eq.1) do
*  wait 'wait a second...' 1
*  shell date
*endwhile
*for/file  20 [file] new 
*close 20
xx=0
if ([xx].eq.1) then
ve/del vcs
exec crsfit load vcs [fname]
exec crsfit vcs 1 c
l=$GRAFINFO('WNXMIN')
r=$GRAFINFO('WNXMAX')
d=$GRAFINFO('WNYMIN')
u=$GRAFINFO('WNYMAX')
l=-17
r=17
null [l] [r] 0 [u]
null $sigma([l]*100+2000) $sigma([r]*100+2000) 0 [u] sab
exec crsfit vcs 1 cs
*exec crsfit vcs 1 ces
*shell rm [file]
endif

shell /work/snd2000/users/konctbel/SNDSimACC/readlifstat [fname].log
ve/del pars0,dpars0
ve/read pars0,dpars0 [fname].log.pars

z0=pars0(58)
dz0=dpars0(58)
sigz0=pars0(59)
dsigz0=dpars0(59)
ag=pars0(60)
dag=dpars0(60)
mg=pars0(61)
dmg=dpars0(61)
sg=pars0(62)
dsg=dpars0(62)
ve/cre  ai(14) r 
ve/cre dai(14) r 
ve/copy  pars0(30:43)  ai(1:14)
ve/copy dpars0(30:43) dai(1:14)
ve/cre  ki(14) r 
ve/cre dki(14) r 
ve/copy  pars0(44:57)  ki(1:14)
ve/copy dpars0(44:57) dki(1:14)
ve/cre  gzi(13) r 
ve/cre dgzi(13) r
ve/copy  pars0(17:29)  gzi(1:13)
ve/copy dpars0(17:29) dgzi(1:13)
ve/cre  alpi(14) r 
ve/cre dalpi(14) r
ve/copy  pars0(58:71)  alpi(1:14)
ve/copy dpars0(58:71) dalpi(1:14)

ve/cre  zi(15) r
ve/cre dzi(15) r 
ve/copy  pars0(2:16)  zi(1:15)
ve/copy dpars0(2:16) dzi(1:15)

kx=100
bx=2000

ve/cre zic(15) r

i1=-1
do i=1,$vdim(zi)
  z=zi([i])
  if (([z].ne.1000).and.([i1].eq.-1)) then
    zx=[kx]*[z]+[bx]
    ve/inp zic([i]) [zx]
    line [zx] 0 [zx] 100
  endif
  if (([z].eq.1000).and.([i1].eq.-1)) then
    i1=[i]-1
    z1=zi([i1])
  endif
  if (([z].ne.1000).and.([i1].ne.-1)) then
    i2=[i]
    z2=zi([i2])
    do j=[i1],[i2]
      z=$sigma([z1]+([z2]-[z1])*([j]-[i1])/([i2]-[i1]))
      zx=[kx]*[z]+[bx]
      ve/inp zic([j]) [zx]
      line [zx] 0 [zx] 100
    enddo
    i1=-1
  endif
enddo

sigma zic = (zic-[bx])/[kx]

*exec $PER/s#tf 0.85 0.9 [nh]_[sm]
exec pl#tf 0.85 0.9 [nh]_[sm]
atitle 'z?r!, cm' 'A, pe'
*exec save [fname].eps f

return


macro zicr vi=zi vo=zic

kx=100
bx=2000

n=$vlen([vi])
ve/cre [vo]([n]) r

i1=-1
do i=1,$vdim([vi])
  z=[vi]([i])
  if (([z].ne.1000).and.([i1].eq.-1)) then
    zx=[kx]*[z]+[bx]
    ve/inp [vo]([i]) [zx]
    line [zx] 0 [zx] 100
  endif
  if (([z].eq.1000).and.([i1].eq.-1)) then
    i1=[i]-1
    z1=[vi]([i1])
  endif
  if (([z].ne.1000).and.([i1].ne.-1)) then
    i2=[i]
    z2=[vi]([i2])
    do j=[i1],[i2]
      z=$sigma([z1]+([z2]-[z1])*([j]-[i1])/([i2]-[i1]))
      zx=[kx]*[z]+[bx]
      ve/inp [vo]([j]) [zx]
      line [zx] 0 [zx] 100
    enddo
    i1=-1
  endif
enddo

sigma [vo] = ([vo]-[bx])/[kx]

return


macro mapcmp nh=201 nh0=! dl=0
n=0
n=[n]+1; ebeam[n]=508.7 ; il[n]=0
n1=[n]
n=[n]+1; ebeam[n]=509.8 ; il[n]=0
n=[n]+1; ebeam[n]=510.5 ; il[n]=0
n=[n]+1; ebeam[n]=525   ; il[n]=363.2
n=[n]+1; ebeam[n]=537.5 ; il[n]=546.1
n=[n]+1; ebeam[n]=550   ; il[n]=485.7
n=[n]+1; ebeam[n]=562.5 ; il[n]=517.3
n=[n]+1; ebeam[n]=575   ; il[n]=419.1
n=[n]+1; ebeam[n]=587.5 ; il[n]=534.4
n=[n]+1; ebeam[n]=600   ; il[n]=489.8999940
n=[n]+1; ebeam[n]=612.5 ; il[n]=546.7000120
n=[n]+1; ebeam[n]=625   ; il[n]=440.5000000
n=[n]+1; ebeam[n]=637.5 ; il[n]=496.2000120
n=[n]+1; ebeam[n]=650   ; il[n]=456.1000060
n=[n]+1; ebeam[n]=662.5 ; il[n]=526.2999880
n=[n]+1; ebeam[n]=675   ; il[n]=559.7999880
n=[n]+1; ebeam[n]=687.5 ; il[n]=576.2000120
n=[n]+1; ebeam[n]=700   ; il[n]=577.0999760
n=[n]+1; ebeam[n]=712.5 ; il[n]=593.2000120
n=[n]+1; ebeam[n]=725   ; il[n]=432.2999880
n=[n]+1; ebeam[n]=737.5 ; il[n]=600.0999760
n=[n]+1; ebeam[n]=750   ; il[n]=694.7000120
n=[n]+1; ebeam[n]=762.5 ; il[n]=487.2999880
n=[n]+1; ebeam[n]=775   ; il[n]=546.5000000
n=[n]+1; ebeam[n]=787.5 ; il[n]=502.3999940
n=[n]+1; ebeam[n]=800   ; il[n]=439.3999940
n=[n]+1; ebeam[n]=812.5 ; il[n]=509.8999940
n=[n]+1; ebeam[n]=825   ; il[n]=475.6000060
n=[n]+1; ebeam[n]=850   ; il[n]=462.2000120
n=[n]+1; ebeam[n]=862.5 ; il[n]=504.6000060
n=[n]+1; ebeam[n]=875   ; il[n]=506.3999940
n=[n]+1; ebeam[n]=887.5 ; il[n]=525.7000120
n=[n]+1; ebeam[n]=900   ; il[n]=384.1000060
n=[n]+1; ebeam[n]=912.5 ; il[n]=481.3999940
n=[n]+1; ebeam[n]=925   ; il[n]=401.2999880
n=[n]+1; ebeam[n]=935   ; il[n]=637.5000000
n=[n]+1; ebeam[n]=945   ; il[n]=577.0999760
n=[n]+1; ebeam[n]=950   ; il[n]=451.2999880
n=[n]+1; ebeam[n]=962.5 ; il[n]=561.9000240
n=[n]+1; ebeam[n]=987.5 ; il[n]=467.5000000
n=[n]+1; ebeam[n]=1000  ; il[n]=538.0999760
n2=[n]

n=0
n=[n]+1; ebeam[n]=750   ; il[n]=694.7000120   ; nfirst[n]=7842  ; nlast[n]=9258 ;
n1=[n]
n=[n]+1; ebeam[n]=525   ; il[n]=363.2         ; nfirst[n]=8368  ; nlast[n]=8537 ;
n=[n]+1; ebeam[n]=550   ; il[n]=485.7         ; nfirst[n]=8538  ; nlast[n]=8643 ;
n=[n]+1; ebeam[n]=575   ; il[n]=419.1         ; nfirst[n]=8650  ; nlast[n]=8712 ;
n=[n]+1; ebeam[n]=600   ; il[n]=489.8999940   ; nfirst[n]=8715  ; nlast[n]=8768 ;
n=[n]+1; ebeam[n]=625   ; il[n]=440.5000000   ; nfirst[n]=8797  ; nlast[n]=8973 ;
n=[n]+1; ebeam[n]=650   ; il[n]=456.1000060   ; nfirst[n]=8978  ; nlast[n]=9032 ;
n=[n]+1; ebeam[n]=675   ; il[n]=559.7999880   ; nfirst[n]=9033  ; nlast[n]=9097 ;
n=[n]+1; ebeam[n]=700   ; il[n]=577.0999760   ; nfirst[n]=9100  ; nlast[n]=9152 ;
n=[n]+1; ebeam[n]=725   ; il[n]=432.2999880   ; nfirst[n]=9153  ; nlast[n]=9193 ;
n=[n]+1; ebeam[n]=775   ; il[n]=546.5000000   ; nfirst[n]=9261  ; nlast[n]=9320 ;
n=[n]+1; ebeam[n]=800   ; il[n]=439.3999940   ; nfirst[n]=9321  ; nlast[n]=9366 ;
n=[n]+1; ebeam[n]=825   ; il[n]=475.6000060   ; nfirst[n]=9367  ; nlast[n]=9435 ;
n=[n]+1; ebeam[n]=850   ; il[n]=462.2000120   ; nfirst[n]=9443  ; nlast[n]=9489 ;
n=[n]+1; ebeam[n]=875   ; il[n]=506.3999940   ; nfirst[n]=9492  ; nlast[n]=9550 ;
n=[n]+1; ebeam[n]=900   ; il[n]=384.1000060   ; nfirst[n]=9566  ; nlast[n]=9595 ;
n=[n]+1; ebeam[n]=925   ; il[n]=401.2999880   ; nfirst[n]=9598  ; nlast[n]=9651 ;
n=[n]+1; ebeam[n]=950   ; il[n]=451.2999880   ; nfirst[n]=9653  ; nlast[n]=9738 ;
n=[n]+1; ebeam[n]=962.5 ; il[n]=561.9000240   ; nfirst[n]=9739  ; nlast[n]=9825 ;
n=[n]+1; ebeam[n]=987.5 ; il[n]=467.5000000   ; nfirst[n]=9932  ; nlast[n]=10033 ;
n=[n]+1; ebeam[n]=1000  ; il[n]=538.0999760   ; nfirst[n]=10035 ; nlast[n]=10127 ;
n=[n]+1; ebeam[n]=945   ; il[n]=577.0999760   ; nfirst[n]=10162 ; nlast[n]=10192 ;
n=[n]+1; ebeam[n]=935   ; il[n]=637.5000000   ; nfirst[n]=10225 ; nlast[n]=10273 ;
n=[n]+1; ebeam[n]=912.5 ; il[n]=481.3999940   ; nfirst[n]=10275 ; nlast[n]=10295 ;
n=[n]+1; ebeam[n]=887.5 ; il[n]=525.7000120   ; nfirst[n]=10296 ; nlast[n]=10315 ;
n=[n]+1; ebeam[n]=862.5 ; il[n]=504.6000060   ; nfirst[n]=10316 ; nlast[n]=10332 ;
n=[n]+1; ebeam[n]=812.5 ; il[n]=509.8999940   ; nfirst[n]=10356 ; nlast[n]=10383 ;
n=[n]+1; ebeam[n]=787.5 ; il[n]=502.3999940   ; nfirst[n]=10384 ; nlast[n]=10424 ;
n=[n]+1; ebeam[n]=712.5 ; il[n]=593.2000120   ; nfirst[n]=10429 ; nlast[n]=10466 ;
n=[n]+1; ebeam[n]=737.5 ; il[n]=600.0999760   ; nfirst[n]=10467 ; nlast[n]=10513 ;
n=[n]+1; ebeam[n]=762.5 ; il[n]=487.2999880   ; nfirst[n]=10514 ; nlast[n]=10555 ;
n=[n]+1; ebeam[n]=687.5 ; il[n]=576.2000120   ; nfirst[n]=10559 ; nlast[n]=10614 ;
n=[n]+1; ebeam[n]=662.5 ; il[n]=526.2999880   ; nfirst[n]=10615 ; nlast[n]=10673 ;
n=[n]+1; ebeam[n]=637.5 ; il[n]=496.2000120   ; nfirst[n]=10703 ; nlast[n]=10745 ;
n=[n]+1; ebeam[n]=612.5 ; il[n]=546.7000120   ; nfirst[n]=10746 ; nlast[n]=10784 ;
n=[n]+1; ebeam[n]=587.5 ; il[n]=534.4         ; nfirst[n]=10785 ; nlast[n]=10825 ;
n=[n]+1; ebeam[n]=562.5 ; il[n]=517.3         ; nfirst[n]=10826 ; nlast[n]=10863 ;
n=[n]+1; ebeam[n]=537.5 ; il[n]=546.1         ; nfirst[n]=10864 ; nlast[n]=10902 ;
n=[n]+1; ebeam[n]=508.7 ; il[n]=0             ; nfirst[n]=10921 ; nlast[n]=10932 ;
n=[n]+1; ebeam[n]=509.8 ; il[n]=0             ; nfirst[n]=10934 ; nlast[n]=10952 ;
n=[n]+1; ebeam[n]=510.5 ; il[n]=0             ; nfirst[n]=10955 ; nlast[n]=10988 ;
n2=[n]

dir=/work/users/konctbel/MinuitTest

if ([nh0].eq.'!') then
  nh0=[nh]
  mess nh0=[nh]
endif

ve/cre prob([n]) r
ve/cre probd([n]) r
ve/cre prob0([n]) r
ve/cre iprb([n]) r
fname=mhad2011_[nh0].txt
if ([dl].eq.1) then
  shell rm [fname]
  shell rm mhad2011_[nh0]_v1.txt
endif
if ($fexist([fname]).eq.1) then
  ve/del prob,probd,prob0
  ve/read prob,probd,prob0 [fname] '3f10.6'
  ve/prin prob0
  fname=mhad2011_[nh0]_v1.txt
  if ($fexist([fname]).eq.1) then
    ve/del prob,probd,prob0
    ve/read prob,probd,prob0 [fname] '3f10.6'
    ve/prin prob0
  endif
endif

close 20

ci=0
hi/del [nh]
ix=0
hi/del 40000
hi/del 50000
do i=[n1],[n2]
  ebeam=[ebeam[i]]
  hfile=hist_[ebeam]_prof.his
  hi/file 20 [dir]/[hfile]
  ci=$sigma(mod([i],5)+1)
  set pmci [ci]
  l=-15.01
  r=15.01
  if ([i].eq.[n1]) then
    hi/pl //lun20/[nh]
    prb=prob0([i])
    if (([prb].gt.0.9).or.([prb].eq.0)) then
      ix=1
      hi/copy [nh] 10000
      ve/inp iprb([i]) 1
    endif
    hi/copy [nh] 10001
    hi/copy [nh] 1000
    hi/del 1001
*    exec hsigma @1001 = @1000*($1000 gt [l] and $1000 lt [r])
*    exec hsigma %1001 = %1000*($1000 gt [l] and $1000 lt [r])
    exec hsigma @1001 = @1000
    exec hsigma %1001 = %1000
    hi/copy 1001 1003
  else
    hi/pl //lun20/[nh] s
    prb=prob0([i])
    if (([prb].gt.0.9).or.([prb].eq.0)) then
      if ([ix].eq.0) then
        ix=1
        hi/copy [nh] 10000
      else
        hi/op/add [nh] 10000 10000
      endif
      ve/inp iprb([i]) 1
      mess [i]
    else
      if ($hexist(40000).eq.0) then
        hi/copy [nh] 40000
      else
        hi/op/add [nh] 40000 40000
      endif
      if ([prb].lt.0.1) then
        if ($hexist(50000).eq.0) then
          hi/copy [nh] 50000
        else
          hi/op/add [nh] 50000 50000
        endif
      endif
    endif
  endif
  set pmci 2
  hi/pl 10001
  set pmci 4
  hi/pl [nh] s
  hi/copy [nh] 10002
  hi/copy [nh] 2000
  hi/del 1002
*  exec hsigma @1002 = @2000*($2000 gt [l] and $2000 lt [r])
*  exec hsigma %1002 = %2000*($2000 gt [l] and $2000 lt [r])
  exec hsigma @1002 = @2000
  exec hsigma %1002 = %2000
  ve/inp prob([i]) $call('mhdiff(1001,1002)')
  ve/inp probd([i]) $call('mhdiff(1003,1002)')
  mess $sigma(prob([i]))
*  read x
  close 20
  hi/copy 1002 1003
enddo

hi/del 1001
exec hsigma @1001 = @10000
exec hsigma %1001 = %10000

ci=0
hi/del [nh]
do i=[n1],[n2]
  ebeam=[ebeam[i]]
  hfile=hist_[ebeam]_prof.his
  hi/file 20 [dir]/[hfile]
  ci=$sigma(mod([i],5)+1)
  set pmci [ci]
  hi/pl //lun20/[nh]
  hi/copy [nh] 2000
  hi/del 1002
  exec hsigma @1002 = @2000
  exec hsigma %1002 = %2000
  ve/inp prob0([i]) $call('mhdiff(1001,1002)')
*  read x
  close 20
enddo

ve/draw prob
ve/cre dprob([n]) r
ve/cre np([n]) r
sigma np = array([n],1#[n])
* exec $PER/s#vpl prob dprob np dprob
* exec $PER/s#vpl probd dprob np dprob o=s iatt=24

if ($fexist([fname]).eq.0) then
  ve/write prob,probd,prob0 [fname] '3f10.6'
endif 

return


macro mapcmpx n1=1 n2=9
do i=[n1],[n2]
  nh=1100+[i]
  exec ../MinuitTest/mapcal#mapcmp [nh] dl=1
  exec ../MinuitTest/mapcal#mapcmp [nh]
  exec ../MinuitTest/mapcal#mapcmp [nh]
  * exec $PER/s#vpl prob0 dprob np dprob
  atitle 'n?beam!' 'prob'
  exec save mhad2011_prob0_[nh].eps f
  set pmci 1
  set mtyp 20
  hi/pl 10000
  set pmci 1
  set mtyp 24
  hi/pl 50000 s
  atitle 'z?r!, cm' 'A, pe.'
  exec save mhad2011_amplitude_vs_zr_counter[nh]_best.eps f
enddo
do i=[n1],[n2]
  nh=2100+[i]
  exec ../MinuitTest/mapcal#mapcmp [nh] dl=1
  exec ../MinuitTest/mapcal#mapcmp [nh]
  exec ../MinuitTest/mapcal#mapcmp [nh]
  * exec $PER/s#vpl prob0 dprob np dprob
  atitle 'n?beam!' 'prob'
  exec save mhad2011_prob0_[nh].eps f
  set pmci 1
  set mtyp 20
  hi/pl 10000
  set pmci 1
  set mtyp 24
  hi/pl 50000 s
  atitle 'z?r!, cm' 'A, pe.'
  exec save mhad2011_amplitude_vs_zr_counter[nh]_best.eps f
enddo
do i=[n1],[n2]
  nh=3100+[i]
  exec ../MinuitTest/mapcal#mapcmp [nh] dl=1
  exec ../MinuitTest/mapcal#mapcmp [nh]
  exec ../MinuitTest/mapcal#mapcmp [nh]
  * exec $PER/s#vpl prob0 dprob np dprob
  atitle 'n?beam!' 'prob'
  exec save mhad2011_prob0_[nh].eps f
  set pmci 1
  set mtyp 20
  hi/pl 10000
  set pmci 1
  set mtyp 24
  hi/pl 50000 s
  atitle 'z?r!, cm' 'A, pe.'
  exec save mhad2011_amplitude_vs_zr_counter[nh]_best.eps f
enddo
do i=[n1],[n2]
  nh=200+[i]
  exec ../MinuitTest/mapcal#mapcmp [nh] dl=1
  exec ../MinuitTest/mapcal#mapcmp [nh]
  exec ../MinuitTest/mapcal#mapcmp [nh]
  * exec $PER/s#vpl prob0 dprob np dprob
  atitle 'n?beam!' 'prob'
  exec save mhad2011_prob0_[nh].eps f
  set pmci 1
  set mtyp 20
  hi/pl 10000
  set pmci 1
  set mtyp 24
  hi/pl 50000 s
  atitle '[f], degree' 'A, pe.'
  exec save mhad2011_amplitude_vs_phi_counter[nh]_best.eps f
enddo
return


macro mapcmpxx n1=1 n2=9
set lwid 1
do i=[n1],-[n2]
  nh=1100+[i]
  fname=mhad2011_[nh]_v1.txt
  ve/del prob1,probd1,prob01
  ve/read prob1,probd1,prob01 [fname] '3f10.6'
  nh=2100+[i]
  fname=mhad2011_[nh]_v1.txt
  ve/del prob2,probd2,prob02
  ve/read prob2,probd2,prob02 [fname] '3f10.6'
  nh=3100+[i]
  fname=mhad2011_[nh]_v1.txt
  ve/del prob3,probd3,prob03
  ve/read prob3,probd3,prob03 [fname] '3f10.6'
  nh=200+[i]
  fname=mhad2011_[nh]_v1.txt
  ve/del prob4,probd4,prob04
  ve/read prob4,probd4,prob04 [fname] '3f10.6'
  nh=[i]
  ve/del prob,probd,prob0
  np=$vdim(prob01)
  ve/cre prob([np]) r
  ve/cre probd([np]) r
  ve/cre prob0([np]) r
  lvl=0.9
  do j=1,[np]
    if (($sigma(prob01([j])).gt.[lvl]).and.($sigma(prob02([j])).gt.[lvl]).and.($sigma(prob03([j])).gt.[lvl]).and.($sigma(prob04([j])).gt.[lvl])) then
      ve/inp prob0([j]) 1
    else
      ve/inp prob0([j]) 0.001
    endif
  enddo
  fname=mhad2011_[nh].txt
  ve/write prob,probd,prob0 [fname] '3f10.6'
  fname=mhad2011_[nh]_v1.txt
  ve/write prob,probd,prob0 [fname] '3f10.6'
enddo    
set pmci 1
do i=[n1],-[n2]
  nh=1100+[i]
  exec ../MinuitTest/mapcal#mapcmp [nh] [i]
  * exec $PER/s#vpl prob0 dprob np dprob iatt=24
  * exec $PER/s#vpl iprb dprob np dprob o=s iatt=20
  atitle 'n?beam!' 'prob'
  exec save mhad2011_prob0_[nh]_all.eps f
  set pmci 1
  set mtyp 20
  hi/pl 10000
  set pmci 1
  set mtyp 24
  hi/pl 50000 s
  atitle 'z?r!, cm' 'A, pe.'
  exec save mhad2011_amplitude_vs_zr_counter[nh]_best_all.eps f
enddo
do i=[n1],-[n2]
  nh=2100+[i]
  exec ../MinuitTest/mapcal#mapcmp [nh] [i]
  * exec $PER/s#vpl prob0 dprob np dprob iatt=24
  * exec $PER/s#vpl iprb dprob np dprob o=s iatt=20
  atitle 'n?beam!' 'prob'
  exec save mhad2011_prob0_[nh]_all.eps f
  set pmci 1
  set mtyp 20
  hi/pl 10000
  set pmci 1
  set mtyp 24
  hi/pl 50000 s
  atitle 'z?r!, cm' 'A, pe.'
  exec save mhad2011_amplitude_vs_zr_counter[nh]_best_all.eps f
enddo
do i=[n1],-[n2]
  nh=3100+[i]
  exec ../MinuitTest/mapcal#mapcmp [nh] [i]
  * exec $PER/s#vpl prob0 dprob np dprob iatt=24
  * exec $PER/s#vpl iprb dprob np dprob o=s iatt=20
  atitle 'n?beam!' 'prob'
  exec save mhad2011_prob0_[nh]_all.eps f
  set pmci 1
  set mtyp 20
  hi/pl 10000
  set pmci 1
  set mtyp 24
  hi/pl 50000 s
  atitle 'z?r!, cm' 'A, pe.'
  exec save mhad2011_amplitude_vs_zr_counter[nh]_best_all.eps f
enddo
do i=[n1],[n2]
  nh=200+[i]
  exec ../MinuitTest/mapcal#mapcmp [nh] [i]
  * exec $PER/s#vpl prob0 dprob np dprob iatt=24
  * exec $PER/s#vpl iprb dprob np dprob o=s iatt=20
  atitle 'n?beam!' 'prob'
  exec save mhad2011_prob0_[nh]_all.eps f
  set pmci 1
  set mtyp 20
  hi/pl 10000
  set pmci 1
  set mtyp 24
  hi/pl 50000 s
  atitle '[f], degree' 'A, pe.'
  exec save mhad2011_amplitude_vs_phi_counter[nh]_best_all.eps f
enddo
return

macro forprof
n=0
n=[n]+1; ebeam[n]=750   ; il[n]=694.7000120   ; nfirst[n]=7842  ; nlast[n]=9258 ;
n1=[n]
n=[n]+1; ebeam[n]=525   ; il[n]=363.2         ; nfirst[n]=8368  ; nlast[n]=8537 ;
n=[n]+1; ebeam[n]=550   ; il[n]=485.7         ; nfirst[n]=8538  ; nlast[n]=8643 ;
n=[n]+1; ebeam[n]=575   ; il[n]=419.1         ; nfirst[n]=8650  ; nlast[n]=8712 ;
n=[n]+1; ebeam[n]=600   ; il[n]=489.8999940   ; nfirst[n]=8715  ; nlast[n]=8768 ;
n=[n]+1; ebeam[n]=625   ; il[n]=440.5000000   ; nfirst[n]=8797  ; nlast[n]=8973 ;
n=[n]+1; ebeam[n]=650   ; il[n]=456.1000060   ; nfirst[n]=8978  ; nlast[n]=9032 ;
n=[n]+1; ebeam[n]=675   ; il[n]=559.7999880   ; nfirst[n]=9033  ; nlast[n]=9097 ;
n=[n]+1; ebeam[n]=700   ; il[n]=577.0999760   ; nfirst[n]=9100  ; nlast[n]=9152 ;
n=[n]+1; ebeam[n]=725   ; il[n]=432.2999880   ; nfirst[n]=9153  ; nlast[n]=9193 ;
n=[n]+1; ebeam[n]=775   ; il[n]=546.5000000   ; nfirst[n]=9261  ; nlast[n]=9320 ;
n=[n]+1; ebeam[n]=800   ; il[n]=439.3999940   ; nfirst[n]=9321  ; nlast[n]=9366 ;
n=[n]+1; ebeam[n]=825   ; il[n]=475.6000060   ; nfirst[n]=9367  ; nlast[n]=9435 ;
n=[n]+1; ebeam[n]=850   ; il[n]=462.2000120   ; nfirst[n]=9443  ; nlast[n]=9489 ;
n=[n]+1; ebeam[n]=875   ; il[n]=506.3999940   ; nfirst[n]=9492  ; nlast[n]=9550 ;
n=[n]+1; ebeam[n]=900   ; il[n]=384.1000060   ; nfirst[n]=9566  ; nlast[n]=9595 ;
n=[n]+1; ebeam[n]=925   ; il[n]=401.2999880   ; nfirst[n]=9598  ; nlast[n]=9651 ;
n=[n]+1; ebeam[n]=950   ; il[n]=451.2999880   ; nfirst[n]=9653  ; nlast[n]=9738 ;
n=[n]+1; ebeam[n]=962.5 ; il[n]=561.9000240   ; nfirst[n]=9739  ; nlast[n]=9825 ;
n=[n]+1; ebeam[n]=987.5 ; il[n]=467.5000000   ; nfirst[n]=9932  ; nlast[n]=10033 ;
n=[n]+1; ebeam[n]=1000  ; il[n]=538.0999760   ; nfirst[n]=10035 ; nlast[n]=10127 ;
n=[n]+1; ebeam[n]=945   ; il[n]=577.0999760   ; nfirst[n]=10162 ; nlast[n]=10192 ;
n=[n]+1; ebeam[n]=935   ; il[n]=637.5000000   ; nfirst[n]=10225 ; nlast[n]=10273 ;
n=[n]+1; ebeam[n]=912.5 ; il[n]=481.3999940   ; nfirst[n]=10275 ; nlast[n]=10295 ;
n=[n]+1; ebeam[n]=887.5 ; il[n]=525.7000120   ; nfirst[n]=10296 ; nlast[n]=10315 ;
n=[n]+1; ebeam[n]=862.5 ; il[n]=504.6000060   ; nfirst[n]=10316 ; nlast[n]=10332 ;
n=[n]+1; ebeam[n]=812.5 ; il[n]=509.8999940   ; nfirst[n]=10356 ; nlast[n]=10383 ;
n=[n]+1; ebeam[n]=787.5 ; il[n]=502.3999940   ; nfirst[n]=10384 ; nlast[n]=10424 ;
n=[n]+1; ebeam[n]=712.5 ; il[n]=593.2000120   ; nfirst[n]=10429 ; nlast[n]=10466 ;
n=[n]+1; ebeam[n]=737.5 ; il[n]=600.0999760   ; nfirst[n]=10467 ; nlast[n]=10513 ;
n=[n]+1; ebeam[n]=762.5 ; il[n]=487.2999880   ; nfirst[n]=10514 ; nlast[n]=10555 ;
n=[n]+1; ebeam[n]=687.5 ; il[n]=576.2000120   ; nfirst[n]=10559 ; nlast[n]=10614 ;
n=[n]+1; ebeam[n]=662.5 ; il[n]=526.2999880   ; nfirst[n]=10615 ; nlast[n]=10673 ;
n=[n]+1; ebeam[n]=637.5 ; il[n]=496.2000120   ; nfirst[n]=10703 ; nlast[n]=10745 ;
n=[n]+1; ebeam[n]=612.5 ; il[n]=546.7000120   ; nfirst[n]=10746 ; nlast[n]=10784 ;
n=[n]+1; ebeam[n]=587.5 ; il[n]=534.4         ; nfirst[n]=10785 ; nlast[n]=10825 ;
n=[n]+1; ebeam[n]=562.5 ; il[n]=517.3         ; nfirst[n]=10826 ; nlast[n]=10863 ;
n=[n]+1; ebeam[n]=537.5 ; il[n]=546.1         ; nfirst[n]=10864 ; nlast[n]=10902 ;
n=[n]+1; ebeam[n]=508.7 ; il[n]=0             ; nfirst[n]=10921 ; nlast[n]=10932 ;
n=[n]+1; ebeam[n]=509.8 ; il[n]=0             ; nfirst[n]=10934 ; nlast[n]=10952 ;
n=[n]+1; ebeam[n]=510.5 ; il[n]=0             ; nfirst[n]=10955 ; nlast[n]=10988 ;
n2=[n]

do i=1,9
  ve/del prob,probd
  ve/read prob,probd,prob[i] mhad2011_[i].txt '3f10.6'
enddo

file=mhad2011_touse.txt

if ($fexist([file]).ne.0) then
  shell rm [file]
endif
for/file  20 [file] new
close 20

do i=[n1],[n2]
  ebeam=[ebeam[i]]
  exec ../SepPar/sp#hists [ebeam]
  gl/imp ntup
  txt=' '
  do j=1,9
    if ($sigma(int(prob[j]([i]))).eq.1) then
      txt=$unquote([txt])$unquote('1 ')
    else
      txt=$unquote([txt])$unquote('0 ')
    endif
  enddo
  do q=1,[ntup]
    fmess [txt] [file]
  enddo
enddo

file=mhad2011_profile.list

if ($fexist([file]).ne.0) then
  shell rm [file]
endif
for/file  20 [file] new
close 20

gl/imp nfiles
dir=/work/snd2000/users/konctbel/exp/MHAD2011-2/col/

fmess [nfiles] [file]
fmess [dir] [file]
fmess [dir] [file]

nfiles=0
do i=[n1],[n2]
  ebeam=[ebeam[i]]
  exec ../SepPar/sp#hists [ebeam]
  gl/imp ntup
  txt=' '
  do q=1,[ntup]
    nfiles=[nfiles]+1
    exec ../SepPar/sp#hists [ebeam] n=[q]
    gl/imp fexpi
    fmess [fexpi] [file]
  enddo
enddo
gl/cre nfiles [nfiles]

return


*=======================================================================


macro puzpl ncnt=1

  npar=25
  fname=mhad2011_amplitude_vs_phi_counter[ncnt]__p1_p41_v2
  ve/read pars,dpars [fname]_fit.par '2f15.10'
  ve/del fitpad
  ve/copy pars(1:19) fitpad
  do i=1,11
    ve/inp fitpad([i]) $sigma(fitpad([i])/(14-10.6))
  enddo
  s1=fitpad(12)
  ss=fitpad(13)
  sm=fitpad(14)
  s2=fitpad(15)
  ve/cre phii(10) r [s1] $sigma([s1]+([ss]-[s1])/3) $sigma([s1]+([ss]-[s1])*2/3) _
                    [ss] $sigma([ss]+([sm]-[ss])/3) $sigma([ss]+([sm]-[ss])*2/3) _
                    [sm] $sigma([sm]+([s2]-[sm])/3) $sigma([sm]+([s2]-[sm])*2/3) [s2]
  ve/cre vpars(186) r
  do i=1,14
    i1=([i]-1)*11+1
    i2=[i]*11
    ve/copy fitpad(1:11) vpars([i1]:[i2])
  enddo
  ve/inp vpars(155) [s1]
  ve/inp vpars(156) [ss]
  ve/inp vpars(157) [sm]
  ve/inp vpars(158) [s2]
  ve/inp vpars(161) $sigma(fitpad(17))
  ve/inp vpars(162) $sigma(fitpad(18))
  ve/inp vpars(163) 0
  ve/inp vpars(164) 0
  ve/inp vpars(165) 2

  ve/cre  vxi(15,3) r
  ve/cre vdxi(13,3) r 
  ve/cre vxmi(14,3) r
  ve/cre valphai(14,3) r
  ve/cre vbeta(1) r
  ve/cre vxsi(15) r

  sm=0
  z1=0
  z2=0
  do il=3,1,-1
  
    nh=110[ncnt]

    fname=mhad2011_amplitude_vs_zr_counter[nh]_[sm]

    shell /work/snd2000/users/konctbel/SNDSimACC/readlifstat [fname].log
    ve/del pars0,dpars0
    ve/read pars0,dpars0 [fname].log.pars
    ve/del parmin
    ve/copy pars0(2:83) parmin

    ve/copy  pars0(2:16)  vxi(1:15,[il])
    ve/copy  pars0(44:57)  vxmi(1:14,[il])
    ve/copy  pars0(58:81)  valphai(1:14,[il])
    
    ve/copy pars0(78) vbeta(1)
    z1=$sigma([z1]+pars0(2))
    z2=$sigma([z2]+pars0(16))
    
  enddo
  ve/inp vpars(159) $sigma([z1]/3)
  ve/inp vpars(160) $sigma([z2]/3)
  xmin=-13
  xmax=12
  do i=1,15
    ii=165+[i]
    ve/inp vpars([ii]) 10
    ve/inp vxsi([i]) $sigma([xmin]+([i]-1)*([xmax]-([xmin]))/14)
  enddo
  ve/inp vpars(181) $sigma([z1]/3-0.5)
  ve/inp vpars(182) $sigma([z2]/3+0.5)
  ve/inp vpars(183) [ss]
  ve/inp vpars(184) 0.9
  ve/inp vpars(185) 20
  ve/inp vpars(186) 5
  
  ve/cre pns(1) r 14
                    
  phi=phii(3); fun/pl puzfitp.f(x,[phi],0.0) -15 15
  phi=phii(2); fun/pl puzfitp.f(x,[phi],0.0) -15 15 s
  phi=phii(4); fun/pl puzfitp.f(x,[phi],0.0) -15 15 s
  phi=phii(1); fun/pl puzfitp.f(x,[phi],0.0) -15 15 s
return





macro prepnew p=parn fout=mapfit.out.old ncnt=1 ebeam=510.2 exp=mhad2011 edata=no ver=0

gl/cre nxg 12
gl/cre nyg 12

chain -exp
chain exp mhad2011_510.5-2_col_p1.hbook
chain exp mhad2011_510.5-2_col_p2.hbook
chain exp mhad2011_510.5-2_col_p3.hbook
chain exp mhad2011_510.5-2_col_p4.hbook
chain exp -P /work/users/konctbel/exp/MHAD2011-2/col/

*chain -exp
*chain exp exp1_0_510.5.hbook
*chain exp exp2_0_510.5.hbook
*chain exp exp3_0_510.5.hbook
*chain exp exp4_0_510.5.hbook

mess shell ./mkcal.sh [exp] [ebeam] [ncnt]
shell ./mkcal.sh [exp] [ebeam] [ncnt]
dir=[exp]_[ebeam]/counter_[ncnt]
datafile=[exp]_[ebeam]_nc[ncnt].dat
datafile=mapfitnew.txt

gl/cre sc $sigma(180/3.1415927)
r1=10.666
r2=13.834
r=([r1]+[r2])/2

*-----------------------------

mnz=300
mzmin=-15
mzmax=15
mnf=320
mfmax=[ncnt]*40+20
mfmin=[ncnt]*40-60
dz=$sigma(([mzmax]-[mzmin])/[mnz])
df=$sigma(([mfmax]-[mfmin])/[mnf])
ve/cre xvi([mnz]) r
ve/cre yvi([mnf]) r
l0=[mzmin]+[dz]/2.
r0=[mzmax]-[dz]/2.
sigma xvi = array([mnz],[l0]#[r0])
l0=[mfmin]+[df]/2.
r0=[mfmax]-[df]/2.
sigma yvi = array([mnf],[l0]#[r0])

if ([edata].eq.'yes') then

hi/file 20 mhad2011_profile_sin.his

ve/cre  avi([mnz]) r
ve/cre davi([mnz]) r
ve/cre  av([mnz],[mnf]) r
ve/cre dav([mnz],[mnf]) r

do i=1,[mnf]

  idh=[i]*100+[ncnt]
  hi/pl //lun20/[idh]
  
  hi/get/cont //lun20/[idh]  avi
  hi/get/err  //lun20/[idh] davi
	
  ve/copy  avi(1:[mnz])  av(1:[mnz],[i])
  ve/copy davi(1:[mnz]) dav(1:[mnz],[i])

enddo
ve/del avi,davi
close 20

ve/write av,dav [dir]/[datafile] 2g15.6

endif

nsc=[mnz]*[mnf]
ve/del av1,dav1
ve/read av1,dav1 [dir]/[datafile] 2g15.6

ve/cre  av([mnz],[mnf]) r
ve/cre dav([mnz],[mnf]) r
ve/read av,dav [dir]/[datafile] 2g15.6


*------------------

r1=10.666
r2=13.834
r=([r1]+[r2])/2
d1=0.5
d2=1
s0=0
mode=1
md1=1
md2=1
pk=3
pn=180
pm=40
zk=3
zn=100
zm=20
kx=3
ky=3
nx=11
ny=11
nx0=[nx]-2
ny0=[ny]-2
nxy0=[nx0]*[ny0]
npar0=21
ind0=[npar0]+[nx0]
nx0st=1
nx0sp=[nx0]
npar=[npar0]+[nx0]+[nx0]*[ny0]

*------------------

ve/cre s1f(9) r 3.9 43.5 82.6 122.6 163.4 203.9 245.2 285.4 324.7
ve/cre s2f(9) r 42.6 81.4 120.2 162.6 202.1 242.5 284.1 323.6 362.1
ve/cre ssf(9) r 17.6 57.6 96.3 136.3 178.0 217.7 259.1 299.5 338.8

ds0=2
*exec mapcal#s1s2 $sigma([ds0]+4.2+40*([ncnt]-1)) $sigma([ds0]+42.2+40*([ncnt]-1)) $sigma([ds0]+17.5+40*([ncnt]-1))
exec mapcal#s1s2 $sigma(s1f([ncnt])) $sigma(s2f([ncnt])) $sigma(ssf([ncnt]))
gl/imp s1g
gl/imp s2g
gl/imp ssg
gl/imp sigs1g
gl/imp sigs2g
gl/imp sigssg
exec mapcal#z1z2 $sigma(ssf([ncnt]))
gl/imp z1g
gl/imp z2g
*z1g=-10
*z2g=10

*------------------

z0=0; dz0=0.1
z1=[z1g]; dz1=0.1
z2=[z2g];  dz2=0.1
s1=[s1g];  ds1=[sigs1g]/10
s2=[s2g];  ds2=[sigs2g]/10
sigz=0.6; dsigz=0.001
sigs=([sigs1g]+[sigs2g])/2; dsigs=([sigs1g]+[sigs2g])/10
sigs0=[sigssg]; dsigs0=0.01
pds=0;    dpds=0.001
rPMT=0.86; drPMT=0.001
sPMT=[ssg]; dsPMT=0.01
apmt=14; daPMT=0.1
zs1=[z1g]; dzs1=0.1
zs2=[z2g]; dzs2=0.1
dzds=0; ddzds=0.00001
ss=[ssg]; dss=0.01
alpha=0; dalpha=0.001
ps1=0.05; dps1=0.01
ps2=0.05; dps2=0.01
ts1=40; dts1=1
ts2=40; dts2=1
ve/cre xi([nx]) r
ve/cre yi([ny]) r
l0=[z1]
r0=[z2]
sigma xi = array([nx],[l0]#[r0])
l0=[s1]
r0=[s2]
sigma yi = array([ny],[l0]#[r0])

ve/cre wls([nx0]) r [nx0]*15
ve/cre dwls([nx0]) r [nx0]*0.3
ve/cre map([nx0],[ny0]) r [nxy0]*2
ve/cre dmap([nx0],[ny0]) r [nxy0]*0.1
exec mapcal#mapc [ssg]
sigma  map =  mape/([r2]-[r1])
sigma dmap = dmape/([r2]-[r1])

ve/cre parn([npar]) r [z1] [z2] [s1] [s2] [sigz] [sigs] [pds] [rpmt] [apmt] [zs1] [zs2] [dzds] [ss] [z0] [sPMT] [sigs0] [ps1] [ps2] [ts1] [ts2]
i1=[ind0]-[nx0]+1
i2=[ind0]
ve/copy wls(1:[nx0]) parn([i1]:[i2])
do j=1,[ny0]
  i1=[i1]+[nx0]
  i2=[i2]+[nx0]
  ve/copy map(1:[nx0],[j]) parn([i1]:[i2])
enddo


if ([p].eq.'parpl') then

  n=0
  n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out
  n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out
  n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out
  n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out
  n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out.cut
  n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out.cut
  n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out.cut
  n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out.cut
  n=[n]+1; file[n]=mhad2011_510.5/counter_[n]/mapfit.out
  fout=[file[ncnt]]
  
*  shell /work/snd2000/users/konctbel/SNDSimACC/readlifstat_new mapfit.out
  shell /work/snd2000/users/konctbel/SNDSimACC/readlifstat [fout]
  ve/del pars0,dpars0
  ve/read pars0,dpars0 [fout].pars

  fout=[dir]/mapfit.out
  dir=.
  ve/del xi,yi
  shell fgrep -e "map0,xi:" [fout] >& tmp.txt
  shell $unquote('cat tmp.txt | sed "s/map0,xi:/ /g" > tmp1.txt')
  ve/read xi tmp1.txt
  shell fgrep -e "map0,yi:" [fout] >& tmp.txt
  shell $unquote('cat tmp.txt | sed "s/map0,yi:/ /g" > tmp1.txt')
  ve/read yi tmp1.txt

  ve/inp xi(1) $sigma(xi(1) - 1000)
  gl/cre nxg $vlen(xi)
  gl/cre nyg $vlen(yi)
  nx=[nxg]
  ny=[nyg]
  nx0=[nx]-2
  ny0=[ny]-2
  ind0=[npar0]+[nx0]
  nx0st=1
  nx0sp=[nx0]
  npar=[npar0]+[nx0]+[nx0]*[ny0]
  ve/cre wls([nx0]) r [nx0]*5
  ve/cre dwls([nx0]) r [nx0]*0.3
  nxy0=[nx0]*[ny0]
  ve/cre map([nx0],[ny0]) r [nxy0]*2
  ve/cre dmap([nx0],[ny0]) r [nxy0]*0.1

  ve/cre parn([npar]) r

  i1=2
  i2=[npar]+1
  ve/copy pars0([i1]:[i2]) parn(1:[npar])
  i=1
  i=[i]+1; z1=pars0([i]);    dz1=dpars0([i])
  i=[i]+1; z2=pars0([i]);    dz2=dpars0([i])
  i=[i]+1; s1=pars0([i]);    ds1=dpars0([i])
  i=[i]+1; s2=pars0([i]);    ds2=dpars0([i])
  i=[i]+1; sigz=pars0([i]);  dsigz=dpars0([i])
  i=[i]+1; sigs=pars0([i]);  dsigs=dpars0([i])
  i=[i]+1; pds=pars0([i]);   dpds=dpars0([i])
  i=[i]+1; rPMT=pars0([i]);  drPMT=dpars0([i])
  i=[i]+1; apmt=pars0([i]);  daPMT=dpars0([i])
  i=[i]+1; zs1=pars0([i]);   dzs1=dpars0([i])
  i=[i]+1; zs2=pars0([i]);   dzs2=dpars0([i])
  i=[i]+1; dzds=pars0([i]);  ddzds=dpars0([i])
  i=[i]+1; ss=pars0([i]);    dss=dpars0([i])
  i=[i]+1; alpha=pars0([i]); dalpha=dpars0([i])
  i=[i]+1; z0=pars0([i]);    dz0=dpars0([i])
  i=[i]+1; sPMT=pars0([i]);    dsPMT=dpars0([i])
  i=[i]+1; sigs0=pars0([i]);  dsigs0=dpars0([i])
  i=[i]+1; ps1=pars0([i]);  dps1=dpars0([i])
  i=[i]+1; ps2=pars0([i]);  dps2=dpars0([i])
  i=[i]+1; ts1=pars0([i]);  dts1=dpars0([i])
  i=[i]+1; ts2=pars0([i]);  dts2=dpars0([i])
  i1=[i]+1;
  i2=[i1]+[nx0]-1
  ve/copy  pars0([i1]:[i2])  wls(1:[nx0])
*  ve/copy dpars0([i1]:[i2]) dwls(1:[nx0])
  do i=1,[ny0]
    i1=[i2]+1
    i2=[i1]+[nx0]-1
    ve/copy  pars0([i1]:[i2])  map(1:[nx0],[i])
*    ve/copy dpars0([i1]:[i2]) dmap(1:[nx0],[i])
  enddo
endif

if ([p].eq.'pars0') then
*  shell /work/snd2000/users/konctbel/SNDSimACC/readlifstat_new mapfit.out
  shell /work/snd2000/users/konctbel/SNDSimACC/readlifstat [fout]
  ve/del pars0,dpars0
  ve/read pars0,dpars0 [fout].pars

  xx=4
  if ([xx].eq.4) then
  gl/imp nxg
  gl/imp nyg
  nx=[nxg]
  ny=[nyg]
  nx0=[nx]-2
  ny0=[ny]-2
  ind0=[npar0]+[nx0]
  nx0st=1
  nx0sp=[nx0]
  npar=[npar0]+[nx0]+[nx0]*[ny0]
  ve/cre wls([nx0]) r [nx0]*5
  ve/cre dwls([nx0]) r [nx0]*0.3
  nxy0=[nx0]*[ny0]
  ve/cre map([nx0],[ny0]) r [nxy0]*2
  ve/cre dmap([nx0],[ny0]) r [nxy0]*0.1
  ve/cre xi([nx]) r
  ve/cre yi([ny]) r
  sigma xi = array([nx],-12#12)
  l0=5+40*([ncnt]-1)
  r0=43+40*([ncnt]-1)
  sigma yi = array([ny],[l0]#[r0])
  ve/cre parn([npar]) r
  endif

  i1=2
  i2=[npar]+1
  ve/copy pars0([i1]:[i2]) parn(1:[npar])
  i=1
  i=[i]+1; z1=pars0([i]);    dz1=dpars0([i])
  i=[i]+1; z2=pars0([i]);    dz2=dpars0([i])
  i=[i]+1; s1=pars0([i]);    ds1=dpars0([i])
  i=[i]+1; s2=pars0([i]);    ds2=dpars0([i])
  i=[i]+1; sigz=pars0([i]);  dsigz=dpars0([i])
  i=[i]+1; sigs=pars0([i]);  dsigs=dpars0([i])
  i=[i]+1; pds=pars0([i]);   dpds=dpars0([i])
  i=[i]+1; rPMT=pars0([i]);  drPMT=dpars0([i])
  i=[i]+1; apmt=pars0([i]);  daPMT=dpars0([i])
  i=[i]+1; zs1=pars0([i]);   dzs1=dpars0([i])
  i=[i]+1; zs2=pars0([i]);   dzs2=dpars0([i])
  i=[i]+1; dzds=pars0([i]);  ddzds=dpars0([i])
  i=[i]+1; ss=pars0([i]);    dss=dpars0([i])
  i=[i]+1; alpha=pars0([i]); dalpha=dpars0([i])
  i=[i]+1; z0=pars0([i]);    dz0=dpars0([i])
  i=[i]+1; sPMT=pars0([i]);    dsPMT=dpars0([i])
  i=[i]+1; sigs0=pars0([i]);  dsigs0=dpars0([i])
  i=[i]+1; ps1=pars0([i]);  dps1=dpars0([i])
  i=[i]+1; ps2=pars0([i]);  dps2=dpars0([i])
  i=[i]+1; ts1=pars0([i]);  dts1=dpars0([i])
  i=[i]+1; ts2=pars0([i]);  dts2=dpars0([i])
  i1=[i]+1;
  i2=[i1]+[nx0]-1
  ve/copy  pars0([i1]:[i2])  wls(1:[nx0])
*  ve/copy dpars0([i1]:[i2]) dwls(1:[nx0])
  do i=1,[ny0]
    i1=[i2]+1
    i2=[i1]+[nx0]-1
    ve/copy  pars0([i1]:[i2])  map(1:[nx0],[i])
*    ve/copy dpars0([i1]:[i2]) dmap(1:[nx0],[i])
  enddo
*  exec mapcal#mapc [ssg]
*  sigma  map =  mape/([r2]-[r1])
*  sigma dmap = dmape/([r2]-[r1])
*  s1=[s1g];  ds1=[sigs1g]/10
*  s2=[s2g];  ds2=[sigs2g]/10
*  sigs=([sigs1g]+[sigs2g])/2; dsigs=([sigs1g]+[sigs2g])/10
*  alpha=0; dalpha=0.0001
*  dzds=0; ddzds=0.0001
*  z1=[z1g];
*  z2=[z2g];
*
  i2=[ind0]
  do i=1,-[ny0]
    i1=[i2]+1
    i2=[i1]+[nx0]-1
*    ve/copy  map(1:[nx0],[i]) parn([i1]:[i2])
  enddo
  
  
*new  

  k0=1
  npar=25
  fname=mhad2011_amplitude_vs_phi_counter[ncnt]__p1_p41_v2
  ve/read pars,dpars [fname]_fit.par '2f15.10'
  ve/del fitpad
  ve/copy pars(1:19) fitpad
  ve/copy dpars(1:19) dfitpad
  do i=1,11
    ve/inp fitpad([i]) $sigma(fitpad([i])/(14-10.6))
    ve/inp dfitpad([i]) $sigma([k0]*dfitpad([i])/(14-10.6))
  enddo
  sig1=fitpad(17); dsig1=dfitpad(17)
  sig2=fitpad(18); dsig2=dfitpad(18)
  s1=fitpad(12); ds1=dfitpad(12)
  ss=fitpad(13); dss=dfitpad(13)
  sm=fitpad(14); dsm=dfitpad(14)
  s2=fitpad(15); ds2=dfitpad(15)
  ve/cre phii(10) r [s1] $sigma([s1]+([ss]-[s1])/3) $sigma([s1]+([ss]-[s1])*2/3) _
                    [ss] $sigma([ss]+([sm]-[ss])/3) $sigma([ss]+([sm]-[ss])*2/3) _
                    [sm] $sigma([sm]+([s2]-[sm])/3) $sigma([sm]+([s2]-[sm])*2/3) [s2]
  ve/cre vyi(12) r [s1] $sigma([s1]+([ss]-[s1])/3) $sigma([s1]+([ss]-[s1])*2/3) [ss] _
                   [ss] $sigma([ss]+([sm]-[ss])/3) $sigma([ss]+([sm]-[ss])*2/3) [sm] _
                   [sm] $sigma([sm]+([s2]-[sm])/3) $sigma([sm]+([s2]-[sm])*2/3) [s2]
  ve/cre  parn(186) r
  ve/cre dparn(186) r
  do i=1,14
    i1=([i]-1)*11+1
    i2=[i]*11
    ve/copy  fitpad(1:11)  parn([i1]:[i2])
    ve/copy dfitpad(1:11) dparn([i1]:[i2])
  enddo
  ve/inp parn(155) [s1]
  ve/inp parn(156) [ss]
  ve/inp parn(157) [sm]
  ve/inp parn(158) [s2]
  ve/inp parn(161) $sigma(fitpad(17))
  ve/inp parn(162) $sigma(fitpad(18))
  ve/inp parn(163) 0
  ve/inp parn(164) 0
  ve/inp parn(165) 2

  ve/cre  vxi(45) r
  ve/cre vdxi(39) r 
  ve/cre vxmi(42) r
  ve/cre valphai(42) r
  ve/cre vbeta(1) r
  ve/cre vxsi(15) r

  smd=0
  z1=0
  z2=0
  do il=3,1,-1
  
    nh=[il]10[ncnt]

    fname=mhad2011_amplitude_vs_zr_counter[nh]_[smd]

    shell /work/snd2000/users/konctbel/SNDSimACC/readlifstat [fname].log
    ve/del pars0,dpars0
    ve/read pars0,dpars0 [fname].log.pars
    ve/del parmin
    ve/copy pars0(2:83) parmin

    i1=15*[il]-14; l2=15*[il]
    ve/copy  pars0(2:16)  vxi([i1]:[i2])
    i1=14*[il]-13; l2=14*[il]
    ve/copy  pars0(44:57)  vxmi([i1]:[i2])
    ve/copy  pars0(58:81)  valphai([i1]:[i2])
    
    ve/copy pars0(78) vbeta(1)
    z1=$sigma([z1]+pars0(2))
    z2=$sigma([z2]+pars0(16))
    
  enddo
  ve/copy valphai(1:14) valphai(15:28)
  ve/copy valphai(1:14) valphai(29:42)
  ve/inp parn(159) $sigma([z1]/3)
  ve/inp parn(160) $sigma([z2]/3)
  z1=parn(159); dz1=0.5*[k0]
  z2=parn(160); dz2=0.5*[k0]
  zs1=parn(159); dzs1=0.5*[k0]
  zs2=parn(160); dzs2=0.5*[k0]
  z0=0; dz0=0
  sigz0=0; dsigz0=0
  xmin=-13
  xmax=12
  ve/cre  asi(15) r 15*10
  ve/cre dasi(15) r 15*$sigma(0.5*[k0])
  do i=1,15
    ii=165+[i]
    ve/inp parn([ii]) 10
    ve/inp vxsi([i]) $sigma([xmin]+([i]-1)*([xmax]-([xmin]))/14)
  enddo
  ve/inp parn(181) $sigma([z1]-0.5)
  ve/inp parn(182) $sigma([z2]+0.5)
  ve/inp parn(183) [ss]
  ve/inp parn(184) 0.9
  ve/inp parn(185) 20
  ve/inp parn(186) 5
  pds=0; dpds=0.02*[k0]
  phiPMT=[ss]; dphiPMT=[dss]*[k0]
  rPMT=0.8; drPMT=0.01*[k0]
  aPMT=10; daPMT=1*[k0]
  a2PMT=0; da2PMT=0.1*[k0]
  ve/inp parn(164) [z0]
  ve/inp parn(165) [sigz0]
  ve/inp parn(183) [phiPMT]
  ve/inp parn(184) [rPMT]
  ve/inp parn(185) [aPMT]
  ve/inp parn(186) [a2PMT]
  
  nb=14
  ve/cre pns(1) r [nb]
  
  fname=[dir]/mapfitnew.out.v0
  shell /work/snd2000/users/konctbel/SNDSimACC/readlifstat [fname]
  ve/del pars0,dpars0
  ve/read pars0,dpars0 [fname].pars
  ve/del parn,dparn
  ve/copy  pars0(2:187)  parn
  ve/copy dpars0(2:187) dparn
  
  s1=parn(155);    ds1=dparn(155)
  ss=parn(156);    dss=dparn(156)
  sm=parn(157);    dsm=dparn(157)
  s2=parn(158);    ds2=dparn(158)
  z1=parn(159);    dz1=dparn(159)
  z2=parn(160);    dz2=dparn(160)
  sig1=parn(161);  dsig1=dparn(161)
  sig2=parn(162);  dsig2=dparn(162)
  pds=parn(163);   dpds=dparn(163)
  z0=parn(164);    dz0=dparn(164)
  sigz0=parn(165); dsigz0=dparn(165)
  ve/copy  parn(166:180)  asi(1:15)
  ve/copy dparn(166:180) dasi(1:15)
  zs1=parn(181);    dzs1=dparn(181)
  zs2=parn(182);    dzs2=dparn(182)
  phiPMT=parn(183);    dphiPMT=dparn(183)
  rPMT=parn(184);    drPMT=dparn(184)
  aPMT=parn(185);    daPMT=dparn(185)
  a2PMT=parn(186);    da2PMT=dparn(186)
  
  
*  ve/inp parn(184) 0.8
*  ve/inp parn(185) 8
*  ve/inp parn(186) 1
*  ve/inp parn(165) 2
                      
endif

if ([p].eq.'chmap') then
*  shell /work/snd2000/users/konctbel/SNDSimACC/readlifstat_new mapfit.out
  shell /work/snd2000/users/konctbel/SNDSimACC/readlifstat [fout]
  ve/del pars0,dpars0
  ve/read pars0,dpars0 [fout].pars

  i1=2
  i2=[npar]+1
  ve/copy pars0([i1]:[i2]) parn(1:[npar])
  i=1
  i=[i]+1; z1=pars0([i]);    dz1=dpars0([i])
  i=[i]+1; z2=pars0([i]);    dz2=dpars0([i])
  i=[i]+1; s1=pars0([i]);    ds1=dpars0([i])
  i=[i]+1; s2=pars0([i]);    ds2=dpars0([i])
  i=[i]+1; sigz=pars0([i]);  dsigz=dpars0([i])
  i=[i]+1; sigs=pars0([i]);  dsigs=dpars0([i])
  i=[i]+1; pds=pars0([i]);   dpds=dpars0([i])
  i=[i]+1; rPMT=pars0([i]);  drPMT=dpars0([i])
  i=[i]+1; apmt=pars0([i]);  daPMT=dpars0([i])
  i=[i]+1; zs1=pars0([i]);   dzs1=dpars0([i])
  i=[i]+1; zs2=pars0([i]);   dzs2=dpars0([i])
  i=[i]+1; dzds=pars0([i]);  ddzds=dpars0([i])
  i=[i]+1; ss=pars0([i]);    dss=dpars0([i])
  i=[i]+1; alpha=pars0([i]); dalpha=dpars0([i])
  i=[i]+1; z0=pars0([i]);    dz0=dpars0([i]); dz0=0.1
  i=[i]+1; sPMT=pars0([i]);    dsPMT=dpars0([i])
  i=[i]+1; sigs0=pars0([i]);  dsigs0=dpars0([i])
  i=[i]+1; ps1=pars0([i]);  dps1=dpars0([i])
  i=[i]+1; ps2=pars0([i]);  dps2=dpars0([i])
  i=[i]+1; ts1=pars0([i]);  dts1=dpars0([i])
  i=[i]+1; ts2=pars0([i]);  dts2=dpars0([i])
  i1=[i]+1;
  i2=[i1]+[nx0]-1
  ve/copy  pars0([i1]:[i2])  wls(1:[nx0])
  ve/copy dpars0([i1]:[i2]) dwls(1:[nx0])
  do i=1,[ny0]
    i1=[i2]+1
    i2=[i1]+[nx0]-1
    ve/copy  pars0([i1]:[i2])  map(1:[nx0],[i])
    ve/copy dpars0([i1]:[i2]) dmap(1:[nx0],[i])
  enddo
  gl/imp nxg
  gl/imp nyg
  nx=[nxg]
  ny=[nyg]
  nx0=[nx]-2
  ny0=[ny]-2
  ind0=[npar0]+[nx0]
  nx0st=1
  nx0sp=[nx0]
  npar=[npar0]+[nx0]+[nx0]*[ny0]
  ve/cre wls([nx0]) r [nx0]*5
  ve/cre dwls([nx0]) r [nx0]*0.3
  nxy0=[nx0]*[ny0]
  ve/cre map([nx0],[ny0]) r [nxy0]*2
  ve/cre dmap([nx0],[ny0]) r [nxy0]*0.1
  ve/cre xi([nx]) r
  ve/cre yi([ny]) r
  sigma xi = array([nx],-12#12)
  l0=5+40*([ncnt]-1)
  r0=43+40*([ncnt]-1)
  sigma yi = array([ny],[l0]#[r0])
  ve/cre xr(1) r
  ve/cre yr(1) r
  do i=1,[nx0]
    ve/inp xr(1) $sigma(xi([i]))
    ve/inp wls([i]) $call('accmappl.sl(xr,yr,2)')
    do j=1,[ny0]
      ve/inp yr(1) $sigma(yi([j]))
      ve/inp map([i],[j]) $call('accmappl.sl(xr,yr,1)')
    enddo
  enddo
  ve/cre parn([npar]) r [z1] [z2] [s1] [s2] [sigz] [sigs] [pds] [rpmt] [apmt] [zs1] [zs2] [dzds] [ss] [z0] [sPMT] [sigs0]
  i1=[ind0]-[nx0]+1
  i2=[ind0]
  ve/copy wls(1:[nx0]) parn([i1]:[i2])
  do j=1,[ny0]
    i1=[i1]+[nx0]
    i2=[i2]+[nx0]
    ve/copy map(1:[nx0],[j]) parn([i1]:[i2])
  enddo
endif

dir=.
fname=mapfitnew

file=[dir]/[fname].inc
if ($fexist([file]).ne.0) then
*  shell rm [file]
  shell mv -v [file] [file].old
endif
for/file  20 [file] new
close 20

nbm=[nb]-1
nbp=[nb]+1
*fmess '      common/pars/nb,xi(3,101),dxi(3,100),yi(3,4),ai(3,4,100),' [file]
txt='      common/pars/'
txt=$unquote([txt])xi(3,[nbp]),dxi(3,[nbm]),yi(3,4),ai(3,4,[nb]),
fmess [txt] [file]
*fmess '     &  xmi(3,100),alphai(3,100),beta,' [file]
txt='     &'
txt=$unquote([txt])  xmi(3,[nb]),alphai(3,[nb]),beta,
fmess [txt] [file]
fmess '     &  z0,sig0,sig1(3),sig2(3),' [file]
fmess '     &  pds,' [file]
fmess '     &  lambda(3), xr1(3), xr2(3), sr1(3), sr2(3),' [file]
*fmess '     &  recalcs, xsi(101), asi(101), xsmin, xsmax,' [file]
txt='     &'
txt=$unquote([txt])  xsi([nbp]), asi([nbp]), xsmin, xsmax,
fmess [txt] [file]
fmess '     &  zPMT, phiPMT, rPMT, aPMT, a2PMT,' [file]
fmess '     &  nb, recalcs' [file]
fmess '      integer nb' [file]
fmess '      double precision xi, dxi, yi, ai, xmi, alphai, beta,' [file]
fmess '     &  z0, sig0, sig1, sig2, pds,' [file]
fmess '     &  lambda, xr1, xr2, sr1, sr2, xsi, asi, xsmin, xsmax,' [file]
fmess '     &  zPMT, phiPMT, rPMT, aPMT, a2PMT' [file]
fmess '      logical recalcs' [file]
fmess '' [file]
txt='      data nb'
txt=$unquote([txt])/[nb]/
fmess [txt] [file]
fmess '' [file]
txt='      double precision '
txt=$unquote([txt])vxi([nbp],3),vdxi([nbm],3),vyi(12)
fmess [txt] [file]                                                                                                              
txt='      double precision '
txt=$unquote([txt])vxmi([nb],3),valphai([nb],3)
fmess [txt] [file]                                                                                                              

fmess '      data vxi /' [file]
exec mapcal#vewrite [file] vxi f10.4 $vdim(vxi) 4
fmess '     &  /' [file]
fmess '      data vdxi /' [file]
exec mapcal#vewrite [file] vdxi f10.4 $vdim(vdxi) 4
fmess '     &  /' [file]
fmess '      data vyi /' [file]
exec mapcal#vewrite [file] vyi f10.4 $vdim(vyi) 4
fmess '     &  /' [file]
fmess '      data vxmi /' [file]
exec mapcal#vewrite [file] vxmi f10.4 $vdim(vxmi) 4
fmess '     &  /' [file]
fmess '      data valphai /' [file]
exec mapcal#vewrite [file] valphai f10.4 $vdim(valphai) 4
fmess '     &  /' [file]
fmess '      data xsi /' [file]
exec mapcal#vewrite [file] vxsi f10.4 $vdim(vxsi) 4
fmess '     &  /' [file]
txt='      data beta'
txt=$unquote([txt])/$sigma(vbeta(1))/
fmess [txt] [file]
fmess '' [file]
txt='      double precision'
txt=$unquote([txt]) r1, r2
fmess [txt] [file]
txt='      data r1,r2'
txt=$unquote([txt])/[r1],[r2]/
fmess [txt] [file]
fmess '      integer kag,kwls,kpds,kpmt,ktail' [file]
fmess '      data kag,kwls,kpds,kpmt,ktail/1,1,1,1,1/' [file]
npar=$vdim(parn)
fmess '      integer npar0' [file]
txt='      data npar0 '
txt=$unquote([txt])/[npar]/
fmess [txt] [file]
txt='      common/parnc/'
txt=$unquote([txt])parn([npar])
fmess [txt] [file]
fmess '      double precision parn' [file]
fmess '      data parn /' [file]
exec mapcal#vewrite [file] parn f10.4 [npar] 4
fmess '     &  /' [file]
fmess '' [file]
fmess '' [file]
fmess '' [file]
fmess '' [file]

*---------------

file=[dir]/[fname]_data.inc
if ($fexist([file]).eq.1) then
*  shell rm -v [file]
  shell mv -v [file] [file].old
endif
for/file  20 [file] new
close 20

fmess '      integer nvx, nvy' [file]
txt='      data nvx,nvy'
txt=$unquote([txt])/[mnz],[mnf]/
fmess [txt] [file]

txt='      common/mapc/'
txt=$unquote([txt])av([mnz],[mnf]),dav([mnz],[mnf]),xvi([mnz]),yvi([mnf])
fmess [txt] [file]
fmess '      double precision av,dav,xvi,yvi' [file]

d=$sigma(([mzmax]-[mzmin])/[mnz]/2)
ve/cre xvi([mnz]) r
l=[mzmin]+[d]
r=[mzmax]-[d]
sigma xvi = array([mnz],[l]#[r])
fmess '      data xvi /' [file]
exec mapcal#vewrite [file] xvi f10.2 [mnz] 4
fmess '     &  /' [file]

d=$sigma(([mfmax]-[mfmin])/[mnf]/2)
ve/cre yvi([mnf]) r
l=[mfmin]+[d]
r=[mfmax]-[d]
sigma yvi = array([mnf],[l]#[r])
fmess '      data yvi /' [file]
exec mapcal#vewrite [file] yvi f10.2 [mnf] 4
fmess '     &  /' [file]

fmess '      data av /' [file]
exec mapcal#vewrite [file] av1 g10.3 [nsc] 4
fmess '     &  /' [file]
fmess '      data dav /' [file]
exec mapcal#vewrite [file] dav1 g10.3 [nsc] 4
fmess '     &  /' [file]

*---------------

file=[dir]/[fname].dat
if ($fexist([file]).eq.1) then
  shell rm -v [file]
endif
for/file  20 [file] new
close 20

fmess 'SET TITLE' [file]
ni=0
fmess 'counter ununiformity fitting' [file]
fmess 'PARAMETERS' [file]
nif=[ni]
do i=1,[nb]
  do j=1,11
    tmp0=a_[i]_[j]
    ni=[ni]+1
    aij=$sigma(parn([ni]))
    aijl=0
    aiju=20
    daij=$sigma(dparn([ni]))
    pname=a_[i]_[j]; tmp1=$format([ni],i-4); tmp2=$format([aij],e13.5); tmp3=$format([daij],e13.5); tmp4=tmp2=$format([aijl],e13.5); tmp5=tmp2=$format([aiju],e13.5);
    tmp=$unquote([tmp1])$quote([pname])   $unquote([tmp2])$unquote([tmp3]) 0 30
    fmess [tmp] [file]
  enddo
enddo
ls=[s1]-1
rs=[s1]+1
ni=[ni]+1; pname=s1; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5); tmp4=$format([ls],e13.5); tmp5=$format([rs],e13.5)
tmp=$unquote([tmp1])$quote([pname])      $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5]) 
fmess [tmp] [file]
ls=[ss]-1
rs=[ss]+1
ni=[ni]+1; pname=ss; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5); tmp4=$format([ls],e13.5); tmp5=$format([rs],e13.5)
tmp=$unquote([tmp1])$quote([pname])      $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
fmess [tmp] [file]
ls=[sm]-1
rs=[sm]+1
ni=[ni]+1; pname=sm; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5); tmp4=$format([ls],e13.5); tmp5=$format([rs],e13.5)
tmp=$unquote([tmp1])$quote([pname])      $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5]) 
fmess [tmp] [file]
ls=[s2]-1
rs=[s2]+1
ni=[ni]+1; pname=s2; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5); tmp4=$format([ls],e13.5); tmp5=$format([rs],e13.5)
tmp=$unquote([tmp1])$quote([pname])      $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
fmess [tmp] [file]
ls=[z1]-1
rs=[z1]+1
ni=[ni]+1; pname=z1; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5); tmp4=$format([ls],e13.5); tmp5=$format([rs],e13.5)
tmp=$unquote([tmp1])$quote([pname])      $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
fmess [tmp] [file]
ls=[z2]-1
rs=[z2]+1
ni=[ni]+1; pname=z2; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5); tmp4=$format([ls],e13.5); tmp5=$format([rs],e13.5)
tmp=$unquote([tmp1])$quote([pname])      $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
fmess [tmp] [file]
ni=[ni]+1; pname=sig1; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5)
tmp=$unquote([tmp1])$quote([pname])    $unquote([tmp2])$unquote([tmp3]) 0 1
fmess [tmp] [file]
ni=[ni]+1; pname=sig2; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5)
tmp=$unquote([tmp1])$quote([pname])    $unquote([tmp2])$unquote([tmp3]) 0 1
fmess [tmp] [file]
ni=[ni]+1; pname=pds; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5)
tmp=$unquote([tmp1])$quote([pname])     $unquote([tmp2])$unquote([tmp3]) -0.01 0.01
fmess [tmp] [file]
ni=[ni]+1; pname=z0; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5)
tmp=$unquote([tmp1])$quote([pname])      $unquote([tmp2])$unquote([tmp3]) -1 1
fmess [tmp] [file]
ni=[ni]+1; pname=sigz0; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5)
tmp=$unquote([tmp1])$quote([pname])   $unquote([tmp2])$unquote([tmp3]) 0 1
fmess [tmp] [file]
do i=1,[nbp]
  asi=$sigma(asi([i]))
  asil=0.8*[asi]
  asiu=1.2*[asi]
  dasi=$sigma(dasi([i]))
  ni=[ni]+1; pname=as_[i]; tmp1=$format([ni],i-4); tmp2=$format([asi],e13.5); tmp3=$format([dasi],e13.5); tmp4=tmp2=$format([asil],e13.5); tmp5=tmp2=$format([asiu],e13.5);
  tmp=$unquote([tmp1])$quote([pname])   $unquote([tmp2])$unquote([tmp3]) 0 30
  fmess [tmp] [file]
enddo
ls=[z1]-2
rs=[z1]+1
ni=[ni]+1; pname=zs1; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5); tmp4=$format([ls],e13.5); tmp5=$format([rs],e13.5)
tmp=$unquote([tmp1])$quote([pname])     $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
fmess [tmp] [file]
ls=[z2]-1
rs=[z2]+2
ni=[ni]+1; pname=zs2; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5); tmp4=$format([ls],e13.5); tmp5=$format([rs],e13.5)
tmp=$unquote([tmp1])$quote([pname])     $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
fmess [tmp] [file]
ls=[phiPMT]-1
rs=[phiPMT]+1
ni=[ni]+1; pname=phiPMT; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5); tmp4=$format([ls],e13.5); tmp5=$format([rs],e13.5)
tmp=$unquote([tmp1])$quote([pname])  $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
fmess [tmp] [file]
ni=[ni]+1; pname=rPMT; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5)
tmp=$unquote([tmp1])$quote([pname])    $unquote([tmp2])$unquote([tmp3]) 0 1
fmess [tmp] [file]
ni=[ni]+1; pname=aPMT; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5)
tmp=$unquote([tmp1])$quote([pname])    $unquote([tmp2])$unquote([tmp3]) 0 30
fmess [tmp] [file]
ni=[ni]+1; pname=a2PMT; tmp1=$format([ni],i-4); tmp2=$format([[pname]],e13.5); tmp3=$format([d[pname]],e13.5)
tmp=$unquote([tmp1])$quote([pname])   $unquote([tmp2])$unquote([tmp3]) -30 30
fmess [tmp] [file]
fmess '' [file]
fmess '' [file]
*fmess 'set err 0.5' [file]
*fmess 'set eps 1e-9' [file]
*fmess 'set strategy 2' [file]
*fmess '' [file]

do i=1,154
  txt=fix [i]
  fmess [txt] [file]
enddo
fmess 'mini' [file]
fmess 'mini' [file]
fmess 'mini' [file]
do i=1,154
  txt=release [i]
  fmess [txt] [file]
enddo
do i=155,186
  txt=fix [i]
  fmess [txt] [file]
enddo
fmess 'mini' [file]
fmess 'mini' [file]
fmess 'mini' [file]
do i=1,186
  txt=release [i]
  fmess [txt] [file]
enddo
fmess 'mini' [file]
fmess 'mini' [file]
fmess 'mini' [file]
fmess '' [file]
fmess 'ret' [file]

*---------------------------------
if ([p].ne.'parpl') then
*  shell ./mknew.sh [dir]
endif

return



macro accmapb xm=0 ym=17.5
np=200
ve/cre yr(1) r
ve/cre xr(1) r [xm]
ve/cre yci([np]) r
ve/cre zci([np]) r
pmin=0
pmax=50
dp=([pmax]-[pmin])/([np]-1)
do i=1,[np]
  ve/inp yr(1) $sigma([pmin]+[dp]*([i]-1))
  ve/inp yci([i]) yr(1)
  ve/inp zci([i]) $call('accmappl_blocks.sl(xr,yr,0)')
enddo
ve/cre dzc([np]) r
* exec $PER/s#vpl zci dzc yci dzc sz=0.01 ll=-1
read x
ve/inp yr(1) [ym]
pmin=-15
pmax=15
dp=([pmax]-[pmin])/([np]-1)
do i=1,[np]
  ve/inp xr(1) $sigma([pmin]+[dp]*([i]-1))
  ve/inp yci([i]) xr(1)
  ve/inp zci([i]) $call('accmappl_blocks.sl(xr,yr,0)')
enddo
* exec $PER/s#vpl zci dzc yci dzc sz=0.1 ll=-1
return

macro mapplxnew nc=1 y=0
hi/file 20 mhad2011_profile_sin.his
zmin=-15
zmix=15
nz=300
pmin=[nc]*40-60
pmax=[nc]*40+20
nf=320
iy=$sigma(int([nf]*([y]-([pmin]))/([pmax]-([pmin])))+1)
idh=[iy]*100+[nc]
mess [idh]
hi/pl //lun20/[idh]
close 20
return


macro fitx nsl=160 nsr=[nsl] idh=0 ncnt=1

cy=1
ve/cre pix(3) r 0.001 -0.1 0
ve/cre pix(3) r 0.001 1 0
x=2
if ([x].eq.1) then

npar=25
fname=mhad2011_amplitude_vs_phi_counter1__p1_p41_v2
ve/read pars,dpars [fname]_fit.par '2f15.10'
ve/del p16
ve/copy pars(1:25) p16
  s1=p16(12)
  ss=p16(13)
  sm=p16(14)
  s2=p16(15)
  ve/cre xi(10) r [s1] $sigma([s1]+([ss]-[s1])/3) $sigma([s1]+([ss]-[s1])*2/3) _
                  [ss] $sigma([ss]+([sm]-[ss])/3) $sigma([ss]+([sm]-[ss])*2/3) _
                  [sm] $sigma([sm]+([s2]-[sm])/3) $sigma([sm]+([s2]-[sm])*2/3) [s2]

i1=$sigma(int([s1]/15))
i2=$sigma(int([s2]/15))
ni=[i2]-[i1]
i1=[i1]+1
f1=[i1]*15
f2=[i2]*15
ve/cre si0([ni]) r
sigma si0 = array([ni],[f1]#[f2])
ve/cre si(3) r 3*-1000
ve/copy si0(1:[ni]) si(1:[ni])
dni=3-[ni]

*$PER/s#vpl ye dye yvi dyvi
np=$vlen(yvi)
lv=$sigma(vmin(yvi))
rv=$sigma(vmax(yvi))
l=$sigma([lv]-([rv]-[lv])/([np]-1)/2)
r=$sigma([rv]+([rv]-[lv])/([np]-1)/2)
if ([idh].eq.0) then
  1d 1000 ! [np] [l] [r]
  mess [l] [r]
  hi/put/cont 1000 ye
  hi/put/erro 1000 dye
else
  hi/del 1000
  hi/copy [idh] 1000
endif
ve/cre mod(4) r 1 1 1 1
hi/pl 1000
set hcol 2 
fun/pl scphip.f [l] [r] s
read x

ve/copy pars(1:25) s16
do i=12,25
ve/inp s16([i]) 0
enddo
set fcol 4
do i=1,2
  hi/fit 1000 scphi.f sb [npar] p16 s16
enddo
ve/copy pars(1:25) s16
do i=1,11
ve/inp s16([i]) 0
enddo
hi/fit 1000 scphi.f sb [npar] p16 s16
endif





if ([x].eq.2) then
npar=29
ve/cre chi2(2) r
ve/cre paru([npar]) r
ve/cre dparu([npar]) r
ve/cre covu([npar],[npar]) r
*
*
*fname=mhad2011_amplitude_vs_phi_counter[ncnt]__p1_p41_v2
*ve/read pars,dpars [fname]_fit.par '2f15.10'
*fname=mhad2011_amplitude_vs_phi_counter[ncnt]_slx145_fit_95.par
gl/imp fname
if ([fname].eq.'') then
  gl/imp expsim
  if ([expsim].eq.'exp') then
*  fname=mhad2011_amplitude_vs_phi_counter[ncnt]_slx145_fit_95.par
    fname=mhad2011-4_amplitude_vs_phi_counter[ncnt]_slx145_fit_95_et0_chit.par
  else
    fname=mhad2011-4_sim_v1_amplitude_vs_phi_counter[ncnt]_slx145_fit_95_et0_chit.par
  endif
endif
ve/del pars,dpars
ve/read pars,dpars [fname] '2f15.10'
ve/cre p16([npar]) r
ve/cre dp16([npar]) r
ve/copy pars(1:[npar]) p16
*ve/copy pars(8) p16(29)
*ve/inp p16(19) 0
*ve/inp p16(13) $sigma(17.5+40*([ncnt]-1))
*ve/inp p16(14) $sigma(31+40*([ncnt]-1))
ve/cre p16y(11) r
ve/copy p16(1:11) p16y(1:11)
sigma p16y = p16y*[cy]
ve/copy p16y(1:11) p16(1:11)
s1=p16(12)
ss=p16(13)
sm=p16(14)
s2=p16(15)
if ([s2].gt.360) then
  s2=362
endif
gl/imp ssm
if ([ssm].ne.0) then
  sm=([s2]+[ss])/2
  ve/inp p16(14) [sm]
endif
*
gl/cre s1g [s1]
gl/cre s2g [s2]
gl/cre ssg [ss]
gl/cre sigs1g 0.2
gl/cre sigs2g 0.2
gl/cre sigssg 0.2
*exec mapcal#s1s2 [s1] [s2] [ss] i1=[nsl] i2=[nsr] idh=[idh]
gl/imp s1g
gl/imp s2g
gl/imp ssg
gl/imp sigs1g
gl/imp sigs2g
gl/imp sigssg
s1=[s1g]
if ([ssf].ne.0) then
  ss=[ssg]
endif
s2=[s2g]
gl/imp expsim
if (([expsim].eq.'sim').or.([expsim].eq.'sim0')) then
*  s1=[s1g]
*  ss=[ssg]
*  s2=[s2g]
  gl/imp dphisi
*  ss=[dphisi]
*  ve/inp p16(13) [ss]
  gl/imp dphisf
  ssf=[dphisf]
  sig=[dphisf]
  sigs=[dphisf]
else
  ssf=0.01
  sig=0.02
  gl/imp dphisf
  ssf=[dphisf]
  sigs=[dphisf]
endif
ss1=0.02
ssf=0.02
ss2=0.02
if ($vexist(vss1)) then
  ss1=vss1(1)
endif
if ($vexist(vsss)) then
  ssf=vsss(1)
endif
if ($vexist(vss2)) then
  ss2=vss2(1)
endif
*
*
np=$vlen(yvi)
lv=$sigma(vmin(yvi))
rv=$sigma(vmax(yvi))
l=$sigma([lv]-([rv]-[lv])/([np]-1)/2)
r=$sigma([rv]+([rv]-[lv])/([np]-1)/2)
if ([idh].eq.0) then
  1d 1000 ! [np] [l] [r]
  mess [l] [r]
  hi/put/cont 1000 ye
  hi/put/erro 1000 dye
else
  hi/del 1000
  hi/copy [idh] 1000
endif

*gl/imp nbrn
*if ([nbrn].ne.0) then
*  hi/copy 1000 2000
*  exec ../SepPar/sp#brn 2000 [nbrn] 100 ! 1
*  hi/copy 2100 1000
*endif


ve/cre mod(4) r 1 1 1 1
*hi/pl 1000
*
ve/cre xi(10) r [s1] $sigma([s1]+([ss]-[s1])/3) $sigma([s1]+([ss]-[s1])*2/3) _
                [ss] $sigma([ss]+([sm]-[ss])/3) $sigma([ss]+([sm]-[ss])*2/3) _
                [sm] $sigma([sm]+([s2]-[sm])/3) $sigma([sm]+([s2]-[sm])*2/3) [s2]

i1=$sigma(int([s1]/15))
i2=$sigma(int([s2]/15))
ni=$sigma(min(3,[i2]-[i1]))
mess [s1] [s2] [i1] [i2]
i1=[i1]+1
f1=[i1]*15
f2=[i2]*15
ve/cre si0([ni]) r
sigma si0 = array([ni],[f1]#[f2])
ve/cre si(3) r 3*-1000
ve/copy si0(1:[ni]) si(1:[ni])
dni=$sigma(max(0,3-[ni]))

nfit=3

if ([dni].eq.0) then
  ve/cre s16([npar]) r 11*$sigma(0.5*[cy]) 0 0 0 0 1*1 0 0  0.0 [ni]*0.00 [ni]*0
else
  ve/cre s16([npar]) r 11*$sigma(0.5*[cy]) 0 0 0 0 1*1 0 0  0.0 [ni]*0.00 [dni]*0 [ni]*0 [dni]*0
endif
ve/inp s16(29) 0.5
ve/cre pmin([npar]) r 11*0 $sigma([s1]-1) $sigma([ss]-2) $sigma([sm]-1) $sigma([s2]-2) 0 3*0.1 3*0 3*-2 3*0 0
ve/cre pmax([npar]) r 11*50 $sigma([s1]+2) $sigma([ss]+2) $sigma([sm]+2) $sigma([s2]+1) 100 3*1.0 3*1 3*2 3*2 50
lf=$sigma([s1]-1.5*[sigs1g])
rf=$sigma([s2]+1.5*[sigs2g])
set ksiz 0.05
set fcol 2
hi/pl 1000 e
ve/inp p16(12) [s1]
ve/inp p16(13) [ss]
ve/inp p16(15) [s2]
ve/inp p16(20) 0
ve/inp p16(21) 0
ve/inp p16(22) 0
ve/inp p16(23) 0
ve/inp p16(24) 0
ve/inp p16(25) 0
ve/inp p16(26) 0.2
ve/inp p16(27) 0.2
ve/inp p16(28) 0.2
do i=1,[nfit]
  hi/fit 1000([lf]:[rf]) scphi.f qsb [npar] p16 s16 pmin pmax dp16
enddo
ve/inp s16(1) 0.0
ve/inp s16(2) 0.0
ve/inp s16(3) 0.0
ve/inp s16(4) $sigma(0.1*[cy])
ve/inp s16(5) $sigma(0.1*[cy])
ve/inp s16(6) 0.0
ve/inp s16(7) 0.0
ve/inp s16(8) $sigma(0.1*[cy])
ve/inp s16(29) $sigma(0.1*[cy])
ve/inp s16(9) 0.0
ve/inp s16(10) 0.0
ve/inp s16(11) 0.0
ve/inp s16(26) $sigma(0.01*[cy])
ve/inp s16(27) $sigma(0.01*[cy])
ve/inp s16(28) $sigma(0.01*[cy])
do i=1,[nfit]
  hi/fit 1000([lf]:[rf]) scphi.f qsb [npar] p16 s16 pmin pmax dp16
enddo
ve/inp s16(1) $sigma(0.1*[cy])
ve/inp s16(2) $sigma(0.1*[cy])
ve/inp s16(3) $sigma(0.1*[cy])
ve/inp s16(4) $sigma(0.1*[cy])
ve/inp s16(5) $sigma(0.1*[cy])
ve/inp s16(6) $sigma(0.1*[cy])
ve/inp s16(7) $sigma(0.1*[cy])
ve/inp s16(8) $sigma(0.1*[cy])
ve/inp s16(9) $sigma(0.1*[cy])
ve/inp s16(10) $sigma(0.1*[cy])
ve/inp s16(11) $sigma(0.1*[cy])
g=0.01
g=0.00
ve/inp s16(23) $sigma([g]*[cy])
ve/inp s16(24) $sigma([g]*[cy])
ve/inp s16(25) $sigma([g]*[cy])
do i=1,[nfit]
  hi/fit 1000([lf]:[rf]) scphi.f qsb [npar] p16 s16 pmin pmax dp16
enddo
ve/inp s16(19) 0.0
ve/inp s16(20) $sigma([g]*[cy])
ve/inp s16(21) $sigma([g]*[cy])
ve/inp s16(22) $sigma([g]*[cy])
do i=1,[nfit]
  hi/fit 1000([lf]:[rf]) scphi.f qsb [npar] p16 s16 pmin pmax dp16
enddo
ve/inp s16(20) 0.0
ve/inp s16(21) 0.0
ve/inp s16(22) 0.0
ve/inp s16(23) 0.0
ve/inp s16(24) 0.0
ve/inp s16(25) 0.0
*
ve/inp s16(12) $sigma([ss1]*[cy])
ve/inp s16(13) $sigma([ssf])
ve/inp s16(14) $sigma([ssm]*[cy])
ve/inp s16(15) $sigma([ss2]*[cy])
do i=1,[nfit]
  hi/fit 1000([lf]:[rf]) scphi.f qsb [npar] p16 s16 pmin pmax dp16
enddo
*hi/pl 1000
i=0
if ([i].gt.1) then
  sg=$sigma(abs(p16(19)))
  dg=$sigma(abs(p16(20)))
  ag=$sigma([dg]/sqrt(2*3.1415927)/[sg]*exp(-([dg]/[sg])**2/8))
  if (([ag].gt.0.1).and.([dg].gt.0.05)) then
    ve/inp s16(23) $sigma([g]*[cy])
  endif
  dg=$sigma(abs(p16(21)))
  ag=$sigma([dg]/sqrt(2*3.1415927)/[sg]*exp(-([dg]/[sg])**2/8))
  if (([ag].gt.0.1).and.([dg].gt.0.05)) then
    ve/inp s16(24) $sigma([g]*[cy])
  endif
  dg=$sigma(abs(p16(22)))
  ag=$sigma([dg]/sqrt(2*3.1415927)/[sg]*exp(-([dg]/[sg])**2/8))
  if (([ag].gt.0.1).and.([dg].gt.0.05)) then
    ve/inp s16(25) $sigma([g]*[cy])
  endif
  ve/inp s16(17) $sigma([sig]*[cy])
  ve/inp s16(18) [sigs]
  ve/inp s16(19) 0.0
  do i=1,[nfit]
    hi/fit 1000([lf]:[rf]) scphi.f qsb [npar] p16 s16 pmin pmax dp16
  enddo
endif
if ([cy].eq.0) then
  goto 10
endif
sigr=$sigma(abs(p16(17)))
do i=1,[nfit]
  sigr=$sigma(abs(p16(17)))
  l=$sigma([s1]-0*[sigr]+0.1)
  r=$sigma(min([s2]+0*[sigr]+0.1,363.1))
  s1=p16(12)
  ss=p16(13)
  sm=p16(14)
  s2=p16(15)
  ve/cre xi(10) r [s1] $sigma([s1]+([ss]-[s1])/3) $sigma([s1]+([ss]-[s1])*2/3) _
                  [ss] $sigma([ss]+([sm]-[ss])/3) $sigma([ss]+([sm]-[ss])*2/3) _
                  [sm] $sigma([sm]+([s2]-[sm])/3) $sigma([sm]+([s2]-[sm])*2/3) [s2]
  ve/inp s16(17) 0.0
  ve/inp s16(18) 0.0
  ve/inp s16(19) 0.0
  ve/inp s16(23) 0.0
  ve/inp s16(24) 0.0
  ve/inp s16(25) 0.0
  ve/inp s16(19) 0.0
  ve/inp s16(20) 0.0
  ve/inp s16(21) 0.0
  ve/inp s16(22) 0.0
  hi/fit 1000([l]:[r]) scphi.f qsb [npar] p16 s16 pmin pmax dp16
  ve/inp s16(17) $sigma([sig]*[cy])
  ve/inp s16(18) [sigs]
  ve/inp s16(19) 0.0
  ve/inp s16(23) 0.0
  ve/inp s16(24) 0.0
  ve/inp s16(25) 0.0
  ve/inp s16(19) 0.0
  ve/inp s16(20) 0.0
  ve/inp s16(21) 0.0
  ve/inp s16(22) 0.0
  hi/fit 1000([l]:[r]) scphi.f qsb [npar] p16 s16 pmin pmax dp16
enddo
do i=1,[nfit]
  hi/fit 1000([l]:[r]) scphi.f qsb [npar] p16 s16 pmin pmax dp16
enddo
nf=0
while (($sigma(dp16(16)).gt.3).and.([nf].lt.5)) do
  nf=[nf]+1
  hi/fit 1000([l]:[r]) scphi.f qsb [npar] p16 s16 pmin pmax dp16
endwhile
hi/fit 1000([l]:[r]) scphi.f b [npar] p16 s16 pmin pmax dp16
10:
call covm.f(1)
call covmpen(chi2,[npar],paru,dparu)
call covmcov([npar],covu)
*ve/cre ndf(1) r $sigma(chi2(2)-vsum(dparu/dparu))
ve/cre ndf(1) r $sigma(chi2(2))
endif



if ([x].eq.3) then
npar=29
ve/cre chi2(2) r
ve/cre paru([npar]) r
ve/cre dparu([npar]) r
ve/cre covu([npar],[npar]) r
*
gl/imp fname
ve/read pars,dpars [fname] '2f15.10'
ve/cre s16([npar]) r
ve/cre p16([npar]) r
ve/cre dp16([npar]) r
ve/copy pars(1:[npar]) p16
ve/copy dpars(1:[npar]) s16
sigma s16=abs(s16)
ve/cre xi(10) r
ve/copy pars(30:39) xi(1:10)
s1=p16(12)
ss=p16(13)
sm=p16(14)
s2=p16(15)
*
exec mapcal#s1s2 [s1] [s2] [ss] i1=[nsl] i2=[nsr] idh=[idh]
*
gl/imp expsim
if ([expsim].eq.'sim') then
  gl/imp dphisi
  ss=[dphisi]
  ve/inp p16(13) [ss]
  gl/imp dphisf
  ve/inp s16(13) [dphisf]
endif
*
np=$vlen(yvi)
lv=$sigma(vmin(yvi))
rv=$sigma(vmax(yvi))
l=$sigma([lv]-([rv]-[lv])/([np]-1)/2)
r=$sigma([rv]+([rv]-[lv])/([np]-1)/2)
if ([idh].eq.0) then
  1d 1000 ! [np] [l] [r]
  mess [l] [r]
  hi/put/cont 1000 ye
  hi/put/erro 1000 dye
else
  hi/del 1000
  hi/copy [idh] 1000
endif
ve/cre mod(4) r 1 1 1 1
hi/pl 1000
*
i1=$sigma(int([s1]/15))
i2=$sigma(int([s2]/15))
ni=[i2]-[i1]
i1=[i1]+1
f1=[i1]*15
f2=[i2]*15
ve/cre si0([ni]) r
sigma si0 = array([ni],[f1]#[f2])
ve/cre si(3) r 3*-1000
ve/copy si0(1:[ni]) si(1:[ni])
dni=3-[ni]

ve/inp s16(29) 0.5
ds=0.5
ve/cre pmin([npar]) r 11*0 $sigma([s1]-[ds]) $sigma([ss]-[ds]) $sigma([sm]-[ds]) $sigma([s2]-[ds]) 0 3*0.2 3*0 3*-2 3*0 0
ve/cre pmax([npar]) r 11*50 $sigma([s1]+[ds]) $sigma([ss]+[ds]) $sigma([sm]+[ds]) $sigma([s2]+[ds]) 100 3*1.0 3*1 3*2 3*2 50

set ksiz 0.05
set fcol 2
hi/pl 1000 e

sig=$sigma(abs(p16(17)))
do i=1,5
  l=$sigma([s1]-0*[sig]-0.01)
  r=$sigma([s2]+0*[sig]+0.01)
  s1=p16(12)
  ss=p16(13)
  sm=p16(14)
  s2=p16(15)
  ve/cre xi(10) r [s1] $sigma([s1]+([ss]-[s1])/3) $sigma([s1]+([ss]-[s1])*2/3) _
                  [ss] $sigma([ss]+([sm]-[ss])/3) $sigma([ss]+([sm]-[ss])*2/3) _
                  [sm] $sigma([sm]+([s2]-[sm])/3) $sigma([sm]+([s2]-[sm])*2/3) [s2]
  ve/inp s16(17) 0.0
  ve/inp s16(18) 0.0
  ve/inp s16(19) 0.0
  ve/inp s16(23) 0.0
  ve/inp s16(24) 0.0
  ve/inp s16(25) 0.0
  ve/inp s16(19) 0.0
  ve/inp s16(20) 0.0
  ve/inp s16(21) 0.0
  ve/inp s16(22) 0.0
  hi/fit 1000([l]:[r]) scphi.f sb [npar] p16 s16 pmin pmax dp16
  ve/inp s16(17) 0.01
  ve/inp s16(18) 0.01
  ve/inp s16(19) 0.0
  ve/inp s16(23) 0.0
  ve/inp s16(24) 0.0
  ve/inp s16(25) 0.0
  ve/inp s16(19) 0.0
  ve/inp s16(20) 0.0
  ve/inp s16(21) 0.0
  ve/inp s16(22) 0.0
  hi/fit 1000([l]:[r]) scphi.f sb [npar] p16 s16 pmin pmax dp16
enddo
do i=1,5
  hi/fit 1000([l]:[r]) scphi.f sb [npar] p16 s16 pmin pmax dp16
enddo
nf=0
while (($sigma(dp16(13)).gt.0.1).and.([nf].lt.10)) do
  nf=[nf]+1
  hi/fit 1000([l]:[r]) scphi.f sb [npar] p16 s16 pmin pmax dp16
endwhile
call covm.f(1)
call covmpen(chi2,[npar],paru,dparu)
call covmcov([npar],covu)
*ve/cre ndf(1) r $sigma(chi2(2)-vsum(dparu/dparu))
ve/cre ndf(1) r $sigma(chi2(2))
endif


if ([x].eq.4) then
npar=4
ve/cre chi2(2) r
ve/cre paru([npar]) r
ve/cre dparu([npar]) r
ve/cre covu([npar],[npar]) r
*
gl/imp fname
if ([fname].eq.'') then
  gl/imp expsim
  if ([expsim].eq.'exp') then
    fname=mhad2011-4_amplitude_vs_phi_counter[ncnt]_slx145_fit_95_et0_chit.par
  else
    fname=mhad2011-4_sim_v1_amplitude_vs_phi_counter[ncnt]_slx145_fit_95_et0_chit.par
  endif
endif
ve/read pars,dpars [fname] '2f15.10'
ve/cre p4(4) r $sigma(pars(16)) $sigma(pars(13)) $sigma(pars(18)) 0
ve/cre dp4(4) r
xmin=$sigma(pars(13)-10.)
xmax=$sigma(pars(13)+10.)
hi/copy [idh] 1000
do i=1,5
  hi/fit 1000([xmin]:[xmax]) scphiw.f s 4 p4 ! ! ! dp4
enddo
ve/cre p16(29) r
ve/cre dp16(29) r
ve/copy pars(12:28) p16(12:28)
ve/copy dpars(12:28) dp16(12:28)
ve/inp p16(16) $sigma(p4(1)); ve/inp dp16(16) $sigma(dp4(1))
ve/inp p16(13) $sigma(p4(2)); ve/inp dp16(13) $sigma(dp4(2))
ve/inp p16(18) $sigma(p4(3)); ve/inp dp16(18) $sigma(dp4(3))
  s1=p16(12)
  ss=p16(13)
  sm=p16(14)
  s2=p16(15)
  ve/cre xi(10) r [s1] $sigma([s1]+([ss]-[s1])/3) $sigma([s1]+([ss]-[s1])*2/3) _
                  [ss] $sigma([ss]+([sm]-[ss])/3) $sigma([ss]+([sm]-[ss])*2/3) _
                  [sm] $sigma([sm]+([s2]-[sm])/3) $sigma([sm]+([s2]-[sm])*2/3) [s2]
ve/cre dxi(10) r
call covm.f(1)
call covmpen(chi2,[npar],paru,dparu)
call covmcov([npar],covu)
ve/cre ndf(1) r $sigma(chi2(2))

endif

do i=1,11
  ve/inp p16([i]) $sigma(abs(p16([i])))
enddo
i=29;  ve/inp p16([i]) $sigma(abs(p16([i])))

return


macro fitxsv
ve/read nfn nfn.txt
do i=1,$vlen(nfn)
  if ($sigma(nfn([i])).eq.1) then
    exec mapcal#fitxs [i] [i] 1 1 f 2
  endif
enddo
return


macro fitxsz0old ncnt=1 z1=-10 z2=9.
*ncnt=1
*exec mapcal#fitpl 110[ncnt] 0
*n=$vlen(zic)
mnz=300
mzmin=-15
mzmax=15
*do i=1,[n]-1
*  z1=zic([i])
  i1=$sigma([mnz]*([z1]-[mzmin])/([mzmax]-[mzmin]))
  j=[i]+1
*  z2=zic([j])
  i2=$sigma([mnz]*([z2]-[mzmin])/([mzmax]-[mzmin]))
  ic=$sigma(int(([i1]+[i2])/2+0.5))
  di=$sigma(int(([i2]-[i1])/2+0.5))

  gl/cre expsim sim
  if ([expsim].eq.'exp') then
    hi/file 30 mhad2011_profile_sin_new.his
  else
    hi/file 30 profilex_ee_sim_corr_sim.his
  endif

  idh=200
  if ($hexist([idh]).eq.1) then
    hi/del [idh]
  endif
  i1=$sigma(int([i1]+0.5))
  i2=$sigma(int([i2]+0.5))
  do j=[i1],[i2]

    idhj=[j]*100+10*3+[ncnt]
    
    if ($hexist([idh]).eq.0) then
      hi/copy [idhj] [idh]
    else
      hi/op/add [idhj] [idh] [idh]
    endif
  enddo
  
  ve/cre z1b(1) r [z1]
  ve/cre z2b(1) r [z2]
  exec mapcal#fitxs [ic] [ic] 1 1 f [di] idh=[idh] ncnt=[ncnt]
  close 30
*enddo
return

macro fitxszpl map=1 i1=1 i2=9 
do i=[i1],[i2]
  exec mapcal#fitxsz [map] [i] 1 14 p
enddo
return

macro fitxszbx map=1 i1=1 i2=9 opt='f'
do i=[i1],[i2]
  exec mapcal#fitxszb [map] [i] [opt]
enddo
return

macro fitxszb map=1 nc=1 opt='f' i1=1 i2=14
do i=[i1],[i2]
  if ([opt].eq.'f') then
    fname=fitxs[map]_counter[nc]_slix[i].kumac
    if ($fexist([fname]).eq.1) then
      shell rm [fname]
    endif
    for/file  20 [fname] new
    close 20
    txt=exec mapcal#fitxsz [map] [nc] [i] [i] [opt]
    fmess [txt] [fname]
    shell .testrelease/.mainrelease/Offline/submit.sh -q clusters,180 pawbigX11 -b [fname]
  endif
  if ([opt].eq.'fl') then
    fname=fitxs[map]_counter[nc]_slix[i].kumac
    if ($fexist([fname]).eq.1) then
      shell rm [fname]
    endif
    for/file  20 [fname] new
    close 20
    txt=exec mapcal#fitxsz [map] [nc] [i] [i] f
    fmess [txt] [fname]
    shell pawbigX11 -b [fname]
*    exec mapcal#fitxsz [map] [nc] [i] [i] f
  endif
  if ([opt].eq.'p') then
    exec mapcal#fitxsz [map] [nc] [i] [i] [opt]
  endif
enddo
return


macro fitxsz map=1 ncnt=1 i1f=1 i2f=[n]-1 opt='f' ism=0
gl/cre ssm 0
gl/cre expsim exp
gl/cre ver 0
*
if ([map].eq.1) then
  gl/cre mname mapcal1
  gl/cre nh1 19
  gl/cre nh2 19
  gl/cre nb 14
  npt=0
  dnb=3
  nbr=[nb]-3
  ve/cre vss1i([nb]) r [nb]*0.02
  ve/cre vsssi([nb]) r [nbr]*0.02 [dnb]*0 
  ve/cre vss2i([nb]) r [nb]*0.02
endif
if ([map].eq.2) then
  gl/cre mname mapcal2
  gl/cre nh1 20
  gl/cre nh2 20
  gl/cre nb 14
  npt=1
  dnb=3
  nbr=[nb]-3
  ve/cre vss1i([nb]) r [nb]*0.02
  ve/cre vsssi([nb]) r [nbr]*0.02 [dnb]*0
  ve/cre vss2i([nb]) r [nb]*0.02
endif
if ([map].eq.3) then
  gl/cre mname mapcal3
  gl/cre nh1 21
  gl/cre nh2 21
  gl/cre nb 10
  npt=2
  dnb=2
  nbr=[nb]-[dnb]
  ve/cre vss1i([nb]) r [nb]*0.0
  ve/cre vsssi([nb]) r [nbr]*0.02 [dnb]*0
  ve/cre vss2i([nb]) r [nb]*0.0
endif
if ([map].eq.4) then
  gl/cre mname mapcal4
  gl/cre nh1 22
  gl/cre nh2 22
  gl/cre nb 14
  npt=3
  dnb=3
  nbr=[nb]-3
  ve/cre vss1i([nb]) r [nb]*0
  ve/cre vsssi([nb]) r [nbr]*0.02 [dnb]*0
  ve/cre vss2i([nb]) r [nb]*0
  if ([ncnt].eq.9) then
    ve/cre vsssi([nb]) r [nb]*0
  endif
endif
*
do l=1,3
  ve/del pars0,dpars0
  flname=[mname]/[mname]_amplitude_vs_zr_counter[l]10[ncnt]_0.log.pars
  if ($fexist([flname]).eq.0) then
    flname=[mname]/[mname]_amplitude_vs_zr_counter[l]10[ncnt]_0.log
    shell /work/snd2000/users/konctbel/SNDSimACC/readlifstat [flname]
  endif
  ve/read pars0,dpars0 [mname]/[mname]_amplitude_vs_zr_counter[l]10[ncnt]_0.log.pars
  nr=[nb]+1
  ve/cre zi([nr]) r
  i1=2; i2=[nr]+1
  ve/copy pars0([i1]:[i2]) zi(1:[nr])
  exec mapcal#zicr zi zic[l]
enddo
sigma zic = (zic1+zic2+zic3)/3
*
ftname=[mname]/[mname]_phi_profiles_blocks_counter[ncnt].tex
if ($fexist([ftname]).eq.1) then
  shell rm [ftname]
endif
for/file  20 [ftname] new
close 20
*
ve/del dphis
*ve/read dphis phi_shift_sim.txt
ve/cre dphis(9) r
ve/cre wp(1) r 1
ve/cre pix(3) r 0.001 1 0
*
i1f0=$sigma(max(1,[i1f]))
i2f0=$sigma(min([nb],[i2f]))
if ([ism].eq.0) then
  i1fs=[i1f0]
  i2fs=[i2f0]
else
  i1fs=[i1f0]
  i2fs=[i1fs]
endif
*
do i=[i1fs],[i2fs]

  if ([ism].eq.0) then
    j=[i]+1
  else
    j=[i2f0]+1
  endif
  z1=zic([i])
  z2=zic([j])
  exec mapcal#histprep z [ncnt] [nh1] [nh2] [z1] [z2]
  gl/imp ic
  gl/imp di

  gl/cre dphisi $sigma(([ncnt]-0.5)*40-2.5+dphis([ncnt]))
  gl/cre dphisf $sigma(vsssi([i]))
  ve/cre vss1(1) r $sigma(vss1i([i]))
  ve/cre vsss(1) r $sigma(vsssi([i]))
  ve/cre vss2(1) r $sigma(vss2i([i]))
  ve/cre z1b(1) r [z1]
  ve/cre z2b(1) r [z2]
  ve/cre lthpl(1) r 1
  
  gl/cre fname [mname]/[mname]_exp_v[ver]_amplitude_vs_phi_counter[ncnt]_slx150_fit_100_et0_chit.par
  exec mapcal#fitxs [ic] [ic] 1 1 [opt] [di] idh=200 ncnt=[ncnt]
  
  gl/imp epsfile
  shell mv [epsfile] [mname]/
  if ([opt].eq.'f') then
    gl/imp fname
    if ($fexist([fname])) then
      shell mv [fname] [mname]/
    else
      mess File [fname] does not exists
    endif
  endif
  
  fmess '\begin{figure}[ht!b]' [ftname]
  fmess '  \begin{minipage}{\textwidth}' [ftname]
  fmess '     {\centering\resizebox*{\textwidth}{!}' [ftname]
  txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
  fmess [txt] [ftname]
  txt=$unquote('     ')\caption{[mname]: counter \No [i]}
  fmess [txt] [ftname]
  fmess '   \end{minipage}' [ftname]
  fmess '\end{figure}' [ftname]
*
enddo
*
  if (([opt].eq.'p').and.([ncnt].eq.9)) then
    fname=[mname]_latex.sh
    if ($fexist([fname]).eq.1) then
      shell rm [fname]
    endif
    for/file  20 [fname] new
    close 20
  
    txt=cd [mname]
    fmess [txt] [fname]
    txt=latex [mname]_phi_profiles_blocks
    fmess [txt] [fname]
    txt=latex [mname]_phi_profiles_blocks
    fmess [txt] [fname]
    txt=dvips [mname]_phi_profiles_blocks
    fmess [txt] [fname]
    
    shell chmod +x [fname]
    shell ./[fname]
  endif
return



macro fitxexp z1=-10 z2=9.
gl/cre ssm 0.01
gl/cre mname mhad2011-4
gl/cre expsim exp
ve/cre z1b(1) r [z1]
ve/cre z2b(1) r [z2]
ve/del dphis
ve/read dphis phi_shift_sim.txt
ve/del phisf
ve/cre phisf(9) r
do i=1,9
  gl/cre dphisi $sigma(([i]-0.5)*40-2.5+dphis([i]))
  gl/cre dphisf $sigma(phisf([i]))
  exec mapcal#expsimf [i] [z1] [z2]
  gl/imp ic
  gl/imp di
  exec mapcal#fitxs [ic] [ic] 1 1 f [di] idh=200 ncnt=[i]
enddo
return


macro fitxexpnb pt=1 z1=-10 z2=10. i1=1 i2=9 opt='f'
do i=[i1],[i2]
  if ([opt].eq.'f') then
    fname=setup[pt]_fitxexpn_counter[i].kumac
    if ($fexist([fname]).eq.1) then
      shell rm [fname]
    endif
    for/file  20 [fname] new
    close 20
    txt=exec mapcal#fitxexpn [pt] [z1] [z2] [i] [i] [opt]
    fmess [txt] [fname]
    shell .testrelease/.mainrelease/Offline/submit.sh -q clusters,180 pawbigX11 -b [fname]
  endif
  if ([opt].eq.'p') then
    exec mapcal#fitxexpn [pt] [z1] [z2] [i] [i] [opt]
  endif
enddo
return



macro fitxexpn pt=1 z1=-10 z2=10. i1=1 i2=9 opt='f'
gl/cre ssm 0.01
ks=0
npt=[pt]-1
if ([pt].eq.0) then
  gl/cre mname setup[pt]
  nh1=1
  nh2=1
  gl/cre expsim exp
endif
if ([pt].eq.1) then
  gl/cre mname setup[pt]
  nh1=2
  nh2=5
  gl/cre expsim exp
endif
if ([pt].eq.2) then
  gl/cre mname setup[pt]
  nh1=6
  nh2=7
  gl/cre expsim exp
endif
if ([pt].eq.3) then
  gl/cre mname setup[pt]
  nh1=8
  nh2=8
  gl/cre expsim exp
endif
if ([pt].eq.4) then
  gl/cre mname setup[pt]
  nh1=9
  nh2=9
  gl/cre expsim exp
endif
if ([pt].eq.15) then
  gl/cre mname mapcal1
  nh1=15
  nh2=15
  npt=0
  gl/cre expsim exp
endif
if ([pt].eq.16) then
  gl/cre mname mapcal2
  nh1=16
  nh2=16
  npt=1
  gl/cre expsim exp
endif
if ([pt].eq.17) then
  gl/cre mname mapcal3
  nh1=17
  nh2=17
  gl/cre expsim exp
  npt=2
endif
if ([pt].eq.18) then
  gl/cre mname mapcal4
  nh1=18
  nh2=18
  gl/cre expsim exp
  npt=3
endif
if ([pt].eq.115) then
  gl/cre mname mapcal1
  nh1=115
  nh2=115
  ks=1
  gl/cre expsim sim
  npt=0
endif
if ([pt].eq.116) then
  gl/cre mname mapcal2
  nh1=116
  nh2=116
  ks=1
  gl/cre expsim sim
  npt=1
endif
if ([pt].eq.117) then
  gl/cre mname mapcal3
  nh1=117
  nh2=117
  ks=1
  gl/cre expsim sim
  npt=2
endif
if ([pt].eq.118) then
  gl/cre mname mapcal4
  nh1=118
  nh2=118
  ks=1
  gl/cre expsim sim
  npt=3
endif
if ([pt].eq.'s') then
  gl/cre mname setup[pt]
  mess [mname]
  nh1=10
  nh2=10
  ks=1
  gl/cre expsim sim
  npt=0
endif

ftname=[mname]_[expsim]_phi_profiles.tex
if ($fexist([ftname]).eq.1) then
  shell rm [ftname]
endif
for/file  20 [ftname] new
close 20
*
gl/cre ver 0
ve/cre wp(1) r 0
ve/cre z1b(1) r [z1]
ve/cre z2b(1) r [z2]
ve/del dphis
ve/read dphis /work/users/konctbel/MinuitTest/phi_shift_sim.txt
ve/del phisf
ve/cre phisf(9) r 9*0.01
ve/del ptl
ve/read ptl /work/users/konctbel/MinuitTest/phicuts_test.txt
do i=[i1],[i2]
  gl/cre dphisi $sigma(([i]-0.5)*40-2.5+dphis([i]))
  gl/cre dphisf $sigma(phisf([i]))
*  exec mapcal#expsimf [i] [z1] [z2]
  exec mapcal#histprep z [i] [nh1] [nh2] [z1] [z2]
  gl/imp ic
  gl/imp di
  ve/cre lthpl(1) r 1
  gl/cre fname [mname]/[mname]_[expsim]_v[ver]_amplitude_vs_phi_counter[i]_slx150_fit_100_et0_chit.par
  if ($fexist([fname]).eq.0) then
    gl/cre fname [mname]/[mname]_exp_v[ver]_amplitude_vs_phi_counter[i]_slx150_fit_100_et0_chit.par
  endif
  exec mapcal#fitxs [ic] [ic] 1 1 [opt] [di] idh=200 ncnt=[i]
*  read x
  gl/imp epsfile
  ind=$sigma(128*[npt]+14*([i]-1)+14)
  dxl=$sigma(ptl([ind]))
  do j=1,5
    ind=$sigma(128*[npt]+14*([i]-1)+8+[j])
    xl=$sigma(ptl([ind])-[dxl]*[ks])
    line [xl] $GRAFINFO('WNYMIN') [xl] $GRAFINFO('WNYMAX')
  enddo
  exec save [epsfile] f
  fmess '\begin{figure}[ht!b]' [ftname]
  fmess '  \begin{minipage}{\textwidth}' [ftname]
  fmess '     {\centering\resizebox*{\textwidth}{!}' [ftname]
  txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
  fmess [txt] [ftname]
  txt=$unquote('     ')\caption{[mname]: counter \No [i]}
  fmess [txt] [ftname]
  fmess '   \end{minipage}' [ftname]
  fmess '\end{figure}' [ftname]
*
  if (([pt].ge.15).and.([pt].le.118)) then
    shell mv [epsfile] [mname]/
    gl/imp fname
    if ($fexist([fname])) then
      shell mv [fname] [mname]/
    else
      mess File [fname] does not exists
    endif
  endif
*  read x
enddo
if (([pt].ge.15).and.([pt].le.118)) then
  shell mv [ftname] [mname]/
  fname=[mname]_latex.sh
  if ([opt].eq.'p') then
    if ($fexist([fname]).eq.1) then
      shell rm [fname]
    endif
    for/file  20 [fname] new
    close 20
  
    txt=cd [mname]
    fmess [txt] [fname]
    txt=latex phi_profiles_cuts
    fmess [txt] [fname]
    txt=latex phi_profiles_cuts
    fmess [txt] [fname]
    txt=dvips phi_profiles_cuts
    fmess [txt] [fname]
    
    shell chmod +x [fname]
    shell ./[fname]
  endif
endif
return


macro histprep par=z ncnt=1 nh1=2 nh2=5 p1=-8.0 p2=8.0 ver=0
nh=0
nh=[nh]+1; nb[nh]=mhad2010_profile.hbook
nh=[nh]+1; nb[nh]=mhad2011-4_profile.hbook
nh=[nh]+1; nb[nh]=mhad2012-2_profile.hbook
nh=[nh]+1; nb[nh]=phi2013-0_profile.hbook
nh=[nh]+1; nb[nh]=omeg2012-0_n1.13_profile.hbook
nh=[nh]+1; nb[nh]=omeg2012-0_n1.05_profile.hbook
nh=[nh]+1; nb[nh]=rho_2012-0_n1.05-1_profile.hbook
nh=[nh]+1; nb[nh]=rho_2012-0_n1.13_profile.hbook
nh=[nh]+1; nb[nh]=rho_2012-0_n1.05-2_profile.hbook
nh=[nh]+1; nb[nh]=sim_profile.hbook
*
nh=[nh]+1; nb[nh]=mapcal1_profilex_exp.his
nh=[nh]+1; nb[nh]=mapcal2_profilex_exp.his
nh=[nh]+1; nb[nh]=mapcal3_profilex_exp.his
nh=[nh]+1; nb[nh]=mapcal4_profilex_exp.his
*
nh=[nh]+1; nb[nh]=mapcal1_profilex_exp_pds.his
nh=[nh]+1; nb[nh]=mapcal2_profilex_exp_pds.his
nh=[nh]+1; nb[nh]=mapcal3_profilex_exp_pds.his
nh=[nh]+1; nb[nh]=mapcal4_profilex_exp_pds.his
*
nh=[nh]+1; nb[nh]=mapcal1_profilex_exp_ltr.his
nh=[nh]+1; nb[nh]=mapcal2_profilex_exp_ltr.his
nh=[nh]+1; nb[nh]=mapcal3_profilex_exp_ltr.his
nh=[nh]+1; nb[nh]=mapcal4_profilex_exp_ltr.his
*-- sim
nh=110
nh=[nh]+1; nb[nh]=mapcal1_profilex_sim_v[ver].his
nh=[nh]+1; nb[nh]=mapcal2_profilex_sim_v[ver].his
nh=[nh]+1; nb[nh]=mapcal3_profilex_sim_v[ver].his
nh=[nh]+1; nb[nh]=mapcal4_profilex_sim_v[ver].his
*
nh=[nh]+1; nb[nh]=mapcal1_profilex_sim_v[vet]_pds.his
nh=[nh]+1; nb[nh]=mapcal2_profilex_sim_v[ver]_pds.his
nh=[nh]+1; nb[nh]=mapcal3_profilex_sim_v[ver]_pds.his
nh=[nh]+1; nb[nh]=mapcal4_profilex_sim_v[ver]_pds.his

if ([par].eq.'z') then
  pzf=3
  mnz=300
  mzmin=-15
  mzmax=15
  i1=$sigma([mnz]*([p1]-[mzmin])/([mzmax]-[mzmin]))
  i2=$sigma([mnz]*([p2]-[mzmin])/([mzmax]-[mzmin]))
else
  pzf=0
  mnf=320
  mfmax=[ncnt]*40+20
  mfmin=[ncnt]*40-60
  i1=$sigma(int([mnf]*([p1]-[mfmin])/([mfmax]-[mfmin])+0.5))
  i2=$sigma(int([mnf]*([p2]-[mfmin])/([mfmax]-[mfmin])+0.5))
endif
ic=$sigma(int(([i1]+[i2])/2+0.5))
di=$sigma(int(([i2]-[i1])/2+0.5))
mess [ic] [di]
gl/cre ic [ic]
gl/cre di [di]


idh=200
if ($hexist([idh]).eq.1) then
  hi/del [idh]
endif

i1=$sigma(int([i1]+0.5))
i2=$sigma(int([i2]+0.5))
do i=[nh1],[nh2]

  hi/file 30 [nb[i]]
  mess [nb[i]] is opend

  do j=[i1],[i2]

    idhj=[j]*100+10*[pzf]+[ncnt]
    
    if ($hexist([idh]).eq.0) then
      hrin [idhj]
*      hi/pl //lun30/[idhj]
      nx=$hinfo([idhj],'xbins')
      xmin=$hinfo([idhj],'xmin')
      xmax=$hinfo([idhj],'xmax')
      if ([nx].eq.0) then
        mess [idhj]
        hi/pl [idhj]
        read x
      endif
      1d [idh] ! [nx] [xmin] [xmax]
      ve/del xis,xi2s,si2s,xi,si,xi2,si2,sis,ni
      ve/cre ni([nx]) r
      ve/cre xis([nx]) r
      ve/cre xi2s([nx]) r
      ve/cre si2s([nx]) r
      ve/cre xi([nx]) r
      ve/cre xi2([nx]) r
      ve/cre si([nx]) r
      ve/cre si2([nx]) r
      ve/cre sis([nx]) r
    endif
    hi/get/cont [idhj] xi
    hi/get/err  [idhj] si
    sigma ni = ni + xi/xi
    sigma xi2 = xi**2
    sigma si2 = si**2
    sigma xis = xis + xi
    sigma xi2s = xi2s + xi2
    sigma si2s = si2s + si2
    hi/del [idhj]
  enddo
  close 30
enddo
sigma xis = xis/ni
sigma sis = sqrt(((si2s+xi2s)/ni-xis**2)/ni)
hi/put/cont [idh] xis
hi/put/erro [idh] sis
return




macro histprepz ncnt=1 nh1=2 nh2=5 z1=-8.0 z2=8.0
nh=0
nh=[nh]+1; nb[nh]=mhad2010_profile.hbook
nh=[nh]+1; nb[nh]=mhad2011-4_profile.hbook
nh=[nh]+1; nb[nh]=mhad2012-2_profile.hbook
nh=[nh]+1; nb[nh]=phi2013-0_profile.hbook
nh=[nh]+1; nb[nh]=omeg2012-0_n1.13_profile.hbook
nh=[nh]+1; nb[nh]=omeg2012-0_n1.05_profile.hbook
nh=[nh]+1; nb[nh]=rho_2012-0_n1.05-1_profile.hbook
nh=[nh]+1; nb[nh]=rho_2012-0_n1.13_profile.hbook
nh=[nh]+1; nb[nh]=rho_2012-0_n1.05-2_profile.hbook
nh=[nh]+1; nb[nh]=sim_profile.hbook
*
nh=[nh]+1; nb[nh]=mapcal1_profilex_exp.his
nh=[nh]+1; nb[nh]=mapcal2_profilex_exp.his
nh=[nh]+1; nb[nh]=mapcal3_profilex_exp.his
nh=[nh]+1; nb[nh]=mapcal4_profilex_exp.his
*
nh=[nh]+1; nb[nh]=mapcal1_profilex_exp_pds.his
nh=[nh]+1; nb[nh]=mapcal2_profilex_exp_pds.his
nh=[nh]+1; nb[nh]=mapcal3_profilex_exp_pds.his
nh=[nh]+1; nb[nh]=mapcal4_profilex_exp_pds.his

mnz=300
mzmin=-15
mzmax=15
i1=$sigma([mnz]*([z1]-[mzmin])/([mzmax]-[mzmin]))
i2=$sigma([mnz]*([z2]-[mzmin])/([mzmax]-[mzmin]))
ic=$sigma(int(([i1]+[i2])/2+0.5))
di=$sigma(int(([i2]-[i1])/2+0.5))
mess [ic] [di]
gl/cre ic [ic]
gl/cre di [di]

idh=200
if ($hexist([idh]).eq.1) then
  hi/del [idh]
endif

i1=$sigma(int([i1]+0.5))
i2=$sigma(int([i2]+0.5))
do i=[nh1],[nh2]

  hi/file 30 [nb[i]]
  mess [nb[i]] is opend

  do j=[i1],[i2]

    idhj=[j]*100+10*3+[ncnt]
    
    if ($hexist([idh]).eq.0) then
      hi/pl //lun30/[idhj]
      nx=$hinfo([idhj],'xbins')
      xmin=$hinfo([idhj],'xmin')
      xmax=$hinfo([idhj],'xmax')
      if ([nx].eq.0) then
        hi/pl [idhj]
        read x
      endif
      1d [idh] ! [nx] [xmin] [xmax]
      ve/del xis,xi2s,si2s,xi,si,xi2,si2,sis,ni
      ve/cre ni([nx]) r
      ve/cre xis([nx]) r
      ve/cre xi2s([nx]) r
      ve/cre si2s([nx]) r
      ve/cre xi([nx]) r
      ve/cre xi2([nx]) r
      ve/cre si([nx]) r
      ve/cre si2([nx]) r
      ve/cre sis([nx]) r
    endif
    hi/get/cont [idhj] xi
    hi/get/err  [idhj] si
    sigma ni = ni + xi/xi
    sigma xi2 = xi**2
    sigma si2 = si**2
    sigma xis = xis + xi
    sigma xi2s = xi2s + xi2
    sigma si2s = si2s + si2
    hi/del [idhj]
  enddo
  close 30
enddo
sigma xis = xis/ni
sigma sis = sqrt(((si2s+xi2s)/ni-xis**2)/ni)
hi/put/cont [idh] xis
hi/put/erro [idh] sis
return



macro histprepf ncnt=1 nh1=2 nh2=5 f1=5.5 f2=15.5
nh=0
nh=[nh]+1; nb[nh]=mhad2010_profile.hbook
nh=[nh]+1; nb[nh]=mhad2011-4_profile.hbook
nh=[nh]+1; nb[nh]=mhad2012-2_profile.hbook
nh=[nh]+1; nb[nh]=phi2013-0_profile.hbook
nh=[nh]+1; nb[nh]=omeg2012-0_n1.13_profile.hbook
nh=[nh]+1; nb[nh]=omeg2012-0_n1.05_profile.hbook
nh=[nh]+1; nb[nh]=rho_2012-0_n1.05-1_profile.hbook
nh=[nh]+1; nb[nh]=rho_2012-0_n1.13_profile.hbook
nh=[nh]+1; nb[nh]=rho_2012-0_n1.05-2_profile.hbook
nh=[nh]+1; nb[nh]=sim_profile.hbook
*
nh=[nh]+1; nb[nh]=mapcal1_profilex_exp.his
nh=[nh]+1; nb[nh]=mapcal2_profilex_exp.his
nh=[nh]+1; nb[nh]=mapcal3_profilex_exp.his
nh=[nh]+1; nb[nh]=mapcal4_profilex_exp.his
*
nh=[nh]+1; nb[nh]=mapcal1_profilex_exp_pds.his
nh=[nh]+1; nb[nh]=mapcal2_profilex_exp_pds.his
nh=[nh]+1; nb[nh]=mapcal3_profilex_exp_pds.his
nh=[nh]+1; nb[nh]=mapcal4_profilex_exp_pds.his

mnf=320
mfmax=[ncnt]*40+20
mfmin=[ncnt]*40-60
i1=$sigma(int([mnf]*([f1]-[mfmin])/([mfmax]-[mfmin])+0.5))
i2=$sigma(int([mnf]*([f2]-[mfmin])/([mfmax]-[mfmin])+0.5))

idh=200
if ($hexist([idh]).eq.1) then
  hi/del [idh]
endif

do i=[nh1],[nh2]

  hi/file 30 [nb[i]]
  mess [nb[i]] is opend

  do j=[i1],[i2]

    idhj=[j]*100+[ncnt]
    
    if ($hexist([idh]).eq.0) then
	  mess 12345
      hi/pl //lun30/[idhj]
      nx=$hinfo([idhj],'xbins')
      xmin=$hinfo([idhj],'xmin')
      xmax=$hinfo([idhj],'xmax')
      if ([nx].eq.0) then
        hi/pl [idhj]
        read x
      endif
      1d [idh] ! [nx] [xmin] [xmax]
      ve/del xis,xi2s,si2s,xi,si,xi2,si2,sis,ni
      ve/cre ni([nx]) r
      ve/cre xis([nx]) r
      ve/cre xi2s([nx]) r
      ve/cre si2s([nx]) r
      ve/cre xi([nx]) r
      ve/cre xi2([nx]) r
      ve/cre si([nx]) r
      ve/cre si2([nx]) r
      ve/cre sis([nx]) r
    endif
    hi/get/cont [idhj] xi
    hi/get/err  [idhj] si
    sigma ni = ni + xi/xi
    sigma xi2 = xi**2
    sigma si2 = si**2
    sigma xis = xis + xi
    sigma xi2s = xi2s + xi2
    sigma si2s = si2s + si2
    hi/del [idhj]
  enddo
  close 30
enddo
sigma xis = xis/ni
sigma sis = sqrt(((si2s+xi2s)/ni-xis**2)/ni)
hi/put/cont [idh] xis
hi/put/erro [idh] sis
return



macro fitxsimz i1=1 i2=[ncz]
zmin=-10
zmax=9
ncz=5
dz=$sigma(([zmax]-[zmin])/[ncz])
do i=[i1],[i2]
  z2=$sigma([zmin]+[i]*[dz])
  z1=[z2]-[dz]
  exec mapcal#fitxsim z1=[z1] z2=[z2]
enddo
return

macro fitxsim z1=-10 z2=9 i1=1 i2=9
gl/cre ssm 0.01
gl/cre mname mhad2011-4
gl/cre expsim sim
gl/cre ver 102
ve/cre z1b(1) r [z1]
ve/cre z2b(1) r [z2]
ve/del dphis
ve/read dphis phi_shift_sim.txt
ve/del phisf
ve/cre phisf(9) r 9*0.01
do i=[i1],[i2]
  gl/cre dphisi $sigma(([i]-0.5)*40-2.5+dphis([i]))
  gl/cre dphisf $sigma(phisf([i]))
  exec mapcal#expsimf [i] [z1] [z2]
  gl/imp ic
  gl/imp di
  exec mapcal#fitxs [ic] [ic] 1 1 f [di] idh=300 ncnt=[i]
enddo
return

macro fitxsim0 z1=-10 z2=9
gl/cre ssm 0.01
gl/cre mname mhad2011-4
gl/cre expsim sim0
ve/cre z1b(1) r [z1]
ve/cre z2b(1) r [z2]
ve/del dphis
ve/cre dphis(9) r
*ve/read dphis phi_shift_sim.txt
ve/del phisf
ve/cre phisf(9) r 
do i=1,9
  gl/cre dphisi $sigma(([i]-0.5)*40-2.5+dphis([i]))
  gl/cre dphisf $sigma(phisf([i]))
  exec mapcal#expsimf [i] [z1] [z2]
  gl/imp ic
  gl/imp di
  exec mapcal#fitxs [ic] [ic] 1 1 f [di] idh=400 ncnt=[i]
enddo
return



macro fitxs ns1=100 ns2=200 dns=1 rd=0 pf='f' di=0 idh=0 ncnt=1

gl/imp mname
  
if ([rd].eq.1) then

  mess [ns1] [ns2] [dns] [rd]

  exec mapcal#cuts mylibinit='yes'

  mess [ns1] [ns2] [dns] [rd]

*  ncnt=1
  mnz=300
  mzmin=-15
  mzmax=15
  mnf=320
  mfmax=[ncnt]*40+20
  mfmin=[ncnt]*40-60
  dz=$sigma(([mzmax]-[mzmin])/[mnz])
  df=$sigma(([mfmax]-[mfmin])/[mnf])
  ve/cre xvi([mnz]) r
  ve/cre yvi([mnf]) r
  l0=[mzmin]+[dz]/2.
  r0=[mzmax]-[dz]/2.
  sigma xvi = array([mnz],[l0]#[r0])
  l0=[mfmin]+[df]/2.
  r0=[mfmax]-[df]/2.
  sigma yvi = array([mnf],[l0]#[r0])

  exp=mhad2011
*  exp=[mname]
  ebeam=all
  
  dir=[exp]_[ebeam]/counter_[ncnt]
  datafile=[exp]_[ebeam]_nc[ncnt].dat
  datafile=mapfitnew.txt

  nsc=[mnz]*[mnf]
  ve/del av1,dav1
*  ve/read av1,dav1 [dir]/[datafile] 2g15.6

  ve/cre  av([mnz],[mnf]) r
  ve/cre dav([mnz],[mnf]) r
*  ve/read av,dav [dir]/[datafile] 2g15.6

endif

ve/cre dxi(10) r
gl/imp expsim
gl/imp ver
*pref=mhad2011
pref=[mname]
if ([expsim].eq.'exp') then
*  pref=mhad2011
  pref=[mname]_exp_v[ver]
else
  if ([expsim].eq.'sim0') then
    pref=[mname]_sim0
  else
    pref=[mname]_sim_v[ver]
  endif
endif
opt nstat
do i=[ns1],[ns2],[dns]
  if ([pf].eq.'f') then
    i1=[i]-[di]
    i2=[i]+[di]
    if ([expsim].eq.'exp') then
      pref0=[pref]
    else
      pref0=[pref]
    endif
*    gl/cre fname [pref0]_amplitude_vs_phi_counter[ncnt]_slx[i]_fit_[di]_et0_chit.par
*    if ($fexist([fname]).eq.0) then
*      txt=File <[fname]> is not exists
*      gl/cre fname ''
*      mess [txt] 
*      read x
*    endif
*    gl/cre fname ''
    exec mapcal#fitx [i1] [i2] idh=[idh] ncnt=[ncnt]
    exec vappend  p16  xi
    exec vappend dp16 dxi
    ve/cre chi2r(1) r $sigma(chi2(1))
    exec vappend  p16 chi2r
    exec vappend dp16 ndf
    exec vappend  p16 z1b
    exec vappend dp16 z2b
    gl/cre fname [pref0]_amplitude_vs_phi_counter[ncnt]_slx[i]_fit_[di]_et0_chit.par
    ve/write p16,dp16 [fname] '2f15.10'
    set hwid 1
    set ksiz 0.05
    point 1000
    hi/pl 1000
*    goto 1
    l=$sigma(p16(12)-p16(17))
    r=$sigma(p16(15)+p16(17))
*    set ltyp 14
*    set basl 0.005
    set hcol 5
    ve/cre mod(4) r 0 0 0 1
    fun/pl scphip.f [l] [r] s
*    set ltyp 13
*    set basl 0.01
    set hcol 4
    ve/cre mod(4) r 1 0 0 0
    fun/pl scphip.f [l] [r] s
    ve/cre mod(4) r 0 1 0 0
    fun/pl scphip.f [l] [r] s
    ve/cre mod(4) r 0 0 1 0
    fun/pl scphip.f [l] [r] s
*    set ltyp 1
    set hcol 2
    ve/cre mod(4) r 1 1 1 1
    fun/pl scphip.f [l] [r] s
    set hcol 1
*    exec /personal/l3-1-205-1/konctbel/kskl/s#tf 0.1 0.9 [i]
1:
    exec pl#tf 0.1 0.9 [i] 
    atitle '[f], degree' 'A, pe.'
    exec save [pref]_amplitude_vs_phi_counter[ncnt]_slix[i]_fit_[di]_et0_chit.eps f
    gl/cre epsfile [pref]_amplitude_vs_phi_counter[ncnt]_slix[i]_fit_[di]_et0_chit.eps
  endif
  if ([pf].eq.'p') then
    ve/del p16,dp16
    ve/read p16,dp16 [mname]/[pref]_amplitude_vs_phi_counter[ncnt]_slx[i]_fit_[di]_et0_chit.par '2f15.10'
    mess [mname]/[pref]_amplitude_vs_phi_counter[ncnt]_slx[i]_fit_[di]_et0_chit.par '2f15.10'
    s1=p16(12)
    ss=p16(13)
    sm=p16(14)
    s2=p16(15)
    ve/cre xi(10) r [s1] $sigma([s1]+([ss]-[s1])/3) $sigma([s1]+([ss]-[s1])*2/3) _
                    [ss] $sigma([ss]+([sm]-[ss])/3) $sigma([ss]+([sm]-[ss])*2/3) _
                    [sm] $sigma([sm]+([s2]-[sm])/3) $sigma([sm]+([s2]-[sm])*2/3) [s2]
    i1=$sigma(int([s1]/15))
    i2=$sigma(int([s2]/15))
    ni=[i2]-[i1]
    i1=[i1]+1
    f1=[i1]*15
    f2=[i2]*15
    ve/cre si0([ni]) r
    sigma si0 = array([ni],[f1]#[f2])
    ve/cre si(3) r 3*-1000
    ve/copy si0(1:[ni]) si(1:[ni])
    set hwid 1
    set ksiz 0.05
    point 1000
    exec hsigma u = vmax( (@[idh]+%[idh]) * ( $[idh] gt [s1] and $[idh] lt [s2] ))
    u=u(1)
*    ve/cre mod(4) r 1 1 1 1
*    ve/cre xf(1) r [ss]
*    u=$call('scphip.f(xf)')
    xmin=$hinfo([idh],'xmin')
    xmax=$hinfo([idh],'xmax')
    null [xmin] [xmax] 0 $sigma(1.2*[u])
    hi/pl [idh] s
    l=$sigma(p16(12)-p16(17))
    r=$sigma(p16(15)+p16(17))
*    exec mapcal#vecut p16 29
*    set ltyp 14
*    set basl 0.005
    set hcol 3
    ve/cre mod(4) r 0 0 0 1
    fun/pl scphip.f [l] [r] s
*    set ltyp 13
*    set basl 0.01
    set hcol 4
    ve/cre mod(4) r 1 0 0 0
    fun/pl scphip.f [l] [r] s
    ve/cre mod(4) r 0 1 0 0
    fun/pl scphip.f [l] [r] s
    ve/cre mod(4) r 0 0 1 0
    fun/pl scphip.f [l] [r] s
*    set ltyp 1
    set hcol 2
    ve/cre mod(4) r 1 1 1 1
    fun/pl scphip.f [l] [r] s
    mess fun/pl scphip.f [l] [r] s
    if ($fexist(lthpl).eq.0) then
      ve/cre lthpl(1) r
    endif
    lthplx=lthpl(1)
    ve/cre lthpl(1) r 0
    set dmod 5
    fun/pl scphip.f [l] [r] s
    ve/cre lthpl(1) r [lthplx]
    set dmod 1
    set hcol 1
*    exec /personal/l3-1-205-1/konctbel/kskl/s#tf 0.1 0.9 [i] 
    exec pl#tf 0.1 0.9 [i] 
    atitle '[f], degree' 'A, pe.'
    exec save [pref]_amplitude_vs_phi_counter[ncnt]_slix[i]_fit_[di]_et0_chit.eps f
    gl/cre epsfile [pref]_amplitude_vs_phi_counter[ncnt]_slix[i]_fit_[di]_et0_chit.eps
  endif
enddo
return


macro fitxss
do i=1,30
  fname=fitxs[i].kumac
  if ($fexist([fname]).eq.1) then
    shell rm [fname]
  endif
  for/file 20 [fname] new
  close 20
  txt=exec mapcal#fitxs [i] [i] 1 1
  fmess [txt] [fname]
  shell .testrelease/.mainrelease/Offline/submit.sh -q clusters,180 pawbigX11 -b [fname]
enddo
return

macro wlsdtime ncnt=1 i1=1 i2=[n]
exec mapcal#slxparn ncnt=[ncnt] npar=16
ve/del pno,dpno,xvio,dxvio
sigma  pno = order(pn,-pn)
sigma dpno = order(dpn,-pn)
sigma  xvio = order(xvi,-pn)
sigma dxvio = order(dxvi,-pn)
n=$vlen(pno)
ve/del pnr,dpnr,xvir,dxvir
ve/copy  pno(1:[n])  pnr
ve/copy dpno(1:[n]) dpnr
ve/copy  xvio(1:[n])  xvir
ve/copy dxvio(1:[n]) dxvir
sigma  pnr = order( pnr,xvir)
sigma dpnr = order(dpnr,xvir)
sigma dxvir = order(dxvir,xvir)
sigma  xvir = order(xvir,xvir)
ve/del pnf,dpnf,xvif,dxvif
ve/copy  pnr([i1]:[i2])  pnf
ve/copy dpnr([i1]:[i2]) dpnf
ve/copy  xvir([i1]:[i2])  xvif
ve/copy dxvir([i1]:[i2]) dxvif
* exec $PER/s#vpl pn dpn xvi dxvi iatt=20 sz=0.1
*
*npar=2
*ve/cre chi2(2) r
*ve/cre paru([npar]) r
*ve/cre dparu([npar]) r
*ve/cre covu([npar],[npar]) r
*
ve/cre p2(2) r 3 0.02
ve/fit xvif pnf dpnf e s 2 p2
*call covm.f(1)
*call covmpen(chi2,[npar],paru,dparu)
*call covmcov([npar],covu)
*
return


macro slxparn ncnt=1 npar=1 pref=mhad2011 suff=.par
mnz=300
mzmin=-15
mzmax=15
mnf=320
mfmax=[ncnt]*40+20
mfmin=[ncnt]*40-60
dz=$sigma(([mzmax]-[mzmin])/[mnz])
df=$sigma(([mfmax]-[mfmin])/[mnf])
ve/cre xvi([mnz]) r
ve/cre yvi([mnf]) r
l0=[mzmin]+[dz]/2.
r0=[mzmax]-[dz]/2.
sigma xvi = array([mnz],[l0]#[r0])
l0=[mfmin]+[df]/2.
r0=[mfmax]-[df]/2.
sigma yvi = array([mnf],[l0]#[r0])

ve/cre dxvi([mnz]) r

ve/cre pn(300) r 
ve/cre dpn(300) r
do i=1,300
*  fname=mhad2011_amplitude_vs_phi_counter[ncnt]_slx[i]_fit_8.par
  fname=[pref]_amplitude_vs_phi_counter[ncnt]_slx[i]_fit_8[suff]
  if ($fexist([fname]).eq.1) then
    ve/del p16,dp16
    ve/read p16,dp16 [fname] '2f15.10'
    ve/inp  pn([i])  p16([npar])
    ve/inp dpn([i]) dp16([npar])
  endif
enddo

* exec $PER/s#vpl pn dpn xvi dxvi sz=0.1

return


macro slzsim npar=1
n=0
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter1_slx69_fit_19_et0.par     
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter1_slx107_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter1_slx145_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter1_slx183_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter1_slx221_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter2_slx69_fit_19_et0.par     
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter2_slx107_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter2_slx145_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter2_slx183_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter2_slx221_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter3_slx69_fit_19_et0.par     
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter3_slx107_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter3_slx145_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter3_slx183_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter3_slx221_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter4_slx69_fit_19_et0.par     
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter4_slx107_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter4_slx145_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter4_slx183_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter4_slx221_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter5_slx69_fit_19_et0.par     
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter5_slx107_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter5_slx145_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter5_slx183_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter5_slx221_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter6_slx69_fit_19_et0.par     
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter6_slx107_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter6_slx145_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter6_slx183_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter6_slx221_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter7_slx69_fit_19_et0.par     
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter7_slx107_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter7_slx145_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter7_slx183_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter7_slx221_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter8_slx69_fit_19_et0.par     
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter8_slx107_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter8_slx145_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter8_slx183_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter8_slx221_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter9_slx69_fit_19_et0.par     
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter9_slx107_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter9_slx145_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter9_slx183_fit_19_et0.par    
n=[n]+1; fpar[n]=mhad2011_sim_amplitude_vs_phi_counter9_slx221_fit_19_et0.par    

mnz=300
mzmin=-15
mzmax=15
mnf=320
mfmax=[ncnt]*40+20
mfmin=[ncnt]*40-60
dz=$sigma(([mzmax]-[mzmin])/[mnz])
df=$sigma(([mfmax]-[mfmin])/[mnf])
ve/cre xvi([mnz]) r
ve/cre yvi([mnf]) r
l0=[mzmin]+[dz]/2.
r0=[mzmax]-[dz]/2.
sigma xvi = array([mnz],[l0]#[r0])
l0=[mfmin]+[df]/2.
r0=[mfmax]-[df]/2.
sigma yvi = array([mnf],[l0]#[r0])

ve/cre dxvi([mnz]) r

ve/cre pn([n]) r 
ve/cre dpn([n]) r
do i=1,[n]
  fname=[fpar[i]]
  if ($fexist([fname]).eq.1) then
    ve/del p16,dp16
    ve/read p16,dp16 [fname] '2f15.10'
    ve/inp  pn([i])  p16([npar])
    ve/inp dpn([i]) dp16([npar])
  endif
enddo

ve/cre an([n]) r
ve/cre dan([n]) r
sigma an = array([n],1#[n])

* exec $PER/s#vpl pn dpn an dan sz=0.1

return


macro errn npar=1
exec mapcal#slxparn [npar]
if ($vlen(dpn).ne.0) then
ns=$vdim(dpn)
ve/cre ni([ns]) r
sigma ni = array([ns],1#[ns])
sigma dpn = abs(dpn)
sigma dpno = order(dpn,-dpn)
sigma nio = order(ni,-dpn)
n=$vlen(dpno)
ve/del dpnor, nior
ve/copy dpno(1:[n]) dpnor
ve/copy nio(1:[n]) nior
sigma nior = order(nior,dpnor)
sigma dpnor = order(dpnor,dpnor)
m=$sigma(vsum(dpnor)/[n])
s=$sigma(sqrt(vsum((dpnor-[m])**2))/[n])
r=$sigma([m]+3*[s])
mess [n] [m] [s]
exec hsigma dpnorr = dpnor*(dpnor lt [r])
nr=$vlen(dpnorr)
mess [n] [nr]
ve/del niorr[npar]
ve/copy nior([nr]:[n]) niorr[npar]
ve/prin niorr[npar]
endif
return

macro errnn
do i=1,-28
  exec mapcal#errn [i]
enddo
ve/cre nfn(300) r
ve/cre gni(7) r 2 3 6 79 10 13 16
do i=1,300
  fn=0
  do j=1,$vlen(gni)
    ng=gni([j])
    vp=niorr[ng]
    if ($vexist([vp]).eq.1) then
      do l=1,$vlen([vp])
        if ($sigma([vp]([l])).eq.[i]) then
          fn=1
        endif
      enddo
    endif
  enddo
  ve/inp nfn([i]) [fn]
enddo
return

macro corra icorra=1
do nc=1,9
  exec mapcal#effz [nc] 6 1
  ve/inp  corraf([nc]) $sigma(corra(1)) 
  ve/inp dcorraf([nc]) $sigma(dcorra(1))
  exec mapcal#effz [nc] 7 1
  ve/inp  corraz([nc]) $sigma(corra(1)) 
  ve/inp dcorraz([nc]) $sigma(dcorra(1))
enddo
ve/write corraf,dcorraf,corraz,dcorraz corra_v[icorra].txt '(4f15.6)'
do v=$sigma([icorra]-1),1,-1
  ve/del icorraf,idcorraf,icorraz,idcorraz
  ve/read icorraf,idcorraf,icorraz,idcorraz corra_v[v].txt '(4f15.6)'
  sigma corraf = corraf*icorraf
  sigma corraz = corraz*icorraz
enddo
ve/write corraf,dcorraf,corraz,dcorraz corra_v[icorra]_tot.txt '(4f15.6)'
return


macro rfdir mname=mapcal1 nc=1
s='"'
l='|'
pref=[mname]_exp_v0
sufix=et0_chit.par
shell rm tmp[0-9].txt
shell ls -l [mname]/[mname]*counter[nc]*.par > tmp0.txt
shell $unquote('cat tmp0.txt | sed "s/\(^.*slx\)/ /g" > tmp1.txt')
txt=cat tmp1.txt [l] sed [s]s/\([sufix]*$\)/ /g[s] > tmp2.txt
shell $unquote([txt])
shell $unquote('cat tmp2.txt | sed "s/[^0-9]/ /g" > tmp.txt')
ve/del is,dis
ve/read is,dis tmp.txt
isi=is(8)
disi=dis(8)
shell cp [mname]/[mname]*counter[nc]*[isi]*[disi]*.par tmp3.txt
shell diff [mname]/[mname]*counter[nc]*[isi]*[disi]*.par tmp3.txt
return



macro mapfitprerp map=1 ver=0
*
exec mapcal#cuts
*
if ([map].eq.1) then
  gl/cre mname mapcal1
  nb=14
endif
if ([map].eq.2) then
  gl/cre mname mapcal2
  nb=14
endif
if ([map].eq.3) then
  gl/cre mname mapcal3
  nb=10
endif
if ([map].eq.4) then
  gl/cre mname mapcal4
  nb=14
  ve/del parc,dparc
  ve/read parc,dparc mapcal4/mapcal4_exp_v0_amplitude_vs_phi_counter9_slx109_fit_32_et0_chit.par '2f15.10'
endif
*
sh=0
gl/cre ver [ver]
fname=[mname]_v[ver].cal
if ($fexist([fname]).eq.1) then
  shell rm [fname]
endif
for/file 20 [fname] new
close 20
*
ve/cre pix(3) r 0.001 1 0
ve/cre si(3) r 3*-1000
ve/cre mod(4) r 4*1
point 1000
*
ve/cre vi1(9) r 1 1 3 1 1 1 1 1 1
ve/cre vi2(9) r 9 9 9 8 10 10 9 9 7
ve/cre vif(9) r 10 12 8 10 10 10 10 10 9
ve/cre vim(9) r 12 14 12 14 14 14 14 14 13
*
ve/cre xis(10) r 
ve/cre xii(10) r
ve/cre yip(4) r
ve/cre xic(4) r
ve/cre yic(4) r
ve/cre xt(1) r
*
pref=[mname]_exp_v0
sufix=et0_chit.par
*
ve/del acorr
ve/read acorr p0_amplitude.dat
icorra=0
ve/cre corraf(9) r 9*1
ve/cre corraz(9) r 9*1
ve/cre dcorraf(9) r 
ve/cre dcorraz(9) r
gl/cre corrz $sigma(min([ver],1))
*
ve/del corra
if ([ver].eq.2) then
  ve/read corra AmplitudeCorrection_mapcal[map].txt
else
  ve/read corra AmplitudeCorrection_mapcal0.txt
endif
*
do nc=1,9
*
  if ([corrz].ne.0) then
*    mapcal#zslixcorrn [map] [nc] 
  endif
  txt=Counter [nc]
  fmess [txt] [fname]
*  
  if ([corrz].eq.10) then
    tmp=[corrz]-1
    tmp=0
    ve/read p3,dp3 corrwls_counter[nc]_v[tmp].txt '2e15.6'
    txt=Shifter $sigma(p3(1)) $sigma(p3(2)) $sigma(p3(3))
  else
    ve/del p3,dp3
    ve/read p3,dp3 [mname]/[mname]_corrwls_counter[nc]_v0.txt '2e15.6'
    txt=Shifter $sigma(p3(1)*corra([nc])) $sigma(p3(2)*corra([nc])) $sigma(p3(3)*corra([nc]))
  endif
  fmess [txt] [fname]
*  
  if ([icorra].ne.0) then
    exec mapcal#effz [nc] 6 1
    ve/inp  corraf([nc]) $sigma(corra(1)) 
    ve/inp dcorraf([nc]) $sigma(dcorra(1))
    corrph=corra(1)
    exec mapcal#effz [nc] 7 1
    ve/inp  corraz([nc]) $sigma(corra(1)) 
    ve/inp dcorraz([nc]) $sigma(dcorra(1)) 
    corrzd=corra(1)
    corra=([corrph]+[corrzd])/2.0
    ve/write corraf,dcorraf,corraz,dcorraz corra_v[icorra].txt '(4f15.6)'
  else
    corra=1
  endif
  txt=AmplitudeCorrection [corra]
*  txt=AmplitudeCorrection 1
  fmess [txt] [fname]
*  
  txt=Aerogel
  fmess [txt] [fname]
*    
  s='"'
  l='|'
  shell rm tmp[0-9].txt
  shell ls -l [mname]/[mname]*counter[nc]*.par > tmp0.txt
  shell $unquote('cat tmp0.txt | sed "s/\(^.*slx\)/ /g" > tmp1.txt')
  txt=cat tmp1.txt [l] sed [s]s/\([sufix]*$\)/ /g[s] > tmp2.txt
  shell $unquote([txt])
  shell $unquote('cat tmp2.txt | sed "s/[^0-9]/ /g" > tmp.txt')
  ve/del is,dis
  ve/read is,dis tmp.txt
*    
  nbi=$vlen(is)
  if ([nbi].ge.[nb]) then
    sigma is=order(is,dis)
    sigma dis=order(dis,dis)
    exec mapcal#vecut is [nb]
    exec mapcal#vecut dis [nb]
    sigma dis=order(dis,is)
    sigma is=order(is,is)
    do i=1,[nb]
      isi=is([i])
      disi=dis([i])
      shell cp [mname]/[mname]*counter[nc]*slx[isi]_fit_[disi]*.par block[i].txt
    enddo
    n=[nb]
  else
    mess Number of blocks files is less than expected!!!
    read x
  endif
*    
  do i=1,[nb]
    ve/del par[i],dpar[i]
    ve/read par[i],dpar[i] block[i].txt '2f15.10' 
  enddo
*
  ve/del xis
  ve/cre xis(10) r 
  do i=1,[nb]
    ve/copy par[i](30:39) xii
    sigma xis = xis + xii
  enddo
  sigma xis = xis/[nb]
*    
  do i=1,10
    dx[i]=$sigma(xis([i])-[sh]*xis(4))
  enddo
*
  do nl=1,3
*  
    if ([nl].eq.1) then
      txt=[n] [dx1] [dx2] [dx3] [dx4]      
    endif
    if ([nl].eq.2) then
      txt=[n] [dx4] [dx5] [dx6] [dx7]
    endif
    if ([nl].eq.3) then
      txt=[n] [dx7] [dx8] [dx9] [dx10]
    endif
    fmess [txt] [fname]
*    
    do i=1,[nb]
    
      if (([map].eq.4).and.([nc].eq.9).and.([i].ge.4).and.([i].le.7)) then
        vpari=parc
      else
        vpari=par[i]
      endif
      
      ve/copy [vpari](30:39) xii
      do j=1,10
        x[j]=xii([j])
      enddo
      do ii=1,11
        ve/inp [vpari]([ii]) $sigma(abs([vpari]([ii])))
      enddo
      ii=29;  ve/inp [vpari]([ii]) $sigma(abs([vpari]([ii])))
      if ([nl].eq.1) then
        ve/copy [vpari](1:4) yip(1:4)
        ve/del xip
        ve/cre xip(4) r [x1] [x2] [x3] [x4]
        do j=1,4
          ve/inp xt(1) $sigma(xis([j]))
          ve/inp xic([j]) $sigma(xis([j]))
          y[j]=$call('f3p.f(xt)')
          mess [x[j]] [y[j]] $sigma(xis([j])) $sigma(yip([j]))
          ve/inp yic([j]) [y[j]]
        enddo
      endif
      if ([nl].eq.2) then
        ve/copy [vpari](5:8) yip(1:4)
        ve/del xip
        ve/cre xip(4) r [x4] [x5] [x6] [x7]
        do j=1,4
          jj=[j]+3
          ve/inp xt(1) $sigma(xis([jj]))
          ve/inp xic([j]) $sigma(xis([jj])) 
          y[j]=$call('f3p.f(xt)')
          mess [x[jj]] [y[j]] $sigma(xis([jj])) $sigma(yip([j]))
          ve/inp yic([j]) [y[j]]
        enddo
      endif
      if ([nl].eq.3) then
        ve/copy [vpari](9:11) yip(2:4)
        ve/copy [vpari](29) yip(1)
        ve/del xip
        ve/cre xip(4) r [x7] [x8] [x9] [x10]
        do j=1,4
          jj=[j]+6
          ve/inp xt(1) $sigma(xis([jj]))
          ve/inp xic([j]) $sigma(xis([jj])) 
          y[j]=$call('f3p.f(xt)')
          mess [x[jj]] [y[j]] $sigma(xis([jj])) $sigma(yip([j]))
          ve/inp yic([j]) [y[j]]
        enddo
      endif
      z1=$sigma(par[i](41))
      z2=$sigma(dpar[i](41))
      if ([corrz].ne.0) then
*        corrslx=corrslix[nc]([i])
*        y1=[y1]*[corrslx]
*        y2=[y2]*[corrslx]
*        y3=[y3]*[corrslx]
*        y4=[y4]*[corrslx]
        do v=1,[corrz]
          ve/del corrslix,dcorrslix,cslxp0,cslxp1,cslxp2
          ve/read corrslix,dcorrslix,cslxp0,cslxp1,cslxp2 [mname]/[mname]_corrslix_counter[nc]_v[v].txt '5e15.6'
          acsx=cslxp0([i])
          bcsx=cslxp1([i])
          ccsx=cslxp2([i])
          corrslx=$sigma(([acsx])+([bcsx])*xip(1)+([ccsx])*xip(1)**2)
          y1=[y1]*[corrslx]
          corrslx=$sigma(([acsx])+([bcsx])*xip(2)+([ccsx])*xip(2)**2)
          y2=[y2]*[corrslx]
          corrslx=$sigma(([acsx])+([bcsx])*xip(3)+([ccsx])*xip(3)**2)
          y3=[y3]*[corrslx]
          corrslx=$sigma(([acsx])+([bcsx])*xip(4)+([ccsx])*xip(4)**2)
          y4=[y4]*[corrslx]
        enddo
      endif
      do v=1,-[icorra]
        ve/del icorraf,idcorraf,icorraz,idcorraz
        ve/read icorraf,idcorraf,icorraz,idcorraz corra_v[v].txt '(4f15.6)'
        corra=$sigma((icorraf([nc])+icorraz([nc]))/2)
        y1=[y1]*[corra]
        y2=[y2]*[corra]
        y3=[y3]*[corra]
        y4=[y4]*[corra]
      enddo
*      y1=0
*      y2=0
*      y3=0
*      y4=0
      y1=$sigma([y1]*corra([nc]))
      y2=$sigma([y2]*corra([nc]))
      y3=$sigma([y3]*corra([nc]))
      y4=$sigma([y4]*corra([nc]))
      txt=[z1] [y1] [y2] [y3] [y4]
      fmess [txt] [fname]
*      
      ve/copy par[i](1:29) p16
      ve/copy par[i](30:39) xi
      set hcol 1
*      fun/pl scphip.f(x) $sigma(vmin(xi)-3) $sigma(vmax(xi)+3)
      fun/pl scphip.f(x) $sigma([nc]*40-50) $sigma([nc]*40+10)
      set hcol 2
      l=$sigma(xip(1)-5)
      r=$sigma(xip(4)+5)
      fun/pl f3p.f(x) [l] [r] s
      graph 4 xic yic *
*      read x
    enddo
  enddo
enddo
exec mapcal#readmap [fname]
return



macro readmap fname=mapcal2011-4shifter0_v2.cal
ve/cre data(2054) r
n=0
n=[n]+1; par[n]=Counter
n=[n]+1; par[n]=Shifter
n=[n]+1; par[n]=AmplitudeCorrection
ve/cre par(100) r
ve/del counter,shifter,ampcorr
if ($fexist([fname]).eq.1) then
  exec mapcal#readspectp [fname] Counter
  ve/copy par(1:9) counter
  exec mapcal#readspectp [fname] Shifter
  ve/copy par(1:27) shifter
  exec mapcal#readspectp [fname] AmplitudeCorrection
  ve/copy par(1:9) ampcorr
  ve/inp data(1) 1
  ve/inp data(2) 1.13
  ind0=0
  do i=1,9
  
    do j=1,3
      ind=$sigma(2+228*([i]-1)+[j])
      if ([j].eq.3) then
        ve/inp data([ind]) $sigma(shifter(3*([i]-1)+[j]))
      else
        ve/inp data([ind]) $sigma(shifter(3*([i]-1)+[j])*ampcorr([i]))
      endif
    enddo

    ind0=[ind0]+4
    
    do p=1,3
      ind0=[ind0]+1
      exec mapcal#readspectn [fname] [ind0]
      do j=1,5
        ind=$sigma(2+228*([i]-1)+3+5*15*([p]-1)+[j])
        ve/inp data([ind]) $sigma(par([j]))
      enddo
      nz=par(1)
      do l=1,[nz]
        ind0=[ind0]+1
        exec mapcal#readspectn [fname] [ind0]
        do j=1,5
          ind=$sigma(2+228*([i]-1)+3+5*15*([p]-1)+5*[l]+[j])
          if ([j].eq.1) then
            ve/inp data([ind]) $sigma(par([j]))
          else
            ve/inp data([ind]) $sigma(par([j])*ampcorr([i]))
          endif
        enddo
      enddo
    enddo
    
  enddo
endif
*
ve/write data /work/users/konctbel/snd2k/R005-999/[fname].v
return


macro readspectp fname=x par=Counter
s='"'
l='|'
*ve/del par
sigma par = par*0
if ($fexist([fname]).eq.1) then
  shell fgrep -e [s][par] [s] [fname] > out.txt
  txt=cat out.txt [l] sed [s]s/[par]/ /g[s] > out1.txt
  shell $unquote([txt])
*  mess shell $unquote([txt])
  ve/read par out1.txt
endif
return

macro readspectn fname=x line=1
s='"'
l='|'
*ve/del par
sigma par = par*0
if ($fexist([fname]).eq.1) then
  txt=cat [fname] [l] sed -n [s][line]p[s] > out.txt
  shell $unquote([txt])
*  mess shell $unquote([txt])
  ve/read par out.txt
endif
return


macro wlscorrm
ve/cre  ab(9) r
ve/cre dab(9) r
ve/cre  af(9) r
ve/cre daf(9) r
ve/cre  tau(9) r
ve/cre dtau(9) r
corrz=1
xmin=-10
xmax=10
null [xmin] [xmax] 0 100
do i=1,9
  ve/del p3,dp3
  ve/read p3,dp3 corrwls_counter[i]_v[corrz].txt '2e15.6'
  ve/inp  ab([i]) $sigma(p3(1))
  ve/inp dab([i]) $sigma(dp3(1))
  ve/inp  af([i]) $sigma(p3(2))
  ve/inp daf([i]) $sigma(dp3(2))
  ve/inp  tau([i]) $sigma(abs(p3(3)))
  ve/inp dtau([i]) $sigma(dp3(3))
*  
  a2=p3(1)
  b2=p3(2)
  c2=p3(3)
  set hcol 4
  fun/pl ([a2])*exp(-x/([c2]))+([b2])*exp(x/([c2])) [xmin] [xmax] s
enddo
read x
ve/cre ni(9) r 1 2 3 4 5 6 7 8 9
ve/cre dni(9) r
ve/cre p0(1) r
* exec $PER/s#vpl ab dab ni dni
ve/fit ni ab dab p0 s 1 p0
abm=p0(1)
read x
* exec $PER/s#vpl af daf ni dni
ve/fit ni af daf p0 s 1 p0
afm=p0(1)
read x
* exec $PER/s#vpl tau dtau ni dni
ve/fit ni tau dtau p0 s 1 p0
taum=p0(1)
mess ab=[abm] af=[afm] tau=[taum]
return


macro zslixcorrx
do i=1,9
  exec mapcal#zslixcorr [i]
enddo
return

macro zslixcorr ncnt=1
gl/imp corrz
*
ver=[corrz]-1
shell fgrep -e Shifter mapcal2011-4shifter_v[ver].cal >& tmp.txt
shell $unquote('cat tmp.txt | sed "s/Shifter/ /g" > tmp1.txt')
ve/del ab,af,tau
ve/read ab,af,tau tmp1.txt
ve/cre wls(3) r $sigma(ab([ncnt])) $sigma(af([ncnt])) $sigma(tau([ncnt]))
*
fname=mhad2011-4_amplitude_vs_phi_counter[ncnt]_slx145_fit_95_et0_chit.par
ve/del pars,dpars
ve/read pars,dpars [fname] '2f15.10'
ve/cre  cwls(14) r
ve/cre dcwls(14) r
ve/cre  zwls(14) r
ve/cre dzwls(14) r
*
ve/del ic,if1,if2,if3,if4,if5
ve/read ic,if1,if2,if3,if4,if5 phi_cuts.txt
f1=if1([ncnt])
f2=if2([ncnt])
f3=if3([ncnt])
f4=if5([ncnt])
mnf=320
mfmax=[ncnt]*40+20
mfmin=[ncnt]*40-60
i1=$sigma(int([mnf]*([f1]-[mfmin])/([mfmax]-[mfmin])+0.5))
i2=$sigma(int([mnf]*([f2]-[mfmin])/([mfmax]-[mfmin])+0.5))
i3=$sigma(int([mnf]*([f3]-[mfmin])/([mfmax]-[mfmin])+0.5))
i4=$sigma(int([mnf]*([f4]-[mfmin])/([mfmax]-[mfmin])+0.5))
ve/cre cf([mnf]) r
do i=[i1],[i2]
  ve/inp cf([i]) 1
enddo
do i=[i3],[i4]
  ve/inp cf([i]) 1
enddo
*ve/cre cf([mnf]) r
*do i=[i1],[i4]
*  ve/inp cf([i]) 1
*enddo
*
gl/imp ncnto [ncnt]
gl/cre mname mhad2011-4
fzic=zic[ncnt].txt
if ($fexist([fzic]).eq.0) then
  exec mapcal#fitpl 110[ncnt] 1
  ve/copy zic zic1
  exec mapcal#fitpl 210[ncnt] 1
  ve/copy zic zic2
  exec mapcal#fitpl 310[ncnt] 1
  ve/copy zic zic3
  sigma zic = (zic1+zic2+zic3)/3
  ve/write zic [fzic] 
else
  ve/del zic
  ve/read zic [fzic]
endif
n=$vlen(zic)
n=[n]-1
ve/cre  corrslix[ncnt]([n]) r
ve/cre dcorrslix[ncnt]([n]) r
ve/cre cslx[ncnt]p0([n]) r
ve/cre cslx[ncnt]p1([n]) r
ve/cre cslx[ncnt]p2([n]) r
ve/cre p4(4) r 1 0 $sigma(pars(13)) $sigma(pars(17))
ve/cre dp4(4) r
ve/cre s4(4) r 0.01 0 0 0
ve/cre p6(6) r 1 0 0 0 $sigma(pars(13)) $sigma(pars(17))
ve/cre s6(6) r 0.01 0.0001 0.00001 0 0 0
ve/cre dp6(6) r
ve/cre ashe(14) r
ve/cre dashe(14) r
ve/cre ashs(14) r
ve/cre dashs(14) r
zz=1
do i=1,[n]
  if ([zz].eq.1) then
  exec mapcal#expsimf ncnt=[ncnt] ! ! slix=[i]
*  read x
  nx=$hinfo(200,'xbins')
  xmin=$hinfo(200,'xmin')
  xmax=$hinfo(200,'xmax')
  1d 400 ! [nx] [xmin] [xmax]
  ve/del v1,v2,dv1,dv2
  ve/cre v1([nx]) r
  ve/cre v2([nx]) r
  ve/cre dv1([nx]) r
  ve/cre dv2([nx]) r
  hi/get/cont 200 v1
  hi/get/cont 300 v2
  hi/get/err 200 dv1
  hi/get/err 300 dv2
  sigma vr = cf*v1/v2
  sigma dvr = cf*vr*sqrt((dv1/v1)**2+(dv2/v2)**2)
  hi/put/cont 400 vr
  hi/put/err 400 dvr
  hi/pl 400
  hi/fit 400 p0+g sb 4 p4 s4 ! ! dp4
  ve/inp  corrslix[ncnt]([i]) $sigma(p4(1))
  ve/inp dcorrslix[ncnt]([i]) $sigma(dp4(1))
  do j=1,5
    hi/fit 400 p2+g sb 6 p6 s6 ! ! dp6
  enddo
  ve/inp cslx[ncnt]p0([i]) $sigma(p6(1))
  ve/inp cslx[ncnt]p1([i]) $sigma(p6(2))
  ve/inp cslx[ncnt]p2([i]) $sigma(p6(3))
  a=p6(1)
  b=p6(2)
  c=p6(3)
  x=p6(5)
  f0=$sigma([a]+([b])*([x])+([c])*([x])**2)
  ve/inp cwls([i]) $sigma(p6(4)+[f0])
  ve/inp dcwls([i]) $sigma(dp6(4))
  endif
  if ([zz].eq.2) then
  j=[i]+1
  ve/inp zwls([i]) $sigma((zic([i])+zic([j]))/2)
*  
*  exec mapcal#slixi [ncnt] [i] mhad2011-4
  exec mapcal#slixi [ncnt] [i] mhad2011-4_exp_v1
  gl/imp fslix
  ve/del pars,dpars
  ve/read pars,dpars [fslix] '2f15.10'
  aexp=pars(16)
  daexp=dpars(16)
*  exec mapcal#slixi [ncnt] [i] mhad2011-4_sim_v[ver]
  exec mapcal#slixi [ncnt] [i] mhad2011-4_sim_v101
  gl/imp fslix
  ve/del pars,dpars
  ve/read pars,dpars [fslix] '2f15.10'
  asim=pars(16)
  dasim=dpars(16)
  rt=$sigma([aexp]/[asim])
  drt=$sigma([rt]*sqrt(([daexp]/[aexp])**2+([dasim]/[asim])**2))
  ve/inp cwls([i]) [rt]
  ve/inp dcwls([i]) [drt]
  ve/inp ashe([i]) [aexp]
  ve/inp dashe([i]) [daexp]
  ve/inp ashs([i]) [asim]
  ve/inp dashs([i]) [dasim]
  endif
  read x
enddo
gl/imp corrz
ve/write corrslix[ncnt],dcorrslix[ncnt],cslx[ncnt]p0,cslx[ncnt]p1,cslx[ncnt]p2 corrslix_counter[ncnt]_v[corrz].txt '5e15.6'
if ([zz].eq.2) then
*
xmin=$sigma(vmin(zic))
xmax=$sigma(vmax(zic))
*
null -12 12 0 $sigma(1.1*max(vmax(ashe),vmax(ashs)))
* exec $PER/s#vpl ashe dashe zwls dzwls o=s iatt=20
* exec $PER/s#vpl ashs dashs zwls dzwls o=s iatt=24
read x
*
* exec $PER/s#vpl cwls dcwls zwls dzwls
ve/del p3(3)
ve/copy wls p3
ve/cre dp3(3) r
n1=2
n2=9
ve/del cwlsr,dcwlsr,zwlsr
ve/copy  cwls([n1]:[n2])  cwlsr
ve/copy dcwls([n1]:[n2]) dcwlsr
ve/copy  zwls([n1]:[n2])  zwlsr
ve/fit zwlsr cwlsr dcwlsr wlscorr.f s 3 p3 ! ! ! dp3
sigma p3 = abs(p3)
read x
set hcol 2
*
a1=wls(1)
b1=wls(2)
c1=wls(3)
x=[xmin]
y1=$sigma(([a1])*exp(-([x])/([c1]))+([b1])*exp([x]/([c1])))
x=[xmax]
y2=$sigma(([a1])*exp(-([x])/([c1]))+([b1])*exp([x]/([c1])))
*
a2=p3(1)
b2=p3(2)
c2=p3(3)
x=[xmin]
y3=$sigma(([a2])*exp(-([x])/([c2]))+([b2])*exp([x]/([c2])))
x=[xmax]
y4=$sigma(([a2])*exp(-([x])/([c2]))+([b2])*exp([x]/([c2])))
*
ymax=$sigma(1.1*max(max([y1],[y2]),max([y3],[y4])))
null [xmin] [xmax] 0 [ymax]
set hcol 2
fun/pl ([a1])*exp(-x/([c1]))+([b1])*exp(x/([c1])) [xmin] [xmax] s
set hcol 4
fun/pl ([a2])*exp(-x/([c2]))+([b2])*exp(x/([c2])) [xmin] [xmax] s
set hcol 1
*
* exec $PER/s#vpl ashe dashe zwls dzwls o=s iatt=20
* exec $PER/s#vpl ashs dashs zwls dzwls o=s iatt=24
*
ve/write p3,dp3 corrwls_counter[ncnt]_v[corrz].txt '2e15.6'
read x
*
ve/cre ai([n]) r
ve/cre dai([n]) r
ve/cre ae([n]) r
ve/cre dae([n]) r
ve/cre as([n]) r
ve/cre das([n]) r
a=wls(1)
b=wls(2)
c=wls(3)
do i=1,[n]
  zi=zwls([i])
  fi=$sigma(([a])*exp(-([zi])/([c])))
  fi=$sigma([fi]+([b])*exp([zi]/([c])))
  sti=$sigma(12/sqrt(12**2+([zi])**2))
  ve/inp ai([i]) $sigma([fi]/[sti])
  ve/inp ae([i]) $sigma(ashe([i])/[sti])
  ve/inp dae([i]) $sigma(dashe([i])/[sti])
  ve/inp as([i]) $sigma(ashs([i])/[sti])
  ve/inp das([i]) $sigma(dashs([i])/[sti])
enddo
amax = $sigma(1.2*vmax(ai))
null 0 [amax] 0 [amax]
* exec $PER/s#vpl as das ai dai o=s sz=0.05
set plci 4
line 0 0 [amax] [amax]
ve/cre p3(3) r 25 0.3 5
ve/cre p3min(3) r 0 0.1 3 
ve/cre p3max(3) r 100 1 100
ve/cre p3stp(3) r 0.1 0.1 0.1
ve/fit ai as das ampsat.f sb 3 p3 p3stp p3min p3max
n1=1
n2=$sigma(int([n]/2))
ve/del yr,dyr,xr,dxr
ve/copy as([n1]:[n2]) yr
ve/copy das([n1]:[n2]) dyr
ve/copy ai([n1]:[n2]) xr
ve/copy dai([n1]:[n2]) dxr
set pmci 4
* exec $PER/s#vpl yr dyr xr dxr o=s sz=0.05
set pmci 1
set hcol 2
fun/pl ampsatp.f 0 [amax] s
*
ve/del aif,daif,asf,dasf
ve/read aif,daif,asf,dasf ampsatsim_counter2.txt '4f15.6'
ve/del aix,daix,asx,dasx
ve/read aix,daix,asx,dasx ampsatsim_counter[ncnt].txt '4f15.6'
set pmci 1
set hcol 1
* exec $PER/s#vpl asx dasx aix daix o=s sz=0.01
read x
*
sigma p3=abs(p3)
ar=p3(1)
am=p3(2)
at=p3(3)
ve/cre ci1(14) r 9*1 5*2
ve/cre ci2(14) r 11*1 3*2
ve/cre ci3(14) r 1 1 1 1 1 1 1 1 1 2 2 2 2 2
ve/cre ci4(14) r 10*1 4*2
ve/cre ci5(14) r 14*1
ve/cre ci6(14) r 10*1 4*2
ve/cre ci7(14) r 12*1 2*2
ve/cre ci8(14) r 10*1 4*2
ve/cre ci9(14) r 9*1 5*2
ve/cre aec([n]) r
ve/cre daec([n]) r
ve/cre ax(1) r
do i=1,[n]
  aei=ae([i])
  daei=dae([i])
  if ($sigma(ci[ncnt]([i])).eq.1) then
    aeci = [aei]
    daeci = [daei]
  else
    aeci = $sigma([ar]-[at]*log(([aei]/[ar]-[am])/(1-[am])))
    daeci = $sigma([at]/abs([aei]-[am]*[ar])*[daei])
  endif
  ve/inp aec([i]) [aeci]
  ve/inp daec([i]) [daeci]
  ve/inp ax(1) [aeci]
  mess [i] [aei] [aeci] $call('ampsatp.f(ax)')
enddo
*sigma  aec = [ar]-[at]*log((ae-[am])/([ar]-[am]))
*sigma daec = [at]/abs(ae-[am])*dae
* exec $PER/s#vpl ae dae aec daec o=s iatt=24 sz=0.05
read x 123
*
* exec $PER/s#vpl ai dai zwls dzwls iatt=23
* exec $PER/s#vpl ae dae zwls dzwls iatt=20 o=s
* exec $PER/s#vpl aec daec zwls dzwls iatt=24 o=s
*
do i=1,[n]
  aeci=aec([i])
  daeci=daec([i])
  zi=zwls([i])
  sti=$sigma(12/sqrt(12**2+([zi])**2))
  ve/inp aec([i]) $sigma([aeci]*[sti])
  ve/inp daec([i]) $sigma([daeci]*[sti])
enddo
set hcol 1
set plci 1
set fcol 1
ve/write ai,dai,ashs,dashs ampsat_counter[ncnt]_v[corrz].txt '4e15.6'
*
read x
* exec $PER/s#vpl ashe dashe zwls dzwls iatt=20
* exec $PER/s#vpl aec daec zwls dzwls iatt=24 o=s
n1=2
n11=2; n21=[n]-5
n12=2; n22=[n]-3
n13=3; n23=[n]-7
n14=2; n24=[n]-4
n15=2; n25=[n]-4
n16=2; n26=[n]-4
n17=2; n27=[n]-5
n18=2; n28=[n]-5
n19=2; n29=[n]-5 
n1=[n1[ncnt]]
n2=[n2[ncnt]]
ve/del aecr,daecr,zwlsr,dzwlsr
ve/copy  aec([n1]:[n2])  aecr
ve/copy daec([n1]:[n2]) daecr
ve/copy  zwls([n1]:[n2])  zwlsr
ve/copy dzwls([n1]:[n2]) dzwlsr
ve/cre wls(3) r 0.5 0.5 1000000
ve/cre daec([n]) r [n]*0.5
do i=1,5
  ve/fit zwlsr aecr daecr wlscorr.f s 3 p3 ! ! ! dp3
enddo
sigma p3=abs(p3)
ve/write p3,dp3 corrwls_counter[ncnt]_v[corrz].txt '2e15.6'
*
ve/cre p6(6) r 0 0 0 0 3 0
ve/copy p3(1:3) p6(1:3)
ve/cre dp6(6) r
n1=1
n2=14
if ([ncnt].eq.7) then
  ve/inp dashe(10) 10
endif
ve/del asher,dasher,zwlsr,dzwlsr
ve/copy  ashe([n1]:[n2])  asher
ve/copy dashe([n1]:[n2]) dasher
ve/copy  zwls([n1]:[n2])  zwlsr
ve/copy dzwls([n1]:[n2]) dzwlsr
ve/cre s6(6) r 0 0 0 0.1 0 0
ve/cre pmin(6) r 0 0 0 -20 -10 0
ve/cre pmax(6) r 100 100 10000 20 10 10
ve/fit zwlsr asher dasher wlssat.f sb 6 p6 s6 pmin pmax dp6
ve/cre s6(6) r 0 0 0 0 0.1 0
ve/fit zwlsr asher dasher wlssat.f sb 6 p6 s6 pmin pmax dp6
read x
ve/cre s6(6) r 0.1 0.1 0.1 0.1 0.1 0
do i=1,5
  ve/fit zwlsr asher dasher wlssat.f sb 6 p6 s6 pmin pmax dp6
enddo
* exec $PER/s#vpl ashe dashe zwls dzwls iatt=20
ve/fit zwlsr asher dasher wlssat.f sb 6 p6 s6 pmin pmax dp6
*ve/fit zwlsr asher dasher wlssat.f s 6 p6 ! ! ! dp6
a=$sigma(abs(p6(1)))
b=$sigma(abs(p6(2)))
c=$sigma(abs(p6(3)))
x=11
null -12 11 0 $sigma(([a])*exp(-[x]/([c]))+([b])*exp([x]/([c])))
fun/pl ([a])*exp(-x/([c]))+([b])*exp(x/([c])) -12 11 s
* exec $PER/s#vpl ashe dashe zwls dzwls iatt=20 o=s
ve/fit zwlsr asher dasher wlssat.f sb 6 p6 s6 pmin pmax dp6
ve/copy p6(1:3) p3(1:3)
ve/copy dp6(1:3) dp3(1:3)
sigma p3 = abs(p3)
ve/write p3,dp3 corrwls_counter[ncnt]_v[corrz].txt '2e15.6'
*
read x
ve/copy p6(4:6) p3(1:3)
fun/pl tadcutp.f 0 45
set pmci 1
set hcol 1
* exec $PER/s#vpl asx dasx aix daix o=s sz=0.01
read x
*
n1=[n1[ncnt]]
n2=[n2[ncnt]]
ve/del ashsr,dashsr,zwlsr,dzwlsr
ve/copy  ashs([n1]:[n2])  ashsr
ve/copy dashs([n1]:[n2]) dashsr
ve/copy  zwls([n1]:[n2])  zwlsr
ve/copy dzwls([n1]:[n2]) dzwlsr
ve/copy p6(1:3) p3(1:3)
ve/cre wls(3) r 0.5 0.5 1000000
do i=1,5
  ve/fit zwlsr ashsr dashsr wlscorr.f ! 3 p3 ! ! ! dp3
enddo
ve/cre p6(6) r 0 0 0 0 0 0.05
ve/inp p6(1) $sigma(ab([ncnt]))
ve/inp p6(2) $sigma(af([ncnt]))
ve/inp p6(3) $sigma(tau([ncnt]))
*sigma p3=abs(p3)
*ve/copy p3(1:3) p6(1:3)
mess '>>>>'
read x
n1=1
n2=14
ve/inp dashs(1) 10
ve/del ashsr,dashsr,zwlsr,dzwlsr
ve/copy  ashs([n1]:[n2])  ashsr
ve/copy dashs([n1]:[n2]) dashsr
ve/copy  zwls([n1]:[n2])  zwlsr
ve/copy dzwls([n1]:[n2]) dzwlsr
ve/cre s6(6) r 0 0 0 0.1 0.1 0.
ve/cre pmin(6) r 0 0 0 -100 -10 0
ve/cre pmax(6) r 100 100 10000 100 10 10
ve/fit zwlsr ashsr dashsr wlssat.f sb 6 p6 s6 pmin pmax dp6
ve/cre s6(6) r 0.1 0.1 0 0.1 0.1 0. 
ve/fit zwlsr ashsr dashsr wlssat.f sb 6 p6 s6 pmin pmax dp6
ve/cre s6(6) r 0.1 0.1 0.1 0.1 0.1 0.
do i=1,10
  ve/fit zwlsr ashsr dashsr wlssat.f sb 6 p6 s6 pmin pmax dp6
enddo
a=ab([ncnt])
b=af([ncnt])
c=tau([ncnt])
set hcol 2
r=12.3
fun/pl ([a])*exp(-x/([c]))+([b])*exp(x/([c])) -12 11
set hcol 1
set plci 1
* exec $PER/s#vpl ashs dashs zwls dzwls iatt=20 o=s
ve/fit zwlsr ashsr dashsr wlssat.f sb 6 p6 s6 pmin pmax dp6
a=p6(1)
b=p6(2)
c=p6(3)
set hcol 1
fun/pl ([a])*exp(-x/([c]))+([b])*exp(x/([c])) -12 11 s
set hcol 1
*
read x
ve/cre aif(14) r
do i=1,14
  zi=zwls([i])
  sti=$sigma(12/sqrt(12**2+([zi])**2))
  fi=$sigma(([a])*exp(-([zi])/([c])))
  fi=$sigma([fi]+([b])*exp([zi]/([c])))
  ve/inp aif([i]) $sigma([fi]/[sti]*(1+p6(6)))
enddo
null 0 $sigma(1.2*vmax(ai)) 0 $sigma(1.2*vmax(as))
set pmci 2
* exec $PER/s#vpl as das aif dai sz=0.05 o=s
set pmci 1
* exec $PER/s#vpl as das ai dai sz=0.05 o=s
* exec $PER/s#vpl asx dasx aix daix o=s sz=0.01
ve/copy p6(4:6) p3(1:3)
fun/pl tadcutp.f 0 90 s
set pmci 1
set hcol 1
*
endif
return


macro zslixcorrnx map=1
gl/cre corrz 1
do nc=1,9
  if (([map].eq.4).and.([nc].eq.9)) then
    ib=4
  else
    ib=1
  endif
  exec mapcal#zslixcorrn [map] [nc] [ib] 1
enddo
return


macro zslixcorrn map=1 ncnt=1 ib=1 auto=1 ver=0
*
if ([map].eq.1) then
  gl/cre mname mapcal1
  nh1=11
  nh2=11
  npt=0
  nb=14
endif
if ([map].eq.2) then
  gl/cre mname mapcal2
  nh1=12
  nh2=12
  npt=1
  nb=14
endif
if ([map].eq.3) then
  gl/cre mname mapcal3
  nh1=13
  nh2=13
  npt=2
  nb=10
endif
if ([map].eq.4) then
  gl/cre mname mapcal4
  nh1=14
  nh2=14
  npt=3
  nb=14
endif
nh1s=[nh1]+100
nh2s=[nh2]+100
*
ve/del ptl
ve/read ptl phicuts_test.txt
ks=0
ind=$sigma(128*[npt]+14*([ncnt]-1)+14)
dxl=$sigma(ptl([ind]))
do i=1,5
  ind=$sigma(128*[npt]+14*([ncnt]-1)+8+[i])
  f[i]s=$sigma(ptl([ind])-[dxl]*[ks])
enddo
f1=[f1s]
f2=[f2s]
f3=[f3s]
f4=[f5s]
*
mnf=320
mfmax=[ncnt]*40+20
mfmin=[ncnt]*40-60
i1=$sigma(int([mnf]*([f1]-[mfmin])/([mfmax]-[mfmin])+0.5))
i2=$sigma(int([mnf]*([f2]-[mfmin])/([mfmax]-[mfmin])+0.5))
i3=$sigma(int([mnf]*([f3]-[mfmin])/([mfmax]-[mfmin])+0.5))
i4=$sigma(int([mnf]*([f4]-[mfmin])/([mfmax]-[mfmin])+0.5))
i1=$sigma(int(([i1]-1)/[ib])*[ib]+1)
i2=$sigma(int(([i2]-1)/[ib])*[ib]+1)
i3=$sigma(int(([i3]-1)/[ib])*[ib]+1)
i4=$sigma(int(([i4]-1)/[ib])*[ib]+1)
ve/cre cf0([mnf]) r
do i=[i1],$sigma([i2]-1)
  ve/inp cf0([i]) 1
enddo
do i=[i3],$sigma([i4]-1)
  ve/inp cf0([i]) 1
enddo
*ve/cre cf0([mnf]) r
*do i=[i1],[i4]
*  ve/inp cf0([i]) 1
*enddo
*
do l=1,3
  ve/del pars0,dpars0
  flname=[mname]/[mname]_amplitude_vs_zr_counter[l]10[ncnt]_0.log.pars
  if ($fexist([flname]).eq.0) then
    flname=[mname]/[mname]_amplitude_vs_zr_counter[l]10[ncnt]_0.log
    shell /work/snd2000/users/konctbel/SNDSimACC/readlifstat [flname]
  endif
  ve/read pars0,dpars0 [mname]/[mname]_amplitude_vs_zr_counter[l]10[ncnt]_0.log.pars
  nr=[nb]+1
  ve/cre zi([nr]) r
  i1=2; i2=[nr]+1
  ve/copy pars0([i1]:[i2]) zi(1:[nr])
  exec mapcal#zicr zi zic[l]
enddo
sigma zic = (zic1+zic2+zic3)/3
n=$vlen(zic)
n=[n]-1
*
ve/cre  corrslix[ncnt]([n]) r
ve/cre dcorrslix[ncnt]([n]) r
ve/cre cslx[ncnt]p0([n]) r
ve/cre cslx[ncnt]p1([n]) r
ve/cre cslx[ncnt]p2([n]) r
*
fname=[mname]/[mname]_exp_v0_amplitude_vs_phi_counter[ncnt]_slx150_fit_100_et0_chit.par
ve/del pars,dpars
ve/read pars,dpars [fname] '2f15.10'
*
ve/cre p4(4) r 1 0 $sigma(pars(13)) $sigma(pars(17))
ve/cre dp4(4) r
ve/cre s4(4) r 0.01 0 0 0
ve/cre p6(6) r 1 0 0 0 $sigma(pars(13)) $sigma(pars(17))
ve/cre s6(6) r 0.01 0.0001 0.00001 0 0 0
ve/cre dp6(6) r
*
do i=1,[n]
*
  j=[i]+1
  z1=zic([i])
  z2=zic([j])
  exec mapcal#histprep z [ncnt] [nh1s] [nh2s] [z1] [z2] [ver]
  hi/copy 200 300
  exec mapcal#histprep z [ncnt] [nh1] [nh2] [z1] [z2] [ver]
  set ksiz 0.05
  set pmci 4
  hi/pl 200
  set pmci 2
  hi/pl 300 s
  line $sigma(([f2s]+[f3s])/2) 0 $sigma(([f2s]+[f3s])/2) 100
  if ([auto].eq.0) then
    read x
  endif
*
  ve/del cf
  ve/copy cf0 cf
  if ([ib].ne.1) then
    exec ../SepPar/sp#brn 200 [ib] 1000
    exec ../SepPar/sp#brn 300 [ib] 1000
    set pmci 4
    hi/pl 1200
    set pmci 2
    hi/pl 1300 s
    line $sigma(([f2s]+[f3s])/2) 0 $sigma(([f2s]+[f3s])/2) 100
    if ([auto].eq.0) then
      read x
    endif
    nx=$hinfo(200,'xbins')
    xmin=$hinfo(200,'xmin')
    xmax=$hinfo(200,'xmax')
    1d 500 ! [nx] [xmin] [xmax]
    hi/put/cont 500 cf
    exec ../SepPar/sp#brn 500 [ib] 1000
    nx=$hinfo(1500,'xbins')
    ve/cre cf([nx]) r
    hi/get/cont 1500 cf
  else
    hi/copy 200 1200
    hi/copy 300 1300
  endif
*
  nx=$hinfo(1200,'xbins')
  xmin=$hinfo(1200,'xmin')
  xmax=$hinfo(1200,'xmax')
  ve/cre xr([nx]) r 
  l=$sigma([xmin]+([xmax]-[xmin])/[nx]/2)
  r=$sigma([xmax]-([xmax]-[xmin])/[nx]/2)
  sigma xr=array([nx],[l]#[r])
  1d 1400 ! [nx] [xmin] [xmax]
  ve/del v1,v2,dv1,dv2
  ve/cre v1([nx]) r
  ve/cre v2([nx]) r
  ve/cre dv1([nx]) r
  ve/cre dv2([nx]) r
  hi/get/cont 1200 v1
  hi/get/cont 1300 v2
  hi/get/err 1200 dv1
  hi/get/err 1300 dv2
  sigma vr = cf*v1/v2
  sigma dvr = cf*vr*sqrt((dv1/v1)**2+(dv2/v2)**2)
  hi/put/cont 1400 vr
  hi/put/err 1400 dvr
*
  sigma xr=order(xr,-vr)
  sigma dvr=order(dvr,-vr)
  sigma vr=order(vr,-vr)
  exec mapcal#vecut vr
  nx=$vlen(vr)
  exec mapcal#vecut dvr [nx]
  exec mapcal#vecut xr [nx]
  sigma vr=order(vr,xr)
  sigma dvr=order(dvr,xr)
  sigma xr=order(xr,xr)
  ve/cre dxr([nx]) r
*  
  hi/pl 1400
  * exec $PER/s#vpl vr dvr xr dxr sz=0.05 o=s
  hi/fit 1400 p0+g sb 4 p4 s4 ! ! dp4
  ve/inp  corrslix[ncnt]([i]) $sigma(p4(1))
  ve/inp dcorrslix[ncnt]([i]) $sigma(dp4(1))
  do j=1,5
    hi/fit 1400 p2+g sb 6 p6 s6 ! ! dp6
  enddo
  set plci 4
  set lwid 5
  ve/cre p3(3) r [ib] 0 0
  do j=1,5
    ve/fit xr vr dvr p2 s 3 p3
  enddo
  set lwid
  set plci 1
*  ve/inp cslx[ncnt]p0([i]) $sigma(p6(1)/[ib])
*  ve/inp cslx[ncnt]p1([i]) $sigma(p6(2)/[ib])
*  ve/inp cslx[ncnt]p2([i]) $sigma(p6(3)/[ib])
  ve/inp cslx[ncnt]p0([i]) $sigma(p3(1)/[ib])
  ve/inp cslx[ncnt]p1([i]) $sigma(p3(2)/[ib])
  ve/inp cslx[ncnt]p2([i]) $sigma(p3(3)/[ib])
  if ([auto].eq.0) then
    read x
  endif
enddo
gl/imp corrz
ve/write corrslix[ncnt],dcorrslix[ncnt],cslx[ncnt]p0,cslx[ncnt]p1,cslx[ncnt]p2 [mname]/[mname]_corrslix_counter[ncnt]_v[corrz].txt '5e15.6'
return


macro expsimfn map=1 ncnt=1 ver=0
*
if ([map].eq.1) then
  gl/cre mname mapcal1
  nh1=11
  nh2=11
  npt=0
  nb=14
endif
if ([map].eq.2) then
  gl/cre mname mapcal2
  nh1=12
  nh2=12
  npt=1
  nb=14
endif
if ([map].eq.3) then
  gl/cre mname mapcal3
  nh1=13
  nh2=13
  npt=2
  nb=10
endif
if ([map].eq.4) then
  gl/cre mname mapcal4
  nh1=14
  nh2=14
  npt=3
  nb=14
endif
nh1s=[nh1]+100
nh2s=[nh2]+100
*
  z1=-10
  z2=9
  exec mapcal#histprep z [ncnt] [nh1s] [nh2s] [z1] [z2] [ver]
  hi/copy 200 300
  exec mapcal#histprep z [ncnt] [nh1] [nh2] [z1] [z2] [ver]
  set ksiz 0.05
  set pmci 4
  hi/pl 200
  set pmci 2
  hi/pl 300 s
*  
return


macro expsimzn map=1 ncnt=1 layer=1 ver=0
*
if ([map].eq.1) then
  gl/cre mname mapcal1
  nh1=11
  nh2=11
  npt=0
  nb=14
endif
if ([map].eq.2) then
  gl/cre mname mapcal2
  nh1=12
  nh2=12
  npt=1
  nb=14
endif
if ([map].eq.3) then
  gl/cre mname mapcal3
  nh1=13
  nh2=13
  npt=2
  nb=10
endif
if ([map].eq.4) then
  gl/cre mname mapcal4
  nh1=14
  nh2=14
  npt=3
  nb=14
endif
nh1s=[nh1]+100
nh2s=[nh2]+100
*
ve/del ptl
ve/read ptl phicuts_test.txt
ks=0
ind=$sigma(128*[npt]+14*([ncnt]-1)+14)
dxl=$sigma(ptl([ind]))
do i=1,5
  ind=$sigma(128*[npt]+14*([ncnt]-1)+8+[i])
  f[i]s=$sigma(ptl([ind])-[dxl]*[ks])
enddo

if ([layer].eq.1) then
  f1=[f1s]
  f2=[f2s]
endif
if ([layer].eq.2) then
  f1=[f3s]
  f2=[f4s]
endif
if ([layer].eq.3) then
  f1=[f4s]
  f2=[f5s]
endif
*
  exec mapcal#histprep f [ncnt] [nh1s] [nh2s] [f1] [f2] [ver]
  hi/copy 200 300
  exec mapcal#histprep f [ncnt] [nh1] [nh2] [f1] [f2] [ver]
  set ksiz 0.05
  set pmci 4
  hi/pl 200
  set pmci 2
  hi/pl 300 s
*  
return



macro ampsatpl
null 0 100 0 50
line 0 0 50 50
amax=35
do i=1,9
  ve/del ai,dai,as,das
  ve/read ai,dai,as,das ampsatsim_counter[i].txt '4f15.6'
  sigma x = order(ai,-as)
  sigma y = order(as,-as)
  mess I=[i] Ai=$sigma(x(1)) As=$sigma(y(1))
  da=$sigma([amax]-y(1))
*  sigma as=as+[da]
*  sigma ai=ai+[da]
  * exec $PER/s#vpl as das ai dai sz=0.01 o=s
enddo
return



macro ampsatsimx
do i=1,9
  exec mapcal#ampsatsim [i]
enddo
return

macro ampsatsim ncnt=1
name=spect0
ve/cre ai(101) r
ve/cre as(101) r
ve/cre dai(101) r
ve/cre das(101) r
1d 100 ! 2200 -20 200
do i=0,100
  fname=/work/users/konctbel/snd2k/R005-999/[name]_[i]pe.txt
  if ($fexist([fname]).eq.1) then
    nt/cre 1 ! 9 ! ! a1 a2 a3 a4 a5 a6 a7 a8 a9
    nt/read 1 [fname]
    nt/pl 1.a[ncnt] idh=100
    j=[i]+1
    mean=$hinfo(100,'mean')
    rms=$hinfo(100,'rms')
    nevt=$hinfo(100,'events')
    ve/inp ai([j]) [i]
    ve/inp dai([j]) 0
    ve/inp as([j]) [mean]
    ve/inp das([j]) $sigma([rms]/sqrt([nevt]))
  endif
enddo
exec mapcal#vecut ai
exec mapcal#vecut dai $vlen(ai)
exec mapcal#vecut as $vlen(ai)
exec mapcal#vecut das $vlen(ai)
amax = $sigma(1.2*vmax(ai))
null 0 [amax] 0 [amax]
* exec $PER/s#vpl as das ai dai o=s sz=0.05
set plci 4
line 0 0 [amax] [amax]
ve/cre p3(3) r 25 0.3 5
ve/cre p3min(3) r 0 0.1 3 
ve/cre p3max(3) r 100 1 100
ve/cre p3stp(3) r 0.1 0.1 0.1
*ve/fit ai as das ampsat.f sb 3 p3 p3stp p3min p3max
set pmci 1
set hcol 1
*fun/pl ampsatp.f 0 [amax] s
atitle 'A?I!, pe' 'A?O!, pe'
exec save ampsatsim0_counter[ncnt].eps f
sigma asai = as - ai
* exec $PER/s#vpl asai das ai dai sz=0.05 ll=-1
ve/cre p3(3) r 0.1 -0.1 30
ve/fit ai asai das expp0.f s 3 p3
set pmci 1
set hcol 1
*fun/pl ampsatp.f 0 [amax] s
atitle 'A?I!, pe' 'A?O!, pe'
exec save ampsatsim0d_counter[ncnt].eps f
sigma x = order(ai,-as)
sigma y = order(as,-as)
mess Ai=$sigma(x(1)) As=$sigma(y(1))
ve/write ai,dai,as,das ampsatsim0_counter[ncnt].txt '4f15.6'
return

macro ampsatfit ncnt=1
ve/del ai,dai,as,das
ve/read ai,dai,as,das ampsatsim_counter[ncnt].txt '4f15.6'
n1=1
n2=42
ve/del air,dair,asr,dasr
ve/copy ai([n1]:[n2]) air
ve/copy dai([n1]:[n2]) dair
ve/copy as([n1]:[n2]) asr
ve/copy das([n1]:[n2]) dasr
* exec $PER/s#vpl as das ai dai sz=0.05
ve/cre p3(3) r 40 20 -1.5
ve/fit air asr dasr tadcut.f s 3 p3
return

macro ampsat
corrz=2
ve/cre  aix(1) r 0
ve/cre daix(1) r 0
ve/cre  ashsx(1) r 0
ve/cre dashsx(1) r 0
do i=1,9
  ve/del ai,dai,ashs,dashs
  ve/read ai,dai,ashs,dashs ampsat_counter[i]_v[corrz].txt '4e15.6'
  exec vappend  aix  ai
  exec vappend daix dai
  exec vappend  ashsx  ashs
  exec vappend dashsx dashs
enddo
sigma  ashsx = order(ashsx,aix)
sigma dashsx = order(dashsx,aix)
sigma   daix = order(daix,aix)
sigma    aix = order(aix,aix)
amax = $sigma(1.2*vmax(aix))
null 0 [amax] 0 [amax]
* exec $PER/s#vpl ashsx dashsx aix daix o=s
set plci 4
line 0 0 [amax] [amax]
ve/cre p3(3) r 25 15 5
ve/cre p3min(3) r 0 10 3
ve/cre p3max(3) r 100 100 100
ve/cre p3stp(3) r 0.1 0.1 0.1
ve/fit aix ashsx dashsx ampsat.f sb 3 p3 p3stp p3min p3max
set hcol 2
fun/pl ampsatp.f 0 [amax] s
*
amax = $sigma(1.2*vmax(aix))
null 0 [amax] 0 [amax]
set plci 4
line 0 0 [amax] [amax]
ve/cre p3(3) r 25 15 5
ve/cre p3min(3) r 0 10 3
ve/cre p3max(3) r 100 100 100
ve/cre p3stp(3) r 0.1 0.1 0.1
do i=1,9
  ve/del ai,dai,ashs,dashs
  ve/read ai,dai,ashs,dashs ampsat_counter[i]_v[corrz].txt '4e15.6'
  * exec $PER/s#vpl ashsx dashsx aix daix o=s sz=0.05
  ve/fit ai ashs dashs ampsat.f sb 3 p3 p3stp p3min p3max
  set hcol 2
  fun/pl ampsatp.f 0 [amax] s
enddo
set hcol 1
set plci 1
set fcol 1
return

macro expsimf ncnt=1 z1=-8.0 z2=8.0 slix=0
if ([slix].ne.0) then
gl/imp ncnto
if ([ncnt].ne.[ncnto]) then
  gl/cre mname mhad2011-4
  exec mapcal#fitpl 110[ncnt] 1
  ve/copy zic zic1
  exec mapcal#fitpl 210[ncnt] 1
  ve/copy zic zic2
  exec mapcal#fitpl 310[ncnt] 1
  ve/copy zic zic3
  sigma zic = (zic1+zic2+zic3)/3
  gl/cre ncnto [ncnt]
endif
  i=[slix]
  z1=zic([i])
  i=[i]+1
  z2=zic([i])
endif
mnz=300
mzmin=-15
mzmax=15
  i1=$sigma([mnz]*([z1]-[mzmin])/([mzmax]-[mzmin]))
  j=[i]+1
  i2=$sigma([mnz]*([z2]-[mzmin])/([mzmax]-[mzmin]))
  ic=$sigma(int(([i1]+[i2])/2+0.5))
  di=$sigma(int(([i2]-[i1])/2+0.5))
  mess [ic] [di]
  gl/cre ic [ic]
  gl/cre di [di]

*  hi/file 30 mhad2011_profile_sin_new.his.old
*  hi/file 30 mhad2011_profile_sin_new_et.his
*   hi/file 30 mhad2011_profile_sin_new.his
*  hi/file 30 mhad2011_profile_sin_new_et_chit.his
*  hi/file 30 mhad2011_profile_sin_new_et_chit_lcorr.his
*  hi/file 30 mhad2011-4_profile_sin_new_et_chit.his
*  hi/file 30 mhad2011-4_profile_sin_new_et_chit_acorr_675.his
 hi/file 30 mhad2011-4_profilex_ee_sim_shifter7_chit.his

  idh=200
  if ($hexist([idh]).eq.1) then
    hi/del [idh]
  endif
  i1=$sigma(int([i1]+0.5))
  i2=$sigma(int([i2]+0.5))
  do j=[i1],[i2]

    idhj=[j]*100+10*3+[ncnt]
    
*    if ($hexist([idh]).eq.0) then
*      hi/copy [idhj] [idh]
*    else
*      hi/op/add [idhj] [idh] [idh]
*    endif
*  enddo
    if ($hexist([idh]).eq.0) then
      hi/pl //lun30/[idhj]
      nx=$hinfo([idhj],'xbins')
      xmin=$hinfo([idhj],'xmin')
      xmax=$hinfo([idhj],'xmax')
      if ([nx].eq.0) then
        hi/pl [idhj]
        read x
      endif
      1d [idh] ! [nx] [xmin] [xmax]
      ve/del xis,xi2s,si2s,xi,si,xi2,si2,sis,ni
      ve/cre ni([nx]) r
      ve/cre xis([nx]) r
      ve/cre xi2s([nx]) r
      ve/cre si2s([nx]) r
      ve/cre xi([nx]) r
      ve/cre xi2([nx]) r
      ve/cre si([nx]) r
      ve/cre si2([nx]) r
      ve/cre sis([nx]) r
    endif
    hi/get/cont [idhj] xi
    hi/get/err  [idhj] si
    sigma ni = ni + xi/xi
    sigma xi2 = xi**2
    sigma si2 = si**2
    sigma xis = xis + xi
    sigma xi2s = xi2s + xi2
    sigma si2s = si2s + si2
    hi/del [idhj]
  enddo
  sigma xis = xis/ni
  sigma sis = sqrt(((si2s+xi2s)/ni-xis**2)/ni)
  hi/put/cont [idh] xis
  hi/put/erro [idh] sis
  
  close 30

*  hi/file 30 profilex_ee_sim_corr_shift.his
*  hi/file 30 profilex_ee_sim_corr_sim.his
*  hi/file 30 profilex_ee_sim_corr_sim_shift.his
*  hi/file 30 mhad2011_profile_sin_new.his
*  hi/file 30 profilex_ee_sim_corr_sim_shift_new.his
*  hi/file 30 profilex_ee_sim_corr_sim_shifttest_new.his
*  hi/file 30 profilex_ee_sim_corr_sim_shiftreal_new.his
*  hi/file 30 profilex_ee_sim_shifter0_chit.his
*  hi/file 30 profilex_ee_sim_shifter1_chit.his
*  hi/file 30 profilex_ee_sim_shifter2_chit.his
*  hi/file 30 mhad2011-4_profile_sin_new_et_chit.his
*  hi/file 30 profilex_ee_sim_shifter_chit.his
* hi/file 30 mhad2011-4_profilex_ee_sim_shifter101_chit.his
* hi/file 30 mhad2011-4_profilex_ee_sim_shifter102_chit.his
* hi/file 30 mhad2011-4_profilex_ee_sim_shifter_chit.his
* hi/file 30 mhad2011-4_profilex_ee_sim_shifter0_chit.his
* hi/file 30 mhad2011-4_profilex_ee_sim_shifter1_chit.his
* hi/file 30 mhad2011-4_profilex_ee_sim_shifter2_chit.his
* hi/file 30 mhad2011-4_profilex_ee_sim_shifter3_chit.his
* hi/file 30 mhad2011-4_profilex_ee_sim_shifter4_chit.his
 hi/file 30 mhad2011-4_profilex_ee_sim_shifter5_chit.his

  idh=300
  if ($hexist([idh]).eq.1) then
    hi/del [idh]
  endif
  i1=$sigma(int([i1]+0.5))
  i2=$sigma(int([i2]+0.5))
  do j=[i1],[i2]

    idhj=[j]*100+10*3+[ncnt]
    
*    if ($hexist([idh]).eq.0) then
*      hi/copy [idhj] [idh]
*    else
*      hi/op/add [idhj] [idh] [idh]
*    endif
*  enddo
    if ($hexist([idh]).eq.0) then
      hi/pl //lun30/[idhj]
      nx=$hinfo([idhj],'xbins')
      xmin=$hinfo([idhj],'xmin')
      xmax=$hinfo([idhj],'xmax')
      if ([nx].eq.0) then
        hi/pl [idhj]
        read x
      endif
      1d [idh] ! [nx] [xmin] [xmax]
      ve/del xis,xi2s,si2s,xi,si,xi2,si2,sis,ni
      ve/cre ni([nx]) r
      ve/cre xis([nx]) r
      ve/cre xi2s([nx]) r
      ve/cre si2s([nx]) r
      ve/cre xi([nx]) r
      ve/cre xi2([nx]) r
      ve/cre si([nx]) r
      ve/cre si2([nx]) r
      ve/cre sis([nx]) r
    endif
    hi/get/cont [idhj] xi
    hi/get/err  [idhj] si
    sigma ni = ni + xi/xi
    sigma xi2 = xi**2
    sigma si2 = si**2
    sigma xis = xis + xi
    sigma xi2s = xi2s + xi2
    sigma si2s = si2s + si2
    hi/del [idhj]
  enddo
  sigma xis = xis/ni
  sigma sis = sqrt(((si2s+xi2s)/ni-xis**2)/ni)
  hi/put/cont [idh] xis
  hi/put/erro [idh] sis
  
  close 30

  hi/file 30 profilex_ee_sim_shifter0_chit.his

  idh=400
  if ($hexist([idh]).eq.1) then
    hi/del [idh]
  endif
  i1=$sigma(int([i1]+0.5))
  i2=$sigma(int([i2]+0.5))
  do j=[i1],[i2]

    idhj=[j]*100+10*3+[ncnt]
    
    if ($hexist([idh]).eq.0) then
      hi/copy [idhj] [idh]
    else
      hi/op/add [idhj] [idh] [idh]
    endif
  enddo
  
  close 30

set ksiz 0.05
set pmci 2
hi/pl 200
set pmci 4
hi/pl 300 s

ve/del ic,if1,if2,if3,if4,if5
ve/read ic,if1,if2,if3,if4,if5 phi_cuts.txt
do i=1,5
  line $sigma(if[i]([ncnt])) 0 $sigma(if[i]([ncnt])) 100
enddo

return

macro expsimz map=1 ncnt=1 f1=5.0 f2=16.0 mode=0
*
if ([map].eq.1) then
  gl/cre mname mapcal1
  nh1=11
  nh2=11
  npt=0
endif
if ([map].eq.2) then
  gl/cre mname mapcal2
  nh1=12
  nh2=12
  npt=1
endif
if ([map].eq.3) then
  gl/cre mname mapcal3
  nh1=13
  nh2=13
  npt=2
endif
if ([map].eq.4) then
  gl/cre mname mapcal4
  nh1=14
  nh2=14
  npt=3
endif
*
ve/del ptl
ve/read ptl phicuts_test.txt
ks=0
ind=$sigma(128*[npt]+14*([ncnt]-1)+14)
dxl=$sigma(ptl([ind]))
do i=1,5
  ind=$sigma(128*[npt]+14*([ncnt]-1)+8+[i])
  f[i]=$sigma(ptl([ind])-[dxl]*[ks])
enddo
*
ve/del ic,if1,if2,if3,if4,if5
ve/read ic,if1,if2,if3,if4,if5 phi_cuts.txt
if ([mode].eq.1) then
  f1=if1([ncnt])
  f2=if2([ncnt])
endif
if ([mode].eq.2) then
  f1=if3([ncnt])
  f2=if4([ncnt])
endif
if ([mode].eq.3) then
  f1=if4([ncnt])
  f2=if5([ncnt])
endif
mnz=300
mzmin=-15
mzmax=15
mnf=320
mfmax=[ncnt]*40+20
mfmin=[ncnt]*40-60
  i1=$sigma([mnf]*([f1]-[mfmin])/([mfmax]-[mfmin]))
  j=[i]+1
  i2=$sigma([mnf]*([f2]-[mfmin])/([mfmax]-[mfmin]))
  ic=$sigma(int(([i1]+[i2])/2+0.5))
  di=$sigma(int(([i2]-[i1])/2+0.5))
  mess [ic] [di]

*  hi/file 30 mhad2011_profile_sin_new.his
*  hi/file 30 mhad2011_profile_sin_new_et.his
*  hi/file 30 mhad2011_profile_sin_new_et_chit_lcorr.his
*  hi/file 30 mhad2011_profile_sin_new_et_chit.his
*  hi/file 30 mhad2011-4_profile_sin_new_et_chit_lcorr.his
*  hi/file 30 mhad2011-4_profile_sin_new_et_chit.his
  hi/file 30 mhad2011-4_profile_sin_new_et_chit_acorr_675.his

  idh=210
  if ($hexist([idh]).eq.1) then
    hi/del [idh]
  endif
  i1=$sigma(int([i1]+0.5))
  i2=$sigma(int([i2]+0.5))
  do j=[i1],[i2]

    idhj=[j]*100+10*0+[ncnt]
    
    if ($hexist([idh]).eq.0) then
      hi/copy [idhj] [idh]
    else
      hi/op/add [idhj] [idh] [idh]
    endif
  enddo
  
  close 30

*  hi/file 30 profilex_ee_sim_corr_sim.his
*  hi/file 30 profilex_ee_sim_corr_sim_shift.his
*  hi/file 30 profilex_ee_sim_corr_sim_shift_new.his
*  hi/file 30 profilex_ee_sim_corr_sim_shifttest_new.his
*  hi/file 30 profilex_ee_sim_corr_sim_shiftreal_new.his
*  hi/file 30 profilex_ee_sim_shifter0_chit.his
*  hi/file 30 profilex_ee_sim_shifter1_chit.his
*  hi/file 30 profilex_ee_sim_shifter_chit.his
*  hi/file 30 mhad2011-4_profile_sin_new_et_chit_lcorr.his
*  hi/file 30 mhad2011-4_profilex_ee_sim_shifter0_chit.his
*  hi/file 30 mhad2011-4_profilex_ee_sim_shifter1_chit.his
*  hi/file 30 mhad2011-4_profilex_ee_sim_shifter2_chit.his
*  hi/file 30 mhad2011-4_profilex_ee_sim_shifter3_chit.his
  hi/file 30 mhad2011-4_profilex_ee_sim_shifter4_chit.his

  idh=310
  if ($hexist([idh]).eq.1) then
    hi/del [idh]
  endif
  i1=$sigma(int([i1]+0.5))
  i2=$sigma(int([i2]+0.5))
  do j=[i1],[i2]

    idhj=[j]*100+10*0+[ncnt]
    
    if ($hexist([idh]).eq.0) then
      hi/copy [idhj] [idh]
    else
      hi/op/add [idhj] [idh] [idh]
    endif
  enddo
  
  close 30

set ksiz 0.05
set pmci 2
hi/pl 210
set pmci 4
hi/pl 310 s

read x

hi/op/div 310 210 410 ! ! E
ve/cre p0(1) r 1
ve/cre dp0(1) r 
*null -15 15 0 2
*hi/fit 410(-8.:8.) p0 s 1 p0 ! ! ! dp0

return


macro expsimc idh=100001 ib=1
gl/imp mname
mess [mname]
if ([mname].eq.'mhad2011-4') then
  hi/file 30 mhad2011-4_profile_sin_new_et_chit.his
  mess mhad2011-4_profile_sin_new_et_chit.his
endif
if ([mname].eq.'mhad2011-4_v0') then
  hi/file 30 mhad2011-4_profile_sin_new_et_chit.his
  mess mhad2011-4_profile_sin_new_et_chit.his
endif
if ([mname].eq.'mhad2011-4_v1') then
  hi/file 30 mhad2011-4_profile_sin_new_et_chit.his
  mess mhad2011-4_profile_sin_new_et_chit.his
endif
if ([mname].eq.'mhad2011-4_v2') then
  hi/file 30 mhad2011-4_profile_sin_new_et_chit.his
  mess mhad2011-4_profile_sin_new_et_chit.his
endif
if ([mname].eq.'mhad2011-4_v3') then
  hi/file 30 mhad2011-4_profile_sin_new_et_chit_acorr_675.his
  mess mhad2011-4_profile_sin_new_et_chit_acorr_675.his
endif
if ([mname].eq.'mhad2011-4_v4') then
  hi/file 30 mhad2011-4_profile_sin_new_et_chit_acorr_925.his
  mess mhad2011-4_profile_sin_new_et_chit_acorr_925.his
endif
if ([mname].eq.'mhad2011-4_v5') then
  hi/file 30 mhad2011-4_profile_sin_new_et_chit_acorr_537.5.his
  mess mhad2011-4_profile_sin_new_et_chit_acorr_537.5.his
endif
if ([mname].eq.'mhad2011-4_v6') then
  hi/file 30 mhad2011-4_profilex_ee_sim_shifter6_chit.his
  mess mhad2011-4_profilex_ee_sim_shifter6_chit.his
endif
if ([mname].eq.'mhad2011-4_v7') then
  hi/file 30 mhad2011-4_profilex_ee_sim_shifter7_chit.his
  mess mhad2011-4_profilex_ee_sim_shifter7_chit.his
endif
if ([mname].eq.'mhad2011-4_v8') then
  hi/file 30 mhad2011-4_profilex_ee_sim_shifter8_chit.his
  mess mhad2011-4_profilex_ee_sim_shifter8_chit.his
endif
*  hi/file 30 profilex_ee_sim_corr_sim_shift_nonunif.his
*  hi/file 30 mhad2011_profile_sin_new.his
*  hi/file 30 mhad2011_profile_sin_new_et.his
*  hi/file 30 mhad2011_profile_sin_new_et_chit.his
*  hi/file 30 mhad2011-4_profile_sin_new_et_chit.his
*  hi/file 30 mhad2011-4_profile_sin_new_et_chit_acorr_675.his
*  hi/file 30 mhad2011-4_profile_sin_new_et_chit_acorr_925.his
*  hi/file 30 mhad2011-4_profile_sin_new_et_chit_acorr_537.5.his
*  hi/file 30 mhad2011-4_profile_sin_new_et_chit_lcorr.his
  hrin [idh]
  nx=$hinfo([idh],'events')
  mue=$hinfo([idh],'mean')
  rmse=$hinfo([idh],'rms')
  dmue=$sigma([rmse]/sqrt([nx]))
  hi/copy [idh] 100
  idh0=[idh]+10000
  hrin [idh0]
  n0=$hinfo([idh0],'events')
  if ([n0].eq.0) then
    n0=1
  endif
  effe=$sigma([n0]/[nx])
  deffe=$sigma(sqrt([effe]*(1-[effe])/[nx]))
  mess [nx] [n0]
  close 30
if ([mname].eq.'mhad2011-4') then
  hi/file 30 mhad2011-4_profilex_ee_sim_shifter0_chit.his
  mess mhad2011-4_profilex_ee_sim_shifter0_chit.his
endif
if ([mname].eq.'mhad2011-4_v0') then
  hi/file 30 mhad2011-4_profilex_ee_sim_shifter0_chit.his
  mess mhad2011-4_profilex_ee_sim_shifter0_chit.his
endif
if ([mname].eq.'mhad2011-4_v1') then
  hi/file 30 mhad2011-4_profilex_ee_sim_shifter1_chit.his
  mess mhad2011-4_profilex_ee_sim_shifter1_chit.his
endif
if ([mname].eq.'mhad2011-4_v2') then
  hi/file 30 mhad2011-4_profilex_ee_sim_shifter2_chit.his
  mess mhad2011-4_profilex_ee_sim_shifter2_chit.his
endif
if ([mname].eq.'mhad2011-4_v3') then
  hi/file 30 mhad2011-4_profilex_ee_sim_shifter3_chit.his
  mess mhad2011-4_profilex_ee_sim_shifter3_chit.his
endif
if ([mname].eq.'mhad2011-4_v4') then
  hi/file 30 mhad2011-4_profilex_ee_sim_shifter4_chit.his
  mess mhad2011-4_profilex_ee_sim_shifter4_chit.his
endif
if ([mname].eq.'mhad2011-4_v5') then
  hi/file 30 mhad2011-4_profilex_ee_sim_shifter5_chit.his
  mess mhad2011-4_profilex_ee_sim_shifter5_chit.his
endif
if ([mname].eq.'mhad2011-4_v6') then
  hi/file 30 mhad2011-4_profilex_ee_sim_shifter5_chit.his
  mess mhad2011-4_profilex_ee_sim_shifter5_chit.his
endif
if ([mname].eq.'mhad2011-4_v7') then
  hi/file 30 mhad2011-4_profilex_ee_sim_shifter5_chit.his
  mess mhad2011-4_profilex_ee_sim_shifter5_chit.his
endif
if ([mname].eq.'mhad2011-4_v8') then
  hi/file 30 mhad2011-4_profilex_ee_sim_shifter5_chit.his
  mess mhad2011-4_profilex_ee_sim_shifter5_chit.his
endif
*  hi/file 30 profilex_ee_sim_corr_sim_shift.his
*  hi/file 30 profilex_ee_sim_corr_sim_shift_nonunif.his
*  hi/file 30 profilex_ee_sim_corr_sim_shift_new.his
*  hi/file 30 profilex_ee_sim_corr_sim_shiftreal_new_spectcorr.his
*  hi/file 30 profilex_ee_sim_corr_sim_shiftreal_new.his
*  hi/file 30 profilex_ee_sim_shifter1_chit.his
*  hi/file 30 profilex_ee_sim_shifter2_chit.his
*  hi/file 30 profilex_ee_sim_shifter3a_chit.his
*  hi/file 30 profilex_ee_sim_shifter_chit.his
*  hi/file 30 mhad2011-4_profilex_ee_sim_shifter0_chit.his
*  hi/file 30 mhad2011-4_profilex_ee_sim_shifter1_chit.his
*  hi/file 30 mhad2011-4_profilex_ee_sim_shifter2_chit.his
*  hi/file 30 mhad2011-4_profilex_ee_sim_shifter3_chit.his
*  hi/file 30 mhad2011-4_profilex_ee_sim_shifter4_chit.his
*  hi/file 30 mhad2011-4_profilex_ee_sim_shifter5_chit.his
  hrin [idh]
  nx=$hinfo([idh],'events')
  mus=$hinfo([idh],'mean')
  rmss=$hinfo([idh],'mean')
  dmus=$sigma([rmss]/sqrt([nx]))
  hi/copy [idh] 200
  idh0=[idh]+10000
  hrin [idh0]
  n0=$hinfo([idh0],'events')
  if ([n0].eq.0) then
    n0=1
  endif
  effs=$sigma([n0]/[nx])
  deffs=$sigma(sqrt([effs]*(1-[effs])/[nx]))
  mess [nx] [n0]
  close 30
  mess exp=[effe]([deffe]) sim=[effs]([deffs])
  gl/cre  effe  [effe]
  gl/cre deffe [deffe]
  gl/cre  effs  [effs]
  gl/cre deffs [deffs]
  gl/cre  mue    [mue] 
  gl/cre  mus    [mus]
  gl/cre  dmue  [dmue] 
  gl/cre  dmus  [dmus]
  gl/cre  rmse  [rmse] 
  gl/cre  rmss  [rmss]
  n100=$hinfo(100,'events')
  n200=$hinfo(200,'events')
  hi/copy 200 300
  hi/op/res 300
  xbins=$hinfo(200,'xbins')
  ve/cre v200([xbins]) r
  ve/cre dv200([xbins]) r
  hi/get/cont 200 v200
  hi/get/err 200 dv200
  sigma  v200 =  v200*[n100]/[n200]
  sigma dv200 = dv200*[n100]/[n200]
  hi/put/cont 300 v200
  hi/put/err 300 dv200
*  hi/op/add 100 100 300 $sigma([n200]/[n100]) 0
*  exec hsigma @300 = @100/vsum(@100)*vsum(@200)
*  exec hsigma %300 = %100/vsum(@100)*vsum(@200)
  me = $hinfo(100,'mean')
  re = $hinfo(100,'rms')
  ms = $hinfo(200,'mean')
  rs = $hinfo(200,'rms')
  se = $sigma(sqrt([re]**2/[me]-1))
  ss = $sigma(sqrt([rs]**2/[ms]-1))
  corr = $sigma([se]/[ss])
  mess SigE = [se] SigS = [ss] Corr = [corr]
  xdiff=$call('mhdiff(100,200)')
  mess Ae = [me] As = [ms] As/Ae = $sigma([ms]/[me]) Diff = [xdiff]
  mess P0e = [effe] P0s = [effs]
  if ([ib].ne.1) then
    exec ../SepPar/sp#brn 100 [ib] 1000
    exec ../SepPar/sp#brn 300 [ib] 1000
    hi/del 1301
    exec hsigma @1301 = @1300
    exec hsigma u = 1.15*max(vmax(@1301),vmax(@1100))
    null -1 35 0 $sigma(u)
    hi/pl 1301 s
    hi/pl 1100 s
  else
    exec hsigma u = 1.2*max(@300+%300,@100+%100)
    null -1 35 0 $sigma(u)
    hi/pl 300 s
    hi/pl 100 s
  endif
  nc=$sigma([idh]-10*int([idh]/10))
  dir=$sigma(([idh]-100*int([idh]/100)-[nc])/10)
  np=$sigma(([idh]-1000*int([idh]/1000)-10*[dir]-[nc])/100)
  if ([dir].eq.6) then
    pos='phi'
  else
    pos='z'
  endif
  set chhe 0.2
  text=diff = [xdiff]
  exec $PER/s#tf 0.6 0.8 [text]
  k=3; p=$sigma(int([mue]*10**[k]+0.5)/10**[k])
  k=3; dp=$sigma(int([dmue]*10**[k]+0.5)/10**[k])
  text=A?exp! = [p] "A# [dp]
*"
  exec $PER/s#tf 0.6 0.7 [text]
  k=3; p=$sigma(int([mus]*10**[k]+0.5)/10**[k])
  k=3; dp=$sigma(int([dmus]*10**[k]+0.5)/10**[k])
  text=A?sim! = [p] "A# [dp]
*"
  exec $PER/s#tf 0.6 0.6 [text]
  dk=2
  k=5; p=$sigma(int([effe]*10**[k]+0.5)/10**([k]-[dk]))
  k=5; dp=$sigma(int([deffe]*10**[k]+0.5)/10**([k]-[dk]))
  text=P?0,exp! = ([p] "A# [dp])[\327]10^-[dk]!
*"
  exec $PER/s#tf 0.6 0.5 [text]
  k=5; p=$sigma(int([effs]*10**[k]+0.5)/10**([k]-[dk]))
  k=5; dp=$sigma(int([deffs]*10**[k]+0.5)/10**([k]-[dk]))
  text=P?0,sim! = ([p] "A# [dp])[\327]10^-[dk]!
*"
  exec $PER/s#tf 0.6 0.4 [text]
  atitle 'pusle height, pe' 'events'
  set chhe 0.309
  exec save [mname]_counter[nc]_amplitude_vs_$unquote([pos])_pos[np]_v1.eps f
*
*  hi/copy 300000 400
*  exec hsigma tot = @400
*  exec hsigma dtot = %400
return


macro expsimcn map=1 idh=100001 ib=1 ver=0
*
if ([ver].le.1) then
  suffix=exp.his
endif
if ([ver].eq.2) then
  suffix=exp_full.his
endif
if ([ver].eq.3) then
  suffix=exp_final.his
endif
if ([map].eq.1) then
  gl/cre mname mapcal1
  hi/file 30 mapcal1_profilex_[suffix]
  mess mapcal1_profilex_[suffix]
endif
if ([map].eq.2) then
  gl/cre mname mapcal2
  hi/file 30 mapcal2_profilex_[suffix]
  mess mapcal2_profilex_[suffix]
endif
if ([map].eq.3) then
  gl/cre mname mapcal3
  hi/file 30 mapcal3_profilex_[suffix]
  mess mapcal3_profilex_[suffix]
endif
if ([map].eq.4) then
  gl/cre mname mapcal4
  hi/file 30 mapcal4_profilex_[suffix]
  mess mapcal4_profilex_[suffix]
endif
*
  hrin [idh]
  nx=$hinfo([idh],'events')
  mue=$hinfo([idh],'mean')
  rmse=$hinfo([idh],'rms')
  dmue=$sigma([rmse]/sqrt([nx]))
  hi/copy [idh] 100
  idh0=[idh]+10000
  hrin [idh0]
  n0=$hinfo([idh0],'events')
  if ([n0].eq.0) then
    n0=1
  endif
  effe=$sigma([n0]/[nx])
  deffe=$sigma(sqrt([effe]*(1-[effe])/[nx]))
  mess [nx] [n0]
  close 30
*  
*vers=$sigma(min(1,[ver]))
vers=[ver]
if ([map].eq.1) then
  gl/cre mname mapcal1
  hi/file 30 mapcal1_profilex_sim_v[vers].his
  mess mapcal1_profilex_sim_v[vers].his
endif
if ([map].eq.2) then
  gl/cre mname mapcal2
  hi/file 30 mapcal2_profilex_sim_v[vers].his
  mess mapcal2_profilex_sim_v[vers].his
endif
if ([map].eq.3) then
  gl/cre mname mapcal3
  hi/file 30 mapcal3_profilex_sim_v[vers].his
  mess mapcal3_profilex_sim_v[vers].his
endif
if ([map].eq.4) then
  gl/cre mname mapcal4
  hi/file 30 mapcal4_profilex_sim_v[vers].his
  mess mapcal4_profilex_sim_v[vers].his
endif
*
  hrin [idh]
  nx=$hinfo([idh],'events')
  mus=$hinfo([idh],'mean')
  rmss=$hinfo([idh],'mean')
  dmus=$sigma([rmss]/sqrt([nx]))
  hi/copy [idh] 200
  idh0=[idh]+10000
  hrin [idh0]
  n0=$hinfo([idh0],'events')
  if ([n0].eq.0) then
    n0=1
  endif
  effs=$sigma([n0]/[nx])
  deffs=$sigma(sqrt([effs]*(1-[effs])/[nx]))
  mess [nx] [n0]
  close 30
  mess exp=[effe]([deffe]) sim=[effs]([deffs])
  gl/cre  effe  [effe]
  gl/cre deffe [deffe]
  gl/cre  effs  [effs]
  gl/cre deffs [deffs]
  gl/cre  mue    [mue] 
  gl/cre  mus    [mus]
  gl/cre  dmue  [dmue] 
  gl/cre  dmus  [dmus]
  gl/cre  rmse  [rmse] 
  gl/cre  rmss  [rmss]
  n100=$hinfo(100,'events')
  n200=$hinfo(200,'events')
  hi/copy 200 300
  hi/op/res 300
  xbins=$hinfo(200,'xbins')
  ve/cre v200([xbins]) r
  ve/cre dv200([xbins]) r
  hi/get/cont 200 v200
  hi/get/err 200 dv200
  sigma  v200 =  v200*[n100]/[n200]
  sigma dv200 = dv200*[n100]/[n200]
  hi/put/cont 300 v200
  hi/put/err 300 dv200
*  hi/op/add 100 100 300 $sigma([n200]/[n100]) 0
*  exec hsigma @300 = @100/vsum(@100)*vsum(@200)
*  exec hsigma %300 = %100/vsum(@100)*vsum(@200)
  me = $hinfo(100,'mean')
  re = $hinfo(100,'rms')
  ms = $hinfo(200,'mean')
  rs = $hinfo(200,'rms')
  se = $sigma(sqrt([re]**2/[me]-1))
  ss = $sigma(sqrt([rs]**2/[ms]-1))
  corr = $sigma([se]/[ss])
  mess SigE = [se] SigS = [ss] Corr = [corr]
  xdiff=$call('mhdiff(100,200)')
  mess Ae = [me] As = [ms] As/Ae = $sigma([ms]/[me]) Diff = [xdiff]
  mess P0e = [effe] P0s = [effs] 
  if ([ib].ne.1) then
    exec ../SepPar/sp#brn 100 [ib] 1000
    exec ../SepPar/sp#brn 300 [ib] 1000
    hi/del 1301
    exec hsigma @1301 = @1300
    exec hsigma u = 1.15*max(vmax(@1301),vmax(@1100))
    null -1 35 0 $sigma(u)
    hi/pl 1301 s
    hi/pl 1100 s
  else
    exec hsigma u1 = @300+%300
    exec hsigma u2 = @100+%100
    u = $sigma(1.2*max(vmax(u1),vmax(u2)))
    null -1 35 0 [u]
    hi/pl 300 s
    hi/pl 100 s
  endif
  nc=$sigma([idh]-10*int([idh]/10))
  dir=$sigma(([idh]-100*int([idh]/100)-[nc])/10)
  np=$sigma(([idh]-1000*int([idh]/1000)-10*[dir]-[nc])/100)
  if ([dir].eq.6) then
    pos='phi'
  else
    pos='z'
  endif
  set chhe 0.2
  text=diff = [xdiff]
  exec $PER/s#tf 0.6 0.8 [text]
  k=3; p=$sigma(int([mue]*10**[k]+0.5)/10**[k])
  k=3; dp=$sigma(int([dmue]*10**[k]+0.5)/10**[k])
  text=A?exp! = [p] "A# [dp]
*"
  exec $PER/s#tf 0.6 0.7 [text]
  k=3; p=$sigma(int([mus]*10**[k]+0.5)/10**[k])
  k=3; dp=$sigma(int([dmus]*10**[k]+0.5)/10**[k])
  text=A?sim! = [p] "A# [dp]
*"
  exec $PER/s#tf 0.6 0.6 [text]
  dk=2
  k=5; p=$sigma(int([effe]*10**[k]+0.5)/10**([k]-[dk]))
  k=5; dp=$sigma(int([deffe]*10**[k]+0.5)/10**([k]-[dk]))
  text=P?0,exp! = ([p] "A# [dp])[\327]10^-[dk]!
*"
  exec $PER/s#tf 0.6 0.5 [text]
  k=5; p=$sigma(int([effs]*10**[k]+0.5)/10**([k]-[dk]))
  k=5; dp=$sigma(int([deffs]*10**[k]+0.5)/10**([k]-[dk]))
  text=P?0,sim! = ([p] "A# [dp])[\327]10^-[dk]!
*"
  exec $PER/s#tf 0.6 0.4 [text]
  atitle 'pusle height, pe' 'events'
  set chhe 0.309
  exec save [mname]/[mname]_counter[nc]_amplitude_vs_$unquote([pos])_pos[np]_v[ver].eps f
*
*  hi/copy 300000 400
*  exec hsigma tot = @400
*  exec hsigma dtot = %400
return



macro effzxpl map=1 ver=1
mname=mapcal[map]
fname0=[mname]_results_v[ver]
fname=[mname]/[fname0]_0.tex
if ($fexist([fname])) then
  shell rm [fname]
endif
for/file 20 [fname]
close 20
*
do nc=1,9
  exec mapcal#effzpl [map] [nc] [ver] [fname]
enddo
*
  flname=latex.sh
  if ($fexist([flname]).eq.1) then
    shell rm [flname]
  endif
  for/file  20 [flname] new
  close 20
  
  txt=cd [mname]
  fmess [txt] [flname]
  txt=latex [fname0]
  fmess [txt] [flname]
  txt=latex [fname0]
  fmess [txt] [flname]
  txt=dvips [fname0]
  fmess [txt] [flname]
  
  shell chmod +x [flname]
  shell ./[flname]
*
return


macro effzpl map=1 nc=1 ver=1 fname=tmp.tex
*
mname=mapcal[map]
*
pos=phi
*
do np=1,3
  fmess '\begin{figure}[ht!b]' [fname]
  fmess '  \begin{minipage}{\textwidth}' [fname]
  fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
  epsfile=[mname]_counter[nc]_amplitude_vs_$unquote([pos])_pos[np]_v[ver].eps
  txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
  fmess [txt] [fname]
  txt=$unquote('     ')\caption{[mname]: counter \No [nc]}
  fmess [txt] [fname]
  fmess '   \end{minipage}' [fname]
  fmess '\end{figure}' [fname]
enddo
*
fmess '\begin{figure}[ht!b]' [fname]
fmess '  \begin{minipage}{\textwidth}' [fname]
fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
epsfile=[mname]_counter[nc]_p0_vs_$unquote([pos])_v[ver].eps
txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
fmess [txt] [fname]
txt=$unquote('     ')\caption{[mname]: counter \No [nc]}
fmess [txt] [fname]
fmess '   \end{minipage}' [fname]
fmess '\end{figure}' [fname]
*
fmess '\begin{figure}[ht!b]' [fname]
fmess '  \begin{minipage}{\textwidth}' [fname]
fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
epsfile=[mname]_counter[nc]_mu_vs_$unquote([pos])_v[ver].eps
txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
fmess [txt] [fname]
txt=$unquote('     ')\caption{[mname]: counter \No [nc]}
fmess [txt] [fname]
fmess '   \end{minipage}' [fname]
fmess '\end{figure}' [fname]
*
fmess '\begin{figure}[ht!b]' [fname]
fmess '  \begin{minipage}{\textwidth}' [fname]
fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
epsfile=[mname]_counter[nc]_mucorr_vs_$unquote([pos])_v[ver].eps
txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
fmess [txt] [fname]
txt=$unquote('     ')\caption{[mname]: counter \No [nc]}
fmess [txt] [fname]
fmess '   \end{minipage}' [fname]
fmess '\end{figure}' [fname]
*
*
pos=z
*
do np=1,5
  fmess '\begin{figure}[ht!b]' [fname]
  fmess '  \begin{minipage}{\textwidth}' [fname]
  fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
  epsfile=[mname]_counter[nc]_amplitude_vs_$unquote([pos])_pos[np]_v[ver].eps
  txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
  fmess [txt] [fname]
  txt=$unquote('     ')\caption{[mname]: counter \No [nc]}
  fmess [txt] [fname]
  fmess '   \end{minipage}' [fname]
  fmess '\end{figure}' [fname]
enddo
*
fmess '\begin{figure}[ht!b]' [fname]
fmess '  \begin{minipage}{\textwidth}' [fname]
fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
epsfile=[mname]_counter[nc]_p0_vs_$unquote([pos])_v[ver].eps
txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
fmess [txt] [fname]
txt=$unquote('     ')\caption{[mname]: counter \No [nc]}
fmess [txt] [fname]
fmess '   \end{minipage}' [fname]
fmess '\end{figure}' [fname]
*
fmess '\begin{figure}[ht!b]' [fname]
fmess '  \begin{minipage}{\textwidth}' [fname]
fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
epsfile=[mname]_counter[nc]_mu_vs_$unquote([pos])_v[ver].eps
txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
fmess [txt] [fname]
txt=$unquote('     ')\caption{[mname]: counter \No [nc]}
fmess [txt] [fname]
fmess '   \end{minipage}' [fname]
fmess '\end{figure}' [fname]
*
fmess '\begin{figure}[ht!b]' [fname]
fmess '  \begin{minipage}{\textwidth}' [fname]
fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
epsfile=[mname]_counter[nc]_mucorr_vs_$unquote([pos])_v[ver].eps
txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
fmess [txt] [fname]
txt=$unquote('     ')\caption{[mname]: counter \No [nc]}
fmess [txt] [fname]
fmess '   \end{minipage}' [fname]
fmess '\end{figure}' [fname]
*
fmess '\begin{figure}[ht!b]' [fname]
fmess '  \begin{minipage}{\textwidth}' [fname]
fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
epsfile=[mname]_counter[nc]_q_vs_ksi_v[ver].eps
txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
fmess [txt] [fname]
txt=$unquote('     ')\caption{[mname]: counter \No [nc]}
fmess [txt] [fname]
fmess '   \end{minipage}' [fname]
fmess '\end{figure}' [fname]
*
fmess '\clearpage' [fname]
return


macro effzxx ver=1 st=1
exec mapcal#effzx 1 1 9 5 [ver] [st]
exec mapcal#effzx 2 1 9 5 [ver] [st]
exec mapcal#effzx 3 1 9 5 [ver] [st]
exec mapcal#effzx 4 1 9 5 [ver] [st]
return


macro effzx map=1 i1=1 i2=9 ib=5 ver=1 st=1
do i=[i1],[i2]
  ve/del crx,dcrx,crsx,dcrsx
  exec mapcal#effz [map] [i] 6 1 [ib] [ver]
  exec vappend   crx   cr
  exec vappend  dcrx  dcr
  exec vappend  crsx  crs
  exec vappend dcrsx dcrs
  exec mapcal#effz [map] [i] 7 1 [ib] [ver]
  exec vappend   crx   cr
  exec vappend  dcrx  dcr
  exec vappend  crsx  crs
  exec vappend dcrsx dcrs
*  
  sigma crx = order(crx,crsx)
  sigma dcrx = order(dcrx,crsx)
  sigma dcrsx = order(dcrsx,crsx)
  sigma crsx = order(crsx,crsx)
  set pmci 4
  * exec $PER/s#vpl crx dcrx crsx dcrsx sz=0.1 ll=1
  d0=$GRAFINFO('WNYMIN')
  u0=$GRAFINFO('WNYMAX')
  d=$sigma(([u0]+[d0])/2-3*([u0]-[d0])/2)
  u=$sigma(([u0]+[d0])/2+3*([u0]-[d0])/2)
  null 1 $GRAFINFO('WNXMAX') [d] [u]
  * exec $PER/s#vpl crx dcrx crsx dcrsx sz=0.1 o=s
  set pmci 2
  * exec $PER/s#vpl cr dcr crs dcrs sz=0.1 o=s
*
  ve/cre p1(1) r 1
  ve/fit crsx crx dcrx p0 s 1 p1
  l=p1(1)
  line $GRAFINFO('WNXMIN') [l] $GRAFINFO('WNXMAX') [l]
  ve/cre p1(1) r 1
  ve/cre dp1(1) r
  ve/fit crsx crx dcrx qx.f s 1 p1 ! ! ! dp1
*  
  fun/pl qxp.f 0 $GRAFINFO('WNXMAX') s
  set basl 0.01
  set ltyp 15
  l=p1(1)
  line $GRAFINFO('WNXMIN') [l] $GRAFINFO('WNXMAX') [l]
  set ltyp 1
*
  gl/imp mname
  exec save [mname]/[mname]_counter[i]_q_vs_ksi_v[ver].eps f
  ve/write p1,dp1 [mname]/[mname]_counter[i]_amplitude_correction_st[st].txt '(2f15.6)'
*  read x
enddo
return

macro ampcorrmapx st=1
do i=0,4
  exec mapcal#ampcorrmap [i] [st]
enddo
return

macro ampcorrmap map=1 st=1
ve/cre acorr(9) r 9*1
ve/del p1,dp1,p2,dp2
if ([map].gt.0) then
  mname=mapcal[map]
  do i=1,9
    acorri=1
    do j=1,[st]
      ve/read p1,dp1 [mname]/[mname]_counter[i]_amplitude_correction_st[j].txt '(2f15.6)'
      acorri=$sigma([acorri]*p1(1))
    enddo
    ve/inp acorr([i]) [acorri]
  enddo
endif
ve/write acorr AmplitudeCorrection_mapcal[map].txt
return

macro effz map=1 nc=1 np=3 auto=0 ib=1 ver=0
opt nstat
set mscf 0.3
set pmci 1
set hcol 1
set plci 1
*
ni=5
pos='z'
if ([np].eq.6) then
  ni=3
  pos='phi'
endif
ve/cre   mue([ni]) r
ve/cre   mus([ni]) r
ve/cre   mue1([ni]) r
ve/cre   mus1([ni]) r
ve/cre   mue2([ni]) r
ve/cre   mus2([ni]) r
ve/cre  dmus2([ni]) r
ve/cre  dmue([ni]) r
ve/cre  dmus([ni]) r
ve/cre  effe([ni]) r
ve/cre deffe([ni]) r
ve/cre  effs([ni]) r
ve/cre deffs([ni]) r
do i=1,[ni]
  exec mapcal#expsimcn map=[map] idh=$sigma(100000+[i]*100+[np]*10+[nc]) ib=[ib] ver=[ver]
*  exec mapcal#expsimc $sigma(100000+[i]*100+[np]*10+[nc]) ib=[ib]
*  mess exec mapcal#expsimc $sigma(100000+[i]*100+[np]*10+[nc])
  gl/imp  effe
  gl/imp deffe
  gl/imp  effs
  gl/imp deffs
  ve/inp  effe([i])  [effe]
  ve/inp deffe([i]) [deffe]
  ve/inp  effs([i])  [effs]
  ve/inp deffs([i]) [deffs]
  gl/imp mue
  gl/imp mus
  ve/inp mue([i]) [mue]
  ve/inp mus([i]) [mus]
  gl/imp dmue
  gl/imp dmus
  ve/inp dmue([i]) [dmue]
  ve/inp dmus([i]) [dmus]
  exec hsigma mean = vsum( ( @100 * ($100 gt 0.2) ) * $100 )/vsum( @100 * ($100 gt 0.2) )
  ve/inp mue1([i]) $sigma(mean(1))
  exec hsigma mean = vsum( ( @200 * ($200 gt 0.2) ) * $200 )/vsum( @200 * ($200 gt 0.2) )
  ve/inp mus1([i]) $sigma(mean(1))
  if ([np].eq.0) then
    ns=[i]
  endif
  if ([np].eq.3) then
    ns=[i]+3
  endif
  if ([np].eq.6) then
    ns=[i]+9
  endif
  if ([np].eq.7) then
    ns=[i]+12
  endif
  ind=([nc]-1)*18+[ns]+1
*  ve/inp mus2([i]) $sigma(tot([ind]))
*  ve/inp dmus2([i]) $sigma(dtot([ind]))
  if ([auto].eq.0) then
    read x
  endif
enddo
*
gl/imp mname
*
ve/cre nz([ni]) r
sigma nz = array([ni],1#[ni])
ve/cre dnz([ni]) r
sigma ceffs = mus*exp(-mus)/10
sigma effsc = effs + ceffs
set pmci 4
set hcol 1
sz=0.2
sigma nzs = nz+0.05
if ($sigma(vmax(effe+deffe)).gt.$sigma(vmax(effs+deffs))) then
  * exec $PER/s#vpl effe deffe nz dnz iatt=20 sz=[sz]
  * exec $PER/s#vpl effs deffs nzs dnz iatt=24 sz=[sz] o=s
*  * exec $PER/s#vpl effsc dnz nz dnz iatt=21 sz=[sz] o=s
*  * exec $PER/s#vpl ceffs dnz nzs dnz iatt=23 sz=[sz] o=s
else
  * exec $PER/s#vpl effs deffs nz dnz iatt=24 sz=[sz] 
  * exec $PER/s#vpl effe deffe nzs dnz iatt=20 sz=[sz] o=s
*  * exec $PER/s#vpl effsc dnz nz dnz iatt=21 sz=[sz] o=s
*  * exec $PER/s#vpl ceffs dnz nzs dnz iatt=23 sz=[sz] o=s
endif
txt=$unquote([pos])-position
atitle [txt] 'P?0!'
exec save [mname]/[mname]_counter[nc]_p0_vs_$unquote([pos])_v[ver].eps f
  if ([auto].eq.0) then
    read x
  endif
if ($sigma(vmax(mue+deffe)).gt.$sigma(vmax(mus+deffs))) then
  * exec $PER/s#vpl mue dmue nz dnz iatt=20 sz=[sz]
*  * exec $PER/s#vpl mus2 dmus2 nz dnz iatt=21 sz=[sz] o=s
*  set pmci 2
  * exec $PER/s#vpl mus dmus nzs dnz iatt=24 sz=[sz] o=s
else
  * exec $PER/s#vpl mus dmus nz dnz iatt=24 sz=[sz] 
  * exec $PER/s#vpl mue dmue nzs dnz iatt=20 sz=[sz] o=s
*  set pmci 2
*  * exec $PER/s#vpl mus2 dmus2 nz dnz iatt=20 sz=[sz] o=s
*  set pmci 4
endif
txt=$unquote([pos])-position
atitle [txt] '[m]'
exec save [mname]/[mname]_counter[nc]_mu_vs_$unquote([pos])_v[ver].eps f
  if ([auto].eq.0) then
    read x
  endif
sigma rmu = mue/mus
sigma drmu = rmu*sqrt((dmue/mue)**2+(dmus/mus)**2)
* exec $PER/s#vpl rmu drmu nz dnz iatt=20 sz=[sz]
atitle [txt] '[m]?exp!/[m]?sim!'
npar=1
ve/cre chi2(2) r
ve/cre paru([npar]) r
ve/cre dparu([npar]) r
ve/cre covu([npar],[npar]) r
ve/fit nz rmu drmu p0 s
call covm.f(1)
call covmpen(chi2,[npar],paru,dparu)
call covmcov([npar],covu)
rt=$sigma(int(paru(1)*10000+0.5)/10000)
drt=$sigma(int(dparu(1)*10000+0.5)/10000)
text=[m]?exp!/[m]?sim! = [rt] "A# [drt]
*"
exec $PER/s#tf 0.5 0.5 [text]
ndf=$sigma(chi2(2))
chi=$sigma(int(chi2(1)*[ndf]*1000+0.5)/1000)
txt=[h]^2!/ndf = [chi] / [ndf]
exec $PER/s#tf 0.5 0.6 [txt]
ve/prin chi2
  if ([auto].eq.0) then
    read x
  endif
*
* exec $PER/s#vpl mus dmus nz dnz iatt=24 sz=[sz]
* exec $PER/s#vpl mue1 dmue nzs dnz iatt=20 sz=[sz] o=s
  if ([auto].eq.0) then
    read x
  endif
*
sigma cr = log(effe)/log(effs)
sigma crc = log(effe)/log(effsc)
sigma dcr = cr*sqrt((deffe/(effe*log(effe)))**2+(deffs/(effs*log(effs)))**2)
** exec $PER/s#vpl crc dcr nz dnz iatt=24 sz=[sz]
* exec $PER/s#vpl cr dcr nzs dnz iatt=20 sz=[sz]

txt=$unquote([pos])-position
atitle [txt] '[m]?exp!/[m]?sim!'

npar=1
ve/cre chi2(2) r
ve/cre paru([npar]) r
ve/cre dparu([npar]) r
ve/cre covu([npar],[npar]) r

ve/cre corra(1) r 1
ve/cre dcorra(1) r
ve/fit nz cr dcr p0 s 1 corra ! ! ! dcorra

call covm.f(1)
call covmpen(chi2,[npar],paru,dparu)
call covmcov([npar],covu)
rt=$sigma(int(paru(1)*1000+0.5)/1000)
drt=$sigma(int(dparu(1)*1000+0.5)/1000)
text=[m]?exp!/[m]?sim! = [rt] "A# [drt]
*"
exec $PER/s#tf 0.5 0.5 [text]
ndf=$sigma(chi2(2))
chi=$sigma(int(chi2(1)*[ndf]*100+0.5)/100)
txt=[h]^2!/ndf = [chi] / [ndf]
exec $PER/s#tf 0.5 0.6 [txt]

exec save [mname]/[mname]_counter[nc]_mucorr_vs_$unquote([pos])_v[ver].eps f

*
mess corr = $sigma(vsum(cr)/[ni]) $sigma(sqrt(vsum(dcr**2))/[ni])

  if ([auto].eq.0) then
    read x
  endif
*
sigma uf = 2*(mus+log(effs))/mus**2
ufm = $sigma(vsum(uf)/[ni])
sigma muec = (1-sqrt(1+2*log(effe)*[ufm]))/[ufm]
null 0 15 0 15
line 0 0 15 15
* exec $PER/s#vpl mue dmue muec dmue o=s sz=[sz]
ve/cre p2(2) r 1 0
ve/fit muec mue dmue musat.f s 2 p2
*
ve/del crs,dcrs
sigma  crs = -mus/log(effs)
sigma dcrs = 0*crs*sqrt((dmus/mus)**2+(deffs/effs)**2)
sigma cr = order(cr,crs)
sigma dcr = order(dcr,crs)
sigma dcrs = order(dcrs,crs)
sigma crs = order(crs,crs)
set pmci 4
* exec $PER/s#vpl cr dcr crs dcrs sz=0.1
ve/cre p1(1) r 1
ve/fit crs cr dcr  qx.f s 1 p1
return



macro effzn map=1 nc=1 np=3 auto=0 ib=1
opt nstat
set mscf 0.3
set pmci 1
*
if ([map].eq.1) then
  gl/cre mname mapcal1
  nh1=11
  nh2=11
  npt=0
endif
if ([map].eq.2) then
  gl/cre mname mapcal2
  nh1=12
  nh2=12
  npt=1
endif
if ([map].eq.3) then
  gl/cre mname mapcal3
  nh1=13
  nh2=13
  npt=2
endif
if ([map].eq.4) then
  gl/cre mname mapcal4
  nh1=14
  nh2=14
  npt=3
endif
*
ni=5
pos='phi'
if ([np].eq.6) then
  ni=3
  pos='z'
endif
ve/cre   mue([ni]) r
ve/cre   mus([ni]) r
ve/cre   mue1([ni]) r
ve/cre   mus1([ni]) r
ve/cre   mue2([ni]) r
ve/cre   mus2([ni]) r
ve/cre  dmus2([ni]) r
ve/cre  dmue([ni]) r
ve/cre  dmus([ni]) r
ve/cre  effe([ni]) r
ve/cre deffe([ni]) r
ve/cre  effs([ni]) r
ve/cre deffs([ni]) r
do i=1,[ni]
  exec mapcal#expsimc $sigma(100000+[i]*100+[np]*10+[nc]) ib=[ib]
*  mess exec mapcal#expsimc $sigma(100000+[i]*100+[np]*10+[nc])
  gl/imp  effe
  gl/imp deffe
  gl/imp  effs
  gl/imp deffs
  ve/inp  effe([i])  [effe]
  ve/inp deffe([i]) [deffe]
  ve/inp  effs([i])  [effs]
  ve/inp deffs([i]) [deffs]
  gl/imp mue
  gl/imp mus
  ve/inp mue([i]) [mue]
  ve/inp mus([i]) [mus]
  gl/imp dmue
  gl/imp dmus
  ve/inp dmue([i]) [dmue]
  ve/inp dmus([i]) [dmus]
  exec hsigma mean = vsum( ( @100 * ($100 gt 0.2) ) * $100 )/vsum( @100 * ($100 gt 0.2) )
  ve/inp mue1([i]) $sigma(mean(1))
  exec hsigma mean = vsum( ( @200 * ($200 gt 0.2) ) * $200 )/vsum( @200 * ($200 gt 0.2) )
  ve/inp mus1([i]) $sigma(mean(1))
  if ([np].eq.0) then
    ns=[i]
  endif
  if ([np].eq.3) then
    ns=[i]+3
  endif
  if ([np].eq.6) then
    ns=[i]+9
  endif
  if ([np].eq.7) then
    ns=[i]+12
  endif
  ind=([nc]-1)*18+[ns]+1
*  ve/inp mus2([i]) $sigma(tot([ind]))
*  ve/inp dmus2([i]) $sigma(dtot([ind]))
  if ([auto].eq.0) then
    read x
  endif
enddo
ve/cre nz([ni]) r
sigma nz = array([ni],1#[ni])
ve/cre dnz([ni]) r
sigma ceffs = mus*exp(-mus)/10
sigma effsc = effs + ceffs
set pmci 4
set hcol 1
sz=0.2
sigma nzs = nz+0.05
if ($sigma(vmax(effe+deffe)).gt.$sigma(vmax(effs+deffs))) then
  * exec $PER/s#vpl effe deffe nz dnz iatt=20 sz=[sz]
  * exec $PER/s#vpl effs deffs nzs dnz iatt=24 sz=[sz] o=s
*  * exec $PER/s#vpl effsc dnz nz dnz iatt=21 sz=[sz] o=s
*  * exec $PER/s#vpl ceffs dnz nzs dnz iatt=23 sz=[sz] o=s
else
  * exec $PER/s#vpl effs deffs nz dnz iatt=24 sz=[sz] 
  * exec $PER/s#vpl effe deffe nzs dnz iatt=20 sz=[sz] o=s
*  * exec $PER/s#vpl effsc dnz nz dnz iatt=21 sz=[sz] o=s
*  * exec $PER/s#vpl ceffs dnz nzs dnz iatt=23 sz=[sz] o=s
endif
txt=$unquote([pos])-position
atitle [txt] 'P?0!'
exec save [mname]/[mname]_counter[nc]_p0_vs_$unquote([pos])_v[ver].eps f
  if ([auto].eq.0) then
    read x
  endif
if ($sigma(vmax(mue+deffe)).gt.$sigma(vmax(mus+deffs))) then
  * exec $PER/s#vpl mue dmue nz dnz iatt=20 sz=[sz]
*  * exec $PER/s#vpl mus2 dmus2 nz dnz iatt=21 sz=[sz] o=s
*  set pmci 2
  * exec $PER/s#vpl mus dmus nzs dnz iatt=24 sz=[sz] o=s
else
  * exec $PER/s#vpl mus dmus nz dnz iatt=24 sz=[sz] 
  * exec $PER/s#vpl mue dmue nzs dnz iatt=20 sz=[sz] o=s
*  set pmci 2
*  * exec $PER/s#vpl mus2 dmus2 nz dnz iatt=20 sz=[sz] o=s
*  set pmci 4
endif
txt=$unquote([pos])-position
atitle [txt] '[m]'
exec save [mname]/[mname]_counter[nc]_mu_vs_$unquote([pos])_v[ver].eps f
  if ([auto].eq.0) then
    read x
  endif
sigma rmu = mue/mus
sigma drmu = rmu*sqrt((dmue/mue)**2+(dmus/mus)**2)
* exec $PER/s#vpl rmu drmu nz dnz iatt=20 sz=[sz]
atitle [txt] '[m]?exp!/[m]?sim!'
npar=1
ve/cre chi2(2) r
ve/cre paru([npar]) r
ve/cre dparu([npar]) r
ve/cre covu([npar],[npar]) r
ve/fit nz rmu drmu p0 s
call covm.f(1)
call covmpen(chi2,[npar],paru,dparu)
call covmcov([npar],covu)
rt=$sigma(int(paru(1)*10000+0.5)/10000)
drt=$sigma(int(dparu(1)*10000+0.5)/10000)
text=[m]?exp!/[m]?sim! = [rt] "A# [drt]
*"
exec $PER/s#tf 0.5 0.5 [text]
ndf=$sigma(chi2(2))
chi=$sigma(int(chi2(1)*[ndf]*1000+0.5)/1000)
txt=[h]^2!/ndf = [chi] / [ndf]
exec $PER/s#tf 0.5 0.6 [txt]
ve/prin chi2
  if ([auto].eq.0) then
    read x
  endif
*
* exec $PER/s#vpl mus dmus nz dnz iatt=24 sz=[sz]
* exec $PER/s#vpl mue1 dmue nzs dnz iatt=20 sz=[sz] o=s
  if ([auto].eq.0) then
    read x
  endif
*
sigma cr = log(effe)/log(effs)
sigma crc = log(effe)/log(effsc)
sigma dcr = cr*sqrt((deffe/(effe*log(effe)))**2+(deffs/(effs*log(effs)))**2)
** exec $PER/s#vpl crc dcr nz dnz iatt=24 sz=[sz]
* exec $PER/s#vpl cr dcr nzs dnz iatt=20 sz=[sz]

txt=$unquote([pos])-position
atitle [txt] '[m]?exp!/[m]?sim!'

npar=1
ve/cre chi2(2) r
ve/cre paru([npar]) r
ve/cre dparu([npar]) r
ve/cre covu([npar],[npar]) r

ve/cre corra(1) r 1
ve/cre dcorra(1) r
ve/fit nz cr dcr p0 s 1 corra ! ! ! dcorra

call covm.f(1)
call covmpen(chi2,[npar],paru,dparu)
call covmcov([npar],covu)
rt=$sigma(int(paru(1)*1000+0.5)/1000)
drt=$sigma(int(dparu(1)*1000+0.5)/1000)
text=[m]?exp!/[m]?sim! = [rt] "A# [drt]
*"
exec $PER/s#tf 0.5 0.5 [text]
ndf=$sigma(chi2(2))
chi=$sigma(int(chi2(1)*[ndf]*100+0.5)/100)
txt=[h]^2!/ndf = [chi] / [ndf]
exec $PER/s#tf 0.5 0.6 [txt]

exec save [mname]/[mname]_counter[nc]_mucorr_vs_$unquote([pos])_v[ver].eps f

*
mess corr = $sigma(vsum(cr)/[ni]) $sigma(sqrt(vsum(dcr**2))/[ni])

  if ([auto].eq.0) then
    read x
  endif
*
sigma uf = 2*(mus+log(effs))/mus**2
ufm = $sigma(vsum(uf)/[ni])
sigma muec = (1-sqrt(1+2*log(effe)*[ufm]))/[ufm]
null 0 15 0 15
line 0 0 15 15
* exec $PER/s#vpl mue dmue muec dmue o=s sz=[sz]
ve/cre p2(2) r 1 0
ve/fit muec mue dmue musat.f s 2 p2
return



macro mutest
ve/cre mus(30) r 
ve/cre pds(9) r 478.142 385.294 448.022 429.255 520.526 591.312 421.836 419.288 660.439
ve/cre a1q(9) r 90.043 77.241 105.399 76.985 97.676 114.848 94.784 94.174  98.95
pds1=pds(1)
amp1=a1q(1)
do i=0,29
  nt/cre 1 ! 2 ! ! npe amp
  nt/read 1 /work/users/konctbel/snd2k/R005-999/mu[i].txt
  exec ntproj $sigma(10000+[i]) 1.npe ! 102 -1.5 100.5 
  exec ntproj $sigma(20000+[i]) 1.amp ! 4000 0 4000 
  hi/del 1
  j=[i]+1
  idh=$sigma(20000+[i])
  ampx=$hinfo([idh],'mean')
  ve/inp mus([j]) $sigma(([ampx]-[pds1])/[amp1])
enddo
ve/cre mus0(30) r 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
return


macro mutest0 cnt=1
ve/cre mus(30) r 
*ve/cre pds(9) r 478.142 385.294 448.022 429.255 520.526 591.312 421.836 419.288 660.439
*ve/cre a1q(9) r 90.043 77.241 105.399 76.985 97.676 114.848 94.784 94.174  98.95
ve/cre pds(9) r 478.142 385.294 448.022 429.255 520.526 591.312 421.836 419.288 660.439
ve/cre a1q(9) r 85.553 73.093 96.013 72.988 94.227 92.113 82.458 136.194 102.672
pds1=pds([cnt])
amp1=a1q([cnt])
do i=0,29
  nt/cre 1 ! 3 ! ! mui npe amp
  s="
*"  
  txt=[s]mumu[i] [s]
  mess [txt]
  shell fgrep -e [txt] /work/users/konctbel/snd2k/R005-999/cnt[cnt].log > tmp.txt
  shell $unquote('cat tmp.txt | sed "s/mumu/ /g" > tmp1.txt')
  nt/read 1 tmp1.txt
  exec ntproj $sigma(10000+[i]) 1.npe ! 102 -1.5 100.5 
  exec ntproj $sigma(20000+[i]) 1.amp ! 4000 0 4000 
  hi/del 1
  j=[i]+1
  idh=$sigma(20000+[i])
  ampx=$hinfo([idh],'mean')
  ve/inp mus([j]) $sigma([ampx]/[amp1])
enddo
ve/cre mus0(30) r 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
ve/cre dmu(30) r
ve/fit mus0 mus dmu p1
graph 30 mus0 mus *
return


macro modx
chain -mod
chain mod ee_t35e10-510-10943-2516-100000.hbook
chain mod ee_t35e10-510-10943-2501-100000.hbook
chain mod ee_t35e10-510-10943-2518-100000.hbook
chain mod ee_t35e10-510-10943-2514-100000.hbook
chain mod ee_t35e10-510-10943-2504-100000.hbook
chain mod ee_t35e10-510-10943-2502-100000.hbook
chain mod ee_t35e10-510-10943-2511-100000.hbook
chain mod ee_t35e10-510-10943-2515-100000.hbook
chain mod ee_t35e10-510-10943-2513-100000.hbook
chain mod ee_t35e10-510-10943-2508-100000.hbook
chain mod ee_t35e10-510-10943-2503-100000.hbook
chain mod ee_t35e10-510-10943-2509-100000.hbook
chain mod ee_t35e10-510-10943-2512-100000.hbook
chain mod ee_t35e10-510-10943-2506-100000.hbook
chain mod ee_t35e10-510-10943-2519-100000.hbook
chain mod ee_t35e10-510-10943-2505-100000.hbook
chain mod ee_t35e10-510-10943-2507-100000.hbook
chain mod ee_t35e10-510-10943-2510-100000.hbook
chain mod ee_t35e10-510-10943-2500-100000.hbook
chain mod ee_t35e10-510-10943-2517-100000.hbook
chain mod ee_t10e10-510-10943-2505-50000.hbook
chain mod ee_t10e10-510-10943-2500-50000.hbook
chain mod ee_t10e10-510-10943-2510-50000.hbook
chain mod ee_t10e10-510-10943-2515-50000.hbook
chain mod ee_t10e10-510-10943-2506-50000.hbook
chain mod ee_t10e10-510-10943-2501-50000.hbook
chain mod ee_t10e10-510-10943-2511-50000.hbook
chain mod ee_t10e10-510-10943-2516-50000.hbook
chain mod ee_t10e10-510-10943-2507-50000.hbook
chain mod ee_t10e10-510-10943-2502-50000.hbook
chain mod ee_t10e10-510-10943-2512-50000.hbook
chain mod ee_t10e10-510-10943-2517-50000.hbook
chain mod ee_t10e10-510-10943-2508-50000.hbook
chain mod ee_t10e10-510-10943-2503-50000.hbook
chain mod ee_t10e10-510-10943-2513-50000.hbook
chain mod ee_t10e10-510-10943-2518-50000.hbook
chain mod ee_t10e10-510-10943-2509-50000.hbook
chain mod ee_t10e10-510-10943-2504-50000.hbook
chain mod ee_t10e10-510-10943-2514-50000.hbook
chain mod ee_t10e10-510-10943-2519-50000.hbook
chain mod -P /work/users/konctbel/snd2k/R005-999/
return


macro fit8
    n=0
*    if ([nc].eq.1) then
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter1_slix41_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter1_slix57_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter1_slix73_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter1_slix89_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter1_slix105_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter1_slix121_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter1_slix136_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter1_slix152_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter1_slix168_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter1_slix184_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter1_slix200_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter1_slix216_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter1_slix232_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter1_slix247_fit_8.eps
*    endif
*    if ([nc].eq.2) then
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter2_slix43_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter2_slix58_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter2_slix74_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter2_slix90_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter2_slix106_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter2_slix121_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter2_slix137_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter2_slix153_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter2_slix169_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter2_slix184_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter2_slix200_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter2_slix216_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter2_slix232_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter2_slix248_fit_8.eps
*    endif
*    if ([nc].eq.3) then
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter3_slix40_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter3_slix55_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter3_slix71_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter3_slix87_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter3_slix103_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter3_slix118_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter3_slix134_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter3_slix150_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter3_slix166_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter3_slix181_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter3_slix197_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter3_slix213_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter3_slix229_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter3_slix244_fit_8.eps
*    endif
*    if ([nc].eq.4) then
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter4_slix43_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter4_slix59_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter4_slix74_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter4_slix90_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter4_slix106_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter4_slix121_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter4_slix137_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter4_slix153_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter4_slix168_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter4_slix184_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter4_slix199_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter4_slix215_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter4_slix231_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter4_slix246_fit_8.eps
*    endif
*    if ([nc].eq.5) then
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter5_slix42_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter5_slix58_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter5_slix74_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter5_slix90_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter5_slix106_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter5_slix121_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter5_slix137_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter5_slix153_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter5_slix169_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter5_slix185_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter5_slix200_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter5_slix216_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter5_slix232_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter5_slix248_fit_8.eps
*    endif
*    if ([nc].eq.6) then
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter6_slix42_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter6_slix57_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter6_slix73_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter6_slix89_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter6_slix105_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter6_slix120_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter6_slix136_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter6_slix152_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter6_slix168_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter6_slix183_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter6_slix199_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter6_slix215_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter6_slix231_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter6_slix246_fit_8.eps
*    endif
*    if ([nc].eq.7) then
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter7_slix44_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter7_slix59_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter7_slix75_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter7_slix91_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter7_slix106_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter7_slix122_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter7_slix138_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter7_slix153_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter7_slix169_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter7_slix185_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter7_slix200_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter7_slix216_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter7_slix232_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter7_slix248_fit_8.eps
*    endif
*    if ([nc].eq.8) then
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter8_slix41_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter8_slix57_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter8_slix72_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter8_slix88_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter8_slix103_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter8_slix119_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter8_slix134_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter8_slix150_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter8_slix166_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter8_slix181_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter8_slix197_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter8_slix212_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter8_slix228_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter8_slix243_fit_8.eps
*    endif
*    if ([nc].eq.9) then
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter9_slix42_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter9_slix58_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter9_slix73_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter9_slix89_fit_8.eps  
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter9_slix105_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter9_slix120_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter9_slix136_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter9_slix151_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter9_slix167_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter9_slix182_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter9_slix198_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter9_slix213_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter9_slix229_fit_8.eps 
      n=[n]+1; file[n]=mhad2011_amplitude_vs_phi_counter9_slix245_fit_8.eps 
*    endif
*
fname=mhad2011_amplitude_vs_phi_counters_slix_fit_8.tex
if ($fexist([fname]).eq.1) then
  shell rm [fname]
endif
for/file 20 [fname] new
close 20
do i=1,[n]
fmess '\begin{figure}[ht!b]' [fname]
fmess '\begin{minipage}{\textwidth}' [fname]
fmess '{\centering\resizebox*{\textwidth}{!}' [fname]
txt={\includegraphics{[file[i]]}}\par}
fmess [txt] [fname]
ncnt=$sigma(int(([i]-1)/14)+1)
txt=\caption{Counter \No [ncnt], slise $sigma([i]-14*([ncnt]-1))}
fmess [txt] [fname]
fmess '\end{minipage}' [fname]
fmess '\end{figure}' [fname]
fmess '' [fname]
if ($sigma(mod([i],14)).eq.0) then
fmess '\clearpage' [fname]
fmess '' [fname]
endif
enddo
return


macro fit8n
    n=0
*    if ([nc].eq.1) then
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter1_slix40_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter1_slix55_fit_7_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter1_slix70_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter1_slix86_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter1_slix104_fit_9_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter1_slix121_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter1_slix137_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter1_slix153_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter1_slix170_fit_9_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter1_slix187_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter1_slix203_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter1_slix219_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter1_slix235_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter1_slix249_fit_6_et0_chit.eps
*    endif
*    if ([nc].eq.2) then
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter2_slix40_fit_7_et0_chit.eps   
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter2_slix56_fit_8_et0_chit.eps   
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter2_slix72_fit_8_et0_chit.eps   
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter2_slix87_fit_8_et0_chit.eps   
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter2_slix103_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter2_slix120_fit_9_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter2_slix138_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter2_slix153_fit_7_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter2_slix168_fit_9_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter2_slix184_fit_7_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter2_slix199_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter2_slix217_fit_10_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter2_slix234_fit_7_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter2_slix248_fit_8_et0_chit.eps 
*    endif
*    if ([nc].eq.3) then
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix38_fit_7_et0_chit.eps   
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix53_fit_8_et0_chit.eps   
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix69_fit_8_et0_chit.eps   
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix87_fit_10_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix105_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix119_fit_6_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix133_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix149_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix166_fit_9_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix181_fit_7_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix196_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix214_fit_10_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix232_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix247_fit_7_et0_chit.eps 
*    endif
*    if ([nc].eq.4) then
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix38_fit_7_et0_chit.eps   
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix53_fit_8_et0_chit.eps   
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix69_fit_8_et0_chit.eps   
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix87_fit_10_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix105_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix119_fit_6_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix133_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix149_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix166_fit_9_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix181_fit_7_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix196_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix214_fit_10_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix232_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix247_fit_7_et0_chit.eps 
*    endif
*    if ([nc].eq.5) then
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix38_fit_7_et0_chit.eps   
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix53_fit_8_et0_chit.eps   
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix69_fit_8_et0_chit.eps   
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix87_fit_10_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix105_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix119_fit_6_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix133_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix149_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix166_fit_9_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix181_fit_7_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix196_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix214_fit_10_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix232_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix247_fit_7_et0_chit.eps 
*    endif
*    if ([nc].eq.6) then
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter6_slix41_fit_10_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter6_slix58_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter6_slix74_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter6_slix88_fit_6_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter6_slix103_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter6_slix121_fit_9_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter6_slix136_fit_6_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter6_slix152_fit_9_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter6_slix170_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter6_slix186_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter6_slix201_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter6_slix217_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter6_slix232_fit_7_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter6_slix246_fit_7_et0_chit.eps
*    endif
*    if ([nc].eq.7) then
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter7_slix41_fit_8_et0_chit.eps   
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter7_slix55_fit_6_et0_chit.eps   
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter7_slix69_fit_8_et0_chit.eps   
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter7_slix86_fit_9_et0_chit.eps   
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter7_slix103_fit_9_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter7_slix118_fit_6_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter7_slix131_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter7_slix149_fit_10_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter7_slix164_fit_6_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter7_slix178_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter7_slix193_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter7_slix210_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter7_slix227_fit_9_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter7_slix244_fit_8_et0_chit.eps 
*    endif
*    if ([nc].eq.8) then
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter8_slix40_fit_10_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter8_slix58_fit_8_et0_chit.eps   
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter8_slix73_fit_7_et0_chit.eps   
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter8_slix88_fit_9_et0_chit.eps   
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter8_slix105_fit_7_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter8_slix118_fit_6_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter8_slix134_fit_9_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter8_slix150_fit_7_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter8_slix165_fit_7_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter8_slix180_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter8_slix197_fit_10_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter8_slix214_fit_6_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter8_slix229_fit_9_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter8_slix245_fit_7_et0_chit.eps 
*    endif
*    if ([nc].eq.9) then
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter9_slix39_fit_10_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter9_slix57_fit_8_et0_chit.eps   
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter9_slix71_fit_6_et0_chit.eps   
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter9_slix85_fit_9_et0_chit.eps   
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter9_slix103_fit_9_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter9_slix121_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter9_slix135_fit_6_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter9_slix149_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter9_slix165_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter9_slix183_fit_10_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter9_slix200_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter9_slix216_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter9_slix232_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter9_slix246_fit_6_et0_chit.eps 
*    endif
*
fname=mhad2011-4_amplitude_vs_phi_counters_slix_fit_8.tex
if ($fexist([fname]).eq.1) then
  shell rm [fname]
endif
for/file 20 [fname] new
close 20
do i=1,[n]
fmess '\begin{figure}[ht!b]' [fname]
fmess '\begin{minipage}{\textwidth}' [fname]
fmess '{\centering\resizebox*{\textwidth}{!}' [fname]
txt={\includegraphics{[file[i]]}}\par}
fmess [txt] [fname]
ncnt=$sigma(int(([i]-1)/14)+1)
txt=\caption{Counter \No [ncnt], slise $sigma([i]-14*([ncnt]-1))}
fmess [txt] [fname]
fmess '\end{minipage}' [fname]
fmess '\end{figure}' [fname]
fmess '' [fname]
if ($sigma(mod([i],14)).eq.0) then
fmess '\clearpage' [fname]
fmess '' [fname]
endif
enddo
return


macro fit8nm
    n=0
*    if ([nc].eq.1) then
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter1_slix41_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter1_slix56_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter1_slix71_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter1_slix87_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter1_slix104_fit_9_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter1_slix121_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter1_slix136_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter1_slix152_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter1_slix168_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter1_slix184_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter1_slix200_fit_7_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter1_slix215_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter1_slix232_fit_9_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter1_slix248_fit_7_et0_chit.eps
*    endif
*    if ([nc].eq.2) then
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter2_slix42_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter2_slix58_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter2_slix73_fit_7_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter2_slix89_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter2_slix105_fit_7_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter2_slix120_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter2_slix137_fit_9_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter2_slix153_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter2_slix170_fit_9_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter2_slix185_fit_7_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter2_slix200_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter2_slix216_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter2_slix233_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter2_slix248_fit_7_et0_chit.eps
*    endif
*    if ([nc].eq.3) then
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix40_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix55_fit_7_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix70_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix87_fit_9_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix104_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix120_fit_7_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix135_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix151_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix168_fit_9_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix183_fit_7_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix198_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix215_fit_9_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix231_fit_7_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter3_slix247_fit_8_et0_chit.eps
*    endif
*    if ([nc].eq.4) then
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter4_slix41_fit_9_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter4_slix56_fit_6_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter4_slix70_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter4_slix87_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter4_slix104_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter4_slix120_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter4_slix136_fit_7_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter4_slix152_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter4_slix168_fit_9_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter4_slix184_fit_7_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter4_slix200_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter4_slix214_fit_7_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter4_slix231_fit_9_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter4_slix248_fit_7_et0_chit.eps
*    endif
*    if ([nc].eq.5) then
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter5_slix42_fit_10_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter5_slix58_fit_6_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter5_slix73_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter5_slix89_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter5_slix105_fit_9_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter5_slix122_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter5_slix138_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter5_slix154_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter5_slix170_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter5_slix186_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter5_slix202_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter5_slix218_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter5_slix233_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter5_slix249_fit_8_et0_chit.eps
*    endif
*    if ([nc].eq.6) then
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter6_slix40_fit_9_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter6_slix57_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter6_slix73_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter6_slix88_fit_7_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter6_slix103_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter6_slix120_fit_9_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter6_slix135_fit_7_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter6_slix151_fit_9_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter6_slix169_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter6_slix183_fit_7_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter6_slix198_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter6_slix215_fit_9_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter6_slix230_fit_7_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter6_slix245_fit_8_et0_chit.eps
*    endif
*    if ([nc].eq.7) then
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter7_slix41_fit_9_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter7_slix57_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter7_slix73_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter7_slix89_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter7_slix104_fit_7_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter7_slix119_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter7_slix134_fit_7_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter7_slix150_fit_9_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter7_slix166_fit_7_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter7_slix181_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter7_slix196_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter7_slix213_fit_9_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter7_slix231_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter7_slix246_fit_7_et0_chit.eps
*    endif
*    if ([nc].eq.8) then
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter8_slix40_fit_9_et0_chit.eps   
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter8_slix57_fit_8_et0_chit.eps   
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter8_slix73_fit_8_et0_chit.eps   
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter8_slix88_fit_7_et0_chit.eps   
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter8_slix103_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter8_slix119_fit_7_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter8_slix135_fit_9_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter8_slix152_fit_7_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter8_slix166_fit_7_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter8_slix182_fit_9_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter8_slix199_fit_9_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter8_slix214_fit_6_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter8_slix230_fit_10_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter8_slix247_fit_8_et0_chit.eps 
*    endif
*    if ([nc].eq.9) then
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter9_slix39_fit_10_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter9_slix56_fit_7_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter9_slix71_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter9_slix87_fit_8_et0_chit.eps  
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter9_slix104_fit_9_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter9_slix120_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter9_slix135_fit_7_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter9_slix150_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter9_slix166_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter9_slix182_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter9_slix198_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter9_slix215_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter9_slix232_fit_8_et0_chit.eps 
      n=[n]+1; file[n]=mhad2011-4_amplitude_vs_phi_counter9_slix247_fit_7_et0_chit.eps
*    endif
*
fname=mhad2011-4_amplitude_vs_phi_counters_slix_fit_8m.tex
if ($fexist([fname]).eq.1) then
  shell rm [fname]
endif
for/file 20 [fname] new
close 20
do i=1,[n]
fmess '\begin{figure}[ht!b]' [fname]
fmess '\begin{minipage}{\textwidth}' [fname]
fmess '{\centering\resizebox*{\textwidth}{!}' [fname]
txt={\includegraphics{[file[i]]}}\par}
fmess [txt] [fname]
ncnt=$sigma(int(([i]-1)/14)+1)
txt=\caption{Counter \No [ncnt], slise $sigma([i]-14*([ncnt]-1))}
fmess [txt] [fname]
fmess '\end{minipage}' [fname]
fmess '\end{figure}' [fname]
fmess '' [fname]
if ($sigma(mod([i],14)).eq.0) then
fmess '\clearpage' [fname]
fmess '' [fname]
endif
enddo
return

macro slixi nc=1 slix=1 pref=mhad2011-4 
    n=0
    if ([nc].eq.1) then
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter1_slx41_fit_8_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter1_slx56_fit_8_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter1_slx71_fit_8_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter1_slx87_fit_8_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter1_slx104_fit_9_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter1_slx121_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter1_slx136_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter1_slx152_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter1_slx168_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter1_slx184_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter1_slx200_fit_7_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter1_slx215_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter1_slx232_fit_9_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter1_slx248_fit_7_et0_chit.par
    endif
    if ([nc].eq.2) then
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter2_slx42_fit_8_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter2_slx58_fit_8_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter2_slx73_fit_7_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter2_slx89_fit_8_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter2_slx105_fit_7_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter2_slx120_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter2_slx137_fit_9_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter2_slx153_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter2_slx170_fit_9_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter2_slx185_fit_7_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter2_slx200_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter2_slx216_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter2_slx233_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter2_slx248_fit_7_et0_chit.par
    endif
    if ([nc].eq.3) then
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter3_slx40_fit_8_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter3_slx55_fit_7_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter3_slx70_fit_8_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter3_slx87_fit_9_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter3_slx104_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter3_slx120_fit_7_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter3_slx135_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter3_slx151_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter3_slx168_fit_9_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter3_slx183_fit_7_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter3_slx198_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter3_slx215_fit_9_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter3_slx231_fit_7_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter3_slx247_fit_8_et0_chit.par
    endif
    if ([nc].eq.4) then
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter4_slx41_fit_9_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter4_slx56_fit_6_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter4_slx70_fit_8_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter4_slx87_fit_8_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter4_slx104_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter4_slx120_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter4_slx136_fit_7_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter4_slx152_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter4_slx168_fit_9_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter4_slx184_fit_7_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter4_slx200_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter4_slx214_fit_7_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter4_slx231_fit_9_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter4_slx248_fit_7_et0_chit.par
    endif
    if ([nc].eq.5) then
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter5_slx42_fit_10_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter5_slx58_fit_6_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter5_slx73_fit_8_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter5_slx89_fit_8_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter5_slx105_fit_9_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter5_slx122_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter5_slx138_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter5_slx154_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter5_slx170_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter5_slx186_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter5_slx202_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter5_slx218_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter5_slx233_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter5_slx249_fit_8_et0_chit.par
    endif
    if ([nc].eq.6) then
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter6_slx40_fit_9_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter6_slx57_fit_8_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter6_slx73_fit_8_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter6_slx88_fit_7_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter6_slx103_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter6_slx120_fit_9_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter6_slx135_fit_7_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter6_slx151_fit_9_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter6_slx169_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter6_slx183_fit_7_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter6_slx198_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter6_slx215_fit_9_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter6_slx230_fit_7_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter6_slx245_fit_8_et0_chit.par
    endif
    if ([nc].eq.7) then
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter7_slx41_fit_9_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter7_slx57_fit_8_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter7_slx73_fit_8_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter7_slx89_fit_8_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter7_slx104_fit_7_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter7_slx119_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter7_slx134_fit_7_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter7_slx150_fit_9_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter7_slx166_fit_7_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter7_slx181_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter7_slx196_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter7_slx213_fit_9_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter7_slx231_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter7_slx246_fit_7_et0_chit.par
    endif
    if ([nc].eq.8) then
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter8_slx40_fit_9_et0_chit.par   
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter8_slx57_fit_8_et0_chit.par   
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter8_slx73_fit_8_et0_chit.par   
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter8_slx88_fit_7_et0_chit.par   
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter8_slx103_fit_8_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter8_slx119_fit_7_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter8_slx135_fit_9_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter8_slx152_fit_7_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter8_slx166_fit_7_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter8_slx182_fit_9_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter8_slx199_fit_9_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter8_slx214_fit_6_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter8_slx230_fit_10_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter8_slx247_fit_8_et0_chit.par 
    endif
    if ([nc].eq.9) then
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter9_slx39_fit_10_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter9_slx56_fit_7_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter9_slx71_fit_8_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter9_slx87_fit_8_et0_chit.par  
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter9_slx104_fit_9_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter9_slx120_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter9_slx135_fit_7_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter9_slx150_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter9_slx166_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter9_slx182_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter9_slx198_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter9_slx215_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter9_slx232_fit_8_et0_chit.par 
      n=[n]+1; file[n]=[pref]_amplitude_vs_phi_counter9_slx247_fit_7_et0_chit.par
    endif
    gl/cre fslix [file[slix]]
return



macro simshift
ve/cre  sse(9) r
ve/cre dsse(9) r
ve/cre  sss(9) r
ve/cre dsss(9) r
do i=1,9
*  fname=mhad2011_amplitude_vs_phi_counter[i]_slx145_fit_95.par
  fname=mhad2011-4_amplitude_vs_phi_counter[i]_slx145_fit_95_et0_chit.par
  ve/read pare,dpare [fname] '2f15.10'
*  fname=mhad2011_sim_amplitude_vs_phi_counter[i]_slx145_fit_95.par
  fname=mhad2011-4_sim0_amplitude_vs_phi_counter[i]_slx145_fit_95_et0_chit.par
  ve/read pars,dpars [fname] '2f15.10'
  ve/inp  sse([i])  pare(13)
  ve/inp dsse([i]) dpare(13)
  ve/inp  sss([i])  pars(13)
  ve/inp dsss([i]) dpars(13)
enddo  
ve/cre nc(9) r 1 2 3 4 5 6 7 8 9
ve/cre dnc(9) r
* exec $PER/s#vpl sss dsss nc dnc iatt=24 sz=0.1
* exec $PER/s#vpl sse dsse nc dnc iatt=20 sz=0.1 o=s
sigma sses = sse-sss
sigma dsses = sqrt(dsse**2+dsss**2)
* exec $PER/s#vpl sses dsses nc dnc iatt=20 sz=0.1 ll=-1
return

macro corrm nc=1
exec mapcal#expsimf [nc] -8 8
*fname=mhad2011_amplitude_vs_phi_counter[nc]_slx145_fit_95.par
fname=mhad2011-4_amplitude_vs_phi_counter[nc]_slx145_fit_95_et0_chit.par
ve/read pare,dpare [fname] '2f15.10'
*fname=mhad2011_sim_amplitude_vs_phi_counter[nc]_slx145_fit_95.par
fname=mhad2011-4_sim0_amplitude_vs_phi_counter[nc]_slx145_fit_95_et0_chit.par
ve/read pars,dpars [fname] '2f15.10'
dfe=$sigma(pare(13)-pare(12)-3*abs(pare(17)))
dfs=$sigma(pars(13)-pars(12)-3*abs(pars(17)))
l=$sigma(pare(13)-min([dfe],[dfs]))
dfe=$sigma(pare(15)-3*abs(pare(17))-pare(13))
dfs=$sigma(pars(15)-3*abs(pars(17))-pars(13))
r=$sigma(pare(13)+min([dfe],[dfs]))
exec hsigma @400 = @300*($300 gt [l] and $300 lt [r])
exec hsigma @400 = @400/@400
gl/cre f1 [l] 
gl/cre f5 [r]
l=$sigma(pare(13)-1.5/123*180/3.1415926535-3*abs(pare(18)))
r=$sigma(pare(13)+1.5/123*180/3.1415926535+3*abs(pare(18)))
gl/cre f2 [l]
gl/cre f3 [r]
gl/cre f4 $sigma(pare(14))
do i=1,5
  line [f[i]] 0 [f[i]] 100
enddo
read x
*
exec hsigma @500 = @300*($300 gt [l] and $300 lt [r])
exec hsigma @500 = abs(1-@500/@500)
exec hsigma @201 = @200*@400*@500
exec hsigma @301 = @300*@400*@500
hi/op/div 301 201 401 ! ! E
ve/cre p0(1) r 7
ve/cre dp0(1) r 
hi/fit 401 p0 ! 1 p0 ! ! ! dp0
*
mess [nc] [f1] [f2] [f3] [f4] [f5]
return

macro corrmip0
ve/cre  p0i(9) r
ve/cre dp0i(9) r
do i=1,9
  exec mapcal#corrm [i]
  ve/inp p0i([i]) $sigma(p0(1))
  ve/inp dp0i([i]) $sigma(dp0(1))
enddo
ve/cre nc(9) r 1 2 3 4 5 6 7 8 9
ve/cre dnc(9) r
set hcol 1
set pmci 4
* exec $PER/s#vpl p0i dp0i nc dnc sz=0.1 ll=1
ve/del p0ir,dp0ir,ncr,dncr
ve/copy p0i(1:7) p0ir
ve/copy dp0i(1:7) dp0ir
ve/copy nc(1:7) ncr
ve/copy dnc(1:7) dncr
ve/fit ncr p0ir dp0ir p0 s
return


macro corrmi
*
exec mapcal#simshift
ve/write sses phi_shift_sim.txt
*
fname=phi_cuts.txt
if ($fexist([fname]).eq.1) then
  shell rm [fname]
endif
for/file 20 [fname] new
close 20
do i=1,9
  exec mapcal#corrm [i]
  gl/imp f1
  gl/imp f2
  gl/imp f3
  gl/imp f4
  gl/imp f5
  txt=[i] [f1] [f2] [f3] [f4] [f5]
  fmess [txt] [fname]
enddo
return

macro empstabx i1=7 i2=17 dir=v2
do i=[i1],[i2]
  fname=empstab[i].kumac
  if ($fexist([fname]).eq.1) then
    shell rm [fname]
  endif
  for/file 20 [fname] new
  close 20
  txt=exec mapcal#empstab [dir] [i]
  fmess [txt] [fname]
  shell pawbigX11 -b [fname] 
*  exec mapcal#empstab [dir] [i]
enddo
return


macro empstab dir=v2 part=1
dir1=/work/users/konctbel/MinuitTest/
gl/cre tl [part]
nmax=20000
ve/cre runs([nmax]) r
ve/cre days([nmax]) r
ve/cre  eff(9,[nmax]) r
ve/cre deff(9,[nmax]) r
ve/cre  effa(9,[nmax]) r
ve/cre deffa(9,[nmax]) r
ve/cre  effb(9,[nmax]) r
ve/cre deffb(9,[nmax]) r
ve/cre  n0b(9,[nmax]) r
ve/cre  nxb(9,[nmax]) r
ve/cre  amp(9,[nmax]) r
ve/cre damp(9,[nmax]) r
ind=0
do nc=1,9
  hi/del 10[nc],11[nc],20[nc],21[nc],30[nc],31[nc]
enddo
if ([dir].eq.'v1') then
  suffh=.his
  sufft=.txt
else
  suffh=_[dir].his
  sufft=_[dir].txt
endif
do i=$sigma([part]*1000),$sigma([part]*1000+999)
  fname=[dir1][dir]/run_[i]_spects[suffh]
  if ($fexist([fname]).eq.1) then
    ind=[ind]+1
    exec mapcal#ndays [i]
    gl/imp ndays
    ve/inp days([ind]) [ndays]
    mess [fname]
    ve/inp runs([ind]) [i]
    hi/file 20 [fname]
    do nc=1,9
      hi/del 100,200
      idh=10+[nc]
      hi/copy //lun20/[idh] 100
      idhc=10[nc]
      if ($hexist([idhc]).eq.0) then
        hi/copy 100 [idhc]
      else
*        hi/op/add 100 [idhc] [idhc]
*        exec hsigma @[idhc] = @[idhc] + @100
      endif
      nx=$hinfo(100,'events')
      idh=[idh]+100000
      hi/copy //lun20/[idh] 200
      idhc=11[nc]
      if ($hexist([idhc]).eq.0) then
        hi/copy 200 [idhc]
      else
*        hi/op/add 200 [idhc] [idhc]
*        exec hsigma @[idhc] = @[idhc] + @200
      endif
      n0=$hinfo(200,'events')
      if ([nx].ne.0) then
        effi=[n0]/([nx]+[n0])
        deffi=$sigma(sqrt([effi]*(1-[effi])/[nx]))
      else
        effi=1
        deffi=1
      endif
      ve/inp eff([nc],[ind]) [effi]
      ve/inp deff([nc],[ind]) [deffi]
      ve/inp amp([nc],[ind]) $hinfo(100,'mean')
      dampi=$hinfo(100,'rms')
      dampi=$sigma([dampi]/sqrt([nx]))
      ve/inp damp([nc],[ind]) [dampi]
      if (([nx].lt.3).or.([dampi].lt.0.1)) then
        ve/inp damp([nc],[ind]) 10
      endif
*
      hi/del 100,200
      idh=20+[nc]
      hi/copy //lun20/[idh] 100
      idhc=20[nc]
      if ($hexist([idhc]).eq.0) then
        hi/copy 100 [idhc]
      else
*        hi/op/add 100 [idhc] [idhc]
*        exec hsigma @[idhc] = @[idhc] + @100
      endif
      nx=$hinfo(100,'events')
      idh=[idh]+100000
      hi/copy //lun20/[idh] 200
      idhc=21[nc]
      if ($hexist([idhc]).eq.0) then
        hi/copy 200 [idhc]
      else
*        hi/op/add 200 [idhc] [idhc]
*        exec hsigma @[idhc] = @[idhc] + @200
      endif
      n0=$hinfo(200,'events')
      if ([nx].ne.0) then
        effi=[n0]/[nx]
        deffi=$sigma(sqrt([effi]*(1-[effi])/[nx]))
      else
        effi=1
        deffi=1
      endif
      ve/inp effa([nc],[ind]) [effi]
      ve/inp deffa([nc],[ind]) [deffi]
*
      hi/del 100,200
      idh=30+[nc]
      hi/copy //lun20/[idh] 100
      idhc=30[nc]
      if ($hexist([idhc]).eq.0) then
        hi/copy 100 [idhc]
      else
*        hi/op/add 100 [idhc] [idhc]
*        exec hsigma @[idhc] = @[idhc] + @100
      endif
      nx=$hinfo(100,'events')
      idh=[idh]+100000
      hi/copy //lun20/[idh] 200
      idhc=31[nc]
      if ($hexist([idhc]).eq.0) then
        hi/copy 200 [idhc]
      else
*        hi/op/add 200 [idhc] [idhc]
*        exec hsigma @[idhc] = @[idhc] + @200
      endif
      n0=$hinfo(200,'events')
      if ([nx].ne.0) then
        effi=[n0]/[nx]
        deffi=$sigma(sqrt([effi]*(1-[effi])/[nx]))
      else
        effi=1
        deffi=1
      endif
      ve/inp effb([nc],[ind]) [effi]
      ve/inp deffb([nc],[ind]) [deffi]
      ve/inp n0b([nc],[ind]) [n0]
      ve/inp nxb([nc],[ind]) [nx]
*      
      ve/inp amp([nc],[ind]) $hinfo(100,'mean')
      dampi=$hinfo(100,'rms')
      dampi=$sigma([dampi]/sqrt([nx]))
      ve/inp damp([nc],[ind]) [dampi]
      if (([nx].lt.3).or.([dampi].lt.0.1)) then
        ve/inp damp([nc],[ind]) 3
      endif
    enddo
    close 20
  endif
enddo
ve/del tmp
ve/copy runs tmp
ve/cre runs([ind]) r
ve/cre druns([ind]) r
ve/copy tmp(1:[ind]) runs(1:[ind])
ve/del tmp
ve/copy days tmp
ve/cre days([ind]) r
ve/cre ddays([ind]) r
ve/copy tmp(1:[ind]) days(1:[ind])
do nc=1,9
  ve/cre  eff[nc]([ind]) r
  ve/cre deff[nc]([ind]) r
  ve/copy eff([nc],1:[ind]) eff[nc](1:[ind])
  ve/copy deff([nc],1:[ind]) deff[nc](1:[ind])
  ve/cre  amp[nc]([ind]) r
  ve/cre damp[nc]([ind]) r
  ve/copy amp([nc],1:[ind]) amp[nc](1:[ind])
  ve/copy damp([nc],1:[ind]) damp[nc](1:[ind])
  ve/write runs,days,amp[nc],damp[nc] [dir]/amplitude_vs_runs_counter[nc]_p[part].txt 4f15.6
*
  ve/cre  eff[nc]a([ind]) r
  ve/cre deff[nc]a([ind]) r
  ve/copy effa([nc],1:[ind]) eff[nc]a(1:[ind])
  ve/copy deffa([nc],1:[ind]) deff[nc]a(1:[ind])
*
  ve/cre  eff[nc]b([ind]) r
  ve/cre deff[nc]b([ind]) r
  ve/copy effb([nc],1:[ind]) eff[nc]b(1:[ind])
  ve/copy deffb([nc],1:[ind]) deff[nc]b(1:[ind])
*
  ve/cre n0[nc]b([ind]) r
  ve/cre nx[nc]b([ind]) r
  ve/copy n0b([nc],1:[ind]) n0[nc]b(1:[ind])
  ve/copy nxb([nc],1:[ind]) nx[nc]b(1:[ind])
  ve/write runs,days,n0[nc]b,nx[nc]b [dir]/n0_nx_vs_runs_counter[nc]_p[part].txt 4f15.6
enddo
* exec $PER/s#vpl eff1 deff1 runs druns sz=0.1
*exec vpl#pl eff1 deff1 runs druns sz=0.1
n=$vlen(days)
do i=2,$sigma([n]-1)
  dd=days([i])
  if ([dd].eq.0) then
    ii=[i]-1
    d1=days([ii])
    ii=[i]+1
    d2=days([ii])
    ve/inp days([i]) $sigma(([d2]+[d1])/2)
  endif
enddo
*
if ($fexist([dir1]amp1pe.txt).eq.1) then
  n=$vlen(runs)
  ve/cre a1pe(9,[n]) r
  ve/read a1pe amp1pe.txt
else
  exec mapcal#reada1
endif
return


macro a1per
if ($fexist(amp1pe.txt).eq.1) then
  n=$vlen(runs)
  ve/cre a1pe(9,[n]) r
  ve/read a1pe amp1pe.txt
else
  exec mapcal#reada1
endif
return


macro ndays run=10000 year0=2011 tl=x
*gl/imp tl
ve/del mlen
ve/cre mlen(12) r 31 28 31 30 31 30 31 31 30 31 30 31
  match=[run] & 20
  quote="
  quote="
  shell fgrep -m 1 -e [quote]$unquote([match])[quote] /work/users/konctbel/SepPar/runlist.tex0 > tmp[tl].txt
  ndig=$sigma(int(log10([run]))+1)
  ve/del year,month,day,hour,minute,second,beami
  nd=3+[ndig];  fmt=([nd]x,f4.0); ve/cre year(1)   r; ve/read year tmp[tl].txt [fmt]
  nd=8+[ndig];  fmt=([nd]x,f2.0); ve/cre month(1)  r; ve/read month tmp[tl].txt [fmt]
  nd=11+[ndig]; fmt=([nd]x,f2.0); ve/cre day(1)    r; ve/read day tmp[tl].txt [fmt]
  nd=14+[ndig]; fmt=([nd]x,f2.0); ve/cre hour(1)   r; ve/read hour tmp[tl].txt [fmt]
  nd=17+[ndig]; fmt=([nd]x,f2.0); ve/cre minute(1) r; ve/read minute tmp[tl].txt [fmt]
  nd=20+[ndig]; fmt=([nd]x,f2.0); ve/cre second(1) r; ve/read second tmp[tl].txt [fmt]
  nd=38+[ndig]; fmt=([nd]x,f10.3); ve/cre beami(1)  r; ve/read beami tmp[tl].txt [fmt]
  m=$sigma(month(1)-1)
  days=0
  do j=1,[m]
    days=$sigma([days]+mlen([j]))
  enddo
  days=$sigma([days]+day(1))
  dday=$sigma((hour(1)+(minute(1)+second(1)/60)/60)/24)
  nday=[days]+[dday]
  if ($sigma(year(1)).eq.2010) then
    nday=$sigma([nday]-365)
  endif
  if ($sigma(year(1)).eq.2012) then
    nday=$sigma([nday]+365)
  endif
  if ($sigma(year(1)).eq.2013) then
    nday=$sigma([nday]+365*2)
  endif
  shell /work/users/konctbel/MinuitTest/gtime [year0].01.01.00.00.00 $sigma(year(1)).$sigma(month(1)).$sigma(day(1)).$sigma(hour(1)).$sigma(minute(1)).$sigma(second(1)) gtime[tl].txt
  ve/read gtime gtime[tl].txt
  gl/cre ndays $sigma(gtime(1))
  gl/cre ebeam $sigma(beami(1))
*  mess [nday]
return


macro wlscorrz ncnt=1
*
shell fgrep -e Shifter mapcal2011-4shifter_v0.cal >& tmp.txt
shell $unquote('cat tmp.txt | sed "s/Shifter/ /g" > tmp1.txt')
ve/del ab,af,tau
ve/read ab,af,tau tmp1.txt
ve/cre wls(3) r $sigma(ab([ncnt])) $sigma(af([ncnt])) $sigma(tau([ncnt]))
*
fname=mhad2011-4_amplitude_vs_phi_counter[ncnt]_slx145_fit_95_et0_chit.par
ve/read pars,dpars [fname] '2f15.10'
ve/cre  cwls(14) r
ve/cre dcwls(14) r
ve/cre  zwls(14) r
ve/cre dzwls(14) r
do i=1,14
  exec mapcal#wlscorri [ncnt] [i]
  ve/inp cwls([i]) $sigma(p4(2)+p4(1))
  ve/inp dcwls([i]) $sigma(dp4(2))
  j=[i]+1
  ve/inp zwls([i]) $sigma((zic([i])+zic([j]))/2)
  read x
enddo
* exec $PER/s#vpl cwls dcwls zwls dzwls ll=-1
ve/del p3(3)
ve/copy wls p3
n1=1
n2=9
ve/del cwlsr,dcwlsr,zwlsr
ve/copy  cwls([n1]:[n2])  cwlsr
ve/copy dcwls([n1]:[n2]) dcwlsr
ve/copy  zwls([n1]:[n2])  zwlsr
ve/fit zwlsr cwlsr dcwlsr wlscorr.f s 3 p3
read x
set hcol 2
xmin=-10
xmax=10
*
a1=wls(1)
b1=wls(2)
c1=wls(3)
x=[xmin]
y1=$sigma(([a1])*exp(-([x])/([c1]))+([b1])*exp([x]/([c1])))
x=[xmax]
y2=$sigma(([a1])*exp(-([x])/([c1]))+([b1])*exp([x]/([c1])))
*
a2=p3(1)
b2=p3(2)
c2=p3(3)
x=[xmin]
y3=$sigma(([a2])*exp(-([x])/([c2]))+([b2])*exp([x]/([c2])))
x=[xmax]
y4=$sigma(([a2])*exp(-([x])/([c2]))+([b2])*exp([x]/([c2])))
*
ymax=$sigma(1.1*max(max([y1],[y2]),max([y3],[y4])))
null [xmin] [xmax] 0 [ymax]
set hcol 2
fun/pl ([a1])*exp(-x/([c1]))+([b1])*exp(x/([c1])) [xmin] [xmax] s
set hcol 4
fun/pl ([a2])*exp(-x/([c2]))+([b2])*exp(x/([c2])) [xmin] [xmax] s
set hcol 1
return


macro wlscorri ncnt=1 slx=0
exec mapcal#expsimf [ncnt] ! ! [slx]
nx=$hinfo(200,'xbins')
xmin=$hinfo(200,'xmin')
xmax=$hinfo(200,'xmax')
ve/cre  nexp([nx]) r
ve/cre dnexp([nx]) r
ve/cre  nsim([nx]) r
ve/cre dnsim([nx]) r
hi/get/cont 200 nexp
hi/get/err 200 dnexp
hi/get/cont 300 nsim
hi/get/err 300 dnsim
sigma res = nexp/nsim
sigma dres = res*sqrt((dnexp/nexp)**2+(dnsim/nsim)**2)
1d 400 ! [nx] [xmin] [xmax]
hi/put/cont 400 res
hi/put/err 400 dres
*
ve/del ic,if1,if2,if3,if4,if5
ve/read ic,if1,if2,if3,if4,if5 phi_cuts.txt
f1=if1([ncnt])
f5=if5([ncnt])
*
hi/pl 400([f1]:[f5])
*ve/cre p4(4) r 1 0.1 $sigma((if2([ncnt])+if3([ncnt]))/2) 0.5
ve/cre p4(4) r 1 0.1 $sigma(pars(13)) $sigma(pars(17))
ve/cre s4(4) r 0.01 0.001 0 0
ve/cre dp4(4) r
hi/fit 400([f1]:[f5]) p0+g sb 4 p4 s4 ! ! dp4
return


macro wlscorr ncnt=1 nr=10 nfe=12
*
exec mapcal#slxparn [ncnt] 16 mhad2011 _et0_chit.par
ve/del pno,dpno,xvio,dxvio
sigma  pno = order(pn,-pn)
sigma dpno = order(dpn,-pn)
sigma  xvio = order(xvi,-pn)
sigma dxvio = order(dxvi,-pn)
n=$vlen(pno)
ve/del pnr,dpnr,xvir,dxvir
ve/copy  pno(1:[n])  pnr
ve/copy dpno(1:[n]) dpnr
ve/copy  xvio(1:[n])  xvir
ve/copy dxvio(1:[n]) dxvir
i1=1
i2=[n]
sigma  pnr = order( pnr,xvir)
sigma dpnr = order(dpnr,xvir)
sigma dxvir = order(dxvir,xvir)
sigma  xvir = order(xvir,xvir)
ve/del pne,dpne,xvie,dxvie
ve/copy  pnr([i1]:[i2])  pne
ve/copy dpnr([i1]:[i2]) dpne
ve/copy  xvir([i1]:[i2])  xvie
ve/copy dxvir([i1]:[i2]) dxvie
i1=2
i2=[nfe]
ve/del pnef,dpnef,xvief,dxvief
ve/copy  pnr([i1]:[i2])  pnef
ve/copy dpnr([i1]:[i2]) dpnef
ve/copy  xvir([i1]:[i2])  xvief
ve/copy dxvir([i1]:[i2]) dxvief
*
exec mapcal#slxparn [ncnt] 16 mhad2011_sim _et0.par
ve/del pno,dpno,xvio,dxvio
sigma  pno = order(pn,-pn)
sigma dpno = order(dpn,-pn)
sigma  xvio = order(xvi,-pn)
sigma dxvio = order(dxvi,-pn)
n=$vlen(pno)
ve/del pnr,dpnr,xvir,dxvir
ve/copy  pno(1:[n])  pnr
ve/copy dpno(1:[n]) dpnr
ve/copy  xvio(1:[n])  xvir
ve/copy dxvio(1:[n]) dxvir
i1=1
i2=[n]
sigma  pnr = order( pnr,xvir)
sigma dpnr = order(dpnr,xvir)
sigma dxvir = order(dxvir,xvir)
sigma  xvir = order(xvir,xvir)
ve/del pns,dpns,xvis,dxvis
ve/copy  pnr([i1]:[i2])  pns
ve/copy dpnr([i1]:[i2]) dpns
ve/copy  xvir([i1]:[i2])  xvis
ve/copy dxvir([i1]:[i2]) dxvis
*
ve/del ai,bi 
ve/read ai,bi shifter.par
*
a=ai([ncnt])
b=bi([ncnt])
*
set pmci 4
* exec $PER/s#vpl pne dpne xvie dxvie sz=0.1 iatt=20
* exec $PER/s#vpl pns dpns xvis dxvis sz=0.1 iatt=24 o=s
fun/pl exp([a]+[b]*x) -15 15 s
*
ve/cre s5(5) r 0 1 1 1 1
ve/cre p5e(5) r 0 30 25 5 7
ve/cre dp5e(5) r
ve/fit xvief pnef dpnef gexp.f sb 5 p5e s5 ! ! dp5e
*
ve/cre p5s(5) r 0 30 25 5 7
ve/cre dp5s(5) r
ve/fit xvis pns dpns gexp.f sb 5 p5s s5 ! ! dp5s
*
ve/cre zr(1) r
ve/copy p5e par
fl=0
xi=0
dft=1
while ( [dft].gt.0 ) do
  xi=[xi]+0.1
  ve/inp zr(1) [xi]
  dft=$call('gexpp.f(zr,2)')
endwhile
ve/inp zr(1) [xi]
ft=$call('gexpp.f(zr,1)')
ns=0
while (($sigma(abs([ft]-[fl])).gt.0.0001).and.([ns].lt.100)) do
  ns=[ns]+1
  dft=$call('gexpp.f(zr,2)')
  dzr=$sigma(([fl]-[ft])/[dft])
  dzr=$sigma(min(abs([dzr]),0.1)*abs([dzr])/([dzr]))
  ve/inp zr(1) $sigma(zr(1)+([dzr]))
  ft=$call('gexpp.f(zr,1)')
*  mess [fl] [ft] $sigma(zr(1)) [dft] 
endwhile
xme=zr(1)
fme=$call('gexpp.f(zr,0)')
mess MaxE: [xme] [fme]
*
ve/copy p5s par
fl=0
xi=0
dft=1
while ( [dft].gt.0 ) do
  xi=[xi]+0.1
  ve/inp zr(1) [xi]
  dft=$call('gexpp.f(zr,2)')
endwhile
ve/inp zr(1) [xi]
ft=$call('gexpp.f(zr,1)')
ns=0
while (($sigma(abs([ft]-[fl])).gt.0.0001).and.([ns].lt.100)) do
  ns=[ns]+1
  dft=$call('gexpp.f(zr,2)')
  dzr=$sigma(([fl]-[ft])/[dft])
  dzr=$sigma(min(abs([dzr]),0.1)*abs([dzr])/([dzr]))
  ve/inp zr(1) $sigma(zr(1)+([dzr]))
  ft=$call('gexpp.f(zr,1)')
*  mess [fl] [ft] $sigma(zr(1)) [dft] 
endwhile
xms=zr(1)
fms=$call('gexpp.f(zr,0)')
mess MaxS: [xms] [fms]
*
ve/copy pne pnec
do i=1,$sigma([nr]-1)
  ve/inp pnec([i]) $sigma(pne([i])*exp([a]+[b]*xvie([i]))/pns([i]))
enddo
*
ve/cre zr(1) r
do i = [nr], 14
xi=xvie([i])
ve/copy p5e par
ve/inp zr(1) [xi]
fl=$call('gexpp.f(zr,0)')
fl=[fl]/[fme]*[fms]
ve/copy p5s par
ve/inp zr(1) $sigma([xms]+1)
ft=$call('gexpp.f(zr,0)')
ns=0
while (($sigma(abs([ft]-[fl])).gt.0.001).and.([ns].lt.1000)) do
  ns=[ns]+1
  dft=$call('gexpp.f(zr,1)')
  dzr=$sigma(([fl]-[ft])/[dft])
  dzr=$sigma(min(abs([dzr]),0.1)*abs([dzr])/([dzr]))
  ve/inp zr(1) $sigma(zr(1)+([dzr]))
  ft=$call('gexpp.f(zr,0)')
*  mess [fl] [ft] $sigma(zr(1)) [dft] [ns]
endwhile
mess '>>>'
fe=$sigma(exp([a]+[b]*zr(1))*[fme]/[fms])
mess $sigma(pne([i])) [fe] [fl] [ns]
ve/inp pnec([i]) [fe]
enddo
kk=0
nk=0
do i=5,10
  nk=[nk]+1
  kk=$sigma([kk]+exp([a]+[b]*xvis([i]))/pns([i]))
enddo
mess [kk] [nk]
kk=$sigma(1+2*(1-[kk]/[nk]))
mess [kk]
sigma pnec = pnec*[kk]
*
*read x
*
null -13 13 0 100
set pmci 2
* exec $PER/s#vpl pnec dpne xvie dxvie sz=0.1 iatt=20 o=s
*
set pmci 4
* exec $PER/s#vpl pne dpne xvie dxvie sz=0.1 iatt=20 o=s
* exec $PER/s#vpl pns dpns xvis dxvis sz=0.1 iatt=24 o=s
fun/pl exp([a]+[b]*x) -15 15 s
*
ve/cre p3(3) r 7.6580 20.842 15.719
ve/del xvier,pnecr,dpnecr
ve/copy xvie(2:14) xvier
ve/copy pnec(2:14) pnecr
ve/copy dpne(2:14) dpnecr
ve/fit xvier pnecr dpnecr dexp.f s 3 p3
*
return



macro p0
a=1.1
p0=0
do i=1,10
  p0=$sigma(exp(-[a]*(1-[p0])))
  mess [p0]
enddo
return


macro zminmax sm=0
ve/cre  zmin(27) r
ve/cre dzmin(27) r
ve/cre  zmax(27) r
ve/cre dzmax(27) r
i=0
do nc=1,9
  do nl=1,3
    i=[i]+1
    nh=$sigma([nl]*1000+100+[nc])
    fname=mhad2011-4_amplitude_vs_zr_counter[nh]_[sm]
    shell /work/snd2000/users/konctbel/SNDSimACC/readlifstat [fname].log
    ve/del pars0,dpars0
    ve/read pars0,dpars0 [fname].log.pars
    ve/copy  pars0(2)  zmin([i])
    ve/copy dpars0(2) dzmin([i])
    ve/copy  pars0(16)  zmax([i])
    ve/copy dpars0(16) dzmax([i])
  enddo
enddo
ve/cre  ni(27) r
ve/cre dni(27) r
sigma ni=array(27,1#27)
*
* exec $PER/s#vpl zmin dzmin ni dni sz=0.1 ll=-1
ve/cre p0(1) r -10
ve/cre dp0(1) r
ve/fit ni zmin dzmin p0 s 1 p0 ! ! ! dp0
zmin=$sigma(vsum(zmin)/27)
dzmin=$sigma(sqrt(vsum((zmin-[zmin])**2)/27))
mess zmin=[zmin] ([dzmin])
read x
* exec $PER/s#vpl zmax dzmax ni dni sz=0.1
ve/cre p0(1) r 10
ve/cre dp0(1) r
ve/fit ni zmax dzmax p0 s 1 p0 ! ! ! dp0
zmax=$sigma(vsum(zmax)/27)
dzmax=$sigma(sqrt(vsum((zmax-[zmax])**2)/27))
mess zmax=[zmax] ([dzmax])
return


macro reada1 fname=accled_amp1pe.txt
n=$vlen(runs)
ve/cre a1pe(9,[n]) r
do i=1,[n]
  mess [i]  $sigma(runs([i]))
  exec mapcal#runaccled amp1pe $sigma(runs([i]))
  ve/copy pars(1:9) a1pe(1:9,[i])
enddo
ve/write a1pe(1),a1pe(2),a1pe(3),a1pe(4),a1pe(5),a1pe(6),a1pe(7),a1pe(8),a1pe(9) [fname] 9f10.3
return

macro readpar vr=runs par=a v=0
fname=accled_[par]_v[v].txt
n=$vlen([vr])
ve/cre vx(9,[n]) r
do i=1,[n]
  mess [i]  $sigma([vr]([i]))
  exec mapcal#runaccled [par] $sigma([vr]([i]))
  ve/copy pars(1:9) vx(1:9,[i])
enddo
ve/write vx(1),vx(2),vx(3),vx(4),vx(5),vx(6),vx(7),vx(8),vx(9) [fname] 9f10.3
return

macro runaccled par=amp1pe run=9000
shell /work/users/konctbel/snd2k/R005-999/i386-SL5-opt/bin/clbarray accled [par] CURRENT [run] > out.txt
shell $unquote('fgrep -e "0000:" out.txt > out0.txt')
shell $unquote('cat out0.txt | sed "s/0000:/ /g" > out1.txt')
*exec cal#readpar
ve/del pars
ve/read pars out1.txt
return

macro readpds
n=$vlen(runs)
ve/cre pds(9,[n]) r
do i=1,[n]
  mess [i]  $sigma(runs([i]))
  exec mapcal#runaccpds pds $sigma(runs([i]))
  ve/copy pars(1:9) pds(1:9,[i])
enddo
ve/write runs,pds(1),pds(2),pds(3),pds(4),pds(5),pds(6),pds(7),pds(8),pds(9) runspds.txt 10f10.3
return

macro runaccpds par=pds run=9000
shell /work/users/konctbel/snd2k/R005-999/i386-SL5-opt/bin/clbarray accpds [par] CURRENT [run] > out.txt
shell $unquote('fgrep -e "0000:" out.txt > out0.txt')
shell $unquote('cat out0.txt | sed "s/0000:/ /g" > out1.txt')
*exec cal#readpar
ve/del pars
ve/read pars out1.txt
return

macro clblist v=amp2 dv=damp2 mode=0 nc=2
if ([mode].eq.0) then
*shell /work/users/konctbel/snd2k/R005-999/i386-SL5-opt/bin/clbixlist accled CURRENT > out.txt
ve/del nr1,nr2
ve/read nr1,nr2 accled.list
sigma nr2=order(nr2,nr1)
sigma nr1=order(nr1,nr1)
n=$vlen(nr2)
ve/inp nr2([n]) 1000000
ve/cre  xi(1000) r
ve/cre dxi(1000) r
ve/cre  yi(1000) r
ve/cre dyi(1000) r
ve/cre if1(1000) r
ve/cre if2(1000) r
ind=1
ind0=0
nm=0
ym=0
dym=0
rmin=runs(1)
rmax=runs(1)
i1=1
i2=1
do i=1,$vlen([v])
  r=runs([i])
  r1=nr1([ind])
  r2=nr2([ind])
  if (([r].ge.[r1]).and.([r].le.[r2])) then
    nm=[nm]+1
    rmax=[r]
    i2=[i]
  else
    if (([i].gt.1).and.([nm].ne.0)) then
      ind0=[ind0]+1
*      ve/inp  yi([ind0]) $sigma([ym]/[nm])
*      ve/inp dyi([ind0]) $sigma(sqrt(vsum([dym]))/[nm])
      ve/del [v]r,[dv]r
      ve/copy [v]([i1]:[i2]) [v]r
      ve/copy [dv]([i1]:[i2]) [dv]r
      ve/inp  yi([ind0]) $sigma(vsum([v]r/[dv]r**2)/vsum(1/[dv]r**2))
      ve/inp dyi([ind0]) $sigma(1/sqrt(vsum(1/[dv]r**2)))
      ve/inp  xi([ind0]) $sigma((days([i2])+days([i1]))/2)
      ve/inp dxi([ind0]) $sigma((days([i2])-days([i1]))/2)
      ve/inp if1([ind0]) [i1]
      ve/inp if2([ind0]) [i2]
      mess [i1] [i2] $sigma(runs([i1])) $sigma(runs([i2]))
      xi0=xi([ind0])
      if ([xi0].eq.0) then
        mess [i1] [i2] $sigma(runs([i1])) $sigma(runs([i2])) $sigma(days([i2])) $sigma(days([i1]))
      endif
    endif
    while ([r].gt.$sigma(nr2([ind]))) do
      ind=[ind]+1
    endwhile
*    mess $sigma(nr1([ind])) [r] $sigma(nr2([ind]))
    nm=0
    ym=0
    dym=0
    r1=nr1([ind])
    r2=nr2([ind])
    if (([r].ge.[r1]).and.([r].le.[r2])) then
      nm=[nm]+1
      rmin=[r]
      rmax=[r]
      i1=[i]
      i2=[i]
    endif
  endif
  r1=nr1([ind])
  r2=nr2([ind])
  if (([r].ge.[r1]).and.([r].le.[r2])) then
    ym=$sigma([ym]+[v]([i]))
    dym=$sigma([dym]+[dv]([i])**2)
  else
    mess run = [r]
  endif
enddo
ind0=[ind0]+1
ve/inp  yi([ind0]) $sigma([ym]/[nm])
ve/inp dyi([ind0]) $sigma(sqrt(vsum([dym]))/[nm])
ve/inp  xi([ind0]) $sigma((days([i2])+days([i1]))/2)
ve/inp dxi([ind0]) $sigma((days([i2])-days([i1]))/2)
ve/inp if1([ind0]) [i1]
ve/inp if2([ind0]) [i2]

exec mapcal#vecut xi
exec mapcal#vecut dxi
exec mapcal#vecut yi
exec mapcal#vecut dyi
exec mapcal#vecut if1
exec mapcal#vecut if2

endif

if ([mode].eq.1) then
n=$vlen(if1)
ve/cre  xi([n]) r
ve/cre dxi([n]) r
ve/cre  yi([n]) r
ve/cre dyi([n]) r
ve/cre  a1([n]) r
ve/cre da1([n]) r
do i=1,[n]
  i1=if1([i])
  i2=if2([i])
  di=[i2]-[i1]+1
  ve/del [v]r,[dv]r
  ve/copy [v]([i1]:[i2]) [v]r
  ve/copy [dv]([i1]:[i2]) [dv]r
*  ve/inp yi([i]) $sigma(vsum([v]r)/[di])
*  ve/inp dyi([i]) $sigma(sqrt(vsum([dv]r**2))/[di])
  ve/inp  xi([i]) $sigma((days([i2])+days([i1]))/2)
  ve/inp dxi([i]) $sigma((days([i2])-days([i1]))/2)
  ve/inp  yi([i]) $sigma(vsum([v]r/[dv]r**2)/vsum(1/[dv]r**2))
  ve/inp dyi([i]) $sigma(1/sqrt(vsum(1/[dv]r**2)))
  if ($substring([v],1,1).eq.'n') then
    ve/inp  yi([i]) $sigma(vsum([v]r)/vsum([dv]r))
    ve/inp dyi([i]) $sigma(sqrt(yi([i])*(1-yi([i]))/vsum([dv]r)))
  endif
  if ($vexist(a1pe).eq.1) then
    a1i1=a1pe([nc],[i1])
    a1i2=a1pe([nc],[i2])
    ve/inp a1([i]) [a1i1]
    if ([a1i1].ne.[a1i2]) then
      mess '>>>>>' [i] [i1] [i2] [a1i1] [a1i2]
    endif
  endif
enddo
endif

** exec $PER/s#vpl yi dyi xi dxi sz=0.1
exec vpl#pl yi dyi xi dxi sz=0.1
ve/fit xi yi dyi e s
ve/cre p3(3) r 0 10 4000
set plci 1
ve/fit xi yi dyi expp0.f s 3 p3
ve/del dyif
ve/copy dyi dyif
ve/cre xic(1) r
do i=1,$vlen(xi)
  ve/copy xi([i]) xic(1)
  yic=$call('expp0p.f(xic)')
  if ($sigma(abs([yic]-yi([i]))).gt.1) then
    ve/inp dyif([i]) 1000
  endif
enddo

set plci 4
ve/fit xi yi dyif expp0.f s 3 p3

*
null 0 700 0 0.02
** exec $PER/s#vpl yi dyi xi dxi sz=0.1 o=s
exec vpl#pl yi dyi xi dxi sz=0.1 o=s
*ve/fit xi yi dyi expp0e.f s 3 p3e
*
i1=1
i2=52
ve/del xip,dxip,yip,dyip
ve/copy  yi([i1]:[i2]) yip
ve/copy dyi([i1]:[i2]) dyip
ve/copy  xi([i1]:[i2]) xip
ve/copy dxi([i1]:[i2]) dxip
ve/fit xip yip dyip expp0e.f s 3 p3e
return

macro vecut v=xi n=0
ve/del tmp
ve/copy [v] tmp
if ([n].eq.0) then
*  n=$vlen([v])
  n=$sigma(vsum([v]/[v]))
endif
ve/del [v]
ve/cre [v]([n]) r
nc=$vlen(tmp)
if ([nc].ne.0) then
  nc=$sigma(min([n],[nc]))
else
  nc=[n]
endif
if ([nc].ne.0) then
  ve/copy tmp(1:[nc]) [v]
endif
return



macro pardayc nc=1 dir=v2
*
ve/del id1,id2
ve/read id1,id2 daycutsx.txt 2f10.2
*
ve/del pei,dpei,effi,deffi,mueic,dmueic,r,dr,xei,dxei,kei
ve/read pei,dpei,effi,deffi,mueic,dmueic,r,dr,xei,dxei,kei [dir]/parday_counter[nc]_[dir].txt 11e15.6
*goto 1
exec mapcal#parday n0[nc]b nx[nc]b [nc]
ve/del effi,deffi
ve/copy pei effi
ve/copy dpei deffi
sigma mueic =-log(pei)
sigma dmueic = dpei/pei
exec mapcal#parday amp[nc] damp[nc] [nc]
set pmci 4
** exec $PER/s#vpl pei dpei xei dxei sz=0.05 
exec vpl#pl pei dpei xei dxei sz=0.05 
set pmci 2
** exec $PER/s#vpl mueic dmueic xei dxei sz=0.05 o=s
exec vpl#pl mueic dmueic xei dxei sz=0.05 o=s
ve/del y1,dy1,y2,dy2
ve/copy  pei  y1
ve/copy dpei dy1
ve/copy  effi  y2
ve/copy deffi dy2
sigma r = sqrt(2*(y1+log(y2))/y1**2)
*  sigma dr = sqrt((dy2/y2/y1)**2+(log(y2)*dy1/y1**2)**2)
sigma dr = sqrt((dy1*(1+2*log(y2)/y1))**2+(dy2/y2)**2)/r/y1**2
*
1:
ve/cre dp3(3) r
ve/cre dp3a(3) r
ve/cre p3(3) r 3 3 1000
do i=1,5
  ve/fit xei mueic dmueic expp0.f s 3 p3
enddo
np=$vlen(xei)
ve/cre kei([np]) r [np]*1
do i=1,$vlen(xei)
  fi=$sigma(p3(1)+p3(2)*exp(-xei([i])/p3(3)))
  if ($sigma(abs([fi]-mueic([i]))/dmueic([i])).gt.3)  then
    ve/inp kei([i]) 10000
  endif
enddo
sigma ek = dmueic*kei
do i=1,5
  ve/fit xei mueic ek expp0.f s 3 p3
enddo
ve/del p3a
ve/copy p3 p3a
ve/cre s3(3) r 0.1 0.1 0
sigma ek = dpei*kei
do i=1,5
  ve/fit xei pei ek expp0.f sb 3 p3a s3 
enddo
ve/del keia
ve/copy kei keia
do i=1,$vlen(xei)
  fi=$sigma(p3a(1)+p3a(2)*exp(-xei([i])/p3a(3)))
  if ($sigma(abs([fi]-pei([i]))).gt.1)  then
    ve/inp keia([i]) 10000
  endif
enddo
sigma ek = dpei*keia
do i=1,5
  ve/fit xei pei ek expp0.f sb 3 p3a s3 
enddo
*
set lwid 1
set fwid 1
set hwid 1
null -40 720 0 12
set pmci 2
** exec $PER/s#vpl mueic dmueic xei dxei sz=0.05 o=s
exec vpl#pl mueic dmueic xei dxei sz=0.05 o=s
set pmci 4
** exec $PER/s#vpl pei dpei xei dxei sz=0.05 o=s
exec vpl#pl pei dpei xei dxei sz=0.05 o=s
set plci 2
sigma ek = dmueic*kei
ve/fit xei mueic ek expp0.f s 3 p3 ! ! ! dp3
sigma ek = dpei*keia
ve/fit xei pei ek expp0.f sb 3 p3a s3 ! ! sp3a
*
txt=[t]?[nc]! = $sigma(int(p3(3)*10+0.5)/10) [\261] $sigma(int(dp3(3)*10+0.5)/10)
*exec $PER/s#tf 0.6 0.1 [txt]
*
set plci 4
ve/fit xei pei ek expp0.f s 3 p3a
txt=A?[nc]!, pe
atitle 'T, days' [txt]
*
set basl 0.01
set ltyp 15
do i=1,$vlen(id1)
  exec mapcal#ndays $sigma(id1([i]))
  gl/imp ndays
  line [ndays] $GRAFINFO('WNYMIN') [ndays] $GRAFINFO('WNYMAX')
enddo
set ltyp 1
*
exec save [dir]/amplitude_vs_days_counter[nc]_[dir].eps f
*
*
set lwid 1
set fwid 1
set hwid 1
null -40 720 0 12
set pmci 2
** exec $PER/s#vpl mueic dmueic xei dxei sz=0.05 o=s
exec vpl#pl mueic dmueic xei dxei sz=0.05 o=s
set pmci 4
** exec $PER/s#vpl pei dpei xei dxei sz=0.05 o=s
exec vpl#pl pei dpei xei dxei sz=0.05 o=s
set plci 2
*
ve/del p14
ve/read p14 [dir]/amplitude_vs_times_fourerf_counter[nc]_[dir].par
ve/del p12
ve/copy p14(1:12) p12
fun/pl fourerfp.f -100 800 s
fun/pl fourerfpa.f -100 800 s
*
txt=A?[nc]!, pe
atitle 'T, days' [txt]
*
set basl 0.01
set ltyp 15
do i=1,$vlen(id1)
  exec mapcal#ndays $sigma(id1([i]))
  gl/imp ndays
  line [ndays] $GRAFINFO('WNYMIN') [ndays] $GRAFINFO('WNYMAX')
enddo
set ltyp 1
*
exec save [dir]/amplitude_vs_days_fourerf_counter[nc]_[dir].eps f
*
ve/cre p2(2) r 0.25 0
ve/del keia
ve/copy kei keia
do j=1,3
a=p2(1)
b=p2(2)
sigma fi1 = abs(r - ([a]+([b])*xei))
sigma fi2 = fi1/dr
do i=1,$vlen(xei)
  if (($sigma(fi1([i])).gt.0.15).or.($sigma(fi2([i])).gt.5).or.($sigma(abs(dr([i]))).lt.0.001))  then
    ve/inp keia([i]) 10000
  else
    ve/inp keia([i]) 1
  endif
enddo
sigma ek = dpei*keia
do i=1,1
  ve/fit xei r ek p1 sb 2 p2 
enddo
enddo
sigma rk = r*keia
*
null -40 720 0 1
set pmci 4
** exec $PER/s#vpl r dr xei dxei sz=0.05 o=s
sigma r = min(r,1)
sigma dr = min(dr,1)
exec vpl#pl r dr xei dxei sz=0.05 o=s
ve/fit xei r ek p1 s 2 p2
ve/fit xei r ek p0 s 2 p2
ve/cre p4(4) r 0.2 0.1 250 300
ve/fit xei r ek erf0.f s 4 p4
ve/cre p3(3) r 0.3 -0.1 300
do i=1,5
  ve/fit xei r ek expp0a.f s 3 p3
enddo
set pmci 2
exec vpl#pl rk dr xei dxei sz=0.05 o=s
txt=[s]?A,[nc]! / A?[nc]!
atitle 'T, days' [txt]
*
set basl 0.01
set ltyp 15
do i=1,$vlen(id1)
  exec mapcal#ndays $sigma(id1([i]))
  gl/imp ndays
  line [ndays] $GRAFINFO('WNYMIN') [ndays] $GRAFINFO('WNYMAX')
enddo
set ltyp 1
*
exec save [dir]/uniformity_vs_days_counter[nc]_[dir].eps f
*
ve/write pei,dpei,effi,deffi,mueic,dmueic,r,dr,xei,dxei,kei [dir]/parday_counter[nc]_[dir].txt 11e15.6
return



macro pardaycn nc=1 dir=v2
exec mapcal#ampeff n0[nc]b nx[nc]b [nc]
ve/del effi,deffi
ve/copy yi effi
ve/copy dyi deffi
sigma mueic =-log(yi)
sigma dmueic = dyi/yi
exec mapcal#ampeff amp[nc] damp[nc] [nc]
set pmci 4
** exec $PER/s#vpl pei dpei xei dxei sz=0.05 
exec vpl#pl yi dyi xi dxi sz=0.05 
set pmci 2
** exec $PER/s#vpl mueic dmueic xei dxei sz=0.05 o=s
exec vpl#pl mueic dmueic xi dxi sz=0.05 o=s
ve/del y1,dy1,y2,dy2,pei,dpei,xei,dxei
ve/copy xi xei
ve/copy dxi dxei
ve/copy yi pei
ve/copy dyi dpei
ve/copy  pei  y1
ve/copy dpei dy1
ve/copy  effi  y2
ve/copy deffi dy2
sigma r = sqrt(2*(y1+log(y2))/y1**2)
*  sigma dr = sqrt((dy2/y2/y1)**2+(log(y2)*dy1/y1**2)**2)
sigma dr = sqrt((dy1*(1+2*log(y2)/y1))**2+(dy2/y2)**2)/r/y1**2
*
1:
*
ve/cre p2(2) r 0.22 0.0002
sigma ek = dr
  iek=0; ve/cre ivek(1) r; ve/del vek
  iek=[iek]+1; ve/inp ivek(1)  1; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1)  2; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1)  3; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1)  4; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1) 21; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1) 22; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1) 23; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1) 38; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1) 39; exec vappend vek ivek
  do i=1,[iek]
    j=vek([i])
    ve/inp ek([j]) 1000
  enddo
ve/fit xei r ek fp1.f ! 2 p2
np=$vlen(xei)
ve/cre keia([np]) r [np]*1
  if ([nc].eq.0) then
    drm=0.001
  else
    drm=0.005
  endif
do j=1,3
  a=p2(1)
  b=$sigma(abs(p2(2)))
  sigma fi1 = abs(r - ([a]+[b]*xei))
  sigma fi2 = fi1/dr
  n=$vlen(xei)
  do i=1,[n]
    if ((($sigma(fi1([i])).gt.0.05).or.($sigma(fi2([i])).gt.3)).or.($sigma(abs(dr([i]))).lt.[drm]))  then
      ve/inp keia([i]) 10000
    else
      ve/inp keia([i]) 1
    endif
  enddo
  ve/inp keia([n]) 10000
  n=[n]-1
  ve/inp keia([n]) 10000
  sigma ek = dr*keia
  do i=1,1
    ve/fit xei r ek fp1.f sb 2 p2 
  enddo
enddo
sigma rk = r*keia
*
  sigma ek = dr*keia
  ve/fit xei r ek p1 sb 2 p2 
*
null -40 720 0 0.4
set pmci 4
** exec $PER/s#vpl r dr xei dxei sz=0.05 o=s
sigma r = min(r,0.4)
sigma dr = min(dr,0.4)
exec vpl#pl r dr xei dxei sz=0.05 o=s
do i=1,5
  ve/fit xei r ek fp1.f s 2 p2
enddo
ve/copy p2 par2
ve/copy keia keiar
*ve/fit xei r ek p0 s 2 p2
*ve/fit xei r ek rt.f s 5 par
*ve/cre p4(4) r 0.2 0.1 250 300
*ve/fit xei r ek erf0.f s 4 p4
*ve/cre p3(3) r 0.3 -0.1 300
do i=1,5
*  ve/fit xei r ek expp0a.f s 3 p3
enddo
*set pmci 2
*exec vpl#pl rk dr xei dxei sz=0.05 o=s
txt=[s]?A,[nc]! / A?[nc]!
atitle 'T, days' [txt]
*
*
if ($vlen(id1).lt.10) then
set basl 0.01
set ltyp 15
do i=1,$vlen(id1)
  exec mapcal#ndays $sigma(id1([i]))
  gl/imp ndays
  line [ndays] $GRAFINFO('WNYMIN') [ndays] $GRAFINFO('WNYMAX')
enddo
endif
set ltyp 1
*
exec save [dir]/uniformity_vs_times_counter[nc]_[dir]_p$vlen(ti).eps f
*read x
if ([nc].eq.0) then
  sigma ek = dpei
  iek=0; ve/cre ivek(1) r; ve/del vek
  iek=[iek]+1; ve/inp ivek(1)  1; exec vappend vek ivek
*  iek=[iek]+1; ve/inp ivek(1)  2; exec vappend vek ivek
*  iek=[iek]+1; ve/inp ivek(1)  3; exec vappend vek ivek
*  iek=[iek]+1; ve/inp ivek(1)  4; exec vappend vek ivek
*  iek=[iek]+1; ve/inp ivek(1) 21; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1) 22; exec vappend vek ivek
*  iek=[iek]+1; ve/inp ivek(1) 23; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1) 39; exec vappend vek ivek
  do i=1,[iek]
    j=vek([i])
    ve/inp ek([j]) 1000
  enddo
  ve/cre par(5) r 0.233 0 0.04 0 100
endif  
*if (([nc].ge.1) then
  sigma ek = dpei
  iek=0; ve/cre ivek(1) r; ve/del vek
  iek=[iek]+1; ve/inp ivek(1)  1; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1)  2; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1)  3; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1)  4; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1) 21; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1) 22; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1) 23; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1) 39; exec vappend vek ivek
  do i=1,[iek]
    j=vek([i])
    ve/inp ek([j]) 1000
  enddo
  ve/cre par(5) r 0.233 0 0.04 0 100
*endif  
if ([nc].eq.1) then
  sigma ek = dpei
  iek=0; ve/cre ivek(1) r; ve/del vek
  iek=[iek]+1; ve/inp ivek(1)  1; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1)  2; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1)  3; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1)  4; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1) 21; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1) 22; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1) 23; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1) 39; exec vappend vek ivek
  do i=1,[iek]
    j=vek([i])
    ve/inp ek([j]) 1000
  enddo
  ve/cre par(5) r 0.233 0 0.04 0 100
endif  
if ([nc].eq.2) then
  sigma ek = dpei
  iek=0; ve/cre ivek(1) r; ve/del vek
  iek=[iek]+1; ve/inp ivek(1)  1; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1)  2; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1)  3; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1)  4; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1) 21; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1) 22; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1) 23; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1) 39; exec vappend vek ivek
  do i=1,[iek]
    j=vek([i])
    ve/inp ek([j]) 1000
  enddo
  ve/cre par(5) r 0.233 0 0.04 0 100
endif  
if ([nc].eq.3) then
  sigma ek = dpei
  iek=0; ve/cre ivek(1) r; ve/del vek
  iek=[iek]+1; ve/inp ivek(1)  1; exec vappend vek ivek
*  iek=[iek]+1; ve/inp ivek(1)  2; exec vappend vek ivek
*  iek=[iek]+1; ve/inp ivek(1)  3; exec vappend vek ivek
*  iek=[iek]+1; ve/inp ivek(1)  4; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1)  6; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1) 21; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1) 22; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1) 23; exec vappend vek ivek
*  iek=[iek]+1; ve/inp ivek(1) 38; exec vappend vek ivek
  iek=[iek]+1; ve/inp ivek(1) 39; exec vappend vek ivek
  do i=1,[iek]
    j=vek([i])
    ve/inp ek([j]) 1000
  enddo
  ve/cre par(5) r 0.233 0 0.04 0 100
endif  
*
set pmci 2
exec vpl#pl pei dpei xei dxei sz=0.05
set pmci 4
exec vpl#pl mueic dmueic xei dxei sz=0.05 o=s
ve/cre dp3(3) r
ve/cre dp3a(3) r
ve/cre p3(3) r 3 3 1000
do i=1,5
  ve/fit xei mueic dmueic expp0.f s 3 p3
enddo
np=$vlen(xei)
ve/cre kei([np]) r [np]*1
do i=1,$vlen(xei)
  fi=$sigma(p3(1)+p3(2)*exp(-xei([i])/p3(3)))
  if ($sigma(abs([fi]-mueic([i]))/dmueic([i])).gt.3)  then
    ve/inp kei([i]) 10000
  endif
enddo
ve/inp kei([np]) 10000
sigma ek = dmueic*kei
do i=1,5
  ve/fit xei mueic ek expp0.f s 3 p3
enddo
ve/del p3a
ve/copy p3 p3a
ve/cre s3(3) r 0.1 0.1 0
sigma ek = dpei*kei
ve/cre ek([np]) r [np]*0.2
do j=1,3
ve/del keia
ve/copy kei keia
do i=1,1
  ve/fit xei pei ek expp0.f sb 3 p3a s3 
enddo
do i=1,$vlen(xei)
  fi=$sigma(abs(p3a(1))+abs(p3a(2))*exp(-xei([i])/p3a(3)))
  if ($sigma(abs([fi]-pei([i]))).gt.0.5)  then
    ve/inp keia([i]) 10000
  endif
enddo
sigma ek = ek*keia
enddo
sigma ek = dpei*keia
do i=1,1
  ve/fit xei pei ek expp0.f sb 3 p3a s3 
enddo
*
set lwid 1
set fwid 1
set hwid 1
null -40 720 0 12
set pmci 2
** exec $PER/s#vpl mueic dmueic xei dxei sz=0.05 o=s
exec vpl#pl mueic dmueic xei dxei sz=0.05 o=s
set pmci 4
** exec $PER/s#vpl pei dpei xei dxei sz=0.05 o=s
exec vpl#pl pei dpei xei dxei sz=0.05 o=s
set plci 2
sigma ek = dmueic*kei
ve/fit xei mueic ek expp0.f s 3 p3 ! ! ! dp3
sigma ek = dpei*keia
ve/fit xei pei ek expp0.f sb 3 p3a s3 ! ! sp3a
if (([nc].ge.0).and.([nc].le.9)) then
  ve/cre s7(7) r 0.1 0.1 0.1 10 1 5 0
  ve/cre p7(7) r 0 6.4147 3.7943 -2281.5 11.770 460.51 0
  sigma ek = dmueic
  do i=1,[iek]
    j=vek([i])
    ve/inp ek([j]) 1000
  enddo
*  ve/fit xei mueic ek gexp0.f sb 7 p7 s7
  sigma ek = dpei
  do i=1,[iek]
    j=vek([i])
    ve/inp ek([j]) 1000
  enddo
*  ve/cre s7(7) r 6*0 0.01
  set plci 3
  do i=1,3
    ve/fit xei pei ek gexp0.f sb0 7 p7 s7
  enddo
*  ve/fit xei pei ek gexp0.f sb 7 p7 s7
*  ve/write p7 [dir]/amplitude_vs_times_counter[nc]_[dir].pars
  sigma ek = dmueic
  do i=1,[iek]
    j=vek([i])
    ve/inp ek([j]) 1000
  enddo
  ve/copy s7 s7t
  ve/inp s7(4) 0
  ve/inp s7(5) 0
  ve/inp s7(6) 0
  ve/fit xei mueic ek gexp0.f sb0 7 p7 s7
  s7(4)=s7t(4)
  do i=1,3
    ve/fit xei mueic ek gexp0.f sb0 7 p7 s7
  enddo
  ve/fit xei mueic ek gexp0.f sb 7 p7 s7
*
  sigma ek = dpei
  do i=1,[iek]
    j=vek([i])
    ve/inp ek([j]) 1000
  enddo
  ve/inp s7(4) 0
  s7(5)=s7t(5)
  s7(6)=s7t(6)
  set plci 3
  do i=1,3
    ve/fit xei pei ek gexp0.f sb0 7 p7 s7
  enddo
  ve/fit xei pei ek gexp0.f sb 7 p7 s7
  ve/write p7 [dir]/amplitude_vs_times_counter[nc]_[dir].pars
  set plci 6
*  fun/pl gexp0s.f -90 800 s
  ve/cre p10(10) r 7.317316 0.2921743 -1.254233 0.3240522 0.3240522 0.1741663 82.14809 458.2414 1.929329 142.6761 4.625659
  ve/del p10
  ve/read p10 [dir]/fourerf.par
  ve/inp p10(5) 0
  ve/del p12,p14
  ve/copy p10 p12
  exec vappend p12 par2
  ve/copy p12 p14
  exec vappend p14 par2
  ve/cre s12(12) r 5*0.01
  ks=3
  ve/cre p12min(12) r 5*-100 0 0 $sigma(p12(8)/[ks]) $sigma(p12(9)/[ks]) $sigma(p12(10)/[ks]) 
  ve/cre p12max(12) r 5*100 1000 10000 $sigma(p12(8)*[ks]) $sigma(p12(9)*[ks]) $sigma(p12(10)*[ks]) 
  sigma ek = dmueic
  do i=1,[iek]
    j=vek([i])
    ve/inp ek([j]) 1000
  enddo
  do i=1,5
    ve/fit xei mueic ek fourerf.f sb 12 p12 s12 p12min p12max
  enddo
  np=$vlen(xei)
  ve/cre keia([np]) r [np]*1
  ve/cre xi(1) r
  do j=1,3
    do i=1,[np]
      ve/inp xi(1) xei([i])
      dfi1=$call('fourerfp.f(xi)')
      dfi1=$sigma(abs(mueic([i])-[dfi1]))
      dmui=$sigma(abs(dmueic([i])))
      dfi2=$sigma([dfi1]/[dmui])
      if (([dfi1].gt.1).or.([dfi2].gt.3).or.([dmui].lt.0.01))  then
        ve/inp keia([i]) 10000
      else
        ve/inp keia([i]) 1
      endif
    enddo
    ve/inp keia(1) 10000
    ve/inp keia(2) 10000
    ve/inp keia(3) 10000
    ve/inp keia(4) 10000
    ve/inp keia([np]) 10000
    sigma ek = dmueic*keia
    do i=1,2
      ve/fit xei mueic ek fourerf.f sb 12 p12 s12 p12min p12max
    enddo
  enddo
  sigma ek = dmueic*keia
  sigma eka = dpei*keiar
  ve/cre s14(12) r 12*0 0.01 1e-6
  do i=1,5
    n=$vlen(xei)
    ve/cre keiar([n]) r [n]*1
    mean=0
    rms=0
    do j=1,[n]
      ve/inp xi(1) xei([j])
      dfi1=$call('fourerfpa.f(xi)')
      dfi1=$sigma(pei([j])-([dfi1]))
      mean=[mean]+([dfi1])
      rms=$sigma([rms]+([dfi1])**2)
    enddo
    mean=[mean]/[n]
    rms=[rms]/[n]
    rms0=$sigma(sqrt([rms]-([mean])**2))
*    
    mean=0
    rms=0
    nevt=0
    do j=1,[n]
      ve/inp xi(1) xei([j])
      dfi1=$call('fourerfpa.f(xi)')
      dfi1=$sigma(pei([j])-([dfi1]))
      if ($sigma(abs([dfi1])/2).lt.[rms0]) then
        mean=[mean]+([dfi1])
        rms=$sigma([rms]+([dfi1])**2)
        nevt=[nevt]+1
      endif
    enddo
    mean=[mean]/[nevt]
    rms=[rms]/[nevt]
    rms0=$sigma(sqrt([rms]-([mean])**2))
*    
    mean=0
    rms=0
    nevt=0
    do j=1,[n]
      ve/inp xi(1) xei([j])
      dfi1=$call('fourerfpa.f(xi)')
      dfi1=$sigma(pei([j])-([dfi1]))
      if ($sigma(abs([dfi1])/2).lt.[rms0]) then
        mean=[mean]+([dfi1])
        rms=$sigma([rms]+([dfi1])**2)
        nevt=[nevt]+1
      endif
    enddo
    mean=[mean]/[nevt]
    rms=[rms]/[nevt]
    rms=$sigma(min(sqrt([rms]-([mean])**2),1))
    rms3=[rms]*3
    do j=1,[n]
      ve/inp xi(1) xei([j])
      dfi1=$call('fourerfpa.f(xi)')
      dfi1=$sigma(abs(pei([j])-([dfi1])))
      if ([rms3].lt.[dfi1])  then
        ve/inp keiar([j]) 10000
      else
        ve/inp keiar([j]) 1
      endif
    enddo
    ve/inp keiar(1) 10000
    ve/inp keiar(2) 10000
    ve/inp keiar(3) 10000
    ve/inp keiar(4) 10000
    ve/inp keiar([n]) 10000
    n=[n]-1
    ve/inp keiar([n]) 1
    ve/cre eka([np]) r [np]*[rms]
    sigma eka = eka*keiar
    ve/fit xei pei eka fourerfa.f sb 14 p14 s14
*    
    ve/copy p14(13:14) p12(11:12)
    ve/fit xei mueic ek fourerf.f sb 12 p12 s12 p12min p12max
    ve/copy p12(1:12) p14(1:12)
*    
  enddo
  ve/write p14 [dir]/amplitude_vs_times_fourerf_counter[nc]_[dir].par
  ve/del p12
  ve/copy p14(1:12) p12
endif
if ([nc].eq.-3) then
  ve/cre s7(7) r 0.1 0.1 0.1 10 1 5 0.01
  ve/cre p7(7) r 0 6.4147 3.7943 -1281.5 11.770 460.51 0.1
  sigma ek = dmueic
  ve/inp ek(1) 1000
  ve/inp ek(2) 1000
  ve/inp ek(3) 1000
  ve/inp ek(4) 1000
  ve/inp ek(6) 1000
  ve/inp ek(21) 1000
  ve/inp ek(22) 1000
  ve/inp ek(23) 1000
  ve/inp ek(39) 1000
*  ve/fit xei mueic ek gexp0.f sb 7 p7 s7
  sigma ek = dpei
  ve/inp ek(1) 1000
  ve/inp ek(2) 1000
  ve/inp ek(3) 1000
  ve/inp ek(4) 1000
  ve/inp ek(6) 1000
  ve/inp ek(21) 1000
  ve/inp ek(22) 1000
  ve/inp ek(23) 1000
  ve/inp ek(39) 1000
*  ve/cre s7(7) r 6*0 0.01
  set plci 4
  ve/fit xei pei ek gexp0.f sb 7 p7 s7
  ve/write p7 [dir]/amplitude_vs_times_counter[nc]_[dir].pars
  sigma ek = dmueic
  ve/inp ek(1) 1000
  ve/inp ek(2) 1000
  ve/inp ek(3) 1000
  ve/inp ek(4) 1000
  ve/inp ek(6) 1000
  ve/inp ek(21) 1000
  ve/inp ek(22) 1000
  ve/inp ek(23) 1000
  ve/inp ek(39) 1000
  ve/inp s7(5) 0
  ve/inp s7(6) 0
  ve/fit xei mueic ek gexp0.f sb 7 p7 s7
  set plci 6
*  fun/pl gexp0s.f -90 800 s
endif
if ([nc].eq.-7) then
  ve/cre s7(7) r 0.1 0.1 0.1 10 1 5 0
  ve/cre p7(7) r 0 6.4147 3.7943 -1281.5 11.770 460.51 1
  sigma ek = dmueic
  ve/inp ek(37) 1000
  ve/inp ek(39) 1000
  ve/fit xei mueic ek gexp0.f sb 7 p7 s7
  sigma ek = dpei
  ve/inp ek(1) 1000
  ve/inp ek(2) 1000
  ve/inp ek(6) 1000
  ve/inp ek(7) 1000
  ve/inp ek(18) 1000
  ve/inp ek(19) 1000
  ve/inp ek(20) 1000
  ve/inp ek(21) 1000
  ve/inp ek(22) 1000
  ve/inp ek(24) 1000
  ve/inp ek(39) 1000
  ve/cre s7(7) r 6*0 0.01
  set plci 4
  ve/fit xei pei ek gexp0.f sb 7 p7 s7
endif
if ([nc].eq.-9) then
  ve/cre s7(7) r 0.1 0.1 0.1 10 1 5 0
  ve/cre p7(7) r 0 6.4 5.5 -12810 11.770 460.51 1
  sigma ek = dmueic
  ve/inp ek(38) 1000
  ve/inp ek(39) 1000
  ve/fit xei mueic ek gexp0.f sb 7 p7 s7
  sigma ek = dpei
  ve/inp ek(1) 1000
  ve/inp ek(2) 1000
  ve/inp ek(6) 1000
  ve/inp ek(7) 1000
  ve/inp ek(18) 1000
  ve/inp ek(19) 1000
  ve/inp ek(20) 1000
  ve/inp ek(21) 1000
  ve/inp ek(22) 1000
  ve/inp ek(24) 1000
  sigma ek = dpei
  ve/inp ek(15) 1000
  ve/inp ek(19) 1000
  ve/inp ek(20) 1000
  ve/inp ek(22) 1000
  ve/inp ek(37) 1000
  ve/inp ek(38) 1000
  ve/inp ek(39) 1000
  ve/cre s7(7) r 6*0 0.01
  set plci 4
  ve/fit xei pei ek gexp0.f sb 7 p7 s7
endif
*
txt=[t]?[nc]! = $sigma(int(p3(3)*10+0.5)/10) [\261] $sigma(int(dp3(3)*10+0.5)/10)
*exec $PER/s#tf 0.6 0.1 [txt]
*
set plci 4
ve/fit xei pei ek expp0.f s 3 p3a
txt=A?[nc]!, pe
atitle 'T, days' [txt]
*
if ($vlen(id1).lt.10) then
set basl 0.01
set ltyp 15
do i=1,$vlen(id1)
  exec mapcal#ndays $sigma(id1([i]))
  gl/imp ndays
  line [ndays] $GRAFINFO('WNYMIN') [ndays] $GRAFINFO('WNYMAX')
enddo
endif
set ltyp 1
*
exec save [dir]/amplitude_vs_times_counter[nc]_[dir]_p$vlen(ti).eps f
*
set plci 1
exec vpl#pl mueic dmueic xei dxei sz=0.05
fun/pl fourerfp.f -100 800 s
exec save [dir]/amplitude_vs_times_fourerf0_counter[nc]_[dir]_p$vlen(ti).eps f
*
set lwid 1
set fwid 1
set hwid 1
null -40 720 0 $sigma(p14(1)*1.2)
set pmci 2
** exec $PER/s#vpl mueic dmueic xei dxei sz=0.05 o=s
exec vpl#pl mueic dmueic xei dxei sz=0.05 o=s
set pmci 4
** exec $PER/s#vpl pei dpei xei dxei sz=0.05 o=s
exec vpl#pl pei dpei xei dxei sz=0.05 o=s
set plci 2
fun/pl fourerfp.f -100 800 s
fun/pl fourerfpa.f -100 800 s
exec save [dir]/amplitude_vs_times_fourerf_counter[nc]_[dir]_p$vlen(ti).eps f
*ve/write pei,dpei,effi,deffi,mueic,dmueic,r,dr,xei,dxei,kei [dir]/parday_counter[nc]_[dir].txt 11e15.6
*
null -40 720 0 0.4
set pmci 4
sigma r = min(r,0.4)
sigma dr = min(dr,0.4)
exec vpl#pl r dr xei dxei sz=0.05 o=s
ve/copy p12(11:12) p2(1:2)
ve/cre s2(2) r
do i=1,1
  ve/fit xei r dr fp1.f sb 2 p2 s2
enddo
txt=[s]?A,[nc]! / A?[nc]!
atitle 'T, days' [txt]
*
*
if ($vlen(id1).lt.10) then
set basl 0.01
set ltyp 15
do i=1,$vlen(id1)
  exec mapcal#ndays $sigma(id1([i]))
  gl/imp ndays
  line [ndays] $GRAFINFO('WNYMIN') [ndays] $GRAFINFO('WNYMAX')
enddo
endif
set ltyp 1
*
exec save [dir]/uniformity_vs_times_counter[nc]_[dir]_p$vlen(ti).eps f
*
ip=0
ip=[ip]+1; it1[ip]=-30; it2[ip]=180;
ip=[ip]+1; it1[ip]=380; it2[ip]=490; 
ip=[ip]+1; it1[ip]=680; it2[ip]=710;
do i=1,[ip]
  null [it1[i]] [it2[i]] 0 0.4
  set pmci 4
  sigma r = min(r,0.4)
  sigma dr = min(dr,0.4)
  exec vpl#pl r dr xei dxei sz=0.05 o=s
  ve/cre s2(2) r
  ve/fit xei r dr fp1.f sb 2 p2 s2
  txt=[s]?A,[nc]! / A?[nc]!
  atitle 'T, days' [txt]
*
  exec save [dir]/uniformity_vs_times_counter[nc]_[dir]_p$vlen(ti)_t[i].eps f
enddo
*
ip=0
ip=[ip]+1; it1[ip]=-30; it2[ip]=180;
ip=[ip]+1; it1[ip]=380; it2[ip]=490; 
ip=[ip]+1; it1[ip]=680; it2[ip]=710;
do i=1,[ip]
  set lwid 1  
  set fwid 1
  set hwid 1
  null [it1[i]] [it2[i]] 0 $sigma(p14(1)*1.2)
  set pmci 2
  exec vpl#pl mueic dmueic xei dxei sz=0.05 o=s
  set pmci 4
  exec vpl#pl pei dpei xei dxei sz=0.05 o=s
  set plci 2
  fun/pl fourerfp.f -100 800 s
  fun/pl fourerfpa.f -100 800 s
  exec save [dir]/amplitude_vs_times_fourerf_counter[nc]_[dir]_p$vlen(ti)_t[i].eps f
enddo
return


macro tix
ve/del id1,id2
ve/read id1,id2 daycuts.txt 2f10.2
n=$vlen(id1)
np=[n]+1
ve/cre ti([np]) r
*
  exec mapcal#ndays $sigma(id1(1))
  gl/imp ndays
  t1=[ndays]-0.007
  ve/inp ti(1) [t1]
*
do i=2,[n]
  j=[i]-1
  exec mapcal#ndays $sigma(id2([j]))
  gl/imp ndays
  t1=[ndays]
  exec mapcal#ndays $sigma(id1([i]))
  gl/imp ndays
  t2=[ndays]
  ve/inp ti([i]) $sigma(([t1]+[t2])/2)
enddo
*
  exec mapcal#ndays $sigma(id2([n]))
  gl/imp ndays
  t2=[ndays]+0.007
  ve/inp ti([np]) [t2]
*
return


macro timecuts dd=7 drun=200
dir=v2
ve/read runs,days,amp1,damp1 [dir]/amplitude_vs_runs_counter1.txt 4f15.6
ve/del id1,id2
ve/read id1,id2 daycutsx.txt 2f10.2
n=$vlen(id1)
ve/cre in1([n]) r
ve/cre in2([n]) r
m=$vlen(runs)
do j=1,[m]
  runj=runs([j])
  do i=1,[n]
    idi=id1([i])
    if ([idi].eq.[runj]) then
      ve/inp in1([i]) [j]
    endif
    if ([idi].eq.[runj]) then
      ve/inp in2([i]) [j]
    endif
  enddo
enddo
*
ve/draw days
*
set basl 0.01
set ltyp 15
n=$vlen(runs)
n=[n]-1
nmax=100
ve/cre is1([nmax]) r
ind=0
do i=1,[n]
*  exec mapcal#ndays $sigma(id1([i]))
*  gl/imp ndays
  d1=days([i])
  j=[i]+1
  d2=days([j])
  if ([d2]-[d1].gt.7) then
    line [j] $GRAFINFO('WNYMIN') [j] $GRAFINFO('WNYMAX')
    mess [i] [j]
    ind=[ind]+1
    ve/inp is1([ind]) [j]
  endif
enddo
exec mapcal#vecut is1
set ltyp 1
*
set basl 0.01
set ltyp 15
set plci 2
do i=1,$vlen(in1)
  j=in1([i])
  line [j] $GRAFINFO('WNYMIN') [j] $GRAFINFO('WNYMAX')
enddo
set ltyp 1
set plci 1
*
exec vappend is1 in1
sigma ic = order(is1,is1)
*
nmax=100
ve/cre is1([nmax]) r
ind=0
n=$vlen(ic)
n=[n]-1
set basl 0.01
set ltyp 15
set plci 4
do i=1,[n]
  n1=ic([i])
  j=[i]+1
  n2=ic([j])
  if ([n2]-[n1].gt.2*[drun]) then
    m=$sigma(int(([n2]-[n1])/[drun])-1)
    dn=$sigma(int(([n2]-[n1])/([m]+1)))
    do j=1,[m]
      ind=[ind]+1
      n1=[n1]+[dn]
      ve/inp is1([ind]) [n1]
      line [n1] $GRAFINFO('WNYMIN') [n1] $GRAFINFO('WNYMAX')
    enddo
  endif
enddo
set ltyp 1
set plci 1
exec mapcal#vecut is1
*
exec vappend ic is1
sigma ic = order(ic,ic)
do i=1,$vlen(ic)
  ic1=ic([i])
  j=[i]+1
  do l=[j],$vlen(ic)
    ic2=ic([l])
    if ([ic1].eq.[ic2]) then
      mess [i] [l]
      ve/inp ic([l]) 0
    endif
  enddo
enddo
sigma ic = order(ic,-ic)
exec mapcal#vecut ic
sigma ic = order(ic,ic)
*
n=$vlen(ic)
ve/cre ti([n]) r
do i=1,[n]
  j=ic([i])
  exec mapcal#ndays $sigma(runs([j]))
  gl/imp ndays
  ndays=[ndays]-0.001
  ve/inp ti([i]) [ndays]  
enddo
n=$sigma(vmax(runs))
exec mapcal#ndays [n]
gtime(1)=gtime(1)+0.001
exec vappend ti gtime
return



macro pardayi dir=v2 nc1=1 nc2=9
gl/cre dir [dir]
*exec mapcal#cuts
*exec mapcal#empstab [dir]
*exec mapcal#daycut
*ve/read in1,in2 daycuts.txt 2f10.2
exec mapcal#timecuts 7 100 
*exec mapcal#tix
ve/read runs,days,amp1,damp1 [dir]/amplitude_vs_runs_counter1.txt 4f15.6
sigma ddays = 0*days
n=$vlen(runs)
ve/cre a1pe(9,[n]) r
ve/read a1pe amp1pe.txt
ve/cre st(9) r 0 8*1
gl/cre dir [dir]
ve/cre dshft(12) r -23 -21 10 46 66 128 137 139 158 169 161 162
ve/cre shft(12) r 0.69 0.4 0.9 0.2 0.3 0.79 0.97 0.9 0.99 0.25 0.2 0.9
do nc=[nc1],[nc2]
  if ([nc].eq.0) then
    do i=1,8
      ve/del runs,days,amp[i],damp[i]
      ve/read runs,days,amp[i],damp[i] [dir]/amplitude_vs_runs_counter[i].txt 4f15.6
      ve/del runs,days,n0[i]b,nx[i]b
      ve/read runs,days,n0[i]b,nx[i]b [dir]/n0_nx_vs_runs_counter[i].txt 4f15.6
    enddo
    np=$vlen(runs)
    ve/cre  amp0([np]) r
    ve/cre damp0([np]) r
    ve/cre n00b([np]) r
    ve/cre nx0b([np]) r
    do i=1,8
      sigma amp0 = amp0 + amp[i]/damp[i]**2
      sigma damp0 = damp0 + 1.0/damp[i]**2 
      sigma n00b = n00b + n0[i]b
      sigma nx0b = nx0b + nx[i]b
    enddo
    sigma amp0 = amp0/damp0
    sigma damp0 = 1.0/sqrt(damp0)
  else
    ve/del runs,days,amp[nc],damp[nc]
    ve/read runs,days,amp[nc],damp[nc] [dir]/amplitude_vs_runs_counter[nc].txt 4f15.6
    ve/del runs,days,n0[nc]b,nx[nc]b
    ve/read runs,days,n0[nc]b,nx[nc]b [dir]/n0_nx_vs_runs_counter[nc].txt 4f15.6
  endif
  exec mapcal#clblist amp[nc] damp[nc] $sigma(0*st([nc]))
  exec mapcal#pardayc [nc] dir=[dir]
*  exec mapcal#pardaycn [nc] dir=[dir]
enddo
return

macro parday par=amp1 dpar=damp1 nc=1
ve/del pei,dpei,xei,dxei,aeib,aeie
gl/imp dir
do i=1,$vlen(if1)
  exec mapcal#vpli np=[i] v=[par] dv=[dpar] ! nc=[nc]
  exec save [dir]/[par]_vs_time[i]_counter[nc]_[dir].eps f
  if ($vexist(p0i).eq.1) then
    if ([i].eq.1) then
      ve/copy  p0i  pei
      ve/copy dp0i dpei
      ve/copy  x0i  xei
      ve/copy dx0i dxei
      ve/copy a0ib aeib
      ve/copy a0ie aeie
    else
      exec vappend  pei  p0i
      exec vappend dpei dp0i
      exec vappend  xei  x0i
      exec vappend dxei dx0i
      exec vappend aeib a0ib
      exec vappend aeie a0ie
    endif
  endif
enddo
return

macro vpli np=1 v=amp1 dv=damp1 o=0 nc=1
i1=if1([np])
i2=if2([np])
ltr=$substring([v],1,1)
ve/del xir,dxir,yir,dyir,air
ve/copy  [v]([i1]:[i2]) yir
ve/copy [dv]([i1]:[i2]) dyir
ve/copy days([i1]:[i2]) xir
ve/copy ddays([i1]:[i2]) dxir
ve/copy a1pe([nc],[i1]:[i2]) air
set pmci 4
*** exec $PER/s#vpl yir dyir xir dxir sz=0.1 o=[o]
if ([ltr].ne.'a') then
  ve/del effx,deffx
  sigma effx = yir/dyir
  sigma deffx = sqrt(effx*(1-effx)/dyir)
  exec vpl#pl effx deffx xir dxir sz=0.1 o=[o]
else
  ve/cre p0(1) r 8
  ve/fit xir yir dyir p0 ! 1 p0
  exec hsigma yirm = yir
  exec hsigma dyirm = dyir
  do i=1,$vlen(xir)
    if ($sigma(dyirm([i])).gt.1) then
      ve/inp yirm([i]) $sigma(p0(1))
      ve/inp dyirm([i]) 1
    endif
  enddo
  exec vpl#pl yirm dyirm xir dxir sz=0.1 o=[o]
  exec vpl#pl yir dyir xir dxir sz=0.1 o=x
endif
ve/del xip,dxip,yip,dyip
ve/copy  yi([np]:[np]) yip
ve/copy dyi([np]:[np]) dyip
ve/copy  xi([np]:[np]) xip
ve/copy dxi([np]:[np]) dxip
set pmci 2
*** exec $PER/s#vpl yip dyip xip dxip sz=0.1 o=s
exec vpl#pl yip dyip xip dxip sz=0.1 o=s
xmin=$sigma(vmin(xir))
xmax=$sigma(vmax(xir))
id1=$sigma(int([xmin]))
id2=$sigma(int([xmax])+1)
ve/del p0i,dp0i,x0i,dx0i,a0ib,a0ie
ve/cre p0(1) r 
ve/cre dp0(1) r
j1=1
do i=[id1],[id2]
  dx=0
  do j=1,$vlen(dshft)
    dsh=dshft([j])
    if (([dsh]).eq.([i])) then
      dx=shft([j])
    endif
  enddo
  if ([dx].eq.0) then
    dx=8.5/24
  endif
  x=$sigma([i]+[dx])
  line [x] $GRAFINFO('WNYMIN') [x] $GRAFINFO('WNYMAX')
  dxmin=10000
  j2=0
  do j=1,$vlen(xir)
    dxmini=$sigma(abs(xir([j])-[x]))
*    if (([dxmini].lt.[dxmin])) then
    if ($sigma(xir([j])).lt.[x]) then
      dxmin=[dxmini]
      j2=[j]
    endif
  enddo
  if ([j1].lt.[j2]) then
    mess [j1] [j2]
    if ([ltr].eq.'a') then
    ve/del xirf,dxirf,yirf,dyirf,airf
    ve/copy  air([j1]:[j2]) airf
    ve/copy  xir([j1]:[j2]) xirf
    ve/copy dxir([j1]:[j2]) dxirf
    ve/copy  yir([j1]:[j2]) yirf
    ve/copy dyir([j1]:[j2]) dyirf
    ve/fit xirf yirf dyirf p0 s 1 p0 ! ! ! dp0
*    ve/fit xirf yirf dyirf p1 s
    if ($vexist(p0i).eq.0) then
      ve/cre p0i(1) r 
      ve/copy p0(1) p0i(1)
      ve/cre dp0i(1) r 
      ve/copy dp0(1) dp0i(1)
      ve/cre x0i(1) r $sigma((vmax(xirf)+vmin(xirf))/2)
      ve/cre dx0i(1) r $sigma((vmax(xirf)-vmin(xirf))/2)
      ve/cre a0ib(1) r $sigma(vmin(airf))
      ve/cre a0ie(1) r $sigma(vmax(airf))
    else
      exec vappend p0i p0
      exec vappend dp0i dp0
      ve/cre x0(1) r $sigma((vmax(xirf)+vmin(xirf))/2)
      exec vappend x0i x0
      ve/cre dx0(1) r $sigma((vmax(xirf)-vmin(xirf))/2)
      exec vappend dx0i dx0
      ve/cre a0i(1) r $sigma(vmin(airf))
      exec vappend a0ib a0i
      ve/cre a0i(1) r $sigma(vmax(airf))
      exec vappend a0ie a0i
    endif
    else
      ve/del xirf,dxirf,yirf,dyirf,airf
      ve/copy  air([j1]:[j2]) airf
      ve/copy  xir([j1]:[j2]) xirf
      ve/copy dxir([j1]:[j2]) dxirf
      ve/copy  yir([j1]:[j2]) yirf
      ve/copy dyir([j1]:[j2]) dyirf
      n0=$sigma(vsum(yirf))
      nx=$sigma(vsum(dyirf))
      ve/inp p0(1) $sigma([n0]/[nx])
      ve/inp dp0(1) $sigma(sqrt(p0(1)*(1-p0(1))/[nx]))
      if ($vexist(p0i).eq.0) then
        ve/cre p0i(1) r 
        ve/copy p0(1) p0i(1)
        ve/cre dp0i(1) r 
        ve/copy dp0(1) dp0i(1)
        ve/cre x0i(1) r $sigma((vmax(xirf)+vmin(xirf))/2)
        ve/cre dx0i(1) r $sigma((vmax(xirf)-vmin(xirf))/2)
        ve/cre a0ib(1) r $sigma(vmin(airf))
        ve/cre a0ie(1) r $sigma(vmax(airf))
      else
        exec vappend p0i p0
        exec vappend dp0i dp0
        ve/cre x0(1) r $sigma((vmax(xirf)+vmin(xirf))/2)
        exec vappend x0i x0
        ve/cre dx0(1) r $sigma((vmax(xirf)-vmin(xirf))/2)
        exec vappend dx0i dx0
        ve/cre a0i(1) r $sigma(vmin(airf))
        exec vappend a0ib a0i
        ve/cre a0i(1) r $sigma(vmax(airf))
        exec vappend a0ie a0i
      endif
    endif
  endif
  mess [i]: [j1] [j2]
  j1=[j2]+1
enddo
return



macro caltest
n=$vlen(nr1)
im=[n]-1
do i=1,[im]
  ii=[i]+1
  mess [i] $sigma(nr1([ii])-nr2([i]))
enddo
return


macro ampeff v=amp1 dv=damp1 nc=1
*ve/cre ti(16) r -50 20 40 60 80 100 120 135 150 300 420 435 450 460 600 800
*exec mapcal#timecuts
n=$vlen(ti)
n=[n]-1
ve/cre  xi([n]) r
ve/cre dxi([n]) r
ve/cre  yi([n]) r
ve/cre dyi([n]) r
ve/cre  a1([n]) r
ve/cre da1([n]) r
do i=1,[n]
  ii=[i]+1
  t1=ti([i])
  t2=ti([ii])
  j1=10000
  j2=0
  do j=1,$vlen([v])
    t=days([j])
    if (([t].gt.[t1]).and.([t].lt.[t2])) then
      if ([j1].gt.[j]) then
        j1=[j]
      endif
      if ([j2].lt.[j]) then
        j2=[j]
      endif
    endif
  enddo
  mess [j1] [j2] [t1] [t2]
  dj=[j2]-[j1]+1
  ve/del [v]r,[dv]r
  ve/copy [v]([j1]:[j2]) [v]r
  ve/copy [dv]([j1]:[j2]) [dv]r
*  ve/inp  yi([i]) $sigma(vsum([v]r)/[dj])
  ve/inp  yi([i]) $sigma(vsum([v]r/[dv]r**2)/vsum(1/[dv]r**2))
  ve/inp dyi([i]) $sigma(1/sqrt(vsum(1/[dv]r**2)))
*  ve/inp dyi([i]) $sigma(sqrt(vsum([dv]r**2))/[dj])
  ve/inp  xi([i]) $sigma((days([j2])+days([j1]))/2)
  ve/inp dxi([i]) $sigma((days([j2])-days([j1]))/2)
  if ($substring([v],1,1).eq.'n') then
    ve/inp  yi([i]) $sigma(vsum([v]r)/vsum([dv]r))
    ve/inp dyi([i]) $sigma(sqrt(yi([i])*(1-yi([i]))/vsum([dv]r)))
  endif
  if ($vexist(a1pe).eq.1) then
    a1j1=a1pe([nc],[j1])
    a1j2=a1pe([nc],[j2])
    ve/copy a1([i]) a1[j1]
    if ([a1j1].ne.[a1j2]) then
      mess '>>>>>' [i] [j1] [j2] [a1j1] [a1j2]
    endif
  endif
enddo
* exec $PER/s#vpl yi dyi xi dxi sz=0.1
return

macro ampcr nc=1 mode=1
*exec mapcal#ampeff amp[nc] damp[nc]
exec mapcal#clblist amp[nc] damp[nc] 1
ve/del x1,dx1,y1,dy1
ve/copy  xi  x1
ve/copy dxi dx1
ve/copy  yi  y1
ve/copy dyi dy1
*exec mapcal#ampeff n0[nc]b nx[nc]b
exec mapcal#clblist n0[nc]b nx[nc]b 1
ve/del x2,dx2,y2,dy2
ve/copy  xi  x2
ve/copy dxi dx2
ve/copy  yi  y2
ve/copy dyi dy2
if ([mode].eq.1) then
  sigma r = sqrt(2*(y1+log(y2))/y1**2)
*  sigma dr = sqrt((dy2/y2/y1)**2+(log(y2)*dy1/y1**2)**2)
  sigma dr = sqrt((dy1*(1+2*log(y2)/y1))**2+(dy2/y2)**2)/r/y1**2
  set pmci 2
  * exec $PER/s#vpl r dr x1 dx1 sz=0.05
endif
if ([mode].eq.2) then
  set pmci 2
  * exec $PER/s#vpl y1 dy1 x1 dx1 sz=0.1
  sigma y3 = -log(y2)
  sigma dy3 = dy2/y2
  set pmci 4
  * exec $PER/s#vpl y3 dy3 x2 dx2 sz=0.1 o=s
endif
return


macro zdays
do i=1,$vlen(days)
  iday=days([i])
  if ([iday].eq.0) then
    mess [i] $sigma(runs([i]))
  endif
enddo
n=$vlen(days)
do i=2,$sigma([n]-1)
  dd=days([i])
  if ([dd].eq.0) then
    ii=[i]-1
    d1=days([ii])
    ii=[i]+1
    d2=days([ii])
    ve/inp days([i]) $sigma(([d2]+[d1])/2)
    mess [i] [d1] [dd] [d2]
  endif
enddo
return


macro newcal nc=1
n=0
*n=[n]+1; cal[n]=2011_Jan_12_0041          
*n=[n]+1; cal[n]=2011_Jan_21_0042          
*n=[n]+1; cal[n]=2011_Jan_24_0043          
*n=[n]+1; cal[n]=2011_Feb_09_0044          
*n=[n]+1; cal[n]=2011_Feb_11_0045          
*n=[n]+1; cal[n]=2011_Feb_13_0046          
*n=[n]+1; cal[n]=2011_Feb_16_0047          
*n=[n]+1; cal[n]=2011_Feb_21_0048          
*n=[n]+1; cal[n]=2011_Feb_25_0049          
*n=[n]+1; cal[n]=2011_Feb_28_0050          
*n=[n]+1; cal[n]=2011_Mar_09_0051          
*n=[n]+1; cal[n]=2011_Mar_14_0052          
*n=[n]+1; cal[n]=2011_Mar_16_0053          
*n=[n]+1; cal[n]=2011_Mar_21_0054          
*n=[n]+1; cal[n]=2011_Mar_25_0055          
*n=[n]+1; cal[n]=2011_Mar_28_0056          
*n=[n]+1; cal[n]=2011_Mar_30_0057          
*n=[n]+1; cal[n]=2011_Apr_04_0058          
*n=[n]+1; cal[n]=2011_Apr_10_0059          
*n=[n]+1; cal[n]=2011_Apr_14_0060          
*n=[n]+1; cal[n]=2011_Apr_27_0061          
*n=[n]+1; cal[n]=2011_Apr_27_0062          
*n=[n]+1; cal[n]=2011_May_05_0063          
*n=[n]+1; cal[n]=2011_May_05_0064          
*n=[n]+1; cal[n]=2011_May_11_0065          
*n=[n]+1; cal[n]=2011_May_14_0066          
*n=[n]+1; cal[n]=2011_May_17_0067          
*n=[n]+1; cal[n]=2011_May_20_0068          
*n=[n]+1; cal[n]=2011_May_25_0069          
*n=[n]+1; cal[n]=2011_May_25_0070          
*n=[n]+1; cal[n]=2011_May_31_0071          
*n=[n]+1; cal[n]=2011_Jun_02_0072          
*n=[n]+1; cal[n]=2011_Jun_02_0073          
*n=[n]+1; cal[n]=2011_Jun_07_0074          
*n=[n]+1; cal[n]=2011_Jun_09_0075          
*n=[n]+1; cal[n]=2011_Jun_15_0076          
*n=[n]+1; cal[n]=2011_Jun_17_0077          
*n=[n]+1; cal[n]=2011_Oct_24_0080          
*n=[n]+1; cal[n]=2011_Oct_24_0081          
*n=[n]+1; cal[n]=2011_Oct_25_0082          
*n=[n]+1; cal[n]=2011_Oct_26_0083          
*n=[n]+1; cal[n]=2011_Oct_27_0084          
*n=[n]+1; cal[n]=2011_Oct_28_0085          
*n=[n]+1; cal[n]=2011_Nov_03_0086          
*n=[n]+1; cal[n]=2011_Nov_03_0087          
*n=[n]+1; cal[n]=2011_Nov_21_0090          
*n=[n]+1; cal[n]=2011_Nov_25_0091          
*n=[n]+1; cal[n]=2011_Nov_25_0092          
*n=[n]+1; cal[n]=2011_Nov_25_0093          
*n=[n]+1; cal[n]=2011_Nov_29_0094          
*n=[n]+1; cal[n]=2011_Nov_30_0095          
*n=[n]+1; cal[n]=2011_Dec_05_0096          
*n=[n]+1; cal[n]=2011_Dec_13_0097          
*n=[n]+1; cal[n]=2011_Dec_19_0098          
*n=[n]+1; cal[n]=2012_Jun_10_0099          
*n=[n]+1; cal[n]=2012_Jun_20_0100          
*n=[n]+1; cal[n]=2012_Jun_25_0101          
*n=[n]+1; cal[n]=2012_Jun_25_0102          
*n=[n]+1; cal[n]=2012_Jun_25_0103          
*n=[n]+1; cal[n]=2012_Feb_01_0104          
*n=[n]+1; cal[n]=2012_Feb_01_0105          
*n=[n]+1; cal[n]=2012_Feb_03_0106          
*n=[n]+1; cal[n]=2012_Feb_03_0107          
*n=[n]+1; cal[n]=2012_Feb_07_0108          
*n=[n]+1; cal[n]=2012_Feb_07_0109          
*n=[n]+1; cal[n]=2012_Feb_13_0110          
*n=[n]+1; cal[n]=2012_Feb_28_0111          
*n=[n]+1; cal[n]=2012_Mar_05_0112          
*n=[n]+1; cal[n]=2012_Mar_12_0113          
*n=[n]+1; cal[n]=2012_Mar_18_0114          
*n=[n]+1; cal[n]=2012_Mar_23_0115          
*n=[n]+1; cal[n]=2012_Mar_30_0116          
*n=[n]+1; cal[n]=2012_Mar_30_0117          
*n=[n]+1; cal[n]=2012_Apr_10_0118          
*n=[n]+1; cal[n]=2012_Apr_16_0119          
*n=[n]+1; cal[n]=2012_Apr_23_0120          
*n=[n]+1; cal[n]=2012_Apr_26_0122          
*n=[n]+1; cal[n]=2012_Oct_03_0123          
*n=[n]+1; cal[n]=2012_Oct_10_0124          
*n=[n]+1; cal[n]=2012_Oct_10_0125          
*n=[n]+1; cal[n]=2012_Oct_15_0126          
*n=[n]+1; cal[n]=2012_Nov_19_0127          
*n=[n]+1; cal[n]=2012_Nov_20_0128          
*n=[n]+1; cal[n]=2012_Nov_30_0131          
*n=[n]+1; cal[n]=2012_Dec_03_0135          
*n=[n]+1; cal[n]=2012_Dec_10_0136          
*n=[n]+1; cal[n]=2012_Dec_18_0137          
*n=[n]+1; cal[n]=2012_Dec_18_0138          
*n=[n]+1; cal[n]=2012_Dec_24_0139          
*n=[n]+1; cal[n]=2013_Jan_02_0140          
*n=[n]+1; cal[n]=2013_Jan_04_0141          
*n=[n]+1; cal[n]=2013_Jan_09_0142          
*n=[n]+1; cal[n]=2013_Jan_14_0143          
*n=[n]+1; cal[n]=2013_Jan_17_0144          
*
*n=[n]+1; cal[n]=2010_dec_06_0037 
*n=[n]+1; cal[n]=2010_dec_07_0038 
*n=[n]+1; cal[n]=2010_dec_07_0039 
*n=[n]+1; cal[n]=2010_dec_07_0040 
n=[n]+1; cal[n]=2011_Jan_12_0041 
n1=[n]
n=[n]+1; cal[n]=2011_Jan_21_0042 
n=[n]+1; cal[n]=2011_Jan_24_0043 
n=[n]+1; cal[n]=2011_Feb_09_0044 
n=[n]+1; cal[n]=2011_Feb_11_0045 
n=[n]+1; cal[n]=2011_Feb_13_0046 
n=[n]+1; cal[n]=2011_Feb_16_0047 
n=[n]+1; cal[n]=2011_Feb_21_0048 
n=[n]+1; cal[n]=2011_Feb_25_0049 
n=[n]+1; cal[n]=2011_Feb_28_0050 
n=[n]+1; cal[n]=2011_Mar_09_0051 
n=[n]+1; cal[n]=2011_Mar_14_0052 
n=[n]+1; cal[n]=2011_Mar_16_0053 
n=[n]+1; cal[n]=2011_Mar_21_0054 
n=[n]+1; cal[n]=2011_Mar_25_0055 
*n=[n]+1; cal[n]=2011_Mar_28_0056 
n=[n]+1; cal[n]=2011_Mar_30_0057 
n=[n]+1; cal[n]=2011_Apr_04_0058 
n=[n]+1; cal[n]=2011_Apr_10_0059 
n=[n]+1; cal[n]=2011_Apr_14_0060 
n=[n]+1; cal[n]=2011_Apr_27_0061 
n=[n]+1; cal[n]=2011_Apr_27_0062 
n=[n]+1; cal[n]=2011_May_05_0063 
n=[n]+1; cal[n]=2011_May_05_0064 
n=[n]+1; cal[n]=2011_May_11_0065 
n=[n]+1; cal[n]=2011_May_14_0066 
n=[n]+1; cal[n]=2011_May_17_0067 
n=[n]+1; cal[n]=2011_May_20_0068 
n=[n]+1; cal[n]=2011_May_25_0069 
n=[n]+1; cal[n]=2011_May_25_0070 
*n=[n]+1; cal[n]=2011_May_31_0071 
n=[n]+1; cal[n]=2011_Jun_02_0072 
n=[n]+1; cal[n]=2011_Jun_02_0073 
*n=[n]+1; cal[n]=2011_Jun_07_0074 
n=[n]+1; cal[n]=2011_Jun_09_0075 
n=[n]+1; cal[n]=2011_Jun_15_0076 
n=[n]+1; cal[n]=2011_Jun_17_0077 
*n=[n]+1; cal[n]=2011_Oct_24_0080 
*n=[n]+1; cal[n]=2011_Oct_24_0081 
n=[n]+1; cal[n]=2011_Oct_25_0082 
n=[n]+1; cal[n]=2011_Oct_26_0083 
n=[n]+1; cal[n]=2011_Oct_27_0084 
n=[n]+1; cal[n]=2011_Oct_28_0085 
n=[n]+1; cal[n]=2011_Nov_03_0086 
n=[n]+1; cal[n]=2011_Nov_03_0087 
n=[n]+1; cal[n]=2011_Nov_21_0090 
*n=[n]+1; cal[n]=2011_Nov_25_0091 
n=[n]+1; cal[n]=2011_Nov_25_0092 
n=[n]+1; cal[n]=2011_Nov_25_0093 
n=[n]+1; cal[n]=2011_Nov_29_0094 
*n=[n]+1; cal[n]=2011_Nov_30_0095 
*n=[n]+1; cal[n]=2011_Dec_05_0096 
*n=[n]+1; cal[n]=2011_Dec_13_0097 
*n=[n]+1; cal[n]=2011_Dec_19_0098 
*n=[n]+1; cal[n]=2012_Jun_10_0099 
*n=[n]+1; cal[n]=2012_Jun_20_0100 
*n=[n]+1; cal[n]=2012_Jun_25_0101 
*n=[n]+1; cal[n]=2012_Jun_25_0102 
*n=[n]+1; cal[n]=2012_Jun_25_0103 
*n=[n]+1; cal[n]=2012_Feb_01_0104 
*n=[n]+1; cal[n]=2012_Feb_01_0105 
*n=[n]+1; cal[n]=2012_Feb_03_0106 
*n=[n]+1; cal[n]=2012_Feb_03_0107 
*n=[n]+1; cal[n]=2012_Feb_07_0108 
*n=[n]+1; cal[n]=2012_Feb_07_0109 
*n=[n]+1; cal[n]=2012_Feb_13_0110 
*n=[n]+1; cal[n]=2012_Feb_28_0111 
*n=[n]+1; cal[n]=2012_Mar_05_0112 
*n=[n]+1; cal[n]=2012_Mar_12_0113 
*n=[n]+1; cal[n]=2012_Mar_18_0114 
n=[n]+1; cal[n]=2012_Mar_23_0115 
n=[n]+1; cal[n]=2012_Mar_30_0116 
n=[n]+1; cal[n]=2012_Mar_30_0117 
*n=[n]+1; cal[n]=2012_Apr_10_0118 
n=[n]+1; cal[n]=2012_Apr_16_0119 
n=[n]+1; cal[n]=2012_Apr_23_0120 
n=[n]+1; cal[n]=2012_Apr_26_0122 
n=[n]+1; cal[n]=2012_Oct_03_0123 
*n=[n]+1; cal[n]=2012_Oct_10_0124 
n=[n]+1; cal[n]=2012_Oct_10_0125 
*n=[n]+1; cal[n]=2012_Oct_15_0126 
n=[n]+1; cal[n]=2012_Nov_19_0127 
*n=[n]+1; cal[n]=2012_Nov_20_0128 
n=[n]+1; cal[n]=2012_Nov_26_0129 
*n=[n]+1; cal[n]=2012_Nov_30_0131 
*n=[n]+1; cal[n]=2012_Nov_30_0132 
n=[n]+1; cal[n]=2012_Nov_30_0133 
n2=[n]
n=[n]+1; cal[n]=2012_Nov_30_0134 
*n=[n]+1; cal[n]=2012_Nov_30b_0134
n=[n]+1; cal[n]=2012_Dec_03_0135 
*n1=[n]
n=[n]+1; cal[n]=2012_Dec_10_0136 
n=[n]+1; cal[n]=2012_Dec_18_0137 
n=[n]+1; cal[n]=2012_Dec_18_0138 
n=[n]+1; cal[n]=2012_Dec_24_0139 
n=[n]+1; cal[n]=2013_Jan_02_0140 
*n=[n]+1; cal[n]=2013_Jan_31_0140 
n=[n]+1; cal[n]=2013_Jan_04_0141 
*n=[n]+1; cal[n]=2013_Jan_31_0141 
n=[n]+1; cal[n]=2013_Jan_09_0142 
*n=[n]+1; cal[n]=2013_Jan_17_0142 
n=[n]+1; cal[n]=2013_Jan_14_0143 
n=[n]+1; cal[n]=2013_Jan_17_0144 
n=[n]+1; cal[n]=2013_Jan_21_0145 
*n=[n]+1; cal[n]=2013_Feb_01_0146 
n=[n]+1; cal[n]=2013_Jan_27_0146 
n=[n]+1; cal[n]=2013_Jan_31_0147 
n=[n]+1; cal[n]=2013_Feb_01_0148 
n=[n]+1; cal[n]=2013_Feb_01_0149 
n=[n]+1; cal[n]=2013_Feb_06_0150 
n=[n]+1; cal[n]=2013_Feb_12_0151 
n=[n]+1; cal[n]=2013_Feb_15_0152 
n=[n]+1; cal[n]=2013_Feb_18_0153 
n=[n]+1; cal[n]=2013_Feb_21_0154 
n=[n]+1; cal[n]=2013_Apr_12_0155 
*n=[n]+1; cal[n]=2013_Apr_15_0155 
*n=[n]+1; cal[n]=2013_Feb_25_0156 
*n2=[n]
n=[n]+1; cal[n]=2013_Feb_26_0157 
n=[n]+1; cal[n]=2013_Mar_04_0158 
n=[n]+1; cal[n]=2013_Mar_06_0160 
n=[n]+1; cal[n]=2013_Mar_12_0161 
n=[n]+1; cal[n]=2013_Mar_18_0162 
n=[n]+1; cal[n]=2013_Mar_29_0163 
n=[n]+1; cal[n]=2013_Apr_01_0164 
n=[n]+1; cal[n]=2013_Apr_08_0165 
n=[n]+1; cal[n]=2013_Apr_15_0167 
n=[n]+1; cal[n]=2013_Apr_24_0168 
n=[n]+1; cal[n]=2013_Apr_25_0169 
n=[n]+1; cal[n]=2013_Apr_25_0170 
n=[n]+1; cal[n]=2013_Apr_29_0171 
n=[n]+1; cal[n]=2013_May_06_0172 
n=[n]+1; cal[n]=2013_May_13_0173 
n=[n]+1; cal[n]=2013_May_16_0174 
n=[n]+1; cal[n]=2013_May_23_0175 
n=[n]+1; cal[n]=2013_May_28_0176 
n=[n]+1; cal[n]=2013_Jun_06_0177 
n=[n]+1; cal[n]=2013_Jun_17_0178 
*n=[n]+1; cal[n]=2013_Jul_09_0179 
n=[n]+1; cal[n]=2013_Jun_23_0179 
n=[n]+1; cal[n]=2013_Jul_04_0180 
n=[n]+1; cal[n]=2013_Jul_10_0180 
*
ve/cre q1i([n]) r
nx=[n2]-[n1]+1
ve/cre u0ax([nx]) r
ve/cre a1x([nx]) r
ve/cre da1x([nx]) r
ve/cre u0bx([nx]) r
ve/cre b1x([nx]) r
ve/cre db1x([nx]) r
ind=0
do i=[n1],[n2]
  mess [i] [cal[i]] [n1] [n2]
  exec mapcal#newsat [nc] [cal[i]]
  if ([xx].eq.1) then
*  gl/imp q1
*  ve/inp q1i([i]) [q1]
  a1=p3(2)
  k=p3(1)
  a0=p3(3)
  ur=uw(1)
  l0=$sigma(vmin(u0))
  r0=$sigma(vmax(u0))
  l=$sigma([l0]-([r0]-[l0])/20)
  r=$sigma([r0]+([r0]-[l0])/20)
  if ([i].eq.[n1]) then
    fun/pl ([a1]-([a0]))*(x/[ur])**([k])+([a0]) [l] [r]
  else
    fun/pl ([a1]-([a0]))*(x/[ur])**([k])+([a0]) [l] [r] s
  endif
  endif
    ind=[ind]+1
    gl/imp a1xi 
    gl/imp da1xi
    ve/inp a1x([ind]) [a1xi]
    ve/inp da1x([ind]) [da1xi]
    gl/imp u0ai
    ve/inp u0ax([ind]) [u0ai]
    gl/imp b1xi 
    gl/imp db1xi
    ve/inp b1x([ind]) [b1xi]
    ve/inp db1x([ind]) [db1xi]
    gl/imp u0bi
    ve/inp u0bx([ind]) [u0bi]
*  ve/read x
enddo
sigma a1x = abs(a1x)
sigma b1x = abs(b1x)
mean=$sigma(vsum(a1x)/[nx])
mean2=$sigma(vsum(a1x**2)/[nx])
rms=$sigma(sqrt([mean2]-([mean])**2))
mess [mean] [rms]
ve/cre  u0axs([nx]) r
ve/cre  u0bxs([nx]) r
ve/cre  a1xs([nx]) r
ve/cre da1xs([nx]) r
ve/cre  b1xs([nx]) r
ve/cre db1xs([nx]) r
ind=0
do i=1,[nx]
  if ($sigma(abs(a1x([i])-[mean])/10).lt.3) then
    ind=[ind]+1
    ve/inp  u0axs([ind]) $sigma(u0ax([i]))
    ve/inp  u0bxs([ind]) $sigma(u0bx([i]))
    ve/inp  a1xs([ind]) $sigma(a1x([i]))
    ve/inp da1xs([ind]) $sigma(da1x([i]))
    ve/inp  b1xs([ind]) $sigma(b1x([i]))
    ve/inp db1xs([ind]) $sigma(db1x([i]))
  endif
enddo
exec mapcal#vecut a1xs
nx=$vlen(a1xs)
exec mapcal#vecut da1xs [nx]
exec mapcal#vecut  b1xs [nx]
exec mapcal#vecut db1xs [nx]
exec mapcal#vecut u0axs [nx]
exec mapcal#vecut u0bxs [nx]
ve/cre tx([nx]) r
ve/cre dtx([nx]) r
sigma tx=array([nx],1#[nx])
* exec $PER/s#vpl a1xs da1xs tx dtx
ve/cre p4(4) r 200 -50 50 50
ve/fit tx a1xs da1xs erf0.f s 4 p4
return                      



macro newsat ncnt=1 pref=2012_Nov_20_0128
ncal=$substring([pref],13,16)
dir=/work/users/konctbel/Calibr/[pref]
hfile=[dir]/[pref]_[ncnt].res

if ($fexist([hfile]).eq.1) then

exec mapcal#uw [dir]/[pref]_[ncnt].ps

np=7
ve/cre c0([np]) r 20 40 50 80 100 150 200

vname=u0;   ve/del [vname]; len=$len([vname])+1; fmt=([len]X,f15.6); ve/read [vname] [hfile] [fmt] ! /[vname]:/(1)
vname=du0;  ve/del [vname]; len=$len([vname])+1; fmt=([len]X,f15.6); ve/read [vname] [hfile] [fmt] ! /[vname]:/(1)
vname=kxrn; ve/del [vname]; len=$len([vname])+1; fmt=([len]X,f15.6); ve/read [vname] [hfile] [fmt] ! /[vname]:/(1)
vname=pn_200; ve/del [vname]; len=$len([vname])+1; fmt=([len]X,f15.6); ve/read [vname] [hfile] [fmt] ! /[vname]:/(1)
vname=dpn_200; ve/del [vname]; len=$len([vname])+1; fmt=([len]X,f15.6); ve/read [vname] [hfile] [fmt] ! /[vname]:/(1)
ax=pn_200(4);  gl/cre  a1xi $sigma([ax]*kxrn(7)); mess [a1xi]
ax=dpn_200(4); gl/cre da1xi $sigma([ax]*kxrn(7)); mess [da1xi]
vname=pn_150; ve/del [vname]; len=$len([vname])+1; fmt=([len]X,f15.6); ve/read [vname] [hfile] [fmt] ! /[vname]:/(1)
vname=dpn_150; ve/del [vname]; len=$len([vname])+1; fmt=([len]X,f15.6); ve/read [vname] [hfile] [fmt] ! /[vname]:/(1)
ax=pn_150(4);  gl/cre  b1xi $sigma([ax]*kxrn(6)); mess [b1xi]
ax=dpn_150(4); gl/cre db1xi $sigma([ax]*kxrn(6)); mess [db1xi]
gl/cre u0ai $sigma(u0(7))
gl/cre u0bi $sigma(u0(6))

goto 1

ve/cre  muc([np]) r
ve/cre smuc([np]) r
ve/cre  mu([np]) r
ve/cre smu([np]) r
ve/cre  pds([np]) r
ve/cre dpds([np]) r
ve/cre  rim([np]) r

do i=1,[np]
*
  c0t=c0([i])
*
  vname=pn_[c0t]; 
  ve/del [vname];
  len=$len([vname])+1;
  fmt=([len]X,f15.6);
  ve/read [vname] [hfile] [fmt] ! /[vname]:/(1)
  ve/copy [vname](1) pds([i])
*
  vname=dpn_[c0t]; 
  ve/del [vname];
  len=$len([vname])+1;
  fmt=([len]X,f15.6);
  ve/read [vname] [hfile] [fmt] ! /[vname]:/(1)
  ve/copy [vname](1) dpds([i])
*
  vname=ps_[c0t]; 
  ve/del [vname];
  len=$len([vname])+1;
  fmt=([len]X,f15.6);
  ve/read [vname] [hfile] [fmt] ! /[vname]:/(1)
*
  N0=$sigma(ps_[c0t](1))
  m0=$sigma(ps_[c0t](4))
  N1=$sigma(ps_[c0t](7))
  m1=$sigma(ps_[c0t](10))
  ve/inp  muc([i]) $sigma(log(1-[m0]/[N0])-log(1-[m1]/[N1]))
  ve/inp smuc([i]) $sigma(sqrt([m0]/[N0]/([N0]-[m0])+[m1]/[N1]/([N1]-[m1])))
  ve/inp  mu([i]) $sigma(-log(1-[m1]/[N1]))
  ve/inp smu([i]) $sigma(sqrt([m1]/[N1]/([N1]-[m1])))
*
*  exec mapcal#rimr [dir]/picts/cal_pmt_[ncal]_[c0t]_ch[ncnt]_RDeff.eps
*  gl/imp rimRD
*  ve/inp rim([i]) [rimRD]
enddo

exec mapcal#rimr1 [dir]/[pref]_[ncnt].ps

sigma pds=pds*kxrn
sigma dpds=dpds*kxrn
ve/cre di(1) r $sigma(vsum(rim-pds)/[np])

set pmci 4
set hcol 1
set plci 4
* exec $PER/s#vpl muc smuc u0 du0 iatt=20 sz=0.1
** exec $PER/s#vpl mu smu u0 du0 iatt=23 sz=0.1

*npar=5
*ve/cre chi2(2) r
*ve/cre paru([npar]) r
*ve/cre dparu([npar]) r
*ve/cre covu([npar],[npar]) r

ve/cre pmin(5) r 1 0.5 1e-12 50 0.1
ve/cre pmax(5) r 1 0.99 1e-8 200 10.
ve/cre p5(5) r 1 0.8 1e-10 100 $sigma(vmax(muc))
ve/cre s5(5) r 0 0 1e-12 0 0.01
ve/fit u0 muc smuc muv.f sb 5 p5 s5 pmin pmax
ve/cre s5(5) r 0 0 1e-12 1 0
ve/fit u0 muc smuc muv.f sb 5 p5 s5 pmin pmax
ve/cre s5(5) r 0 0.1 1e-12 1 0.01
ve/fit u0 muc smuc muv.f sb 5 p5 s5 pmin pmax
do i=1,1
  ve/fit u0 muc smuc muv.f sb 5 p5 s5 pmin pmax
enddo

ve/cre ui(1) r
im=0
dmuim=0
do i=1,[np]
  ve/inp ui(1) $sigma(u0([i]))
  mui=$call('muvp.f(ui)')
  dmui=$sigma(abs([mui]-muc([i]))/smuc([i]))
  if (([dmui].gt.[dmuim]).and.([dmui].gt.3)) then
    im=[i]
    dmuim=[dmui]
  endif
enddo

if ([im].ne.0) then
np1=[np]-1
ve/cre  u0r([np1]) r
ve/cre du0r([np1]) r
ve/cre  mucr([np1]) r
ve/cre smucr([np1]) r
ind=0
do i=1,[np]
  if ([i].ne.[im]) then
    ind=[ind]+1
    ve/inp  u0r([ind])  u0([i])
    ve/inp du0r([ind]) du0([i])
    ve/inp  mucr([ind])  muc([i])
    ve/inp smucr([ind]) smuc([i])
  endif
enddo

* exec $PER/s#vpl muc smuc u0 du0 iatt=20 sz=0.1
ve/fit u0r mucr smucr muv.f sb 5 p5 s5 pmin pmax

else

* exec $PER/s#vpl muc smuc u0 du0 iatt=20 sz=0.1
ve/fit u0 muc smuc muv.f sb 5 p5 s5 pmin pmax
mess '<<<<there is no out point>>>>'

endif

*call covm.f(1)
*call covmpen(chi2,[npar],paru,dparu)
*call covmcov([npar],covu)

*ve/cre dxx([npar]) r 
*do i=1,-[npar]
*  dx=covu([i],[i])
*  ve/inp dxx([i]) $sigma(sqrt([dx]))
*enddo
*ve/cre corr([npar],[npar]) r
*nparf=4
*do i=1,-[nparf]
*  do j=1,-[nparf]
*    dxy=covu([i],[j])
*    dxy=[dxy]/dxx([i])/dxx([j])
*    ve/inp corr([i],[j]) [dxy]
*  enddo
*enddo
*ve/write covu(1),covu(2),covu(3),covu(4),covu(5) ! 5f10.5
*ve/write corr(1),corr(2),corr(3),corr(4),corr(5) ! 5f10.5

ve/cre p4(4) r 0.000000 0.510029 2813.002686 -92.598007
ve/cre s4(4) r 0 0.01 10 10
set plci 2
do i=1,1
  ve/fit u0 muc smuc erf0.f sb 4 p4 s4
enddo

gl/cre q1 $call('muvp.f(uw)')

else
  gl/cre q1 0
endif

*mess R=$sigma(p5(5)/p4(2))

read x

1:
*
*vname=q1RDc;
*ve/del [vname];
*len=$len([vname])+1;
*fmt=([len]X,f15.6);
*ve/read [vname] [hfile]1 [fmt] ! /[vname]:/(1)
*ve/copy [vname] p3

return


macro rimr fname=/work/users/konctbel/Calibr/2012_Nov_20_0128/picts/cal_pmt_0128_80_ch3_RDeff.eps
if ($fexist([fname]).eq.1) then
  shell fgrep -e "\040)" [fname] >& tmp.txt
  shell $unquote('cat tmp.txt | sed "s/\\\040)/ /g" > tmp1.txt')
  shell $unquote('cat tmp1.txt | sed "s/(/ /g" > tmp2.txt')
  shell $unquote('cat tmp2.txt | sed "s/=/ /g" > tmp3.txt')
  shell $unquote('cat tmp3.txt | sed "s/*/ /g" > tmp4.txt')
  ve/del tmp
  ve/read tmp tmp4.txt
  gl/cre rimRD $sigma(tmp(2))
endif
return

macro rimr1 fname=/work/users/konctbel/Calibr/2011_Jan_12_0041/2011_Jan_12_0041_1.ps
if ($fexist([fname]).eq.1) then
  shell fgrep -e "\040)" [fname] >& tmp.txt
  shell $unquote('cat tmp.txt | sed "s/\\\040)/ /g" > tmp1.txt')
  shell $unquote('cat tmp1.txt | sed "s/(/ /g" > tmp2.txt')
  shell $unquote('cat tmp2.txt | sed "s/=/ /g" > tmp3.txt')
  shell $unquote('cat tmp3.txt | sed "s/*/ /g" > tmp4.txt')
  ve/del tmp
  ve/read tmp tmp4.txt
  ve/cre rim(7) r
  do i=1,7
    ind=54+[i]*4
    ve/inp rim([i]) tmp([ind])
  enddo
endif
return


macro uw fname=/work/users/konctbel/Calibr/2011_Jan_12_0041/2011_Jan_12_0041_1.ps
if ($fexist([fname]).eq.1) then
  shell fgrep -e "\(U=" [fname] >& tmp.txt
  shell $unquote('cat tmp.txt | sed "s/(1, \\\(U=/ /g" > tmp1.txt')
  shell $unquote('cat tmp1.txt | sed "s/ V\\\))/ /g" > tmp2.txt')
  ve/del tmp
  ve/read tmp tmp2.txt
  gl/cre uw $sigma(tmp(1))
else
  gl/cre uw 0
endif
ve/cre uw(1) r [uw]
return


macro phicutpl nc=1
ve/del ic,if1,if2,if3,if4,if5
ve/read ic,if1,if2,if3,if4,if5 phi_cuts.txt
ve/cre y1b(6) r
ve/cre y2b(6) r
ve/cre y1r(8) r
ve/cre y2r(8) r
ve/inp y1b(1) $sigma(if1([nc]))
ve/inp y1b(2) $sigma(if3([nc]))
ve/inp y1b(3) $sigma(if1([nc])-180)
ve/inp y1b(4) $sigma(if3([nc])-180)
ve/inp y1b(5) $sigma(if1([nc])+180)
ve/inp y1b(6) $sigma(if3([nc])+180)
ve/inp y2b(1) $sigma(if2([nc]))
ve/inp y2b(2) $sigma(if5([nc]))
ve/inp y2b(3) $sigma(if2([nc])-180)
ve/inp y2b(4) $sigma(if5([nc])-180)
ve/inp y2b(5) $sigma(if2([nc])+180)
ve/inp y2b(6) $sigma(if5([nc])+180)
*
nc1=[nc]+4
if ([nc1].gt.9) then
  nc1=[nc1]-9
endif
nc2=[nc]+5
if ([nc2].gt.9) then
  nc2=[nc2]-9
endif
ve/inp y1r(1) $sigma(if1([nc1])-180)
ve/inp y1r(2) $sigma(if3([nc1])-180)
ve/inp y1r(3) $sigma(if1([nc1])+180)
ve/inp y1r(4) $sigma(if3([nc1])+180)
ve/inp y1r(5) $sigma(if1([nc2])-180)
ve/inp y1r(6) $sigma(if3([nc2])-180)
ve/inp y1r(7) $sigma(if1([nc2])+180)
ve/inp y1r(8) $sigma(if3([nc2])+180)
*
ve/inp y2r(1) $sigma(if2([nc1])-180)
ve/inp y2r(2) $sigma(if5([nc1])-180)
ve/inp y2r(3) $sigma(if2([nc1])+180)
ve/inp y2r(4) $sigma(if5([nc1])+180)
ve/inp y2r(5) $sigma(if2([nc2])-180)
ve/inp y2r(6) $sigma(if5([nc2])-180)
ve/inp y2r(7) $sigma(if2([nc2])+180)
ve/inp y2r(8) $sigma(if5([nc2])+180)
*
exec seteps 
opt ngrid
op nstat
y1bx=$sigma(40*[nc]-80)
y2bx=$sigma(40*[nc]+40)
null -12 12 [y1bx] [y2bx]
*
x2=10
x1=-[x2]
*
set plci 4
set lwid 5
set fais 0
do i=1,6
  ymin=$sigma(max(y1b([i]),[y1bx]))
  ymax=$sigma(min(y2b([i]),[y2bx]))
  if ([ymin].lt.[ymax]) then
    dbox [x1] [x2] [ymin] [ymax]
  endif
enddo
*
x1=[x1]+0.2
x2=[x2]-0.2
set plci 2
set lwid 1
do i=1,8
  ymin=$sigma(max(y1r([i]),[y1bx]))
  ymax=$sigma(min(y2r([i]),[y2bx]))
  if ([ymin].lt.[ymax]) then
    dbox [x1] [x2] [ymin] [ymax]
  endif
enddo
*
set fais 3
set fasi 245
do i=1,6
  do j=1,8
    ymin=$sigma(max(y1b([i]),y1r([j])))
    ymax=$sigma(min(y2b([i]),y2r([j])))
    if ([ymin].lt.[ymax]) then
      dbox [x1] [x2] [ymin] [ymax]
    endif
  enddo
enddo
*
atitle 'z?R!, cm' '[f], degree'
exec save phicuts_counter[nc].eps f
return




macro goodspect scan=2011 nst=1
if ([scan].eq.2011) then
  i1=7700
  i2=11000
  d1=-150
  d2=200
  fname=mhad2011_goodspect_[nst].his
  gname=mhad2011_goodspect_[nst].txt
  nst1=[nst]-1
  fname0=mhad2011_goodspect_[nst1].his
  gname0=mhad2011_goodspect_[nst1].txt
endif
if ([scan].eq.2012) then
  i1=11000
  i2=14000
  d1=200
  d2=556
endif
do nc=1,9
  hi/del 10[nc]
  hi/del 20[nc]
enddo
if ($fexist([fname0]).eq.1) then
  hi/file 20 [fname0]
  do nc=1,9
    idh=10[nc]
    hi/copy //lun20/[idh] 20[nc]
  enddo
  close 20
endif
*
nmax=[i2]-[i1]+1
ve/cre gruns([nmax]) r
ve/cre prob(9,[nmax]) r
ind=0
do i=[i1],[i2]
  fnamei=run_[i]_spects_v2.his
  if ($fexist([fnamei]).eq.1) then
    exec mapcal#ndays [i]
    gl/imp ndays
    if (([ndays].gt.[d1]).and.([ndays].lt.[d2])) then
      ind=[ind]+1
      ve/inp gruns([ind]) [i]
      hi/file 20 [fnamei]
      mess [fnamei]
      do nc=1,9
        hi/del 100
        idh=10+[nc]
        hi/copy //lun20/[idh] 100
        if ([nst].gt.1) then
          hi/del 1001,1002
          hi/copy 20[nc] 1002
          exec hsigma @1001 = @100*($100 gt 0.2)
          exec hsigma @1002 = @1002*($1002 gt 0.2)
          hdiff=$call('mhdiff(1001,1002)')
        else
          hdiff=1
        endif
        mess [nc] [hdiff]
        ve/inp prob([nc],[ind]) [hdiff]
        mean=$hinfo(100,'mean')
        if (([mean].gt.1).and.([mean].lt.15).and.([hdiff].gt.0.9)) then
          idhc=10[nc]
          if ($hexist([idhc]).eq.0) then
            hi/copy 100 [idhc]
          else
            hi/op/add 100 [idhc] [idhc]
          endif
        endif
      enddo
      close 20
    endif
  endif
enddo
hi/file 20 [fname] ! N
do nc=1,9
  hrout 10[nc]
enddo
close 20
*
exec mapcal#vecut gruns
ve/write gruns,prob(1),prob(2),prob(3),prob(4),prob(5),prob(6),prob(7),prob(8),prob(9) [gname] 10f15.6
return



macro ettest ix=10 iy=10
nmax=10000
ve/cre runs([nmax]) r
ve/cre  et([nmax]) r
ve/cre det([nmax]) r
ve/read id1,id2 daycuts.txt 2f10.2
ind=0
tup=0
do i=6000,15000
  if ([tup].eq.1) then
    fnamei=run_[i]_spects_v2.his
  else
    fnamei=run_[i]_eventtime_v2.txt
  endif
  if ($fexist([fnamei]).eq.1) then
    exec mapcal#ndays [i]
    gl/imp ndays
      ind=[ind]+1
      ve/inp runs([ind]) [i]
      mess [fnamei]
      if ([tup].eq.1) then
        hi/file 20 [fnamei]
        hi/del 1000
        hi/copy //lun20/100 1000
        exec hsigma tx = @1000
        hi/del 1000
        hi/copy //lun20/200 1000
        exec hsigma t1 = @1000
        close 20
      else
        ve/del tx,t1
        ve/cre tx(10,20) r
        ve/cre t1(10,20) r
        ve/read tx,t1 [fnamei]
      endif
      nx=tx([ix],[iy])
      n1=t1([ix],[iy])
      if ([nx].ne.0) then
      eff=1-[n1]/[nx]
      deff=$sigma(sqrt([eff]*(1-[eff])/[nx]))
      ve/inp  et([ind])  [eff]
      ve/inp det([ind]) [deff]
      endif
  endif
enddo
*
exec mapcal#vecut runs
n=$vlen(runs)
exec mapcal#vecut et [n]
exec mapcal#vecut det [n]
sigma druns = runs*0
* exec $PER/s#vpl et det runs druns sz=0.05
return

macro ettest00 
fname=Eventtime/eventtime_00.tex
if ($fexist([fname]).eq.1) then
  shell rm [fname]
endif
for/file  20 [fname] new
close 20
ve/del runsp,drunsp,daysp,ddaysp,beamp,dbeamp,nevtp,meanp,dmeanp,rmsp,drmsp,wp,dwp,vp,dvp
ve/read runsp,drunsp,daysp,ddaysp,beamp,dbeamp,nevtp,meanp,dmeanp,rmsp,drmsp,wp,dwp,vp,dvp v2/dphi_vs_beams_new.txt 15e15.6
n=$vlen(runsp)
do i=1,[n]
  dirp=Eventtime/d[i]_e$sigma(beamp([i]))
  txt='null'
  fi=d[i]_e$sigma(beamp([i]))/eventtime_runs_t0_f1.eps
  dfi=Eventtime/[fi]
  if ($fexist([dfi]).eq.1) then
    txt=$unquote('     '){\includegraphics{[fi]}}\par}
  endif
  fi=d[i]_e$sigma(beamp([i]))/eventtime_runs_t0_f0.eps
  dfi=Eventtime/[fi]
  if ($fexist([dfi]).eq.1) then
    txt=$unquote('     '){\includegraphics{[fi]}}\par}
  endif
  if ([txt].ne.'null') then
    fmess '\begin{figure}[ht!b]' [fname]
    fmess '  \begin{minipage}{0.65\textwidth}' [fname]
    fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
    fmess [txt] [fname]
    txt=$unquote('     ')\caption{d[i]:e$sigma(beamp([i]))}
    fmess [txt] [fname]
    fmess '   \end{minipage}' [fname]
    fmess '\end{figure}' [fname]
  endif
  if ($sigma(mod([i],2)).eq.0) then
    fmess '\clearpage' [fname]
  endif
enddo
return


macro listtest0 
suffix=.col.gz
fo=absentruns.txt
if ($fexist([fo]).eq.1) then
  shell rm [fo]
endif
for/file  20 [fo] new
close 20
ns=0
*ns=[ns]+1; sc[ns]=MHAD2010; dir[ns]=
ns=[ns]+1; sc[ns]=MHAD2011; dir[ns]=/online/gridspool/proc/MHAD2011-4/
ns=[ns]+1; sc[ns]=MHAD2012; dir[ns]=/online/gridspool/proc/MHAD2012-2/
ns=[ns]+1; sc[ns]=OMEG2012; dir[ns]=/online/gridspool/proc/OMEG2012-0/
ns=[ns]+1; sc[ns]=PHI_2013; dir[ns]=/online/gridspool/proc/PHI_2013-0/
ns=[ns]+1; sc[ns]=RHO_2012; dir[ns]=/online/gridspool/proc/RHO_2012-0/
do i=1,[ns]
  shell TAKE_ALL_RUNS=1 .testrelease/.mainrelease/Offline/recodb.py points [sc[i]] > tmp.list
  shell $unquote('cat tmp.list | sed "s/\// \n/g" > out.txt')
  ve/del pnts
  ve/read pnts out.txt
  nl=$vlen(pnts)
  j=0
  while ([j].lt.[nl]) do
    j=[j]+1
    pnt=$format(pnts([j]),f6.1)
    if ([j].lt.[nl]) then
      j=[j]+1
      pnt0=pnts([j])
      if ([pnt0].lt.10) then
        pnt=$unquote([pnt])$unquote('/')[pnt0]
      else
        j=[j]-1
      endif
    endif
    shell TAKE_ALL_RUNS=1 .testrelease/.mainrelease/Offline/recodb.py list [sc[i]] [pnt] > tmp.list
    ve/del runs
    ve/read runs tmp.list
    do l=1,$vlen(runs)
      prt=0
      run=runs([l])
      if ([run].eq.13846) then
        prt=1
      endif
      txt=[sc[i]] $unquote([pnt]) [run]
      fn1=v2/run_[run]_spects_v2.his
      if ($fexist([fn1]).eq.0) then
        txt=$unquote([txt]): [fn1]
        prt=1
      endif
      fn2=[dir[i]]/exp$format([run],i8.8)[suffix]
      if ($fexist([fn2]).eq.0) then
        txt=$unquote([txt]) and [fn2] are not exist
        prt=1
      else
        txt=$unquote([txt]) is not exists
      endif
      if ([prt].eq.1) then
        fmess [txt] [fo]
      endif
    enddo
  endwhile
enddo
return


macro listtest
ve/del runsp,drunsp,daysp,ddaysp,beamp,dbeamp,nevtp,meanp,dmeanp,rmsp,drmsp,wp,dwp,vp,dvp
ve/read runsp,drunsp,daysp,ddaysp,beamp,dbeamp,nevtp,meanp,dmeanp,rmsp,drmsp,wp,dwp,vp,dvp v2/dphi_vs_beams_new.txt 15e15.6
ns=0
ns=[ns]+1; sc[ns]=MHAD2010
ns=[ns]+1; sc[ns]=MHAD2011
ns=[ns]+1; sc[ns]=MHAD2012
ns=[ns]+1; sc[ns]=OMEG2012
ns=[ns]+1; sc[ns]=PHI_2013
ns=[ns]+1; sc[ns]=RHO_2012
ve/cre tmp(1) r
do i=1,[ns]
  shell TAKE_ALL_RUNS=1 .testrelease/.mainrelease/Offline/recodb.py first [sc[i]] > tmp.list
  ve/read tmp tmp.list
  rf[i]=tmp(1)
  shell TAKE_ALL_RUNS=1 .testrelease/.mainrelease/Offline/recodb.py last [sc[i]] > tmp.list
  ve/read tmp tmp.list
  rl[i]=tmp(1)
  mess [sc[i]]={[rf[i]],[rl[i]]}
enddo
read x
n=$vlen(runsp)
ntot1=0
ntot2=0
ve/cre runst(30000) r
do i=1,[n]
  run=$sigma(int(runsp([i])))
  do j=1,[ns]
    if (([run].ge.[rf[j]]).and.([run].le.[rl[j]])) then
      scn=[sc[j]]
    endif
  enddo
  beam=beamp([i])
  shell TAKE_ALL_RUNS=1 .testrelease/.mainrelease/Offline/recodb.py list [scn] [beam] > tmp.list
  ve/del lruns
  ve/read lruns tmp.list
  do j=1,$vlen(lruns)
    runx=lruns([j])
    ve/inp runst([runx]) [runx]
  enddo
  run1=$sigma(runsp([i])-drunsp([i]))
  run2=$sigma(runsp([i])+drunsp([i]))
  mess [run1] [run2]
  ind=0
  do j=[run1],[run2]
    fnamei=v2/run_[j]_spects_v2.his
    if ($fexist([fnamei]).eq.1) then
      ind=[ind]+1
      ve/del tmp
      sigma tmp = abs(lruns-[j])
      sigma tmp1 = order(tmp,tmp)
      ntfd=tmp1(1)
      if ([ntfd].ne.0) then
        txt=point [i]: scan [scn]: beam [beam]: run [j] is not identified is recodb list
        mess [txt]
      else
*        txt=Runs [j] is identified
*        mess [txt]
      endif
    endif
  enddo
  mess [ind]($vlen(lruns)) files were tested
  ntot1=[ntot1]+[ind]
  ntot2=[ntot2]+$vlen(lruns)
*  read x
enddo
sigma runst=order(runst,-runst)
exec mapcal#vecut runst
mess total [ntot1] $vlen(runst) [ntot2]
return

macro sometest
ve/del runsp,drunsp,daysp,ddaysp,beamp,dbeamp,nevtp,meanp,dmeanp,rmsp,drmsp,wp,dwp,vp,dvp
ve/read runsp,drunsp,daysp,ddaysp,beamp,dbeamp,nevtp,meanp,dmeanp,rmsp,drmsp,wp,dwp,vp,dvp v2/dphi_vs_beams_new.txt 15e15.6
n=$vlen(runsp)
fname=Eventtime/eff.lt.0.1.tex
if ($fexist([fname]).eq.1) then
  shell rm [fname]
endif
for/file  20 [fname] new
close 20
ind=0
do i=1,[n]
  dirp=Eventtime/d[i]_e$sigma(beamp([i]))
  ve/del pars
  ve/read pars [dirp]/eventtime_runs_t0_f0.txt
  q=pars(15)
  prb=pars(18)
  dst=pars(19)
*
  if ((([q].gt.0.1).and.([q].lt.0.9)).and.([prb].gt.0.005).and.([dst].gt.25)) then
    txt='null'
    fi=d[i]_e$sigma(beamp([i]))/eventtime_runs_t0_f0.eps
    dfi=Eventtime/[fi]
    if ($fexist([dfi]).eq.1) then
      txt=$unquote('     '){\includegraphics{[fi]}}\par}
    endif
    if ([txt].ne.'null') then
      fmess '\begin{figure}[ht!b]' [fname]
      fmess '  \begin{minipage}{0.65\textwidth}' [fname]
      fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
      fmess [txt] [fname]
      txt=$unquote('     ')\caption{d[i]:e$sigma(beamp([i])): ix=0 iy=0}
      fmess [txt] [fname]
      fmess '   \end{minipage}' [fname]
      fmess '\end{figure}' [fname]
      ind=[ind]+1
      if ($sigma(mod([ind],2)).eq.0) then
        fmess '\clearpage' [fname]
      endif
    endif
  endif
  do ix=2,9
    do iy=1,20
      ve/del pars
      ve/read pars [dirp]/eventtime_runs_t[ix]_f[iy].txt
      q=pars(15)
      prb=pars(18)
      dst=pars(19)
*
      if ((([q].gt.0.1).and.([q].lt.0.9)).and.([prb].gt.0.005).and.([dst].gt.25)) then
        txt='null'
        fi=d[i]_e$sigma(beamp([i]))/eventtime_runs_t[ix]_f[iy].eps
        dfi=Eventtime/[fi]
        if ($fexist([dfi]).eq.1) then
          txt=$unquote('     '){\includegraphics{[fi]}}\par}
        endif
        if ([txt].ne.'null') then
          fmess '\begin{figure}[ht!b]' [fname]
          fmess '  \begin{minipage}{0.65\textwidth}' [fname]
          fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
          fmess [txt] [fname]
          txt=$unquote('     ')\caption{d[i]:e$sigma(beamp([i])): ix=[ix] iy=[iy]}
          fmess [txt] [fname]
          fmess '   \end{minipage}' [fname]
          fmess '\end{figure}' [fname]
          ind=[ind]+1
          if ($sigma(mod([ind],2)).eq.0) then
            fmess '\clearpage' [fname]
          endif
        endif
      endif
    enddo
  enddo
enddo
return

macro ettestxyq i1=1 i2=[n] opt=q
ve/del runsp,drunsp,daysp,ddaysp,beamp,dbeamp,nevtp,meanp,dmeanp,rmsp,drmsp,wp,dwp,vp,dvp
ve/read runsp,drunsp,daysp,ddaysp,beamp,dbeamp,nevtp,meanp,dmeanp,rmsp,drmsp,wp,dwp,vp,dvp v2/dphi_vs_beams_new.txt 15e15.6
n=$vlen(runsp)
i2=$sigma([i2])
fo=fails.txt
if ($fexist([fo]).eq.1) then
  shell rm [fo]
endif
for/file  20 [fo] new
close 20
do i=[i1],[i2]
  dirp=Eventtime/d[i]_e$sigma(beamp([i]))
  if ([opt].eq.'q') then
  fname=et_d[i].kumac
  if ($fexist([fname]).eq.1) then
    shell rm [fname]
  endif
  for/file  20 [fname] new
  close 20
*  txt=shell rm -v [dirp]/
*  fmess [txt] [fname]
  txt=gl/cre tl [i]
  fmess [txt] [fname]
  txt=exec mapcal#ettestxy 1 20 [i]
  fmess [txt] [fname]
  fnames=et_d[i].sh
  if ($fexist([fnames]).eq.1) then
    shell rm [fnames]
  endif
  for/file  20 [fnames] new
  close 20
  txt=#!/bin/bash
  fmess [txt] [fnames]
  txt=pwd 
  fmess [txt] [fnames]
  txt=cd [dirp] 
  fmess [txt] [fnames]
  txt=pwd 
  fmess [txt] [fnames]
  txt=ls ../eventtimei.tex
  fmess [txt] [fnames]
  txt=cp ../eventtimei.tex ./
  fmess [txt] [fnames]
  txt=latex eventtimei
  fmess [txt] [fnames]
  txt=latex eventtimei
  fmess [txt] [fnames]
  txt=dvips eventtimei
  fmess [txt] [fnames]
  shell chmod +x [fnames]
  txt=shell ./[fnames]
  fmess [txt] [fname]
  shell .testrelease/.mainrelease/Offline/submit.sh -q clusters,180 pawbigX11 -b [fname]
  else
    f=[dirp]/eventtime_runs_t1_f1.eps
    if ($fexist([f]).eq.1) then
      fnames=et_d[i].sh
      shell ./[fnames]
    else
      fmess [f] [fo]
    endif
  endif  
enddo
return



macro ettestxy iy1=1 iy2=20 tl=1
*ve/del beam,r1,r2
*ve/read beam,r1,r2 beamcuts.txt '3f10.2'
*ve/del runs,days,beams
*ve/read runs,days,beams runpars.txt 3f15.6
ve/del runsp,drunsp,daysp,ddaysp,beamp,dbeamp,nevtp,meanp,dmeanp,rmsp,drmsp,wp,dwp,vp,dvp
ve/read runsp,drunsp,daysp,ddaysp,beamp,dbeamp,nevtp,meanp,dmeanp,rmsp,drmsp,wp,dwp,vp,dvp v2/dphi_vs_beams_new.txt 15e15.6
imin=$sigma(runsp([tl])-drunsp([tl]))
imax=$sigma(runsp([tl])+drunsp([tl]))
gl/cre dirp Eventtime/d[tl]_e$sigma(beamp([tl]))
if ($fexist([dirp]).eq.0) then
  shell mkdir [dirp]
endif
shell rm -v [dirp]/*
mess [dirp]
do iy=[iy1],[iy2]
  exec mapcal#ettestx [iy] [imin] [imax]
*  exec mapcal#ettestxd [iy]
enddo
exec mapcal#ettestabl 0 v2 [imin] [imax]
return



macro ettestx iy=10 imin=6000 imax=15000
nmax=10000
ve/cre runs([nmax]) r
do ix=0,10
  ve/cre  et[ix]([nmax]) r
  ve/cre det[ix]([nmax]) r
enddo
ve/cre nx0([nmax]) r
dir=v2
if ([dir].eq.'v1') then
  suffh=.his
  sufft=.txt
else
  suffh=_[dir].his
  sufft=_[dir].txt
endif
dire=Eventtime
dire=[dir]
ind=0
tup=0
do i=[imin],[imax]
  if ([tup].eq.1) then
    fnamei=[dir]/run_[i]_spects[suffh]
  else
    fnamei=[dir]/run_[i]_eventtime[sufft]
  endif
  if ($fexist([fnamei]).eq.1) then
    exec mapcal#ndays [i]
    gl/imp ndays
      ind=[ind]+1
      ve/inp runs([ind]) [i]
      mess [fnamei]
      if ([tup].eq.1) then
        hi/file 20 [fnamei]
        hi/del 1000
        hi/copy //lun20/100 1000
        exec hsigma tx = @1000
        hi/del 1000
        hi/copy //lun20/200 1000
        exec hsigma t1 = @1000
        close 20
      else
        ve/del txa,t1a,txb,t1b
        ve/cre txa(10,20) r
        ve/cre t1a(10,20) r
        ve/cre txb(10,20) r
        ve/cre t1b(10,20) r
        ve/read txa,t1a,txb,t1b [fnamei]
        ve/del tx,t1
        sigma tx = txa + txb
        sigma t1 = t1a + t1b
      endif
      do ix=1,10
        nx=tx([ix],[iy])
        n1=t1([ix],[iy])
        if ([nx].ne.0) then
          eff=1-[n1]/[nx]
          deff=$sigma(sqrt([eff]*(1-[eff])/[nx]))
          ve/inp  et[ix]([ind])  [eff]
          ve/inp det[ix]([ind]) [deff]
        endif
      enddo
      nx=$sigma(vsum(tx))
      n1=$sigma(vsum(t1))
      if ([nx].ne.0) then
        eff=1-[n1]/[nx]
        deff=$sigma(sqrt([eff]*(1-[eff])/[nx]))
        ve/inp  et0([ind])  [eff]
        ve/inp det0([ind]) [deff]
        ve/inp  nx0([ind]) [nx]
      else
        ind=[ind]-1
      endif
  endif
enddo
*
exec mapcal#vecut runs
sigma druns = runs*0
n=$vlen(runs)
gl/imp dirp
point 1000
exec mapcal#vecut nx0 [n]
ve/cre  etx(1000) r
ve/cre detx(1000) r
ve/cre res(8) r
ve/cre nx(1) i [n]
do ix=0,10
  exec mapcal#vecut et[ix] [n]
  exec mapcal#vecut det[ix] [n]
*
  ve/del etxt,detxt
  ve/copy  et[ix](1:[n])  etxt
  ve/copy det[ix](1:[n]) detxt
  amean=$sigma(vsum(etxt/detxt**2)/vsum(1/detxt**2))
  mean=$sigma(vsum(etxt/detxt)/vsum(1/detxt))
  sigma svetx = sumv((etxt-[mean])/detxt)
  ve/draw svetx
  mess [mean]
*  read x
  dst=$sigma((vmax(svetx)-vmin(svetx)))
  sigma detxt=order(detxt,etxt)
  sigma etxt=order(etxt,etxt)
  ve/copy  etxt(1:[n])  etx(1:[n])
  ve/copy detxt(1:[n]) detx(1:[n])
  xxx=$call('twodelta.f(nx)')
*    a=$sigma(vsum(et[ix])/[n])
*    b=$sigma(vsum(et[ix]**2)/[n])
*    c=$sigma(vsum(et[ix]**3)/[n])
*    d=0
*    d=$sigma([d]-3*([a]*[b])**2)
*    d=$sigma([d]+4*[b]**3)
*    d=$sigma([d]+4*[a]**3*[c])
*    d=$sigma([d]-6*[a]*[b]*[c])
*    d=$sigma([d]+[c]**2)
*    if ([d].gt.0) then
*      d=$sigma(sqrt([d]))
*    else
*      d=0
*    endif
*    dd=$sigma(2*([a]**2-[b]))
*    y1=$sigma([a]*[b]-[c]+[d])
*    y1=$sigma([y1]/([dd]))
*    y2=$sigma([a]*[b]-[c]-[d])
*    y2=$sigma([y2]/([dd]))
*    q=$sigma(2*[a]**3-3*[a]*[b]+[c]+[d])
*    q=$sigma([q]/2/[d])
    q=res(1)
    y1=res(2)
    y2=res(3)
    a=res(4)
    b=res(5)
    c=res(6)
    d=res(7)
    g=res(8)
    a=[amean]
    sigma chi=((et[ix]-[a])/det[ix])**2
**    sigma chi=((et[ix]-[y1])/det[ix])**2+((et[ix]-[y2])/det[ix])**2
    ve/cre ch(1) r $sigma(vsum(chi))
    ve/cre ndf(1) i [n] 
    ns=$sigma(sqrt(vsum(chi)/[n]/2-1))
    prb=$call('prob(ch,ndf)')
*    s21=0
*    s21=$sigma([s21]+[c]-3*[b]*[y2])
*    s21=$sigma([s21]+2*[y2]**3-[y1]**3*[q])
*    s21=$sigma([s21]+3*[y1]**2*[y2]*[q])
*    s21=$sigma([s21]-2*[y2]**3*[q])
*    s21=$sigma([s21]/([y1]-[y2])/3/[q])
*    s22=0
*    s22=$sigma([s22]+[c]-3*[b]*[y1])
*    s22=$sigma([s22]-[y2]**3*(1-[q]))
*    s22=$sigma([s22]+3*[y2]**2*[y1]*(1-[q]))
*    s22=$sigma([s22]+2*[y1]**3*[q])
*    s22=$sigma([s22]/([y1]-[y2])/3/(1-[q]))
*    mess s21=[s21] s22=[s22]
*
  if ([ix].eq.0) then
    exec mapcal#step x=runs y=et[ix] dy=det[ix] rf=0    
    rf=p3(3)
    ng=0
    nb=0
    do i=1,[n]
      if ($sigma(p3(1)).lt.$sigma(p3(2))) then
        pb=2
        if ($sigma(runs([i])).lt.[rf]) then
          ng=[ng]+nx0([i])
        else
          nb=[nb]+nx0([i])
        endif
      else
        pb=1
        if ($sigma(runs([i])).ge.[rf]) then
          ng=[ng]+nx0([i])
        else
          nb=[ng]+nx0([i])
        endif
      endif
    enddo
    eb=[nb]+[ng]
    if ($sigma([eb]).ne.0) then
      eb=[nb]/([nb]+[ng])
    else
      eb=0
    endif
  else
    exec mapcal#step x=runs y=et[ix] dy=det[ix] rf=[rf]
  endif
  l=runs(1)
  r=runs([n])
  mess [ix] [iy] [pb] [eb] [l] [r] [q] [y1] [y2] [prb]
  if ([ix].ne.0) then
    ve/cre pars(20) r [ix] [iy] [pb] [eb] [l] [r] _
                      $sigma(p3(1)) $sigma(p3(2)) $sigma(p3(3)) _
                      $sigma(dp3(1)) $sigma(dp3(2)) $sigma(dp3(3)) _
                      $sigma(chi2(1)) $sigma(chi2(2)) _
                      [q] [y1] [y2] [prb] [dst] [g]
    ve/write pars [dirp]/eventtime_runs_t[ix]_f[iy].txt
  else
    ve/cre pars(20) r 0 0 [pb] [eb] [l] [r] _
                      $sigma(p3(1)) $sigma(p3(2)) $sigma(p3(3)) _
                      $sigma(dp3(1)) $sigma(dp3(2)) $sigma(dp3(3)) _
                      $sigma(chi2(1)) $sigma(chi2(2)) _
                      [q] [y1] [y2] [prb] [dst] [g]
    ve/write pars [dirp]/eventtime_runs_t0_f0.txt
  endif
*  * exec $PER/s#vpl et[ix] det[ix] runs druns sz=0.05
  exec vpl#pl0 et[ix] det[ix] runs druns sz=0.1
  wl=$GRAFINFO('WNXMIN')
  wr=$GRAFINFO('WNXMAX')
  wd=$GRAFINFO('WNYMIN')
  wu=$GRAFINFO('WNYMAX')
  wu=1.3*[wu]
  null [wl] [wr] [wd] [wu]
  exec vpl#pl0 et[ix] det[ix] runs druns sz=0.1 o=s
*  exec vpl#pl et[ix] det[ix] runs druns sz=0.05
  atitle 'runs' '[s]'
  mess $sigma((p3(1)-p3(2))/sqrt(dp3(1)**2+dp3(2)**2))
  l=runs(1)
  r=runs([n])
  fun/pl stepp.f [l] [r] s
*  read x
*  l=$GRAFINFO('WNXMIN')
*  r=$GRAFINFO('WNXMAX')
*  d=$GRAFINFO('WNYMIN')
*  u=$GRAFINFO('WNYMAX')
*  do i=1,-$vlen(r1)
*    c=r1([i])
*    if (([c].ge.[l]).and.([c].le.[r])) then
*      line [c] [d] [c] [u]
*    endif
*  enddo
  if ([y1].lt.[y2]) then
    ym=[y2]
    yl=[y1]
    p=1-[q]
  else
    ym=[y1]
    yl=[y2]
    p=[q]
  endif
  txt=prob = $sigma(int([prb]*10000+0.5)/10000)
  exec pl#tf 0.05 0.9 [txt]
  txt=n?[s]! = $sigma(int([g]*10000+0.5)/10000)
  exec pl#tf 0.05 0.8 [txt]
  txt=[e] = $sigma(int([p]*10000+0.5)/10000)
  exec pl#tf 0.35 0.9 [txt]
  txt=[y] = $sigma(int([dst]*10000+0.5)/10000)
  exec pl#tf 0.35 0.8 [txt]
  txt=l?max! = $sigma(int([ym]*10000+0.5)/10000)
  exec pl#tf 0.65 0.9 [txt]
  txt=l?mid! = $sigma(int([a]*10000+0.5)/10000)
  exec pl#tf 0.65 0.8 [txt]
  if ([ix].ne.0) then
    exec save [dirp]/eventtime_runs_t[ix]_f[iy].eps f
  else
    line [l] [y1] [r] [y1]
    line [l] [y2] [r] [y2]
    exec save [dirp]/eventtime_runs_t0_f0.eps f
  endif
  mess [ix] [iy] q=[q] y1=[y1] y2=[y2] prob=[prb] ns=[ns] d=[d] a=[a] b=[b] c=[c]
*  read x
enddo
return


macro step x=runs y=et0 dy=det0 rf=0
n=$vlen([x])
l=[x](1)
r=[x]([n])
mess [l] [r] [rf]
ve/cre p3(3) r
ve/cre dp3(3) r
npar=3
ve/cre chi2(2) r
ve/cre paru([npar]) r
ve/cre dparu([npar]) r
ve/cre covu([npar],[npar]) r
if (([n].gt.10).and.($sigma(vsum([y])).ne.0)) then
  sigma xo = order([x],[y])
  sigma yo = order([y],[y])
  i1=1
  i2=$sigma(int([n]/10))
  ve/del xol,yol
  ve/copy xo([i1]:[i2]) xol
  ve/copy yo([i1]:[i2]) yol
  xolm=$sigma(vsum(xol)/([i2]-[i1]+1))
  yolm=$sigma(vsum(yol)/([i2]-[i1]+1))
  i1=$sigma(int(9*[n]/10))
  i2=[n]
  ve/del xor,yor
  ve/copy xo([i1]:[i2]) xor
  ve/copy yo([i1]:[i2]) yor
  xorm=$sigma(vsum(xor)/([i2]-[i1]+1))
  yorm=$sigma(vsum(yor)/([i2]-[i1]+1))
  s=$sigma([xolm]+[xorm]-([l]+[r])/2)
  if ([rf].eq.0) then
    ve/cre s3(3) r 0.001 0.001 0.1
  else
    s=[rf]
    ve/cre s3(3) r 0.001 0.001 0
  endif
  if ([xolm].lt.[xorm]) then
    ve/cre p3(3) r [yolm] [yorm] [s]
  else
    ve/cre p3(3) r [yolm] [yorm] [s]
  endif
  ve/cre dp3(3) r
  ve/cre st(1) r 1
  do i=1,3
    ve/fit [x] [y] [dy] step.f sb 3 p3 s3 ! ! dp3
  enddo
  ve/cre st(1) r 0.3
  do i=1,3
    ve/fit [x] [y] [dy] step.f sb 3 p3 s3 ! ! dp3
  enddo
  call covm.f(1)
  call covmpen(chi2,[npar],paru,dparu)
  call covmcov([npar],covu)
endif
sigma p3 = abs(p3)
return



macro ttt
nmax=100
ve/cre beam([nmax]) r
ve/cre r1([nmax]) r
ve/cre r2([nmax]) r
ve/cre par(100) r
ind=0
do i=1,100
exec mapcal#readspectn /work/users/konctbel/snd2k/R005-999/exp1012.tbl [i]
ibeam=par(1)
if ([ibeam].ne.0) then
  ind=[ind]+1
  ve/inp beam([ind]) [ibeam]
  exec mapcal#vecut par
  n=$vlen(par)
  ve/inp r1([ind]) $sigma(par(2))
  ve/inp r2([ind]) $sigma(par([n]))
endif
enddo
exec mapcal#vecut beam
exec mapcal#vecut r1 $vlen(beam)
exec mapcal#vecut r2 $vlen(beam)
ve/write beam,r1,r2 beamcuts.txt '3f10.2'
return

macro readspectn fname=x line=1
s='"'
l='|'
*ve/del par
sigma par = par*0
if ($fexist([fname]).eq.1) then
  txt=cat [fname] [l] sed -n [s][line]p[s] > out.txt
  shell $unquote([txt])
  shell $unquote('cat out.txt | sed "s/\t/ /g" > out1.txt')
  shell $unquote('cat out1.txt | sed "s/ /\n/g" > out2.txt')
  ve/read par out2.txt
endif
return

macro ettestxd iy=10 imin=6000 imax=15000
*ve/del id1,id2
*ve/read id1,id2 daycuts.txt 2f10.2
*nmax=$vlen(id1)
nmax=10000
ve/cre runs([nmax]) r
ve/cre druns([nmax]) r
do ix=1,10
  ve/cre  et[ix]([nmax]) r
  ve/cre det[ix]([nmax]) r
enddo
dir=v2
if ([dir].eq.'v1') then
  suffh=.his
  sufft=.txt
else
  suffh=_[dir].his
  sufft=_[dir].txt
endif
dire=Eventtime
dire=[dir]
ind=0
tup=0
do d=1,[nmax]
ve/del tx,t1
ve/cre tx(10,20) r
ve/cre t1(10,20) r
i1=id1([d])
i2=id2([d])
do i=[i1],[i2]
  if ([tup].eq.1) then
    fnamei=[dir]/run_[i]_spects[suffh]
  else
    fnamei=[dir]/run_[i]_eventtime[sufft]
  endif
  if ($fexist([fnamei]).eq.1) then
    exec mapcal#ndays [i]
    gl/imp ndays
      mess [fnamei]
      if ([tup].eq.1) then
        hi/file 20 [fnamei]
        hi/del 1000
        hi/copy //lun20/100 1000
        exec hsigma tx = @1000
        hi/del 1000
        hi/copy //lun20/200 1000
        exec hsigma t1 = @1000
        close 20
      else
        ve/del txa,t1a,txb,t1b
        ve/cre txa(10,20) r
        ve/cre t1a(10,20) r
        ve/cre txb(10,20) r
        ve/cre t1b(10,20) r
        ve/read txa,t1a,txb,t1b [fnamei]
        sigma tx = tx + txa + txb
        sigma t1 = t1 + t1a + t1b
      endif
  endif
enddo
      ind=[ind]+1
      ve/inp  runs([ind]) $sigma(([i2]+[i1])/2)
      ve/inp druns([ind]) $sigma(([i2]-[i1])/2)
      do ix=1,10
        nx=tx([ix],[iy])
        n1=t1([ix],[iy])
        if ([nx].ne.0) then
          eff=1-[n1]/[nx]
          deff=$sigma(sqrt([eff]*(1-[eff])/[nx]))
          ve/inp  et[ix]([ind])  [eff]
          ve/inp det[ix]([ind]) [deff]
        endif
      enddo
enddo
*
do ix=1,10
  * exec $PER/s#vpl et[ix] det[ix] runs druns sz=0.05
  atitle 'runs' '[s]'
  exec save eventtime_days_t[ix]_f[iy].eps f 
enddo
return



macro ettestabl d=1 dir=v2 i1=0 i2=0
*goto 2
ve/del id1,id2
ve/read id1,id2 daycutsx.txt 2f10.2
ind=0
tup=0
*
ve/del txa,t1a,txb,t1b
ve/cre txa(10,20) r
ve/cre t1a(10,20) r
ve/cre txb(10,20) r
ve/cre t1b(10,20) r
if ([i1].eq.0) then
  i1=id1([d])
endif
if ([i2].eq.0) then
  i2=id2([d])
endif
if ([dir].eq.'v1') then
  suffh=.his
  sufft=.txt
else
  suffh=_[dir].his
  sufft=_[dir].txt
endif
dire=Eventtime
dire=[dir]
do i=[i1],[i2]
  if ([tup].eq.1) then
    fnamei=[dir]/run_[i]_spects[suffh]
  else
    fnamei=[dir]/run_[i]_eventtime[sufft]
  endif
  if ($fexist([fnamei]).eq.1) then
    exec mapcal#ndays [i]
    gl/imp ndays
      mess [fnamei]
      if ([tup].eq.1) then
        hi/file 20 [fnamei]
        hi/del 1000
        hi/copy //lun20/100 1000
        exec hsigma tx = @1000
        hi/del 1000
        hi/copy //lun20/200 1000
        exec hsigma t1 = @1000
        close 20
      else
        ve/del txai,t1ai,txbi,t1bi
        ve/cre txai(10,20) r
        ve/cre t1ai(10,20) r
        ve/cre txbi(10,20) r
        ve/cre t1bi(10,20) r
        ve/read txai,t1ai,txbi,t1bi [fnamei]
        sigma txa = txa + txai
        sigma t1a = t1a + t1ai
        sigma txb = txb + txbi
        sigma t1b = t1b + t1bi
      endif
  endif
enddo
*          eff=1-[n1]/[nx]
*          deff=$sigma(sqrt([eff]*(1-[eff])/[nx]))
*
*
sigma tx = txa + txb
sigma t1 = t1a + t1b
1:
gl/imp dirp
*
exec seteps 0
set lwid 1
opt ngrid
ve/cre txr(10,20) r 200*1
ve/cre t1r(10,20) r 200*1
do i=1,20
  ve/copy tx(2:9,[i]) txr(2:9,[i])
  ve/copy t1(2:9,[i]) t1r(2:9,[i])
enddo
2d 100 ! 10 -18 162 20 0 360
2d 200 ! 10 -18 162 20 0 360
hi/put/cont 100 txr
hi/put/cont 200 t1r
*exec hsigma @300 = 1 - @200/@100
*hi/pl 300 box
atitle '[q], degree' '[f], degree'
*exec save [dirp]/eventtime_table_[i1]_[i2]_v0.eps f
*
2:
*exec hsigma rt = @300
exec hsigma rt = 1 - t1/tx
null 18 162 0 360
do i=2,9
  do j=1,20
    dx=rt([i],[j])
    dx=$sigma(sqrt([dx])*9)
    xc=$sigma(([i]-0.5)*18)
    x1=$sigma([xc]-[dx])
    x2=$sigma([xc]+[dx])
    yc=$sigma(([j]-0.5)*18)
    y1=$sigma([yc]-[dx])
    y2=$sigma([yc]+[dx])
    box [x1] [x2] [y1] [y2]
  enddo
enddo
set basl 0.01
set ltyp 15
*set lwid 0
do i=2,8
  line $sigma(18*[i]) 0 $sigma(18*[i]) 360
enddo
do j=1,19
  line 18 $sigma(18*[j]) 162 $sigma(18*[j])
enddo
set ltyp 1
atitle '[q], degree' '[f], degree'
exec save [dirp]/eventtime_table_[i1]_[i2]_v1.eps f
*
exec hsigma rt = 1 - t1a/txa
null 18 162 0 360
do i=2,9
  do j=1,20
    dx=rt([i],[j])
    dx=$sigma(sqrt([dx])*9)
    xc=$sigma(([i]-0.5)*18)
    x1=$sigma([xc]-[dx])
    x2=$sigma([xc]+[dx])
    yc=$sigma(([j]-0.5)*18)
    y1=$sigma([yc]-[dx])
    y2=$sigma([yc]+[dx])
    box [x1] [x2] [y1] [y2]
  enddo
enddo
set basl 0.01
set ltyp 15
*set lwid 0
do i=2,8
  line $sigma(18*[i]) 0 $sigma(18*[i]) 360
enddo
do j=1,19
  line 18 $sigma(18*[j]) 162 $sigma(18*[j])
enddo
set ltyp 1
atitle '[q], degree' '[f], degree'
*exec save [dirp]/eventtime_table_[i1]_[i2]_v1a.eps f
*
exec hsigma rt = 1 - t1b/txb
null 18 162 0 360
do i=2,9
  do j=1,20
    dx=rt([i],[j])
    dx=$sigma(sqrt([dx])*9)
    xc=$sigma(([i]-0.5)*18)
    x1=$sigma([xc]-[dx])
    x2=$sigma([xc]+[dx])
    yc=$sigma(([j]-0.5)*18)
    y1=$sigma([yc]-[dx])
    y2=$sigma([yc]+[dx])
    box [x1] [x2] [y1] [y2]
  enddo
enddo
set basl 0.01
set ltyp 15
*set lwid 0
do i=2,8
  line $sigma(18*[i]) 0 $sigma(18*[i]) 360
enddo
do j=1,19
  line 18 $sigma(18*[j]) 162 $sigma(18*[j])
enddo
set ltyp 1
atitle '[q], degree' '[f], degree'
*exec save [dirp]/eventtime_table_[i1]_[i2]_v1b.eps f
*
*exec hsigma rt = @300
exec hsigma rt = 1 - t1/tx
set txal 23
set chhe 0.1
null 18 162 0 360
do i=2,9
  do j=1,20
    dx=rt([i],[j])
    xc=$sigma(([i]-0.5)*18)
    yc=$sigma(([j]-0.5)*18)
    txt=$sigma(int([dx]*1000+0.5)/10) %
    if ([dx].lt.0.05) then
      set txfp -10
    else
      set txfp -20
    endif
    itx [xc] [yc] [txt]
  enddo
enddo
set txfp -10
set basl 0.01
set ltyp 15
*set lwid 0.1
do i=2,8
  line $sigma(18*[i]) 0 $sigma(18*[i]) 360
enddo
do j=1,19
  line 18 $sigma(18*[j]) 162 $sigma(18*[j])
enddo
set ltyp 1
atitle '[q], degree' '[f], degree'
exec save [dirp]/eventtime_table_[i1]_[i2]_v2.eps f
*
exec hsigma rt = 1 - t1a/txa
set txal 23
set chhe 0.1
null 18 162 0 360
do i=2,9
  do j=1,20
    dx=rt([i],[j])
    xc=$sigma(([i]-0.5)*18)
    yc=$sigma(([j]-0.5)*18)
    txt=$sigma(int([dx]*1000+0.5)/10) %
    if ([dx].lt.0.05) then
      set txfp -10
    else
      set txfp -20
    endif
    itx [xc] [yc] [txt]
  enddo
enddo
set txfp -10
set basl 0.01
set ltyp 15
*set lwid 0.1
do i=2,8
  line $sigma(18*[i]) 0 $sigma(18*[i]) 360
enddo
do j=1,19
  line 18 $sigma(18*[j]) 162 $sigma(18*[j])
enddo
set ltyp 1
atitle '[q], degree' '[f], degree'
*exec save [dirp]/eventtime_table_[i1]_[i2]_v2a.eps f
*
exec hsigma rt = 1 - t1b/txb
set txal 23
set chhe 0.1
null 18 162 0 360
do i=2,9
  do j=1,20
    dx=rt([i],[j])
    xc=$sigma(([i]-0.5)*18)
    yc=$sigma(([j]-0.5)*18)
    txt=$sigma(int([dx]*1000+0.5)/10) %
    if ([dx].lt.0.05) then
      set txfp -10
    else
      set txfp -20
    endif
    itx [xc] [yc] [txt]
  enddo
enddo
set txfp -10
set basl 0.01
set ltyp 15
*set lwid 0.1
do i=2,8
  line $sigma(18*[i]) 0 $sigma(18*[i]) 360
enddo
do j=1,19
  line 18 $sigma(18*[j]) 162 $sigma(18*[j])
enddo
set ltyp 1
atitle '[q], degree' '[f], degree'
*exec save [dirp]/eventtime_table_[i1]_[i2]_v2b.eps f
*
exec hsigma rt = 1 - t1/tx
exec hsigma drt = sqrt(rt*(1-rt)/tx)
ve/write rt,drt [dirp]/eventtime_table_[i1]_[i2].txt 2f15.6
*
ve/del rtm,drtm
ve/copy rt rtm
ve/copy drt drtm
ve/del pars
ve/read pars [dirp]/eventtime_runs_t0_f0.txt
eb=pars(4)
q=pars(15)
y1=pars(16)
y2=pars(17)
if ([y1].gt.[y2]) then
  ym=[y1]
else
  ym=[y2]
  q=1-[q]
endif
*if ([eb].gt.0.1) then
q=1
if ([q].gt.0.1) then
  do ix=1,10
    do iy=1,20
      ve/del pars
      ve/read pars [dirp]/eventtime_runs_t[ix]_f[iy].txt
      q=pars(15)
      prb=pars(18)
      dst=pars(19)
      if (([q].gt.0.1).and.([dst].gt.25)) then
        y1=pars(16)
        y2=pars(17)
        if ([y1].gt.[y2]) then
          ym=[y1]
        else
          ym=[y2]
        endif
        ve/inp rtm([ix],[iy]) [ym]
        ve/inp drtm([ix],[iy]) 0
      endif
*      lb=$sigma(abs(pars(7)-pars(8))/sqrt(pars(10)**2+pars(11)**2))
*      mess [ix] [iy] $sigma(pars(7)) $sigma(pars(8)) $sigma(pars(10)) $sigma(pars(11)) [lb]
*      if ([lb].gt.3) then
*        ib=pars(3)+6
*        ve/inp rtm([ix],[iy]) $sigma(pars([ib]))
*        ib=pars(3)+9
*        ve/inp drtm([ix],[iy]) $sigma(pars([ib]))
*      endif
    enddo
  enddo
endif
ve/write rtm,drtm [dirp]/eventtime_table_max_[i1]_[i2].txt 2f15.6
*
null 18 162 0 360
do i=2,9
  do j=1,20
    dx=rtm([i],[j])
    dx=$sigma(sqrt([dx])*9)
    xc=$sigma(([i]-0.5)*18)
    x1=$sigma([xc]-[dx])
    x2=$sigma([xc]+[dx])
    yc=$sigma(([j]-0.5)*18)
    y1=$sigma([yc]-[dx])
    y2=$sigma([yc]+[dx])
    box [x1] [x2] [y1] [y2]
  enddo
enddo
set basl 0.01
set ltyp 15
*set lwid 0
do i=2,8
  line $sigma(18*[i]) 0 $sigma(18*[i]) 360
enddo
do j=1,19
  line 18 $sigma(18*[j]) 162 $sigma(18*[j])
enddo
set ltyp 1
atitle '[q], degree' '[f], degree'
exec save [dirp]/eventtime_table_max_[i1]_[i2]_v1.eps f
*
set txal 23
set chhe 0.1
null 18 162 0 360
do i=2,9
  do j=1,20
    dx=rtm([i],[j])
    xc=$sigma(([i]-0.5)*18)
    yc=$sigma(([j]-0.5)*18)
    txt=$sigma(int([dx]*1000+0.5)/10) %
    if ([dx].lt.0.05) then
      set txfp -10
    else
      set txfp -20
    endif
    itx [xc] [yc] [txt]
  enddo
enddo
set txfp -10
set basl 0.01
set ltyp 15
*set lwid 0.1
do i=2,8
  line $sigma(18*[i]) 0 $sigma(18*[i]) 360
enddo
do j=1,19
  line 18 $sigma(18*[j]) 162 $sigma(18*[j])
enddo
set ltyp 1
atitle '[q], degree' '[f], degree'
exec save [dirp]/eventtime_table_max_[i1]_[i2]_v2.eps f
*
fname=[dirp]/eventtime_table.tex
if ($fexist([fname]).eq.1) then
  shell rm [fname]
endif
for/file  20 [fname] new
close 20
fmess '\begin{figure}[ht!b]' [fname]
fmess '  \begin{minipage}{0.65\textwidth}' [fname]
fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
txt=$unquote('     '){\includegraphics{eventtime_table_[i1]_[i2]_v1.eps}}\par}
fmess [txt] [fname]
txt=$unquote('     ')\caption{calorimeter map(avg) runs=([i1],[i2])}
fmess [txt] [fname]
fmess '   \end{minipage}' [fname]
fmess '\end{figure}' [fname]
fmess '\begin{figure}[ht!b]' [fname]
fmess '  \begin{minipage}{0.65\textwidth}' [fname]
fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
txt=$unquote('     '){\includegraphics{eventtime_table_[i1]_[i2]_v2.eps}}\par}
fmess [txt] [fname]
txt=$unquote('     ')\caption{calorimeter map(avg) runs=([i1],[i2])}
fmess [txt] [fname]
fmess '   \end{minipage}' [fname]
fmess '\end{figure}' [fname]
fmess '\begin{figure}[ht!b]' [fname]
fmess '  \begin{minipage}{0.65\textwidth}' [fname]
fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
txt=$unquote('     '){\includegraphics{eventtime_table_max_[i1]_[i2]_v1.eps}}\par}
fmess [txt] [fname]
txt=$unquote('     ')\caption{calorimeter map(max) runs=([i1],[i2])}
fmess [txt] [fname]
fmess '   \end{minipage}' [fname]
fmess '\end{figure}' [fname]
fmess '\begin{figure}[ht!b]' [fname]
fmess '  \begin{minipage}{0.65\textwidth}' [fname]
fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
txt=$unquote('     '){\includegraphics{eventtime_table_max_[i1]_[i2]_v2.eps}}\par}
fmess [txt] [fname]
txt=$unquote('     ')\caption{calorimeter map(max) runs=([i1],[i2])}
fmess [txt] [fname]
fmess '   \end{minipage}' [fname]
fmess '\end{figure}' [fname]
return

macro daycut
ve/cre id1(10000) r
ve/cre id2(10000) r
dr=$sigma(-300+8.5/24)
ind=0
do i=5000,20000
  fname=v2/run_[i]_spects_v2.his
  if ($fexist([fname]).eq.1) then
    exec mapcal#ndays [i]
    gl/imp ndays
    if ([ndays].gt.[dr]) then
      dr=$sigma(int([ndays])+8.5/24)
      if ([ndays].gt.[dr]) then
        dr=[dr]+1
      endif
      if ([ind].gt.0) then
        ve/inp id2([ind]) [io]
        ind=[ind]+1
        ve/inp id1([ind]) [i]
      else
        ind=[ind]+1
        ve/inp id1([ind]) [i]
      endif
    endif
    io=[i]
  endif
enddo
ve/inp id2([ind]) [io]
exec mapcal#vecut id1
exec mapcal#vecut id2
ve/write id1,id2 daycuts.txt 2f10.2
return


macro filetest exp=mhad2011
n=0
if ([exp].eq.'mhad2011') then
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_1000_p1.fwi          
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_1000_p2.fwi          
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_1000_p3.fwi          
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_508.7.fwi            
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_509.8.fwi            
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_510.5-2_p1.fwi       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_510.5-2_p2.fwi       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_510.5-2_p3.fwi       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_510.5-2_p4.fwi       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_525_p1.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_525_p2.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_525_p3.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_525_p4.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_537.5_p1.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_537.5_p2.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_537.5_p3.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_550_p1.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_550_p2.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_550_p3.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_562.5_p1.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_562.5_p2.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_575_p1.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_575_p2.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_587.5_p1.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_587.5_p2.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_600_p1.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_600_p2.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_612.5_p1.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_612.5_p2.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_625_p1.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_625_p2.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_637.5_p1.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_637.5_p2.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_650_9000.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_650_p1.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_650_p2.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_650_p3.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_662.5_p1.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_662.5_p2.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_675_p1.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_675_p2.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_687.5_p1.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_687.5_p2.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_700_p1.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_700_p2.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_712.5_p1.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_712.5_p2.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_725_p1.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_725_p2.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_737.5_p1.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_737.5_p2.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_750-1_p1.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_750-1_p2.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_750_p1.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_750_p2.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_762.5_p1.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_762.5_p2.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_775_p1.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_775_p2.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_787.5_p1.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_787.5_p2.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_800_p1.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_800_p2.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_812.5.fwi            
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_825_p1.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_825_p2.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_825_p3.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_837.5.fwi            
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_850_p1.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_850_p2.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_862.5.fwi            
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_875_p1.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_875_p2.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_887.5.fwi            
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_900.fwi              
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_912.5.fwi            
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_925_p1.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_925_p2.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_935_p1.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_935_p2.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_945_p1.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_945_p2.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_950_p1.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_950_p2.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_950_p3.fwi           
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_962.5_p1.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_962.5_p2.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_962.5_p3.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_987.5_p1.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_987.5_p2.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_987.5_p3.fwi         
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_975_p1.fwi
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_975_p2.fwi
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011/MHAD2011_975_p3.fwi
endif
if ([exp].eq.'mhad2012') then
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_505.3.runlist          
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_508.0_p1.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_508.0_p2.runlist       
*n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_508.0.runlist          
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_509.0_p1.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_509.0_p2.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_509.0_p3.runlist       
*n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_509.0.runlist          
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_510.0_p1.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_510.0_p2.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_510.0_p3.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_510.0_p4.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_510.0_p5.runlist       
*n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_510.0.runlist          
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_511.0_p1.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_511.0_p2.runlist       
*n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_511.0.runlist          
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_512.0.runlist          
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_512.3.runlist          
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_513.0.runlist          
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_640.0_p1.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_640.0_p2.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_640.0_p3.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_640.0_p4.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_640.0_p5.runlist       
*n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_640.0.runlist          
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_680.0_p1.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_680.0_p2.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_680.0_p3.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_680.0_p4.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_680.0_p5.runlist       
*n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_680.0.runlist          
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_720.0_p1.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_720.0_p2.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_720.0_p3.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_720.0_p4.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_720.0_p5.runlist       
*n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_720.0.runlist          
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_760.0_p1.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_760.0_p2.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_760.0_p3.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_760.0_p4.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_760.0_p5.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_760.0_p6.runlist       
*n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_760.0.runlist          
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_800.0_p1.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_800.0_p2.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_800.0_p3.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_800.0_p4.runlist       
*n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_800.0.runlist          
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_840.0_p1.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_840.0_p2.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_840.0_p3.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_840.0_p4.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_840.0_p5.runlist       
*n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_840.0.runlist          
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_860.0_p1.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_860.0_p2.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_860.0_p3.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_860.0_p4.runlist       
*n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_860.0.runlist          
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_880.0_p1.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_880.0_p2.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_880.0_p3.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_880.0_p4.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_880.0_p5.runlist       
*n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_880.0.runlist          
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_900.0_p1.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_900.0_p2.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_900.0_p3.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_900.0_p4.runlist       
*n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_900.0.runlist          
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_920.0_p1.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_920.0_p2.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_920.0_p3.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_920.0_p4.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_920.0_p5.runlist       
*n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_920.0.runlist          
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_936.0_p1.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_936.0_p2.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_936.0_p3.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_936.0_p4.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_936.0_p5.runlist       
*n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_936.0.runlist          
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_950.0_p1.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_950.0_p2.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_950.0_p3.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_950.0_p4.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_950.0_p5.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_950.0_p6.runlist       
*n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_950.0.runlist          
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_960.0_p1.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_960.0_p2.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_960.0_p3.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_960.0_p4.runlist       
*n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_960.0.runlist          
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_970.0_p1.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_970.0_p2.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_970.0_p3.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_970.0_p4.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_970.0_p5.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_970.0_p6.runlist       
*n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_970.0.runlist          
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_980.0_p1.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_980.0_p2.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_980.0_p3.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_980.0_p4.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_980.0_p5.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_980.0_p6.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_980.0_p7.runlist       
*n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_980.0.runlist          
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_990.0_p1.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_990.0_p2.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_990.0_p3.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_990.0_p4.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_990.0_p5.runlist       
n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_990.0_p6.runlist       
*n=[n]+1; f[n]=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2012/mhad2012_990.0.runlist          
endif
*
do i=1,[n]
  shell fgrep -e "runlist" [f[i]] > tmp.txt
  shell $unquote('cat tmp.txt | sed "s/runlist = \[/ /g" > tmp1.txt')
  shell $unquote('cat tmp1.txt | sed "s/,/\n/g" > tmp2.txt')
  shell $unquote('cat tmp2.txt | sed "s/\]/ /g" > tmp3.txt')
  ve/del rn
  ve/read rn tmp3.txt
  if ([i].eq.1) then
    ve/del runs11
    ve/copy rn runs11
  else
    if ($vlen(rn).ne.0) then
      exec vappend runs11 rn
    else
      mess [f[i]]
    endif
  endif
enddo
sigma runs11 = order(runs11,runs11)
return

macro filetest1
*exec mapcal#filetest mhad2012
ve/del runs11a
*ve/read runs11a mhad2011-4.list
ve/read runs11a mhad2012-2.list
nmax=$sigma(max(vmax(runs11),vmax(runs11a)))
nmin=$sigma(min(vmin(runs11),vmin(runs11a)))
ve/cre r1([nmax]) r
ve/cre r2([nmax]) r
do i=1,$vlen(runs11)
  ind=runs11([i])
  if ($sigma(r1([ind])).eq.1) then
    mess 1:[ind]
  endif
  ve/inp r1([ind]) 1
enddo
do i=1,$vlen(runs11a)
  ind=runs11a([i])
  if ($sigma(r2([ind])).eq.1) then
    mess 2:[ind]
  endif
  ve/inp r2([ind]) 1
enddo
ind=0
do i=[nmin],[nmax]
  r1i=r1([i])
  r2i=r2([i])
  if (([r1i].eq.0).and.([r2i].eq.1)) then
    exec mapcal#ndays [i]
    ind=[ind]+1
    mess 1,[ind]:[i] $sigma(beam(1))
  endif
enddo
do i=[nmin],[nmax]
  r1i=r1([i])
  r2i=r2([i])
  if (([r1i].eq.1).and.([r2i].eq.0)) then
    exec mapcal#ndays [i]
    mess 2:[i] $sigma(beam(1))
  endif
enddo
return



macro runsprep expn=PHI2013
ve/del rn
if ([expn].eq.'MHAD2011') then
  ve/read rn mhad2011-4.list
endif
if ([expn].eq.'MHAD2012') then
  ve/read rn mhad2012-2.list
endif
if ([expn].eq.'OMEGA2012') then
  ve/read rn omeg2012-0.list
endif
if ([expn].eq.'PHI2013') then
  ve/read rn phi_2013-0.list
endif
if ([expn].eq.'RHO2012') then
  ve/read rn rho_2012-0.list
endif
nmax=$vlen(rn)
ve/cre en([nmax]) r
ve/cre ep([nmax]) r
ve/cre ip1([nmax]) r
ve/cre ip2([nmax]) r
ind=0
ei=0
indi=0
ve/inp ip1(1) 1
do i=1,[nmax]
  run=rn([i])
  exec mapcal#ndays [run]
  beam=beam(1)
  mess [i]:[run] [beam]
  ve/inp en([i]) [beam]
  if ([ei].ne.[beam]) then
    indi=[indi]+1
    ve/inp ep([indi]) [beam]
    ve/inp ip1([indi]) [i]
    if ([indi].gt.1) then
      j=[indi]-1
      ve/inp ip2([j]) $sigma([i]-1)
    endif
    ei=[beam]
  endif
enddo
ve/inp ip2([indi]) [nmax]
exec mapcal#vecut ep
exec mapcal#vecut ip1
exec mapcal#vecut ip2
return


macro preprunlists expn=PHI2013
exec mapcal#runsprep [expn]
do l=1,$vlen(ep)
  e=ep([l])
  n1=ip1([l])
  n2=ip2([l])
  ve/del r
  ve/copy rn([n1]:[n2]) r
  n=$vlen(r)
  dn=10
  nf=$sigma(int([n]/[dn]))
  do j=0,[nf]
    list='runlist = ['
    i1=[j]*[dn]+1
    i2=$sigma(min([i1]+[dn]-1,[n]))
    do i=[i1],[i2]
      if ([i].lt.[i2]) then
        list=$unquote([list]) $sigma(r([i])),
      else
        list=$unquote([list]) $sigma(r([i]))
      endif
    enddo
    list=$unquote([list])$unquote(' ]')
    fo=[expn]_$format([e],f5.1)_p[j].runlist
    if ($fexist([fo]).eq.1) then
      shell rm [fo]
    endif
    txt=# [expn], точка $format([e],f5.1) МэВ, список файлов, часть [j]
    fmess [txt] [fo]
    fmess [list] [fo]
    mess File [fo] is prepared
  enddo
enddo
return



macro mvs
do i=6000,14000
  fname=run_[i]_eventtime_full_v2.
  fnameo=run_[i]_eventtime_full_v2.txt
  if ($fexist([fname]).eq.1) then
    shell mv -v [fname] [fnameo]
  endif
enddo
return

macro phicut
nmax=10000
ve/cre  dphi([nmax]) r
ve/cre ddphi([nmax]) r
ve/cre runs([nmax]) r
ve/cre druns([nmax]) r
ve/cre p3(3) r 100 0 1
ve/cre dp3(3) r
ind=0
do i=6000,14000
  fname=run_[i]_spects_v2.his
  if ($fexist([fname]).eq.1) then
    ind=[ind]+1
    ve/inp runs([ind]) [i]
    mess [fname]
    hi/del 6000
    hi/file 20 [fname]
    hrin 6000
    hi/pl 6000
    hi/fit 6000 g s 3 p3 ! ! ! dp3
    ve/inp dphi([ind]) $sigma(p3(3))
    ve/inp ddphi([ind]) $sigma(dp3(3))
    close 20
  endif
enddo
exec mapcal#vecut runs
n=$vlen(runs)
exec mapcal#vecut druns [n]
exec mapcal#vecut dphi
exec mapcal#vecut ddphi
* exec $PER/s#vpl dphi ddphi runs druns sz=0.05 
return

macro parcut par=dphi
imin=6000
imax=18000
nmax=[imax]-[imin]+1
ve/cre nevt([nmax]) r
ve/cre mean([nmax]) r
ve/cre rms([nmax]) r
ve/cre  w([nmax]) r
ve/cre dw([nmax]) r
ve/cre  v([nmax]) r
ve/cre dv([nmax]) r
ve/cre runs([nmax]) r
ve/cre druns([nmax]) r
ve/cre p3(3) r 100 0 1
ve/cre dp3(3) r
ve/cre step(3) r 1 0.1 0.1
ve/cre pmin(3) r 0 -5 0
ve/cre pmax(3) r 10000 5 5
ind=0
if ([par].eq.'dphi') then
  nh=6000
endif
if ([par].eq.'dtheta') then
  nh=7000
endif
if ([par].eq.'dz0') then
  nh=3000
endif
if ([par].eq.'dd0') then
  nh=5000
endif
if ([par].eq.'z0') then
  nh=2000
endif
set fcol 2
xmin=-10000
do i=[imin],[imax]
  fname=v2/run_[i]_spects_v2.his
  if ($fexist([fname]).eq.1) then
    ind=[ind]+1
    ve/inp runs([ind]) [i]
    mess [fname]
    hi/del [nh]
    hi/file 20 [fname]
    hrin [nh]
    hi/pl [nh]
    if ([xmin].eq.-10000) then
      xmin=$hinfo([nh],'xmin')
      xmax=$hinfo([nh],'xmax')
      xbin=$hinfo([nh],'xbins')
      ve/cre x([xbin]) r
      ve/cre y([xbin]) r
      xminp=$sigma([xmin]+([xmax]-[xmin])/[xbin]/2)
      xmaxp=$sigma([xmax]-([xmax]-[xmin])/[xbin]/2)
      sigma x = array([xbin],[xminp]#[xmaxp])
    endif
    hi/get/cont [nh] y
    nev=$sigma(vsum(y))
    mean=$sigma(vsum(x*y)/[nev])
    mean2=$sigma(vsum(x*x*y)/[nev])
    rms=$sigma(sqrt([mean2]-[mean]**2))
    ve/cre p3(3) r 50 [mean] [rms]
    ve/cre step(3) r 1 0 0
    hi/fit [nh] g sb 3 p3 step pmin pmax dp3
    ve/cre step(3) r 1 0.001 0.001
    do j=1,3
      l=$sigma(p3(2)-2*abs(p3(3)))
      r=$sigma(p3(2)+2*abs(p3(3)))
      hi/fit [nh]([l]:[r]) g sb 3 p3 step pmin pmax dp3
    enddo
    ve/inp nevt([ind]) [nev]
    ve/inp mean([ind]) [mean]
    ve/inp rms([ind]) [rms]
    ve/inp w([ind]) $sigma(p3(2))
    ve/inp dw([ind]) $sigma(abs(dp3(2)))
    ve/inp v([ind]) $sigma(abs(p3(3)))
    ve/inp dv([ind]) $sigma(abs(dp3(3)))
    close 20
  endif
enddo
exec mapcal#vecut runs
n=$vlen(runs)
exec mapcal#vecut druns [n]
exec mapcal#vecut w [n]
exec mapcal#vecut dw [n]
exec mapcal#vecut v [n]
exec mapcal#vecut dv [n]
exec mapcal#vecut nevt [n]
exec mapcal#vecut mean [n]
exec mapcal#vecut rms [n]
* exec $PER/s#vpl v dv runs druns sz=0.05
ve/write runs,nevt,mean,rms,w,dw,v,dv [par]_vs_runs.txt 8f15.6
return


macro parcutp par=dphi dir=v2 ver=0
ve/del runs,days,beams
if ([dir].eq.'v2') then
  ve/read runs,days,beams runpars.txt 3f15.6
else  
  ve/read runs,days,beams [dir]/runpars.txt 3f15.6
endif
ve/cre vnull(1) r 1
exec vappend runs vnull
exec vappend days vnull
exec vappend beams vnull
*
n=$vlen(runs)
ve/cre beamp([n]) r
ve/cre nevtp([n]) r
ve/cre meanp([n]) r
ve/cre rmsp([n]) r
ve/cre wp([n]) r
ve/cre dwp([n]) r
ve/cre vp([n]) r
ve/cre dvp([n]) r
ve/cre daysp([n]) r
ve/cre ddaysp([n]) r
ve/cre runsp([n]) r
ve/cre drunsp([n]) r
ind=0
if ([par].eq.'dphi') then
  nh=6000+[ver]
endif
if ([par].eq.'dtheta') then
  nh=7000+[ver]
endif
if ([par].eq.'dz0') then
  nh=3000+[ver]
endif
if ([par].eq.'dd0') then
  nh=5000+[ver]
endif
if ([par].eq.'z0') then
  nh=2000+[ver]
endif
set fcol 2
ind=0
bo=0
xmin=-10000
nhs=100000
ve/cre y(1) r
do i=1,[n]
  b=beams([i])
  if ([b].ne.[bo]) then
    if ([bo].ne.0) then
      i2=[i]-1
      if ([i1].le.[i2]) then
        ind=[ind]+1
        nev=$sigma(vsum(y))
        mean=$sigma(vsum(x*y)/[nev])
        mean2=$sigma(vsum(x*x*y)/[nev])
        rms=$sigma(sqrt([mean2]-[mean]**2))
        st=$sigma([rms]/100)
        ym=$sigma(vmax(y))
        ve/cre p3(3) r $sigma(vmax(y)) [mean] [rms]
        ve/cre dp3(3) r
        ve/cre pmin(3) r $sigma([ym]/10) $sigma([mean]-[rms]) $sigma([rms]/3)
        ve/cre pmax(3) r $sigma(10*[ym]) $sigma([mean]+[rms]) $sigma(3*[rms])
        hi/pl [nhs]
*        hi/fit [nhs] g s 3 p3 ! ! ! dp3
        ve/cre step(3) r 1 [st] [st]
        do j=1,3
          l=$sigma(p3(2)-(5-[j])*abs(p3(3)))
          r=$sigma(p3(2)+(5-[j])*abs(p3(3)))
          hi/fit [nhs]([l]:[r]) g sb 3 p3 step pmin pmax dp3
        enddo
        do j=1,3
          l=$sigma(p3(2)-2*abs(p3(3)))
          r=$sigma(p3(2)+2*abs(p3(3)))
          hi/fit [nhs]([l]:[r]) g sb 3 p3 step pmin pmax dp3
        enddo
*        hi/fit [nhs]([l]:[r]) g s 3 p3 ! ! ! dp3
        hi/pl [nhs]
        ve/inp nevtp([ind]) [nev]
        ve/inp meanp([ind]) [mean]
        ve/inp rmsp([ind]) [rms]
        ve/inp wp([ind]) $sigma(p3(2))
        ve/inp dwp([ind]) $sigma(abs(dp3(2)))
        ve/inp vp([ind]) $sigma(abs(p3(3)))
        ve/inp dvp([ind]) $sigma(abs(dp3(3)))
        ve/inp beamp([ind]) [bo]
        ve/inp daysp([ind]) $sigma(([t2]+[t1])/2)
        ve/inp ddaysp([ind]) $sigma(([t2]-[t1])/2)
        ve/inp runsp([ind]) $sigma(([r2]+[r1])/2)
        ve/inp drunsp([ind]) $sigma(([r2]-[r1])/2)
      endif
    endif
    bo=[b]
    i1=[i]
    t1=1000000
    t2=-1000000
    r1=1000000
    r2=-1000000
    sigma y = y*0
    if ($hexist([nhs]).eq.1) then
      hi/op/res [nhs]
    endif
  endif
  r=runs([i])
  fname=[dir]/run_[r]_spects_v2.his
  mess [b]: [fname]
  if ($fexist([fname]).eq.1) then
  hi/del [nh]
  hi/file 20 [fname]
  hrin [nh]
*  hi/pl [nh]
  if ([xmin].eq.-10000) then
    xmin=$hinfo([nh],'xmin')
    xmax=$hinfo([nh],'xmax')
    xbin=$hinfo([nh],'xbins')
    ve/cre x([xbin]) r
    ve/cre y([xbin]) r
    ve/cre yi([xbin]) r
    xminp=$sigma([xmin]+([xmax]-[xmin])/[xbin]/2)
    xmaxp=$sigma([xmax]-([xmax]-[xmin])/[xbin]/2)
    sigma x = array([xbin],[xminp]#[xmaxp])
    hi/copy [nh] [nhs]
    hi/op/res [nhs]
  endif
  hi/op/add [nhs] [nh] [nhs]
  hi/get/cont [nh] yi
  close 20
  sigma y = y + yi
  daysi=days([i])
  if ([t1].gt.[daysi]) then
    t1=[daysi]
  endif
  if ([t2].lt.[daysi]) then
    t2=[daysi]
  endif
*  t1=$sigma(min([t1],[daysi]))
*  t2=$sigma(max([t2],[daysi]))
  runsi=runs([i])
  if ([r1].gt.[runsi]) then
    r1=[runsi]
  endif
  if ([r2].lt.[runsi]) then
    r2=[runsi]
  endif
*  r1=$sigma(min([r1],[runsi]))
*  r2=$sigma(max([r2],[runsi]))
  endif
enddo
exec mapcal#vecut beamp
n=$vlen(beamp)
exec mapcal#vecut nevtp [n]
exec mapcal#vecut meanp [n]
exec mapcal#vecut rmsp [n]
exec mapcal#vecut wp [n]
exec mapcal#vecut dwp [n]
exec mapcal#vecut vp [n]
exec mapcal#vecut dvp [n]
exec mapcal#vecut daysp [n]
exec mapcal#vecut ddaysp [n]
exec mapcal#vecut runsp [n]
exec mapcal#vecut drunsp [n]
sigma dbeams = beams*0
sigma dbeamp = beamp*0
sigma dmeanp = rmsp/sqrt(nevtp)
sigma drmsp = rmsp/sqrt(nevtp*2)
*
ve/write runsp,drunsp,daysp,ddaysp,beamp,dbeamp,nevtp,meanp,dmeanp,rmsp,drmsp,wp,dwp,vp,dvp [dir]/[par]_vs_beams_new_cut[ver].txt 15e15.6
return


macro runpar
ve/del runs,nevt,mean,rms,w,dw,v,dv
ve/read runs,nevt,mean,rms,w,dw,v,dv dphi_vs_runs.txt 8f15.6
n=$vlen(runs)
ve/cre days([n]) r
ve/cre beams([n]) r
do i=1,[n]
  r=runs([i])
  exec mapcal#ndays [r]
  gl/imp ndays
  gl/imp ebeam 
  mess $sigma(int([i]/[n]*10000+0.5)/100) %    [ndays] [ebeam]
  ve/inp days([i]) [ndays]
  ve/inp beams([i]) [ebeam]
enddo
*ve/write runs,days,beams runpars.txt 3f15.6
return

macro runpars
nmax=10000
ve/cre runs([nmax]) r
ve/cre days([nmax]) r
ve/cre beams([nmax]) r
ind=0
do i=1,10000
  fname=vs/run_[i]_spects_v2.his
  if ($fexist([fname]).eq.1) then
    fname=vs/run_[i]_eventtime_beam_v2.txt
    ve/read p3 [fname]
    ind=[ind]+1
    ve/inp runs([ind]) [i]
    ve/inp days([ind]) 0
    ve/inp beams([ind]) $sigma(int(p3(2)*10+0.5)/10)
    mess $sigma(int(p3(2)*10+0.5)/10) $sigma(p3(2))
  endif
enddo
exec mapcal#vecut runs [ind]
exec mapcal#vecut days [ind]
exec mapcal#vecut beams [ind]
ve/write runs,days,beams vs/runpars.txt 3f15.6
return


macro parcutpl par=dphi auto=0 ver=0
ve/del runs,nevt,mean,rms,w,dw,v,dv
ve/read runs,nevt,mean,rms,w,dw,v,dv [par]_vs_runs.txt 8f15.6
ve/del runs,days,beams
ve/read runs,days,beams runpars.txt 3f15.6
*
n=$vlen(runs)
ve/cre druns([n]) r
sigma dmean = rms/sqrt(nevt)
sigma drms = rms/sqrt(nevt*2)
*
goto 1
*
ve/cre beamp([n]) r
ve/cre nevtp([n]) r
ve/cre meanp([n]) r
ve/cre rmsp([n]) r
ve/cre wp([n]) r
ve/cre dwp([n]) r
ve/cre vp([n]) r
ve/cre dvp([n]) r
ve/cre daysp([n]) r
ve/cre ddaysp([n]) r
ve/cre runsp([n]) r
ve/cre drunsp([n]) r
*
ind=0
bo=0
do i=1,[n]
  b=beams([i])
  if ([b].ne.[bo]) then
    if ([bo].ne.0) then
      i2=[i]-1
      if ([i1].le.[i2]) then
        ind=[ind]+1
        ve/del nevtr,meanr,rmsr,wr,dwr,vr,dvr
        ve/copy nevt([i1]:[i2]) nevtr
        ve/copy mean([i1]:[i2]) meanr
        ve/copy rms([i1]:[i2]) rmsr
        ve/copy w([i1]:[i2]) wr
        ve/copy dw([i1]:[i2]) dwr
        ve/copy v([i1]:[i2]) vr
        ve/copy dv([i1]:[i2]) dvr
        ve/inp beamp([ind]) [bo]
        nevti=$sigma(vsum(nevtr))
        ve/inp nevtp([ind]) [nevti]
        meani=$sigma(vsum(meanr*nevtr)/[nevti])
        ve/inp meanp([ind]) [meani]
        mean2i=$sigma(vsum((rmsr**2+([meani])**2)*nevtr)/[nevti])
        ve/inp rmsp([ind]) $sigma(sqrt([mean2i]-[meani]**2))
        dd=$sigma(vsum(1/dwr**2))
        ve/inp dwp([ind]) $sigma(1/sqrt([dd]))
        ve/inp wp([ind]) $sigma(vsum(wr/dwr**2)/[dd])
        dd=$sigma(vsum(1/dvr**2))
        ve/inp dvp([ind]) $sigma(1/sqrt([dd]))
        ve/inp vp([ind]) $sigma(vsum(vr/dvr**2)/[dd])
        ve/inp daysp([ind]) $sigma(([t2]+[t1])/2)
        ve/inp ddaysp([ind]) $sigma(([t2]-[t1])/2)
        ve/inp runsp([ind]) $sigma(([r2]+[r1])/2)
        ve/inp drunsp([ind]) $sigma(([r2]-[r1])/2)
      endif
    endif
    bo=[b]
    i1=[i]
    t1=1000000
    t2=-1000000
    r1=1000000
    r2=-1000000
  endif
  daysi=days([i])
  if ([t1].gt.[daysi]) then
    t1=[daysi]
  endif
  if ([t2].lt.[daysi]) then
    t2=[daysi]
  endif
*  t1=$sigma(min([t1],[daysi]))
*  t2=$sigma(max([t2],[daysi]))
  runsi=runs([i])
  if ([r1].gt.[runsi]) then
    r1=[runsi]
  endif
  if ([r2].lt.[runsi]) then
    r2=[runsi]
  endif
*  r1=$sigma(min([r1],[runsi]))
*  r2=$sigma(max([r2],[runsi]))
enddo
exec mapcal#vecut beamp
n=$vlen(beamp)
exec mapcal#vecut nevtp [n]
exec mapcal#vecut meanp [n]
exec mapcal#vecut rmsp [n]
exec mapcal#vecut wp [n]
exec mapcal#vecut dwp [n]
exec mapcal#vecut vp [n]
exec mapcal#vecut dvp [n]
exec mapcal#vecut daysp [n]
exec mapcal#vecut ddaysp [n]
exec mapcal#vecut runsp [n]
exec mapcal#vecut drunsp [n]
sigma dbeams = beams*0
sigma dbeamp = beamp*0
sigma dmeanp = rmsp/sqrt(nevtp)
sigma drmsp = rmsp/sqrt(nevtp*2)
*
ve/write runsp,drunsp,daysp,ddaysp,beamp,dbeamp,nevtp,meanp,dmeanp,rmsp,drmsp,wp,dwp,vp,dvp [par]_vs_beams.txt 15e15.6
*
1:
sigma dbeams = beams*0
ve/del runsp,drunsp,daysp,ddaysp,beamp,dbeamp,nevtp,meanp,dmeanp,rmsp,drmsp,wp,dwp,vp,dvp
*ve/read runsp,drunsp,daysp,ddaysp,beamp,dbeamp,nevtp,meanp,dmeanp,rmsp,drmsp,wp,dwp,vp,dvp [par]_vs_beams.txt 15e15.6
ve/read runsp,drunsp,daysp,ddaysp,beamp,dbeamp,nevtp,meanp,dmeanp,rmsp,drmsp,wp,dwp,vp,dvp v2/[par]_vs_beams_new_cut[ver].txt 15e15.6
*
if ([par].eq.'dphi') then
  tm=[Df], degree
  tr=[s]?[Df]!, degree
endif
if ([par].eq.'dtheta') then
  tm=[Dq], dergee
  tr=[s]?[Dq]!, degree
endif
if ([par].eq.'dz0') then
  tm=[Dz], cm
  tr=[s]?[Dz]!, cm
endif
if ([par].eq.'dd0') then
  tm=[Dr], cm
  tr=[s]?[Dr]!, cm
endif
if ([par].eq.'z0') then
  tm=[z]0, cm
  tr=[s]?[z]0!, cm
endif
*
ks=0
sigma runsm = runs + 50*[ks]
sigma runspm = runsp + 50*[ks]
sigma beampm = beamp + 5*[ks]
sigma dayspm = daysp + 5*[ks]
*
ve/cre runsx(10) r 7842 17868 7842 10988 11285 13846 13847 16698 16699 17868
exec mapcal#dcuts runsx daysx
do i=0,4
  set pmci 1
*
  j=[i]*2+1
  tx=runsx([j])
  exec mapcal#ixndl [tx] runs
  gl/imp inx
  in1=[inx]
  j=[i]*2+2
  tx=runsx([j])
  exec mapcal#ixndr [tx] runs
  gl/imp inx
  in2=[inx]

  ve/del vt,dvt,runst,drunst,rmst,drmst,runsmt
  ve/copy  v([in1]:[in2])  vt
  ve/copy dv([in1]:[in2]) dvt
  ve/copy  runs([in1]:[in2])  runst
  ve/copy druns([in1]:[in2]) drunst
  ve/copy  rms([in1]:[in2])  rmst
  ve/copy drms([in1]:[in2]) drmst
  ve/copy  runsm([in1]:[in2])  runsmt
  u0=$sigma(max(vmax(vt),vmax(rmst)))
  d0=$sigma(min(vmin(vt),vmin(rmst)))
  l0=$sigma(vmin(runst))
  r0=$sigma(vmax(runst))
  d=$sigma([d0]-0.2*([u0]-[d0]))
  u=$sigma([u0]+0.2*([u0]-[d0]))
  l=$sigma([l0]-0.05*([r0]-[l0]))
  r=$sigma([r0]+0.05*([r0]-[l0]))
  null [l] [r] [d] [u]
*set pmci 4
  * exec $PER/s#vpl vt dvt runst drunst sz=0.1 iatt=20 o=s
*set pmci 2
  * exec $PER/s#vpl rmst drmst runsmt drunst sz=0.1 iatt=24 o=s
  atitle 'runs' [tr]
  exec save v2/mapcal[i]_[par]_rms_vs_runs0_cut[ver].eps f
  if ([auto].eq.0) then
    read x
  endif
*
  ve/del wt,dwt,rmst,drmst
  ve/copy  w([in1]:[in2])  wt
  ve/copy dw([in1]:[in2]) dwt
  ve/copy  mean([in1]:[in2])  meant
  ve/copy dmean([in1]:[in2]) dmeant
  u0=$sigma(max(vmax(wt),vmax(meant)))
  d0=$sigma(min(vmin(wt),vmin(meant)))
  l0=$sigma(vmin(runst))
  r0=$sigma(vmax(runst))
  d=$sigma([d0]-0.2*([u0]-[d0]))
  u=$sigma([u0]+0.2*([u0]-[d0]))
  l=$sigma([l0]-0.05*([r0]-[l0]))
  r=$sigma([r0]+0.05*([r0]-[l0]))
  null [l] [r] [d] [u]
*set pmci 4
  * exec $PER/s#vpl wt dwt runst drunst sz=0.1 iatt=20 ll=-1 o=s
*set pmci 2
  * exec $PER/s#vpl meant dmeant runsmt drunst sz=0.1 iatt=24 ll=-1 o=s
  atitle 'runs' [tm]
  exec save v2/mapcal[i]_[par]_mean_vs_runs0_cut[ver].eps f
  if ([auto].eq.0) then
    read x
  endif
*
  j=[i]*2+1
  tx=daysx([j])
  exec mapcal#ixndl [tx] daysp
  gl/imp inx
  in1=[inx]
  j=[i]*2+2
  tx=daysx([j])
  exec mapcal#ixndr [tx] daysp
  gl/imp inx
  in2=[inx]

  ve/del vpt,dvpt,beampt,dbeampt,rmspt,drmspt,beampmt
  ve/copy  vp([in1]:[in2])  vpt
  ve/copy dvp([in1]:[in2]) dvpt
  ve/copy  beamp([in1]:[in2])  beampt
  ve/copy dbeamp([in1]:[in2]) dbeampt
  ve/copy  rmsp([in1]:[in2])  rmspt
  ve/copy drmsp([in1]:[in2]) drmspt
  ve/copy  beampm([in1]:[in2])  beampmt
  u0=$sigma(max(vmax(vpt),vmax(rmspt)))
  d0=$sigma(min(vmin(vpt),vmin(rmspt)))
  l0=$sigma(vmin(beampt))
  r0=$sigma(vmax(beampt))
  d=$sigma([d0]-0.2*([u0]-[d0]))
  u=$sigma([u0]+0.2*([u0]-[d0]))
  l=$sigma([l0]-0.05*([r0]-[l0]))
  r=$sigma([r0]+0.05*([r0]-[l0]))
  null [l] [r] [d] [u]
**set pmci 1
** exec $PER/s#vpl rms drms beams dbeams sz=0.01
**set pmci 2
*set pmci 2
  * exec $PER/s#vpl rmspt drmspt beampt dbeampt sz=0.1 iatt=24 o=s
*set pmci 4
  * exec $PER/s#vpl vpt dvpt beampt dbeampt sz=0.1 iatt=20 o=s
  atitle 'beam, MeV' [tr]
  exec mapcal#partnamepl
  exec save v2/mapcal[i]_[par]_rms_vs_beam_cut[ver].eps f
  if ([auto].eq.0) then
    read x
  endif
*
  ve/del wpt,dwpt,meanpt,dmeanpt
  ve/copy  wp([in1]:[in2])  wpt
  ve/copy dwp([in1]:[in2]) dwpt
  ve/copy  meanp([in1]:[in2])  meanpt
  ve/copy dmeanp([in1]:[in2]) dmeanpt
  u0=$sigma(max(vmax(wpt),vmax(meanpt)))
  d0=$sigma(min(vmin(wpt),vmin(meanpt)))
  l0=$sigma(vmin(beampt))
  r0=$sigma(vmax(beampt))
  d=$sigma([d0]-0.2*([u0]-[d0]))
  u=$sigma([u0]+0.2*([u0]-[d0]))
  l=$sigma([l0]-0.05*([r0]-[l0]))
  r=$sigma([r0]+0.05*([r0]-[l0]))
  null [l] [r] [d] [u]
**set pmci 1
** exec $PER/s#vpl mean dmean beams dbeams sz=0.01 ll=-1
**set pmci 2
*set pmci 2
  * exec $PER/s#vpl meanpt dmeanpt beampt dbeampt sz=0.1 iatt=24 ll=-1 o=s
*set pmci 4
  * exec $PER/s#vpl wpt dwpt beampt dbeampt sz=0.1 iatt=20 ll=-1 o=s
  atitle 'beam, MeV' [tm]
  exec mapcal#partnamepl
  exec save v2/mapcal[i]_[par]_mean_vs_beam_cut[ver].eps f
  if ([auto].eq.0) then
    read x
  endif
*
  ve/del dayspt,ddayspt,dayspmt
  ve/copy  daysp([in1]:[in2])  dayspt
  ve/copy ddaysp([in1]:[in2]) ddayspt
  ve/copy  dayspm([in1]:[in2])  dayspmt
  u0=$sigma(max(vmax(vpt),vmax(rmspt)))
  d0=$sigma(min(vmin(vpt),vmin(rmspt)))
  l0=$sigma(vmin(dayspt))
  r0=$sigma(vmax(dayspt))
  d=$sigma([d0]-0.2*([u0]-[d0]))
  u=$sigma([u0]+0.2*([u0]-[d0]))
  l=$sigma([l0]-0.05*([r0]-[l0]))
  r=$sigma([r0]+0.05*([r0]-[l0]))
  null [l] [r] [d] [u]
*set pmci 2
  * exec $PER/s#vpl rmspt drmspt dayspt ddayspt sz=0.1 iatt=24 o=s
*set pmci 4
  * exec $PER/s#vpl vpt dvpt dayspt ddayspt sz=0.1 iatt=20 o=s
  atitle 'time, days' [tr]
  exec save v2/mapcal[i]_[par]_rms_vs_days_cut[ver].eps f
  if ([auto].eq.0) then
    read x
  endif
*
  u0=$sigma(max(vmax(wpt),vmax(meanpt)))
  d0=$sigma(min(vmin(wpt),vmin(meanpt)))
  l0=$sigma(vmin(dayspt))
  r0=$sigma(vmax(dayspt))
  d=$sigma([d0]-0.2*([u0]-[d0]))
  u=$sigma([u0]+0.2*([u0]-[d0]))
  l=$sigma([l0]-0.05*([r0]-[l0]))
  r=$sigma([r0]+0.05*([r0]-[l0]))
  null [l] [r] [d] [u]
*set pmci 2
  * exec $PER/s#vpl meanpt dmeanpt dayspt ddayspt sz=0.1 iatt=24 ll=-1 o=s
*set pmci 4
  * exec $PER/s#vpl wpt dwpt dayspt ddayspt sz=0.1 iatt=20 ll=-1 o=s
  atitle 'time, days' [tm]
  exec save v2/mapcal[i]_[par]_mean_vs_days_cut[ver].eps f
  if ([auto].eq.0) then
    read x
  endif  
*
  ve/del runspt,drunspt,runspmt
  ve/copy  runsp([in1]:[in2])  runspt
  ve/copy drunsp([in1]:[in2]) drunspt
  ve/copy  runspm([in1]:[in2])  runspmt
  u0=$sigma(max(vmax(vpt),vmax(rmspt)))
  d0=$sigma(min(vmin(vpt),vmin(rmspt)))
  l0=$sigma(vmin(runspt))
  r0=$sigma(vmax(runspt))
  d=$sigma([d0]-0.2*([u0]-[d0]))
  u=$sigma([u0]+0.2*([u0]-[d0]))
  l=$sigma([l0]-0.05*([r0]-[l0]))
  r=$sigma([r0]+0.05*([r0]-[l0]))
  null [l] [r] [d] [u]
*set pmci 2
  * exec $PER/s#vpl rmspt drmspt runspt drunspt sz=0.1 iatt=24 o=s
*set pmci 4
  * exec $PER/s#vpl vpt dvpt runspt drunspt sz=0.1 iatt=20 o=s
  atitle 'runs' [tr]
  exec save v2/mapcal[i]_[par]_rms_vs_runs_cut[ver].eps f
  if ([auto].eq.0) then
    read x
  endif
*
  u0=$sigma(max(vmax(wpt),vmax(meanpt)))
  d0=$sigma(min(vmin(wpt),vmin(meanpt)))
  l0=$sigma(vmin(runspt))
  r0=$sigma(vmax(runspt))
  d=$sigma([d0]-0.2*([u0]-[d0]))
  u=$sigma([u0]+0.2*([u0]-[d0]))
  l=$sigma([l0]-0.05*([r0]-[l0]))
  r=$sigma([r0]+0.05*([r0]-[l0]))
  null [l] [r] [d] [u]
*set pmci 2
  * exec $PER/s#vpl meanpt dmeanpt runspt drunspt sz=0.1 iatt=24 ll=-1 o=s
*set pmci 4
  * exec $PER/s#vpl wpt dwpt runspt drunspt sz=0.1 iatt=20 ll=-1 o=s
  atitle 'runs' [tm]
  exec save v2/mapcal[i]_[par]_mean_vs_runs_cut[ver].eps f
  if ([auto].eq.0) then
    read x
  endif
*
enddo
return


macro partnamepl
npt=0
npt=[npt]+1; ms[npt]=775.49/2; nm[npt]='[r]'; dy[npt]=0
npt=[npt]+1; ms[npt]=782.65/2; nm[npt]='[w]'; dy[npt]=0.1
npt=[npt]+1; ms[npt]=1019.455/2; nm[npt]='[f]'; dy[npt]=0
lw=$GRAFINFO('WNXMIN')
rw=$GRAFINFO('WNXMAX')
dw=$GRAFINFO('WNYMIN')
uw=$GRAFINFO('WNYMAX')
mess [lw] [rw] [uw] [dw]
ud=[uw]-[dw]
do i=1,[npt]
  if (([ms[i]].ge.[lw]).and.([ms[i]].le.[rw])) then
    arrow [ms[i]] [ms[i]] $sigma([dw]+0.15*[ud]+[dy[i]]) $sigma([dw]+0.05*[ud]+[dy[i]]) $sigma(0.1*[ud])
    itx $sigma([ms[i]]+0.05) $sigma([dw]+0.05*[ud]+[dy[i]]) [nm[i]]
  endif
enddo
return



macro runcorr dir=v2
hi/file 30 mhad2011-4_profilex_ee_sim_shifter2_chit.his
ve/cre runs(20000) r
ve/cre ampx(9) r
do i=1,9
  hrin 10000[i]
  ve/inp ampx([i]) $hinfo(10000[i],'mean')
  ve/cre qecorr[i](20000) r
  ve/cre corr[i](20000) r
  ve/cre sigs[i](20000) r
  ve/del p14
  ve/read p14 [dir]/amplitude_vs_times_fourerf_counter[i]_[dir].par
  ve/del p12c[i]
  ve/copy p14(1:12) p12c[i]
  ve/del pei,dpei,effi,deffi,mueic,dmueic,r,dr,xei,dxei,kei
  ve/read pei,dpei,effi,deffi,mueic,dmueic,r,dr,xei,dxei,kei [dir]/parday_counter[i]_[dir].txt 11e15.6
  ve/del a[i],da[i],time,dtime
  ve/copy  pei     a[i]
  ve/copy dpei    da[i]
  ve/copy  xei  time
  ve/copy dxei dtime
  ve/del runs0,days,amp[i],damp[i]
  ve/read runs0,days,amp[i],damp[i] [dir]/amplitude_vs_runs_counter[i].txt 4f15.6
*
  ve/copy p12c[i] p12
  fun/pl fourerfpa.f -100 800
  ve/cre xi(1) r 0
  axi=$call('fourerfpa.f(xi)')
  if ([axi].gt.$sigma(ampx([i]))) then
  dxi=100
  set pmci 4
  while ($sigma(abs([dxi])).gt.0.001) do
    key $sigma(xi(1)) [axi] 20 ! 0.05
    ve/inp xi(1) $sigma(xi(1)+[dxi])
    axi=$call('fourerfpa.f(xi)')
    if ($sigma([dxi]*([axi]-ampx([i]))).lt.0) then
      dxi=-[dxi]/2
    endif
*    mess $sigma(ampx([i])) [axi] $sigma(xi(1)) [dxi]
  endwhile
  set pmci 2
  key $sigma(xi(1)) [axi] 20 ! 0.1
  else
    mess $sigma(ampx([i])) [axi]
    set pmci 4
    key $sigma(xi(1)) $sigma(ampx([i])) 20 ! 0.3
    set pmci 2
    key $sigma(xi(1)) [axi] 20 ! 0.1
    ve/inp ampx([i]) [axi]
  endif
*  read x
enddo
close 30
*read x
tmin=$sigma(vmin(time-dtime)-1.0/24)
tmax=$sigma(vmax(time+dtime)+1.0/24)
nd=$sigma(int(([tmax]-([tmin]))*24*6))
ve/cre tn([nd]) r
dt=$sigma(([tmax]-[tmin])/[nd])
do i=1,$vlen(time)
  t=$sigma(time([i])-dtime([i]))
  l=$sigma(int(([t]-[tmin])/[dt])+1)
  t=$sigma(time([i])+dtime([i]))
  r=$sigma(int(([t]-[tmin])/[dt])+1)
  do j=[l],[r]
    ve/inp tn([j]) [i]
  enddo
enddo
ind=0
*fun/pl fourerfp.f -100 800 
*fun/pl fourerfpa.f -100 800 s
ve/cre xi(1) r
ve/cre aexp(1) r
ve/cre daexp(1) r
do i=5000,20000
  fname=v2/run_[i]_spects_v2.his
  if ($fexist([fname]).eq.1) then
    ind=[ind]+1
    ve/inp runs([ind]) [i]
    exec mapcal#ndays [i]
    gl/imp ndays
    ve/inp xi(1) [ndays]
    nj=$sigma(int(([ndays]-[tmin])/[dt])+1)
    j=tn([nj])
    if ([j].eq.0) then
      njp=[nj]
      while (([j].eq.0).and.([njp].lt.[nd])) do
        njp=[njp]+1
        j=tn([njp])
      endwhile
      j=0
      njm=[nj]
      while (([j].eq.0).and.([njm].gt.1)) do
        njm=[njm]-1
        j=tn([njm])
      endwhile
      if ($sigma([njp]-[nj]).lt.$sigma([nj]-[njm])) then
        j=tn([njp])
      else
        j=tn([njm])
      endif
    endif
*    nj=$vlen(time)
*    do j=1,[nj]
*      l=$sigma(time([j])-dtime([j])-0.007)
*      r=$sigma(time([j])+dtime([j])+0.007)
*      if (([ndays].ge.[l]).and.([ndays].le.[r])) then
        do nc=1,9
          ve/del p12
          ve/copy p12c[nc] p12
          af=$call('fourerfpa.f(xi)')
          ve/copy a[nc]([j]) aexp(1)
          ae=aexp(1)
          ve/copy amp[nc]([ind]) aexp(1)
          ve/copy damp[nc]([ind]) daexp(1)
          aex=aexp(1)
          daex=daexp(1)
          mess [i]($sigma(runs0([ind]))) [ind] [nc] [j] [af] [ae] $sigma(ampx([nc]))
          ve/inp corr[nc]([ind]) $sigma(min(abs([af]/[ae]),1000))
          ds=([ae]-([aex]))/([daex])
          ve/inp sigs[nc]([ind]) $sigma(min(abs([ds]),1000))
          ve/inp qecorr[nc]([ind]) $sigma([af]/ampx([nc]))
        enddo
*      endif
*    enddo
  endif
enddo
exec mapcal#vecut runs
n=$vlen(runs)
do i=1,9
  exec mapcal#vecut corr[i] [n]
  exec mapcal#vecut sigs[i] [n]
  exec mapcal#vecut qecorr[i] [n]
enddo
ve/write runs,corr1,corr2,corr3,corr4,corr5,corr6,corr7,corr8,corr9 runscorr.txt 10f15.6
ve/write runs,sigs1,sigs2,sigs3,sigs4,sigs5,sigs6,sigs7,sigs8,sigs9 runssigs.txt 10f15.6
ve/write runs,qecorr1,qecorr2,qecorr3,qecorr4,qecorr5,qecorr6,qecorr7,qecorr8,qecorr9 runsqecorr.txt 10f15.6
return


macro qereduce
ve/del runs,corr1,corr2,corr3,corr4,corr5,corr6,corr7,corr8,corr9
ve/read runs,corr1,corr2,corr3,corr4,corr5,corr6,corr7,corr8,corr9 runscorr.txt 10f15.6
ve/del runs,sigs1,sigs2,sigs3,sigs4,sigs5,sigs6,sigs7,sigs8,sigs9
ve/read runs,sigs1,sigs2,sigs3,sigs4,sigs5,sigs6,sigs7,sigs8,sigs9 runssigs.txt 10f15.6
ve/del runs,qecorr1,qecorr2,qecorr3,qecorr4,qecorr5,qecorr6,qecorr7,qecorr8,qecorr9
ve/read runs,qecorr1,qecorr2,qecorr3,qecorr4,qecorr5,qecorr6,qecorr7,qecorr8,qecorr9 runsqecorr.txt 10f15.6
ve/del r1,r2
ve/read r1,r2 daycuts.txt 2f10.2
n=$vlen(r1)
ve/cre runsx([n]) r
do i=1,9
  ve/cre qecorr[i]x([n]) r
enddo
do i=1,[n]
  r1i=r1([i])
  r2i=r2([i])
  do c=1,9
    c[c]=0
  enddo
  ni=0
  do j=1,$vlen(runs)
    rj=runs([j])
    if (([rj].ge.[r1i]).and.([rj].le.[r2i])) then
      do c=1,9
        c[c]=$sigma([c[c]]+qecorr[c]([j]))
      enddo
      ni=[ni]+1
    endif
  enddo
  ve/inp runsx([i]) [r1i]
  do c=1,9
    c[c]=[c[c]]/[ni]
    ve/inp qecorr[c]x([i]) [c[c]]
  enddo
enddo  
ve/write runsx,qecorr1x,qecorr2x,qecorr3x,qecorr4x,qecorr5x,qecorr6x,qecorr7x,qecorr8x,qecorr9x runsqecorrx.txt 10f15.6
return



macro qereduceold
ve/del runs,corr1,corr2,corr3,corr4,corr5,corr6,corr7,corr8,corr9
ve/read runs,corr1,corr2,corr3,corr4,corr5,corr6,corr7,corr8,corr9 runscorr.txt 10f15.6
ve/del runs,sigs1,sigs2,sigs3,sigs4,sigs5,sigs6,sigs7,sigs8,sigs9
ve/read runs,sigs1,sigs2,sigs3,sigs4,sigs5,sigs6,sigs7,sigs8,sigs9 runssigs.txt 10f15.6
ve/del runs,qecorr1,qecorr2,qecorr3,qecorr4,qecorr5,qecorr6,qecorr7,qecorr8,qecorr9
ve/read runs,qecorr1,qecorr2,qecorr3,qecorr4,qecorr5,qecorr6,qecorr7,qecorr8,qecorr9 runsqecorr.txt 10f15.6
n=$vlen(runs)
ve/cre runsx([n]) r
do i=1,9
  ve/cre qecorr[i]x([n]) r
enddo
p1=0
ni=0
ind=0
do i=1,[n]
  p2=corr1([i])
  if ([p2].ne.[p1]) then
    if ([ni].ne.0) then
      ind=[ind]+1
      do c=1,9
        c[c]=[c[c]]/[ni]
        ve/inp qecorr[c]x([ind]) [c[c]]
      enddo
      ve/inp runsx([ind]) [run]
      mess [run] [c1] [c2] [c3] [c4] [c5] [c6] [c7] [c8] [c9] [ni]
    endif
    p1=[p2]
    run=runs([i])
    do c=1,9
      c[c]=0
    enddo
    ni=0
  endif
  do c=1,9
    c[c]=$sigma([c[c]]+qecorr[c]([i]))
  enddo
  ni=[ni]+1
enddo
exec mapcal#vecut runsx
n=$vlen(runsx)
do i=1,9
  exec mapcal#vecut qecorr[i]x [n]
enddo
ve/write runsx,qecorr1x,qecorr2x,qecorr3x,qecorr4x,qecorr5x,qecorr6x,qecorr7x,qecorr8x,qecorr9x runsqecorrx.txt 10f15.6
return



macro wlsspectzx 
do i=1,9
  exec mapcal#wlsspectz [i]
enddo
return


macro wlsspectz nc=1
gl/imp corrz
ver=[corrz]-1
ncz=10
zmin=-10
zmax=9
*goto 1
ve/cre  asz([ncz]) r 
ve/cre dasz([ncz]) r
ve/cre  aszs([ncz]) r 
ve/cre daszs([ncz]) r
ve/cre  zsz([ncz]) r 
ve/cre dzsz([ncz]) r
do i=1,[ncz]
  exec mapcal#wlsspect [nc] [i]
  ve/inp  asz([i]) $sigma(p3(2))
  ve/inp dasz([i]) $sigma(dp3(2))
  ve/inp  aszs([i]) $sigma(p3s(2))
  ve/inp daszs([i]) $sigma(dp3s(2))
  ve/inp  zsz([i]) $sigma([zmin]+([i]-0.5)*([zmax]-[zmin])/[ncz])
enddo
* exec $PER/s#vpl asz dasz zsz dzsz sz=0.1
* exec $PER/s#vpl aszs daszs zsz dzsz sz=0.1 iatt=24 o=s
n1=1
n2=$sigma(min(7,[ncz]))
ve/del aszr,daszr,zszr,dzszr
ve/copy  asz([n1]:[n2])  aszr
ve/copy dasz([n1]:[n2]) daszr
ve/copy  zsz([n1]:[n2])  zszr
ve/copy dzsz([n1]:[n2]) dzszr
ve/cre wls(3) r 0.5 0.5 1000000
ve/cre p3(3) r 10 10 50
ve/cre dp3(3) r
do i=1,2
  ve/fit zszr aszr daszr wlscorr.f s 3 p3 ! ! ! dp3
enddo
sigma p3=abs(p3)
ve/write p3,dp3 corrwls_counter[nc]_v0.txt '2e15.6'
*
read x
*
1:
*
ver=[corrz]-1
shell fgrep -e Shifter mapcal2011-4shifter_v[ver].cal >& tmp.txt
shell $unquote('cat tmp.txt | sed "s/Shifter/ /g" > tmp1.txt')
ve/del ab,af,tau
ve/read ab,af,tau tmp1.txt
ve/cre wls(3) r $sigma(ab([nc])) $sigma(af([nc])) $sigma(tau([nc]))
*
ve/cre cwls([ncz]) r
ve/cre dcwls([ncz]) r
ve/copy zsz zwls
ve/copy dzsz dzwls
a=wls(1)
b=wls(2)
c=wls(3)
ve/cre aszm([ncz]) r
ve/cre daszm([ncz]) r
do i=1,[ncz]
  z=zwls([i])
  f=$sigma(([a])*exp(-([z])/([c]))+([b])*exp(([z])/([c])))
  ve/inp aszm([i]) [f]
enddo
sigma cwls = aszm*asz/aszs
sigma dcwls = cwls*sqrt((dasz/asz)**2+(daszs/aszs)**2+(daszm/aszm))
null -10 10 0 $sigma(1.5*min(30,vmax(cwls)))
* exec $PER/s#vpl cwls dcwls zwls dzwls sz=0.1 o=s
ve/del p3(3)
ve/copy wls p3
ve/cre dp3(3) r
n1=1
n2=7
ve/del cwlsr,dcwlsr,zwlsr
ve/copy  cwls([n1]:[n2])  cwlsr
ve/copy dcwls([n1]:[n2]) dcwlsr
ve/copy  zwls([n1]:[n2])  zwlsr
ve/cre wls(3) r 0.5 0.5 1000000
ve/fit zwlsr cwlsr dcwlsr wlscorr.f s 3 p3 ! ! ! dp3
fun/pl [a]*exp(-x/[c])+[b]*exp(x/[c]) -10 10 s
a=p3(1)
b=p3(2)
c=p3(3)
fun/pl [a]*exp(-x/[c])+[b]*exp(x/[c]) -8 8 s
sigma p3 = abs(p3)
ve/write p3,dp3 corrwls_counter[nc]_v[corrz].txt '2e15.6'
read x
return



macro wlsspect nc=1 ncz=0
*
idh=[nc]+10*8+100*[ncz]+100000
*
hi/del [idh]
hi/file 20 mhad2011-4_profilex_ee_sim_shifter2_chit.his
hrin [idh]
close 20
hi/copy [idh] 1000
exec ../SepPar/sp#brn 1000 10 100
exec ../SepPar/sp#brn 1000 20 900
*
hi/del [idh]
hi/file 20 mhad2011-4_profile_sin_new_et_chit.his
hrin [idh]
close 20
hi/copy [idh] 1000
exec ../SepPar/sp#brn 1000 10 200
exec ../SepPar/sp#brn 1000 20 800
*
nsim=$hinfo(1100,'events')
nexp=$hinfo(1200,'events')
*exec hsigma @1300 = @1200/vsum(@1200)*vsum(@1100)
hi/op/add 1200 1200 1300 $sigma([nsim]/[nexp]) 0
hi/pl 1100 e
hi/pl 1300 s
hi/copy 1300 1400
*
l=20.01
r=40.01
ve/cre   p3(3) r  100 25 5
ve/cre pmin(3) r    0 20 3
ve/cre pmax(3) r 1000 35 8
ve/cre   s3(3) r    1  1 0
hi/fit 1800([l]:[r]) g b 3 p3 s3 pmin pmax dp3
ve/cre dp3(3) r
do i=1,3
  hi/fit 1800([l]:[r]) g b 3 p3 ! pmin pmax dp3
  l=$sigma(p3(2)-2.5*abs(p3(3)))
  r=$sigma(p3(2)+2.5*abs(p3(3)))
enddo
l=$sigma(p3(2)-1.5*abs(p3(3)))
r=$sigma(p3(2)+1.5*abs(p3(3)))
hi/fit 1800([l]:[r]) g b 3 p3 ! pmin pmax dp3
ve/copy p3 p3r
do i=1,2
  l=$sigma(p3r(2)-0.5*abs(p3r(3)))
  r=$sigma(p3r(2)+2.0*abs(p3r(3)))
  hi/fit 1800([l]:[r]) g ! 3 p3r
enddo
*
l=20.01
r=40.01
ve/cre  p3s(3) r 100 25 5
ve/cre  s3s(3) r   1  1 0
hi/fit 1900([l]:[r]) g b 3 p3s s3s pmin pmax dp3s
ve/cre dp3s(3) r
do i=1,3
  hi/fit 1900([l]:[r]) g b 3 p3s ! pmin pmax dp3s
  l=$sigma(p3s(2)-2.5*abs(p3s(3)))
  r=$sigma(p3s(2)+2.5*abs(p3s(3)))
enddo
l=$sigma(p3s(2)-1.5*abs(p3s(3)))
r=$sigma(p3s(2)+1.5*abs(p3s(3)))
hi/fit 1900([l]:[r]) g b 3 p3s ! pmin pmax dp3s
ve/copy p3s p3sr
do i=1,2
  l=$sigma(p3sr(2)-0.5*abs(p3sr(3)))
  r=$sigma(p3sr(2)+2.0*abs(p3sr(3)))
  hi/fit 1900([l]:[r]) g ! 3 p3sr
enddo
*
ae=$sigma(p3r(2)+abs(p3r(3)))
as=$sigma(p3sr(2)+abs(p3sr(3)))
mess Ae=[ae] As=[as] r=$sigma([as]/[ae]) dA=Ae-As=$sigma([ae]-[as])
*
goto 1
name=spect
hdiff=-1
ve/cre p3s(3) r 100 25 5
dp3=100
do i=20,40
  fname=/work/users/konctbel/snd2k/R005-999/[name]_[i]pe.txt
  if ($fexist([fname]).eq.1) then
    nt/cre 1 ! 9 ! ! a1 a2 a3 a4 a5 a6 a7 a8 a9
    nt/read 1 [fname]
    nt/pl 1.a[nc] idh=1400
    l=15.01
    r=50.01
    ve/cre p3s(3) r 100 25 5
    do j=1,3
      hi/fit 1400([l]:[r]) g ! 3 p3s
      l=$sigma(p3s(2)-1.5*abs(p3s(3)))
      r=$sigma(p3s(2)+1*abs(p3s(3)))
    enddo
    dp3i=$sigma(abs(p3(2)-p3s(2)))
    if ([dp3i].lt.[dp3]) then
      dp3=[dp3i]
      hi/copy 1400 1500
      ve/copy p3s p3sb
      ib=[i]
    endif
*    hdiffi=$call('mhdiff(1200,1400)')
*    diff 1200 1400
*    mess [hdiffi]
*    if ([hdiffi].gt.[hdiff]) then
*      hi/copy 1400 1500
*      hdiff=[hdiffi]
*      mess [hdiff]
*    endif
  endif
enddo
nsim=$hinfo(1100,'events')
nexp=$hinfo(1500,'events')
hi/op/add 1500 1500 1600 $sigma([nsim]/[nexp]) 0
set hcol 1 
hi/pl 1100 e
set hcol 2
hi/pl 1300 s
set hcol 4
hi/pl 1600 s
mess iBest=[ib]
1:
set hcol 1
hi/pl 1100 e
hi/del 1700
exec hsigma @1700 = @1300
hi/pl 1700 s 
return


macro wlsspectexpxx map=1
do nc=1,9
  exec mapcal#wlsspectexpx [map] [nc]
enddo
return


macro wlsspectexpx map=1 nc=1
*
ncz=10
zmin=-10
zmax=9
*
ve/cre  asz([ncz]) r 
ve/cre dasz([ncz]) r
ve/cre  aszs([ncz]) r 
ve/cre daszs([ncz]) r
ve/cre  zsz([ncz]) r 
ve/cre dzsz([ncz]) r
do i=1,[ncz]
  exec mapcal#wlsspectexp [map] [nc] [i]
  ve/inp  asz([i]) $sigma(p3(2))
  ve/inp dasz([i]) $sigma(dp3(2))
  ve/inp  zsz([i]) $sigma([zmin]+([i]-0.5)*([zmax]-[zmin])/[ncz])
  exec mapcal#wlsspectsim [map] [nc] [i] 2
  ve/inp  aszs([i]) $sigma(p3(2))
  ve/inp daszs([i]) $sigma(dp3(2))
enddo
set pmci 1
* exec $PER/s#vpl asz dasz zsz dzsz sz=0.1 iatt=20
n1=2
n2=$sigma(min(8,[ncz]))
ve/del aszr,daszr,zszr,dzszr
ve/copy  asz([n1]:[n2])  aszr
ve/copy dasz([n1]:[n2]) daszr
ve/copy  zsz([n1]:[n2])  zszr
ve/copy dzsz([n1]:[n2]) dzszr
ve/cre wls(3) r 0.5 0.5 1000000
ve/cre p3(3) r 10 10 50
ve/cre dp3(3) r
do i=1,2
  ve/fit zszr aszr daszr wlscorr.f s 3 p3 ! ! ! dp3
enddo
* exec $PER/s#vpl asz dasz zsz dzsz sz=0.1
* exec $PER/s#vpl aszs daszs zsz dzsz sz=0.1 iatt=24 o=s
ve/fit zszr aszr daszr wlscorr.f s 3 p3 ! ! ! dp3
gl/imp mname
sigma p3=abs(p3)
*ve/write p3,dp3 [mname]/[mname]_corrwls_counter[nc]_v[ver].txt '2e15.6'
*
txt=[mname]
exec $PER/s#tf 0.05 0.9 [txt]
txt=Counter [nc]
exec $PER/s#tf 0.05 0.8 [txt]
*
atitle 'z?r!,  cm' 'A?wls!, pe'
*exec save [mname]/[mname]_corrwls_counter[nc]_v[ver].eps f
exec save [mname]/[mname]_corrwls_counter[nc]_full.eps f
return


macro wlsspectexp map=1 nc=1 ncz=0
*
if ([map].eq.1) then
  gl/cre mname mapcal1
  gl/cre nh1 15
  gl/cre nh2 15
  gl/cre nb 14
endif
if ([map].eq.2) then
  gl/cre mname mapcal2
  gl/cre nh1 16
  gl/cre nh2 16
  gl/cre nb 14
endif
if ([map].eq.3) then
  gl/cre mname mapcal3
  gl/cre nh1 17
  gl/cre nh2 17
  gl/cre nb 10
endif
if ([map].eq.4) then
  gl/cre mname mapcal4
  gl/cre nh1 18
  gl/cre nh2 18
  gl/cre nb 14
endif
*
idh=[nc]+10*8+100*[ncz]+100000
*
hi/del [idh]
*hi/file 20 [mname]/[mname]_profilex_exp_pds.his
hi/file 20 [mname]_profilex_exp_full.his
hrin [idh]
close 20
hi/copy [idh] 1000
exec ../SepPar/sp#brn 1000 10 200
exec ../SepPar/sp#brn 1000 20 800
*
l=20.01
r=40.01
ve/cre   p3(3) r  100 25 5
ve/cre pmin(3) r    0 10 2
ve/cre pmax(3) r 10000 65 15
ve/cre   s3(3) r    1  1 0
hi/fit 1800([l]:[r]) g b 3 p3 s3 pmin pmax dp3
ve/cre dp3(3) r
lmin=10.01
lmin=8.01
ve/cre   s3(3) r    1  1 0.3
do i=1,3
  hi/fit 1800([l]:[r]) g b 3 p3 ! pmin pmax dp3
  l=$sigma(max([lmin],p3(2)-2.5*abs(p3(3))))
  r=$sigma(p3(2)+2.5*abs(p3(3)))
enddo
l=$sigma(max([lmin],p3(2)-1.5*abs(p3(3))))
r=$sigma(p3(2)+1.5*abs(p3(3)))
*hi/fit 1800([l]:[r]) g b 3 p3 ! pmin pmax dp3
ve/copy p3 p3r
do i=1,3
  l=$sigma(max([lmin],p3(2)-1.5*abs(p3(3))))
  r=$sigma(p3(2)+2.*abs(p3(3)))
  hi/fit 1800([l]:[r]) g b 3 p3 ! ! ! dp3
enddo
*
nevt=$hinfo(1800,'events')
nevt=[nevt]/2
ve/cre p8(8) r [nevt] 4 2 0.1 [nevt] 20 10 0.1
ve/cre pmin(8) r 0 1 1 0 0 15 5 0
ve/cre pmax(8) r 1000 7 5 0.2 1000 50 20 0.2
ve/cre s8(8) r 1 1 1 0 1 1 1 0
*hi/fit 1200 fun32.for b 8 p8 s8 pmin pmax
*
gl/cre rime $sigma(p3(2)+1.17741*abs(p3(3)))
return


macro wlsspectsim map=1 nc=1 ncz=0 ver=0
*
if ([map].eq.1) then
  gl/cre mname mapcal1
  gl/cre nh1 15
  gl/cre nh2 15
  gl/cre nb 14
endif
if ([map].eq.2) then
  gl/cre mname mapcal2
  gl/cre nh1 16
  gl/cre nh2 16
  gl/cre nb 14
endif
if ([map].eq.3) then
  gl/cre mname mapcal3
  gl/cre nh1 17
  gl/cre nh2 17
  gl/cre nb 10
endif
if ([map].eq.4) then
  gl/cre mname mapcal4
  gl/cre nh1 18
  gl/cre nh2 18
  gl/cre nb 14
endif
*
idh=[nc]+10*8+100*[ncz]+100000
*
hi/del [idh]
hi/file 20 [mname]_profilex_sim_v[ver].his
hrin [idh]
close 20
hi/copy [idh] 1000
exec ../SepPar/sp#brn 1000 10 200
exec ../SepPar/sp#brn 1000 20 800
*
l=20.01
r=40.01
ve/cre   p3(3) r  100 25 5
ve/cre pmin(3) r    0 10 2
ve/cre pmax(3) r 1000 35 15
ve/cre   s3(3) r    1  1 0
hi/fit 1800([l]:[r]) g b 3 p3 s3 pmin pmax dp3
ve/cre dp3(3) r
lmin=10.01
lmin=8.01
ve/cre   s3(3) r    1  1 0.3
do i=1,3
  hi/fit 1800([l]:[r]) g b 3 p3 ! pmin pmax dp3
  l=$sigma(max([lmin],p3(2)-2.5*abs(p3(3))))
  r=$sigma(p3(2)+2.5*abs(p3(3)))
enddo
l=$sigma(max([lmin],p3(2)-1.5*abs(p3(3))))
r=$sigma(p3(2)+1.5*abs(p3(3)))
*hi/fit 1800([l]:[r]) g b 3 p3 ! pmin pmax dp3
ve/copy p3 p3r
do i=1,3
  l=$sigma(max([lmin],p3(2)-1.*abs(p3(3))))
  r=$sigma(p3(2)+2.5*abs(p3(3)))
  hi/fit 1800([l]:[r]) g b 3 p3 ! ! ! dp3
enddo
*
nevt=$hinfo(1800,'events')
nevt=[nevt]/2
ve/cre p8(8) r [nevt] 4 2 0.1 [nevt] 20 10 0.1
ve/cre pmin(8) r 0 1 1 0 0 15 5 0
ve/cre pmax(8) r 1000 7 5 0.2 1000 50 20 0.2
ve/cre s8(8) r 1 1 1 0 1 1 1 0
*hi/fit 1200 fun32.for b 8 p8 s8 pmin pmax
*
gl/cre rims $sigma(p3(2)+1.17741*abs(p3(3)))
return


macro wlscomp map=1 nc=1 ncz=0
exec mapcal#wlsspectsim [map] [nc] [ncz] 1
hi/del 1801,1802
hi/copy 1800 1801
exec mapcal#wlsspectexp [map] [nc] [ncz]
exec hsigma @1802 = @1800/vsum(@1800)*vsum(@1801)
exec hsigma %1802 = %1800/vsum(@1800)*vsum(@1801)
set hcol 1
set ksiz 0.15
set pmci 2
hi/pl 1802(:100.)
set pmci 4
hi/pl 1801 s
gl/imp rime
gl/imp rims
line [rime] $GRAFINFO('WNYMIN') [rime] $GRAFINFO('WNYMAX')
line [rims] $GRAFINFO('WNYMIN') [rims] $GRAFINFO('WNYMAX')
*
read x
n=$hinfo(1800,'xbins')
ve/cre nx([n]) r
ve/cre n0([n]) r
hi/get/cont 1800 nx
hi/get/cont 1801 n0
hi/copy 1802 100
sigma qx = n0/nx
sigma dqx=sqrt(n0*(1-n0/nx))/nx
hi/put/cont 100 qx
hi/put/err  100 dqx
hi/pl 100(10.:)
ve/cre p3(3) r 0.1 60 5
hi/fit 100(10.:) erf0q.f s 3 p3
return


macro a1shift
ve/cre a1s(9) r 90.043 77.241 105.399 76.985 97.676 114.848 94.784 94.174  98.95
fname=accled_amp1pe_mc.txt
ve/del zvo,runo,aeo1,aeo2,aeo3,aeo4,aeo5,aeo6,aeo7,aeo8,aeo9
ve/read zvo,runo,aeo1,aeo2,aeo3,aeo4,aeo5,aeo6,aeo7,aeo8,aeo9 [fname] '(2f10.1,9f15.6)'
*
n=0
n=[n]+1; nrf[n]=7842; 
n=[n]+1; nrf[n]=11285;
n=[n]+1; nrf[n]=13847;
n=[n]+1; nrf[n]=16699;
*n=[n]+1; nrf[n]=17869;
*
do nc=1,9
  ve/cre ae[nc]([n]) r
  ve/cre zv[nc]([n]) r
enddo
ve/cre run([n]) r 
*
do i=1,[n]
  ve/del corra
  fcorra=AmplitudeCorrection_mapcal[i].txt
  if ($fexist([fcorra]).eq.1) then
    ve/read corra [fcorra]
  else
    ve/cre corra(9) r 9*1
  endif
  ve/inp run([i]) [nrf[i]]
  do nc=1,-9
    if ([i].le.2) then
      exec mapcal#wlscomp [i] [nc] 0
      gl/imp rime
      gl/imp rims
      ve/inp ae[nc]([i]) $sigma(aeo[nc]([i])*[rims]/[rime]/corra([nc]))
    else
      ve/inp ae[nc]([i]) $sigma(aeo[nc]([i])/corra([nc]))
    endif
  enddo
enddo
ve/cre zv([n]) r
*
n=0
n=[n]+1; v[n]=ae; f[n]=amp1pe
n=[n]+1; v[n]=zv; f[n]=u 
n=[n]+1; v[n]=zv; f[n]=a 
n=[n]+1; v[n]=zv; f[n]=b 
n=[n]+1; v[n]=zv; f[n]=c 
n=[n]+1; v[n]=zv; f[n]=eff1pe
n=[n]+1; v[n]=zv; f[n]=u0
n=[n]+1; v[n]=zv; f[n]=s0
n=[n]+1; v[n]=zv; f[n]=pds
n=[n]+1; v[n]=zv; f[n]=tmin
n=[n]+1; v[n]=zv; f[n]=tmax
*
do v=2,[n]
*do v=5,5
  fname=accled_[f[v]]_mc_v1.txt
  if ($fexist([fname])) then
    shell rm [fname]
  endif
  for/file 20 [fname] ! n
  close 20
  ve/write zv,run,[v[v]]1,[v[v]]2,[v[v]]3,[v[v]]4,[v[v]]5,[v[v]]6,[v[v]]7,[v[v]]8,[v[v]]9 [fname] '(2f10.1,9f15.6)'
enddo

return

macro rho_2012_list
shell ls /online/gridspool/proc/RHO_2012-0/*.col.gz > rho_2012.list
ve/del runs
ve/read runs rho_2012.list '37X,f9.0'
n=$vlen(runs)
ve/cre beams([n]) r
do i=1,[n]
  exec mapcal#ndays $sigma(runs([i]))
  gl/imp ebeam
  ve/inp beams([i]) [ebeam]
enddo
ve/write runs,beams rho_2012_beam.list '2f10.3'
return



macro omeg2012list
n=0
n=[n]+1; file[n]=omega2012_385.0_p0.hbook
n=[n]+1; file[n]=omega2012_385.0_p1.hbook
n=[n]+1; file[n]=omega2012_385.0_p2.hbook
n=[n]+1; file[n]=omega2012_385.0_p3.hbook
n=[n]+1; file[n]=omega2012_387.0_p0.hbook
n=[n]+1; file[n]=omega2012_387.0_p1.hbook
n=[n]+1; file[n]=omega2012_387.0_p2.hbook
n=[n]+1; file[n]=omega2012_387.0_p3.hbook
n=[n]+1; file[n]=omega2012_388.0_p0.hbook
n=[n]+1; file[n]=omega2012_388.0_p1.hbook
n=[n]+1; file[n]=omega2012_388.0_p2.hbook
n=[n]+1; file[n]=omega2012_388.0_p3.hbook
n=[n]+1; file[n]=omega2012_388.0_p4.hbook
n=[n]+1; file[n]=omega2012_389.0_p0.hbook
n=[n]+1; file[n]=omega2012_389.0_p1.hbook
n=[n]+1; file[n]=omega2012_389.0_p2.hbook
n=[n]+1; file[n]=omega2012_389.0_p3.hbook
n=[n]+1; file[n]=omega2012_389.0_p4.hbook
n=[n]+1; file[n]=omega2012_389.0_p5.hbook
n=[n]+1; file[n]=omega2012_389.5_p0.hbook
n=[n]+1; file[n]=omega2012_390.0_p0.hbook
n=[n]+1; file[n]=omega2012_390.0_p1.hbook
n=[n]+1; file[n]=omega2012_390.0_p2.hbook
n=[n]+1; file[n]=omega2012_390.0_p3.hbook
n=[n]+1; file[n]=omega2012_391.0_p0.hbook
n=[n]+1; file[n]=omega2012_391.0_p1.hbook
n=[n]+1; file[n]=omega2012_391.0_p2.hbook
n=[n]+1; file[n]=omega2012_392.0_p0.hbook
n=[n]+1; file[n]=omega2012_392.0_p1.hbook
n=[n]+1; file[n]=omega2012_392.0_p2.hbook
n=[n]+1; file[n]=omega2012_394.0_p0.hbook
n=[n]+1; file[n]=omega2012_394.0_p1.hbook
n=[n]+1; file[n]=omega2012_394.0_p2.hbook
n=[n]+1; file[n]=omega2012_394.0_p3.hbook
n=[n]+1; file[n]=omega2012_394.0_p4.hbook
n=[n]+1; file[n]=omega2012_395.0_p0.hbook
dir=/work/users/konctbel/exp/OMEG2012-0/col/
f1=omeg2012-0_n1.13_profile.list
if ($fexist([f1]).eq.1) then
  shell rm [f1]
endif
for/file 20 [f1] n
close 20
fmess 0 [f1]
fmess [dir] [f1]
fmess [dir] [f1]
f2=omeg2012-0_n1.05_profile.list
if ($fexist([f2]).eq.1) then
  shell rm [f2]
endif
for/file 20 [f2] n
close 20
fmess 0 [f2]
fmess [dir] [f2]
fmess [dir] [f2]
runx=13847
1d 100 ! 10000 10000.5 20000.5
ve/cre irun0(10000) r
ve/cre ind0(10000) r
sigma ind0 = array(10000,10001#20000)
do i=1,[n]
  hi/file 20 [dir][file[i]]
  nt/pl 1.run idh=100
  hi/get/cont 100 irun0
  ve/del ind,irun
  sigma ind=order(ind0,-irun0)
  sigma irun=order(irun0,-irun0)
  exec mapcal#vecut irun
  m=$vlen(irun)
  exec mapcal#vecut ind [m]
  rmin=$sigma(vmin(ind))
  rmax=$sigma(vmax(ind))
  mean=$hinfo(100,'mean')
  rms=$hinfo(100,'rms')
  l=$sigma([mean]-3*max([rms],3))
  r=$sigma([mean]+3*max([rms],3))
  if ($sigma(int([l])).eq.[l]) then
    l=[l]-0.5
  endif
  if ($sigma(int([r])).eq.[r]) then
    r=[r]+0.5
  endif
  hi/pl 100([l]:[r])
  line [rmin] $GRAFINFO('WNYMIN') [rmin] $GRAFINFO('WNYMAX')
  line [rmax] $GRAFINFO('WNYMIN') [rmax] $GRAFINFO('WNYMAX')
  if (([runx].gt.[rmin]).and.([runx].le.[rmax])) then
    mess wery strange!!! -> [file[i]]: [runx] in {[rmin],[rmax]}
  endif
  if ([runx].gt.[rmax]) then
    fmess [file[i]] [f1]
    mess [rmin] [rmax] [f1]
  endif
  if ([runx].le.[rmin]) then
    fmess [file[i]] [f2]
    mess [rmin] [rmax] [f2]
  endif
  close 20
*  read x
enddo
return


macro rho_2012list
n=0
n=[n]+1; file[n]=rho2012_160.0_p0.hbook         
n=[n]+1; file[n]=rho2012_160.0_p10.hbook        
n=[n]+1; file[n]=rho2012_160.0_p11.hbook        
n=[n]+1; file[n]=rho2012_160.0_p12.hbook        
n=[n]+1; file[n]=rho2012_160.0_p13.hbook        
n=[n]+1; file[n]=rho2012_160.0_p14.hbook        
n=[n]+1; file[n]=rho2012_160.0_p15.hbook        
n=[n]+1; file[n]=rho2012_160.0_p16.hbook        
n=[n]+1; file[n]=rho2012_160.0_p17.hbook        
n=[n]+1; file[n]=rho2012_160.0_p18.hbook        
n=[n]+1; file[n]=rho2012_160.0_p1.hbook         
n=[n]+1; file[n]=rho2012_160.0_p2.hbook         
n=[n]+1; file[n]=rho2012_160.0_p3.hbook         
n=[n]+1; file[n]=rho2012_160.0_p4.hbook         
n=[n]+1; file[n]=rho2012_160.0_p5.hbook         
n=[n]+1; file[n]=rho2012_160.0_p6.hbook         
n=[n]+1; file[n]=rho2012_160.0_p7.hbook         
n=[n]+1; file[n]=rho2012_160.0_p8.hbook         
n=[n]+1; file[n]=rho2012_160.0_p9.hbook         
n=[n]+1; file[n]=rho2012_170.0_p0.hbook         
n=[n]+1; file[n]=rho2012_170.0_p10.hbook        
n=[n]+1; file[n]=rho2012_170.0_p11.hbook        
n=[n]+1; file[n]=rho2012_170.0_p12.hbook        
n=[n]+1; file[n]=rho2012_170.0_p13.hbook        
n=[n]+1; file[n]=rho2012_170.0_p14.hbook        
n=[n]+1; file[n]=rho2012_170.0_p15.hbook        
n=[n]+1; file[n]=rho2012_170.0_p1.hbook         
n=[n]+1; file[n]=rho2012_170.0_p2.hbook         
n=[n]+1; file[n]=rho2012_170.0_p3.hbook         
n=[n]+1; file[n]=rho2012_170.0_p4.hbook         
n=[n]+1; file[n]=rho2012_170.0_p5.hbook         
n=[n]+1; file[n]=rho2012_170.0_p6.hbook         
n=[n]+1; file[n]=rho2012_170.0_p7.hbook         
n=[n]+1; file[n]=rho2012_170.0_p8.hbook         
n=[n]+1; file[n]=rho2012_170.0_p9.hbook         
n=[n]+1; file[n]=rho2012_180.0_p0.hbook         
n=[n]+1; file[n]=rho2012_180.0_p10.hbook        
n=[n]+1; file[n]=rho2012_180.0_p11.hbook        
n=[n]+1; file[n]=rho2012_180.0_p12.hbook        
n=[n]+1; file[n]=rho2012_180.0_p13.hbook        
n=[n]+1; file[n]=rho2012_180.0_p14.hbook        
n=[n]+1; file[n]=rho2012_180.0_p15.hbook        
n=[n]+1; file[n]=rho2012_180.0_p16.hbook        
n=[n]+1; file[n]=rho2012_180.0_p17.hbook        
n=[n]+1; file[n]=rho2012_180.0_p18.hbook        
n=[n]+1; file[n]=rho2012_180.0_p19.hbook        
n=[n]+1; file[n]=rho2012_180.0_p1.hbook         
n=[n]+1; file[n]=rho2012_180.0_p20.hbook        
n=[n]+1; file[n]=rho2012_180.0_p2.hbook         
n=[n]+1; file[n]=rho2012_180.0_p3.hbook         
n=[n]+1; file[n]=rho2012_180.0_p4.hbook         
n=[n]+1; file[n]=rho2012_180.0_p5.hbook         
n=[n]+1; file[n]=rho2012_180.0_p6.hbook         
n=[n]+1; file[n]=rho2012_180.0_p7.hbook         
n=[n]+1; file[n]=rho2012_180.0_p8.hbook         
n=[n]+1; file[n]=rho2012_180.0_p9.hbook         
n=[n]+1; file[n]=rho2012_190.0_p0.hbook         
n=[n]+1; file[n]=rho2012_190.0_p1.hbook         
n=[n]+1; file[n]=rho2012_190.0_p2.hbook         
n=[n]+1; file[n]=rho2012_190.0_p3.hbook         
n=[n]+1; file[n]=rho2012_190.0_p4.hbook         
n=[n]+1; file[n]=rho2012_200.0_p0.hbook         
n=[n]+1; file[n]=rho2012_200.0_p1.hbook         
n=[n]+1; file[n]=rho2012_200.0_p2.hbook         
n=[n]+1; file[n]=rho2012_200.0_p3.hbook         
n=[n]+1; file[n]=rho2012_200.0_p4.hbook         
n=[n]+1; file[n]=rho2012_200.0_p5.hbook         
n=[n]+1; file[n]=rho2012_200.0_p6.hbook         
n=[n]+1; file[n]=rho2012_210.0_p0.hbook         
n=[n]+1; file[n]=rho2012_210.0_p1.hbook         
n=[n]+1; file[n]=rho2012_210.0_p2.hbook         
n=[n]+1; file[n]=rho2012_210.0_p3.hbook         
n=[n]+1; file[n]=rho2012_210.0_p4.hbook         
n=[n]+1; file[n]=rho2012_210.0_p5.hbook         
n=[n]+1; file[n]=rho2012_220.0_p0.hbook         
n=[n]+1; file[n]=rho2012_220.0_p1.hbook         
n=[n]+1; file[n]=rho2012_220.0_p2.hbook         
n=[n]+1; file[n]=rho2012_220.0_p3.hbook         
n=[n]+1; file[n]=rho2012_220.0_p4.hbook         
n=[n]+1; file[n]=rho2012_220.0_p5.hbook         
n=[n]+1; file[n]=rho2012_220.0_p6.hbook         
n=[n]+1; file[n]=rho2012_220.0_p7.hbook         
n=[n]+1; file[n]=rho2012_230.0_p0.hbook         
n=[n]+1; file[n]=rho2012_230.0_p1.hbook         
n=[n]+1; file[n]=rho2012_230.0_p2.hbook         
n=[n]+1; file[n]=rho2012_230.0_p3.hbook         
n=[n]+1; file[n]=rho2012_230.0_p4.hbook         
n=[n]+1; file[n]=rho2012_230.0_p5.hbook         
n=[n]+1; file[n]=rho2012_230.0_p6.hbook         
n=[n]+1; file[n]=rho2012_230.0_p7.hbook         
n=[n]+1; file[n]=rho2012_240.0_p0.hbook         
n=[n]+1; file[n]=rho2012_240.0_p1.hbook         
n=[n]+1; file[n]=rho2012_240.0_p2.hbook         
n=[n]+1; file[n]=rho2012_240.0_p3.hbook         
n=[n]+1; file[n]=rho2012_240.0_p4.hbook         
n=[n]+1; file[n]=rho2012_240.0_p5.hbook         
n=[n]+1; file[n]=rho2012_240.0_p6.hbook         
n=[n]+1; file[n]=rho2012_240.0_p7.hbook         
n=[n]+1; file[n]=rho2012_240.0_p8.hbook         
n=[n]+1; file[n]=rho2012_250.0_p0.hbook         
n=[n]+1; file[n]=rho2012_250.0_p10.hbook        
n=[n]+1; file[n]=rho2012_250.0_p11.hbook        
n=[n]+1; file[n]=rho2012_250.0_p12.hbook        
n=[n]+1; file[n]=rho2012_250.0_p13.hbook        
n=[n]+1; file[n]=rho2012_250.0_p14.hbook        
n=[n]+1; file[n]=rho2012_250.0_p15.hbook        
n=[n]+1; file[n]=rho2012_250.0_p16.hbook        
n=[n]+1; file[n]=rho2012_250.0_p1.hbook         
n=[n]+1; file[n]=rho2012_250.0_p2.hbook         
n=[n]+1; file[n]=rho2012_250.0_p3.hbook         
n=[n]+1; file[n]=rho2012_250.0_p4.hbook         
n=[n]+1; file[n]=rho2012_250.0_p5.hbook         
n=[n]+1; file[n]=rho2012_250.0_p6.hbook         
n=[n]+1; file[n]=rho2012_250.0_p7.hbook         
n=[n]+1; file[n]=rho2012_250.0_p8.hbook         
n=[n]+1; file[n]=rho2012_250.0_p9.hbook         
n=[n]+1; file[n]=rho2012_260.0_p0.hbook         
n=[n]+1; file[n]=rho2012_260.0_p1.hbook         
n=[n]+1; file[n]=rho2012_260.0_p2.hbook         
n=[n]+1; file[n]=rho2012_270.0_p0.hbook         
n=[n]+1; file[n]=rho2012_270.0_p1.hbook         
n=[n]+1; file[n]=rho2012_270.0_p2.hbook         
n=[n]+1; file[n]=rho2012_270.0_p3.hbook         
n=[n]+1; file[n]=rho2012_270.0_p4.hbook         
n=[n]+1; file[n]=rho2012_270.0_p5.hbook         
n=[n]+1; file[n]=rho2012_280.0_p0.hbook         
n=[n]+1; file[n]=rho2012_280.0_p1.hbook         
n=[n]+1; file[n]=rho2012_280.0_p2.hbook         
n=[n]+1; file[n]=rho2012_280.0_p3.hbook         
n=[n]+1; file[n]=rho2012_290.0_p0.hbook         
n=[n]+1; file[n]=rho2012_290.0_p1.hbook         
n=[n]+1; file[n]=rho2012_290.0_p2.hbook         
n=[n]+1; file[n]=rho2012_300.0_p0.hbook         
n=[n]+1; file[n]=rho2012_300.0_p1.hbook         
n=[n]+1; file[n]=rho2012_300.0_p2.hbook         
n=[n]+1; file[n]=rho2012_300.0_p3.hbook         
n=[n]+1; file[n]=rho2012_310.0_p0.hbook         
n=[n]+1; file[n]=rho2012_310.0_p1.hbook         
n=[n]+1; file[n]=rho2012_310.0_p2.hbook         
n=[n]+1; file[n]=rho2012_320.0_p0.hbook         
n=[n]+1; file[n]=rho2012_320.0_p1.hbook         
n=[n]+1; file[n]=rho2012_320.0_p2.hbook         
n=[n]+1; file[n]=rho2012_330.0_p0.hbook         
n=[n]+1; file[n]=rho2012_330.0_p10.hbook        
n=[n]+1; file[n]=rho2012_330.0_p11.hbook        
n=[n]+1; file[n]=rho2012_330.0_p12.hbook        
n=[n]+1; file[n]=rho2012_330.0_p13.hbook        
n=[n]+1; file[n]=rho2012_330.0_p14.hbook        
n=[n]+1; file[n]=rho2012_330.0_p15.hbook        
n=[n]+1; file[n]=rho2012_330.0_p16.hbook        
n=[n]+1; file[n]=rho2012_330.0_p17.hbook        
n=[n]+1; file[n]=rho2012_330.0_p18.hbook        
n=[n]+1; file[n]=rho2012_330.0_p19.hbook        
n=[n]+1; file[n]=rho2012_330.0_p1.hbook         
n=[n]+1; file[n]=rho2012_330.0_p20.hbook        
n=[n]+1; file[n]=rho2012_330.0_p21.hbook        
n=[n]+1; file[n]=rho2012_330.0_p22.hbook        
n=[n]+1; file[n]=rho2012_330.0_p23.hbook        
n=[n]+1; file[n]=rho2012_330.0_p24.hbook        
n=[n]+1; file[n]=rho2012_330.0_p25.hbook        
n=[n]+1; file[n]=rho2012_330.0_p26.hbook        
n=[n]+1; file[n]=rho2012_330.0_p2.hbook         
n=[n]+1; file[n]=rho2012_330.0_p3.hbook         
n=[n]+1; file[n]=rho2012_330.0_p4.hbook         
n=[n]+1; file[n]=rho2012_330.0_p5.hbook         
n=[n]+1; file[n]=rho2012_330.0_p6.hbook         
n=[n]+1; file[n]=rho2012_330.0_p7.hbook         
n=[n]+1; file[n]=rho2012_330.0_p8.hbook         
n=[n]+1; file[n]=rho2012_330.0_p9.hbook         
n=[n]+1; file[n]=rho2012_340.0_p0.hbook         
n=[n]+1; file[n]=rho2012_340.0_p1.hbook         
n=[n]+1; file[n]=rho2012_340.0_p2.hbook         
n=[n]+1; file[n]=rho2012_340.0_p3.hbook         
n=[n]+1; file[n]=rho2012_340.0_p4.hbook         
n=[n]+1; file[n]=rho2012_350.0_p0.hbook         
n=[n]+1; file[n]=rho2012_350.0_p1.hbook         
n=[n]+1; file[n]=rho2012_350.0_p2.hbook         
n=[n]+1; file[n]=rho2012_350.0_p3.hbook         
n=[n]+1; file[n]=rho2012_360.0_p0.hbook         
n=[n]+1; file[n]=rho2012_360.0_p1.hbook         
n=[n]+1; file[n]=rho2012_360.0_p2.hbook         
n=[n]+1; file[n]=rho2012_360.0_p3.hbook         
n=[n]+1; file[n]=rho2012_360.0_p4.hbook         
n=[n]+1; file[n]=rho2012_364.0_p0.hbook         
n=[n]+1; file[n]=rho2012_364.0_p1.hbook         
n=[n]+1; file[n]=rho2012_364.0_p2.hbook         
n=[n]+1; file[n]=rho2012_364.0_p3.hbook         
n=[n]+1; file[n]=rho2012_368.0_p0.hbook         
n=[n]+1; file[n]=rho2012_368.0_p10.hbook        
n=[n]+1; file[n]=rho2012_368.0_p11.hbook        
n=[n]+1; file[n]=rho2012_368.0_p1.hbook         
n=[n]+1; file[n]=rho2012_368.0_p2.hbook         
n=[n]+1; file[n]=rho2012_368.0_p3.hbook         
n=[n]+1; file[n]=rho2012_368.0_p4.hbook         
n=[n]+1; file[n]=rho2012_368.0_p5.hbook         
n=[n]+1; file[n]=rho2012_368.0_p6.hbook         
n=[n]+1; file[n]=rho2012_368.0_p7.hbook         
n=[n]+1; file[n]=rho2012_368.0_p8.hbook         
n=[n]+1; file[n]=rho2012_368.0_p9.hbook         
n=[n]+1; file[n]=rho2012_370.0_p0.hbook         
n=[n]+1; file[n]=rho2012_370.0_p1.hbook         
n=[n]+1; file[n]=rho2012_370.0_p2.hbook         
n=[n]+1; file[n]=rho2012_372.0_p0.hbook         
n=[n]+1; file[n]=rho2012_372.0_p1.hbook         
n=[n]+1; file[n]=rho2012_374.0_p0.hbook         
n=[n]+1; file[n]=rho2012_374.0_p1.hbook         
n=[n]+1; file[n]=rho2012_374.0_p2.hbook         
n=[n]+1; file[n]=rho2012_376.0_p0.hbook         
n=[n]+1; file[n]=rho2012_376.0_p1.hbook         
n=[n]+1; file[n]=rho2012_376.0_p2.hbook         
n=[n]+1; file[n]=rho2012_376.0_p3.hbook         
n=[n]+1; file[n]=rho2012_378.0_p0.hbook         
n=[n]+1; file[n]=rho2012_378.0_p1.hbook         
n=[n]+1; file[n]=rho2012_378.0_p2.hbook         
n=[n]+1; file[n]=rho2012_380.0_p0.hbook         
n=[n]+1; file[n]=rho2012_380.0_p10.hbook        
n=[n]+1; file[n]=rho2012_380.0_p11.hbook        
n=[n]+1; file[n]=rho2012_380.0_p12.hbook        
n=[n]+1; file[n]=rho2012_380.0_p13.hbook        
n=[n]+1; file[n]=rho2012_380.0_p14.hbook        
n=[n]+1; file[n]=rho2012_380.0_p15.hbook        
n=[n]+1; file[n]=rho2012_380.0_p16.hbook        
n=[n]+1; file[n]=rho2012_380.0_p1.hbook         
n=[n]+1; file[n]=rho2012_380.0_p2.hbook         
n=[n]+1; file[n]=rho2012_380.0_p3.hbook         
n=[n]+1; file[n]=rho2012_380.0_p4.hbook         
n=[n]+1; file[n]=rho2012_380.0_p5.hbook         
n=[n]+1; file[n]=rho2012_380.0_p6.hbook         
n=[n]+1; file[n]=rho2012_380.0_p7.hbook         
n=[n]+1; file[n]=rho2012_380.0_p8.hbook         
n=[n]+1; file[n]=rho2012_380.0_p9.hbook         
n=[n]+1; file[n]=rho2012_382.0_p0.hbook         
n=[n]+1; file[n]=rho2012_382.0_p10.hbook        
n=[n]+1; file[n]=rho2012_382.0_p1.hbook         
n=[n]+1; file[n]=rho2012_382.0_p2.hbook         
n=[n]+1; file[n]=rho2012_382.0_p3.hbook         
n=[n]+1; file[n]=rho2012_382.0_p4.hbook         
n=[n]+1; file[n]=rho2012_382.0_p5.hbook         
n=[n]+1; file[n]=rho2012_382.0_p6.hbook         
n=[n]+1; file[n]=rho2012_382.0_p7.hbook         
n=[n]+1; file[n]=rho2012_382.0_p8.hbook         
n=[n]+1; file[n]=rho2012_382.0_p9.hbook         
n=[n]+1; file[n]=rho2012_384.0_p0.hbook         
n=[n]+1; file[n]=rho2012_384.0_p1.hbook         
n=[n]+1; file[n]=rho2012_384.0_p2.hbook         
n=[n]+1; file[n]=rho2012_386.0_p0.hbook         
n=[n]+1; file[n]=rho2012_386.0_p1.hbook         
n=[n]+1; file[n]=rho2012_387.5_p0.hbook         
n=[n]+1; file[n]=rho2012_387.5_p10.hbook        
n=[n]+1; file[n]=rho2012_387.5_p11.hbook        
n=[n]+1; file[n]=rho2012_387.5_p12.hbook        
n=[n]+1; file[n]=rho2012_387.5_p13.hbook        
n=[n]+1; file[n]=rho2012_387.5_p14.hbook        
n=[n]+1; file[n]=rho2012_387.5_p1.hbook         
n=[n]+1; file[n]=rho2012_387.5_p2.hbook         
n=[n]+1; file[n]=rho2012_387.5_p3.hbook         
n=[n]+1; file[n]=rho2012_387.5_p4.hbook         
n=[n]+1; file[n]=rho2012_387.5_p5.hbook         
n=[n]+1; file[n]=rho2012_387.5_p6.hbook         
n=[n]+1; file[n]=rho2012_387.5_p7.hbook         
n=[n]+1; file[n]=rho2012_387.5_p8.hbook         
n=[n]+1; file[n]=rho2012_387.5_p9.hbook         
n=[n]+1; file[n]=rho2012_388.5_p0.hbook         
n=[n]+1; file[n]=rho2012_388.5_p1.hbook         
n=[n]+1; file[n]=rho2012_388.5_p2.hbook         
n=[n]+1; file[n]=rho2012_389.0_p0.hbook         
n=[n]+1; file[n]=rho2012_389.0_p1.hbook         
n=[n]+1; file[n]=rho2012_389.5_p0.hbook         
n=[n]+1; file[n]=rho2012_389.5_p10.hbook        
n=[n]+1; file[n]=rho2012_389.5_p11.hbook        
n=[n]+1; file[n]=rho2012_389.5_p12.hbook        
n=[n]+1; file[n]=rho2012_389.5_p13.hbook        
n=[n]+1; file[n]=rho2012_389.5_p14.hbook        
n=[n]+1; file[n]=rho2012_389.5_p15.hbook        
n=[n]+1; file[n]=rho2012_389.5_p16.hbook        
n=[n]+1; file[n]=rho2012_389.5_p1.hbook         
n=[n]+1; file[n]=rho2012_389.5_p2.hbook         
n=[n]+1; file[n]=rho2012_389.5_p3.hbook         
n=[n]+1; file[n]=rho2012_389.5_p4.hbook         
n=[n]+1; file[n]=rho2012_389.5_p5.hbook         
n=[n]+1; file[n]=rho2012_389.5_p6.hbook         
n=[n]+1; file[n]=rho2012_389.5_p7.hbook         
n=[n]+1; file[n]=rho2012_389.5_p8.hbook         
n=[n]+1; file[n]=rho2012_389.5_p9.hbook         
n=[n]+1; file[n]=rho2012_390.0_p0.hbook         
n=[n]+1; file[n]=rho2012_390.0_p1.hbook         
n=[n]+1; file[n]=rho2012_390.5_p0.hbook         
n=[n]+1; file[n]=rho2012_390.5_p1.hbook         
n=[n]+1; file[n]=rho2012_390.5_p2.hbook         
n=[n]+1; file[n]=rho2012_391.5_p0.hbook         
n=[n]+1; file[n]=rho2012_391.5_p10.hbook        
n=[n]+1; file[n]=rho2012_391.5_p11.hbook        
n=[n]+1; file[n]=rho2012_391.5_p12.hbook        
n=[n]+1; file[n]=rho2012_391.5_p13.hbook        
n=[n]+1; file[n]=rho2012_391.5_p14.hbook        
n=[n]+1; file[n]=rho2012_391.5_p1.hbook         
n=[n]+1; file[n]=rho2012_391.5_p2.hbook         
n=[n]+1; file[n]=rho2012_391.5_p3.hbook         
n=[n]+1; file[n]=rho2012_391.5_p4.hbook         
n=[n]+1; file[n]=rho2012_391.5_p5.hbook         
n=[n]+1; file[n]=rho2012_391.5_p6.hbook         
n=[n]+1; file[n]=rho2012_391.5_p7.hbook         
n=[n]+1; file[n]=rho2012_391.5_p8.hbook         
n=[n]+1; file[n]=rho2012_391.5_p9.hbook         
n=[n]+1; file[n]=rho2012_393.0_p0.hbook         
n=[n]+1; file[n]=rho2012_393.0_p1.hbook         
n=[n]+1; file[n]=rho2012_395.0_p0.hbook         
n=[n]+1; file[n]=rho2012_395.0_p1.hbook         
n=[n]+1; file[n]=rho2012_395.0_p2.hbook         
n=[n]+1; file[n]=rho2012_397.0_p0.hbook         
n=[n]+1; file[n]=rho2012_397.0_p1.hbook         
n=[n]+1; file[n]=rho2012_397.0_p2.hbook         
n=[n]+1; file[n]=rho2012_400.0_p0.hbook         
n=[n]+1; file[n]=rho2012_400.0_p1.hbook         
n=[n]+1; file[n]=rho2012_400.0_p2.hbook         
n=[n]+1; file[n]=rho2012_400.0_p3.hbook         
n=[n]+1; file[n]=rho2012_409.0_p0.hbook         
n=[n]+1; file[n]=rho2012_409.0_p10.hbook        
n=[n]+1; file[n]=rho2012_409.0_p11.hbook        
n=[n]+1; file[n]=rho2012_409.0_p1.hbook         
n=[n]+1; file[n]=rho2012_409.0_p2.hbook         
n=[n]+1; file[n]=rho2012_409.0_p3.hbook         
n=[n]+1; file[n]=rho2012_409.0_p4.hbook         
n=[n]+1; file[n]=rho2012_409.0_p5.hbook         
n=[n]+1; file[n]=rho2012_409.0_p6.hbook         
n=[n]+1; file[n]=rho2012_409.0_p7.hbook         
n=[n]+1; file[n]=rho2012_409.0_p8.hbook         
n=[n]+1; file[n]=rho2012_409.0_p9.hbook         
n=[n]+1; file[n]=rho2012_420.0_p0.hbook         
n=[n]+1; file[n]=rho2012_420.0_p1.hbook         
n=[n]+1; file[n]=rho2012_420.0_p2.hbook         
n=[n]+1; file[n]=rho2012_420.0_p3.hbook         
n=[n]+1; file[n]=rho2012_420.0_p4.hbook         
n=[n]+1; file[n]=rho2012_420.0_p5.hbook         
n=[n]+1; file[n]=rho2012_420.0_p6.hbook         
n=[n]+1; file[n]=rho2012_420.0_p7.hbook         
n=[n]+1; file[n]=rho2012_430.0_p0.hbook         
n=[n]+1; file[n]=rho2012_430.0_p1.hbook         
n=[n]+1; file[n]=rho2012_430.0_p2.hbook         
n=[n]+1; file[n]=rho2012_430.0_p3.hbook         
n=[n]+1; file[n]=rho2012_430.0_p4.hbook         
n=[n]+1; file[n]=rho2012_440.0_p0.hbook         
n=[n]+1; file[n]=rho2012_440.0_p1.hbook         
n=[n]+1; file[n]=rho2012_440.0_p2.hbook         
n=[n]+1; file[n]=rho2012_450.0_p0.hbook         
n=[n]+1; file[n]=rho2012_450.0_p1.hbook         
n=[n]+1; file[n]=rho2012_450.0_p2.hbook         
n=[n]+1; file[n]=rho2012_460.0_p0.hbook         
n=[n]+1; file[n]=rho2012_460.0_p10.hbook        
n=[n]+1; file[n]=rho2012_460.0_p11.hbook        
n=[n]+1; file[n]=rho2012_460.0_p12.hbook        
n=[n]+1; file[n]=rho2012_460.0_p13.hbook        
n=[n]+1; file[n]=rho2012_460.0_p14.hbook        
n=[n]+1; file[n]=rho2012_460.0_p15.hbook        
n=[n]+1; file[n]=rho2012_460.0_p1.hbook         
n=[n]+1; file[n]=rho2012_460.0_p2.hbook         
n=[n]+1; file[n]=rho2012_460.0_p3.hbook         
n=[n]+1; file[n]=rho2012_460.0_p4.hbook         
n=[n]+1; file[n]=rho2012_460.0_p5.hbook         
n=[n]+1; file[n]=rho2012_460.0_p6.hbook         
n=[n]+1; file[n]=rho2012_460.0_p7.hbook         
n=[n]+1; file[n]=rho2012_460.0_p8.hbook         
n=[n]+1; file[n]=rho2012_460.0_p9.hbook         
n=[n]+1; file[n]=rho2012_470.0_p0.hbook         
n=[n]+1; file[n]=rho2012_470.0_p1.hbook         
n=[n]+1; file[n]=rho2012_470.0_p2.hbook         
n=[n]+1; file[n]=rho2012_470.0_p3.hbook         
n=[n]+1; file[n]=rho2012_470.0_p4.hbook         
n=[n]+1; file[n]=rho2012_480.0_p0.hbook         
n=[n]+1; file[n]=rho2012_480.0_p1.hbook         
n=[n]+1; file[n]=rho2012_480.0_p2.hbook         
n=[n]+1; file[n]=rho2012_480.0_p3.hbook         
n=[n]+1; file[n]=rho2012_480.0_p4.hbook         
n=[n]+1; file[n]=rho2012_480.0_p5.hbook         
n=[n]+1; file[n]=rho2012_490.0_p0.hbook         
n=[n]+1; file[n]=rho2012_490.0_p1.hbook         
n=[n]+1; file[n]=rho2012_490.0_p2.hbook         
n=[n]+1; file[n]=rho2012_490.0_p3.hbook         
n=[n]+1; file[n]=rho2012_490.0_p4.hbook         
n=[n]+1; file[n]=rho2012_490.0_p5.hbook         
dir=/work/users/konctbel/exp/RHO_2012-0/col/
f1=rho_2012-0_n1.05-1_profile.list
if ($fexist([f1]).eq.1) then
  shell rm [f1]
endif
for/file 20 [f1] n
close 20
fmess 0 [f1]
fmess [dir] [f1]
fmess [dir] [f1]
f2=rho_2012-0_n1.13_profile.list
if ($fexist([f2]).eq.1) then
  shell rm [f2]
endif
for/file 20 [f2] n
close 20
fmess 0 [f2]
fmess [dir] [f2]
fmess [dir] [f2]
f3=rho_2012-0_n1.05-2_profile.list
if ($fexist([f3]).eq.1) then
  shell rm [f3]
endif
for/file 20 [f3] n
close 20
fmess 0 [f3]
fmess [dir] [f3]
fmess [dir] [f3]
run1=16699
run2=17869
1d 100 ! 10000 10000.5 20000.5
ve/cre irun0(10000) r
ve/cre ind0(10000) r
sigma ind0 = array(10000,10001#20000)
do i=1,[n]
  hi/file 20 [dir][file[i]]
  nt/pl 1.run idh=100
  hi/get/cont 100 irun0
  ve/del ind,irun
  sigma ind=order(ind0,-irun0)
  sigma irun=order(irun0,-irun0)
  exec mapcal#vecut irun
  m=$vlen(irun)
  exec mapcal#vecut ind [m]
  rmin=$sigma(vmin(ind))
  rmax=$sigma(vmax(ind))
  mean=$hinfo(100,'mean')
  rms=$hinfo(100,'rms')
  l=$sigma([mean]-3*max([rms],3))
  r=$sigma([mean]+3*max([rms],3))
  if ($sigma(int([l])).eq.[l]) then
    l=[l]-0.5
  endif
  if ($sigma(int([r])).eq.[r]) then
    r=[r]+0.5
  endif
  hi/pl 100([l]:[r])
  line [rmin] $GRAFINFO('WNYMIN') [rmin] $GRAFINFO('WNYMAX')
  line [rmax] $GRAFINFO('WNYMIN') [rmax] $GRAFINFO('WNYMAX')
  if (([run1].gt.[rmin]).and.([run1].le.[rmax])) then
    mess wery strange!!! -> [file[i]]: [run1] in {[rmin],[rmax]}
  endif
  if (([run2].gt.[rmin]).and.([run2].le.[rmax])) then
    mess wery strange!!! -> [file[i]]: [run2] in {[rmin],[rmax]}
  endif
  if ([run1].gt.[rmax]) then
    fmess [file[i]] [f1]
    mess [rmin] [rmax] [f1]
  endif
  if (([run1].le.[rmin]).and.([run2].gt.[rmax])) then
    fmess [file[i]] [f2]
    mess [rmin] [rmax] [f2]
  endif
  if ([run2].le.[rmin]) then
    fmess [file[i]] [f3]
    mess [rmin] [rmax] [f3]
  endif
  close 20
*  read x
enddo
return

macro wlspos pt=1
ve/cre  fw(9) r
ve/cre dfw(9) r 9*0.3
do i=1,9
  ve/del p16,dp16
  ve/read p16,dp16 setup[pt]_exp_v0_amplitude_vs_phi_counter[i]_slx150_fit_100_et0_chit.par '2f15.10'
  ve/inp fw([i]) p16(13)
  ve/inp dfw([i]) $sigma(dp16(13)*sqrt(p16(40)))
enddo
ve/cre nw(9) r 1 2 3 4 5 6 7 8 9
ve/cre dnw(9) r 
* exec $PER/s#vpl fw dfw nw dnw 
ve/cre p2(2) r -25 40
ve/cre dp2(2) r
ve/cre s2(2) r 1 0
ve/fit nw fw dfw p1 sb 2 p2 s2 ! ! dp2
a=p2(1)
b=p2(2)
sigma lfw = fw-(([a])+([b])*nw)
* exec $PER/s#vpl lfw dfw nw dnw ll=-1
mess $sigma([a]+25) ($sigma(dp2(1)))
return


macro phicutsprep
n=0
*n=[n]+1; pt[n]=s; nrf[n]=1;       mod[n]=sim; pref[n]=setup[pt[n]]
*n=[n]+1; pt[n]=0; nrf[n]=5781;  mod[n]=exp
n=[n]+1; pt[n]=11; nrf[n]=7842;   mod[n]=sim; pref[n]=mapcal1
n=[n]+1; pt[n]=12; nrf[n]=11285;  mod[n]=sim; pref[n]=mapcal2
n=[n]+1; pt[n]=13; nrf[n]=13847;  mod[n]=sim; pref[n]=mapcal3
n=[n]+1; pt[n]=14; nrf[n]=16699;  mod[n]=sim; pref[n]=mapcal4
*
n=[n]+1; pt[n]=11; nrf[n]=7842;   mod[n]=exp; pref[n]=mapcal1
n=[n]+1; pt[n]=12; nrf[n]=11285;  mod[n]=exp; pref[n]=mapcal2
n=[n]+1; pt[n]=13; nrf[n]=13847;  mod[n]=exp; pref[n]=mapcal3
n=[n]+1; pt[n]=14; nrf[n]=16699;  mod[n]=exp; pref[n]=mapcal4
*n=[n]+1; pt[n]=4; nrf[n]=17869;  mod[n]=exp
*
di=[n]/2
ve/cre nrx([di]) r
ve/cre ix([di]) r 
sigma ix=array([di],1#[di])
do i=1,[n]
  if ([mod[i]].eq.'sim') then
    shell rm tmp[0-9].txt
    shell ls [pref[i]]/sim_v0/*.hbook > tmp0.txt
    shell $unquote('cat tmp0.txt | sed "s/\(^.*kb\)/ /g" > tmp1.txt')
    shell $unquote('cat tmp1.txt | sed "s/-/ /g" > tmp2.txt')
    shell $unquote('cat tmp2.txt | sed "s/[^0-9]/ /g" > tmp3.txt')
    ve/del vb,vf,vr,vn
    ve/read vb,vf,vr,vn tmp3.txt
    mrf[i]=$sigma(vmin(vr))
    ie=[i]+[di]
    mrf[ie]=$sigma(vmin(vr))
    ve/inp nrx([i]) $sigma(vmin(vr))
  endif
enddo
sigma ix=order(ix,nrx)
*
fname=setup_all_phi_cuts.txt
if ($fexist([fname]).eq.1) then
  shell rm [fname]
endif
for/file 20 [fname] n
close 20
*
fnamem=setup_sim_all_phi_cuts.txt
if ($fexist([fnamem]).eq.1) then
  shell rm [fnamem]
endif
for/file 20 [fnamem] n
close 20
*
do j=1,[n]
  ve/cre  fw(9) r
  ve/cre dfw(9) r 9*0.3
*  fname[j]=setup[pt[j]]_phi_cuts.txt
  fname[j]=[pref[j]]_[mod[j]]_phi_cuts.txt
  if ($fexist([fname[j]]).eq.1) then
    shell rm [fname[j]]
  endif
  for/file 20 [fname[j]] n
  close 20
  do i=1,9
    ve/del p16,dp16
*    ve/read p16,dp16 setup[pt[j]]_exp_v0_amplitude_vs_phi_counter[i]_slx150_fit_100_et0_chit.par '2f15.10'
    ve/read p16,dp16 [pref[j]]/[pref[j]]_[mod[j]]_v0_amplitude_vs_phi_counter[i]_slx150_fit_100_et0_chit.par '2f15.10'
    ve/inp fw([i]) p16(13)
    ve/inp dfw([i]) $sigma(dp16(13)*sqrt(p16(40)))
    txt=$sigma(p16(12)) $sigma(p16(13)) $sigma(p16(14)) $sigma(p16(15))
    fmess [txt] [fname[j]]
  enddo
  txt=[mod[j]] [nrf[j]] [fname[j]]
  fmess [txt] [fname]
  txt[j]=[mod[j]] [mrf[j]] [fname[j]]
enddo
*
fname=setup_sim_all_phi_cuts.txt
if ($fexist([fname]).eq.1) then
  shell rm [fname]
endif
for/file 20 [fname] n
close 20
*
do i=1,[n]
  ii=[i]
  if ([i].gt.[di]) then
    ii=[ii]-[di]
  endif
  j=ix([ii])
  if ([i].gt.[di]) then
    j=[j]+[di]
  endif
  fmess [txt[j]] [fname]
enddo
return

macro plxxxx i1=1 i2=9 prep=0 dir=v2
if ([dir].eq.'v3') then
  corr=0
else
  corr=1
endif
do i=[i1],[i2]
  fname=plxxx[i].kumac
  if ($fexist([fname]).eq.1) then
    shell rm [fname]
  endif
  for/file 20 [fname] new
  close 20
  txt=exec mapcal#plxxx [i] ! ! [prep] days [corr] [dir]
  fmess [txt] [fname]
*  shell pawbigX11 -b [fname] >& log[i].log &
  shell pawbigX11 -b [fname]
*  exec mapcal#plxxx [i] 1 100000 [prep] days [corr] [dir]
enddo
return


macro plxxx nc=1 run1=1 run2=100000 prep=0 msum=days corr=0 dir=v2
*
if ([prep].eq.1) then
*
  dirc=/work/users/konctbel/AGPARS
  ve/del runc,ac1,ac2,ac3,ac4,ac5,ac6,ac7,ac8,ac9
  ve/read runc,ac1,ac2,ac3,ac4,ac5,ac6,ac7,ac8,ac9 [dirc]/AmpCorrHV/amplitude_correction_hv.txt '(10f15.6)'  
  rb=runc(1)
*  
  ve/del runs,days,amp[nc],damp[nc],amp[nc]c,damp[nc]c,n0[nc]b,nx[nc]b,rms
*  
  if ([nc].ne.0) then
    ncx1=[nc]
    ncx2=[nc]
  else
    ncx1=1
    ncx2=9
  endif
*  
  do ncx=[ncx1],[ncx2]
    ve/del runs,days,amp[ncx]i,damp[ncx]i,n0[ncx]bi,nx[ncx]bi,eff[ncx]i,deff[ncx]i
    do i=7,17
      ve/read runs[i],days[i],amp[ncx]p[i],damp[ncx]p[i] [dir]/amplitude_vs_runs_counter[ncx]_p[i].txt 4f15.6
      ve/read runs[i],days[i],n0[ncx]b[i],nx[ncx]b[i] [dir]/n0_nx_vs_runs_counter[ncx]_p[i].txt 4f15.6
      exec vappend runs runs[i]
      exec vappend days days[i]
      exec vappend  amp[ncx]i  amp[ncx]p[i]
      exec vappend damp[ncx]i damp[ncx]p[i]
      exec vappend n0[ncx]bi n0[ncx]b[i]
      exec vappend nx[ncx]bi nx[ncx]b[i]
      ve/del runs[i],days[i],amp[ncx]p[i],damp[ncx]p[i],n0[ncx]b[i],nx[ncx]b[i]
    enddo  
    n=$vlen(runs)
    ve/cre  ac[ncx]x([n]) r
    do i=1,[n]
      ind=$sigma(runs([i])-[rb]+1)
      ve/inp ac[ncx]x([i]) ac[ncx]([ind])
    enddo
    if ([ncx].eq.[ncx1]) then
      sigma amp[nc] = amp[ncx]i/damp[ncx]i**2
      sigma rms = 1.0/damp[ncx]i**2
      sigma amp[nc]c = amp[ncx]i/damp[ncx]i**2/ac[ncx]x
      sigma rmsc = 1.0/damp[ncx]i**2/ac[ncx]x**2
      sigma n0[nc]b = n0[ncx]bi
      sigma nx[nc]b = nx[ncx]bi
    else
      sigma amp[nc] = amp[nc] + amp[ncx]i/damp[ncx]i**2
      sigma rms = rms + 1.0/damp[ncx]i**2
      sigma amp[nc]c = amp[nc]c + amp[ncx]i/damp[ncx]i**2/ac[ncx]x
      sigma rmsc = rmsc + 1.0/damp[ncx]i**2/ac[ncx]x**2
      sigma n0[nc]b = n0[nc]b + n0[ncx]bi
      sigma nx[nc]b = nx[nc]b + nx[ncx]bi
    endif
  enddo
*  
  sigma amp[nc] = amp[nc]/rms
  sigma damp[nc] = 1.0/sqrt(rms)
  sigma amp[nc]c = amp[nc]c/rmsc
  sigma damp[nc]c = 1.0/sqrt(rmsc)
  ve/del rms,rmsc
*
  sigma druns = runs*0
  rmin=$sigma(max([run1],vmin(runs)))
  rmax=$sigma(min([run2],vmax(runs)))
  null [rmin] [rmax] 0 12
  set pmci 2
  * exec $PER/s#vpl amp[nc] damp[nc] runs druns sz=0.01 o=s
  sigma eff[nc]b = n0[nc]b/nx[nc]b
  sigma deff[nc]b = sqrt(eff[nc]b*(1-eff[nc]b)/nx[nc]b)
  sigma eff[nc]p = -log(eff[nc]b)
  sigma deff[nc]p = deff[nc]b/eff[nc]b
  set pmci 4
  * exec $PER/s#vpl eff[nc]p deff[nc]p runs druns sz=0.01 o=s
*  read x
  ve/write runs,days,amp[nc],damp[nc],n0[nc]b,nx[nc]b [dir]/amp_event_vs_runs_counter[nc].txt '6f15.6'
*
  ve/cre dcuts(1) r -22.31
*
  exec mapcal#daysum amp[nc] damp[nc] runs days dcuts
  exec mapcal#daysumn nx[nc]b runs days dcuts
  exec mapcal#daysumn n0[nc]b runs days dcuts
  ve/write runsz,drunsz,daysz,ddaysz,amp[nc]z,damp[nc]z,n0[nc]bz,nx[nc]bz [dir]/amp_event_vs_days_counter[nc].txt '8f15.6'
*
  exec mapcal#daysum amp[nc]c damp[nc]c runs days dcuts
  ve/write runsz,drunsz,daysz,ddaysz,amp[nc]cz,damp[nc]cz,n0[nc]bz,nx[nc]bz [dir]/ampc_event_vs_days_counter[nc].txt '8f15.6'
*
  exec mapcal#calsum amp[nc] damp[nc] runs days
  exec mapcal#calsumn nx[nc]b runs days
  exec mapcal#calsumn n0[nc]b runs days
  ve/write runsz,drunsz,daysz,ddaysz,amp[nc]z,damp[nc]z,n0[nc]bz,nx[nc]bz [dir]/amp_event_vs_cals_counter[nc].txt '8f15.6'
*
  exec mapcal#calsum amp[nc]c damp[nc]c runs days
  ve/write runsz,drunsz,daysz,ddaysz,amp[nc]cz,damp[nc]cz,n0[nc]bz,nx[nc]bz [dir]/ampc_event_vs_cals_counter[nc].txt '8f15.6'
*
endif
*
if ([dir].eq.'v3') then
  goto 1
endif
*
if ([corr].eq.0) then
  ve/del runsz,drunsz,daysz,ddaysz,amp[nc]z,damp[nc]z,n0[nc]bz,nx[nc]bz
  ve/read runsz,drunsz,daysz,ddaysz,amp[nc]z,damp[nc]z,n0[nc]bz,nx[nc]bz [dir]/amp_event_vs_[msum]_counter[nc].txt '8f15.6'
else
  ve/del runsz,drunsz,daysz,ddaysz,amp[nc]z,damp[nc]z,n0[nc]bz,nx[nc]bz
  ve/read runsz,drunsz,daysz,ddaysz,amp[nc]z,damp[nc]z,n0[nc]bz,nx[nc]bz [dir]/ampc_event_vs_[msum]_counter[nc].txt '8f15.6'
endif

set pmci 2
* exec $PER/s#vpl amp[nc]z damp[nc]z daysz ddaysz sz=0.01
sigma eff[nc]bz = n0[nc]bz/nx[nc]bz
sigma deff[nc]bz = sqrt(eff[nc]bz*(1-eff[nc]bz)/nx[nc]bz)
sigma amp[nc]pz = -log(eff[nc]bz) 
sigma damp[nc]pz = deff[nc]bz/eff[nc]bz

*sigma smu[nc]z = sqrt(2*(amp[nc]z-amp[nc]pz))/amp[nc]z
*sigma sigmu1 = ((amp[nc]z-2*amp[nc]pz)*damp[nc]z/amp[nc]z)**2
*sigma sigmu2 = (damp[nc]pz/amp[nc]pz)**2
*sigma dsmu[nc]z = sqrt(sigmu1+sigmu2)/smu[nc]z/amp[nc]z**2

sigma smu[nc]z = 2*(amp[nc]z-amp[nc]pz)/amp[nc]z**2
sigma sigmu1 = ((2*amp[nc]pz/amp[nc]z-1)*damp[nc]z)**2
sigma sigmu2 = (damp[nc]pz)**2
sigma dsmu[nc]z = 2*sqrt(sigmu1+sigmu2)/amp[nc]z**2
  
set pmci 4
* exec $PER/s#vpl amp[nc]pz damp[nc]pz daysz ddaysz sz=0.01 o=s
*read x
*
*ve/cre runsx(16) r 7842 10988 11285 13219 13515 13724 13725 13846 13847 14047 14048 16698 16699 17868 17869 21075
ve/cre runsx(8) r 7842 10988 11285 13846 13847 16698 16699 17868
exec mapcal#dcuts runsx daysx
ve/cre daysx1(2) r -30 180
ve/cre daysx2(2) r 410 710
ve/cre daysx3(4) r 695 750 760 790
ve/cre daysx4(2) r 786 835
*
ve/del act,act0,nact,nact0,dst,aft,aft0
n=$vlen(daysz)
*
do i=1,4
*
  j=[i]*2-1
  tx=daysx([j])
  exec mapcal#ixndl [tx] daysz
  gl/imp inx
  in1=[inx]
  j=[i]*2
  tx=daysx([j])
  exec mapcal#ixndr [tx] daysz
  gl/imp inx
  in2=[inx]
*  
  ve/del ampt,dampt,ampq,dampq,dayst,ddayst,runst,drunst,smut,dsmut
  ve/copy  smu[nc]z([in1]:[in2])  smut
  ve/copy dsmu[nc]z([in1]:[in2]) dsmut
  ve/copy  amp[nc]z([in1]:[in2])  ampt
  ve/copy damp[nc]z([in1]:[in2]) dampt
  ve/copy  amp[nc]pz([in1]:[in2])  ampq
  ve/copy damp[nc]pz([in1]:[in2]) dampq
  ve/copy  daysz([in1]:[in2])  dayst
  ve/copy ddaysz([in1]:[in2]) ddayst
  ve/copy  runsz([in1]:[in2])  runst
  ve/copy drunsz([in1]:[in2]) drunst
  
*  sigma smut = 2*(ampt-ampq)/ampt**2
*  sigma sigmu1 = ((2*ampq/ampt-1)*dampt)**2
*  sigma sigmu2 = (dampq)**2
*  sigma dsmut = 2*sqrt(sigmu1+sigmu2)/ampt**2
  
  sigma smut = sqrt(2*(ampt-ampq))/ampt
  sigma sigmu1 = ((ampt-2*ampq)*dampt/ampt)**2
  sigma sigmu2 = dampq**2
  sigma dsmut = sqrt(sigmu1+sigmu2)/smut/ampt**2
  
  set pmci 1
  set lwid 1
  set hwid 1
  * exec $PER/s#vpl ampt dampt dayst ddayst sz=0.1 iatt=20
  * exec $PER/s#vpl ampq dampq dayst ddayst sz=0.1 iatt=24 o=s
*  
  nz=$vlen(daysx[i])
  nz=[nz]/2
  ve/del ampt0,dampt0,ampq0,dampq0,dayst0,ddayst0,runst0,drunst0,smut0,dsmut0
  do k=1,[nz]
  
    j=[k]*2-1
    tx=daysx[i]([j])
    exec mapcal#ixndl [tx] dayst
    gl/imp inx
    in1=[inx]
    j=[k]*2
    tx=daysx[i]([j])
    exec mapcal#ixndr [tx] dayst
    gl/imp inx
    in2=[inx]
    
    ve/del amptt,damptt,ampqt,dampqt,daystt,ddaystt,runstt,drunstt,smutt,dsmutt
    ve/copy  smut([in1]:[in2])  smutt
    ve/copy dsmut([in1]:[in2]) dsmutt
    ve/copy  ampt([in1]:[in2])  amptt
    ve/copy dampt([in1]:[in2]) damptt
    ve/copy  ampq([in1]:[in2])  ampqt
    ve/copy dampq([in1]:[in2]) dampqt
    ve/copy  dayst([in1]:[in2])  daystt
    ve/copy ddayst([in1]:[in2]) ddaystt
    ve/copy  runst([in1]:[in2])  runstt
    ve/copy drunst([in1]:[in2]) drunstt
    
    exec vappend  smut0  smutt
    exec vappend dsmut0 dsmutt
    exec vappend  ampt0  amptt
    exec vappend dampt0 damptt
    exec vappend  ampq0  ampqt
    exec vappend dampq0 dampqt
    exec vappend dayst0 daystt
    exec vappend ddayst0 ddaystt
    exec vappend runst0 runstt
    exec vappend drunst0 drunstt

  enddo
*------------  
  nz=$vlen(dayst0)
  ve/cre is([nz]) r
  sigma is = array([nz],1#[nz])
  nz0=[nz]
  do l=1,[nz]
    ampl=ampq0([l])
    if (([ampl].lt.2).or.([ampl].gt.15)) then
      ve/inp is([l]) $sigma(is([l])+[nz])
      nz0=[nz0]-1
    endif
  enddo
  ve/del amptf,damptf,ampqf,dampqf,daystf,ddaystf,runstf,drunstf,smutf,dsmutf
  sigma  smutf = order( smut0,is)
  sigma dsmutf = order(dsmut0,is)
  sigma  amptf = order( ampt0,is)
  sigma damptf = order(dampt0,is)
  sigma  ampqf = order( ampq0,is)
  sigma dampqf = order(dampq0,is)
  sigma  daystf = order( dayst0,is)
  sigma ddaystf = order(ddayst0,is)
  sigma  runstf = order( runst0,is)
  sigma drunstf = order(drunst0,is)
  exec mapcal#vecut  smutf  [nz0]
  exec mapcal#vecut dsmutf  [nz0]
  exec mapcal#vecut  amptf  [nz0]
  exec mapcal#vecut damptf  [nz0]
  exec mapcal#vecut  ampqf  [nz0]
  exec mapcal#vecut dampqf  [nz0]
  exec mapcal#vecut  daystf [nz0]
  exec mapcal#vecut ddaystf [nz0]
  exec mapcal#vecut  runstf [nz0]
  exec mapcal#vecut drunstf [nz0]
  mess nz: [nz] [nz0]
* 
  ve/del dampqfx
  ve/copy dampqf dampqfx
  nz=$vlen(daystf)
  ve/cre is([nz]) r
  sigma is = array([nz],1#[nz])
  nz0=[nz]
  do ll=1,1
    mean=$sigma(vsum(ampqf)/[nz])
    mean2=$sigma(vsum(ampqf**2)/[nz])
    rms=$sigma(sqrt([mean2]-[mean]**2))
    mess mean=[mean] rms=[rms]
    do l=1,[nz]
      if ($sigma(abs(ampqf([l])-[mean])/[rms]).gt.3) then
        ve/inp dampqfx([l]) 10
        ve/inp is([l]) $sigma(is([l])+[nz])
        mess ll=[ll] l=[l] 
        nz0=[nz0]-1
      endif
    enddo
  enddo
  sigma  smutf = order( smutf,is)
  sigma dsmutf = order(dsmutf,is)
  sigma  amptf = order( amptf,is)
  sigma damptf = order(damptf,is)
  sigma  ampqf = order( ampqf,is)
  sigma dampqf = order(dampqf,is)
  sigma  daystf = order( daystf,is)
  sigma ddaystf = order(ddaystf,is)
  sigma  runstf = order( runstf,is)
  sigma drunstf = order(drunstf,is)
  exec mapcal#vecut  smutf  [nz0]
  exec mapcal#vecut dsmutf  [nz0]
  exec mapcal#vecut  amptf  [nz0]
  exec mapcal#vecut damptf  [nz0]
  exec mapcal#vecut  ampqf  [nz0]
  exec mapcal#vecut dampqf  [nz0]
  exec mapcal#vecut  daystf [nz0]
  exec mapcal#vecut ddaystf [nz0]
  exec mapcal#vecut  runstf [nz0]
  exec mapcal#vecut drunstf [nz0]
  mess nz: [nz] [nz0]
*
  if ([i].eq.1) then
*  
    npar=7
    ve/cre chi2(2) r
    ve/cre paru([npar]) r
    ve/cre dparu([npar]) r
    ve/cre covu([npar],[npar]) r
*    
    w=$sigma(vmax(daystf)-vmin(daystf))
*    
    fname=[dir]/mapcal[i]_fit_ampt_vs_[msum]_counter0.txt
    if ($fexist([fname]).eq.1) then
      ve/read p7,dp7 [fname] '2f15.6'
    else
      ve/cre p7(7) r $sigma([mean]-0.5) +1 $sigma(vsum(daystf)/[nz0]-[w]/3) $sigma([w]/8) _
                                        -1 $sigma(vsum(daystf)/[nz0]+[w]/3) $sigma([w]/4)
      ve/cre dp7(7) r                  
    endif
*    
    if ([nc].eq.0) then
      ve/cre s7(7) r $sigma([mean]/100) 0.01 $sigma(vsum(daystf)/100) $sigma([w]/100) _
                                        0.01 $sigma(vsum(daystf)/100) $sigma([w]/100)
      ve/cre pmin(7) r $sigma([mean]/2) $sigma(-[mean]/2) $sigma(vmin(daystf)) $sigma(30) _
                                        $sigma(-[mean]/2) $sigma(vmin(daystf)) $sigma(30)
      ve/cre pmax(7) r $sigma([mean]*2) $sigma(+[mean]/2) $sigma(vmax(daystf)) $sigma([w]) _
                                        $sigma(+[mean]/2) $sigma(vmax(daystf)) $sigma([w])
    else
      ve/cre s7(7) r $sigma([mean]/100) 0.01 0 0 0.01 0 0
      ve/del pmin,pmax
      ve/cre st7(7) r 2 1 0 0 1 0 0
      sigma pmin = p7-abs(dp7)-st7
      sigma pmax = p7+abs(dp7)+st7
    endif
*    
    ve/cre xf(1) r
    do nf=1,3
    
      do l=1,3
        ve/fit daystf ampqf dampqf derf0.f sb 7 p7 s7 pmin pmax dp7
      enddo
    
      nz=$vlen(dayst0)
      ve/cre is([nz]) r
      sigma is = array([nz],1#[nz])
      nz0=[nz]
      do l=1,[nz]
        ve/inp xf(1) $sigma(dayst0([l]))
        derf0pl=$call('derf0p.f(xf)')
        if ($sigma(abs(ampq0([l])-[derf0pl])/dampq0([l])).gt.3) then
          ve/inp is([l]) $sigma(is([l])+[nz])
          nz0=[nz0]-1
        endif
      enddo
      ve/del amptf,damptf,ampqf,dampqf,daystf,ddaystf,runstf,drunstf,smutf,dsmutf
      sigma  smutf = order( smut0,is)
      sigma dsmutf = order(dsmut0,is)
      sigma  amptf = order( ampt0,is)
      sigma damptf = order(dampt0,is)
      sigma  ampqf = order( ampq0,is)
      sigma dampqf = order(dampq0,is)
      sigma  daystf = order( dayst0,is)
      sigma ddaystf = order(ddayst0,is)
      sigma  runstf = order( runst0,is)
      sigma drunstf = order(drunst0,is)
      exec mapcal#vecut  smutf  [nz0]
      exec mapcal#vecut dsmutf  [nz0]
      exec mapcal#vecut  amptf  [nz0]
      exec mapcal#vecut damptf  [nz0]
      exec mapcal#vecut  ampqf  [nz0]
      exec mapcal#vecut dampqf  [nz0]
      exec mapcal#vecut  daystf [nz0]
      exec mapcal#vecut ddaystf [nz0]
      exec mapcal#vecut  runstf [nz0]
      exec mapcal#vecut drunstf [nz0]
      
    enddo
*    
    set pmci 1
    set lwid 1
    set hwid 1
    * exec $PER/s#vpl ampt dampt dayst ddayst sz=0.1 iatt=20
    * exec $PER/s#vpl ampq dampq dayst ddayst sz=0.1 iatt=24 o=s
    if ([nc].eq.0) then
      opt=s
    else
      opt=sb
    endif
    do l=1,3
      ve/fit daystf ampqf dampqf derf0.f sb 7 p7 s7 pmin pmax dp7
    enddo
*    
    fname=[dir]/mapcal[i]_fit_ampt_vs_[msum]_counter[nc].txt
    ve/write p7,dp7 [fname] '2f15.6'
    ve/write p7,dp7 ! '2f15.6'
*    
  endif
*  
  if ([i].ne.1) then
*  
    npar=4
    ve/cre chi2(2) r
    ve/cre paru([npar]) r
    ve/cre dparu([npar]) r
    ve/cre covu([npar],[npar]) r
*    
    w=$sigma(vmax(daystf)-vmin(daystf))
*    
    fname=[dir]/mapcal[i]_fit_ampt_vs_[msum]_counter0.txt
    if ($fexist([fname]).eq.1) then
      ve/read p4,dp4 [fname] '2f15.6'
    else
      ve/cre p4(4) r $sigma([mean]+0.5) -1 $sigma(vsum(daystf)/[nz0]) $sigma((vmax(daystf)-vmin(daystf))/4)
      ve/cre dp4(4) r
    endif
*    
    if ([nc].eq.0) then
      ve/cre s4(4) r $sigma([mean]/100) 0.01 $sigma(vsum(daystf)/100) $sigma([w]/100)
      ve/cre pmin(4) r $sigma([mean]/2) $sigma(-[mean]/2) $sigma(vmin(daystf)) $sigma(5)
      ve/cre pmax(4) r $sigma([mean]*2) $sigma(+[mean]/2) $sigma(vmax(daystf)) $sigma(vmax(daystf)-vmin(daystf))
    else
      ve/cre s4(4) r $sigma([mean]/100) 0.01 0 0
      ve/del pmin,pmax
      ve/cre st4(4) r 2 1 $sigma(abs(p4(4))) $sigma(abs(p4(4))/2)
      sigma pmin = p4-st4
      sigma pmax = p4+st4
    endif
*
    ve/cre xf(1) r
    do nf=1,3
    
      do l=1,3
        ve/fit daystf ampqf dampqf erf0.f sb 4 p4 s4 pmin pmax dp4
      enddo
    
      nz=$vlen(dayst0)
      ve/cre is([nz]) r
      sigma is = array([nz],1#[nz])
      nz0=[nz]
      do l=1,[nz]
        ve/inp xf(1) $sigma(dayst0([l]))
        erf0pl=$call('erf0p.f(xf)')
        if ($sigma(abs(ampq0([l])-[erf0pl])/dampq0([l])).gt.5) then
          ve/inp is([l]) $sigma(is([l])+[nz])
          nz0=[nz0]-1
        endif
      enddo
      ve/del amptf,damptf,ampqf,dampqf,daystf,ddaystf,runstf,drunstf,smutf,dsmutf
      sigma  smutf = order( smut0,is)
      sigma dsmutf = order(dsmut0,is)
      sigma  amptf = order( ampt0,is)
      sigma damptf = order(dampt0,is)
      sigma  ampqf = order( ampq0,is)
      sigma dampqf = order(dampq0,is)
      sigma  daystf = order( dayst0,is)
      sigma ddaystf = order(ddayst0,is)
      sigma  runstf = order( runst0,is)
      sigma drunstf = order(drunst0,is)
      exec mapcal#vecut  smutf  [nz0]
      exec mapcal#vecut dsmutf  [nz0]
      exec mapcal#vecut  amptf  [nz0]
      exec mapcal#vecut damptf  [nz0]
      exec mapcal#vecut  ampqf  [nz0]
      exec mapcal#vecut dampqf  [nz0]
      exec mapcal#vecut  daystf [nz0]
      exec mapcal#vecut ddaystf [nz0]
      exec mapcal#vecut  runstf [nz0]
      exec mapcal#vecut drunstf [nz0]

    enddo
*    
    set pmci 1
    set lwid 1
    set hwid 1
    * exec $PER/s#vpl ampt dampt dayst ddayst sz=0.1 iatt=20
    * exec $PER/s#vpl ampq dampq dayst ddayst sz=0.1 iatt=24 o=s
    if ([nc].eq.0) then
      opt=s
    else
      opt=sb
    endif
    do l=1,3
      ve/fit daystf ampqf dampqf erf0.f sb 4 p4 s4 pmin pmax dp4
    enddo
    
    fname=[dir]/mapcal[i]_fit_ampt_vs_[msum]_counter[nc].txt
    ve/write p4,dp4 [fname] '2f15.6'
        
  endif
*  
  set pmci 2
  * exec $PER/s#vpl ampqf dampqf daystf ddaystf sz=0.1 iatt=24 o=s
*------------------------  
  nz=$vlen(dayst0)
  ve/cre is([nz]) r
  sigma is = array([nz],1#[nz])
  nz0=[nz]
  do l=1,[nz]
    ampl=ampt0([l])
    if (([ampl].lt.2).or.([ampl].gt.15)) then
      ve/inp is([l]) $sigma(is([l])+[nz])
      nz0=[nz0]-1
    endif
  enddo
  ve/del ampth,dampth,ampqh,dampqh,daysth,ddaysth,runsth,drunsth,smuth,dsmuth
  sigma  smuth = order( smut0,is)
  sigma dsmuth = order(dsmut0,is)
  sigma  ampth = order( ampt0,is)
  sigma dampth = order(dampt0,is)
  sigma  ampqh = order( ampq0,is)
  sigma dampqh = order(dampq0,is)
  sigma  daysth = order( dayst0,is)
  sigma ddaysth = order(ddayst0,is)
  sigma  runsth = order( runst0,is)
  sigma drunsth = order(drunst0,is)
  exec mapcal#vecut  smuth  [nz0]
  exec mapcal#vecut dsmuth  [nz0]
  exec mapcal#vecut  ampth  [nz0]
  exec mapcal#vecut dampth  [nz0]
  exec mapcal#vecut  ampqh  [nz0]
  exec mapcal#vecut dampqh  [nz0]
  exec mapcal#vecut  daysth [nz0]
  exec mapcal#vecut ddaysth [nz0]
  exec mapcal#vecut  runsth [nz0]
  exec mapcal#vecut drunsth [nz0]
  mess nz: [nz] [nz0]
* 
  ve/del dampqhx
  ve/copy dampqh dampqhx
  nz=$vlen(daysth)
  ve/cre is([nz]) r
  sigma is = array([nz],1#[nz])
  nz0=[nz]
  do ll=1,1
    mean=$sigma(vsum(ampqh)/[nz])
    mean2=$sigma(vsum(ampqh**2)/[nz])
    rms=$sigma(sqrt([mean2]-[mean]**2))
    mess mean=[mean] rms=[rms]
    do l=1,[nz]
      if ($sigma(abs(ampqh([l])-[mean])/[rms]).gt.3) then
        ve/inp dampqhx([l]) 10
        ve/inp is([l]) $sigma(is([l])+[nz])
        mess ll=[ll] l=[l] 
        nz0=[nz0]-1
      endif
    enddo
  enddo
  sigma  smuth = order( smuth,is)
  sigma dsmuth = order(dsmuth,is)
  sigma  ampth = order( ampth,is)
  sigma dampth = order(dampth,is)
  sigma  ampqh = order( ampqh,is)
  sigma dampqh = order(dampqh,is)
  sigma  daysth = order( daysth,is)
  sigma ddaysth = order(ddaysth,is)
  sigma  runsth = order( runsth,is)
  sigma drunsth = order(drunsth,is)
  exec mapcal#vecut  smuth  [nz0]
  exec mapcal#vecut dsmuth  [nz0]
  exec mapcal#vecut  ampth  [nz0]
  exec mapcal#vecut dampth  [nz0]
  exec mapcal#vecut  ampqh  [nz0]
  exec mapcal#vecut dampqh  [nz0]
  exec mapcal#vecut  daysth [nz0]
  exec mapcal#vecut ddaysth [nz0]
  exec mapcal#vecut  runsth [nz0]
  exec mapcal#vecut drunsth [nz0]
  mess nz: [nz] [nz0]
*
  ve/cre t0(1) r $sigma(vmin(dayst0))
*
  if ([i].eq.1) then
*  
    npar=2
    ve/cre chi2(2) r
    ve/cre paru([npar]) r
    ve/cre dparu([npar]) r
    ve/cre covu([npar],[npar]) r
*    
    w=$sigma(vmax(daystf)-vmin(daystf))
*    
    fname=[dir]/mapcal[i]_fit_ampt_vs_[msum]_counter[nc].txt
    if ($fexist([fname]).eq.1) then
      ve/read p7,dp7 [fname] '2f15.6'
    endif
    fname=[dir]/mapcal[i]_fit_amp_vs_[msum]_counter0.txt
    if ($fexist([fname]).eq.1) then
      ve/read p3,dp3 [fname] '2f15.6'
    endif
*    
    ve/cre pmin(3) r 0.1 0 5
    ve/cre pmax(3) r 0.5 0.5 100
    if ([nc].eq.0) then
      tau=30
      dtau=1
    else
      tau=p3(3)
      dtau=0
    endif
    ve/cre p3(3) r 0.25 0.1 [tau]
    ve/cre dp3(3) r 
    ve/cre s3(3) r 0.001 0.001 [dtau]
*    
    ve/cre xf(1) r
    do nf=1,3
    
      do l=1,3
        ve/fit daysth ampth dampth derf0u.f sb 3 p3 s3 pmin pmax dp3
      enddo
    
      nz=$vlen(dayst0)
      ve/cre is([nz]) r
      sigma is = array([nz],1#[nz])
      nz0=[nz]
      do l=1,[nz]
        ve/inp xf(1) $sigma(dayst0([l]))
        derf0pl=$call('derf0up.f(xf)')
        if ($sigma(abs(ampt0([l])-[derf0pl])/dampt0([l])).gt.10) then
          ve/inp is([l]) $sigma(is([l])+[nz])
          nz0=[nz0]-1
        endif
      enddo
      ve/del ampth,dampth,ampqh,dampqh,daysth,ddaysth,runsth,drunsth,smuth,dsmuth
      sigma  smuth = order( smut0,is)
      sigma dsmuth = order(dsmut0,is)
      sigma  ampth = order( ampt0,is)
      sigma dampth = order(dampt0,is)
      sigma  ampqh = order( ampq0,is)
      sigma dampqh = order(dampq0,is)
      sigma  daysth = order( dayst0,is)
      sigma ddaysth = order(ddayst0,is)
      sigma  runsth = order( runst0,is)
      sigma drunsth = order(drunst0,is)
      exec mapcal#vecut  smuth  [nz0]
      exec mapcal#vecut dsmuth  [nz0]
      exec mapcal#vecut  ampth  [nz0]
      exec mapcal#vecut dampth  [nz0]
      exec mapcal#vecut  ampqh  [nz0]
      exec mapcal#vecut dampqh  [nz0]
      exec mapcal#vecut  daysth [nz0]
      exec mapcal#vecut ddaysth [nz0]
      exec mapcal#vecut  runsth [nz0]
      exec mapcal#vecut drunsth [nz0]
      
    enddo
*    
    set pmci 1
    set lwid 1
    set hwid 1
    * exec $PER/s#vpl ampt dampt dayst ddayst sz=0.1 iatt=20
    * exec $PER/s#vpl ampq dampq dayst ddayst sz=0.1 iatt=24 o=s

    do l=1,3
      ve/fit daysth ampth dampth derf0u.f sb 3 p3 s3 pmin pmax dp3
    enddo
*    
    fname=[dir]/mapcal[i]_fit_amp_vs_[msum]_counter[nc].txt
    ve/write p3,dp3 [fname] '2f15.6'
    ve/write p3,dp3 ! '2f15.6'
*    
  endif
*  
  if ([i].ne.1) then
*  
    npar=4
    ve/cre chi2(2) r
    ve/cre paru([npar]) r
    ve/cre dparu([npar]) r
    ve/cre covu([npar],[npar]) r
*    
    w=$sigma(vmax(daystf)-vmin(daystf))
*    
    fname=[dir]/mapcal[i]_fit_ampt_vs_[msum]_counter[nc].txt
    if ($fexist([fname]).eq.1) then
      ve/read p4,dp4 [fname] '2f15.6'
    endif
    fname=[dir]/mapcal[i]_fit_amp_vs_[msum]_counter0.txt
    if ($fexist([fname]).eq.1) then
      ve/read p4a,dp4a [fname] '2f15.6'
    endif
*    
    ve/cre pmin(4) r 0.1 0 10 0.1
    ve/cre pmax(4) r $sigma(1/sqrt(2*p4(1))) 0.5 100 2
    if ([nc].eq.0) then
      tau=30
      dtau=1
    else
      tau=p4a(3)
      dtau=0
    endif
    if (([i].eq.2).or.([i].eq.4)) then
      ve/cre p4a(4) r 0.25 0 10000 1
      ve/cre dp4a(4) r 
      ve/cre s4a(4) r 0.001 0 0 0
      ve/cre tr(1) r 100000
    endif
    if ([i].eq.3) then
      ve/cre p4a(4) r 0.25 0.1 [tau] 1
      ve/cre dp4a(4) r 
      ve/cre s4a(4) r 0.001 0.001 [dtau] 0.01
      ve/cre tr(1) r 751
      ks=0.01
      ve/cre s4a(4) r 0.001 0 0 [ks]
      ve/fit daysth ampth dampth erf0u.f sb 4 p4a s4a pmin pmax dp4a
      ve/cre s4a(4) r 0.001 0.001 0 [ks]
      ve/fit daysth ampth dampth erf0u.f sb 4 p4a s4a pmin pmax dp4a
      ve/cre s4a(4) r 0 0.001 0 0
      ve/fit daysth ampth dampth erf0u.f sb 4 p4a s4a pmin pmax dp4a
      ve/cre s4a(4) r 0.001 0.001 [dtau] [ks]
      ve/fit daysth ampth dampth erf0u.f sb 4 p4a s4a pmin pmax dp4a
    endif
*
    ve/cre xf(1) r
    do nf=1,3
    
      do l=1,3
        ve/fit daysth ampth dampth erf0u.f sb 4 p4a s4a pmin pmax dp4a
      enddo
    
      nz=$vlen(dayst0)
      ve/cre is([nz]) r
      sigma is = array([nz],1#[nz])
      nz0=[nz]
      do l=1,[nz]
        ve/inp xf(1) $sigma(dayst0([l]))
        erf0pl=$call('erf0up.f(xf)')
        if ($sigma(abs(ampt0([l])-[erf0pl])/dampt0([l])).gt.10) then
          ve/inp is([l]) $sigma(is([l])+[nz])
          nz0=[nz0]-1
        endif
      enddo
      ve/del ampth,dampth,ampqh,dampqh,daysth,ddaysth,runsth,drunsth,smuth,dsmuth
      sigma  smuth = order( smut0,is)
      sigma dsmuth = order(dsmut0,is)
      sigma  ampth = order( ampt0,is)
      sigma dampth = order(dampt0,is)
      sigma  ampqh = order( ampq0,is)
      sigma dampqh = order(dampq0,is)
      sigma  daysth = order( dayst0,is)
      sigma ddaysth = order(ddayst0,is)
      sigma  runsth = order( runst0,is)
      sigma drunsth = order(drunst0,is)
      exec mapcal#vecut  smuth  [nz0]
      exec mapcal#vecut dsmuth  [nz0]
      exec mapcal#vecut  ampth  [nz0]
      exec mapcal#vecut dampth  [nz0]
      exec mapcal#vecut  ampqh  [nz0]
      exec mapcal#vecut dampqh  [nz0]
      exec mapcal#vecut  daysth [nz0]
      exec mapcal#vecut ddaysth [nz0]
      exec mapcal#vecut  runsth [nz0]
      exec mapcal#vecut drunsth [nz0]

    enddo
*    
    fname=[dir]/mapcal[i]_fit_amp_vs_[msum]_counter[nc].txt
    ve/write p4a,dp4a [fname] '2f15.6'
    
    set pmci 1
    set lwid 1
    set hwid 1
    * exec $PER/s#vpl ampt dampt dayst ddayst sz=0.1 iatt=20
    * exec $PER/s#vpl ampq dampq dayst ddayst sz=0.1 iatt=24 o=s
    if ([nc].eq.0) then
      opt=s
    else
      opt=sb
    endif
    do l=1,3
      ve/fit daysth ampth dampth erf0u.f sb 4 p4a s4a pmin pmax dp4a
    enddo
    
  endif
*  
  set pmci 2
  * exec $PER/s#vpl ampth dampth daysth ddaysth sz=0.1 iatt=24 o=s
*--------------  
  if ([i].eq.1) then
  
    fname=[dir]/mapcal[i]_fit_ampt_vs_[msum]_counter[nc].txt
    ve/read p7,dp7 [fname] '2f15.6'
    npi=$vlen(dayst)
*    
    fname=[dir]/mapcal[i]_fit_amp_vs_[msum]_counter[nc].txt
    ve/read p3,dp3 [fname] '2f15.6'
    ve/inp p3(2) 0
    ve/cre  aci0([npi]) r 
    ve/cre naci0([npi]) r 
    ve/cre afti0([npi]) r 
    ve/cre xpi(1) r
    do ip=1,[npi]
      ve/inp xpi(1) $sigma(dayst([ip]))
      ai=$call('derf0up.f(xpi)')
      ve/inp  aci0([ip]) $sigma([ai]/ampt([ip]))
      ve/inp naci0([ip]) $sigma(abs([ai]-ampt([ip]))/dampt([ip]))
      ve/inp afti0([ip]) [ai]
    enddo
    exec vappend  act0  aci0
    exec vappend nact0 naci0
    exec vappend  aft0 afti0
*    
    ve/read p3,dp3 [fname] '2f15.6'
    fun/pl derf0up.f $GRAFINFO('WNXMIN') $GRAFINFO('WNXMAX') s
    ve/cre  aci([npi]) r 
    ve/cre naci([npi]) r 
    ve/cre afti([npi]) r
    do ip=1,[npi]
      ve/inp xpi(1) $sigma(dayst([ip]))
      ai=$call('derf0up.f(xpi)')
      ve/inp  aci([ip]) $sigma([ai]/ampt([ip]))
      ve/inp naci([ip]) $sigma(abs([ai]-ampt([ip]))/dampt([ip]))
      ve/inp afti([ip]) [ai]
    enddo
    exec vappend  act  aci
    exec vappend nact naci
    exec vappend  aft afti
    exec vappend dct dayst
    
  endif
*  
  if ([i].ge.2) then
  
    fname=[dir]/mapcal[i]_fit_ampt_vs_[msum]_counter[nc].txt
    ve/read p4,dp4 [fname] '2f15.6'
    npi=$vlen(dayst)
*    
    fname=[dir]/mapcal[i]_fit_amp_vs_[msum]_counter[nc].txt
    ve/read p4a,dp4a [fname] '2f15.6'
    ve/inp p4a(2) 0
    ve/inp p4a(4) 1
    ve/cre  aci0([npi]) r 
    ve/cre naci0([npi]) r 
    ve/cre afti0([npi]) r
    ve/cre xpi(1) r
    do ip=1,[npi]
      ve/inp xpi(1) $sigma(dayst([ip]))
      ai=$call('erf0up.f(xpi)')
      ve/inp  aci0([ip]) $sigma([ai]/ampt([ip]))
      ve/inp naci0([ip]) $sigma(abs([ai]-ampt([ip]))/dampt([ip]))
      ve/inp afti0([ip]) [ai]
    enddo
    exec vappend  act0  aci0
    exec vappend nact0 naci0
    exec vappend  aft0 afti0
*    
    ve/read p4a,dp4a [fname] '2f15.6'
    ve/inp p4a(4) 1
    fun/pl erf0up.f $GRAFINFO('WNXMIN') $GRAFINFO('WNXMAX') s
    ve/cre  aci([npi]) r
    ve/cre naci([npi]) r
    ve/cre afti([npi]) r
    do ip=1,[npi]
      ve/inp xpi(1) $sigma(dayst([ip]))
      ai=$call('erf0up.f(xpi)')
      ve/inp  aci([ip]) $sigma([ai]/ampt([ip]))
      ve/inp naci([ip]) $sigma(abs([ai]-ampt([ip]))/dampt([ip]))
      ve/inp afti([ip]) [ai]
    enddo
    exec vappend  act  aci
    exec vappend nact naci
    exec vappend  aft afti
    exec vappend dct dayst
    ve/read p4a,dp4a [fname] '2f15.6'
    
    if ([i].eq.3) then
    
      tx=daysx[i](2)
      exec mapcal#ixndl [tx] dayst
      gl/imp inx
      in1=[inx]
      tx=daysx[i](3)
      exec mapcal#ixndr [tx] dayst
      gl/imp inx
      in2=[inx]
      
      do j=[in1],[in2]
        ve/inp aci0([j]) 1
        ve/inp aci([j]) 1
        ve/inp naci0([j]) 0
        ve/inp naci([j]) 0
      enddo
      
    endif
    
  endif
*
  ve/write runst,drunst,dayst,ddayst,aci0,naci0,afti0 [dir]/mapcal[i]_ampcorr0_[msum]_counter[nc].txt '(7f15.6)'
  ve/write runst,drunst,dayst,ddayst,aci,naci,afti [dir]/mapcal[i]_ampcorr_[msum]_counter[nc].txt '(7f15.6)'
*  
*    
  set pmci 1
  call covm.f(1)
  call covmpen(chi2,[npar],paru,dparu)
  call covmcov([npar],covu)
  ndf=$sigma(chi2(2))
  chi=$sigma(int(chi2(1)*[ndf]*100+0.5)/100)
  txt=[h]^2!/ndf = [chi] / [ndf]
  exec $PER/s#tf 0.1 0.9 [txt]
*  
  l0=$sigma(vmin(dayst))
  r0=$sigma(vmax(dayst))
  l=$sigma([l0]-0.05*([r0]-[l0]))
  r=$sigma([r0]+0.05*([r0]-[l0]))
  if ([i].eq.3) then
    ve/cre vum(10) r 6 6 6 6 6 6 6 6 6 6
  else
    ve/cre vum(10) r 10 10 13 13 10 10 10 10 10 10 
  endif
  inc=[nc]+1
  um=vum([inc])
  null [l] [r] 0 [um]
  set pmci 1
  set lwid 1
  set hwid 1
  * exec $PER/s#vpl ampt dampt dayst ddayst sz=0.1 iatt=20 o=s
  * exec $PER/s#vpl ampq dampq dayst ddayst sz=0.1 iatt=24 o=s
  set pmci 2
  * exec $PER/s#vpl ampqf dampqf daystf ddaystf sz=0.1 iatt=24 o=s
  set pmci 6
  * exec $PER/s#vpl ampth dampth daysth ddaysth sz=0.1 iatt=24 o=s
  if ([i].eq.1) then
    fun/pl derf0p.f [l0] [r0] s
    fun/pl derf0up.f [l0] [r0] s
    ve/del p3t
    ve/copy p3 p3t
    ve/inp p3(2) 0
    set plci 2
    set ltyp 12
    set basl 0.005
    fun/pl derf0up.f [l0] [r0] s
    set plci 1
    set ltyp 1
    ve/copy p3t p3
  endif
  if ([i].ne.1) then
    fun/pl erf0p.f [l0] [r0] s
    fun/pl erf0up.f [l0] [r0] s
    ve/del p4at
    ve/copy p4a p4at
    ve/inp p4a(2) 0
    ve/inp p4a(4) 1
    set plci 2
    set ltyp 14
    set basl 0.005
    fun/pl erf0up.f [l0] [r0] s
    set plci 1
    set ltyp 1
    ve/copy p4at p4a
  endif
  ve/del amptc,damptc
  sigma amptc = ampt*aci
  sigma damptc = 0*aci
  set pmci 4
  * exec $PER/s#vpl amptc damptc dayst ddayst sz=0.03 iatt=20 o=s
  atitle 't, days' '[m], pe.' 
  exec save [dir]/mapcal[i]_amp_vs_[msum]_counter[nc].eps f
*  read x
*  
  null [l] [r] 0 [um]
  set pmci 1
  set lwid 1
  set hwid 1
  * exec $PER/s#vpl ampt dampt dayst ddayst sz=0.1 iatt=20 o=s
  * exec $PER/s#vpl ampq dampq dayst ddayst sz=0.1 iatt=24 o=s
  atitle 't, days' '[m], pe.' 
  exec save [dir]/mapcal[i]_amp_vs_[msum]_0_counter[nc].eps f
*
  null [l] [r] 0 0.4
  set pmci 1
  set lwid 1
  set hwid 1
  * exec $PER/s#vpl smut0 dsmut0 dayst ddayst sz=0.1 iatt=20 o=s
  if ([i].eq.1) then
    v=p3
  else
    v=p4a
  endif
  umax=$sigma(abs([v](1)))
  umin=$sigma(abs([v](2)))
  tau=$sigma(abs([v](3)))
  t0i=t0(1)
  fun/pl [umax]-[umin]*exp(-(x-([t0i]))/[tau]) [l] [r] s
  atitle 't, days' '[s]?[m]!/[m]'
  exec save [dir]/mapcal[i]_sigmu_vs_[msum]_counter[nc].eps f
*  read x
*
  ve/del amptc,damptc,smutc,sigmu1,sigmu2,dsmutc
  sigma amptc = ampt*aci
  sigma damptc = dampt*aci
  sigma smutc = sqrt(2*(amptc-ampq))/amptc
  sigma sigmu1 = ((amptc-2*ampq)*damptc/amptc)**2
  sigma sigmu2 = dampq**2
  sigma dsmutc = sqrt(sigmu1+sigmu2)/smutc/amptc**2
  null [l] [r] 0 0.4
  set pmci 1
  set lwid 1
  set hwid 1
  * exec $PER/s#vpl smutc dsmutc dayst ddayst sz=0.1 iatt=20 o=s
  if ([i].eq.1) then
    v=p3
  else
    v=p4a
  endif
  umax=$sigma(abs([v](1)))
  umin=$sigma(abs([v](2)))
  tau=$sigma(abs([v](3)))
  t0i=t0(1)
  fun/pl [umax]-[umin]*exp(-(x-([t0i]))/[tau]) [l] [r] s
  atitle 't, days' '[s]?[m]!/[m]'
  exec save [dir]/mapcal[i]_sigmuc_vs_[msum]_counter[nc].eps f
*  read x
enddo
*ve/write dct,act [dir]/mapcal_ampcorr_[msum]_counter[nc].txt '(2f15.6)'
shell cat [dir]/mapcal?_ampcorr0_[msum]_counter[nc].txt > [dir]/mapcal_ampcorr0_[msum]_counter[nc].txt
shell cat [dir]/mapcal?_ampcorr_[msum]_counter[nc].txt > [dir]/mapcal_ampcorr_[msum]_counter[nc].txt
*
1:
return


macro ampcorrprx
do i=0,7
  shell mv ampcorrpr[i].log ampcorrpr[i].log.old
  shell mv ampcorr_full_new_p[i].txt ampcorr_full_new_p[i].txt.old
  shell mv ampcorr_nsigs_full_new_p[i].txt ampcorr_nsigs_full_new_p[i].txt.old
  fname=ampcorrpr[i].kumac
  if ($fexist([fname]).eq.1) then
    shell rm [fname]
  endif
  for/file 20 [fname] new
  close 20
  txt=exec mapcal#ampcorrpr [i]
  fmess [txt] [fname]
  shell pawbigX11 -b [fname]
enddo
*
shell cat ampcorrpr*.log > ampcorrpr.log
shell cat ampcorr_full_new_p*.txt > ampcorr_full_new.txt
shell cat ampcorr_nsigs_full_new_p*.txt > ampcorr_nsigs_full_new.txt
*
dir=v2
do nc=1,9
  ve/del runsi,daysi,amp[nc]i,damp[nc]i,n0[nc]bi,nx[nc]bi
  do i=7,17
    ve/del runs[i],days[i],amp[nc]p[i],damp[nc]p[i]
    ve/read runs[i],days[i],amp[nc]p[i],damp[nc]p[i] [dir]/amplitude_vs_runs_counter[nc]_p[i].txt 4f15.6
    exec vappend runsi runs[i]
    exec vappend daysi days[i]
    exec vappend amp[nc]i amp[nc]p[i]
    exec vappend damp[nc]i damp[nc]p[i]
    ve/del runs[i],days[i],n0[nc]b[i],nx[nc]b[i]
    ve/read runs[i],days[i],n0[nc]b[i],nx[nc]b[i] [dir]/n0_nx_vs_runs_counter[nc]_p[i].txt 4f15.6
    exec vappend n0[nc]bi n0[nc]b[i]
    exec vappend nx[nc]bi nx[nc]b[i]
  enddo
  ve/write runsi,daysi,amp[nc]i,damp[nc]i [dir]/amplitude_vs_runs_counter[nc].txt 4f15.6
  ve/write runsi,daysi,n0[nc]bi,nx[nc]bi  [dir]/n0_nx_vs_runs_counter[nc].txt 4f15.6
enddo
return



macro ampcorrpr ipart=0
*
dir=v2
do nc=1,9
  ve/del runs,days,amp[nc],damp[nc],n0[nc]b,nx[nc]b
  ve/read runs,days,amp[nc],damp[nc],n0[nc]b,nx[nc]b [dir]/amp_event_vs_runs_counter[nc].txt '6f15.6'
enddo
nruns=$vlen(runs)
*
dirc=/work/users/konctbel/AGPARS
ve/del runv,av1,av2,av3,av4,av5,av6,av7,av8,av9
ve/read runv,av1,av2,av3,av4,av5,av6,av7,av8,av9 [dirc]/AmpCorrHV/amplitude_correction_hv.txt '(10f15.6)'  
*
do nc=1,9
  ve/del runst,drunst,dayst,ddayst,aci0[nc],naci0[nc],afti0[nc]
  ve/read runst,drunst,dayst,ddayst,aci0[nc],naci0[nc],afti0[nc] mapcal_ampcorr0_days_counter[nc].txt '(7f15.6)'
enddo
sigma rb = runst-drunst
sigma re = runst+drunst
*
rmin=$sigma(vmin(rb))
rmax=$sigma(vmax(re))
*
ind1=$sigma([ipart]*1000+1)
ind2=$sigma(min(([ipart]+1)*1000,[nruns]))
if ([ind1].gt.[ind2]) then
  goto 1
endif
rmin=$sigma(runs([ind1]))
rmax=$sigma(runs([ind2]))
*
n=$sigma([rmax]-[rmin]+1)
ve/cre runc([n]) r
sigma runc=array([n],[rmin]#[rmax])
*
do nc=1,9
  ve/cre  ac[nc]([n]) r [n]*1
  ve/cre nac[nc]([n]) r [n]*0
enddo
*
fname=ampcorrpr[ipart].log
if ($fexist([fname])) then
  shell rm [fname]
endif
for/file 20 [fname] new
close 20
*
id=0
rei=0
rvmin=$sigma(vmin(runv)-1)
do i=[ind1],[ind2]
  ri=runs([i])
  while ([rei].lt.[ri]) do
    id=[id]+1
    if ([id].le.$vlen(re)) then
      rei=re([id])
      rbi=rb([id])
    else
      rei=100000
      rbi=100000
    endif
  endwhile
  ind=$sigma([ri]-[rvmin])
  do nc=1,9
    chv=av[nc]([ind])
    if ([rbi].le.[ri]) then
      cft=aci0[nc]([id])
      sft=$sigma(abs(amp[nc]([i])-afti0[nc]([id]))/damp[nc]([i]))
      txt=[i] [rbi] [ri] [rei] $sigma(runv([ind]))
    else
      cft=1
      sft=0
      txt=[i] [rbi] [ri] [rei] $sigma(runv([ind])) *
    endif
    if ([nc].eq.1) then
      fmess [txt] [fname]
    endif
    j=[ri]-[rmin]+1
    ve/inp  ac[nc]([j]) $sigma([chv]*[cft])
    ve/inp nac[nc]([j]) [sft]
  enddo
enddo
*
ve/write runc,ac1,ac2,ac3,ac4,ac5,ac6,ac7,ac8,ac9 ampcorr_full_new_p[ipart].txt '(10f15.6)'
ve/write runc,nac1,nac2,nac3,nac4,nac5,nac6,nac7,nac8,nac9 ampcorr_nsigs_full_new_p[ipart].txt '(10f15.6)'
*
1:
return




macro amp1pecorr
*
ve/del runc,ac1,ac2,ac3,ac4,ac5,ac6,ac7,ac8,ac9
ve/read runc,ac1,ac2,ac3,ac4,ac5,ac6,ac7,ac8,ac9 ampcorr_full_new.txt '(10f15.6)'
ve/del clb,run,amp1,amp2,amp3,amp4,amp5,amp6,amp7,amp8,amp9
ve/read clb,run,amp1,amp2,amp3,amp4,amp5,amp6,amp7,amp8,amp9 accled_amp1pe_work.txt '(2f10.1,9f15.6)'
*
n=$vlen(run)
do nc=1,9
  ve/cre kc[nc]([n]) r [n]*1
enddo
*
do i=1,[n]
*
  rb=run([i])
  if ([i].lt.[n]) then
    j=[i]+1
    re=run([j])
  else
    re=$sigma(vmax(runc))
  endif
*  
  exec mapcal#ixndl [rb] runc
  gl/imp inx
  in1=[inx]
  exec mapcal#ixndr [re] runc
  gl/imp inx
  in2=[inx]
*
  if ([in1].lt.[in2]) then
    nk=[in2]-[in1]+1
    ve/del act
    ve/cre act([nk]) r
    do nc=1,9
      ve/copy ac[nc]([in1]:[in2]) act    
      ve/inp kc[nc]([i]) $sigma(vsum(act)/[nk])
    enddo
  endif
enddo
*
ve/cre dn([n]) r
*
n=0
n=[n]+1; nrf[n]=7842; 
n=[n]+1; nrf[n]=11285;
n=[n]+1; nrf[n]=13847;
n=[n]+1; nrf[n]=16699;
n=[n]+1; nrf[n]=17869;
*
do nc=1,9
  ve/cre ae[nc]([n]) r
  ve/cre zv[nc]([n]) r
enddo
ve/cre runf([n]) r 
*
do i=1,[n]
  ve/inp runf([i]) [nrf[i]]
*
  rb=[nrf[i]]
  if ([i].lt.[n]) then
    j=[i]+1
    re=[nrf[j]]
  else
    re=$sigma(vmax(run))
  endif
*  
  exec mapcal#ixndl [rb] run
  gl/imp inx
  in1=[inx]
  exec mapcal#ixndr [re] run
  gl/imp inx
  in2=[inx]
*
  if ([in1].lt.[in2]) then
    nk=[in2]-[in1]+1
    nz=$sigma(int([nk]/10+0.5))
    no=[nk]-2*[nz]
    ve/del act
    ve/cre act([nk]) r
    ve/cre ix([nk]) r
    do nc=1,9
      ve/copy amp[nc]([in1]:[in2]) act    
      sigma ix=array([nk],1#[nk])
      sigma ix=order(ix,act)
      ve/cre fx([nk]) r [nz]*0 [no]*1 [nz]*0
      sigma fx=order(fx,ix)
      ve/inp ae[nc]([i]) $sigma(vsum(act*fx)/[no])
    enddo
  endif
enddo
*
ve/cre a1s(9) r 85.553 73.093 96.013 72.988 94.227 92.113 82.458 136.194 102.672
ve/cre a1s(9) r 90.043 77.241 105.399 76.985 97.676 114.848 94.784 94.174  98.95
*
do nc=1,9
  set pmci 2
  * exec $PER/s#vpl amp[nc] dn run dn sz=0.1
  sigma amp[nc]c = amp1/kc[nc]
  set pmci 4
  * exec $PER/s#vpl amp[nc]c dn run dn sz=0.05 o=s
*
  do i=1,[n]
    line [nrf[i]] $GRAFINFO('WNYMIN') [nrf[i]] $GRAFINFO('WNYMAX')
*    
    rb=[nrf[i]]
    if ([i].lt.[n]) then
      j=[i]+1
      re=[nrf[j]]
    else
      re=$sigma(vmax(run))
    endif
    l=ae[nc]([i])
    line [rb] [l] [re] [l]
  enddo
  l=a1s([nc])
  line $GRAFINFO('WNXMIN') [l] $GRAFINFO('WNXMAX') [l]
  read x
enddo
*
ve/cre zv([n]) r
*
n=0
n=[n]+1; v[n]=zv; f[n]=u 
n=[n]+1; v[n]=zv; f[n]=a 
n=[n]+1; v[n]=zv; f[n]=b 
n=[n]+1; v[n]=zv; f[n]=c 
n=[n]+1; v[n]=ae; f[n]=amp1pe
n=[n]+1; v[n]=zv; f[n]=eff1pe
n=[n]+1; v[n]=zv; f[n]=u0
n=[n]+1; v[n]=zv; f[n]=s0
n=[n]+1; v[n]=zv; f[n]=pds
n=[n]+1; v[n]=zv; f[n]=tmin
n=[n]+1; v[n]=zv; f[n]=tmax
*
do v=1,[n]
  fname=accled_[f[v]]_mc.txt
  if ($fexist([fname])) then
    shell rm [fname]
  endif
  for/file 20 [fname] ! n
  close 20
  ve/write zv,runf,[v[v]]1,[v[v]]2,[v[v]]3,[v[v]]4,[v[v]]5,[v[v]]6,[v[v]]7,[v[v]]8,[v[v]]9 [fname] '(2f10.1,9f15.6)'
enddo
return



macro acccorrprep
*
ve/del runc,ac1,ac2,ac3,ac4,ac5,ac6,ac7,ac8,ac9
ve/read runc,ac1,ac2,ac3,ac4,ac5,ac6,ac7,ac8,ac9 ampcorr_full_new.txt '(10f15.6)'
*
n=0
n=[n]+1; nrf[n]=7842; 
n=[n]+1; nrf[n]=11285;
n=[n]+1; nrf[n]=13847;
n=[n]+1; nrf[n]=16699;
*n=[n]+1; nrf[n]=17869;
*
ve/del acr
do nc=1,9
  ve/del acr[nc]
enddo
*
do i=1,[n]
*
  rb=[nrf[i]]
  if ([i].eq.[n]) then
    re=$sigma(vmax(runc))
  else
    j=[i]+1
    re=[nrf[j]]-1
  endif
*  
  exec mapcal#ixndl [rb] runc
  gl/imp inx
  in1=[inx]
  exec mapcal#ixndr [re] runc
  gl/imp inx
  in2=[inx]
  
  mess $sigma(runc([in1])) $sigma(runc([in2])) [rb] [re]
*
  ve/read acr AmplitudeCorrection_mapcal[i].txt
*
  nx=[in2]-[in1]+1
  do nc=1,9
    ve/cre acrx([nx]) r [nx]*$sigma(acr([nc]))
    exec vappend acr[nc] acrx
  enddo
enddo
*
do nc=1,9
  sigma ac[nc]=ac[nc]*acr[nc]
  sigma ac[nc]=ac[nc]+(1-ac[nc]/ac[nc])
enddo
*
ve/write runc,ac1,ac2,ac3,ac4,ac5,ac6,ac7,ac8,ac9 acccorr_vs_run_to_db.txt '(10f15.6)'  
return



macro plxxxold nc=1 run1=1 run2=100000 prep=0 msum=days corr=0
dir=v2
if ([prep].eq.1) then
dirc=/work/users/konctbel/AGPARS
ve/del runc,ac1,ac2,ac3,ac4,ac5,ac6,ac7,ac8,ac9
ve/read runc,ac1,ac2,ac3,ac4,ac5,ac6,ac7,ac8,ac9 [dirc]/AmpCorrHV/amplitude_correction_hv.txt '(10f15.6)'
rb=runc(1)
ve/del runs,days,amp[nc],damp[nc],amp[nc]c,damp[nc]c,n0[nc]b,nx[nc]b,rms
if ([nc].ne.0) then
  ncx1=[nc]
  ncx2=[nc]
else
  ncx1=1
  ncx2=9
endif
do ncx=[ncx1],[ncx2]
  ve/del runs,days,amp[ncx]i,damp[ncx]i,n0[ncx]bi,nx[ncx]bi,eff[ncx]i,deff[ncx]i
  do i=7,17
    ve/read runs[i],days[i],amp[ncx]p[i],damp[ncx]p[i] [dir]/amplitude_vs_runs_counter[ncx]_p[i].txt 4f15.6
    ve/read runs[i],days[i],n0[ncx]b[i],nx[ncx]b[i] [dir]/n0_nx_vs_runs_counter[ncx]_p[i].txt 4f15.6
    exec vappend runs runs[i]
    exec vappend days days[i]
    exec vappend  amp[ncx]i  amp[ncx]p[i]
    exec vappend damp[ncx]i damp[ncx]p[i]
    exec vappend n0[ncx]bi n0[ncx]b[i]
    exec vappend nx[ncx]bi nx[ncx]b[i]
    ve/del runs[i],days[i],amp[ncx]p[i],damp[ncx]p[i],n0[ncx]b[i],nx[ncx]b[i]
  enddo  
  n=$vlen(runs)
  ve/cre  ac[ncx]x([n]) r
  do i=1,[n]
    ind=$sigma(runs([i])-[rb]+1)
    ve/inp ac[ncx]x([i]) ac[ncx]([ind])
  enddo
  if ([ncx].eq.[ncx1]) then
    sigma amp[nc] = amp[ncx]i/damp[ncx]i**2
    sigma rms = 1.0/damp[ncx]i**2
    sigma amp[nc]c = amp[ncx]i/damp[ncx]i**2/ac[ncx]x
    sigma rmsc = 1.0/damp[ncx]i**2/ac[ncx]x**2
    sigma n0[nc]b = n0[ncx]bi
    sigma nx[nc]b = nx[ncx]bi
  else
    sigma amp[nc] = amp[nc] + amp[ncx]i/damp[ncx]i**2
    sigma rms = rms + 1.0/damp[ncx]i**2
    sigma amp[nc]c = amp[nc]c + amp[ncx]i/damp[ncx]i**2/ac[ncx]x
    sigma rmsc = rmsc + 1.0/damp[ncx]i**2/ac[ncx]x**2
    sigma n0[nc]b = n0[nc]b + n0[ncx]bi
    sigma nx[nc]b = nx[nc]b + nx[ncx]bi
  endif
  sigma 
enddo
sigma amp[nc] = amp[nc]/rms
sigma damp[nc] = 1.0/sqrt(rms)
sigma amp[nc]c = amp[nc]c/rmsc
sigma damp[nc]c = 1.0/sqrt(rmsc)
ve/del rms,rmsc
*
sigma druns = runs*0
rmin=$sigma(max([run1],vmin(runs)))
rmax=$sigma(min([run2],vmax(runs)))
null [rmin] [rmax] 0 10
set pmci 2
* exec $PER/s#vpl amp[nc] damp[nc] runs druns sz=0.01 o=s
sigma eff[nc]b = n0[nc]b/nx[nc]b
sigma deff[nc]b = sqrt(eff[nc]b*(1-eff[nc]b)/nx[nc]b)
sigma eff[nc]p = -log(eff[nc]b)
sigma deff[nc]p = deff[nc]b/eff[nc]b
set pmci 4
* exec $PER/s#vpl eff[nc]p deff[nc]p runs druns sz=0.01 o=s
*read x
ve/write runs,days,amp[nc],damp[nc],n0[nc]b,nx[nc]b [dir]/amp_event_vs_runs_counter[nc].txt '6f15.6'
*
exec mapcal#daysum amp[nc] damp[nc] days
exec mapcal#daysumn nx[nc]b days
exec mapcal#daysumn n0[nc]b days
ve/write daysz,ddaysz,amp[nc]z,damp[nc]z,n0[nc]bz,nx[nc]bz [dir]/amp_event_vs_days_counter[nc].txt '6f15.6'
*
exec mapcal#daysum amp[nc]c damp[nc]c days
ve/write daysz,ddaysz,amp[nc]cz,damp[nc]cz,n0[nc]bz,nx[nc]bz [dir]/ampc_event_vs_days_counter[nc].txt '6f15.6'
*
exec mapcal#calsum amp[nc] damp[nc] runs days
exec mapcal#calsumn nx[nc]b runs days
exec mapcal#calsumn n0[nc]b runs days
ve/write daysz,ddaysz,amp[nc]z,damp[nc]z,n0[nc]bz,nx[nc]bz [dir]/amp_event_vs_cals_counter[nc].txt '6f15.6'
*
exec mapcal#calsum amp[nc]c damp[nc]c runs days
ve/write daysz,ddaysz,amp[nc]cz,damp[nc]cz,n0[nc]bz,nx[nc]bz [dir]/ampc_event_vs_cals_counter[nc].txt '6f15.6'
*
endif
*
if ([corr].eq.0) then
  ve/del daysz,ddaysz,amp[nc]z,damp[nc]z,n0[nc]bz,nx[nc]bz
  ve/read daysz,ddaysz,amp[nc]z,damp[nc]z,n0[nc]bz,nx[nc]bz [dir]/amp_event_vs_[msum]_counter[nc].txt '6f15.6'
else
  ve/del daysz,ddaysz,amp[nc]z,damp[nc]z,n0[nc]bz,nx[nc]bz
  ve/read daysz,ddaysz,amp[nc]z,damp[nc]z,n0[nc]bz,nx[nc]bz [dir]/ampc_event_vs_[msum]_counter[nc].txt '6f15.6'
endif

set pmci 2
* exec $PER/s#vpl amp[nc]z damp[nc]z daysz ddaysz sz=0.01
sigma eff[nc]bz = n0[nc]bz/nx[nc]bz
sigma deff[nc]bz = sqrt(eff[nc]bz*(1-eff[nc]bz)/nx[nc]bz)
sigma amp[nc]pz = -log(eff[nc]bz) 
sigma damp[nc]pz = deff[nc]bz/eff[nc]bz

*sigma smu[nc]z = sqrt(2*(amp[nc]z-amp[nc]pz))/amp[nc]z
*sigma sigmu1 = ((amp[nc]z-2*amp[nc]pz)*damp[nc]z/amp[nc]z)**2
*sigma sigmu2 = (damp[nc]pz/amp[nc]pz)**2
*sigma dsmu[nc]z = sqrt(sigmu1+sigmu2)/smu[nc]z/amp[nc]z**2

sigma smu[nc]z = 2*(amp[nc]z-amp[nc]pz)/amp[nc]z**2
sigma sigmu1 = ((2*amp[nc]pz/amp[nc]z-1)*damp[nc]z)**2
sigma sigmu2 = (damp[nc]pz)**2
sigma dsmu[nc]z = 2*sqrt(sigmu1+sigmu2)/amp[nc]z**2
  
set pmci 4
* exec $PER/s#vpl amp[nc]pz damp[nc]pz daysz ddaysz sz=0.01 o=s
*read x
*
*ve/cre runsx(16) r 7842 10988 11285 13219 13515 13724 13725 13846 13847 14047 14048 16698 16699 17868 17869 21075
ve/cre runsx(8) r 7842 10988 11285 13846 13847 16698 16699 17868
exec mapcal#dcuts runsx daysx
ve/cre daysx1(2) r -30 180
ve/cre daysx2(2) r 410 710
ve/cre daysx3(4) r 695 750 760 790
ve/cre daysx4(2) r 786 835
if ([dir].eq.'vr') then
  ve/cre kt(4) r 1 1 2.5 1
else
  ve/cre kt(4) r 1 1 1 1
endif
ve/del act,dst
n=$vlen(daysz)
do i=1,4
  j=[i]*2-1
  tx=daysx([j])
  exec mapcal#ixndl [tx] daysz
  gl/imp inx
  in1=[inx]
  j=[i]*2
  tx=daysx([j])
  exec mapcal#ixndr [tx] daysz
  gl/imp inx
  in2=[inx]
  
  ve/del ampt,dampt,ampq,dampq,dayst,ddayst,smut,dsmut
  ve/copy  smu[nc]z([in1]:[in2])  smut
  ve/copy dsmu[nc]z([in1]:[in2]) dsmut
  ve/copy  amp[nc]z([in1]:[in2])  ampt
  ve/copy damp[nc]z([in1]:[in2]) dampt
  ve/copy  amp[nc]pz([in1]:[in2])  ampq
  ve/copy damp[nc]pz([in1]:[in2]) dampq
  ve/copy  daysz([in1]:[in2])  dayst
  ve/copy ddaysz([in1]:[in2]) ddayst
  kti=kt([i])
  sigma  ampt =  ampt * [kti]
  sigma dampt = dampt * [kti]
  sigma smut = 2*(ampt-ampq)/ampt**2
  sigma sigmu1 = ((2*ampq/ampt-1)*dampt)**2
  sigma sigmu2 = (dampq)**2
  sigma dsmut = 2*sqrt(sigmu1+sigmu2)/ampt**2
  set pmci 1
  set lwid 1
  set hwid 1
  * exec $PER/s#vpl ampt dampt dayst ddayst sz=0.1 iatt=20
  * exec $PER/s#vpl ampq dampq dayst ddayst sz=0.1 iatt=24 o=s
  
  nz=$vlen(daysx[i])
  nz=[nz]/2
  ve/del ampt0,dampt0,ampq0,dampq0,dayst0,ddayst0,smut0,dsmut0
  do k=1,[nz]
  
    j=[k]*2-1
    tx=daysx[i]([j])
    exec mapcal#ixndl [tx] dayst
    gl/imp inx
    in1=[inx]
    j=[k]*2
    tx=daysx[i]([j])
    exec mapcal#ixndr [tx] dayst
    gl/imp inx
    in2=[inx]
    
    ve/del amptt,damptt,ampqt,dampqt,daystt,ddaystt,smutt,dsmutt
    ve/copy  smut([in1]:[in2])  smutt
    ve/copy dsmut([in1]:[in2]) dsmutt
    ve/copy  ampt([in1]:[in2])  amptt
    ve/copy dampt([in1]:[in2]) damptt
    ve/copy  ampq([in1]:[in2])  ampqt
    ve/copy dampq([in1]:[in2]) dampqt
    ve/copy  dayst([in1]:[in2])  daystt
    ve/copy ddayst([in1]:[in2]) ddaystt
    
    exec vappend  smut0  smutt
    exec vappend dsmut0 dsmutt
    exec vappend  ampt0  amptt
    exec vappend dampt0 damptt
    exec vappend  ampq0  ampqt
    exec vappend dampq0 dampqt
    exec vappend dayst0 daystt
    exec vappend ddayst0 ddaystt

  enddo
  
  nz=$vlen(dayst0)
  ve/cre is([nz]) r
  sigma is = array([nz],1#[nz])
  nz0=[nz]
  do l=1,[nz]
    ampl=ampq0([l])
    if (([ampl].lt.2).or.([ampl].gt.15)) then
      ve/inp is([l]) $sigma(is([l])+[nz])
      nz0=[nz0]-1
    endif
  enddo
  ve/del amptf,damptf,ampqf,dampqf,daystf,ddaystf,smutf,dsmutf
  sigma  smutf = order( smut0,is)
  sigma dsmutf = order(dsmut0,is)
  sigma  amptf = order( ampt0,is)
  sigma damptf = order(dampt0,is)
  sigma  ampqf = order( ampq0,is)
  sigma dampqf = order(dampq0,is)
  sigma  daystf = order( dayst0,is)
  sigma ddaystf = order(ddayst0,is)
  exec mapcal#vecut  smutf  [nz0]
  exec mapcal#vecut dsmutf  [nz0]
  exec mapcal#vecut  amptf  [nz0]
  exec mapcal#vecut damptf  [nz0]
  exec mapcal#vecut  ampqf  [nz0]
  exec mapcal#vecut dampqf  [nz0]
  exec mapcal#vecut  daystf [nz0]
  exec mapcal#vecut ddaystf [nz0]
  mess nz: [nz] [nz0]
 
  ve/del dampqfx
  ve/copy dampqf dampqfx
  nz=$vlen(daystf)
  ve/cre is([nz]) r
  sigma is = array([nz],1#[nz])
  nz0=[nz]
  do ll=1,1
    mean=$sigma(vsum(ampqf)/[nz])
    mean2=$sigma(vsum(ampqf**2)/[nz])
    rms=$sigma(sqrt([mean2]-[mean]**2))
    mess mean=[mean] rms=[rms]
    do l=1,[nz]
      if ($sigma(abs(ampqf([l])-[mean])/[rms]).gt.3) then
        ve/inp dampqfx([l]) 10
        ve/inp is([l]) $sigma(is([l])+[nz])
        mess ll=[ll] l=[l] 
        nz0=[nz0]-1
      endif
    enddo
  enddo
  sigma  smutf = order( smutf,is)
  sigma dsmutf = order(dsmutf,is)
  sigma  amptf = order( amptf,is)
  sigma damptf = order(damptf,is)
  sigma  ampqf = order( ampqf,is)
  sigma dampqf = order(dampqf,is)
  sigma  daystf = order( daystf,is)
  sigma ddaystf = order(ddaystf,is)
  exec mapcal#vecut  smutf  [nz0]
  exec mapcal#vecut dsmutf  [nz0]
  exec mapcal#vecut  amptf  [nz0]
  exec mapcal#vecut damptf  [nz0]
  exec mapcal#vecut  ampqf  [nz0]
  exec mapcal#vecut dampqf  [nz0]
  exec mapcal#vecut  daystf [nz0]
  exec mapcal#vecut ddaystf [nz0]
  mess nz: [nz] [nz0]

  do l=1,-3
    ve/fit daystf ampqf dampqf e s
  enddo
  if ([i].eq.1) then
    npar=7
    ve/cre chi2(2) r
    ve/cre paru([npar]) r
    ve/cre dparu([npar]) r
    ve/cre covu([npar],[npar]) r
    fname=[dir]/mapcal[i]_fit_amp_vs_[msum]_counter0.txt
    w=$sigma(vmax(daystf)-vmin(daystf))
    ve/cre dp7(7) r                  
    if ($fexist([fname]).eq.1) then
      ve/read p7,dp7 [fname] '2f15.6'
    else
      ve/cre p7(7) r $sigma([mean]-0.5) +1 $sigma(vsum(daystf)/[nz0]-[w]/3) $sigma([w]/8) _
                                        -1 $sigma(vsum(daystf)/[nz0]+[w]/3) $sigma([w]/4)
    endif
    if ([nc].eq.0) then
      ve/cre s7(7) r $sigma([mean]/100) 0.01 $sigma(vsum(daystf)/100) $sigma([w]/100) _
                                        0.01 $sigma(vsum(daystf)/100) $sigma([w]/100)
      ve/cre pmin(7) r $sigma([mean]/2) $sigma(-[mean]/2) $sigma(vmin(daystf)) $sigma(30) _
                                        $sigma(-[mean]/2) $sigma(vmin(daystf)) $sigma(30)
      ve/cre pmax(7) r $sigma([mean]*2) $sigma(+[mean]/2) $sigma(vmax(daystf)) $sigma([w]) _
                                        $sigma(+[mean]/2) $sigma(vmax(daystf)) $sigma([w])
    else
      ve/cre s7(7) r $sigma([mean]/100) 0.01 0 0 0.01 0 0
      ve/del pmin,pmax
      ve/cre st7(7) r 2 1 0 0 1 0 0
      sigma pmin = p7-abs(dp7)-st7
      sigma pmax = p7+abs(dp7)+st7
      sigma s7 = abs(dp7)/10
    endif
    
    ve/cre s7t(7) r 0.1 0.1 $sigma(abs(p7(4))) $sigma(abs(p7(4))/2) 0.1 $sigma(abs(p7(7))) $sigma(abs(p7(7))/2)
    do l=1,3
      ve/fit daystf ampqf dampqf derf0.f sb 7 p7 s7t pmin pmax dp7
    enddo
    
    ve/cre xf(1) r
    do nf=1,3
    
    do l=1,3
      ve/fit daystf ampqf dampqf derf0.f sb 7 p7 s7 pmin pmax dp7
    enddo
    
    nz=$vlen(dayst0)
    ve/cre is([nz]) r
    sigma is = array([nz],1#[nz])
    nz0=[nz]
    do ll=1,1
      do l=1,[nz]
        ve/inp xf(1) $sigma(dayst0([l]))
        derf0pl=$call('derf0p.f(xf)')
        if ($sigma(abs(ampq0([l])-[derf0pl])/dampq0([l])).gt.3) then
          ve/inp is([l]) $sigma(is([l])+[nz])
          nz0=[nz0]-1
        endif
      enddo
    enddo
    ve/del amptf,damptf,ampqf,dampqf,daystf,ddaystf,smutf,dsmutf
    sigma  smutf = order( smut0,is)
    sigma dsmutf = order(dsmut0,is)
    sigma  amptf = order( ampt0,is)
    sigma damptf = order(dampt0,is)
    sigma  ampqf = order( ampq0,is)
    sigma dampqf = order(dampq0,is)
    sigma  daystf = order( dayst0,is)
    sigma ddaystf = order(ddayst0,is)
    exec mapcal#vecut  smutf  [nz0]
    exec mapcal#vecut dsmutf  [nz0]
    exec mapcal#vecut  amptf  [nz0]
    exec mapcal#vecut damptf  [nz0]
    exec mapcal#vecut  ampqf  [nz0]
    exec mapcal#vecut dampqf  [nz0]
    exec mapcal#vecut  daystf [nz0]
    exec mapcal#vecut ddaystf [nz0]
    enddo
    
    set pmci 1
    set lwid 1
    set hwid 1
    * exec $PER/s#vpl ampt dampt dayst ddayst sz=0.1 iatt=20
    * exec $PER/s#vpl ampq dampq dayst ddayst sz=0.1 iatt=24 o=s
    if ([nc].eq.0) then
      opt=s
    else
      opt=sb
    endif
    do l=1,3
      ve/fit daystf ampqf dampqf derf0.f [opt] 7 p7 s7 pmin pmax dp7
    enddo
    
    fname=[dir]/mapcal[i]_fit_amp_vs_[msum]_counter[nc].txt
    ve/write p7,dp7 [fname] '2f15.6'
    ve/write p7,dp7 ! '2f15.6'
    
  endif
  if ([i].ne.1) then
    npar=4
    ve/cre chi2(2) r
    ve/cre paru([npar]) r
    ve/cre dparu([npar]) r
    ve/cre covu([npar],[npar]) r
    
    fname=[dir]/mapcal[i]_fit_amp_vs_[msum]_counter0.txt
    w=$sigma(vmax(daystf)-vmin(daystf))
    ve/cre dp4(4) r
    if ($fexist([fname]).eq.1) then
      ve/read p4,dp4 [fname] '2f15.6'
    else
      ve/cre p4(4) r $sigma([mean]+0.5) -1 $sigma(vsum(daystf)/[nz0]) $sigma((vmax(daystf)-vmin(daystf))/4)
    endif
    if ([nc].eq.0) then
      ve/cre s4(4) r $sigma([mean]/100) 0.01 $sigma(vsum(daystf)/100) $sigma([w]/100)
      ve/cre pmin(4) r $sigma([mean]/2) $sigma(-[mean]/2) $sigma(vmin(daystf)) $sigma(5)
      ve/cre pmax(4) r $sigma([mean]*2) $sigma(+[mean]/2) $sigma(vmax(daystf)) $sigma(vmax(daystf)-vmin(daystf))
    else
      ve/cre s4(4) r $sigma([mean]/100) 0.01 0 0
      ve/del pmin,pmax
      ve/cre st4(4) r 2 1 $sigma(abs(p4(4))) $sigma(abs(p4(4))/2)
      sigma pmin = p4-st4
      sigma pmax = p4+st4
      sigma s4 = abs(dp4)/10
    endif

    ve/cre s4t(4) r 0.1 0.1 0 0
    do l=1,3
      ve/fit daystf ampqf dampqf erf0.f sb 4 p4 s4t pmin pmax dp4
    enddo

    ve/cre xf(1) r
    do nf=1,3
    
    do l=1,3
      ve/fit daystf ampqf dampqf erf0.f sb 4 p4 s4 pmin pmax dp4
    enddo
    
    nz=$vlen(dayst0)
    ve/cre is([nz]) r
    sigma is = array([nz],1#[nz])
    nz0=[nz]
    do ll=1,1
      do l=1,[nz]
        ve/inp xf(1) $sigma(dayst0([l]))
        erf0pl=$call('erf0p.f(xf)')
        if ($sigma(abs(ampq0([l])-[erf0pl])/dampq0([l])).gt.5) then
          ve/inp is([l]) $sigma(is([l])+[nz])
          nz0=[nz0]-1
        endif
      enddo
    enddo
    ve/del amptf,damptf,ampqf,dampqf,daystf,ddaystf,smutf,dsmutf
    sigma  smutf = order( smut0,is)
    sigma dsmutf = order(dsmut0,is)
    sigma  amptf = order( ampt0,is)
    sigma damptf = order(dampt0,is)
    sigma  ampqf = order( ampq0,is)
    sigma dampqf = order(dampq0,is)
    sigma  daystf = order( dayst0,is)
    sigma ddaystf = order(ddayst0,is)
    exec mapcal#vecut  smutf  [nz0]
    exec mapcal#vecut dsmutf  [nz0]
    exec mapcal#vecut  amptf  [nz0]
    exec mapcal#vecut damptf  [nz0]
    exec mapcal#vecut  ampqf  [nz0]
    exec mapcal#vecut dampqf  [nz0]
    exec mapcal#vecut  daystf [nz0]
    exec mapcal#vecut ddaystf [nz0]
    enddo
    
    fname=[dir]/mapcal[i]_fit_amp_vs_[msum]_counter[nc].txt
    ve/write p4,dp4 [fname] '2f15.6'
    
    set pmci 1
    set lwid 1
    set hwid 1
    * exec $PER/s#vpl ampt dampt dayst ddayst sz=0.1 iatt=20
    * exec $PER/s#vpl ampq dampq dayst ddayst sz=0.1 iatt=24 o=s
    if ([nc].eq.0) then
      opt=s
    else
      opt=sb
    endif
    do l=1,3
      ve/fit daystf ampqf dampqf erf0.f [opt] 4 p4 s4 pmin pmax
    enddo
    
  endif
  
  set pmci 2
  * exec $PER/s#vpl ampqf dampqf daystf ddaystf sz=0.1 iatt=24 o=s
*  
  ve/del p7t,dp7t,p4t,dp4t,p2t,dp2t
  ve/copy p7 p7t
  ve/copy p4 p4t
  ve/copy p2 p2t
  ve/copy dp7 dp7t
  ve/copy dp4 dp4t
  ve/copy dp2 dp2t
    if ([i].eq.1) then
      fname=[dir]/mapcal[i]_fit_amp_vs_[msum]_counter[nc].txt
      ve/read p7,dp7 [fname] '2f15.6'
*      fun/pl derf0p.f $GRAFINFO('WNXMIN') $GRAFINFO('WNXMAX') s
      ve/read p4,dp4 mapcal[i]_sigmu2_fit_[msum]_counter[nc].txt '(2f15.6)'
      ve/inp p4(2) 0
      ve/inp p4(4) $sigma(vmin(dayst0))
      fun/pl mucorr.f $GRAFINFO('WNXMIN') $GRAFINFO('WNXMAX') s
      npi=$vlen(dayst)
      ve/cre aci([npi]) r 
      ve/cre xpi(1) r
      do ip=1,[npi]
        ve/inp xpi(1) $sigma(dayst([ip]))
        ai=$call('mucorr.f(xpi)')
        ve/inp aci([ip]) $sigma([ai]/ampt([ip]))
      enddo
      exec vappend act aci
      exec vappend dct dayst
    endif
    if ([i].eq.3) then
      ncx=[nc]
      if ([nc].eq.9) then
        ncx=0
      endif
      fname=[dir]/mapcal[i]_fit_amp_vs_[msum]_counter[ncx].txt
      ve/read p4,dp4 [fname] '2f15.6'
*      fun/pl derf0p.f $GRAFINFO('WNXMIN') $GRAFINFO('WNXMAX') s
      ve/read p4s,dp4s mapcal[i]_sigmu2_fit_[msum]_counter[ncx].txt '(2f15.6)'
      ve/inp p4s(2) 0
      ve/inp p4s(4) $sigma(vmin(dayst0))
      fun/pl mucorr3.f $GRAFINFO('WNXMIN') $GRAFINFO('WNXMAX') s
      npi=$vlen(dayst)
      ve/cre aci([npi]) r 
      ve/cre xpi(1) r
      do ip=1,[npi]
        ve/inp xpi(1) $sigma(dayst([ip]))
        ai=$call('mucorr3.f(xpi)')
        ve/inp aci([ip]) $sigma([ai]/ampt([ip]))
      enddo
      exec vappend act aci
      exec vappend dct dayst
    endif
    if ([i].eq.2) then
      fname=[dir]/mapcal[i]_fit_amp_vs_[msum]_counter[nc].txt
      ve/read p4,dp4 [fname] '2f15.6'
*      fun/pl erf0p.f $GRAFINFO('WNXMIN') $GRAFINFO('WNXMAX') s
      ve/read p2,dp2 mapcal[i]_sigmu2_fit_[msum]_counter[nc].txt '(2f15.6)'
      fun/pl mucorr2.f $GRAFINFO('WNXMIN') $GRAFINFO('WNXMAX') s
      npi=$vlen(dayst)
      ve/cre aci([npi]) r 
      ve/cre xpi(1) r
      do ip=1,[npi]
        ve/inp xpi(1) $sigma(dayst([ip]))
        ai=$call('mucorr2.f(xpi)')
        ve/inp aci([ip]) $sigma([ai]/ampt([ip]))
      enddo
      exec vappend act aci
      exec vappend dct dayst
    endif
    if ([i].eq.4) then
      fname=[dir]/mapcal[i]_fit_amp_vs_[msum]_counter[nc].txt
      ve/read p4,dp4 [fname] '2f15.6'
*      fun/pl erf0p.f $GRAFINFO('WNXMIN') $GRAFINFO('WNXMAX') s
      ve/read p2,dp2 mapcal[i]_sigmu2_fit_[msum]_counter[nc].txt '(2f15.6)'
      fun/pl mucorr2.f $GRAFINFO('WNXMIN') $GRAFINFO('WNXMAX') s
      npi=$vlen(dayst)
      ve/cre aci([npi]) r 
      ve/cre xpi(1) r
      do ip=1,[npi]
        ve/inp xpi(1) $sigma(dayst([ip]))
        ai=$call('mucorr2.f(xpi)')
        ve/inp aci([ip]) $sigma([ai]/ampt([ip]))
      enddo
      exec vappend act aci
      exec vappend dct dayst
    endif
    ve/write dayst,aci mapcal[i]_ampcorr_[msum]_counter[nc].txt '(2f15.6)'
    sigma amptc = ampt*aci
    sigma damptc = 0*aci
    set pmci 4
    * exec $PER/s#vpl amptc damptc dayst ddayst sz=0.03 iatt=20 o=s
*    read x
  ve/copy p7t p7
  ve/copy p4t p4
  ve/copy p2t p2
  ve/copy dp7t dp7
  ve/copy dp4t dp4
  ve/copy dp2t dp2
*    
  set pmci 1
  call covm.f(1)
  call covmpen(chi2,[npar],paru,dparu)
  call covmcov([npar],covu)
  ndf=$sigma(chi2(2))
  chi=$sigma(int(chi2(1)*[ndf]*100+0.5)/100)
  txt=[h]^2!/ndf = [chi] / [ndf]
  exec $PER/s#tf 0.1 0.9 [txt]
  atitle 't, days' '[m], pe.' 
  exec save mapcal[i]_amp_vs_[msum]_counter[nc].eps f
  
  set pmci 1
  set lwid 1
  set hwid 1
  * exec $PER/s#vpl ampt dampt dayst ddayst sz=0.1 iatt=20
  * exec $PER/s#vpl ampq dampq dayst ddayst sz=0.1 iatt=24 o=s
  atitle 't, days' '[m], pe.' 
  exec save mapcal[i]_amp_vs_[msum]_0_counter[nc].eps f

  ve/del p4e
  ve/copy p4 p4e
  ve/cre tr(1) r 751
*  read x
  
*  sigma sigmu0 = sqrt(2*(ampt0-ampq0))/ampt0
*  sigma sigmu1 = ((ampt0-2*ampq0)*dampt0/ampt0)**2
*  sigma sigmu2 = (dampq0/ampq0)**2
*  sigma dsigmu0 = sqrt(sigmu1+sigmu2)/sigmu0/ampt0**2
  ve/del sigmu0,dsigmu0
  ve/copy smut0 sigmu0
  ve/copy dsmut0 dsigmu0
  ve/del sigmuz,dsigmuz,daystz,ddaystz
  ve/copy sigmu0 sigmuz
  ve/copy dsigmu0 dsigmuz
  ve/copy dayst0 daystz
  ve/copy ddayst0 ddaystz
  
  do nf=1,2
    nz=$vlen(daystz)
    ve/cre is([nz]) r
    sigma is = array([nz],1#[nz])
    nz0=[nz]
    
    mean=0.05
    rms=0.02
    do ll=1,1
      mess mean=[mean] rms=[rms]
      do l=1,[nz]
        if ($sigma(sigmuz([l])-dsigmuz([l])).gt.0.1) then
          ve/inp is([l]) $sigma(is([l])+[nz])
          nz0=[nz0]-1
        endif
        if (($sigma(abs(sigmuz([l])-[mean])/[rms]).gt.3000).or.($sigma(dsigmuz([l])).eq.0)) then
          ve/inp is([l]) $sigma(is([l])+[nz])
          nz0=[nz0]-1
        endif
        if (([i].eq.1).and.($sigma(daystz([l])).lt.-10)) then
          ve/inp is([l]) $sigma(is([l])+[nz])
          nz0=[nz0]-1
        endif
        if (([i].eq.3).and.($sigma(daystz([l])).lt.705)) then
          ve/inp is([l]) $sigma(is([l])+[nz])
          nz0=[nz0]-1
        endif
        if (([i].eq.3).and.($sigma(daystz([l])).gt.$sigma(tr(1)))) then
*          ve/inp is([l]) $sigma(is([l])+[nz])
*          nz0=[nz0]-1
        endif
      enddo
    enddo
*    ve/del sigmuz,dsigmuz,daystz,ddaystz
    sigma  sigmuz = order( sigmuz,is)
    sigma dsigmuz = order(dsigmuz,is)
    sigma  daystz = order( daystz,is)
    sigma ddaystz = order(ddaystz,is)
    exec mapcal#vecut  sigmuz [nz0]
    exec mapcal#vecut dsigmuz [nz0]
    exec mapcal#vecut  daystz [nz0]
    exec mapcal#vecut ddaystz [nz0]
    nz=$vlen(daystz)
    mean=$sigma(vsum(sigmuz)/[nz])
    mean2=$sigma(vsum(sigmuz**2)/[nz])
    rms=$sigma(sqrt([mean2]-[mean]**2))
  enddo
  nz=$vlen(daystz)
  mean=$sigma(vsum(sigmuz)/[nz])
  mean2=$sigma(vsum(sigmuz**2)/[nz])
  rms=$sigma(sqrt([mean2]-[mean]**2))
  mess mean=[mean] rms=[rms]
  ve/cre p3(3) r [mean] [rms] [nz]
  ve/write p3 mapcal[i]_sigmu_pars_vs_[msum]_counter[nc].txt  
  rms=[rms]/2
*  sigma dsigmuz = array([nz],[rms]#[rms])
  
  * exec $PER/s#vpl sigmuz dsigmuz daystz ddaystz sz=0.1 iatt=20
  ve/cre xf(1) r
  ve/cre p2(2) r 0.25 0
  ve/cre dp2(2) r
  ve/del sigmuf,dsigmuf,daystf,ddaystf
  ve/copy  sigmuz  sigmuf
  ve/copy dsigmuz dsigmuf
  ve/copy  daystz  daystf
  ve/copy ddaystz ddaystf
  ve/cre p4(4) r 0.25 0.1 30 0.9
  ve/cre dp4(4) r
  sigma dsigmux = dsigmu0 + abs(1-dsigmu0/dsigmu0)
  ve/cre t0r(1) r $sigma(vmin(dayst0))
  ve/write t0r mapcal[i]_t0r_counter[nc].txt
  do nf=1,3
    
    set pmci 1
    * exec $PER/s#vpl sigmuz dsigmuz daystz ddaystz sz=0.1 iatt=20
    do l=1,3
      if ([i].eq.1) then
        null $GRAFINFO('WNXMIN') $GRAFINFO('WNXMAX') $sigma([mean]-5*[rms]) $sigma([mean]+5*[rms])
        * exec $PER/s#vpl sigmuz dsigmuz daystz ddaystz sz=0.1 iatt=20 ll=-1
        fname=[dir]/mapcal[i]_fit_amp_vs_[msum]_counter[nc].txt
        ve/read p7e,dp7e [fname] '2f15.6'
        if ([nc].eq.0) then
          if ([l].eq.1) then
          ve/cre p4(4) r 0.25 0.1 30 1
          ve/cre dp4(4) r
          ve/cre s4(4) r 0.01 0 0 0
          ve/fit daystf sigmuf dsigmuf sigmu2d.f sb 4 p4 s4 
          ve/cre s4(4) r 0.01 0 0 0
          ve/fit daystf sigmuf dsigmuf sigmu2d.f sb 4 p4 s4
          ve/cre s4(4) r 0.01 0.01 0 0
          ve/fit daystf sigmuf dsigmuf sigmu2d.f sb 4 p4 s4
          ve/cre s4(4) r 0.01 0.01 1 0
          ve/fit daystf sigmuf dsigmuf sigmu2d.f sb 4 p4 s4
          endif
          ve/fit daystf sigmuf dsigmuf sigmu2d.f sb 4 p4 s4
*          read x
        else
          if ([l].eq.1) then
          ve/del p4,dp4
          ve/read p4,dp4 mapcal[i]_sigmu2_fit_[msum]_counter0.txt '(2f15.6)'
          ve/cre s4(4) r 0.01 0 0 0
          ve/fit daystf sigmuf dsigmuf sigmu2d.f sb 4 p4 s4
          ve/cre s4(4) r 0.01 0 0 0
          ve/fit daystf sigmuf dsigmuf sigmu2d.f sb 4 p4 s4
          ve/cre s4(4) r 0.01 0.01 0 0
          ve/fit daystf sigmuf dsigmuf sigmu2d.f sb 4 p4 s4
          endif
          ve/fit daystf sigmuf dsigmuf sigmu2d.f sb 4 p4 s4
*          ve/cre pmin(4) r $sigma([mean]-3*[rms]) $sigma(0) 0 0.5
*          ve/cre pmax(4) r $sigma([mean]+5*[rms]) $sigma(3*[rms]) 1000 1.0
*          ve/fit daystf sigmuf dsigmuf sigmu2d.f sb 4 p4 s4 pmin pmax dp4
        endif
        ve/write p4,dp4 mapcal[i]_sigmu2_fit_[msum]_counter[nc].txt '(2f15.6)'
      endif
      if ([i].eq.3) then
        null $GRAFINFO('WNXMIN') $GRAFINFO('WNXMAX') $sigma([mean]-5*[rms]) $sigma([mean]+5*[rms])
        * exec $PER/s#vpl sigmuz dsigmuz daystz ddaystz sz=0.1 iatt=20 ll=-1 o=s
        fname=[dir]/mapcal[i]_fit_amp_vs_[msum]_counter[nc].txt
        ve/write p4e,dp4e [fname] '2f15.6'
        if ([nc].eq.0) then
          if ([l].eq.1) then
          ve/cre p4(4) r 0.25 0.1 200 0.9
          ve/cre dp4(4) r
          ve/cre s4(4) r 0.01 0 0 0
          ve/fit daystf sigmuf dsigmuf sigmu2.f sb 4 p4 s4
          ve/cre s4(4) r 0.01 0 0 0.01
          ve/fit daystf sigmuf dsigmuf sigmu2.f sb 4 p4 s4
          ve/cre s4(4) r 0.01 0.01 0 0.01
          ve/fit daystf sigmuf dsigmuf sigmu2.f sb 4 p4 s4
          ve/cre s4(4) r 0 0 1 0
          ve/fit daystf sigmuf dsigmuf sigmu2.f sb 4 p4 s4
          endif
          ve/fit daystf sigmuf dsigmuf sigmu2.f s 4 p4 ! ! ! dp4
*          read x
        else
          if ([l].eq.1) then
          ve/del p4,dp4
          ve/read p4,dp4 mapcal[i]_sigmu2_fit_[msum]_counter0.txt '(2f15.6)'
          ve/cre s4(4) r 0.01 0 0 0
          ve/fit daystf sigmuf dsigmuf sigmu2.f sb 4 p4 s4
          ve/cre s4(4) r 0.01 0 0 0.01
          ve/fit daystf sigmuf dsigmuf sigmu2.f sb 4 p4 s4
          ve/cre s4(4) r 0.01 0.01 0 0.01
          ve/fit daystf sigmuf dsigmuf sigmu2.f sb 4 p4 s4
          endif
          ve/fit daystf sigmuf dsigmuf sigmu2.f sb 4 p4 s4
*          ve/cre pmin(4) r $sigma([mean]-3*[rms]) $sigma(0) 0 0.5
*          ve/cre pmax(4) r $sigma([mean]+5*[rms]) $sigma(3*[rms]) 1000 1.0
*          ve/fit daystf sigmuf dsigmuf sigmu2.f sb 4 p4 s4 pmin pmax dp4
        endif
        ve/write p4,dp4 mapcal[i]_sigmu2_fit_[msum]_counter[nc].txt '(2f15.6)'
      endif
      if (([i].eq.2).or.([i].eq.4)) then
*        if ([nf].eq.1) then
*          ve/fit daystf sigmuf dsigmuf p0 s 1 p2
*        else
          ve/cre s2(2) r 0.001 0
          ve/fit daystf sigmuf dsigmuf p1.f sb 2 p2 s2
*        endif
        ve/write p2,dp2 mapcal[i]_sigmu2_fit_[msum]_counter[nc].txt '(2f15.6)'
      endif
    enddo
*    read x
    
    nz=$vlen(daystz)
    ve/cre is([nz]) r
    sigma is = array([nz],1#[nz])
    nz0=[nz]
    if (([nc].eq.0).and.([i].eq.1)) then
      nsig=3
    else
      nsig=10
    endif
    do ll=1,1
      do l=1,[nz]
        ve/inp xf(1) $sigma(daystz([l]))
        a=p2(1)
        b=p2(2)
        p2l=$sigma(([a])+xf(1)*([b]))
        p2l=$call('p1p.f(xf)')
        if (([i].eq.3).or.([i].eq.1)) then
          p2l=$call('sigmu2p.f(xf)')
        endif
*        if ($sigma(abs(sigmuz([l])-[p2l])/[rms]).gt.3) then
        if ($sigma(abs(sigmuz([l])-[p2l])/dsigmuz([l])).gt.[nsig]) then
          ve/inp is([l]) $sigma(is([l])+[nz])
          nz0=[nz0]-1
        endif
      enddo
    enddo
    ve/del sigmuf,dsigmuf,daystf,ddaystf
    sigma  sigmuf = order( sigmuz,is)
    sigma dsigmuf = order(dsigmuz,is)
    sigma  daystf = order( daystz,is)
    sigma ddaystf = order(ddaystz,is)
    exec mapcal#vecut  sigmuf [nz0]
    exec mapcal#vecut dsigmuf [nz0]
    exec mapcal#vecut  daystf [nz0]
    exec mapcal#vecut ddaystf [nz0]
    
    set pmci 2
    * exec $PER/s#vpl sigmuf dsigmuf daystf ddaystf sz=0.1 iatt=20 o=s
*    read x
  enddo
  atitle 't, days' '[s]?[m]!/[m]'
  exec save mapcal[i]_sigmu_vs_[msum]_counter[nc].eps f
*  read x
enddo
*ve/write dct,act mapcal_ampcorr_[msum]_counter[nc].txt '(2f15.6)'
shell cat mapcal?_ampcorr_[msum]_counter[nc].txt > mapcal_ampcorr_[msum]_counter[nc].txt
return




macro ixndl tx=0 v=daysz
  ve/del va,ix
  sigma va = abs([v]-([tx]))
  n=$vlen([v])
  ve/cre ix([n]) r
  sigma ix = array([n],1#[n])
  sigma ixa = order(ix,va)
  ind=ixa(1)
  if ($sigma([v]([ind])).lt.[tx]) then
    ind=$sigma(min([ind]+1,[n]))
  endif
  gl/cre inx [ind]
return

macro ixndr tx=0 v=daysz
  ve/del va,ix
  sigma va = abs([v]-([tx]))
  n=$vlen([v])
  ve/cre ix([n]) r
  sigma ix = array([n],1#[n])
  sigma ixa = order(ix,va)
  ind=ixa(1)
  if ($sigma([v]([ind])).gt.[tx]) then
    ind=$sigma(max([ind]-1,1))
  endif
  gl/cre inx [ind]
return

macro calsum v=amp1 dv=damp1 r=runs d=days 
ve/del r1,r2
*ve/read r1,r2 accled_corr.ixlist
ve/read r1,r2 accled_ixlist_actcals.txt
n=$vlen(r1)
ve/inp r2([n]) 30000
ve/cre  [v]z([n]) r
ve/cre [dv]z([n]) r
ve/cre  [r]z([n]) r
ve/cre d[r]z([n]) r
ve/cre  [d]z([n]) r
ve/cre d[d]z([n]) r
ve/cre id1([n]) r
ve/cre id2([n]) r
dr=$sigma(-300+8.5/24)
ind=0
do i=1,[n]
  rb=r1([i])
  exec mapcal#ixndl [rb] [r]
  gl/imp inx
  in1=[inx]
  re=r2([i])
  exec mapcal#ixndr [re] [r]
  gl/imp inx
  in2=[inx]
  if ([in1].lt.[in2]) then
    ind=[ind]+1
    ve/inp id1([ind]) [in1]
    ve/inp id2([ind]) [in2]
    mess [in1] [in2]
*  
    j1=id1([ind])
    j2=id2([ind])
    ve/del [r]r,[v]r,[dv]r,[d]r
    ve/copy  [r]([j1]:[j2])  [r]r
    ve/copy  [d]([j1]:[j2])  [d]r
    ve/copy  [v]([j1]:[j2])  [v]r
    ve/copy [dv]([j1]:[j2]) [dv]r
    sigma  vm = vsum([v]r/([dv]r)**2)
    sigma dvm = vsum(1.0/([dv]r)**2)
    mean=$sigma(vm(1)/dvm(1))
    sig=$sigma(sqrt(1.0/dvm(1)))
    ve/inp  [v]z([ind]) [mean]
    ve/inp [dv]z([ind]) [sig]
    dmin=$sigma(vmin([d]r))
    dmax=$sigma(vmax([d]r))
    ve/inp  [d]z([ind]) $sigma(([dmax]+[dmin])/2)
    ve/inp d[d]z([ind]) $sigma(([dmax]-[dmin])/2)
    rmin=$sigma(vmin([r]r))
    rmax=$sigma(vmax([r]r))
    ve/inp  [r]z([ind]) $sigma(([rmax]+[rmin])/2)
    ve/inp d[r]z([ind]) $sigma(([rmax]-[rmin])/2)
  endif
enddo
exec mapcal#vecut [d]z
n=$vlen([d]z)
exec mapcal#vecut  [v]z [n]
exec mapcal#vecut [dv]z [n]
exec mapcal#vecut  [r]z [n]
exec mapcal#vecut d[r]z [n]
exec mapcal#vecut d[d]z [n]
return


macro calsumn v=nx r=runs d=days 
ve/del r1,r2
*ve/read r1,r2 accled_corr.ixlist
ve/read r1,r2 accled_ixlist_actcals.txt
n=$vlen(r1)
ve/inp r2([n]) 30000
ve/cre  [v]z([n]) r
ve/cre  [r]z([n]) r
ve/cre d[r]z([n]) r
ve/cre  [d]z([n]) r
ve/cre d[d]z([n]) r
ve/cre id1([n]) r
ve/cre id2([n]) r
ind=0
do i=1,[n]
  rb=r1([i])
  exec mapcal#ixndl [rb] [r]
  gl/imp inx
  in1=[inx]
  re=r2([i])
  exec mapcal#ixndr [re] [r]
  gl/imp inx
  in2=[inx]
  if ([in1].lt.[in2]) then
    ind=[ind]+1
    ve/inp id1([ind]) [in1]
    ve/inp id2([ind]) [in2]
*  
    j1=id1([ind])
    j2=id2([ind])
    ve/del [r]r,[v]r,[d]r
    ve/copy  [r]([j1]:[j2])  [r]r
    ve/copy  [d]([j1]:[j2])  [d]r
    ve/copy  [v]([j1]:[j2])  [v]r
    sigma  vm = vsum([v]r)
    mean=$sigma(vm(1))
    ve/inp  [v]z([ind]) [mean]
    dmin=$sigma(vmin([d]r))
    dmax=$sigma(vmax([d]r))
    ve/inp  [d]z([ind]) $sigma(([dmax]+[dmin])/2)
    ve/inp d[d]z([ind]) $sigma(([dmax]-[dmin])/2)
    rmin=$sigma(vmin([r]r))
    rmax=$sigma(vmax([r]r))
    ve/inp  [r]z([ind]) $sigma(([rmax]+[rmin])/2)
    ve/inp d[r]z([ind]) $sigma(([rmax]-[rmin])/2)
  endif
enddo
exec mapcal#vecut [d]z
n=$vlen([d]z)
exec mapcal#vecut  [v]z [n]
exec mapcal#vecut  [r]z [n]
exec mapcal#vecut d[r]z [n]
exec mapcal#vecut d[d]z [n]
return


macro daysum v=amp1 dv=damp1 r=runs d=days dc=dcuts
n=$vlen(days)
ve/cre  [v]z([n]) r
ve/cre [dv]z([n]) r
ve/cre  [r]z([n]) r
ve/cre d[r]z([n]) r
ve/cre  [d]z([n]) r
ve/cre d[d]z([n]) r
ve/cre id1([n]) r
ve/cre id2([n]) r
dr=$sigma(-300+8.5/24)
ind=0
do i=1,[n]
  ndays=[d]([i])
  if ([i].eq.[n]) then
    ndays=[ndays]+1
    io=[n]
  endif
  if ([ndays].gt.[dr]) then
    dr=$sigma(int([ndays])+8.5/24)
    if ([ndays].gt.[dr]) then
      dr=[dr]+1
    endif
    if ([ind].gt.0) then
      ve/inp id2([ind]) [io]
*
      j1=id1([ind])
      j2=id2([ind])
      ve/del [v]r,[dv]r,[d]r,[r]r
      ve/copy  [d]([j1]:[j2])  [d]r
      ve/copy  [r]([j1]:[j2])  [r]r
      ve/copy  [v]([j1]:[j2])  [v]r
      ve/copy [dv]([j1]:[j2]) [dv]r
*
      in1=0
      if ($vexist([dc])) then
        dmin=$sigma(vmin([d]r))
        dmax=$sigma(vmax([d]r))
        exec mapcal#ixndl [dmin] [dc]
        gl/imp inx
        inx1=[inx]
        exec mapcal#ixndr [dmax] [dc]
        gl/imp inx
        inx2=[inx]
        mess [j1] [j2] [dmin] [dmax] [inx1] [inx2]
        dcmin=$sigma([dc]([inx1]))
        dcmax=$sigma([dc]([inx2]))
        if (([dmin].le.[dcmin]).and.([dcmax].le.[dmax])) then
          mess $vlen([d]r)
          in1=1
          do j=[inx1],$sigma([inx2]+1)
            if ([j].le.[inx2]) then
              exec mapcal#ixndr $sigma([dc]([j])) [d]r
              gl/imp inx
              in2=[inx]
              dcmax=$sigma([dc]([j]))
            else
              in2=$vlen([d]r)
              dcmax=[dmax]
            endif
            if ([in1].le.[in2]) then
              mess --> [in1]:[in2] $sigma([d]r([in1])) $sigma([d]r([in2])) [dcmax]
*             
              ve/inp id1([ind]) $sigma([in1]+[j1]-1)
              ve/inp id2([ind]) $sigma([in2]+[j1]-1)
              ve/del [v]t,[dv]t,[d]t,[r]t
              ve/copy  [d]r([in1]:[in2])  [d]t
              ve/copy  [r]r([in1]:[in2])  [r]t
              ve/copy  [v]r([in1]:[in2])  [v]t
              ve/copy [dv]r([in1]:[in2]) [dv]t
*              
              sigma  vm = vsum([v]t/([dv]t)**2)
              sigma dvm = vsum(1.0/([dv]t)**2)
              mean=$sigma(vm(1)/dvm(1))
              sig=$sigma(sqrt(1.0/dvm(1)))
              ve/inp  [v]z([ind]) [mean]
              ve/inp [dv]z([ind]) [sig]
              dminr=$sigma(vmin([d]t))
              dmaxr=$sigma(vmax([d]t))
              ve/inp  [d]z([ind]) $sigma(([dmaxr]+[dminr])/2)
              ve/inp d[d]z([ind]) $sigma(([dmaxr]-[dminr])/2)
              rminr=$sigma(vmin([r]t))
              rmaxr=$sigma(vmax([r]t))
              ve/inp  [r]z([ind]) $sigma(([rmaxr]+[rminr])/2)
              ve/inp d[r]z([ind]) $sigma(([rmaxr]-[rminr])/2)
*      
              in1=[in2]+1
              ind=[ind]+1
            endif
          enddo
        endif
      endif
*      
      if ([in1].eq.0) then
        sigma  vm = vsum([v]r/([dv]r)**2)
        sigma dvm = vsum(1.0/([dv]r)**2)
        mean=$sigma(vm(1)/dvm(1))
        sig=$sigma(sqrt(1.0/dvm(1)))
        ve/inp  [v]z([ind]) [mean]
        ve/inp [dv]z([ind]) [sig]
        dmin=$sigma(vmin([d]r))
        dmax=$sigma(vmax([d]r))
        ve/inp  [d]z([ind]) $sigma(([dmax]+[dmin])/2)
        ve/inp d[d]z([ind]) $sigma(([dmax]-[dmin])/2)
        rmin=$sigma(vmin([r]r))
        rmax=$sigma(vmax([r]r))
        ve/inp  [r]z([ind]) $sigma(([rmax]+[rmin])/2)
        ve/inp d[r]z([ind]) $sigma(([rmax]-[rmin])/2)
*        
        ind=[ind]+1
      endif
*      read x
*
*      sigma  vm = vsum([v]r/([dv]r)**2)
*      sigma dvm = vsum(1.0/([dv]r)**2)
*      mean=$sigma(vm(1)/dvm(1))
*      sig=$sigma(sqrt(1.0/dvm(1)))
*      ve/inp  [v]z([ind]) [mean]
*      ve/inp [dv]z([ind]) [sig]
*      dmin=$sigma(vmin([d]r))
*      dmax=$sigma(vmax([d]r))
*      ve/inp  [d]z([ind]) $sigma(([dmax]+[dmin])/2)
*      ve/inp d[d]z([ind]) $sigma(([dmax]-[dmin])/2)
*      rmin=$sigma(vmin([r]r))
*      rmax=$sigma(vmax([r]r))
*      ve/inp  [r]z([ind]) $sigma(([rmax]+[rmin])/2)
*      ve/inp d[r]z([ind]) $sigma(([rmax]-[rmin])/2)
*
      ve/inp id1([ind]) [i]
    else
      ind=[ind]+1
      ve/inp id1([ind]) [i]
    endif
  endif
  io=[i]
enddo
exec mapcal#vecut [d]z
n=$vlen([d]z)
exec mapcal#vecut  [v]z [n]
exec mapcal#vecut [dv]z [n]
exec mapcal#vecut  [r]z [n]
exec mapcal#vecut d[r]z [n]
exec mapcal#vecut d[d]z [n]
return


macro daysumn v=nx r=runs d=days dc=dcuts
n=$vlen(days)
ve/cre  [v]z([n]) r
ve/cre  [d]z([n]) r
ve/cre d[d]z([n]) r
ve/cre  [r]z([n]) r
ve/cre d[r]z([n]) r
dr=$sigma(-300+8.5/24)
ind=0
do i=1,[n]
  ndays=[d]([i])
  if ([i].eq.[n]) then
    ndays=[ndays]+1
    io=[n]
  endif
  if ([ndays].gt.[dr]) then
    dr=$sigma(int([ndays])+8.5/24)
    if ([ndays].gt.[dr]) then
      dr=[dr]+1
    endif
    if ([ind].gt.0) then
      ve/inp id2([ind]) [io]
*
      j1=id1([ind])
      j2=id2([ind])
      ve/del [v]r,[d]r,[r]r
      ve/copy  [d]([j1]:[j2])  [d]r
      ve/copy  [r]([j1]:[j2])  [r]r
      ve/copy  [v]([j1]:[j2])  [v]r
*
      in1=0
      if ($vexist([dc])) then
        dmin=$sigma(vmin([d]r))
        dmax=$sigma(vmax([d]r))
        exec mapcal#ixndl [dmin] [dc]
        gl/imp inx
        inx1=[inx]
        exec mapcal#ixndr [dmax] [dc]
        gl/imp inx
        inx2=[inx]
        mess [j1] [j2] [dmin] [dmax] [inx1] [inx2]
        dcmin=$sigma([dc]([inx1]))
        dcmax=$sigma([dc]([inx2]))
        if (([dmin].le.[dcmin]).and.([dcmax].le.[dmax])) then
          mess $vlen([d]r)
          in1=1
          do j=[inx1],$sigma([inx2]+1)
            if ([j].le.[inx2]) then
              exec mapcal#ixndr $sigma([dc]([j])) [d]r
              gl/imp inx
              in2=[inx]
              dcmax=$sigma([dc]([j]))
            else
              in2=$vlen([d]r)
              dcmax=[dmax]
            endif
            if ([in1].le.[in2]) then
              mess --> [in1]:[in2] $sigma([d]r([in1])) $sigma([d]r([in2])) [dcmax]
*             
              ve/inp id1([ind]) $sigma([in1]+[j1]-1)
              ve/inp id2([ind]) $sigma([in2]+[j1]-1)
              ve/del [v]t,[d]t,[r]t
              ve/copy  [d]r([in1]:[in2])  [d]t
              ve/copy  [r]r([in1]:[in2])  [r]t
              ve/copy  [v]r([in1]:[in2])  [v]t
*              
              sigma  vm = vsum([v]t)
              mean=$sigma(vm(1))
              ve/inp  [v]z([ind]) [mean]
              dminr=$sigma(vmin([d]t))
              dmaxr=$sigma(vmax([d]t))
              ve/inp  [d]z([ind]) $sigma(([dmaxr]+[dminr])/2)
              ve/inp d[d]z([ind]) $sigma(([dmaxr]-[dminr])/2)
              rminr=$sigma(vmin([r]t))
              rmaxr=$sigma(vmax([r]t))
              ve/inp  [r]z([ind]) $sigma(([rmaxr]+[rminr])/2)
              ve/inp d[r]z([ind]) $sigma(([rmaxr]-[rminr])/2)
*      
              in1=[in2]+1
              ind=[ind]+1
            endif
          enddo
        endif
      endif
*      
      if ([in1].eq.0) then
        sigma  vm = vsum([v]r)
        mean=$sigma(vm(1))
        ve/inp  [v]z([ind]) [mean]
        dmin=$sigma(vmin([d]r))
        dmax=$sigma(vmax([d]r))
        ve/inp  [d]z([ind]) $sigma(([dmax]+[dmin])/2)
        ve/inp d[d]z([ind]) $sigma(([dmax]-[dmin])/2)
        rmin=$sigma(vmin([r]r))
        rmax=$sigma(vmax([r]r))
        ve/inp  [r]z([ind]) $sigma(([rmax]+[rmin])/2)
        ve/inp d[r]z([ind]) $sigma(([rmax]-[rmin])/2)
*        
        ind=[ind]+1
      endif
*
*      sigma  vm = vsum([v]r)
*      mean=$sigma(vm(1))
*      ve/inp  [v]z([ind]) [mean]
*      dmin=$sigma(vmin([d]r))
*      dmax=$sigma(vmax([d]r))
*      ve/inp  [d]z([ind]) $sigma(([dmax]+[dmin])/2)
*      ve/inp d[d]z([ind]) $sigma(([dmax]-[dmin])/2)
*      rmin=$sigma(vmin([r]r))
*      rmax=$sigma(vmax([r]r))
*      ve/inp  [r]z([ind]) $sigma(([rmax]+[rmin])/2)
*      ve/inp d[r]z([ind]) $sigma(([rmax]-[rmin])/2)
*
      ve/inp id1([ind]) [i]
    else
      ind=[ind]+1
      ve/inp id1([ind]) [i]
    endif
  endif
  io=[i]
enddo
exec mapcal#vecut [d]z
n=$vlen([d]z)
exec mapcal#vecut  [v]z [n]
exec mapcal#vecut d[d]z [n]
exec mapcal#vecut  [r]z [n]
exec mapcal#vecut d[r]z [n]
return


macro dcuts r=runsx d=daysx
n=$vlen([r])
ve/cre [d]([n])
do i=1,[n]
  exec mapcal#ndays $sigma([r]([i]))
  gl/imp ndays
  ve/inp [d]([i]) [ndays]
*  line [ndays] 0 [ndays] 15
enddo
return


macro sigmu map=1
ve/cre smu(9) r
ve/cre rmu(9) r
ve/cre nx(9) r
do i=1,9
  ve/read p3 mapcal[map]_sigmu_pars_vs_days_counter[i].txt
  ve/copy p3(1) smu([i])
  ve/copy p3(2) rmu([i])
  ve/copy p3(3)  nx([i])
enddo
ve/cre icnt(9) r 1 2 3 4 5 6 7 8 9
ve/cre dicnt(9) r
set pmci 1
sigma dsmu = rmu/sqrt(nx)
* exec $PER/s#vpl smu rmu icnt dicnt sz=0.1
atitle 'counter' '[s]?[m]!/[m]'
ve/cre p1(1) r 0.25
ve/fit icnt smu rmu p0 s 1 p1
txt=[s]?[m]!/[m] = $sigma(int(p1(1)*1000+0.5)/1000)
exec $PER/s#tf 0.05 0.9 [txt]
exec save mapcal[map]_sigmu_vs_counters.eps f
return



macro nevtsum dir=v2
ve/del runs,days,beams
ve/read runs,days,beams [dir]/runpars.txt 3f15.6
n=$vlen(runs)
ve/cre runf([n]) r
ve/cre run1([n]) r
ve/cre run2([n]) r
ve/cre day1([n]) r
ve/cre day2([n]) r
ve/cre beam([n]) r
ve/cre nevt([n]) r
ve/cre nevtf([n]) r
if ([dir].eq.'v1') then
  suffh=.his
  sufft=.txt
else
  suffh=_[dir].his
  sufft=_[dir].txt
endif
if ([dir].eq.'vs') then
  suffh=_v2.his
  sufft=_v2.txt
endif
ind=0
ebo=0
do i=1,[n]
  run=runs([i])
  fname=[dir]/run_[run]_spects[suffh]
  if ($fexist([fname]).eq.1) then
    hi/file 20 [fname]
    ebi=beams([i])
    if ([ebi].ne.[ebo]) then
      if ([ebo].ne.0) then
        ind=[ind]+1
        ve/inp beam([ind]) [ebo]
        ve/inp run1([ind]) [r1]
        ve/inp run2([ind]) [r2]
        ve/inp runf([ind]) [rf]
        ve/inp nevt([ind]) [nevts]
        ve/inp nevtf([ind]) [nmax]
        mess [ind]: [ebo] [r1] [r2] [nevts]
      endif
      r1=[run]
      nevts=0
      nmax=0
      rf=0
      ebo=[ebi]
    endif
    nxf=0
    do nc=1,9
      idh=30+[nc]
      hi/del 100
      hi/copy //lun20/[idh] 100
      nx=$hinfo(100,'events')
      nxf=[nxf]+[nx]
    enddo
    nevts=[nevts]+[nxf]
    if ([nxf].gt.[nmax]) then
      nmax=[nxf]
      rf=[run]
    endif
    r2=[run]
    close 20
  endif
enddo
ind=[ind]+1
ve/inp beam([ind]) [ebo]
ve/inp run1([ind]) [r1]
ve/inp run2([ind]) [r2]
ve/inp nevt([ind]) [nevts]
ve/inp runf([ind]) [rf]
ve/inp nevtf([ind]) [nmax]
mess [ind]: [ebo] [r1] [r2] [nevts]
exec mapcal#vecut run1
exec mapcal#vecut run2
exec mapcal#vecut runf $vlen(run1)
exec mapcal#vecut beam
exec mapcal#vecut nevt $vlen(run1)
exec mapcal#vecut nevtf $vlen(run1)
ve/write run1,run2,runf,beam,nevt,nevtf [dir]/nevt_selected_vs_beams.txt '(6f15.6)'
return

macro nevtsim
dir=v2
ve/del run1,run2,runf,beam,nevt,nevtf
ve/read run1,run2,runf,beam,nevt,nevtf [dir]/nevt_selected_vs_beams.txt '(6f15.6)'
ve/cre ind1(4) r 1 45 89 137
ve/cre ind2(4) r 44 89 135 145
ve/read beami,frun flist.fwi.txt
sigma beami = order(beami,frun)
sigma frun  = order(frun,frun)
fn1=map1_n1.13.fwi
fn2=map2_n1.13.fwi
fn3=map3_n1.05.fwi
fn4=map4_n1.13.fwi
k=2000000/390140
do i=1,4
  i1=ind1([i])
  i2=ind2([i])
  ve/del beamr,nevtr,run1r,run2r,runfr
  ve/copy beam([i1]:[i2]) beamr
  ve/copy nevt([i1]:[i2]) nevtr
  ve/copy run1([i1]:[i2]) run1r
  ve/copy run2([i1]:[i2]) run2r
  ve/copy runf([i1]:[i2]) runfr
  sigma nevtx = int(2000*nevtr/vsum(nevtr)+0.5)*1000
  sigma qevt = (10000*nevtx/nevtr+0.5)/100
*  nx=[i2]-[i1]+1
*  ve/cre frunx([nx]) r
*  ve/cre beamx([nx]) r
*  do j=1,$vlen(beamr)
*    exec mapcal#ixndl $sigma(run1r([j])) frun
*    gl/imp inx
*    inx1=[inx]
*    exec mapcal#ixndr $sigma(run2r([j])) frun
*    gl/imp inx
*    inx2=[inx]
*    ve/inp frunx([j]) frun([inx])
*    ve/inp beamx([j]) beami([inx])
*    mess $sigma(run1r([j])) $sigma(run2r([j])) $sigma(frun([inx1])) $sigma(frun([inx2]))
*  enddo
*  sigma qbeam=beamx/beamr
  ve/write beamr,nevtx,nevtr,qevt,runfr ! '(5f10.2)'
  if ($fexist([fn[i]]).eq.1) then
    shell rm [fn[i]]
  endif
  for/file 20 [fn[i]] n
  close 20
  fmess 'flist=[\' [fn[i]]
  do j=1,$vlen(beamr)
    en=beamr([j])
    fr=runfr([j])
    nv=nevtx([j])
    txt1=$format([en],'f6.1')
    txt2=$format([fr],'i6')
    txt3=$format([nv],'i6')
    txt=($unquote([txt1]),$unquote([txt2]),$unquote([txt3])),
    fmess [txt] [fn[i]]
  enddo
  fmess ']' [fn[i]]
  read x
enddo
return

macro convcsvx
do i=1,9
  exec mapcal#convcsv ag[i]hv.csv
enddo
return


macro convcsv fname=ag1hv.csv year0=2011
*shell cp [fname] tmp1.txt 
*shell $unquote('cat tmp1.txt | sed "s/\"from_unixtime(stamp\/1e6)\",\"value\"/ /g" > tmp2.txt')
*shell $unquote('cat tmp2.txt | sed "s/\"/ /g" > tmp3.txt')
*shell $unquote('cat tmp3.txt | sed "s/:/ /g" > tmp4.txt')
*shell $unquote('cat tmp4.txt | sed "s/-/ /g" > tmp5.txt')
*ve/del year,month,day,hour,minute,second,hv
*ve/read year,month,day,hour,minute,second,hv tmp5.txt
*fnameo=[fname].txt
*if ($fexist([fnameo])) then
*  shell rm [fnameo]
*endif
*for/file 20 [fnameo] n
*close 20
*do i=1,$vlen(year)
*  shell /work/users/konctbel/MinuitTest/gtime [year0].01.01.00.00.00 $sigma(year([i])).$sigma(month([i])).$sigma(day([i])).$sigma(hour([i])).$sigma(minute([i])).$sigma(second([i])) gtime.txt
*  ve/read gtime gtime.txt
*  txt=$sigma(gtime(1)) $sigma(hv([i]))
*  fmess [txt] [fnameo]
*enddo
shell ./cvstxt [fname] [fname].txt
return


macro qetest nc=1 i1=1 i2=145
ve/del runsp,drunsp,daysp,ddaysp,beamp,dbeamp,nevtp,meanp,dmeanp,rmsp,drmsp,wp,dwp,vp,dvp
ve/read runsp,drunsp,daysp,ddaysp,beamp,dbeamp,nevtp,meanp,dmeanp,rmsp,drmsp,wp,dwp,vp,dvp v2/dphi_vs_beams_new.txt 15e15.6
dir=v2
if ([dir].eq.'v1') then
  suffh=.his
  sufft=.txt
else
  suffh=_[dir].his
  sufft=_[dir].txt
endif
n=[i2]-[i1]+1
ve/cre ai([n]) r
ve/cre dai([n]) r
ve/cre ti([n]) r
ve/cre dti([n]) r
ve/cre asi([n]) r
ve/cre dasi([n]) r
ve/cre dp8(8) r
do i=[i1],[i2]
  idhc=100
  hi/del [idhc]
  r1=$sigma(runsp([i])-drunsp([i]))
  r2=$sigma(runsp([i])+drunsp([i]))
  do ri=[r1],[r2]
    fname=[dir]/run_[ri]_spects[suffh]
    if ($fexist([fname])) then
      hi/file 20 [fname]
      idh=40+[nc]
      if ($hexist([idhc])) then
*        hi/op/add [idhc] //lun20/[idh] [idhc]
        hi/get/cont [idh] vdhci
      else
        hi/copy //lun20/[idh] [idhc]
        nx=$hinfo([idhc],'xbins')
        xmin=$hinfo([idhc],'xmin')
        xmax=$hinfo([idhc],'xmax')
        ve/cre vdhci([nx]) r
        hi/get/cont [idhc] vdhci
        ve/cre vdhc([nx]) r
      endif
      sigma vdhc = vdhc + vdhci
      close 20
    endif
  enddo
  hi/put/cont [idhc] vdhc
  exec ../SepPar/sp#brn 100 10 100
  ve/cre p4(4) r 20 35 12 0.2
  ve/cre p8(8) r 20 31 12 0.25266 150 3 2.5 0.25266
  ve/cre s4(4) r 1 1 0 0
  ve/cre s8(8) r 1 1 0 0 1 1 0 0
*  hi/fit 200(15.:) ../SepPar/fun3.for b 4 p4 s4 ! ! dp4
*  hi/fit 200(15.:) ../SepPar/fun3.for ! 4 p4 ! ! ! dp4
  hi/fit 200 ../SepPar/fun32.for b 8 p8 s8 ! ! dp8
  ve/cre s8(8) r 1 1 1 0 1 1 0.1 0.01
  hi/fit 200 ../SepPar/fun32.for b 8 p8 s8 ! ! dp8
  ind=[i]-[i1]+1
  ve/inp  ai([ind]) $sigma(p8(2))
  ve/inp dai([ind]) $sigma(dp8(2))
  ve/inp  asi([ind]) $sigma(p8(4))
  ve/inp dasi([ind]) $sigma(dp8(4))
  ve/inp  ti([ind]) $sigma(daysp([i]))
  ve/inp dti([ind]) $sigma(ddaysp([i]))
*  read x
enddo
ve/cre ki([n]) r [n]*1
ve/cre kix(5) r 5*10
ve/copy kix(1:5) ki(13:17)
sigma daix = dai*ki
* exec $PER/s#vpl ai dai ti dti ll=1
ve/cre p3(3) r 23 5 1000
ve/fit ti ai daix expp0.f s 3 p3
return


macro accledsave
ve/del r1,r2
*ve/read r1,r2 accled_corr.ixlist
ve/read r1,r2 accled_ixlist_actcals.txt
v=3
exec mapcal#readpar r1 a [v]
exec mapcal#readpar r1 b [v]
exec mapcal#readpar r1 c [v]
exec mapcal#readpar r1 amp1pe [v]
exec mapcal#readpar r1 eff1pe [v]
exec mapcal#readpar r1 pds [v]
exec mapcal#readpar r1 s0 [v]
exec mapcal#readpar r1 u0 [v]
exec mapcal#readpar r1 tmin [v]
exec mapcal#readpar r1 tmax [v]
return



macro a1corr
ve/del r,e,t1,t2,v1,v2,v3,v4,v5,v6,v7,v8,v9,k1,k2,k3,k4,k5,k6,k7,k8,k9
ve/read r,e,t1,t2,v1,v2,v3,v4,v5,v6,v7,v8,v9,k1,k2,k3,k4,k5,k6,k7,k8,k9 /work/users/konctbel/AGPARS/run_hv.txt
*exec mapcal#readpar r a
*exec mapcal#readpar r b
*exec mapcal#readpar r c
*exec mapcal#readpar r amp1pe
p=a
ve/del v[p]1,v[p]2,v[p]3,v[p]4,v[p]5,v[p]6,v[p]7,v[p]8,v[p]9
ve/read v[p]1,v[p]2,v[p]3,v[p]4,v[p]5,v[p]6,v[p]7,v[p]8,v[p]9 accled_[p].txt 9f10.3
p=b
ve/del v[p]1,v[p]2,v[p]3,v[p]4,v[p]5,v[p]6,v[p]7,v[p]8,v[p]9
ve/read v[p]1,v[p]2,v[p]3,v[p]4,v[p]5,v[p]6,v[p]7,v[p]8,v[p]9 accled_[p].txt 9f10.3
p=c
ve/del v[p]1,v[p]2,v[p]3,v[p]4,v[p]5,v[p]6,v[p]7,v[p]8,v[p]9
ve/read v[p]1,v[p]2,v[p]3,v[p]4,v[p]5,v[p]6,v[p]7,v[p]8,v[p]9 accled_[p].txt 9f10.3
p=ae
ve/del v[p]1,v[p]2,v[p]3,v[p]4,v[p]5,v[p]6,v[p]7,v[p]8,v[p]9
ve/read v[p]1,v[p]2,v[p]3,v[p]4,v[p]5,v[p]6,v[p]7,v[p]8,v[p]9 accled_amp1pe.txt 9f10.3
*
sigma u1 = va1*0
sigma u1 = -3000*((vae1-vb1)/va1)**(1/vc1)
sigma du = u1*0
null 690 790 -2900 -2800
set pmci 2
* exec $PER/s#vpl v1 du t1 du sz=0.01 ll=-1 o=s
set pmci 4
* exec $PER/s#vpl u1 du t1 du sz=0.01 ll=-1 o=s
return


macro accledread
ve/del r1,r2
ve/read r1,r2 accled.ixlist
exec mapcal#readpar r1 a 1
exec mapcal#readpar r1 b 1
exec mapcal#readpar r1 c 1
exec mapcal#readpar r1 amp1pe 1
p=a
ve/del v[p]1,v[p]2,v[p]3,v[p]4,v[p]5,v[p]6,v[p]7,v[p]8,v[p]9
ve/read v[p]1,v[p]2,v[p]3,v[p]4,v[p]5,v[p]6,v[p]7,v[p]8,v[p]9 accled_[p]_v1.txt 9f10.3
p=b
ve/del v[p]1,v[p]2,v[p]3,v[p]4,v[p]5,v[p]6,v[p]7,v[p]8,v[p]9
ve/read v[p]1,v[p]2,v[p]3,v[p]4,v[p]5,v[p]6,v[p]7,v[p]8,v[p]9 accled_[p]_v1.txt 9f10.3
p=c
ve/del v[p]1,v[p]2,v[p]3,v[p]4,v[p]5,v[p]6,v[p]7,v[p]8,v[p]9
ve/read v[p]1,v[p]2,v[p]3,v[p]4,v[p]5,v[p]6,v[p]7,v[p]8,v[p]9 accled_[p]_v1.txt 9f10.3
p=ae
ve/del v[p]1,v[p]2,v[p]3,v[p]4,v[p]5,v[p]6,v[p]7,v[p]8,v[p]9
ve/read v[p]1,v[p]2,v[p]3,v[p]4,v[p]5,v[p]6,v[p]7,v[p]8,v[p]9 accled_amp1pe_v1.txt 9f10.3
ve/write r1,r2,v[p]1,v[p]2,v[p]3,v[p]4,v[p]5,v[p]6,v[p]7,v[p]8,v[p]9 accled_amp1pe_db.txt 11f15.6
return


macro runtest
ve/del runs,nevt,mean,rms,w,dw,v,dv
ve/read runs,nevt,mean,rms,w,dw,v,dv dphi_vs_runs.txt 8f15.6
n=$vlen(runs)
ve/cre runs0([n]) r
ind=0
fo=tmp.txt
if ($fexist([fo])) then
  shell rm [fo]
endif
for/file 20 [fo] n
close 20
fmess '-->' [fo]
do i=1,[n]
  r=runs([i])
  fname=v2/run_[r]_spects_v2.his
  if ($fexist([fname]).eq.0) then
    ind=[ind]+1
    ve/inp runs0([ind]) [r]
    shell fgrep -e [r] /work/users/konctbel/snd2k/R005-999/fwk/*/*.fwi > tmp1.txt
    shell cat tmp.txt tmp1.txt > tmp2.txt
    shell mv tmp2.txt tmp.txt
  endif
enddo
exec mapcal#vecut runs0
return


macro ampstep p=4
ve/cre st(10) r
ve/cre dst(10) r
ve/cre ni(10) r 0 1 2 3 4 5 6 7 8 9
ve/cre dni(10) r
do i=0,9
  ve/read p4,dp4 mapcal3_sigmu2_fit_counter[i].txt '(2f15.6)'
  sigma p4=abs(p4)
  ind=[i]+1
  ve/copy p4([p]) st([ind])
  ve/copy dp4([p]) dst([ind])
enddo
* exec $PER/s#vpl st dst ni dni
return


macro hvtestx p1=7 p2=21 per=r
exec seteps 0
opt ngrid
do p=[p1],[p2]
  do nc=1,9
    exec mapcal#hvtest [p] hv [nc] [per]
    exec mapcal#hvtest [p] i [nc] [per]
    exec mapcal#hvtest [p] f [nc] [per]
    exec mapcal#hvtest [p] r [nc] [per]
  enddo
enddo
return


macro hvtest part=10 p=hv nc=1 per=r
dir=/work/users/konctbel/AGPARS
fname=[dir]/run_part[part].txt_[p].txt
if ([p].ne.'r') then
  shell cp [fname] tmp.txt
  shell $unquote('cat tmp.txt | sed "s/nan/0/g" > tmp1.txt')
  shell $unquote('cat tmp1.txt | sed "s/inf/0/g" > tmp.txt')
*  shell cp tmp.txt [fname]
  fname=tmp.txt
endif
if ([p].eq.'hv') then
  ve/del run,beam,ts,tf,[p]1,[p]2,[p]3,[p]4,[p]5,[p]6,[p]7,[p]8,[p]9,k[p]1,k[p]2,k[p]3,k[p]4,k[p]5,k[p]6,k[p]7,k[p]8,k[p]9
  ve/read run,beam,ts,tf,[p]1,[p]2,[p]3,[p]4,[p]5,[p]6,[p]7,[p]8,[p]9,k[p]1,k[p]2,k[p]3,k[p]4,k[p]5,k[p]6,k[p]7,k[p]8,k[p]9 [fname]
  s=+
  tl=U?[nc]! , V
  E=U?[nc]!
  M=V
  k=10
*
  ve/del r1,r2
  ve/read r1,r2 accled_corr.ixlist
  q=a
  ve/del v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9
  ve/read v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9 accled_[q]_v1.txt 9f10.3
  q=b
  ve/del v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9
  ve/read v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9 accled_[q]_v1.txt 9f10.3
  q=c
  ve/del v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9
  ve/read v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9 accled_[q]_v1.txt 9f10.3
  q=ae
  ve/del v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9
  ve/read v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9 accled_amp1pe_v1.txt 9f10.3
*
  n=$vlen(run)
  ve/cre vpx0([n]) r
  ind=1
  do i=1,[n]
    ri=run([i])
    while (([ri].gt.$sigma(r2([ind]))).and.([ind].lt.$vlen(r2))) do
      ind=[ind]+1
    endwhile
    a=va[nc]([ind])
    b=vb[nc]([ind])
    c=vc[nc]([ind])
    ae=vae[nc]([ind])
    if (([a].ne.0).and.([b].ne.0).and.([c].ne.0).and.([ae].ne.0)) then
      ve/inp vpx0([i]) $sigma(-3000*(([ae]-([b]))/([a]))**(1/([c])))
    endif
  enddo
*
  ve/del r1db,r2db,u1db,u2db,u3db,u4db,u5db,u6db,u7db,u8db,u9db
  ve/read r1db,r2db,u1db,u2db,u3db,u4db,u5db,u6db,u7db,u8db,u9db hv_cals_corr_set.txt
  ve/cre vpy0([n]) r
  nr=$vlen(r2db)
  ve/inp r2db([nr]) 1000000
  ind=1
  do i=1,[n]
    ri=run([i])
    while (([ri].gt.$sigma(r2db([ind]))).and.([ind].lt.$vlen(r2db))) do
      ind=[ind]+1
    endwhile
    ve/inp vpy0([i]) $sigma(-u[nc]db([ind]))
  enddo
*  
  ve/del tic
  ve/read tic calstime.txt
  ve/del ric
  ve/read ric calsrun.txt
*
endif
if ([p].eq.'i') then
  ve/del run,beam,ts,tf,[p]1,[p]2,[p]3,[p]4,[p]5,[p]6,[p]7,[p]8,[p]9,k[p]1,k[p]2,k[p]3,k[p]4,k[p]5,k[p]6,k[p]7,k[p]8,k[p]9
  ve/read run,beam,ts,tf,[p]1,[p]2,[p]3,[p]4,[p]5,[p]6,[p]7,[p]8,[p]9,k[p]1,k[p]2,k[p]3,k[p]4,k[p]5,k[p]6,k[p]7,k[p]8,k[p]9 [fname]
  s=-
  tl=I?[nc]! , [m]A
  E=I?[nc]!
  M=mkA
  k=10
endif
if ([p].eq.'f') then
  ve/del run,beam,dt,ts,tf,[p]1,[p]2,[p]3,[p]4,[p]5,[p]6,[p]7,[p]8,[p]9
  ve/read run,beam,dt,ts,tf,[p]1,[p]2,[p]3,[p]4,[p]5,[p]6,[p]7,[p]8,[p]9 [fname]
  s=-
  tl=f?[nc]! , Hz
  E=f?[nc]!
  M=Hz
  k=100
endif
if ([p].eq.'r') then
  sigma r[nc] = abs(hv[nc]/i[nc])
  s=-
  tl=R?[nc]! , MOm
  E=R?[nc]!
  M=MOm
  k=10
endif
ve/del vx,vp
if ([per].eq.r) then
  sigma vx = run
  perx=run
  tx='runs'
else
  sigma vx = (ts+tf)/2
  perx=time
  tx='time, days'
endif
n=$vlen(vx)
ve/cre ix([n]) r 
sigma ix = array([n],1#[n])
sigma vx = order(vx,[s][p][nc])
sigma vpx = order(vpx0,[s][p][nc])
sigma vpy = order(vpy0,[s][p][nc])
sigma ix = order(ix,[s][p][nc])
sigma vp = order([p][nc],[s][p][nc])
exec mapcal#vecut vp
n=$vlen(vp)
if ([n].ne.0) then
  ve/cre dn([n]) r
  exec mapcal#vecut vx [n]
  exec mapcal#vecut vpx [n]
  exec mapcal#vecut vpy [n]
  exec mapcal#vecut ix [n]
  sigma vx = order(vx,ix)
  sigma vpx = order(vpx,ix)
  sigma vpy = order(vpy,ix)
  sigma vp = order(vp,ix)
  if ([p].eq.'hv') then
    l0=$sigma(vmin(vx))
    r0=$sigma(vmax(vx))
    l=$sigma([l0]-0.05*([r0]-[l0]))
    r=$sigma([r0]+0.05*([r0]-[l0]))
    d0=$sigma(min(vmin(vp),vmin(vpx)))
    u0=$sigma(max(vmax(vp),vmax(vpx)))
    if ([u0].eq.0) then
      u0=$sigma(vmax(vp))
    endif
    d=$sigma([d0]-0.1*([u0]-[d0]))
    u=$sigma([u0]+0.1*([u0]-[d0]))
    null [l] [r] [d] [u]
    set pmci 1
    * exec $PER/s#vpl vp dn vx dn sz=0.05 ll=-1 o=s
    set pmci 2
    * exec $PER/s#vpl vpx dn vx dn sz=0.1 ll=-1 o=s
    set pmci 4
    * exec $PER/s#vpl vpy dn vx dn sz=0.05 ll=-1 o=s
    set pmci 1
  else
    set pmci 1
    * exec $PER/s#vpl vp dn vx dn sz=0.05 ll=-1
  endif
  atitle [tx] [tl]
*
  mean0=$sigma(vsum(vp)/[n])
  ve/del vt
  sigma vt = vp-([mean0])
  mean=$sigma(vsum(vt)/[n])
  mean2=$sigma(vsum(vt**2)/[n])
  rms=$sigma(sqrt([mean2]-[mean]**2))
  mean=[mean]+([mean0])
  mess [mean] [rms] [mean2]
  meanr=$sigma(int([mean]*[k]+0.5)/[k]) 
  sig=$sigma(int([rms]*[k]+0.5)/[k])
  txt=[E] = [meanr] [\261] [sig] [M]
  exec $PER/s#tf 0.05 0.9 [txt]
*
*
  if ([p].eq.'hv') then
    mean0=$sigma(vsum(vpx)/[n])
    ve/del vt
    sigma vt = vpx-([mean0])
    mean=$sigma(vsum(vt)/[n])
    mean2=$sigma(vsum(vt**2)/[n])
    rms=$sigma(sqrt([mean2]-[mean]**2))
    mean=[mean]+([mean0])
    mess [mean] [rms] [mean2]
    meanr=$sigma(int([mean]*[k]+0.5)/[k]) 
    sig=$sigma(int([rms]*[k]+0.5)/[k])
    txt=[E] = [meanr] [M] (cal)
    exec $PER/s#tf 0.55 0.9 [txt]
*    
    l0=$GRAFINFO('WNXMIN')
    r0=$GRAFINFO('WNXMAX')
    d0=$GRAFINFO('WNYMIN')
    u0=$GRAFINFO('WNYMAX')
    if ([per].eq.r) then
      vp=ric
    else
      vp=tic
    endif
    set hcol 3
    set dmod 2
    set basl 0.01
    set ltyp 12
    set lwid 1
    do i=1,$vlen([vp])
      pi=[vp]([i])
      if (([pi].ge.[l0]).and.([pi].lt.[r0])) then
        line [pi] [d0] [pi] [u0]
      endif
    enddo
    set hcol 1
    set dmod 1
    set ltyp 1
    set pmci 4
*    
  endif
*
else
  null 0 1 0 1
endif
exec save [dir]/[p]_vs_[perx]_part[part]_counter[nc].eps f 
return


macro hvsteps
ve/del r1,r2
ve/read r1,r2 accled_corr.ixlist
q=a
ve/del v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9
ve/read v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9 accled_[q]_v1.txt 9f10.3
q=b
ve/del v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9
ve/read v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9 accled_[q]_v1.txt 9f10.3
q=c
ve/del v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9
ve/read v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9 accled_[q]_v1.txt 9f10.3
q=ae
ve/del v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9
ve/read v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9 accled_amp1pe_v1.txt 9f10.3
*
n=$vlen(r1)
n=[n]-1
do nc=1,9
  ve/cre hvst[nc]([n]) r
enddo
do i=1,[n]
  do nc=1,9
    ind=[i]
    hv1=$sigma(-3000*((vae[nc]([ind])-vb[nc]([ind]))/va[nc]([ind]))**(1/vc[nc]([ind])))
    ind=[i]+1
    hv2=$sigma(-3000*((vae[nc]([ind])-vb[nc]([ind]))/va[nc]([ind]))**(1/vc[nc]([ind])))
    ve/inp hvst[nc]([i]) $sigma([hv2]-[hv1])
  enddo
enddo
return


macro hvscans sc=1 nc=1 u0=3000
ve/del r1,r2
ve/read r1,r2 accled_corr.ixlist
q=a
ve/del v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9
ve/read v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9 accled_[q]_v1.txt 9f10.3
q=b
ve/del v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9
ve/read v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9 accled_[q]_v1.txt 9f10.3
q=c
ve/del v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9
ve/read v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9 accled_[q]_v1.txt 9f10.3
q=ae
ve/del v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9
ve/read v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9 accled_amp1pe_v1.txt 9f10.3
*
read x
ns=0
ns=[ns]+1; sc[ns]=MHAD2010; dir[ns]=null
ns=[ns]+1; sc[ns]=MHAD2011; dir[ns]=/online/gridspool/proc/MHAD2011-4/
ns=[ns]+1; sc[ns]=MHAD2012; dir[ns]=/online/gridspool/proc/MHAD2012-2/
ns=[ns]+1; sc[ns]=OMEG2012; dir[ns]=/online/gridspool/proc/OMEG2012-0/
ns=[ns]+1; sc[ns]=PHI_2013; dir[ns]=/online/gridspool/proc/PHI_2013-0/
ns=[ns]+1; sc[ns]=RHO_2012; dir[ns]=/online/gridspool/proc/RHO_2012-0/
*
ns=0
ns=[ns]+1; stp[ns]=cn113; ns[ns]=5000;  nf[ns]=7841;
ns=[ns]+1; stp[ns]=cn113; ns[ns]=7842;  nf[ns]=13847;
ns=[ns]+1; stp[ns]=cn105; ns[ns]=13847; nf[ns]=16698;
ns=[ns]+1; stp[ns]=cn113; ns[ns]=16699; nf[ns]=17868;
ns=[ns]+1; stp[ns]=cn105; ns[ns]=17869; nf[ns]=21000;
scx=4
*
fst=1
n=$vlen(r1)
ve/cre xae[nc]([n]) r
ve/cre va[nc]r([n]) r
ve/cre vb[nc]r([n]) r
ve/cre vc[nc]r([n]) r
ve/cre ti([n]) r
ind=0
p=hv
dir=/work/users/konctbel/AGPARS
fnameo=figvam
do i=1,[n]
  ri=r1([i])
  if ((([ri].ge.[ns[sc]]).and.([ri].le.[nf[sc]])).or.(([ri].ge.[ns[scx]]).and.([ri].le.[nf[scx]]))) then
    part=$sigma(int([ri]/1000))
    fname=[dir]/run_part[part].txt_[p].txt
    if ([fname].ne.[fnameo]) then
      fnameo=[fname]
      shell cp [fname] tmp.txt
      shell $unquote('cat tmp.txt | sed "s/nan/0/g" > tmp1.txt')
      shell $unquote('cat tmp1.txt | sed "s/inf/0/g" > tmp.txt')
      fname=tmp.txt
      ve/del run,beam,ts,tf,[p]1,[p]2,[p]3,[p]4,[p]5,[p]6,[p]7,[p]8,[p]9,k[p]1,k[p]2,k[p]3,k[p]4,k[p]5,k[p]6,k[p]7,k[p]8,k[p]9
      ve/read run,beam,ts,tf,[p]1,[p]2,[p]3,[p]4,[p]5,[p]6,[p]7,[p]8,[p]9,k[p]1,k[p]2,k[p]3,k[p]4,k[p]5,k[p]6,k[p]7,k[p]8,k[p]9 [fname]
    endif
    ig=$sigma([ri]-1000*[part])
*    
    a=$sigma(va[nc]([i]))
    b=$sigma(vb[nc]([i]))
    c=$sigma(vc[nc]([i]))
    if ([c].ne.0) then
      ind=[ind]+1
      ve/inp va[nc]r([ind]) [a]
      ve/inp vb[nc]r([ind]) [b]
      ve/inp vc[nc]r([ind]) [c]
      ve/inp ti([ind]) $sigma(ts([ig]))
    endif
    a1=$sigma(vae[nc]([i]))
    hv=$sigma(3000*((([a1])-([b]))/([a]))**(1/([c])))
*    c=$sigma(vcm([nc]))
    fu0=$sigma(([a])*([u0]/3000)**([c])+([b]))
    if ([fu0].ne.0) then
      a=[a]*100/[fu0]
      b=[b]*100/[fu0]
    endif
*    mess [a] [b] [c] [a1] [hv]
    if ([fst].eq.1) then
      fun/pl ([a])*(x/3000)**([c])+([b]) $sigma([hv]-100) $sigma([hv]+100)
      fst=0
    else
      fun/pl ([a])*(x/3000)**([c])+([b]) $sigma([hv]-100) $sigma([hv]+100) s
    endif
    key [hv] [a1] 20 ! 0.1
    ve/inp xae[nc]([i]) [fu0]
  endif
enddo
read x
exec mapcal#vecut va[nc]r [ind]
exec mapcal#vecut vb[nc]r [ind]
exec mapcal#vecut vc[nc]r [ind]
exec mapcal#vecut ti [ind]
ve/cre dy([ind]) r [ind]*10
ve/cre xi(1) r
*
npar=3
ve/cre chi2(2) r
ve/cre paru([npar]) r
ve/cre dparu([npar]) r
ve/cre covu([npar],[npar]) r
ve/cre daysx(6) r 0 250 250 550 642 850
do ii=1,3
  j=[ii]*2-1
  tx=daysx([j])
  exec mapcal#ixndl [tx] ti
  gl/imp inx
  in1=[inx]
  j=[ii]*2
  tx=daysx([j])
  exec mapcal#ixndr [tx] ti
  gl/imp inx
  in2=[inx]
  mess [in1] [in2]
*
  ve/cre p3(3) r 150 -2 -100
  sig=100000
  do k=1,2
    ve/del vy,dvy,vx,dvx
    ve/copy va[nc]r([in1]:[in2]) vy
    ve/copy dy([in1]:[in2]) dvy
    ve/copy ti([in1]:[in2]) vx
    ve/copy dy([in1]:[in2]) dvx
    id=0
    do i=1,$vlen(vx)
      ve/inp xi(1) $sigma(vx([i]))
      vyi=$call('expp0p.f(xi)')
      if ($sigma(abs(vy([i])-([vyi]))/[sig]).lt.3) then
        id=[id]+1
        ve/copy  vy([i])  vy([id])
        ve/copy dvy([i]) dvy([id])
        ve/copy  vx([i])  vx([id])
        ve/copy dvx([i]) dvx([id])
      endif
    enddo
    exec mapcal#vecut  vy [id]
    exec mapcal#vecut dvy [id]
    exec mapcal#vecut  vx [id]
    exec mapcal#vecut dvx [id]
*    ve/cre p3(3) r 150 -2 -100
    * exec $PER/s#vpl vy dvy vx dvx
    ve/fit vx vy dvy expp0.f s 3 p3
    call covm.f(1)
    call covmpen(chi2,[npar],paru,dparu)
    call covmcov([npar],covu)
    sig=$sigma(dparu(1)*sqrt(chi2(1)))
    sig=100
*    read x
  enddo
  ve/copy  vy  vy[ii]
  ve/copy dvy dvy[ii]
  ve/copy  vx  vx[ii]
  ve/copy dvx dvx[ii]
  ve/copy p3 p3x[ii]
enddo
*
* exec $PER/s#vpl va[nc]r dy ti dy sz=0.1 iatt=24
do ii=1,3
  * exec $PER/s#vpl vy[ii] dvy[ii] vx[ii] dvx[ii] sz=0.1 iarr=20 o=s
  ve/fit vx[ii] vy[ii] dvy[ii] expp0.f s 3 p3x[ii]
enddo
atitle 'a?i!' 'time, days'
exec save setup[ns]_a_vs_days_counter[nc].eps f
read x
*
ve/cre dy([ind]) r [ind]*1
npar=1
ve/cre chi2(2) r
ve/cre paru([npar]) r
ve/cre dparu([npar]) r
ve/cre covu([npar],[npar]) r
sig=1000
*
do k=1,2
  ve/del vy,dvy,vx,dvx
  ve/copy vb[nc]r vy
  ve/copy dy dvy
  ve/copy ti vx
  ve/copy dy dvx
  id=0
  do i=1,$vlen(vx)
    if ($sigma(abs(vy([i])-paru(1))/[sig]).lt.3) then
      id=[id]+1
      ve/copy  vy([i])  vy([id])
      ve/copy dvy([i]) dvy([id])
      ve/copy  vx([i])  vx([id])
      ve/copy dvx([i]) dvx([id])
    endif
  enddo
  exec mapcal#vecut  vy [id]
  exec mapcal#vecut dvy [id]
  exec mapcal#vecut  vx [id]
  exec mapcal#vecut dvx [id]
  ve/fit vx vy dvy p0 s
  call covm.f(1)
  call covmpen(chi2,[npar],paru,dparu)
  call covmcov([npar],covu)
  sig=$sigma(dparu(1)*sqrt(chi2(1)))
  sig=1
  mess $sigma(paru(1)) [sig]
enddo
*
* exec $PER/s#vpl vb[nc]r dy ti dy sz=0.1 ll=-1 iatt=24
* exec $PER/s#vpl vy dvy vx dvx sz=0.1 ll=-1 iatt=20 o=s
ve/fit vx vy dvy p0 s
*
atitle 'b?i!' 'time, days'
exec save setup[ns]_b_vs_days_counter[nc].eps f
read x
*
npar=1
ve/cre chi2(2) r
ve/cre paru([npar]) r
ve/cre dparu([npar]) r
ve/cre covu([npar],[npar]) r
sig=1000
*
do k=1,2
  ve/del vy,dvy,vx,dvx
  ve/copy vc[nc]r vy
  ve/copy dy dvy
  ve/copy ti vx
  ve/copy dy dvx
  id=0
  do i=1,$vlen(vx)
    if ($sigma(abs(vy([i])-paru(1))/[sig]).lt.3) then
      id=[id]+1
      ve/copy  vy([i])  vy([id])
      ve/copy dvy([i]) dvy([id])
      ve/copy  vx([i])  vx([id])
      ve/copy dvx([i]) dvx([id])
    endif
  enddo
  exec mapcal#vecut  vy [id]
  exec mapcal#vecut dvy [id]
  exec mapcal#vecut  vx [id]
  exec mapcal#vecut dvx [id]
  ve/fit vx vy dvy p0 s
  call covm.f(1)
  call covmpen(chi2,[npar],paru,dparu)
  call covmcov([npar],covu)
  sig=$sigma(dparu(1)*sqrt(chi2(1)))
  sig=1
  mess $sigma(paru(1)) [sig]
enddo
*
* exec $PER/s#vpl vc[nc]r dy ti dy sz=0.1 iatt=24
* exec $PER/s#vpl vy dvy vx dvx sz=0.1 iatt=20 o=s
ve/fit vx vy dvy p0 s
*
atitle 'c?i!' 'time, days'
exec save setup[ns]_c_vs_days_counter[nc].eps f
read x
*
return


macro hvcals  i1=1 i2=10000
ve/del r1,r2
*ve/read r1,r2 accled_corr.ixlist
ve/read r1,r2 accled_ixlist_actcals.txt
dir1 = /work/users/konctbel/MinuitTest/
v=3
q=a
ve/del v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9
ve/read v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9 [dir1]accled_[q]_v[v].txt 9f10.3
q=b
ve/del v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9
ve/read v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9 [dir1]accled_[q]_v[v].txt 9f10.3
q=c
ve/del v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9
ve/read v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9 [dir1]accled_[q]_v[v].txt 9f10.3
q=ae
ve/del v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9
ve/read v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9 [dir1]accled_amp1pe_v[v].txt 9f10.3
*
p=hv
dir=/work/users/konctbel/AGPARS
n=$vlen(r1)
ve/del runc,ac1,ac2,ac3,ac4,ac5,ac6,ac7,ac8,ac9
do i=$sigma(max(1,[i1])),$sigma(min([n],[i2]))
  rb=r1([i])
  re=r2([i])
  if ([re].eq.-1) then
    re=25000
  endif
  partmin=$sigma(int([rb]/1000)) 
  partmax=$sigma(int([re]/1000))
*  
  ve/del runz,beamz,tsz,tfz
  do j=1,9
    ve/del [p][j]z,k[p][j]z
  enddo
  do part=[partmin],[partmax]
    fname=[dir]/run_part[part].txt_[p].txt
    if ($fexist([fname])) then
      shell cp [fname] tmp.txt
      shell $unquote('cat tmp.txt | sed "s/nan/0/g" > tmp1.txt')
      shell $unquote('cat tmp1.txt | sed "s/inf/0/g" > tmp.txt')
      fname=tmp.txt
      ve/del run,beam,ts,tf,[p]1,[p]2,[p]3,[p]4,[p]5,[p]6,[p]7,[p]8,[p]9,k[p]1,k[p]2,k[p]3,k[p]4,k[p]5,k[p]6,k[p]7,k[p]8,k[p]9
      ve/read run,beam,ts,tf,[p]1,[p]2,[p]3,[p]4,[p]5,[p]6,[p]7,[p]8,[p]9,k[p]1,k[p]2,k[p]3,k[p]4,k[p]5,k[p]6,k[p]7,k[p]8,k[p]9 [fname]
*
      exec mapcal#ixndl [rb] run
      gl/imp inx
      in1=[inx]
      exec mapcal#ixndr [re] run
      gl/imp inx
      in2=[inx]
      mess [part],[in1],[in2]
*  
      vpar=run;  ve/del [vpar]p; ve/copy [vpar]([in1]:[in2]) [vpar]p; exec vappend [vpar]z [vpar]p
      vpar=beam; ve/del [vpar]p; ve/copy [vpar]([in1]:[in2]) [vpar]p; exec vappend [vpar]z [vpar]p
      vpar=ts;   ve/del [vpar]p; ve/copy [vpar]([in1]:[in2]) [vpar]p; exec vappend [vpar]z [vpar]p
      vpar=tf;   ve/del [vpar]p; ve/copy [vpar]([in1]:[in2]) [vpar]p; exec vappend [vpar]z [vpar]p
      do j=1,9
        vpar=[p][j];   ve/del [vpar]p; ve/copy [vpar]([in1]:[in2]) [vpar]p; exec vappend [vpar]z [vpar]p
        vpar=k[p][j];  ve/del [vpar]p; ve/copy [vpar]([in1]:[in2]) [vpar]p; exec vappend [vpar]z [vpar]p
      enddo
    endif
  enddo
*
  do nc=1,9
    ae=vae[nc]([i])
    a=va[nc]([i])
    b=vb[nc]([i])
    c=vc[nc]([i])
    ve/del aei,dn,ti
    vhv=[p][nc]z
    if ($vlen([vhv]).ne.0) then
      sigma  aei = ([a])*(-[p][nc]z/3000)**([c])+([b])
    else
      npi=$vlen(runz)
      ve/cre aei([npi]) r [npi]*[ae]
    endif
    sigma raei = aei/[ae]
    sigma dn = aei*0
    sigma ti = (tsz+tfz)/2
*   
    if ([a].ne.0) then
    ni=$vlen(runz)
    ve/cre ix([ni]) r
    sigma ix = array([ni],1#[ni])
    ve/del runzs,dns,raeis,ixs
    sigma runzs = order(runz,raei) 
    sigma dns   = order(dn,raei) 
    sigma ixs   = order(ix,raei) 
    sigma raeis = order(raei,raei)
    raeimin=0.7
    raeimax=1.3
    exec mapcal#ixndl [raeimin] raeis
    gl/imp inx
    in1=[inx]
    exec mapcal#ixndr [raeimax] raeis
    gl/imp inx
    in2=[inx]
    mess [in1] [in2] [re] [rb]
    if ([in1].lt.[in2]) then
      ve/del runzsr,dnsr,raeisr,ixsr
      ve/copy runzs([in1]:[in2]) runzsr
      ve/copy dns([in1]:[in2]) dnsr
      ve/copy raeis([in1]:[in2]) raeisr
      ve/copy ixs([in1]:[in2]) ixsr
      ve/del runzs,dns,raeis,ixs
      sigma runzs = order(runzsr,ixsr)
      sigma dns   = order(dnsr,ixsr)
      sigma raeis = order(raeisr,ixsr)
      sigma ixs   = order(ixsr,ixsr)
*    
      set pmci 1
      * exec $PER/s#vpl raeis dns runzs dns sz=0.1 ll=-1
      u0=$GRAFINFO('WNYMAX')
      d0=$GRAFINFO('WNYMIN')
      d=$sigma(min(1-0.05*([u0]-[d0]),[d0]))
      u=$sigma(max(1+0.20*([u0]-[d0]),[u0]))
      null $GRAFINFO('WNXMIN') $GRAFINFO('WNXMAX') [d] [u]
      * exec $PER/s#vpl raeis dns runzs dns sz=0.1 ll=-1 o=s
      set plci 2
      line $GRAFINFO('WNXMIN') 1 $GRAFINFO('WNXMAX') 1
      set plci 1
      txt=calibr [i]: ([rb],[re])
      exec $PER/s#tf 0.05 0.9 [txt]
      txt=Counter [nc]
      exec $PER/s#tf 0.65 0.9 [txt]
      atitle 'N, runs' 'A?1!(U?real!)/A?1!(U?set!)'
      exec save [dir]/AmpCorrHV/amplitude_correction_hv_cal[i]_counter[nc].eps f
*
    endif
    if ($vlen(runz).ne.$vlen(raei)) then
      mess $vlen(runz) $vlen(raei) [i] [nc]
      read x
    endif
    endif
    sigma aci = 1/raei
    sigma uni = [p][nc]z/[p][nc]z
    sigma acix = (aci-1)*uni+1
*    ve/write aci,acix ! 2f15.6
    exec vappend ac[nc] acix
*    read x
  enddo
  exec vappend runc runz
*  read x
enddo
nf=$vlen(runc)
rb=runc(1)
re=runc([nf])
n=[re]-[rb]+1
ve/cre runcf([n]) r
sigma runcf=array([n],[rb]#[re])
do i=1,9
  ve/cre ac[i]f([n]) r [n]*1
enddo
do j=1,[nf]
  ind=$sigma(runc([j])-[rb]+1)
  ve/inp runcf([ind]) $sigma(runc([j]))
  do i=1,9 
    ve/inp ac[i]f([ind]) $sigma(ac[i]([j]))
  enddo
enddo
ve/write runcf,ac1f,ac2f,ac3f,ac4f,ac5f,ac6f,ac7f,ac8f,ac9f [dir]/AmpCorrHV/amplitude_correction_hv.txt '(10f15.6)'
return


macro hvcalstex
ve/del r1,r2
ve/read r1,r2 accled_ixlist_actcals.txt
n=$vlen(r1)
*
dir=/work/users/konctbel/AGPARS/AmpCorrHV
*
fname=[dir]/AmpCorrHV0.tex
if ($fexist([fname])) then
  shell rm [fname]
endif
for/file 20 [fname] ! N
close 20
*
do i=1,[n]
  do nc=1,9
    epsfile=[dir]/amplitude_correction_hv_cal[i]_counter[nc].eps
    if ($fexist([epsfile]).eq.0) then
      epsfile=/work/users/konctbel/Calibr/null.eps
    endif
    fmess '\begin{figure}[ht!b]' [fname]
    fmess '  \begin{minipage}{0.7\textwidth}' [fname]
    fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
    txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
    fmess [txt] [fname]
    txt=$unquote('     ')\caption{cal=[i], counter=[nc]}
*    fmess [txt] [fname]
    fmess '   \end{minipage}' [fname]
    fmess '\end{figure}' [fname]
    if ($sigma(mod([nc],9)).eq.0) then
      fmess '\clearpage' [fname]
    endif
  enddo
enddo
return


macro hvcalsth
dir=/work/users/konctbel/AGPARS
ve/del runcf,ac1f,ac2f,ac3f,ac4f,ac5f,ac6f,ac7f,ac8f,ac9f
ve/read runcf,ac1f,ac2f,ac3f,ac4f,ac5f,ac6f,ac7f,ac8f,ac9f [dir]/AmpCorrHV/amplitude_correction_hv.txt '(10f15.6)'
imin=$sigma(int((vmin(runcf)-999)/1000))
imax=$sigma(int((vmax(runcf)+999)/1000))
n=[imax]-[imin]
do i=[imin],[imax]
  rb=1000*[i]
  re=[rb]+999
  exec mapcal#ixndl [rb] runcf
  gl/imp inx
  in1=[inx]
  exec mapcal#ixndr [re] runcf
  gl/imp inx
  in2=[inx]
  if ([in1].lt.[in2]) then
    do nc=1,9
      ve/del runcfr,ac[nc]fr
      mess [part],[in1],[in2]
      ve/copy runcf([in1]:[in2]) runcfr
      ve/copy ac[nc]f([in1]:[in2]) ac[nc]fr
      ni=$vlen(runcfr)
      ve/cre ix([ni]) r
      sigma ix = array([ni],1#[ni])
      ve/del runzs,dns,raeis,ixs
      sigma runcfrs   = order(runcfr,-ac[nc]fr) 
      sigma ac[nc]frs = order(ac[nc]fr,-ac[nc]fr) 
      sigma ixs       = order(ix,-ac[nc]fr)
      exec mapcal#vecut ac[nc]frs
      ni=$vlen(ac[nc]frs)
      exec mapcal#vecut runcfrs [ni]
      exec mapcal#vecut ixs [ni]
      sigma runcfrs = order(runcfrs,ixs)
      sigma ac[nc]frs = order(ac[nc]frs,ixs)
      np=$vlen(runcfrs)
      if ([np].ne.0) then
        ve/cre dns([np]) r
        set pmci 1
        * exec $PER/s#vpl ac[nc]frs dns runcfrs dns sz=0.05 ll=-1
        set plci 2
        line $GRAFINFO('WNXMIN') 1 $GRAFINFO('WNXMAX') 1
        set plci 1
        txt=part [i]
        exec $PER/s#tf 0.05 0.9 [txt]
        txt=Counter [nc]
        exec $PER/s#tf 0.65 0.9 [txt]
        atitle 'N, runs' 'A?1!(U?real!)/A?1!(U?set!)'
        exec save [dir]/AmpCorrHV/amplitude_correction_hv_part[i]_counter[nc].eps f
*        read x
      endif
    enddo
  endif
enddo
return


macro vecorr
dir=/work/users/konctbel/AGPARS
nf=$vlen(runc)
rb=runc(1)
re=runc([nf])
n=[re]-[rb]+1
ve/cre runcf([n]) r
sigma runcf=array([n],[rb]#[re])
do i=1,9
  ve/cre ac[i]f([n]) r [n]*1
enddo
do j=1,[nf]
  ind=$sigma(runc([j])-[rb]+1)
  ve/inp runcf([ind]) $sigma(runc([j]))
  do i=1,9 
    ve/inp ac[i]f([ind]) $sigma(ac[i]([j]))
  enddo
enddo
ve/write runcf,ac1f,ac2f,ac3f,ac4f,ac5f,ac6f,ac7f,ac8f,ac9f [dir]/AmpCorrHV/amplitude_correction_hv.txt '(10f15.6)'
return

macro corrtest
npf=$vlen(ci)
ve/cre dni([npf]) r
exec vpl#pl0 ci dni bi dni sz=0.1 ll=-1
*
bim=$sigma(vsum(bi)/[npf])
cim=$sigma(vsum(ci)/[npf])
dxx=$sigma(vsum((bi-[bim])*(bi-[bim]))/[npf])
dxy=$sigma(vsum((bi-[bim])*(ci-[cim]))/[npf])
dyy=$sigma(vsum((ci-[cim])*(ci-[cim]))/[npf])
sbi=$sigma(sqrt([dxx]))
sci=$sigma(sqrt([dyy]))
ve/cre bir([npf]) r
ve/cre cir([npf]) r
ind=0
do i=1,[npf]
  bii=bi([i])
  cii=ci([i])
  if ((abs(([bim]-([bii]))/[sbi]).lt.2).and.(abs(([cim]-([cii]))/[sci]).lt.2)) then
    ind=[ind]+1
    ve/inp bir([ind]) [bii]
    ve/inp cir([ind]) [cii]
  endif
enddo
exec mapcal#vecut bir [ind]
exec mapcal#vecut cir [ind]
bim=$sigma(vsum(bir)/[ind])
cim=$sigma(vsum(cir)/[ind])
dxx=$sigma(vsum((bir-[bim])*(bir-[bim]))/[ind])
dxy=$sigma(vsum((bir-[bim])*(cir-[cim]))/[ind])
dyy=$sigma(vsum((cir-[cim])*(cir-[cim]))/[ind])
*
tg2phi=$sigma(2*([dxy])/(([dyy])-([dxx])))
phi=$sigma(atan([tg2phi]))
mess [phi] [tg2phi] [dxx] [dxy] [dyy]
if ([phi].lt.0) then
  phi=[phi]+3.1415927/2
endif
phi=[phi]/2
mess phi=$sigma([phi]*180/3.1415927) [phi] $sigma(tan([phi]))
line $sigma(([bim])-10) $sigma(([cim])-10*tan([phi])) $sigma(([bim])+10) $sigma(([cim])+10*tan([phi]))
*
dxx=$sigma(vsum(di4)/[npf])
dxy=$sigma(vsum(di5)/[npf])
dyy=$sigma(vsum(di6)/[npf])
tg2phi=$sigma(2*([dxy])/(([dyy])-([dxx])))
phi=$sigma(atan([tg2phi]))
mess [phi] [tg2phi] [dxx] [dxy] [dyy]
if ([phi].lt.0) then
  phi=[phi]+3.1415927/2
endif
phi=[phi]/2
mess phi=$sigma([phi]*180/3.1415927) [phi] $sigma(tan([phi]))
line $sigma(([bim])-10) $sigma(([cim])-10*tan([phi])) $sigma(([bim])+10) $sigma(([cim])+10*tan([phi]))
*
return


macro totcals
ve/cre vnopt(2) r 113 105

do mode=0,2

fname=/work/users/konctbel/Calibr/totcals0_m[mode].tex
if ($fexist([fname])) then
  shell rm [fname]
endif
for/file 20 [fname] ! N
close 20

do i=1,$vlen(vnopt)

nopt=vnopt([i])

do nc=1,9
  fkname=totcals_n[nopt]_counter[nc].kumac
  if ($fexist([fkname]).eq.1) then
    shell rm [fkname]
  endif
  for/file 20 [fkname] new
  close 20
  txt=exec mapcal#newcalnx [nc] 1 3 1 180 [nopt] 1 [mode]
  fmess [txt] [fkname]
  txt=exec mapcal#newcalnx [nc] 1 4 1 180 [nopt] 1 [mode]
  fmess [txt] [fkname]
  shell pawbigX11 -b [fkname]
*  
  fmess '\begin{figure}[ht!b]' [fname]
  fmess '  \begin{minipage}{\textwidth}' [fname]
  fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
  epsfile=/work/users/konctbel/MinuitTest/calfit/c_vs_a_correlation_mcal_m[mode]_n[nopt]_counter[nc].eps
  txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
  fmess [txt] [fname]
  txt=$unquote('     ')\caption{n=[nopt] counter [nc]}
  fmess [txt] [fname]
  fmess '   \end{minipage}' [fname]
  fmess '\end{figure}' [fname]
  
  fmess '\begin{figure}[ht!b]' [fname]
  fmess '  \begin{minipage}{\textwidth}' [fname]
  fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
  epsfile=/work/users/konctbel/MinuitTest/calfit/b_vs_a_correlation_mcal_m[mode]_n[nopt]_counter[nc].eps
  txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
  fmess [txt] [fname]
  txt=$unquote('     ')\caption{n=[nopt] counter [nc]}
  fmess [txt] [fname]
  fmess '   \end{minipage}' [fname]
  fmess '\end{figure}' [fname]
  
  fmess '\clearpage' [fname]
  
  fmess '\begin{figure}[ht!b]' [fname]
  fmess '  \begin{minipage}{\textwidth}' [fname]
  fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
  epsfile=/work/users/konctbel/MinuitTest/calfit/c_vs_b_correlation_mcal_m[mode]_n[nopt]_counter[nc].eps
  txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
  fmess [txt] [fname]
  txt=$unquote('     ')\caption{n=[nopt] counter [nc]}
  fmess [txt] [fname]
  fmess '   \end{minipage}' [fname]
  fmess '\end{figure}' [fname]
  
  fmess '\begin{figure}[ht!b]' [fname]
  fmess '  \begin{minipage}{\textwidth}' [fname]
  fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
  epsfile=/work/users/konctbel/MinuitTest/calfit/mu_vs_days_correlation_mcal_m[mode]_n[nopt]_counter[nc].eps
  txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
  fmess [txt] [fname]
  txt=$unquote('     ')\caption{n=[nopt] counter [nc]}
  fmess [txt] [fname]
  fmess '   \end{minipage}' [fname]
  fmess '\end{figure}' [fname]
  
  fmess '\clearpage' [fname]
  
  fmess '\begin{figure}[ht!b]' [fname]
  fmess '  \begin{minipage}{\textwidth}' [fname]
  fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
  epsfile=/work/users/konctbel/MinuitTest/calfit/a_vs_days_correlation_0_mcal_m[mode]_n[nopt]_counter[nc].eps
  txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
  fmess [txt] [fname]
  txt=$unquote('     ')\caption{n=[nopt] counter [nc]}
  fmess [txt] [fname]
  fmess '   \end{minipage}' [fname]
  fmess '\end{figure}' [fname]
  
  fmess '\begin{figure}[ht!b]' [fname]
  fmess '  \begin{minipage}{\textwidth}' [fname]
  fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
  epsfile=/work/users/konctbel/MinuitTest/calfit/a_vs_days_correlation_mcal_m[mode]_n[nopt]_counter[nc].eps
  txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
  fmess [txt] [fname]
  txt=$unquote('     ')\caption{n=[nopt] counter [nc]}
  fmess [txt] [fname]
  fmess '   \end{minipage}' [fname]
  fmess '\end{figure}' [fname]
  
  fmess '\clearpage' [fname]
  
  fmess '\begin{figure}[ht!b]' [fname]
  fmess '  \begin{minipage}{\textwidth}' [fname]
  fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
  epsfile=/work/users/konctbel/MinuitTest/calfit/b_vs_days_correlation_0_mcal_m[mode]_n[nopt]_counter[nc].eps
  txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
  fmess [txt] [fname]
  txt=$unquote('     ')\caption{n=[nopt] counter [nc]}
  fmess [txt] [fname]
  fmess '   \end{minipage}' [fname]
  fmess '\end{figure}' [fname]
  
  fmess '\begin{figure}[ht!b]' [fname]
  fmess '  \begin{minipage}{\textwidth}' [fname]
  fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
  epsfile=/work/users/konctbel/MinuitTest/calfit/b_vs_days_correlation_mcal_m[mode]_n[nopt]_counter[nc].eps
  txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
  fmess [txt] [fname]
  txt=$unquote('     ')\caption{n=[nopt] counter [nc]}
  fmess [txt] [fname]
  fmess '   \end{minipage}' [fname]
  fmess '\end{figure}' [fname]
  
  fmess '\clearpage' [fname]
  
  fmess '\begin{figure}[ht!b]' [fname]
  fmess '  \begin{minipage}{\textwidth}' [fname]
  fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
  epsfile=/work/users/konctbel/MinuitTest/calfit/c_vs_days_correlation_0_mcal_m[mode]_n[nopt]_counter[nc].eps
  txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
  fmess [txt] [fname]
  txt=$unquote('     ')\caption{n=[nopt] counter [nc]}
  fmess [txt] [fname]
  fmess '   \end{minipage}' [fname]
  fmess '\end{figure}' [fname]
  fmess '\begin{figure}[ht!b]' [fname]
  
  fmess '  \begin{minipage}{\textwidth}' [fname]
  fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
  epsfile=/work/users/konctbel/MinuitTest/calfit/c_vs_days_correlation_mcal_m[mode]_n[nopt]_counter[nc].eps
  txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
  fmess [txt] [fname]
  txt=$unquote('     ')\caption{n=[nopt] counter [nc]}
  fmess [txt] [fname]
  fmess '   \end{minipage}' [fname]
  fmess '\end{figure}' [fname]
  
  fmess '\clearpage' [fname]
  
  fmess '\begin{figure}[ht!b]' [fname]
  fmess '  \begin{minipage}{\textwidth}' [fname]
  fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
  epsfile=/work/users/konctbel/MinuitTest/calfit/amp_vs_days_correlation_mcal_m[mode]_n[nopt]_counter[nc].eps
  txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
  fmess [txt] [fname]
  txt=$unquote('     ')\caption{n=[nopt] counter [nc]}
  fmess [txt] [fname]
  fmess '   \end{minipage}' [fname]
  fmess '\end{figure}' [fname]
  
  fmess '\begin{figure}[ht!b]' [fname]
  fmess '  \begin{minipage}{\textwidth}' [fname]
  fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
  epsfile=/work/users/konctbel/MinuitTest/calfit/amp_ratio_vs_days_correlation_mcal_m[mode]_n[nopt]_counter[nc].eps
  txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
  fmess [txt] [fname]
  txt=$unquote('     ')\caption{n=[nopt] counter [nc]}
  fmess [txt] [fname]
  fmess '   \end{minipage}' [fname]
  fmess '\end{figure}' [fname]
  
  fmess '\clearpage' [fname]
enddo
enddo
enddo
return



macro newcalnx nc=1 np=1 mode=1 n1m=-1 n2m=10000 nopt=100 auto=0 fmode=0
n=0
n=[n]+1; cal[n]=Cal_2010_03_16_0001; ncal[n]=0; n1=[n]
n=[n]+1; cal[n]=Cal_2010_03_15_0002; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2010_03_15_0003; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2010_03_15_0004; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2010_04_08_0005; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2010_04_13_0006; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2010_04_19_0007; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2010_04_20_0008; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2010_04_20_0009; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2010_04_20_0010; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2010_04_20_0011; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2010_04_20_0012; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2010_04_20_0013; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2010_04_20_0014; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2010_04_20_0015; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2010_04_21_0016; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2010_04_21_0017; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2010_04_22_0018; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2010_04_26_0019; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2010_04_29_0020; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2010_04_29_0021; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2010_04_29_0022; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2010_04_29_0023; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2010_04_29_0024; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2010_04_29_0025; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2010_05_06_0026; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2010_05_07_0027; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2010_05_11_0028; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2010_05_18_0029; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2010_05_18_0030; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2010_05_26_0031; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2010_05_26_0032; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2010_05_26_0033; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2010_06_07_0034; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2010_06_14_0035; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2010_12_06_0036; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2010_12_06_0037; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2010_12_07_0038; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2010_12_07_0039; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2010_12_07_0040; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_01_12_0041; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_01_21_0042; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_01_24_0043; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_02_09_0044; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_02_11_0045; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_02_13_0046; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_02_16_0047; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_02_21_0048; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_02_25_0049; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_02_28_0050; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_03_09_0051; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_03_14_0052; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_03_16_0053; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_03_21_0054; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_03_25_0055; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_03_28_0056; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_03_30_0057; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_04_04_0058; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_04_10_0059; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_04_14_0060; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_04_27_0061; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_04_27_0062; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_05_05_0063; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_05_05_0064; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_05_11_0065; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_05_14_0066; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_05_17_0067; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_05_20_0068; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_05_25_0069; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_05_25_0070; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_05_31_0071; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_06_02_0072; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_06_02_0073; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_06_07_0074; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_06_09_0075; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_06_15_0076; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_06_17_0077; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_10_21_0078; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_10_21_0079; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_10_22_0080; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_10_24_0081; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_10_25_0082; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_10_26_0083; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_10_27_0084; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_10_28_0085; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_11_02_0086; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_11_03_0087; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_11_03_0088; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_11_07_0089; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_11_21_0090; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_11_21_0091; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_11_25_0092; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_11_25_0093; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_11_29_0094; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2011_11_30_0095; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_12_05_0096; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_12_13_0097; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2011_12_19_0098; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_01_10_0099; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_01_20_0100; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_01_25_0101; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_01_25_0102; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_01_25_0103; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_02_01_0104; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_02_01_0105; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_02_03_0106; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_02_03_0107; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_02_07_0108; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_02_07_0109; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_02_13_0110; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_02_28_0111; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_03_05_0112; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_03_12_0113; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_03_18_0114; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_03_23_0115; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_03_27_0116; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_03_30_0117; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_04_04_0118; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_04_16_0119; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_04_23_0120; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_04_24_0121; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_04_26_0122; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_10_03_0123; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_10_10_0124; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_10_10_0125; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_10_15_0126; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_11_19_0127; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_11_20_0128; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_11_26_0129; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_11_29_0130; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_11_30_0131; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_11_30_0132; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_11_30_0133; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_11_30_0134; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_12_03_0135; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_12_10_0136; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_12_18_0137; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_12_18_0138; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2012_12_24_0139; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2013_01_02_0140; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2013_01_04_0141; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2013_01_17_0142; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2013_01_14_0143; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2013_01_17_0144; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2013_01_21_0145; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2013_01_27_0146; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2013_01_31_0147; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2013_02_01_0148; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2013_02_01_0149; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2013_02_06_0150; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2013_02_12_0151; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2013_02_15_0152; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2013_02_18_0153; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2013_02_21_0154; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2013_04_12_0155; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2013_02_25_0156; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2013_02_26_0157; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2013_03_04_0158; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2013_03_06_0159; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2013_03_06_0160; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2013_03_12_0161; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2013_03_18_0162; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2013_03_29_0163; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2013_04_01_0164; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2013_04_08_0165; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2013_04_15_0166; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2013_04_15_0167; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2013_04_24_0168; ncal[n]=0;
n=[n]+1; cal[n]=Cal_2013_04_25_0169; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2013_04_25_0170; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2013_04_29_0171; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2013_05_06_0172; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2013_05_13_0173; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2013_05_16_0174; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2013_05_23_0175; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2013_05_28_0176; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2013_06_06_0177; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2013_06_17_0178; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2013_06_23_0179; ncal[n]=0;      
n=[n]+1; cal[n]=Cal_2013_07_04_0180; ncal[n]=0;
n2=[n]
dir0=/work/users/konctbel/Calibr
nm=0
nm=[nm]+1; m[nm]=jan
nm=[nm]+1; m[nm]=feb
nm=[nm]+1; m[nm]=mar
nm=[nm]+1; m[nm]=apr
nm=[nm]+1; m[nm]=may
nm=[nm]+1; m[nm]=jun
nm=[nm]+1; m[nm]=jul
nm=[nm]+1; m[nm]=agu
nm=[nm]+1; m[nm]=sep
nm=[nm]+1; m[nm]=oct
nm=[nm]+1; m[nm]=nov
nm=[nm]+1; m[nm]=dec
ve/cre n113(4) r 40 129 156 165
ve/cre n105(4) r 131 155 166 180
ve/cre n100(2) r 1 [n]
if ([mode].eq.'0a') then
  do ni=[n2],[n1],-1
    exec mapcal#newcaln [nc] [nc] [cal[ni]] 2
  enddo
endif
if ([mode].eq.0) then
  do ni=[n1],[n2]
    fname=cal[ni]_counter[nc].kumac
    if ($fexist([fname]).eq.1) then
      shell rm [fname]
    endif
    for/file  20 [fname] new
    close 20
    txt=exec mapcal#newcaln [nc] [nc] [cal[ni]] 2
    fmess [txt] [fname]
    shell .testrelease/.mainrelease/Offline/submit.sh -q clusters,180 pawbigX11 -b [fname]
  enddo
endif
if ([mode].eq.'0l') then
  do ni=$sigma(max([n1],[n1m])),$sigma(min([n2],[n2m]))
    fname=cal[ni]_counter[nc].kumac
    if ($fexist([fname]).eq.1) then
      shell rm [fname]
    endif
    for/file  20 [fname] new
    close 20
    txt=exec mapcal#newcaln [nc] [nc] [cal[ni]] 2
    fmess [txt] [fname]
    shell pawbigX11 -b [fname]
  enddo
endif
if ([mode].eq.'0lp') then
  do ni=$sigma(max([n1],[n1m])),$sigma(min([n2],[n2m]))
    fname=cal[ni]_counter[nc].kumac
    if ($fexist([fname]).eq.1) then
      shell rm [fname]
    endif
    for/file  20 [fname] new
    close 20
    txt=exec mapcal#newcalnp [nc] [nc] [cal[ni]] 2
    fmess [txt] [fname]
    shell pawbigX11 -b [fname]
  enddo
endif
if ([mode].eq.1) then
  ind=0
  do ni=[n1],[n2]
    date=$substring([cal[ni]],1,11)
    cal=$substring([cal[ni]],13,4)
    fname=[dir0]/[date]_[cal]/pars_correlation_counter[nc].txt
    if ($fexist([fname]).eq.1) then
      ve/del pars
      ve/read pars [fname]
      if ([ind].eq.0) then
        ve/cre parx(1000) r
      endif
      chi2i=pars(1)
      parxi=pars([np])
      if ([chi2i].lt.10) then
        ind=[ind]+1
        ve/inp parx([ind]) [parxi]
      else
        mess [cal[ni]]: [parxi]
      endif
    endif
  enddo
  exec mapcal#vecut parx
  ii=$vlen(parx)
  ve/cre ix([ii]) r
  sigma ix = array([ii],1#[ii])
  ve/cre dn([ii]) r
  exec vpl#pl0 parx dn ix dn sz=0.1 ll=-1
endif
if ([mode].eq.'1a') then
  fnamed=[dir0]/[date]_[cal]/fit_data_counter[nc].txt
  fnamem=[dir0]/[date]_[cal]/fit_mu_counter[nc].txt
  if (($fexist([fnamem]).eq.1).and.($fexist([fnamed]).eq.1)) then
    ve/del u0,amp1,damp1,mu,dmu,eff,deff
    ve/read u0,amp1,damp1,mu,dmu,eff,deff [fnamed] '(7f15.6)'
    ve/del parm,dparm
    ve/read parm,dparm [fnamem] '(2f15.6)'
  endif
endif
if ([mode].eq.2) then
  ind=0
  nf=[n2]-[n1]+1
  ve/cre ti0([n]) r
  do ni=[n1],[n2]
    date=$substring([cal[ni]],1,11)
    cal=$substring([cal[ni]],13,4)
    year=$substring([cal[ni]],1,4)
    monthx=$substring([cal[ni]],6,3)
    day=$substring([cal[ni]],10,2)
    do j=1,[nm]
      if ($lower([monthx]).eq.[m[j]]) then
        month=[j]
      endif
    enddo
    year0=2011
    txt=$format([year],i4).$format([month],i2.2).$format([day],i2.2).09.00.00
    shell /work/users/konctbel/MinuitTest/gtime [year0].01.01.00.00.00 [txt] gtime.txt
    ve/read gtime gtime.txt
    gl/cre ndays $sigma(gtime(1))
    mess [year] [monthx] [day] [ndays]
    ind=[ni]-[n1]+1
    ve/inp ti0([ind]) [ndays]
  enddo
  exec mapcal#vecut ti0
*  
  ve/del r1,r2
  ve/read r1,r2 accled_corr.ixlist
  nf=$vlen(r1)
  ve/cre ti([nf]) r
  p=hv
  dir=/work/users/konctbel/AGPARS
  fnameo=figvam
  ind=0
  do i=1,[nf]
    ri=r1([i])
    part=$sigma(int([ri]/1000))
    fname=[dir]/run_part[part].txt_[p].txt
    if ($fexist([fname])) then
      if ([fname].ne.[fnameo]) then
        fnameo=[fname]
        shell cp [fname] tmp.txt
        shell $unquote('cat tmp.txt | sed "s/nan/0/g" > tmp1.txt')
        shell $unquote('cat tmp1.txt | sed "s/inf/0/g" > tmp.txt')
        fname=tmp.txt
        ve/del run,beam,ts,tf,[p]1,[p]2,[p]3,[p]4,[p]5,[p]6,[p]7,[p]8,[p]9,k[p]1,k[p]2,k[p]3,k[p]4,k[p]5,k[p]6,k[p]7,k[p]8,k[p]9
        ve/read run,beam,ts,tf,[p]1,[p]2,[p]3,[p]4,[p]5,[p]6,[p]7,[p]8,[p]9,k[p]1,k[p]2,k[p]3,k[p]4,k[p]5,k[p]6,k[p]7,k[p]8,k[p]9 [fname]
      endif
      ig=$sigma([ri]-1000*[part])
      ind=[ind]+1
      ve/inp ti([ind]) $sigma(ts([ig]))
*    mess [ri] $sigma(ts([ig]))
    endif
  enddo
  exec mapcal#vecut ti
*
  nf=$vlen(ti)
  nf0=$vlen(ti0)
  ve/cre ix0([nf0]) r
  sigma ix0 = array([nf0],1#[nf0])
  ve/cre ti0c([nf]) r
  ve/cre ixc([nf]) r
  do i=1,[nf]
    tii=ti([i])
    sigma dti0=abs(ti0-([tii]))
    sigma ix0s = order(ix0,dti0)
    ix0i=ix0s(1)
    ve/inp ixc([i]) [ix0i]
    t0ii=ti0([ix0i])
    ve/inp ti0c([i]) [t0ii]
*    mess [tii] [t0ii] [ix0i]
  enddo
  ixcio=1
  do i=1,[nf]
    ixci=ixc([i])
    do j=[ixcio],$sigma([ixci]-1)
      ni=[j]+[n1]-1
      mess -- ------- $sigma(ti0([j])) [cal[ni]]
    enddo
    ni=[ixci]+[n1]-1
    mess [i] $sigma(ti([i])) $sigma(ti0([ixci])) [cal[ni]]
    ixcio=[ixci]+1
  enddo
endif
if ([mode].eq.3) then
  ind=0
  ve/cre idi(6) r 1 2 3 5 6 9
*  
  fname=calfit/calfit_sigs_m0_n[nopt]_counter[nc].txt
  if (($fexist([fname]).eq.1).and.([fmode].eq.1)) then
    ve/del pars
    ve/read pars [fname]
    amid=pars(1)
    bmid=pars(2)
    cmid=pars(3)
    asig=pars(4)
    bsig=pars(5)
    csig=pars(6)
  else
    amid=100
    bmid=0
    cmid=20
    asig=10000
    bsig=10000
    csig=10000
  endif
  al=$sigma([amid]-sqrt(3.0)*[asig])
  ar=$sigma([amid]+sqrt(3.0)*[asig])
  if ([fmode].eq.2) then
    nsig=0
  else
    nsig=1
  endif
*  
  if ([nopt].ne.100) then
    n1m=$sigma(vmin(n[nopt]))
    n2m=$sigma(vmax(n[nopt]))
  endif
  nit=$vlen(n[nopt])
  do ni=$sigma(max([n1],[n1m])),$sigma(min([n2],[n2m]))
    gni=0
    do is=1,$sigma([nit]/2)
      il=[is]*2-1
      ir=[is]*2
      nl=n[nopt]([il])
      nr=n[nopt]([ir])
      if (([ni].ge.[nl]).and.([ni].le.[nr])) then
        gni=1
      endif
    enddo
    if ([gni].eq.1) then
*    
    dir=[dir0]/[cal[ni]]
    cal=$substring([cal[ni]],16,4)
    fname=[dir]/pars_correlation_counter[nc].txt
    fnamed=[dir]/fit_data_counter[nc].txt
    fnamem=[dir]/fit_mu_counter[nc].txt
    if (($fexist([fname]).eq.1).and.($fexist([fnamed]).eq.1)) then
      ve/del pars
      shell cp [fname] tmp.txt
      shell $unquote('cat tmp.txt | sed "s/NaN/10000/g" > tmp1.txt')
      shell $unquote('cat tmp1.txt | sed "s/nan/10000/g" > tmp2.txt')
      shell $unquote('cat tmp2.txt | sed "s/inf/10000/g" > tmp.txt')
      fname=tmp.txt
      ve/read pars [fname]
      ve/del u0,amp1,damp1,amp,damp,mu,dmu,eff,deff
      shell cp [fnamed] tmp.txt
      shell $unquote('cat tmp.txt | sed "s/nan/10000/g" > tmp1.txt')
      shell $unquote('cat tmp1.txt | sed "s/inf/10000/g" > tmp.txt')
      fnamed=tmp.txt
      ve/read u0,amp1,damp1,amp,damp,mu,dmu,eff,deff [fnamed] '(9f15.6)'
      ve/del parm,dparm
      shell cp [fnamem] tmp.txt
      shell $unquote('cat tmp.txt | sed "s/nan/10000/g" > tmp1.txt')
      shell $unquote('cat tmp1.txt | sed "s/inf/10000/g" > tmp.txt')
      fnamem=tmp.txt
      ve/read parm,dparm [fnamem] '(2f15.6)'
      if ([ind].eq.0) then
        imax=1000
        ve/cre chi2i([imax]) r
        ve/cre ndfi([imax]) r
        ve/cre imi([imax]) r
        ve/cre ki([imax]) r
        ve/cre ai([imax]) r
        ve/cre bi([imax]) r
        ve/cre ci([imax]) r
        ve/cre mi([imax]) r
        ve/cre dmi([imax]) r
        do j=1,6
          ve/cre di[j]([imax]) r
        enddo
        npi=$vlen(u0)
        do j=1,[npi]
          ve/cre xi[j]([imax]) r
          ve/cre yi[j]([imax]) r
          ve/cre si[j]([imax]) r
        enddo
      endif
      chi2ii=pars(1)
      ndfii=pars(2)
      ve/cre uch(1) r $sigma([chi2ii]*[ndfii])
      ve/cre nch(1) i $sigma([ndfii])
      prob=$call('prob(uch,nch)')
      gcal=1
      if ([prob].lt.0.005) then
        gcal=0
      endif
      aii=pars(4)
      saii=[asig]
      if ($sigma(abs(([aii]-([amid]))/[saii])).gt.3) then
        gcal=0
      endif
*      if (([aii].lt.[al]).or.([aii].gt.[ar])) then
*        gcal=0
*      endif
      bii=pars(5)
      sbii=[bsig]
      if ($sigma(abs(([bii]-([bmid]))/[sbii])).gt.3) then
        gcal=0
      endif
      cii=pars(6)
      scii=[csig]
      if ($sigma(abs(([cii]-([cmid]))/[scii])).gt.3) then
        gcal=0
      endif
      if ([gcal].eq.1) then
        ind=[ind]+1
        ve/inp chi2i([ind]) [chi2ii]
        ve/inp ndfi([ind]) [ndfii]
        ve/inp imi([ind]) $sigma(pars(3))
        ve/inp ai([ind]) [aii]
        ve/inp bi([ind]) [bii]
        ve/inp ci([ind]) [cii]
        ve/inp mi([ind]) $sigma(parm(2))
        ve/inp dmi([ind]) $sigma(dparm(2))
        do j=1,6
          jj=$sigma(9+idi([j]))
          ve/inp di[j]([ind]) $sigma(pars([jj]))
        enddo
        do j=1,[npi]
          ve/inp xi[j]([ind]) $sigma(u0([j]))
          ve/inp yi[j]([ind]) $sigma(amp1([j]))
          ve/inp si[j]([ind]) $sigma(damp1([j]))
        enddo
*        
        told=0
        if ([told].eq.1) then
          year=$substring([cal[ni]],1,4)
          monthx=$substring([cal[ni]],6,3)
          day=$substring([cal[ni]],10,2)
          do j=1,[nm]
            if ($lower([monthx]).eq.[m[j]]) then
              month=[j]
            endif
          enddo
          hour=9
          minute=0
          second=0
        else
          fnamet=[dir]/time.txt
          ve/del t
          ve/read t [fnamet]
          year=t(1)
          month=t(2)
          day=t(3)
          hour=t(4)
          minute=t(5)
          second=t(6)
        endif
        txt=$format([year],i4).$format([month],i2.2).$format([day],i2.2).$format([hour],i2.2).$format([minute],i2.2).$format([second],i2.2)
        shell /work/users/konctbel/MinuitTest/gtime 2011.01.01.00.00.00 [txt] gtime.txt
        ve/read gtime gtime.txt
        ve/inp ki([ind]) $sigma(gtime(1))
      else
        mess [cal[ni]]: [parxi]
      endif
    endif
*    
    endif
  enddo
  exec mapcal#vecut chi2i [ind]
  exec mapcal#vecut ndfi [ind]
  exec mapcal#vecut imi [ind]
  exec mapcal#vecut ki [ind]
  exec mapcal#vecut ai [ind]
  exec mapcal#vecut bi [ind]
  exec mapcal#vecut ci [ind]
  exec mapcal#vecut mi [ind]
  exec mapcal#vecut dmi [ind]
  do j=1,6
    exec mapcal#vecut di[j] [ind]
  enddo
  ve/del xi,yi,si
  do j=1,[npi]
    exec mapcal#vecut xi[j] [ind]
    exec mapcal#vecut yi[j] [ind]
    exec mapcal#vecut si[j] [ind]
    exec vappend xi xi[j]
    exec vappend yi yi[j]
    exec vappend si si[j]
  enddo
  npf=$vlen(ai)

  fname=calfit
  dir=calfit

  file=[dir]/[fname].inc
  if ($fexist([file]).ne.0) then
    shell mv -v [file] [file].old
  endif
  for/file  20 [file] new
  close 20

  txt='      common/parfit/'
  txt=$unquote([txt])am,bm,cm,ncal,npnt,mi([npf]),
  fmess [txt] [file]
  txt='     &  '
  txt=$unquote([txt])ai([npf]),bi([npf]),ci([npf]),
  fmess [txt] [file]
  txt='     &  '
  txt=$unquote([txt])Daa([npf]),Dab([npf]),Dac([npf]),Dbb([npf]),Dbc([npf]),Dcc([npf]),
  fmess [txt] [file]
  txt='     &  '
  txt=$unquote([txt])xi([npf],[npi]),yi([npf],[npi]),si([npf],[npi])
  fmess [txt] [file]
  fmess '      double precision ai,bi,ci,Daa,Dab,Dac,Dbb,Dbc,Dcc' [file]
  fmess '      double precision xi,yi,si' [file]
  fmess '      double precision am,bm,cm' [file]
  fmess '      double precision mi' [file]
  fmess '      integer ncal,npnt' [file]
  fmess '' [file]
  txt='      data am,bm,cm,ncal,npnt'
  am=$sigma(vsum(ai)/[npf])
  bm=$sigma(vsum(bi)/[npf])
  cm=$sigma(vsum(ci)/[npf])
  txt=$unquote([txt])/[am],[bm],[cm],[npf],[npi]/
  fmess [txt] [file]
  fmess '' [file]
  fmess '      data ai /' [file]
  exec mapcal#vewrite [file] ai e13.6 [npf] 4
  fmess '     &  /' [file]
  fmess '' [file]
  fmess '      data bi /' [file]
  exec mapcal#vewrite [file] bi e13.6 [npf] 4
  fmess '     &  /' [file]
  fmess '' [file]
  fmess '      data ci /' [file]
  exec mapcal#vewrite [file] ci e13.6 [npf] 4
  fmess '     &  /' [file]
  fmess '' [file]
  fmess '      data Daa /' [file]
  exec mapcal#vewrite [file] di1 e13.6 [npf] 4
  fmess '     &  /' [file]
  fmess '' [file]
  fmess '      data Dab /' [file]
  exec mapcal#vewrite [file] di2 e13.6 [npf] 4
  fmess '     &  /' [file]
  fmess '' [file]
  fmess '      data Dac /' [file]
  exec mapcal#vewrite [file] di3 e13.6 [npf] 4
  fmess '     &  /' [file]
  fmess '' [file]
  fmess '      data Dbb /' [file]
  exec mapcal#vewrite [file] di4 e13.6 [npf] 4
  fmess '     &  /' [file]
  fmess '' [file]
  fmess '      data Dbc /' [file]
  exec mapcal#vewrite [file] di5 e13.6 [npf] 4
  fmess '     &  /' [file]
  fmess '' [file]
  fmess '      data Dcc /' [file]
  exec mapcal#vewrite [file] di6 e13.6 [npf] 4
  fmess '     &  /' [file]
  fmess '' [file]
  fmess '      data xi /' [file]
  exec mapcal#vewrite [file] xi e13.6 $sigma([npf]*[npi]) 4
  fmess '     &  /' [file]
  fmess '' [file]
  fmess '      data yi /' [file]
  exec mapcal#vewrite [file] yi e13.6 $sigma([npf]*[npi]) 4
  fmess '     &  /' [file]
  fmess '' [file]
  fmess '      data si /' [file]
  exec mapcal#vewrite [file] si e13.6 $sigma([npf]*[npi]) 4
  fmess '     &  /' [file]
  fmess '' [file]
  fmess '      data mi /' [file]
  exec mapcal#vewrite [file] imi f5.1 [npf] 6
  fmess '     &  /' [file]
  fmess '' [file]
  
  file=[dir]/[fname].dat
  if ($fexist([file]).eq.1) then
    shell mv -v [file] [file].old
  endif
  for/file  20 [file] new
  close 20

  dam=1*[nsig]
  dbm=1*[nsig]
  dcm=1*[nsig]
  pb=0; dpb=0.01*[nsig]; 
  qb=0; dqb=0.001*[nsig];
  pc=0; dpc=0.01*[nsig];
  qc=0; dqc=0.001*[nsig];
*  
  bim=$sigma(vsum(bi)/[npf])
  cim=$sigma(vsum(ci)/[npf])
  dxx=$sigma(vsum((bi-[bim])*(bi-[bim]))/[npf])
  dxy=$sigma(vsum((bi-[bim])*(ci-[cim]))/[npf])
  dyy=$sigma(vsum((ci-[cim])*(ci-[cim]))/[npf])
  sbi=$sigma(sqrt([dxx]))
  sci=$sigma(sqrt([dyy]))
  ve/cre bir([npf]) r
  ve/cre cir([npf]) r
  ind=0
  do i=1,[npf]
    bii=bi([i])
    cii=ci([i])
    if ((abs(([bim]-([bii]))/[sbi]).lt.3).and.(abs(([cim]-([cii]))/[sci]).lt.3)) then
      ind=[ind]+1
      ve/inp bir([ind]) [bii]
      ve/inp cir([ind]) [cii]
    endif
  enddo
  exec mapcal#vecut bir [ind]
  exec mapcal#vecut cir [ind]
  bim=$sigma(vsum(bir)/[ind])
  cim=$sigma(vsum(cir)/[ind])
  dxx=$sigma(vsum((bir-[bim])*(bir-[bim]))/[ind])
  dxy=$sigma(vsum((bir-[bim])*(cir-[cim]))/[ind])
  dyy=$sigma(vsum((cir-[cim])*(cir-[cim]))/[ind])
*  
*  dxx=$sigma(vsum(di4)/[npf])
*  dxy=$sigma(vsum(di5)/[npf])
*  dyy=$sigma(vsum(di6)/[npf])
  tg2phi=$sigma(2*([dxy])/(([dxx])-([dyy])))
  phi=$sigma(atan([tg2phi])/2)
  mess [phi] [tg2phi] [dxx] [dxy] [dyy]
  if ([tg2phi].lt.0) then
    phi=[phi]+3.1415927/2
  endif
  mess phi=$sigma([phi]*180/3.1415927) [phi] $sigma(tan([phi]))
  pc=[phi]; dpc=0.000*[nsig];
  qc=0; dqc=0.000*[nsig];
*
  if ([fmode].eq.2) then
    fout=calfit/calfit_m1_n[nopt]_counter[nc].out
    shell /work/snd2000/users/konctbel/SNDSimACC/readlifstat [fout]
    ve/del pars0,dpars0
    ve/read pars0,dpars0 [fout].pars
    am=pars0(2)
    bm=pars0(3)
    cm=pars0(4)
    pb=pars0(5)
    qb=pars0(6)
    pc=pars0(7)
    qc=pars0(8)
  endif

  fmess 'SET TITLE' [file]
  fmess 'calibration fit' [file]
  fmess 'PARAMETERS' [file]
  ni=0
*  
  pname=am
  ls=[[pname]]-10
  rs=[[pname]]+10
  ni=[ni]+1;
  tmp1=$format([ni],i-4); 
  tmp2=$format([[pname]],e13.5); 
  tmp3=$format([d[pname]],e13.5); 
  tmp4=$format([ls],e13.5); 
  tmp5=$format([rs],e13.5)
  tmp=$unquote([tmp1])$quote([pname])      $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
  fmess [tmp] [file]
*  
  pname=bm
  ls=[[pname]]-10
  rs=[[pname]]+10
  ni=[ni]+1;
  tmp1=$format([ni],i-4); 
  tmp2=$format([[pname]],e13.5); 
  tmp3=$format([d[pname]],e13.5); 
  tmp4=$format([ls],e13.5); 
  tmp5=$format([rs],e13.5)
  tmp=$unquote([tmp1])$quote([pname])      $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
  fmess [tmp] [file]
*  
  pname=cm
  ls=[[pname]]-10
  rs=[[pname]]+10
  ni=[ni]+1;
  tmp1=$format([ni],i-4); 
  tmp2=$format([[pname]],e13.5); 
  tmp3=$format([d[pname]],e13.5); 
  tmp4=$format([ls],e13.5); 
  tmp5=$format([rs],e13.5)
  tmp=$unquote([tmp1])$quote([pname])      $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
  fmess [tmp] [file]
*  
  ni0=[ni]
  pname=pb
  ls=[[pname]]-1
  rs=[[pname]]+1
  ni=[ni]+1;
  tmp1=$format([ni],i-4); 
  tmp2=$format([[pname]],e13.5); 
  tmp3=$format([d[pname]],e13.5); 
  tmp4=$format([ls],e13.5); 
  tmp5=$format([rs],e13.5)
  tmp=$unquote([tmp1])$quote([pname])      $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
  fmess [tmp] [file]
*  
  pname=qb
  ls=[[pname]]-1
  rs=[[pname]]+1
  ni=[ni]+1;
  tmp1=$format([ni],i-4); 
  tmp2=$format([[pname]],e13.5); 
  tmp3=$format([d[pname]],e13.5); 
  tmp4=$format([ls],e13.5); 
  tmp5=$format([rs],e13.5)
  tmp=$unquote([tmp1])$quote([pname])      $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
  fmess [tmp] [file]
*  
  pname=pc
  ls=[[pname]]-10
  rs=[[pname]]+10
  ni=[ni]+1;
  tmp1=$format([ni],i-4); 
  tmp2=$format([[pname]],e13.5); 
  tmp3=$format([d[pname]],e13.5); 
  tmp4=$format([ls],e13.5); 
  tmp5=$format([rs],e13.5)
  tmp=$unquote([tmp1])$quote([pname])      $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
  fmess [tmp] [file]
*  
  pname=qc
  ls=[[pname]]-1
  rs=[[pname]]+1
  ni=[ni]+1;
  tmp1=$format([ni],i-4); 
  tmp2=$format([[pname]],e13.5); 
  tmp3=$format([d[pname]],e13.5); 
  tmp4=$format([ls],e13.5); 
  tmp5=$format([rs],e13.5)
  tmp=$unquote([tmp1])$quote([pname])      $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
  fmess [tmp] [file]
*  
  nif=[ni]
  do i=1,[npf]
    tmp0=t_[i]
    ti=$sigma(ai([i])-([am]))
    dti=$sigma(sqrt(di1([i])))
*    dti=1
    til=$sigma([ti]-3*[dti])
    tir=$sigma([ti]+3*[dti])
    ni=[ni]+1; 
    pname=t_$format([i],i3.3); 
    tmp1=$format([ni],i-4); 
    tmp2=$format([ti],e13.5); 
    tmp3=$format([dti],e13.5); 
    tmp4=$format([til],e13.5); 
    tmp5=$format([tir],e13.5);
    tmp=$unquote([tmp1])$quote([pname])   $unquote([tmp2])$unquote([tmp3])$unquote([tmp4])$unquote([tmp5])
    fmess [tmp] [file]
  enddo
*    
  fmess '' [file]
  fmess '' [file]
*  fmess 'set err 0.5' [file]
*  fmess 'set eps 1e-9' [file]
*  fmess 'set strategy 2' [file]
*  fmess '' [file]

  do i=1,[npf]
    txt=fix $sigma([i]+[nif])
    fmess [txt] [file]
  enddo
  fmess 'mini' [file]
*  fmess 'ret' [file]
  do i=1,3
    txt=fix $sigma([i])
    fmess [txt] [file]
  enddo
  do i=1,4
    txt=fix $sigma([i]+[ni0])
    fmess [txt] [file]
  enddo
  do i=1,[npf]
    txt=release $sigma([i]+[nif])
    fmess [txt] [file]
  enddo
  fmess 'mini' [file]
  do i=1,-4
    txt=release $sigma([i]+[ni0])
    fmess [txt] [file]
  enddo
  do i=1,[npf]
    txt=fix $sigma([i]+[nif])
    fmess [txt] [file]
  enddo
  fmess 'mini' [file]
  do i=1,[npf]
    txt=release $sigma([i]+[nif])
    fmess [txt] [file]
  enddo
  fmess 'mini' [file]
  do i=1,[npf]
    txt=fix $sigma([i]+[nif])
    fmess [txt] [file]
  enddo
  fmess 'mini' [file]
  fmess 'mini' [file]
  fmess '' [file]
  fmess 'ret' [file]
*
  shell ./make_calfit.sh
*  
  shell cp calfit/calfit.inc calfit/calfit_m[fmode]_n[nopt]_counter[nc].inc
  shell cp calfit/calfit.out calfit/calfit_m[fmode]_n[nopt]_counter[nc].out
*
  npf=$vlen(ci)
  aim=$sigma(vsum(ai)/[npf])
  bim=$sigma(vsum(bi)/[npf])
  cim=$sigma(vsum(ci)/[npf])
  daa=$sigma(vsum((ai-([aim]))*(ai-([aim])))/[npf])
  dbb=$sigma(vsum((bi-([bim]))*(bi-([bim])))/[npf])
  dcc=$sigma(vsum((ci-([cim]))*(ci-([cim])))/[npf])
  sa=$sigma(sqrt([daa]))
  sb=$sigma(sqrt([dbb]))
  sc=$sigma(sqrt([dcc]))
  mess [aim] [bim] [cim] [sa] [sb] [sc]
* 
  ve/del dai,dbi,dci
  sigma dai=sqrt(di1)
  sigma dbi=sqrt(di4)
  sigma dci=sqrt(di6)
  sa=$sigma(vsum(1/di1))
  aim=$sigma(vsum(ai/di1)/[sa])
  sa=$sigma(1/sqrt([sa]))
  sb=$sigma(vsum(1/di4))
  bim=$sigma(vsum(bi/di4)/[sb])
  sb=$sigma(1/sqrt([sb]))
  sc=$sigma(vsum(1/di6))
  cim=$sigma(vsum(ci/di6)/[sc])
  sc=$sigma(1/sqrt([sc]))
  mess [npf] [aim] [bim] [cim] [sa] [sb] [sc]
*  
  al=$sigma([aim]-sqrt(3.0*[npf])*[sa])
  ar=$sigma([aim]+sqrt(3.0*[npf])*[sa])
  ve/cre air([npf]) r
  ve/cre bir([npf]) r
  ve/cre cir([npf]) r
  ve/cre di1r([npf]) r
  ve/cre di4r([npf]) r
  ve/cre di6r([npf]) r
  ind=0
  snpf=$sigma(sqrt([npf]))
  do i=1,[npf]
    aii=ai([i]); di1i=di1([i]); saii=$sigma(sqrt([di1i]))
    bii=bi([i]); di4i=di4([i]); sbii=$sigma(sqrt([di4i]))
    cii=ci([i]); di6i=di6([i]); scii=$sigma(sqrt([di6i]))
    nsa=$sigma(abs(([aii]-([aim]))/[saii]))
    nsb=$sigma(abs(([bii]-([bim]))/[sbii]))
    nsc=$sigma(abs(([cii]-([cim]))/[scii]))
*    if (([nsb].lt.3).and.([nsc].lt.3).and.([aii].ge.[al]).and.([aii].le.[ar])) then
    if (([nsb].lt.3).and.([nsc].lt.3).and.([nsa].lt.3)) then
      ind=[ind]+1
      ve/inp air([ind]) [aii]
      ve/inp bir([ind]) [bii]
      ve/inp cir([ind]) [cii]
      ve/inp di1r([ind]) [di1i]
      ve/inp di4r([ind]) [di4i]
      ve/inp di6r([ind]) [di6i]
    endif
  enddo
  exec mapcal#vecut air [ind]
  exec mapcal#vecut bir [ind]
  exec mapcal#vecut cir [ind]
  exec mapcal#vecut di1r [ind]
  exec mapcal#vecut di4r [ind]
  exec mapcal#vecut di6r [ind]
*  
  aim=$sigma(vsum(air)/[ind])
  bim=$sigma(vsum(bir)/[ind])
  cim=$sigma(vsum(cir)/[ind])
  daa=$sigma(vsum((air-([aim]))*(air-([aim])))/[ind])
  dbb=$sigma(vsum((bir-([bim]))*(bir-([bim])))/[ind])
  dcc=$sigma(vsum((cir-([cim]))*(cir-([cim])))/[ind])
  sa=$sigma(sqrt([daa]))
  sb=$sigma(sqrt([dbb]))
  sc=$sigma(sqrt([dcc]))
*  
  sa=$sigma(vsum(1/di1r))
  aim=$sigma(vsum(air/di1r)/[sa])
  sa=$sigma(1/sqrt([sa]))
  sb=$sigma(vsum(1/di4r))
  bim=$sigma(vsum(bir/di4r)/[sb])
  sb=$sigma(1/sqrt([sb]))
  sc=$sigma(vsum(1/di6r))
  cim=$sigma(vsum(cir/di6r)/[sc])
  sc=$sigma(1/sqrt([sc]))
*  
  n=$vlen(ai)
*  
  ve/cre ixa([n]) r
  ve/cre ixaf([n]) r [n]*1
  sigma ixa=array([n],1#[n])
  sigma ixa=order(ixa,ai)
  nz=$sigma(int(0.1*[n]+0.5))
  do i=1,[nz]
    ve/inp ixaf([i]) 0
    j=[n]-[i]+1
    ve/inp ixaf([j]) 0
  enddo
  sigma ixaf=order(ixaf,ixa)
*  
  ve/cre ixb([n]) r
  ve/cre ixbf([n]) r [n]*1
  sigma ixb=array([n],1#[n])
  sigma ixb=order(ixb,bi)
  nz=$sigma(int(0.1*[n]+0.5))
  do i=1,[nz]
    ve/inp ixbf([i]) 0
    j=[n]-[i]+1
    ve/inp ixbf([j]) 0
  enddo
  sigma ixbf=order(ixbf,ixb)
*  
  ve/cre ixc([n]) r
  ve/cre ixcf([n]) r [n]*1
  sigma ixc=array([n],1#[n])
  sigma ixc=order(ixc,ci)
  nz=$sigma(int(0.1*[n]+0.5))
  do i=1,[nz]
    ve/inp ixcf([i]) 0
    j=[n]-[i]+1
    ve/inp ixcf([j]) 0
  enddo
  sigma ixcf=order(ixcf,ixc)
*
  ve/cre ixf([n]) r
  sigma ixf=ixaf*ixbf*ixcf
  ve/cre air([n]) r
  ve/cre bir([n]) r
  ve/cre cir([n]) r
  ind=0
  do i=1,[n]
    if ($sigma(ixf([i])).eq.1) then
      ind=[ind]+1
      ve/inp air([ind]) $sigma(ai([i]))
      ve/inp bir([ind]) $sigma(bi([i]))
      ve/inp cir([ind]) $sigma(ci([i]))
    endif
  enddo
  exec mapcal#vecut air [ind]
  exec mapcal#vecut bir [ind]
  exec mapcal#vecut cir [ind]
*  
  aim=$sigma(vsum(air)/[ind])
  bim=$sigma(vsum(bir)/[ind])
  cim=$sigma(vsum(cir)/[ind])
  daa=$sigma(vsum((air-([aim]))*(air-([aim])))/[ind])
  dbb=$sigma(vsum((bir-([bim]))*(bir-([bim])))/[ind])
  dcc=$sigma(vsum((cir-([cim]))*(cir-([cim])))/[ind])
  sa=$sigma(sqrt([daa]))
  sb=$sigma(sqrt([dbb]))
  sc=$sigma(sqrt([dcc]))
*  
  mess [ind] [aim] [bim] [cim] [sa] [sb] [sc] [al] [ar]
  ve/cre pars(6) r [aim] [bim] [cim] [sa] [sb] [sc] 
  ve/write pars calfit/calfit_sigs_m[fmode]_n[nopt]_counter[nc].txt
endif
if ([mode].eq.4) then
  fout=calfit/calfit_m[fmode]_n[nopt]_counter[nc].out
  shell /work/snd2000/users/konctbel/SNDSimACC/readlifstat [fout]
  ve/del pars0,dpars0
  ve/read pars0,dpars0 [fout].pars
  n=$vlen(ai)
  ve/cre dn([n]) r
  am=pars0(2)
  bm=pars0(3)
  cm=pars0(4)
  pb=pars0(5)
  qb=pars0(6)
  pc=pars0(7)
  qc=pars0(8)
  n1=9
  n2=[n1]+[n]-1
  ve/del ti
  ve/copy pars0([n1]:[n2]) ti
  ve/del sai
  ve/copy dpars0([n1]:[n2]) sai
  ve/del aif,bif,cif,rif,phi
  ve/cre aif([n]) r
  ve/cre bif([n]) r
  ve/cre cif([n]) r
  sigma aif = [am] + ti
*  sigma bif = [bm] + ([pb])*ti + ([qb])*ti*ti
*  sigma cif = [cm] + ([pc])*ti + ([qc])*ti*ti
  sigma rif = [pb]*ti + ([qb])*ti*ti
  sigma phi = [pc] + ([qc])*ti
  sigma bif = [bm] + rif*cos(phi)
  sigma cif = [cm] + rif*sin(phi)
  sigma dai = sqrt(di1)
  sigma dbi = sqrt(di4)
  sigma dci = sqrt(di6)
*  
szc=0.1 
*  
  gl/imp popt
  fname=calfit/calfit_sigs_m0_n[nopt]_counter[nc].txt
  if ($fexist([fname])) then
    ve/del pars
    ve/read pars [fname]
    amid=pars(1)
    bmid=pars(2)
    cmid=pars(3)
    asig=pars(4)
    bsig=pars(5)
    csig=pars(6)
    nsig=3
    amin=$sigma([amid]-[nsig]*[asig])
    amax=$sigma([amid]+[nsig]*[asig])
    bmin=$sigma([bmid]-[nsig]*[bsig])
    bmax=$sigma([bmid]+[nsig]*[bsig])
    cmin=$sigma([cmid]-[nsig]*[csig])
    cmax=$sigma([cmid]+[nsig]*[csig])
    mess [amin] [amax] [bmin] [bmax] [cmin] [cmax]
  endif
  l=$sigma(vmin(ki))
  r=$sigma(vmax(ki))
  kmin=$sigma([l]-0.05*([r]-([l])))
  kmax=$sigma([r]+0.05*([r]-([l])))
*  
  null [amin] [amax] [cmin] [cmax]
  set pmci 4
  * exec $PER/s#vpl ci dn ai dn sz=[szc] ll=-1 
  set pmci 2
  * exec $PER/s#vpl cif dn aif sai sz=[szc] ll=-1 o=s
  atitle 'a' 'c'
  exec save calfit/c_vs_a_correlation_mcal_m[fmode]_n[nopt]_counter[nc].eps f
  if ([auto].le.0) then
    read x
  endif
*  
  null [amin] [amax] [bmin] [bmax]
  set pmci 4
  * exec $PER/s#vpl bi dn ai dn sz=[szc] ll=-1 
  set pmci 2
  * exec $PER/s#vpl bif dn aif sai sz=[szc] ll=-1 o=s
  atitle 'a' 'b'
  exec save calfit/b_vs_a_correlation_mcal_m[fmode]_n[nopt]_counter[nc].eps f
  if ([auto].le.0) then
    read x
  endif
*  
  null [bmin] [bmax] [cmin] [cmax]
  set pmci 4
  * exec $PER/s#vpl ci dn bi dn sz=[szc] ll=-1 
  set pmci 2
  * exec $PER/s#vpl cif dn bif dn sz=[szc] ll=-1 o=s
  l=[bm]-10
  r=[bm]+10
  tphi=$sigma(tan([pc]))
  line [l] $sigma([cm]-10*([tphi])) [r] $sigma([cm]+10*([tphi]))
  atitle 'b' 'c'
  exec save calfit/c_vs_b_correlation_mcal_m[fmode]_n[nopt]_counter[nc].eps f
  if ([auto].le.0) then
    read x
  endif
*  
  null [kmin] [kmax] [amin] [amax]
  set pmci 4
  * exec $PER/s#vpl ai dai ki dn sz=[szc] ll=-1 
  atitle 'time, days' 'a'
  exec save calfit/a_vs_days_correlation_0_mcal_m[fmode]_n[nopt]_counter[nc].eps f
  null [kmin] [kmax] [amin] [amax]
  set pmci 4
  * exec $PER/s#vpl ai dai ki dn sz=[szc] ll=-1 
  atitle 'time, days' 'a'
  set pmci 2
  * exec $PER/s#vpl aif sai ki dn sz=[szc] ll=-1 o=s
  atitle 'time, days' 'a'
  exec save calfit/a_vs_days_correlation_mcal_m[fmode]_n[nopt]_counter[nc].eps f 
  if ([auto].le.0) then
    read x
  endif
*  
  null [kmin] [kmax] [bmin] [bmax]
  set pmci 4
  * exec $PER/s#vpl bi dbi ki dn sz=[szc] ll=-1 
  atitle 'time, days' 'b'
  exec save calfit/b_vs_days_correlation_0_mcal_m[fmode]_n[nopt]_counter[nc].eps f
  null [kmin] [kmax] [bmin] [bmax]
  set pmci 4
  * exec $PER/s#vpl bi dbi ki dn sz=[szc] ll=-1 
  atitle 'time, days' 'b'
  set pmci 2
  * exec $PER/s#vpl bif dn ki dn sz=[szc] ll=-1 o=s
  atitle 'time, days' 'b'
  exec save calfit/b_vs_days_correlation_mcal_m[fmode]_n[nopt]_counter[nc].eps f 
  if ([auto].le.0) then
    read x
  endif
*  
  null [kmin] [kmax] [cmin] [cmax]
  set pmci 4
  * exec $PER/s#vpl ci dci ki dn sz=[szc] ll=-1 
  atitle 'time, days' 'c'
  exec save calfit/c_vs_days_correlation_0_mcal_m[fmode]_n[nopt]_counter[nc].eps f
  null [kmin] [kmax] [cmin] [cmax]
  set pmci 4
  * exec $PER/s#vpl ci dci ki dn sz=[szc] ll=-1 
  atitle 'time, days' 'c'
  set pmci 2
  * exec $PER/s#vpl cif dn ki dn sz=[szc] ll=-1 o=s
  atitle 'time, days' 'c'
  exec save calfit/c_vs_days_correlation_mcal_m[fmode]_n[nopt]_counter[nc].eps f 
  if ([auto].le.0) then
    read x
  endif
*  
  ve/del abi, abif
  sigma abi = ai + bi
  sigma abif = aif + bif
  null [kmin] [kmax] [amin] [amax]
  set pmci 4
  * exec $PER/s#vpl abi dn ki dn sz=[szc] ll=-1 
  set pmci 2
  * exec $PER/s#vpl abif sai ki dn sz=[szc] ll=-1 o=s
  atitle 'time, days' 'a+b'
  exec save calfit/amp_vs_days_correlation_mcal_m[fmode]_n[nopt]_counter[nc].eps f 
  if ([auto].le.0) then
    read x
  endif
*  
  ve/del rabi, rabif
  sigma rabi = abi/abif
  sigma rabif = abif/abif
  null [kmin] [kmax] 0.9 1.1
  set pmci 4
  * exec $PER/s#vpl rabi dn ki dn sz=[szc] ll=-1 
  set pmci 2
  * exec $PER/s#vpl rabif dn ki dn sz=[szc] ll=-1 o=s
  atitle 'time, days' '(a+b)?cal!/(a+b)?fit!'
  npt=$vlen(rabi)
  mean=$sigma(vsum(rabi)/[npt])
  mean2=$sigma(vsum(rabi**2)/[npt])
  sig=$sigma(sqrt([mean2]-[mean]**2))
  txt=(a+b)?cal!/(a+b)?fit! = $sigma(int([mean]*1000+0.5)/1000) [\261] $sigma(int([sig]*1000+0.5)/1000)
  exec pl#tf 0.1 0.9 [txt]
  exec save calfit/amp_ratio_vs_days_correlation_mcal_m[fmode]_n[nopt]_counter[nc].eps f 
  if ([auto].le.0) then
    read x
  endif
*  
  set pmci 4
  ve/del dmi0
  ve/copy dmi dmi0
  do i=1,$vlen(dmi0)
    if ($sigma(dmi0([i])).gt.0.1) then
      ve/inp dmi0([i]) 0
    endif
  enddo
  * exec $PER/s#vpl mi dmi0 ki dn sz=[szc] ll=-1
  atitle 'time, days' '[m]?max!'
  exec save calfit/mu_vs_days_correlation_mcal_m[fmode]_n[nopt]_counter[nc].eps f
  if ([auto].le.0) then
    read x
  endif
endif
if ([mode].eq.5) then
  fname=/work/users/konctbel/Calibr/total_cals0.tex
  if ($fexist([fname])) then
    shell rm [fname]
  endif
  for/file 20 [fname] ! N
  close 20
  do ni=[n1],[n2]
    year=$substring([cal[ni]],5,4)
    month=$substring([cal[ni]],10,2)
    day=$substring([cal[ni]],13,2)
    cal=$substring([cal[ni]],16,4)
    do nc=1,9
      fmess '\begin{figure}[ht!b]' [fname]
      fmess '  \begin{minipage}{0.49\textwidth}' [fname]
      fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
      epsfile=[cal[ni]]/picts/mu_vs_hv_counter[nc].eps
      fepsfile=/work/users/konctbel/Calibr/[epsfile]
      if ($fexist([fepsfile]).eq.0) then
        epsfile=null.eps
      endif
      txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
      fmess [txt] [fname]
      txt=$unquote('     ')\caption{[year]-[month]-[day] [cal]:[nc]}
      fmess [txt] [fname]
      fmess '   \end{minipage}' [fname]
      fmess '  \begin{minipage}{0.49\textwidth}' [fname]
      fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
      epsfile=[cal[ni]]/picts/amp_vs_hv_counter[nc].eps
      fepsfile=/work/users/konctbel/Calibr/[epsfile]
      if ($fexist([fepsfile]).eq.0) then
        epsfile=null.eps
      endif
      txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
      fmess [txt] [fname]
      txt=$unquote('     ')\caption{[year]-[month]-[day] [cal]:[nc]}
      fmess [txt] [fname]
      fmess '   \end{minipage}' [fname]
      fmess '\end{figure}' [fname]
      if ($sigma(mod([nc],3)).eq.0) then
        fmess '\clearpage' [fname]
      endif
    enddo
  enddo
endif
return


macro newcalnh nc=1 ncal=180
shell rm tmp[0-9].txt
shell ls -d /work/users/konctbel/Calibr/Cal_* > tmp0.txt
shell $unquote('cat tmp0.txt | sed "s/[^0-9]/ /g" > tmp1.txt')
ve/del y,m,d,c
ve/read y,m,d,c tmp1.txt
exec mapcal#ixndl [ncal] c
gl/imp inx
in1=[inx]
exec mapcal#ixndr [ncal] c
gl/imp inx
in2=[inx]
if ([in1].eq.[in2]) then
  year=y([in1])
  month=m([in1])
  day=d([in1])
  cal=$format([ncal],i4.4)
  dir=Cal_$format([year],i4.4)_$format([month],i2.2)_$format([day],i2.2)_[cal]
  exec mapcal#newcaln [nc] [nc] [dir] 1 s 1
endif
return


macro newcalnr im=1
dir0=/work/users/konctbel/Calibr
shell rm tmp[0-9].txt
shell ls -d [dir0]/Cal_* > tmp0.txt
shell $unquote('cat tmp0.txt | sed "s/[^0-9]/ /g" > tmp1.txt')
ve/del y,m,d,c
ve/read y,m,d,c tmp1.txt
do i=[im],$vlen(c)
  year=y([i])
  month=m([i])
  day=d([i])
  ncal=c([i])
  cal=$format([ncal],i4.4)
  dir=Cal_$format([year],i4.4)_$format([month],i2.2)_$format([day],i2.2)_$format([ncal],i4.4)
  do nc=1,9
    ve/del pars
    ve/read pars [dir0]/[dir]/pars_correlation_counter[nc].txt
    im=pars(3)
    if ([im].eq.1) then
      exec mapcal#newcaln [nc] [nc] [dir] 1 s 1
    endif
  enddo
enddo
return


macro newcaln nc1=1 nc2=9 dir=Cal_2013_07_04_0180 auto=0 fopt=sq rtup=1
date=$substring([dir],1,14)
cal=$substring([dir],16,4)
dir0=/work/users/konctbel/Calibr
fname=[dir0]/[dir]/c0.txt
if ($fexist([fname]).eq.0) then
  mess file [fname] das not exists
  goto 1
endif
ve/del c0
ve/read c0 [dir0]/[dir]/c0.txt
np=$vlen(c0)
*np0=[np]
*ve/cre c0([np]) r 20 40 50 80 100 150 200
fname=[dir0]/[dir]/Uset_[cal].txt
if ($fexist([fname]).eq.0) then
  mess file [fname] das not exists
  goto 1
endif
ve/del u01,u02,u03,u04,u05,u06,u07,u08,u09
ve/read u01,u02,u03,u04,u05,u06,u07,u08,u09 [dir0]/[dir]/Uset_[cal].txt
*
ucal=[cal]
if ([cal].eq.155) then
  ucal=165.5
endif
dirw=/work/snd2000/users/martin/snd2k/testrel/results
if ([ucal].le.180) then
  fuwork=[dirw]/U_work.txt
endif
if ([ucal].le.166) then
  fuwork=[dirw]/U_work_before_15_04_13.txt
endif
if ([ucal].le.165.5) then
  fuwork=[dirw]/U_work_before_12_04_13.txt
endif
if ([ucal].le.157) then
  fuwork=[dirw]/U_work_before_04_03_13.txt
endif
if ([ucal].le.155) then
  fuwork=[dirw]/U_work_before_25_02_13.txt
endif
if ([ucal].le.148) then
  fuwork=[dirw]/U_work_before_01_02_13.txt
endif
if ([ucal].le.133) then
  fuwork=[dirw]/U_work_before_30_11_12.txt
endif
if ([ucal].le.90) then
  fuwork=[dirw]/U_work_before_21_11_11.txt
endif
if ([ucal].le.87) then
  fuwork=[dirw]/U_work_before_03_11_11.txt
endif
if ([ucal].le.81) then
  fuwork=[dirw]/U_work_before_25_10_11.txt
endif
if ([ucal].le.64) then
  fuwork=[dirw]/U_work_before_110511.txt
endif
if ([ucal].le.62) then
  fuwork=[dirw]/U_work_before_02_05_11.txt
endif
if ([ucal].le.42) then
  fuwork=[dirw]/U_work_15_12_10.txt
endif
if ([ucal].le.40) then
  fuwork=[dirw]/U_work_07_12_10.txt
endif
if ([ucal].le.35) then
  fuwork=[dirw]/U_work_old.txt
endif
ve/del uw
ve/read uw [fuwork]
*
ve/cre du0([np]) r
ve/cre pds0([np]) r
ve/cre dpds0([np]) r
ve/cre pds0f([np]) r
ve/cre pds([np]) r
ve/cre pdsf([np]) r
ve/cre dpdsf([np]) r
ve/cre ampx([np]) r
ve/cre dampx([np]) r
ve/cre amp([np]) r
ve/cre damp([np]) r
ve/cre amp1([np]) r
ve/cre damp1([np]) r
ve/cre eff([np]) r
ve/cre deff([np]) r
ve/cre mu([np]) r
ve/cre dmu([np]) r
*
ve/cre  mumax(9) r
ve/cre dmumax(9) r
nfit=1
if (([nc1].le.9).and.([nc2].ge.9)) then
  c0i=$sigma(vmax(c0))
  fname=[dir0]/[date]_[cal]/data/cal_pmt_[cal]_[c0i].tup
  if ($fexist([fname]).eq.0) then
    mess file [fname] das not exists
    goto 1
  endif
  hi/file 20 [fname]
  1d 11111 ! 4000 0 4000
  nt/pl //lun20/1.t9 t9>0 idh=11111
  ve/cre vpn(4000) r
  hi/get/cont 11111 vpn
  ve/cre vpx(4000) r
  sigma vpx = array(4000,0.5#3999.5)
  nevt1=$sigma(vsum(vpn))
  mean1=$sigma(vsum(vpx*vpn)/[nevt1])
  rms1=$sigma(vsum((vpx-[mean1])**2*vpn)/[nevt1])
  rms1=$sigma(sqrt([rms1]))
  nt/pl //lun20/1.t10 t10>0 idh=11111
  hi/get/cont 11111 vpn
  nevt2=$sigma(vsum(vpn))
  mean2=$sigma(vsum(vpx*vpn)/[nevt2])
  rms2=$sigma(vsum((vpx-[mean2])**2*vpn)/[nevt2])
  rms2=$sigma(sqrt([rms2]))
  mess nevt1=[nevt1] nevt2=[nevt2] [rms1] [rms2]
  if ([nevt1].gt.[nevt2]) then
    nct0=9
  else
    nct0=10
  endif
  nt/pl //lun20/1.a9 t[nct0]>0 idh=11111
  hi/get/cont 11111 vpn
  nevt1=$sigma(vsum(vpn))
  mean1=$sigma(vsum(vpx*vpn)/[nevt1])
  rms1=$sigma(vsum((vpx-[mean1])**2*vpn)/[nevt1])
  rms1=$sigma(sqrt([rms1]))
  nt/pl //lun20/1.a10 t[nct0]>0 idh=11111
  hi/get/cont 11111 vpn
  nevt2=$sigma(vsum(vpn))
  mean2=$sigma(vsum(vpx*vpn)/[nevt2])
  rms2=$sigma(vsum((vpx-[mean2])**2*vpn)/[nevt2])
  rms2=$sigma(sqrt([rms2]))
  if ([rms1].gt.[rms2]) then
    nca0=9
  else
    nca0=10
  endif
  close 20
endif
do nc=[nc1],[nc2]
*
  if ([rtup].eq.1) then
  do i=1,[np]
    c0i=c0([i])
*    
    if ([nc].ne.9) then
      nca=[nc]
      nct=[nc]
    else
      nca=[nca0]
      nct=[nct0]
    endif
*
    fname=[dir0]/[date]_[cal]/data/cal_pmt_p_[cal]_[c0i].tup
    if ($fexist([fname]).eq.0) then
      mess file [fname] das not exists
      goto 1
    endif
    hi/file 20 [fname]
    idh0=100
    1d [idh0] ! 4000 0 4000
    nt/pl //lun20/1.a[nca] idh=[idh0] option=N
*    exec ntproj [idh0] //lun20/1.a[nca] ! 4000 0 4000 
    idh0t=200
    1d [idh0t] ! 4000 0 4000
    nt/pl //lun20/1.a[nca] t[nct]>0 idh=[idh0t] option=N 
*    exec ntproj [idh0t] //lun20/1.a[nca] t[nct]>0 4000 0 4000
    idh00=201
    1d [idh00] ! 4000 0 4000
    nt/pl //lun20/1.a[nca] t[nct].eq.0 idh=[idh00] option=N 
*    exec ntproj [idh00] //lun20/1.a[nca] [nct].eq.0 4000 0 4000
    ve/inp pds0([i]) $hinfo([idh0],'mean')
    nx=$hinfo([idh0],'events')
    rms=$hinfo([idh0],'rms')
    ve/inp dpds0([i]) $sigma([rms]/sqrt([nx]))
    close 20
*    
    fname=[dir0]/[date]_[cal]/data/cal_pmt_[cal]_[c0i].tup
    if ($fexist([fname]).eq.0) then
      mess file [fname] das not exists
      goto 1
    endif
    hi/file 20 [fname]
    idhx=300
    1d [idhx] ! 4000 0 4000
    nt/pl //lun20/1.a[nca] idh=[idhx] option=N 
*    exec ntproj [idhx] //lun20/1.a[nca] ! 4000 0 4000
    idhxt=400
    1d [idhxt] ! 4000 0 4000
    nt/pl //lun20/1.a[nca] t[nct]>0 idh=[idhxt] option=N 
*    exec ntproj [idhxt] //lun20/1.a[nca] t[nct]>0 4000 0 4000
    idhx0=401
    1d [idhx0] ! 4000 0 4000
    nt/pl //lun20/1.a[nca] t[nct].eq.0 idh=[idhx0] option=N 
*    exec ntproj [idhx0] //lun20/1.a[nca] t[nct].eq.0 4000 0 4000
    ve/inp ampx([i]) $hinfo([idhx],'mean')
    nx=$hinfo([idhx],'events')
    rms=$hinfo([idhx],'rms')
    ve/inp dampx([i]) $sigma([rms]/sqrt([nx]))
    close 20
*    
    nx=$hinfo([idh00],'xbins')
    xmin=$hinfo([idh00],'xmin')
    xmax=$hinfo([idh00],'xmax')
    ve/cre vx([nx]) r
    hi/get/cont [idh00] vx
    ve/cre ix([nx]) r
    sigma ix = array([nx],1#[nx])
    sigma ix = order(ix,-vx)
    sigma vx = order(vx,-vx)
    mean=$sigma((ix(1)-0.5)*([xmax]-[xmin])/[nx])
    ve/cre p3(3) r $sigma(vx(1)) [mean] 2
    do j=1,[nfit]
      l=$sigma(p3(2)-3*abs(p3(3)))
      r=$sigma(p3(2)+1.5*abs(p3(3)))
      hi/fit [idh00]([l]:[r]) g [fopt] 3 p3
    enddo
    ve/inp pds0f([i]) $sigma(p3(2))
    if ([auto].le.0) then
      l=$sigma(p3(2)-10*abs(p3(3)))
      r=$sigma(p3(2)+10*abs(p3(3)))
      hi/pl [idh00]([l]:[r])
      read x
    endif
*    
    nx=$hinfo([idhx0],'xbins')
    xmin=$hinfo([idhx0],'xmin')
    xmax=$hinfo([idhx0],'xmax')
    ve/cre vx([nx]) r
    hi/get/cont [idhx0] vx
    ve/cre ix([nx]) r
    sigma ix = array([nx],1#[nx])
    sigma ix = order(ix,-vx)
    sigma vx = order(vx,-vx)
    mean=$sigma((ix(1)-0.5)*([xmax]-[xmin])/[nx])
    ve/cre p3(3) r $sigma(vx(1)) [mean] 2
    do j=1,[nfit]
      l=$sigma(p3(2)-3*abs(p3(3)))
      r=$sigma(p3(2)+1.5*abs(p3(3)))
      hi/fit [idhx0]([l]:[r]) g [fopt] 3 p3
    enddo
    ve/inp pdsf([i]) $sigma(p3(2))
    if ([auto].le.0) then
      l=$sigma(p3(2)-10*abs(p3(3)))
      r=$sigma(p3(2)+10*abs(p3(3)))
      hi/pl [idhx0]([l]:[r])
      read x
    endif
*
    n0=$hinfo([idh0],'events')
    n0t=$hinfo([idh0t],'events')
    e0=$sigma([n0t]/[n0])
    de0=$sigma(sqrt([e0]*(1-[e0])/[n0]))
    n1=$hinfo([idhx],'events')
    n1t=$hinfo([idhxt],'events')
    e1=$sigma([n1t]/[n1])
    de1=$sigma(sqrt([e1]*(1-[e1])/[n1]))
    effi=$sigma(1-(1-[e1])/(1-[e0]))
    pe0=$sigma(([de0]/(1-[e0]))**2)
    pe1=$sigma(([de1]/(1-[e1]))**2)
    deffi=$sigma((1-[effi])*sqrt([pe0]+[pe1]))
    mess eff0=[e0] [de0]
    mess effx=[e1] [de1]
    mess effi=[effi] [deffi]
    ve/inp eff([i]) [effi]
    ve/inp deff([i]) [deffi]
    if ([auto].le.0) then
      read x
    endif
*    
  enddo
  sigma pds = pds0 - pds0f + pdsf
  sigma dpds = dpds0
  sigma  amp = ampx - pds
  sigma damp = sqrt(dampx**2+dpds**2)
  sigma mu = -log(1-eff)
  sigma dmu = deff/(1-eff)
  else
    ve/del u0[nc],amp1,damp1,amp,damp,mu,dmu,eff,deff
    ve/read u0[nc],amp1,damp1,amp,damp,mu,dmu,eff,deff [dir0]/[date]_[cal]/fit_data_counter[nc].txt '(9f15.6)'
  endif  
*
*  if ($sigma(vmin(u0[nc])).lt.1000) then
*    ind=0
*    do ii=1,[np]
*      if ($sigma(u0[nc]([ii])).gt.1000) then
*        ind=[ind]+1
*        ve/inp u0[nc]([ind]) $sigma(u0[nc]([ii]))
*        ve/inp du0[nc]([ind]) $sigma(du0[nc]([ii]))
*        ve/inp mu([ind]) $sigma(mu([ii]))
*        ve/inp dmu([ind]) $sigma(dmu([ii]))
*        ve/inp amp([ind]) $sigma(amp([ii]))
*        ve/inp damp([ind]) $sigma(damp([ii]))
*      endif
*    enddo
*    exec mapcal#vecut u0[nc] [ind]
*    exec mapcal#vecut du0[nc] [ind]
*    exec mapcal#vecut mu [ind]
*    exec mapcal#vecut dmu [ind]
*    exec mapcal#vecut amp [ind]
*    exec mapcal#vecut damp [ind]
*    np=[ind]
*  endif
*
  exec vpl#pl0 mu dmu u0[nc] du0 
  atitle 'U?0!, V' '[m], pe.'
  mmax=mu([np])
  mmid=[mmax]/2
  sigma rmu = abs(mu-[mmid])
  ve/cre ix([np]) r
  sigma ix = array([np],1#[np])
  sigma ix = order(ix,rmu)
  ic=ix(1)
  il=$sigma(max(1,[ic]-1))
  ir=$sigma(min([np],[ic]+1))
  dmudx=$sigma((mu([ir])-mu([il]))/(u0[nc]([ir])-u0[nc]([il])))
  sig=$sigma(abs([mmax]/sqrt(2*3.1415927)/([dmudx])))
  ve/cre p4(4) r 0 [mmax] $sigma(u0[nc]([ic])) [sig]
  ve/cre s4(4) r 0 0 10 0
  ve/fit u0[nc] mu dmu erfp0.f sb 4 p4 s4
  ve/cre s4(4) r 0 0 10 10
  ve/fit u0[nc] mu dmu erfp0.f sb 4 p4 s4
  ve/cre s4(4) r 1 1 10 10
  ve/fit u0[nc] mu dmu erfp0.f sb 4 p4 s4
*  
  npar=4
  ve/cre chi2(2) r
  ve/cre paru([npar]) r
  ve/cre dparu([npar]) r
  ve/cre covu([npar],[npar]) r
  do j=1,[nfit]
    ve/fit u0[nc] mu dmu erfp0.f s 4 p4
  enddo
  call covm.f(1)
  call covmpen(chi2,[npar],paru,dparu)
*  
  ve/del dmuo
  ve/copy dmu dmuo
  ve/cre uch(1) r $sigma(chi2(1)*chi2(2))
  ve/cre nch(1) i $sigma(chi2(2))
  prob=$call('prob(uch,nch)')
  im=0
  if ([prob].lt.0.005) then
  chi2m=1000000
  do ii=1,[np]
    ve/del p4o
    ve/copy p4 p4o
    dmut=dmu([ii])
    ve/inp dmu([ii]) 0
    ve/fit u0[nc] mu dmu erfp0.f s 4 p4o
    call covm.f(1)
    call covmpen(chi2,[npar],paru,dparu)
    chi2i=chi2(1)
    if ([chi2i].lt.[chi2m]) then
      chi2m=[chi2i]
      im=[ii]
    endif
    ve/inp dmu([ii]) [dmut]
  enddo
  ve/inp dmu([im]) 0
  endif
*  
*  ve/cre p4(4) r 0 [mmax] $sigma(u0[nc]([ic])) [sig]
*  ve/cre s4(4) r 0 0 10 10
  ve/fit u0[nc] mu dmu erfp0.f sb 4 p4 s4
  ve/write mu,dmu,dmuo ! 3f15.6
  exec vpl#pl0 mu dmuo u0[nc] du0
  atitle 'U?0!, V' '[m], pe.'
  do j=1,[nfit]
    ve/fit u0[nc] mu dmu erfp0.f s 4 p4
  enddo
  ve/inp p4(4) $sigma(abs(p4(4)))
  call covm.f(1)
  call covmpen(chi2,[npar],paru,dparu)
  call covmcov([npar],covu)
  ve/write paru,dparu [dir0]/[date]_[cal]/fit_mu_counter[nc].txt '(2f15.6)'
  sig=$sigma(dparu(2)*sqrt(chi2(1)*chi2(2)))
  ve/inp mumax([nc]) $sigma(p4(2))
  ve/inp dmumax([nc]) $sigma([sig])
  mui=mumax([nc])
  dmui=dmumax([nc])
  mess mumax[nc]=[mui] [dmui]
*  
  l=$GRAFINFO('WNXMIN')
  r=$GRAFINFO('WNXMAX')
  d=$GRAFINFO('WNYMIN')
  u=$GRAFINFO('WNYMAX')
  u=$sigma([u]+0.1*([u]-[d]))
  null [l] [r] [d] [u]
  exec vpl#pl0 mu dmu u0[nc] du0 o=s
  atitle 'U?0!, V' '[m], pe.'
  ve/fit u0[nc] mu dmu erfp0.f s 4 p4
  fun/pl erfp0p.f $sigma(vmin(u0[nc])) $sigma(vmax(u0[nc])) s
  set dmod 5
  line [l] [mui] [r] [mui]
  set dmod 1
  txt=[m]?[nc],max! = $sigma(int([mui]*1000+0.5)/1000) [\261] $sigma(int([dmui]*1000+0.5)/1000)
  exec pl#tf 0.55 0.9 [txt]
*  
  ndf=$sigma(chi2(2))
  chi=$sigma(int(chi2(1)*[ndf]*1000+0.5)/1000)
  txt=[h]^2!/ndf = [chi] / [ndf]
*  exec $PER/s#tf 0.05 0.9 [txt]
  exec pl#tf 0.05 0.9 [txt]
*  
  exec save [dir0]/[date]_[cal]/picts/mu_vs_hv_counter[nc].eps f
  ve/del pars
  ve/cre imx(1) r [im]
  exec vappend pars chi2
  exec vappend pars imx
  exec vappend pars paru
  exec vappend pars dparu
  do j=1,[npar]
    ve/copy covu([j]) p4
    exec vappend pars p4
  enddo
  ve/write pars [dir0]/[date]_[cal]/pars_mu_correlation_counter[nc].txt
*  
  if ([auto].le.1) then
    read x
  endif
  sigma amp1 = amp/[mui]
*  sigma damp1 = amp1*sqrt((damp/amp)**2+([dmui]/([mui]))**2)
  sigma damp1 = damp/[mui]
  ve/del damp2
  ve/copy damp1 damp2
  do j=1,[np]
    if ($sigma(amp1([j])).lt.0) then
      ve/inp damp1([j]) 0
    endif
  enddo
  exec vpl#pl0 amp1 damp1 u0[nc] du0
  atitle 'U?0!, V' 'Amplitude, ch.'
  ve/del yi,dyi,xi
  ve/cre xi([np]) r
  ve/cre yi([np]) r
  ve/cre dyi([np]) r
  sigma xi = log(u0[nc]/3000)
  sigma yi = log(amp)
  sigma dyi = damp/amp
  ve/cre p2(2) r 3 23
  ve/fit xi yi dyi p1 s 2 p2 
*
  npar=3
  ve/cre chi2(2) r
  ve/cre paru([npar]) r
  ve/cre dparu([npar]) r
  ve/cre covu([npar],[npar]) r
  do j=1,[nfit]
    ve/fit u0[nc] amp1 damp1 ampcal.f s 3 p3
  enddo
  ve/inp p3(1) $sigma(abs(p3(1)))
  ve/inp p3(3) $sigma(abs(p3(3)))
  call covm.f(1)
  call covmpen(chi2,[npar],paru,dparu)
*  
  ve/del damp1o
  ve/copy damp1 damp1o
  ve/cre uch(1) r $sigma(chi2(1)*chi2(2))
  ve/cre nch(1) i $sigma(chi2(2))
  prob=$call('prob(uch,nch)')
  im=0
  if ([prob].lt.0.005) then
  chi2m=1000000
  ampi1=amp1(1)
  ampi2=amp1(2)
  if ([ampi1].lt.[ampi2]) then
    iimin=2
  else
    iimin=1
  endif
  do ii=[iimin],[np]
    ve/cre p3(3) r $sigma(exp(p2(1))) 0 $sigma(p2(2))
    damp1t=damp1([ii])
    ve/inp damp1([ii]) 0
    ve/fit u0[nc] amp1 damp1 ampcal.f s 3 p3
    call covm.f(1)
    call covmpen(chi2,[npar],paru,dparu)
    chi2i=chi2(1)
    if ([chi2i].lt.[chi2m]) then
      chi2m=[chi2i]
      im=[ii]
    endif
    ve/inp damp1([ii]) [damp1t]
  enddo
  ve/inp damp1([im]) 0
  endif
*
  exec vpl#pl0 amp1 damp1o u0[nc] du0
  atitle 'U?0!, V' 'Amplitude, ch.'
  ve/cre p3(3) r $sigma(exp(p2(1))) 0 $sigma(p2(2))
  do j=1,[nfit]
    ve/fit u0[nc] amp1 damp1 ampcal.f s 3 p3
  enddo
  ve/inp p3(1) $sigma(abs(p3(1)))
  ve/inp p3(3) $sigma(abs(p3(3)))
  fun/pl ampcalp.f $sigma(vmin(u0[nc])) $sigma(vmax(u0[nc])) s
  call covm.f(1)
  call covmpen(chi2,[npar],paru,dparu)
  call covmcov([npar],covu)
  ve/inp paru(1) $sigma(abs(paru(1)))
  ve/inp paru(3) $sigma(abs(paru(3)))
  ve/write paru,dparu [dir0]/[date]_[cal]/fit_amp_counter[nc].txt '(2f15.6)'
  ve/del pars
  ve/cre imx(1) r [im]
  exec vappend pars chi2
  exec vappend pars imx
  exec vappend pars paru
  exec vappend pars dparu
  do j=1,[npar]
    ve/copy covu([j]) p3
    exec vappend pars p3
  enddo
  ve/del damp1
  ve/copy damp2 damp1
  ve/write pars [dir0]/[date]_[cal]/pars_correlation_counter[nc].txt
  ve/write u0[nc],amp1,damp1,amp,damp,mu,dmu,eff,deff [dir0]/[date]_[cal]/fit_data_counter[nc].txt '(9f15.6)'
*
  nva=0
  nvb=0
  nvc=0
  nvd=0
  mode=1
  if ([mode].eq.0) then
  fcal=/work/snd2000/onldat/acc/a_[cal].txt
  if ($fexist([fcal])) then
    nva=1
    shell cp [fcal] tmp.txt
    shell $unquote('cat tmp.txt | sed "s/nan/0/g" > tmp1.txt')
    shell $unquote('cat tmp1.txt | sed "s/inf/0/g" > tmp.txt')
    fcal=tmp.txt
    ve/del vai
    ve/read vai [fcal]
  endif
  fcal=/work/snd2000/onldat/acc/b_[cal].txt
  if ($fexist([fcal])) then
    nvb=1
    shell cp [fcal] tmp.txt
    shell $unquote('cat tmp.txt | sed "s/nan/0/g" > tmp1.txt')
    shell $unquote('cat tmp1.txt | sed "s/inf/0/g" > tmp.txt')
    fcal=tmp.txt
    ve/del vbi
    ve/read vbi [fcal]
  endif
  fcal=/work/snd2000/onldat/acc/c_[cal].txt
  if ($fexist([fcal])) then
    nvc=1
    shell cp [fcal] tmp.txt
    shell $unquote('cat tmp.txt | sed "s/nan/0/g" > tmp1.txt')
    shell $unquote('cat tmp1.txt | sed "s/inf/0/g" > tmp.txt')
    fcal=tmp.txt
    ve/del vci
    ve/read vci [fcal]
  endif
  fcal=/work/snd2000/onldat/acc/amp1pe_[cal].txt
  if ($fexist([fcal])) then
    nvd=1
    shell cp [fcal] tmp.txt
    shell $unquote('cat tmp.txt | sed "s/nan/0/g" > tmp1.txt')
    shell $unquote('cat tmp1.txt | sed "s/inf/0/g" > tmp.txt')
    fcal=tmp.txt
    ve/del vamp1pe
    ve/read vamp1pe [fcal]
  endif
  fcal=/work/snd2000/users/martin/snd2k/testrel/results/a_[cal].txt
  if ($fexist([fcal])) then
    nva=1
    shell cp [fcal] tmp.txt
    shell $unquote('cat tmp.txt | sed "s/nan/0/g" > tmp1.txt')
    shell $unquote('cat tmp1.txt | sed "s/inf/0/g" > tmp.txt')
    fcal=tmp.txt
    ve/del vai
    ve/read vai [fcal]
  endif
  fcal=/work/snd2000/users/martin/snd2k/testrel/results/b_[cal].txt
  if ($fexist([fcal])) then
    nvb=1
    shell cp [fcal] tmp.txt
    shell $unquote('cat tmp.txt | sed "s/nan/0/g" > tmp1.txt')
    shell $unquote('cat tmp1.txt | sed "s/inf/0/g" > tmp.txt')
    fcal=tmp.txt
    ve/del vbi
    ve/read vbi [fcal]
  endif
  fcal=/work/snd2000/users/martin/snd2k/testrel/results/c_[cal].txt
  if ($fexist([fcal])) then
    nvc=1
    shell cp [fcal] tmp.txt
    shell $unquote('cat tmp.txt | sed "s/nan/0/g" > tmp1.txt')
    shell $unquote('cat tmp1.txt | sed "s/inf/0/g" > tmp.txt')
    fcal=tmp.txt
    ve/del vci
    ve/read vci [fcal]
  endif
  fcal=/work/snd2000/users/martin/snd2k/testrel/results/amp1pe_[cal].txt
  if ($fexist([fcal])) then
    nvd=1
    shell cp [fcal] tmp.txt
    shell $unquote('cat tmp.txt | sed "s/nan/0/g" > tmp1.txt')
    shell $unquote('cat tmp1.txt | sed "s/inf/0/g" > tmp.txt')
    fcal=tmp.txt
    ve/del vamp1pe
    ve/read vamp1pe [fcal]
  endif
  fcal=figvam
  fcal0=/work/snd2000/users/martin/snd2k/testrel/results/rez_cal_pmt_[cal].txt
  if ($fexist([fcal0])) then
    fcal=[fcal0]
  endif
  fcal0=/work/snd2000/onldat/acc/rez_cal_pmt_[cal].txt
  if ($fexist([fcal0])) then
    fcal=[fcal0]
  endif
  if ($fexist([fcal])) then
    ve/del vai,vbi,vci
    ve/cre vai(9) r
    ve/cre vbi(9) r
    ve/cre vci(9) r
    shell grep -E "^par3_" [fcal] >& tmp.txt
    shell $unquote('cat tmp.txt | sed "s/nan/1000000/g" > tmp1.txt')
    shell $unquote('cat tmp1.txt | sed "s/inf/1000000/g" > tmp.txt')
    shell $unquote('cat tmp.txt | sed "s/par3_/ /g" > tmp1.txt')
    shell cat tmp0.txt tmp1.txt > tmp2.txt
    ve/del np,vp
    ve/read np,vp tmp2.txt
    if ($vlen(vp).eq.30) then
      nva=1
      nvb=1
      nvc=1
      di=0
      do i=1,9
        i3=$sigma(([i]-[di])*3+1)
        ve/inp vai([i]) $sigma(vp([i3]))
        i3=$sigma(([i]-[di])*3+2)
        ve/inp vbi([i]) $sigma(vp([i3]))
        i3=$sigma(([i]-[di])*3+3)
        ve/inp vci([i]) $sigma(vp([i3]))
      enddo
    endif
  endif
  endif
  if ([nvd].eq.0) then
    exec mapcal#calsort
    exec mapcal#ixndr [cal] vcal
    gl/imp inx
    in1=[inx]
    exec mapcal#ixndl [cal] vcal
    gl/imp inx
    in2=[inx]
    if ([in1].eq.[in2]) then
      nvd=1
      q=ae
      ve/del v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9
      ve/read v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9 accled_amp1pe_v1.txt 9f10.3
      ve/cre vamp1pe(9) r
      do i=1,9
        ve/inp vamp1pe([i]) $sigma(v[q][i]([in1]))
      enddo
    endif
  endif
*
  if ([mode].eq.1) then
    n=0
    n=[n]+1; v[n]=a           ; vn[n]=v[v[n]]
    n=[n]+1; v[n]=b           ; vn[n]=v[v[n]]
    n=[n]+1; v[n]=c           ; vn[n]=v[v[n]]
    n=[n]+1; v[n]=amp1pe      ; vn[n]=[v[n]]
    do v=1,[n]
      shell cp [v[v]]_calslist.txt tmp.txt
      shell $unquote('cat tmp.txt | sed "s/\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*/       0.000000/g" > tmp1.txt')
      ve/del ixc,[vn[v]]1,[vn[v]]2,[vn[v]]3,[vn[v]]4,[vn[v]]5,[vn[v]]6,[vn[v]]7,[vn[v]]8,[vn[v]]9
      ve/read ixc,[vn[v]]1,[vn[v]]2,[vn[v]]3,[vn[v]]4,[vn[v]]5,[vn[v]]6,[vn[v]]7,[vn[v]]8,[vn[v]]9 tmp1.txt '(10f15.6)'
      vv1=[vn[v]]1
      if ($vexist([vv1]).ne.1) then
        wait 'Wait a ten seconds...' 10
        ve/del ixc,[vn[v]]1,[vn[v]]2,[vn[v]]3,[vn[v]]4,[vn[v]]5,[vn[v]]6,[vn[v]]7,[vn[v]]8,[vn[v]]9
        ve/read ixc,[vn[v]]1,[vn[v]]2,[vn[v]]3,[vn[v]]4,[vn[v]]5,[vn[v]]6,[vn[v]]7,[vn[v]]8,[vn[v]]9 tmp1.txt '(10f15.6)'
      endif
    enddo
    ve/del vai,vbi,vci,vamp1pei
    ve/cre vai(9) r
    ve/cre vbi(9) r
    ve/cre vci(9) r
    ve/cre vamp1pei(9) r
    ixcmin=$sigma(vmin(ixc))
    ixcmax=$sigma(vmax(ixc))
    if (([cal].ge.[ixcmin]).and.([cal].le.[ixcmax])) then
      exec mapcal#ixndl [cal] ixc
      gl/imp inx
      in1=[inx]
      exec mapcal#ixndr [cal] ixc
      gl/imp inx
      in2=[inx]
      if ([in1].eq.[in2]) then
        do v=1,[n]
          do i=1,9
            ve/inp v[v[v]]i([i]) $sigma([vn[v]][i]([in1]))
          enddo
        enddo
      endif
    endif
    nva=$sigma(vsum(vai)/vsum(vai))
    nvb=$sigma(vsum(vbi)/vsum(vbi))
    nvc=$sigma(vsum(vci)/vsum(vci))
    ve/copy vamp1pei vxi
    nvd=$sigma(vsum(vxi)/vsum(vxi))
  endif
*
  a=ndf
  b=ndf
  c=ndf
  amp1pe=ndf
  nv=[nva]+[nvb]+[nvc]
  if ([nv].eq.3) then
    a=vai([nc])
    b=vbi([nc])
    c=vci([nc])
    set hcol 2
    set dmod 2
    set basl 0.01
    set ltyp 12
    fun/pl [a]*(x/3000)**([c])+([b]) $sigma(vmin(u0[nc])) $sigma(vmax(u0[nc])) s
    set hcol 1
    set dmod 1
    set ltyp 1
    set pmci 2
    if ([nvd].eq.1) then
      amp1pe=vamp1pei([nc])
      u=$sigma(3000*(([amp1pe]-([b]))/([a]))**(1/([c])))
      set ksiz 0.2
      key [u] $sigma([a]*([u]/3000)**([c])+([b])) 20 ! 1.0
      set pmci 1
      txt=A?1! = $sigma(int(100*[amp1pe]+0.5)/100) (db)
*      exec $PER/s#tf 0.05 0.8 [txt]
      exec pl#tf 0.05 0.7 [txt]
    endif
  endif
*  
  uwi=uw([nc])
  a=paru(1)
  b=paru(2)
  c=paru(3)
  amp1pe=$sigma([a]*([uwi]/3000)**([c])+([b]))
  txt=A?1! = $sigma(int(100*[amp1pe]+0.5)/100) (my)
*  exec $PER/s#tf 0.05 0.7 [txt] 
  exec pl#tf 0.05 0.6 [txt]
  set pmci 4
  set ksiz 0.2
  key [uwi] [amp1pe] 20 ! 1.0
  set pmci 1
*  
  txt=counter [nc]   U?w! = [uwi] V
  exec pl#tf 0.05 0.9 [txt]
  ndf=$sigma(chi2(2))
  chi=$sigma(int(chi2(1)*[ndf]*1000+0.5)/1000)
  txt=[h]^2!/ndf = [chi] / [ndf]
*  exec $PER/s#tf 0.05 0.9 [txt]
  exec pl#tf 0.05 0.8 [txt]
  exec pl#tf 0.05 0.5 a?db!=$sigma(vai([nc]))
  exec pl#tf 0.05 0.4 b?db!=$sigma(vbi([nc]))
  exec pl#tf 0.05 0.3 c?db!=$sigma(vci([nc]))
  exec pl#tf 0.7 0.25 a?my!=$sigma(paru(1))
  exec pl#tf 0.7 0.15 b?my!=$sigma(paru(2))
  exec pl#tf 0.7 0.05 c?my!=$sigma(paru(3))
*
  exec save [dir0]/[date]_[cal]/picts/amp_vs_hv_counter[nc].eps f
  if ([auto].le.1) then
    read x
  endif  
  mode=1
  if ([mode].eq.0) then
  exec mapcal#elr x y paru covu
  atitle 'a' 'b'
  exec save [dir0]/[date]_[cal]/picts/a_vs_b_correlation_counter[nc].eps f
  if ([auto].le.1) then
    read x
  endif  
  exec mapcal#elr x z paru covu
  atitle 'a' 'c'
  exec save [dir0]/[date]_[cal]/picts/a_vs_c_correlation_counter[nc].eps f
  if ([auto].le.1) then
    read x
  endif  
  exec mapcal#elr y z paru covu
  atitle 'b' 'c'
  exec save [dir0]/[date]_[cal]/picts/b_vs_c_correlation_counter[nc].eps f
  if ([auto].le.1) then
    read x
  endif  
  endif
*  np=[np0]
enddo
1:
return


macro newcalnp nc1=1 nc2=9 dir=Cal_2013_07_04_0180 auto=0 fopt=sq 
date=$substring([dir],1,14)
cal=$substring([dir],16,4)
dir0=/work/users/konctbel/Calibr
fname=[dir0]/[dir]/c0.txt
if ($fexist([fname]).eq.0) then
  mess file [fname] das not exists
  goto 1
endif
ve/del c0
ve/read c0 [dir0]/[dir]/c0.txt
np=$vlen(c0)
*np0=[np]
*ve/cre c0([np]) r 20 40 50 80 100 150 200
fname=[dir0]/[dir]/Uset_[cal].txt
if ($fexist([fname]).eq.0) then
  mess file [fname] das not exists
  goto 1
endif
ve/del u01,u02,u03,u04,u05,u06,u07,u08,u09
ve/read u01,u02,u03,u04,u05,u06,u07,u08,u09 [dir0]/[dir]/Uset_[cal].txt
*
ucal=[cal]
if ([cal].eq.155) then
  ucal=165.5
endif
dirw=/work/snd2000/users/martin/snd2k/testrel/results
if ([ucal].le.180) then
  fuwork=[dirw]/U_work.txt
endif
if ([ucal].le.166) then
  fuwork=[dirw]/U_work_before_15_04_13.txt
endif
if ([ucal].le.165.5) then
  fuwork=[dirw]/U_work_before_12_04_13.txt
endif
if ([ucal].le.157) then
  fuwork=[dirw]/U_work_before_04_03_13.txt
endif
if ([ucal].le.155) then
  fuwork=[dirw]/U_work_before_25_02_13.txt
endif
if ([ucal].le.148) then
  fuwork=[dirw]/U_work_before_01_02_13.txt
endif
if ([ucal].le.133) then
  fuwork=[dirw]/U_work_before_30_11_12.txt
endif
if ([ucal].le.90) then
  fuwork=[dirw]/U_work_before_21_11_11.txt
endif
if ([ucal].le.87) then
  fuwork=[dirw]/U_work_before_03_11_11.txt
endif
if ([ucal].le.81) then
  fuwork=[dirw]/U_work_before_25_10_11.txt
endif
if ([ucal].le.64) then
  fuwork=[dirw]/U_work_before_110511.txt
endif
if ([ucal].le.62) then
  fuwork=[dirw]/U_work_before_02_05_11.txt
endif
if ([ucal].le.42) then
  fuwork=[dirw]/U_work_15_12_10.txt
endif
if ([ucal].le.40) then
  fuwork=[dirw]/U_work_07_12_10.txt
endif
if ([ucal].le.35) then
  fuwork=[dirw]/U_work_old.txt
endif
ve/del uw
ve/read uw [fuwork]
*
ve/cre du0([np]) r
ve/cre pds0([np]) r
ve/cre dpds0([np]) r
ve/cre pds0f([np]) r
ve/cre pds([np]) r
ve/cre pdsf([np]) r
ve/cre dpdsf([np]) r
ve/cre ampx([np]) r
ve/cre dampx([np]) r
ve/cre amp([np]) r
ve/cre damp([np]) r
ve/cre amp1([np]) r
ve/cre damp1([np]) r
ve/cre eff([np]) r
ve/cre deff([np]) r
ve/cre mu([np]) r
ve/cre dmu([np]) r
ve/cre tmin([np]) r
ve/cre tmax([np]) r
*
ve/cre  mumax(9) r
ve/cre dmumax(9) r
nfit=1
if (([nc1].le.9).and.([nc2].ge.9)) then
  c0i=$sigma(vmax(c0))
  fname=[dir0]/[date]_[cal]/data/cal_pmt_[cal]_[c0i].tup
  if ($fexist([fname]).eq.0) then
    mess file [fname] das not exists
    goto 1
  endif
  hi/file 20 [fname]
  1d 11111 ! 4000 0 4000
  nt/pl //lun20/1.t9 t9>0 idh=11111
  ve/cre vpn(4000) r
  hi/get/cont 11111 vpn
  ve/cre vpx(4000) r
  sigma vpx = array(4000,0.5#3999.5)
  nevt1=$sigma(vsum(vpn))
  mean1=$sigma(vsum(vpx*vpn)/[nevt1])
  rms1=$sigma(vsum((vpx-[mean1])**2*vpn)/[nevt1])
  rms1=$sigma(sqrt([rms1]))
  nt/pl //lun20/1.t10 t10>0 idh=11111
  hi/get/cont 11111 vpn
  nevt2=$sigma(vsum(vpn))
  mean2=$sigma(vsum(vpx*vpn)/[nevt2])
  rms2=$sigma(vsum((vpx-[mean2])**2*vpn)/[nevt2])
  rms2=$sigma(sqrt([rms2]))
  mess nevt1=[nevt1] nevt2=[nevt2] [rms1] [rms2]
  if ([nevt1].gt.[nevt2]) then
    nct0=9
  else
    nct0=10
  endif
  nt/pl //lun20/1.a9 t[nct0]>0 idh=11111
  hi/get/cont 11111 vpn
  nevt1=$sigma(vsum(vpn))
  mean1=$sigma(vsum(vpx*vpn)/[nevt1])
  rms1=$sigma(vsum((vpx-[mean1])**2*vpn)/[nevt1])
  rms1=$sigma(sqrt([rms1]))
  nt/pl //lun20/1.a10 t[nct0]>0 idh=11111
  hi/get/cont 11111 vpn
  nevt2=$sigma(vsum(vpn))
  mean2=$sigma(vsum(vpx*vpn)/[nevt2])
  rms2=$sigma(vsum((vpx-[mean2])**2*vpn)/[nevt2])
  rms2=$sigma(sqrt([rms2]))
  if ([rms1].gt.[rms2]) then
    nca0=9
  else
    nca0=10
  endif
  close 20
endif
do nc=[nc1],[nc2]
*
  do i=1,[np]
    c0i=c0([i])
*    
    if ([nc].ne.9) then
      nca=[nc]
      nct=[nc]
    else
      nca=[nca0]
      nct=[nct0]
    endif
*
    fname=[dir0]/[date]_[cal]/data/cal_pmt_p_[cal]_[c0i].tup
    if ($fexist([fname]).eq.0) then
      mess file [fname] das not exists
      goto 1
    endif
    hi/file 20 [fname]
    idh0=100
    1d [idh0] ! 4000 0 4000
    nt/pl //lun20/1.a[nca] idh=[idh0] option=N
*    exec ntproj [idh0] //lun20/1.a[nca] ! 4000 0 4000 
    idh0t=200
    1d [idh0t] ! 4000 0 4000
    nt/pl //lun20/1.a[nca] t[nct]>0 idh=[idh0t] option=N 
*    exec ntproj [idh0t] //lun20/1.a[nca] t[nct]>0 4000 0 4000
    idh00=201
    1d [idh00] ! 4000 0 4000
    nt/pl //lun20/1.a[nca] t[nct].eq.0 idh=[idh00] option=N 
*    exec ntproj [idh00] //lun20/1.a[nca] [nct].eq.0 4000 0 4000
    ve/inp pds0([i]) $hinfo([idh0],'mean')
    nx=$hinfo([idh0],'events')
    rms=$hinfo([idh0],'rms')
    ve/inp dpds0([i]) $sigma([rms]/sqrt([nx]))
    close 20
*    
    fname=[dir0]/[date]_[cal]/data/cal_pmt_[cal]_[c0i].tup
    if ($fexist([fname]).eq.0) then
      mess file [fname] das not exists
      goto 1
    endif
    hi/file 20 [fname]
    idhx=300
    1d [idhx] ! 4000 0 4000
    nt/pl //lun20/1.a[nca] idh=[idhx] option=N 
*    exec ntproj [idhx] //lun20/1.a[nca] ! 4000 0 4000
    idhxt=400
    1d [idhxt] ! 4000 0 4000
    nt/pl //lun20/1.a[nca] t[nct]>0 idh=[idhxt] option=N 
*    exec ntproj [idhxt] //lun20/1.a[nca] t[nct]>0 4000 0 4000
    idhx0=401
    1d [idhx0] ! 4000 0 4000
    nt/pl //lun20/1.a[nca] t[nct].eq.0 idh=[idhx0] option=N 
*    exec ntproj [idhx0] //lun20/1.a[nca] t[nct].eq.0 4000 0 4000
    ve/inp ampx([i]) $hinfo([idhx],'mean')
    nx=$hinfo([idhx],'events')
    rms=$hinfo([idhx],'rms')
    ve/inp dampx([i]) $sigma([rms]/sqrt([nx]))
*
    idhxtt=500
    1d [idhxtt] ! 4000 0 4000
    nt/pl //lun20/1.t[nct] t[nct]>0 idh=[idhxtt] option=N
    mean=$hinfo([idhxtt],'mean')
    rms=$hinfo([idhxtt],'rms')
    tl=$sigma([mean]-1.5*sqrt(3.0)*[rms])
    tr=$sigma([mean]+1.5*sqrt(3.0)*[rms])
    nt/pl //lun20/1.t[nct] t[nct]>0.and.[tl]<t[nct].and.t[nct]<[tr] idh=[idhxtt] option=N
    mean=$hinfo([idhxtt],'mean')
    rms=$hinfo([idhxtt],'rms')
    ve/inp tmin([i]) $sigma([mean]-sqrt(3.0)*[rms])
    ve/inp tmax([i]) $sigma([mean]+sqrt(3.0)*[rms])
    close 20
    if ([auto].le.0) then
      l=$sigma(tmin([i])-[rms])
      r=$sigma(tmax([i])+[rms])
      hi/pl [idhxtt]([l]:[r])
      line $sigma(tmin([i])) $GRAFINFO('WNYMIN') $sigma(tmin([i])) $GRAFINFO('WNYMAX')
      line $sigma(tmax([i])) $GRAFINFO('WNYMIN') $sigma(tmax([i])) $GRAFINFO('WNYMAX')
      read x
    endif
*    
    nx=$hinfo([idh00],'xbins')
    xmin=$hinfo([idh00],'xmin')
    xmax=$hinfo([idh00],'xmax')
    ve/cre vx([nx]) r
    hi/get/cont [idh00] vx
    ve/cre ix([nx]) r
    sigma ix = array([nx],1#[nx])
    sigma ix = order(ix,-vx)
    sigma vx = order(vx,-vx)
    mean=$sigma((ix(1)-0.5)*([xmax]-[xmin])/[nx])
    ve/cre p3(3) r $sigma(vx(1)) [mean] 2
    do j=1,[nfit]
      l=$sigma(p3(2)-3*abs(p3(3)))
      r=$sigma(p3(2)+1.5*abs(p3(3)))
      hi/fit [idh00]([l]:[r]) g [fopt] 3 p3
    enddo
    ve/inp pds0f([i]) $sigma(p3(2))
    if ([auto].le.0) then
      l=$sigma(p3(2)-10*abs(p3(3)))
      r=$sigma(p3(2)+10*abs(p3(3)))
      hi/pl [idh00]([l]:[r])
      read x
    endif
*    
    nx=$hinfo([idhx0],'xbins')
    xmin=$hinfo([idhx0],'xmin')
    xmax=$hinfo([idhx0],'xmax')
    ve/cre vx([nx]) r
    hi/get/cont [idhx0] vx
    ve/cre ix([nx]) r
    sigma ix = array([nx],1#[nx])
    sigma ix = order(ix,-vx)
    sigma vx = order(vx,-vx)
    mean=$sigma((ix(1)-0.5)*([xmax]-[xmin])/[nx])
    ve/cre p3(3) r $sigma(vx(1)) [mean] 2
    do j=1,[nfit]
      l=$sigma(p3(2)-3*abs(p3(3)))
      r=$sigma(p3(2)+1.5*abs(p3(3)))
      hi/fit [idhx0]([l]:[r]) g [fopt] 3 p3
    enddo
    ve/inp pdsf([i]) $sigma(p3(2))
    if ([auto].le.0) then
      l=$sigma(p3(2)-10*abs(p3(3)))
      r=$sigma(p3(2)+10*abs(p3(3)))
      hi/pl [idhx0]([l]:[r])
      read x
    endif
*
    n0=$hinfo([idh0],'events')
    n0t=$hinfo([idh0t],'events')
    e0=$sigma([n0t]/[n0])
    de0=$sigma(sqrt([e0]*(1-[e0])/[n0]))
    n1=$hinfo([idhx],'events')
    n1t=$hinfo([idhxt],'events')
    e1=$sigma([n1t]/[n1])
    de1=$sigma(sqrt([e1]*(1-[e1])/[n1]))
    effi=$sigma(1-(1-[e1])/(1-[e0]))
    pe0=$sigma(([de0]/(1-[e0]))**2)
    pe1=$sigma(([de1]/(1-[e1]))**2)
    deffi=$sigma((1-[effi])*sqrt([pe0]+[pe1]))
    mess eff0=[e0] [de0]
    mess effx=[e1] [de1]
    mess effi=[effi] [deffi]
    ve/inp eff([i]) [effi]
    ve/inp deff([i]) [deffi]
    if ([auto].le.0) then
      read x
    endif
*    
  enddo
  sigma pds = pds0 - pds0f + pdsf
  sigma dpds = dpds0
  ve/write u0[nc],pds,dpds,pds0,pds0f,pdsf,tmin,tmax [dir0]/[date]_[cal]/pds_data_counter[nc].txt '(8f15.6)'
enddo
1:
return


macro elr c1=x c2=y vm=paru cm=covu
dxx=[cm](1,1)
dxy=[cm](1,2)
dxz=[cm](1,3)
dyy=[cm](2,2)
dyz=[cm](2,3)
dzz=[cm](3,3)
mess [dxx] [dxy] [dxz] [dyy] [dyz] [dzz]
*
d=$sigma(([d[c1][c1]])*([d[c2][c2]])-([d[c1][c2]])*([d[c1][c2]]))
dxxi=$sigma([d[c2][c2]]/([d]))
dxyi=$sigma(-([d[c1][c2]])/([d]))
dyyi=$sigma([d[c1][c1]]/[d])
mess [c1] [c2]
mess [d] [d[c1][c1]] [d[c1][c2]] [d[c2][c2]]
mess [d] [dxxi] [dxyi] [dyyi]
np=1000
ve/cre xs([np]) r
ve/cre ys([np]) r
do i=1,[np]
  phi=$sigma(([i]-1)*6.28319/([np]-1))
  cphi=$sigma(cos([phi]))
  sphi=$sigma(sin([phi]))
  ri=$sigma([dxxi]*([cphi])**2)
  ri=$sigma([ri]+2*([dxyi])*([cphi])*([sphi]))
  ri=$sigma([ri]+[dyyi]*([sphi])**2)
  ri=$sigma(sqrt([ri]))
  xst=[cphi]/[ri]
  yst=[sphi]/[ri]
  ve/inp xs([i]) [xst]
  ve/inp ys([i]) [yst]
enddo
x0=[vm](1)
y0=[vm](2)
z0=[vm](3)
sigma x = xs+([[c1]0])
sigma y = ys+([[c2]0])
xmax=$sigma(vmax(x))
xmin=$sigma(vmin(x))
ymax=$sigma(vmax(y))
ymin=$sigma(vmin(y))
dx=$sigma([xmax]-[xmin])
dy=$sigma([ymax]-[ymin])
dm=$sigma(max([dx],[dy]))
*dx=[dm]
*dy=[dm]
l=$sigma([xmin]-0.25*[dx])
r=$sigma([xmax]+0.25*[dx])
d=$sigma([ymin]-0.25*[dy])
u=$sigma([ymax]+0.25*[dy])
exec seteps 0
opt ngrid
null [l] [r] [d] [u]
set plci 1
set ltyp 1
graph $vlen(x) x y l
key [[c1]0] [[c2]0] 20 ! 0.1
*
return



macro daycutsf
goto 1
*
p=hv
dir=/work/users/konctbel/AGPARS
ve/del grun,gbeam,gts,gtf
do i=0,30
  fname=[dir]/run_part[i].txt_[p].txt
  if ($fexist([fname])) then
    shell cp [fname] tmp.txt
    shell $unquote('cat tmp.txt | sed "s/nan/0/g" > tmp1.txt')
    shell $unquote('cat tmp1.txt | sed "s/inf/0/g" > tmp.txt')
    fname=tmp.txt
    ve/del run,beam,ts,tf,[p]1,[p]2,[p]3,[p]4,[p]5,[p]6,[p]7,[p]8,[p]9,k[p]1,k[p]2,k[p]3,k[p]4,k[p]5,k[p]6,k[p]7,k[p]8,k[p]9
    ve/read run,beam,ts,tf,[p]1,[p]2,[p]3,[p]4,[p]5,[p]6,[p]7,[p]8,[p]9,k[p]1,k[p]2,k[p]3,k[p]4,k[p]5,k[p]6,k[p]7,k[p]8,k[p]9 [fname]
    exec vappend grun run
    exec vappend gbeam beam
    exec vappend gts ts
    exec vappend gtf tf
  endif
enddo
ve/write grun,gbeam,gts,gtf run_energy_stime_ftime.txt '(4f15.6)'
*
1:
ve/del run,beam,ts,tf
ve/read run,beam,ts,tf run_energy_stime_ftime.txt '(4f15.6)'
nf=$vlen(run)
*
ve/cre id1(20000) r
ve/cre id2(20000) r
dr=$sigma(-500+9.0/24)
ind=0
do i=1,[nf]
  r=run([i])
  ndays=ts([i])
  if ([ndays].gt.[dr]) then
    dr=$sigma(int([ndays])+9.0/24)
    if ([ndays].gt.[dr]) then
      dr=[dr]+1
    endif
    if ([ind].gt.0) then
      ve/inp id2([ind]) [io]
      ind=[ind]+1
      ve/inp id1([ind]) [i]
    else
      ind=[ind]+1
      ve/inp id1([ind]) [i]
    endif
  endif
  io=[i]
enddo
ve/inp id2([ind]) [io]
exec mapcal#vecut id1
exec mapcal#vecut id2
*
nf=$vlen(id1)
ve/cre ts1([nf]) r
ve/cre ts2([nf]) r
ve/cre tf1([nf]) r
ve/cre tf2([nf]) r
ve/cre rn1([nf]) r
ve/cre rn2([nf]) r
do i=1,[nf]
  ind=id1([i])
  ve/inp ts1([i]) ts([ind])
  ve/inp tf1([i]) tf([ind])
  ve/inp rn1([i]) run([ind])
  ind=id2([i])
  ve/inp ts2([i]) ts([ind])
  ve/inp tf2([i]) tf([ind])
  ve/inp rn2([i]) run([ind])
enddo
*
ve/write rn1,rn2,id1,id2,ts1,ts2,tf1,tf2 daycutsfull.txt 8f15.6
return

macro ampcorrx vx=days
do i=1,9
  exec mapcal#ampcorr [i] [vx]
enddo
return


macro ampcorr nc=1 vx=days
set pmci 1
set hcol 1
set fcol 1
set plci 1
*
ve/del r1,r2
ve/read r1,r2 accled_corr.ixlist
ve/del rn1,rn2,id1,id2,ts1,ts2,tf1,tf2
ve/read rn1,rn2,id1,id2,ts1,ts2,tf1,tf2 daycutsfull.txt 8f15.6
ncal=$vlen(r1)
ve/cre runsx(8) r 7842 10988 11285 13846 13847 16698 16699 17868
do i=1,[ncal]
  rb=r1([i])
  re=r2([i])
*  
  exec mapcal#ixndl [rb] rn1
  gl/imp inx
  in1=[inx]
  exec mapcal#ixndl [rb] rn2
  gl/imp inx
  in2=[inx]
  imin=$sigma(min([in1],[in2]))
*  
  exec mapcal#ixndr [re] rn1
  gl/imp inx
  in1=[inx]
  exec mapcal#ixndr [re] rn2
  gl/imp inx
  in2=[inx]
  imax=$sigma(max([in1],[in2]))
*
  rmin=rn1([imin])
  rmax=rn2([imax])
  partmin=$sigma(int([rmin]/1000)) 
  partmax=$sigma(int([rmax]/1000))
*
  mcal=0
  do l=1,4
    j=2*[l]-1
    rx1=$sigma(runsx([j]))
    j=2*[l]
    rx2=$sigma(runsx([j]))
    if ($sigma(max([rb],[rx1])).lt.$sigma(min([re],[rx2]))) then
      mcal=[l]
    endif
  enddo
*
  dir=v2
  ve/del runs,days,ampi,dampi,n0bi,nxbi,effi,deffi
  mess [rb] [re] [rmin] [rmax] [partmin] [partmax]
  do j=[partmin],[partmax]
    fname1=[dir]/amplitude_vs_runs_counter[nc]_p[j].txt
    fname2=[dir]/n0_nx_vs_runs_counter[nc]_p[j].txt
    if (($fexist([fname1]).eq.1).and.($fexist([fname2]).eq.1)) then
      ve/read runs[j],days[j],ampp[j],dampp[j] [fname1] 4f15.6
      ve/read runs[j],days[j],n0b[j],nxb[j] [fname2] 4f15.6
      vrmin=$sigma(vmin(runs[j]))
      vrmax=$sigma(vmax(runs[j]))
      if (([rmin].le.[vrmax]).and.([rmax].ge.[vrmin])) then
        exec mapcal#ixndl [rmin] runs[j]
        gl/imp inx
        in1=[inx]
        exec mapcal#ixndr [rmax] runs[j]
        gl/imp inx
        in2=[inx]
        mess [in1] [in2]
        if ([in1].le.[in2]) then
          ve/del runsp,daysp,amppp,damppp,n0bp,nxbp
          ve/copy runs[j]([in1]:[in2]) runsp
          ve/copy days[j]([in1]:[in2]) daysp
          ve/copy ampp[j]([in1]:[in2]) amppp
          ve/copy dampp[j]([in1]:[in2]) damppp
          ve/copy n0b[j]([in1]:[in2]) n0bp
          ve/copy nxb[j]([in1]:[in2]) nxbp
          exec vappend runs runsp
          exec vappend days daysp
          exec vappend ampi amppp
          exec vappend dampi damppp
          exec vappend n0bi n0bp
          exec vappend nxbi nxbp
        endif
      endif
    endif
  enddo
  ni=$vlen(days)
  if ([ni].ne.0) then
    ve/cre dn([ni]) r
    set pmci 1
    * exec $PER/s#vpl ampi dampi [vx] dn sz=0.1
    sigma effi = 1-max(n0bi,1)/nxbi
    sigma deffi = sqrt(effi*(1-effi)/nxbi)
    sigma ampic = -log(1-effi)
    sigma dampic = deffi/(1-effi)
    set pmci 2
*    * exec $PER/s#vpl ampic dampic [vx] dn sz=0.1 o=s
    set pmci 1
    do j=[imin],[imax]
      rl=rn1([j])
      rr=rn2([j])
      exec mapcal#ixndl [rl] runs
      gl/imp inx
      in1=[inx]
      exec mapcal#ixndr [rr] runs
      gl/imp inx
      in2=[inx]
      mess [in1] [in2]
      if ([in1].lt.[in2]) then
        ve/del runsd,daysd,ampid,dampid,n0bid,nxbid
        ve/copy runs([in1]:[in2]) runsd
        ve/copy days([in1]:[in2]) daysd
        ve/copy ampi([in1]:[in2]) ampid
        ve/copy dampi([in1]:[in2]) dampid
        ve/copy n0bi([in1]:[in2]) n0bid
        ve/copy nxbi([in1]:[in2]) nxbid
        ve/cre p1(1) r 7
        ve/cre dp1(1) r
        set plci 4
        ve/fit daysd ampid dampid p0 s 1 p1 ! ! ! dp1
        set plci 1
        n0bis=$sigma(vsum(n0bid))
        nxbis=$sigma(vsum(nxbid))
        effis=$sigma(1-[n0bis]/[nxbis])
        deffis=$sigma(sqrt([effis]*(1-[effis])/[nxbis]))
        mess [n0bis] [nxbis] [effis] [deffis]
        ve/cre ampics(1) r $sigma(-log(1-[effis]))
        ve/cre dampics(1) r $sigma([deffis]/(1-[effis]))
        ve/cre [vx]s(1) r $sigma((vmax([vx]d)+vmin([vx]d))/2)
        ve/cre d[vx]s(1) r $sigma((vmax([vx]d)-vmin([vx]d))/2)
        set pmci 2
        set hcol 2
        * exec $PER/s#vpl ampics dampics [vx]s d[vx]s sz=0.1 o=s
        set pmci 1
        set hcol 1
      endif
    enddo
    if (([mcal].eq.1).or.([mcal].eq.3)) then
      fname=[dir]/mapcal[mcal]_fit_amp_vs_days_counter[nc].txt
      ve/read p7,dp7 [fname] '2f15.6'
      fun/pl derf0p.f $GRAFINFO('WNXMIN') $GRAFINFO('WNXMAX') s
      ve/read p4,dp4 mapcal[mcal]_sigmu2_fit_days_counter[nc].txt '(2f15.6)'
      ve/inp p4(4) 0
      fun/pl mucorr.f $GRAFINFO('WNXMIN') $GRAFINFO('WNXMAX') s
    endif
    if (([mcal].eq.2).or.([mcal].eq.4)) then
      fname=[dir]/mapcal[mcal]_fit_amp_vs_days_counter[nc].txt
      ve/read p4,dp4 [fname] '2f15.6'
      fun/pl erf0p.f $GRAFINFO('WNXMIN') $GRAFINFO('WNXMAX') s
      ve/read p2,dp2 mapcal[mcal]_sigmu2_fit_days_counter[nc].txt '(2f15.6)'
      fun/pl mucorr2.f $GRAFINFO('WNXMIN') $GRAFINFO('WNXMAX') s
    endif
*    
    if ([vx].eq.'days') then
      atitle 'T, days' 'Signal, pe.'
    endif
    if ([vx].eq.'runs') then
      atitle 'N, runs' 'Signal, pe.'
    endif
    exec save vcal/accled[i]_amp_vs_[vx]_counter[nc].eps f
*    read x
  endif
enddo
return



macro sigmu3 i=1
ve/cre  sm(10) r
ve/cre dsm(10) r
msum=days
ve/del p4s,dp4s
do nc=0,9
  ve/read p4s,dp4s mapcal[i]_sigmu2_fit_[msum]_counter[nc].txt '(2f15.6)'
  ind=[nc]+1
  ve/inp sm([ind]) $sigma(p4s(1))
  ve/inp dsm([ind]) $sigma(dp4s(1))
enddo  
ve/cre vnc(10) r 0 1 2 3 4 5 6 7 8 9
ve/cre dn(10) r
* exec $PER/s#vpl sm dsm vnc dn 
return


macro ampcorrtot
dirc=/work/users/konctbel/AGPARS
ve/del runc,ac1,ac2,ac3,ac4,ac5,ac6,ac7,ac8,ac9
ve/read runc,ac1,ac2,ac3,ac4,ac5,ac6,ac7,ac8,ac9 [dirc]/AmpCorrHV/amplitude_correction_hv.txt.old '(10f15.6)'
do nc=1,9
  ve/del days,aci[nc]
  ve/read days,aci[nc] mapcal_ampcorr_days_counter[nc].txt '(2f15.6)'
enddo
dmin=$sigma(vmin(days)-10)
dt=9/24
sigma daysl=int(days-([dt])-([dmin]))+([dt])+([dmin])
sigma daysr=daysl+1
ve/del runs,beams,tsx,tfx
dir=/work/users/konctbel/AGPARS
p=hv
do part=0,30
  fname=[dir]/run_part[part].txt_[p].txt
  if ($fexist([fname])) then
    shell cp [fname] tmp.txt
    shell $unquote('cat tmp.txt | sed "s/nan/0/g" > tmp1.txt')
    shell $unquote('cat tmp1.txt | sed "s/inf/0/g" > tmp.txt')
    fname=tmp.txt
    if ([p].eq.'hv') then
      ve/del run,beam,ts,tf,[p]1,[p]2,[p]3,[p]4,[p]5,[p]6,[p]7,[p]8,[p]9,k[p]1,k[p]2,k[p]3,k[p]4,k[p]5,k[p]6,k[p]7,k[p]8,k[p]9
      ve/read run,beam,ts,tf,[p]1,[p]2,[p]3,[p]4,[p]5,[p]6,[p]7,[p]8,[p]9,k[p]1,k[p]2,k[p]3,k[p]4,k[p]5,k[p]6,k[p]7,k[p]8,k[p]9 [fname]
      exec vappend runs run
      exec vappend beams beam
      exec vappend tsx ts
      exec vappend tfx tf
    endif
  endif
enddo
do nc=1,9
  ve/del ac[nc]c
  ve/copy ac[nc] ac[nc]c
enddo
rmin=$sigma(min(vmin(runc),vmin(runs)))
rmax=$sigma(max(vmax(runc),vmax(runs)))
n=[rmax]-[rmin]+1
ve/cre ix1([n]) r
ve/cre ix2([n]) r
do i=1,$vlen(runc)
  ind=runc([i])-[rmin]+1
  ve/inp ix1([ind]) [i]
enddo
do i=1,$vlen(runs)
  ind=runs([i])-[rmin]+1
  ve/inp ix2([ind]) [i]
enddo
ind=0
do i=1,[n]
  i1=ix1([i])
  i2=ix2([i])
  if (([i1].ne.0).and.([i2].ne.0)) then
    if ($sigma(runc([i1])).ne.$sigma(runs([i2]))) then
      mess ---> [i] [i1] [i2] ($sigma(runc([i1])) $sigma(runs([i2]))
    endif
    ti=$sigma((tsx([i2])+tfx([i2]))/2)
    exec mapcal#ixndr [ti] daysl
    gl/imp inx
    in1=[inx]
    exec mapcal#ixndl [ti] daysr
    gl/imp inx
    in2=[inx]
    if ([in1].eq.[in2]) then
      do nc=1,9
        ve/inp ac[nc]c([i1]) $sigma(ac[nc]c([i1])*aci[nc]([in1]))
      enddo
    endif
  endif
  if ($sigma(mod([i],100)).eq.0) then
    mess [i]
  endif
enddo
*
nf=$vlen(runc)
do i=1,-[nf]
  runi=runc([i])
  exec mapcal#ixndl [runi] runs
  gl/imp inx
  in0=[inx]
  if ($sigma(runs([in0])).ne.[runi]) then
    mess ---> [runi]
  endif
  ti=$sigma((tsx([in0])+tfx([in0]))/2)
  exec mapcal#ixndr [ti] daysl
  gl/imp inx
  in1=[inx]
  exec mapcal#ixndl [ti] daysr
  gl/imp inx
  in2=[inx]
  if ([in1].eq.[in2]) then
    do nc=1,9
      k=aci[nc]([in1])
      if ([k].le.0) then
        k=1
      endif
      ve/inp ac[nc]c([i]) $sigma(ac[nc]c([i])*[k])
    enddo
  endif
  if ($sigma(mod([i],100)).eq.0) then
    mess [i]
  endif
enddo
ve/write runc,ac1c,ac2c,ac3c,ac4c,ac5c,ac6c,ac7c,ac8c,ac9c ampcorr_full.txt '(10f15.6)'
return


macro ampspecttest nc=1 runi=16500 idh=100
par=dphi
ver=0
ve/del runsp,drunsp,daysp,ddaysp,beamp,dbeamp,nevtp,meanp,dmeanp,rmsp,drmsp,wp,dwp,vp,dvp
ve/read runsp,drunsp,daysp,ddaysp,beamp,dbeamp,nevtp,meanp,dmeanp,rmsp,drmsp,wp,dwp,vp,dvp v2/[par]_vs_beams_new_cut[ver].txt 15e15.6
n=$vlen(runsp)
ve/cre rl([n]) r
ve/cre rr([n]) r
sigma rl=runsp-drunsp
sigma rr=runsp+drunsp
*
exec mapcal#ixndl [runi] rr
gl/imp inx
in1=[inx]
exec mapcal#ixndr [runi] rl
gl/imp inx
in2=[inx]
*
dir=v2
if ([dir].eq.'v1') then
  suffh=.his
  sufft=.txt
else
  suffh=_[dir].his
  sufft=_[dir].txt
endif
*
ind=0
if ([in1].eq.[in2]) then
  mess beam=$sigma(beamp([in1]))
  gl/cre beam $sigma(beamp([in1]))
  rb=rl([in1])
  re=rr([in1])
  do i=[rb],[re]
    fname=[dir]/run_[i]_spects[suffh]
    if ($fexist([fname]).eq.1) then
      hi/file 20 [fname]
      ind=[ind]+1
      if ([ind].eq.1) then
        if ($hexist([idh])) then
          hi/del [idh]
        endif
        hrin 3[nc]
        nx=$hinfo(3[nc],'xbins')
        xmin=$hinfo(3[nc],'xmin')
        xmax=$hinfo(3[nc],'xmax')
        ve/cre nevt([nx]) r
        ve/cre nevti([nx]) r
        hi/get/cont 3[nc] nevti
      else
        hrin 3[nc]
        hi/get/cont 3[nc] nevti
      endif
      sigma nevt=nevt+nevti
      close 20
    endif
  enddo
  1d [idh] ! [nx] [xmin] [xmax]
  hi/put/cont [idh] nevt
endif
return

macro amphistcomptex
par=dphi
ver=0
ve/del runsp,drunsp,daysp,ddaysp,beamp,dbeamp,nevtp,meanp,dmeanp,rmsp,drmsp,wp,dwp,vp,dvp
ve/read runsp,drunsp,daysp,ddaysp,beamp,dbeamp,nevtp,meanp,dmeanp,rmsp,drmsp,wp,dwp,vp,dvp v2/[par]_vs_beams_new_cut[ver].txt 15e15.6
fname=amplitude_spectr_test_0.tex
if ($fexist([fname]).eq.1) then
  shell rm [fname]
endif
for/file  20 [fname] new
close 20
do i=110,134
  beam=$sigma(beamp([i]))
  do nc=1,9
    epsfile=amplitude_spectr_test_e[beam]_counter[nc].eps
    fmess '\begin{figure}[ht!b]' [fname]
    fmess '  \begin{minipage}{\textwidth}' [fname]
    fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
    txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
    fmess [txt] [fname]
    txt=$unquote('     ')\caption{E=[beam] MeV, counter \No [nc]}
    fmess [txt] [fname]
    fmess '   \end{minipage}' [fname]
    fmess '\end{figure}' [fname]
  enddo
  fmess '\clearpage' [fname]
enddo
return


macro amphistcompx
par=dphi
ver=0
ve/del runsp,drunsp,daysp,ddaysp,beamp,dbeamp,nevtp,meanp,dmeanp,rmsp,drmsp,wp,dwp,vp,dvp
ve/read runsp,drunsp,daysp,ddaysp,beamp,dbeamp,nevtp,meanp,dmeanp,rmsp,drmsp,wp,dwp,vp,dvp v2/[par]_vs_beams_new_cut[ver].txt 15e15.6
do i=110,134
  run=$sigma(int(runsp([i])))
  do nc=1,9
    exec mapcal#amphistcomp [nc] [run] 15000
  enddo
enddo
return


macro amphistcomp nc=1 run1=16550 run0=15000
idh1=1111
idh2=2222
idh3=3333
exec mapcal#ampspecttest [nc] [run0] [idh2]
exec mapcal#ampspecttest [nc] [run1] [idh1] 
nx=$hinfo([idh1],'xbins')
xmin=$hinfo([idh1],'xmin')
xmax=$hinfo([idh1],'xmax')
ve/cre nevt1([nx]) r
ve/cre nevt2([nx]) r
ve/cre nevt3([nx]) r
hi/get/cont [idh1] nevt1
hi/get/cont [idh2] nevt2
sigma nevt3 = nevt2/vsum(nevt2)*vsum(nevt1)
1d [idh3] ! [nx] [xmin] [xmax]
hi/put/cont [idh3] nevt3
*exec hsigma @[idh3] = @[idh2]/vsum(@[idh2])*vsum(@[idh1])
exec hsigma u = max(vmax(@[idh3]),vmax(@[idh1]))
*xmin=$hinfo([idh3],'xmin')
*xmax=$hinfo([idh3],'xmax')
null [xmin] [xmax] 0 $sigma(u(1)*1.2)
hi/pl [idh3] s
hi/pl [idh1] se
mess $hinfo([idh1],'mean') $hinfo([idh2],'mean')
*
ve/cre pn(3) r $hinfo([idh1],'xbins') $hinfo([idh1],'xmin') $hinfo([idh1],'xmax')
ve/cre nexp(1) r $hinfo([idh1],'events')
ve/cre hkk(1000) r
hi/get/cont [idh3] hkk
*ve/cre p2(2) r 1 5
*ve/cre dp2(2) r
*ve/cre s2(2) r 0.01 0
*hi/fit [idh1] hifit.f sbn 2 p2 s2 ! ! dp2
*
opt nstat
hi/pl [idh1](-3.:20.) e
hi/pl [idh3] s
*
txt=Counter [nc]
exec pl#tf 0.55 0.9 [txt]
*
gl/imp beam
txt=E?0! = [beam] MeV
exec pl#tf 0.55 0.8 [txt]
*
mean1=$hinfo([idh1],'mean')
rms=$hinfo([idh1],'rms')
nevt=$hinfo([idh1],'events')
sig1=$sigma([rms]/sqrt([nevt]))
txt=[m]?mid! = $sigma(int([mean1]*100+0.5)/100) [\261] $sigma(int([sig1]*100+0.5)/100)
exec pl#tf 0.55 0.7 [txt]
*
mean2=$hinfo([idh2],'mean')
rms=$hinfo([idh2],'rms')
nevt=$hinfo([idh2],'events')
sig2=$sigma([rms]/sqrt([nevt]))
txt=[m]?mid,n! = $sigma(int([mean2]*100+0.5)/100) [\261] $sigma(int([sig2]*100+0.5)/100)
exec pl#tf 0.55 0.6 [txt]
*
*txt=R?fit! = $sigma(int(p2(1)*1000+0.5)/1000) [\261] $sigma(int(dp2(1)*1000+0.5)/1000)
*exec pl#tf 0.55 0.6 [txt]
*
r=$sigma([mean1]/[mean2])
dr=$sigma([r]*sqrt(([sig1]/[mean1])**2+([sig2]/[mean2])**2))
txt=R = $sigma(int([r]*1000+0.5)/1000) [\261] $sigma(int([dr]*1000+0.5)/1000)
exec pl#tf 0.55 0.5 [txt]
*
ve/cre p2(2) r [r] $sigma((1-([r]))*[nx]*(0-([xmin]))/([xmax]-([xmin])))
point 1000
set hcol 2
fun/pl hifitp.f(x) -3 20 s
set hcol 1
*
atitle '[m], pe' 'events'
exec save amplitude_spectr_test_e[beam]_counter[nc].eps f
return

macro calsort
ve/read indc,nndc cal_db.list0 '2f15.0'
n=$sigma(vmax(indc))
ve/cre vcal([n]) r
ve/cre ncal([n]) r
ind=0
do i=1,[n]
  nndcm=-10
  jm=-10
  do j=1,$vlen(indc)
    indcj=indc([j])
    nndcj=nndc([j])
    if (([nndcj].gt.[nndcm]).and.([indcj].eq.[i])) then
      nndcm=nndc([j])
      jm=[j]
    endif
  enddo
  if ([nndcm].lt.10) then
    jm=0
    nndcm=0
  endif
  ind=[ind]+1
  ve/inp vcal([ind]) [jm]
  ve/inp ncal([ind]) [nndcm]
enddo
*exec mapcal#vecut vcal
*exec mapcal#vecut ncal
return


macro calsearchx
n=180
ve/cre indc([n]) r
ve/cre nndc([n]) r
do i=1,[n]
  gl/cre inc -1
  gl/cre nnc -1
  exec mapcal#calsearch $format([i],i4.4) 1
  gl/imp inc
  ve/inp indc([i]) [inc]
  gl/imp nnc
  ve/inp nndc([i]) [nnc]
enddo
ve/write indc,nndc cal_db.list0 '2f15.0'
return


macro calsearch cal=0180 mode=1
*
  ve/del r1,r2
  ve/read r1,r2 accled_corr.ixlist
  q=a
  ve/del v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9
  ve/read v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9 accled_[q]_v1.txt 9f10.3
  q=b
  ve/del v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9
  ve/read v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9 accled_[q]_v1.txt 9f10.3
  q=c
  ve/del v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9
  ve/read v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9 accled_[q]_v1.txt 9f10.3
  q=ae
  ve/del v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9
  ve/read v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9 accled_amp1pe_v1.txt 9f10.3
* 
  if ([mode].eq.0) then
  nva=0
  nvb=0
  nvc=0
  fcal=/work/snd2000/onldat/acc/a_[cal].txt
  if ($fexist([fcal])) then
    nva=1
    shell cp [fcal] tmp.txt
    shell $unquote('cat tmp.txt | sed "s/nan/0/g" > tmp1.txt')
    shell $unquote('cat tmp1.txt | sed "s/inf/0/g" > tmp.txt')
    fcal=tmp.txt
    ve/del vai
    ve/read vai [fcal]
  endif
  fcal=/work/snd2000/onldat/acc/b_[cal].txt
  if ($fexist([fcal])) then
    nvb=1
    shell cp [fcal] tmp.txt
    shell $unquote('cat tmp.txt | sed "s/nan/0/g" > tmp1.txt')
    shell $unquote('cat tmp1.txt | sed "s/inf/0/g" > tmp.txt')
    fcal=tmp.txt
    ve/del vbi
    ve/read vbi [fcal]
  endif
  fcal=/work/snd2000/onldat/acc/c_[cal].txt
  if ($fexist([fcal])) then
    nvc=1
    shell cp [fcal] tmp.txt
    shell $unquote('cat tmp.txt | sed "s/nan/0/g" > tmp1.txt')
    shell $unquote('cat tmp1.txt | sed "s/inf/0/g" > tmp.txt')
    fcal=tmp.txt
    ve/del vci
    ve/read vci [fcal]
  endif
  fcal=/work/snd2000/onldat/acc/amp1pe_[cal].txt
  if ($fexist([fcal])) then
    nvd=1
    shell cp [fcal] tmp.txt
    shell $unquote('cat tmp.txt | sed "s/nan/0/g" > tmp1.txt')
    shell $unquote('cat tmp1.txt | sed "s/inf/0/g" > tmp.txt')
    fcal=tmp.txt
    ve/del vamp1pei
    ve/read vamp1pei [fcal]
  endif
  fcal=/work/snd2000/users/martin/snd2k/testrel/results/a_[cal].txt
  if ($fexist([fcal])) then
    nva=1
    shell cp [fcal] tmp.txt
    shell $unquote('cat tmp.txt | sed "s/nan/0/g" > tmp1.txt')
    shell $unquote('cat tmp1.txt | sed "s/inf/0/g" > tmp.txt')
    fcal=tmp.txt
    ve/del vai
    ve/read vai [fcal]
  endif
  fcal=/work/snd2000/users/martin/snd2k/testrel/results/b_[cal].txt
  if ($fexist([fcal])) then
    nvb=1
    shell cp [fcal] tmp.txt
    shell $unquote('cat tmp.txt | sed "s/nan/0/g" > tmp1.txt')
    shell $unquote('cat tmp1.txt | sed "s/inf/0/g" > tmp.txt')
    fcal=tmp.txt
    ve/del vbi
    ve/read vbi [fcal]
  endif
  fcal=/work/snd2000/users/martin/snd2k/testrel/results/c_[cal].txt
  if ($fexist([fcal])) then
    nvc=1
    shell cp [fcal] tmp.txt
    shell $unquote('cat tmp.txt | sed "s/nan/0/g" > tmp1.txt')
    shell $unquote('cat tmp1.txt | sed "s/inf/0/g" > tmp.txt')
    fcal=tmp.txt
    ve/del vci
    ve/read vci [fcal]
  endif
  fcal=/work/snd2000/users/martin/snd2k/testrel/results/amp1pe_[cal].txt
  if ($fexist([fcal])) then
    nvd=1
    shell cp [fcal] tmp.txt
    shell $unquote('cat tmp.txt | sed "s/nan/0/g" > tmp1.txt')
    shell $unquote('cat tmp1.txt | sed "s/inf/0/g" > tmp.txt')
    fcal=tmp.txt
    ve/del vamp1pei
    ve/read vamp1pei [fcal]
  endif
  fcal=figvam
  fcal0=/work/snd2000/users/martin/snd2k/testrel/results/rez_cal_pmt_[cal].txt
  if ($fexist([fcal0])) then
    fcal=[fcal0]
  endif
  fcal0=/work/snd2000/onldat/acc/rez_cal_pmt_[cal].txt
  if ($fexist([fcal0])) then
    fcal=[fcal0]
  endif
  if ($fexist([fcal])) then
    ve/del vai,vbi,vci
    ve/cre vai(9) r
    ve/cre vbi(9) r
    ve/cre vci(9) r
    shell grep -E "^par3_" [fcal] >& tmp.txt
    shell $unquote('cat tmp.txt | sed "s/nan/1000000/g" > tmp1.txt')
    shell $unquote('cat tmp1.txt | sed "s/inf/1000000/g" > tmp.txt')
    shell $unquote('cat tmp.txt | sed "s/par3_/ /g" > tmp1.txt')
    shell cat tmp0.txt tmp1.txt > tmp2.txt
    ve/del np,vp
    ve/read np,vp tmp2.txt
    if ($vlen(vp).eq.30) then
      nva=1
      nvb=1
      nvc=1
      di=0
      do i=1,9
        i3=$sigma(([i]-[di])*3+1)
        ve/inp vai([i]) $sigma(vp([i3]))
        i3=$sigma(([i]-[di])*3+2)
        ve/inp vbi([i]) $sigma(vp([i3]))
        i3=$sigma(([i]-[di])*3+3)
        ve/inp vci([i]) $sigma(vp([i3]))
      enddo
    endif
  endif
  endif
  if ([mode].eq.1) then
    n=0
    n=[n]+1; v[n]=a           ; vn[n]=[v[n]]
    n=[n]+1; v[n]=b           ; vn[n]=[v[n]]
    n=[n]+1; v[n]=c           ; vn[n]=[v[n]]
    n=[n]+1; v[n]=amp1pe      ; vn[n]=[v[n]]
    do v=1,[n]
      shell cp [v[v]]_calslist.txt tmp.txt
      shell $unquote('cat tmp.txt | sed "s/\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*/       0.000000/g" > tmp1.txt')
      ve/del ixc,[vn[v]]1,[vn[v]]2,[vn[v]]3,[vn[v]]4,[vn[v]]5,[vn[v]]6,[vn[v]]7,[vn[v]]8,[vn[v]]9
      ve/read ixc,[vn[v]]1,[vn[v]]2,[vn[v]]3,[vn[v]]4,[vn[v]]5,[vn[v]]6,[vn[v]]7,[vn[v]]8,[vn[v]]9 tmp1.txt '(10f15.6)'
    enddo
    ve/del vai,vbi,vci,vamp1pei
    ve/cre vai(9) r
    ve/cre vbi(9) r
    ve/cre vci(9) r
    ve/cre vamp1pei(9) r
    if (([cal].ge.$sigma(vmin(ixc))).and.([cal].le.$sigma(vmax(ixc)))) then
    exec mapcal#ixndl [cal] ixc
    gl/imp inx
    in1=[inx]
    exec mapcal#ixndr [cal] ixc
    gl/imp inx
    in2=[inx]
    if ([in1].eq.[in2]) then
      do v=1,[n]
        do nc=1,9
          ve/inp v[v[v]]i([nc]) $sigma([vn[v]][nc]([in1]))
        enddo
      enddo
    endif
    endif
    nva=$sigma(vsum(vai)/vsum(vai))
    nvb=$sigma(vsum(vbi)/vsum(vbi))
    nvc=$sigma(vsum(vci)/vsum(vci))
    ve/copy vamp1pei vxi
    nvd=$sigma(vsum(vxi)/vsum(vxi))
  endif
  nv=[nva]+[nvb]+[nvc]
  if ([nv].eq.3) then
  dp=0.002
  n=$vlen(r1)
  ve/cre iy([n]) r
  ve/cre ix([n]) r
  sigma ix = array([n],1#[n])
  do nc=1,9
*    
    a=$sigma(int(1000*vai([nc]))/1000-[dp])
    exec mapcal#ixndl [a] va[nc]
    gl/imp inx
    ve/inp iy([inx]) $sigma(iy([inx])+1)
    a=$sigma(int(1000*vai([nc]))/1000+[dp])
    exec mapcal#ixndr [a] va[nc]
    gl/imp inx
    ve/inp iy([inx]) $sigma(iy([inx])+1)
*      
    b=$sigma(int(1000*vbi([nc]))/1000-[dp])
    exec mapcal#ixndl [b] vb[nc]
    gl/imp inx
    ve/inp iy([inx]) $sigma(iy([inx])+1)
    b=$sigma(int(1000*vbi([nc]))/1000+[dp])
    exec mapcal#ixndr [b] vb[nc]
    gl/imp inx
    ve/inp iy([inx]) $sigma(iy([inx])+1)
*      
    c=$sigma(int(1000*vci([nc]))/1000-[dp])
    exec mapcal#ixndl [c] vc[nc]
    gl/imp inx
    ve/inp iy([inx]) $sigma(iy([inx])+1)
    c=$sigma(int(1000*vci([nc]))/1000+[dp])
    exec mapcal#ixndr [c] vc[nc]
    gl/imp inx
    ve/inp iy([inx]) $sigma(iy([inx])+1)
*  
  enddo
  sigma ix=order(ix,-iy)
  sigma iy=order(iy,-iy)
  gl/cre inc $sigma(ix(1))
  gl/cre nnc $sigma(iy(1))
  mess [cal] $sigma(ix(1)) $sigma(iy(1))
  endif
return


macro calprepx n1=1 n2=200 mode=i
do i=[n1],[n2]
  if ([mode].eq.'b') then
    fkumac=tmp[i]x.kumac
    if ($fexist([fkumac])) then
      shell rm [fkumac]
    endif
    for/file 20 [fkumac] ! N
    close 20
    txt=exec mapcal#calprep [i]
    fmess [txt] [fkumac]
    shell .testrelease/.mainrelease/Offline/submit.sh -q clusters,180 pawbigX11 -b [fkumac]
  else
    exec mapcal#calprep [i]
  endif
enddo
return



macro calprep ncal=180 
dir1=/work/snd2000/users/martin/snd2k/testrel/results
dir2=/work/snd2000/onldat/acc
dir3=/work/snd2000/users/martin/snd2k/testrel/AccCalibr
dir4=/work/snd2000/users/martin/snd2k/testrel/AccCalibr/*
dir5=/work/snd2000/users/martin/snd2k/testrel/AccCalibr/*/*
dir6=/work/snd2000/users/martin/snd2k/testrel/AccCalibr/*/*/*
cal=$format([ncal],i4.4)
fname=*_[cal]*
shell rm tmp[0-9].txt
shell ls --full-time -R [dir1]/[fname] > tmp1.txt
shell ls --full-time -R [dir2]/[fname] > tmp2.txt
shell ls --full-time -R [dir3]/[fname] > tmp3.txt
shell ls --full-time -R [dir4]/[fname] > tmp4.txt
shell ls --full-time -R [dir5]/[fname] > tmp5.txt
shell ls --full-time -R [dir6]/[fname] > tmp6.txt
shell cat tmp1.txt tmp2.txt tmp3.txt tmp4.txt tmp5.txt tmp6.txt > tmp0.txt
shell grep -E [cal] tmp0.txt > tmp5.txt
shell grep -E "^[-]" tmp5.txt > tmp6.txt
shell $unquote('cat tmp6.txt | sed "s/ d[0-9] / hrenvam /g" > tmp.txt')
shell $unquote('cat tmp.txt | sed "s/\(+0700.*$\)/ /g" > tmp1.txt')
shell $unquote('cat tmp1.txt | sed "s/\(+0600.*$\)/ /g" > tmp2.txt')
shell $unquote('cat tmp2.txt | sed "s/[^0-9]/ /g" > tmp3.txt')
ve/del v1,v2,v3,v4,v5,v6,v7,v8,v9
ve/read v1,v2,v3,v4,v5,v6,v7,v8,v9 tmp3.txt
*  
if ($vexist(v1)) then
  year0=2011
  tmin=100000000
  im=-1
  do i=1,$vlen(v1)
    year=v3([i])   ; txt=$format([year],i4.4)
    month=v4([i])  ; txt=[txt].$format([month],i2.2)
    day=v5([i])    ; txt=[txt].$format([day],i2.2)
    hour=v6([i])   ; txt=[txt].$format([hour],i2.2)
    minute=v7([i]) ; txt=[txt].$format([minute],i2.2)
    second=v8([i]) ; txt=[txt].$format([second],i2.2)
    shell /work/users/konctbel/MinuitTest/gtime [year0].01.01.00.00.00 [txt] gtime.txt
    ve/read gtime gtime.txt
    ti=gtime(1)
    if ([ti].lt.[tmin]) then
      tmin=[ti]
      im=[i]
    endif
  enddo
  year=v3([im])   ; txt=$format([year],i4.4)        
  month=v4([im])  ; txt=[txt]-$format([month],i2.2) 
  day=v5([im])    ; txt=[txt]-$format([day],i2.2)   
  hour=v6([im])   ; txt=[txt] $format([hour],i2.2)  
  minute=v7([im]) ; txt=$unquote([txt]):$format([minute],i2.2)
  second=v8([im]) ; txt=$unquote([txt]):$format([second],i2.2)
  msec=v9([im])
  mess [txt]
  dir0=/work/users/konctbel/Calibr/Cal_$format([year],i4.4)_$format([month],i2.2)_$format([day],i2.2)_[cal]
  mess [dir0]
*
  if ($fexist([dir0]).eq.0) then
    shell mkdir [dir0]
    shell mkdir [dir0]/data
    shell mkdir [dir0]/data_el
    shell mkdir [dir0]/picts
  endif
*  
  logfile=[dir0]/files_origin.txt
  if ($fexist([logfile])) then
    shell rm [logfile]
  endif
  for/file 20 [logfile] ! N
  close 20
  shell cp tmp.txt [dir0]/full_files_list.txt
*
  fname=[dir0]/time.txt
  if ($fexist([fname])) then
    shell rm [fname]
  endif
  for/file 20 [fname] ! N
  close 20
  txt=$format([year],i4.4) $format([month],i2.2) $format([day],i2.2) $format([hour],i2.2) $format([minute],i2.2) $format([second],i2.2) $format([msec],i9.9)
  fmess [txt] [fname]
*
  prefix=cal_pmt
  fname=[prefix]_[cal]*.vec
  shell rm tmp[0-9].txt
  shell ls --full-time -R [dir1]/[fname] > tmp1.txt
  shell ls --full-time -R [dir2]/[fname] > tmp2.txt
  shell ls --full-time -R [dir3]/[fname] > tmp3.txt
  shell ls --full-time -R [dir4]/[fname] > tmp4.txt
  shell ls --full-time -R [dir5]/[fname] > tmp5.txt
  shell ls --full-time -R [dir6]/[fname] > tmp6.txt
  shell cat tmp1.txt tmp2.txt tmp3.txt tmp4.txt tmp5.txt tmp6.txt > tmp0.txt
  shell grep -E [cal] tmp0.txt > tmp6.txt
  shell $unquote('cat tmp6.txt | sed "s/ d[0-9] / hrenvam /g" > tmp5.txt')
  shell $unquote('cat tmp5.txt | sed "s/\(^.*+0600\)/ /g" > tmp6.txt')
  shell $unquote('cat tmp6.txt | sed "s/\(^.*+0700\)/ /g" > tmp7.txt')
  shell $unquote('cat tmp5.txt | sed "s/\(^.*\/\)/ /g" > tmp8.txt')
  shell $unquote('cat tmp8.txt | sed "s/[^0-9]/ /g" > tmp9.txt')
  shell cat [logfile] > tmpx.txt
  shell cat tmpx.txt tmp7.txt > [logfile]
  ve/del xtmp,c0
  ve/read xtmp,c0 tmp9.txt
  nc0=$vlen(c0)
  if ([nc0].ne.0) then
    ve/cre c0r([nc0]) r
    ve/del c0s
    sigma c0s = order(c0,c0)
    co=-1
    ind=0
    do i=1,[nc0]
      c0i=c0s([i])
      if ([co].ne.[c0i]) then
        ind=[ind]+1
        ve/inp c0r([ind]) [c0i]
        co=[c0i]
      endif
    enddo
    exec mapcal#vecut c0r
    ve/write c0r [dir0]/c0.txt
  endif
  s='"'
  l='|'
  t='`'
  do i=1,[nc0]
    fkumac=tmp[cal].kumac
    if ($fexist([fkumac])) then
      shell rm [fkumac]
    endif
    for/file 20 [fkumac] ! N
    close 20
    txt=echo fname=[t]cat tmp7.txt [l] sed -n [s][i]p[s][t] > [fkumac]
    shell $unquote([txt])
    txt=tname=[dir0]/data/[prefix]_[cal]_$sigma(c0([i])).tup
    fmess [txt] [fkumac]
    fmess 'nt/cre 1 ! 20 ! ! a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 t1 t2 t3 t4 t5 t6 t7 t8 t9 t10' [fkumac]
    fmess 'nt/read 1 [fname]' [fkumac]
    fmess 'hi/file 20 [tname] ! N' [fkumac]
    fmess 'hrout 1' [fkumac]
    fmess 'close 20' [fkumac]
    shell cat [fkumac]
    shell pawbigX11 -b [fkumac]
  enddo
*
  prefix=cal_pmt_p
  fname=[prefix]_[cal]*.vec
  shell rm tmp[0-9].txt
  shell ls --full-time -R [dir1]/[fname] > tmp1.txt
  shell ls --full-time -R [dir2]/[fname] > tmp2.txt
  shell ls --full-time -R [dir3]/[fname] > tmp3.txt
  shell ls --full-time -R [dir4]/[fname] > tmp4.txt
  shell ls --full-time -R [dir5]/[fname] > tmp5.txt
  shell ls --full-time -R [dir6]/[fname] > tmp6.txt
  shell cat tmp1.txt tmp2.txt tmp3.txt tmp4.txt tmp5.txt tmp6.txt > tmp0.txt
  shell grep -E [cal] tmp0.txt > tmp6.txt
  shell $unquote('cat tmp6.txt | sed "s/ d[0-9] / hrenvam /g" > tmp5.txt')
  shell $unquote('cat tmp5.txt | sed "s/\(^.*+0600\)/ /g" > tmp6.txt')
  shell $unquote('cat tmp6.txt | sed "s/\(^.*+0700\)/ /g" > tmp7.txt')
  shell $unquote('cat tmp5.txt | sed "s/\(^.*\/\)/ /g" > tmp8.txt')
  shell $unquote('cat tmp8.txt | sed "s/[^0-9]/ /g" > tmp9.txt')
  shell cat [logfile] > tmpx.txt
  shell cat tmpx.txt tmp7.txt > [logfile]
  ve/del xtmp,c0
  ve/read xtmp,c0 tmp9.txt
  nc0=$vlen(c0)
  if ([nc0].ne.0) then
    ve/cre c0r([nc0]) r
    ve/del c0s
    sigma c0s = order(c0,c0)
    co=-1
    ind=0
    do i=1,[nc0]
      c0i=c0s([i])
      if ([co].ne.[c0i]) then
        ind=[ind]+1
        ve/inp c0r([ind]) [c0i]
        co=[c0i]
      endif
    enddo
    exec mapcal#vecut c0r
    ve/write c0r [dir0]/c0p.txt
  endif
  s='"'
  l='|'
  t='`'
  do i=1,[nc0]
    fkumac=tmp[cal].kumac
    if ($fexist([fkumac])) then
      shell rm [fkumac]
    endif
    for/file 20 [fkumac] ! N
    close 20
    txt=echo fname=[t]cat tmp7.txt [l] sed -n [s][i]p[s][t] > [fkumac]
    shell $unquote([txt])
    txt=tname=[dir0]/data/[prefix]_[cal]_$sigma(c0([i])).tup
    fmess [txt] [fkumac]
    fmess 'nt/cre 1 ! 20 ! ! a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 t1 t2 t3 t4 t5 t6 t7 t8 t9 t10' [fkumac]
    fmess 'nt/read 1 [fname]' [fkumac]
    fmess 'hi/file 20 [tname] ! N' [fkumac]
    fmess 'hrout 1' [fkumac]
    fmess 'close 20' [fkumac]
    shell cat [fkumac]
    shell pawbigX11 -b [fkumac]
  enddo
*
  prefix=cal_el
  fname=[prefix]_[cal]*.vec
  shell rm tmp[0-9].txt
  shell ls --full-time -R [dir1]/[fname] > tmp1.txt
  shell ls --full-time -R [dir2]/[fname] > tmp2.txt
  shell ls --full-time -R [dir3]/[fname] > tmp3.txt
  shell ls --full-time -R [dir4]/[fname] > tmp4.txt
  shell ls --full-time -R [dir5]/[fname] > tmp5.txt
  shell ls --full-time -R [dir6]/[fname] > tmp6.txt
  shell cat tmp1.txt tmp2.txt tmp3.txt tmp4.txt tmp5.txt tmp6.txt > tmp0.txt
  shell grep -E [cal] tmp0.txt > tmp6.txt
  shell $unquote('cat tmp6.txt | sed "s/ d[0-9] / hrenvam /g" > tmp5.txt')
  shell $unquote('cat tmp5.txt | sed "s/\(^.*+0600\)/ /g" > tmp6.txt')
  shell $unquote('cat tmp6.txt | sed "s/\(^.*+0700\)/ /g" > tmp7.txt')
  shell $unquote('cat tmp5.txt | sed "s/\(^.*\/\)/ /g" > tmp8.txt')
  shell $unquote('cat tmp8.txt | sed "s/[^0-9]/ /g" > tmp9.txt')
  shell cat [logfile] > tmpx.txt
  shell cat tmpx.txt tmp7.txt > [logfile]
  ve/del xtmp,g0,n0
  ve/read xtmp,g0,n0 tmp9.txt
  ng0=$vlen(g0)
  if ([ng0].ne.0) then
    ve/cre n0r([ng0]) r
    ve/del n0s
    sigma n0s = order(n0,n0)
    no=-1
    ind=0
    do i=1,[ng0]
      n0i=n0s([i])
      if ([no].ne.[n0i]) then
        ind=[ind]+1
        ve/inp n0r([ind]) [n0i]
        no=[n0i]
      endif
    enddo
    exec mapcal#vecut n0r
    ve/write n0r [dir0]/c0e.txt
  endif
  s='"'
  l='|'
  t='`'
  do i=1,[ng0]
    fkumac=tmp[cal].kumac
    if ($fexist([fkumac])) then
      shell rm [fkumac]
    endif
    for/file 20 [fkumac] ! N
    close 20
    txt=echo fname=[t]cat tmp7.txt [l] sed -n [s][i]p[s][t] > [fkumac]
    shell $unquote([txt])
    txt=tname=[dir0]/data_el/[prefix]_[cal]_$sigma(g0([i]))_$sigma(n0([i])).tup
    fmess [txt] [fkumac]
    fmess 'nt/cre 1 ! 20 ! ! a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 t1 t2 t3 t4 t5 t6 t7 t8 t9 t10' [fkumac]
    fmess 'nt/read 1 [fname]' [fkumac]
    fmess 'hi/file 20 [tname] ! N' [fkumac]
    fmess 'hrout 1' [fkumac]
    fmess 'close 20' [fkumac]
    shell cat [fkumac]
    shell pawbigX11 -b [fkumac]
  enddo
*
  prefix=cal_el_all
  fname=[prefix]_[cal]*.vec
  shell rm tmp[0-9].txt
  shell ls --full-time -R [dir1]/[fname] > tmp1.txt
  shell ls --full-time -R [dir2]/[fname] > tmp2.txt
  shell ls --full-time -R [dir3]/[fname] > tmp3.txt
  shell ls --full-time -R [dir4]/[fname] > tmp4.txt
  shell ls --full-time -R [dir5]/[fname] > tmp5.txt
  shell ls --full-time -R [dir6]/[fname] > tmp6.txt
  shell cat tmp1.txt tmp2.txt tmp3.txt tmp4.txt tmp5.txt tmp6.txt > tmp0.txt
  shell grep -E [cal] tmp0.txt > tmp6.txt
  shell $unquote('cat tmp6.txt | sed "s/ d[0-9] / hrenvam /g" > tmp5.txt')
  shell $unquote('cat tmp5.txt | sed "s/\(^.*+0600\)/ /g" > tmp6.txt')
  shell $unquote('cat tmp6.txt | sed "s/\(^.*+0700\)/ /g" > tmp7.txt')
  shell $unquote('cat tmp5.txt | sed "s/\(^.*\/\)/ /g" > tmp8.txt')
  shell $unquote('cat tmp8.txt | sed "s/[^0-9]/ /g" > tmp9.txt')
  shell cat [logfile] > tmpx.txt
  shell cat tmpx.txt tmp7.txt > [logfile]
  ve/del xtmp,g0
  ve/read xtmp,g0 tmp9.txt
  ng0=$vlen(g0)
  s='"'
  l='|'
  t='`'
  do i=1,[ng0]
    fkumac=tmp[cal].kumac
    if ($fexist([fkumac])) then
      shell rm [fkumac]
    endif
    for/file 20 [fkumac] ! N
    close 20
    txt=echo fname=[t]cat tmp7.txt [l] sed -n [s][i]p[s][t] > [fkumac]
    shell $unquote([txt])
    txt=tname=[dir0]/data_el/[prefix]_[cal]_$sigma(g0([i])).tup
    fmess [txt] [fkumac]
    fmess 'nt/cre 1 ! 20 ! ! a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 t1 t2 t3 t4 t5 t6 t7 t8 t9 t10' [fkumac]
    fmess 'nt/read 1 [fname]' [fkumac]
    fmess 'hi/file 20 [tname] ! N' [fkumac]
    fmess 'hrout 1' [fkumac]
    fmess 'close 20' [fkumac]
    shell cat [fkumac]
    shell pawbigX11 -b [fkumac]
  enddo
*
  prefix=cal_pmt
  fname=[prefix]_[cal]*.vec.bz2
  shell rm tmp[0-9].txt
  shell ls --full-time -R [dir1]/[fname] > tmp1.txt
  shell ls --full-time -R [dir2]/[fname] > tmp2.txt
  shell ls --full-time -R [dir3]/[fname] > tmp3.txt
  shell ls --full-time -R [dir4]/[fname] > tmp4.txt
  shell ls --full-time -R [dir5]/[fname] > tmp5.txt
  shell ls --full-time -R [dir6]/[fname] > tmp6.txt
  shell cat tmp1.txt tmp2.txt tmp3.txt tmp4.txt tmp5.txt tmp6.txt > tmp0.txt
  shell grep -E [cal] tmp0.txt > tmp6.txt
  shell $unquote('cat tmp6.txt | sed "s/ d[0-9] / hrenvam /g" > tmp5.txt')
  shell $unquote('cat tmp5.txt | sed "s/\(^.*+0600\)/ /g" > tmp6.txt')
  shell $unquote('cat tmp6.txt | sed "s/\(^.*+0700\)/ /g" > tmp7.txt')
  shell $unquote('cat tmp5.txt | sed "s/\(^.*\/\)/ /g" > tmp8.txt')
  shell $unquote('cat tmp8.txt | sed "s/[^0-9]/ /g" > tmp9.txt')
  shell cat [logfile] > tmpx.txt
  shell cat tmpx.txt tmp7.txt > [logfile]
  ve/del xtmp,c0,xt
  ve/read xtmp,c0,xt tmp9.txt
  nc0=$vlen(c0)
  if ([nc0].ne.0) then
    ve/cre c0r([nc0]) r
    ve/del c0s
    sigma c0s = order(c0,c0)
    co=-1
    ind=0
    do i=1,[nc0]
      c0i=c0s([i])
      if ([co].ne.[c0i]) then
        ind=[ind]+1
        ve/inp c0r([ind]) [c0i]
        co=[c0i]
      endif
    enddo
    exec mapcal#vecut c0r
    ve/write c0r [dir0]/c0.txt
  endif
  s='"'
  l='|'
  t='`'
  do i=1,[nc0]
    fkumac=tmp[cal].kumac
    if ($fexist([fkumac])) then
      shell rm [fkumac]
    endif
    for/file 20 [fkumac] ! N
    close 20
    txt=echo fname=[t]cat tmp7.txt [l] sed -n [s][i]p[s][t] > [fkumac]
    shell $unquote([txt])
    fmess 'shell bunzip2 -d -c [fname] > tmp.vec' [fkumac]
    txt=fname=tmp.vec
    fmess [txt] [fkumac]
    txt=tname=[dir0]/data/[prefix]_[cal]_$sigma(c0([i])).tup
    fmess [txt] [fkumac]
    fmess 'nt/cre 1 ! 20 ! ! a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 t1 t2 t3 t4 t5 t6 t7 t8 t9 t10' [fkumac]
    fmess 'nt/read 1 [fname]' [fkumac]
    fmess 'hi/file 20 [tname] ! N' [fkumac]
    fmess 'hrout 1' [fkumac]
    fmess 'close 20' [fkumac]
    shell cat [fkumac]
    shell pawbigX11 -b [fkumac]
  enddo
*
  prefix=cal_pmt_p
  fname=[prefix]_[cal]*.vec.bz2
  shell rm tmp[0-9].txt
  shell ls --full-time -R [dir1]/[fname] > tmp1.txt
  shell ls --full-time -R [dir2]/[fname] > tmp2.txt
  shell ls --full-time -R [dir3]/[fname] > tmp3.txt
  shell ls --full-time -R [dir4]/[fname] > tmp4.txt
  shell ls --full-time -R [dir5]/[fname] > tmp5.txt
  shell ls --full-time -R [dir6]/[fname] > tmp6.txt
  shell cat tmp1.txt tmp2.txt tmp3.txt tmp4.txt tmp5.txt tmp6.txt > tmp0.txt
  shell grep -E [cal] tmp0.txt > tmp6.txt
  shell $unquote('cat tmp6.txt | sed "s/ d[0-9] / hrenvam /g" > tmp5.txt')
  shell $unquote('cat tmp5.txt | sed "s/\(^.*+0600\)/ /g" > tmp6.txt')
  shell $unquote('cat tmp6.txt | sed "s/\(^.*+0700\)/ /g" > tmp7.txt')
  shell $unquote('cat tmp5.txt | sed "s/\(^.*\/\)/ /g" > tmp8.txt')
  shell $unquote('cat tmp8.txt | sed "s/[^0-9]/ /g" > tmp9.txt')
  shell cat [logfile] > tmpx.txt
  shell cat tmpx.txt tmp7.txt > [logfile]
  ve/del xtmp,c0,xt
  ve/read xtmp,c0,xt tmp9.txt
  nc0=$vlen(c0)
  if ([nc0].ne.0) then
    ve/cre c0r([nc0]) r
    ve/del c0s
    sigma c0s = order(c0,c0)
    co=-1
    ind=0
    do i=1,[nc0]
      c0i=c0s([i])
      if ([co].ne.[c0i]) then
        ind=[ind]+1
        ve/inp c0r([ind]) [c0i]
        co=[c0i]
      endif
    enddo
    exec mapcal#vecut c0r
    ve/write c0r [dir0]/c0p.txt
  endif
  s='"'
  l='|'
  t='`'
  do i=1,[nc0]
    fkumac=tmp[cal].kumac
    if ($fexist([fkumac])) then
      shell rm [fkumac]
    endif
    for/file 20 [fkumac] ! N
    close 20
    txt=echo fname=[t]cat tmp7.txt [l] sed -n [s][i]p[s][t] > [fkumac]
    shell $unquote([txt])
    fmess 'shell bunzip2 -d -c [fname] > tmp.vec' [fkumac]
    txt=fname=tmp.vec
    fmess [txt] [fkumac]
    txt=tname=[dir0]/data/[prefix]_[cal]_$sigma(c0([i])).tup
    fmess [txt] [fkumac]
    fmess 'nt/cre 1 ! 20 ! ! a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 t1 t2 t3 t4 t5 t6 t7 t8 t9 t10' [fkumac]
    fmess 'nt/read 1 [fname]' [fkumac]
    fmess 'hi/file 20 [tname] ! N' [fkumac]
    fmess 'hrout 1' [fkumac]
    fmess 'close 20' [fkumac]
    shell cat [fkumac]
    shell pawbigX11 -b [fkumac]
  enddo
*
  prefix=cal_el
  fname=[prefix]_[cal]*.vec.bz2
  shell rm tmp[0-9].txt
  shell ls --full-time -R [dir1]/[fname] > tmp1.txt
  shell ls --full-time -R [dir2]/[fname] > tmp2.txt
  shell ls --full-time -R [dir3]/[fname] > tmp3.txt
  shell ls --full-time -R [dir4]/[fname] > tmp4.txt
  shell ls --full-time -R [dir5]/[fname] > tmp5.txt
  shell ls --full-time -R [dir6]/[fname] > tmp6.txt
  shell cat tmp1.txt tmp2.txt tmp3.txt tmp4.txt tmp5.txt tmp6.txt > tmp0.txt
  shell grep -E [cal] tmp0.txt > tmp6.txt
  shell $unquote('cat tmp6.txt | sed "s/ d[0-9] / hrenvam /g" > tmp5.txt')
  shell $unquote('cat tmp5.txt | sed "s/\(^.*+0600\)/ /g" > tmp6.txt')
  shell $unquote('cat tmp6.txt | sed "s/\(^.*+0700\)/ /g" > tmp7.txt')
  shell $unquote('cat tmp5.txt | sed "s/\(^.*\/\)/ /g" > tmp8.txt')
  shell $unquote('cat tmp8.txt | sed "s/[^0-9]/ /g" > tmp9.txt')
  shell cat [logfile] > tmpx.txt
  shell cat tmpx.txt tmp7.txt > [logfile]
  ve/del xtmp,g0,n0,xt
  ve/read xtmp,g0,n0,xt tmp9.txt
  ng0=$vlen(g0)
  if ([ng0].ne.0) then
    ve/cre n0r([ng0]) r
    ve/del n0s
    sigma n0s = order(n0,n0)
    no=-1
    ind=0
    do i=1,[ng0]
      n0i=n0s([i])
      if ([no].ne.[n0i]) then
        ind=[ind]+1
        ve/inp n0r([ind]) [n0i]
        no=[n0i]
      endif
    enddo
    exec mapcal#vecut n0r
    ve/write n0r [dir0]/c0e.txt
  endif
  s='"'
  l='|'
  t='`'
  do i=1,[ng0]
    fkumac=tmp[cal].kumac
    if ($fexist([fkumac])) then
      shell rm [fkumac]
    endif
    for/file 20 [fkumac] ! N
    close 20
    txt=echo fname=[t]cat tmp7.txt [l] sed -n [s][i]p[s][t] > [fkumac]
    shell $unquote([txt])
    fmess 'shell bunzip2 -d -c [fname] > tmp.vec' [fkumac]
    txt=fname=tmp.vec
    fmess [txt] [fkumac]
    txt=tname=[dir0]/data_el/[prefix]_[cal]_$sigma(g0([i]))_$sigma(n0([i])).tup
    fmess [txt] [fkumac]
    fmess 'nt/cre 1 ! 20 ! ! a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 t1 t2 t3 t4 t5 t6 t7 t8 t9 t10' [fkumac]
    fmess 'nt/read 1 [fname]' [fkumac]
    fmess 'hi/file 20 [tname] ! N' [fkumac]
    fmess 'hrout 1' [fkumac]
    fmess 'close 20' [fkumac]
    shell cat [fkumac]
    shell pawbigX11 -b [fkumac]
  enddo
*
  prefix=cal_el_all
  fname=[prefix]_[cal]*.vec.bz2
  shell rm tmp[0-9].txt
  shell ls --full-time -R [dir1]/[fname] > tmp1.txt
  shell ls --full-time -R [dir2]/[fname] > tmp2.txt
  shell ls --full-time -R [dir3]/[fname] > tmp3.txt
  shell ls --full-time -R [dir4]/[fname] > tmp4.txt
  shell ls --full-time -R [dir5]/[fname] > tmp5.txt
  shell ls --full-time -R [dir6]/[fname] > tmp6.txt
  shell cat tmp1.txt tmp2.txt tmp3.txt tmp4.txt tmp5.txt tmp6.txt > tmp0.txt
  shell grep -E [cal] tmp0.txt > tmp6.txt
  shell $unquote('cat tmp6.txt | sed "s/ d[0-9] / hrenvam /g" > tmp5.txt')
  shell $unquote('cat tmp5.txt | sed "s/\(^.*+0600\)/ /g" > tmp6.txt')
  shell $unquote('cat tmp6.txt | sed "s/\(^.*+0700\)/ /g" > tmp7.txt')
  shell $unquote('cat tmp5.txt | sed "s/\(^.*\/\)/ /g" > tmp8.txt')
  shell $unquote('cat tmp8.txt | sed "s/[^0-9]/ /g" > tmp9.txt')
  shell cat [logfile] > tmpx.txt
  shell cat tmpx.txt tmp7.txt > [logfile]
  ve/del xtmp,g0
  ve/read xtmp,g0 tmp9.txt
  ng0=$vlen(g0)
  s='"'
  l='|'
  t='`'
  do i=1,[ng0]
    fkumac=tmp[cal].kumac
    if ($fexist([fkumac])) then
      shell rm [fkumac]
    endif
    for/file 20 [fkumac] ! N
    close 20
    txt=echo fname=[t]cat tmp7.txt [l] sed -n [s][i]p[s][t] > [fkumac]
    shell $unquote([txt])
    fmess 'shell bunzip2 -d -c [fname] > tmp.vec' [fkumac]
    txt=fname=tmp.vec
    fmess [txt] [fkumac]
    txt=tname=[dir0]/data_el/[prefix]_[cal]_$sigma(g0([i])).tup
    fmess [txt] [fkumac]
    fmess 'nt/cre 1 ! 20 ! ! a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 t1 t2 t3 t4 t5 t6 t7 t8 t9 t10' [fkumac]
    fmess 'nt/read 1 [fname]' [fkumac]
    fmess 'hi/file 20 [tname] ! N' [fkumac]
    fmess 'hrout 1' [fkumac]
    fmess 'close 20' [fkumac]
    shell cat [fkumac]
    shell pawbigX11 -b [fkumac]
  enddo
*
  prefix=Uset
  fname=[prefix]_[cal]*
  shell rm tmp[0-9].txt
  shell ls --full-time -R [dir1]/[fname] > tmp1.txt
  shell ls --full-time -R [dir2]/[fname] > tmp2.txt
  shell ls --full-time -R [dir3]/[fname] > tmp3.txt
  shell ls --full-time -R [dir4]/[fname] > tmp4.txt
  shell ls --full-time -R [dir5]/[fname] > tmp5.txt
  shell ls --full-time -R [dir6]/[fname] > tmp6.txt
  shell cat tmp1.txt tmp2.txt tmp3.txt tmp4.txt tmp5.txt tmp6.txt > tmp0.txt
  shell grep -E [cal] tmp0.txt > tmp6.txt
  shell $unquote('cat tmp6.txt | sed "s/ d[0-9] / hrenvam /g" > tmp5.txt')
  shell $unquote('cat tmp5.txt | sed "s/\(^.*+0600\)/ /g" > tmp6.txt')
  shell $unquote('cat tmp6.txt | sed "s/\(^.*+0700\)/ /g" > tmp7.txt')
  shell $unquote('cat tmp5.txt | sed "s/\(^.*\/\)/ /g" > tmp8.txt')
  shell $unquote('cat tmp8.txt | sed "s/[^0-9]/ /g" > tmp9.txt')
  shell cat [logfile] > tmpx.txt
  shell cat tmpx.txt tmp7.txt > [logfile]
  ve/del xtmp
  ve/read xtmp tmp9.txt
  nu=$vlen(xtmp)
  s='"'
  l='|'
  t='`'
  do i=1,[nu]
    txt=cp -v [t]cat tmp7.txt [l] sed -n [s][i]p[s][t] [dir0]
    shell $unquote([txt])
  enddo
*
  prefix=code
  fname=[prefix]_[cal]*
  shell rm tmp[0-9].txt
  shell ls --full-time -R [dir1]/[fname] > tmp1.txt
  shell ls --full-time -R [dir2]/[fname] > tmp2.txt
  shell ls --full-time -R [dir3]/[fname] > tmp3.txt
  shell ls --full-time -R [dir4]/[fname] > tmp4.txt
  shell ls --full-time -R [dir5]/[fname] > tmp5.txt
  shell ls --full-time -R [dir6]/[fname] > tmp6.txt
  shell cat tmp1.txt tmp2.txt tmp3.txt tmp4.txt tmp5.txt tmp6.txt > tmp0.txt
  shell grep -E [cal] tmp0.txt > tmp6.txt
  shell $unquote('cat tmp6.txt | sed "s/ d[0-9] / hrenvam /g" > tmp5.txt')
  shell $unquote('cat tmp5.txt | sed "s/\(^.*+0600\)/ /g" > tmp6.txt')
  shell $unquote('cat tmp6.txt | sed "s/\(^.*+0700\)/ /g" > tmp7.txt')
  shell $unquote('cat tmp5.txt | sed "s/\(^.*\/\)/ /g" > tmp8.txt')
  shell $unquote('cat tmp8.txt | sed "s/[^0-9]/ /g" > tmp9.txt')
  shell cat [logfile] > tmpx.txt
  shell cat tmpx.txt tmp7.txt > [logfile]
  ve/del xtmp
  ve/read xtmp tmp9.txt
  nu=$vlen(xtmp)
  s='"'
  l='|'
  t='`'
  do i=1,[nu]
    txt=cp -v [t]cat tmp7.txt [l] sed -n [s][i]p[s][t] [dir0]
    shell $unquote([txt])
  enddo
endif
return


macro abcgetx
nmax=180
n=0
n=[n]+1; v[n]=a           ; vn[n]=[v[n]]
n=[n]+1; v[n]=b           ; vn[n]=[v[n]]
n=[n]+1; v[n]=c           ; vn[n]=[v[n]]
n=[n]+1; v[n]=amp1pe      ; vn[n]=[v[n]]
n=[n]+1; v[n]=rez_cal_pmt ; vn[n]=rez
goto 1
do v=1,[n]
  do nc=1,9
    ve/cre [vn[v]][nc]([nmax]) r
  enddo
enddo
do i=1,[nmax]
  do v=1,[n]
    exec mapcal#abcget [v[v]] [i]
    gl/imp nfc
    if ([nfc].ne.0) then
      do nc=1,9
        if ([v].lt.[n]) then
          ve/inp [vn[v]][nc]([i]) $sigma([v[v]]([nc]))
        else
          ve/inp [vn1][nc]([i]) $sigma(vai([nc]))
          ve/inp [vn2][nc]([i]) $sigma(vbi([nc]))
          ve/inp [vn3][nc]([i]) $sigma(vci([nc]))
        endif
      enddo
    endif
  enddo
enddo
1:
ve/cre ix([nmax]) r
sigma ix=array([nmax],1#[nmax])
do v=1,[n]
  ve/write ix,[vn[v]]1,[vn[v]]2,[vn[v]]3,[vn[v]]4,[vn[v]]5,[vn[v]]6,[vn[v]]7,[vn[v]]8,[vn[v]]9 [v[v]]_calslist.txt '(10f15.6)'
enddo
return


macro abcget v=a ncal=180
shell rm tmp[0-9].txt
shell ls -d /work/users/konctbel/Calibr/Cal_* > tmp0.txt
shell $unquote('cat tmp0.txt | sed "s/[^0-9]/ /g" > tmp1.txt')
ve/del y,m,d,c
ve/read y,m,d,c tmp1.txt
exec mapcal#ixndl [ncal] c
gl/imp inx
in1=[inx]
exec mapcal#ixndr [ncal] c
gl/imp inx
in2=[inx]
if ([in1].eq.[in2]) then
  year=y([in1])
  month=m([in1])
  day=d([in1])
  cal=$format([ncal],i4.4)
  dir=/work/users/konctbel/Calibr/Cal_$format([year],i4.4)_$format([month],i2.2)_$format([day],i2.2)_[cal]
  fname=[dir]/full_files_list.txt
  file=[v]_[cal].txt
  shell grep -e [file] [fname] > tmp2.txt
  shell cat tmp2.txt
  shell $unquote('cat tmp2.txt | sed "s/\(^.*+0600\)/ /g" > tmp3.txt')
  shell $unquote('cat tmp3.txt | sed "s/\(^.*+0700\)/ /g" > tmp4.txt')
  shell $unquote('cat tmp4.txt | sed "s/\(^.*\/\)/ /g" > tmp5.txt')
  shell $unquote('cat tmp5.txt | sed "s/amp1pe/ /g" > tmp6.txt')
  shell $unquote('cat tmp6.txt | sed "s/[^0-9]/ /g" > tmp7.txt')
  shell cat tmp4.txt
  ve/del cx
  ve/read cx tmp7.txt
  nf=$vlen(cx)
  gl/cre nfc [nf]
  s='"'
  l='|'
  t='`'
  do i=1,[nf]
    fkumac=tmp[cal].kumac
    if ($fexist([fkumac])) then
      shell rm [fkumac]
    endif
    for/file 20 [fkumac] ! N
    close 20
    txt=cp [t]cat tmp4.txt [l] sed -n [s][i]p[s][t] tmp8.txt
    shell $unquote([txt])
    shell $unquote('cat tmp8.txt | sed "s/nan/0/g" > tmp9.txt')
    shell $unquote('cat tmp9.txt | sed "s/inf/0/g" > tmp_cal.txt')
    if ([v].ne.'rez_cal_pmt') then
      ve/del [v]
      ve/read [v] tmp_cal.txt
      ve/prin [v]
    else
      ve/del vai,vbi,vci
      ve/cre vai(9) r
      ve/cre vbi(9) r
      ve/cre vci(9) r
      shell grep -E "^par3_" tmp_cal.txt >& tmp.txt
      shell $unquote('cat tmp.txt | sed "s/par3_/ /g" > tmp_cal.txt')
      ve/cre v3(3) r 
      ve/write v3,v3 tmp0.txt '(2f10.2)'
      shell cat tmp0.txt tmp_cal.txt > tmp1.txt
      ve/del np,vp
      ve/read np,vp tmp1.txt
      if ($vlen(vp).eq.30) then
        di=0
        do i=1,9
          i3=$sigma(([i]-[di])*3+1)
          ve/inp vai([i]) $sigma(vp([i3]))
          i3=$sigma(([i]-[di])*3+2)
          ve/inp vbi([i]) $sigma(vp([i3]))
          i3=$sigma(([i]-[di])*3+3)
          ve/inp vci([i]) $sigma(vp([i3]))
        enddo
      endif
    endif
  enddo
endif
return


macro calrungetx
p=hv
dir=/work/users/konctbel/AGPARS
ve/del runs,times
do part=1,30
  fname=[dir]/run_part[part].txt_[p].txt
  if ($fexist([fname])) then
    ve/del run,beam,ts,tf,[p]1,[p]2,[p]3,[p]4,[p]5,[p]6,[p]7,[p]8,[p]9,k[p]1,k[p]2,k[p]3,k[p]4,k[p]5,k[p]6,k[p]7,k[p]8,k[p]9
    ve/read run,beam,ts,tf,[p]1,[p]2,[p]3,[p]4,[p]5,[p]6,[p]7,[p]8,[p]9,k[p]1,k[p]2,k[p]3,k[p]4,k[p]5,k[p]6,k[p]7,k[p]8,k[p]9 [fname]
    exec vappend runs run
    exec vappend times ts
  endif
enddo
ve/del run,beam,ts,tf,[p]1,[p]2,[p]3,[p]4,[p]5,[p]6,[p]7,[p]8,[p]9,k[p]1,k[p]2,k[p]3,k[p]4,k[p]5,k[p]6,k[p]7,k[p]8,k[p]9
ve/del ti
ve/read ti calstime.txt
n=$vlen(ti)
ve/cre ri([n]) r
do i=1,[n]
  tii=ti([i])
  exec mapcal#ixndl [tii] times
  gl/imp inx
  ve/inp ri([i]) $sigma(runs([inx]))
  mess [i] [tii] $sigma(ri([i]))
enddo
ve/write ri calsrun.txt
return


macro caltimegetx
nmax=180
ve/cre ti([nmax]) r
do i=1,[nmax]
  exec mapcal#caltimeget [i]
  gl/imp ndays
  ve/inp ti([i]) [ndays]
enddo
ve/write ti calstime.txt
return


macro caltimeget ncal=180
shell rm tmp[0-9].txt
shell ls -d /work/users/konctbel/Calibr/Cal_* > tmp0.txt
shell $unquote('cat tmp0.txt | sed "s/[^0-9]/ /g" > tmp1.txt')
ve/del y,m,d,c
ve/read y,m,d,c tmp1.txt
exec mapcal#ixndl [ncal] c
gl/imp inx
in1=[inx]
exec mapcal#ixndr [ncal] c
gl/imp inx
in2=[inx]
if ([in1].eq.[in2]) then
  year=y([in1])
  month=m([in1])
  day=d([in1])
  cal=$format([ncal],i4.4)
  dir=/work/users/konctbel/Calibr/Cal_$format([year],i4.4)_$format([month],i2.2)_$format([day],i2.2)_[cal]
  fname=[dir]/time.txt
  ve/del t
  ve/read t [fname]
  year0=2011
  year=t(1)
  month=t(2)
  day=t(3)
  hour=t(4)
  minute=t(5)
  second=t(6)
  txt=$format([year],i4).$format([month],i2.2).$format([day],i2.2).$format([hour],i2.2).$format([minute],i2.2).$format([second],i2.2)
  shell /work/users/konctbel/MinuitTest/gtime [year0].01.01.00.00.00 [txt] gtime.txt
  ve/read gtime gtime.txt
  gl/cre ndays $sigma(gtime(1))
endif
return


macro callistpl
shell rm tmp[0-9].txt
shell ls -d /work/users/konctbel/Calibr/Cal_* > tmp0.txt
shell $unquote('cat tmp0.txt | sed "s/[^0-9]/ /g" > tmp1.txt')
ve/del y,m,d,c
ve/read y,m,d,c tmp1.txt
ve/del ri,ti
ve/read ri calsrun.txt
ve/read ti calstime.txt
fname=calslist.txt
if ($fexist([fname])) then
  shell rm [fname]
endif
for/file 20 [fname] ! n
close 20
do ncal=1,180
  exec mapcal#ixndl [ncal] c
  gl/imp inx
  in1=[inx]
  exec mapcal#ixndr [ncal] c
  gl/imp inx
  in2=[inx]
  if ([in1].eq.[in2]) then
    year=y([in1])
    month=m([in1])
    day=d([in1])
    cal=$format([ncal],i4.4)
    dir=Cal_$format([year],i4.4)_$format([month],i2.2)_$format([day],i2.2)_[cal]
    rii=ri([ncal])
    trii=$format([rii],i06)
    tcal=$format([cal],i04)
    tii=$sigma(ti([ncal]))
    ttii=$format([tii],f10.3)
    txt=$unquote([tcal])     [dir]     $unquote([trii])     $unquote([ttii])
    fmess [txt] [fname]
  endif
enddo
return



macro calparstab ncal1=36 ncal2=129 nc=1
*
dir0=/work/users/konctbel/Calibr
shell ls -d [dir0]/Cal_* > tmp0.txt
shell $unquote('cat tmp0.txt | sed "s/[^0-9]/ /g" > tmp1.txt')
ve/del y,m,d,c
ve/read y,m,d,c tmp1.txt
*
n=[ncal2]-[ncal1]+1
ve/cre chi2([n]) r
ve/cre ndf([n]) r
ve/cre af([n]) r
ve/cre daf([n]) r
ve/cre bf([n]) r
ve/cre dbf([n]) r
ve/cre cf([n]) r
ve/cre dcf([n]) r
ve/cre ndays([n]) r
ve/cre dndays([n]) r
*
ind=0
do ncal=[ncal1],[ncal2]
  exec mapcal#ixndl [ncal] c
  gl/imp inx
  in1=[inx]
  exec mapcal#ixndr [ncal] c
  gl/imp inx
  in2=[inx]
  if ([in1].eq.[in2]) then
    year=y([in1])
    month=m([in1])
    day=d([in1])
    cal=$format([ncal],i4.4)
    dir=[dir0]/Cal_$format([year],i4.4)_$format([month],i2.2)_$format([day],i2.2)_[cal]
    fname=[dir]/pars_correlation_counter[nc].txt
    shell cp [fname] tmp.txt
    shell $unquote('cat tmp.txt | sed "s/nan/10000/g" > tmp1.txt')
    shell $unquote('cat tmp1.txt | sed "s/inf/10000/g" > tmp.txt')
    ve/del pars
    ve/read pars tmp.txt
    if ($sigma(pars(1)).lt.100) then
      ind=[ind]+1
      ve/inp chi2([ind]) $sigma(pars(1))
      ve/inp ndf([ind]) $sigma(pars(2))
      ve/inp af([ind]) $sigma(pars(3))
      ve/inp daf([ind]) $sigma(pars(6))
      ve/inp bf([ind]) $sigma(pars(4))
      ve/inp dbf([ind]) $sigma(pars(7))
      ve/inp cf([ind]) $sigma(pars(5))
      ve/inp dcf([ind]) $sigma(pars(8))
      ve/del ymd
      ve/read ymd [dir]/time.txt
      year0=2011
      year=ymd(1)   ; txt=$format([year],i4.4)
      month=ymd(2)  ; txt=[txt].$format([month],i2.2)
      day=ymd(3)    ; txt=[txt].$format([day],i2.2)
      hour=ymd(4)   ; txt=[txt].$format([hour],i2.2)
      minute=ymd(5) ; txt=[txt].$format([minute],i2.2)
      second=ymd(6) ; txt=[txt].$format([second],i2.2)
      shell /work/users/konctbel/MinuitTest/gtime [year0].01.01.00.00.00 [txt] gtime.txt
      ve/read gtime gtime.txt
      ve/inp ndays([ind]) $sigma(gtime(1))
    endif
  endif
enddo
exec vpl#pl0 af daf ndays dndays sz=0.1
return



macro accled_cals_table
shell rm tmp[0-9].txt
shell ls -d /work/users/konctbel/Calibr/Cal_* > tmp0.txt
shell $unquote('cat tmp0.txt | sed "s/[^0-9]/ /g" > tmp1.txt')
ve/del y,m,d,c
ve/read y,m,d,c tmp1.txt
*
fname=accled_cals_table0.tex
if ($fexist([fname])) then
  shell rm [fname]
endif
for/file 20 [fname] ! N
close 20
*
do i=1,$vlen(c)
  year=y([i])
  month=m([i])
  day=d([i])
  ncal=c([i])
  dir=Cal\_$format([year],i4.4)\_$format([month],i2.2)\_$format([day],i2.2)\_$format([ncal],i4.4)
  txt=[dir] & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & \\
  fmess [txt] [fname]
  if ($sigma(mod([i],5)).eq.0) then
    txt=\hline
    fmess [txt] [fname]
  endif
enddo
*
fname=accled_cals_table_problems.txt
if ($fexist([fname])) then
  shell rm [fname]
endif
for/file 20 [fname] ! N
close 20
*
do i=1,$vlen(c)
  ncal=c([i])
  txt=$format([ncal],i4.4) 0 0 0 0 0 0 0 0 0
  fmess [txt] [fname]
enddo
return


macro accled_problems
*  
shell rm tmp[0-9].txt
shell ls -d /work/users/konctbel/Calibr/Cal_* > tmp0.txt
shell $unquote('cat tmp0.txt | sed "s/[^0-9]/ /g" > tmp1.txt')
ve/del y,m,d,c
ve/read y,m,d,c tmp1.txt
*
fname=accled_cals_table_problems.txt
ve/del ncals,p1,p2,p3,p4,p5,p6,p7,p8,p9
ve/read ncals,p1,p2,p3,p4,p5,p6,p7,p8,p9 accled_cals_table_problems.txt
*
fname=/work/users/konctbel/Calibr/problem_cals0.tex
if ($fexist([fname])) then
  shell rm [fname]
endif
for/file 20 [fname] ! N
close 20
*
n=$vlen(ncals)
np=0
do nc=1,9
  i123=0
  i12=0
  i1=0
  i2=0
  do i=1,[n]
    p=p[nc]([i])
    if (([p].ne.0).and.([i].eq.1)) then
      i1=-1
    endif
    if (([p].eq.0).and.([i123].ne.0)) then
      i2=[i]
    endif
    if (([p].eq.0).and.([i2].eq.0)) then
      i1=[i]
    endif
    if ([p].ne.0) then
      i123=[i]
    endif
    if (([p].ge.1).and.([p].le.2)) then
      i12=[i]
    endif
    if (([i1].ne.0).and.([i2].ne.0)) then
      if ([i12].ne.0) then
        np=[np]+1
        txt=Problem [np]: channel [nc]: { 
        do j=[i1],$sigma([i2]-1)
          txt=$unquote([txt])[j],
        enddo
        txt=$unquote([txt])[i2]}
        mess [txt]
*        
        ni=[i2]-([i1])+1
        ve/cre ji([ni]) r [i1] [i2]
        do j=$sigma([i1]+1),$sigma([i2]-1)
          jj=[j]-([i1])+2
          ve/inp ji([jj]) [j]
        enddo
*        
        do j=1,[ni]
          jj=ji([j])
          if (([jj].ge.1).and.([jj].le.[n])) then
            year=y([jj])
            month=m([jj])
            day=d([jj])
            ncal=c([jj])
            cal=$format([ncal],i4.4)
            dir=Cal_$format([year],i4.4)_$format([month],i2.2)_$format([day],i2.2)_$format([ncal],i4.4)
            fmess '\begin{figure}[ht!b]' [fname]
            fmess '  \begin{minipage}{0.49\textwidth}' [fname]
            fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
            epsfile=[dir]/picts/mu_vs_hv_counter[nc].eps
            fepsfile=/work/users/konctbel/Calibr/[epsfile]
            if ($fexist([fepsfile]).eq.0) then
              epsfile=null.eps
            endif
            txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
            fmess [txt] [fname]
            txt=$unquote('     ')\caption{[year]-[month]-[day] [cal]:[nc]}
            fmess [txt] [fname]
            fmess '   \end{minipage}' [fname]
            fmess '  \begin{minipage}{0.49\textwidth}' [fname]
            fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
            epsfile=[dir]/picts/amp_vs_hv_counter[nc].eps
            fepsfile=/work/users/konctbel/Calibr/[epsfile]
            if ($fexist([fepsfile]).eq.0) then
              epsfile=null.eps
            endif
            txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
            fmess [txt] [fname]
            txt=$unquote('     ')\caption{[year]-[month]-[day] [cal]:[nc]}
            fmess [txt] [fname]
            fmess '   \end{minipage}' [fname]
            fmess '\end{figure}' [fname]
          endif
        enddo
        fmess '\clearpage' [fname]
      endif
      i1=[i2]
      i2=0
      i12=0
      i123=0
    endif
  enddo
enddo
return


macro hv_cals
*
  ve/del r1,r2
  ve/read r1,r2 accled_corr.ixlist
  q=a
  ve/del v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9
  ve/read v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9 accled_[q]_v1.txt 9f10.3
  q=b
  ve/del v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9
  ve/read v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9 accled_[q]_v1.txt 9f10.3
  q=c
  ve/del v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9
  ve/read v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9 accled_[q]_v1.txt 9f10.3
  q=ae
  ve/del v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9
  ve/read v[q]1,v[q]2,v[q]3,v[q]4,v[q]5,v[q]6,v[q]7,v[q]8,v[q]9 accled_amp1pe_v1.txt 9f10.3
*
  n=$vlen(r1)
  do nc=1,9
    ve/cre u[nc]([n]) r
  enddo
  do i=1,[n]
    do nc=1,9
      a=va[nc]([i])
      b=vb[nc]([i])
      c=vc[nc]([i])
      a1=vae[nc]([i])
      if (([a].ne.0).and.([b].ne.0).and.([c].ne.0).and.([a1].ne.0)) then
        ui=$sigma(-3000*(([a1]-([b]))/([a]))**(1/([c])))
        ui=$sigma(int(abs([ui])+0.5))
        ve/inp u[nc]([i]) [ui]
      endif
    enddo
  enddo
*
  ve/write r1,r2,u1,u2,u3,u4,u5,u6,u7,u8,u9 hv_cals.txt '(11f10.1)'
return


macro cals_change
n=0
n=[n]+1; ic[n]=1; bc[n]=131; gc[n]=133
n=[n]+1; ic[n]=1; bc[n]=132; gc[n]=133
n=[n]+1; ic[n]=2; bc[n]=123; gc[n]=124
n=[n]+1; ic[n]=2; bc[n]=126; gc[n]=125
n=[n]+1; ic[n]=3; bc[n]=041; gc[n]=042
n=[n]+1; ic[n]=3; bc[n]=082; gc[n]=081
n=[n]+1; ic[n]=3; bc[n]=100; gc[n]=097
n=[n]+1; ic[n]=3; bc[n]=101; gc[n]=102
n=[n]+1; ic[n]=3; bc[n]=104; gc[n]=108
n=[n]+1; ic[n]=3; bc[n]=105; gc[n]=108
n=[n]+1; ic[n]=3; bc[n]=106; gc[n]=108
n=[n]+1; ic[n]=3; bc[n]=107; gc[n]=108
n=[n]+1; ic[n]=4; bc[n]=022; gc[n]=023
n=[n]+1; ic[n]=7; bc[n]=034; gc[n]=033
n=[n]+1; ic[n]=7; bc[n]=035; gc[n]=033
n=[n]+1; ic[n]=7; bc[n]=036; gc[n]=039
n=[n]+1; ic[n]=7; bc[n]=037; gc[n]=039
n=[n]+1; ic[n]=7; bc[n]=038; gc[n]=039
n=[n]+1; ic[n]=7; bc[n]=098; gc[n]=097
n=[n]+1; ic[n]=7; bc[n]=099; gc[n]=097
n=[n]+1; ic[n]=7; bc[n]=104; gc[n]=105
n=[n]+1; ic[n]=7; bc[n]=106; gc[n]=107
n=[n]+1; ic[n]=7; bc[n]=123; gc[n]=129
n=[n]+1; ic[n]=7; bc[n]=124; gc[n]=129
n=[n]+1; ic[n]=7; bc[n]=125; gc[n]=129
n=[n]+1; ic[n]=7; bc[n]=126; gc[n]=129
n=[n]+1; ic[n]=7; bc[n]=127; gc[n]=129
n=[n]+1; ic[n]=7; bc[n]=128; gc[n]=129
n=[n]+1; ic[n]=8; bc[n]=098; gc[n]=097
n=[n]+1; ic[n]=9; bc[n]=007; gc[n]=009
n=[n]+1; ic[n]=9; bc[n]=008; gc[n]=009
n=[n]+1; ic[n]=9; bc[n]=016; gc[n]=015
n=[n]+1; ic[n]=9; bc[n]=017; gc[n]=018
n=[n]+1; ic[n]=9; bc[n]=019; gc[n]=022
n=[n]+1; ic[n]=9; bc[n]=025; gc[n]=024
n=[n]+1; ic[n]=9; bc[n]=026; gc[n]=030
n=[n]+1; ic[n]=9; bc[n]=027; gc[n]=030
n=[n]+1; ic[n]=9; bc[n]=028; gc[n]=030
n=[n]+1; ic[n]=9; bc[n]=032; gc[n]=031
n=[n]+1; ic[n]=9; bc[n]=034; gc[n]=033
n=[n]+1; ic[n]=9; bc[n]=035; gc[n]=033
n=[n]+1; ic[n]=9; bc[n]=038; gc[n]=037
n=[n]+1; ic[n]=9; bc[n]=041; gc[n]=042
n=[n]+1; ic[n]=9; bc[n]=089; gc[n]=088
n=[n]+1; ic[n]=9; bc[n]=105; gc[n]=104
n=[n]+1; ic[n]=9; bc[n]=123; gc[n]=129
n=[n]+1; ic[n]=9; bc[n]=124; gc[n]=129
n=[n]+1; ic[n]=9; bc[n]=125; gc[n]=129
n=[n]+1; ic[n]=9; bc[n]=126; gc[n]=129
n=[n]+1; ic[n]=9; bc[n]=127; gc[n]=129
n=[n]+1; ic[n]=9; bc[n]=128; gc[n]=129
n=[n]+1; ic[n]=9; bc[n]=131; gc[n]=132
*
n=[n]+1; ic[n]=1; bc[n]=166; gc[n]=155
n=[n]+1; ic[n]=2; bc[n]=166; gc[n]=155
n=[n]+1; ic[n]=3; bc[n]=166; gc[n]=155
n=[n]+1; ic[n]=4; bc[n]=166; gc[n]=155
n=[n]+1; ic[n]=5; bc[n]=166; gc[n]=155
n=[n]+1; ic[n]=6; bc[n]=166; gc[n]=155
n=[n]+1; ic[n]=7; bc[n]=166; gc[n]=155
n=[n]+1; ic[n]=8; bc[n]=166; gc[n]=155
n=[n]+1; ic[n]=9; bc[n]=166; gc[n]=155
*
shell rm tmp[0-9].txt
shell ls -d /work/users/konctbel/Calibr/Cal_* > tmp0.txt
shell $unquote('cat tmp0.txt | sed "s/[^0-9]/ /g" > tmp1.txt')
ve/del y,m,d,c
ve/read y,m,d,c tmp1.txt
do i=1,[n]
  mess [i] counter=[ic[i]] change:  [bc[i]]  -->  [gc[i]]
*
  ncal=[bc[i]]
  exec mapcal#ixndl [ncal] c
  gl/imp inx
  in1=[inx]
  exec mapcal#ixndr [ncal] c
  gl/imp inx
  in2=[inx]
  if ([in1].eq.[in2]) then
    year=y([in1])
    month=m([in1])
    day=d([in1])
    dir=/work/users/konctbel/Calibr/Cal_$format([year],i4.4)_$format([month],i2.2)_$format([day],i2.2)_$format([ncal],i4.4)
    fname=[dir]/change_counter[ic[i]].txt
    if ($fexist([fname])) then
      shell rm [fname]
    endif
    for/file 20 [fname] ! n
    close 20
    fmess $format([gc[i]],i4.4) [fname]
  endif
enddo
return


macro cals_to_db
n=0
*n=[n]+1; clb[n]=0001; rb[n]= 6407;
*n=[n]+1; clb[n]=0002; rb[n]= 6384;
n=[n]+1; clb[n]=0003; rb[n]= 6384; x=1
*n=[n]+1; clb[n]=0004; rb[n]= 6384;
n=[n]+1; clb[n]=0005; rb[n]= 6858;
n=[n]+1; clb[n]=0006; rb[n]= 6886;
*n=[n]+1; clb[n]=0007; rb[n]= 6975;
*n=[n]+1; clb[n]=0008; rb[n]= 6975;
*n=[n]+1; clb[n]=0009; rb[n]= 6975;
*n=[n]+1; clb[n]=0010; rb[n]= 6975;
*n=[n]+1; clb[n]=0011; rb[n]= 6975;
*n=[n]+1; clb[n]=0012; rb[n]= 6975;
*n=[n]+1; clb[n]=0013; rb[n]= 6975;
*n=[n]+1; clb[n]=0014; rb[n]= 6975;
n=[n]+1; clb[n]=0015; rb[n]= 6975; x=1
n=[n]+1; clb[n]=0016; rb[n]= 6991; x=1
*n=[n]+1; clb[n]=0017; rb[n]= 6992;
n=[n]+1; clb[n]=0018; rb[n]= 7018;
n=[n]+1; clb[n]=0019; rb[n]= 7089;
*n=[n]+1; clb[n]=0020; rb[n]= 7116;
*n=[n]+1; clb[n]=0021; rb[n]= 7116;
*n=[n]+1; clb[n]=0022; rb[n]= 7116;
*n=[n]+1; clb[n]=0023; rb[n]= 7116;
*n=[n]+1; clb[n]=0024; rb[n]= 7116;
n=[n]+1; clb[n]=0025; rb[n]= 7102; x=1
n=[n]+1; clb[n]=0026; rb[n]= 7198; x=1
n=[n]+1; clb[n]=0027; rb[n]= 7206;
n=[n]+1; clb[n]=0028; rb[n]= 7264;
*n=[n]+1; clb[n]=0029; rb[n]= 7340;
n=[n]+1; clb[n]=0030; rb[n]= 7340;
*n=[n]+1; clb[n]=0031; rb[n]= 7340;
*n=[n]+1; clb[n]=0032; rb[n]= 7340;
*n=[n]+1; clb[n]=0033; rb[n]= 7340;
n=[n]+1; clb[n]=0034; rb[n]= 7372;
*n=[n]+1; clb[n]=0035; rb[n]= 7405;
*n=[n]+1; clb[n]=0036; rb[n]= 7828;
n=[n]+1; clb[n]=0037; rb[n]= 7828;
*n=[n]+1; clb[n]=0038; rb[n]= 7840;
n=[n]+1; clb[n]=0039; rb[n]= 7840; x=1
n=[n]+1; clb[n]=0039; rb[n]= 7933; x=1
*n=[n]+1; clb[n]=0040; rb[n]= 7840;
n=[n]+1; clb[n]=0041; rb[n]= 8604;
n=[n]+1; clb[n]=0042; rb[n]= 8816;
n=[n]+1; clb[n]=0043; rb[n]= 8849; x=1
n=[n]+1; clb[n]=0044; rb[n]= 8945;
n=[n]+1; clb[n]=0045; rb[n]= 8974;
n=[n]+1; clb[n]=0046; rb[n]= 9028;
n=[n]+1; clb[n]=0047; rb[n]= 9098;
n=[n]+1; clb[n]=0048; rb[n]= 9222;
n=[n]+1; clb[n]=0049; rb[n]= 9294;
n=[n]+1; clb[n]=0050; rb[n]= 9367;
n=[n]+1; clb[n]=0051; rb[n]= 9508;
n=[n]+1; clb[n]=0052; rb[n]= 9556;
n=[n]+1; clb[n]=0053; rb[n]= 9579;
n=[n]+1; clb[n]=0054; rb[n]= 9640;
n=[n]+1; clb[n]=0055; rb[n]= 9703;
n=[n]+1; clb[n]=0056; rb[n]= 9767;
n=[n]+1; clb[n]=0057; rb[n]= 9799;
n=[n]+1; clb[n]=0058; rb[n]= 9932;
n=[n]+1; clb[n]=0059; rb[n]=10051;
n=[n]+1; clb[n]=0060; rb[n]=10103;
n=[n]+1; clb[n]=0061; rb[n]=10138;
*n=[n]+1; clb[n]=0062; rb[n]=10138;
n=[n]+1; clb[n]=0063; rb[n]=10167; x=1
n=[n]+1; clb[n]=0064; rb[n]=10205;
n=[n]+1; clb[n]=0065; rb[n]=10316; x=1
n=[n]+1; clb[n]=0066; rb[n]=10374;
n=[n]+1; clb[n]=0067; rb[n]=10429;
n=[n]+1; clb[n]=0068; rb[n]=10467;
*n=[n]+1; clb[n]=0069; rb[n]=10539;
n=[n]+1; clb[n]=0070; rb[n]=10539;
*n=[n]+1; clb[n]=0071; rb[n]=10675;
*n=[n]+1; clb[n]=0072; rb[n]=10675;
n=[n]+1; clb[n]=0073; rb[n]=10675;
n=[n]+1; clb[n]=0074; rb[n]=10703;
n=[n]+1; clb[n]=0075; rb[n]=10746;
n=[n]+1; clb[n]=0076; rb[n]=10864;
n=[n]+1; clb[n]=0077; rb[n]=10905;
*n=[n]+1; clb[n]=0078; rb[n]=10998;
*n=[n]+1; clb[n]=0079; rb[n]=10998;
*n=[n]+1; clb[n]=0080; rb[n]=10998;
*n=[n]+1; clb[n]=0081; rb[n]=10998;
*n=[n]+1; clb[n]=0082; rb[n]=10998;
*n=[n]+1; clb[n]=0083; rb[n]=10998;
*n=[n]+1; clb[n]=0084; rb[n]=10998;
*n=[n]+1; clb[n]=0085; rb[n]=10998;
*n=[n]+1; clb[n]=0086; rb[n]=10998;
*n=[n]+1; clb[n]=0087; rb[n]=10998;
n=[n]+1; clb[n]=0088; rb[n]=10998; x=1
*n=[n]+1; clb[n]=0089; rb[n]=10998;
*n=[n]+1; clb[n]=0090; rb[n]=11022;
n=[n]+1; clb[n]=0091; rb[n]=11022;
*n=[n]+1; clb[n]=0092; rb[n]=11046;
n=[n]+1; clb[n]=0093; rb[n]=11046;
n=[n]+1; clb[n]=0094; rb[n]=11065;
n=[n]+1; clb[n]=0095; rb[n]=11092;
n=[n]+1; clb[n]=0096; rb[n]=11159;
n=[n]+1; clb[n]=0097; rb[n]=11243;
*n=[n]+1; clb[n]=0098; rb[n]=11243;
*n=[n]+1; clb[n]=0099; rb[n]=11278;
n=[n]+1; clb[n]=0100; rb[n]=11340;
*n=[n]+1; clb[n]=0101; rb[n]=11345;
*n=[n]+1; clb[n]=0102; rb[n]=11345;
n=[n]+1; clb[n]=0103; rb[n]=11345;
*n=[n]+1; clb[n]=0104; rb[n]=11453;
n=[n]+1; clb[n]=0105; rb[n]=11453;
n=[n]+1; clb[n]=0106; rb[n]=11522;
*n=[n]+1; clb[n]=0107; rb[n]=11522;
*n=[n]+1; clb[n]=0108; rb[n]=11570;
*n=[n]+1; clb[n]=0109; rb[n]=11570;
n=[n]+1; clb[n]=0110; rb[n]=11570;
n=[n]+1; clb[n]=0111; rb[n]=11750;
n=[n]+1; clb[n]=0112; rb[n]=11915;
n=[n]+1; clb[n]=0113; rb[n]=12108;
n=[n]+1; clb[n]=0114; rb[n]=12256;
n=[n]+1; clb[n]=0115; rb[n]=12396;
n=[n]+1; clb[n]=0116; rb[n]=12492;
n=[n]+1; clb[n]=0117; rb[n]=12574; x=1
n=[n]+1; clb[n]=0118; rb[n]=12706;
n=[n]+1; clb[n]=0119; rb[n]=12971;
*n=[n]+1; clb[n]=0120; rb[n]=13097;
*n=[n]+1; clb[n]=0121; rb[n]=13115;
n=[n]+1; clb[n]=0122; rb[n]=13154;
n=[n]+1; clb[n]=0123; rb[n]=13265;
*n=[n]+1; clb[n]=0124; rb[n]=13305;
n=[n]+1; clb[n]=0125; rb[n]=13305;
n=[n]+1; clb[n]=0126; rb[n]=13386;
*
n=[n]+1; clb[n]=0127; rb[n]=13504; x=1
n=[n]+1; clb[n]=0127; rb[n]=13589; x=1
n=[n]+1; clb[n]=0127; rb[n]=13601; x=1
n=[n]+1; clb[n]=0127; rb[n]=13606; x=1
n=[n]+1; clb[n]=0127; rb[n]=13624; x=1
*
n=[n]+1; clb[n]=0128; rb[n]=13648; x=1
n=[n]+1; clb[n]=0129; rb[n]=13845;
*n=[n]+1; clb[n]=0130; rb[n]=13892;
*n=[n]+1; clb[n]=0131; rb[n]=13935;
*n=[n]+1; clb[n]=0132; rb[n]=13935;
*n=[n]+1; clb[n]=0133; rb[n]=13935;
n=[n]+1; clb[n]=0134; rb[n]=13935; x=1
n=[n]+1; clb[n]=0135; rb[n]=14031;
n=[n]+1; clb[n]=0136; rb[n]=14171;
*n=[n]+1; clb[n]=0137; rb[n]=14451;
n=[n]+1; clb[n]=0138; rb[n]=14451;
n=[n]+1; clb[n]=0139; rb[n]=14647;
n=[n]+1; clb[n]=0140; rb[n]=14936;
n=[n]+1; clb[n]=0141; rb[n]=14992;
n=[n]+1; clb[n]=0142; rb[n]=15535;
n=[n]+1; clb[n]=0143; rb[n]=15437;
n=[n]+1; clb[n]=0144; rb[n]=15535;
n=[n]+1; clb[n]=0145; rb[n]=15731;
n=[n]+1; clb[n]=0146; rb[n]=15940;
n=[n]+1; clb[n]=0147; rb[n]=15951;
*n=[n]+1; clb[n]=0148; rb[n]=15955;
n=[n]+1; clb[n]=0149; rb[n]=15955; x=1
n=[n]+1; clb[n]=0150; rb[n]=15990;
n=[n]+1; clb[n]=0151; rb[n]=16319;
n=[n]+1; clb[n]=0152; rb[n]=16436;
n=[n]+1; clb[n]=0153; rb[n]=16549;
n=[n]+1; clb[n]=0154; rb[n]=16599;
*n=[n]+1; clb[n]=0155; rb[n]=17869;
n=[n]+1; clb[n]=0156; rb[n]=16697; x=1
n=[n]+1; clb[n]=0157; rb[n]=16701; x=1
n=[n]+1; clb[n]=0158; rb[n]=16835; x=1
*n=[n]+1; clb[n]=0159; rb[n]=16882;
n=[n]+1; clb[n]=0160; rb[n]=16882;
n=[n]+1; clb[n]=0161; rb[n]=16980;
n=[n]+1; clb[n]=0162; rb[n]=17094;
n=[n]+1; clb[n]=0163; rb[n]=17275;
n=[n]+1; clb[n]=0164; rb[n]=17355;
n=[n]+1; clb[n]=0165; rb[n]=17670;
n=[n]+1; clb[n]=0166; rb[n]=17869; x=1
n=[n]+1; clb[n]=0167; rb[n]=17975; x=1
n=[n]+1; clb[n]=0168; rb[n]=18325;
*n=[n]+1; clb[n]=0169; rb[n]=18384;
n=[n]+1; clb[n]=0170; rb[n]=18384;
n=[n]+1; clb[n]=0171; rb[n]=18610;
n=[n]+1; clb[n]=0172; rb[n]=18830;
n=[n]+1; clb[n]=0173; rb[n]=19134;
n=[n]+1; clb[n]=0174; rb[n]=19243;
n=[n]+1; clb[n]=0175; rb[n]=19565;
n=[n]+1; clb[n]=0176; rb[n]=19802;
n=[n]+1; clb[n]=0177; rb[n]=20279;
n=[n]+1; clb[n]=0178; rb[n]=20776;
n=[n]+1; clb[n]=0179; rb[n]=21041;
n=[n]+1; clb[n]=0180; rb[n]=21558;
*
shell rm tmp[0-9].txt
shell ls -d /work/users/konctbel/Calibr/Cal_* > tmp0.txt
shell $unquote('cat tmp0.txt | sed "s/[^0-9]/ /g" > tmp1.txt')
ve/del y,m,d,c
ve/read y,m,d,c tmp1.txt
*
ve/del r1db,r2db,u1db,u2db,u3db,u4db,u5db,u6db,u7db,u8db,u9db
ve/read r1db,r2db,u1db,u2db,u3db,u4db,u5db,u6db,u7db,u8db,u9db hv_cals_corr_set.txt
do nc=1,9
  ve/cre uw[nc]([n]) r
  ve/cre a[nc]([n]) r
  ve/cre b[nc]([n]) r
  ve/cre c[nc]([n]) r
  ve/cre ae[nc]([n]) r
  ve/cre ee[nc]([n]) r
  ve/cre u0[nc]([n]) r
  ve/cre s0[nc]([n]) r
  ve/cre pds[nc]([n]) r
  ve/cre tmin[nc]([n]) r
  ve/cre tmax[nc]([n]) r
enddo
ve/cre run([n]) r
ve/cre clb([n]) r
nr=$vlen(r2db)
ve/inp r2db([nr]) 1000000
ve/cre uwx(1) r
ind=1
do i=1,[n]
  ri=[rb[i]]
  ve/inp run([i]) [ri]
  while (([ri].gt.$sigma(r2db([ind]))).and.([ind].lt.$vlen(r2db))) do
    ind=[ind]+1
  endwhile
  do nc=1,9
    ve/inp uw[nc]([i]) $sigma(u[nc]db([ind]))
  enddo
* 
  ncal=[clb[i]]
  ve/inp clb([i]) [ncal]
  exec mapcal#ixndl [ncal] c
  gl/imp inx
  in1=[inx]
  exec mapcal#ixndr [ncal] c
  gl/imp inx
  in2=[inx]
  if ([in1].eq.[in2]) then
    year=y([in1])
    month=m([in1])
    day=d([in1])
    dir0=/work/users/konctbel/Calibr/Cal_$format([year],i4.4)_$format([month],i2.2)_$format([day],i2.2)_$format([ncal],i4.4)
    do nc=1,9
      dir=[dir0]
      fname=[dir]/change_counter[nc].txt
      if ($fexist([fname])) then
        ve/read vt [fname]
        ncalch=vt(1)
        exec mapcal#ixndl [ncalch] c
        gl/imp inx
        in1=[inx]
        exec mapcal#ixndr [ncalch] c
        gl/imp inx
        in2=[inx]
        if ([in1].eq.[in2]) then
          year=y([in1])
          month=m([in1])
          day=d([in1])
          dir=/work/users/konctbel/Calibr/Cal_$format([year],i4.4)_$format([month],i2.2)_$format([day],i2.2)_$format([ncalch],i4.4)
        else
          dir=figvam
        endif
      endif
      fname=[dir]/pars_correlation_counter[nc].txt
      fnamed=[dir]/fit_data_counter[nc].txt
      fnamem=[dir]/pars_mu_correlation_counter[nc].txt
      fnamep=[dir]/pds_data_counter[nc].txt
      if (($fexist([fname]).eq.1).and.($fexist([fnamem]).eq.1)) then
        ve/read apars [fname]
        ve/read mpars [fnamem]
        ve/read u0i,pdsi,dpdsi,pds0i,pds0fi,pdsfi,tmini,tmaxi [fnamep]
        npi=$vlen(u0i)
        a=apars(4)
        b=apars(5)
        c=apars(6)
        u=$sigma(uw[nc]([i]))
        ae=$sigma([a]*([u]/3000)**([c])+([b]))
        ve/inp a[nc]([i]) [a]
        ve/inp b[nc]([i]) [b]
        ve/inp c[nc]([i]) [c]
        ve/inp ae[nc]([i]) [ae]
        ve/copy mpars(4:7) p4
        ve/inp uwx(1) [u]
        mux=$call('erfp0p.f(uwx)')
        ve/inp ee[nc]([i]) $sigma([mux]/p4(2))
        ve/inp u0[nc]([i]) $sigma(p4(3))
        ve/inp s0[nc]([i]) $sigma(abs(p4(4)))
        ve/inp pds[nc]([i]) $sigma(vsum(pdsi)/[npi])
        ve/inp tmin[nc]([i]) $sigma(vsum(tmini)/[npi])
        ve/inp tmax[nc]([i]) $sigma(vsum(tmaxi)/[npi])
      endif
    enddo
  endif
enddo
*
n=0
n=[n]+1; v[n]=uw; f[n]=u 
n=[n]+1; v[n]=a; f[n]=a 
n=[n]+1; v[n]=b; f[n]=b 
n=[n]+1; v[n]=c; f[n]=c 
n=[n]+1; v[n]=ae; f[n]=amp1pe
n=[n]+1; v[n]=ee; f[n]=eff1pe
n=[n]+1; v[n]=u0; f[n]=u0
n=[n]+1; v[n]=s0; f[n]=s0
n=[n]+1; v[n]=pds; f[n]=pds
n=[n]+1; v[n]=tmin; f[n]=tmin
n=[n]+1; v[n]=tmax; f[n]=tmax
*
do v=1,[n]
  fname=accled_[f[v]]_work.txt
  if ($fexist([fname])) then
    shell rm [fname]
  endif
  for/file 20 [fname] ! n
  close 20
  ve/write clb,run,[v[v]]1,[v[v]]2,[v[v]]3,[v[v]]4,[v[v]]5,[v[v]]6,[v[v]]7,[v[v]]8,[v[v]]9 [fname] '(2f10.1,9f15.6)'
enddo
return


macro mapcal1
dir0=/work/snd2000/users/konctbel/exp/MHAD2011-4/col
dir1=mapcal1
n=0
n=[n]+1; fname[n]=mhad2011_550_col_p1_a.hbook          
n=[n]+1; fname[n]=mhad2011_750-1_col_p1_a.hbook         
n=[n]+1; fname[n]=mhad2011_525_col_p4.hbook             
n=[n]+1; fname[n]=mhad2011_550_col_p1_b.hbook           
n=[n]+1; fname[n]=mhad2011_525_col_p3.hbook             
n=[n]+1; fname[n]=mhad2011_525_col_p2.hbook             
n=[n]+1; fname[n]=mhad2011_525_col_p1.hbook             
n=[n]+1; fname[n]=mhad2011_600_col_p1_a.hbook           
n=[n]+1; fname[n]=mhad2011_600_col_p2_a.hbook           
n=[n]+1; fname[n]=mhad2011_650_col_p2.hbook             
n=[n]+1; fname[n]=mhad2011_750-1_col_p1_b.hbook         
n=[n]+1; fname[n]=mhad2011_575_col_p3.hbook             
n=[n]+1; fname[n]=mhad2011_600_col_p1_b.hbook           
n=[n]+1; fname[n]=mhad2011_650_col_p1.hbook             
n=[n]+1; fname[n]=mhad2011_600_col_p2_b.hbook           
n=[n]+1; fname[n]=mhad2011_575_col_p2.hbook             
n=[n]+1; fname[n]=mhad2011_675_col_p2_a.hbook           
n=[n]+1; fname[n]=mhad2011_550_col_p3.hbook             
n=[n]+1; fname[n]=mhad2011_650_col_p3.hbook             
n=[n]+1; fname[n]=mhad2011_625_col_p2.hbook             
n=[n]+1; fname[n]=mhad2011_675_col_p2_b.hbook           
n=[n]+1; fname[n]=mhad2011_550_col_p2.hbook             
n=[n]+1; fname[n]=mhad2011_750-1_col_p2.hbook           
n=[n]+1; fname[n]=mhad2011_625_col_p1.hbook             
n=[n]+1; fname[n]=mhad2011_825_col_p2.hbook             
n=[n]+1; fname[n]=mhad2011_575_col_p1.hbook             
n=[n]+1; fname[n]=mhad2011_825_col_p1.hbook             
n=[n]+1; fname[n]=mhad2011_750_col_p1_a.hbook           
n=[n]+1; fname[n]=mhad2011_675_col_p1.hbook             
n=[n]+1; fname[n]=mhad2011_750_col_p1_b.hbook           
n=[n]+1; fname[n]=mhad2011_750_col_p2_b.hbook           
n=[n]+1; fname[n]=mhad2011_750_col_p2_a.hbook           
n=[n]+1; fname[n]=mhad2011_725_col_p1.hbook             
n=[n]+1; fname[n]=mhad2011_800_col_p1.hbook             
n=[n]+1; fname[n]=mhad2011_725_col_p2.hbook             
n=[n]+1; fname[n]=mhad2011_700_col_p1.hbook             
n=[n]+1; fname[n]=mhad2011_800_col_p2.hbook             
n=[n]+1; fname[n]=mhad2011_850_col_p1.hbook             
n=[n]+1; fname[n]=mhad2011_950_col_p1.hbook             
n=[n]+1; fname[n]=mhad2011_925_col_p2.hbook             
n=[n]+1; fname[n]=mhad2011_950_col_p2.hbook             
n=[n]+1; fname[n]=mhad2011_700_col_p2.hbook             
n=[n]+1; fname[n]=mhad2011_975_col_p2.hbook             
n=[n]+1; fname[n]=mhad2011_775_col_p1.hbook             
n=[n]+1; fname[n]=mhad2011_925_col_p1.hbook             
n=[n]+1; fname[n]=mhad2011_775_col_p2.hbook             
n=[n]+1; fname[n]=mhad2011_875_col_p1.hbook             
n=[n]+1; fname[n]=mhad2011_962.5_col_p2.hbook           
n=[n]+1; fname[n]=mhad2011_850_col_p2.hbook             
n=[n]+1; fname[n]=mhad2011_975_col_p1.hbook             
n=[n]+1; fname[n]=mhad2011_950_col_p3.hbook             
n=[n]+1; fname[n]=mhad2011_825_col_p3.hbook             
n=[n]+1; fname[n]=mhad2011_962.5_col_p3.hbook           
n=[n]+1; fname[n]=mhad2011_975_col_p3.hbook             
n=[n]+1; fname[n]=mhad2011_987.5_col_p1.hbook           
n=[n]+1; fname[n]=mhad2011_962.5_col_p1.hbook           
n=[n]+1; fname[n]=mhad2011_875_col_p2.hbook             
n=[n]+1; fname[n]=mhad2011_1000_col_p3.hbook            
n=[n]+1; fname[n]=mhad2011_900_col.hbook                
n=[n]+1; fname[n]=mhad2011_987.5_col_p2.hbook           
n=[n]+1; fname[n]=mhad2011_987.5_col_p3.hbook           
n=[n]+1; fname[n]=mhad2011_1000_col_p1.hbook            
n=[n]+1; fname[n]=mhad2011_1000_col_p2.hbook            
n=[n]+1; fname[n]=mhad2011_935_col_p2.hbook             
n=[n]+1; fname[n]=mhad2011_862.5_col_a.hbook            
n=[n]+1; fname[n]=mhad2011_945_col_p1.hbook             
n=[n]+1; fname[n]=mhad2011_812.5_col_p2.hbook           
n=[n]+1; fname[n]=mhad2011_935_col_p1.hbook             
n=[n]+1; fname[n]=mhad2011_945_col_p2.hbook             
n=[n]+1; fname[n]=mhad2011_862.5_col_b.hbook            
n=[n]+1; fname[n]=mhad2011_812.5_col_p1.hbook           
n=[n]+1; fname[n]=mhad2011_762.5_col_p1.hbook           
n=[n]+1; fname[n]=mhad2011_787.5_col_p2.hbook           
n=[n]+1; fname[n]=mhad2011_712.5_col_p1.hbook           
n=[n]+1; fname[n]=mhad2011_787.5_col_p1.hbook           
n=[n]+1; fname[n]=mhad2011_762.5_col_p2.hbook           
n=[n]+1; fname[n]=mhad2011_737.5_col_p2.hbook           
n=[n]+1; fname[n]=mhad2011_637.5_col_p2_a.hbook         
n=[n]+1; fname[n]=mhad2011_737.5_col_p1.hbook           
n=[n]+1; fname[n]=mhad2011_887.5_col.hbook              
n=[n]+1; fname[n]=mhad2011_712.5_col_p2.hbook           
n=[n]+1; fname[n]=mhad2011_587.5_col_p2_a.hbook         
n=[n]+1; fname[n]=mhad2011_612.5_col_p1_a.hbook         
n=[n]+1; fname[n]=mhad2011_912.5_col.hbook              
n=[n]+1; fname[n]=mhad2011_837.5_col.hbook              
n=[n]+1; fname[n]=mhad2011_587.5_col_p1_a.hbook         
n=[n]+1; fname[n]=mhad2011_637.5_col_p2_b.hbook         
n=[n]+1; fname[n]=mhad2011_637.5_col_p1.hbook           
n=[n]+1; fname[n]=mhad2011_562.5_col_p1_a.hbook         
n=[n]+1; fname[n]=mhad2011_687.5_col_p1.hbook           
n=[n]+1; fname[n]=mhad2011_612.5_col_p1_b.hbook         
n=[n]+1; fname[n]=mhad2011_587.5_col_p1_b.hbook         
n=[n]+1; fname[n]=mhad2011_562.5_col_p1_b.hbook         
n=[n]+1; fname[n]=mhad2011_662.5_col_p1.hbook           
n=[n]+1; fname[n]=mhad2011_687.5_col_p2.hbook          
n=[n]+1; fname[n]=mhad2011_587.5_col_p2_b.hbook        
n=[n]+1; fname[n]=mhad2011_537.5_col_p1_a.hbook        
n=[n]+1; fname[n]=mhad2011_537.5_col_p1_b.hbook        
n=[n]+1; fname[n]=mhad2011_662.5_col_p2.hbook          
n=[n]+1; fname[n]=mhad2011_508.7_col.hbook             
n=[n]+1; fname[n]=mhad2011_537.5_col_p2.hbook          
n=[n]+1; fname[n]=mhad2011_510.5-2_col_p1.hbook        
n=[n]+1; fname[n]=mhad2011_612.5_col_p2.hbook          
n=[n]+1; fname[n]=mhad2011_510.5-2_col_p4.hbook        
n=[n]+1; fname[n]=mhad2011_510.5-2_col_p2.hbook        
n=[n]+1; fname[n]=mhad2011_537.5_col_p3.hbook          
n=[n]+1; fname[n]=mhad2011_562.5_col_p2.hbook          
n=[n]+1; fname[n]=mhad2011_510.5-2_col_p3.hbook        
n=[n]+1; fname[n]=mhad2011_509.8_col.hbook             
do i=1,[n]
  shell ln -s [dir0]/[fname[i]] [dir1]/[fname[i]]
enddo
mess [n]
return



macro mapcal2
dir0=/work/users/konctbel/exp/MHAD2012-2/col
dir1=mapcal2
n=0
n=[n]+1; fname[n]=mhad2012_505.3.hbook             
n=[n]+1; fname[n]=mhad2012_508.0_p2.hbook          
n=[n]+1; fname[n]=mhad2012_508.0_p1.hbook          
n=[n]+1; fname[n]=mhad2012_510.0_p5.hbook          
n=[n]+1; fname[n]=mhad2012_509.0_p3.hbook          
n=[n]+1; fname[n]=mhad2012_512.3.hbook             
n=[n]+1; fname[n]=mhad2012_509.0_p1.hbook          
n=[n]+1; fname[n]=mhad2012_509.0_p2.hbook          
n=[n]+1; fname[n]=mhad2012_511.0_p1.hbook          
n=[n]+1; fname[n]=mhad2012_510.0_p4.hbook          
n=[n]+1; fname[n]=mhad2012_513.0.hbook             
n=[n]+1; fname[n]=mhad2012_640.0_p5.hbook          
n=[n]+1; fname[n]=mhad2012_511.0_p2.hbook          
n=[n]+1; fname[n]=mhad2012_510.0_p1.hbook          
n=[n]+1; fname[n]=mhad2012_512.0.hbook             
n=[n]+1; fname[n]=mhad2012_760.0_p6.hbook          
n=[n]+1; fname[n]=mhad2012_510.0_p3.hbook          
n=[n]+1; fname[n]=mhad2012_640.0_p3.hbook          
n=[n]+1; fname[n]=mhad2012_510.0_p2.hbook          
n=[n]+1; fname[n]=mhad2012_720.0_p1.hbook          
n=[n]+1; fname[n]=mhad2012_640.0_p4.hbook          
n=[n]+1; fname[n]=mhad2012_760.0_p4.hbook          
n=[n]+1; fname[n]=mhad2012_680.0_p1.hbook          
n=[n]+1; fname[n]=mhad2012_760.0_p1.hbook          
n=[n]+1; fname[n]=mhad2012_720.0_p5.hbook          
n=[n]+1; fname[n]=mhad2012_640.0_p2.hbook          
n=[n]+1; fname[n]=mhad2012_640.0_p1.hbook          
n=[n]+1; fname[n]=mhad2012_760.0_p2.hbook          
n=[n]+1; fname[n]=mhad2012_680.0_p4.hbook          
n=[n]+1; fname[n]=mhad2012_860.0_p4.hbook          
n=[n]+1; fname[n]=mhad2012_720.0_p4.hbook          
n=[n]+1; fname[n]=mhad2012_760.0_p3.hbook          
n=[n]+1; fname[n]=mhad2012_680.0_p2.hbook          
n=[n]+1; fname[n]=mhad2012_680.0_p3.hbook          
n=[n]+1; fname[n]=mhad2012_680.0_p5.hbook          
n=[n]+1; fname[n]=mhad2012_760.0_p5.hbook          
n=[n]+1; fname[n]=mhad2012_720.0_p2.hbook          
n=[n]+1; fname[n]=mhad2012_840.0_p3.hbook          
n=[n]+1; fname[n]=mhad2012_840.0_p1.hbook          
n=[n]+1; fname[n]=mhad2012_720.0_p3.hbook          
n=[n]+1; fname[n]=mhad2012_800.0_p1.hbook          
n=[n]+1; fname[n]=mhad2012_860.0_p2.hbook          
n=[n]+1; fname[n]=mhad2012_880.0_p5.hbook          
n=[n]+1; fname[n]=mhad2012_800.0_p2.hbook          
n=[n]+1; fname[n]=mhad2012_800.0_p4.hbook          
n=[n]+1; fname[n]=mhad2012_840.0_p5.hbook          
n=[n]+1; fname[n]=mhad2012_840.0_p2.hbook          
n=[n]+1; fname[n]=mhad2012_860.0_p3.hbook          
n=[n]+1; fname[n]=mhad2012_860.0_p1.hbook          
n=[n]+1; fname[n]=mhad2012_900.0_p4.hbook          
n=[n]+1; fname[n]=mhad2012_920.0_p5.hbook          
n=[n]+1; fname[n]=mhad2012_800.0_p3.hbook          
n=[n]+1; fname[n]=mhad2012_880.0_p2.hbook          
n=[n]+1; fname[n]=mhad2012_880.0_p1.hbook          
n=[n]+1; fname[n]=mhad2012_880.0_p3.hbook          
n=[n]+1; fname[n]=mhad2012_936.0_p5.hbook          
n=[n]+1; fname[n]=mhad2012_840.0_p4.hbook          
n=[n]+1; fname[n]=mhad2012_880.0_p4.hbook          
n=[n]+1; fname[n]=mhad2012_920.0_p2.hbook          
n=[n]+1; fname[n]=mhad2012_920.0_p3.hbook          
n=[n]+1; fname[n]=mhad2012_936.0_p2.hbook          
n=[n]+1; fname[n]=mhad2012_920.0_p4.hbook          
n=[n]+1; fname[n]=mhad2012_920.0_p1.hbook          
n=[n]+1; fname[n]=mhad2012_960.0_p4.hbook          
n=[n]+1; fname[n]=mhad2012_950.0_p3.hbook          
n=[n]+1; fname[n]=mhad2012_950.0_p1.hbook          
n=[n]+1; fname[n]=mhad2012_950.0_p6.hbook          
n=[n]+1; fname[n]=mhad2012_936.0_p3.hbook          
n=[n]+1; fname[n]=mhad2012_900.0_p2.hbook          
n=[n]+1; fname[n]=mhad2012_900.0_p3.hbook          
n=[n]+1; fname[n]=mhad2012_900.0_p1.hbook          
n=[n]+1; fname[n]=mhad2012_936.0_p1.hbook          
n=[n]+1; fname[n]=mhad2012_950.0_p2.hbook          
n=[n]+1; fname[n]=mhad2012_960.0_p3.hbook          
n=[n]+1; fname[n]=mhad2012_950.0_p4.hbook          
n=[n]+1; fname[n]=mhad2012_960.0_p1.hbook          
n=[n]+1; fname[n]=mhad2012_950.0_p5.hbook          
n=[n]+1; fname[n]=mhad2012_980.0_p7.hbook          
n=[n]+1; fname[n]=mhad2012_970.0_p1.hbook          
n=[n]+1; fname[n]=mhad2012_936.0_p4.hbook          
n=[n]+1; fname[n]=mhad2012_960.0_p2.hbook          
n=[n]+1; fname[n]=mhad2012_980.0_p1.hbook          
n=[n]+1; fname[n]=mhad2012_970.0_p6.hbook          
n=[n]+1; fname[n]=mhad2012_970.0_p5.hbook          
n=[n]+1; fname[n]=mhad2012_980.0_p5.hbook          
n=[n]+1; fname[n]=mhad2012_980.0_p4.hbook          
n=[n]+1; fname[n]=mhad2012_980.0_p2.hbook          
n=[n]+1; fname[n]=mhad2012_970.0_p3.hbook          
n=[n]+1; fname[n]=mhad2012_970.0_p4.hbook          
n=[n]+1; fname[n]=mhad2012_970.0_p2.hbook          
n=[n]+1; fname[n]=mhad2012_990.0_p2.hbook          
n=[n]+1; fname[n]=mhad2012_990.0_p1.hbook          
n=[n]+1; fname[n]=mhad2012_980.0_p6.hbook          
n=[n]+1; fname[n]=mhad2012_980.0_p3.hbook          
n=[n]+1; fname[n]=mhad2012_990.0_p4.hbook          
n=[n]+1; fname[n]=mhad2012_990.0_p5.hbook          
n=[n]+1; fname[n]=mhad2012_990.0_p3.hbook          
n=[n]+1; fname[n]=mhad2012_990.0_p6.hbook          
do i=1,[n]
  shell ln -s [dir0]/[fname[i]] [dir1]/[fname[i]]
enddo
mess [n]
*
dir0=/work/users/konctbel/exp/OMEG2012-0/col
dir1=mapcal2
n=0
n=[n]+1; fname[n]=omega2012_389.0_p0.hbook
n=[n]+1; fname[n]=omega2012_389.0_p1.hbook
n=[n]+1; fname[n]=omega2012_389.0_p2.hbook
n=[n]+1; fname[n]=omega2012_389.0_p3.hbook
n=[n]+1; fname[n]=omega2012_389.0_p4.hbook
n=[n]+1; fname[n]=omega2012_389.0_p5.hbook
n=[n]+1; fname[n]=omega2012_389.5_p0.hbook
n=[n]+1; fname[n]=omega2012_390.0_p0.hbook
n=[n]+1; fname[n]=omega2012_390.0_p1.hbook
n=[n]+1; fname[n]=omega2012_390.0_p2.hbook
n=[n]+1; fname[n]=omega2012_390.0_p3.hbook
n=[n]+1; fname[n]=omega2012_391.0_p0.hbook
n=[n]+1; fname[n]=omega2012_391.0_p1.hbook
n=[n]+1; fname[n]=omega2012_391.0_p2.hbook
n=[n]+1; fname[n]=omega2012_395.0_p0.hbook
do i=1,[n]
  shell ln -s [dir0]/[fname[i]] [dir1]/[fname[i]]
enddo
mess [n]
*
dir0=/work/users/konctbel/exp/PHI2013-0/col
dir1=mapcal2
n=0
n=[n]+1; fname[n]=phi2013_507.0_p0.hbook
n=[n]+1; fname[n]=phi2013_512.0_p0.hbook
n=[n]+1; fname[n]=phi2013_512.0_p1.hbook
n=[n]+1; fname[n]=phi2013_512.0_p2.hbook
n=[n]+1; fname[n]=phi2013_508.0_p0.hbook
n=[n]+1; fname[n]=phi2013_514.0_p0.hbook
n=[n]+1; fname[n]=phi2013_508.0_p1.hbook
n=[n]+1; fname[n]=phi2013_514.0_p1.hbook
n=[n]+1; fname[n]=phi2013_509.0_p0.hbook
n=[n]+1; fname[n]=phi2013_525.0_p0.hbook
n=[n]+1; fname[n]=phi2013_509.0_p1.hbook
n=[n]+1; fname[n]=phi2013_525.0_p1.hbook
n=[n]+1; fname[n]=phi2013_525.0_p2.hbook
n=[n]+1; fname[n]=phi2013_509.1_p0.hbook
n=[n]+1; fname[n]=phi2013_525.0_p3.hbook
n=[n]+1; fname[n]=phi2013_509.1_p1.hbook
n=[n]+1; fname[n]=phi2013_525.0_p4.hbook
n=[n]+1; fname[n]=phi2013_525.0_p5.hbook
n=[n]+1; fname[n]=phi2013_509.1_p2.hbook
n=[n]+1; fname[n]=phi2013_510.0_p0.hbook
n=[n]+1; fname[n]=phi2013_510.0_p1.hbook
n=[n]+1; fname[n]=phi2013_511.0_p0.hbook
n=[n]+1; fname[n]=phi2013_511.0_p1.hbook
do i=1,[n]
  shell ln -s [dir0]/[fname[i]] [dir1]/[fname[i]]
enddo
mess [n]
return



macro mapcal3
dir0=/work/users/konctbel/exp/RHO_2012-0/col
dir1=mapcal3
n=0
n=[n]+1; fname[n]=rho2012_250.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_250.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_250.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_250.0_p3.hbook               
n=[n]+1; fname[n]=rho2012_260.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_260.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_260.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_270.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_270.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_270.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_270.0_p3.hbook               
n=[n]+1; fname[n]=rho2012_270.0_p4.hbook               
n=[n]+1; fname[n]=rho2012_270.0_p5.hbook               
n=[n]+1; fname[n]=rho2012_280.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_280.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_280.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_280.0_p3.hbook               
n=[n]+1; fname[n]=rho2012_290.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_290.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_290.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_300.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_300.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_300.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_300.0_p3.hbook               
n=[n]+1; fname[n]=rho2012_310.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_310.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_310.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_320.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_320.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_320.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_330.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_330.0_p10.hbook              
n=[n]+1; fname[n]=rho2012_330.0_p11.hbook              
n=[n]+1; fname[n]=rho2012_330.0_p12.hbook              
n=[n]+1; fname[n]=rho2012_330.0_p13.hbook              
n=[n]+1; fname[n]=rho2012_330.0_p14.hbook              
n=[n]+1; fname[n]=rho2012_330.0_p15.hbook              
n=[n]+1; fname[n]=rho2012_330.0_p16.hbook              
n=[n]+1; fname[n]=rho2012_330.0_p17.hbook              
n=[n]+1; fname[n]=rho2012_330.0_p18.hbook              
n=[n]+1; fname[n]=rho2012_330.0_p19.hbook              
n=[n]+1; fname[n]=rho2012_330.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_330.0_p20.hbook              
n=[n]+1; fname[n]=rho2012_330.0_p21.hbook              
n=[n]+1; fname[n]=rho2012_330.0_p22.hbook              
n=[n]+1; fname[n]=rho2012_330.0_p23.hbook              
n=[n]+1; fname[n]=rho2012_330.0_p24.hbook              
n=[n]+1; fname[n]=rho2012_330.0_p25.hbook              
n=[n]+1; fname[n]=rho2012_330.0_p26.hbook              
n=[n]+1; fname[n]=rho2012_330.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_330.0_p3.hbook               
n=[n]+1; fname[n]=rho2012_330.0_p4.hbook               
n=[n]+1; fname[n]=rho2012_330.0_p5.hbook               
n=[n]+1; fname[n]=rho2012_330.0_p6.hbook               
n=[n]+1; fname[n]=rho2012_330.0_p7.hbook               
n=[n]+1; fname[n]=rho2012_330.0_p8.hbook               
n=[n]+1; fname[n]=rho2012_330.0_p9.hbook               
n=[n]+1; fname[n]=rho2012_340.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_340.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_340.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_340.0_p3.hbook               
n=[n]+1; fname[n]=rho2012_340.0_p4.hbook               
n=[n]+1; fname[n]=rho2012_350.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_350.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_350.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_350.0_p3.hbook               
n=[n]+1; fname[n]=rho2012_360.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_360.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_360.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_360.0_p3.hbook               
n=[n]+1; fname[n]=rho2012_360.0_p4.hbook               
n=[n]+1; fname[n]=rho2012_364.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_364.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_364.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_364.0_p3.hbook               
n=[n]+1; fname[n]=rho2012_368.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_368.0_p10.hbook              
n=[n]+1; fname[n]=rho2012_368.0_p11.hbook              
n=[n]+1; fname[n]=rho2012_368.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_368.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_368.0_p3.hbook               
n=[n]+1; fname[n]=rho2012_368.0_p4.hbook               
n=[n]+1; fname[n]=rho2012_368.0_p5.hbook               
n=[n]+1; fname[n]=rho2012_368.0_p6.hbook               
n=[n]+1; fname[n]=rho2012_368.0_p7.hbook               
n=[n]+1; fname[n]=rho2012_368.0_p8.hbook               
n=[n]+1; fname[n]=rho2012_368.0_p9.hbook               
n=[n]+1; fname[n]=rho2012_370.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_370.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_370.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_372.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_372.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_374.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_374.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_374.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_376.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_376.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_376.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_376.0_p3.hbook               
n=[n]+1; fname[n]=rho2012_378.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_378.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_378.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_380.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_380.0_p10.hbook              
n=[n]+1; fname[n]=rho2012_380.0_p11.hbook              
n=[n]+1; fname[n]=rho2012_380.0_p12.hbook              
n=[n]+1; fname[n]=rho2012_380.0_p13.hbook              
n=[n]+1; fname[n]=rho2012_380.0_p14.hbook              
n=[n]+1; fname[n]=rho2012_380.0_p15.hbook              
n=[n]+1; fname[n]=rho2012_380.0_p16.hbook              
n=[n]+1; fname[n]=rho2012_380.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_380.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_380.0_p3.hbook               
n=[n]+1; fname[n]=rho2012_380.0_p4.hbook               
n=[n]+1; fname[n]=rho2012_380.0_p5.hbook               
n=[n]+1; fname[n]=rho2012_380.0_p6.hbook               
n=[n]+1; fname[n]=rho2012_380.0_p7.hbook               
n=[n]+1; fname[n]=rho2012_380.0_p8.hbook               
n=[n]+1; fname[n]=rho2012_380.0_p9.hbook               
n=[n]+1; fname[n]=rho2012_382.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_382.0_p10.hbook              
n=[n]+1; fname[n]=rho2012_382.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_382.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_382.0_p3.hbook               
n=[n]+1; fname[n]=rho2012_382.0_p4.hbook               
n=[n]+1; fname[n]=rho2012_382.0_p5.hbook               
n=[n]+1; fname[n]=rho2012_382.0_p6.hbook               
n=[n]+1; fname[n]=rho2012_382.0_p7.hbook               
n=[n]+1; fname[n]=rho2012_382.0_p8.hbook               
n=[n]+1; fname[n]=rho2012_382.0_p9.hbook               
n=[n]+1; fname[n]=rho2012_384.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_384.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_384.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_386.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_386.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_387.5_p0.hbook               
n=[n]+1; fname[n]=rho2012_387.5_p10.hbook              
n=[n]+1; fname[n]=rho2012_387.5_p11.hbook              
n=[n]+1; fname[n]=rho2012_387.5_p12.hbook              
n=[n]+1; fname[n]=rho2012_387.5_p13.hbook              
n=[n]+1; fname[n]=rho2012_387.5_p14.hbook              
n=[n]+1; fname[n]=rho2012_387.5_p1.hbook               
n=[n]+1; fname[n]=rho2012_387.5_p2.hbook               
n=[n]+1; fname[n]=rho2012_387.5_p3.hbook               
n=[n]+1; fname[n]=rho2012_387.5_p4.hbook               
n=[n]+1; fname[n]=rho2012_387.5_p5.hbook               
n=[n]+1; fname[n]=rho2012_387.5_p6.hbook               
n=[n]+1; fname[n]=rho2012_387.5_p7.hbook               
n=[n]+1; fname[n]=rho2012_387.5_p8.hbook               
n=[n]+1; fname[n]=rho2012_387.5_p9.hbook               
n=[n]+1; fname[n]=rho2012_388.5_p0.hbook               
n=[n]+1; fname[n]=rho2012_388.5_p1.hbook               
n=[n]+1; fname[n]=rho2012_388.5_p2.hbook               
n=[n]+1; fname[n]=rho2012_389.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_389.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_389.5_p0.hbook               
n=[n]+1; fname[n]=rho2012_389.5_p10.hbook              
n=[n]+1; fname[n]=rho2012_389.5_p11.hbook              
n=[n]+1; fname[n]=rho2012_389.5_p12.hbook              
n=[n]+1; fname[n]=rho2012_389.5_p13.hbook              
n=[n]+1; fname[n]=rho2012_389.5_p14.hbook              
n=[n]+1; fname[n]=rho2012_389.5_p15.hbook              
n=[n]+1; fname[n]=rho2012_389.5_p16.hbook              
n=[n]+1; fname[n]=rho2012_389.5_p1.hbook               
n=[n]+1; fname[n]=rho2012_389.5_p2.hbook               
n=[n]+1; fname[n]=rho2012_389.5_p3.hbook               
n=[n]+1; fname[n]=rho2012_389.5_p4.hbook               
n=[n]+1; fname[n]=rho2012_389.5_p5.hbook               
n=[n]+1; fname[n]=rho2012_389.5_p6.hbook               
n=[n]+1; fname[n]=rho2012_389.5_p7.hbook               
n=[n]+1; fname[n]=rho2012_389.5_p8.hbook               
n=[n]+1; fname[n]=rho2012_389.5_p9.hbook               
n=[n]+1; fname[n]=rho2012_390.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_390.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_390.5_p0.hbook               
n=[n]+1; fname[n]=rho2012_390.5_p1.hbook               
n=[n]+1; fname[n]=rho2012_390.5_p2.hbook               
n=[n]+1; fname[n]=rho2012_391.5_p0.hbook               
n=[n]+1; fname[n]=rho2012_391.5_p10.hbook              
n=[n]+1; fname[n]=rho2012_391.5_p11.hbook              
n=[n]+1; fname[n]=rho2012_391.5_p12.hbook              
n=[n]+1; fname[n]=rho2012_391.5_p13.hbook              
n=[n]+1; fname[n]=rho2012_391.5_p14.hbook              
n=[n]+1; fname[n]=rho2012_391.5_p1.hbook               
n=[n]+1; fname[n]=rho2012_391.5_p2.hbook               
n=[n]+1; fname[n]=rho2012_391.5_p3.hbook               
n=[n]+1; fname[n]=rho2012_391.5_p4.hbook               
n=[n]+1; fname[n]=rho2012_391.5_p5.hbook               
n=[n]+1; fname[n]=rho2012_391.5_p6.hbook               
n=[n]+1; fname[n]=rho2012_391.5_p7.hbook               
n=[n]+1; fname[n]=rho2012_391.5_p8.hbook               
n=[n]+1; fname[n]=rho2012_391.5_p9.hbook               
n=[n]+1; fname[n]=rho2012_393.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_393.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_395.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_395.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_395.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_397.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_397.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_397.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_400.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_400.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_400.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_400.0_p3.hbook               
n=[n]+1; fname[n]=rho2012_409.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_409.0_p10.hbook              
n=[n]+1; fname[n]=rho2012_409.0_p11.hbook              
n=[n]+1; fname[n]=rho2012_409.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_409.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_409.0_p3.hbook               
n=[n]+1; fname[n]=rho2012_409.0_p4.hbook               
n=[n]+1; fname[n]=rho2012_409.0_p5.hbook               
n=[n]+1; fname[n]=rho2012_409.0_p6.hbook               
n=[n]+1; fname[n]=rho2012_409.0_p7.hbook               
n=[n]+1; fname[n]=rho2012_409.0_p8.hbook               
n=[n]+1; fname[n]=rho2012_409.0_p9.hbook               
n=[n]+1; fname[n]=rho2012_420.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_420.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_420.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_420.0_p3.hbook               
n=[n]+1; fname[n]=rho2012_420.0_p4.hbook               
n=[n]+1; fname[n]=rho2012_420.0_p5.hbook               
n=[n]+1; fname[n]=rho2012_420.0_p6.hbook               
n=[n]+1; fname[n]=rho2012_420.0_p7.hbook               
n=[n]+1; fname[n]=rho2012_430.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_430.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_430.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_430.0_p3.hbook               
n=[n]+1; fname[n]=rho2012_430.0_p4.hbook               
n=[n]+1; fname[n]=rho2012_440.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_440.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_440.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_450.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_450.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_450.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_460.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_460.0_p10.hbook              
n=[n]+1; fname[n]=rho2012_460.0_p11.hbook              
n=[n]+1; fname[n]=rho2012_460.0_p12.hbook              
n=[n]+1; fname[n]=rho2012_460.0_p13.hbook              
n=[n]+1; fname[n]=rho2012_460.0_p14.hbook              
n=[n]+1; fname[n]=rho2012_460.0_p15.hbook              
n=[n]+1; fname[n]=rho2012_460.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_460.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_460.0_p3.hbook               
n=[n]+1; fname[n]=rho2012_460.0_p4.hbook               
n=[n]+1; fname[n]=rho2012_460.0_p5.hbook               
n=[n]+1; fname[n]=rho2012_460.0_p6.hbook               
n=[n]+1; fname[n]=rho2012_460.0_p7.hbook               
n=[n]+1; fname[n]=rho2012_460.0_p8.hbook               
n=[n]+1; fname[n]=rho2012_460.0_p9.hbook               
n=[n]+1; fname[n]=rho2012_470.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_470.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_470.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_470.0_p3.hbook               
n=[n]+1; fname[n]=rho2012_470.0_p4.hbook               
n=[n]+1; fname[n]=rho2012_480.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_480.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_480.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_480.0_p3.hbook               
n=[n]+1; fname[n]=rho2012_480.0_p4.hbook               
n=[n]+1; fname[n]=rho2012_480.0_p5.hbook               
n=[n]+1; fname[n]=rho2012_490.0_p0.hbook               
n=[n]+1; fname[n]=rho2012_490.0_p1.hbook               
n=[n]+1; fname[n]=rho2012_490.0_p2.hbook               
n=[n]+1; fname[n]=rho2012_490.0_p3.hbook               
n=[n]+1; fname[n]=rho2012_490.0_p4.hbook               
n=[n]+1; fname[n]=rho2012_490.0_p5.hbook               
do i=1,[n]
  shell ln -s [dir0]/[fname[i]] [dir1]/[fname[i]]
enddo
mess [n]
*
dir0=/work/users/konctbel/exp/OMEG2012-0/col
dir1=mapcal3
n=0
n=[n]+1; fname[n]=omega2012_385.0_p0.hbook
n=[n]+1; fname[n]=omega2012_385.0_p1.hbook
n=[n]+1; fname[n]=omega2012_385.0_p2.hbook
n=[n]+1; fname[n]=omega2012_385.0_p3.hbook
n=[n]+1; fname[n]=omega2012_387.0_p0.hbook
n=[n]+1; fname[n]=omega2012_387.0_p1.hbook
n=[n]+1; fname[n]=omega2012_387.0_p2.hbook
n=[n]+1; fname[n]=omega2012_387.0_p3.hbook
n=[n]+1; fname[n]=omega2012_388.0_p0.hbook
n=[n]+1; fname[n]=omega2012_388.0_p1.hbook
n=[n]+1; fname[n]=omega2012_388.0_p2.hbook
n=[n]+1; fname[n]=omega2012_388.0_p3.hbook
n=[n]+1; fname[n]=omega2012_388.0_p4.hbook
n=[n]+1; fname[n]=omega2012_392.0_p0.hbook
n=[n]+1; fname[n]=omega2012_392.0_p1.hbook
n=[n]+1; fname[n]=omega2012_392.0_p2.hbook
n=[n]+1; fname[n]=omega2012_394.0_p0.hbook
n=[n]+1; fname[n]=omega2012_394.0_p1.hbook
n=[n]+1; fname[n]=omega2012_394.0_p2.hbook
n=[n]+1; fname[n]=omega2012_394.0_p3.hbook
n=[n]+1; fname[n]=omega2012_394.0_p4.hbook
do i=1,[n]
  shell ln -s [dir0]/[fname[i]] [dir1]/[fname[i]]
enddo
mess [n]
return



macro mapcal4
dir0=/work/users/konctbel/exp/RHO_2012-0/col
dir1=mapcal4
n=0
n=[n]+1; fname[n]=rho2012_160.0_p0.hbook             
n=[n]+1; fname[n]=rho2012_160.0_p10.hbook            
n=[n]+1; fname[n]=rho2012_160.0_p11.hbook            
n=[n]+1; fname[n]=rho2012_160.0_p12.hbook            
n=[n]+1; fname[n]=rho2012_160.0_p13.hbook            
n=[n]+1; fname[n]=rho2012_160.0_p14.hbook            
n=[n]+1; fname[n]=rho2012_160.0_p15.hbook            
n=[n]+1; fname[n]=rho2012_160.0_p16.hbook            
n=[n]+1; fname[n]=rho2012_160.0_p17.hbook            
n=[n]+1; fname[n]=rho2012_160.0_p18.hbook            
n=[n]+1; fname[n]=rho2012_160.0_p1.hbook             
n=[n]+1; fname[n]=rho2012_160.0_p2.hbook             
n=[n]+1; fname[n]=rho2012_160.0_p3.hbook             
n=[n]+1; fname[n]=rho2012_160.0_p4.hbook             
n=[n]+1; fname[n]=rho2012_160.0_p5.hbook             
n=[n]+1; fname[n]=rho2012_160.0_p6.hbook             
n=[n]+1; fname[n]=rho2012_160.0_p7.hbook             
n=[n]+1; fname[n]=rho2012_160.0_p8.hbook             
n=[n]+1; fname[n]=rho2012_160.0_p9.hbook             
n=[n]+1; fname[n]=rho2012_170.0_p0.hbook             
n=[n]+1; fname[n]=rho2012_170.0_p10.hbook            
n=[n]+1; fname[n]=rho2012_170.0_p11.hbook            
n=[n]+1; fname[n]=rho2012_170.0_p12.hbook            
n=[n]+1; fname[n]=rho2012_170.0_p13.hbook            
n=[n]+1; fname[n]=rho2012_170.0_p14.hbook            
n=[n]+1; fname[n]=rho2012_170.0_p15.hbook            
n=[n]+1; fname[n]=rho2012_170.0_p1.hbook             
n=[n]+1; fname[n]=rho2012_170.0_p2.hbook             
n=[n]+1; fname[n]=rho2012_170.0_p3.hbook             
n=[n]+1; fname[n]=rho2012_170.0_p4.hbook             
n=[n]+1; fname[n]=rho2012_170.0_p5.hbook             
n=[n]+1; fname[n]=rho2012_170.0_p6.hbook             
n=[n]+1; fname[n]=rho2012_170.0_p7.hbook             
n=[n]+1; fname[n]=rho2012_170.0_p8.hbook             
n=[n]+1; fname[n]=rho2012_170.0_p9.hbook             
n=[n]+1; fname[n]=rho2012_180.0_p0.hbook             
n=[n]+1; fname[n]=rho2012_180.0_p10.hbook            
n=[n]+1; fname[n]=rho2012_180.0_p11.hbook            
n=[n]+1; fname[n]=rho2012_180.0_p12.hbook            
n=[n]+1; fname[n]=rho2012_180.0_p13.hbook            
n=[n]+1; fname[n]=rho2012_180.0_p14.hbook            
n=[n]+1; fname[n]=rho2012_180.0_p15.hbook            
n=[n]+1; fname[n]=rho2012_180.0_p16.hbook            
n=[n]+1; fname[n]=rho2012_180.0_p17.hbook            
n=[n]+1; fname[n]=rho2012_180.0_p18.hbook            
n=[n]+1; fname[n]=rho2012_180.0_p19.hbook            
n=[n]+1; fname[n]=rho2012_180.0_p1.hbook             
n=[n]+1; fname[n]=rho2012_180.0_p20.hbook            
n=[n]+1; fname[n]=rho2012_180.0_p2.hbook             
n=[n]+1; fname[n]=rho2012_180.0_p3.hbook             
n=[n]+1; fname[n]=rho2012_180.0_p4.hbook             
n=[n]+1; fname[n]=rho2012_180.0_p5.hbook             
n=[n]+1; fname[n]=rho2012_180.0_p6.hbook             
n=[n]+1; fname[n]=rho2012_180.0_p7.hbook             
n=[n]+1; fname[n]=rho2012_180.0_p8.hbook             
n=[n]+1; fname[n]=rho2012_180.0_p9.hbook             
n=[n]+1; fname[n]=rho2012_190.0_p0.hbook             
n=[n]+1; fname[n]=rho2012_190.0_p1.hbook             
n=[n]+1; fname[n]=rho2012_190.0_p2.hbook             
n=[n]+1; fname[n]=rho2012_190.0_p3.hbook             
n=[n]+1; fname[n]=rho2012_190.0_p4.hbook             
n=[n]+1; fname[n]=rho2012_200.0_p0.hbook             
n=[n]+1; fname[n]=rho2012_200.0_p1.hbook             
n=[n]+1; fname[n]=rho2012_200.0_p2.hbook             
n=[n]+1; fname[n]=rho2012_200.0_p3.hbook             
n=[n]+1; fname[n]=rho2012_200.0_p4.hbook             
n=[n]+1; fname[n]=rho2012_200.0_p5.hbook             
n=[n]+1; fname[n]=rho2012_200.0_p6.hbook             
n=[n]+1; fname[n]=rho2012_210.0_p0.hbook             
n=[n]+1; fname[n]=rho2012_210.0_p1.hbook             
n=[n]+1; fname[n]=rho2012_210.0_p2.hbook             
n=[n]+1; fname[n]=rho2012_210.0_p3.hbook             
n=[n]+1; fname[n]=rho2012_210.0_p4.hbook             
n=[n]+1; fname[n]=rho2012_210.0_p5.hbook             
n=[n]+1; fname[n]=rho2012_220.0_p0.hbook             
n=[n]+1; fname[n]=rho2012_220.0_p1.hbook             
n=[n]+1; fname[n]=rho2012_220.0_p2.hbook             
n=[n]+1; fname[n]=rho2012_220.0_p3.hbook             
n=[n]+1; fname[n]=rho2012_220.0_p4.hbook             
n=[n]+1; fname[n]=rho2012_220.0_p5.hbook             
n=[n]+1; fname[n]=rho2012_220.0_p6.hbook             
n=[n]+1; fname[n]=rho2012_220.0_p7.hbook             
n=[n]+1; fname[n]=rho2012_230.0_p0.hbook             
n=[n]+1; fname[n]=rho2012_230.0_p1.hbook             
n=[n]+1; fname[n]=rho2012_230.0_p2.hbook             
n=[n]+1; fname[n]=rho2012_230.0_p3.hbook             
n=[n]+1; fname[n]=rho2012_230.0_p4.hbook             
n=[n]+1; fname[n]=rho2012_230.0_p5.hbook             
n=[n]+1; fname[n]=rho2012_230.0_p6.hbook             
n=[n]+1; fname[n]=rho2012_230.0_p7.hbook             
n=[n]+1; fname[n]=rho2012_240.0_p0.hbook             
n=[n]+1; fname[n]=rho2012_240.0_p1.hbook             
n=[n]+1; fname[n]=rho2012_240.0_p2.hbook             
n=[n]+1; fname[n]=rho2012_240.0_p3.hbook             
n=[n]+1; fname[n]=rho2012_240.0_p4.hbook             
n=[n]+1; fname[n]=rho2012_240.0_p5.hbook             
n=[n]+1; fname[n]=rho2012_240.0_p6.hbook             
n=[n]+1; fname[n]=rho2012_240.0_p7.hbook             
n=[n]+1; fname[n]=rho2012_240.0_p8.hbook             
n=[n]+1; fname[n]=rho2012_250.0_p10.hbook            
n=[n]+1; fname[n]=rho2012_250.0_p11.hbook            
n=[n]+1; fname[n]=rho2012_250.0_p12.hbook            
n=[n]+1; fname[n]=rho2012_250.0_p13.hbook            
n=[n]+1; fname[n]=rho2012_250.0_p14.hbook            
n=[n]+1; fname[n]=rho2012_250.0_p15.hbook            
n=[n]+1; fname[n]=rho2012_250.0_p16.hbook            
n=[n]+1; fname[n]=rho2012_250.0_p5.hbook             
n=[n]+1; fname[n]=rho2012_250.0_p6.hbook             
n=[n]+1; fname[n]=rho2012_250.0_p7.hbook             
n=[n]+1; fname[n]=rho2012_250.0_p8.hbook             
n=[n]+1; fname[n]=rho2012_250.0_p9.hbook             
do i=1,[n]
  shell ln -s [dir0]/[fname[i]] [dir1]/[fname[i]]
enddo
mess [n]
return



macro pdsposx i1=7 i2=17 pdsx=raw
do i=[i1],[i2]
  fname=pdspos[i].kumac
  if ($fexist([fname]).eq.1) then
    shell rm [fname]
  endif
  for/file 20 [fname] new
  close 20
  txt=exec mapcal#pdspos [i] [pdsx]
  fmess [txt] [fname]
*  shell .testrelease/.mainrelease/Offline/submit.sh -q clusters,180 pawbigX11 -b [fname]
  shell paw -b [fname]
enddo
return





macro pdspos part=7 pdsx=raw
*
if ([pdsx].eq.'cal') then
  idh0=10000
  sig0=0.1
endif
if ([pdsx].eq.'raw') then
  idh0=20000
  sig0=5
endif
*
nmax=1000
ve/cre runs([nmax]) r
do nc=1,9
  ve/cre mean[nc]([nmax]) r
  ve/cre rms[nc]([nmax]) r
  ve/cre nevt[nc]([nmax]) r
  ve/cre xp[nc]([nmax]) r
  ve/cre dxp[nc]([nmax]) r
  ve/cre sp[nc]([nmax]) r
  ve/cre dsp[nc]([nmax]) r
enddo
dir=/work/users/konctbel/MinuitTest/v2
dir1=/work/users/kladov/MinuitTest/v2
ind=0
ve/cre p3(3) r 1000 0 [sig0]
ve/cre dp3(3) r
ve/cre s3(3) r 10 1 $sigma([sig0]/10)
ve/cre pmin(3) r 0 0 0
ve/cre pmax(3) r 100000 0 0
do i=$sigma([part]*1000),$sigma([part]*1000+999)
  fname=[dir]/run_[i]_spects_v2.his
  if ($fexist([fname])) then
    do nc=1,9
      idh=[idh0]+[nc]
      if ($hexist([idh])) then
        hi/del [idh]
      endif
    enddo
    hi/file 20 [fname]
    ind=[ind]+1
    ve/inp runs([ind]) [i]
    do nc=1,9
      idh=[idh0]+[nc]
      hrin [idh]
      meani=$hinfo([idh],'mean')
      rmsi=$hinfo([idh],'rms')
      nevti=$hinfo([idh],'events')
      nx=$hinfo([idh],'xbins')
      xmin=$hinfo([idh],'xmin')
      xmax=$hinfo([idh],'xmax')
      ve/cre ic([nx]) r
      hi/get/cont [idh] ic
      ve/cre ix([nx]) r
      sigma ix=array([nx],1#[nx])
      sigma ix=order(ix,-ic)
      sigma ic=order(ic,-ic)
      ve/inp mean[nc]([ind]) [meani]
      ve/inp rms[nc]([ind]) [rmsi]
      ve/inp nevt[nc]([ind]) [nevti]
      ve/inp p3(1) $sigma(0.2*[nevti])
      ve/inp p3(2) $sigma([xmin]+(ix(1)-0.5)*([xmax]-([xmin]))/[nx])
      ve/inp p3(3) [rmsi]
      ve/inp pmin(2) $sigma([meani]-[rmsi])
      ve/inp pmax(2) $sigma([meani]+[rmsi])
      ve/inp pmin(3) $sigma(max([rmsi]/100,1))
      ve/inp pmax(3) $sigma([rmsi])
      if ($sigma(ic(1)).gt.1) then
        ve/copy s3 s3t
        ve/inp s3t(2) 0
        l=$sigma(p3(2)-2*abs(p3(3)))
        r=$sigma(p3(2)+3*abs(p3(3)))
        hi/fit [idh]([l]:[r]) g b 3 p3 s3t pmin pmax dp3
        ve/copy s3 s3t
        ve/inp s3t(1) 0
        ve/inp s3t(3) 0
        l=$sigma(p3(2)-2*abs(p3(3)))
        r=$sigma(p3(2)+3*abs(p3(3)))
        hi/fit [idh]([l]:[r]) g b 3 p3 s3t pmin pmax dp3
        do j=1,3
          l=$sigma(p3(2)-2*abs(p3(3)))
          r=$sigma(p3(2)+3*abs(p3(3)))
          hi/fit [idh]([l]:[r]) g b 3 p3 s3 pmin pmax dp3
        enddo
      endif
      ve/inp  xp[nc]([ind]) $sigma(p3(2))
      ve/inp dxp[nc]([ind]) $sigma(dp3(2))
      ve/inp  sp[nc]([ind]) $sigma(abs(p3(3)))
      ve/inp dsp[nc]([ind]) $sigma(dp3(3))
*      l=$sigma(p3(2)-10*abs(p3(3)))
*      r=$sigma(p3(2)+10*abs(p3(3)))
*      hi/pl [idh]([l]:[r])
*      read x
    enddo
    close 20
  endif
enddo
exec mapcal#vecut runs [ind]
do nc=1,9
  exec mapcal#vecut mean[nc] [ind]
  exec mapcal#vecut rms[nc] [ind]
  exec mapcal#vecut nevt[nc] [ind]
  exec mapcal#vecut xp[nc] [ind]
  exec mapcal#vecut dxp[nc] [ind]
  exec mapcal#vecut sp[nc] [ind]
  exec mapcal#vecut dsp[nc] [ind]
  ve/write runs,mean[nc],rms[nc],nevt[nc],xp[nc],dxp[nc],sp[nc],dsp[nc] [dir1]/pds_position_p[part]_[pdsx]_counter[nc].txt '(8f15.6)'
enddo
return


macro pdsposb pdsx=raw
*
if ([pdsx].eq.'cal') then
  idh0=10000
  sig0=0.1
endif
if ([pdsx].eq.'raw') then
  idh0=20000
  sig0=10
endif
idhs0=40000
*
ve/del runs0,days0,beams0
ve/read runs0,days0,beams0 runpars.txt 3f15.6
rmin=$sigma(vmin(runs0))
rmax=$sigma(vmax(runs0))
nrun=[rmax]-[rmin]+1
ve/cre runi([nrun]) r
ve/cre dayi([nrun]) r
ve/cre beami([nrun]) r
do i=1,$vlen(runs0)
  ri=runs0([i])
  di=days0([i])
  bi=beams0([i])
  ind=[ri]-[rmin]+1
  ve/inp runi([ind]) [ri]
  ve/inp dayi([ind]) [di]
  ve/inp beami([ind]) [bi]
enddo
*
nmax=1000
ve/cre  runs([nmax]) r
ve/cre druns([nmax]) r
ve/cre  days([nmax]) r
ve/cre ddays([nmax]) r
ve/cre  beams([nmax]) r
ve/cre dbeams([nmax]) r
do nc=1,9
  ve/cre mean[nc]([nmax]) r
  ve/cre rms[nc]([nmax]) r
  ve/cre nevt[nc]([nmax]) r
  ve/cre xp[nc]([nmax]) r
  ve/cre dxp[nc]([nmax]) r
  ve/cre sp[nc]([nmax]) r
  ve/cre dsp[nc]([nmax]) r
enddo
dir=v2
ind=0
ve/cre p3(3) r 1000 0 [sig0]
ve/cre dp3(3) r
ve/cre s3(3) r 10 0.001 $sigma([sig0]/10)
ve/cre pmin(3) r 0 0 0
ve/cre pmax(3) r 100000 0 0
*
fname=pdspos/pds_beams_[pdsx]_v0.tex
if ($fexist([fname]).eq.1) then
  shell rm [fname]
endif
for/file  20 [fname] new
close 20
*
bo=0
do i=[rmin],[rmax]
  frname=[dir]/run_[i]_spects_v2.his
  if ($fexist([frname])) then
    j=[i]-[rmin]+1
    bi=beami([j])
*    
    if ([bi].ne.[bo]) then
*
      if ([bo].ne.0) then
        ind=[ind]+1
        fmess '\begin{figure}[ht!b]' [fname]
        do nc=1,9
          idhs=[idhs0]+[nc]
*          exec ../SepPar/sp#brn [idhs] 10 100 ! 1
*          idhs=[idhs]+100
          meani=$hinfo([idhs],'mean')
          rmsi=$hinfo([idhs],'rms')
          nevti=$hinfo([idhs],'events')
          ve/inp mean[nc]([ind]) [meani]
          ve/inp rms[nc]([ind]) [rmsi]
          ve/inp nevt[nc]([ind]) [nevti]
          ve/inp p3(1) $sigma(0.2*[nevti])
          ve/inp p3(2) [meani]
          ve/inp p3(3) [sig0]
          ve/inp pmin(2) $sigma([meani]-[rmsi])
          ve/inp pmax(2) $sigma([meani]+[rmsi])
          ve/inp pmin(3) $sigma([rmsi]/100)
          ve/inp pmax(3) $sigma([rmsi])
          do j=1,3
            l=$sigma(p3(2)-2*abs(p3(3)))
            r=$sigma(p3(2)+2*abs(p3(3)))
            hi/fit [idhs]([l]:[r]) g b 3 p3 s3 pmin pmax dp3
          enddo
          ve/inp  xp[nc]([ind]) $sigma(p3(2))
          ve/inp dxp[nc]([ind]) $sigma(dp3(2))
          ve/inp  sp[nc]([ind]) $sigma(abs(p3(3)))
          ve/inp dsp[nc]([ind]) $sigma(dp3(3))
          ve/inp  runs([ind]) $sigma(([ro]+[rb])/2)
          ve/inp druns([ind]) $sigma(([ro]-[rb])/2)
          ve/inp  days([ind]) $sigma(([to]+[tb])/2)
          ve/inp ddays([ind]) $sigma(([to]-[tb])/2)
          ve/inp  beams([ind]) [bo]
          ve/inp dbeams([ind]) 0
*
          l=$sigma(p3(2)-10*abs(p3(3)))
          r=$sigma(p3(2)+10*abs(p3(3)))
          hi/pl [idhs]([l]:[r])
          atitle 'signal, pe' 'events'
          txt=counter [nc]
          exec $PER/s#tf 0.05 0.9 [txt]
          txt=E?beam! = [bo] MeV
          exec $PER/s#tf 0.05 0.8 [txt]
          txt=Runs = {$sigma(int([rb])),$sigma(int([ro]))}
          exec $PER/s#tf 0.05 0.7 [txt]
          epsfile=pds_r$sigma(int([rb]))-$sigma(int([ro]))_counter[nc].eps
          exec save pdspos/[epsfile] f
*          
          fmess '  \begin{minipage}{0.5\textwidth}' [fname]
          fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
          txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
          fmess [txt] [fname]
          fmess '   \end{minipage}' [fname]
*          txt=$unquote('     ')\caption{}
*          fmess [txt] [fname]
*            
*          nx=$hinfo([idhs],'xbins')
*          xmin=$hinfo([idhs],'xmin')
*          xmax=$hinfo([idhs],'xmax')
*          dx=$sigma(([xmax]-[xmin])/[nx]/2)
*          ve/cre ni([nx]) r
*          hi/get/cont [idhs] ni
*          sigma sni = sumv(ni)/[nevti]
*          ir=$sigma(int([nx]/2))
*          while ($sigma(sni([ir])).lt.0.99) do
*            ir=[ir]+1
*          endwhile
*          xir=$sigma([xmin]+[dx]+2*([ir]-1)*[dx])
*          hi/pl [idhs](-0.5:0.5)
*          mess [bo] {[rb],[ro]} {[tb],[to]} {[l]:[r]} {[xir],$sigma(p3(2)+3*abs(p3(3)))}  
*          read x
        enddo
        fmess '\end{figure}' [fname]
        fmess '\clearpage' [fname]
      endif
*    
      do nc=1,9
        idhs=[idhs0]+[nc]
        if ($hexist([idhs])) then
          hi/del [idhs]
        endif
      enddo
*      
      bo=[bi]
      rb=[i]
      tb=dayi([j])
*      
    endif
*    
    ro=[i]
    to=dayi([j])
*    
    do nc=1,9
      idh=[idh0]+[nc]
      if ($hexist([idh])) then
        hi/del [idh]
      endif
    enddo
*      
    hi/file 20 [frname]
    do nc=1,9
      idh=[idh0]+[nc]
      hrin [idh]
      idhs=[idhs0]+[nc]
      if ($hexist([idhs]).eq.0) then
        hi/copy [idh] [idhs]
        hi/op/res [idhs]
      endif
      hi/op/add [idh] [idhs] [idhs]
    enddo
    close 20
*    
  endif
enddo
*
exec mapcal#vecut  runs [ind]
exec mapcal#vecut druns [ind]
exec mapcal#vecut  days [ind]
exec mapcal#vecut ddays [ind]
exec mapcal#vecut  beams [ind]
exec mapcal#vecut dbeams [ind]
do nc=1,9
  exec mapcal#vecut mean[nc] [ind]
  exec mapcal#vecut rms[nc] [ind]
  exec mapcal#vecut nevt[nc] [ind]
  exec mapcal#vecut xp[nc] [ind]
  exec mapcal#vecut dxp[nc] [ind]
  exec mapcal#vecut sp[nc] [ind]
  exec mapcal#vecut dsp[nc] [ind]
  ve/write runs,druns,days,ddays,beams,dbeams,mean[nc],rms[nc],nevt[nc],xp[nc],dxp[nc],sp[nc],dsp[nc] [dir]/pds_position_vs_beam_[pdsx]_counter[nc].txt '(13f15.6)'
enddo
return



macro pdsposc pdsx=raw
*
if ([pdsx].eq.'cal') then
  idh0=10000
  sig0=0.1
endif
if ([pdsx].eq.'raw') then
  idh0=20000
  sig0=10
endif
idhs0=40000
*
fname=accled_amp1pe_work.txt
ve/read clb,run,p1,p2,p3,p4,p5,p6,p7,p8,p9 [fname] '(2f10.1,9f15.6)'
*
rmin=$sigma(vmin(run))
rmax=$sigma(vmax(run))
nrun=[rmax]-[rmin]+1
ve/cre runi([nrun]) r
sigma runi=array([nrun],[rmin]#[rmax])
ve/cre clbi([nrun]) r
n=$vlen(run)
do i=1,$sigma([n]-1)
  ri=run([i])
  j=[i]+1; rj=run([j])
  ci=clb([i])
  indi=[ri]-[rmin]+1
  indj=[rj]-[rmin]+1
  ni=[indj]-[indi]+1
  ve/cre clbx([ni]) r [ni]*[ci]
  ve/copy clbx(1:[ni]) clbi([indi]:[indj])
enddo
*
nmax=1000
ve/cre  runs(1000) r
ve/cre druns(1000) r
ve/cre  clbs(1000) r
ve/cre dclbs(1000) r
do nc=1,9
  ve/cre mean[nc](1000) r
  ve/cre rms[nc](1000) r
  ve/cre nevt[nc](1000) r
  ve/cre xp[nc](1000) r
  ve/cre dxp[nc](1000) r
  ve/cre sp[nc](1000) r
  ve/cre dsp[nc](1000) r
enddo
dir=/work/users/konctbel/MinuitTest/v2
dir1=/work/users/kladov/MinuitTest/v2
ind=0
ve/cre p3(3) r 1000 0 [sig0]
ve/cre dp3(3) r
ve/cre s3(3) r 100 0.01 $sigma([sig0]/10)
ve/cre pmin(3) r 0 0 0
ve/cre pmax(3) r 10000000 0 0
*
fname=pdspos/pds_clbs_[pdsx]_v0.tex
if ($fexist([fname]).eq.1) then
  shell rm [fname]
endif
for/file  20 [fname] new
close 20
*
co=0
do i=[rmin],[rmax]
  frname=[dir]/run_[i]_spects_v2.his
  if ($fexist([frname])) then
    j=[i]-[rmin]+1
    ci=clbi([j])
*    
    if ([ci].ne.[co]) then
*
      if ([co].ne.0) then
	mess co=[co]
        ind=[ind]+1
        fmess '\begin{figure}[ht!b]' [fname]
        do nc=1,9
          idhs=[idhs0]+[nc]
*          exec ../SepPar/sp#brn [idhs] 10 100 ! 1
*          idhs=[idhs]+100
          meani=$hinfo([idhs],'mean')
          rmsi=$hinfo([idhs],'rms')
          nevti=$hinfo([idhs],'events')
          nx=$hinfo([idhs],'xbins')
          xmin=$hinfo([idhs],'xmin')
          xmax=$hinfo([idhs],'xmax')
          ve/cre ic([nx]) r
          hi/get/cont [idhs] ic
          ve/cre ix([nx]) r
          sigma ix=array([nx],1#[nx])
          sigma ix=order(ix,-ic)
          sigma ic=order(ic,-ic)
          ve/inp mean[nc]([ind]) [meani]
          ve/inp rms[nc]([ind]) [rmsi]
          ve/inp nevt[nc]([ind]) [nevti]
          ve/inp p3(1) $sigma(0.2*[nevti])
          ve/inp p3(2) $sigma([xmin]+(ix(1)-0.5)*([xmax]-([xmin]))/[nx])
          ve/inp p3(3) [rmsi]
          ve/inp pmin(2) $sigma([meani]-[rmsi])
          ve/inp pmax(2) $sigma([meani]+[rmsi])
          ve/inp pmin(3) $sigma(max([rmsi]/100,0.01))
          ve/inp pmax(3) $sigma([rmsi])
          if ($sigma(ic(1)).gt.1) then
            ve/inp s3(1) $sigma([nevti]/100)
            ve/copy s3 s3t
            ve/inp s3t(2) 0
            l=$sigma(p3(2)-2*abs(p3(3)))
            r=$sigma(p3(2)+3*abs(p3(3)))
*            hi/fit [idhs]([l]:[r]) g b 3 p3 s3t pmin pmax dp3
*            ve/inp s3(1) $sigma(p3(1)/10)
            ve/copy s3 s3t
            ve/inp s3t(1) 0
            ve/inp s3t(3) 0
            l=$sigma(p3(2)-2*abs(p3(3)))
            r=$sigma(p3(2)+3*abs(p3(3)))
*            hi/fit [idhs]([l]:[r]) g b 3 p3 s3t pmin pmax dp3
*            ve/inp s3(1) $sigma(p3(1)/10)
            if ($sigma(ic(1)).gt.100) then
              do j=1,3
                l=$sigma(p3(2)-2*abs(p3(3)))
                r=$sigma(p3(2)+3*abs(p3(3)))
*                hi/fit [idhs]([l]:[r]) g b 3 p3 s3 pmin pmax dp3
                hi/fit [idhs]([l]:[r]) g ! 3 p3 s3 pmin pmax dp3
              enddo
            endif
          endif
*          l=$sigma(p3(2)-10*abs(p3(3)))
*          r=$sigma(p3(2)+10*abs(p3(3)))
*          hi/pl [idhs]([l]:[r])
*          read x
          ve/inp  xp[nc]([ind]) $sigma(p3(2))
          ve/inp dxp[nc]([ind]) $sigma(dp3(2))
          ve/inp  sp[nc]([ind]) $sigma(abs(p3(3)))
          ve/inp dsp[nc]([ind]) $sigma(dp3(3))
          ve/inp  runs([ind]) $sigma(([ro]+[rb])/2)
          ve/inp druns([ind]) $sigma(([ro]-[rb])/2)
          ve/inp  clbs([ind]) [co]
          ve/inp dclbs([ind]) 0
*
          l=$sigma(p3(2)-10*abs(p3(3)))
          r=$sigma(p3(2)+10*abs(p3(3)))
          hi/pl [idhs]([l]:[r])
          atitle 'signal, pe' 'events'
          txt=counter [nc]
          exec /sweet/home/konctbel/kumac/pl#tf 0.05 0.9 [txt]
          txt=E?beam! = [co] MeV
          exec /sweet/home/konctbel/kumac/pl#tf 0.05 0.8 [txt]
          txt=Runs = {$sigma(int([rb])),$sigma(int([ro]))}
          exec /sweet/home/konctbel/kumac/pl#tf 0.05 0.7 [txt]
          epsfile=pds_r$sigma(int([rb]))-$sigma(int([ro]))_cal_counter[nc].eps
          exec save pdspos/[epsfile] f
*          
          fmess '  \begin{minipage}{0.5\textwidth}' [fname]
          fmess '     {\centering\resizebox*{\textwidth}{!}' [fname]
          txt=$unquote('     '){\includegraphics{[epsfile]}}\par}
          fmess [txt] [fname]
          fmess '   \end{minipage}' [fname]
*          txt=$unquote('     ')\caption{}
*          fmess [txt] [fname]
*            
*          nx=$hinfo([idhs],'xbins')
*          xmin=$hinfo([idhs],'xmin')
*          xmax=$hinfo([idhs],'xmax')
*          dx=$sigma(([xmax]-[xmin])/[nx]/2)
*          ve/cre ni([nx]) r
*          hi/get/cont [idhs] ni
*          sigma sni = sumv(ni)/[nevti]
*          ir=$sigma(int([nx]/2))
*          while ($sigma(sni([ir])).lt.0.99) do
*            ir=[ir]+1
*          endwhile
*          xir=$sigma([xmin]+[dx]+2*([ir]-1)*[dx])
*          hi/pl [idhs](-0.5:0.5)
*          mess [bo] {[rb],[ro]} {[tb],[to]} {[l]:[r]} {[xir],$sigma(p3(2)+3*abs(p3(3)))}  
*          read x
        enddo
        fmess '\end{figure}' [fname]
        fmess '\clearpage' [fname]
      endif
*    
      do nc=1,9
        idhs=[idhs0]+[nc]
        if ($hexist([idhs])) then
          hi/del [idhs]
        endif
      enddo
*      
      co=[ci]
      rb=[i]
*      
    endif
*    
    ro=[i]
*    
    do nc=1,9
      idh=[idh0]+[nc]
      if ($hexist([idh])) then
        hi/del [idh]
      endif
    enddo
*      
    hi/file 20 [frname]
    do nc=1,9
      idh=[idh0]+[nc]
      hrin [idh]
      idhs=[idhs0]+[nc]
      if ($hexist([idhs]).eq.0) then
        hi/copy [idh] [idhs]
        hi/op/res [idhs]
      endif
      hi/op/add [idh] [idhs] [idhs]
    enddo
    close 20
*    
  else
*    mess file [frname] dosn't exist
  endif
enddo
*
mess [ind] $vlen(runs)
exec mapcal#vecut  runs [ind]
exec mapcal#vecut druns [ind]
exec mapcal#vecut  clbs [ind]
exec mapcal#vecut dclbs [ind]
do nc=1,9
  exec mapcal#vecut mean[nc] [ind]
  exec mapcal#vecut rms[nc] [ind]
  exec mapcal#vecut nevt[nc] [ind]
  exec mapcal#vecut xp[nc] [ind]
  exec mapcal#vecut dxp[nc] [ind]
  exec mapcal#vecut sp[nc] [ind]
  exec mapcal#vecut dsp[nc] [ind]
  ve/write runs,druns,clbs,dclbs,mean[nc],rms[nc],nevt[nc],xp[nc],dxp[nc],sp[nc],dsp[nc] [dir1]/pds_position_vs_clb_[pdsx]_counter[nc].txt '(13f15.6)'
enddo
return


macro pdspocuts
dir=/work/users/konctbel/MinuitTest
dir1=v2
pdsx=cal
do nc=1,9
  ve/del runs,druns,clbs,dclbs,mean[nc],rms[nc],nevt[nc],xp[nc],dxp[nc],sp[nc],dsp[nc]
  ve/read runs,druns,clbs,dclbs,mean[nc],rms[nc],nevt[nc],xp[nc],dxp[nc],sp[nc],dsp[nc] [dir1]/pds_position_vs_clb_[pdsx]_counter[nc].txt '(13f15.6)'
  sigma up[nc]=xp[nc]+3*sp[nc]
enddo
sigma rb=runs-druns
*ve/write rb,up1,up2,up3,up4,up5,up6,up7,up8,up9 pds_cuts_vs_clbs.txt '(10f15.6)'
*
ve/del runsp,drunsp,daysp,ddaysp,beamp,dbeamp,nevtp,meanp,dmeanp,rmsp,drmsp,wp,dwp,vp,dvp
ve/read runsp,drunsp,daysp,ddaysp,beamp,dbeamp,nevtp,meanp,dmeanp,rmsp,drmsp,wp,dwp,vp,dvp [dir]/v2/dphi_vs_beams_new.txt 15e15.6
*
n=0
n=[n]+1; fname[n]=[dir]/map1_n1.13.fwi
n=[n]+1; fname[n]=[dir]/map2_n1.13.fwi
n=[n]+1; fname[n]=[dir]/map3_n1.05.fwi
n=[n]+1; fname[n]=[dir]/map4_n1.13.fwi
*
ve/del vrun,vfrun
do f=1,[n]
*
  shell cp [fname[f]] tmp0.txt
  shell $unquote('cat tmp0.txt | sed "s/[^0-9]/ /g" > tmp1.txt')
  ve/del beam,dbeam,frun,nevt
  ve/read beam,dbeam,frun,nevt tmp1.txt
  exec vappend vfrun frun
*
  do i=1,$vlen(frun)
    fakerun=frun([i])
    mess fakerun = [fakerun]
    shell ls /online/gridback/MC/R006-001/ee/output/ee_bhwide_*9979*.psy.gz > tmp.txt
*    if (shell grep "No such" tmp.txt) then
*    else

    shell $unquote('cat tmp.txt | sed "s/\(^.*kb\)/ /g" > tmp1.txt')
    shell $unquote('cat tmp1.txt | sed "s/.mod.gz/ /g" > tmp2.txt')
    shell $unquote('cat tmp2.txt | sed "s/[^0-9]/ /g" > tmp3.txt')
    ve/del ven,vfr,vrn,vnt
    ve/read ven,vfr,vrn,vnt tmp3.txt
    mess [vrn]
    exec vappend vrun vrn
*    endif
  enddo
*  
enddo
*
mess $vlen(vfrun) 
mess $vlen(vrun)
sigma vfrun=order(vfrun,vrun)
sigma vrun=order(vrun,vrun)
*
n=$vlen(vfrun)
do nc=1,9
  ve/cre ups[nc]([n]) r
enddo
*
do i=1,[n]
  frn=vfrun([i])
  exec mapcal#ixndr [frn] rb
  gl/imp inx
  do nc=1,9
    ve/inp ups[nc]([i]) $sigma(up[nc]([inx]))
  enddo
enddo
ve/write vrun,ups1,ups2,ups3,ups4,ups5,ups6,ups7,ups8,ups9 pds_cuts_vs_run_sim.txt '(10f15.6)'
return



macro pdspospl nc=1 pdsx=raw
dir=v2
ve/del runs,mean,rms,nevt,xp,dxp,sp,dsp
ve/read runs,mean,rms,nevt,xp,dxp,sp,dsp [dir]/pds_position_[pdsx]_counter[nc].txt '(8f15.6)'
*
dmax=0.001
sigma dr=abs(dxp-([dmax]))-(dxp-([dmax]))
sigma dr=dr/dr
sigma dxpr=dxp*dr+abs([dmax])*(1-dr)
dmax=0.001
sigma dr=abs(dsp-([dmax]))-(dsp-([dmax]))
sigma dr=dr/dr
sigma dspr=dsp*dr+abs([dmax])*(1-dr)
*
n=$vlen(runs)
ve/cre dn([n]) r 
rmin=$sigma(vmin(runs))
rmax=$sigma(vmax(runs))
l=$sigma([rmin]-0.05*([rmax]-[rmin]))
r=$sigma([rmax]+0.05*([rmax]-[rmin]))
*
v=mean
xm=$sigma(vsum([v])/[n])
xm2=$sigma(vsum([v]**2)/[n])
rm=$sigma(sqrt([xm2]-[xm]**2))
d=$sigma([xm]-2*[rm])
u=$sigma([xm]+2*[rm])
null [l] [r] [d] [u]
sigma dmean=rms/sqrt(nevt)
if ([pdsx].eq.'raw') then
  ve/del runs0,pds1,pds2,pds3,pds4,pds5,pds6,pds7,pds8,pds9
  ve/read runs0,pds1,pds2,pds3,pds4,pds5,pds6,pds7,pds8,pds9 v2/runspdsnew.txt 10f10.3
  set pmci 2
  * exec $PER/s#vpl pds[nc] dn runs0 dn sz=0.05 ll=-1 o=s
endif
set pmci 4
* exec $PER/s#vpl mean dmean runs dn sz=0.05 ll=-1 o=s
read x
*
v=rms
xm=$sigma(vsum([v])/[n])
xm2=$sigma(vsum([v]**2)/[n])
rm=$sigma(sqrt([xm2]-[xm]**2))
d=$sigma([xm]-2*[rm])
u=$sigma([xm]+10*[rm])
null [l] [r] 0 [u]
sigma drms=rms/sqrt(2*nevt)
* exec $PER/s#vpl rms drms runs dn sz=0.05 ll=-1 o=s
read x
*
v=xp
xm=$sigma(vsum([v])/[n])
xm2=$sigma(vsum([v]**2)/[n])
rm=$sigma(sqrt([xm2]-[xm]**2))
d=$sigma([xm]-2*[rm])
u=$sigma([xm]+2*[rm])
null [l] [r] [d] [u]
* exec $PER/s#vpl xp dxpr runs dn sz=0.05 ll=-1 o=s
read x
*
v=sp
xm=$sigma(vsum([v])/[n])
xm2=$sigma(vsum([v]**2)/[n])
rm=$sigma(sqrt([xm2]-[xm]**2))
d=$sigma([xm]-2*[rm])
u=$sigma([xm]+7*[rm])
null [l] [r] 0 [u]
* exec $PER/s#vpl sp dspr runs dn sz=0.05 ll=-1 o=s
read x
*
sigma up = xp+3*sp
v=up
xm=$sigma(vsum([v])/[n])
xm2=$sigma(vsum([v]**2)/[n])
rm=$sigma(sqrt([xm2]-[xm]**2))
d=$sigma([xm]-2*[rm])
u=$sigma([xm]+7*[rm])
null [l] [r] [d] [u]
sigma dup=sqrt(dxpr**2+9*dspr**2)
* exec $PER/s#vpl up dup runs dn sz=0.05 ll=-1 o=s
line [rmin] 0.2 [rmax] 0.2
read x
return



macro runt
exec mapcal#mylib
ve/cre runs(10000) r
ve/cre runc(10000) r
ve/cre prob(10000) r
ve/cre rnevt(10000) r
dir=v2
ind=0
ind1=0
shell rm tmp.txt
do i=6000,18000
  frname=[dir]/run_[i]_spects_v2.his
  if ($fexist([frname])) then
    if ($hexist(1000)) then
      hi/del 1000
    endif
    if ($hexist(100)) then
      hi/del 100
    endif
    hi/file 20 [frname]
    hrin 1000
    hi/copy 1000 100
    close 20
    frname=run_[i]_spects_v2.his
    if ($fexist([frname]).eq.0) then
      ind=[ind]+1
      ve/inp runs([ind]) [i]
      shell fgrep -e [i] /work/users/konctbel/snd2k/R005-999/fwk/*/*fwi > tmpi.txt
      if ([ind].eq.1) then
        shell mv tmpi.txt tmp.txt
      else
        shell cat tmp.txt tmpi.txt > tmpii.txt
        shell mv tmpii.txt tmp.txt
      endif
    else
      ind1=[ind1]+1
      if ($hexist(1000)) then
        hi/del 1000
      endif
      if ($hexist(200)) then
        hi/del 200
      endif
      hi/file 20 [frname]
      hrin 1000
      hi/copy 1000 200
      close 20
      ve/inp prob([ind1]) $call('mhdiff(100,200)')
      ve/inp runc([ind1]) [i]
      nevt1=$hinfo(100,'events')
      nevt2=$hinfo(200,'events')
      if (([nevt1].eq.0).and.([nevt2].eq.0)) then
        ve/inp rnevt([ind1]) 1
      endif
    endif
  endif
enddo
exec mapcal#vecut runs [ind]
exec mapcal#vecut runc [ind1]
exec mapcal#vecut prob [ind1]
exec mapcal#vecut rnevt [ind1]
return


macro fwidiv
dir=/work/users/konctbel/snd2k/R005-999/fwk/MHAD2011
n=0
n=[n]+1; fwi[n]=MHAD2011_750-1_p1
n=[n]+1; fwi[n]=MHAD2011_550_p1
n=[n]+1; fwi[n]=MHAD2011_600_p1
n=[n]+1; fwi[n]=MHAD2011_600_p2
n=[n]+1; fwi[n]=MHAD2011_675_p2
n=[n]+1; fwi[n]=MHAD2011_750_p1
n=[n]+1; fwi[n]=MHAD2011_750_p2
n=[n]+1; fwi[n]=MHAD2011_862.5
n=[n]+1; fwi[n]=MHAD2011_637.5_p2
n=[n]+1; fwi[n]=MHAD2011_612.5_p1
n=[n]+1; fwi[n]=MHAD2011_587.5_p1
n=[n]+1; fwi[n]=MHAD2011_587.5_p2
n=[n]+1; fwi[n]=MHAD2011_562.5_p1
n=[n]+1; fwi[n]=MHAD2011_537.5_p1
do i=1,[n]
  shell fgrep -e "runlist" [dir]/[fwi[i]].fwi >& tmp0.txt
  shell $unquote('cat tmp0.txt | sed "s/runlist = \[/ /g" > tmp1.txt')
  shell $unquote('cat tmp1.txt | sed "s/,/\n/g" > tmp2.txt')
  shell $unquote('cat tmp2.txt | sed "s/\]/ /g" > tmp.txt')
  ve/del runs
  ve/read runs tmp.txt
  nf=$vlen(runs)
  ve/prin runs
  mess [dir]/[fwi[i]].fwi
  read x
*  
  fname=[dir]/[fwi[i]]_a.fwi
  if ($fexist([fname])) then
    shell rm [fname]
  endif
  for/file 20 [fname] new
  close 20
*  
  nd=$sigma(int([nf]/2)-1)
  txt='runlist = ['
  do j=1,[nd]
    txt=$unquote([txt]) $sigma(runs([j])),
  enddo
  nd=[nd]+1
  txt=$unquote([txt]) $sigma(runs([nd])) $unquote(']')
  fmess [txt] [fname]
*  
  fname=[dir]/[fwi[i]]_b.fwi
  if ($fexist([fname])) then
    shell rm [fname]
  endif
  for/file 20 [fname] new
  close 20
*  
  nd=[nd]+1
  nf=[nf]-1
  txt='runlist = ['
  do j=[nd],[nf]
    txt=$unquote([txt]) $sigma(runs([j])),
  enddo
  nf=[nf]+1
  txt=$unquote([txt]) $sigma(runs([nf])) $unquote(']')
  fmess [txt] [fname]
*  
enddo
return

macro readpdsall
ve/cre runs(10000) r
ind=0
dir=v2
do i=6000,18000
  frname=[dir]/run_[i]_spects_v2.his
  if ($fexist([frname])) then
    ind=[ind]+1
    ve/inp runs([ind]) [i]
  endif
enddo
exec mapcal#vecut runs [ind]
*exec mapcal#readpds
ve/del pds1,pds2,pds3,pds4,pds5,pds6,pds7,pds8,pds9
ve/read pds1,pds2,pds3,pds4,pds5,pds6,pds7,pds8,pds9 runspds.txt 10f10.3
ve/write runs,pds1,pds2,pds3,pds4,pds5,pds6,pds7,pds8,pds9 v2/runspdsnew.txt 10f10.3
return


macro pdsposort pdsx=raw
dir=v2
do nc=1,9
  if ([pdsx].eq.'raw') then
    shell cat [dir]/pds_position_p*_[pdsx]_counter[nc].txt > [dir]/pds_position_[pdsx]_counter[nc].txt
  endif
  if ([pdsx].eq.'cal') then
    shell cat [dir]/pds_position_p{[0-9],[0-9][0-9]}_counter[nc].txt > [dir]/pds_position_[pdsx]_counter[nc].txt
  endif
  ve/del runs,mean,rms,nevt,xp,dxp,sp,dsp
  ve/read runs,mean,rms,nevt,xp,dxp,sp,dsp [dir]/pds_position_[pdsx]_counter[nc].txt '(8f15.6)'
  sigma mean=order(mean,runs)
  sigma rms=order(rms,runs)
  sigma nevt=order(nevt,runs)
  sigma xp=order(xp,runs)
  sigma dxp=order(dxp,runs)
  sigma sp=order(sp,runs)
  sigma dsp=order(dsp,runs)
  sigma runs=order(runs,runs)
  ve/write runs,mean,rms,nevt,xp,dxp,sp,dsp [dir]/pds_position_[pdsx]_counter[nc].txt '(8f15.6)'
enddo
return


macro pdsposbeam
dir=v2
*
ve/del runs,pds1,pds2,pds3,pds4,pds5,pds6,pds7,pds8,pds9
ve/read runs,pds1,pds2,pds3,pds4,pds5,pds6,pds7,pds8,pds9 [dir]/runspdsnew.txt 10f10.3
*
ve/del runs0,days,beams
ve/read runs0,days,beams [dir]/runpars.txt 3f15.6
*
ve/cre vnull(1) r 1
exec vappend runs vnull
exec vappend days vnull
exec vappend beams vnull
*
n=$vlen(runs)
ve/cre beamp([n]) r
ve/cre daysp([n]) r
ve/cre ddaysp([n]) r
ve/cre runsp([n]) r
ve/cre drunsp([n]) r
do nc=1,9
  ve/cre  pds[nc]p([n]) r
  ve/cre dpds[nc]p([n]) r
enddo
*
ind=0
bo=0
do i=1,[n]
  b=beams([i])
  if ([b].ne.[bo]) then
    if ([bo].ne.0) then
      i2=[i]-1
      if ([i1].le.[i2]) then
        ind=[ind]+1
        rb=runs0([i1])
        re=runs0([i2])
        exec mapcal#ixndl [rb] runs
        gl/imp inx
        in1=[inx]
        exec mapcal#ixndr [re] runs
        gl/imp inx
        in2=[inx]
        if ([in1].le.[in2]) then
          ve/inp beamp([ind]) [bo]
          t1=days([in1])
          t2=days([in2])
          r1=runs([in1])
          r2=runs([in2])
          ve/inp daysp([ind]) $sigma(([t2]+[t1])/2)
          ve/inp ddaysp([ind]) $sigma(([t2]-[t1])/2)
          ve/inp runsp([ind]) $sigma(([r2]+[r1])/2)
          ve/inp drunsp([ind]) $sigma(([r2]-[r1])/2)
          do nc=1,9
            ve/del pdst
            ve/copy pds[nc]([in1]:[in2]) pdst
            np=[in2]-[in1]+1
            mean=$sigma(vsum(pdst)/[np])
            mean2=$sigma(vsum(pdst**2)/[np])
            ve/inp  pds[nc]p([ind]) [mean]
            ve/inp dpds[nc]p([ind]) $sigma(sqrt([mean2]-[mean]**2))
          enddo
        endif
      endif
    endif
    bo=[b]
    i1=[i]
  endif
enddo
*
exec mapcal#vecut beamp
n=$vlen(beamp)
exec mapcal#vecut  daysp [n]
exec mapcal#vecut ddaysp [n]
exec mapcal#vecut  runsp [n]
exec mapcal#vecut drunsp [n]
do nc=1,9
  exec mapcal#vecut  pds[nc]p [n]
  exec mapcal#vecut dpds[nc]p [n]
enddo
ve/write beamp,daysp,ddaysp,runsp,drunsp,pds1p,pds2p,pds3p,pds4p,pds5p,pds6p,pds7p,pds8p,pds9p [dir]/pds_beam.txt 14f10.3
return


macro pdsposcomp
dir=v2
*
ve/del beamp,daysp,ddaysp,runsp,drunsp,pds1p,pds2p,pds3p,pds4p,pds5p,pds6p,pds7p,pds8p,pds9p
ve/read beamp,daysp,ddaysp,runsp,drunsp,pds1p,pds2p,pds3p,pds4p,pds5p,pds6p,pds7p,pds8p,pds9p [dir]/pds_beam.txt 14f10.3
n=$vlen(beamp)
ve/cre dp([n]) r
*
ve/del runs,pds1,pds2,pds3,pds4,pds5,pds6,pds7,pds8,pds9
ve/read runs,pds1,pds2,pds3,pds4,pds5,pds6,pds7,pds8,pds9 v2/runspdsnew.txt 10f10.3
*
ve/del runm,mn1,mn2,mn3,mn4,mn5,mn6,mn7,mn8,mn9
ve/read runm,mn1,mn2,mn3,mn4,mn5,mn6,mn7,mn8,mn9 v2/runspdscorr.txt 10f10.3
n=$vlen(runm)
ve/cre dm([n]) r
*
pdsx=raw
rmin=$sigma(vmin(runs))
rmax=$sigma(vmax(runs))
l=$sigma([rmin]-0.05*([rmax]-[rmin]))
r=$sigma([rmax]+0.05*([rmax]-[rmin]))
do nc=1,9
  ve/del runs,mean,rms,nevt,xp,dxp,sp,dsp
  ve/read runs,mean,rms,nevt,xp,dxp,sp,dsp [dir]/pds_position_[pdsx]_counter[nc].txt '(8f15.6)'
*
dmax=0.1
sigma dr=abs(dxp-([dmax]))-(dxp-([dmax]))
sigma dr=dr/dr
sigma dxpr=dxp*dr+abs([dmax])*(1-dr)
dmax=0.1
sigma dr=abs(dsp-([dmax]))-(dsp-([dmax]))
sigma dr=dr/dr
sigma dspr=dsp*dr+abs([dmax])*(1-dr)
*
  n=$vlen(runs)
  ve/cre dn([n]) r
  m=$sigma(vsum(mean)/[n])
  u=[m]+100
  d=[m]-100
  null [l] [r] [d] [u]
  sigma dmean = rms/sqrt(nevt)
  set pmci 2
  * exec $PER/s#vpl pds[nc] dn runs dn sz=0.05 o=s
  set pmci 4
*  * exec $PER/s#vpl mean dmean runs dn sz=0.05 o=s
  * exec $PER/s#vpl xp dxpr runs dn sz=0.05 o=s
  set pmci 6
*  * exec $PER/s#vpl mean dmean runs dn sz=0.05 o=s
  * exec $PER/s#vpl mn[nc] dm runm dm sz=0.05 o=s
  set pmci 3
*  * exec $PER/s#vpl mean dmean runs dn sz=0.05 o=s
  * exec $PER/s#vpl pds[nc]p dp runsp drunsp sz=0.05 o=s
  read x
  sigma dpeak = pds[nc] - xp
  ve/pl dpeak
  1d 100 ! 100 -10 10
  ve/hfill dpeak 100
  hi/pl 100
  hi/fit 100(-1.:1.) g s
  read x
enddo
return




macro pdsposprep
dir=v2
ve/del runs,pds1,pds2,pds3,pds4,pds5,pds6,pds7,pds8,pds9
ve/read runs,pds1,pds2,pds3,pds4,pds5,pds6,pds7,pds8,pds9 v2/runspdsnew.txt 10f10.3
pdsx=raw
rmin=$sigma(vmin(runs))
rmax=$sigma(vmax(runs))
l=$sigma([rmin]-0.05*([rmax]-[rmin]))
r=$sigma([rmax]+0.05*([rmax]-[rmin]))
do nc=1,9
  ve/del runs,mean,rms,nevt,xp,dxp,sp,dsp
  ve/read runs,mean,rms,nevt,xp,dxp,sp,dsp [dir]/pds_position_[pdsx]_counter[nc].txt '(8f15.6)'
*
  n=$vlen(runs)
  m=$sigma(vsum(pds[nc])/[n])
  u=[m]+100
  d=[m]-100
  1d 100 ! 10000 [d] [u]
  ve/hfill pds[nc] 100
*  
  nx=$hinfo(100,'xbins')
  xmin=$hinfo(100,'xmin')
  xmax=$hinfo(100,'xmax')
  ve/cre ic([nx]) r
  hi/get/cont 100 ic
  ve/cre ix([nx]) r
  sigma ix=array([nx],1#[nx])
  sigma ix=order(ix,-ic)
  sigma ic=order(ic,-ic)
  mch=$sigma([xmin]+(ix(1)-0.5)*([xmax]-([xmin]))/[nx])
  mess [mch] $sigma(ic(1))
  l=$sigma([mch]-1)
  r=$sigma([mch]+1)
  hi/pl 100([l]:[r])
*  
  dmax=$sigma(([u]-[d])/10000)
  sigma dpds=abs(pds[nc]-([mch]))
  sigma dr=abs(dpds-([dmax]))-(dpds-([dmax]))
  sigma dr=dr+(1-pds[nc]/pds[nc])
  sigma ind[nc]=dr/dr
*  sigma pds[nc]r=pds[nc]*ind[nc]
  ve/del mn[nc]
  ve/copy mean mn[nc]
*  
*  read x
enddo
sigma ind=ind1
do nc=2,9
  sigma ind=ind+ind[nc]
enddo
sigma ind=ind/ind
ve/cre ind0(53) r 53*1
ve/copy ind0(1:53) ind(1:53)
np=$sigma(vsum(ind))
sigma runs=order(runs,-ind)
exec mapcal#vecut runs [np]
do nc=1,9
  sigma mn[nc]=order(mn[nc],-ind)
  exec mapcal#vecut mn[nc] [np]
enddo
n=$vlen(runs)
do i=1,[n]
  mn[nc]i=0
enddo
do i=1,[n]
  do nc=1,9
    mn[nc]x=mn[nc]([i])
    if ([mn[nc]x].eq.0) then
      ve/inp mn[nc]([i]) [mn[nc]i]
    else
      mn[nc]i=[mn[nc]x]
    endif
  enddo
enddo
ve/write runs,mn1,mn2,mn3,mn4,mn5,mn6,mn7,mn8,mn9 v2/runspdscorr.txt 10f10.3
return




macro qsubprep map=1 ver=0 prep=0
*
if ([map].eq.1) then
  fname=map1_n1.13.fwi
endif
if ([map].eq.2) then
  fname=map2_n1.13.fwi
endif
if ([map].eq.3) then
  fname=map3_n1.05.fwi
endif
if ([map].eq.4) then
  fname=map4_n1.13.fwi
endif
*
shell cp [fname] tmp0.txt
shell $unquote('cat tmp0.txt | sed "s/[^0-9]/ /g" > tmp1.txt')
ve/del beam,dbeam,frun,nevt
ve/read beam,dbeam,frun,nevt tmp1.txt
*
ve/cre a1s(9) r 90.043 77.241 105.399 76.985 97.676 114.848 94.784 94.174  98.95
ve/del zv,runf,amp1,amp2,amp3,amp4,amp5,amp6,amp7,amp8,amp9
ve/read zv,runf,amp1,amp2,amp3,amp4,amp5,amp6,amp7,amp8,amp9 accled_amp1pe_mc_v1.txt '(2f10.1,9f15.6)'
ve/cre amp(9) r
do nc=1,9
  ve/inp amp([nc]) $sigma(amp[nc]([map]))
enddo
fa1=mapcal[map]_amp1pe.txt
if ([map].le.-2) then
  ve/write a1s /work/users/kladov/snd2k/R006-003/[fa1]
*  ve/write a1s /work/users/konctbel/snd2k/R005-999/[fa1]
else
  ve/write amp /work/users/kladov/snd2k/R006-003/[fa1]
*  ve/write amp /work/users/konctbel/snd2k/R005-999/[fa1]
endif
*
if ([prep].eq.1) then
*
ve/del run,beam,ts,tf
ve/read run,beam,ts,tf run_energy_stime_ftime.txt '(4f15.6)'
*
dir=v2
msum=days
n=$vlen(frun)
do nc=1,9
  ve/cre a[nc]([n]) r
enddo
ve/cre tx([n]) r
ve/cre xpi(1) r
*
do i=1,[n]
  fakerun=frun([i])
  exec mapcal#ixndl [fakerun] run
  gl/imp inx
  in1=[inx]
  exec mapcal#ixndr [fakerun] run
  gl/imp inx
  in2=[inx]
  if ([in1].eq.[in2]) then
    ti=ts([in1])
    ve/inp tx([i]) [ti]
    mess [i] [ti]
*    
    do nc=1,9
    
      if ([map].eq.1) then
  
        fname=[dir]/mapcal[map]_fit_ampt_vs_[msum]_counter[nc].txt
        ve/read p7,dp7 [fname] '2f15.6'
        fname=[dir]/mapcal[map]_fit_amp_vs_[msum]_counter[nc].txt
        ve/read p3,dp3 [fname] '2f15.6'
        ve/inp p3(2) 0
        ve/inp xpi(1) [ti]
        ai=$call('derf0up.f(xpi)')
        ve/inp a[nc]([i]) [ai]

      endif
*  
      if ([map].ge.2) then
  
        fname=[dir]/mapcal[map]_fit_ampt_vs_[msum]_counter[nc].txt
        ve/read p4,dp4 [fname] '2f15.6'
        fname=[dir]/mapcal[map]_fit_amp_vs_[msum]_counter[nc].txt
        ve/read p4a,dp4a [fname] '2f15.6'
        ve/inp p4a(2) 0
        ve/inp p4a(4) 1
        ve/inp xpi(1) [ti]
        ai=$call('erf0up.f(xpi)')
        ve/inp a[nc]([i]) [ai]
      
      endif
    
    enddo

  endif
enddo
*
ve/cre am(9) r
do nc=1,9
  ami=$sigma(vsum(a[nc]*nevt)/vsum(nevt))
  ve/cre qe[nc]([n]) r
  sigma qe[nc] = a[nc]/[ami]
  ve/inp am([nc]) [ami]
enddo
*
ve/cre dn([n]) r
do nc=1,9
  * exec $PER/s#vpl a[nc] dn tx dn sz=0.05
  if ([map].eq.1) then
    fname=[dir]/mapcal[map]_fit_ampt_vs_[msum]_counter[nc].txt
    ve/read p7,dp7 [fname] '2f15.6'
    fname=[dir]/mapcal[map]_fit_amp_vs_[msum]_counter[nc].txt
    ve/read p3,dp3 [fname] '2f15.6'
    ve/inp p3(2) 0
    fun/pl derf0up.f $GRAFINFO('WNXMIN') $GRAFINFO('WNXMAX') s
  else
    fname=[dir]/mapcal[map]_fit_ampt_vs_[msum]_counter[nc].txt
    ve/read p4,dp4 [fname] '2f15.6'
    fname=[dir]/mapcal[map]_fit_amp_vs_[msum]_counter[nc].txt
    ve/read p4a,dp4a [fname] '2f15.6'
    ve/inp p4a(2) 0
    ve/inp p4a(4) 1
    fun/pl erf0up.f $GRAFINFO('WNXMIN') $GRAFINFO('WNXMAX') s
  endif
  ami=am([nc])
  line $GRAFINFO('WNXMIN') [ami] $GRAFINFO('WNXMAX') [ami]
  read x
enddo
*
ve/cre qe(9) r
do i=1,[n]
  fakerun=frun([i])
  fqe=mapcal[map]_v0_qe_f[fakerun].txt
  do nc=1,9
    ve/inp qe([nc]) $sigma(qe[nc]([i]))
  enddo
  ve/write qe /work/users/kladov/snd2k/R006-003/[fqe]
*  ve/write qe /work/users/konctbel/snd2k/R005-999/[fqe]
  mess [fqe] wrote
enddo  
*fqe=mapcal[map]_v[ver]_qe.txt
*ve/write frun,qe1,qe2,qe3,qe4,qe5,qe6,qe7,qe8,qe9 /work/users/konctbel/snd2k/R005-999/[fqe]  '(10f10.6)'
endif
*
cmd=.mainrelease/Offline/submit.sh -q clusters,180
fmap=mapcal[map]_v[ver].cal.v
dir=/online/gridspool/MC/bhabha/

hdir=/work/users/konctbel/MinuitTest/mapcal[map]/sim_v[ver]/
hdir1=/work/users/kladov/MinuitTest/mapcal[map]/sim_v[ver]/
if ($fexist([hdir1]).eq.0) then
  shell mkdir [hdir1]
endif
*
shell rm tmp*.sh
do i=1,$vlen(frun)
  fakerun=frun([i])
  fakerunx=[fakerun]
  if ([fakerun].gt.14044) then
*    fakerunx=14044
  endif
  fqe=mapcal[map]_v0_qe_f[fakerun].txt
  shell ls /online/gridspool/MC/bhabha/ee_bhwide_*[fakerun]*.mod.gz > tmp.txt
  shell $unquote('cat tmp.txt | sed "s/\/online\/gridspool\/MC\/bhabha\// zzz/g" > tmp1.txt')
  shell $unquote('cat tmp1.txt | sed "s/ zzz/=/g" > tmp2.txt')
  shell $unquote('cat tmp2.txt | sed "s/.mod.gz/ /g" > tmp3.txt')
*  shell echo [cmd] FAKERUN=[fakerunx] CLBMAP=[fmap] CLBQE=[fqe] CLBA1=[fa1] MODFILENAME`cat tmp3.txt` MODFILEDIR=[dir] HBOOKDIR=[hdir1] SimRecApp fwk/simreco_col_point.fwk > tmp[i].sh
  shell echo [cmd] FAKERUN=[fakerunx] MODFILENAME`cat tmp3.txt` MODFILEDIR=[dir] HBOOKDIR=[hdir1] SimRecApp fwk/simrecapp_col_point.fwk > tmp[i].sh
enddo
shell cat tmp*.sh > /work/users/konctbel/snd2k/R005-999/mapcal[map]_v[ver].sh
shell chmod +x /work/users/konctbel/snd2k/R005-999/mapcal[map]_v[ver].sh
return



macro accqeprep
*
ve/del run,beam,ts,tf
ve/read run,beam,ts,tf run_energy_stime_ftime.txt '(4f15.6)'
*
ve/del rb
do nc=1,9
  ve/del a[nc]
enddo
*
do map=1,4
*
if ([map].eq.1) then
  fname=map1_n1.13.fwi
endif
if ([map].eq.2) then
  fname=map2_n1.13.fwi
endif
if ([map].eq.3) then
  fname=map3_n1.05.fwi
endif
if ([map].eq.4) then
  fname=map4_n1.13.fwi
endif
*
shell cp [fname] tmp0.txt
shell $unquote('cat tmp0.txt | sed "s/[^0-9]/ /g" > tmp1.txt')
ve/del beami,dbeam,frun,nevt
ve/read beami,dbeam,frun,nevt tmp1.txt
*
n=$vlen(frun)
do nc=1,9
  ve/cre am[nc]([n]) r
enddo
ve/cre rbm([n]) r
*
ve/del qe
do i=1,[n]
  fakerun=frun([i])
*
  shell rm tmp*.txt
  shell fgrep -e [fakerun] /work/users/konctbel/snd2k/R005-999/fwk/*/*.* > tmp0.txt
  shell $unquote('cat tmp0.txt | sed "s/\(^.*=\)/ /g" > tmp1.txt')
  shell $unquote('cat tmp1.txt | sed "s/,/\n/g" > tmp2.txt')
  shell $unquote('cat tmp2.txt | sed "s/[^0-9]/ /g" > tmp3.txt')
  shell ls -l tmp3.txt > tmp4.txt
  shell $unquote('cat tmp4.txt | sed "s/\(^.*konctbel\)/ /g" > tmp5.txt')
  shell $unquote('cat tmp5.txt | sed "s/\(Dec.*$\)/ /g" > tmp6.txt')
  ve/read lf tmp6.txt
  if ($sigma(lf(1)).ne.0) then
    ve/del rp
    ve/read rp tmp3.txt
    ri=$sigma(vmin(rp))
    mess [ri] [fakerun]
  else
    ri=[fakerun]
  endif
*
  ve/read qe /work/users/konctbel/snd2k/R005-999/mapcal[map]_v0_qe_f[fakerun].txt
  ve/inp rbm([i]) [ri]
  do nc=1,9
    ve/inp am[nc]([i]) $sigma(qe([nc]))
  enddo
*  read x
enddo
ve/write rb,a1,a2,a3,a4,a5,a6,a7,a8,a9 ! 10f15.6
exec vappend rb rbm
do nc=1,9
  exec vappend a[nc] am[nc]
enddo
enddo
do nc=1,9
  sigma a[nc]=order(a[nc],rb)
enddo
sigma rb=order(rb,rb)
ve/write rb,a1,a2,a3,a4,a5,a6,a7,a8,a9 accqe_vs_beams.txt 10f15.6
return
