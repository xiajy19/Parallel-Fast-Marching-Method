!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM slice
IMPLICIT NONE
INTEGER :: i,j,k,l,m,n,i1,j1,k1,extrds,extrnss,extrews,idp,iew,ins
INTEGER :: checkstat,isw,conp,cont,conr,id1,id2
INTEGER :: nnr,nnt,nnp,rnnr,rnnt,rnnp,pvgd,pvgns,pvgew,prpd,prpns,prpew,nre,nr
INTEGER :: idm1,idm2,idm3,nnx,nnz,stp,stt,str
INTEGER :: ptfd,ptfns,ptfew,rnode,ri,xi,zi
INTEGER :: ddt,ddp,dnsr,dnst,dewr,dewp
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
INTEGER, PARAMETER :: i5=SELECTED_REAL_KIND(5,10)
REAL(KIND=i5) :: rdep,rlat,rlon
REAL(KIND=i10) :: gor,got,gop,rgor,rgot,rgop,rgsr,rgst,rgsp,rrgsr,rrgst,rrgsp,sldep,slns,slew
REAL(KIND=i10) :: tt,ttt,ttb,sumi,sumj,sumk
REAL(KIND=i10) :: lft,rgt,btm,top
REAL(KIND=i10), SAVE :: gox,goz,dvx,dvz,rgox,rgoz,dnx,dnz,distz,distx,distr
REAL(KIND=i10) :: rdm,rd1,rd2,rd3,u,v,w
REAL(KIND=i10), DIMENSION(:) :: wi(4)
REAL(KIND=i10), DIMENSION (4) :: ai,bi,ci
REAL(KIND=i10), DIMENSION(:,:), ALLOCATABLE :: ui,vi,vela
REAL(KIND=i10), DIMENSION(:,:,:), ALLOCATABLE :: ttn,veln,velv
CHARACTER (LEN=30) :: ifilet,ifilev,irfile,ifilevr,sep
CHARACTER (LEN=30) :: ofiledb,ofilensb,ofileewb
CHARACTER (LEN=30) :: ofiledv,ofilensv,ofileewv
CHARACTER (LEN=30) :: ofiledt,ofilenst,ofileewt
CHARACTER (LEN=30) :: ofiledr,ofilensr,ofileewr
REAL(KIND=i10), PARAMETER :: pi=3.1415926535898
!
! sldep = slice depth
! slns = NS slice
! slew = EW slice
! ifilet = input 3-D traveltime grid file
! ifilev = input velocity grid file
! ifilevr = input reference velocity grid file
! irfile = input ray path file
! ofiledt = output depth slice file for traveltimes
! ofilenst = output N-S slice file for traveltimes
! ofileewt = output E-W slice file for traveltimes
! ofiledb = bounds for output depth slice file
! ofilensb = bounds for output N-S slice file
! ofileewb = bounds for output E-W slice file
! ofiledv = output velocity file for depth slice
! ofilensv = output velocity file for N-S slice
! ofileewv = output velocity file for E-W slice
! ofiledr = output ray path file in depth
! ofilensr = output ray path file in N-S
! ofileewr = output ray path file in E-W
! nnr,nnt,nnp = number of diced nodes in r,theta,phi
! gor = grid origin in radius
! got = grid origin in theta (N-S)
! gop = grid origin in phi (E-W)
! rgsr,rgst,rgsp = Refined node spacing in r,theta,phi
! idp = radial index of slice depth
! ins = NS index for latitude slice
! iew = EW index for longitude slice
! tt = traveltime at or near a node
! ttn = traveltime grid values
! ttt = traveltime of node above slice 
! ttb =  traveltime of node below slice
! lft,rgt,btm,top = plotting bounds
! pvgd = plot velocity depth slice (0=no,1=yes)
! pvgns = plot velocity N-S slice (0=no, 1=yes)
! pvgew = plot velocity E-W slice (0=no, 1=yes)
! ptfd = plot traveltime field depth slice? (0=no,1=yes)
! ptfns = plot traveltime field N-S slice? (0=no, 1=yes)
! ptfew = plot traveltime field E-W slice? (0=no,1=yes)
! veln = velocity grid values
! prpd = plot ray path in depth? (0=no,1=yes)
! prpns = plot ray path in N-S? (0=no,1=yes)
! prpew = plot ray path in E-W? (0=no,1=yes)
! nre = number of ray elements
! rdep,rlat,rlon = ray point depth, latitude, longitude
! nr = number of receivers
! sep = character marker for separating rays
! extrds = extract depth slice? (0=no,1=yes)
! extrnss = extract N-S slice? (0=no,1=yes)
! extrews = extract E-W slice? (0=no,1=yes)
! ddt,ddp = dicing of velocity grid for depth slice
! dnsr,dnst = dicing of velocity grid for N-S slice
! dewr,dewp = dicing of velocity grid for E-W slice
! u = bspline independent coordinate
! ui,vi = bspline basis functions
! vela = diced velocity values
! nnx,nnz = dimensions of vela
! conr,conp,cont = variables for edge of bspline grid
! str,stp,stt = counters for vela grid points
! rnode = reference node for slice
!
OPEN(UNIT=10,FILE='gmtslicet.in',STATUS='old')

READ(10,*)
READ(10,*)
READ(10,*)
READ(10,'(a26)')ifilev
READ(10,'(a26)')ifilevr
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)extrds
READ(10,*)sldep
READ(10,'(a22)')ofiledb
READ(10,*)extrnss
READ(10,*)slns
READ(10,'(a22)')ofilensb
READ(10,*)extrews
READ(10,*)slew
READ(10,'(a22)')ofileewb

OPEN(UNIT=20,FILE=ifilev,status='old')
READ(20,*)nnr,nnt,nnp
READ(20,*)gor,got,gop
READ(20,*)rgsr,rgst,rgsp
READ(20,*)rnnr,rnnt,rnnp
READ(20,*)rgor,rgot,rgop
READ(20,*)rrgsr,rrgst,rrgsp

ALLOCATE(veln(0:nnr+1,0:nnt+1,0:nnp+1), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM tslice: REAL veln'
ENDIF
DO i=0,nnp+1
   DO j=0,nnt+1
      DO k=0,nnr+1
         READ(20,*)veln(k,j,i)
      ENDDO
   ENDDO
ENDDO
CLOSE(20)

READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)pvgd
READ(10,*)ddt,ddp
READ(10,'(a22)')ofiledv
READ(10,*)pvgns
READ(10,*)dnsr,dnst  
READ(10,'(a22)')ofilensv 
READ(10,*)pvgew
READ(10,*)dewr,dewp   
READ(10,'(a22)')ofileewv 
CLOSE(10)


IF(extrds.EQ.1)THEN
   lft=rgop
   rgt=rgop+(rnnp-1)*rrgsp
   btm=rgot-(rnnt-1)*rrgst
   top=rgot
   OPEN(UNIT=50,FILE=ofiledb,STATUS='unknown')
   WRITE(50,'(f16.11)')lft
   WRITE(50,'(f16.11)')rgt
   WRITE(50,'(f16.11)')btm
   WRITE(50,'(f16.11)')top
   id1=(rnnp-1)*ddp+1
   id2=(rnnt-1)*ddt+1
   WRITE(50,*)id1
   WRITE(50,*)id2
ENDIF


IF(extrnss.EQ.1)THEN
   rgt=rgot-(rnnt-1)*rrgst
   lft=rgot
   btm=rgor-(rnnr-1)*rrgsr
   top=rgor
   OPEN(UNIT=60,FILE=ofilensb,STATUS='unknown')
   WRITE(60,'(f16.11)')rgt
   WRITE(60,'(f16.11)')lft
   WRITE(60,'(f16.11)')btm
   WRITE(60,'(f16.11)')top
   id1=(rnnt-1)*dnst+1
   id2=(rnnr-3)*dnsr+1
   WRITE(60,*)id1
   WRITE(60,*)id2
ENDIF

IF(extrews.EQ.1)THEN
   lft=rgop
   rgt=rgop+(rnnp-1)*rrgsp
   btm=rgor-(rnnr-1)*rrgsr
   top=rgor
   OPEN(UNIT=70,FILE=ofileewb,STATUS='unknown')
   WRITE(70,'(f16.11)')lft
   WRITE(70,'(f16.11)')rgt
   WRITE(70,'(f16.11)')btm
   WRITE(70,'(f16.11)')top
   id1=(rnnp-1)*dewp+1
   id2=(rnnr-3)*dewr+1
   WRITE(70,*)id1
   WRITE(70,*)id2
ENDIF

IF(pvgd.EQ.1.OR.pvgns.EQ.1.OR.pvgew.EQ.1)THEN
   isw=0
   OPEN(UNIT=20,FILE=ifilevr,status='old')
   READ(20,*)idm1,idm2,idm3
   IF(idm1.NE.nnr.OR.idm2.NE.nnt.OR.idm3.NE.nnp)isw=1
   READ(20,*)rd1,rd2,rd3
   IF(rd1.NE.gor.OR.rd2.NE.got.OR.rd3.NE.gop)isw=1
   READ(20,*)rd1,rd2,rd3
   IF(rd1.NE.rgsr.OR.rd2.NE.rgst.OR.rd3.NE.rgsp)isw=1
   
   IF(isw.EQ.1)THEN
      WRITE(6,*)'ERROR! Actual velocity grid and reference'
      WRITE(6,*)'velocity grid have different dimensions or'
      WRITE(6,*)'different numbers of grid points!'
      WRITE(6,*)'TERMINATING PROGRAM!!!'
   ENDIF
   READ(20,*)
   READ(20,*)
   READ(20,*)
   
   DO i=0,nnp+1
      DO j=0,nnt+1
         DO k=0,nnr+1
            READ(20,*)rd1
            veln(k,j,i)=veln(k,j,i)-rd1
         ENDDO
      ENDDO
   ENDDO
   CLOSE(20)
ENDIF

ALLOCATE(velv(rnnr,0:rnnt+1,0:rnnp+1))
dvx=rgst*pi/180.0
dvz=rgsp*pi/180.0
gox=(90.0-got)*pi/180.0
goz=gop*pi/180.0

dnx=rrgst*pi/180.0
dnz=rrgsp*pi/180.0
rgox=(90.0-rgot)*pi/180.0
rgoz=rgop*pi/180.0

DO k=1,rnnr
   ri=INT((gor-rgor+(k-1)*rrgsr)/rgsr)+1
   distr=-(rgor-(k-1)*rrgsr)+(gor-(ri-1)*rgsr)
   u=distr/rgsr
   ai(1)=(1.0-u)**3/6.0
   ai(2)=(4.0-6.0*u**2+3.0*u**3)/6.0
   ai(3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
   ai(4)=u**3/6.0
   DO j=0,rnnt+1
      xi=INT((rgox+(j-1)*dnx-gox)/dvx)+1
      distx=rgox+(j-1)*dnx-gox-(xi-1)*dvx
      v=distx/dvx
      bi(1)=(1.0-v)**3/6.0
      bi(2)=(4.0-6.0*v**2+3.0*v**3)/6.0
      bi(3)=(1.0+3.0*v+3.0*v**2-3.0*v**3)/6.0
      bi(4)=v**3/6.0
      DO i=0,rnnp+1
         zi=INT((rgoz+(i-1)*dnz-goz)/dvz)+1
         distz=rgoz+(i-1)*dnz-goz-(zi-1)*dvz
         w=distz/dvz
         ci(1)=(1.0-w)**3/6.0
         ci(2)=(4.0-6.0*w**2+3.0*w**3)/6.0
         ci(3)=(1.0+3.0*w+3.0*w**2-3.0*w**3)/6.0
         ci(4)=w**3/6.0
         sumi=0.0
         DO i1=1,4 !z
            sumj=0.0 
            DO j1=1,4 !x
               sumk=0.0
               DO k1=1,4 !r
                  sumk=sumk+ai(k1)*veln(ri-2+k1,xi-2+j1,zi-2+i1)
               ENDDO
               sumj=sumj+bi(j1)*sumk
            ENDDO
            sumi=sumi+ci(i1)*sumj
         ENDDO
         velv(k,j,i)=sumi
      ENDDO
   ENDDO
ENDDO

IF(pvgd.EQ.1)THEN

  nnx=(rnnt-1)*ddt+1
  nnz=(rnnp-1)*ddp+1
  ALLOCATE(vela(nnx,nnz), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL vela'
  ENDIF

  ALLOCATE(ui(ddt+1,4), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL ui'
  ENDIF
  DO i=1,ddt+1
     u=ddt
     u=(i-1)/u
     ui(i,1)=(1.0-u)**3/6.0
     ui(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
     ui(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
     ui(i,4)=u**3/6.0
  ENDDO
  ALLOCATE(vi(ddp+1,4), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL vi'
  ENDIF
  DO i=1,ddp+1
     u=ddp
     u=(i-1)/u
     vi(i,1)=(1.0-u)**3/6.0
     vi(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
     vi(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
     vi(i,4)=u**3/6.0
  ENDDO

  rnode=INT((rgor-sldep)/rrgsr)+1
  u=ABS(rgor-sldep-(rnode-1)*rrgsr)/rrgsr
  wi(1)=(1.0-u)**3/6.0
  wi(2)=(4.0-6.0*u**2+3.0*u**3)/6.0
  wi(3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
  wi(4)=u**3/6.0
  DO i=1,rnnp-1
     conp=ddp
     IF(i==rnnp-1)conp=ddp+1
     DO j=1,rnnt-1
        cont=ddt
        IF(j==rnnt-1)cont=ddt+1
        DO l=1,conp
           stp=ddp*(i-1)+l
           DO m=1,cont
              stt=ddt*(j-1)+m
              sumi=0.0
              DO i1=1,4
                 sumj=0.0
                 DO j1=1,4
                    sumk=0.0
                    DO k1=1,4
                       rdm=wi(k1)*velv(rnode-2+k1,j-2+j1,i-2+i1)
                       sumk=sumk+rdm
                    ENDDO
                    sumj=sumj+ui(m,j1)*sumk
                 ENDDO
                 sumi=sumi+vi(l,i1)*sumj
              ENDDO
              vela(stt,stp)=sumi
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  OPEN(UNIT=30,FILE=ofiledv,STATUS='unknown')
  DO i=1,nnz
     DO j=nnx,1,-1
         WRITE(30,*)vela(j,i)
     ENDDO
  ENDDO
  CLOSE(30)
  DEALLOCATE(vela,ui,vi, STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with DEALLOCATE: PROGRAM slice: REAL vela,ui,vi'
  ENDIF
ENDIF

IF(pvgns.EQ.1)THEN
  nnx=(rnnr-3)*dnsr+1
  nnz=(rnnt-1)*dnst+1
  ALLOCATE(vela(nnx,nnz), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL vela'
  ENDIF

  ALLOCATE(ui(dnsr+1,4), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL ui'
  ENDIF
  DO i=1,dnsr+1
     u=dnsr
     u=(i-1)/u
     ui(i,1)=(1.0-u)**3/6.0
     ui(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
     ui(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
     ui(i,4)=u**3/6.0
  ENDDO
  ALLOCATE(vi(dnst+1,4), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL vi'
  ENDIF
  DO i=1,dnst+1
     u=dnst
     u=(i-1)/u
     vi(i,1)=(1.0-u)**3/6.0
     vi(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
     vi(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
     vi(i,4)=u**3/6.0
  ENDDO

  rnode=INT((slns-rgop)/rrgsp)+1
  u=ABS(slns-rgop-(rnode-1)*rrgsp)/rrgsp
  wi(1)=(1.0-u)**3/6.0
  wi(2)=(4.0-6.0*u**2+3.0*u**3)/6.0
  wi(3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
  wi(4)=u**3/6.0
  DO j=1,rnnt-1
     cont=dnst
     IF(j==rnnt-1)cont=dnst+1
     DO k=2,rnnr-2
        conr=dnsr
        IF(k==rnnr-2)conr=dnsr+1
        DO m=1,cont
           stt=dnst*(j-1)+m
           DO n=1,conr
              str=dnsr*(k-2)+n
              sumi=0.0
              DO i1=1,4
                 sumj=0.0
                 DO j1=1,4
                    sumk=0.0
                    DO k1=1,4
                       rdm=ui(n,k1)*velv(k-2+k1,j-2+j1,rnode-2+i1)
                       sumk=sumk+rdm
                    ENDDO
                    sumj=sumj+vi(m,j1)*sumk
                 ENDDO
                 sumi=sumi+wi(i1)*sumj
              ENDDO
              vela(str,stt)=sumi
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  OPEN(UNIT=30,FILE=ofilensv,STATUS='unknown')
  DO i=1,nnz
     DO j=nnx,1,-1
         WRITE(30,*)vela(j,i)
     ENDDO
  ENDDO
  CLOSE(30)
  DEALLOCATE(vela,ui,vi, STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with DEALLOCATE: PROGRAM slice: REAL vela,ui,vi'
  ENDIF
ENDIF

IF(pvgew.EQ.1)THEN

  nnx=(rnnr-3)*dewr+1
  nnz=(rnnp-1)*dewp+1
  ALLOCATE(vela(nnx,nnz), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL vela'
  ENDIF

  ALLOCATE(ui(dewr+1,4), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL ui'
  ENDIF
  DO i=1,dewr+1
     u=dewr
     u=(i-1)/u
     ui(i,1)=(1.0-u)**3/6.0
     ui(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
     ui(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
     ui(i,4)=u**3/6.0
  ENDDO
  ALLOCATE(vi(dewp+1,4), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL vi'
  ENDIF
  DO i=1,dewp+1
     u=dewp
     u=(i-1)/u
     vi(i,1)=(1.0-u)**3/6.0
     vi(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
     vi(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
     vi(i,4)=u**3/6.0
  ENDDO

  rnode=INT((rgot-slew)/rrgst)+1
  u=ABS(rgot-slew-(rnode-1)*rrgst)/rrgst
  wi(1)=(1.0-u)**3/6.0
  wi(2)=(4.0-6.0*u**2+3.0*u**3)/6.0
  wi(3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
  wi(4)=u**3/6.0
  DO i=1,rnnp-1
     conp=dewp
     IF(i==rnnp-1)conp=dewp+1
     DO k=2,rnnr-2
        conr=dewr
        IF(k==rnnr-2)conr=dewr+1
        DO l=1,conp
           stp=dewp*(i-1)+l
           DO n=1,conr
              str=dewr*(k-2)+n
              sumi=0.0
              DO i1=1,4
                 sumj=0.0
                 DO j1=1,4
                    sumk=0.0
                    DO k1=1,4
                       rdm=ui(n,k1)*velv(k-2+k1,rnode-2+j1,i-2+i1)
                       sumk=sumk+rdm
                    ENDDO
                    sumj=sumj+wi(j1)*sumk
                 ENDDO
                 sumi=sumi+vi(l,i1)*sumj
              ENDDO
              vela(str,stp)=sumi
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  OPEN(UNIT=30,FILE=ofileewv,STATUS='unknown')
  DO i=1,nnz
     DO j=nnx,1,-1
         WRITE(30,*)vela(j,i)
     ENDDO
  ENDDO
  CLOSE(30)
  DEALLOCATE(vela,ui,vi, STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with DEALLOCATE: PROGRAM slice: REAL vela,ui,vi'
  ENDIF
ENDIF
DEALLOCATE(veln, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: PROGRAM slice: REAL veln'
ENDIF
DEALLOCATE(velv, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: PROGRAM slice: REAL veln'
ENDIF

WRITE(50,*)nnp
WRITE(50,*)nnt
WRITE(60,*)nnt
WRITE(60,*)nnr
WRITE(70,*)nnp
WRITE(70,*)nnr
CLOSE(50)
CLOSE(60)
CLOSE(70)

STOP
END PROGRAM slice
