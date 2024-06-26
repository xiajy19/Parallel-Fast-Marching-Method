!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NAME:     PARALLEL FAST MARCHING METHOD
! CODE:     FORTRAN 90, MPI
! AUTHOR:   JUNYI XIA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE globalp
   IMPLICIT NONE
   INTEGER, PARAMETER :: i5=SELECTED_REAL_KIND(5,10)
   INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
   INTEGER :: checkstat
   INTEGER, SAVE :: s1,s2
   INTEGER, SAVE :: nnr,nnx,nnz,nrc,fom,ltfr
   INTEGER, SAVE :: nvr,nvx,nvz
   INTEGER, SAVE :: nsnn,nsns,nsne,nsnw
   INTEGER, SAVE :: nprocs,myrank
   INTEGER, SAVE :: nr,nx,nz
   INTEGER, SAVE :: rank_r,rank_x,rank_z
   INTEGER, SAVE :: left,right,in,out,up,down
   INTEGER, SAVE :: rloca,xloca,zloca
   INTEGER, SAVE :: roffset,xoffset,zoffset
   INTEGER, SAVE :: xnum,znum,rnum
   INTEGER, SAVE :: comm_cart,comm_start
   REAL(KIND=i10), SAVE :: mini,maxi
   INTEGER, SAVE :: ntr,outsize 
   INTEGER, SAVE :: sw_fmm,cal_already
   REAL(KIND=i10), SAVE :: error_to_stop,pre_error_back_propagation,cur_error_back_propagation
   INTEGER, DIMENSION (6), SAVE :: direction_send,send_direction
   INTEGER, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: nsts
   REAL(KIND=i10), DIMENSION (:,:,:), ALLOCATABLE, SAVE :: ttn,ttn_previous_step
   REAL(KIND=i10), DIMENSION (:,:), ALLOCATABLE, SAVE :: ttn1,ttn2
   CHARACTER, DIMENSION (:), ALLOCATABLE :: BUF,LBUF1,LBUF2,RBUF1,RBUF2
   INTEGER, SAVE :: type_left,type_right,type_up,type_down
   INTEGER, SAVE :: type_recv1, type_recv2
   INTEGER, SAVE :: type_in_up_down, type_out_up_down
   INTEGER, SAVE :: type_left_in_out, type_right_in_out
   REAL(KIND=i10), DIMENSION (:,:,:), ALLOCATABLE, SAVE :: ttn_ghost_left,ttn_ghost_right
   REAL(KIND=i10), DIMENSION (:,:,:), ALLOCATABLE, SAVE :: ttn_ghost_in,ttn_ghost_out
   REAL(KIND=i10), DIMENSION (:,:,:), ALLOCATABLE, SAVE :: ttn_ghost_up,ttn_ghost_down
   INTEGER, SAVE :: srtimes_myrank,srtimes_nprocs,comm_srtimes,comm_srtimes_to_file
   INTEGER, SAVE :: number_of_trr
   INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE :: rec_location
   INTEGER, DIMENSION (:), ALLOCATABLE, SAVE :: number_file,temp_number_file,number
   INTEGER, DIMENSION (:), ALLOCATABLE, SAVE :: displs_srtimes,rcounts_srtimes
   REAL, DIMENSION (:), ALLOCATABLE, SAVE :: trr_file
   REAL(KIND=i5), DIMENSION (:), ALLOCATABLE, SAVE :: trr
   INTEGER, SAVE :: rpaths_stop_criterion,communicate_rpaths_length
   INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: raypath_rec_location
   REAL(KIND=i10), DIMENSION (:,:), ALLOCATABLE ,SAVE :: total_communicate_rpaths,communicate_rpaths,communicate_only
   INTEGER, DIMENSION (:), ALLOCATABLE, SAVE :: rcounts_rpaths, displs_rpaths
   INTEGER, SAVE :: temp_fdm,rcounts_process
   INTEGER, DIMENSION (:), ALLOCATABLE, SAVE :: displs_fdm,rcounts_fdm
   INTEGER, DIMENSION (:), ALLOCATABLE, SAVE :: coln_process
   REAL(KIND=i5), DIMENSION (:), ALLOCATABLE, SAVE :: fdm_process
   REAL(KIND=i10), SAVE :: gor,gox,goz,dnr,dnx,dnz,earth,rgor,rgox,rgoz
   REAL(KIND=i10), SAVE :: dvr,dvx,dvz
   REAL(KIND=i10), SAVE :: goxd,gozd,rgoxd,rgozd,dnxd,dnzd
   REAL(KIND=i10), DIMENSION (:,:,:), ALLOCATABLE, SAVE :: veln,velv
   REAL(KIND=i10), DIMENSION (:), ALLOCATABLE, SAVE :: rcr,rcx,rcz
   REAL(KIND=i10), PARAMETER :: pi=3.1415926535898
END MODULE globalp


MODULE traveltime
   USE globalp
   IMPLICIT NONE
   TYPE backpointer
      INTEGER(KIND=2) :: pr,px,pz
   END TYPE backpointer
   TYPE temppointer
      INTEGER(KIND=2) :: pa,pb,pwhere
      REAL(KIND=i10) :: pttn
   END TYPE temppointer
   TYPE(backpointer), DIMENSION (:), ALLOCATABLE :: btg
   CONTAINS

   SUBROUTINE local_fmm
      IMPLICIT NONE
      INTEGER :: i,j,ir,ix,iz,sw,ittmin,jttmin
      INTEGER :: isum,tnrn
      REAL(KIND=i10) :: old_value
      REAL(KIND=i10) :: rd1,ttmin
      IF(ltfr.EQ.1)THEN
         tnrn=(nsns+1-nsnn)*(nsne+1-nsnw)
         isum=0
      ENDIF
      DO WHILE(ntr.gt.0)
         IF(cal_already.EQ.0)THEN
            cal_already=1
         ENDIF
         ir=btg(1)%pr
         ix=btg(1)%px
         iz=btg(1)%pz
         nsts(iz,ix,ir)=0
         CALL downtree_all
         DO i=ir-1,ir+1,2
            IF(i.GE.1.and.i.LE.rnum)THEN
               IF(nsts(iz,ix,i).EQ.-1)THEN
                  IF(fom.EQ.0)THEN
                     CALL fouds1_all(iz,ix,i)
                  ELSE
                     CALL fouds2_all(iz,ix,i)
                  ENDIF
                  CALL addtree_all(iz,ix,i)
               ELSE IF(nsts(iz,ix,i).GT.0)THEN
                  IF(fom.EQ.0)THEN
                     CALL fouds1_all(iz,ix,i)
                  ELSE
                     CALL fouds2_all(iz,ix,i)
                  ENDIF
                  CALL updtree_all(iz,ix,i)
               ENDIF
            ENDIF
         ENDDO
         DO i=ix-1,ix+1,2
            IF(i.GE.1.and.i.LE.xnum)THEN
               IF(nsts(iz,i,ir).EQ.-1)THEN
                  IF(fom.EQ.0)THEN
                     CALL fouds1_all(iz,i,ir)
                  ELSE
                     CALL fouds2_all(iz,i,ir)
                  ENDIF
                  CALL addtree_all(iz,i,ir)
               ELSE IF(nsts(iz,i,ir).GT.0)THEN
                  IF(fom.EQ.0)THEN
                     CALL fouds1_all(iz,i,ir)
                  ELSE
                     CALL fouds2_all(iz,i,ir)
                  ENDIF
                  CALL updtree_all(iz,i,ir)
               ENDIF
            ENDIF
         ENDDO
         DO i=iz-1,iz+1,2
            IF(i.GE.1.and.i.LE.znum)THEN
               IF(nsts(i,ix,ir).EQ.-1)THEN
                  IF(fom.EQ.0)THEN
                     CALL fouds1_all(i,ix,ir)
                  ELSE
                     CALL fouds2_all(i,ix,ir)
                  ENDIF
                  CALL addtree_all(i,ix,ir)
               ELSE IF(nsts(i,ix,ir).GT.0)THEN
                  IF(fom.EQ.0)THEN
                     CALL fouds1_all(i,ix,ir)
                  ELSE
                     CALL fouds2_all(i,ix,ir)
                  ENDIF
                  CALL updtree_all(i,ix,ir)
               ENDIF
            ENDIF
         ENDDO
      ENDDO
   END SUBROUTINE local_fmm

   SUBROUTINE order_overlap
      IMPLICIT NONE
      INTEGER :: i,j,k,temp,sw
      INTEGER :: ix,iz,ir,tpp,tpc,maxbt
      REAL(KIND=i10) :: old_value,ttn_temp,temp1,temp2
      TYPE(temppointer) :: exchh
      TYPE(temppointer), DIMENSION (:), ALLOCATABLE :: btgg
      maxbt=nint(0.5*znum*xnum*rnum)
      IF(send_direction(1).EQ.1)THEN
         ALLOCATE(btgg(maxbt))
         temp=0
         DO i=1,rnum
            DO j=1,znum
               IF(ttn_ghost_left(j,2,i).LT.ttn(j,1,i).OR.ttn(j,1,i).EQ.0)THEN
                  temp=temp+1
                  btgg(temp)%pa=j
                  btgg(temp)%pb=i
                  btgg(temp)%pwhere=1
                  btgg(temp)%pttn=ttn_ghost_left(j,2,i)
                  tpc=temp
                  tpp=tpc/2
                  DO WHILE(tpp.gt.0)
                     IF(ttn_ghost_left(j,2,i).LT.btgg(tpp)%pttn)THEN
                        exchh=btgg(tpc)
                        btgg(tpc)=btgg(tpp)
                        btgg(tpp)=exchh
                        tpc=tpp
                        tpp=tpc/2
                     ELSE
                        tpp=0
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO       
         DO WHILE(temp.GT.0)
            iz=btgg(1)%pa
            ix=1
            ir=btgg(1)%pb
            old_value = ttn(iz,ix,ir)
            CALL fouds2_left(iz,ix,ir,old_value)
            CALL downtree_left(btgg,temp)
         ENDDO
         DEALLOCATE(btgg)
      ENDIF
      IF(send_direction(2).EQ.1)THEN
         ALLOCATE(btgg(maxbt))
         temp=0 
         DO i=1,rnum
            DO j=1,znum
               IF(ttn_ghost_right(j,1,i).LT.ttn(j,xnum,i).OR.ttn(j,xnum,i).EQ.0)THEN
                  temp=temp+1
                  btgg(temp)%pa=j
                  btgg(temp)%pb=i
                  btgg(temp)%pwhere=2
                  btgg(temp)%pttn=ttn_ghost_right(j,1,i)
                  tpc=temp
                  tpp=tpc/2
                  DO WHILE(tpp.gt.0)
                     IF(ttn_ghost_right(j,1,i).LT.btgg(tpp)%pttn)THEN
                        exchh=btgg(tpc)
                        btgg(tpc)=btgg(tpp)
                        btgg(tpp)=exchh
                        tpc=tpp
                        tpp=tpc/2
                     ELSE
                        tpp=0
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
         DO WHILE(temp.GT.0)
            iz=btgg(1)%pa
            ix=xnum
            ir=btgg(1)%pb
            old_value = ttn(iz,ix,ir)
            CALL fouds2_right(iz,ix,ir,old_value)
            CALL downtree_right(btgg,temp)
         ENDDO
         DEALLOCATE(btgg)
      ENDIF
      IF(send_direction(3).EQ.1)THEN         
         ALLOCATE(btgg(maxbt))
         temp=0 
         DO i=1,rnum
            DO j=1,xnum
               IF(ttn_ghost_in(2,j,i).LT.ttn(1,j,i).OR.ttn(1,j,i).EQ.0)THEN
                  temp=temp+1
                  btgg(temp)%pa=j
                  btgg(temp)%pb=i
                  btgg(temp)%pwhere=3
                  btgg(temp)%pttn=ttn_ghost_in(2,j,i)
                  tpc=temp
                  tpp=tpc/2
                  DO WHILE(tpp.GT.0)
                     IF(ttn_ghost_in(2,j,i).LT.btgg(tpp)%pttn)THEN
                        exchh=btgg(tpc)
                        btgg(tpc)=btgg(tpp)
                        btgg(tpp)=exchh
                        tpc=tpp
                        tpp=tpc/2
                     ELSE
                        tpp=0
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
         DO WHILE(temp.GT.0)
            iz=1            
            ix=btgg(1)%pa   
            ir=btgg(1)%pb   
            old_value = ttn(iz,ix,ir)
            CALL fouds2_in(iz,ix,ir,old_value)
            CALL downtree_in(btgg,temp)
         ENDDO
         DEALLOCATE(btgg)
      ENDIF
      IF(send_direction(4).EQ.1)THEN 
         ALLOCATE(btgg(maxbt))
         temp=0
         DO i=1,rnum
            DO j=1,xnum
               IF(ttn_ghost_out(1,j,i).LT.ttn(znum,j,i).OR.ttn(znum,j,i).EQ.0)THEN
                  temp=temp+1
                  btgg(temp)%pa=j
                  btgg(temp)%pb=i
                  btgg(temp)%pwhere=4
                  btgg(temp)%pttn=ttn_ghost_out(1,j,i)
                  tpc=temp
                  tpp=tpc/2
                  DO WHILE(tpp.gt.0)
                     IF(ttn_ghost_out(1,j,i).LT.btgg(tpp)%pttn)THEN
                        exchh=btgg(tpc)
                        btgg(tpc)=btgg(tpp)
                        btgg(tpp)=exchh
                        tpc=tpp
                        tpp=tpc/2
                     ELSE
                        tpp=0
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
         DO WHILE(temp.GT.0)
            iz=znum         
            ix=btgg(1)%pa   
            ir=btgg(1)%pb   
            old_value = ttn(iz,ix,ir)
            CALL fouds2_out(iz,ix,ir,old_value)
            CALL downtree_out(btgg,temp)
         ENDDO
         DEALLOCATE(btgg)
      ENDIF
      IF(send_direction(5).EQ.1)THEN
         ALLOCATE(btgg(maxbt))
         temp=0
         DO i=1,xnum
            DO j=1,znum
               IF(ttn_ghost_up(j,i,2).LT.ttn(j,i,1).OR.ttn(j,i,1).EQ.0)THEN
                  temp=temp+1
                  btgg(temp)%pa=j
                  btgg(temp)%pb=i
                  btgg(temp)%pwhere=5
                  btgg(temp)%pttn=ttn_ghost_up(j,i,2)
                  tpc=temp
                  tpp=tpc/2
                  DO WHILE(tpp.gt.0)
                     IF(ttn_ghost_up(j,i,2).LT.btgg(tpp)%pttn)THEN
                        exchh=btgg(tpc)
                        btgg(tpc)=btgg(tpp)
                        btgg(tpp)=exchh
                        tpc=tpp
                        tpp=tpc/2
                     ELSE
                        tpp=0
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
         DO WHILE(temp.GT.0)
            iz=btgg(1)%pa   
            ix=btgg(1)%pb   
            ir=1            
            old_value = ttn(iz,ix,ir)
            CALL fouds2_up(iz,ix,ir,old_value)
            CALL downtree_up(btgg,temp)
         ENDDO
         DEALLOCATE(btgg)
      ENDIF
      IF(send_direction(6).EQ.1)THEN
         ALLOCATE(btgg(maxbt))
         temp=0
         DO i=1,xnum
            DO j=1,znum
               IF(ttn_ghost_down(j,i,1).LT.ttn(j,i,rnum).OR.ttn(j,i,rnum).EQ.0)THEN
                  temp=temp+1
                  btgg(temp)%pa=j
                  btgg(temp)%pb=i
                  btgg(temp)%pwhere=6
                  btgg(temp)%pttn=ttn_ghost_down(j,i,1)
                  tpc=temp
                  tpp=tpc/2
                  DO WHILE(tpp.gt.0)           
                     IF(ttn_ghost_down(j,i,1).LT.btgg(tpp)%pttn)THEN
                        exchh=btgg(tpc)
                        btgg(tpc)=btgg(tpp)
                        btgg(tpp)=exchh
                        tpc=tpp
                        tpp=tpc/2
                     ELSE
                        tpp=0
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO  
         DO WHILE(temp.GT.0)
            iz=btgg(1)%pa   
            ix=btgg(1)%pb   
            ir=rnum         
            old_value = ttn(iz,ix,ir)
            CALL fouds2_down(iz,ix,ir,old_value)
            CALL downtree_down(btgg,temp)
         ENDDO
         DEALLOCATE(btgg)
      ENDIF
      IF(ntr.GE.1)THEN
         sw_fmm=1
         mini=ttn(btg(1)%pz,btg(1)%px,btg(1)%pr)
         maxi=0.0
         DO i=1,ntr
            IF(ttn(btg(i)%pz,btg(i)%px,btg(i)%pr).GT.maxi)THEN   
               maxi=ttn(btg(i)%pz,btg(i)%px,btg(i)%pr)
            ENDIF 
         ENDDO 
         DO i=1,rnum
            DO j=1,xnum
               DO k=1,znum
                  IF(nsts(k,j,i).LE.0)THEN
                     IF(ttn(k,j,i).GT.0)THEN
                        IF(ttn(k,j,i).LE.mini)THEN
                           nsts(k,j,i)=0
                        ELSEIF(ttn(k,j,i).GE.maxi)THEN
                           nsts(k,j,i)=-1             
                        ELSE
                           CALL addtree_all(k,j,i)
                        ENDIF
                     ELSE
                        nsts(k,j,i)=-1
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDIF
   END SUBROUTINE order_overlap

   SUBROUTINE direction_mixed_case(iz,ix,ir)
      IMPLICIT NONE
      INTEGER :: iz,ix,ir
      IF(rank_x.NE.0.AND.ix==1)THEN
         IF(send_direction(1).EQ.1.AND.ttn_ghost_left(iz,2,ir).LT.ttn(iz,ix,ir))THEN
            IF(rank_z.NE.0.AND.iz==1)THEN
               IF(send_direction(3).EQ.1.AND.ttn_ghost_in(2,ix,ir).LT.ttn(iz,ix,ir))THEN
                  CALL fouds2_left_in(iz,ix,ir)
               ENDIF
            ELSEIF(rank_z.NE.nz-1.AND.iz==znum)THEN
               IF(send_direction(4).EQ.1.AND.ttn_ghost_out(1,ix,ir).LT.ttn(iz,ix,ir))THEN
                  CALL fouds2_left_out(iz,ix,ir)
               ENDIF
            ENDIF
            IF(rank_r.NE.0.AND.ir==1)THEN
               IF(send_direction(5).EQ.1.AND.ttn_ghost_up(iz,ix,2).LT.ttn(iz,ix,ir))THEN
                  CALL fouds2_left_up(iz,ix,ir)
               ENDIF
            ELSEIF(rank_r.NE.nr-1.AND.ir==rnum)THEN
               IF(send_direction(6).EQ.1.AND.ttn_ghost_down(iz,ix,1).LT.ttn(iz,ix,ir))THEN
                  CALL fouds2_left_down(iz,ix,ir)
               ENDIF
            ENDIF
         ENDIF
      ELSEIF(rank_x.NE.nx-1.AND.ix==xnum)THEN
         IF(send_direction(2).EQ.1.AND.ttn_ghost_right(iz,1,ir).LT.ttn(iz,ix,ir))THEN
            IF(rank_z.NE.0.AND.iz==1)THEN
               IF(send_direction(3).EQ.1.AND.ttn_ghost_in(2,ix,ir).LT.ttn(iz,ix,ir))THEN
                  CALL fouds2_right_in(iz,ix,ir)
               ENDIF
            ELSEIF(rank_z.NE.nz-1.AND.iz==znum)THEN
               IF(send_direction(4).EQ.1.AND.ttn_ghost_out(1,ix,ir).LT.ttn(iz,ix,ir))THEN
                  CALL fouds2_right_out(iz,ix,ir)
               ENDIF
            ENDIF
            IF(rank_r.NE.0.AND.ir==1)THEN
               IF(send_direction(5).EQ.1.AND.ttn_ghost_up(iz,ix,2).LT.ttn(iz,ix,ir))THEN
                  CALL fouds2_right_up(iz,ix,ir)
               ENDIF
            ELSEIF(rank_r.NE.nr-1.AND.ir==rnum)THEN
               IF(send_direction(6).EQ.1.AND.ttn_ghost_down(iz,ix,1).LT.ttn(iz,ix,ir))THEN
                  CALL fouds2_right_down(iz,ix,ir)
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      IF(rank_z.NE.0.AND.iz==1)THEN
         IF(send_direction(3).EQ.1.AND.ttn_ghost_in(2,ix,ir).LT.ttn(iz,ix,ir))THEN
            IF(rank_r.NE.0.AND.ir==1)THEN
               IF(send_direction(5).EQ.1.AND.ttn_ghost_up(iz,ix,2).LT.ttn(iz,ix,ir))THEN
                  CALL fouds2_in_up(iz,ix,ir)
               ENDIF
            ELSEIF(rank_r.NE.nr-1.AND.ir==rnum)THEN
               IF(send_direction(6).EQ.1.AND.ttn_ghost_down(iz,ix,1).LT.ttn(iz,ix,ir))THEN
                  CALL fouds2_in_down(iz,ix,ir)
               ENDIF
            ENDIF
         ENDIF
      ELSEIF(rank_z.NE.nz-1.AND.iz==znum)THEN
         IF(send_direction(4).EQ.1.AND.ttn_ghost_out(1,ix,ir).LT.ttn(iz,ix,ir))THEN
            IF(rank_r.NE.0.AND.ir==1)THEN
               IF(send_direction(5).EQ.1.AND.ttn_ghost_up(iz,ix,2).LT.ttn(iz,ix,ir))THEN
                  CALL fouds2_out_up(iz,ix,ir)
               ENDIF
            ELSEIF(rank_r.NE.nr-1.AND.ir==rnum)THEN
               IF(send_direction(6).EQ.1.AND.ttn_ghost_down(iz,ix,1).LT.ttn(iz,ix,ir))THEN
                  CALL fouds2_out_down(iz,ix,ir)
               ENDIF
            ENDIF
         ENDIF
      ENDIF
   END SUBROUTINE direction_mixed_case

   SUBROUTINE downtree_ghost_point(btgg,temp)
      IMPLICIT NONE
      INTEGER :: tpp,tpc
      INTEGER :: temp
      REAL(KIND=i10) :: rd1,rd2
      TYPE(temppointer) :: exchh
      TYPE(temppointer), DIMENSION (:) :: btgg
      IF(temp.EQ.1)THEN
         temp=temp-1
         RETURN
      ENDIF
      btgg(1)=btgg(temp)
      temp=temp-1
      tpp=1
      tpc=2*tpp
      DO WHILE(tpc.lt.temp)
         rd1=btgg(tpc)%pttn
         rd2=btgg(tpc+1)%pttn
         IF(rd1.gt.rd2)THEN
            tpc=tpc+1
         ENDIF
         rd1=btgg(tpc)%pttn
         rd2=btgg(tpp)%pttn
         IF(rd1.lt.rd2)THEN
            exchh=btgg(tpc)
            btgg(tpc)=btgg(tpp)
            btgg(tpp)=exchh
            tpp=tpc
            tpc=2*tpp
         ELSE
            tpc=temp+1
         ENDIF
      ENDDO
      IF(tpc.eq.temp)THEN
         rd1=btgg(tpc)%pttn
         rd2=btgg(tpp)%pttn
         IF(rd1.lt.rd2)THEN
            exchh=btgg(tpc)
            btgg(tpc)=btgg(tpp)
            btgg(tpp)=exchh
         ENDIF
      ENDIF
   END SUBROUTINE downtree_ghost_point

   SUBROUTINE fouds1_all(iz,ix,ir)
      IMPLICIT NONE
      INTEGER :: i,j,k,ir,ix,iz,tsw1,swsol
      REAL(KIND=i10) :: trav,travm,slown,tdsh,tref
      REAL(KIND=i10) :: a,b,c,u,v,w,em,en,ri,risti
      REAL(KIND=i10) :: rd1,rd2,rd3
      tsw1=0
      slown=1.0/veln(iz,ix,ir)
      ri=rgor-(roffset-1)*dnr-(ir-1)*dnr+earth
      risti=ri*sin(rgox+(xoffset-1)*dnx+(ix-1)*dnx)
      DO i=ir-1,ir+1,2
         DO j=ix-1,ix+1,2
            DO k=iz-1,iz+1,2 
               IF(i.GE.1.AND.i.LE.rnum)THEN
                  IF(j.GE.1.AND.j.LE.xnum)THEN
                     IF(k.GE.1.AND.k.LE.znum)THEN
                        swsol=0
                        IF(nsts(iz,ix,i).EQ.0)THEN
                           swsol=1
                           IF(nsts(iz,j,ir).EQ.0.AND.nsts(k,ix,ir).EQ.0)THEN
                              u=dnr
                              v=ri*dnx
                              w=risti*dnz
                              em=ttn(iz,j,ir)-ttn(iz,ix,i)
                              en=ttn(k,ix,ir)-ttn(iz,ix,i)
                              rd1=v**2*w**2
                              rd2=u**2*w**2
                              rd3=u**2*v**2
                              a=rd1+rd2+rd3
                              b=-2.0*(em*rd2+en*rd3)
                              c=em**2*rd2+en**2*rd3-slown**2*u**2*rd1
                              tref=ttn(iz,ix,i)
                           ELSE IF(nsts(iz,j,ir).EQ.0)THEN
                              u=dnr
                              v=ri*dnx
                              em=ttn(iz,j,ir)-ttn(iz,ix,i)
                              a=u**2+v**2
                              b=-2.0*u**2*em
                              c=u**2*(em**2-v**2*slown**2)
                              tref=ttn(iz,ix,i)
                           ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                              u=dnr
                              v=risti*dnz
                              em=ttn(k,ix,ir)-ttn(iz,ix,i)
                              a=u**2+v**2
                              b=-2.0*u**2*em
                              c=u**2*(em**2-v**2*slown**2)
                              tref=ttn(iz,ix,i)
                           ELSE
                              a=1.0
                              b=0.0
                              c=-slown**2*dnr**2
                              tref=ttn(iz,ix,i)
                           ENDIF
                        ELSE IF(nsts(iz,j,ir).EQ.0)THEN
                           swsol=1
                           IF(nsts(k,ix,ir).EQ.0)THEN
                              u=ri*dnx
                              v=risti*dnz
                              em=ttn(k,ix,ir)-ttn(iz,j,ir)
                              a=u**2+v**2
                              b=-2.0*u**2*em
                              c=u**2*(em**2-v**2*slown**2)
                              tref=ttn(iz,j,ir)
                           ELSE
                              a=1.0
                              b=0.0
                              c=-slown**2*ri**2*dnx**2
                              tref=ttn(iz,j,ir)
                           ENDIF
                        ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                           swsol=1
                           a=1.0
                           b=0.0
                           c=-(slown*risti*dnz)**2
                           tref=ttn(k,ix,ir)
                        ENDIF
                        IF(swsol.EQ.1)THEN
                           rd1=b**2-4.0*a*c
                           IF(rd1.LT.0.0)rd1=0.0
                           tdsh=(-b+sqrt(rd1))/(2.0*a)
                           trav=tref+tdsh
                           IF(tsw1.EQ.1)THEN
                              travm=MIN(trav,travm)
                           ELSE
                              travm=trav
                              tsw1=1
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      IF(ttn(iz,ix,ir).EQ.0.OR.ttn(iz,ix,ir).GT.travm)THEN
         ttn(iz,ix,ir)=travm
      ENDIF
   END SUBROUTINE fouds1_all

   SUBROUTINE fouds2_all(iz,ix,ir)
      IMPLICIT NONE
      INTEGER :: i,j,k,i2,j2,k2,ir,ix,iz,tsw1
      INTEGER :: swi,swj,swk,swsol
      REAL(KIND=i10) :: trav,travm,slown,tdsh,tref,tdiv
      REAL(KIND=i10) :: a,b,c,u,v,w,em,en,ri,risti,rd1
      tsw1=0
      slown=1.0/veln(iz,ix,ir)
      ri=rgor-(roffset-1)*dnr-(ir-1)*dnr+earth
      risti=ri*sin(rgox+(xoffset-1)*dnx+(ix-1)*dnx)
      DO i=ir-1,ir+1,2
         swi=-1
         IF(i.EQ.ir-1)THEN
            i2=i-1
            IF(i2.GE.1)THEN
               IF(nsts(iz,ix,i2).EQ.0)swi=0
            ENDIF
         ELSE
            i2=i+1
            IF(i2.LE.rnum)THEN
               IF(nsts(iz,ix,i2).EQ.0)swi=0
            ENDIF
         ENDIF
         IF(nsts(iz,ix,i).EQ.0.AND.swi.EQ.0)THEN
            swi=-1
            IF(ttn(iz,ix,i).GT.ttn(iz,ix,i2))THEN
               swi=0
            ENDIF
         ELSE
            swi=-1
         ENDIF
         DO j=ix-1,ix+1,2
            swj=-1
            IF(j.eq.ix-1)THEN
               j2=j-1
               IF(j2.GE.1)THEN
                  IF(nsts(iz,j2,ir).EQ.0)swj=0
               ENDIF
            ELSE
               j2=j+1
               IF(j2.LE.xnum)THEN
                  IF(nsts(iz,j2,ir).EQ.0)swj=0
               ENDIF
            ENDIF
            IF(nsts(iz,j,ir).EQ.0.AND.swj.EQ.0)THEN
               swj=-1
               IF(ttn(iz,j,ir).GT.ttn(iz,j2,ir))THEN
                  swj=0
               ENDIF
            ELSE
               swj=-1
            ENDIF
            DO k=iz-1,iz+1,2
               swk=-1
               IF(k.eq.iz-1)THEN
                  k2=k-1
                  IF(k2.GE.1)THEN
                     IF(nsts(k2,ix,ir).EQ.0)swk=0
                  ENDIF
               ELSE
                  k2=k+1
                  IF(k2.LE.znum)THEN
                     IF(nsts(k2,ix,ir).EQ.0)swk=0
                  ENDIF
               ENDIF
               IF(nsts(k,ix,ir).EQ.0.AND.swk.EQ.0)THEN
                  swk=-1
                  IF(ttn(k,ix,ir).GT.ttn(k2,ix,ir))THEN
                     swk=0
                  ENDIF
               ELSE
                  swk=-1
               ENDIF
               IF(i.GE.1.AND.i.LE.rnum)THEN
                  IF(j.GE.1.AND.j.LE.xnum)THEN
                     IF(k.GE.1.AND.k.LE.znum)THEN
                        swsol=0
                        IF(swi.EQ.0)THEN
                           swsol=1
                           IF(swj.EQ.0)THEN
                              IF(swk.EQ.0)THEN
                                 u=2.0*dnr
                                 v=2.0*ri*dnx
                                 w=2.0*risti*dnz
                                 em=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn(iz,j,ir)
                                 em=em+ttn(iz,j2,ir)
                                 en=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn(k,ix,ir)
                                 en=en+ttn(k2,ix,ir)
                                 a=v**2*w**2+u**2*w**2+u**2*v**2
                                 b=2.0*u**2*(em*w**2+en*v**2)
                                 c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                                 tref=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)
                                 tdiv=3.0
                              ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                                 u=risti*dnz
                                 v=2.0*ri*dnx
                                 w=2.0*dnr
                                 em=3.0*ttn(k,ix,ir)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                                 en=3.0*ttn(k,ix,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                                 a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                                 b=6.0*u**2*(em*w**2+en*v**2)
                                 c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                                 tref=ttn(k,ix,ir)
                                 tdiv=1.0
                              ELSE
                                 u=2.0*dnr
                                 v=2.0*ri*dnx
                                 em=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn(iz,j,ir)
                                 em=em+ttn(iz,j2,ir)
                                 a=v**2+u**2
                                 b=2.0*em*u**2
                                 c=u**2*(em**2-slown**2*v**2) 
                                 tref=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)
                                 tdiv=3.0
                              ENDIF
                           ELSE IF(nsts(iz,j,ir).EQ.0)THEN
                              IF(swk.EQ.0)THEN
                                 u=ri*dnx
                                 v=2.0*dnr
                                 w=2.0*risti*dnz
                                 em=3.0*ttn(iz,j,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                                 en=3.0*ttn(iz,j,ir)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                                 a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                                 b=6.0*u**2*(em*w**2+en*v**2)
                                 c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                                 tref=ttn(iz,j,ir)
                                 tdiv=1.0
                              ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                                 u=ri*dnx
                                 v=2.0*dnr
                                 w=risti*dnz
                                 em=3.0*ttn(iz,j,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                                 en=ttn(iz,j,ir)-ttn(k,ix,ir)
                                 a=w**2*v**2+9.0*u**2*w**2+u**2*v**2
                                 b=u**2*(6.0*em*w**2+2.0*en*v**2)
                                 c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                                 tref=ttn(iz,j,ir)
                                 tdiv=1.0
                              ELSE
                                 u=ri*dnx
                                 v=2.0*dnr
                                 em=3.0*ttn(iz,j,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                                 a=v**2+9.0*u**2
                                 b=6.0*em*u**2
                                 c=u**2*(em**2-slown**2*v**2)
                                 tref=ttn(iz,j,ir)
                                 tdiv=1.0
                              ENDIF
                           ELSE
                              IF(swk.EQ.0)THEN
                                 u=2.0*dnr
                                 v=2.0*risti*dnz
                                 em=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn(k,ix,ir)
                                 em=em+ttn(k2,ix,ir)
                                 a=u**2+v**2
                                 b=2.0*em*u**2
                                 c=u**2*(em**2-v**2*slown**2)
                                 tref=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)
                                 tdiv=3.0
                              ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                                 u=risti*dnz
                                 v=2.0*dnr
                                 em=3.0*ttn(k,ix,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                                 a=v**2+9.0*u**2
                                 b=6.0*em*u**2
                                 c=u**2*(em**2-v**2*slown**2)
                                 tref=ttn(k,ix,ir)
                                 tdiv=1.0
                              ELSE
                                 u=2.0*dnr
                                 a=1.0
                                 b=0.0
                                 c=-u**2*slown**2
                                 tref=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)
                                 tdiv=3.0
                              ENDIF
                           ENDIF
                        ELSE IF(nsts(iz,ix,i).EQ.0)THEN
                           swsol=1
                           IF(swj.EQ.0)THEN
                              IF(swk.EQ.0)THEN
                                 u=dnr
                                 v=2.0*ri*dnx
                                 w=2.0*risti*dnz
                                 em=3.0*ttn(iz,ix,i)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                                 en=3.0*ttn(iz,ix,i)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                                 a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                                 b=6.0*u**2*(em*w**2+en*v**2)
                                 c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                                 tref=ttn(iz,ix,i)
                                 tdiv=1.0
                              ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                                 u=dnr
                                 v=2.0*ri*dnx
                                 w=risti*dnz
                                 em=3.0*ttn(iz,ix,i)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                                 en=ttn(iz,ix,i)-ttn(k,ix,ir)
                                 a=v**2*w**2+9.0*u**2*w**2+u**2*v**2
                                 b=u**2*(6.0*em*w**2+2.0*en*v**2)
                                 c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                                 tref=ttn(iz,ix,i)
                                 tdiv=1.0
                              ELSE
                                 u=dnr
                                 v=2.0*ri*dnx
                                 em=3.0*ttn(iz,ix,i)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                                 a=v**2+9.0*u**2
                                 b=6.0*em*u**2
                                 c=u**2*(em**2-v**2*slown**2)
                                 tref=ttn(iz,ix,i)
                                 tdiv=1.0
                              ENDIF
                           ELSE IF(nsts(iz,j,ir).EQ.0)THEN
                              IF(swk.EQ.0)THEN
                                 u=dnr
                                 v=ri*dnx
                                 w=2.0*risti*dnz
                                 em=ttn(iz,ix,i)-ttn(iz,j,ir)
                                 en=3.0*ttn(iz,ix,i)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                                 a=v**2*w**2+u**2*w**2+9.0*u**2*v**2
                                 b=u**2*(2.0*em*w**2+6.0*en*v**2)
                                 c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                                 tref=ttn(iz,ix,i)
                                 tdiv=1.0
                              ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                                 u=dnr
                                 v=ri*dnx
                                 w=risti*dnz
                                 em=ttn(iz,j,ir)-ttn(iz,ix,i)
                                 en=ttn(k,ix,ir)-ttn(iz,ix,i)
                                 a=v**2*w**2+u**2*w**2+u**2*v**2
                                 b=-2.0*u**2*(em*w**2+en*v**2)
                                 c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                                 tref=ttn(iz,ix,i)
                                 tdiv=1.0
                              ELSE
                                 u=dnr
                                 v=ri*dnx
                                 em=ttn(iz,j,ir)-ttn(iz,ix,i)
                                 a=u**2+v**2
                                 b=-2.0*u**2*em
                                 c=u**2*(em**2-v**2*slown**2)
                                 tref=ttn(iz,ix,i)
                                 tdiv=1.0
                              ENDIF
                           ELSE
                              IF(swk.EQ.0)THEN
                                 u=dnr
                                 v=2.0*risti*dnz
                                 em=3.0*ttn(iz,ix,i)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                                 a=v**2+9.0*u**2
                                 b=6.0*em*u**2
                                 c=u**2*(em**2-v**2*slown**2)
                                 tref=ttn(iz,ix,i)
                                 tdiv=1.0
                              ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                                 u=dnr
                                 v=risti*dnz
                                 em=ttn(k,ix,ir)-ttn(iz,ix,i)
                                 a=u**2+v**2
                                 b=-2.0*u**2*em
                                 c=u**2*(em**2-v**2*slown**2)
                                 tref=ttn(iz,ix,i)
                                 tdiv=1.0
                              ELSE
                                 a=1.0
                                 b=0.0
                                 c=-slown**2*dnr**2
                                 tref=ttn(iz,ix,i)
                                 tdiv=1.0
                              ENDIF
                           ENDIF
                        ELSE
                           IF(swj.EQ.0)THEN
                              swsol=1
                              IF(swk.EQ.0)THEN
                                 u=2.0*ri*dnx
                                 v=2.0*risti*dnz
                                 em=4.0*ttn(iz,j,ir)-ttn(iz,j2,ir)-4.0*ttn(k,ix,ir)
                                 em=em+ttn(k2,ix,ir)
                                 a=v**2+u**2
                                 b=2.0*em*u**2
                                 c=u**2*(em**2-slown**2*v**2)
                                 tref=4.0*ttn(iz,j,ir)-ttn(iz,j2,ir)
                                 tdiv=3.0
                              ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                                 u=risti*dnz
                                 v=2.0*ri*dnx
                                 em=3.0*ttn(k,ix,ir)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                                 a=v**2+9.0*u**2
                                 b=6.0*em*u**2
                                 c=u**2*(em**2-slown**2*v**2)
                                 tref=ttn(k,ix,ir)
                                 tdiv=1.0
                              ELSE
                                 u=2.0*ri*dnx
                                 a=1.0
                                 b=0.0
                                 c=-u**2*slown**2
                                 tref=4.0*ttn(iz,j,ir)-ttn(iz,j2,ir)
                                 tdiv=3.0
                              ENDIF
                           ELSE IF(nsts(iz,j,ir).EQ.0)THEN
                              swsol=1
                              IF(swk.EQ.0)THEN
                                 u=ri*dnx
                                 v=2.0*risti*dnz
                                 em=3.0*ttn(iz,j,ir)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                                 a=v**2+9.0*u**2
                                 b=6.0*em*u**2
                                 c=u**2*(em**2-v**2*slown**2)
                                 tref=ttn(iz,j,ir)
                                 tdiv=1.0
                              ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                                 u=ri*dnx
                                 v=risti*dnz
                                 em=ttn(k,ix,ir)-ttn(iz,j,ir)
                                 a=u**2+v**2
                                 b=-2.0*u**2*em
                                 c=u**2*(em**2-v**2*slown**2)
                                 tref=ttn(iz,j,ir)
                                 tdiv=1.0
                              ELSE
                                 a=1.0
                                 b=0.0
                                 c=-slown**2*ri**2*dnx**2
                                 tref=ttn(iz,j,ir)
                                 tdiv=1.0
                              ENDIF
                           ELSE
                              IF(swk.EQ.0)THEN
                                 swsol=1
                                 u=2.0*risti*dnz
                                 a=1.0
                                 b=0.0
                                 c=-u**2*slown**2
                                 tref=4.0*ttn(k,ix,ir)-ttn(k2,ix,ir)
                                 tdiv=3.0
                              ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                                 swsol=1
                                 a=1.0
                                 b=0.0
                                 c=-slown**2*risti**2*dnz**2
                                 tref=ttn(k,ix,ir)
                                 tdiv=1.0
                              ENDIF
                           ENDIF
                        ENDIF
                        IF(swsol.EQ.1)THEN
                           rd1=b**2-4.0*a*c
                           IF(rd1.LT.0.0)rd1=0.0
                           tdsh=(-b+sqrt(rd1))/(2.0*a)
                           trav=(tref+tdsh)/tdiv
                           IF(tsw1.EQ.1)THEN
                              travm=MIN(trav,travm)
                           ELSE
                              travm=trav
                              tsw1=1
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      IF(ttn(iz,ix,ir).EQ.0.OR.ttn(iz,ix,ir).GT.travm)THEN
         ttn(iz,ix,ir)=travm
      ENDIF
   END SUBROUTINE fouds2_all

   SUBROUTINE addtree_all(iz,ix,ir)
      IMPLICIT NONE
      INTEGER :: ix,iz,ir,tpp,tpc
      TYPE(backpointer) :: exch
      ntr=ntr+1
      nsts(iz,ix,ir)=ntr
      btg(ntr)%pr=ir
      btg(ntr)%px=ix
      btg(ntr)%pz=iz
      tpc=ntr
      tpp=tpc/2
      DO WHILE(tpp.gt.0)
         IF(ttn(iz,ix,ir).lt.ttn(btg(tpp)%pz,btg(tpp)%px,btg(tpp)%pr))THEN
            nsts(iz,ix,ir)=tpp
            nsts(btg(tpp)%pz,btg(tpp)%px,btg(tpp)%pr)=tpc
            exch=btg(tpc)
            btg(tpc)=btg(tpp)
            btg(tpp)=exch
            tpc=tpp
            tpp=tpc/2
         ELSE
            tpp=0
         ENDIF
      ENDDO
   END SUBROUTINE addtree_all

   SUBROUTINE downtree_all
      IMPLICIT NONE
      INTEGER :: tpp,tpc
      REAL(KIND=i10) :: rd1,rd2
      TYPE(backpointer) :: exch
      IF(ntr.EQ.1)THEN
         ntr=ntr-1
         RETURN
      ENDIF
      nsts(btg(ntr)%pz,btg(ntr)%px,btg(ntr)%pr)=1
      btg(1)=btg(ntr)
      ntr=ntr-1
      tpp=1
      tpc=2*tpp
      DO WHILE(tpc.lt.ntr)
         rd1=ttn(btg(tpc)%pz,btg(tpc)%px,btg(tpc)%pr)
         rd2=ttn(btg(tpc+1)%pz,btg(tpc+1)%px,btg(tpc+1)%pr)
         IF(rd1.gt.rd2)THEN
            tpc=tpc+1
         ENDIF
         rd1=ttn(btg(tpc)%pz,btg(tpc)%px,btg(tpc)%pr)
         rd2=ttn(btg(tpp)%pz,btg(tpp)%px,btg(tpp)%pr)
         IF(rd1.lt.rd2)THEN
            nsts(btg(tpp)%pz,btg(tpp)%px,btg(tpp)%pr)=tpc
            nsts(btg(tpc)%pz,btg(tpc)%px,btg(tpc)%pr)=tpp
            exch=btg(tpc)
            btg(tpc)=btg(tpp)
            btg(tpp)=exch
            tpp=tpc
            tpc=2*tpp
         ELSE
            tpc=ntr+1
         ENDIF
      ENDDO
      IF(tpc.eq.ntr)THEN
         rd1=ttn(btg(tpc)%pz,btg(tpc)%px,btg(tpc)%pr)
         rd2=ttn(btg(tpp)%pz,btg(tpp)%px,btg(tpp)%pr)
         IF(rd1.lt.rd2)THEN
            nsts(btg(tpp)%pz,btg(tpp)%px,btg(tpp)%pr)=tpc
            nsts(btg(tpc)%pz,btg(tpc)%px,btg(tpc)%pr)=tpp
            exch=btg(tpc)
            btg(tpc)=btg(tpp)
            btg(tpp)=exch
         ENDIF
      ENDIF
   END SUBROUTINE downtree_all

   SUBROUTINE updtree_all(iz,ix,ir)
      IMPLICIT NONE
      INTEGER :: ir,ix,iz,tpp,tpc
      TYPE(backpointer) :: exch
      tpc=nsts(iz,ix,ir)
      tpp=tpc/2
      DO WHILE(tpp.gt.0)
         IF(ttn(iz,ix,ir).lt.ttn(btg(tpp)%pz,btg(tpp)%px,btg(tpp)%pr))THEN
            nsts(iz,ix,ir)=tpp
            nsts(btg(tpp)%pz,btg(tpp)%px,btg(tpp)%pr)=tpc
            exch=btg(tpc)
            btg(tpc)=btg(tpp)
            btg(tpp)=exch
            tpc=tpp
            tpp=tpc/2
         ELSE
            tpp=0
         ENDIF
      ENDDO
   END SUBROUTINE updtree_all

   SUBROUTINE fouds2_left(iz,ix,ir,old_value)
      IMPLICIT NONE
      INTEGER :: i,j,k,i2,j2,k2,ir,ix,iz,tsw1
      INTEGER :: swi,swj,swk,swsol
      REAL(KIND=i10) :: old_value
      REAL(KIND=i10) :: trav,travm,slown,tdsh,tref,tdiv
      REAL(KIND=i10) :: a,b,c,u,v,w,em,en,ri,risti,rd1
      tsw1=0
      slown=1.0/veln(iz,ix,ir)
      ri=rgor-(roffset-1)*dnr-(ir-1)*dnr+earth
      risti=ri*sin(rgox+(xoffset-1)*dnx+(ix-1)*dnx)
      DO i=ir-1,ir+1,2
         swi=-1
         IF(i.eq.ir-1)THEN
            i2=i-1
            IF(i2.GE.1.AND.ttn(iz,ix,i2).NE.0)THEN
               swi=0
            ENDIF
         ELSE
            i2=i+1
            IF(i2.LE.rnum.AND.ttn(iz,ix,i2).NE.0)THEN
               swi=0
            ENDIF
         ENDIF
         IF((ttn(iz,ix,i).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(iz,ix,i).NE.0.AND.swi.EQ.0)THEN
            swi=-1
            IF(ttn(iz,ix,i).GT.ttn(iz,ix,i2))THEN
               swi=0
            ENDIF
         ELSE
            swi=-1
         ENDIF
         swj=-1
         IF(ttn_ghost_left(iz,2,ir).GT.ttn_ghost_left(iz,1,ir))THEN
            swj=0
         ENDIF
         DO k=iz-1,iz+1,2
            swk=-1
            IF(k.eq.iz-1)THEN
               k2=k-1
               IF(k2.GE.1.AND.ttn(k2,ix,ir).NE.0)THEN
                  swk=0
               ENDIF
            ELSE
               k2=k+1
               IF(k2.LE.znum.AND.ttn(k2,ix,ir).NE.0)THEN
                  swk=0
               ENDIF
            ENDIF
            IF((ttn(k,ix,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(k,ix,ir).NE.0.AND.swk.EQ.0)THEN
               swk=-1
               IF(ttn(k,ix,ir).GT.ttn(k2,ix,ir))THEN
                  swk=0
               ENDIF
            ELSE
               swk=-1
            ENDIF
            IF(i.GE.1.AND.i.LE.rnum)THEN
               IF(k.GE.1.AND.k.LE.znum)THEN
                  swsol=0
                  IF(swi.EQ.0)THEN
                     swsol=1
                     IF(swj.EQ.0)THEN
                        IF(swk.EQ.0)THEN
                           u=2.0*dnr
                           v=2.0*ri*dnx
                           w=2.0*risti*dnz
                           em=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn_ghost_left(iz,2,ir)
                           em=em+ttn_ghost_left(iz,1,ir)
                           en=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn(k,ix,ir)
                           en=en+ttn(k2,ix,ir)
                           a=v**2*w**2+u**2*w**2+u**2*v**2
                           b=2.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                           tref=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)
                           tdiv=3.0
                        ELSEIF((ttn(k,ix,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(k,ix,ir).NE.0)THEN
                           u=risti*dnz
                           v=2.0*ri*dnx
                           w=2.0*dnr
                           em=3.0*ttn(k,ix,ir)-4.0*ttn_ghost_left(iz,2,ir)+ttn_ghost_left(iz,1,ir)
                           en=3.0*ttn(k,ix,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                           a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                           b=6.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn(k,ix,ir)
                           tdiv=1.0
                        ELSE
                           u=2.0*dnr
                           v=2.0*ri*dnx
                           em=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn_ghost_left(iz,2,ir)
                           em=em+ttn_ghost_left(iz,1,ir)
                           a=v**2+u**2
                           b=2.0*em*u**2
                           c=u**2*(em**2-slown**2*v**2) 
                           tref=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)
                           tdiv=3.0
                        ENDIF
                     ELSE
                        IF(swk.EQ.0)THEN
                           u=ri*dnx
                           v=2.0*dnr
                           w=2.0*risti*dnz
                           em=3.0*ttn_ghost_left(iz,2,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                           en=3.0*ttn_ghost_left(iz,2,ir)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                           a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                           b=6.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn_ghost_left(iz,2,ir)
                           tdiv=1.0
                        ELSEIF((ttn(k,ix,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(k,ix,ir).NE.0)THEN
                           u=ri*dnx
                           v=2.0*dnr
                           w=risti*dnz
                           em=3.0*ttn_ghost_left(iz,2,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                           en=ttn_ghost_left(iz,2,ir)-ttn(k,ix,ir)
                           a=w**2*v**2+9.0*u**2*w**2+u**2*v**2
                           b=u**2*(6.0*em*w**2+2.0*en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn_ghost_left(iz,2,ir)
                           tdiv=1.0
                        ELSE
                           u=ri*dnx
                           v=2.0*dnr
                           em=3.0*ttn_ghost_left(iz,2,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                           a=v**2+9.0*u**2
                           b=6.0*em*u**2
                           c=u**2*(em**2-slown**2*v**2)
                           tref=ttn_ghost_left(iz,2,ir)
                           tdiv=1.0
                        ENDIF
                     ENDIF
                  ELSEIF((ttn(iz,ix,i).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(iz,ix,i).NE.0)THEN
                     swsol=1
                     IF(swj.EQ.0)THEN
                        IF(swk.EQ.0)THEN
                           u=dnr
                           v=2.0*ri*dnx
                           w=2.0*risti*dnz
                           em=3.0*ttn(iz,ix,i)-4.0*ttn_ghost_left(iz,2,ir)+ttn_ghost_left(iz,1,ir)
                           en=3.0*ttn(iz,ix,i)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                           a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                           b=6.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ELSEIF((ttn(k,ix,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(k,ix,ir).NE.0)THEN
                           u=dnr
                           v=2.0*ri*dnx
                           w=risti*dnz
                           em=3.0*ttn(iz,ix,i)-4.0*ttn_ghost_left(iz,2,ir)+ttn_ghost_left(iz,1,ir)
                           en=ttn(iz,ix,i)-ttn(k,ix,ir)
                           a=v**2*w**2+9.0*u**2*w**2+u**2*v**2
                           b=u**2*(6.0*em*w**2+2.0*en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ELSE
                           u=dnr
                           v=2.0*ri*dnx
                           em=3.0*ttn(iz,ix,i)-4.0*ttn_ghost_left(iz,2,ir)+ttn_ghost_left(iz,1,ir)
                           a=v**2+9.0*u**2
                           b=6.0*em*u**2
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ENDIF
                     ELSE
                        IF(swk.EQ.0)THEN
                           u=dnr
                           v=ri*dnx
                           w=2.0*risti*dnz
                           em=ttn(iz,ix,i)-ttn_ghost_left(iz,2,ir)
                           en=3.0*ttn(iz,ix,i)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                           a=v**2*w**2+u**2*w**2+9.0*u**2*v**2
                           b=u**2*(2.0*em*w**2+6.0*en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ELSEIF((ttn(k,ix,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(k,ix,ir).NE.0)THEN
                           u=dnr
                           v=ri*dnx
                           w=risti*dnz
                           em=ttn_ghost_left(iz,2,ir)-ttn(iz,ix,i)
                           en=ttn(k,ix,ir)-ttn(iz,ix,i)
                           a=v**2*w**2+u**2*w**2+u**2*v**2
                           b=-2.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ELSE
                           u=dnr
                           v=ri*dnx
                           em=ttn_ghost_left(iz,2,ir)-ttn(iz,ix,i)
                           a=u**2+v**2
                           b=-2.0*u**2*em
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ENDIF
                     ENDIF
                  ELSE
                     IF(swj.EQ.0)THEN
                        swsol=1
                        IF(swk.EQ.0)THEN
                           u=2.0*ri*dnx
                           v=2.0*risti*dnz
                           em=4.0*ttn_ghost_left(iz,2,ir)-ttn_ghost_left(iz,1,ir)-4.0*ttn(k,ix,ir)
                           em=em+ttn(k2,ix,ir)
                           a=v**2+u**2
                           b=2.0*em*u**2
                           c=u**2*(em**2-slown**2*v**2)
                           tref=4.0*ttn_ghost_left(iz,2,ir)-ttn_ghost_left(iz,1,ir)
                           tdiv=3.0
                        ELSEIF((ttn(k,ix,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(k,ix,ir).NE.0)THEN
                           u=risti*dnz
                           v=2.0*ri*dnx
                           em=3.0*ttn(k,ix,ir)-4.0*ttn_ghost_left(iz,2,ir)+ttn_ghost_left(iz,1,ir)
                           a=v**2+9.0*u**2
                           b=6.0*em*u**2
                           c=u**2*(em**2-slown**2*v**2)
                           tref=ttn(k,ix,ir)
                           tdiv=1.0
                        ELSE
                           u=2.0*ri*dnx
                           a=1.0
                           b=0.0
                           c=-u**2*slown**2
                           tref=4.0*ttn_ghost_left(iz,2,ir)-ttn_ghost_left(iz,1,ir)
                           tdiv=3.0
                        ENDIF
                     ELSE
                        swsol=1
                        IF(swk.EQ.0)THEN
                           u=ri*dnx
                           v=2.0*risti*dnz
                           em=3.0*ttn_ghost_left(iz,2,ir)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                           a=v**2+9.0*u**2
                           b=6.0*em*u**2
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn_ghost_left(iz,2,ir)
                           tdiv=1.0
                        ELSEIF((ttn(k,ix,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(k,ix,ir).NE.0)THEN
                           u=ri*dnx
                           v=risti*dnz
                           em=ttn(k,ix,ir)-ttn_ghost_left(iz,2,ir)
                           a=u**2+v**2
                           b=-2.0*u**2*em
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn_ghost_left(iz,2,ir)
                           tdiv=1.0
                        ELSE
                           a=1.0
                           b=0.0
                           c=-slown**2*ri**2*dnx**2
                           tref=ttn_ghost_left(iz,2,ir)
                           tdiv=1.0
                        ENDIF
                     ENDIF
                  ENDIF
   
                  IF(swsol.EQ.1)THEN
                     rd1=b**2-4.0*a*c
                     IF(rd1.LT.0.0)rd1=0.0
                     tdsh=(-b+sqrt(rd1))/(2.0*a)
                     trav=(tref+tdsh)/tdiv
                     IF(tsw1.EQ.1)THEN
                        travm=MIN(trav,travm)
                     ELSE
                        travm=trav
                        tsw1=1
                     ENDIF
                  ENDIF      
               ENDIF
            ENDIF
         ENDDO
         
      ENDDO
      IF(old_value.GE.travm.OR.old_value.EQ.0)THEN
         ttn(iz,ix,ir)=travm
         IF(nsts(iz,ix,ir).LE.0)THEN
            CALL addtree_all(iz,ix,ir)
         ELSE
            CALL updtree_all(iz,ix,ir)
         ENDIF
      ENDIF
   END SUBROUTINE fouds2_left

   SUBROUTINE fouds2_right(iz,ix,ir,old_value)
      IMPLICIT NONE
      INTEGER :: i,j,k,i2,j2,k2,ir,ix,iz,tsw1
      INTEGER :: swi,swj,swk,swsol
      REAL(KIND=i10) :: old_value
      REAL(KIND=i10) :: trav,travm,slown,tdsh,tref,tdiv
      REAL(KIND=i10) :: a,b,c,u,v,w,em,en,ri,risti,rd1
      tsw1=0
      slown=1.0/veln(iz,ix,ir)
      ri=rgor-(roffset-1)*dnr-(ir-1)*dnr+earth
      risti=ri*sin(rgox+(xoffset-1)*dnx+(ix-1)*dnx)
      DO i=ir-1,ir+1,2
         swi=-1
         IF(i.eq.ir-1)THEN
            i2=i-1
            IF(i2.GE.1.AND.ttn(iz,ix,i2).NE.0)THEN
               swi=0
            ENDIF
         ELSE
            i2=i+1
            IF(i2.LE.rnum.AND.ttn(iz,ix,i2).NE.0)THEN
               swi=0
            ENDIF
         ENDIF
         IF((ttn(iz,ix,i).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(iz,ix,i).NE.0.AND.swi.EQ.0)THEN
            swi=-1
            IF(ttn(iz,ix,i).GT.ttn(iz,ix,i2))THEN
               swi=0
            ENDIF
         ELSE
            swi=-1
         ENDIF
         swj=-1
         IF(ttn_ghost_right(iz,1,ir).GT.ttn_ghost_right(iz,2,ir))THEN
            swj=0
         ENDIF
         DO k=iz-1,iz+1,2
            swk=-1
            IF(k.eq.iz-1)THEN
               k2=k-1
               IF(k2.GE.1.AND.ttn(k2,ix,ir).NE.0)THEN
                  swk=0
               ENDIF
            ELSE
               k2=k+1
               IF(k2.LE.znum.AND.ttn(k2,ix,ir).NE.0)THEN
                  swk=0
               ENDIF
            ENDIF
            IF((ttn(k,ix,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(k,ix,ir).NE.0.AND.swk.EQ.0)THEN
               swk=-1
               IF(ttn(k,ix,ir).GT.ttn(k2,ix,ir))THEN
                  swk=0
               ENDIF
            ELSE
               swk=-1
            ENDIF
            IF(i.GE.1.AND.i.LE.rnum)THEN
               IF(k.GE.1.AND.k.LE.znum)THEN
                  swsol=0
                  IF(swi.EQ.0)THEN
                     swsol=1
                     IF(swj.EQ.0)THEN
                        IF(swk.EQ.0)THEN
                           u=2.0*dnr
                           v=2.0*ri*dnx
                           w=2.0*risti*dnz
                           em=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn_ghost_right(iz,1,ir)
                           em=em+ttn_ghost_right(iz,2,ir)
                           en=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn(k,ix,ir)
                           en=en+ttn(k2,ix,ir)
                           a=v**2*w**2+u**2*w**2+u**2*v**2
                           b=2.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                           tref=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)
                           tdiv=3.0
                        ELSEIF((ttn(k,ix,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(k,ix,ir).NE.0)THEN
                           u=risti*dnz
                           v=2.0*ri*dnx
                           w=2.0*dnr
                           em=3.0*ttn(k,ix,ir)-4.0*ttn_ghost_right(iz,1,ir)+ttn_ghost_right(iz,2,ir)
                           en=3.0*ttn(k,ix,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                           a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                           b=6.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn(k,ix,ir)
                           tdiv=1.0
                        ELSE
                           u=2.0*dnr
                           v=2.0*ri*dnx
                           em=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn_ghost_right(iz,1,ir)
                           em=em+ttn_ghost_right(iz,2,ir)
                           a=v**2+u**2
                           b=2.0*em*u**2
                           c=u**2*(em**2-slown**2*v**2) 
                           tref=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)
                           tdiv=3.0
                        ENDIF
                     ELSE
                        IF(swk.EQ.0)THEN
                           u=ri*dnx
                           v=2.0*dnr
                           w=2.0*risti*dnz
                           em=3.0*ttn_ghost_right(iz,1,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                           en=3.0*ttn_ghost_right(iz,1,ir)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                           a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                           b=6.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn_ghost_right(iz,1,ir)
                           tdiv=1.0
                        ELSEIF((ttn(k,ix,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(k,ix,ir).NE.0)THEN
                           u=ri*dnx
                           v=2.0*dnr
                           w=risti*dnz
                           em=3.0*ttn_ghost_right(iz,1,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                           en=ttn_ghost_right(iz,1,ir)-ttn(k,ix,ir)
                           a=w**2*v**2+9.0*u**2*w**2+u**2*v**2
                           b=u**2*(6.0*em*w**2+2.0*en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn_ghost_right(iz,1,ir)
                           tdiv=1.0
                        ELSE
                           u=ri*dnx
                           v=2.0*dnr
                           em=3.0*ttn_ghost_right(iz,1,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                           a=v**2+9.0*u**2
                           b=6.0*em*u**2
                           c=u**2*(em**2-slown**2*v**2)
                           tref=ttn_ghost_right(iz,1,ir)
                           tdiv=1.0
                        ENDIF
                     ENDIF
                  ELSEIF((ttn(iz,ix,i).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(iz,ix,i).NE.0)THEN
                     swsol=1
                     IF(swj.EQ.0)THEN
                        IF(swk.EQ.0)THEN
                           u=dnr
                           v=2.0*ri*dnx
                           w=2.0*risti*dnz
                           em=3.0*ttn(iz,ix,i)-4.0*ttn_ghost_right(iz,1,ir)+ttn_ghost_right(iz,2,ir)
                           en=3.0*ttn(iz,ix,i)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                           a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                           b=6.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ELSEIF((ttn(k,ix,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(k,ix,ir).NE.0)THEN
                           u=dnr
                           v=2.0*ri*dnx
                           w=risti*dnz
                           em=3.0*ttn(iz,ix,i)-4.0*ttn_ghost_right(iz,1,ir)+ttn_ghost_right(iz,2,ir)
                           en=ttn(iz,ix,i)-ttn(k,ix,ir)
                           a=v**2*w**2+9.0*u**2*w**2+u**2*v**2
                           b=u**2*(6.0*em*w**2+2.0*en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ELSE
                           u=dnr
                           v=2.0*ri*dnx
                           em=3.0*ttn(iz,ix,i)-4.0*ttn_ghost_right(iz,1,ir)+ttn_ghost_right(iz,2,ir)
                           a=v**2+9.0*u**2
                           b=6.0*em*u**2
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ENDIF
                     ELSE
                        IF(swk.EQ.0)THEN
                           u=dnr
                           v=ri*dnx
                           w=2.0*risti*dnz
                           em=ttn(iz,ix,i)-ttn_ghost_right(iz,1,ir)
                           en=3.0*ttn(iz,ix,i)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                           a=v**2*w**2+u**2*w**2+9.0*u**2*v**2
                           b=u**2*(2.0*em*w**2+6.0*en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ELSEIF((ttn(k,ix,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(k,ix,ir).NE.0)THEN
                           u=dnr
                           v=ri*dnx
                           w=risti*dnz
                           em=ttn_ghost_right(iz,1,ir)-ttn(iz,ix,i)
                           en=ttn(k,ix,ir)-ttn(iz,ix,i)
                           a=v**2*w**2+u**2*w**2+u**2*v**2
                           b=-2.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ELSE
                           u=dnr
                           v=ri*dnx
                           em=ttn_ghost_right(iz,1,ir)-ttn(iz,ix,i)
                           a=u**2+v**2
                           b=-2.0*u**2*em
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ENDIF
                     ENDIF
                  ELSE
                     IF(swj.EQ.0)THEN
                        swsol=1
                        IF(swk.EQ.0)THEN
                           u=2.0*ri*dnx
                           v=2.0*risti*dnz
                           em=4.0*ttn_ghost_right(iz,1,ir)-ttn_ghost_right(iz,2,ir)-4.0*ttn(k,ix,ir)
                           em=em+ttn(k2,ix,ir)
                           a=v**2+u**2
                           b=2.0*em*u**2
                           c=u**2*(em**2-slown**2*v**2)
                           tref=4.0*ttn_ghost_right(iz,1,ir)-ttn_ghost_right(iz,2,ir)
                           tdiv=3.0
                        ELSEIF((ttn(k,ix,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(k,ix,ir).NE.0)THEN
                           u=risti*dnz
                           v=2.0*ri*dnx
                           em=3.0*ttn(k,ix,ir)-4.0*ttn_ghost_right(iz,1,ir)+ttn_ghost_right(iz,2,ir)
                           a=v**2+9.0*u**2
                           b=6.0*em*u**2
                           c=u**2*(em**2-slown**2*v**2)
                           tref=ttn(k,ix,ir)
                           tdiv=1.0
                        ELSE
                           u=2.0*ri*dnx
                           a=1.0
                           b=0.0
                           c=-u**2*slown**2
                           tref=4.0*ttn_ghost_right(iz,1,ir)-ttn_ghost_right(iz,2,ir)
                           tdiv=3.0
                        ENDIF
                     ELSE
                        swsol=1
                        IF(swk.EQ.0)THEN
                           u=ri*dnx
                           v=2.0*risti*dnz
                           em=3.0*ttn_ghost_right(iz,1,ir)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                           a=v**2+9.0*u**2
                           b=6.0*em*u**2
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn_ghost_right(iz,1,ir)
                           tdiv=1.0
                        ELSEIF((ttn(k,ix,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(k,ix,ir).NE.0)THEN
                           u=ri*dnx
                           v=risti*dnz
                           em=ttn(k,ix,ir)-ttn_ghost_right(iz,1,ir)
                           a=u**2+v**2
                           b=-2.0*u**2*em
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn_ghost_right(iz,1,ir)
                           tdiv=1.0
                        ELSE
                           a=1.0
                           b=0.0
                           c=-slown**2*ri**2*dnx**2
                           tref=ttn_ghost_right(iz,1,ir)
                           tdiv=1.0
                        ENDIF
                     ENDIF
                  ENDIF
                  IF(swsol.EQ.1)THEN
                     rd1=b**2-4.0*a*c
                     IF(rd1.LT.0.0)rd1=0.0
                     tdsh=(-b+sqrt(rd1))/(2.0*a)
                     trav=(tref+tdsh)/tdiv
                     IF(tsw1.EQ.1)THEN
                        travm=MIN(trav,travm)
                     ELSE
                        travm=trav
                        tsw1=1
                     ENDIF
                  ENDIF      
               ENDIF
            ENDIF
         ENDDO 
      ENDDO
      IF(old_value.GE.travm.OR.old_value.EQ.0)THEN
         ttn(iz,ix,ir)=travm
         IF(nsts(iz,ix,ir).LE.0)THEN
            CALL addtree_all(iz,ix,ir)
         ELSE
            CALL updtree_all(iz,ix,ir)
         ENDIF
      ENDIF
   END SUBROUTINE fouds2_right
   
   SUBROUTINE fouds2_in(iz,ix,ir,old_value)
      IMPLICIT NONE
      INTEGER :: i,j,k,i2,j2,k2,ir,ix,iz,tsw1
      INTEGER :: swi,swj,swk,swsol
      INTEGER :: rtmp,xtmp,ztmp
      REAL(KIND=i10) :: old_value
      REAL(KIND=i10) :: trav,travm,slown,tdsh,tref,tdiv
      REAL(KIND=i10) :: a,b,c,u,v,w,em,en,ri,risti,rd1
      tsw1=0
      slown=1.0/veln(iz,ix,ir)
      ri=rgor-(roffset-1)*dnr-(ir-1)*dnr+earth
      risti=ri*sin(rgox+(xoffset-1)*dnx+(ix-1)*dnx)
      DO i=ir-1,ir+1,2
         swi=-1
         IF(i.eq.ir-1)THEN
            i2=i-1
            IF(i2.GE.1.AND.ttn(iz,ix,i2).NE.0)THEN
               swi=0
            ENDIF
         ELSE
            i2=i+1
            IF(i2.LE.rnum.AND.ttn(iz,ix,i2).NE.0)THEN
              swi=0
            ENDIF
         ENDIF
         IF((ttn(iz,ix,i).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(iz,ix,i).NE.0.AND.swi.EQ.0)THEN
            swi=-1
            IF(ttn(iz,ix,i).GT.ttn(iz,ix,i2))THEN
               swi=0
            ENDIF
         ELSE
            swi=-1
         ENDIF
         DO j=ix-1,ix+1,2
            swj=-1
            IF(j.eq.ix-1)THEN
               j2=j-1
               IF(j2.GE.1.AND.ttn(iz,j2,ir).NE.0)THEN
                  swj=0
               ENDIF
            ELSE
               j2=j+1
               IF(j2.LE.xnum.AND.ttn(iz,j2,ir).NE.0)THEN
                  swj=0
               ENDIF
            ENDIF
            IF((ttn(iz,j,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(iz,j,ir).NE.0.AND.swj.EQ.0)THEN
               swj=-1
               IF(ttn(iz,j,ir).GT.ttn(iz,j2,ir))THEN
                  swj=0
               ENDIF
            ELSE
               swj=-1
            ENDIF
            swk=-1
            IF(ttn_ghost_in(2,ix,ir).GT.ttn_ghost_in(1,ix,ir))THEN
               swk=0
            ENDIF
            IF(i.GE.1.AND.i.LE.rnum)THEN
               IF(j.GE.1.AND.j.LE.xnum)THEN
                  swsol=0
                  IF(swi.EQ.0)THEN
                     swsol=1
                     IF(swj.EQ.0)THEN
                        IF(swk.EQ.0)THEN
                           u=2.0*dnr
                           v=2.0*ri*dnx
                           w=2.0*risti*dnz
                           em=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn(iz,j,ir)
                           em=em+ttn(iz,j2,ir)
                           en=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn_ghost_in(2,ix,ir)
                           en=en+ttn_ghost_in(1,ix,ir)
                           a=v**2*w**2+u**2*w**2+u**2*v**2
                           b=2.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                           tref=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)
                           tdiv=3.0
                        ELSE
                           u=risti*dnz
                           v=2.0*ri*dnx
                           w=2.0*dnr
                           em=3.0*ttn_ghost_in(2,ix,ir)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                           en=3.0*ttn_ghost_in(2,ix,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                           a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                           b=6.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn_ghost_in(2,ix,ir)
                           tdiv=1.0
                        ENDIF
                     ELSEIF((ttn(iz,j,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(iz,j,ir).NE.0)THEN
                        IF(swk.EQ.0)THEN
                           u=ri*dnx
                           v=2.0*dnr
                           w=2.0*risti*dnz
                           em=3.0*ttn(iz,j,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                           en=3.0*ttn(iz,j,ir)-4.0*ttn_ghost_in(2,ix,ir)+ttn_ghost_in(1,ix,ir)
                           a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                           b=6.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn(iz,j,ir)
                           tdiv=1.0
                        ELSE
                           u=ri*dnx
                           v=2.0*dnr
                           w=risti*dnz
                           em=3.0*ttn(iz,j,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                           en=ttn(iz,j,ir)-ttn_ghost_in(2,ix,ir)
                           a=w**2*v**2+9.0*u**2*w**2+u**2*v**2
                           b=u**2*(6.0*em*w**2+2.0*en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn(iz,j,ir)
                           tdiv=1.0
                        ENDIF
                     ELSE
                        IF(swk.EQ.0)THEN
                           u=2.0*dnr
                           v=2.0*risti*dnz
                           em=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn_ghost_in(2,ix,ir)
                           em=em+ttn_ghost_in(1,ix,ir)
                           a=u**2+v**2
                           b=2.0*em*u**2
                           c=u**2*(em**2-v**2*slown**2)
                           tref=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)
                           tdiv=3.0
                        ELSE
                           u=risti*dnz
                           v=2.0*dnr
                           em=3.0*ttn_ghost_in(2,ix,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                           a=v**2+9.0*u**2
                           b=6.0*em*u**2
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn_ghost_in(2,ix,ir)
                           tdiv=1.0
                        ENDIF
                     ENDIF
                  ELSEIF((ttn(iz,ix,i).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(iz,ix,i).NE.0)THEN
                     swsol=1
                     IF(swj.EQ.0)THEN
                        IF(swk.EQ.0)THEN
                           u=dnr
                           v=2.0*ri*dnx
                           w=2.0*risti*dnz
                           em=3.0*ttn(iz,ix,i)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                           en=3.0*ttn(iz,ix,i)-4.0*ttn_ghost_in(2,ix,ir)+ttn_ghost_in(1,ix,ir)
                           a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                           b=6.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ELSE
                           u=dnr
                           v=2.0*ri*dnx
                           w=risti*dnz
                           em=3.0*ttn(iz,ix,i)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                           en=ttn(iz,ix,i)-ttn_ghost_in(2,ix,ir)
                           a=v**2*w**2+9.0*u**2*w**2+u**2*v**2
                           b=u**2*(6.0*em*w**2+2.0*en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ENDIF
                     ELSEIF((ttn(iz,j,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(iz,j,ir).NE.0)THEN
                        IF(swk.EQ.0)THEN
                           u=dnr
                           v=ri*dnx
                           w=2.0*risti*dnz
                           em=ttn(iz,ix,i)-ttn(iz,j,ir)
                           en=3.0*ttn(iz,ix,i)-4.0*ttn_ghost_in(2,ix,ir)+ttn_ghost_in(1,ix,ir)
                           a=v**2*w**2+u**2*w**2+9.0*u**2*v**2
                           b=u**2*(2.0*em*w**2+6.0*en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ELSE
                           u=dnr
                           v=ri*dnx
                           w=risti*dnz
                           em=ttn(iz,j,ir)-ttn(iz,ix,i)
                           en=ttn_ghost_in(2,ix,ir)-ttn(iz,ix,i)
                           a=v**2*w**2+u**2*w**2+u**2*v**2
                           b=-2.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ENDIF
                     ELSE
                        IF(swk.EQ.0)THEN
                           u=dnr
                           v=2.0*risti*dnz
                           em=3.0*ttn(iz,ix,i)-4.0*ttn_ghost_in(2,ix,ir)+ttn_ghost_in(1,ix,ir)
                           a=v**2+9.0*u**2
                           b=6.0*em*u**2
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ELSE 
                           u=dnr
                           v=risti*dnz
                           em=ttn_ghost_in(2,ix,ir)-ttn(iz,ix,i)
                           a=u**2+v**2
                           b=-2.0*u**2*em
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ENDIF
                     ENDIF
                  ELSE
                     IF(swj.EQ.0)THEN
                        swsol=1
                        IF(swk.EQ.0)THEN
                           u=2.0*ri*dnx
                           v=2.0*risti*dnz
                           em=4.0*ttn(iz,j,ir)-ttn(iz,j2,ir)-4.0*ttn_ghost_in(2,ix,ir)
                           em=em+ttn_ghost_in(1,ix,ir)
                           a=v**2+u**2
                           b=2.0*em*u**2
                           c=u**2*(em**2-slown**2*v**2)
                           tref=4.0*ttn(iz,j,ir)-ttn(iz,j2,ir)
                           tdiv=3.0
                        ELSE
                           u=risti*dnz
                           v=2.0*ri*dnx
                           em=3.0*ttn_ghost_in(2,ix,ir)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                           a=v**2+9.0*u**2
                           b=6.0*em*u**2
                           c=u**2*(em**2-slown**2*v**2)
                           tref=ttn_ghost_in(2,ix,ir)
                           tdiv=1.0
                        ENDIF
                     ELSEIF((ttn(iz,j,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(iz,j,ir).NE.0)THEN
                        swsol=1
                        IF(swk.EQ.0)THEN
                           u=ri*dnx
                           v=2.0*risti*dnz
                           em=3.0*ttn(iz,j,ir)-4.0*ttn_ghost_in(2,ix,ir)+ttn_ghost_in(1,ix,ir)
                           a=v**2+9.0*u**2
                           b=6.0*em*u**2
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn(iz,j,ir)
                           tdiv=1.0
                        ELSE
                           u=ri*dnx
                           v=risti*dnz
                           em=ttn_ghost_in(2,ix,ir)-ttn(iz,j,ir)
                           a=u**2+v**2
                           b=-2.0*u**2*em
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn(iz,j,ir)
                           tdiv=1.0
                        ENDIF
                     ELSE
                        IF(swk.EQ.0)THEN
                           swsol=1
                           u=2.0*risti*dnz
                           a=1.0
                           b=0.0
                           c=-u**2*slown**2
                           tref=4.0*ttn_ghost_in(2,ix,ir)-ttn_ghost_in(1,ix,ir)
                           tdiv=3.0
                        ELSE
                           swsol=1
                           a=1.0
                           b=0.0
                           c=-slown**2*risti**2*dnz**2
                           tref=ttn_ghost_in(2,ix,ir)
                           tdiv=1.0
                        ENDIF
                     ENDIF
                  ENDIF
                  IF(swsol.EQ.1)THEN
                     rd1=b**2-4.0*a*c
                     IF(rd1.LT.0.0)rd1=0.0
                     tdsh=(-b+sqrt(rd1))/(2.0*a)
                     trav=(tref+tdsh)/tdiv
                     IF(tsw1.EQ.1)THEN
                        travm=MIN(trav,travm)
                     ELSE
                        travm=trav
                        tsw1=1
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      IF(old_value.GE.travm.OR.old_value.EQ.0)THEN
         ttn(iz,ix,ir)=travm
         IF(nsts(iz,ix,ir).LE.0)THEN
            CALL addtree_all(iz,ix,ir)
         ELSE
            CALL updtree_all(iz,ix,ir)
         ENDIF
      ENDIF
   END SUBROUTINE fouds2_in

   SUBROUTINE fouds2_out(iz,ix,ir,old_value)
      IMPLICIT NONE
      INTEGER :: i,j,k,i2,j2,k2,ir,ix,iz,tsw1
      INTEGER :: swi,swj,swk,swsol
      REAL(KIND=i10) :: old_value
      REAL(KIND=i10) :: trav,travm,slown,tdsh,tref,tdiv
      REAL(KIND=i10) :: a,b,c,u,v,w,em,en,ri,risti,rd1
      tsw1=0
      slown=1.0/veln(iz,ix,ir)  
      ri=rgor-(roffset-1)*dnr-(ir-1)*dnr+earth
      risti=ri*sin(rgox+(xoffset-1)*dnx+(ix-1)*dnx)
      DO i=ir-1,ir+1,2
         swi=-1
         IF(i.eq.ir-1)THEN
            i2=i-1
            IF(i2.GE.1.AND.ttn(iz,ix,i2).NE.0)THEN
               swi=0
            ENDIF
         ELSE
            i2=i+1
            IF(i2.LE.rnum.AND.ttn(iz,ix,i2).NE.0)THEN
               swi=0
            ENDIF
         ENDIF
         IF((ttn(iz,ix,i).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(iz,ix,i).NE.0.AND.swi.EQ.0)THEN
            swi=-1
            IF(ttn(iz,ix,i).GT.ttn(iz,ix,i2))THEN
               swi=0
            ENDIF
         ELSE
            swi=-1
         ENDIF
         DO j=ix-1,ix+1,2
            swj=-1
            IF(j.eq.ix-1)THEN
               j2=j-1
               IF(j2.GE.1.AND.ttn(iz,j2,ir).NE.0)THEN
                  swj=0
               ENDIF
            ELSE
               j2=j+1
               IF(j2.LE.xnum.AND.ttn(iz,j2,ir).NE.0)THEN
                  swj=0
               ENDIF
            ENDIF
            IF((ttn(iz,j,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(iz,j,ir).NE.0.AND.swj.EQ.0)THEN
               swj=-1
               IF(ttn(iz,j,ir).GT.ttn(iz,j2,ir))THEN
                  swj=0
               ENDIF
            ELSE
               swj=-1
            ENDIF
            swk=-1
            IF(ttn_ghost_out(1,ix,ir).GT.ttn_ghost_out(2,ix,ir))THEN
               swk=0
            ENDIF   
            IF(i.GE.1.AND.i.LE.rnum)THEN
               IF(j.GE.1.AND.j.LE.xnum)THEN
                  swsol=0
                  IF(swi.EQ.0)THEN
                     swsol=1
                     IF(swj.EQ.0)THEN
                        IF(swk.EQ.0)THEN
                           u=2.0*dnr
                           v=2.0*ri*dnx
                           w=2.0*risti*dnz
                           em=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn(iz,j,ir)
                           em=em+ttn(iz,j2,ir)
                           en=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn_ghost_out(1,ix,ir)
                           en=en+ttn_ghost_out(2,ix,ir)
                           a=v**2*w**2+u**2*w**2+u**2*v**2
                           b=2.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                           tref=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)
                           tdiv=3.0
                        ELSE
                           u=risti*dnz
                           v=2.0*ri*dnx
                           w=2.0*dnr
                           em=3.0*ttn_ghost_out(1,ix,ir)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                           en=3.0*ttn_ghost_out(1,ix,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                           a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                           b=6.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn_ghost_out(1,ix,ir)
                           tdiv=1.0
                        ENDIF
                     ELSEIF((ttn(iz,j,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(iz,j,ir).NE.0)THEN
                        IF(swk.EQ.0)THEN
                           u=ri*dnx
                           v=2.0*dnr
                           w=2.0*risti*dnz
                           em=3.0*ttn(iz,j,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                           en=3.0*ttn(iz,j,ir)-4.0*ttn_ghost_out(1,ix,ir)+ttn_ghost_out(2,ix,ir)
                           a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                           b=6.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn(iz,j,ir)
                           tdiv=1.0
                        ELSE
                           u=ri*dnx
                           v=2.0*dnr
                           w=risti*dnz
                           em=3.0*ttn(iz,j,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                           en=ttn(iz,j,ir)-ttn_ghost_out(1,ix,ir)
                           a=w**2*v**2+9.0*u**2*w**2+u**2*v**2
                           b=u**2*(6.0*em*w**2+2.0*en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn(iz,j,ir)
                           tdiv=1.0
                        ENDIF
                     ELSE
                        IF(swk.EQ.0)THEN
                           u=2.0*dnr
                           v=2.0*risti*dnz
                           em=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn_ghost_out(1,ix,ir)
                           em=em+ttn_ghost_out(2,ix,ir)
                           a=u**2+v**2
                           b=2.0*em*u**2
                           c=u**2*(em**2-v**2*slown**2)
                           tref=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)
                           tdiv=3.0
                        ELSE
                           u=risti*dnz
                           v=2.0*dnr
                           em=3.0*ttn_ghost_out(1,ix,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                           a=v**2+9.0*u**2
                           b=6.0*em*u**2
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn_ghost_out(1,ix,ir)
                           tdiv=1.0
                        ENDIF
                     ENDIF
                  ELSEIF((ttn(iz,ix,i).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(iz,ix,i).NE.0)THEN
                     swsol=1
                     IF(swj.EQ.0)THEN
                        IF(swk.EQ.0)THEN
                           u=dnr
                           v=2.0*ri*dnx
                           w=2.0*risti*dnz
                           em=3.0*ttn(iz,ix,i)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                           en=3.0*ttn(iz,ix,i)-4.0*ttn_ghost_out(1,ix,ir)+ttn_ghost_out(2,ix,ir)
                           a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                           b=6.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ELSE
                           u=dnr
                           v=2.0*ri*dnx
                           w=risti*dnz
                           em=3.0*ttn(iz,ix,i)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                           en=ttn(iz,ix,i)-ttn_ghost_out(1,ix,ir)
                           a=v**2*w**2+9.0*u**2*w**2+u**2*v**2
                           b=u**2*(6.0*em*w**2+2.0*en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ENDIF
                     ELSEIF((ttn(iz,j,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(iz,j,ir).NE.0)THEN
                        IF(swk.EQ.0)THEN
                           u=dnr
                           v=ri*dnx
                           w=2.0*risti*dnz
                           em=ttn(iz,ix,i)-ttn(iz,j,ir)
                           en=3.0*ttn(iz,ix,i)-4.0*ttn_ghost_out(1,ix,ir)+ttn_ghost_out(2,ix,ir)
                           a=v**2*w**2+u**2*w**2+9.0*u**2*v**2
                           b=u**2*(2.0*em*w**2+6.0*en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ELSE
                           u=dnr
                           v=ri*dnx
                           w=risti*dnz
                           em=ttn(iz,j,ir)-ttn(iz,ix,i)
                           en=ttn_ghost_out(1,ix,ir)-ttn(iz,ix,i)
                           a=v**2*w**2+u**2*w**2+u**2*v**2
                           b=-2.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ENDIF
                     ELSE
                        IF(swk.EQ.0)THEN
                           u=dnr
                           v=2.0*risti*dnz
                           em=3.0*ttn(iz,ix,i)-4.0*ttn_ghost_out(1,ix,ir)+ttn_ghost_out(2,ix,ir)
                           a=v**2+9.0*u**2
                           b=6.0*em*u**2
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ELSE 
                           u=dnr
                           v=risti*dnz
                           em=ttn_ghost_out(1,ix,ir)-ttn(iz,ix,i)
                           a=u**2+v**2
                           b=-2.0*u**2*em
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ENDIF
                     ENDIF
                  ELSE
                     IF(swj.EQ.0)THEN
                        swsol=1
                        IF(swk.EQ.0)THEN
                           u=2.0*ri*dnx
                           v=2.0*risti*dnz
                           em=4.0*ttn(iz,j,ir)-ttn(iz,j2,ir)-4.0*ttn_ghost_out(1,ix,ir)
                           em=em+ttn_ghost_out(2,ix,ir)
                           a=v**2+u**2
                           b=2.0*em*u**2
                           c=u**2*(em**2-slown**2*v**2)
                           tref=4.0*ttn(iz,j,ir)-ttn(iz,j2,ir)
                           tdiv=3.0
                        ELSE
                           u=risti*dnz
                           v=2.0*ri*dnx
                           em=3.0*ttn_ghost_out(1,ix,ir)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                           a=v**2+9.0*u**2
                           b=6.0*em*u**2
                           c=u**2*(em**2-slown**2*v**2)
                           tref=ttn_ghost_out(1,ix,ir)
                           tdiv=1.0
                        ENDIF
                     ELSEIF((ttn(iz,j,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(iz,j,ir).NE.0)THEN
                        swsol=1
                        IF(swk.EQ.0)THEN
                           u=ri*dnx
                           v=2.0*risti*dnz
                           em=3.0*ttn(iz,j,ir)-4.0*ttn_ghost_out(1,ix,ir)+ttn_ghost_out(2,ix,ir)
                           a=v**2+9.0*u**2
                           b=6.0*em*u**2
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn(iz,j,ir)
                           tdiv=1.0
                        ELSE
                           u=ri*dnx
                           v=risti*dnz
                           em=ttn_ghost_out(1,ix,ir)-ttn(iz,j,ir)
                           a=u**2+v**2
                           b=-2.0*u**2*em
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn(iz,j,ir)
                           tdiv=1.0
                        ENDIF
                     ELSE
                        IF(swk.EQ.0)THEN
                           swsol=1
                           u=2.0*risti*dnz
                           a=1.0
                           b=0.0
                           c=-u**2*slown**2
                           tref=4.0*ttn_ghost_out(1,ix,ir)-ttn_ghost_out(2,ix,ir)
                           tdiv=3.0
                        ELSE
                           swsol=1
                           a=1.0
                           b=0.0
                           c=-slown**2*risti**2*dnz**2
                           tref=ttn_ghost_out(1,ix,ir)
                           tdiv=1.0
                        ENDIF
                     ENDIF
                  ENDIF
                  IF(swsol.EQ.1)THEN
                     rd1=b**2-4.0*a*c
                     IF(rd1.LT.0.0)rd1=0.0
                     tdsh=(-b+sqrt(rd1))/(2.0*a)
                     trav=(tref+tdsh)/tdiv
                     IF(tsw1.EQ.1)THEN
                        travm=MIN(trav,travm)
                     ELSE
                        travm=trav
                        tsw1=1
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      IF(old_value.GE.travm.OR.old_value.EQ.0)THEN
         ttn(iz,ix,ir)=travm
         IF(nsts(iz,ix,ir).LE.0)THEN
            CALL addtree_all(iz,ix,ir)
         ELSE
            CALL updtree_all(iz,ix,ir)
         ENDIF
      ENDIF
   END SUBROUTINE fouds2_out

   SUBROUTINE fouds2_up(iz,ix,ir,old_value)
      IMPLICIT NONE
      INTEGER :: i,j,k,i2,j2,k2,ir,ix,iz,tsw1
      INTEGER :: swi,swj,swk,swsol
      REAL(KIND=i10) :: old_value
      REAL(KIND=i10) :: trav,travm,slown,tdsh,tref,tdiv
      REAL(KIND=i10) :: a,b,c,u,v,w,em,en,ri,risti,rd1
      tsw1=0
      slown=1.0/veln(iz,ix,ir) 
      ri=rgor-(roffset-1)*dnr-(ir-1)*dnr+earth
      risti=ri*sin(rgox+(xoffset-1)*dnx+(ix-1)*dnx)
      swi=-1
      IF(ttn_ghost_up(iz,ix,2).GT.ttn_ghost_up(iz,ix,1))THEN
         swi=0
      ENDIF
      DO j=ix-1,ix+1,2
         swj=-1
         IF(j.EQ.ix-1)THEN
            j2=j-1
            IF(j2.GE.1.AND.ttn(iz,j2,ir).NE.0)THEN
               swj=0
            ENDIF
         ELSE
            j2=j+1
            IF(j2.LE.xnum.AND.ttn(iz,j2,ir).NE.0)THEN
               swj=0
            ENDIF
         ENDIF
         IF((ttn(iz,j,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(iz,j,ir).NE.0.AND.swj.EQ.0)THEN
            swj=-1
            IF(ttn(iz,j,ir).GT.ttn(iz,j2,ir))THEN
               swj=0
            ENDIF
         ELSE
            swj=-1
         ENDIF
         DO k=iz-1,iz+1,2
            swk=-1
            IF(k.eq.iz-1)THEN
               k2=k-1
               IF(k2.GE.1.AND.ttn(k2,ix,ir).NE.0)THEN
                  swk=0
               ENDIF
            ELSE
               k2=k+1
               IF(k2.LE.znum.AND.ttn(k2,ix,ir).NE.0)THEN
                  swk=0
               ENDIF
            ENDIF
            IF((ttn(k,ix,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(k,ix,ir).NE.0.AND.swk.EQ.0)THEN
               swk=-1
               IF(ttn(k,ix,ir).GT.ttn(k2,ix,ir))THEN
                  swk=0
               ENDIF
            ELSE
               swk=-1
            ENDIF
            IF(j.GE.1.AND.j.LE.xnum)THEN
               IF(k.GE.1.AND.k.LE.znum)THEN
                  swsol=0
                  IF(swi.EQ.0)THEN
                     swsol=1
                     IF(swj.EQ.0)THEN
                        IF(swk.EQ.0)THEN
                           u=2.0*dnr
                           v=2.0*ri*dnx
                           w=2.0*risti*dnz
                           em=4.0*ttn_ghost_up(iz,ix,2)-ttn_ghost_up(iz,ix,1)-4.0*ttn(iz,j,ir)
                           em=em+ttn(iz,j2,ir)
                           en=4.0*ttn_ghost_up(iz,ix,2)-ttn_ghost_up(iz,ix,1)-4.0*ttn(k,ix,ir)
                           en=en+ttn(k2,ix,ir)
                           a=v**2*w**2+u**2*w**2+u**2*v**2
                           b=2.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                           tref=4.0*ttn_ghost_up(iz,ix,2)-ttn_ghost_up(iz,ix,1)
                           tdiv=3.0
                        ELSEIF((ttn(k,ix,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(k,ix,ir).NE.0)THEN
                           u=risti*dnz
                           v=2.0*ri*dnx
                           w=2.0*dnr
                           em=3.0*ttn(k,ix,ir)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                           en=3.0*ttn(k,ix,ir)-4.0*ttn_ghost_up(iz,ix,2)+ttn_ghost_up(iz,ix,1)
                           a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                           b=6.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn(k,ix,ir)
                           tdiv=1.0
                        ELSE
                           u=2.0*dnr
                           v=2.0*ri*dnx
                           em=4.0*ttn_ghost_up(iz,ix,2)-ttn_ghost_up(iz,ix,1)-4.0*ttn(iz,j,ir)
                           em=em+ttn(iz,j2,ir)
                           a=v**2+u**2
                           b=2.0*em*u**2
                           c=u**2*(em**2-slown**2*v**2) 
                           tref=4.0*ttn_ghost_up(iz,ix,2)-ttn_ghost_up(iz,ix,1)
                           tdiv=3.0
                        ENDIF
                     ELSEIF((ttn(iz,j,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(iz,j,ir).NE.0)THEN
                        IF(swk.EQ.0)THEN
                           u=ri*dnx
                           v=2.0*dnr
                           w=2.0*risti*dnz
                           em=3.0*ttn(iz,j,ir)-4.0*ttn_ghost_up(iz,ix,2)+ttn_ghost_up(iz,ix,1)
                           en=3.0*ttn(iz,j,ir)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                           a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                           b=6.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn(iz,j,ir)
                           tdiv=1.0
                        ELSEIF((ttn(k,ix,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(k,ix,ir).NE.0)THEN
                           u=ri*dnx
                           v=2.0*dnr
                           w=risti*dnz
                           em=3.0*ttn(iz,j,ir)-4.0*ttn_ghost_up(iz,ix,2)+ttn_ghost_up(iz,ix,1)
                           en=ttn(iz,j,ir)-ttn(k,ix,ir)
                           a=w**2*v**2+9.0*u**2*w**2+u**2*v**2
                           b=u**2*(6.0*em*w**2+2.0*en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn(iz,j,ir)
                           tdiv=1.0
                        ELSE
                           u=ri*dnx
                           v=2.0*dnr
                           em=3.0*ttn(iz,j,ir)-4.0*ttn_ghost_up(iz,ix,2)+ttn_ghost_up(iz,ix,1)
                           a=v**2+9.0*u**2
                           b=6.0*em*u**2
                           c=u**2*(em**2-slown**2*v**2)
                           tref=ttn(iz,j,ir)
                           tdiv=1.0
                        ENDIF
                     ELSE
                        IF(swk.EQ.0)THEN
                           u=2.0*dnr
                           v=2.0*risti*dnz
                           em=4.0*ttn_ghost_up(iz,ix,2)-ttn_ghost_up(iz,ix,1)-4.0*ttn(k,ix,ir)
                           em=em+ttn(k2,ix,ir)
                           a=u**2+v**2
                           b=2.0*em*u**2
                           c=u**2*(em**2-v**2*slown**2)
                           tref=4.0*ttn_ghost_up(iz,ix,2)-ttn_ghost_up(iz,ix,1)
                           tdiv=3.0
                        ELSEIF((ttn(k,ix,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(k,ix,ir).NE.0)THEN
                           u=risti*dnz
                           v=2.0*dnr
                           em=3.0*ttn(k,ix,ir)-4.0*ttn_ghost_up(iz,ix,2)+ttn_ghost_up(iz,ix,1)
                           a=v**2+9.0*u**2
                           b=6.0*em*u**2
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn(k,ix,ir)
                           tdiv=1.0
                        ELSE
                           u=2.0*dnr
                           a=1.0
                           b=0.0
                           c=-u**2*slown**2
                           tref=4.0*ttn_ghost_up(iz,ix,2)-ttn_ghost_up(iz,ix,1)
                           tdiv=3.0
                        ENDIF
                     ENDIF
                  ELSE
                     swsol=1
                     IF(swj.EQ.0)THEN
                        IF(swk.EQ.0)THEN
                           u=dnr
                           v=2.0*ri*dnx
                           w=2.0*risti*dnz
                           em=3.0*ttn_ghost_up(iz,ix,2)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                           en=3.0*ttn_ghost_up(iz,ix,2)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                           a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                           b=6.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn_ghost_up(iz,ix,2)
                           tdiv=1.0
                        ELSEIF((ttn(k,ix,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(k,ix,ir).NE.0)THEN
                           u=dnr
                           v=2.0*ri*dnx
                           w=risti*dnz
                           em=3.0*ttn_ghost_up(iz,ix,2)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                           en=ttn_ghost_up(iz,ix,2)-ttn(k,ix,ir)
                           a=v**2*w**2+9.0*u**2*w**2+u**2*v**2
                           b=u**2*(6.0*em*w**2+2.0*en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                           tref=ttn_ghost_up(iz,ix,2)
                           tdiv=1.0
                        ELSE
                           u=dnr
                           v=2.0*ri*dnx
                           em=3.0*ttn_ghost_up(iz,ix,2)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                           a=v**2+9.0*u**2
                           b=6.0*em*u**2
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn_ghost_up(iz,ix,2)
                           tdiv=1.0
                        ENDIF
                     ELSEIF((ttn(iz,j,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(iz,j,ir).NE.0)THEN
                        IF(swk.EQ.0)THEN
                           u=dnr
                           v=ri*dnx
                           w=2.0*risti*dnz
                           em=ttn_ghost_up(iz,ix,2)-ttn(iz,j,ir)
                           en=3.0*ttn_ghost_up(iz,ix,2)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                           a=v**2*w**2+u**2*w**2+9.0*u**2*v**2
                           b=u**2*(2.0*em*w**2+6.0*en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                           tref=ttn_ghost_up(iz,ix,2)
                           tdiv=1.0
                        ELSEIF((ttn(k,ix,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(k,ix,ir).NE.0)THEN
                           u=dnr
                           v=ri*dnx
                           w=risti*dnz
                           em=ttn(iz,j,ir)-ttn_ghost_up(iz,ix,2)
                           en=ttn(k,ix,ir)-ttn_ghost_up(iz,ix,2)
                           a=v**2*w**2+u**2*w**2+u**2*v**2
                           b=-2.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn_ghost_up(iz,ix,2)
                           tdiv=1.0
                        ELSE
                           u=dnr
                           v=ri*dnx
                           em=ttn(iz,j,ir)-ttn_ghost_up(iz,ix,2)
                           a=u**2+v**2
                           b=-2.0*u**2*em
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn_ghost_up(iz,ix,2)
                           tdiv=1.0
                        ENDIF
                     ELSE
                        IF(swk.EQ.0)THEN
                           u=dnr
                           v=2.0*risti*dnz
                           em=3.0*ttn_ghost_up(iz,ix,2)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                           a=v**2+9.0*u**2
                           b=6.0*em*u**2
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn_ghost_up(iz,ix,2)
                           tdiv=1.0
                        ELSEIF((ttn(k,ix,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(k,ix,ir).NE.0)THEN
                           u=dnr
                           v=risti*dnz
                           em=ttn(k,ix,ir)-ttn_ghost_up(iz,ix,2)
                           a=u**2+v**2
                           b=-2.0*u**2*em
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn_ghost_up(iz,ix,2)
                           tdiv=1.0
                        ELSE
                           a=1.0
                           b=0.0
                           c=-slown**2*dnr**2
                           tref=ttn_ghost_up(iz,ix,2)
                           tdiv=1.0
                        ENDIF
                     ENDIF
                  ENDIF
                  IF(swsol.EQ.1)THEN
                     rd1=b**2-4.0*a*c
                     IF(rd1.LT.0.0)rd1=0.0
                     tdsh=(-b+sqrt(rd1))/(2.0*a)
                     trav=(tref+tdsh)/tdiv
                     IF(tsw1.EQ.1)THEN
                        travm=MIN(trav,travm)
                     ELSE
                        travm=trav
                        tsw1=1
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      IF(old_value.GE.travm.OR.old_value.EQ.0)THEN
         ttn(iz,ix,ir)=travm
         IF(nsts(iz,ix,ir).LE.0)THEN
            CALL addtree_all(iz,ix,ir)
         ELSE
            CALL updtree_all(iz,ix,ir)
         ENDIF
      ENDIF
   END SUBROUTINE fouds2_up

   SUBROUTINE fouds2_down(iz,ix,ir,old_value)
      IMPLICIT NONE
      INTEGER :: i,j,k,i2,j2,k2,ir,ix,iz,tsw1
      INTEGER :: swi,swj,swk,swsol
      REAL(KIND=i10) :: old_value
      REAL(KIND=i10) :: trav,travm,slown,tdsh,tref,tdiv
      REAL(KIND=i10) :: a,b,c,u,v,w,em,en,ri,risti,rd1
      tsw1=0
      slown=1.0/veln(iz,ix,ir)  
      ri=rgor-(roffset-1)*dnr-(ir-1)*dnr+earth
      risti=ri*sin(rgox+(xoffset-1)*dnx+(ix-1)*dnx)
      swi=-1
      IF(ttn_ghost_down(iz,ix,1).GT.ttn_ghost_down(iz,ix,2))THEN
         swi=0
      ENDIF
      DO j=ix-1,ix+1,2
         swj=-1
         IF(j.EQ.ix-1)THEN
            j2=j-1
            IF(j2.GE.1.AND.ttn(iz,j2,ir).NE.0)THEN
               swj=0
            ENDIF
         ELSE
            j2=j+1
            IF(j2.LE.xnum.AND.ttn(iz,j2,ir).NE.0)THEN
               swj=0
            ENDIF
         ENDIF
         IF((ttn(iz,j,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(iz,j,ir).NE.0.AND.swj.EQ.0)THEN
            swj=-1
            IF(ttn(iz,j,ir).GT.ttn(iz,j2,ir))THEN
               swj=0
            ENDIF
         ELSE
            swj=-1
         ENDIF
         DO k=iz-1,iz+1,2
            swk=-1
            IF(k.EQ.iz-1)THEN
               k2=k-1
               IF(k2.GE.1.AND.ttn(k2,ix,ir).NE.0)THEN
                  swk=0
               ENDIF
            ELSE
               k2=k+1
               IF(k2.LE.znum.AND.ttn(k2,ix,ir).NE.0)THEN
                  swk=0
               ENDIF
            ENDIF
            IF((ttn(k,ix,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(k,ix,ir).NE.0.AND.swk.EQ.0)THEN
               swk=-1
               IF(ttn(k,ix,ir).GT.ttn(k2,ix,ir))THEN
                  swk=0
               ENDIF
            ELSE
               swk=-1
            ENDIF
            IF(j.GE.1.AND.j.LE.xnum)THEN
               IF(k.GE.1.AND.k.LE.znum)THEN
                  swsol=0
                  IF(swi.EQ.0)THEN
                     swsol=1
                     IF(swj.EQ.0)THEN
                        IF(swk.EQ.0)THEN
                           u=2.0*dnr
                           v=2.0*ri*dnx
                           w=2.0*risti*dnz
                           em=4.0*ttn_ghost_down(iz,ix,1)-ttn_ghost_down(iz,ix,2)-4.0*ttn(iz,j,ir)
                           em=em+ttn(iz,j2,ir)
                           en=4.0*ttn_ghost_down(iz,ix,1)-ttn_ghost_down(iz,ix,2)-4.0*ttn(k,ix,ir)
                           en=en+ttn(k2,ix,ir)
                           a=v**2*w**2+u**2*w**2+u**2*v**2
                           b=2.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                           tref=4.0*ttn_ghost_down(iz,ix,1)-ttn_ghost_down(iz,ix,2)
                           tdiv=3.0
                        ELSEIF((ttn(k,ix,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(k,ix,ir).NE.0)THEN
                           u=risti*dnz
                           v=2.0*ri*dnx
                           w=2.0*dnr
                           em=3.0*ttn(k,ix,ir)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                           en=3.0*ttn(k,ix,ir)-4.0*ttn_ghost_down(iz,ix,1)+ttn_ghost_down(iz,ix,2)
                           a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                           b=6.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn(k,ix,ir)
                           tdiv=1.0
                        ELSE
                           u=2.0*dnr
                           v=2.0*ri*dnx
                           em=4.0*ttn_ghost_down(iz,ix,1)-ttn_ghost_down(iz,ix,2)-4.0*ttn(iz,j,ir)
                           em=em+ttn(iz,j2,ir)
                           a=v**2+u**2
                           b=2.0*em*u**2
                           c=u**2*(em**2-slown**2*v**2) 
                           tref=4.0*ttn_ghost_down(iz,ix,1)-ttn_ghost_down(iz,ix,2)
                           tdiv=3.0
                        ENDIF
                     ELSEIF((ttn(iz,j,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(iz,j,ir).NE.0)THEN
                        IF(swk.EQ.0)THEN
                           u=ri*dnx
                           v=2.0*dnr
                           w=2.0*risti*dnz
                           em=3.0*ttn(iz,j,ir)-4.0*ttn_ghost_down(iz,ix,1)+ttn_ghost_down(iz,ix,2)
                           en=3.0*ttn(iz,j,ir)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                           a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                           b=6.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn(iz,j,ir)
                           tdiv=1.0
                        ELSEIF((ttn(k,ix,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(k,ix,ir).NE.0)THEN
                           u=ri*dnx
                           v=2.0*dnr
                           w=risti*dnz
                           em=3.0*ttn(iz,j,ir)-4.0*ttn_ghost_down(iz,ix,1)+ttn_ghost_down(iz,ix,2)
                           en=ttn(iz,j,ir)-ttn(k,ix,ir)
                           a=w**2*v**2+9.0*u**2*w**2+u**2*v**2
                           b=u**2*(6.0*em*w**2+2.0*en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn(iz,j,ir)
                           tdiv=1.0
                        ELSE
                           u=ri*dnx
                           v=2.0*dnr
                           em=3.0*ttn(iz,j,ir)-4.0*ttn_ghost_down(iz,ix,1)+ttn_ghost_down(iz,ix,2)
                           a=v**2+9.0*u**2
                           b=6.0*em*u**2
                           c=u**2*(em**2-slown**2*v**2)
                           tref=ttn(iz,j,ir)
                           tdiv=1.0
                        ENDIF
                     ELSE
                        IF(swk.EQ.0)THEN
                           u=2.0*dnr
                           v=2.0*risti*dnz
                           em=4.0*ttn_ghost_down(iz,ix,1)-ttn_ghost_down(iz,ix,2)-4.0*ttn(k,ix,ir)
                           em=em+ttn(k2,ix,ir)
                           a=u**2+v**2
                           b=2.0*em*u**2
                           c=u**2*(em**2-v**2*slown**2)
                           tref=4.0*ttn_ghost_down(iz,ix,1)-ttn_ghost_down(iz,ix,2)
                           tdiv=3.0
                        ELSEIF((ttn(k,ix,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(k,ix,ir).NE.0)THEN
                           u=risti*dnz
                           v=2.0*dnr
                           em=3.0*ttn(k,ix,ir)-4.0*ttn_ghost_down(iz,ix,1)+ttn_ghost_down(iz,ix,2)
                           a=v**2+9.0*u**2
                           b=6.0*em*u**2
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn(k,ix,ir)
                           tdiv=1.0
                        ELSE
                           u=2.0*dnr
                           a=1.0
                           b=0.0
                           c=-u**2*slown**2
                           tref=4.0*ttn_ghost_down(iz,ix,1)-ttn_ghost_down(iz,ix,2)
                           tdiv=3.0
                        ENDIF
                     ENDIF
                  ELSE
                     swsol=1
                     IF(swj.EQ.0)THEN
                        IF(swk.EQ.0)THEN
                           u=dnr
                           v=2.0*ri*dnx
                           w=2.0*risti*dnz
                           em=3.0*ttn_ghost_down(iz,ix,1)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                           en=3.0*ttn_ghost_down(iz,ix,1)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                           a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                           b=6.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn_ghost_down(iz,ix,1)
                           tdiv=1.0
                        ELSEIF((ttn(k,ix,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(k,ix,ir).NE.0)THEN
                           u=dnr
                           v=2.0*ri*dnx
                           w=risti*dnz
                           em=3.0*ttn_ghost_down(iz,ix,1)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                           en=ttn_ghost_down(iz,ix,1)-ttn(k,ix,ir)
                           a=v**2*w**2+9.0*u**2*w**2+u**2*v**2
                           b=u**2*(6.0*em*w**2+2.0*en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                           tref=ttn_ghost_down(iz,ix,1)
                           tdiv=1.0
                        ELSE
                           u=dnr
                           v=2.0*ri*dnx
                           em=3.0*ttn_ghost_down(iz,ix,1)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                           a=v**2+9.0*u**2
                           b=6.0*em*u**2
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn_ghost_down(iz,ix,1)
                           tdiv=1.0
                        ENDIF
                     ELSEIF((ttn(iz,j,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(iz,j,ir).NE.0)THEN
                        IF(swk.EQ.0)THEN
                           u=dnr
                           v=ri*dnx
                           w=2.0*risti*dnz
                           em=ttn_ghost_down(iz,ix,1)-ttn(iz,j,ir)
                           en=3.0*ttn_ghost_down(iz,ix,1)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                           a=v**2*w**2+u**2*w**2+9.0*u**2*v**2
                           b=u**2*(2.0*em*w**2+6.0*en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                           tref=ttn_ghost_down(iz,ix,1)
                           tdiv=1.0
                        ELSEIF((ttn(k,ix,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(k,ix,ir).NE.0)THEN
                           u=dnr
                           v=ri*dnx
                           w=risti*dnz
                           em=ttn(iz,j,ir)-ttn_ghost_down(iz,ix,1)
                           en=ttn(k,ix,ir)-ttn_ghost_down(iz,ix,1)
                           a=v**2*w**2+u**2*w**2+u**2*v**2
                           b=-2.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn_ghost_down(iz,ix,1)
                           tdiv=1.0
                        ELSE
                           u=dnr
                           v=ri*dnx
                           em=ttn(iz,j,ir)-ttn_ghost_down(iz,ix,1)
                           a=u**2+v**2
                           b=-2.0*u**2*em
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn_ghost_down(iz,ix,1)
                           tdiv=1.0
                        ENDIF
                     ELSE
                        IF(swk.EQ.0)THEN
                           u=dnr
                           v=2.0*risti*dnz
                           em=3.0*ttn_ghost_down(iz,ix,1)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                           a=v**2+9.0*u**2
                           b=6.0*em*u**2
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn_ghost_down(iz,ix,1)
                           tdiv=1.0
                        ELSEIF((ttn(k,ix,ir).LT.ttn(iz,ix,ir).OR.ttn(iz,ix,ir).EQ.0).AND.ttn(k,ix,ir).NE.0)THEN
                           u=dnr
                           v=risti*dnz
                           em=ttn(k,ix,ir)-ttn_ghost_down(iz,ix,1)
                           a=u**2+v**2
                           b=-2.0*u**2*em
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn_ghost_down(iz,ix,1)
                           tdiv=1.0
                        ELSE
                           a=1.0
                           b=0.0
                           c=-slown**2*dnr**2
                           tref=ttn_ghost_down(iz,ix,1)
                           tdiv=1.0
                        ENDIF
                     ENDIF
                  ENDIF
                  IF(swsol.EQ.1)THEN
                     rd1=b**2-4.0*a*c
                     IF(rd1.LT.0.0)rd1=0.0
                     tdsh=(-b+sqrt(rd1))/(2.0*a)
                     trav=(tref+tdsh)/tdiv
                     IF(tsw1.EQ.1)THEN
                        travm=MIN(trav,travm)
                     ELSE
                        travm=trav
                        tsw1=1
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      IF(old_value.GT.travm.OR.old_value.EQ.0)THEN
         ttn(iz,ix,ir)=travm
         IF(nsts(iz,ix,ir).LE.0)THEN
            CALL addtree_all(iz,ix,ir)
         ELSE
            CALL updtree_all(iz,ix,ir)
         ENDIF
      ENDIF
   END SUBROUTINE fouds2_down

   SUBROUTINE fouds2_interior_left(iz,ix,ir)
      IMPLICIT NONE
      INTEGER :: i,j,k,i2,j2,k2,ir,ix,iz,tsw1
      INTEGER :: swi,swj,swk,swsol
      REAL(KIND=i10) :: trav,travm,slown,tdsh,tref,tdiv
      REAL(KIND=i10) :: a,b,c,u,v,w,em,en,ri,risti,rd1
      tsw1=0
      slown=1.0/veln(iz,ix,ir)
      ri=rgor-(roffset-1)*dnr-(ir-1)*dnr+earth
      risti=ri*sin(rgox+(xoffset-1)*dnx+(ix-1)*dnx)
      DO i=ir-1,ir+1,2
         swi=-1
         IF(i.EQ.ir-1)THEN
            i2=i-1
            IF(i2.GE.1)THEN
               IF(nsts(iz,ix,i2).EQ.0)swi=0
            ENDIF
         ELSE
            i2=i+1
            IF(i2.LE.rnum)THEN
               IF(nsts(iz,ix,i2).EQ.0)swi=0
            ENDIF
         ENDIF
         IF(nsts(iz,ix,i).EQ.0.AND.swi.EQ.0)THEN
            swi=-1
            IF(ttn(iz,ix,i).GT.ttn(iz,ix,i2))THEN
               swi=0
            ENDIF
         ELSE
            swi=-1
         ENDIF
         DO j=ix-1,ix+1,2
            swj=-1
            IF(j.eq.ix-1)THEN
               j2=j-1
               IF(j2.GE.1)THEN
                  IF(nsts(iz,j2,ir).EQ.0)swj=0
               ENDIF
            ELSE
               j2=j+1
               IF(j2.LE.xnum)THEN
                  IF(nsts(iz,j2,ir).EQ.0)swj=0
               ENDIF
            ENDIF
            IF(nsts(iz,j,ir).EQ.0.AND.swj.EQ.0)THEN
               swj=-1
               IF(ttn(iz,j,ir).GT.ttn(iz,j2,ir))THEN
                  swj=0
               ENDIF
            ELSE
               swj=-1
            ENDIF
            DO k=iz-1,iz+1,2
               swk=-1
               IF(k.eq.iz-1)THEN
                  k2=k-1
                  IF(k2.GE.1)THEN
                     IF(nsts(k2,ix,ir).EQ.0)swk=0
                  ENDIF
               ELSE
                  k2=k+1
                  IF(k2.LE.znum)THEN
                     IF(nsts(k2,ix,ir).EQ.0)swk=0
                  ENDIF
               ENDIF
               IF(nsts(k,ix,ir).EQ.0.AND.swk.EQ.0)THEN
                  swk=-1
                  IF(ttn(k,ix,ir).GT.ttn(k2,ix,ir))THEN
                     swk=0
                  ENDIF
               ELSE
                  swk=-1
               ENDIF
               IF(i.GE.1.AND.i.LE.rnum)THEN
                  IF(j.GE.1.AND.j.LE.xnum)THEN
                     IF(k.GE.1.AND.k.LE.znum)THEN
                        swsol=0
                        IF(swi.EQ.0)THEN
                           swsol=1
                           IF(swj.EQ.0)THEN
                              IF(swk.EQ.0)THEN
                                 u=2.0*dnr
                                 v=2.0*ri*dnx
                                 w=2.0*risti*dnz
                                 em=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn(iz,j,ir)
                                 em=em+ttn(iz,j2,ir)
                                 en=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn(k,ix,ir)
                                 en=en+ttn(k2,ix,ir)
                                 a=v**2*w**2+u**2*w**2+u**2*v**2
                                 b=2.0*u**2*(em*w**2+en*v**2)
                                 c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                                 tref=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)
                                 tdiv=3.0
                              ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                                 u=risti*dnz
                                 v=2.0*ri*dnx
                                 w=2.0*dnr
                                 em=3.0*ttn(k,ix,ir)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                                 en=3.0*ttn(k,ix,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                                 a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                                 b=6.0*u**2*(em*w**2+en*v**2)
                                 c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                                 tref=ttn(k,ix,ir)
                                 tdiv=1.0
                              ELSE
                                 u=2.0*dnr
                                 v=2.0*ri*dnx
                                 em=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn(iz,j,ir)
                                 em=em+ttn(iz,j2,ir)
                                 a=v**2+u**2
                                 b=2.0*em*u**2
                                 c=u**2*(em**2-slown**2*v**2) 
                                 tref=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)
                                 tdiv=3.0
                              ENDIF
                           ELSE IF(nsts(iz,j,ir).EQ.0)THEN
                              IF(swk.EQ.0)THEN
                                 u=ri*dnx
                                 v=2.0*dnr
                                 w=2.0*risti*dnz
                                 em=3.0*ttn(iz,j,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                                 en=3.0*ttn(iz,j,ir)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                                 a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                                 b=6.0*u**2*(em*w**2+en*v**2)
                                 c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                                 tref=ttn(iz,j,ir)
                                 tdiv=1.0
                              ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                                 u=ri*dnx
                                 v=2.0*dnr
                                 w=risti*dnz
                                 em=3.0*ttn(iz,j,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                                 en=ttn(iz,j,ir)-ttn(k,ix,ir)
                                 a=w**2*v**2+9.0*u**2*w**2+u**2*v**2
                                 b=u**2*(6.0*em*w**2+2.0*en*v**2)
                                 c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                                 tref=ttn(iz,j,ir)
                                 tdiv=1.0
                              ELSE
                                 u=ri*dnx
                                 v=2.0*dnr
                                 em=3.0*ttn(iz,j,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                                 a=v**2+9.0*u**2
                                 b=6.0*em*u**2
                                 c=u**2*(em**2-slown**2*v**2)
                                 tref=ttn(iz,j,ir)
                                 tdiv=1.0
                              ENDIF
                           ELSE
                              IF(swk.EQ.0)THEN
                                 u=2.0*dnr
                                 v=2.0*risti*dnz
                                 em=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn(k,ix,ir)
                                 em=em+ttn(k2,ix,ir)
                                 a=u**2+v**2
                                 b=2.0*em*u**2
                                 c=u**2*(em**2-v**2*slown**2)
                                 tref=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)
                                 tdiv=3.0
                              ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                                 u=risti*dnz
                                 v=2.0*dnr
                                 em=3.0*ttn(k,ix,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                                 a=v**2+9.0*u**2
                                 b=6.0*em*u**2
                                 c=u**2*(em**2-v**2*slown**2)
                                 tref=ttn(k,ix,ir)
                                 tdiv=1.0
                              ELSE
                                 u=2.0*dnr
                                 a=1.0
                                 b=0.0
                                 c=-u**2*slown**2
                                 tref=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)
                                 tdiv=3.0
                              ENDIF
                           ENDIF
                        ELSE IF(nsts(iz,ix,i).EQ.0)THEN
                           swsol=1
                           IF(swj.EQ.0)THEN
                              IF(swk.EQ.0)THEN
                                 u=dnr
                                 v=2.0*ri*dnx
                                 w=2.0*risti*dnz
                                 em=3.0*ttn(iz,ix,i)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                                 en=3.0*ttn(iz,ix,i)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                                 a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                                 b=6.0*u**2*(em*w**2+en*v**2)
                                 c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                                 tref=ttn(iz,ix,i)
                                 tdiv=1.0
                              ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                                 u=dnr
                                 v=2.0*ri*dnx
                                 w=risti*dnz
                                 em=3.0*ttn(iz,ix,i)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                                 en=ttn(iz,ix,i)-ttn(k,ix,ir)
                                 a=v**2*w**2+9.0*u**2*w**2+u**2*v**2
                                 b=u**2*(6.0*em*w**2+2.0*en*v**2)
                                 c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                                 tref=ttn(iz,ix,i)
                                 tdiv=1.0
                              ELSE
                                 u=dnr
                                 v=2.0*ri*dnx
                                 em=3.0*ttn(iz,ix,i)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                                 a=v**2+9.0*u**2
                                 b=6.0*em*u**2
                                 c=u**2*(em**2-v**2*slown**2)
                                 tref=ttn(iz,ix,i)
                                 tdiv=1.0
                              ENDIF
                           ELSE IF(nsts(iz,j,ir).EQ.0)THEN
                              IF(swk.EQ.0)THEN
                                 u=dnr
                                 v=ri*dnx
                                 w=2.0*risti*dnz
                                 em=ttn(iz,ix,i)-ttn(iz,j,ir)
                                 en=3.0*ttn(iz,ix,i)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                                 a=v**2*w**2+u**2*w**2+9.0*u**2*v**2
                                 b=u**2*(2.0*em*w**2+6.0*en*v**2)
                                 c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                                 tref=ttn(iz,ix,i)
                                 tdiv=1.0
                              ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                                 u=dnr
                                 v=ri*dnx
                                 w=risti*dnz
                                 em=ttn(iz,j,ir)-ttn(iz,ix,i)
                                 en=ttn(k,ix,ir)-ttn(iz,ix,i)
                                 a=v**2*w**2+u**2*w**2+u**2*v**2
                                 b=-2.0*u**2*(em*w**2+en*v**2)
                                 c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                                 tref=ttn(iz,ix,i)
                                 tdiv=1.0
                              ELSE
                                 u=dnr
                                 v=ri*dnx
                                 em=ttn(iz,j,ir)-ttn(iz,ix,i)
                                 a=u**2+v**2
                                 b=-2.0*u**2*em
                                 c=u**2*(em**2-v**2*slown**2)
                                 tref=ttn(iz,ix,i)
                                 tdiv=1.0
                              ENDIF
                           ELSE
                              IF(swk.EQ.0)THEN
                                 u=dnr
                                 v=2.0*risti*dnz
                                 em=3.0*ttn(iz,ix,i)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                                 a=v**2+9.0*u**2
                                 b=6.0*em*u**2
                                 c=u**2*(em**2-v**2*slown**2)
                                 tref=ttn(iz,ix,i)
                                 tdiv=1.0
                              ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                                 u=dnr
                                 v=risti*dnz
                                 em=ttn(k,ix,ir)-ttn(iz,ix,i)
                                 a=u**2+v**2
                                 b=-2.0*u**2*em
                                 c=u**2*(em**2-v**2*slown**2)
                                 tref=ttn(iz,ix,i)
                                 tdiv=1.0
                              ELSE
                                 a=1.0
                                 b=0.0
                                 c=-slown**2*dnr**2
                                 tref=ttn(iz,ix,i)
                                 tdiv=1.0
                              ENDIF
                           ENDIF
                        ELSE
                           IF(swj.EQ.0)THEN
                              swsol=1
                              IF(swk.EQ.0)THEN
                                 u=2.0*ri*dnx
                                 v=2.0*risti*dnz
                                 em=4.0*ttn(iz,j,ir)-ttn(iz,j2,ir)-4.0*ttn(k,ix,ir)
                                 em=em+ttn(k2,ix,ir)
                                 a=v**2+u**2
                                 b=2.0*em*u**2
                                 c=u**2*(em**2-slown**2*v**2)
                                 tref=4.0*ttn(iz,j,ir)-ttn(iz,j2,ir)
                                 tdiv=3.0
                              ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                                 u=risti*dnz
                                 v=2.0*ri*dnx
                                 em=3.0*ttn(k,ix,ir)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                                 a=v**2+9.0*u**2
                                 b=6.0*em*u**2
                                 c=u**2*(em**2-slown**2*v**2)
                                 tref=ttn(k,ix,ir)
                                 tdiv=1.0
                              ELSE
                                 u=2.0*ri*dnx
                                 a=1.0
                                 b=0.0
                                 c=-u**2*slown**2
                                 tref=4.0*ttn(iz,j,ir)-ttn(iz,j2,ir)
                                 tdiv=3.0
                              ENDIF
                           ELSE IF(nsts(iz,j,ir).EQ.0)THEN
                              swsol=1
                              IF(swk.EQ.0)THEN
                                 u=ri*dnx
                                 v=2.0*risti*dnz
                                 em=3.0*ttn(iz,j,ir)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                                 a=v**2+9.0*u**2
                                 b=6.0*em*u**2
                                 c=u**2*(em**2-v**2*slown**2)
                                 tref=ttn(iz,j,ir)
                                 tdiv=1.0
                              ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                                 u=ri*dnx
                                 v=risti*dnz
                                 em=ttn(k,ix,ir)-ttn(iz,j,ir)
                                 a=u**2+v**2
                                 b=-2.0*u**2*em
                                 c=u**2*(em**2-v**2*slown**2)
                                 tref=ttn(iz,j,ir)
                                 tdiv=1.0
                              ELSE
                                 a=1.0
                                 b=0.0
                                 c=-slown**2*ri**2*dnx**2
                                 tref=ttn(iz,j,ir)
                                 tdiv=1.0
                              ENDIF
                           ELSE
                              IF(swk.EQ.0)THEN
                                 swsol=1
                                 u=2.0*risti*dnz
                                 a=1.0
                                 b=0.0
                                 c=-u**2*slown**2
                                 tref=4.0*ttn(k,ix,ir)-ttn(k2,ix,ir)
                                 tdiv=3.0
                              ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                                 swsol=1
                                 a=1.0
                                 b=0.0
                                 c=-slown**2*risti**2*dnz**2
                                 tref=ttn(k,ix,ir)
                                 tdiv=1.0
                              ENDIF
                           ENDIF
                        ENDIF
                        IF(swsol.EQ.1)THEN
                           rd1=b**2-4.0*a*c
                           IF(rd1.LT.0.0)rd1=0.0
                           tdsh=(-b+sqrt(rd1))/(2.0*a)
                           trav=(tref+tdsh)/tdiv
                           IF(tsw1.EQ.1)THEN
                              travm=MIN(trav,travm)
                           ELSE
                              travm=trav
                              tsw1=1
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      IF(ttn(iz,ix,ir).EQ.0.OR.ttn(iz,ix,ir).GT.travm)THEN
         ttn(iz,ix,ir)=travm
      ENDIF
   END SUBROUTINE fouds2_interior_left   

   SUBROUTINE fouds2_left_in(iz,ix,ir)
      IMPLICIT NONE
      INTEGER :: i,j,k,i2,j2,k2,ir,ix,iz,tsw1
      INTEGER :: swi,swj,swk,swsol
      REAL(KIND=i10) :: old_value
      REAL(KIND=i10) :: trav,travm,slown,tdsh,tref,tdiv
      REAL(KIND=i10) :: a,b,c,u,v,w,em,en,ri,risti,rd1
      tsw1=0
      slown=1.0/veln(iz,ix,ir)  
      ri=rgor-(roffset-1)*dnr-(ir-1)*dnr+earth
      risti=ri*sin(rgox+(xoffset-1)*dnx+(ix-1)*dnx)
      DO i=ir-1,ir+1,2
         swi=-1
         IF(i.eq.ir-1)THEN
            i2=i-1
            IF(i2.GE.1.AND.ttn(iz,ix,i2).NE.0)THEN
               swi=0
            ENDIF
         ELSE
            i2=i+1
            IF(i2.LE.rnum.AND.ttn(iz,ix,i2).NE.0)THEN
               swi=0
            ENDIF
         ENDIF
         IF(ttn(iz,ix,i).LT.ttn(iz,ix,ir).AND.swi.EQ.0.AND.ttn(iz,ix,i).NE.0)THEN
            swi=-1
            IF(ttn(iz,ix,i).GT.ttn(iz,ix,i2))THEN
               swi=0
            ENDIF
         ELSE
            swi=-1
         ENDIF
         swj=-1
         IF(ttn_ghost_left(iz,2,ir).GT.ttn_ghost_left(iz,1,ir))THEN
            swj=0
         ENDIF
         swk=-1
         IF(ttn_ghost_in(2,ix,ir).GT.ttn_ghost_in(1,ix,ir))THEN
            swk=0
         ENDIF
         IF(i.GE.1.AND.i.LE.rnum)THEN
            swsol=0
            IF(swi.EQ.0)THEN
               swsol=1
               IF(swj.EQ.0)THEN
                  IF(swk.EQ.0)THEN
                     u=2.0*dnr
                     v=2.0*ri*dnx
                     w=2.0*risti*dnz
                     em=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn_ghost_left(iz,2,ir)
                     em=em+ttn_ghost_left(iz,1,ir)
                     en=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn_ghost_in(2,ix,ir)
                     en=en+ttn_ghost_in(1,ix,ir)
                     a=v**2*w**2+u**2*w**2+u**2*v**2
                     b=2.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)
                     tdiv=3.0
                  ELSE
                     u=risti*dnz
                     v=2.0*ri*dnx
                     w=2.0*dnr
                     em=3.0*ttn_ghost_in(2,ix,ir)-4.0*ttn_ghost_left(iz,2,ir)+ttn_ghost_left(iz,1,ir)
                     en=3.0*ttn_ghost_in(2,ix,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_in(2,ix,ir)
                     tdiv=1.0
                  ENDIF
               ELSE
                  IF(swk.EQ.0)THEN
                     u=ri*dnx
                     v=2.0*dnr
                     w=2.0*risti*dnz
                     em=3.0*ttn_ghost_left(iz,2,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                     en=3.0*ttn_ghost_left(iz,2,ir)-4.0*ttn_ghost_in(2,ix,ir)+ttn_ghost_in(1,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_left(iz,2,ir)
                     tdiv=1.0
                  ELSE
                     u=ri*dnx
                     v=2.0*dnr
                     w=risti*dnz
                     em=3.0*ttn_ghost_left(iz,2,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                     en=ttn_ghost_left(iz,2,ir)-ttn_ghost_in(2,ix,ir)
                     a=w**2*v**2+9.0*u**2*w**2+u**2*v**2
                     b=u**2*(6.0*em*w**2+2.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_left(iz,2,ir)
                     tdiv=1.0
                  ENDIF
               ENDIF
            ELSE IF(ttn(iz,ix,i).LT.ttn(iz,ix,ir).AND.ttn(iz,ix,i).NE.0)THEN
               swsol=1
               IF(swj.EQ.0)THEN
                  IF(swk.EQ.0)THEN
                     u=dnr
                     v=2.0*ri*dnx
                     w=2.0*risti*dnz
                     em=3.0*ttn(iz,ix,i)-4.0*ttn_ghost_left(iz,2,ir)+ttn_ghost_left(iz,1,ir)
                     en=3.0*ttn(iz,ix,i)-4.0*ttn_ghost_in(2,ix,ir)+ttn_ghost_in(1,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn(iz,ix,i)
                     tdiv=1.0
                  ELSE
                     u=dnr
                     v=2.0*ri*dnx
                     w=risti*dnz
                     em=3.0*ttn(iz,ix,i)-4.0*ttn_ghost_left(iz,2,ir)+ttn_ghost_left(iz,1,ir)
                     en=ttn(iz,ix,i)-ttn_ghost_in(2,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+u**2*v**2
                     b=u**2*(6.0*em*w**2+2.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=ttn(iz,ix,i)
                     tdiv=1.0
                  ENDIF
               ELSE
                  IF(swk.EQ.0)THEN
                     u=dnr
                     v=ri*dnx
                     w=2.0*risti*dnz
                     em=ttn(iz,ix,i)-ttn_ghost_left(iz,2,ir)
                     en=3.0*ttn(iz,ix,i)-4.0*ttn_ghost_in(2,ix,ir)+ttn_ghost_in(1,ix,ir)
                     a=v**2*w**2+u**2*w**2+9.0*u**2*v**2
                     b=u**2*(2.0*em*w**2+6.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=ttn(iz,ix,i)
                     tdiv=1.0
                  ELSE
                     u=dnr
                     v=ri*dnx
                     w=risti*dnz
                     em=ttn_ghost_left(iz,2,ir)-ttn(iz,ix,i)
                     en=ttn_ghost_in(2,ix,ir)-ttn(iz,ix,i)
                     a=v**2*w**2+u**2*w**2+u**2*v**2
                     b=-2.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn(iz,ix,i)
                     tdiv=1.0
                  ENDIF
               ENDIF
            ELSE
               IF(swj.EQ.0)THEN
                  swsol=1
                  IF(swk.EQ.0)THEN
                     u=2.0*ri*dnx
                     v=2.0*risti*dnz
                     em=4.0*ttn_ghost_left(iz,2,ir)-ttn_ghost_left(iz,1,ir)-4.0*ttn_ghost_in(2,ix,ir)
                     em=em+ttn_ghost_in(1,ix,ir)
                     a=v**2+u**2
                     b=2.0*em*u**2
                     c=u**2*(em**2-slown**2*v**2)
                     tref=4.0*ttn_ghost_left(iz,2,ir)-ttn_ghost_left(iz,1,ir)
                     tdiv=3.0
                  ELSE
                     u=risti*dnz
                     v=2.0*ri*dnx
                     em=3.0*ttn_ghost_in(2,ix,ir)-4.0*ttn_ghost_left(iz,2,ir)+ttn_ghost_left(iz,1,ir)
                     a=v**2+9.0*u**2
                     b=6.0*em*u**2
                     c=u**2*(em**2-slown**2*v**2)
                     tref=ttn_ghost_in(2,ix,ir)
                     tdiv=1.0
                  ENDIF
               ELSE 
                  swsol=1
                  IF(swk.EQ.0)THEN
                     u=ri*dnx
                     v=2.0*risti*dnz
                     em=3.0*ttn_ghost_left(iz,2,ir)-4.0*ttn_ghost_in(2,ix,ir)+ttn_ghost_in(1,ix,ir)
                     a=v**2+9.0*u**2
                     b=6.0*em*u**2
                     c=u**2*(em**2-v**2*slown**2)
                     tref=ttn_ghost_left(iz,2,ir)
                     tdiv=1.0
                  ELSE
                     u=ri*dnx
                     v=risti*dnz
                     em=ttn_ghost_in(2,ix,ir)-ttn_ghost_left(iz,2,ir)
                     a=u**2+v**2
                     b=-2.0*u**2*em
                     c=u**2*(em**2-v**2*slown**2)
                     tref=ttn_ghost_left(iz,2,ir)
                     tdiv=1.0
                  ENDIF
               ENDIF
            ENDIF
            IF(swsol.EQ.1)THEN
               rd1=b**2-4.0*a*c
               IF(rd1.LT.0.0)rd1=0.0
               tdsh=(-b+sqrt(rd1))/(2.0*a)
               trav=(tref+tdsh)/tdiv
               IF(tsw1.EQ.1)THEN
                  travm=MIN(trav,travm)
               ELSE
                  travm=trav
                  tsw1=1
               ENDIF
            ENDIF      
         ENDIF
      ENDDO
      IF(travm.LT.ttn(iz,ix,ir))THEN
         ttn(iz,ix,ir)=travm
         IF(nsts(iz,ix,ir).EQ.0)THEN
            CALL addtree_all(iz,ix,ir)
         ELSE
            CALL updtree_all(iz,ix,ir)
         ENDIF
      ENDIF
   END SUBROUTINE fouds2_left_in

   SUBROUTINE fouds2_left_out(iz,ix,ir)
      IMPLICIT NONE
      INTEGER :: i,j,k,i2,j2,k2,ir,ix,iz,tsw1
      INTEGER :: swi,swj,swk,swsol
      REAL(KIND=i10) :: old_value
      REAL(KIND=i10) :: trav,travm,slown,tdsh,tref,tdiv
      REAL(KIND=i10) :: a,b,c,u,v,w,em,en,ri,risti,rd1
      tsw1=0
      slown=1.0/veln(iz,ix,ir)
      ri=rgor-(roffset-1)*dnr-(ir-1)*dnr+earth
      risti=ri*sin(rgox+(xoffset-1)*dnx+(ix-1)*dnx)
      DO i=ir-1,ir+1,2
         swi=-1
         IF(i.eq.ir-1)THEN
            i2=i-1
            IF(i2.GE.1.AND.ttn(iz,ix,i2).NE.0)THEN
               swi=0
            ENDIF
         ELSE
            i2=i+1
            IF(i2.LE.rnum.AND.ttn(iz,ix,i2).NE.0)THEN
               swi=0
            ENDIF
         ENDIF
         IF(ttn(iz,ix,i).LT.ttn(iz,ix,ir).AND.swi.EQ.0.AND.ttn(iz,ix,i).NE.0)THEN
            swi=-1
            IF(ttn(iz,ix,i).GT.ttn(iz,ix,i2))THEN
               swi=0
            ENDIF
         ELSE
            swi=-1
         ENDIF
         swj=-1
         IF(ttn_ghost_left(iz,2,ir).GT.ttn_ghost_left(iz,1,ir))THEN
            swj=0
         ENDIF
         swk=-1
         IF(ttn_ghost_out(1,ix,ir).GT.ttn_ghost_out(2,ix,ir))THEN
            swk=0
         ENDIF
         IF(i.GE.1.AND.i.LE.rnum)THEN
            swsol=0
            IF(swi.EQ.0)THEN
               swsol=1
               IF(swj.EQ.0)THEN
                  IF(swk.EQ.0)THEN
                     u=2.0*dnr
                     v=2.0*ri*dnx
                     w=2.0*risti*dnz
                     em=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn_ghost_left(iz,2,ir)
                     em=em+ttn_ghost_left(iz,1,ir)
                     en=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn_ghost_out(1,ix,ir)
                     en=en+ttn_ghost_out(2,ix,ir)
                     a=v**2*w**2+u**2*w**2+u**2*v**2
                     b=2.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)
                     tdiv=3.0
                  ELSE
                     u=risti*dnz
                     v=2.0*ri*dnx
                     w=2.0*dnr
                     em=3.0*ttn_ghost_out(1,ix,ir)-4.0*ttn_ghost_left(iz,2,ir)+ttn_ghost_left(iz,1,ir)
                     en=3.0*ttn_ghost_out(1,ix,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_out(1,ix,ir)
                     tdiv=1.0
                  ENDIF
               ELSE
                  IF(swk.EQ.0)THEN
                     u=ri*dnx
                     v=2.0*dnr
                     w=2.0*risti*dnz
                     em=3.0*ttn_ghost_left(iz,2,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                     en=3.0*ttn_ghost_left(iz,2,ir)-4.0*ttn_ghost_out(1,ix,ir)+ttn_ghost_out(2,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_left(iz,2,ir)
                     tdiv=1.0
                  ELSE
                     u=ri*dnx
                     v=2.0*dnr
                     w=risti*dnz
                     em=3.0*ttn_ghost_left(iz,2,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                     en=ttn_ghost_left(iz,2,ir)-ttn_ghost_out(1,ix,ir)
                     a=w**2*v**2+9.0*u**2*w**2+u**2*v**2
                     b=u**2*(6.0*em*w**2+2.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_left(iz,2,ir)
                     tdiv=1.0
                  ENDIF
               ENDIF
            ELSE IF(ttn(iz,ix,i).LT.ttn(iz,ix,ir).AND.ttn(iz,ix,i).NE.0)THEN
               swsol=1
               IF(swj.EQ.0)THEN
                  IF(swk.EQ.0)THEN
                     u=dnr
                     v=2.0*ri*dnx
                     w=2.0*risti*dnz
                     em=3.0*ttn(iz,ix,i)-4.0*ttn_ghost_left(iz,2,ir)+ttn_ghost_left(iz,1,ir)
                     en=3.0*ttn(iz,ix,i)-4.0*ttn_ghost_out(1,ix,ir)+ttn_ghost_out(2,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn(iz,ix,i)
                     tdiv=1.0
                  ELSE
                     u=dnr
                     v=2.0*ri*dnx
                     w=risti*dnz
                     em=3.0*ttn(iz,ix,i)-4.0*ttn_ghost_left(iz,2,ir)+ttn_ghost_left(iz,1,ir)
                     en=ttn(iz,ix,i)-ttn_ghost_out(1,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+u**2*v**2
                     b=u**2*(6.0*em*w**2+2.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=ttn(iz,ix,i)
                     tdiv=1.0
                  ENDIF
               ELSE 
                  IF(swk.EQ.0)THEN
                     u=dnr
                     v=ri*dnx
                     w=2.0*risti*dnz
                     em=ttn(iz,ix,i)-ttn_ghost_left(iz,2,ir)
                     en=3.0*ttn(iz,ix,i)-4.0*ttn_ghost_out(1,ix,ir)+ttn_ghost_out(2,ix,ir)
                     a=v**2*w**2+u**2*w**2+9.0*u**2*v**2
                     b=u**2*(2.0*em*w**2+6.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=ttn(iz,ix,i)
                     tdiv=1.0
                  ELSE
                     u=dnr
                     v=ri*dnx
                     w=risti*dnz
                     em=ttn_ghost_left(iz,2,ir)-ttn(iz,ix,i)
                     en=ttn_ghost_out(1,ix,ir)-ttn(iz,ix,i)
                     a=v**2*w**2+u**2*w**2+u**2*v**2
                     b=-2.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn(iz,ix,i)
                     tdiv=1.0
                  ENDIF
               ENDIF
            ELSE
               IF(swj.EQ.0)THEN
                  swsol=1
                  IF(swk.EQ.0)THEN
                     u=2.0*ri*dnx
                     v=2.0*risti*dnz
                     em=4.0*ttn_ghost_left(iz,2,ir)-ttn_ghost_left(iz,1,ir)-4.0*ttn_ghost_out(1,ix,ir)
                     em=em+ttn_ghost_out(2,ix,ir)
                     a=v**2+u**2
                     b=2.0*em*u**2
                     c=u**2*(em**2-slown**2*v**2)
                     tref=4.0*ttn_ghost_left(iz,2,ir)-ttn_ghost_left(iz,1,ir)
                     tdiv=3.0
                  ELSE
                     u=risti*dnz
                     v=2.0*ri*dnx
                     em=3.0*ttn_ghost_out(1,ix,ir)-4.0*ttn_ghost_left(iz,2,ir)+ttn_ghost_left(iz,1,ir)
                     a=v**2+9.0*u**2
                     b=6.0*em*u**2
                     c=u**2*(em**2-slown**2*v**2)
                     tref=ttn_ghost_out(1,ix,ir)
                     tdiv=1.0
                  ENDIF
               ELSE
                  swsol=1
                  IF(swk.EQ.0)THEN
                     u=ri*dnx
                     v=2.0*risti*dnz
                     em=3.0*ttn_ghost_left(iz,2,ir)-4.0*ttn_ghost_out(1,ix,ir)+ttn_ghost_out(2,ix,ir)
                     a=v**2+9.0*u**2
                     b=6.0*em*u**2
                     c=u**2*(em**2-v**2*slown**2)
                     tref=ttn_ghost_left(iz,2,ir)
                     tdiv=1.0
                  ELSE
                     u=ri*dnx
                     v=risti*dnz
                     em=ttn_ghost_out(1,ix,ir)-ttn_ghost_left(iz,2,ir)
                     a=u**2+v**2
                     b=-2.0*u**2*em
                     c=u**2*(em**2-v**2*slown**2)
                     tref=ttn_ghost_left(iz,2,ir)
                     tdiv=1.0
                  ENDIF
               ENDIF
            ENDIF
            IF(swsol.EQ.1)THEN
               rd1=b**2-4.0*a*c
               IF(rd1.LT.0.0)rd1=0.0
               tdsh=(-b+sqrt(rd1))/(2.0*a)
               trav=(tref+tdsh)/tdiv
               IF(tsw1.EQ.1)THEN
                  travm=MIN(trav,travm)
               ELSE
                  travm=trav
                  tsw1=1
               ENDIF
            ENDIF      
         ENDIF
      ENDDO
      IF(travm.LT.ttn(iz,ix,ir))THEN
         ttn(iz,ix,ir)=travm
         IF(nsts(iz,ix,ir).EQ.0)THEN
            CALL addtree_all(iz,ix,ir)
         ELSE
            CALL updtree_all(iz,ix,ir)
         ENDIF
      ENDIF
   END SUBROUTINE fouds2_left_out
   
   SUBROUTINE fouds2_left_up(iz,ix,ir)
      IMPLICIT NONE
      INTEGER :: i,j,k,i2,j2,k2,ir,ix,iz,tsw1
      INTEGER :: swi,swj,swk,swsol
      REAL(KIND=i10) :: old_value
      REAL(KIND=i10) :: trav,travm,slown,tdsh,tref,tdiv
      REAL(KIND=i10) :: a,b,c,u,v,w,em,en,ri,risti,rd1
      tsw1=0
      slown=1.0/veln(iz,ix,ir)
      ri=rgor-(roffset-1)*dnr-(ir-1)*dnr+earth
      risti=ri*sin(rgox+(xoffset-1)*dnx+(ix-1)*dnx)
      swi=-1
      IF(ttn_ghost_up(iz,ix,2).GT.ttn_ghost_up(iz,ix,1))THEN
         swi=0
      ENDIF
      swj=-1
      IF(ttn_ghost_left(iz,2,ir).GT.ttn_ghost_left(iz,1,ir))THEN
         swj=0
      ENDIF
      DO k=iz-1,iz+1,2
         swk=-1
         IF(k.eq.iz-1)THEN
            k2=k-1
            IF(k2.GE.1.AND.ttn(k2,ix,ir).NE.0)THEN
               swk=0
            ENDIF
         ELSE
            k2=k+1
            IF(k2.LE.znum.AND.ttn(k2,ix,ir).NE.0)THEN
               swk=0
            ENDIF
         ENDIF
         IF(ttn(k,ix,ir).LT.ttn(iz,ix,ir).AND.swk.EQ.0.AND.ttn(k,ix,ir).NE.0)THEN
            swk=-1
            IF(ttn(k,ix,ir).GT.ttn(k2,ix,ir))THEN
               swk=0
            ENDIF
         ELSE
            swk=-1
         ENDIF
         IF(k.GE.1.AND.k.LE.znum)THEN
            swsol=0
            IF(swi.EQ.0)THEN
               swsol=1
               IF(swj.EQ.0)THEN
                  IF(swk.EQ.0)THEN
                     u=2.0*dnr
                     v=2.0*ri*dnx
                     w=2.0*risti*dnz
                     em=4.0*ttn_ghost_up(iz,ix,2)-ttn_ghost_up(iz,ix,1)-4.0*ttn_ghost_left(iz,2,ir)
                     em=em+ttn_ghost_left(iz,1,ir)
                     en=4.0*ttn_ghost_up(iz,ix,2)-ttn_ghost_up(iz,ix,1)-4.0*ttn(k,ix,ir)
                     en=en+ttn(k2,ix,ir)
                     a=v**2*w**2+u**2*w**2+u**2*v**2
                     b=2.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=4.0*ttn_ghost_up(iz,ix,2)-ttn_ghost_up(iz,ix,1)
                     tdiv=3.0
                  ELSE IF(ttn(k,ix,ir).LT.ttn(iz,ix,ir).AND.ttn(k,ix,ir).NE.0)THEN
                     u=risti*dnz
                     v=2.0*ri*dnx
                     w=2.0*dnr
                     em=3.0*ttn(k,ix,ir)-4.0*ttn_ghost_left(iz,2,ir)+ttn_ghost_left(iz,1,ir)
                     en=3.0*ttn(k,ix,ir)-4.0*ttn_ghost_up(iz,ix,2)+ttn_ghost_up(iz,ix,1)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn(k,ix,ir)
                     tdiv=1.0
                  ELSE
                     u=2.0*dnr
                     v=2.0*ri*dnx
                     em=4.0*ttn_ghost_up(iz,ix,2)-ttn_ghost_up(iz,ix,1)-4.0*ttn_ghost_left(iz,2,ir)
                     em=em+ttn_ghost_left(iz,1,ir)
                     a=v**2+u**2
                     b=2.0*em*u**2
                     c=u**2*(em**2-slown**2*v**2) 
                     tref=4.0*ttn_ghost_up(iz,ix,2)-ttn_ghost_up(iz,ix,1)
                     tdiv=3.0
                  ENDIF
               ELSE
                  IF(swk.EQ.0)THEN
                     u=ri*dnx
                     v=2.0*dnr
                     w=2.0*risti*dnz
                     em=3.0*ttn_ghost_left(iz,2,ir)-4.0*ttn_ghost_up(iz,ix,2)+ttn_ghost_up(iz,ix,1)
                     en=3.0*ttn_ghost_left(iz,2,ir)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_left(iz,2,ir)
                     tdiv=1.0
                  ELSE IF(ttn(k,ix,ir).LT.ttn(iz,ix,ir).AND.ttn(k,ix,ir).NE.0)THEN
                     u=ri*dnx
                     v=2.0*dnr
                     w=risti*dnz
                     em=3.0*ttn_ghost_left(iz,2,ir)-4.0*ttn_ghost_up(iz,ix,2)+ttn_ghost_up(iz,ix,1)
                     en=ttn_ghost_left(iz,2,ir)-ttn(k,ix,ir)
                     a=w**2*v**2+9.0*u**2*w**2+u**2*v**2
                     b=u**2*(6.0*em*w**2+2.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_left(iz,2,ir)
                     tdiv=1.0
                  ELSE
                     u=ri*dnx
                     v=2.0*dnr
                     em=3.0*ttn_ghost_left(iz,2,ir)-4.0*ttn_ghost_up(iz,ix,2)+ttn_ghost_up(iz,ix,1)
                     a=v**2+9.0*u**2
                     b=6.0*em*u**2
                     c=u**2*(em**2-slown**2*v**2)
                     tref=ttn_ghost_left(iz,2,ir)
                     tdiv=1.0
                  ENDIF
               ENDIF
            ELSE
               swsol=1
               IF(swj.EQ.0)THEN
                  IF(swk.EQ.0)THEN
                     u=dnr
                     v=2.0*ri*dnx
                     w=2.0*risti*dnz
                     em=3.0*ttn_ghost_up(iz,ix,2)-4.0*ttn_ghost_left(iz,2,ir)+ttn_ghost_left(iz,1,ir)
                     en=3.0*ttn_ghost_up(iz,ix,2)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_up(iz,ix,2)
                     tdiv=1.0
                  ELSE IF(ttn(k,ix,ir).LT.ttn(iz,ix,ir).AND.ttn(k,ix,ir).NE.0)THEN
                     u=dnr
                     v=2.0*ri*dnx
                     w=risti*dnz
                     em=3.0*ttn_ghost_up(iz,ix,2)-4.0*ttn_ghost_left(iz,2,ir)+ttn_ghost_left(iz,1,ir)
                     en=ttn_ghost_up(iz,ix,2)-ttn(k,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+u**2*v**2
                     b=u**2*(6.0*em*w**2+2.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=ttn_ghost_up(iz,ix,2)
                     tdiv=1.0
                  ELSE
                     u=dnr
                     v=2.0*ri*dnx
                     em=3.0*ttn_ghost_up(iz,ix,2)-4.0*ttn_ghost_left(iz,2,ir)+ttn_ghost_left(iz,1,ir)
                     a=v**2+9.0*u**2
                     b=6.0*em*u**2
                     c=u**2*(em**2-v**2*slown**2)
                     tref=ttn_ghost_up(iz,ix,2)
                     tdiv=1.0
                  ENDIF
               ELSE
                  IF(swk.EQ.0)THEN
                     u=dnr
                     v=ri*dnx
                     w=2.0*risti*dnz
                     em=ttn_ghost_up(iz,ix,2)-ttn_ghost_left(iz,2,ir)
                     en=3.0*ttn_ghost_up(iz,ix,2)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                     a=v**2*w**2+u**2*w**2+9.0*u**2*v**2
                     b=u**2*(2.0*em*w**2+6.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=ttn_ghost_up(iz,ix,2)
                     tdiv=1.0
                  ELSE IF(ttn(k,ix,ir).LT.ttn(iz,ix,ir).AND.ttn(k,ix,ir).NE.0)THEN
                     u=dnr
                     v=ri*dnx
                     w=risti*dnz
                     em=ttn_ghost_left(iz,2,ir)-ttn_ghost_up(iz,ix,2)
                     en=ttn(k,ix,ir)-ttn_ghost_up(iz,ix,2)
                     a=v**2*w**2+u**2*w**2+u**2*v**2
                     b=-2.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_up(iz,ix,2)
                     tdiv=1.0
                  ELSE
                     u=dnr
                     v=ri*dnx
                     em=ttn_ghost_left(iz,2,ir)-ttn_ghost_up(iz,ix,2)
                     a=u**2+v**2
                     b=-2.0*u**2*em
                     c=u**2*(em**2-v**2*slown**2)
                     tref=ttn_ghost_up(iz,ix,2)
                     tdiv=1.0
                  ENDIF
               ENDIF
            ENDIF
            IF(swsol.EQ.1)THEN
               rd1=b**2-4.0*a*c
               IF(rd1.LT.0.0)rd1=0.0
               tdsh=(-b+sqrt(rd1))/(2.0*a)
               trav=(tref+tdsh)/tdiv
               IF(tsw1.EQ.1)THEN
                  travm=MIN(trav,travm)
               ELSE
                  travm=trav
                  tsw1=1
               ENDIF
            ENDIF      
         ENDIF
      ENDDO
      IF(travm.LT.ttn(iz,ix,ir))THEN
         ttn(iz,ix,ir)=travm
         IF(nsts(iz,ix,ir).EQ.0)THEN  
            CALL addtree_all(iz,ix,ir)
         ELSE
            CALL updtree_all(iz,ix,ir)
         ENDIF
      ENDIF
   END SUBROUTINE fouds2_left_up

   SUBROUTINE fouds2_left_down(iz,ix,ir)
      IMPLICIT NONE
      INTEGER :: i,j,k,i2,j2,k2,ir,ix,iz,tsw1
      INTEGER :: swi,swj,swk,swsol
      REAL(KIND=i10) :: old_value
      REAL(KIND=i10) :: trav,travm,slown,tdsh,tref,tdiv
      REAL(KIND=i10) :: a,b,c,u,v,w,em,en,ri,risti,rd1
      tsw1=0
      slown=1.0/veln(iz,ix,ir) 
      ri=rgor-(roffset-1)*dnr-(ir-1)*dnr+earth
      risti=ri*sin(rgox+(xoffset-1)*dnx+(ix-1)*dnx)
      swi=-1
      IF(ttn_ghost_down(iz,ix,1).GT.ttn_ghost_down(iz,ix,2))THEN
         swi=0
      ENDIF
      swj=-1
      IF(ttn_ghost_left(iz,2,ir).GT.ttn_ghost_left(iz,1,ir))THEN
         swj=0
      ENDIF
      DO k=iz-1,iz+1,2
         swk=-1
         IF(k.eq.iz-1)THEN
            k2=k-1
            IF(k2.GE.1.AND.ttn(k2,ix,ir).NE.0)THEN
               swk=0
            ENDIF
         ELSE
            k2=k+1
            IF(k2.LE.znum.AND.ttn(k2,ix,ir).NE.0)THEN
               swk=0
            ENDIF
         ENDIF
         IF(ttn(k,ix,ir).LT.ttn(iz,ix,ir).AND.swk.EQ.0.AND.ttn(k,ix,ir).NE.0)THEN
            swk=-1
            IF(ttn(k,ix,ir).GT.ttn(k2,ix,ir))THEN
               swk=0
            ENDIF
         ELSE
            swk=-1
         ENDIF
         IF(k.GE.1.AND.k.LE.znum)THEN
            swsol=0
            IF(swi.EQ.0)THEN
               swsol=1
               IF(swj.EQ.0)THEN
                  IF(swk.EQ.0)THEN
                     u=2.0*dnr
                     v=2.0*ri*dnx
                     w=2.0*risti*dnz
                     em=4.0*ttn_ghost_down(iz,ix,1)-ttn_ghost_down(iz,ix,2)-4.0*ttn_ghost_left(iz,2,ir)
                     em=em+ttn_ghost_left(iz,1,ir)
                     en=4.0*ttn_ghost_down(iz,ix,1)-ttn_ghost_down(iz,ix,2)-4.0*ttn(k,ix,ir)
                     en=en+ttn(k2,ix,ir)
                     a=v**2*w**2+u**2*w**2+u**2*v**2
                     b=2.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=4.0*ttn_ghost_down(iz,ix,1)-ttn_ghost_down(iz,ix,2)
                     tdiv=3.0
                  ELSE IF(ttn(k,ix,ir).LT.ttn(iz,ix,ir).AND.ttn(k,ix,ir).NE.0)THEN
                     u=risti*dnz
                     v=2.0*ri*dnx
                     w=2.0*dnr
                     em=3.0*ttn(k,ix,ir)-4.0*ttn_ghost_left(iz,2,ir)+ttn_ghost_left(iz,1,ir)
                     en=3.0*ttn(k,ix,ir)-4.0*ttn_ghost_down(iz,ix,1)+ttn_ghost_down(iz,ix,2)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn(k,ix,ir)
                     tdiv=1.0
                  ELSE
                     u=2.0*dnr
                     v=2.0*ri*dnx
                     em=4.0*ttn_ghost_down(iz,ix,1)-ttn_ghost_down(iz,ix,2)-4.0*ttn_ghost_left(iz,2,ir)
                     em=em+ttn_ghost_left(iz,1,ir)
                     a=v**2+u**2
                     b=2.0*em*u**2
                     c=u**2*(em**2-slown**2*v**2) 
                     tref=4.0*ttn_ghost_down(iz,ix,1)-ttn_ghost_down(iz,ix,2)
                     tdiv=3.0
                  ENDIF
               ELSE
                  IF(swk.EQ.0)THEN
                     u=ri*dnx
                     v=2.0*dnr
                     w=2.0*risti*dnz
                     em=3.0*ttn_ghost_left(iz,2,ir)-4.0*ttn_ghost_down(iz,ix,1)+ttn_ghost_down(iz,ix,2)
                     en=3.0*ttn_ghost_left(iz,2,ir)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_left(iz,2,ir)
                     tdiv=1.0
                  ELSE IF(ttn(k,ix,ir).LT.ttn(iz,ix,ir).AND.ttn(k,ix,ir).NE.0)THEN
                     u=ri*dnx
                     v=2.0*dnr
                     w=risti*dnz
                     em=3.0*ttn_ghost_left(iz,2,ir)-4.0*ttn_ghost_down(iz,ix,1)+ttn_ghost_down(iz,ix,2)
                     en=ttn_ghost_left(iz,2,ir)-ttn(k,ix,ir)
                     a=w**2*v**2+9.0*u**2*w**2+u**2*v**2
                     b=u**2*(6.0*em*w**2+2.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_left(iz,2,ir)
                     tdiv=1.0
                  ELSE
                     u=ri*dnx
                     v=2.0*dnr
                     em=3.0*ttn_ghost_left(iz,2,ir)-4.0*ttn_ghost_down(iz,ix,1)+ttn_ghost_down(iz,ix,2)
                     a=v**2+9.0*u**2
                     b=6.0*em*u**2
                     c=u**2*(em**2-slown**2*v**2)
                     tref=ttn_ghost_left(iz,2,ir)
                     tdiv=1.0
                  ENDIF
               ENDIF
            ELSE
               swsol=1
               IF(swj.EQ.0)THEN
                  IF(swk.EQ.0)THEN
                     u=dnr
                     v=2.0*ri*dnx
                     w=2.0*risti*dnz
                     em=3.0*ttn_ghost_down(iz,ix,1)-4.0*ttn_ghost_left(iz,2,ir)+ttn_ghost_left(iz,1,ir)
                     en=3.0*ttn_ghost_down(iz,ix,1)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_down(iz,ix,1)
                     tdiv=1.0
                  ELSE IF(ttn(k,ix,ir).LT.ttn(iz,ix,ir).AND.ttn(k,ix,ir).NE.0)THEN
                     u=dnr
                     v=2.0*ri*dnx
                     w=risti*dnz
                     em=3.0*ttn_ghost_down(iz,ix,1)-4.0*ttn_ghost_left(iz,2,ir)+ttn_ghost_left(iz,1,ir)
                     en=ttn_ghost_down(iz,ix,1)-ttn(k,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+u**2*v**2
                     b=u**2*(6.0*em*w**2+2.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=ttn_ghost_down(iz,ix,1)
                     tdiv=1.0
                  ELSE
                     u=dnr
                     v=2.0*ri*dnx
                     em=3.0*ttn_ghost_down(iz,ix,1)-4.0*ttn_ghost_left(iz,2,ir)+ttn_ghost_left(iz,1,ir)
                     a=v**2+9.0*u**2
                     b=6.0*em*u**2
                     c=u**2*(em**2-v**2*slown**2)
                     tref=ttn_ghost_down(iz,ix,1)
                     tdiv=1.0
                  ENDIF
               ELSE
                  IF(swk.EQ.0)THEN
                     u=dnr
                     v=ri*dnx
                     w=2.0*risti*dnz
                     em=ttn_ghost_down(iz,ix,1)-ttn_ghost_left(iz,2,ir)
                     en=3.0*ttn_ghost_down(iz,ix,1)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                     a=v**2*w**2+u**2*w**2+9.0*u**2*v**2
                     b=u**2*(2.0*em*w**2+6.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=ttn_ghost_down(iz,ix,1)
                     tdiv=1.0
                  ELSE IF(ttn(k,ix,ir).LT.ttn(iz,ix,ir).AND.ttn(k,ix,ir).NE.0)THEN
                     u=dnr
                     v=ri*dnx
                     w=risti*dnz
                     em=ttn_ghost_left(iz,2,ir)-ttn_ghost_down(iz,ix,1)
                     en=ttn(k,ix,ir)-ttn_ghost_down(iz,ix,1)
                     a=v**2*w**2+u**2*w**2+u**2*v**2
                     b=-2.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_down(iz,ix,1)
                     tdiv=1.0
                  ELSE
                     u=dnr
                     v=ri*dnx
                     em=ttn_ghost_left(iz,2,ir)-ttn_ghost_down(iz,ix,1)
                     a=u**2+v**2
                     b=-2.0*u**2*em
                     c=u**2*(em**2-v**2*slown**2)
                     tref=ttn_ghost_down(iz,ix,1)
                     tdiv=1.0
                  ENDIF
               ENDIF
            ENDIF
            IF(swsol.EQ.1)THEN
               rd1=b**2-4.0*a*c
               IF(rd1.LT.0.0)rd1=0.0
               tdsh=(-b+sqrt(rd1))/(2.0*a)
               trav=(tref+tdsh)/tdiv
               IF(tsw1.EQ.1)THEN
                  travm=MIN(trav,travm)
               ELSE
                  travm=trav
                  tsw1=1
               ENDIF
            ENDIF      
         ENDIF
      ENDDO
      IF(travm.LT.ttn(iz,ix,ir))THEN
         ttn(iz,ix,ir)=travm
         IF(nsts(iz,ix,ir).EQ.0)THEN 
            CALL addtree_all(iz,ix,ir)
         ELSE 
            CALL updtree_all(iz,ix,ir)
         ENDIF
      ENDIF
   END SUBROUTINE fouds2_left_down

   SUBROUTINE fouds2_right_in(iz,ix,ir)
      IMPLICIT NONE
      INTEGER :: i,j,k,i2,j2,k2,ir,ix,iz,tsw1
      INTEGER :: swi,swj,swk,swsol
      REAL(KIND=i10) :: old_value
      REAL(KIND=i10) :: trav,travm,slown,tdsh,tref,tdiv
      REAL(KIND=i10) :: a,b,c,u,v,w,em,en,ri,risti,rd1      
      tsw1=0
      slown=1.0/veln(iz,ix,ir)      
      ri=rgor-(roffset-1)*dnr-(ir-1)*dnr+earth
      risti=ri*sin(rgox+(xoffset-1)*dnx+(ix-1)*dnx)
      DO i=ir-1,ir+1,2
         swi=-1
         IF(i.eq.ir-1)THEN
            i2=i-1
            IF(i2.GE.1.AND.ttn(iz,ix,i2).NE.0)THEN
               swi=0
            ENDIF
         ELSE
            i2=i+1
            IF(i2.LE.rnum.AND.ttn(iz,ix,i2).NE.0)THEN
               swi=0
            ENDIF
         ENDIF
         IF(ttn(iz,ix,i).LT.ttn(iz,ix,ir).AND.swi.EQ.0.AND.ttn(iz,ix,i).NE.0)THEN
            swi=-1
            IF(ttn(iz,ix,i).GT.ttn(iz,ix,i2))THEN
               swi=0
            ENDIF
         ELSE
            swi=-1
         ENDIF         
         swj=-1
         IF(ttn_ghost_right(iz,1,ir).GT.ttn_ghost_right(iz,2,ir))THEN
            swj=0
         ENDIF
         swk=-1
         IF(ttn_ghost_in(2,ix,ir).GT.ttn_ghost_in(1,ix,ir))THEN
            swk=0
         ENDIF
         IF(i.GE.1.AND.i.LE.rnum)THEN
            swsol=0
            IF(swi.EQ.0)THEN
               swsol=1
               IF(swj.EQ.0)THEN 
                  IF(swk.EQ.0)THEN
                     u=2.0*dnr
                     v=2.0*ri*dnx
                     w=2.0*risti*dnz
                     em=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn_ghost_right(iz,1,ir)
                     em=em+ttn_ghost_right(iz,2,ir)
                     en=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn_ghost_in(2,ix,ir)
                     en=en+ttn_ghost_in(1,ix,ir)
                     a=v**2*w**2+u**2*w**2+u**2*v**2
                     b=2.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)
                     tdiv=3.0
                  ELSE
                     u=risti*dnz
                     v=2.0*ri*dnx
                     w=2.0*dnr
                     em=3.0*ttn_ghost_in(2,ix,ir)-4.0*ttn_ghost_right(iz,1,ir)+ttn_ghost_right(iz,2,ir)
                     en=3.0*ttn_ghost_in(2,ix,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_in(2,ix,ir)
                     tdiv=1.0
                  ENDIF
               ELSE 
                  IF(swk.EQ.0)THEN
                     u=ri*dnx
                     v=2.0*dnr
                     w=2.0*risti*dnz
                     em=3.0*ttn_ghost_right(iz,1,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                     en=3.0*ttn_ghost_right(iz,1,ir)-4.0*ttn_ghost_in(2,ix,ir)+ttn_ghost_in(1,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_right(iz,1,ir)
                     tdiv=1.0
                  ELSE
                     u=ri*dnx
                     v=2.0*dnr
                     w=risti*dnz
                     em=3.0*ttn_ghost_right(iz,1,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                     en=ttn_ghost_right(iz,1,ir)-ttn_ghost_in(2,ix,ir)
                     a=w**2*v**2+9.0*u**2*w**2+u**2*v**2
                     b=u**2*(6.0*em*w**2+2.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_right(iz,1,ir)
                     tdiv=1.0
                  ENDIF               
               ENDIF
            ELSE IF(ttn(iz,ix,i).LT.ttn(iz,ix,ir).AND.ttn(iz,ix,i).NE.0)THEN
               swsol=1
               IF(swj.EQ.0)THEN
                  IF(swk.EQ.0)THEN
                     u=dnr
                     v=2.0*ri*dnx
                     w=2.0*risti*dnz
                     em=3.0*ttn(iz,ix,i)-4.0*ttn_ghost_right(iz,1,ir)+ttn_ghost_right(iz,2,ir)
                     en=3.0*ttn(iz,ix,i)-4.0*ttn_ghost_in(2,ix,ir)+ttn_ghost_in(1,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn(iz,ix,i)
                     tdiv=1.0
                  ELSE
                     u=dnr
                     v=2.0*ri*dnx
                     w=risti*dnz
                     em=3.0*ttn(iz,ix,i)-4.0*ttn_ghost_right(iz,1,ir)+ttn_ghost_right(iz,2,ir)
                     en=ttn(iz,ix,i)-ttn_ghost_in(2,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+u**2*v**2
                     b=u**2*(6.0*em*w**2+2.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=ttn(iz,ix,i)
                     tdiv=1.0   
                  ENDIF
               ELSE
                  IF(swk.EQ.0)THEN
                     u=dnr
                     v=ri*dnx
                     w=2.0*risti*dnz
                     em=ttn(iz,ix,i)-ttn_ghost_right(iz,1,ir)
                     en=3.0*ttn(iz,ix,i)-4.0*ttn_ghost_in(2,ix,ir)+ttn_ghost_in(1,ix,ir)
                     a=v**2*w**2+u**2*w**2+9.0*u**2*v**2
                     b=u**2*(2.0*em*w**2+6.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=ttn(iz,ix,i)
                     tdiv=1.0
                  ELSE
                     u=dnr
                     v=ri*dnx
                     w=risti*dnz
                     em=ttn_ghost_right(iz,1,ir)-ttn(iz,ix,i)
                     en=ttn_ghost_in(2,ix,ir)-ttn(iz,ix,i)
                     a=v**2*w**2+u**2*w**2+u**2*v**2
                     b=-2.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn(iz,ix,i)
                     tdiv=1.0
                  ENDIF
               ENDIF
            ELSE
               IF(swj.EQ.0)THEN
                  swsol=1
                  IF(swk.EQ.0)THEN
                     u=2.0*ri*dnx
                     v=2.0*risti*dnz
                     em=4.0*ttn_ghost_right(iz,1,ir)-ttn_ghost_right(iz,2,ir)-4.0*ttn_ghost_in(2,ix,ir)
                     em=em+ttn_ghost_in(1,ix,ir)
                     a=v**2+u**2
                     b=2.0*em*u**2
                     c=u**2*(em**2-slown**2*v**2)
                     tref=4.0*ttn_ghost_right(iz,1,ir)-ttn_ghost_right(iz,2,ir)
                     tdiv=3.0
                  ELSE
                     u=risti*dnz
                     v=2.0*ri*dnx
                     em=3.0*ttn_ghost_in(2,ix,ir)-4.0*ttn_ghost_right(iz,1,ir)+ttn_ghost_right(iz,2,ir)
                     a=v**2+9.0*u**2
                     b=6.0*em*u**2
                     c=u**2*(em**2-slown**2*v**2)
                     tref=ttn_ghost_in(2,ix,ir)
                     tdiv=1.0
                  ENDIF
               ELSE
                  swsol=1
                  IF(swk.EQ.0)THEN
                     u=ri*dnx
                     v=2.0*risti*dnz
                     em=3.0*ttn_ghost_right(iz,1,ir)-4.0*ttn_ghost_in(2,ix,ir)+ttn_ghost_in(1,ix,ir)
                     a=v**2+9.0*u**2
                     b=6.0*em*u**2
                     c=u**2*(em**2-v**2*slown**2)
                     tref=ttn_ghost_right(iz,1,ir)
                     tdiv=1.0
                  ELSE
                     u=ri*dnx
                     v=risti*dnz
                     em=ttn_ghost_in(2,ix,ir)-ttn_ghost_right(iz,1,ir)
                     a=u**2+v**2
                     b=-2.0*u**2*em
                     c=u**2*(em**2-v**2*slown**2)
                     tref=ttn_ghost_right(iz,2,ir)
                     tdiv=1.0                 
                  ENDIF
               ENDIF
            ENDIF
            IF(swsol.EQ.1)THEN
               rd1=b**2-4.0*a*c
               IF(rd1.LT.0.0)rd1=0.0
               tdsh=(-b+sqrt(rd1))/(2.0*a)
               trav=(tref+tdsh)/tdiv
               IF(tsw1.EQ.1)THEN
                  travm=MIN(trav,travm)
               ELSE
                  travm=trav
                  tsw1=1
               ENDIF
            ENDIF      
         ENDIF         
      ENDDO     
      IF(travm.LT.ttn(iz,ix,ir))THEN
         ttn(iz,ix,ir)=travm
         IF(nsts(iz,ix,ir).EQ.0)THEN
            CALL addtree_all(iz,ix,ir)
         ELSE 
            CALL updtree_all(iz,ix,ir)
         ENDIF
      ENDIF
   END SUBROUTINE fouds2_right_in

   SUBROUTINE fouds2_right_out(iz,ix,ir)
      IMPLICIT NONE
      INTEGER :: i,j,k,i2,j2,k2,ir,ix,iz,tsw1
      INTEGER :: swi,swj,swk,swsol
      REAL(KIND=i10) :: old_value
      REAL(KIND=i10) :: trav,travm,slown,tdsh,tref,tdiv
      REAL(KIND=i10) :: a,b,c,u,v,w,em,en,ri,risti,rd1
      tsw1=0
      slown=1.0/veln(iz,ix,ir)      
      ri=rgor-(roffset-1)*dnr-(ir-1)*dnr+earth
      risti=ri*sin(rgox+(xoffset-1)*dnx+(ix-1)*dnx)
      DO i=ir-1,ir+1,2
         swi=-1
         IF(i.eq.ir-1)THEN
            i2=i-1
            IF(i2.GE.1.AND.ttn(iz,ix,i2).NE.0)THEN
               swi=0
            ENDIF
         ELSE
            i2=i+1
            IF(i2.LE.rnum.AND.ttn(iz,ix,i2).NE.0)THEN
               swi=0
            ENDIF
         ENDIF
         IF(ttn(iz,ix,i).LT.ttn(iz,ix,ir).AND.swi.EQ.0.AND.ttn(iz,ix,i).NE.0)THEN
            swi=-1
            IF(ttn(iz,ix,i).GT.ttn(iz,ix,i2))THEN
               swi=0
            ENDIF
         ELSE
            swi=-1
         ENDIF         
         swj=-1
         IF(ttn_ghost_right(iz,1,ir).GT.ttn_ghost_right(iz,2,ir))THEN
            swj=0
         ENDIF
         swk=-1
         IF(ttn_ghost_out(1,ix,ir).GT.ttn_ghost_out(2,ix,ir))THEN
            swk=0
         ENDIF
         IF(i.GE.1.AND.i.LE.rnum)THEN
            swsol=0
            IF(swi.EQ.0)THEN
               swsol=1
               IF(swj.EQ.0)THEN
                  IF(swk.EQ.0)THEN
                     u=2.0*dnr
                     v=2.0*ri*dnx
                     w=2.0*risti*dnz
                     em=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn_ghost_right(iz,1,ir)
                     em=em+ttn_ghost_right(iz,2,ir)
                     en=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn_ghost_out(1,ix,ir)
                     en=en+ttn_ghost_out(2,ix,ir)
                     a=v**2*w**2+u**2*w**2+u**2*v**2
                     b=2.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)
                     tdiv=3.0
                  ELSE
                     u=risti*dnz
                     v=2.0*ri*dnx
                     w=2.0*dnr
                     em=3.0*ttn_ghost_out(1,ix,ir)-4.0*ttn_ghost_right(iz,1,ir)+ttn_ghost_right(iz,2,ir)
                     en=3.0*ttn_ghost_out(1,ix,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_out(1,ix,ir)
                     tdiv=1.0
                  ENDIF
               ELSE
                  IF(swk.EQ.0)THEN
                     u=ri*dnx
                     v=2.0*dnr
                     w=2.0*risti*dnz
                     em=3.0*ttn_ghost_right(iz,1,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                     en=3.0*ttn_ghost_right(iz,1,ir)-4.0*ttn_ghost_out(1,ix,ir)+ttn_ghost_out(2,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_right(iz,1,ir)
                     tdiv=1.0
                  ELSE
                     u=ri*dnx
                     v=2.0*dnr
                     w=risti*dnz
                     em=3.0*ttn_ghost_right(iz,1,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                     en=ttn_ghost_right(iz,1,ir)-ttn_ghost_out(1,ix,ir)
                     a=w**2*v**2+9.0*u**2*w**2+u**2*v**2
                     b=u**2*(6.0*em*w**2+2.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_right(iz,1,ir)
                     tdiv=1.0
                  ENDIF               
               ENDIF
            ELSE IF(ttn(iz,ix,i).LT.ttn(iz,ix,ir).AND.ttn(iz,ix,i).NE.0)THEN
               swsol=1
               IF(swj.EQ.0)THEN
                  IF(swk.EQ.0)THEN
                     u=dnr
                     v=2.0*ri*dnx
                     w=2.0*risti*dnz
                     em=3.0*ttn(iz,ix,i)-4.0*ttn_ghost_right(iz,1,ir)+ttn_ghost_right(iz,2,ir)
                     en=3.0*ttn(iz,ix,i)-4.0*ttn_ghost_out(1,ix,ir)+ttn_ghost_out(2,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn(iz,ix,i)
                     tdiv=1.0
                  ELSE
                     u=dnr
                     v=2.0*ri*dnx
                     w=risti*dnz
                     em=3.0*ttn(iz,ix,i)-4.0*ttn_ghost_right(iz,1,ir)+ttn_ghost_right(iz,2,ir)
                     en=ttn(iz,ix,i)-ttn_ghost_out(1,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+u**2*v**2
                     b=u**2*(6.0*em*w**2+2.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=ttn(iz,ix,i)
                     tdiv=1.0   
                  ENDIF
               ELSE 
                  IF(swk.EQ.0)THEN
                     u=dnr
                     v=ri*dnx
                     w=2.0*risti*dnz
                     em=ttn(iz,ix,i)-ttn_ghost_right(iz,1,ir)
                     en=3.0*ttn(iz,ix,i)-4.0*ttn_ghost_out(1,ix,ir)+ttn_ghost_out(2,ix,ir)
                     a=v**2*w**2+u**2*w**2+9.0*u**2*v**2
                     b=u**2*(2.0*em*w**2+6.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=ttn(iz,ix,i)
                     tdiv=1.0
                  ELSE
                     u=dnr
                     v=ri*dnx
                     w=risti*dnz
                     em=ttn_ghost_right(iz,1,ir)-ttn(iz,ix,i)
                     en=ttn_ghost_out(1,ix,ir)-ttn(iz,ix,i)
                     a=v**2*w**2+u**2*w**2+u**2*v**2
                     b=-2.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn(iz,ix,i)
                     tdiv=1.0
                  ENDIF
               ENDIF
            ELSE
               IF(swj.EQ.0)THEN
                  swsol=1
                  IF(swk.EQ.0)THEN
                     u=2.0*ri*dnx
                     v=2.0*risti*dnz
                     em=4.0*ttn_ghost_right(iz,1,ir)-ttn_ghost_right(iz,2,ir)-4.0*ttn_ghost_out(1,ix,ir)
                     em=em+ttn_ghost_out(2,ix,ir)
                     a=v**2+u**2
                     b=2.0*em*u**2
                     c=u**2*(em**2-slown**2*v**2)
                     tref=4.0*ttn_ghost_right(iz,1,ir)-ttn_ghost_right(iz,2,ir)
                     tdiv=3.0
                  ELSE
                     u=risti*dnz
                     v=2.0*ri*dnx
                     em=3.0*ttn_ghost_out(1,ix,ir)-4.0*ttn_ghost_right(iz,1,ir)+ttn_ghost_right(iz,2,ir)
                     a=v**2+9.0*u**2
                     b=6.0*em*u**2
                     c=u**2*(em**2-slown**2*v**2)
                     tref=ttn_ghost_out(1,ix,ir)
                     tdiv=1.0   
                  ENDIF
               ELSE 
                  swsol=1
                  IF(swk.EQ.0)THEN
                     u=ri*dnx
                     v=2.0*risti*dnz
                     em=3.0*ttn_ghost_right(iz,1,ir)-4.0*ttn_ghost_out(1,ix,ir)+ttn_ghost_out(2,ix,ir)
                     a=v**2+9.0*u**2
                     b=6.0*em*u**2
                     c=u**2*(em**2-v**2*slown**2)
                     tref=ttn_ghost_right(iz,1,ir)
                     tdiv=1.0
                  ELSE
                     u=ri*dnx
                     v=risti*dnz
                     em=ttn_ghost_out(1,ix,ir)-ttn_ghost_right(iz,1,ir)
                     a=u**2+v**2
                     b=-2.0*u**2*em
                     c=u**2*(em**2-v**2*slown**2)
                     tref=ttn_ghost_right(iz,2,ir)
                     tdiv=1.0
                  ENDIF
               ENDIF
            ENDIF

            IF(swsol.EQ.1)THEN
               rd1=b**2-4.0*a*c
               IF(rd1.LT.0.0)rd1=0.0
               tdsh=(-b+sqrt(rd1))/(2.0*a)
               trav=(tref+tdsh)/tdiv
               IF(tsw1.EQ.1)THEN
                  travm=MIN(trav,travm)
               ELSE
                  travm=trav
                  tsw1=1
               ENDIF
            ENDIF      
         ENDIF         
      ENDDO
      IF(travm.LT.ttn(iz,ix,ir))THEN
         ttn(iz,ix,ir)=travm
         IF(nsts(iz,ix,ir).EQ.0)THEN
            CALL addtree_all(iz,ix,ir)
         ELSE
            CALL updtree_all(iz,ix,ir)
         ENDIF
      ENDIF
   END SUBROUTINE fouds2_right_out

   SUBROUTINE fouds2_right_up(iz,ix,ir)
      IMPLICIT NONE
      INTEGER :: i,j,k,i2,j2,k2,ir,ix,iz,tsw1
      INTEGER :: swi,swj,swk,swsol
      REAL(KIND=i10) :: old_value
      REAL(KIND=i10) :: trav,travm,slown,tdsh,tref,tdiv
      REAL(KIND=i10) :: a,b,c,u,v,w,em,en,ri,risti,rd1
      tsw1=0
      slown=1.0/veln(iz,ix,ir)      
      ri=rgor-(roffset-1)*dnr-(ir-1)*dnr+earth
      risti=ri*sin(rgox+(xoffset-1)*dnx+(ix-1)*dnx)
      swi=-1
      IF(ttn_ghost_up(iz,ix,2).GT.ttn_ghost_up(iz,ix,1))THEN
         swi=0
      ENDIF         
      swj=-1
      IF(ttn_ghost_right(iz,1,ir).GT.ttn_ghost_right(iz,2,ir))THEN
         swj=0
      ENDIF
      DO k=iz-1,iz+1,2
         swk=-1
         IF(k.eq.iz-1)THEN
            k2=k-1
            IF(k2.GE.1.AND.ttn(k2,ix,ir).NE.0)THEN
               swk=0
            ENDIF
         ELSE
            k2=k+1
            IF(k2.LE.znum.AND.ttn(k2,ix,ir).NE.0)THEN
               swk=0
            ENDIF
         ENDIF
         IF(ttn(k,ix,ir).LT.ttn(iz,ix,ir).AND.swk.EQ.0.AND.ttn(k,ix,ir).NE.0)THEN
            swk=-1
            IF(ttn(k,ix,ir).GT.ttn(k2,ix,ir))THEN
               swk=0
            ENDIF
         ELSE
            swk=-1
         ENDIF
         IF(k.GE.1.AND.k.LE.znum)THEN
            swsol=0
            IF(swi.EQ.0)THEN
               swsol=1
               IF(swj.EQ.0)THEN
                  IF(swk.EQ.0)THEN
                     u=2.0*dnr
                     v=2.0*ri*dnx
                     w=2.0*risti*dnz
                     em=4.0*ttn_ghost_up(iz,ix,2)-ttn_ghost_up(iz,ix,1)-4.0*ttn_ghost_right(iz,1,ir)
                     em=em+ttn_ghost_right(iz,2,ir)
                     en=4.0*ttn_ghost_up(iz,ix,2)-ttn_ghost_up(iz,ix,1)-4.0*ttn(k,ix,ir)
                     en=en+ttn(k2,ix,ir)
                     a=v**2*w**2+u**2*w**2+u**2*v**2
                     b=2.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=4.0*ttn_ghost_up(iz,ix,2)-ttn_ghost_up(iz,ix,1)
                     tdiv=3.0
                  ELSE IF(ttn(k,ix,ir).LT.ttn(iz,ix,ir).AND.ttn(k,ix,ir).NE.0)THEN
                     u=risti*dnz
                     v=2.0*ri*dnx
                     w=2.0*dnr
                     em=3.0*ttn(k,ix,ir)-4.0*ttn_ghost_right(iz,1,ir)+ttn_ghost_right(iz,2,ir)
                     en=3.0*ttn(k,ix,ir)-4.0*ttn_ghost_up(iz,ix,2)+ttn_ghost_up(iz,ix,1)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn(k,ix,ir)
                     tdiv=1.0
                  ELSE
                     u=2.0*dnr
                     v=2.0*ri*dnx
                     em=4.0*ttn_ghost_up(iz,ix,2)-ttn_ghost_up(iz,ix,1)-4.0*ttn_ghost_right(iz,1,ir)
                     em=em+ttn_ghost_right(iz,2,ir)
                     a=v**2+u**2
                     b=2.0*em*u**2
                     c=u**2*(em**2-slown**2*v**2) 
                     tref=4.0*ttn_ghost_up(iz,ix,2)-ttn_ghost_up(iz,ix,1)
                     tdiv=3.0
                  ENDIF
               ELSE
                  IF(swk.EQ.0)THEN
                     u=ri*dnx
                     v=2.0*dnr
                     w=2.0*risti*dnz
                     em=3.0*ttn_ghost_right(iz,1,ir)-4.0*ttn_ghost_up(iz,ix,2)+ttn_ghost_up(iz,ix,1)
                     en=3.0*ttn_ghost_right(iz,1,ir)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_right(iz,1,ir)
                     tdiv=1.0
                  ELSE IF(ttn(k,ix,ir).LT.ttn(iz,ix,ir).AND.ttn(k,ix,ir).NE.0)THEN
                     u=ri*dnx
                     v=2.0*dnr
                     w=risti*dnz
                     em=3.0*ttn_ghost_right(iz,1,ir)-4.0*ttn_ghost_up(iz,ix,2)+ttn_ghost_up(iz,ix,1)
                     en=ttn_ghost_right(iz,1,ir)-ttn(k,ix,ir)
                     a=w**2*v**2+9.0*u**2*w**2+u**2*v**2
                     b=u**2*(6.0*em*w**2+2.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_right(iz,1,ir)
                     tdiv=1.0
                  ELSE
                     u=ri*dnx
                     v=2.0*dnr
                     em=3.0*ttn_ghost_right(iz,1,ir)-4.0*ttn_ghost_up(iz,ix,2)+ttn_ghost_up(iz,ix,1)
                     a=v**2+9.0*u**2
                     b=6.0*em*u**2
                     c=u**2*(em**2-slown**2*v**2)
                     tref=ttn_ghost_right(iz,1,ir)
                     tdiv=1.0
                  ENDIF
               ENDIF
            ELSE
               swsol=1
               IF(swj.EQ.0)THEN
                  IF(swk.EQ.0)THEN
                     u=dnr
                     v=2.0*ri*dnx
                     w=2.0*risti*dnz
                     em=3.0*ttn_ghost_up(iz,ix,2)-4.0*ttn_ghost_right(iz,1,ir)+ttn_ghost_right(iz,2,ir)
                     en=3.0*ttn_ghost_up(iz,ix,2)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_up(iz,ix,2)
                     tdiv=1.0
                  ELSE IF(ttn(k,ix,ir).LT.ttn(iz,ix,ir).AND.ttn(k,ix,ir).NE.0)THEN
                     u=dnr
                     v=2.0*ri*dnx
                     w=risti*dnz
                     em=3.0*ttn_ghost_up(iz,ix,2)-4.0*ttn_ghost_right(iz,1,ir)+ttn_ghost_right(iz,2,ir)
                     en=ttn_ghost_up(iz,ix,2)-ttn(k,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+u**2*v**2
                     b=u**2*(6.0*em*w**2+2.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=ttn_ghost_up(iz,ix,2)
                     tdiv=1.0
                  ELSE
                     u=dnr
                     v=2.0*ri*dnx
                     em=3.0*ttn_ghost_up(iz,ix,2)-4.0*ttn_ghost_right(iz,1,ir)+ttn_ghost_right(iz,2,ir)
                     a=v**2+9.0*u**2
                     b=6.0*em*u**2
                     c=u**2*(em**2-v**2*slown**2)
                     tref=ttn_ghost_up(iz,ix,2)
                     tdiv=1.0
                  ENDIF
               ELSE
                  IF(swk.EQ.0)THEN
                     u=dnr
                     v=ri*dnx
                     w=2.0*risti*dnz
                     em=ttn_ghost_up(iz,ix,2)-ttn_ghost_right(iz,1,ir)
                     en=3.0*ttn_ghost_up(iz,ix,2)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                     a=v**2*w**2+u**2*w**2+9.0*u**2*v**2
                     b=u**2*(2.0*em*w**2+6.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=ttn_ghost_up(iz,ix,2)
                     tdiv=1.0
                  ELSE IF(ttn(k,ix,ir).LT.ttn(iz,ix,ir).AND.ttn(k,ix,ir).NE.0)THEN
                     u=dnr
                     v=ri*dnx
                     w=risti*dnz
                     em=ttn_ghost_right(iz,1,ir)-ttn_ghost_up(iz,ix,2)
                     en=ttn(k,ix,ir)-ttn_ghost_up(iz,ix,2)
                     a=v**2*w**2+u**2*w**2+u**2*v**2
                     b=-2.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_up(iz,ix,2)
                     tdiv=1.0
                  ELSE
                     u=dnr
                     v=ri*dnx
                     em=ttn_ghost_right(iz,1,ir)-ttn_ghost_up(iz,ix,2)
                     a=u**2+v**2
                     b=-2.0*u**2*em
                     c=u**2*(em**2-v**2*slown**2)
                     tref=ttn_ghost_up(iz,ix,2)
                     tdiv=1.0
                  ENDIF
               ENDIF          
            ENDIF
            IF(swsol.EQ.1)THEN
               rd1=b**2-4.0*a*c
               IF(rd1.LT.0.0)rd1=0.0
               tdsh=(-b+sqrt(rd1))/(2.0*a)
               trav=(tref+tdsh)/tdiv
               IF(tsw1.EQ.1)THEN
                  travm=MIN(trav,travm)
               ELSE
                  travm=trav
                  tsw1=1
               ENDIF
            ENDIF      
         ENDIF
      ENDDO         
      IF(travm.LT.ttn(iz,ix,ir))THEN
         ttn(iz,ix,ir)=travm
         IF(nsts(iz,ix,ir).EQ.0)THEN
            CALL addtree_all(iz,ix,ir)
         ELSE
            CALL updtree_all(iz,ix,ir)
         ENDIF
      ENDIF
   END SUBROUTINE fouds2_right_up

   SUBROUTINE fouds2_right_down(iz,ix,ir)
      IMPLICIT NONE
      INTEGER :: i,j,k,i2,j2,k2,ir,ix,iz,tsw1
      INTEGER :: swi,swj,swk,swsol
      REAL(KIND=i10) :: old_value
      REAL(KIND=i10) :: trav,travm,slown,tdsh,tref,tdiv
      REAL(KIND=i10) :: a,b,c,u,v,w,em,en,ri,risti,rd1
      tsw1=0
      slown=1.0/veln(iz,ix,ir)      
      ri=rgor-(roffset-1)*dnr-(ir-1)*dnr+earth
      risti=ri*sin(rgox+(xoffset-1)*dnx+(ix-1)*dnx)
      swi=-1
      IF(ttn_ghost_down(iz,ix,1).GT.ttn_ghost_down(iz,ix,2))THEN
         swi=0
      ENDIF         
      swj=-1
      IF(ttn_ghost_right(iz,1,ir).GT.ttn_ghost_right(iz,2,ir))THEN
         swj=0
      ENDIF
      DO k=iz-1,iz+1,2
         swk=-1
         IF(k.eq.iz-1)THEN
            k2=k-1
            IF(k2.GE.1.AND.ttn(k2,ix,ir).NE.0)THEN
               swk=0
            ENDIF
         ELSE
            k2=k+1
            IF(k2.LE.znum.AND.ttn(k2,ix,ir).NE.0)THEN
               swk=0
            ENDIF
         ENDIF
         IF(ttn(k,ix,ir).LT.ttn(iz,ix,ir).AND.swk.EQ.0.AND.ttn(k,ix,ir).NE.0)THEN
            swk=-1
            IF(ttn(k,ix,ir).GT.ttn(k2,ix,ir))THEN
               swk=0
            ENDIF
         ELSE
            swk=-1
         ENDIF
         IF(k.GE.1.AND.k.LE.znum)THEN
            swsol=0
            IF(swi.EQ.0)THEN
               swsol=1
               IF(swj.EQ.0)THEN
                  IF(swk.EQ.0)THEN
                     u=2.0*dnr
                     v=2.0*ri*dnx
                     w=2.0*risti*dnz
                     em=4.0*ttn_ghost_down(iz,ix,1)-ttn_ghost_down(iz,ix,2)-4.0*ttn_ghost_right(iz,1,ir)
                     em=em+ttn_ghost_right(iz,2,ir)
                     en=4.0*ttn_ghost_down(iz,ix,1)-ttn_ghost_down(iz,ix,2)-4.0*ttn(k,ix,ir)
                     en=en+ttn(k2,ix,ir)
                     a=v**2*w**2+u**2*w**2+u**2*v**2
                     b=2.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=4.0*ttn_ghost_down(iz,ix,1)-ttn_ghost_down(iz,ix,2)
                     tdiv=3.0
                  ELSE IF(ttn(k,ix,ir).LT.ttn(iz,ix,ir).AND.ttn(k,ix,ir).NE.0)THEN
                     u=risti*dnz
                     v=2.0*ri*dnx
                     w=2.0*dnr
                     em=3.0*ttn(k,ix,ir)-4.0*ttn_ghost_right(iz,1,ir)+ttn_ghost_right(iz,2,ir)
                     en=3.0*ttn(k,ix,ir)-4.0*ttn_ghost_down(iz,ix,1)+ttn_ghost_down(iz,ix,2)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn(k,ix,ir)
                     tdiv=1.0
                  ELSE
                     u=2.0*dnr
                     v=2.0*ri*dnx
                     em=4.0*ttn_ghost_down(iz,ix,1)-ttn_ghost_down(iz,ix,2)-4.0*ttn_ghost_right(iz,1,ir)
                     em=em+ttn_ghost_right(iz,2,ir)
                     a=v**2+u**2
                     b=2.0*em*u**2
                     c=u**2*(em**2-slown**2*v**2) 
                     tref=4.0*ttn_ghost_down(iz,ix,1)-ttn_ghost_down(iz,ix,2)
                     tdiv=3.0
                  ENDIF
               ELSE
                  IF(swk.EQ.0)THEN
                     u=ri*dnx
                     v=2.0*dnr
                     w=2.0*risti*dnz
                     em=3.0*ttn_ghost_right(iz,1,ir)-4.0*ttn_ghost_down(iz,ix,1)+ttn_ghost_down(iz,ix,2)
                     en=3.0*ttn_ghost_right(iz,1,ir)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_right(iz,1,ir)
                     tdiv=1.0
                  ELSE IF(ttn(k,ix,ir).LT.ttn(iz,ix,ir).AND.ttn(k,ix,ir).NE.0)THEN
                     u=ri*dnx
                     v=2.0*dnr
                     w=risti*dnz
                     em=3.0*ttn_ghost_right(iz,1,ir)-4.0*ttn_ghost_down(iz,ix,1)+ttn_ghost_down(iz,ix,2)
                     en=ttn_ghost_right(iz,1,ir)-ttn(k,ix,ir)
                     a=w**2*v**2+9.0*u**2*w**2+u**2*v**2
                     b=u**2*(6.0*em*w**2+2.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_right(iz,1,ir)
                     tdiv=1.0
                  ELSE
                     u=ri*dnx
                     v=2.0*dnr
                     em=3.0*ttn_ghost_right(iz,1,ir)-4.0*ttn_ghost_down(iz,ix,1)+ttn_ghost_down(iz,ix,2)
                     a=v**2+9.0*u**2
                     b=6.0*em*u**2
                     c=u**2*(em**2-slown**2*v**2)
                     tref=ttn_ghost_right(iz,1,ir)
                     tdiv=1.0
                  ENDIF               
               ENDIF
            ELSE
               swsol=1
               IF(swj.EQ.0)THEN
                  IF(swk.EQ.0)THEN
                     u=dnr
                     v=2.0*ri*dnx
                     w=2.0*risti*dnz
                     em=3.0*ttn_ghost_down(iz,ix,1)-4.0*ttn_ghost_right(iz,1,ir)+ttn_ghost_right(iz,2,ir)
                     en=3.0*ttn_ghost_down(iz,ix,1)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_down(iz,ix,1)
                     tdiv=1.0
                  ELSE IF(ttn(k,ix,ir).LT.ttn(iz,ix,ir).AND.ttn(k,ix,ir).NE.0)THEN
                     u=dnr
                     v=2.0*ri*dnx
                     w=risti*dnz
                     em=3.0*ttn_ghost_down(iz,ix,1)-4.0*ttn_ghost_right(iz,1,ir)+ttn_ghost_right(iz,2,ir)
                     en=ttn_ghost_down(iz,ix,1)-ttn(k,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+u**2*v**2
                     b=u**2*(6.0*em*w**2+2.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=ttn_ghost_down(iz,ix,1)
                     tdiv=1.0
                  ELSE
                     u=dnr
                     v=2.0*ri*dnx
                     em=3.0*ttn_ghost_down(iz,ix,1)-4.0*ttn_ghost_right(iz,1,ir)+ttn_ghost_right(iz,2,ir)
                     a=v**2+9.0*u**2
                     b=6.0*em*u**2
                     c=u**2*(em**2-v**2*slown**2)
                     tref=ttn_ghost_down(iz,ix,1)
                     tdiv=1.0
                  ENDIF
               ELSE
                  IF(swk.EQ.0)THEN
                     u=dnr
                     v=ri*dnx
                     w=2.0*risti*dnz
                     em=ttn_ghost_down(iz,ix,1)-ttn_ghost_right(iz,1,ir)
                     en=3.0*ttn_ghost_down(iz,ix,1)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                     a=v**2*w**2+u**2*w**2+9.0*u**2*v**2
                     b=u**2*(2.0*em*w**2+6.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=ttn_ghost_down(iz,ix,1)
                     tdiv=1.0
                  ELSE IF(ttn(k,ix,ir).LT.ttn(iz,ix,ir).AND.ttn(k,ix,ir).NE.0)THEN
                     u=dnr
                     v=ri*dnx
                     w=risti*dnz
                     em=ttn_ghost_right(iz,1,ir)-ttn_ghost_down(iz,ix,1)
                     en=ttn(k,ix,ir)-ttn_ghost_down(iz,ix,1)
                     a=v**2*w**2+u**2*w**2+u**2*v**2
                     b=-2.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_down(iz,ix,1)
                     tdiv=1.0
                  ELSE
                     u=dnr
                     v=ri*dnx
                     em=ttn_ghost_right(iz,1,ir)-ttn_ghost_down(iz,ix,1)
                     a=u**2+v**2
                     b=-2.0*u**2*em
                     c=u**2*(em**2-v**2*slown**2)
                     tref=ttn_ghost_down(iz,ix,1)
                     tdiv=1.0
                  ENDIF
               ENDIF          
            ENDIF
            IF(swsol.EQ.1)THEN
               rd1=b**2-4.0*a*c
               IF(rd1.LT.0.0)rd1=0.0
               tdsh=(-b+sqrt(rd1))/(2.0*a)
               trav=(tref+tdsh)/tdiv
               IF(tsw1.EQ.1)THEN
                  travm=MIN(trav,travm)
               ELSE
                  travm=trav
                  tsw1=1
               ENDIF
            ENDIF      
         ENDIF
      ENDDO               
      IF(travm.LT.ttn(iz,ix,ir))THEN
         ttn(iz,ix,ir)=travm
         IF(nsts(iz,ix,ir).EQ.0)THEN
            CALL addtree_all(iz,ix,ir)
         ELSE  
            CALL updtree_all(iz,ix,ir)
         ENDIF
      ENDIF
   END SUBROUTINE fouds2_right_down

   SUBROUTINE fouds2_in_up(iz,ix,ir)
      IMPLICIT NONE
      INTEGER :: i,j,k,i2,j2,k2,ir,ix,iz,tsw1
      INTEGER :: swi,swj,swk,swsol
      INTEGER :: rtmp,xtmp,ztmp
      REAL(KIND=i10) :: old_value
      REAL(KIND=i10) :: trav,travm,slown,tdsh,tref,tdiv
      REAL(KIND=i10) :: a,b,c,u,v,w,em,en,ri,risti,rd1
      tsw1=0
      slown=1.0/veln(iz,ix,ir)
      ri=rgor-(roffset-1)*dnr-(ir-1)*dnr+earth
      risti=ri*sin(rgox+(xoffset-1)*dnx+(ix-1)*dnx)      
      swi=-1
      IF(ttn_ghost_up(iz,ix,2).GT.ttn_ghost_up(iz,ix,1))THEN
         swi=0
      ENDIF
      DO j=ix-1,ix+1,2
         swj=-1
         IF(j.eq.ix-1)THEN
            j2=j-1
            IF(j2.GE.1.AND.ttn(iz,j2,ir).NE.0)THEN
              swj=0
            ENDIF
         ELSE
            j2=j+1
            IF(j2.LE.xnum.AND.ttn(iz,j2,ir).NE.0)THEN
               swj=0
            ENDIF
         ENDIF
         IF(ttn(iz,j,ir).LT.ttn(iz,ix,ir).AND.swj.EQ.0.AND.ttn(iz,j,ir).NE.0)THEN
            swj=-1
            IF(ttn(iz,j,ir).GT.ttn(iz,j2,ir))THEN
               swj=0
            ENDIF
         ELSE
            swj=-1
         ENDIF         
         swk=-1
         IF(ttn_ghost_in(2,ix,ir).GT.ttn_ghost_in(1,ix,ir))THEN
            swk=0
         ENDIF         
         IF(j.GE.1.AND.j.LE.xnum)THEN
            swsol=0
            IF(swi.EQ.0)THEN
               swsol=1
               IF(swj.EQ.0)THEN
                  IF(swk.EQ.0)THEN
                     u=2.0*dnr
                     v=2.0*ri*dnx
                     w=2.0*risti*dnz
                     em=4.0*ttn_ghost_up(iz,ix,2)-ttn_ghost_up(iz,ix,1)-4.0*ttn(iz,j,ir)
                     em=em+ttn(iz,j2,ir)
                     en=4.0*ttn_ghost_up(iz,ix,2)-ttn_ghost_up(iz,ix,1)-4.0*ttn_ghost_in(2,ix,ir)
                     en=en+ttn_ghost_in(1,ix,ir)
                     a=v**2*w**2+u**2*w**2+u**2*v**2
                     b=2.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=4.0*ttn_ghost_up(iz,ix,2)-ttn_ghost_up(iz,ix,1)
                     tdiv=3.0
                  ELSE
                     u=risti*dnz
                     v=2.0*ri*dnx
                     w=2.0*dnr
                     em=3.0*ttn_ghost_in(2,ix,ir)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                     en=3.0*ttn_ghost_in(2,ix,ir)-4.0*ttn_ghost_up(iz,ix,2)+ttn_ghost_up(iz,ix,1)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_in(2,ix,ir)
                     tdiv=1.0
                  ENDIF
               ELSEIF(ttn(iz,j,ir).LT.ttn(iz,ix,ir).AND.ttn(iz,j,ir).NE.0)THEN
                  IF(swk.EQ.0)THEN
                     u=ri*dnx
                     v=2.0*dnr
                     w=2.0*risti*dnz
                     em=3.0*ttn(iz,j,ir)-4.0*ttn_ghost_up(iz,ix,2)+ttn_ghost_up(iz,ix,1)
                     en=3.0*ttn(iz,j,ir)-4.0*ttn_ghost_in(2,ix,ir)+ttn_ghost_in(1,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn(iz,j,ir)
                     tdiv=1.0
                  ELSE
                     u=ri*dnx
                     v=2.0*dnr
                     w=risti*dnz
                     em=3.0*ttn(iz,j,ir)-4.0*ttn_ghost_up(iz,ix,2)+ttn_ghost_up(iz,ix,1)
                     en=ttn(iz,j,ir)-ttn_ghost_in(2,ix,ir)
                     a=w**2*v**2+9.0*u**2*w**2+u**2*v**2
                     b=u**2*(6.0*em*w**2+2.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn(iz,j,ir)
                     tdiv=1.0
                  ENDIF
               ELSE
                  IF(swk.EQ.0)THEN
                     u=2.0*dnr
                     v=2.0*risti*dnz
                     em=4.0*ttn_ghost_up(iz,ix,2)-ttn_ghost_up(iz,ix,1)-4.0*ttn_ghost_in(2,ix,ir)
                     em=em+ttn_ghost_in(1,ix,ir)
                     a=u**2+v**2
                     b=2.0*em*u**2
                     c=u**2*(em**2-v**2*slown**2)
                     tref=4.0*ttn_ghost_up(iz,ix,2)-ttn_ghost_up(iz,ix,1)
                     tdiv=3.0
                  ELSE
                     u=risti*dnz
                     v=2.0*dnr
                     em=3.0*ttn_ghost_in(2,ix,ir)-4.0*ttn_ghost_up(iz,ix,2)+ttn_ghost_up(iz,ix,1)
                     a=v**2+9.0*u**2
                     b=6.0*em*u**2
                     c=u**2*(em**2-v**2*slown**2)
                     tref=ttn_ghost_in(2,ix,ir)
                     tdiv=1.0
                  ENDIF
               ENDIF
            ELSE
               swsol=1
               IF(swj.EQ.0)THEN
                  IF(swk.EQ.0)THEN
                     u=dnr
                     v=2.0*ri*dnx
                     w=2.0*risti*dnz
                     em=3.0*ttn_ghost_up(iz,ix,2)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                     en=3.0*ttn_ghost_up(iz,ix,2)-4.0*ttn_ghost_in(2,ix,ir)+ttn_ghost_in(1,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_up(iz,ix,2)
                     tdiv=1.0
                  ELSE
                     u=dnr
                     v=2.0*ri*dnx
                     w=risti*dnz
                     em=3.0*ttn_ghost_up(iz,ix,2)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                     en=ttn_ghost_up(iz,ix,2)-ttn_ghost_in(2,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+u**2*v**2
                     b=u**2*(6.0*em*w**2+2.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=ttn_ghost_up(iz,ix,2)
                     tdiv=1.0
                  ENDIF
               ELSE IF(ttn(iz,j,ir).LT.ttn(iz,ix,ir).AND.ttn(iz,j,ir).NE.0)THEN
                  IF(swk.EQ.0)THEN
                     u=dnr
                     v=ri*dnx
                     w=2.0*risti*dnz
                     em=ttn_ghost_up(iz,ix,2)-ttn(iz,j,ir)
                     en=3.0*ttn_ghost_up(iz,ix,2)-4.0*ttn_ghost_in(2,ix,ir)+ttn_ghost_in(1,ix,ir)
                     a=v**2*w**2+u**2*w**2+9.0*u**2*v**2
                     b=u**2*(2.0*em*w**2+6.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=ttn_ghost_up(iz,ix,2)
                     tdiv=1.0
                  ELSE
                     u=dnr
                     v=ri*dnx
                     w=risti*dnz
                     em=ttn(iz,j,ir)-ttn_ghost_up(iz,ix,2)
                     en=ttn_ghost_in(2,ix,ir)-ttn_ghost_up(iz,ix,2)
                     a=v**2*w**2+u**2*w**2+u**2*v**2
                     b=-2.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_up(iz,ix,2)
                     tdiv=1.0
                  ENDIF
               ELSE
                  IF(swk.EQ.0)THEN
                     u=dnr
                     v=2.0*risti*dnz
                     em=3.0*ttn_ghost_up(iz,ix,2)-4.0*ttn_ghost_in(2,ix,ir)+ttn_ghost_in(1,ix,ir)
                     a=v**2+9.0*u**2
                     b=6.0*em*u**2
                     c=u**2*(em**2-v**2*slown**2)
                     tref=ttn_ghost_up(iz,ix,2)
                     tdiv=1.0
                  ELSE 
                     u=dnr
                     v=risti*dnz
                     em=ttn_ghost_in(2,ix,ir)-ttn_ghost_up(iz,ix,2)
                     a=u**2+v**2
                     b=-2.0*u**2*em
                     c=u**2*(em**2-v**2*slown**2)
                     tref=ttn_ghost_up(iz,ix,2)
                     tdiv=1.0
                  ENDIF
               ENDIF            
            ENDIF
            IF(swsol.EQ.1)THEN
               rd1=b**2-4.0*a*c
               IF(rd1.LT.0.0)rd1=0.0
               tdsh=(-b+sqrt(rd1))/(2.0*a)
               trav=(tref+tdsh)/tdiv
               IF(tsw1.EQ.1)THEN
                  travm=MIN(trav,travm)
               ELSE
                  travm=trav
                  tsw1=1
               ENDIF
            ENDIF
         ENDIF         
      ENDDO
      IF(travm.LT.ttn(iz,ix,ir))THEN
         ttn(iz,ix,ir)=travm
         IF(nsts(iz,ix,ir).EQ.0)THEN 
            CALL addtree_all(iz,ix,ir)
         ELSE
            CALL updtree_all(iz,ix,ir)
         ENDIF
      ENDIF
   END SUBROUTINE fouds2_in_up

   SUBROUTINE fouds2_in_down(iz,ix,ir)
      IMPLICIT NONE
      INTEGER :: i,j,k,i2,j2,k2,ir,ix,iz,tsw1
      INTEGER :: swi,swj,swk,swsol
      INTEGER :: rtmp,xtmp,ztmp
      REAL(KIND=i10) :: old_value
      REAL(KIND=i10) :: trav,travm,slown,tdsh,tref,tdiv
      REAL(KIND=i10) :: a,b,c,u,v,w,em,en,ri,risti,rd1
      tsw1=0
      slown=1.0/veln(iz,ix,ir)
      ri=rgor-(roffset-1)*dnr-(ir-1)*dnr+earth
      risti=ri*sin(rgox+(xoffset-1)*dnx+(ix-1)*dnx)      
      swi=-1
      IF(ttn_ghost_down(iz,ix,1).GT.ttn_ghost_down(iz,ix,2))THEN
         swi=0
      ENDIF
      DO j=ix-1,ix+1,2
         swj=-1
         IF(j.eq.ix-1)THEN
            j2=j-1
            IF(j2.GE.1.AND.ttn(iz,j2,ir).NE.0)THEN
              swj=0
            ENDIF
         ELSE
            j2=j+1
            IF(j2.LE.xnum.AND.ttn(iz,j2,ir).NE.0)THEN
               swj=0
            ENDIF
         ENDIF
         IF(ttn(iz,j,ir).LT.ttn(iz,ix,ir).AND.swj.EQ.0.AND.ttn(iz,j,ir).NE.0)THEN
            swj=-1
            IF(ttn(iz,j,ir).GT.ttn(iz,j2,ir))THEN
               swj=0
            ENDIF
         ELSE
            swj=-1
         ENDIF         
         swk=-1
         IF(ttn_ghost_in(2,ix,ir).GT.ttn_ghost_in(1,ix,ir))THEN
            swk=0
         ENDIF         
         IF(j.GE.1.AND.j.LE.xnum)THEN
            swsol=0
            IF(swi.EQ.0)THEN
               swsol=1
               IF(swj.EQ.0)THEN
                  IF(swk.EQ.0)THEN
                     u=2.0*dnr
                     v=2.0*ri*dnx
                     w=2.0*risti*dnz
                     em=4.0*ttn_ghost_down(iz,ix,2)-ttn_ghost_down(iz,ix,2)-4.0*ttn(iz,j,ir)
                     em=em+ttn(iz,j2,ir)
                     en=4.0*ttn_ghost_down(iz,ix,1)-ttn_ghost_down(iz,ix,2)-4.0*ttn_ghost_in(2,ix,ir)
                     en=en+ttn_ghost_in(1,ix,ir)
                     a=v**2*w**2+u**2*w**2+u**2*v**2
                     b=2.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=4.0*ttn_ghost_down(iz,ix,1)-ttn_ghost_down(iz,ix,2)
                     tdiv=3.0
                  ELSE
                     u=risti*dnz
                     v=2.0*ri*dnx
                     w=2.0*dnr
                     em=3.0*ttn_ghost_in(2,ix,ir)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                     en=3.0*ttn_ghost_in(2,ix,ir)-4.0*ttn_ghost_down(iz,ix,1)+ttn_ghost_down(iz,ix,2)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_in(2,ix,ir)
                     tdiv=1.0
                  ENDIF
               ELSEIF(ttn(iz,j,ir).LT.ttn(iz,ix,ir).AND.ttn(iz,j,ir).NE.0)THEN
                  IF(swk.EQ.0)THEN
                     u=ri*dnx
                     v=2.0*dnr
                     w=2.0*risti*dnz
                     em=3.0*ttn(iz,j,ir)-4.0*ttn_ghost_down(iz,ix,1)+ttn_ghost_down(iz,ix,2)
                     en=3.0*ttn(iz,j,ir)-4.0*ttn_ghost_in(2,ix,ir)+ttn_ghost_in(1,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn(iz,j,ir)
                     tdiv=1.0
                  ELSE
                     u=ri*dnx
                     v=2.0*dnr
                     w=risti*dnz
                     em=3.0*ttn(iz,j,ir)-4.0*ttn_ghost_down(iz,ix,1)+ttn_ghost_down(iz,ix,2)
                     en=ttn(iz,j,ir)-ttn_ghost_in(2,ix,ir)
                     a=w**2*v**2+9.0*u**2*w**2+u**2*v**2
                     b=u**2*(6.0*em*w**2+2.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn(iz,j,ir)
                     tdiv=1.0
                  ENDIF
               ELSE
                  IF(swk.EQ.0)THEN
                     u=2.0*dnr
                     v=2.0*risti*dnz
                     em=4.0*ttn_ghost_down(iz,ix,1)-ttn_ghost_down(iz,ix,2)-4.0*ttn_ghost_in(2,ix,ir)
                     em=em+ttn_ghost_in(1,ix,ir)
                     a=u**2+v**2
                     b=2.0*em*u**2
                     c=u**2*(em**2-v**2*slown**2)
                     tref=4.0*ttn_ghost_down(iz,ix,1)-ttn_ghost_down(iz,ix,2)
                     tdiv=3.0
                  ELSE
                     u=risti*dnz
                     v=2.0*dnr
                     em=3.0*ttn_ghost_in(2,ix,ir)-4.0*ttn_ghost_down(iz,ix,1)+ttn_ghost_down(iz,ix,2)
                     a=v**2+9.0*u**2
                     b=6.0*em*u**2
                     c=u**2*(em**2-v**2*slown**2)
                     tref=ttn_ghost_in(2,ix,ir)
                     tdiv=1.0
                  ENDIF
               ENDIF
            ELSE
               swsol=1
               IF(swj.EQ.0)THEN
                  IF(swk.EQ.0)THEN
                     u=dnr
                     v=2.0*ri*dnx
                     w=2.0*risti*dnz
                     em=3.0*ttn_ghost_down(iz,ix,1)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                     en=3.0*ttn_ghost_down(iz,ix,1)-4.0*ttn_ghost_in(2,ix,ir)+ttn_ghost_in(1,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_down(iz,ix,1)
                     tdiv=1.0
                  ELSE
                     u=dnr
                     v=2.0*ri*dnx
                     w=risti*dnz
                     em=3.0*ttn_ghost_down(iz,ix,1)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                     en=ttn_ghost_down(iz,ix,1)-ttn_ghost_in(2,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+u**2*v**2
                     b=u**2*(6.0*em*w**2+2.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=ttn_ghost_down(iz,ix,1)
                     tdiv=1.0
                  ENDIF
               ELSE IF(ttn(iz,j,ir).LT.ttn(iz,ix,ir).AND.ttn(iz,j,ir).NE.0)THEN
                  IF(swk.EQ.0)THEN
                     u=dnr
                     v=ri*dnx
                     w=2.0*risti*dnz
                     em=ttn_ghost_down(iz,ix,1)-ttn(iz,j,ir)
                     en=3.0*ttn_ghost_down(iz,ix,1)-4.0*ttn_ghost_in(2,ix,ir)+ttn_ghost_in(1,ix,ir)
                     a=v**2*w**2+u**2*w**2+9.0*u**2*v**2
                     b=u**2*(2.0*em*w**2+6.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=ttn_ghost_down(iz,ix,1)
                     tdiv=1.0
                  ELSE
                     u=dnr
                     v=ri*dnx
                     w=risti*dnz
                     em=ttn(iz,j,ir)-ttn_ghost_down(iz,ix,1)
                     en=ttn_ghost_in(2,ix,ir)-ttn_ghost_down(iz,ix,1)
                     a=v**2*w**2+u**2*w**2+u**2*v**2
                     b=-2.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_down(iz,ix,1)
                     tdiv=1.0
                  ENDIF
               ELSE
                  IF(swk.EQ.0)THEN
                     u=dnr
                     v=2.0*risti*dnz
                     em=3.0*ttn_ghost_down(iz,ix,1)-4.0*ttn_ghost_in(2,ix,ir)+ttn_ghost_in(1,ix,ir)
                     a=v**2+9.0*u**2
                     b=6.0*em*u**2
                     c=u**2*(em**2-v**2*slown**2)
                     tref=ttn_ghost_down(iz,ix,1)
                     tdiv=1.0
                  ELSE 
                     u=dnr
                     v=risti*dnz
                     em=ttn_ghost_in(2,ix,ir)-ttn_ghost_down(iz,ix,1)
                     a=u**2+v**2
                     b=-2.0*u**2*em
                     c=u**2*(em**2-v**2*slown**2)
                     tref=ttn_ghost_down(iz,ix,1)
                     tdiv=1.0
                  ENDIF
               ENDIF            
            ENDIF
            IF(swsol.EQ.1)THEN
               rd1=b**2-4.0*a*c
               IF(rd1.LT.0.0)rd1=0.0
               tdsh=(-b+sqrt(rd1))/(2.0*a)
               trav=(tref+tdsh)/tdiv
               IF(tsw1.EQ.1)THEN
                  travm=MIN(trav,travm)
               ELSE
                  travm=trav
                  tsw1=1
               ENDIF
            ENDIF
         ENDIF         
      ENDDO      
      IF(travm.LT.ttn(iz,ix,ir))THEN
         ttn(iz,ix,ir)=travm
         IF(nsts(iz,ix,ir).EQ.0)THEN
            CALL addtree_all(iz,ix,ir)
         ELSE
            CALL updtree_all(iz,ix,ir)
         ENDIF
      ENDIF
   END SUBROUTINE fouds2_in_down

   SUBROUTINE fouds2_out_up(iz,ix,ir)
      IMPLICIT NONE
      INTEGER :: i,j,k,i2,j2,k2,ir,ix,iz,tsw1
      INTEGER :: swi,swj,swk,swsol
      INTEGER :: rtmp,xtmp,ztmp
      REAL(KIND=i10) :: old_value
      REAL(KIND=i10) :: trav,travm,slown,tdsh,tref,tdiv
      REAL(KIND=i10) :: a,b,c,u,v,w,em,en,ri,risti,rd1
      tsw1=0
      slown=1.0/veln(iz,ix,ir)
      ri=rgor-(roffset-1)*dnr-(ir-1)*dnr+earth
      risti=ri*sin(rgox+(xoffset-1)*dnx+(ix-1)*dnx)      
      swi=-1
      IF(ttn_ghost_up(iz,ix,2).GT.ttn_ghost_up(iz,ix,1))THEN
         swi=0
      ENDIF
      DO j=ix-1,ix+1,2
         swj=-1
         IF(j.eq.ix-1)THEN
            j2=j-1
            IF(j2.GE.1.AND.ttn(iz,j2,ir).NE.0)THEN
               swj=0
            ENDIF
         ELSE
            j2=j+1
            IF(j2.LE.xnum.AND.ttn(iz,j2,ir).NE.0)THEN
              swj=0
            ENDIF
         ENDIF
         IF(ttn(iz,j,ir).LT.ttn(iz,ix,ir).AND.swj.EQ.0.AND.ttn(iz,j,ir).NE.0)THEN
            swj=-1
            IF(ttn(iz,j,ir).GT.ttn(iz,j2,ir))THEN
               swj=0
            ENDIF
         ELSE
            swj=-1
         ENDIF         
         swk=-1
         IF(ttn_ghost_out(1,ix,ir).GT.ttn_ghost_out(2,ix,ir))THEN
            swk=0
         ENDIF         
         IF(j.GE.1.AND.j.LE.xnum)THEN
            swsol=0
            IF(swi.EQ.0)THEN
               swsol=1
               IF(swj.EQ.0)THEN
                  IF(swk.EQ.0)THEN
                     u=2.0*dnr
                     v=2.0*ri*dnx
                     w=2.0*risti*dnz
                     em=4.0*ttn_ghost_up(iz,ix,2)-ttn_ghost_up(iz,ix,1)-4.0*ttn(iz,j,ir)
                     em=em+ttn(iz,j2,ir)
                     en=4.0*ttn_ghost_up(iz,ix,2)-ttn_ghost_up(iz,ix,1)-4.0*ttn_ghost_out(1,ix,ir)
                     en=en+ttn_ghost_out(2,ix,ir)
                     a=v**2*w**2+u**2*w**2+u**2*v**2
                     b=2.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=4.0*ttn_ghost_up(iz,ix,2)-ttn_ghost_up(iz,ix,1)
                     tdiv=3.0
                  ELSE
                     u=risti*dnz
                     v=2.0*ri*dnx
                     w=2.0*dnr
                     em=3.0*ttn_ghost_out(1,ix,ir)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                     en=3.0*ttn_ghost_out(1,ix,ir)-4.0*ttn_ghost_up(iz,ix,2)+ttn_ghost_up(iz,ix,1)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_out(1,ix,ir)
                     tdiv=1.0
                  ENDIF
               ELSEIF(ttn(iz,j,ir).LT.ttn(iz,ix,ir).AND.ttn(iz,j,ir).NE.0)THEN
                  IF(swk.EQ.0)THEN
                     u=ri*dnx
                     v=2.0*dnr
                     w=2.0*risti*dnz
                     em=3.0*ttn(iz,j,ir)-4.0*ttn_ghost_up(iz,ix,2)+ttn_ghost_up(iz,ix,1)
                     en=3.0*ttn(iz,j,ir)-4.0*ttn_ghost_out(1,ix,ir)+ttn_ghost_out(2,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn(iz,j,ir)
                     tdiv=1.0
                  ELSE
                     u=ri*dnx
                     v=2.0*dnr
                     w=risti*dnz
                     em=3.0*ttn(iz,j,ir)-4.0*ttn_ghost_up(iz,ix,2)+ttn_ghost_up(iz,ix,1)
                     en=ttn(iz,j,ir)-ttn_ghost_out(1,ix,ir)
                     a=w**2*v**2+9.0*u**2*w**2+u**2*v**2
                     b=u**2*(6.0*em*w**2+2.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn(iz,j,ir)
                     tdiv=1.0
                  ENDIF
               ELSE
                  IF(swk.EQ.0)THEN
                     u=2.0*dnr
                     v=2.0*risti*dnz
                     em=4.0*ttn_ghost_up(iz,ix,2)-ttn_ghost_up(iz,ix,1)-4.0*ttn_ghost_out(1,ix,ir)
                     em=em+ttn_ghost_out(2,ix,ir)
                     a=u**2+v**2
                     b=2.0*em*u**2
                     c=u**2*(em**2-v**2*slown**2)
                     tref=4.0*ttn_ghost_up(iz,ix,2)-ttn_ghost_up(iz,ix,1)
                     tdiv=3.0
                  ELSE
                     u=risti*dnz
                     v=2.0*dnr
                     em=3.0*ttn_ghost_out(1,ix,ir)-4.0*ttn_ghost_up(iz,ix,2)+ttn_ghost_up(iz,ix,1)
                     a=v**2+9.0*u**2
                     b=6.0*em*u**2
                     c=u**2*(em**2-v**2*slown**2)
                     tref=ttn_ghost_out(1,ix,ir)
                     tdiv=1.0
                  ENDIF
               ENDIF
            ELSE
               swsol=1
               IF(swj.EQ.0)THEN
                  IF(swk.EQ.0)THEN
                     u=dnr
                     v=2.0*ri*dnx
                     w=2.0*risti*dnz
                     em=3.0*ttn_ghost_up(iz,ix,2)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                     en=3.0*ttn_ghost_up(iz,ix,2)-4.0*ttn_ghost_out(1,ix,ir)+ttn_ghost_out(2,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_up(iz,ix,2)
                     tdiv=1.0
                  ELSE
                     u=dnr
                     v=2.0*ri*dnx
                     w=risti*dnz
                     em=3.0*ttn_ghost_up(iz,ix,2)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                     en=ttn_ghost_up(iz,ix,2)-ttn_ghost_out(1,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+u**2*v**2
                     b=u**2*(6.0*em*w**2+2.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=ttn_ghost_up(iz,ix,2)
                     tdiv=1.0
                  ENDIF
               ELSE IF(ttn(iz,j,ir).LT.ttn(iz,ix,ir).AND.ttn(iz,j,ir).NE.0)THEN
                  IF(swk.EQ.0)THEN
                     u=dnr
                     v=ri*dnx
                     w=2.0*risti*dnz
                     em=ttn_ghost_up(iz,ix,2)-ttn(iz,j,ir)
                     en=3.0*ttn_ghost_up(iz,ix,2)-4.0*ttn_ghost_out(1,ix,ir)+ttn_ghost_out(2,ix,ir)
                     a=v**2*w**2+u**2*w**2+9.0*u**2*v**2
                     b=u**2*(2.0*em*w**2+6.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=ttn_ghost_up(iz,ix,2)
                     tdiv=1.0
                  ELSE
                     u=dnr
                     v=ri*dnx
                     w=risti*dnz
                     em=ttn(iz,j,ir)-ttn_ghost_up(iz,ix,2)
                     en=ttn_ghost_out(1,ix,ir)-ttn_ghost_up(iz,ix,2)
                     a=v**2*w**2+u**2*w**2+u**2*v**2
                     b=-2.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_up(iz,ix,2)
                     tdiv=1.0
                  ENDIF
               ELSE
                  IF(swk.EQ.0)THEN
                     u=dnr
                     v=2.0*risti*dnz
                     em=3.0*ttn_ghost_up(iz,ix,2)-4.0*ttn_ghost_out(1,ix,ir)+ttn_ghost_out(2,ix,ir)
                     a=v**2+9.0*u**2
                     b=6.0*em*u**2
                     c=u**2*(em**2-v**2*slown**2)
                     tref=ttn_ghost_up(iz,ix,2)
                     tdiv=1.0
                  ELSE 
                     u=dnr
                     v=risti*dnz
                     em=ttn_ghost_out(1,ix,ir)-ttn_ghost_up(iz,ix,2)
                     a=u**2+v**2
                     b=-2.0*u**2*em
                     c=u**2*(em**2-v**2*slown**2)
                     tref=ttn_ghost_up(iz,ix,2)
                     tdiv=1.0
                  ENDIF
               ENDIF            
            ENDIF
            IF(swsol.EQ.1)THEN
               rd1=b**2-4.0*a*c
               IF(rd1.LT.0.0)rd1=0.0
               tdsh=(-b+sqrt(rd1))/(2.0*a)
               trav=(tref+tdsh)/tdiv
               IF(tsw1.EQ.1)THEN
                  travm=MIN(trav,travm)
               ELSE
                  travm=trav
                  tsw1=1
               ENDIF
            ENDIF
         ENDIF         
      ENDDO      
      IF(travm.LT.ttn(iz,ix,ir))THEN
         ttn(iz,ix,ir)=travm
         IF(nsts(iz,ix,ir).EQ.0)THEN
            CALL addtree_all(iz,ix,ir)
         ELSE
            CALL updtree_all(iz,ix,ir)
         ENDIF
      ENDIF
   END SUBROUTINE fouds2_out_up

   SUBROUTINE fouds2_out_down(iz,ix,ir)
      IMPLICIT NONE
      INTEGER :: i,j,k,i2,j2,k2,ir,ix,iz,tsw1
      INTEGER :: swi,swj,swk,swsol
      INTEGER :: rtmp,xtmp,ztmp
      REAL(KIND=i10) :: old_value
      REAL(KIND=i10) :: trav,travm,slown,tdsh,tref,tdiv
      REAL(KIND=i10) :: a,b,c,u,v,w,em,en,ri,risti,rd1
      tsw1=0
      slown=1.0/veln(iz,ix,ir)
      ri=rgor-(roffset-1)*dnr-(ir-1)*dnr+earth
      risti=ri*sin(rgox+(xoffset-1)*dnx+(ix-1)*dnx)      
      swi=-1
      IF(ttn_ghost_down(iz,ix,1).GT.ttn_ghost_down(iz,ix,2))THEN
         swi=0
      ENDIF
      DO j=ix-1,ix+1,2
         swj=-1
         IF(j.eq.ix-1)THEN
            j2=j-1
            IF(j2.GE.1.AND.ttn(iz,j2,ir).NE.0)THEN
              swj=0
            ENDIF
         ELSE
            j2=j+1
            IF(j2.LE.xnum.AND.ttn(iz,j2,ir).NE.0)THEN
               swj=0
            ENDIF
         ENDIF
         IF(ttn(iz,j,ir).LT.ttn(iz,ix,ir).AND.swj.EQ.0.AND.ttn(iz,j,ir).NE.0)THEN
            swj=-1
            IF(ttn(iz,j,ir).GT.ttn(iz,j2,ir))THEN
               swj=0
            ENDIF
         ELSE
            swj=-1
         ENDIF         
         swk=-1
         IF(ttn_ghost_out(1,ix,ir).GT.ttn_ghost_out(2,ix,ir))THEN
            swk=0
         ENDIF         
         IF(j.GE.1.AND.j.LE.xnum)THEN
            swsol=0
            IF(swi.EQ.0)THEN
               swsol=1
               IF(swj.EQ.0)THEN
                  IF(swk.EQ.0)THEN
                     u=2.0*dnr
                     v=2.0*ri*dnx
                     w=2.0*risti*dnz
                     em=4.0*ttn_ghost_down(iz,ix,1)-ttn_ghost_down(iz,ix,2)-4.0*ttn(iz,j,ir)
                     em=em+ttn(iz,j2,ir)
                     en=4.0*ttn_ghost_down(iz,ix,1)-ttn_ghost_down(iz,ix,2)-4.0*ttn_ghost_out(1,ix,ir)
                     en=en+ttn_ghost_out(2,ix,ir)
                     a=v**2*w**2+u**2*w**2+u**2*v**2
                     b=2.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=4.0*ttn_ghost_down(iz,ix,1)-ttn_ghost_down(iz,ix,2)
                     tdiv=3.0
                  ELSE
                     u=risti*dnz
                     v=2.0*ri*dnx
                     w=2.0*dnr
                     em=3.0*ttn_ghost_out(1,ix,ir)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                     en=3.0*ttn_ghost_out(1,ix,ir)-4.0*ttn_ghost_down(iz,ix,1)+ttn_ghost_down(iz,ix,2)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_out(1,ix,ir)
                     tdiv=1.0
                  ENDIF
               ELSEIF(ttn(iz,j,ir).LT.ttn(iz,ix,ir).AND.ttn(iz,j,ir).NE.0)THEN
                  IF(swk.EQ.0)THEN
                     u=ri*dnx
                     v=2.0*dnr
                     w=2.0*risti*dnz
                     em=3.0*ttn(iz,j,ir)-4.0*ttn_ghost_down(iz,ix,1)+ttn_ghost_down(iz,ix,2)
                     en=3.0*ttn(iz,j,ir)-4.0*ttn_ghost_out(1,ix,ir)+ttn_ghost_out(2,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn(iz,j,ir)
                     tdiv=1.0
                  ELSE
                     u=ri*dnx
                     v=2.0*dnr
                     w=risti*dnz
                     em=3.0*ttn(iz,j,ir)-4.0*ttn_ghost_down(iz,ix,1)+ttn_ghost_down(iz,ix,2)
                     en=ttn(iz,j,ir)-ttn_ghost_out(1,ix,ir)
                     a=w**2*v**2+9.0*u**2*w**2+u**2*v**2
                     b=u**2*(6.0*em*w**2+2.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn(iz,j,ir)
                     tdiv=1.0
                  ENDIF
               ELSE
                  IF(swk.EQ.0)THEN
                     u=2.0*dnr
                     v=2.0*risti*dnz
                     em=4.0*ttn_ghost_down(iz,ix,1)-ttn_ghost_down(iz,ix,2)-4.0*ttn_ghost_out(1,ix,ir)
                     em=em+ttn_ghost_out(2,ix,ir)
                     a=u**2+v**2
                     b=2.0*em*u**2
                     c=u**2*(em**2-v**2*slown**2)
                     tref=4.0*ttn_ghost_down(iz,ix,1)-ttn_ghost_down(iz,ix,2)
                     tdiv=3.0
                  ELSE
                     u=risti*dnz
                     v=2.0*dnr
                     em=3.0*ttn_ghost_out(1,ix,ir)-4.0*ttn_ghost_down(iz,ix,1)+ttn_ghost_down(iz,ix,2)
                     a=v**2+9.0*u**2
                     b=6.0*em*u**2
                     c=u**2*(em**2-v**2*slown**2)
                     tref=ttn_ghost_out(1,ix,ir)
                     tdiv=1.0
                  ENDIF
               ENDIF
            ELSE
               swsol=1
               IF(swj.EQ.0)THEN
                  IF(swk.EQ.0)THEN
                     u=dnr
                     v=2.0*ri*dnx
                     w=2.0*risti*dnz
                     em=3.0*ttn_ghost_down(iz,ix,1)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                     en=3.0*ttn_ghost_down(iz,ix,1)-4.0*ttn_ghost_out(1,ix,ir)+ttn_ghost_out(2,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                     b=6.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_down(iz,ix,1)
                     tdiv=1.0
                  ELSE
                     u=dnr
                     v=2.0*ri*dnx
                     w=risti*dnz
                     em=3.0*ttn_ghost_down(iz,ix,1)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                     en=ttn_ghost_down(iz,ix,1)-ttn_ghost_out(1,ix,ir)
                     a=v**2*w**2+9.0*u**2*w**2+u**2*v**2
                     b=u**2*(6.0*em*w**2+2.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=ttn_ghost_down(iz,ix,1)
                     tdiv=1.0
                  ENDIF
               ELSE IF(ttn(iz,j,ir).LT.ttn(iz,ix,ir).AND.ttn(iz,j,ir).NE.0)THEN
                  IF(swk.EQ.0)THEN
                     u=dnr
                     v=ri*dnx
                     w=2.0*risti*dnz
                     em=ttn_ghost_down(iz,ix,1)-ttn(iz,j,ir)
                     en=3.0*ttn_ghost_down(iz,ix,1)-4.0*ttn_ghost_out(1,ix,ir)+ttn_ghost_out(2,ix,ir)
                     a=v**2*w**2+u**2*w**2+9.0*u**2*v**2
                     b=u**2*(2.0*em*w**2+6.0*en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                     tref=ttn_ghost_down(iz,ix,1)
                     tdiv=1.0
                  ELSE
                     u=dnr
                     v=ri*dnx
                     w=risti*dnz
                     em=ttn(iz,j,ir)-ttn_ghost_down(iz,ix,1)
                     en=ttn_ghost_out(1,ix,ir)-ttn_ghost_down(iz,ix,1)
                     a=v**2*w**2+u**2*w**2+u**2*v**2
                     b=-2.0*u**2*(em*w**2+en*v**2)
                     c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                     tref=ttn_ghost_down(iz,ix,1)
                     tdiv=1.0
                  ENDIF
               ELSE
                  IF(swk.EQ.0)THEN
                     u=dnr
                     v=2.0*risti*dnz
                     em=3.0*ttn_ghost_down(iz,ix,1)-4.0*ttn_ghost_out(1,ix,ir)+ttn_ghost_out(2,ix,ir)
                     a=v**2+9.0*u**2
                     b=6.0*em*u**2
                     c=u**2*(em**2-v**2*slown**2)
                     tref=ttn_ghost_down(iz,ix,1)
                     tdiv=1.0
                  ELSE 
                     u=dnr
                     v=risti*dnz
                     em=ttn_ghost_out(1,ix,ir)-ttn_ghost_down(iz,ix,1)
                     a=u**2+v**2
                     b=-2.0*u**2*em
                     c=u**2*(em**2-v**2*slown**2)
                     tref=ttn_ghost_down(iz,ix,1)
                     tdiv=1.0
                  ENDIF
               ENDIF            
            ENDIF
            IF(swsol.EQ.1)THEN
               rd1=b**2-4.0*a*c
               IF(rd1.LT.0.0)rd1=0.0
               tdsh=(-b+sqrt(rd1))/(2.0*a)
               trav=(tref+tdsh)/tdiv
               IF(tsw1.EQ.1)THEN
                  travm=MIN(trav,travm)
               ELSE
                  travm=trav
                  tsw1=1
               ENDIF
            ENDIF
         ENDIF         
      ENDDO      
      IF(travm.LT.ttn(iz,ix,ir))THEN
         ttn(iz,ix,ir)=travm
         IF(nsts(iz,ix,ir).EQ.0)THEN
            CALL addtree_all(iz,ix,ir)
         ELSE  
            CALL updtree_all(iz,ix,ir)
         ENDIF
      ENDIF
   END SUBROUTINE fouds2_out_down

   SUBROUTINE downtree_in(btgg,temp)
      IMPLICIT NONE
      INTEGER :: tpp,tpc
      INTEGER :: temp
      REAL(KIND=i10) :: rd1,rd2
      TYPE(temppointer) :: exchh
      TYPE(temppointer), DIMENSION (:) :: btgg
      IF(temp.EQ.1)THEN
         temp=temp-1
         RETURN
      ENDIF
      btgg(1)=btgg(temp)
      temp=temp-1
      tpp=1
      tpc=2*tpp
      DO WHILE(tpc.lt.temp)
         rd1=ttn_ghost_in(2,btgg(tpc)%pa,btgg(tpc)%pb)
         rd2=ttn_ghost_in(2,btgg(tpc+1)%pa,btgg(tpc+1)%pb)
         IF(rd1.gt.rd2)THEN
            tpc=tpc+1
         ENDIF
         rd1=ttn_ghost_in(2,btgg(tpc)%pa,btgg(tpc)%pb)
         rd2=ttn_ghost_in(2,btgg(tpp)%pa,btgg(tpp)%pb)
         IF(rd1.lt.rd2)THEN
            exchh=btgg(tpc)
            btgg(tpc)=btgg(tpp)
            btgg(tpp)=exchh
            tpp=tpc
            tpc=2*tpp
         ELSE
            tpc=temp+1
         ENDIF
      ENDDO
      IF(tpc.eq.temp)THEN
         rd1=ttn_ghost_in(2,btgg(tpc)%pa,btgg(tpc)%pb)
         rd2=ttn_ghost_in(2,btgg(tpp)%pa,btgg(tpp)%pb)
         IF(rd1.lt.rd2)THEN
            exchh=btgg(tpc)
            btgg(tpc)=btgg(tpp)
            btgg(tpp)=exchh
         ENDIF
      ENDIF
   END SUBROUTINE downtree_in

   SUBROUTINE downtree_out(btgg,temp)
      IMPLICIT NONE
      INTEGER :: tpp,tpc
      INTEGER :: temp
      REAL(KIND=i10) :: rd1,rd2
      TYPE(temppointer) :: exchh
      TYPE(temppointer), DIMENSION (:):: btgg
      IF(temp.EQ.1)THEN
         temp=temp-1
         RETURN
      ENDIF
      btgg(1)=btgg(temp)
      temp=temp-1
      tpp=1
      tpc=2*tpp      
      DO WHILE(tpc.lt.temp)
         rd1=ttn_ghost_out(1,btgg(tpc)%pa,btgg(tpc)%pb)
         rd2=ttn_ghost_out(1,btgg(tpc+1)%pa,btgg(tpc+1)%pb)
         IF(rd1.gt.rd2)THEN
            tpc=tpc+1
         ENDIF
         rd1=ttn_ghost_out(1,btgg(tpc)%pa,btgg(tpc)%pb)
         rd2=ttn_ghost_out(1,btgg(tpp)%pa,btgg(tpp)%pb)
         IF(rd1.lt.rd2)THEN
            exchh=btgg(tpc)
            btgg(tpc)=btgg(tpp)
            btgg(tpp)=exchh
            tpp=tpc
            tpc=2*tpp
         ELSE
            tpc=temp+1
         ENDIF
      ENDDO
      IF(tpc.eq.temp)THEN
         rd1=ttn_ghost_out(1,btgg(tpc)%pa,btgg(tpc)%pb)
         rd2=ttn_ghost_out(1,btgg(tpp)%pa,btgg(tpp)%pb)
         IF(rd1.lt.rd2)THEN
            exchh=btgg(tpc)
            btgg(tpc)=btgg(tpp)
            btgg(tpp)=exchh
         ENDIF
      ENDIF
   END SUBROUTINE downtree_out

   SUBROUTINE downtree_left(btgg,temp)
      IMPLICIT NONE
      INTEGER :: tpp,tpc
      INTEGER :: temp
      REAL(KIND=i10) :: rd1,rd2
      TYPE(temppointer) :: exchh
      TYPE(temppointer), DIMENSION (:) :: btgg
      IF(temp.EQ.1)THEN
         temp=temp-1
         RETURN
      ENDIF
      btgg(1)=btgg(temp)
      temp=temp-1
      tpp=1
      tpc=2*tpp
      DO WHILE(tpc.lt.temp)
         rd1=ttn_ghost_left(btgg(tpc)%pa,2,btgg(tpc)%pb)
         rd2=ttn_ghost_left(btgg(tpc+1)%pa,2,btgg(tpc+1)%pb)
         IF(rd1.gt.rd2)THEN
            tpc=tpc+1
         ENDIF
         rd1=ttn_ghost_left(btgg(tpc)%pa,2,btgg(tpc)%pb)
         rd2=ttn_ghost_left(btgg(tpp)%pa,2,btgg(tpp)%pb)
         IF(rd1.lt.rd2)THEN
            exchh=btgg(tpc)
            btgg(tpc)=btgg(tpp)
            btgg(tpp)=exchh
            tpp=tpc
            tpc=2*tpp
         ELSE
            tpc=temp+1
         ENDIF
      ENDDO
      IF(tpc.eq.temp)THEN
         rd1=ttn_ghost_left(btgg(tpc)%pa,2,btgg(tpc)%pb)
         rd2=ttn_ghost_left(btgg(tpp)%pa,2,btgg(tpp)%pb)
         IF(rd1.lt.rd2)THEN
            exchh=btgg(tpc)
            btgg(tpc)=btgg(tpp)
            btgg(tpp)=exchh
         ENDIF
      ENDIF
   END SUBROUTINE downtree_left

   SUBROUTINE downtree_right(btgg,temp)
      IMPLICIT NONE
      INTEGER :: tpp,tpc
      INTEGER :: temp
      REAL(KIND=i10) :: rd1,rd2
      TYPE(temppointer) :: exchh
      TYPE(temppointer), DIMENSION (:) :: btgg
      IF(temp.EQ.1)THEN
         temp=temp-1
         RETURN
      ENDIF
      btgg(1)=btgg(temp)
      temp=temp-1
      tpp=1
      tpc=2*tpp
      DO WHILE(tpc.lt.temp)
         rd1=ttn_ghost_right(btgg(tpc)%pa,1,btgg(tpc)%pb)
         rd2=ttn_ghost_right(btgg(tpc+1)%pa,1,btgg(tpc+1)%pb)
         IF(rd1.gt.rd2)THEN
            tpc=tpc+1
         ENDIF
         rd1=ttn_ghost_right(btgg(tpc)%pa,1,btgg(tpc)%pb)
         rd2=ttn_ghost_right(btgg(tpp)%pa,1,btgg(tpp)%pb)
         IF(rd1.lt.rd2)THEN
            exchh=btgg(tpc)
            btgg(tpc)=btgg(tpp)
            btgg(tpp)=exchh
            tpp=tpc
            tpc=2*tpp
         ELSE
            tpc=temp+1
         ENDIF
      ENDDO
      IF(tpc.eq.temp)THEN
         rd1=ttn_ghost_right(btgg(tpc)%pa,1,btgg(tpc)%pb)
         rd2=ttn_ghost_right(btgg(tpp)%pa,1,btgg(tpp)%pb)
         IF(rd1.lt.rd2)THEN
            exchh=btgg(tpc)
            btgg(tpc)=btgg(tpp)
            btgg(tpp)=exchh
         ENDIF
      ENDIF
   END SUBROUTINE downtree_right

   SUBROUTINE downtree_up(btgg,temp)
      IMPLICIT NONE
      INTEGER :: tpp,tpc
      INTEGER :: temp
      REAL(KIND=i10) :: rd1,rd2
      TYPE(temppointer) :: exchh
      TYPE(temppointer), DIMENSION (:) :: btgg
      IF(temp.EQ.1)THEN
         temp=temp-1
         RETURN
      ENDIF
      btgg(1)=btgg(temp)
      temp=temp-1
      tpp=1
      tpc=2*tpp
      DO WHILE(tpc.lt.temp)
         rd1=ttn_ghost_up(btgg(tpc)%pa,btgg(tpc)%pb,2)
         rd2=ttn_ghost_up(btgg(tpc+1)%pa,btgg(tpc+1)%pb,2)
         IF(rd1.gt.rd2)THEN
            tpc=tpc+1
         ENDIF
         rd1=ttn_ghost_up(btgg(tpc)%pa,btgg(tpc)%pb,2)
         rd2=ttn_ghost_up(btgg(tpp)%pa,btgg(tpp)%pb,2)
         IF(rd1.lt.rd2)THEN
            exchh=btgg(tpc)
            btgg(tpc)=btgg(tpp)
            btgg(tpp)=exchh
            tpp=tpc
            tpc=2*tpp
         ELSE
            tpc=temp+1
         ENDIF
      ENDDO
      IF(tpc.eq.temp)THEN
         rd1=ttn_ghost_up(btgg(tpc)%pa,btgg(tpc)%pb,2)
         rd2=ttn_ghost_up(btgg(tpp)%pa,btgg(tpp)%pb,2)
         IF(rd1.lt.rd2)THEN
            exchh=btgg(tpc)
            btgg(tpc)=btgg(tpp)
            btgg(tpp)=exchh
         ENDIF
      ENDIF
   END SUBROUTINE downtree_up

   SUBROUTINE downtree_down(btgg,temp)
      IMPLICIT NONE
      INTEGER :: tpp,tpc
      INTEGER :: temp
      REAL(KIND=i10) :: rd1,rd2
      TYPE(temppointer) :: exchh
      TYPE(temppointer), DIMENSION (:) :: btgg
      IF(temp.EQ.1)THEN
         temp=temp-1
         RETURN
      ENDIF
      btgg(1)=btgg(temp)
      temp=temp-1
      tpp=1
      tpc=2*tpp
      DO WHILE(tpc.lt.temp)
         rd1=ttn_ghost_down(btgg(tpc)%pa,btgg(tpc)%pb,1)
         rd2=ttn_ghost_down(btgg(tpc+1)%pa,btgg(tpc+1)%pb,1)
         IF(rd1.gt.rd2)THEN
            tpc=tpc+1
         ENDIF
         rd1=ttn_ghost_down(btgg(tpc)%pa,btgg(tpc)%pb,1)
         rd2=ttn_ghost_down(btgg(tpp)%pa,btgg(tpp)%pb,1)
         IF(rd1.lt.rd2)THEN
            exchh=btgg(tpc)
            btgg(tpc)=btgg(tpp)
            btgg(tpp)=exchh
            tpp=tpc
            tpc=2*tpp
         ELSE
            tpc=temp+1
         ENDIF
      ENDDO
      IF(tpc.eq.temp)THEN
         rd1=ttn_ghost_down(btgg(tpc)%pa,btgg(tpc)%pb,1)
         rd2=ttn_ghost_down(btgg(tpp)%pa,btgg(tpp)%pb,1)
         IF(rd1.lt.rd2)THEN
            exchh=btgg(tpc)
            btgg(tpc)=btgg(tpp)
            btgg(tpp)=exchh
         ENDIF
      ENDIF
   END SUBROUTINE downtree_down
END MODULE traveltime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM fmmin3d
   USE globalp
   USE traveltime
   IMPLICIT NONE
   include 'mpif.h'
   EXTERNAL sgn
   CHARACTER (LEN=20) :: sources,itimes,receivers,gridv
   CHARACTER (LEN=20) :: rtravel,wrays,cdum
   CHARACTER (LEN=20) :: frechet,pht
   INTEGER :: ii,i,temp,j,k,l,p,ll,nsrc,fsrt,cfd,tnr,nra,idum,count
   INTEGER :: awttf,maxbt
   INTEGER :: ierr,status(mpi_status_size),srequest,root,nmp,tempp,rtmp,xtmp,ztmp
   INTEGER :: array_of_status(mpi_status_size,6)
   INTEGER :: cal_all_process
   INTEGER :: temp_comm
   INTEGER :: myfile,filetype
   INTEGER :: color,key
   INTEGER :: root_wfile,fdm_start,fdm_end
   INTEGER :: sgn
   REAL :: fracg
   REAL(KIND=i10) :: cslat,cslong,abw,abe,abn,abs
   REAL(KIND=i10) :: tol,start_time,end_time,total_time
   REAL(KIND=i10), DIMENSION (:), ALLOCATABLE :: scr,scx,scz
   INTEGER :: sw_coln
   REAL(KIND=i5) :: sw_fdm
   REAL(KIND=i5), DIMENSION (:), ALLOCATABLE :: fdm_file,temp_fdm_file,ttemp_fdm_file
   REAL(KIND=i10), DIMENSION (4) :: sww
   INTEGER, DIMENSION (:), ALLOCATABLE :: coln_file,temp_coln_file,ttemp_coln_file
   INTEGER, DIMENSION (2) :: array_of_size,array_of_subsize,array_of_starts
   INTEGER, DIMENSION (6) :: rrequest
   LOGICAL, DIMENSION (3) :: periods
   LOGICAL ::reorder
   INTEGER, DIMENSION (3) :: dims
   INTEGER :: comm_x, comm_z, comm_r
   INTEGER(KIND=MPI_OFFSET_KIND), parameter :: zero_off = 0
   CALL MPI_INIT(ierr)
   start_time = MPI_Wtime()
   nr=2
   nx=2
   nz=2   
   root_wfile=nr*nx*nz
   nprocs=root_wfile
   CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
   periods(1) = .FALSE.
   periods(2) = .FALSE.
   periods(3) = .FALSE.
   reorder = .TRUE.
   tol=0.0001
   dims(:) = [nr,nx,nz]
   CALL MPI_CART_CREATE(MPI_COMM_WORLD,3,dims,periods,reorder,comm_cart,ierr)
   IF(comm_cart.NE.MPI_COMM_NULL)THEN 
      CALL MPI_COMM_SIZE(comm_cart,nprocs,ierr)
      periods(1) = .TRUE.
      periods(2) = .FALSE.
      periods(3) = .FALSE.
      CALL MPI_CART_SUB(comm_cart,periods,comm_r,ierr)
      CALL MPI_COMM_RANK(comm_r,rank_r,ierr)
      periods(1) = .FALSE.
      periods(2) = .TRUE.
      periods(3) = .FALSE.
      CALL MPI_CART_SUB(comm_cart,periods,comm_x,ierr)
      CALL MPI_COMM_RANK(comm_x,rank_x,ierr)
      periods(1) = .FALSE.
      periods(2) = .FALSE.
      periods(3) = .TRUE.
      CALL MPI_CART_SUB(comm_cart,periods,comm_z,ierr)
      CALL MPI_COMM_RANK(comm_z,rank_z,ierr)
      OPEN(UNIT=10,FILE='fm3dt_parallel.in',STATUS='old')
      READ(10,1)cdum
      READ(10,1)cdum
      READ(10,1)cdum
      READ(10,1)sources
      READ(10,1)itimes
      READ(10,1)receivers
      READ(10,1)gridv
      READ(10,*)earth
      READ(10,*)fom
      READ(10,*)ltfr
      READ(10,*)cslat,cslong
      READ(10,1)cdum
      READ(10,1)cdum
      READ(10,1)cdum
      READ(10,*)fsrt
      READ(10,1)rtravel           
      READ(10,*)cfd
      READ(10,1)frechet
      1   FORMAT(a20)
      CLOSE(10)
      OPEN(UNIT=10,FILE=gridv,STATUS='old')
      READ(10,*)nvr,nvx,nvz
      READ(10,*)gor,goxd,gozd
      READ(10,*)dvr,dvx,dvz
      READ(10,*)nnr,nnx,nnz
      READ(10,*)rgor,rgoxd,rgozd
      READ(10,*)dnr,dnxd,dnzd
      dvx=dvx*pi/180.0
      dvz=dvz*pi/180.0
      gox=(90.0-goxd)*pi/180.0
      goz=gozd*pi/180.0
      dnx=dnxd*pi/180.0
      dnz=dnzd*pi/180.0
      rgox=(90.0-rgoxd)*pi/180.0
      rgoz=rgozd*pi/180.0
      ALLOCATE(velv(0:nvr+1,0:nvx+1,0:nvz+1))
      DO i=0,nvz+1
         DO j=0,nvx+1
            DO k=0,nvr+1
               READ(10,*)velv(k,j,i)
            ENDDO
         ENDDO
      ENDDO
      CLOSE(10)
      tempp = mod(nnr,nr)
      IF(tempp.EQ.0)THEN
         rloca = nnr/nr
      ELSEIF(rank_r.LT.tempp)THEN
         rloca = nnr/nr + 1
      ELSE
         rloca = nnr/nr
      ENDIF
      tempp = mod(nnx,nx)
      IF(tempp.EQ.0)THEN
         xloca = nnx/nx
      ELSEIF(rank_x.LT.tempp)THEN
         xloca = nnx/nx + 1
      ELSE
         xloca = nnx/nx
      ENDIF
      tempp = mod(nnz,nz)
      IF(tempp.EQ.0)THEN
         zloca = nnz/nz
      ELSEIF(rank_z.LT.tempp)THEN
         zloca = nnz/nz + 1
      ELSE
         zloca = nnz/nz
      ENDIF
      xnum = xloca
      znum = zloca
      rnum = rloca
      xloca = xloca + 2 + 2*sgn(INT(mod(rank_x,nx-1)))
      zloca = zloca + 2 + 2*sgn(INT(mod(rank_z,nz-1)))
      rloca = rloca + 2 + 2*sgn(INT(mod(rank_r,nr-1)))
      rtmp=rloca - 2 - 2
      CALL MPI_SCAN(rtmp,roffset,1,MPI_INTEGER,MPI_SUM,comm_r,ierr)
      roffset = roffset - rtmp
      xtmp=xloca - 2 - 2
      CALL MPI_SCAN(xtmp,xoffset,1,MPI_INTEGER,MPI_SUM,comm_x,ierr)
      xoffset = xoffset - xtmp
      ztmp=zloca - 2 - 2
      CALL MPI_SCAN(ztmp,zoffset,1,MPI_INTEGER,MPI_SUM,comm_z,ierr)
      zoffset = zoffset - ztmp     
      roffset=roffset+1+2*sgn(rank_r)
      xoffset=xoffset+1+2*sgn(rank_x)
      zoffset=zoffset+1+2*sgn(rank_z)
      CALL MPI_CART_SHIFT(comm_cart,0,1,up,down,ierr)
      CALL MPI_CART_SHIFT(comm_cart,1,1,left,right,ierr)
      CALL MPI_CART_SHIFT(comm_cart,2,1,in,out,ierr)      
      CALL gridder(gridv)
      DEALLOCATE(velv)
      fracg=0.1
      nmp=(nvr+2)*(nvx+2)*(nvz+2)
      Open(UNIT=60,FILE=sources,STATUS='old')
      READ(60,*)nra
      IF(fsrt.eq.1)THEN
         OPEN(UNIT=70,FILE=receivers,STATUS='old')
         READ(70,*)idum
         IF(idum.NE.nra)THEN
            WRITE(6,*)'ERROR!!!'
            WRITE(6,*)'Source and receiver files are'
            WRITE(6,*)'inconsistent!!!!'
            WRITE(6,*)'First line of each file should'
            WRITE(6,*)'be identical!!!'
            WRITE(6,*)'TERMINATING PROGRAM!!!'
            STOP
         ENDIF
      ENDIF   
      ALLOCATE(ttn(0:znum+1,0:xnum+1,0:rnum+1),ttn_previous_step(znum,xnum,rnum),&
               & nsts(0:znum+1,0:xnum+1,0:rnum+1), STAT=checkstat)
      IF(checkstat > 0)THEN
         WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin3d: REAL ttn,nsts'
      ENDIF
      ALLOCATE(ttn_ghost_left(znum,2,rnum),ttn_ghost_right(znum,2,rnum))
      ALLOCATE(ttn_ghost_in(2,xnum,rnum),ttn_ghost_out(2,xnum,rnum))
      ALLOCATE(ttn_ghost_up(znum,xnum,2),ttn_ghost_down(znum,xnum,2))
      ALLOCATE(ttn1(znum,xnum),ttn2(znum,xnum))
      periods(1) = .FALSE.
      periods(2) = .TRUE.
      periods(3) = .TRUE.
      CALL MPI_CART_SUB(comm_cart,periods,comm_srtimes,ierr)   
      periods(1) = .FALSE.
      periods(2) = .TRUE.
      periods(3) = .TRUE.
      CALL MPI_CART_SUB(comm_cart,periods,comm_start,ierr)
      IF(rank_r.EQ.nr-1)THEN
         array_of_size(1) = nnz
         array_of_size(2) = nnx
         array_of_subsize(1) = znum
         array_of_subsize(2) = xnum
         array_of_starts(1) = zoffset-1
         array_of_starts(2) = xoffset-1
         CALL MPI_FILE_OPEN(comm_start,"itimes.out",MPI_MODE_RDONLY,MPI_INFO_NULL,myfile,ierr)
         CALL MPI_TYPE_CREATE_SUBARRAY(2,array_of_size,array_of_subsize,array_of_starts,MPI_ORDER_FORTRAN,&
                                    & MPI_DOUBLE_PRECISION,filetype,ierr)
         CALL MPI_TYPE_COMMIT(filetype,ierr)
         CALL MPI_FILE_SET_VIEW(myfile,zero_off,MPI_DOUBLE_PRECISION,filetype,"native",MPI_INFO_NULL,ierr)
      ENDIF
      ALLOCATE(rcounts_rpaths(nprocs),displs_rpaths(nprocs),STAT=checkstat)
      IF(checkstat>0)THEN
         WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin3d: REAL rcounts_rpaths,displs_rpaths'
      ENDIF
      CALL type_create_for_communication(ierr)
      CALL MPI_PACK_SIZE(2*xnum*rnum,MPI_DOUBLE_PRECISION,comm_cart,outsize,ierr)
      ALLOCATE(BUF(outsize))
      ALLOCATE(LBUF1(outsize),LBUF2(outsize), STAT=checkstat)
      IF(checkstat > 0)THEN
         WRITE(6,*)'Error with ALLOCATE: SUBROUTINE travel: TYPE BUF'
      ENDIF
      ALLOCATE(RBUF1(outsize),RBUF2(outsize), STAT=checkstat)
      IF(checkstat > 0)THEN
         WRITE(6,*)'Error with ALLOCATE: SUBROUTINE travel: TYPE BUF'
      ENDIF
   ELSE      
      OPEN(UNIT=10,FILE='fm3dt_parallel.in',STATUS='old')
      READ(10,1)cdum
      READ(10,1)cdum
      READ(10,1)cdum
      READ(10,1)cdum
      READ(10,1)cdum
      READ(10,1)cdum
      READ(10,1)gridv
      READ(10,*)cdum
      READ(10,*)cdum
      READ(10,*)cdum
      READ(10,*)cdum
      READ(10,1)cdum
      READ(10,1)cdum
      READ(10,1)cdum
      READ(10,*)cdum
      READ(10,1)cdum           
      READ(10,*)cfd
      READ(10,1)frechet
      CLOSE(10)    
      IF(cfd.EQ.1)THEN
         OPEN(UNIT=50,FILE=frechet,FORM='unformatted',STATUS='unknown')
      ENDIF      
      ALLOCATE(displs_fdm(nprocs+1),rcounts_fdm(nprocs+1), STAT=checkstat)
      IF(checkstat > 0)THEN
         WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin3d: REAL displs_fdm,rcounts_fdm'
      ENDIF
      rcounts_process=0
   ENDIF
   IF(myrank==0)THEN
      CALL MPI_SEND(nra,1,MPI_INTEGER,root_wfile,111,MPI_COMM_WORLD,ierr)
      CALL MPI_SEND(nmp,1,MPI_INTEGER,root_wfile,111,MPI_COMM_WORLD,ierr)
   ELSEIF(myrank==root_wfile)THEN
      CALL MPI_RECV(nra,1,MPI_INTEGER,0,111,MPI_COMM_WORLD,status,ierr)
      CALL MPI_RECV(nmp,1,MPI_INTEGER,0,111,MPI_COMM_WORLD,status,ierr)
   ENDIF  
   DO ii=1,nra
      IF(comm_cart.NE.MPI_COMM_NULL)THEN 
         READ(60,*)nsrc
         ALLOCATE(scr(nsrc),scx(nsrc),scz(nsrc), STAT=checkstat)
         IF(checkstat > 0)THEN
            WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin3d: REAL scr,scx,scz'
         ENDIF       
         DO i=1,nsrc
            READ(60,*)scx(i),scz(i),scr(i)
            READ(60,'(a8)')pht
            scx(i)=(90.0-scx(i))*pi/180.0
            scz(i)=scz(i)*pi/180.0
         ENDDO      
         IF(fsrt.EQ.1.or.ltfr.EQ.1)THEN
            READ(70,*)nrc
            ALLOCATE(rcr(nrc),rcx(nrc),rcz(nrc), STAT=checkstat)
            IF(checkstat > 0)THEN
               WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin3d: REAL rcr,rcx,rcz'
            ENDIF
            DO i=1,nrc
               READ(70,*)rcr(i),rcx(i),rcz(i)
               IF(ltfr.EQ.1)THEN
                  IF(i.EQ.1)THEN
                     abe=rcz(i)
                     abw=abe
                     abn=rcx(i)
                     abs=abn
                  ELSE
                     IF(rcz(i).LT.abw)THEN
                        abw=rcz(i)
                     ELSE IF(rcz(i).GT.abe)THEN
                        abe=rcz(i)
                     ENDIF
                     IF(rcx(i).LT.abs)THEN
                        abs=rcx(i)
                     ELSE IF(rcx(i).GT.abn)THEN
                        abn=rcx(i)
                     ENDIF
                  ENDIF
               ENDIF
               rcx(i)=(90.0-rcx(i))*pi/180.0
               rcz(i)=rcz(i)*pi/180.0
            ENDDO
            IF(ltfr.EQ.1)THEN
               abw=(abw-cslong)*pi/180.0
               abe=(abe+cslong)*pi/180.0
               abs=abs-cslat
               abn=abn+cslat
               abs=(90.0-abs)*pi/180.0
               abn=(90.0-abn)*pi/180.0            
               nsnn=INT((abn-rgox)/dnx)
               IF(nsnn.LT.1)nsnn=1               
               nsns=INT((abs-rgox)/dnx)+2
               IF(nsns.GT.nnx)nsns=nnx               
               nsne=INT((abe-rgoz)/dnz)+2
               IF(nsne.GT.nnz)nsne=nnz               
               nsnw=INT((abw-rgoz)/dnz)
               IF(nsnw.LT.1)nsnw=1            
            ENDIF
         ENDIF      
         IF(rank_r.EQ.0)THEN
            i=nrc/2             
            ALLOCATE(rec_location(3,i),number(i),STAT=checkstat)
            IF(checkstat > 0)THEN
               WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin3d: REAL rec_location,number'
            ENDIF
         ENDIF
         ALLOCATE(total_communicate_rpaths(4,nrc),communicate_only(4,nrc),STAT=checkstat)
         IF(checkstat > 0)THEN
            WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin3d: REAL total_communicate_rpaths, communicate_only'
         ENDIF
         DO i=1,nrc
            total_communicate_rpaths(1:3,i)=[rcr(i),rcx(i),rcz(i)]
            total_communicate_rpaths(4,i)=i
         ENDDO      
         i=fracg*nmp*nrc
         ALLOCATE(fdm_process(i),STAT=checkstat)
         IF(checkstat > 0)THEN
            WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin3d: REAL fdm_process'
         ENDIF   
         ALLOCATE(coln_process(i),STAT=checkstat)
         IF(checkstat > 0)THEN
            WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin3d: INTEGER coln_process'
         ENDIF 
      ELSE
         i=0.8*nmp
         ALLOCATE(ttemp_fdm_file(i),ttemp_coln_file(i), STAT=checkstat)
         IF(checkstat > 0)THEN
            WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin3d: REAL ttemp_fdm_file'
         ENDIF
      ENDIF
      IF(myrank==0)THEN
         CALL MPI_SEND(nsrc,1,MPI_INTEGER,root_wfile,111,MPI_COMM_WORLD,ierr)
         CALL MPI_SEND(nrc,1,MPI_INTEGER,root_wfile,112,MPI_COMM_WORLD,ierr)
      ELSEIF(myrank==root_wfile)THEN
         CALL MPI_RECV(nsrc,1,MPI_INTEGER,0,111,MPI_COMM_WORLD,status,ierr)
         CALL MPI_RECV(nrc,1,MPI_INTEGER,0,112,MPI_COMM_WORLD,status,ierr)
      ENDIF  
      DO i=1,nsrc
         IF(comm_cart.NE.MPI_COMM_NULL)THEN
            nsts=-1
            ntr=0
            ttn=0.0
            ttn_previous_step=0.0
            pre_error_back_propagation=0.0
            cur_error_back_propagation=0.0
            sw_fmm=0
            cal_already=0
            cal_all_process=0
            error_to_stop=100
            send_direction=0
            temp_fdm=1
            communicate_rpaths_length=nrc
            rcounts_process=0
            ttn_ghost_left=0.0
            ttn_ghost_right=0.0
            ttn_ghost_in=0.0
            ttn_ghost_out=0.0
            ttn_ghost_up=0.0
            ttn_ghost_down=0.0        
            maxbt=rnum*xnum*znum
            ALLOCATE(btg(maxbt), STAT=checkstat)
            IF(checkstat > 0)THEN
               WRITE(6,*)'Error with ALLOCATE: SUBROUTINE travel: TYPE btg'
            ENDIF            
            IF(rank_r.EQ.nr-1)THEN
               CALL start_fmm(myfile,i,ierr)
            ENDIF          
            CALL error_calculate(ierr)
            count=0          
            DO WHILE(error_to_stop.GT.0.00001)
               ttn_previous_step=ttn(1:znum,1:xnum,1:rnum)
               IF(sw_fmm.EQ.0)THEN
                  direction_send=0
               ELSE
                  direction_send=1
               ENDIF               
               sw_fmm=0             
               CALL MPI_ISEND(direction_send(1),1,MPI_INTEGER,left,11,comm_cart,srequest,ierr)
               CALL MPI_ISEND(direction_send(2),1,MPI_INTEGER,right,112,comm_cart,srequest,ierr)
               CALL MPI_ISEND(direction_send(3),1,MPI_INTEGER,in,113,comm_cart,srequest,ierr)
               CALL MPI_ISEND(direction_send(4),1,MPI_INTEGER,out,114,comm_cart,srequest,ierr)
               CALL MPI_ISEND(direction_send(5),1,MPI_INTEGER,up,115,comm_cart,srequest,ierr)
               CALL MPI_ISEND(direction_send(6),1,MPI_INTEGER,down,116,comm_cart,srequest,ierr)      
               
               CALL MPI_IRECV(send_direction(2),1,MPI_INTEGER,right,11,comm_cart,rrequest(2),ierr)
               CALL MPI_IRECV(send_direction(1),1,MPI_INTEGER,left,112,comm_cart,rrequest(1),ierr)
               CALL MPI_IRECV(send_direction(4),1,MPI_INTEGER,out,113,comm_cart,rrequest(4),ierr)
               CALL MPI_IRECV(send_direction(3),1,MPI_INTEGER,in,114,comm_cart,rrequest(3),ierr)
               CALL MPI_IRECV(send_direction(6),1,MPI_INTEGER,down,115,comm_cart,rrequest(6),ierr)
               CALL MPI_IRECV(send_direction(5),1,MPI_INTEGER,up,116,comm_cart,rrequest(5),ierr)    

               CALL MPI_WAITALL(6,rrequest,array_of_status,ierr)
               CALL MPI_BARRIER(comm_cart,ierr)
               CALL boundary_communication(ierr)
               CALL MPI_BARRIER(comm_cart,ierr)              
               CALL order_overlap     
               CALL local_fmm              
               CALL MPI_BARRIER(comm_cart,ierr)            
               count=count+1      
               CALL error_calculate(ierr)
               CALL MPI_BARRIER(comm_cart,ierr)      
            ENDDO
            IF(myrank==0)THEN
               print*,"Finish Source:",i,"Final Error",error_to_stop
            ENDIF         
            direction_send=1
            send_direction=1
            CALL boundary_communication(ierr)
            CALL MPI_BARRIER(comm_cart,ierr)            
            ttn(1:znum,0,1:rnum)=ttn_ghost_left(1:znum,2,1:rnum)
            ttn(1:znum,xnum+1,1:rnum)=ttn_ghost_right(1:znum,1,1:rnum)
            ttn(0,1:xnum,1:rnum)=ttn_ghost_in(2,1:xnum,1:rnum)
            ttn(znum+1,1:xnum,1:rnum)=ttn_ghost_out(1,1:xnum,1:rnum)
            ttn(1:znum,1:xnum,0)=ttn_ghost_up(1:znum,1:xnum,2)
            ttn(1:znum,1:xnum,rnum+1)=ttn_ghost_down(1:znum,1:xnum,1)
            CALL MPI_BARRIER(comm_cart,ierr)
            CALL cross_data_communication(ierr)
            CALL MPI_BARRIER(comm_cart,ierr)            
            IF(fsrt.eq.1.AND.rank_r.EQ.0)THEN
               CALL srtimes(i,ii,nra,fsrt,nsrc,rtravel,ierr)
            ENDIF       
         ENDIF       
         IF(cfd.EQ.1.OR.myrank.EQ.root_wfile)THEN            
            IF(comm_cart.NE.MPI_COMM_NULL)THEN
               ALLOCATE(communicate_rpaths(4,nrc), STAT=checkstat)
               IF(checkstat > 0)THEN
                  WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin3d: REAL communicate_rpaths'
               ENDIF
               communicate_rpaths=total_communicate_rpaths
               count=0               
               DO                  
                  CALL parallel_rpaths(i,ii,cfd,count,ierr)
                  temp=0                  
                  DO j=1,communicate_rpaths_length
                     IF(communicate_rpaths(3,j).NE.0)THEN  
                        temp=temp+1
                        communicate_only(:,temp)=communicate_rpaths(:,j)
                     ENDIF
                  ENDDO              
                  CALL MPI_ALLGATHER(temp,1,MPI_INTEGER,rcounts_rpaths,1,MPI_INTEGER,comm_cart,ierr)
                  CALL MPI_BARRIER(comm_cart,ierr)                  
                  tempp=SUM(rcounts_rpaths)   
                  DEALLOCATE(communicate_rpaths)                  
                  IF(tempp.GT.0)THEN
                     ALLOCATE(communicate_rpaths(4,tempp))
                     communicate_rpaths_length=tempp
                  ELSE
                     EXIT
                  ENDIF
                  rcounts_rpaths=4*rcounts_rpaths
                  displs_rpaths=0
                  DO l=2,nprocs
                     displs_rpaths(l)=displs_rpaths(l-1)+rcounts_rpaths(l-1)
                  ENDDO
                  CALL MPI_ALLGATHERV(communicate_only,4*temp,MPI_DOUBLE_PRECISION,communicate_rpaths,&
                                    & rcounts_rpaths,displs_rpaths,MPI_DOUBLE_PRECISION,comm_cart,ierr)
                  CALL MPI_BARRIER(comm_cart,ierr)     
                  DO l=1,communicate_rpaths_length-1
                     DO j=1,communicate_rpaths_length-l-1
                        IF(communicate_rpaths(4,j).GT.communicate_rpaths(4,j+1))THEN
                           sww=communicate_rpaths(:,j)
                           communicate_rpaths(:,j)=communicate_rpaths(:,j+1)
                           communicate_rpaths(:,j+1)=sww
                        ENDIF
                     ENDDO
                  ENDDO
                  count=count+1
               ENDDO
            ENDIF              
            CALL MPI_GATHER(rcounts_process,1,MPI_INTEGER,rcounts_fdm,1,MPI_INTEGER,root_wfile,MPI_COMM_WORLD,ierr)
            CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)  
            IF(comm_cart.EQ.MPI_COMM_NULL)THEN
               displs_fdm=0
               DO l=2,nprocs+1
                  displs_fdm(l)=displs_fdm(l-1)+rcounts_fdm(l-1)
               ENDDO     
               ALLOCATE(fdm_process(1),coln_process(1),STAT=checkstat)
               IF(checkstat > 0)THEN
                  WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin3d: REAL fdm_process,coln_process'
               ENDIF   
               fdm_process=0.0
               coln_process=0
               j=SUM(rcounts_fdm)
               ALLOCATE(fdm_file(j),temp_fdm_file(j),STAT=checkstat)
               IF(checkstat > 0)THEN
                  WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin3d: REAL fdm_file'
               ENDIF
               ALLOCATE(coln_file(j),temp_coln_file(j),STAT=checkstat)
               IF(checkstat > 0)THEN
                  WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin3d: INTEGER coln_file'
               ENDIF
            ENDIF                     
            CALL MPI_GATHERV(fdm_process,rcounts_process,MPI_FLOAT,temp_fdm_file,&
                           & rcounts_fdm,displs_fdm,MPI_FLOAT,root_wfile,MPI_COMM_WORLD,ierr)                 
            CALL MPI_GATHERV(coln_process,rcounts_process,MPI_INTEGER,temp_coln_file,&
                           & rcounts_fdm,displs_fdm,MPI_INTEGER,root_wfile,MPI_COMM_WORLD,ierr)           
            CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)            
            IF(comm_cart.EQ.MPI_COMM_NULL)THEN            
               j=sum(rcounts_fdm)
               fdm_start=1             
               DO k=1,nrc                  
                  fdm_end=fdm_start
                  l=1
                  idum=1
                  ttemp_coln_file=0
                  ttemp_fdm_file=0.0                 
                  DO WHILE(l.LT.j)
                     DO WHILE(NINT(temp_fdm_file(l)).NE.k)
                        l=l+temp_coln_file(l)+1
                        IF(l.GE.j)THEN
                           EXIT
                        ENDIF
                     ENDDO
                     IF(l.LT.j)THEN
                        IF(idum.EQ.1)THEN
                           tempp=temp_coln_file(l)
                           ttemp_coln_file(idum:idum+tempp)=temp_coln_file(l:l+tempp)
                           ttemp_fdm_file(idum:idum+tempp)=temp_fdm_file(l:l+tempp)
                           idum=idum+tempp+1
                           l=l+tempp+1
                        ELSE
                           tempp=temp_coln_file(l)
                           ttemp_coln_file(1)=ttemp_coln_file(1)+tempp
                           ttemp_coln_file(idum:idum+tempp-1)=temp_coln_file(l+1:l+tempp)
                           ttemp_fdm_file(idum:idum+tempp-1)=temp_fdm_file(l+1:l+tempp)
                           idum=idum+tempp
                           l=l+tempp+1
                        ENDIF
                     ENDIF
                  ENDDO
                  DO p=2,idum-1-1
                     DO ll=2,idum-p-1
                        IF(ttemp_coln_file(ll).GT.ttemp_coln_file(ll+1))THEN
                           sw_fdm=ttemp_fdm_file(ll)
                           ttemp_fdm_file(ll)=ttemp_fdm_file(ll+1)
                           ttemp_fdm_file(ll+1)=sw_fdm
                           sw_coln=ttemp_coln_file(ll)
                           ttemp_coln_file(ll)=ttemp_coln_file(ll+1)
                           ttemp_coln_file(ll+1)=sw_coln
                        ENDIF
                     ENDDO
                  ENDDO                
                  p=2
                  DO WHILE(p.LE.idum-1)
                     ll=p+1
                     DO WHILE(ttemp_coln_file(ll).EQ.ttemp_coln_file(p).AND.ll.LE.idum-1)
                        ttemp_fdm_file(p)=ttemp_fdm_file(p)+ttemp_fdm_file(ll)
                        ttemp_coln_file(1)=ttemp_coln_file(1)-1
                        ll=ll+1
                     ENDDO
                     fdm_end=fdm_end+1
                     fdm_file(fdm_end)=ttemp_fdm_file(p)
                     coln_file(fdm_end)=ttemp_coln_file(p)
                     p=ll
                  ENDDO
                  coln_file(fdm_start)=ttemp_coln_file(1)
                  fdm_file(fdm_start)=ttemp_fdm_file(1)    
                  fdm_start=fdm_end+1
               ENDDO
               DO k=1,fdm_end
                  WRITE(50)coln_file(k),fdm_file(k)
               ENDDO               
               DEALLOCATE(fdm_process,coln_process)
               DEALLOCATE(fdm_file,coln_file)
               DEALLOCATE(temp_fdm_file,temp_coln_file)             
            ENDIF              
         ENDIF        
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)         
         IF(comm_cart.NE.MPI_COMM_NULL)THEN          
            DEALLOCATE(btg)
         ENDIF        
      ENDDO
      IF(comm_cart.NE.MPI_COMM_NULL)THEN
         DEALLOCATE(scr,scx,scz)
         DEALLOCATE(communicate_only)
         DEALLOCATE(fdm_process,coln_process)
         DEALLOCATE(total_communicate_rpaths)
         IF(fsrt.EQ.1.or.ltfr.EQ.1)THEN
            DEALLOCATE(rcr,rcx,rcz)
         ENDIF
         IF(rank_r.EQ.0)THEN
            DEALLOCATE(rec_location,number)
         ENDIF         
      ELSE
         DEALLOCATE(ttemp_fdm_file,ttemp_coln_file)
      ENDIF
   ENDDO
   IF(comm_cart.NE.MPI_COMM_NULL)THEN
      CLOSE(60)
      IF(fsrt.eq.1)THEN
         CLOSE(70)
      ENDIF
      DEALLOCATE(veln,ttn,ttn_previous_step,nsts, STAT=checkstat)
      IF(checkstat > 0)THEN
         WRITE(6,*)'Error with DEALLOCATE: PROGRAM fmmin3d: final deallocate'
      ENDIF
      DEALLOCATE(ttn_ghost_left,ttn_ghost_right)
      DEALLOCATE(ttn_ghost_in,ttn_ghost_out)
      DEALLOCATE(ttn_ghost_up,ttn_ghost_down)     
      DEALLOCATE(rcounts_rpaths,displs_rpaths, STAT=checkstat)
      IF(checkstat > 0)THEN
         WRITE(6,*)'Error with DEALLOCATE: PROGRAM fmmin3d: final deallocate'
      ENDIF    
      CALL MPI_TYPE_FREE(type_left,ierr)
      CALL MPI_TYPE_FREE(type_right,ierr)
      CALL MPI_TYPE_FREE(type_up,ierr)
      CALL MPI_TYPE_FREE(type_down,ierr)
      CALL MPI_TYPE_FREE(type_recv1,ierr)
      CALL MPI_TYPE_FREE(type_recv2,ierr)
      CALL MPI_COMM_FREE(comm_r,ierr)
      CALL MPI_COMM_FREE(comm_x,ierr)
      CALL MPI_COMM_FREE(comm_z,ierr)
      DEALLOCATE(LBUF1,LBUF2,RBUF1,RBUF2)
      DEALLOCATE(BUF)
      IF(rank_r.EQ.nr-1)THEN
         CALL MPI_FILE_CLOSE(myfile,ierr)
         CALL MPI_TYPE_FREE(filetype,ierr)
      ENDIF
   ENDIF
   IF(comm_cart.EQ.MPI_COMM_NULL)THEN
      IF(cfd.EQ.1)THEN
         CLOSE(50)
      ENDIF
      DEALLOCATE(displs_fdm,rcounts_fdm)
   ENDIF
   CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)   
   CALL MPI_FINALIZE(ierr)
   STOP
END PROGRAM fmmin3d

SUBROUTINE type_create_for_communication(ierr)
   USE globalp
   IMPLICIT NONE
   include 'mpif.h'
   INTEGER :: ierr
   CALL MPI_TYPE_VECTOR(rnum,znum,(xnum+2)*(znum+2),MPI_DOUBLE_PRECISION,type_left,ierr)
   CALL MPI_TYPE_COMMIT(type_left,ierr)
   CALL MPI_TYPE_VECTOR(rnum,znum,(xnum+2)*(znum+2),MPI_DOUBLE_PRECISION,type_right,ierr)
   CALL MPI_TYPE_COMMIT(type_right,ierr)
   CALL MPI_TYPE_VECTOR(rnum,znum,2*znum,MPI_DOUBLE_PRECISION,type_recv1,ierr)
   CALL MPI_TYPE_COMMIT(type_recv1,ierr)   
   CALL MPI_TYPE_VECTOR(xnum,znum,znum+2,MPI_DOUBLE_PRECISION,type_up,ierr)
   CALL MPI_TYPE_COMMIT(type_up,ierr)
   CALL MPI_TYPE_VECTOR(xnum,znum,znum+2,MPI_DOUBLE_PRECISION,type_down,ierr)
   CALL MPI_TYPE_COMMIT(type_down,ierr)
   CALL MPI_TYPE_CONTIGUOUS(znum*xnum,MPI_DOUBLE_PRECISION,type_recv2,ierr)
   CALL MPI_TYPE_COMMIT(type_recv2,ierr)
   CALL MPI_TYPE_VECTOR(xnum,1,znum+2,MPI_DOUBLE_PRECISION,type_in_up_down,ierr)
   CALL MPI_TYPE_COMMIT(type_in_up_down,ierr)
   CALL MPI_TYPE_VECTOR(xnum,1,znum+2,MPI_DOUBLE_PRECISION,type_out_up_down,ierr)
   CALL MPI_TYPE_COMMIT(type_out_up_down,ierr)
   CALL MPI_TYPE_VECTOR(rnum,1,(xnum+2)*(znum+2),MPI_DOUBLE_PRECISION,type_left_in_out,ierr)
   CALL MPI_TYPE_COMMIT(type_left_in_out,ierr)
   CALL MPI_TYPE_VECTOR(rnum,1,(xnum+2)*(znum+2),MPI_DOUBLE_PRECISION,type_right_in_out,ierr)
   CALL MPI_TYPE_COMMIT(type_right_in_out,ierr)
END SUBROUTINE type_create_for_communication

SUBROUTINE boundary_communication(ierr)
   USE globalp
   IMPLICIT NONE
   include 'mpif.h'
   INTEGER :: i,j,k,count
   INTEGER :: ierr
   INTEGER :: status(mpi_status_size),position
   INTEGER :: srequest,rrequest1,rrequest2
   INTEGER :: array_of_status(mpi_status_size,8)
   INTEGER, DIMENSION(8) :: rrequest
      count=0    
      IF(direction_send(3).EQ.1.OR.send_direction(3).EQ.1)THEN        
         IF(direction_send(3).EQ.1)THEN
            position=0
            DO k=1,rnum
               DO j=1,xnum
                  CALL MPI_PACK(ttn(1,j,k),2,MPI_DOUBLE_PRECISION,LBUF1,outsize,position,comm_cart,ierr)
               ENDDO
            ENDDO
            CALL MPI_ISEND(LBUF1,outsize,MPI_PACKED,in,33,comm_cart,srequest,ierr)
         ENDIF
         IF(send_direction(3).EQ.1)THEN
            CALL MPI_IRECV(LBUF2,outsize,MPI_PACKED,in,111,comm_cart,rrequest1,ierr)
         ENDIF        
      ENDIF   
      IF(direction_send(4).EQ.1.OR.send_direction(4).EQ.1)THEN
         IF(direction_send(4).EQ.1)THEN
            position=0
            DO k=1,rnum
               DO j=1,xnum
                  CALL MPI_PACK(ttn(znum-1,j,k),2,MPI_DOUBLE_PRECISION,RBUF1,outsize,position,comm_cart,ierr)
               ENDDO
            ENDDO
            CALL MPI_ISEND(RBUF1,outsize,MPI_PACKED,out,111,comm_cart,srequest,ierr)
         ENDIF
         IF(send_direction(4).EQ.1)THEN
            CALL MPI_IRECV(RBUF2,outsize,MPI_PACKED,out,33,comm_cart,rrequest2,ierr)
         ENDIF 
      ENDIF    
      IF(send_direction(3).EQ.1)THEN
         CALL MPI_WAIT(rrequest1,status,ierr)
         position=0
         DO k=1,rnum
            DO j=1,xnum
               CALL MPI_UNPACK(LBUF2,outsize,position,ttn_ghost_in(1,j,k),2,MPI_DOUBLE_PRECISION,comm_cart,ierr)
            ENDDO
         ENDDO
      ENDIF
      IF(send_direction(4).EQ.1)THEN
         CALL MPI_WAIT(rrequest2,status,ierr)
         position=0
         DO k=1,rnum
            DO j=1,xnum
               CALL MPI_UNPACK(RBUF2,outsize,position,ttn_ghost_out(1,j,k),2,MPI_DOUBLE_PRECISION,comm_cart,ierr)
            ENDDO
         ENDDO
      ENDIF      
      IF(direction_send(1).EQ.1.OR.send_direction(1).EQ.1)THEN      
         IF(direction_send(1).EQ.1)THEN
            CALL MPI_ISEND(ttn(1,1,1),1,type_left,left,11,comm_cart,srequest,ierr)
            CALL MPI_ISEND(ttn(1,2,1),1,type_left,left,22,comm_cart,srequest,ierr)
         ENDIF
         IF(send_direction(1).EQ.1)THEN
            count=count+2
            CALL MPI_IRECV(ttn_ghost_left(1,1,1),1,type_recv1,left,11,comm_cart,rrequest(count-1),ierr)
            CALL MPI_IRECV(ttn_ghost_left(1,2,1),1,type_recv1,left,22,comm_cart,rrequest(count),ierr)
         ENDIF         
      ENDIF
      IF(direction_send(2).EQ.1.OR.send_direction(2).EQ.1)THEN
         IF(direction_send(2).EQ.1)THEN
            CALL MPI_ISEND(ttn(1,xnum-1,1),1,type_right,right,11,comm_cart,srequest,ierr)
            CALL MPI_ISEND(ttn(1,xnum,1),1,type_right,right,22,comm_cart,srequest,ierr)
         ENDIF
         IF(send_direction(2).EQ.1)THEN
            count=count+2
            CALL MPI_IRECV(ttn_ghost_right(1,1,1),1,type_recv1,right,11,comm_cart,rrequest(count-1),ierr)
            CALL MPI_IRECV(ttn_ghost_right(1,2,1),1,type_recv1,right,22,comm_cart,rrequest(count),ierr)
         ENDIF    
      ENDIF  
      IF(direction_send(5).EQ.1.OR.send_direction(5).EQ.1)THEN         
         IF(direction_send(5).EQ.1)THEN
            CALL MPI_ISEND(ttn(1,1,1),1,type_up,up,11,comm_cart,srequest,ierr)
            CALL MPI_ISEND(ttn(1,1,2),1,type_up,up,22,comm_cart,srequest,ierr)
         ENDIF       
         IF(send_direction(5).EQ.1)THEN
            count=count+2
            CALL MPI_IRECV(ttn_ghost_up(1,1,1),1,type_recv2,up,11,comm_cart,rrequest(count-1),ierr)
            CALL MPI_IRECV(ttn_ghost_up(1,1,2),1,type_recv2,up,22,comm_cart,rrequest(count),ierr)
         ENDIF
      ENDIF   
      IF(direction_send(6).EQ.1.OR.send_direction(6).EQ.1)THEN
         IF(direction_send(6).EQ.1)THEN
            CALL MPI_ISEND(ttn(1,1,rnum-1),1,type_down,down,11,comm_cart,srequest,ierr)
            CALL MPI_ISEND(ttn(1,1,rnum),1,type_down,down,22,comm_cart,srequest,ierr)
         ENDIF
         IF(send_direction(6).EQ.1)THEN
            count=count+2
            CALL MPI_IRECV(ttn_ghost_down(1,1,1),1,type_recv2,down,11,comm_cart,rrequest(count-1),ierr)
            CALL MPI_IRECV(ttn_ghost_down(1,1,2),1,type_recv2,down,22,comm_cart,rrequest(count),ierr)
         ENDIF
      ENDIF
      IF(count.GT.0)THEN
         CALL MPI_WAITALL(count,rrequest(1:count),array_of_status,ierr)
      ENDIF   
END SUBROUTINE boundary_communication

SUBROUTINE cross_data_communication(ierr)
   USE globalp
   IMPLICIT NONE
   include 'mpif.h'
   INTEGER :: i,j,k,count
   INTEGER :: ierr
   INTEGER :: srequest
   INTEGER :: array_of_status(mpi_status_size,12)
   INTEGER, DIMENSION(12) :: rrequest
   INTEGER :: left_up,right_up,in_up,out_up,left_down,right_down,in_down,out_down
   INTEGER :: left_in,left_out,right_in,right_out
   count=0
   IF(rank_r.NE.0)THEN
      IF(rank_x.NE.0)THEN         
         count=count+1
         left_up = myrank - nx*nz - nz
         CALL MPI_IRECV(ttn(1,0,0),znum,MPI_DOUBLE_PRECISION,left_up,11,comm_cart,rrequest(count),ierr)
         CALL MPI_ISEND(ttn(1,1,1),znum,MPI_DOUBLE_PRECISION,left_up,12,comm_cart,srequest,ierr)
      ENDIF
      IF(rank_x.NE.nx-1)THEN        
         count=count+1
         right_up = myrank - nx*nz + nz
         CALL MPI_IRECV(ttn(1,xnum+1,0),znum,MPI_DOUBLE_PRECISION,right_up,13,comm_cart,rrequest(count),ierr)
         CALL MPI_ISEND(ttn(1,xnum,1),znum,MPI_DOUBLE_PRECISION,right_up,14,comm_cart,srequest,ierr)
      ENDIF
      IF(rank_z.NE.0)THEN       
         count=count+1
         in_up = myrank - nx*nz - 1
         CALL MPI_IRECV(ttn(0,1,0),1,type_in_up_down,in_up,15,comm_cart,rrequest(count),ierr)
         CALL MPI_ISEND(ttn(1,1,1),1,type_in_up_down,in_up,16,comm_cart,srequest,ierr)      
      ENDIF
      IF(rank_z.NE.nz-1)THEN        
         count=count+1
         out_up = myrank - nx*nz + 1
         CALL MPI_IRECV(ttn(znum+1,1,0),1,type_out_up_down,out_up,17,comm_cart,rrequest(count),ierr)
         CALL MPI_ISEND(ttn(znum,1,1),1,type_out_up_down,out_up,18,comm_cart,srequest,ierr)
      ENDIF
   ENDIF
   IF(rank_r.NE.nr-1)THEN
      IF(rank_x.NE.0)THEN         
         count=count+1
         left_down = myrank + nx*nz -nz
         CALL MPI_IRECV(ttn(1,0,rnum+1),znum,MPI_DOUBLE_PRECISION,left_down,14,comm_cart,rrequest(count),ierr)
         CALL MPI_ISEND(ttn(1,1,rnum),znum,MPI_DOUBLE_PRECISION,left_down,13,comm_cart,srequest,ierr)      
      ENDIF
      IF(rank_x.NE.nx-1)THEN     
         count=count+1
         right_down = myrank + nx*nz + nz
         CALL MPI_IRECV(ttn(1,xnum+1,rnum+1),znum,MPI_DOUBLE_PRECISION,right_down,12,comm_cart,rrequest(count),ierr)
         CALL MPI_ISEND(ttn(1,xnum,rnum),znum,MPI_DOUBLE_PRECISION,right_down,11,comm_cart,srequest,ierr)
      ENDIF
      IF(rank_z.NE.0)THEN         
         count=count+1
         in_down = myrank + nx*nz - 1
         CALL MPI_IRECV(ttn(0,1,rnum+1),1,type_in_up_down,in_down,18,comm_cart,rrequest(count),ierr)
         CALL MPI_ISEND(ttn(1,1,rnum),1,type_in_up_down,in_down,17,comm_cart,srequest,ierr)      
      ENDIF      
      IF(rank_z.NE.nz-1)THEN        
         count=count+1
         out_down = myrank + nx*nz + 1
         CALL MPI_IRECV(ttn(znum+1,1,rnum+1),1,type_out_up_down,out_down,16,comm_cart,rrequest(count),ierr)
         CALL MPI_ISEND(ttn(znum,1,rnum),1,type_out_up_down,out_down,15,comm_cart,srequest,ierr)
      ENDIF
   ENDIF    
   IF(rank_x.NE.0)THEN
      IF(rank_z.NE.0)THEN
         count=count+1
         left_in = myrank - nz - 1
         CALL MPI_IRECV(ttn(0,0,1),1,type_left_in_out,left_in,19,comm_cart,rrequest(count),ierr)
         CALL MPI_ISEND(ttn(1,1,1),1,type_left_in_out,left_in,20,comm_cart,srequest,ierr)
      ENDIF
      IF(rank_z.NE.nz-1)THEN
         count=count+1
         left_out = myrank - nz + 1
         CALL MPI_IRECV(ttn(znum+1,0,1),1,type_left_in_out,left_out,21,comm_cart,rrequest(count),ierr)
         CALL MPI_ISEND(ttn(znum,1,1),1,type_left_in_out,left_out,22,comm_cart,srequest,ierr)
      ENDIF
   ENDIF
   IF(rank_x.NE.nx-1)THEN
      IF(rank_z.NE.0)THEN
         count=count+1
         right_in = myrank + nz - 1
         CALL MPI_IRECV(ttn(0,xnum+1,1),1,type_right_in_out,right_in,22,comm_cart,rrequest(count),ierr)
         CALL MPI_ISEND(ttn(1,xnum,1),1,type_right_in_out,right_in,21,comm_cart,srequest,ierr)
      ENDIF
      IF(rank_z.NE.nz-1)THEN
         count=count+1
         right_out = myrank + nz + 1
         CALL MPI_IRECV(ttn(znum+1,xnum+1,1),1,type_right_in_out,right_out,20,comm_cart,rrequest(count),ierr)
         CALL MPI_ISEND(ttn(znum,xnum,1),1,type_right_in_out,right_out,19,comm_cart,srequest,ierr)
      ENDIF
   ENDIF
   IF(count.GT.0)THEN
      CALL MPI_WAITALL(count,rrequest(1:count),array_of_status,ierr)
   ENDIF
END SUBROUTINE cross_data_communication

SUBROUTINE start_fmm(myfile,srcid,ierr)
   USE globalp
   USE traveltime
   IMPLICIT NONE
   include 'mpif.h'
   INTEGER :: i,j,k,myfile,ierr
   INTEGER :: sw,srcid,ittmin,jttmin
   REAL(KIND=i10) :: rd1,ttmin
   REAL(KIND=i10), DIMENSION (:), ALLOCATABLE :: rcounts_ttnmin
   REAL(KIND=i10), DIMENSION (:,:), ALLOCATABLE :: ttn_start
   ALLOCATE(ttn_start(znum,xnum), STAT=checkstat)
   IF(checkstat > 0)THEN
      WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin3d: REAL ttn'
   ENDIF   
   CALL MPI_FILE_READ_ALL(myfile,ttn_start,znum*xnum,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
   DO i=1,xnum
      DO j=1,znum
         ttn(j,i,rnum)=ttn_start(j,i)
      ENDDO
   ENDDO 
   DEALLOCATE(ttn_start, STAT=checkstat)
   IF(checkstat > 0)THEN
      WRITE(6,*)'Error with DEALLOCATE: SUBROUTINE start_fmm: ttn_start'
   ENDIF
   sw=0
   ttmin = -100.0
   DO i=1,xnum
      DO j=1,znum
         rd1 = ttn(j,i,rnum)
         IF(rd1.GT.0.0)THEN
            IF(sw.eq.0)THEN
               ttmin=rd1
               ittmin=i
               jttmin=j
               sw=1
            ELSE IF(rd1.LT.ttmin)THEN
               ttmin=rd1
               ittmin=i
               jttmin=j
            ENDIF
            nsts(j,i,rnum)=0
         ELSE
            ttn(j,i,rnum)=0.0
         ENDIF
      ENDDO
   ENDDO
   ALLOCATE(rcounts_ttnmin(nx*nz))
   CALL MPI_ALLGATHER(ttmin,1,MPI_DOUBLE_PRECISION,rcounts_ttnmin,1,MPI_DOUBLE_PRECISION,comm_start,ierr)
   sw=0
   IF(ttmin.EQ.-100)THEN
      sw=1
   ELSE
      DO i=1,nx*nz
         IF(rcounts_ttnmin(i).NE.-100)THEN   
            IF(ttmin.GT.rcounts_ttnmin(i))THEN
               sw=1
               EXIT 
            ENDIF
         ENDIF
      ENDDO 
   ENDIF
   DEALLOCATE(rcounts_ttnmin)
   IF(sw.EQ.0)THEN
      CALL addtree_all(jttmin,ittmin,rnum)      
      CALL local_fmm 
      sw_fmm=1       
   ENDIF
END SUBROUTINE start_fmm

SUBROUTINE error_calculate(ierr)
   USE globalp
   IMPLICIT NONE
   include 'mpif.h'
   INTEGER :: i,j,k,ierr
   REAL(KIND=i10) :: temp,tempp
   REAL(KIND=i10) :: dtr,dtx,dtz,delta
   temp=0.0
   tempp=0.0
   error_to_stop=0.0
   DO i=1,rnum
      DO j=1,xnum
         DO k=1,znum
            IF(temp.EQ.0)THEN       
               temp=abs(ttn(k,j,i)-ttn_previous_step(k,j,i))
            ELSE
               IF(temp.LT.abs(ttn(k,j,i)-ttn_previous_step(k,j,i)))then
                  temp=abs(ttn(k,j,i)-ttn_previous_step(k,j,i))
               ENDIF
            ENDIF
         ENDDO
      ENDDO
   ENDDO   
   CALL MPI_ALLREDUCE(temp,error_to_stop,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm_cart,ierr)
END SUBROUTINE error_calculate

SUBROUTINE gridder(gridv)
   USE globalp
   IMPLICIT NONE
   INTEGER :: i,j,k,l,m,n,i1,j1,k1,ri,xi,zi
   INTEGER :: conr,conx,conz,str,stx,stz
   REAL(KIND=i10) :: u,v,w,sumi,sumj,sumk,distz,distx,distr
   REAL(KIND=i10), DIMENSION (4) :: ui,vi,wi
   CHARACTER (LEN=20) :: gridv
   ALLOCATE(veln(-1:znum+1,-1:xnum+1,-1:rnum+1))
   DO i=-1,znum+1
      zi=FLOOR((rgoz+(zoffset-1)*dnz+(i-1)*dnz-goz)/dvz)+1
      distz=rgoz+(zoffset-1)*dnz+(i-1)*dnz-goz-(zi-1)*dvz
      w=distz/dvz
      wi(1)=(1.0-w)**3/6.0
      wi(2)=(4.0-6.0*w**2+3.0*w**3)/6.0
      wi(3)=(1.0+3.0*w+3.0*w**2-3.0*w**3)/6.0
      wi(4)=w**3/6.0      
      DO j=-1,xnum+1
         xi=FLOOR((rgox+(xoffset-1)*dnx+(j-1)*dnx-gox)/dvx)+1
         distx=rgox+(xoffset-1)*dnx+(j-1)*dnx-gox-(xi-1)*dvx
         v=distx/dvx
         vi(1)=(1.0-v)**3/6.0
         vi(2)=(4.0-6.0*v**2+3.0*v**3)/6.0
         vi(3)=(1.0+3.0*v+3.0*v**2-3.0*v**3)/6.0
         vi(4)=v**3/6.0      
         DO k=-1,rnum+1
            ri=FLOOR((gor-rgor+(roffset-1)*dnr+(k-1)*dnr)/dvr)+1
            distr=-(rgor-(roffset-1)*dnr-(k-1)*dnr)+(gor-(ri-1)*dvr)
            u=distr/dvr
            ui(1)=(1.0-u)**3/6.0
            ui(2)=(4.0-6.0*u**2+3.0*u**3)/6.0
            ui(3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
            ui(4)=u**3/6.0               
            sumi=0.0
            DO i1=1,4                  
               sumj=0.0 
               DO j1=1,4                     
                  sumk=0.0	
                  DO k1=1,4
                     sumk=sumk+ui(k1)*velv(ri-2+k1,xi-2+j1,zi-2+i1)
                  ENDDO
                  sumj=sumj+vi(j1)*sumk
               ENDDO
               sumi=sumi+wi(i1)*sumj
            ENDDO
            veln(i,j,k)=sumi               
         ENDDO
      ENDDO
   ENDDO 
END SUBROUTINE gridder

SUBROUTINE srtimes(srcid,caid,nra,fsrt,nsrc,rtravel,ierr)
   USE globalp
   IMPLICIT NONE
   include 'mpif.h'
   INTEGER :: i,j,k,l,m,irr,irx,irz,sw,temp,srcid,caid,nra,fsrt,nsrc
   INTEGER :: temp_myrank,temp_nprocs,ierr
   INTEGER :: group,newgroup
   REAL :: sww
   REAL(KIND=i10) :: drr,drx,drz,produ
   CHARACTER (LEN=20) :: rtravel
   INTEGER, DIMENSION (:), ALLOCATABLE :: RANKS
   INTEGER, DIMENSION (:,:), ALLOCATABLE :: TEMPP   
   IF(srcid.EQ.1)THEN
      CALL MPI_COMM_SIZE(comm_srtimes,temp_nprocs,ierr)
      ALLOCATE(rcounts_srtimes(temp_nprocs),TEMPP(temp_nprocs,2),STAT=checkstat)
      IF(checkstat > 0)THEN
         WRITE(6,*)'Error with ALLOCATE: PROGRAM srtimes: INTEGER rcounts_srtimes,TEMPP'
      ENDIF      
      temp = nsrc*nrc/6
      ALLOCATE(trr(temp),STAT=checkstat)
      IF(checkstat > 0)THEN
         WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin3d: REAL trr'
      ENDIF
      trr=0.0
      number_of_trr=0
      DO i=1,nrc      
         irr=FLOOR((rgor-(roffset-1)*dnr-rcr(i))/dnr)+1
         irx=FLOOR((rcx(i)-rgox-(xoffset-1)*dnx)/dnx)+1
         irz=FLOOR((rcz(i)-rgoz-(zoffset-1)*dnz)/dnz)+1         
         IF(rank_x.EQ.nx-1.AND.irx.EQ.xnum)irx=irx-1
         IF(rank_z.EQ.nz-1.AND.irz.EQ.znum)irz=irz-1         
         sw=0
         IF(irr.LT.1.OR.irr.GE.rnum)sw=1   
         IF(rank_x.EQ.0)THEN
            IF(irx.LT.1.OR.irx.GE.xnum)THEN
               sw=1
            ENDIF
         ELSE
            IF(irx.LT.0.OR.irx.GE.xnum)THEN
               sw=1
            ENDIF
         ENDIF   
         IF(rank_z.EQ.0)THEN
            IF(irz.LT.1.OR.irz.GE.znum)THEN
               sw=1
            ENDIF
         ELSE
            IF(irz.LT.0.OR.irz.GE.znum)THEN
               sw=1
            ENDIF
         ENDIF
         IF(sw.EQ.0)THEN           
            number_of_trr = number_of_trr + 1           
            rec_location(:,number_of_trr)=[irr,irx,irz]
            number(number_of_trr)=i
            drx=(rcx(i)-rgox-(xoffset-1)*dnx)-(irx-1)*dnx
            drz=(rcz(i)-rgoz-(zoffset-1)*dnz)-(irz-1)*dnz
            drr=(rgor-(roffset-1)*dnr-rcr(i))-(irr-1)*dnr           
            DO j=1,2
               DO k=1,2
                  DO l=1,2         
                     produ=(1.0-ABS(((l-1)*dnz-drz)/dnz))*(1.0-ABS(((k-1)*dnx-drx)/dnx))
                     produ=produ*(1.0-ABS(((j-1)*dnr-drr)/dnr))
                     trr(number_of_trr)=trr(number_of_trr)+ttn(irz-1+l,irx-1+k,irr-1+j)*produ
                  ENDDO
               ENDDO
            ENDDO            
         ENDIF
      ENDDO
      CALL MPI_ALLGATHER(number_of_trr,1,MPI_INTEGER,rcounts_srtimes,1,MPI_INTEGER,comm_srtimes,ierr)
      CALL MPI_BARRIER(comm_srtimes,ierr)
      temp=0
      DO i=1,temp_nprocs
         IF(rcounts_srtimes(i).GT.0)THEN
            temp=temp+1
            TEMPP(temp,1)=i-1
            TEMPP(temp,2)=rcounts_srtimes(i)
         ENDIF
      ENDDO
      DEALLOCATE(rcounts_srtimes)     
      ALLOCATE(RANKS(temp),STAT=checkstat)
      IF(checkstat > 0)THEN
         WRITE(6,*)'Error with ALLOCATE: PROGRAM srtimes: INTEGER RANKS'
      ENDIF
      CALL MPI_BARRIER(comm_srtimes,ierr)
      CALL MPI_COMM_GROUP(comm_srtimes,group,ierr)
      RANKS=TEMPP(1:temp,1)
      CALL MPI_GROUP_INCL(group,temp,RANKS,newgroup,ierr)
      comm_srtimes_to_file=MPI_COMM_NULL
      CALL MPI_COMM_CREATE(comm_srtimes,newgroup,comm_srtimes_to_file,ierr)     
      IF(comm_srtimes_to_file.NE.MPI_COMM_NULL)THEN
         CALL MPI_COMM_RANK(comm_srtimes_to_file,srtimes_myrank,ierr)
         CALL MPI_COMM_SIZE(comm_srtimes_to_file,srtimes_nprocs,ierr)
         IF(srtimes_myrank.EQ.0)THEN
            IF(fsrt.EQ.1)THEN
               IF(caid==1)THEN
                  OPEN(UNIT=10,FILE=rtravel,STATUS='unknown')
               ELSE
                  OPEN(UNIT=10,FILE=rtravel,ACCESS='append')
               ENDIF
            ENDIF
            ALLOCATE(rcounts_srtimes(srtimes_nprocs),displs_srtimes(srtimes_nprocs))
            ALLOCATE(trr_file(nrc),number_file(nrc),temp_number_file(nrc),STAT=checkstat)
            IF(checkstat > 0)THEN
               WRITE(6,*)'Error with ALLOCATE: PROGRAM srtimes: INTEGER rcounts_srtimes,displs_srtimes,trr_file,number_file'
            ENDIF
            rcounts_srtimes=TEMPP(1:srtimes_nprocs,2)
            displs_srtimes=0
            DO i=2,srtimes_nprocs
               displs_srtimes(i)=displs_srtimes(i-1)+rcounts_srtimes(i-1)
            ENDDO
         ENDIF
         CALL MPI_GATHERV(trr,number_of_trr,MPI_FLOAT,trr_file,& 
                        & rcounts_srtimes,displs_srtimes,MPI_FLOAT,0,comm_srtimes_to_file,ierr)
         CALL MPI_GATHERV(number,number_of_trr,MPI_INTEGER,number_file,& 
                        & rcounts_srtimes,displs_srtimes,MPI_INTEGER,0,comm_srtimes_to_file,ierr)              
         CALL MPI_BARRIER(comm_srtimes_to_file,ierr)
         IF(srtimes_myrank.EQ.0)THEN
            temp_number_file=number_file
            temp=1
            DO WHILE(temp.NE.0)
               temp=0
               DO i=1,nrc
                  IF(temp_number_file(i).NE.i)THEN                     
                     temp=1
                     sww=trr_file(temp_number_file(i))
                     trr_file(temp_number_file(i))=trr_file(i)
                     trr_file(i)=sww
                     sw=temp_number_file(temp_number_file(i))
                     temp_number_file(temp_number_file(i))=temp_number_file(i)
                     temp_number_file(i)=sw                     
                  ENDIF
               ENDDO
            ENDDO
            DO j=1,nrc
               WRITE(10,*)trr_file(j)
            ENDDO
         ENDIF         
      ENDIF     
      DEALLOCATE(trr,TEMPP,RANKS)     
   ELSE
      IF(number_of_trr.GT.0)THEN         
         IF(srcid.EQ.2)THEN
            ALLOCATE(trr(number_of_trr),STAT=checkstat)
            IF(checkstat > 0)THEN
               WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin3d: REAL trr'
            ENDIF 
         ENDIF
         trr=0.0
         DO m=1,number_of_trr
            i=number(m)
            irr=rec_location(1,m)
            irx=rec_location(2,m)
            irz=rec_location(3,m)
            drx=(rcx(i)-rgox-(xoffset-1)*dnx)-(irx-1)*dnx
            drz=(rcz(i)-rgoz-(zoffset-1)*dnz)-(irz-1)*dnz
            drr=(rgor-(roffset-1)*dnr-rcr(i))-(irr-1)*dnr
            DO j=1,2
               DO k=1,2
                  DO l=1,2
                     produ=(1.0-ABS(((l-1)*dnz-drz)/dnz))*(1.0-ABS(((k-1)*dnx-drx)/dnx))
                     produ=produ*(1.0-ABS(((j-1)*dnr-drr)/dnr))
                     trr(m)=trr(m)+ttn(irz-1+l,irx-1+k,irr-1+j)*produ                    
                  ENDDO
               ENDDO
            ENDDO           
         ENDDO         
         IF(comm_srtimes_to_file.NE.MPI_COMM_NULL)THEN
            CALL MPI_GATHERV(trr,number_of_trr,MPI_FLOAT,trr_file,&
            & rcounts_srtimes,displs_srtimes,MPI_FLOAT,0,comm_srtimes_to_file,ierr)            
            CALL MPI_BARRIER(comm_srtimes_to_file,ierr)
            IF(srtimes_myrank.EQ.0)THEN
               temp_number_file=number_file
               temp=1
               DO WHILE(temp.NE.0)
                  temp=0
                  DO i=1,nrc
                     IF(temp_number_file(i).NE.i)THEN
                        temp=1
                        sww=trr_file(temp_number_file(i))
                        trr_file(temp_number_file(i))=trr_file(i)
                        trr_file(i)=sww      
                        sw=temp_number_file(temp_number_file(i))
                        temp_number_file(temp_number_file(i))=temp_number_file(i)
                        temp_number_file(i)=sw                           
                     ENDIF
                  ENDDO
               ENDDO
               DO j=1,nrc
                  WRITE(10,*)trr_file(j)
               ENDDO                  
            ENDIF
         ENDIF
      ENDIF
   ENDIF
   IF(srcid.EQ.nsrc)THEN
      IF(comm_srtimes_to_file.NE.MPI_COMM_NULL)THEN
         DEALLOCATE(trr)
         IF(srtimes_myrank.EQ.0)THEN
            DEALLOCATE(rcounts_srtimes,displs_srtimes)
            DEALLOCATE(trr_file,number_file,temp_number_file)
            IF(fsrt.EQ.1)THEN
               CLOSE(10)
            ENDIF
         ENDIF
         CALL MPI_COMM_FREE(comm_srtimes_to_file,ierr)
      ENDIF
   ENDIF
END SUBROUTINE srtimes

SUBROUTINE parallel_rpaths(csid,caid,cfd,count,ierr)
   USE globalp
   IMPLICIT NONE
   include 'mpif.h'
   INTEGER :: i,j,k,l,m,n,ipr,ipx,ipz,nrp
   INTEGER :: sw,swr,swx,swz
   INTEGER :: cfd,csid,ipro,ipxo,ipzo,caid
   INTEGER :: ivr,ivx,ivz,ivro,ivxo,ivzo,nhp
   INTEGER :: ivrt,ivxt,ivzt,iprt,ipxt,ipzt,isum
   INTEGER :: temp,tempp
   INTEGER :: ierr,count
   INTEGER, DIMENSION (4) :: chp
   INTEGER, PARAMETER :: maxrp=100000
   REAL(KIND=i5), PARAMETER :: ftol=1.0e-6
   REAL(KIND=i5) :: rayr,rayx,rayz
   REAL(KIND=i10) :: dpl,crad,rd1,rd2,atio,ri,xi,zi,vel,velo
   REAL(KIND=i10) :: u,v,w,rigz,rigx,rigr,dinc
   REAL(KIND=i10) :: dtr,dtx,dtz,drr,drx,drz,produ
   REAL(KIND=i10) :: dtr1,dtr2,dtr3,dtr4
   REAL(KIND=i10), DIMENSION (:), ALLOCATABLE :: rgr,rgx,rgz
   REAL(KIND=i5), DIMENSION (:,:,:), ALLOCATABLE :: fdm
   REAL(KIND=i10), DIMENSION (4) :: vrat,ui,vi,wi,uio,vio,wio
   ALLOCATE(rgr(xnum*rnum), STAT=checkstat)
   IF(checkstat > 0)THEN
      WRITE(6,*)'Error with ALLOCATE: SUBROUTINE parallel_rpaths: REAL rgr'
   ENDIF
   ALLOCATE(rgx(xnum*rnum), STAT=checkstat)
   IF(checkstat > 0)THEN
      WRITE(6,*)'Error with ALLOCATE: SUBROUTINE parallel_rpaths: REAL rgx'
   ENDIF
   ALLOCATE(rgz(xnum*rnum), STAT=checkstat)
   IF(checkstat > 0)THEN
      WRITE(6,*)'Error with ALLOCATE: SUBROUTINE parallel_rpaths: REAL rgz'
   ENDIF
   IF(cfd.EQ.1)THEN
      ALLOCATE(fdm(0:nvz+1,0:nvx+1,0:nvr+1), STAT=checkstat)
      IF(checkstat > 0)THEN
         WRITE(6,*)'Error with ALLOCATE: SUBROUTINE parallel_rpaths: REAL fdm'
      ENDIF
   ENDIF 
   dpl=dnr
   DO i=1,communicate_rpaths_length
      IF(cfd.EQ.1)THEN
         fdm=0.0
      ENDIF
      ipr=FLOOR((rgor-(roffset-1)*dnr-communicate_rpaths(1,i))/dnr)+1
      ipx=FLOOR((communicate_rpaths(2,i)-rgox-(xoffset-1)*dnx)/dnx)+1
      ipz=FLOOR((communicate_rpaths(3,i)-rgoz-(zoffset-1)*dnz)/dnz)+1
      sw=0
      IF(rank_r.EQ.0)THEN
         IF(ipr.LT.1.OR.ipr.GE.rnum)THEN
            sw=1
         ENDIF
      ELSE
         IF(ipr.LT.0.OR.ipr.GE.rnum)THEN
            sw=1
         ENDIF
      ENDIF
      IF(rank_x.EQ.0)THEN
         IF(ipx.LT.1.OR.ipx.GE.xnum)THEN
            sw=1
         ENDIF
      ELSE
         IF(ipx.LT.0.OR.ipx.GE.xnum)THEN
            sw=1
         ENDIF
      ENDIF
      IF(rank_z.EQ.0)THEN
         IF(ipz.LT.1.OR.ipz.GE.znum)THEN
            sw=1
         ENDIF
      ELSE
         IF(ipz.LT.0.OR.ipz.GE.znum)THEN
            sw=1
         ENDIF
      ENDIF
      IF(count.EQ.0)THEN
         IF(ipr.EQ.rnum.AND.rank_r.EQ.nr-1)THEN
            ipr=ipr-1
            sw=0
         ENDIF
         IF(ipx.EQ.xnum.AND.rank_x.EQ.nx-1)THEN
            ipx=ipx-1
            sw=0
         ENDIF
         IF(ipz.EQ.znum.AND.rank_z.EQ.nz-1)THEN
            ipz=ipz-1
            sw=0
         ENDIF
      ENDIF    
      IF(sw.NE.0)THEN
         communicate_rpaths(1:3,i)=[0.0,0.0,0.0]
         CYCLE
      ENDIF
      rgr(1)=communicate_rpaths(1,i)
      rgx(1)=communicate_rpaths(2,i)
      rgz(1)=communicate_rpaths(3,i)
      DO j=1,maxrp
         dtr=0.0
         temp=0
         IF(ttn(ipz,ipx,ipr+1).LT.ttn(ipz,ipx,ipr))THEN
            dtr=dtr+ttn(ipz,ipx,ipr+1)-ttn(ipz,ipx,ipr)
            temp=temp+1
         ENDIF
         IF(ttn(ipz+1,ipx,ipr+1).LT.ttn(ipz+1,ipx,ipr))THEN
            dtr=dtr+ttn(ipz+1,ipx,ipr+1)-ttn(ipz+1,ipx,ipr)
            temp=temp+1
         ENDIF
         IF(ttn(ipz+1,ipx+1,ipr+1).LT.ttn(ipz+1,ipx+1,ipr))THEN
            dtr=dtr+ttn(ipz+1,ipx+1,ipr+1)-ttn(ipz+1,ipx+1,ipr)
            temp=temp+1
         ENDIF
         IF(ttn(ipz,ipx+1,ipr+1).LT.ttn(ipz,ipx+1,ipr))THEN
            dtr=dtr+ttn(ipz,ipx+1,ipr+1)-ttn(ipz,ipx+1,ipr)
            temp=temp+1
         ENDIF
         IF(temp.NE.0)THEN
            dtr=dtr/(temp*dnr)
         ENDIF
         dtx=ttn(ipz,ipx+1,ipr)-ttn(ipz,ipx,ipr)
         dtx=dtx+ttn(ipz+1,ipx+1,ipr)-ttn(ipz+1,ipx,ipr)
         dtx=dtx+ttn(ipz+1,ipx+1,ipr+1)-ttn(ipz+1,ipx,ipr+1)
         dtx=dtx+ttn(ipz,ipx+1,ipr+1)-ttn(ipz,ipx,ipr+1)
         dtx=dtx/(4.0*(earth+rgor-(roffset-1)*dnr-(ipr-1)*dnr)*dnx)
         dtz=ttn(ipz+1,ipx,ipr)-ttn(ipz,ipx,ipr)
         dtz=dtz+ttn(ipz+1,ipx+1,ipr+1)-ttn(ipz,ipx+1,ipr+1)
         dtz=dtz+ttn(ipz+1,ipx+1,ipr)-ttn(ipz,ipx+1,ipr)
         dtz=dtz+ttn(ipz+1,ipx,ipr+1)-ttn(ipz,ipx,ipr+1)
         dtz=dtz/(4.0*(earth+rgor-(roffset-1)*dnr-(ipr-1)*dnr)*SIN(rgx(j))*dnz)
         crad=earth+rgr(j)
         rd1=SQRT(dtr**2+dtx**2+dtz**2)
         rgr(j+1)=rgr(j)+dpl*dtr/rd1
         rgx(j+1)=rgx(j)-dpl*dtx/(crad*rd1)
         rgz(j+1)=rgz(j)-dpl*dtz/(crad*SIN(rgx(j))*rd1)        
         ipro=ipr
         ipxo=ipx
         ipzo=ipz
         ipr=FLOOR((rgor-(roffset-1)*dnr-rgr(j+1))/dnr)+1
         ipx=FLOOR((rgx(j+1)-rgox-(xoffset-1)*dnx)/dnx)+1
         ipz=FLOOR((rgz(j+1)-rgoz-(zoffset-1)*dnz)/dnz)+1
         sw=0
         swr=0
         swx=0
         swz=0      
         IF(rank_r.EQ.0)THEN
            IF(ipr.LT.1)THEN
               sw=1
               swr=1
            ELSEIF(ipr.GE.rnum)THEN
               sw=-1
            ENDIF
         ELSEIF(rank_r.EQ.nr-1)THEN
            IF(ipr.LT.0)THEN
               sw=-1
            ELSEIF(ipr.GE.rnum)THEN
               sw=1
               swr=1
            ENDIF
         ELSE 
            IF(ipr.LT.0.OR.ipr.GE.rnum)THEN
               sw=-1
            ENDIF
         ENDIF
         IF(rank_x.EQ.0)THEN
            IF(ipx.LT.1)THEN
               sw=1
               swx=1
            ELSEIF(ipx.GE.xnum)THEN
               IF(sw.NE.1)THEN
                  sw=-1
               ENDIF
            ENDIF
         ELSEIF(rank_x.EQ.nx-1)THEN
            IF(ipx.LT.0)THEN
               IF(sw.NE.1)THEN
                  sw=-1
               ENDIF
            ELSEIF(ipx.GE.xnum)THEN
               sw=1
               swx=1
            ENDIF
         ELSE
            IF(ipx.LT.0.OR.ipx.GE.xnum)THEN
               IF(sw.NE.1)THEN
                  sw=-1
               ENDIF
            ENDIF
         ENDIF
         IF(rank_z.EQ.0)THEN
            IF(ipz.LT.1)THEN
               sw=1
               swz=1
            ELSEIF(ipz.GE.znum)THEN
               IF(sw.NE.1)THEN
                  sw=-1
               ENDIF
            ENDIF
         ELSEIF(rank_z.EQ.nz-1)THEN
            IF(ipz.LT.0)THEN
               IF(sw.NE.1)THEN
                  sw=-1
               ENDIF
            ELSEIF(ipz.GE.znum)THEN
               sw=1
               swz=1
            ENDIF
         ELSE
            IF(ipz.LT.0.OR.ipz.GE.znum)THEN
               IF(sw.NE.1)THEN
                  sw=-1
               ENDIF
            ENDIF
         ENDIF
         IF(sw.EQ.1)THEN
            sw=0
            communicate_rpaths(1:3,i)=[0.0,0.0,0.0]
            IF(swr.EQ.1)THEN
               IF(ipr.LT.1.OR.ipr.GE.rnum)THEN
                  IF(ipr.LT.1)THEN
                     ri=rgor-(roffset-1)*dnr
                     ipr=1
                  ELSE
                     ri=rgor-(roffset-1)*dnr-(rnum-1)*dnr
                     ipr=rnum-1
                  ENDIF
                  atio=(ri-rgr(j))/(rgr(j+1)-rgr(j))
                  sw=1
               ENDIF
            ENDIF            
            IF(swx.EQ.1)THEN
               IF(ipx.LT.1.OR.ipx.GE.xnum)THEN
                  IF(ipx.LT.1)THEN
                     xi=rgox+(xoffset-1)*dnx
                     ipx=1
                  ELSE
                     xi=rgox+(xoffset-1)*dnx+(xnum-1)*dnx
                     ipx=xnum-1
                  ENDIF
                  rd1=(xi-rgx(j))/(rgx(j+1)-rgx(j))
                  IF(sw.eq.1)then
                     IF(rd1.LT.atio)THEN
                        atio=rd1
                        sw=2
                     ENDIF
                  ELSE
                     atio=rd1
                     sw=2
                  ENDIF
               ENDIF
            ENDIF              
            IF(swz.EQ.1)THEN
               IF(ipz.LT.1.OR.ipz.GE.znum)THEN
                  IF(ipz.LT.1)THEN
                     zi=rgoz+(zoffset-1)*dnz
                     ipz=1
                  ELSE
                     zi=rgoz+(zoffset-1)*dnz+(znum-1)*dnz
                     ipz=znum-1
                  ENDIF
                  rd1=(zi-rgz(j))/(rgz(j+1)-rgz(j))
                  IF(sw.NE.0)then
                     IF(rd1.LT.atio)THEN
                        atio=rd1
                        sw=3
                     ENDIF
                  ELSE
                     atio=rd1
                     sw=3
                  ENDIF
               ENDIF
            ENDIF
            IF(sw.EQ.1)THEN
               ipx=ipxo
               ipz=ipzo
            ELSE IF(sw.EQ.2)THEN
               ipr=ipro
               ipz=ipzo
            ELSE IF(sw.EQ.3)THEN
               ipr=ipro
               ipx=ipxo
            ENDIF         
            rgz(j+1)=rgz(j)+atio*(rgz(j+1)-rgz(j))
            rgx(j+1)=rgx(j)+atio*(rgx(j+1)-rgx(j))
            rgr(j+1)=rgr(j)+atio*(rgr(j+1)-rgr(j))         
            IF(cfd.EQ.1)THEN
               sw=1
            ELSE
               EXIT
            ENDIF
         ELSEIF(sw.EQ.-1)THEN
            communicate_rpaths(1:3,i)=[rgr(j+1),rgx(j+1),rgz(j+1)]
         ENDIF 
         IF(cfd.EQ.1)THEN
            ivr=FLOOR((gor-rgr(j+1))/dvr)+1
            ivx=FLOOR((rgx(j+1)-gox)/dvx)+1
            ivz=FLOOR((rgz(j+1)-goz)/dvz)+1
            ivro=FLOOR((gor-rgr(j))/dvr)+1
            ivxo=FLOOR((rgx(j)-gox)/dvx)+1
            ivzo=FLOOR((rgz(j)-goz)/dvz)+1
            nhp=0                       
            IF(ivr.NE.ivro)THEN
               nhp=nhp+1
               IF(ivr.GT.ivro)THEN
                  ri=gor-(ivr-1)*dvr
               ELSE
                  ri=gor-ivr*dvr
               ENDIF
               vrat(nhp)=(ri-rgr(j))/(rgr(j+1)-rgr(j))
               chp(nhp)=1
            ENDIF
            IF(ivx.NE.ivxo)THEN
               nhp=nhp+1
               IF(ivx.GT.ivxo)THEN
                  xi=gox+(ivx-1)*dvx
               ELSE
                  xi=gox+ivx*dvx
               ENDIF
               rd1=(xi-rgx(j))/(rgx(j+1)-rgx(j))
               IF(nhp.EQ.1)THEN
                  vrat(nhp)=rd1
                  chp(nhp)=2
               ELSE
                  IF(rd1.GE.vrat(nhp-1))THEN
                     vrat(nhp)=rd1
                     chp(nhp)=2
                  ELSE
                     vrat(nhp)=vrat(nhp-1)
                     chp(nhp)=chp(nhp-1)
                     vrat(nhp-1)=rd1
                     chp(nhp-1)=2
                  ENDIF
               ENDIF
            ENDIF
            IF(ivz.NE.ivzo)THEN
               nhp=nhp+1
               IF(ivz.GT.ivzo)THEN
                  zi=goz+(ivz-1)*dvz 
               ELSE
                  zi=goz+ivz*dvz
               ENDIF
               rd1=(zi-rgz(j))/(rgz(j+1)-rgz(j))
               IF(nhp.EQ.1)THEN
                  vrat(nhp)=rd1
                  chp(nhp)=3
               ELSE IF(nhp.EQ.2)THEN
                  IF(rd1.GE.vrat(nhp-1))THEN
                     vrat(nhp)=rd1
                     chp(nhp)=3
                  ELSE
                     vrat(nhp)=vrat(nhp-1)
                     chp(nhp)=chp(nhp-1)
                     vrat(nhp-1)=rd1
                     chp(nhp-1)=3
                  ENDIF
               ELSE
                  IF(rd1.GE.vrat(nhp-1))THEN
                     vrat(nhp)=rd1
                     chp(nhp)=3
                  ELSE IF(rd1.GE.vrat(nhp-2))THEN
                     vrat(nhp)=vrat(nhp-1)
                     chp(nhp)=chp(nhp-1)
                     vrat(nhp-1)=rd1
                     chp(nhp-1)=3
                  ELSE
                     vrat(nhp)=vrat(nhp-1)
                     chp(nhp)=chp(nhp-1)
                     vrat(nhp-1)=vrat(nhp-2)
                     chp(nhp-1)=chp(nhp-2)
                     vrat(nhp-2)=rd1
                     chp(nhp-2)=3
                  ENDIF
               ENDIF
            ENDIF
            nhp=nhp+1
            vrat(nhp)=1.0
            chp(nhp)=0           
            drr=(rgor-(roffset-1)*dnr-rgr(j))-(ipro-1)*dnr
            drx=(rgx(j)-rgox-(xoffset-1)*dnx)-(ipxo-1)*dnx
            drz=(rgz(j)-rgoz-(zoffset-1)*dnz)-(ipzo-1)*dnz
            vel=0.0
            DO k=1,2
               DO l=1,2
                  DO m=1,2
                     produ=(1.0-ABS(((m-1)*dnz-drz)/dnz))
                     produ=produ*(1.0-ABS(((l-1)*dnx-drx)/dnx))
                     produ=produ*(1.0-ABS(((k-1)*dnr-drr)/dnr))   
                     IF(ipzo-1+m.LE.znum+1.AND.ipxo-1+l.LE.xnum+1.AND.ipro-1+k.LE.rnum+1)THEN                
                        vel=vel+veln(ipzo-1+m,ipxo-1+l,ipro-1+k)*produ
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO          
            drr=(gor-rgr(j))-(ivro-1)*dvr
            drx=(rgx(j)-gox)-(ivxo-1)*dvx
            drz=(rgz(j)-goz)-(ivzo-1)*dvz
            u=drr/dvr
            v=drx/dvx
            w=drz/dvz
            ui(1)=(1.0-u)**3/6.0
            ui(2)=(4.0-6.0*u**2+3.0*u**3)/6.0
            ui(3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
            ui(4)=u**3/6.0         
            vi(1)=(1.0-v)**3/6.0
            vi(2)=(4.0-6.0*v**2+3.0*v**3)/6.0
            vi(3)=(1.0+3.0*v+3.0*v**2-3.0*v**3)/6.0
            vi(4)=v**3/6.0            
            wi(1)=(1.0-w)**3/6.0
            wi(2)=(4.0-6.0*w**2+3.0*w**3)/6.0
            wi(3)=(1.0+3.0*w+3.0*w**2-3.0*w**3)/6.0
            wi(4)=w**3/6.0            
            ivrt=ivro
            ivxt=ivxo
            ivzt=ivzo           
            DO k=1,nhp
               velo=vel
               uio=ui
               vio=vi
               wio=wi
               IF(k.GT.1)THEN
                  IF(chp(k-1).EQ.1)THEN
                     ivrt=ivr
                  ELSE IF(chp(k-1).EQ.2)THEN
                     ivxt=ivx
                  ELSE IF(chp(k-1).EQ.3)THEN
                     ivzt=ivz
                  ENDIF
               ENDIF
               rigz=rgz(j)+vrat(k)*(rgz(j+1)-rgz(j))
               rigx=rgx(j)+vrat(k)*(rgx(j+1)-rgx(j))
               rigr=rgr(j)+vrat(k)*(rgr(j+1)-rgr(j))
               iprt=FLOOR((rgor-(roffset-1)*dnr-rigr)/dnr)+1
               ipxt=FLOOR((rigx-rgox-(xoffset-1)*dnx)/dnx)+1
               ipzt=FLOOR((rigz-rgoz-(zoffset-1)*dnz)/dnz)+1              
               drr=(rgor-(roffset-1)*dnr-rigr)-(iprt-1)*dnr
               drx=(rigx-rgox-(xoffset-1)*dnx)-(ipxt-1)*dnx
               drz=(rigz-rgoz-(zoffset-1)*dnz)-(ipzt-1)*dnz               
               vel=0.0
               DO l=1,2
                  DO m=1,2
                     DO n=1,2
                        produ=(1.0-ABS(((n-1)*dnz-drz)/dnz))
                        produ=produ*(1.0-ABS(((m-1)*dnx-drx)/dnx))
                        produ=produ*(1.0-ABS(((l-1)*dnr-drr)/dnr))
                        IF(ipzt-1+n.LE.znum+1.AND.ipxt-1+m.LE.xnum+1.AND.iprt-1+l.LE.rnum+1)THEN
                           vel=vel+veln(ipzt-1+n,ipxt-1+m,iprt-1+l)*produ
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
               drr=(gor-rigr)-(ivrt-1)*dvr
               drx=(rigx-gox)-(ivxt-1)*dvx
               drz=(rigz-goz)-(ivzt-1)*dvz
               u=drr/dvr
               v=drx/dvx
               w=drz/dvz         
               ui(1)=(1.0-u)**3/6.0
               ui(2)=(4.0-6.0*u**2+3.0*u**3)/6.0
               ui(3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
               ui(4)=u**3/6.0
               vi(1)=(1.0-v)**3/6.0
               vi(2)=(4.0-6.0*v**2+3.0*v**3)/6.0
               vi(3)=(1.0+3.0*v+3.0*v**2-3.0*v**3)/6.0
               vi(4)=v**3/6.0
               wi(1)=(1.0-w)**3/6.0
               wi(2)=(4.0-6.0*w**2+3.0*w**3)/6.0
               wi(3)=(1.0+3.0*w+3.0*w**2-3.0*w**3)/6.0
               wi(4)=w**3/6.0                 
               IF(k.EQ.1)THEN
                  dinc=vrat(k)*dpl
               ELSE 
                  dinc=(vrat(k)-vrat(k-1))*dpl
               ENDIF
               DO n=1,4
                  DO m=1,4
                     DO l=1,4
                        rd1=ui(n)*vi(m)*wi(l)/vel**2
                        rd2=uio(n)*vio(m)*wio(l)/velo**2
                        rd1=-(rd1+rd2)*dinc/2.0
                        rd2=fdm(ivzt-2+l,ivxt-2+m,ivrt-2+n)
                        fdm(ivzt-2+l,ivxt-2+m,ivrt-2+n)=rd1+rd2
                     ENDDO
                  ENDDO
               ENDDO  
            ENDDO
            IF(sw.NE.0)THEN
               EXIT
            ENDIF            
         ENDIF        
      ENDDO     
      IF(cfd.EQ.1)THEN      
         isum=0
         DO l=0,nvr+1
            DO k=0,nvx+1
               DO j=0,nvz+1
                  IF(ABS(fdm(j,k,l)).GE.ftol)isum=isum+1
               ENDDO
            ENDDO
         ENDDO         
         rcounts_process=rcounts_process+isum+1  
         coln_process(temp_fdm)=isum
         fdm_process(temp_fdm)=communicate_rpaths(4,i)
         temp_fdm=temp_fdm+1
         isum=0         
         DO j=0,nvz+1
            DO k=0,nvx+1
               DO l=0,nvr+1
                  isum=isum+1			   
                  IF(ABS(fdm(j,k,l)).GE.ftol)THEN
                     coln_process(temp_fdm)=isum
                     fdm_process(temp_fdm)=fdm(j,k,l)
                     temp_fdm=temp_fdm+1
                  ENDIF
               ENDDO 
            ENDDO
         ENDDO         
      ENDIF     
   ENDDO   
   IF(cfd.EQ.1)THEN
      DEALLOCATE(fdm, STAT=checkstat)
      IF(checkstat > 0)THEN
         WRITE(6,*)'Error with DEALLOCATE: SUBROUTINE parallel_rpaths: fdm'
      ENDIF
   ENDIF  
   DEALLOCATE(rgr,rgx,rgz, STAT=checkstat)
   IF(checkstat > 0)THEN
      WRITE(6,*)'Error with DEALLOCATE: SUBROUTINE parallel_rpaths: rgr,rgx,rgz'
   ENDIF
END SUBROUTINE parallel_rpaths

FUNCTION sgn(x)
   implicit none
   INTEGER :: x,sgn
   IF(x.EQ.0)THEN
      sgn=0
   ELSE
      sgn=1
   ENDIF
END FUNCTION
