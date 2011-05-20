if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /afs/ihep.ac.cn/bes3/offline/tool/CMT/v1r20
endif
source ${CMTROOT}/mgr/setup.csh
set tempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set tempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt cleanup -csh -pack=JPsi -version=JPsi-00-00-01 -path=/ihepbatch/bes/nikolaev/work655 $* >${tempfile}; source ${tempfile}
/bin/rm -f ${tempfile}

