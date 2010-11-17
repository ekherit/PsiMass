# echo "Setting JPsi JPsi-00-00-01 in /afs/ihep.ac.cn/users/n/nikolaev/cmthome/workarea"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /afs/ihep.ac.cn/bes3/offline/tool/CMT/v1r20
endif
source ${CMTROOT}/mgr/setup.csh

set tempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set tempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt setup -csh -pack=JPsi -version=JPsi-00-00-01 -path=/afs/ihep.ac.cn/users/n/nikolaev/cmthome/workarea  -no_cleanup $* >${tempfile}; source ${tempfile}
/bin/rm -f ${tempfile}

