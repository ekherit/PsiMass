# echo "Setting JPsi JPsi-00-00-01 in /afs/ihep.ac.cn/users/n/nikolaev/cmthome/workarea"

if test "${CMTROOT}" = ""; then
  CMTROOT=/afs/ihep.ac.cn/bes3/offline/tool/CMT/v1r20; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh

tempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then tempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt setup -sh -pack=JPsi -version=JPsi-00-00-01 -path=/afs/ihep.ac.cn/users/n/nikolaev/cmthome/workarea  -no_cleanup $* >${tempfile}; . ${tempfile}
/bin/rm -f ${tempfile}

