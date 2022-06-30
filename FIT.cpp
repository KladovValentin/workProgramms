#!/bin/sh
# (c) 1996 by BukinD
# $id$

usage() {
  echo '' >&2
  echo '�������������: FIT inname [options]' >&2
  echo '' >&2
  echo ' inname  - ��� �������� ����� (���������� .for,.dat)' >&2
  echo ' options - ����� �� ������:' >&2
  echo '   debug  - ���� �������' >&2
  echo '   xdebug - ���� ����� (�� � ddd)' >&2
  echo '   norun  - ��������������, �� �� ���������' >&2
  echo '���������� �����: (���������� �������� advanced features)' >&2
  echo "   SNDROOT - �������� ����������, ���������� ����� ����� (${SNDROOT:-�� ����������})" >&2
  echo " * FITLIBS - ������� � ������������ ����� � makefile (${FITLIBS:-�� ����������})" >&2
  echo " * FITINCS - ������� � include��� ����� � makefile (${FITINCS:-�� ����������})" >&2
  echo " * FITMAKE - ������ ��� makefile-� ������������ (${FITMAKE:-�� ����������})" >&2
  echo " * FITDEPS - �������������� ����������� (${FITDEPS:-�� ����������})" >&2
  echo "   FIT_DEF_DIR - ������� ���������� (${FIT_DEF_DIR:-�� ����������})" >&2
  echo '' >&2
  exit 2
}

echo "!----------------------------------------------------------------!"
echo "!         ### # ###          V.N.Ivanchenko                      !"
echo "!         ##  #  #     by    A.V.Bozhenok                        !"
echo "!         #   #  #           D.A.Bukin                           !"
echo "!----------------------------------------------------------------!"
date

[ $# -lt 1 ] && usage

infile="$1"

# parse options
debugger=
norun=
if [ $# -gt 1 ]; then
  # skip first argument
  shift
  while [ $# != 0 ]
   do \
    case "$1" in
     debug)	debugger=gdb ;;
     xdebug)	debugger=ddd ;;
     norun)	norun=1 ;;
     *)		usage ;;
    esac
    shift
   done
fi

if [ -z "${SNDMKROOT}" -o -z "${SNDROOT}" ]; then
  echo '�� ���������� SNDMKROOT ��� SNDROOT'
  exit 1
fi

if [ -z "${FITMAKE}" ]; then
  export FITMAKE
  FITMAKE=${SNDROOT}/include/fit.mk
fi

if [ ! -f ${FITMAKE} ]; then
  echo '�� ���� ����� ���� fit.mk'
  exit 1
fi

savedir=`pwd`
if [ -z "$FIT_DEF_DIR" ] ; then
  export FIT_DEF_DIR
  FIT_DEF_DIR=${savedir}
else
  cd ${FIT_DEF_DIR}
fi

if [ -z "${FITLIBS}" ]; then
  export FITLIBS
  FITLIBS=
fi

if [ -z "${FITINCS}" ]; then
  export FITINCS
  FITINCS=
fi

if [ -z "${FITDEPS}" ]; then
  export FITDEPS
  FITDEPS=
fi

echo "current directory ${FIT_DEF_DIR}"
echo Host: `uname -n`

# if [ "$ENVIRONMENT" = "BATCH" -a -n "${TMPDIR}" ] ; then
#   exedir="${TMPDIR}/"
#   exediropt="FITEXEDIR=${exedir}"
# else
  exedir=./
  exediropt=
# fi

make -f ${FITMAKE} ${exediropt} FITEXE=${infile} DEBUG=${debugger}

if [ $? -ne 0 ]; then
  echo "Make status $?"
  exit 1
fi

if [ -z "${norun}" ]; then
  if [ -n "${debugger}" ]; then
    ${debugger} ${exedir}${infile}
  else
    time ${exedir}$infile
    if [ $? -eq 0 ]; then
      make -f ${FITMAKE} ${exediropt} FITEXE=${infile} clean
    fi
  fi
else
  echo "������� ${infile} ������������ � ������" >&2
fi
