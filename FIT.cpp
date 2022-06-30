#!/bin/sh
# (c) 1996 by BukinD
# $id$

usage() {
  echo '' >&2
  echo 'Использование: FIT inname [options]' >&2
  echo '' >&2
  echo ' inname  - имя входного файла (расширение .for,.dat)' >&2
  echo ' options - опции из списка:' >&2
  echo '   debug  - мода отладки' >&2
  echo '   xdebug - тоже самое (но с ddd)' >&2
  echo '   norun  - скомпилировать, но не запускать' >&2
  echo 'переменные среды: (звездочкой помечены advanced features)' >&2
  echo "   SNDROOT - корневая директория, содержащая общие файлы (${SNDROOT:-не определена})" >&2
  echo " * FITLIBS - довесок с библиотеками юзера в makefile (${FITLIBS:-не определена})" >&2
  echo " * FITINCS - довесок с includeами юзера в makefile (${FITINCS:-не определена})" >&2
  echo " * FITMAKE - полное имя makefile-а пользователя (${FITMAKE:-не определено})" >&2
  echo " * FITDEPS - дополнительные зависимости (${FITDEPS:-не определена})" >&2
  echo "   FIT_DEF_DIR - рабочая директория (${FIT_DEF_DIR:-не определена})" >&2
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
  echo 'Не определено SNDMKROOT или SNDROOT'
  exit 1
fi

if [ -z "${FITMAKE}" ]; then
  export FITMAKE
  FITMAKE=${SNDROOT}/include/fit.mk
fi

if [ ! -f ${FITMAKE} ]; then
  echo 'Не могу найти файл fit.mk'
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
  echo "Задание ${infile} подготовлено к работе" >&2
fi
