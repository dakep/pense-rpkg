AC_INIT(pense, version-XX)

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

CC=`"${R_HOME}/bin/R" CMD config CC`
CPP=`"${R_HOME}/bin/R" CMD config CPP`
CXX=`"${R_HOME}/bin/R" CMD config CXX`
CXXCPP=`"${R_HOME}/bin/R" CMD config CXXCPP`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`

AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX

####################################################################################################
##
## Check for necessary support by C compiler
##
####################################################################################################
AC_LANG(C)
AC_CHECK_HEADERS(stdint.h inttypes.h limits.h float.h)

#check for -restrict CFLAGS support
oldCflags="$CFLAGS"
CFLAGS="$CFLAGS -restrict"
AC_MSG_CHECKING([whether C supports -restrict])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
    [AC_MSG_RESULT([yes])]
    [RESTRICT_CFLAGS=-restrict],
    [AC_MSG_RESULT([no])]
)
CFLAGS="$oldCflags"
AC_SUBST([RESTRICT_CFLAGS])

PTHREAD_CPPFLAGS=" -pthread"
$oldCPPFLAGS="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS $PTHREAD_CPPFLAGS"
AC_CHECK_HEADER([pthread.h],[],[
    PTHREAD_CPPFLAGS=" -D_REENTRANT"
    CPPFLAGS="$oldCPPFLAGS $PTHREAD_CPPFLAGS"
    AC_CHECK_HEADER([pthread.h],[],[
        AC_MSG_RESULT([no])
        PTHREAD_CPPFLAGS=""
    ])
])
CPPFLAGS="$oldCPPFLAGS $PTHREAD_CPPFLAGS"
AC_SUBST([PTHREAD_CPPFLAGS])

AC_MSG_CHECKING([for pthread library support])
PTHREAD_LIBS=" -pthread"
oldLIBS="$LIBS"
LIBS="$LIBS $PTHREAD_LIBS"
AC_LINK_IFELSE([
    AC_LANG_PROGRAM([[#include <pthread.h>]],
                    [[pthread_t pthread_self();]])
],[
    AC_MSG_RESULT([yes (-pthread)])
],[
    PTHREAD_LIBS=" -lpthread"
    LIBS="$oldLIBS $PTHREAD_LIBS"
    AC_LINK_IFELSE([
        AC_LANG_PROGRAM([[#include <pthread.h>]],
                        [[pthread_t pthread_self();]])
    ],[
        AC_MSG_RESULT([yes (-lpthread)])
    ],[
        PTHREAD_LIBS=""
        AC_MSG_RESULT([no])
    ])
])
AC_SUBST([PTHREAD_LIBS])

LIBS="$oldLIBS"

####################################################################################################
##
## Check for necessary support by C++ compiler
##
####################################################################################################

AC_LANG(C++)
AC_CHECK_HEADERS(stdint.h inttypes.h climits cfloat)

# Check for "restrict" support
restrict_support="no"

AC_MSG_CHECKING([whether C++ supports restrict keyword])
AC_COMPILE_IFELSE([
    AC_LANG_PROGRAM([[]],
                [[double *restrict ptr;]])
],
[
    AC_MSG_RESULT([yes])
    AC_DEFINE([RESTRICT],[restrict])
    restrict_support=yes
],
[AC_MSG_RESULT([no])]
)

if test $restrict_support != "yes"; then
    AC_MSG_CHECKING([whether C++ supports __restrict__ keyword])
    AC_COMPILE_IFELSE([
        AC_LANG_PROGRAM([[]],
                    [[double *__restrict__ ptr;]])
    ],
    [
        AC_MSG_RESULT([yes])
        AC_DEFINE([RESTRICT],[__restrict__])
        restrict_support="yes"
    ],
    [AC_MSG_RESULT([no])]
    )
fi

if test $restrict_support != "yes"; then
    AC_MSG_CHECKING([whether C++ supports __restrict keyword])
    AC_COMPILE_IFELSE([
        AC_LANG_PROGRAM([[]],
                    [[double *__restrict ptr;]])
    ],
    [
        AC_MSG_RESULT([yes])
        AC_DEFINE([RESTRICT],[__restrict])
        restrict_support=yes
    ],
    [AC_MSG_RESULT([no])]
    )
fi

#check for -restrict CXXFLAGS support
oldCXXflags="$CXXFLAGS"
CXXFLAGS="$CXXFLAGS -restrict"
AC_MSG_CHECKING([whether C++ supports -restrict])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
    [AC_MSG_RESULT([yes])]
    [RESTRICT_CXXFLAGS=-restrict],
    [AC_MSG_RESULT([no])]
)
CXXFLAGS="$oldCXXflags"
AC_SUBST([RESTRICT_CXXFLAGS])


AC_CONFIG_HEADERS([src/autoconfig.h])
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
