AC_INIT(pense, 0.6.6)

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

################################################################################
##
## Check for necessary support by C compiler
##
################################################################################
AC_LANG(C)
AC_CHECK_HEADERS(stdint.h inttypes.h limits.h float.h)

##
## check for -restrict CFLAGS support
##
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

################################################################################
##
## Check for necessary support by C++ compiler
##
################################################################################
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

##
## check for -restrict CXXFLAGS support
##
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
