AC_INIT(pense, 2.0.1)
AC_PREREQ(2.62)

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

RBIN="${R_HOME}/bin/R"
CC=`"${RBIN}" CMD config CC`
CPP="$CC -E"

CXX11=`"${RBIN}" CMD config CXX11`
if test -z "$CXX11"; then
  AC_MSG_ERROR([No C++11 compiler is available])
fi
CXX11STD=`"${RBIN}" CMD config CXX11STD`
CXX="${CXX11} ${CXX11STD}"
CXXFLAGS=`"${RBIN}" CMD config CXX11FLAGS`
CPPFLAGS=`"${RBIN}" CMD config CPPFLAGS`
LDFLAGS=`"${RBIN}" CMD config LDFLAGS`

AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX

################################################################################
##
## Check for necessary support by C++ compiler
##
################################################################################
AC_LANG(C++)

## Check for strict aliasing support
STRICT_ALIASING_FLAG="-fstrict-aliasing"
AC_MSG_CHECKING([whether _AC_LANG compiler accepts $STRICT_ALIASING_FLAG])
ax_check_save_flags=$[]_AC_LANG_PREFIX[]FLAGS
_AC_LANG_PREFIX[]FLAGS="$[]_AC_LANG_PREFIX[]FLAGS -Werror $STRICT_ALIASING_FLAG -Wstrict-aliasing"
AC_COMPILE_IFELSE([AC_LANG_PROGRAM()],
                  [AC_MSG_RESULT([yes])],
                  [
                    AC_MSG_RESULT([no])
                    STRICT_ALIASING_FLAG=""
                  ])
_AC_LANG_PREFIX[]FLAGS=$ax_check_save_flags

AC_SUBST(STRICT_ALIASING_FLAG)

## Check for OpenMP support
R_HAVE_OPENMP='FALSE'
AC_OPENMP

# Some systems have broken OMP libraries, so also check that the actual package will work
ac_pkg_openmp_cxx=no
AC_MSG_CHECKING([C++ support for OpenMP])
AC_LANG_CONFTEST([
  AC_LANG_PROGRAM([[#include <omp.h>]], [[ return omp_get_num_threads (); ]])
])
PKG_CFLAGS="${OPENMP_CXXFLAGS}" PKG_LIBS="${OPENMP_CXXFLAGS}" "$RBIN" CMD SHLIB conftest.c 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && \
  "$RBIN" --vanilla -q -e "dyn.load(paste('conftest',.Platform\$dynlib.ext,sep=''))" 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && ac_pkg_openmp_cxx=yes
AC_MSG_RESULT([${ac_pkg_openmp_cxx}])

# Check if OpenMP supports tasks, reductions, and if constant member variables need to be shared
LDFLAGS="$LDFLAGS $OPENMP_CXXFLAGS"
ax_check_save_flags="$LDFLAGS"
if test "${ac_pkg_openmp_cxx}" = yes; then
  AC_MSG_CHECKING([OpenMP support for tasks])

  # Check if tasks are supported with constant member not shared
  AC_LINK_IFELSE([AC_LANG_PROGRAM([[
    #include <omp.h>
    class Foo {
    public:
      Foo(const int x) : x_(x) {}
      int Test(const int y) {
        int sum = 0;
        #pragma omp parallel num_threads(2) default(none) firstprivate(y) shared(sum)
        {
          #pragma omp single nowait
          for (int i = 0; i < y; ++i) {
            #pragma omp task default(none) firstprivate(i) shared(sum)
            {
              sum += x_ * i;
            }
          }
        }
        return sum;
      }
    private:
      const int x_;
    };
  ]],
  [[
    const int mult = omp_get_num_threads();
    Foo f(mult);
    return f.Test(2 * mult);
  ]])],
  [
    ac_pkg_openmp_cxx_tasks=yes
  ],
  [
    ac_pkg_openmp_cxx_tasks=no
  ])

  if test "${ac_pkg_openmp_cxx_tasks}" = no; then
    # Maybe constant member variables need to be shared
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[
      #include <omp.h>
      class Foo {
      public:
        Foo(const int x) : x_(x) {}
        int Test(const int y) {
          int sum = 0;
          #pragma omp parallel num_threads(2) default(none) shared(sum)
          {
            #pragma omp single nowait
            for (int i = 0; i < y; ++i) {
              #pragma omp task default(none) firstprivate(i) shared(sum, x_)
              {
                sum += x_ * i;
              }
            }
          }
          return sum;
        }
      private:
        const int x_;
      };
    ]],
    [[
      const int mult = omp_get_num_threads();
      Foo f(mult);
      return f.Test(2 * mult);
    ]])],
    [
      AC_DEFINE([PENSE_OPENMP_ADD_CONST_SHARED], [1], [
        Set to 1 if OpenMP requires constant class members to be shared
      ])
      ac_pkg_openmp_cxx_tasks="yes, with shared constant members"
    ],
    [
      ac_pkg_openmp_cxx_tasks=no
      ac_pkg_openmp_cxx=no
    ])
  fi

  AC_MSG_RESULT([$ac_pkg_openmp_cxx_tasks])
fi

if test "${ac_pkg_openmp_cxx}" = yes; then
  AC_MSG_CHECKING([OpenMP support for custom reductions])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([[
    #include <omp.h>
    #include <forward_list>
    #include <type_traits>

    namespace omp_reductions {
    template<typename T>
    using FwdList = std::forward_list<T>;

    class Foo {
    public:
      explicit Foo(const int x) : x_(x) {}
      int x() const {
        return x_;
      }
    private:
      int x_;
    };

    void CustomReduction(FwdList<Foo>* single, FwdList<Foo>* combined) {
      combined->splice_after(combined->before_begin(), *single);
    }

    template<typename T, typename = typename std::enable_if<!std::is_arithmetic<T>::value>::type>
    T ComputeWithReduction(const T mult, const Foo & y) {
      return Foo(mult.x() * y.x());
    }

    template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    T ComputeWithReduction(const T mult, const int y) {
      using FooList = FwdList<Foo>;

      #pragma omp declare reduction(c:FwdList<Foo>:CustomReduction(&omp_in, &omp_out))

      FooList combined;
      int num_threads = omp_get_num_threads();
      #pragma omp parallel num_threads(num_threads) default(none) firstprivate(y) shared(combined)
      {
        #pragma omp for reduction(c:combined)
        for (int i = 0; i < y; ++i) {
          combined.emplace_front(i);
        }
      }

      int sum = 0;
      for (auto&& f : combined) {
        sum += f.x();
      }
      return mult * sum;
    }
    }  // namespace omp_reductions
    ]],
    [[
      using omp_reductions::Foo;
      double sum = omp_reductions::ComputeWithReduction(2.3, 12);
      Foo sumfoo = omp_reductions::ComputeWithReduction(Foo(5), Foo(12));
      return sum + sumfoo.x();
    ]])],
    [ac_pkg_openmp_cxx=yes],
    [ac_pkg_openmp_cxx=no])
  AC_MSG_RESULT([$ac_pkg_openmp_cxx])
fi
LDFLAGS=$ax_check_save_flags

# if ${ac_pkg_openmp_cxx} = "yes" then OMP is good, otherwise it will be "no"
if test "${ac_pkg_openmp_cxx}" = no; then
  OPENMP_CXXFLAGS=''
  AC_DEFINE([PENSE_DISABLE_OPENMP], [1], [Set to 1 if the compiler does not support OpenMP])
else
  AC_DEFINE([PENSE_ENABLE_OPENMP], [1], [Set to 1 if the compiler does supports recent version of OpenMP])
  R_HAVE_OPENMP='TRUE'
fi

AC_SUBST(OPENMP_CXXFLAGS)
AC_SUBST(R_HAVE_OPENMP)

##
## Write output
##
AC_CONFIG_HEADERS([src/autoconfig.hpp])
AC_CONFIG_FILES([src/Makevars R/autoconfig.R])
AC_OUTPUT
