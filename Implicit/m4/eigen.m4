# EIGEN M4 Macro 
# Author - Krishna Kumar
AC_DEFUN([EIGEN_PATH], [
# eigen parser library (taken from a static internally bundled version if no --with flag specified)
AC_ARG_WITH([eigen], AS_HELP_STRING([--with-eigen], [indicate the position of the eigen directory]))
if test "$with_eigen"; then
    EIGEN_BUILTIN=no
    EIGEN_CPPFLAGS="-I$with_eigen/include"
else
    AC_MSG_ERROR([eigen path could not be determined])
fi

AC_SUBST(EIGEN_CPPFLAGS)
])
