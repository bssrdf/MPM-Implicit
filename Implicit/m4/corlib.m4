# Corlib Path M4 Macro
# Author - Krishna Kumar
AC_DEFUN([CORLIB_PATH], [
# corlib parser library (taken from a static internally bundled version if no --with flag specified)
AC_ARG_WITH([corlib], AS_HELP_STRING([--with-corlib], [indicate the position of the corlib directory]))
AC_ARG_WITH([corlib-lib], AS_HELP_STRING([--with-corlib-lib], [indicate the position of the corlib lib directory (the one containing libcorlib.so)]))
AC_ARG_WITH([corlib-inc], AS_HELP_STRING([--with-corlib-inc], [indicate the position of the corlib include directory (the one *containing* the 'corlib' directory)]))
if test "$with_corlib" && (test -z "$with_corlib_lib" || test -z "$with_corlib_inc"); then
    CORLIB_BUILTIN=no
    CORLIB_CPPFLAGS="-I$with_corlib/include"
    CORLIB_LDFLAGS="-L$with_corlib/lib"
elif test "$with_corlib_lib" && test "$with_corlib_inc"; then
    CORLIB_BUILTIN=no
    CORLIB_CPPFLAGS="-I$with_corlib_inc"
    CORLIB_LDFLAGS="-L$with_corlib_lib"
elif test -z "$with_corlib" && test -z "$with_corlib_lib" && test -z "$with_corlib_inc"; then
    AC_MSG_NOTICE([Using built-in corlib i.e., /usr/lib (or) system-wide location])
    CORLIB_BUILTIN=yes
    CORLIB_CPPFLAGS=
    CORLIB_LDFLAGS=
else
    AC_MSG_ERROR([corlib path could not be determined])
fi

AM_CONDITIONAL(USE_BUILTIN_CORLIB, [test $CORLIB_BUILTIN = yes])
AC_SUBST(CORLIB_CPPFLAGS)
AC_SUBST(CORLIB_LDFLAGS)
])
