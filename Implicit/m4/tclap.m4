# TCLAP M4 Macro 
# Author - Krishna Kumar
AC_DEFUN([TCLAP_PATH], [
# tclap parser library (taken from a static internally bundled version if no --with flag specified)
AC_ARG_WITH([tclap], AS_HELP_STRING([--with-tclap], [indicate the position of the tclap directory]))
if test "$with_tclap"; then
    TCLAP_BUILTIN=no
    TCLAP_CPPFLAGS="-I$with_tclap/include"
else
    AC_MSG_ERROR([tclap path could not be determined])
fi

AC_SUBST(TCLAP_CPPFLAGS)
])
