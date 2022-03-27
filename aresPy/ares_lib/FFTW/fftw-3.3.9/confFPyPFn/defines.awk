BEGIN {
D["PACKAGE_NAME"]=" \"fftw\""
D["PACKAGE_TARNAME"]=" \"fftw\""
D["PACKAGE_VERSION"]=" \"3.3.9\""
D["PACKAGE_STRING"]=" \"fftw 3.3.9\""
D["PACKAGE_BUGREPORT"]=" \"fftw@fftw.org\""
D["PACKAGE_URL"]=" \"\""
D["PACKAGE"]=" \"fftw\""
D["VERSION"]=" \"3.3.9\""
D["FFTW_ENABLE_ALLOCA"]=" 1"
D["STDC_HEADERS"]=" 1"
D["HAVE_SYS_TYPES_H"]=" 1"
D["HAVE_SYS_STAT_H"]=" 1"
D["HAVE_STDLIB_H"]=" 1"
D["HAVE_STRING_H"]=" 1"
D["HAVE_MEMORY_H"]=" 1"
D["HAVE_STRINGS_H"]=" 1"
D["HAVE_INTTYPES_H"]=" 1"
D["HAVE_STDINT_H"]=" 1"
D["HAVE_UNISTD_H"]=" 1"
D["HAVE_DLFCN_H"]=" 1"
D["LT_OBJDIR"]=" \".libs/\""
D["HAVE_MPI"]=" 1"
D["SIZEOF_MPI_FINT"]=" 4"
D["STDC_HEADERS"]=" 1"
D["HAVE_FCNTL_H"]=" 1"
D["HAVE_FENV_H"]=" 1"
D["HAVE_LIMITS_H"]=" 1"
D["HAVE_MALLOC_H"]=" 1"
D["HAVE_STDDEF_H"]=" 1"
D["HAVE_SYS_TIME_H"]=" 1"
D["TIME_WITH_SYS_TIME"]=" 1"
D["HAVE_LONG_DOUBLE"]=" 1"
D["SIZEOF_INT"]=" 4"
D["SIZEOF_UNSIGNED_INT"]=" 4"
D["SIZEOF_LONG"]=" 8"
D["SIZEOF_UNSIGNED_LONG"]=" 8"
D["SIZEOF_LONG_LONG"]=" 8"
D["SIZEOF_UNSIGNED_LONG_LONG"]=" 8"
D["SIZEOF_SIZE_T"]=" 8"
D["SIZEOF_PTRDIFF_T"]=" 8"
D["HAVE_PTRDIFF_T"]=" 1"
D["HAVE_UINTPTR_T"]=" 1"
D["SIZEOF_FLOAT"]=" 4"
D["SIZEOF_DOUBLE"]=" 8"
D["SIZEOF_FFTW_R2R_KIND"]=" 4"
D["HAVE_ALLOCA_H"]=" 1"
D["HAVE_ALLOCA"]=" 1"
D["HAVE_VPRINTF"]=" 1"
D["HAVE_LIBM"]=" 1"
D["HAVE_GETTIMEOFDAY"]=" 1"
D["HAVE_DRAND48"]=" 1"
D["HAVE_SQRT"]=" 1"
D["HAVE_MEMSET"]=" 1"
D["HAVE_POSIX_MEMALIGN"]=" 1"
D["HAVE_MEMALIGN"]=" 1"
D["HAVE__MM_MALLOC"]=" 1"
D["HAVE__MM_FREE"]=" 1"
D["HAVE_CLOCK_GETTIME"]=" 1"
D["HAVE_SYSCTL"]=" 1"
D["HAVE_ABORT"]=" 1"
D["HAVE_SINL"]=" 1"
D["HAVE_COSL"]=" 1"
D["HAVE_SNPRINTF"]=" 1"
D["HAVE_MEMMOVE"]=" 1"
D["HAVE_STRCHR"]=" 1"
D["HAVE_GETPAGESIZE"]=" 1"
D["HAVE_DECL_SINL"]=" 1"
D["HAVE_DECL_COSL"]=" 1"
D["HAVE_DECL_SINQ"]=" 0"
D["HAVE_DECL_COSQ"]=" 0"
D["HAVE_DECL_MEMALIGN"]=" 1"
D["HAVE_DECL_DRAND48"]=" 1"
D["HAVE_DECL_SRAND48"]=" 1"
D["HAVE_DECL_POSIX_MEMALIGN"]=" 1"
D["HAVE_ISNAN"]=" 1"
P["F77_FUNC"]="(name,NAME)"
D["F77_FUNC"]=" name ## _"
P["F77_FUNC_"]="(name,NAME)"
D["F77_FUNC_"]=" name ## _"
D["F77_FUNC_EQUIV"]=" 1"
D["WITH_G77_WRAPPERS"]=" 1"
D["FFTW_CC"]=" \"icc -std=gnu99 -O3 -ansi-alias -malign-double\""
  for (key in D) D_is_set[key] = 1
  FS = ""
}
/^[\t ]*#[\t ]*(define|undef)[\t ]+[_abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ][_abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789]*([\t (]|$)/ {
  line = $ 0
  split(line, arg, " ")
  if (arg[1] == "#") {
    defundef = arg[2]
    mac1 = arg[3]
  } else {
    defundef = substr(arg[1], 2)
    mac1 = arg[2]
  }
  split(mac1, mac2, "(") #)
  macro = mac2[1]
  prefix = substr(line, 1, index(line, defundef) - 1)
  if (D_is_set[macro]) {
    # Preserve the white space surrounding the "#".
    print prefix "define", macro P[macro] D[macro]
    next
  } else {
    # Replace #undef with comments.  This is necessary, for example,
    # in the case of _POSIX_SOURCE, which is predefined and required
    # on some systems where configure will not decide to define it.
    if (defundef == "undef") {
      print "/*", prefix defundef, macro, "*/"
      next
    }
  }
}
{ print }
