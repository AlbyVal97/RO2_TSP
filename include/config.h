/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE                                           */
/*                                                                          */
/*  (c) Copyright 1995--1999 by David Applegate, Robert Bixby,              */
/*  Vasek Chvatal, and William Cook                                         */
/*                                                                          */
/*  Permission is granted for academic research use.  For other uses,       */
/*  contact the authors for licensing options.                              */
/*                                                                          */
/*  Use at your own risk.  We make no guarantees about the                  */
/*  correctness or usefulness of this code.                                 */
/*                                                                          */
/****************************************************************************/


/* Define if your compiler is missing the appropriate function prototype */

#undef CC_PROTO_PRINTF
#undef CC_PROTO_GETHOSTNAME
#undef CC_PROTO_GETRUSAGE

/* Define if you want to use posix threads */
#undef CC_POSIXTHREADS

/* Define if <signal.h> needs to be included before <pthreads.h> */
#undef CC_SIGNAL_BEFORE_PTHREAD

/* Define to empty if the keyword `const' does not work.  */
#undef const

/* Define to `int' if <sys/types.h> doesn't define.  */
#undef pid_t

/* Define to `unsigned' if <sys/types.h> doesn't define.  */
#undef size_t

/* Define to `unsigned char' if <sys/types.h> doesn't define.  */
#undef u_char

/* Define to `int' if the builtin type `void' does not work.  */
#undef void

/* Define as the return type of signal handlers (int or void).  */
#undef RETSIGTYPE

/* Define one of the following three to specify the type of signal
 * handling to use. */
#undef  CCSIGNAL_SIGACTION /* sigaction(), preferred */
#undef  CCSIGNAL_SIGNAL    /* signal() */
#undef  CCSIGNAL_NONE      /* no signal handling */

/* Define if you have the gethostname function.  */
#undef HAVE_GETHOSTNAME

/* Define if you have the socket function.  */
#undef HAVE_SOCKET

/* Define if you have the strdup function.  */
#undef HAVE_STRDUP

/* Define if you have the getrusage function.  */
#undef HAVE_GETRUSAGE

/* Define if you have the times function.  */
#undef HAVE_TIMES

/* Define if you have the clock function.  */
#undef HAVE_CLOCK

/* Define if you have the sleep function.  */
#undef HAVE_SLEEP

/* Define if you have the <stdlib.h> header file.  */
#define HAVE_STDLIB_H

/* Define if you have the <math.h> header file.  */
#undef HAVE_MATH_H

/* Define if you have the <string.h> header file.  */
#define HAVE_STRING_H

/* Define if you have the <strings.h> header file.  */
#undef HAVE_STRINGS_H

/* Define if you have the <errno.h> header file.  */
#undef HAVE_ERRNO_H

/* Define if you have the <assert.h> header file.  */
#undef HAVE_ASSERT_H

/* Define if you can safely include both <sys/time.h> and <time.h>.  */
#undef TIME_WITH_SYS_TIME

/* Define if you have the <sys/time.h> header file.  */
#undef HAVE_SYS_TIME_H

/* Define if you have the <time.h> header file.  */
#undef HAVE_TIME_H

/* Define if you have the <stddef.h> header file.  */
#undef HAVE_STDDEF_H

/* Define if you have the <unistd.h> header file.  */
#undef HAVE_UNISTD_H

/* Define if you have the <malloc.h> header file.  */
#undef HAVE_MALLOC_H

/* Define if you have the <sys/types.h> header file.  */
#undef HAVE_SYS_TYPES_H

/* Define if you have the <sys/stat.h> header file.  */
#undef HAVE_SYS_STAT_H

/* Define if you have the <sys/resource.h> header file.  */
#undef HAVE_SYS_RESOURCE_H

/* Define if you have the <fcntl.h> header file.  */
#undef HAVE_FCNTL_H

/* Define if you have the <signal.h> header file.  */
#undef HAVE_SIGNAL_H

/* Define if you have the <sys/socket.h> header file.  */
#undef HAVE_SYS_SOCKET_H

/* Define if you have the <netdb.h> header file.  */
#undef HAVE_NETDB_H

/* Define if you have the <netinet/in.h> header file.  */
#undef HAVE_NETINET_IN_H

/* Define if you have the <sys/param.h> header file.  */
#undef HAVE_SYS_PARAM_H

/* Define if you have the <sys/times.h> header file.  */
#undef HAVE_SYS_TIMES_H

/* Define if your compiler supports attribute modifiers  */
/* such as __attribute__ ((unused)) (gcc 2.8.1 does)     */
#undef CC_ATTRIBUTE

/* Define if your header files use non-Ansi casts for SIG_ERR, SIG_IGN, */
/* or SIG_DFL */
#undef CC_BADSIGDEF_CAST

/* Some machine (o/s) specific problems */

/* Define if unistd.h uses __vfork but does not prototype it */
/* This happens under Irix 6 */
#undef CC_PROTO___VFORK

