
/***************************************************************************
 *
 * $Id$
 *
 * CODE DEVELOPER
 *   John Tannahill, LLNL (Original code from Bill Bosl, LLNL)
 *   jrt@llnl.gov
 *
 * FILE
 *   gem_sys_options.h
 *
 * DESCRIPTION
 *   The options specified in the following #define, etc. clauses are used
 *   to build the desired executable program.
 *
 * HISTORY
 *   - July 1, 2004 - Jules Kouatchou
 *      o Added the pre-processing option MPI_2_OPTION to remove/include
 *        MPI-2 calls in the code at compilation. MPI-2 is not installed
 *        on all the platforms. This new option replaces the parameter
 *        "do_mpi2_border".
 *      o The pre-processing option MSG_OPTION can now be set to MSG_NONE 
 *        on any platform.
 ***************************************************************************/

/* Current ARCH_OPTION choices: */
#define ARCH_COMPAQ     1
#define ARCH_CRAY       2
#define ARCH_IBM_SP     3
#define ARCH_INTEL      4
#define ARCH_SGI_ORIG   5
#define ARCH_SUN4       6
#define ARCH_T3E        7

/* Current HOST_MACH choices: */
#define PALM 1
#define DISCOVER 2

/* 
   Choose architecture and host machine for current machine: 
*/

#define ARCH_OPTION ARCH_INTEL
#define HOST_MACH DISCOVER 

/* 
   Define the message passing option
   MSG_NONE allows you to run the code without MPI and
   can be selected for any platform. 
*/
#define MSG_NONE  0
#define MSG_MPI   1

/* 
   Define MPI-2 options
   WITH_MPI_2 can only be considered if MPI is selected.
   NO_MPI_2 freezes all the MPI-2 calls within the code.
   This option replaces the parameter do_mpi2_border 
*/
#define NO_MPI_2   0
#define WITH_MPI_2 1
/*
   Make your selection here
*/
#if ((ARCH_OPTION == ARCH_CRAY) || (ARCH_OPTION == ARCH_SGI_IND) || \
     (ARCH_OPTION == ARCH_SUN4))
#  define MSG_OPTION    MSG_NONE
#  define MPI_2_OPTION  NO_MPI_2
#else
#  define MSG_OPTION    MSG_NONE
#  define MPI_2_OPTION  NO_MPI_2
#endif


#define CODENAME  "gem"

