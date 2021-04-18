SUBDIRS = reference 

#if ENABLE_OPTIMIZED
#SUBDIRS = optimized
#endif

if ENABLE_WRAPPER
SUBDIRS+= lapackwrapper
endif

if ENABLE_FORTRANWRAPPER
SUBDIRS+= fortranwrapper
endif

if ENABLE_TEST
SUBDIRS+= test
endif

