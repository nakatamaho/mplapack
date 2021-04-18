SUBDIRS = optimized

if !DISABLE_REFERENCE
SUBDIRS += reference 
endif

if ENABLE_WRAPPER
SUBDIRS+= lapackwrapper
endif

if ENABLE_FORTRANWRAPPER
SUBDIRS+= fortranwrapper
endif

if ENABLE_TEST
SUBDIRS+= test
endif

