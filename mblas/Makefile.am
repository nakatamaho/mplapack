SUBDIRS = optimized

if !DISABLE_REFERENCE
SUBDIRS += reference 
endif

if ENABLE_WRAPPER
SUBDIRS+= blaswrapper
endif

if ENABLE_FORTRANWRAPPER
SUBDIRS+= fortranwrapper
endif

if ENABLE_DEBUG
SUBDIRS+= debug
endif
