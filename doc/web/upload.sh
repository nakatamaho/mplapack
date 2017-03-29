#!/bin/sh
#cvs -z3 -d:ext:nakatamaho@mplapack.cvs.sourceforge.net:/cvsroot/mplapack commit . 
rsync --delete -av --exclude '.svn*' --exclude 'CVS*' --exclude '*~' --exclude 'upload.sh' -e ssh . nakatamaho,mplapack@web.sourceforge.net:htdocs/
