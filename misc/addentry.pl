#!/usr/bin/perl

#adding entry in the html document.

$MPNAME=$ARGV[0]; 
$INITIAL=substr($MPNAME,0,1);
$MIDDLE=substr($MPNAME,1);
print "$INITIAL\n";
if ($INITIAL eq 'C'){
$DNAME = "z" . "$MIDDLE";
$SNAME = "c" . "$MIDDLE";
} else {
$DNAME = "d" . "$MIDDLE";
$SNAME = "s" . "$MIDDLE";
}

print "<TR VALIGN=TOP>\n";
print "<TD WIDTH=25%>  <P><a href=\"http://mplapack.cvs.sourceforge.net/viewvc/*checkout*/mplapack/mplapack/mpblas/reference/$MPNAME.cpp\">$MPNAME</a> </P> </TD>\n";

print "<TD WIDTH=25%>  <P> x := A*x, or x := A'*x </P> </TD>\n";
print "<TD WIDTH=25%>  <P>                     </P> </TD>\n";
print "<TD WIDTH=25%>  <P><a href=\"http://www.netlib.org/blas/$SNAME.f\">$SNAME</a>,\n";
print "<a href=\"http://www.netlib.org/blas/$DNAME.f\">$DNAME</a></p>\n";
print "</TR>\n";
exit
