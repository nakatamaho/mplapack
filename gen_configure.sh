#!/bin/bash
set -eu

aclocal
autoheader
automake -a -v --add-missing
autoconf
