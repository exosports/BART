#!/bin/tcsh

if ( ! -x $1 ) then
    echo "cproto is not installed and hence, cannot 'make proto'. STOPING..."
    exec false
endif

cproto -V >& cproto.ver
set ver = `cat cproto.ver`
\rm -f cproto.ver

if ( $ver == "4.7c" || $ver == "4.7d" || $ver == "4.7e" || $ver == "4.7f" || $ver == "4.8" ) exec true

echo "  Sorry, but your cproto is not version '4.7c', '4.7d', '4.7e' nor '4.8', cannot 'make proto'."
echo "Note that 4.8 was non-existent at the time this was written:D, but guessed to"
echo "support inline functions. STOPING..."
exec false


