#!/bin/sh

if [ "$1" = "-u" ]; then 
    echo -n " Uni";
else
    echo -n " I";
fi

echo -n "nstalling TLI's magic bit... "
cat /etc/magic |grep -v 'TLI(Transit' >/tmp/magic.tmp

if [ "$1" != "-u" ]; then 
    cat >>/tmp/magic.tmp <./scripts/newmagic 
fi

cat /tmp/magic.tmp >/etc/magic  && echo "Done" || echo "Failed"
