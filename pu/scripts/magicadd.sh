#!/bin/sh

cat /etc/magic |grep -v 'TLI(Transit' >/tmp/magic.tmp

cat >>/tmp/magic.tmp <./scripts/newmagic 

cat /tmp/magic.tmp >/etc/magic
