

set start=0

dd if=../oth/pands/h2ofast.bin bs=8 skip=$start count=1 |od -tx1 >test/pands.$start.hex

more test/pands.$start.hex

