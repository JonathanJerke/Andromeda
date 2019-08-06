#!/bin/csh

foreach r (`cat states`)
foreach c (`cat catalog`)
	grep C${c}C $1-$r.vector > $1-${r}C${c}C.vector
end
end

