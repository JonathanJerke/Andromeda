#!/bin/bash

grep -rl solo .         | xargs sed -i 's/solo/vector/g'
grep -rl canon .         | xargs sed -i 's/canon/lambda/g'
grep -rl basisRank .    | xargs sed -i 's/basisRank/canonRank/g'
grep -rl basisStage .   | xargs sed -i 's/basisStage/canonStage/g'
grep -rl type .         | xargs sed -i 's/type/irrep/g'
