```bash
voms-proxy-init --rfc --voms cms

## Momentum Resolution

echo "#define Y2016">IsData.h
FILES=files/2016/ZToMuMu*.txt
for i in $FILES; do  root -l -b -q "RunMomentumResolution.C(\"$i\", 2)";  done

echo "#define Y2017">IsData.h
FILES=files/UL2017/ZToMuMu*.txt
for i in $FILES; do  root -l -b -q "RunMomentumResolution.C(\"$i\", 2)";  done

echo "#define Y2018">IsData.h
FILES=files/UL2018/ZToMuMu*.txt
for i in $FILES; do  root -l -b -q "RunMomentumResolution.C(\"$i\", 2)";  done


## Z Peak

echo "#define Y2016">IsData.h
FILES=files/UL2016/DY*.txt
for i in $FILES; do  root -l -b -q "RunZPeakResolution.C(\"$i\", 2)";  done

echo "#define Y2017">IsData.h
FILES=files/UL2017/DY*.txt
for i in $FILES; do  root -l -b -q "RunZPeakResolution.C(\"$i\", 2)";  done

echo "#define Y2018">IsData.h
FILES=files/UL2018/DY*.txt
for i in $FILES; do  root -l -b -q "RunZPeakResolution.C(\"$i\", 2)";  done

echo -e "#define Y2016\n#define CMSDATA">IsData.h
FILES=files/UL2016/data/ULSingleMuon.txt
for i in $FILES; do  root -l -b -q "RunZPeakResolution.C(\"$i\", 6)";  done

echo -e "#define Y2017\n#define CMSDATA">IsData.h
FILES=files/UL2017/data/ULSingleMuon.txt
for i in $FILES; do  root -l -b -q "RunZPeakResolution.C(\"$i\", 6)";  done

echo -e "#define Y2018\n#define CMSDATA">IsData.h
FILES=files/UL2018/data/ULSingleMuon.txt
for i in $FILES; do  root -l -b -q "RunZPeakResolution.C(\"$i\", 6)";  done
```