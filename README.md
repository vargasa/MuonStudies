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

#1st Option

echo "#define Y2016">IsData.h
FILES=files/UL2016/DYJetsToMuMu_M-50*.txt
for i in $FILES; do  root -l -b -q "RunZPeakResolution.C(\"$i\", 2)";  done

echo "#define Y2017">IsData.h
FILES=files/UL2017/DYJetsToMuMu_M-50*.txt
for i in $FILES; do  root -l -b -q "RunZPeakResolution.C(\"$i\", 2)";  done

echo "#define Y2018">IsData.h
FILES=files/UL2018/DYJetsToMuMu_M-50*.txt
for i in $FILES; do  root -l -b -q "RunZPeakResolution.C(\"$i\", 2)";  done


## Jet binning

echo "#define Y2016">IsData.h
FILES=files/UL2016/DYJetsToLL_*J*.txt
for i in $FILES; do  root -l -b -q "RunZPeakResolution.C(\"$i\", 2)";  done

echo "#define Y2017">IsData.h
FILES=files/UL2017/DYJetsToLL_*J*.txt
for i in $FILES; do  root -l -b -q "RunZPeakResolution.C(\"$i\", 2)";  done

echo "#define Y2018">IsData.h
FILES=files/UL2018/DYJetsToLL_*J*.txt
for i in $FILES; do  root -l -b -q "RunZPeakResolution.C(\"$i\", 2)";  done



## HT-Binning

echo "#define Y2016">IsData.h
FILES=files/2016/DYJetsToLL_*HT*.txt
for i in $FILES; do  root -l -b -q "RunZPeakResolution.C(\"$i\", 4)";  done

echo "#define Y2017">IsData.h
FILES=files/UL2017/DYJetsToLL_*HT*.txt
for i in $FILES; do  root -l -b -q "RunZPeakResolution.C(\"$i\", 4)";  done

echo "#define Y2018">IsData.h
FILES=files/UL2018/DYJetsToLL_*HT*.txt
for i in $FILES; do  root -l -b -q "RunZPeakResolution.C(\"$i\", 4)";  done

##ZPt

echo "#define Y2016">IsData.h
FILES=files/UL2016/DYJetsToMuMu_M-50_Zpt-150toInf_TuneCP5_13TeV-madgraphMLM_pdfwgt_F-pythia8.txt
for i in $FILES; do  root -l -b -q "RunZPeakResolution.C(\"$i\", 2)";  done

echo "#define Y2017">IsData.h
FILES=files/UL2017/DYJetsToMuMu_M-50_Zpt-150toInf_TuneCP5_13TeV-madgraphMLM_pdfwgt_F-pythia8.txt
for i in $FILES; do  root -l -b -q "RunZPeakResolution.C(\"$i\", 2)";  done

echo "#define Y2018">IsData.h
FILES=files/UL2018/DYJetsToMuMu_M-50_Zpt-150toInf_TuneCP5_13TeV-madgraphMLM_pdfwgt_F-pythia8.txt
for i in $FILES; do  root -l -b -q "RunZPeakResolution.C(\"$i\", 2)";  done



# DATA

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