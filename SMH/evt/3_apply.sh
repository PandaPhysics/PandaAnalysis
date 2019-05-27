for i in /data/t3home000/hbb/zhnn/2018_v8/sr/*root; do

    python /home/bmaier/cms/Hbb/2018/CMSSW_10_3_1/src/PandaAnalysis/SMH/evt/infer_on_root.py  --ifile $i --json /home/bmaier/cms/Hbb/2018/CMSSW_10_3_1/src/PandaAnalysis/SMH/evt/root.json --model /home/bmaier/cms/Hbb/2018/SubTLEnet/train/smh/models/evt/2018_v8/weights.h5

done
