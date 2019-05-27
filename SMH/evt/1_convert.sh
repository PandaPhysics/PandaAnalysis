# 2018 SR
#echo TT SingleTop ZtoNuNu WJets VH Diboson | xargs -n1 -P5 python ./convert.py --out /data/t3home000/hbb/zhnn/2018_v8/sr/train/ --json root.json --name

# 2018 CRs
#echo TT SingleTop ZtoNuNu WJets Diboson | xargs -n1 -P5 python ./convert.py --out /data/t3home000/hbb/zhnn/2018_v8/cr_zhf/ --json root_zhf.json --name
#echo TT SingleTop ZtoNuNu WJets Diboson | xargs -n1 -P5 python ./convert.py --out /data/t3home000/hbb/zhnn/2018_v8/cr_zlf/ --json root_zlf.json --name
#echo TT SingleTop ZtoNuNu WJets Diboson | xargs -n1 -P5 python ./convert.py --out /data/t3home000/hbb/zhnn/2018_v8/cr_ttbar/ --json root_ttbar.json --name

# 2017 SR
echo TT SingleTop ZtoNuNu WJets VH Diboson | xargs -n1 -P5 python ./convert.py --out /data/t3home000/hbb/zhnn/v12/sr/train/with_fj --json root_with_fj_2017.json --name
