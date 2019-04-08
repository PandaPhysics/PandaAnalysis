#echo TT SingleTop ZtoNuNu WJets VH Diboson | xargs -n1 -P5 python ./convert.py --out /data/t3home000/hbb/zhnn/v12/sr/train/ --json root.json --name

# 2018 SR
#echo TT SingleTop ZtoNuNu WJets VH Diboson | xargs -n1 -P5 python ./convert.py --out /data/t3home000/hbb/zhnn/2018_v1/sr/train/ --json root.json --name


# 2018 CRs
echo TT SingleTop ZtoNuNu WJets Diboson | xargs -n1 -P5 python ./convert.py --out /data/t3home000/hbb/zhnn/2018_v7/cr_zhf/ --json root_zhf.json --name
echo TT SingleTop ZtoNuNu WJets Diboson | xargs -n1 -P5 python ./convert.py --out /data/t3home000/hbb/zhnn/2018_v7/cr_zlf/ --json root_zlf.json --name
echo TT SingleTop ZtoNuNu WJets Diboson | xargs -n1 -P5 python ./convert.py --out /data/t3home000/hbb/zhnn/2018_v7/cr_ttbar/ --json root_ttbar.json --name
