# Panda Analysis

Full doxygen-generated documentation can be found [here](http://t3serv001.mit.edu/~snarayan/doxy/Panda).

Throughout this readme, I will refer to several user-defined environment variables. 
These are typically defined in `T3/setup.sh`.

## Installation

```bash
cmsrel CMSSW_9_4_6 # need at least this for Tensorflow C++ API
cd CMSSW_9_4_6/src
cmsenv
git clone https://github.com/PandaPhysics/PandaTree         # data format
git clone https://github.com/sidnarayanan/PandaCore         # core utilities 
git lfs clone https://github.com/sidnarayanan/PandaAnalysis # analysis package
scram b -j8
```

Most of your code development will be in `PandaAnalysis`. 

## Defining a data format

`PandaAnalysis` data formats are semi-flat (i.e. singleton or dynamic array) `TTree`s. 
They are defined using a metalanguage that is converted into C++. 
The language contains a set of keywords.
One can define two types of branches:
- singleton: `branch_name    type    [filled=X] [default=X]`
- array:     `branch_name    type[max_static_size]    [tree_size=X] [filled=X] [default=X] [shifts=X]`

These can be decorated with various keywords, either per-branch or in blocks (see below):
- shift:     `shift:shift_name    [opt1, opt2, ...]`
- conditional:    `cond:cond_name   condition [? true_val : false_val]` 
- constants:          `const:const_name val`

Note that keywords can nest, i.e one could have a conditional that refers to a previously-defined constant. 
There can be no whitespace in definitions!
Here is an example:
```
const:NJET  20
cond:hbb                   is_monohiggs||is_hbb
cond:njot                  (cond:hbb)?"nJotMax":"2"
shift:jetrings             [0,1,2,3,4]
# branches
runNumber                  int filled=cond:hbb
block:shifts=shift:jetrings,default=0
    dummyvar                      float[const:NJET]
endblock
```
The tree can then be generated by calling `Flat/bin/generateTreeClass.py --cfg Flat/config/someThing.cfg`.
This will create a class called `someThing`. 
If you would like to add any custom code to `someThing.[h,cc]`, do so between the `//STARTCUSTOM`...`//ENDCUSTOM` blocks. 
Most users will use the `GeneralTree` class, which is intended for most analyses.

## Defining an analysis
Here, I am going to talk about the structure of `PandaAnalyzer` and its modules, but there are examples of more task-specific analyses in the code. 
`PandaAnalyzer` reads in `panda` files and outputs `GeneralTree`.
`PandaAnalyzer` is essentially a container class that links together several things:
- An instance of `GeneralTree`
- An `Analysis` object, defined externally, which defines the analysis, the input, the output, etc
- `Config` and `Utils` instances, which contain, respectively, many configuration constants and useful utilities. These are initialized as a function of the aforementioned `Analysis` by `ConfigMod`
- Many `AnalysisMod`s, which do physics tasks. These can be chained together in a tree-like structure
- A `Registry`, which holds `shared_ptr`s to objects that should be shared between `AnalysisMods` (or really anywhere else in `PandaAnalyzer`). 

Here is what that looks like in a (clickable!) image:

![click me](http://t3serv001.mit.edu/~snarayan/doxy/PandaAnalysis/classpa_1_1PandaAnalyzer__coll__graph_org.svg) 

### Modules
The inheritance diagram for `Module`s looks like:
```
Module  # abstract class
  \--> ConfigMod
  \--> AnalysisMod = BaseAnalysisMod<GeneralTree> # also abstract
          \--> ContainerMod
          \--> GlobalMod
          \--> BaseJetMod
                 \--> JetMod
                 \--> FatJetMod
          \--> TFInferMod
                 \--> BRegDeepMod
          ...
```
Every realized subclass of `AnalysisMod` must define a `void do_execute()` function, which does the actual per-event execution.
Further virtual functions can be overridden as needed:
- Constructor and destructor, although it is preferred to use smart pointers whereever possible, so destructor should be avoided. 
- `void do_init(Registry& registry)`: call `registry.publish`(`publishConst`) and `registry.access`(`accessConst`) as needed to publish and and access data by name. 
- `void do_readData(TString path)`: read any data to initialize member objects
- `void do_reset()`: called every event before `do_execute()`
- `void do_terminate()`: called at the end of each file
- `bool on()`: whether the mod is on or not, typically a function of `this->analysis`. 

Each `BaseAnalysisMod` instance has a member function `addSubMod<MOD>()`, which adds a `MOD` instance to the calling mod. 
When the calling mod's `do_init` is called, this call is cascaded down to the sub-mods.
This is true for all `do_*` functions, EXCEPT `do_execute`. 
Calls to `do_execute` are left for the user to define in the calling mod's `do_execute` function. 
Note that the caller should actually call `someSubMod->execute()`, not `someSubMod->do_execute()`. 

Here is an example of a mod definition that shows most features:
```C++
class HbbSystemMod : public AnalysisMod {
public:
  HbbSystemMod(panda::EventAnalysis& event_,
               Config& cfg_,
               Utils& utils_,
               GeneralTree& gt_,
               int level_=0) :
    AnalysisMod("hbbsystem", event_, cfg_, utils_, gt_, level_),
    hbbdJet(std::make_shared<JetWrapper*>(nullptr)) {
      deepreg = addSubMod<BRegDeepMod>();
      bdtreg = addSubMod<BRegBDTMod>();
    }   
  virtual ~HbbSystemMod() { } 

  bool on() { return analysis.hbb; }
protected:
  void do_init(Registry& registry) {
    currentJES = registry.access<JESHandler*>("currentJES");
    looseLeps = registry.accessConst<std::vector<panda::Lepton*>>("looseLeps");
    dilep = registry.accessConst<TLorentzVector>("dilep");
    registry.publish("higgsDaughterJet", hbbdJet);
  }   
  void do_execute();
  void do_reset() { btagsorted.clear(); }
private:
  std::shared_ptr<JESHandler*> currentJES{nullptr};
  std::shared_ptr<const std::vector<panda::Lepton*>> looseLeps{nullptr};
  std::shared_ptr<const TLorentzVector> dilep{nullptr};

  std::vector<JetWrapper*> btagsorted;

  std::shared_ptr<JetWrapper*> hbbdJet;
  BRegDeepMod *deepreg{nullptr};
  BRegBDTMod *bdtreg{nullptr};
};  

```

### Chaining together an analysis
All `AnalysisMod`s are added to `PandaAnalyzer` in the constructor, even if your analysis does not want to run them.
Mods are turned on and off by `AnalysisMod::on()`.
This adds a bit of memory overhead, but it's small in the grand scale of things and allows for more run-time flexibility. 
One adds a mod by calling `ADDMOD(ModClass)`. 
The mods are executed in the order they are added.
Mods are called either before or after the preselections are checked (see below).
Depending on the needs of your mod, you may want to put things in one place or another.
Here is an example of what an analysis can look like:
```
INFO    [AnalysisMod::print            ]: -> global (0)
INFO    [AnalysisMod::print            ]: -> pre-sel (0)
INFO    [AnalysisMod::print            ]:       -> map (1)
INFO    [AnalysisMod::print            ]:       -> trigger (1)
INFO    [AnalysisMod::print            ]:       -> simplep (1)
INFO    [AnalysisMod::print            ]:       -> simplepho (1)
INFO    [AnalysisMod::print            ]:       -> recoil (1)
INFO    [AnalysisMod::print            ]:       -> fatjet (1)
INFO    [AnalysisMod::print            ]:       -> jet (1)
INFO    [AnalysisMod::print            ]:             -> jetflavor (2)
INFO    [AnalysisMod::print            ]:             -> isojet (2)
INFO    [AnalysisMod::print            ]:             -> bjetreg (2)
INFO    [AnalysisMod::print            ]:             -> vbfsystem (2)
INFO    [AnalysisMod::print            ]:             -> hbbsystem (2)
INFO    [AnalysisMod::print            ]:                   -> bregdeep (3)
INFO    [AnalysisMod::print            ]:       -> tau (1)
INFO    [AnalysisMod::print            ]: -> post-sel (0)
INFO    [AnalysisMod::print            ]:       -> hbbmisc (1)
INFO    [AnalysisMod::print            ]:       -> inclep (1)
INFO    [AnalysisMod::print            ]:       -> btagsf (1)
INFO    [AnalysisMod::print            ]:       -> btagweight (1)
```

### Pre-selections
The abstract base class for all selections is `Selection`.
All selections make decisions based on the `GeneralTree` instance (a const pointer is provided). 
There are two ways of defining your own selection:
- Define a sub-class of `Selection` that defines `bool do_accept() const`
- Define a sub-class of `LambdaSel` by passing function with signature `bool f(const GeneralTree*)` to the constructor

The latter can be done in python as well:
```python
from PandaAnalysis.Flat.selection import build 
triggersel = build(root.Selection.sReco, 'Trigger', '(gt->isData==0) || (gt->trigger!=0)', anded=True)
skimmer.AddPresel(triggersel)
```
Two things to note from this snippet:
- Pre-selections are added to `PandaAnalyzer` by the `AddPresel` method
- There are two types of pre-selection: `Selection::sReco`, which is called after a series of reconstruction-based mods, but before more computationally-intensive mods, and `Selection::sGen`, which is called at the very end and reduces disk size.

### Testing your analyzer

Inside `Flat/test`, there is a testing script that runs as:
```bash
./test.py /path/to/input/panda.root [DEBUG_LEVEL]
```
You have to open up the script and modify the number of events you want to run, the type of file it is (data, W+jets MC, etc), and what flags are on.


## Running an analysis
The analysis framework is heavily integrated with HTCondor resources at MIT. 
There are three running options:
- Using only the T3 (`export SUBMIT_CONFIG=T3`). This is used by `check` and `submit` since local filesystem access is needed.
- Using the T3 with opportunistic flocking to the T2 (`export SUBMIT_CONFIG=T2`). This is default for analysis jobs.
- Using the entire SubMIT grid (`export SUBMIT_CONFIG=SubMIT`). This is failure-prone due to the heterogenity of the grid, but with aggressive resubmission, can be useful for CPU-intensive tasks. Instructions for this are below.

For what follows, I assume you're in `$CMSSW_BASE/src/PandaAnalysis/T3`.

First, open up `T3/setup.sh`. 
This defines all the environment variables we'll need.
Make sure anything that is a path (`PANDA_FLATDIR`,`SUBMIT_LOGDIR`,etc) is writable by you.
The `SUBMIT*` environment variables are used for running the jobs, and the `PANDA*` variables are for running things on the outputs of the jobs.
`SUBMIT_TMPL` points to a script inside `inputs/`, that will define your analysis.
For a straightforward example, do `export SUBMIT_TMPL=tnp_tmpl.py`.

Now, let's open up `inputs/tnp_tmpl.py` as a concrete example and look at it.
The main function you have to worry about is `fn`.
Here is an example:
```python
def fn(input_name, isData, full_path):
    logger.info(sname+'.fn','Starting to process '+input_name)
    # now we instantiate and configure the analyzer
    a = monotop(True)
    a.recalcECF = True 
    a.mcTriggers = True
    a.inpath = input_name
    a.outpath = utils.input_to_output(input_name)
    a.datapath = data_dir
    a.isData = isData
    utils.set_year(a, 2016)
    a.processType = utils.classify_sample(full_path, isData)    

    skimmer = root.pa.PandaAnalyzer(a)
    skimmer.AddPresel(root.pa.FatJetSel())

    return utils.run_PandaAnalyzer(skimmer, isData, a.outpath)
```

### Cataloging inputs

```bash
./catalogT2Prod.py --outfile /path/to/config.cfg [ --include datasets to include ] [ --exclude datasets to skip ] [ --force ] [--smartcache]
```

I recommend you put the config in a web-accessible place for use in later steps. For example:
```bash
./catalogT2Prod.py --force --outfile ~/public_html/histcatalog/$(date +%Y%m%d).cfg --include TT --exclude TTbarDM --smartcache
```

The above command will do the following things:

- It will only check datasets that contain `TT` and do not contain `TTbarDM` in the dataset's nickname

- If a dataset is not found in `PandaCore.Tools.process`, it will guess the nickname of the dataset, give it a xsec of 1, and write it to the catalog (`--force`)

- If the file does not exist locally on the T3, a smartcache request will be made

- The output will be a timestamped web-facing config file


### Building the work environment

First, make sure the following environment variables are defined (some examples shown):
```bash
export PANDA_CFG="http://snarayan.web.cern.ch/snarayan/eoscatalog/20170127.cfg"  # location of config file from previous section
export SUBMIT_TMPL="skim_merge_tmpl.py"  # name of template script in T3/inputs
export SUBMIT_NAME="v_8024_2_0"  # name for this job
export SUBMIT_WORKDIR="/data/t3serv014/snarayan/condor/"${SUBMIT_NAME}"/work/"  # staging area for submission
export SUBMIT_LOGDIR="/data/t3serv014/snarayan/condor/"${SUBMIT_NAME}"/logs/"  # log directory
export SUBMIT_LOGDIR="/data/t3serv014/snarayan/condor/"${SUBMIT_NAME}"/locks/"  # lock directory
export SUBMIT_OUTDIR="/mnt/hadoop/scratch/snarayan/panda/"${SUBMIT_NAME}"/batch/"  # location of unmerged files
export PANDA_FLATDIR="${HOME}/home000/store/panda/v_8024_2_0/"   # merged output
export SUBMIT_CONFIG=T2  # allow running on T3 or T2. if $SUBMIT_CONFIG==T3, then only run on T3
```

`T3/inputs/$SUBMIT_TMPL` should be the skimming configuration you wish to run your files through. 

### Testing the jobs

If you want, you can test-run a job locally:
```bash
./task.py --build_only --nfiles 1   # just build the job, with one file per job
cd $SUBMIT_WORKDIR                  # go to the working directory
python skim.py 0 0                  # run the first job in the configuration
cd -                                # back to bin
./task.py --clean                   # clean up after yourself
```

### Submitting and re-submitting

To submit jobs, simply do
```bash 
./task.py --submit [--nfiles NFILES] [--clean]
```
where NFILES is the number of files in each job. 
The default is 25 files if that flag is not passed.
The clean argument will make sure to wipe out all staging directories to ensure a clean release (it's optional because sometimes you don't want to do this).

To check the status of your jobs, simply do:
```bash
./task.py --check [--silent] [--force] [--nfiles NFILES] [--monitor NSECONDS] [--submit_only]
```
Note that the above command overwrites `$SUBMIT_WORKDIR/local.cfg`, with the intention of preparing it for resubmission.
The file will be recreated as a configuration to rerun files that are not present in the output and not running.
- `--silent` will skip the per-sample breakdown.
- `--nfiles` will repackage `local.cfg` into a different number of files per job.
- `--force` will re-catalog files that are incomplete, not just missing.
- `--monitor` will capture the terminal screen and refresh the status after `NSECONDS` has elapsed since the last refresh.
- `--submit_only` will silently re-submit failed jobs if `--monitor` is turned on. 

To manually resubmit missing files, simply do
```bash
./task.py --submit
```
In the case that you are using the `--force` option, make sure you have no running jobs before resubmitting, or you may end up with duplicated outputs.
Forcing resubmission is generally discouraged and is largely included for historical reasons.
The job framework is robust enough at this point that forcing should not be necessary.

To be absolutely sure you have no duplicated outputs, you can run:
```bash
./task.py --check_duplicates
# ./task.py --clean_duplicates # same as above, but it also deletes the offenders and prepares for resubmission
```

If you are having lots of failures, it may be interesting to analyze the logs located in `$SUBMIT_LOGDIR`. 
You can do this manually, or:
```bash
./analyzeLogs.py [--dump]
```
This will print to screen a basic analysis of the errors observed, which errors are correlated, and how frequently they occur.
If the `--dump` flag is passed, then a directory `log_dumps/` is created, containing detailed information on each failure class (where it failed and on for what inputs).

### Running on SubMIT

The above steps still apply, modulo a single change: the `SUBMIT_*` directories must be on `/data/t3home000/` and have `777` permissions.


## Merging

Make sure `$PANDA_FLATDIR` exists. Then, go into `T3/merging` and do:
```bash
./merge.py [--cfg CONFIG] [--silent] TTbar_Powheg
```
to merge the Powheg TT sample, for example. 
If provided, `CONFIG` is the module that is imported from `configs/`. 
The default is `common`, but there are others, like `leptonic`.
To merge en-masse (e.g. many many signal outputs), you can do something like:
```bash
submit --exec merge.py --arglist list_of_signals.txt
```
This assumes that you have `PandaCore/bin` in your `$PATH`
The `submit` command will print a cache directory it created to keep track of the jobs.
You can check the status of the jobs by doing
```bash
check --cache <cache_directory> [--resubmit_failed]
```
The last flag is optional and will resubmit anything that failed (exited without code 0).
NB: the `submit` and `check` executables are very generic and don't know anything about the code they are running.
For them, "success" simply means "exited with code 0".
So it is important to check that the output looks sane to you.

