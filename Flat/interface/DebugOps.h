#include "Operator.h"

namespace pa {
    class PFDumpOp : public AnalysisOp {
        public:
            PFDumpOp(panda::Event& event_,
                      Config& cfg_,
                      Utils& utils_,
                      GeneralTree& gt_,
                      int level_=0) : 
            AnalysisOp("pfdump", event_, cfg_, utils_, gt_, level_) { }
            virtual bool on() { return true; }
        protected:
            void do_execute() {
                logger.debug("PFDumpOp::do_execute", 
                             Form("event = %llu\n", gt.eventNumber));
                int i = 0;
            }
    };
}
