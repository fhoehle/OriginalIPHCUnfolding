#ifndef SAMPLELOCATIONS_HPP
#define SAMPLELOCATIONS_HPP

// WARNING there are still many other explicit filenames in the various code files, for example presel files. take care.

#include "usesystematics.hpp" // gives us whether to use systematics or not

const char* lumiString = "5.0";



#ifdef USESYSTEMATICS
// #include "samplelocations_sys_jesplus.hpp"
// #include "samplelocations_sys_jesminus.hpp"
// #include "samplelocations_sys_jrsplus.hpp"
// #include "samplelocations_sys_jrsminus.hpp"

//#include "samplelocations_sys_puplus.hpp"
//#include "samplelocations_sys_puminus.hpp"

//#include "samplelocations_sys_btagplus.hpp"
//#include "samplelocations_sys_btagminus.hpp"

//#include "samplelocations_sys_lepscaleplus.hpp"
//#include "samplelocations_sys_lepscaleminus.hpp"

//#include "samplelocations_sys_q2wplus.hpp"
//#include "samplelocations_sys_q2wminus.hpp"

//#include "samplelocations_sys_q2ttplus.hpp"
//#include "samplelocations_sys_q2ttminus.hpp"

#elif defined DOPDFSYS
  #include "samplelocations_sys_pdf.hpp"
#else

//const char* sl_data_comb = "/storage/8/froscher/outputroot/fall11/data_histos.root";
const char* sl_data_comb = "/afs/cern.ch/work/c/cardaci/Tuples_and_Results_24_May/TupleFileSkimmedDATA2011AB.root";
//const char* sl_data_comb = "/afs/cern.ch/work/c/cardaci/Tuples_and_Results_24_May/TupleFileSkimmed_TTbarAllMCatNLO_PileUpReweighted.root";
//const char* sl_data_comb = "/afs/cern.ch/work/c/cardaci/Tuples_and_Results_13_July_Generator_Systematics_woSC_nobtag_selection_newMWT/TupleFileSkimmed_TTPowheg_woSC_Pythia.root";
//const char* sl_data_comb = "/afs/cern.ch/work/c/cardaci/Tuples_and_Results_13_July_Generator_Systematics_woSC_nobtag_selection_newMWT/TupleFileSkimmed_TTMCatNLO_woSC_Herwig.root";
//const char* sl_data_comb = "/afs/cern.ch/work/c/cardaci/Tuples_and_Results_4_July_Generator_Systematics_no1btag_selection/TupleFileSkimmed_TTPowheg_Pythia.root";
//const char* sl_data_comb = "/afs/cern.ch/work/c/cardaci/Tuples_and_Results_4_July_Generator_Systematics_no1btag_selection/TupleFileSkimmed_TTMCatNLO_Herwig.root";
//const char* sl_data_comb = "/afs/cern.ch/work/c/cardaci/Tuples_and_Results_24_May/systematics/TupleFileSkimmed_TTbarSignalPowheg.root";
//const char* sl_data_comb = "/afs/cern.ch/work/c/cardaci/Tuples_and_Results_24_May/systematics/TupleFileSkimmed_TTbarSignalTuneZ2_7TeV-madgraph-tauola.root";
//const char* sl_data_comb = "/afs/cern.ch/work/c/cardaci/Tuples_and_Results_24_May/systematics/TupleFileSkimmed_TTbarSignalTuneZ2_7TeV-madgraph-tauola.root";

/////////////const char* sl_ttbarele = "/storage/8/froscher/outputroot/fall11/ttbar_ele_histos.root";
/////////////const char* sl_ttbarele = "ttbar_mu_histos.root";

const char* sl_ttbarele = "/afs/cern.ch/work/c/cardaci/Tuples_and_Results_24_May/TupleFileSkimmed_TTbarAllMCatNLO_PileUpReweighted.root";
//const char* sl_ttbarele = "/afs/cern.ch/work/c/cardaci/Tuples_and_Results_24_May/TupleFileSkimmed_TTbarAllMCatNLO.root";
//const char* sl_ttbarele = "/afs/cern.ch/work/c/cardaci/Tuples_and_Results_24_May/systematics/TupleFileSkimmed_TTbarSignalTuneZ2_7TeV-madgraph-tauola_PileUpReweighted.root";
//const char* sl_ttbarele = "/afs/cern.ch/work/c/cardaci/Tuples_and_Results_24_May/TupleFileSkimmed_TTbarSignalPowheg.root";
//const char* sl_ttbarele = "/afs/cern.ch/work/c/cardaci/Tuples_and_Results_24_May/LeptonAsymmetry/TupleFileSkimmed_TTbarAllMCatNLO_PileUpReweighted_LeptonAsymmetry.root";



//const char* sl_ttbarele = "/storage/8/froscher/outputroot/fall11/ttbar_madgraph_ele_histos.root";
// const char* sl_ttbarele = "/storage/8/froscher/outputroot/fall11/ttbar_mcatnlo_ele_histos.root ";
//const char* sl_ttbarele_mg = "/storage/8/froscher/outputroot/fall11/ttbar_madgraph_ele_histos.root";
//const char* sl_wjetscombele = "/storage/8/froscher/outputroot/fall11/wijets_ele.root";
////////////////const char* sl_wjetsposele = "/storage/8/froscher/outputroot/fall11/wijets_pos_ele.root";

const char* sl_wjetsposele = "/afs/cern.ch/work/c/cardaci/Tuples_and_Results_24_May/TupleFileSkimmed_TTbarSemileptonicMCatNLO_PileUpReweighted.root";
//const char* sl_wjetsposele = "/afs/cern.ch/work/c/cardaci/Tuples_and_Results_24_May/TupleFileSkimmed_TTbarSemileptonicMCatNLO.root";
//const char* sl_wjetsposele = "/afs/cern.ch/work/c/cardaci/Tuples_and_Results_24_May/TupleFileSkimmed_TTbarSemileptonicMadgraph_PileUpReweighted.root";


//const char* sl_wjetsnegele = "/storage/8/froscher/outputroot/fall11/wijets_neg_ele.root";
//const char* sl_wjetscombele = "/storage/8/froscher/outputroot/fall11/wjets_ele_histos.root";
//const char* sl_wjetsposele = "/storage/8/froscher/outputroot/fall11/wjets_pos_ele_histos.root";
//const char* sl_wjetsnegele = "/storage/8/froscher/outputroot/fall11/wjets_neg_ele_histos.root";
////////////const char* sl_zjetsele = "/storage/8/froscher/outputroot/fall11/zjets_ele_histos.root";
const char* sl_zjetsele = "/afs/cern.ch/work/c/cardaci/Tuples_and_Results_24_May/TupleFileSkimmed_ZJets_PileUpReweighted.root";
const char* sl_dibosonele = "/afs/cern.ch/work/c/cardaci/Tuples_and_Results_24_May/TupleFileSkimmed_DiBosons_PileUpReweighted.root";
//const char* sl_zjetsele = "/afs/cern.ch/work/c/cardaci/Tuples_and_Results_24_May/TupleFileSkimmed_ZJets.root";
//const char* sl_dibosonele = "/afs/cern.ch/work/c/cardaci/Tuples_and_Results_24_May/TupleFileSkimmed_DiBosons.root";
//const char* sl_singlet_tchan_t_ele = "/storage/8/froscher/outputroot/fall11/singlet_tchan_top_ele_histos.root";
//const char* sl_singlet_tchan_tbar_ele = "/storage/8/froscher/outputroot/fall11/singlet_tchan_tbar_ele_histos.root";
////////const char* sl_singlet_twchan_t_ele = "/storage/8/froscher/outputroot/fall11/singlet_twchan_top_ele_histos.root";
////////const char* sl_singlet_twchan_tbar_ele = "/storage/8/froscher/outputroot/fall11/singlet_twchan_tbar_ele_histos.root";

const char* sl_singlet_twchan_t_ele = "/afs/cern.ch/work/c/cardaci/Tuples_and_Results_24_May/TupleFileSkimmed_SingleTtW_PileUpReweighted.root";
const char* sl_singlet_twchan_tbar_ele = "/afs/cern.ch/work/c/cardaci/Tuples_and_Results_24_May/TupleFileSkimmed_SingleTbartW_PileUpReweighted.root";
//const char* sl_singlet_twchan_t_ele = "/afs/cern.ch/work/c/cardaci/Tuples_and_Results_24_May/TupleFileSkimmed_SingleTtW.root";
//const char* sl_singlet_twchan_tbar_ele = "/afs/cern.ch/work/c/cardaci/Tuples_and_Results_24_May/TupleFileSkimmed_SingleTbartW.root";

//const char* sl_qcdele = "/storage/8/froscher/outputroot/fall11/qcd_ele_histos.root";
//const char* sl_dataele = "/storage/8/froscher/outputroot/fall11/data_ele_histos.root";




/////////////const char* sl_ttbarmu = "/storage/8/froscher/outputroot/fall11/ttbar_mu_histos.root";
////////////const char* sl_ttbarmu = "ttbar_mu_histos.root";
//const char* sl_ttbarmu = "TupleFileSkimmed.root";


//const char* sl_ttbarmu = "/storage/8/froscher/outputroot/fall11/ttbar_madgraph_mu_histos.root";
// const char* sl_ttbarmu = "/storage/8/froscher/outputroot/fall11/ttbar_mcatnlo_mu_histos.root";

/*
const char* sl_ttbarmu_mg = "/storage/8/froscher/outputroot/fall11/ttbar_madgraph_mu_histos.root";
const char* sl_wjetscombmu = "/storage/8/froscher/outputroot/fall11/wijets_mu.root";
const char* sl_wjetsposmu = "/storage/8/froscher/outputroot/fall11/wijets_pos_mu.root";
const char* sl_wjetsnegmu = "/storage/8/froscher/outputroot/fall11/wijets_neg_mu.root";
*/
//const char* sl_wjetscombmu = "/storage/8/froscher/outputroot/fall11/wjets_mu_histos.root";
//const char* sl_wjetsposmu = "/storage/8/froscher/outputroot/fall11/wjets_pos_mu_histos.root";
//const char* sl_wjetsnegmu = "/storage/8/froscher/outputroot/fall11/wjets_neg_mu_histos.root";
/*
const char* sl_zjetsmu = "/storage/8/froscher/outputroot/fall11/zjets_mu_histos.root";
const char* sl_singlet_tchan_t_mu = "/storage/8/froscher/outputroot/fall11/singlet_tchan_top_mu_histos.root";
const char* sl_singlet_tchan_tbar_mu = "/storage/8/froscher/outputroot/fall11/singlet_tchan_tbar_mu_histos.root";
const char* sl_singlet_twchan_t_mu = "/storage/8/froscher/outputroot/fall11/singlet_twchan_top_mu_histos.root";
const char* sl_singlet_twchan_tbar_mu = "/storage/8/froscher/outputroot/fall11/singlet_twchan_tbar_mu_histos.root";

const char* sl_qcdmu = "/storage/8/froscher/outputroot/fall11/qcd_mu_histos.root";
const char* sl_datamu = "/storage/8/froscher/outputroot/fall11/data_mu_histos.root";
*/
// const char* sl_ttbarcomb = "/storage/8/froscher/outputroot/fall11/ttbar_histos.root";
// const char* sl_ttbarcomb = "/storage/8/froscher/outputroot/fall11/ttbar_mcatnlo_histos.root ";

  
/////////////////const char* sl_preselfile0 = "/storage/8/froscher/outputroot/fall11/pre_ttbar0.root";
//const char* sl_preselfile1 = "/storage/8/froscher/outputroot/fall11/pre_ttbar1.root";
//const char* sl_preselfile2 = "/storage/8/froscher/outputroot/fall11/pre_ttbar2.root";
//const char* sl_preselfile3 = "/storage/8/froscher/outputroot/fall11/pre_ttbar3.root";
//const char* sl_preselfile4 = "/storage/8/froscher/outputroot/fall11/pre_ttbar4.root";
//const char* sl_preselfile5 = "/storage/8/froscher/outputroot/fall11/pre_ttbar5.root";
//const char* sl_preselfile6 = "/storage/8/froscher/outputroot/fall11/pre_ttbar6.root";
//const char* sl_preselfile7 = "/storage/8/froscher/outputroot/fall11/pre_ttbar7.root";
//const char* sl_preselfile8 = "/storage/8/froscher/outputroot/fall11/pre_ttbar8.root";
//const char* sl_preselfile9 = "/storage/8/froscher/outputroot/fall11/pre_ttbar9.root";
//const char* sl_preselfile10 = "/storage/8/froscher/outputroot/fall11/pre_ttbar10.root";
//const char* sl_preselfile11 = "/storage/8/froscher/outputroot/fall11/pre_ttbar11.root";
//const char* sl_preselfile12 = "/storage/8/froscher/outputroot/fall11/pre_ttbar12.root";
//const char* sl_preselfile13 = "/storage/8/froscher/outputroot/fall11/pre_ttbar13.root";
//const char* sl_preselfile14 = "/storage/8/froscher/outputroot/fall11/pre_ttbar14.root";
//const char* sl_preselfile15 = "/storage/8/froscher/outputroot/fall11/pre_ttbar15.root";
//////////////////const char* sl_preselfile16 = "/storage/8/froscher/outputroot/fall11/pre_ttbar16.root";

/////////////const char* sl_preselfile0 = "ttbar_mu_histos.root";
//const char* sl_preselfile1 = "ttbar_mu_histos.root";
//const char* sl_preselfile2 = "ttbar_mu_histos.root";
//const char* sl_preselfile3 = "ttbar_mu_histos.root";
//const char* sl_preselfile4 = "ttbar_mu_histos.root";
//const char* sl_preselfile5 = "ttbar_mu_histos.root";
//const char* sl_preselfile6 = "ttbar_mu_histos.root";
//const char* sl_preselfile7 = "ttbar_mu_histos.root";
//const char* sl_preselfile8 = "ttbar_mu_histos.root";
//const char* sl_preselfile9 = "ttbar_mu_histos.root";
//const char* sl_preselfile10 = "ttbar_mu_histos.root";
//const char* sl_preselfile11 = "ttbar_mu_histos.root";
//const char* sl_preselfile12 = "ttbar_mu_histos.root";
//const char* sl_preselfile13 = "ttbar_mu_histos.root";
//const char* sl_preselfile14 = "ttbar_mu_histos.root";
//const char* sl_preselfile15 = "ttbar_mu_histos.root";
///////////const char* sl_preselfile16 = "ttbar_mu_histos.root";



const char* sl_preselfile0 = "/afs/cern.ch/work/c/cardaci/Tuples_and_Results_24_May/TupleFileNonSkimmed_TTbarAllMCatNLO.root";

//const char* sl_preselfile0 = "/afs/cern.ch/work/c/cardaci/precut_merged_nlo.root";
//const char* sl_preselfile0 = "/afs/cern.ch/work/c/cardaci/Tuples_and_Results_24_May/LeptonAsymmetry/TupleFileNonSkimmed_TTbarAllMCatNLO_LeptonAsymmetry.root";



//const char* sl_preselfile0 = "/afs/cern.ch/work/c/cardaci/precut_merged_nlo.root";

//const char* sl_preselfile1 = "TupleFileNonSkimmed.root";
//const char* sl_preselfile2 = "TupleFileNonSkimmed.root";
//const char* sl_preselfile3 = "TupleFileNonSkimmed.root";
//const char* sl_preselfile4 = "TupleFileNonSkimmed.root";
//const char* sl_preselfile5 = "TupleFileNonSkimmed.root";
//const char* sl_preselfile6 = "TupleFileNonSkimmed.root";
//const char* sl_preselfile7 = "TupleFileNonSkimmed.root";
//const char* sl_preselfile8 = "TupleFileNonSkimmed.root";
//const char* sl_preselfile9 = "TupleFileNonSkimmed.root";
//const char* sl_preselfile10 = "TupleFileNonSkimmed.root";
//const char* sl_preselfile11 = "TupleFileNonSkimmed.root";
//const char* sl_preselfile12 = "TupleFileNonSkimmed.root";
//const char* sl_preselfile13 = "TupleFileNonSkimmed.root";
//const char* sl_preselfile14 = "TupleFileNonSkimmed.root";
//const char* sl_preselfile15 = "TupleFileNonSkimmed.root";
//const char* sl_preselfile16 = "TupleFileNonSkimmed.root";









// const char* sl_preselfile0 = "/storage/8/froscher/outputroot/fall11/precut_0_pre_ttbar_mg.root";
// const char* sl_preselfile1 = "/storage/8/froscher/outputroot/fall11/precut_1_pre_ttbar_mg.root";
// const char* sl_preselfile2 = "/storage/8/froscher/outputroot/fall11/precut_2_pre_ttbar_mg.root";
// const char* sl_preselfile3 = "/storage/8/froscher/outputroot/fall11/precut_3_pre_ttbar_mg.root";
// const char* sl_preselfile4 = "/storage/8/froscher/outputroot/fall11/precut_4_pre_ttbar_mg.root";
// const char* sl_preselfile5 = "/storage/8/froscher/outputroot/fall11/precut_5_pre_ttbar_mg.root";
// const char* sl_preselfile6 = "/storage/8/froscher/outputroot/fall11/precut_6_pre_ttbar_mg.root";
// const char* sl_preselfile7 = "/storage/8/froscher/outputroot/fall11/precut_7_pre_ttbar_mg.root";
// const char* sl_preselfile8 = "/storage/8/froscher/outputroot/fall11/precut_8_pre_ttbar_mg.root";
// const char* sl_preselfile9 = "/storage/8/froscher/outputroot/fall11/precut_9_pre_ttbar_mg.root";
// const char* sl_preselfile10 = "/storage/8/froscher/outputroot/fall11/precut_10_pre_ttbar_mg.root";
// const char* sl_preselfile11 = "/storage/8/froscher/outputroot/fall11/precut_11_pre_ttbar_mg.root";
// const char* sl_preselfile12 = "/storage/8/froscher/outputroot/fall11/precut_12_pre_ttbar_mg.root";
// const char* sl_preselfile13 = "/storage/8/froscher/outputroot/fall11/precut_13_pre_ttbar_mg.root";
// const char* sl_preselfile14 = "/storage/8/froscher/outputroot/fall11/precut_14_pre_ttbar_mg.root";
// const char* sl_preselfile15 = "/storage/8/froscher/outputroot/fall11/precut_15_pre_ttbar_mg.root";
// const char* sl_preselfile16 = "/storage/8/froscher/outputroot/fall11/precut_16_pre_ttbar_mg.root";


// const char* sl_preselfile0 = "/storage/8/froscher/outputroot/fall11/Precut_TTbar_MCatNLO_Fall11/precut_0_pre_ttbar_nlo.root";
// const char* sl_preselfile1 = "/storage/8/froscher/outputroot/fall11/Precut_TTbar_MCatNLO_Fall11/precut_1_pre_ttbar_nlo.root";
// const char* sl_preselfile2 = "/storage/8/froscher/outputroot/fall11/Precut_TTbar_MCatNLO_Fall11/precut_2_pre_ttbar_nlo.root";
// const char* sl_preselfile3 = "/storage/8/froscher/outputroot/fall11/Precut_TTbar_MCatNLO_Fall11/precut_3_pre_ttbar_nlo.root";
// const char* sl_preselfile4 = "/storage/8/froscher/outputroot/fall11/Precut_TTbar_MCatNLO_Fall11/precut_4_pre_ttbar_nlo.root";
// const char* sl_preselfile5 = "/storage/8/froscher/outputroot/fall11/Precut_TTbar_MCatNLO_Fall11/precut_5_pre_ttbar_nlo.root";
// const char* sl_preselfile6 = "/storage/8/froscher/outputroot/fall11/Precut_TTbar_MCatNLO_Fall11/precut_6_pre_ttbar_nlo.root";
// const char* sl_preselfile7 = "/storage/8/froscher/outputroot/fall11/Precut_TTbar_MCatNLO_Fall11/precut_7_pre_ttbar_nlo.root";
// const char* sl_preselfile8 = "/storage/8/froscher/outputroot/fall11/Precut_TTbar_MCatNLO_Fall11/precut_8_pre_ttbar_nlo.root";
// const char* sl_preselfile9 = "/storage/8/froscher/outputroot/fall11/Precut_TTbar_MCatNLO_Fall11/precut_9_pre_ttbar_nlo.root";
// const char* sl_preselfile10 = "/storage/8/froscher/outputroot/fall11/Precut_TTbar_MCatNLO_Fall11/precut_10_pre_ttbar_nlo.root";
// const char* sl_preselfile11 = "/storage/8/froscher/outputroot/fall11/Precut_TTbar_MCatNLO_Fall11/precut_11_pre_ttbar_nlo.root";
// const char* sl_preselfile12 = "/storage/8/froscher/outputroot/fall11/Precut_TTbar_MCatNLO_Fall11/precut_12_pre_ttbar_nlo.root";
// const char* sl_preselfile13 = "/storage/8/froscher/outputroot/fall11/Precut_TTbar_MCatNLO_Fall11/precut_13_pre_ttbar_nlo.root";
// const char* sl_preselfile14 = "/storage/8/froscher/outputroot/fall11/Precut_TTbar_MCatNLO_Fall11/precut_14_pre_ttbar_nlo.root";
// const char* sl_preselfile15 = "/storage/8/froscher/outputroot/fall11/Precut_TTbar_MCatNLO_Fall11/precut_15_pre_ttbar_nlo.root";
// const char* sl_preselfile16 = "/storage/8/froscher/outputroot/fall11/Precut_TTbar_MCatNLO_Fall11/precut_16_pre_ttbar_nlo.root";

#endif
  
#endif
