#include "MUSK_CH.h"

#include "utils_math.h"
#include "text.h"
#include "ChannelRoutingCommon.h"

using namespace utils_math;

MUSK_CH::MUSK_CH() :
    m_dt(-1), m_inputSubbsnID(-1), m_nreach(-1), m_outletID(-1),
    m_Epch(NODATA_VALUE), m_Bnk0(NODATA_VALUE), m_Chs0_perc(NODATA_VALUE),
    m_aBank(NODATA_VALUE), m_bBank(NODATA_VALUE), m_subbsnID(nullptr),
    m_mskX(NODATA_VALUE), m_mskCoef1(NODATA_VALUE), m_mskCoef2(NODATA_VALUE),
    m_chWth(nullptr), m_chDepth(nullptr), m_chLen(nullptr), m_chArea(nullptr),
    m_chSideSlope(nullptr), m_chSlope(nullptr), m_chMan(nullptr),
    m_Kchb(nullptr), m_Kbank(nullptr), m_reachDownStream(nullptr),
    // Inputs from other modules
    m_petSubbsn(nullptr), m_gwSto(nullptr),
    m_olQ2Rch(nullptr), m_ifluQ2Rch(nullptr), m_gndQ2Rch(nullptr),
    // Temporary variables
    m_ptSub(nullptr), m_flowIn(nullptr), m_flowOut(nullptr), m_seepage(nullptr),
    // Outputs
    m_qRchOut(nullptr), m_qsRchOut(nullptr), m_qiRchOut(nullptr), m_qgRchOut(nullptr),
    m_chSto(nullptr), m_rteWtrIn(nullptr), m_rteWtrOut(nullptr), m_bankSto(nullptr),
    m_chWtrDepth(nullptr), m_chWtrWth(nullptr), m_chBtmWth(nullptr), m_chCrossArea(nullptr),
    //ljj++
    m_GWMAX(NODATA_VALUE),m_Kg(NODATA_VALUE), m_Base_ex(NODATA_VALUE),m_ispermafrost(nullptr),
    gw_height(nullptr), m_rch_ht(nullptr), m_gw_sh(nullptr),m_qgsRchOut(nullptr),
    m_nCells(-1),m_subbasinsInfo(nullptr),m_prec(nullptr),curBasinArea(nullptr),m_area(nullptr),m_netPcp(nullptr),
    m_islake(nullptr),m_lakearea(nullptr),m_evlake(NODATA_VALUE),m_lakeseep(NODATA_VALUE),m_petFactor(NODATA_VALUE),
    m_minvol(NODATA_VALUE),m_lakedpini(nullptr),m_lakedp(nullptr),m_lakevol(nullptr),
    m_qin1(nullptr),m_qout1(nullptr),m_lakealpha(nullptr),m_natural_flow(nullptr),
    m_isres(nullptr),m_ResLc(nullptr),m_ResLn(nullptr),m_ResLf(nullptr),m_ResAdjust(nullptr)
{
}

MUSK_CH::~MUSK_CH() {
    /// 1. reaches related variables will be released in ~clsReaches(). By lj, 2017-12-26.
    /// 2. m_ptSrcFactory will be released by DataCenter->Scenario. lj

    if (nullptr != m_ptSub) Release1DArray(m_ptSub);
    if (nullptr != m_flowIn) Release1DArray(m_flowIn);
    if (nullptr != m_flowOut) Release1DArray(m_flowOut);
    if (nullptr != m_seepage) Release1DArray(m_seepage);

    if (nullptr != m_qRchOut) Release1DArray(m_qRchOut);
    if (nullptr != m_qsRchOut) Release1DArray(m_qsRchOut);
    if (nullptr != m_qiRchOut) Release1DArray(m_qiRchOut);
    if (nullptr != m_qgRchOut) Release1DArray(m_qgRchOut);

    if (nullptr != m_chSto) Release1DArray(m_chSto);
    if (nullptr != m_rteWtrIn) Release1DArray(m_rteWtrIn);
    if (nullptr != m_rteWtrOut) Release1DArray(m_rteWtrOut);
    if (nullptr != m_bankSto) Release1DArray(m_bankSto);

    if (nullptr != m_chWtrDepth) Release1DArray(m_chWtrDepth);
    if (nullptr != m_chWtrWth) Release1DArray(m_chWtrWth);
    if (nullptr != m_chBtmWth) Release1DArray(m_chBtmWth);
    if (nullptr != m_chCrossArea) Release1DArray(m_chCrossArea);

    //ljj++
    if (nullptr != m_rch_ht) Release1DArray(m_rch_ht);
    if (nullptr != m_qgsRchOut) Release1DArray(m_qgsRchOut);
    if (nullptr != curBasinArea) Release1DArray(curBasinArea);
    if (nullptr != m_lakedp) Release1DArray(m_lakedp);
    if (nullptr != m_qin1) Release1DArray(m_qin1);
    if (nullptr != m_qout1) Release1DArray(m_qout1);
}

bool MUSK_CH::CheckInputData() {
    CHECK_POSITIVE(MID_MUSK_CH, m_dt);
    CHECK_NONNEGATIVE(MID_MUSK_CH, m_inputSubbsnID);
    CHECK_POSITIVE(MID_MUSK_CH, m_nreach);
    CHECK_POSITIVE(MID_MUSK_CH, m_outletID);
    CHECK_NODATA(MID_MUSK_CH, m_Epch);
    CHECK_NODATA(MID_MUSK_CH, m_Bnk0);
    CHECK_NODATA(MID_MUSK_CH, m_Chs0_perc);
    CHECK_NODATA(MID_MUSK_CH, m_aBank);
    CHECK_NODATA(MID_MUSK_CH, m_bBank);
    //CHECK_NODATA(MID_MUSK_CH, m_mskX); // Do not throw exception, since they have default values.
    //CHECK_NODATA(MID_MUSK_CH, m_mskCoef1);
    //CHECK_NODATA(MID_MUSK_CH, m_mskCoef2);
    CHECK_POINTER(MID_MUSK_CH, m_subbsnID);

    CHECK_POINTER(MID_MUSK_CH, m_petSubbsn);
    CHECK_POINTER(MID_MUSK_CH, m_gwSto);
    CHECK_POINTER(MID_MUSK_CH, m_olQ2Rch);
    CHECK_POINTER(MID_MUSK_CH, m_ifluQ2Rch);
    CHECK_POINTER(MID_MUSK_CH, m_gndQ2Rch);
    return true;
}

void MUSK_CH::InitialOutputs() {
    CHECK_POSITIVE(MID_MUSK_CH, m_nreach);
    if (nullptr != m_qRchOut) return; // DO NOT Initial Outputs repeatedly.
    if (m_mskX < 0.f) m_mskX = 0.2f;
    if (m_mskCoef1 < 0.f || m_mskCoef1 > 1.f) {
        m_mskCoef1 = 0.75f;
        m_mskCoef2 = 0.25f;
    } else {
        // There is no need to use mskCoef2 as input parameter.
        // Make sure m_mskCoef1 + m_mskCoef2 = 1.
        //float msk1 = m_mskCoef1 / (m_mskCoef1 + m_mskCoef2);
        //float msk2 = m_mskCoef2 / (m_mskCoef1 + m_mskCoef2);
        //m_mskCoef1 = msk1;
        //m_mskCoef2 = msk2;
    }
    m_mskCoef2 = 1.f - m_mskCoef1;

    m_flowIn = new(nothrow) float[m_nreach + 1];
    m_flowOut = new(nothrow) float[m_nreach + 1];
    m_seepage = new(nothrow) float[m_nreach + 1];

    m_qRchOut = new(nothrow) float[m_nreach + 1];
    m_qsRchOut = new(nothrow) float[m_nreach + 1];
    m_qiRchOut = new(nothrow) float[m_nreach + 1];
    m_qgRchOut = new(nothrow) float[m_nreach + 1];

    m_chSto = new(nothrow) float[m_nreach + 1];
    m_rteWtrIn = new(nothrow) float[m_nreach + 1];
    m_rteWtrOut = new(nothrow) float[m_nreach + 1];
    m_bankSto = new(nothrow) float[m_nreach + 1];

    m_chWtrDepth = new(nothrow) float[m_nreach + 1];
    m_chWtrWth = new(nothrow) float[m_nreach + 1];
    m_chBtmWth = new(nothrow) float[m_nreach + 1];
    m_chCrossArea = new(nothrow) float[m_nreach + 1];

    //ljj++
    m_rch_ht = new(nothrow) float[m_nreach + 1];
    m_qgsRchOut = new(nothrow) float[m_nreach + 1];
    curBasinArea = new(nothrow) float[m_nreach + 1];
    m_prec = new(nothrow) float[m_nreach + 1];
    m_lakedp = new(nothrow) float[m_nreach + 1];
    m_qin1 = new(nothrow) float[m_nreach + 1];
    m_qout1 = new(nothrow) float[m_nreach + 1];
    for (int i = 1; i <= m_nreach; i++) {
        m_qRchOut[i] = m_olQ2Rch[i];
        m_qsRchOut[i] = m_olQ2Rch[i];
        if (nullptr != m_ifluQ2Rch) {
            m_qRchOut[i] += m_ifluQ2Rch[i];
            m_qiRchOut[i] = m_ifluQ2Rch[i];
        } else {
            m_qiRchOut[i] = 0.f;
        }
        if (nullptr != m_gndQ2Rch) {
            m_qRchOut[i] += m_gndQ2Rch[i];
            m_qgRchOut[i] = m_gndQ2Rch[i];
        } else {
            m_qgRchOut[i] = 0.f;
        }
        m_seepage[i] = 0.f;
        m_bankSto[i] = m_Bnk0 * m_chLen[i];
        m_chBtmWth[i] = ChannleBottomWidth(m_chWth[i], m_chSideSlope[i], m_chDepth[i]);
        m_chCrossArea[i] = ChannelCrossSectionalArea(m_chBtmWth[i], m_chDepth[i], m_chSideSlope[i]);
        m_chWtrDepth[i] = m_chDepth[i] * m_Chs0_perc;
        m_chWtrWth[i] = m_chBtmWth[i] + 2.f * m_chSideSlope[i] * m_chWtrDepth[i];
        m_chSto[i] = m_chLen[i] * m_chWtrDepth[i] * (m_chBtmWth[i] + m_chSideSlope[i] * m_chWtrDepth[i]);
        m_flowIn[i] = m_chSto[i];
        m_flowOut[i] = m_chSto[i];
        m_rteWtrIn[i] = 0.f;
        m_rteWtrOut[i] = 0.f;
        //ljj++
        m_rch_ht[i] = m_chWtrDepth[i];
        m_lakedp[i] = m_lakedpini[i]; //初值改了
        curBasinArea[i] = 0.f;
        m_prec[i] = 0.f;
        m_qin1[i]=0.f;
        m_qout1[i]=0.f;
        if(m_islake[i]==1) m_chSto[i] = m_lakedpini[i]*m_lakearea[i];
        if(m_isres[i]==1) m_chSto[i] = m_lakevol[i];
    }
    /// initialize point source loadings
    if (nullptr == m_ptSub) {
        Initialize1DArray(m_nreach + 1, m_ptSub, 0.f);
    }
}

void MUSK_CH::PointSourceLoading() {
    /// load point source water discharge (m3/s) on current day from Scenario
    for (auto it = m_ptSrcFactory.begin(); it != m_ptSrcFactory.end(); ++it) {
        /// reset point source loading water to 0.f
        for (int i = 0; i <= m_nreach; i++) {
            m_ptSub[i] = 0.f;
        }
        //cout<<"unique Point Source Factory ID: "<<it->first<<endl;
        vector<int>& ptSrcMgtSeqs = it->second->GetPointSrcMgtSeqs();
        map<int, PointSourceMgtParams *>& pointSrcMgtMap = it->second->GetPointSrcMgtMap();
        vector<int>& ptSrcIDs = it->second->GetPointSrcIDs();
        map<int, PointSourceLocations *>& pointSrcLocsMap = it->second->GetPointSrcLocsMap();
        // 1. looking for management operations from m_pointSrcMgtMap
        for (auto seqIter = ptSrcMgtSeqs.begin(); seqIter != ptSrcMgtSeqs.end(); ++seqIter) {
            PointSourceMgtParams* curPtMgt = pointSrcMgtMap.at(*seqIter);
            // 1.1 If current day is beyond the date range, then continue to next management
            if (curPtMgt->GetStartDate() != 0 && curPtMgt->GetEndDate() != 0) {
                if (m_date < curPtMgt->GetStartDate() || m_date > curPtMgt->GetEndDate()) {
                    continue;
                }
            }
            // 1.2 Otherwise, get the water volume
            float per_wtrVol = curPtMgt->GetWaterVolume(); /// m3/'size'/day
            // 1.3 Sum up all point sources
            for (auto locIter = ptSrcIDs.begin(); locIter != ptSrcIDs.end(); ++locIter) {
                if (pointSrcLocsMap.find(*locIter) != pointSrcLocsMap.end()) {
                    PointSourceLocations* curPtLoc = pointSrcLocsMap.at(*locIter);
                    int curSubID = curPtLoc->GetSubbasinID();
                    m_ptSub[curSubID] += per_wtrVol * curPtLoc->GetSize() / 86400.f; /// m3/'size'/day ==> m3/s
                }
            }
        }
    }
}

int MUSK_CH::Execute() {
    InitialOutputs();
    /// load point source water volume from m_ptSrcFactory
    PointSourceLoading();
    //ljj++
    for (auto id = m_subbasinIDs.begin(); id != m_subbasinIDs.end(); ++id) {
        Subbasin* sub = m_subbasinsInfo->GetSubbasinByID(*id);
        int curCellsNum = sub->GetCellCount();
        int* curCells = sub->GetCells();
        m_prec[*id] = 0.f;
        
        float total_area=0.f;
        if(m_islake[*id] == 1 || m_isres[*id] == 1){
            for (int i = 0; i < curCellsNum; i++) {
                int index = curCells[i];
                total_area += m_area[index];
                m_prec[*id]+= m_netPcp[index]/1000* m_area[index]; //m3
                curBasinArea[*id] += m_area[index];
            }
            m_prec[*id] = m_prec[*id] / total_area * m_lakearea[*id]/ m_dt;  //m3/s
        }
    }
    for (auto it = m_rteLyrs.begin(); it != m_rteLyrs.end(); ++it) {
        // There are not any flow relationship within each routing layer.
        // So parallelization can be done here.
        int reachNum = CVT_INT(it->second.size());
        size_t errCount = 0;
        // the size of m_rteLyrs (map) is equal to the maximum stream order
#pragma omp parallel for reduction(+:errCount)
        for (int i = 0; i < reachNum; i++) {
            int reachIndex = it->second[i]; // index in the array, i.e., subbasinID
            if (m_inputSubbsnID == 0 || m_inputSubbsnID == reachIndex) {
                // for OpenMP version, all reaches will be executed,
                // for MPI version, only the current reach will be executed.
                if(m_islake[reachIndex] == 1){ 
                    if (!LakeBudget(reachIndex)) {
                        errCount++;
                    }
                }
                else if (m_isres[reachIndex] == 1){
                    if (!ResBudget(reachIndex)) {
                      errCount++;
                    }
                }
                else{
                    if (!ChannelFlow(reachIndex)) {
                      errCount++;
                    }
                }
            }
        }
        if (errCount > 0) {
            throw ModelException(MID_MUSK_CH, "Execute", "Error occurred!");
        }
    }
    return 0;
}

void MUSK_CH::SetValue(const char* key, const float value) {
    string sk(key);
    if (StringMatch(sk, Tag_ChannelTimeStep)) m_dt = CVT_INT(value);
    else if (StringMatch(sk, Tag_SubbasinId)) m_inputSubbsnID = CVT_INT(value);
    else if (StringMatch(sk, VAR_OUTLETID)) m_outletID = CVT_INT(value);
    else if (StringMatch(sk, VAR_EP_CH)) m_Epch = value;
    else if (StringMatch(sk, VAR_BNK0)) m_Bnk0 = value;
    else if (StringMatch(sk, VAR_CHS0_PERC)) m_Chs0_perc = value;
    else if (StringMatch(sk, VAR_A_BNK)) m_aBank = value;
    else if (StringMatch(sk, VAR_B_BNK)) m_bBank = value;
    else if (StringMatch(sk, VAR_MSK_X)) m_mskX = value;
    else if (StringMatch(sk, VAR_MSK_CO1)) m_mskCoef1 = value;
    //ljj++
    else if (StringMatch(sk, VAR_GWMAX)) m_GWMAX = value;
    else if (StringMatch(sk, VAR_KG)) m_Kg = value;
    else if (StringMatch(sk, VAR_Base_ex)) m_Base_ex = value;
    else if (StringMatch(sk, VAR_LAKE_EVP)) m_evlake = value;
    else if (StringMatch(sk, VAR_LAKE_SEEP)) m_lakeseep = value;
    else if (StringMatch(sk, VAR_K_PET)) m_petFactor = value;
    else if (StringMatch(sk, VAR_LAKE_MNVOL))  m_minvol = value;
    else {
        throw ModelException(MID_MUSK_CH, "SetValue", "Parameter " + sk + " does not exist.");
    }
}

void MUSK_CH::SetValueByIndex(const char* key, const int index, const float value) {
    if (m_inputSubbsnID == 0) return;           // Not for omp version
    if (index <= 0 || index > m_nreach) return; // index should belong 1 ~ m_nreach
    if (nullptr == m_qRchOut) InitialOutputs();
    string sk(key);
    /// Set single value of array1D of current subbasin
    /// IN/OUTPUT variables
    if (StringMatch(sk, VAR_QRECH)) m_qRchOut[index] = value;
    else if (StringMatch(sk, VAR_QS)) m_qsRchOut[index] = value;
    else if (StringMatch(sk, VAR_QI)) m_qiRchOut[index] = value;
    else if (StringMatch(sk, VAR_QG)) m_qgRchOut[index] = value;
    else if (StringMatch(sk, VAR_QGS)) m_qgsRchOut[index] = value;
    else {
        throw ModelException(MID_MUSK_CH, "SetValueByIndex", "Parameter " + sk + " does not exist.");
    }
}

void MUSK_CH::Set1DData(const char* key, const int n, float* data) {
    string sk(key);
    //check the input data
    if (StringMatch(sk, VAR_SUBBSN)) {
        m_subbsnID = data;
    } else if (StringMatch(sk, VAR_SBPET)) {
        CheckInputSize(MID_MUSK_CH, key, n - 1, m_nreach);
        m_petSubbsn = data;
    } else if (StringMatch(sk, VAR_SBGS)) {
        CheckInputSize(MID_MUSK_CH, key, n - 1, m_nreach);
        m_gwSto = data;
    } else if (StringMatch(sk, VAR_SBOF)) {
        CheckInputSize(MID_MUSK_CH, key, n - 1, m_nreach);
        m_olQ2Rch = data;
    } else if (StringMatch(sk, VAR_SBIF)) {
        CheckInputSize(MID_MUSK_CH, key, n - 1, m_nreach);
        m_ifluQ2Rch = data;
    } else if (StringMatch(sk, VAR_SBQG)) {
        CheckInputSize(MID_MUSK_CH, key, n - 1, m_nreach);
        m_gndQ2Rch = data;
    }
    //ljj++
    else if (StringMatch(sk, VAR_GWH)) {
        CheckInputSize(MID_MUSK_CH, key, n - 1, m_nreach);
        gw_height = data;
    }
    else if (StringMatch(sk, VAR_GW_SH)) {
        CheckInputSize(MID_MUSK_CH, key, n - 1, m_nreach);
        m_gw_sh = data;
    } else if (StringMatch(sk, VAR_PCP)) {
        CheckInputSize(MID_MUSK_CH, key, n, m_nCells);
        m_netPcp = data;
    }
    else if (StringMatch(sk, VAR_AHRU)) {
        CheckInputSize(MID_MUSK_CH, key, n, m_nCells);
        m_area = data;
    }
    else {
        throw ModelException(MID_MUSK_CH, "Set1DData", "Parameter " + sk + " does not exist.");
    }
}

void MUSK_CH::GetValue(const char* key, float* value) {
    InitialOutputs();
    string sk(key);
    /// IN/OUTPUT variables
    if (StringMatch(sk, VAR_QRECH) && m_inputSubbsnID > 0) *value = m_qRchOut[m_inputSubbsnID];
    else if (StringMatch(sk, VAR_QS) && m_inputSubbsnID > 0) *value = m_qsRchOut[m_inputSubbsnID];
    else if (StringMatch(sk, VAR_QI) && m_inputSubbsnID > 0) *value = m_qiRchOut[m_inputSubbsnID];
    else if (StringMatch(sk, VAR_QG) && m_inputSubbsnID > 0) *value = m_qgRchOut[m_inputSubbsnID];
     else if (StringMatch(sk, VAR_QGS) && m_inputSubbsnID > 0) *value = m_qgsRchOut[m_inputSubbsnID];
    else {
        throw ModelException(MID_MUSK_CH, "GetValue", "Parameter " + sk + " does not exist.");
    }
}

void MUSK_CH::Get1DData(const char* key, int* n, float** data) {
    InitialOutputs();
    string sk(key);
    *n = m_nreach + 1;
    if (StringMatch(sk, VAR_QRECH)) {
        m_qRchOut[0] = m_qRchOut[m_outletID];
        *data = m_qRchOut;
    } else if (StringMatch(sk, VAR_QS)) {
        m_qsRchOut[0] = m_qsRchOut[m_outletID];
        *data = m_qsRchOut;
    } else if (StringMatch(sk, VAR_QI)) {
        m_qiRchOut[0] = m_qiRchOut[m_outletID];
        *data = m_qiRchOut;
    } else if (StringMatch(sk, VAR_QG)) {
        m_qgRchOut[0] = m_qgRchOut[m_outletID];
        *data = m_qgRchOut;
    } else if (StringMatch(sk, VAR_CHST)) {
        m_chSto[0] = m_chSto[m_outletID];
        *data = m_chSto;
    } else if (StringMatch(sk, VAR_RTE_WTRIN)) {
        m_rteWtrIn[0] = m_rteWtrIn[m_outletID];
        *data = m_rteWtrIn;
    } else if (StringMatch(sk, VAR_RTE_WTROUT)) {
        m_rteWtrOut[0] = m_rteWtrOut[m_outletID];
        *data = m_rteWtrOut;
    } else if (StringMatch(sk, VAR_BKST)) {
        m_bankSto[0] = m_bankSto[m_outletID];
        *data = m_bankSto;
    } else if (StringMatch(sk, VAR_CHWTRDEPTH)) {
        m_chWtrDepth[0] = m_chWtrDepth[m_outletID];
        *data = m_chWtrDepth;
    } else if (StringMatch(sk, VAR_CHWTRWIDTH)) {
        m_chWtrWth[0] = m_chWtrWth[m_outletID];
        *data = m_chWtrWth;
    } else if (StringMatch(sk, VAR_CHBTMWIDTH)) {
        m_chBtmWth[0] = m_chBtmWth[m_outletID];
        *data = m_chBtmWth;
    } else if (StringMatch(sk, VAR_CHCROSSAREA)) {
        m_chCrossArea[0] = m_chCrossArea[m_outletID];
        *data = m_chCrossArea;
    } 
    //ljj++
    else if (StringMatch(sk, VAR_CHSEEPAGE)) {
        m_seepage[0] = m_seepage[m_outletID];
        *data = m_seepage;
    }else if (StringMatch(sk, VAR_QGS)) {
        m_qgsRchOut[0] = m_qgsRchOut[m_outletID];
        *data = m_qgsRchOut;
    }
    else {
        throw ModelException(MID_MUSK_CH, "Get1DData", "Output " + sk + " does not exist.");
    }
}

void MUSK_CH::SetScenario(Scenario* sce) {
    if (nullptr != sce) {
        map<int, BMPFactory *>& tmpBMPFactories = sce->GetBMPFactories();
        for (auto it = tmpBMPFactories.begin(); it != tmpBMPFactories.end(); ++it) {
            /// Key is uniqueBMPID, which is calculated by BMP_ID * 100000 + subScenario;
            if (it->first / 100000 == BMP_TYPE_POINTSOURCE) {
#ifdef HAS_VARIADIC_TEMPLATES
                m_ptSrcFactory.emplace(it->first, static_cast<BMPPointSrcFactory*>(it->second));
#else
                m_ptSrcFactory.insert(make_pair(it->first, static_cast<BMPPointSrcFactory*>(it->second)));
#endif
            }
        }
    } else {
        throw ModelException(MID_MUSK_CH, "SetScenario", "The scenario can not to be NULL.");
    }
}

void MUSK_CH::SetReaches(clsReaches* reaches) {
    if (nullptr == reaches) {
        throw ModelException(MID_MUSK_CH, "SetReaches", "The reaches input can not to be NULL.");
    }
    m_nreach = reaches->GetReachNumber();

    if (nullptr == m_chWth) reaches->GetReachesSingleProperty(REACH_WIDTH, &m_chWth);
    if (nullptr == m_chDepth) reaches->GetReachesSingleProperty(REACH_DEPTH, &m_chDepth);
    if (nullptr == m_chLen) reaches->GetReachesSingleProperty(REACH_LENGTH, &m_chLen);
    if (nullptr == m_chArea) reaches->GetReachesSingleProperty(REACH_AREA, &m_chArea);
    if (nullptr == m_chSideSlope) reaches->GetReachesSingleProperty(REACH_SIDESLP, &m_chSideSlope);
    if (nullptr == m_chSlope) reaches->GetReachesSingleProperty(REACH_SLOPE, &m_chSlope);
    if (nullptr == m_chMan) reaches->GetReachesSingleProperty(REACH_MANNING, &m_chMan);
    if (nullptr == m_Kbank) reaches->GetReachesSingleProperty(REACH_BNKK, &m_Kbank);
    if (nullptr == m_Kchb) reaches->GetReachesSingleProperty(REACH_BEDK, &m_Kchb);
    if (nullptr == m_reachDownStream) reaches->GetReachesSingleProperty(REACH_DOWNSTREAM, &m_reachDownStream);
    if (nullptr == m_ispermafrost) reaches->GetReachesSingleProperty(REACH_PERMAFORST, &m_ispermafrost);
    //ljj++
    if (nullptr == m_islake) reaches->GetReachesSingleProperty(REACH_ISLAKE, &m_islake);
    if (nullptr == m_lakearea) reaches->GetReachesSingleProperty(REACH_LAKEAREA, &m_lakearea);
    if (nullptr == m_lakevol) reaches->GetReachesSingleProperty(REACH_LAKEVOL, &m_lakevol);
    if (nullptr == m_lakedpini) reaches->GetReachesSingleProperty(REACH_LAKEDPINI, &m_lakedpini);
    if (nullptr == m_lakealpha) reaches->GetReachesSingleProperty(REACH_LAKEALPHA, &m_lakealpha);
    if (nullptr == m_isres) reaches->GetReachesSingleProperty(REACH_ISRES, &m_isres);
    if (nullptr == m_natural_flow) reaches->GetReachesSingleProperty(REACH_NATURAL_FLOW, &m_natural_flow);
    if (nullptr == m_ResLc) reaches->GetReachesSingleProperty(REACH_RES_LC, &m_ResLc);
    if (nullptr == m_ResLn) reaches->GetReachesSingleProperty(REACH_RES_LN, &m_ResLn);
    if (nullptr == m_ResLf) reaches->GetReachesSingleProperty(REACH_RES_LF, &m_ResLf);
    if (nullptr == m_ResAdjust) reaches->GetReachesSingleProperty(REACH_RES_ADJUST, &m_ResAdjust);

    m_reachUpStream = reaches->GetUpStreamIDs();
    m_rteLyrs = reaches->GetReachLayers();
}

bool MUSK_CH::ChannelFlow(const int i) {
    // 1. first add all the inflow water
    float qIn = 0.f; /// Water entering reach on current day from both current subbasin and upstreams
    // 1.1. water from this subbasin
    qIn += m_olQ2Rch[i]; /// surface flow
    float qiSub = 0.f;   /// interflow flow
    if (nullptr != m_ifluQ2Rch && m_ifluQ2Rch[i] >= 0.f) {
        qiSub = m_ifluQ2Rch[i];
        qIn += qiSub;
    }
    float qgSub = 0.f; /// groundwater flow
    if (nullptr != m_gndQ2Rch && m_gndQ2Rch[i] >= 0.f) {
        qgSub = m_gndQ2Rch[i];
        qIn += qgSub;
    }
    float ptSub = 0.f; /// point sources flow
    if (nullptr != m_ptSub && m_ptSub[i] >= 0.f) {
        ptSub = m_ptSub[i];
        qIn += ptSub;
    }
    float qgsSub = 0.f; /// shallow groundwater flow
    // 1.2. water from upstream reaches
    float qsUp = 0.f;
    float qiUp = 0.f;
    float qgUp = 0.f;
    float qgsUp = 0.f;
    for (auto upRchID = m_reachUpStream.at(i).begin(); upRchID != m_reachUpStream.at(i).end(); ++upRchID) {
        if (m_qsRchOut[*upRchID] != m_qsRchOut[*upRchID]) {
            cout << "DayOfYear: " << m_dayOfYear << ", rchID: " << i << ", upRchID: " << *upRchID <<
                    ", surface part illegal!" << endl;
            return false;
        }
        if (m_qiRchOut[*upRchID] != m_qiRchOut[*upRchID]) {
            cout << "DayOfYear: " << m_dayOfYear << ", rchID: " << i << ", upRchID: " << *upRchID <<
                    ", subsurface part illegal!" << endl;
            return false;
        }
        if (m_qgRchOut[*upRchID] != m_qgRchOut[*upRchID]) {
            cout << "DayOfYear: " << m_dayOfYear << ", rchID: " << i << ", upRchID: " << *upRchID <<
                    ", groundwater part illegal!" << endl;
            return false;
        }
        if (m_qsRchOut[*upRchID] > 0.f) qsUp += m_qsRchOut[*upRchID];
        if (m_qiRchOut[*upRchID] > 0.f) qiUp += m_qiRchOut[*upRchID];
        if (m_qgRchOut[*upRchID] > 0.f) qgUp += m_qgRchOut[*upRchID];
        if (m_qgsRchOut[*upRchID] > 0.f) qgsUp += m_qgsRchOut[*upRchID];
    }
    qIn += qsUp + qiUp + qgUp + qgsUp;
#ifdef PRINT_DEBUG
    cout << "ID: " << i << ", surfaceQ: " << m_qsSub[i] << ", subsurfaceQ: " << qiSub <<
        ", groundQ: " << qgSub << ", pointQ: " << ptSub <<
        ", UPsurfaceQ: " << qsUp << ", UPsubsurface: " << qiUp << ", UPground: " << qgUp << endl;
#endif
    // 1.3. water from bank storage
    float bankOut = m_bankSto[i] * (1.f - exp(-m_aBank));
    m_bankSto[i] -= bankOut;
    qIn += bankOut / m_dt;

    // loss the water from bank storage to the adjacent unsaturated zone and groundwater storage
    float bankOutGw = m_bankSto[i] * (1.f - exp(-m_bBank));
    m_bankSto[i] -= bankOutGw;
    if (nullptr != m_gwSto) {
        m_gwSto[i] += bankOutGw / m_chArea[i] * 1000.f; // updated groundwater storage
    }

    // Compute storage time constant (ratio of storage to discharge) for reach
    // Wetting perimeter at bankfull
    float wet_perimeter = ChannelWettingPerimeter(m_chBtmWth[i], m_chDepth[i], m_chSideSlope[i]);
    // Cross-sectional area at bankfull
    float cross_area = ChannelCrossSectionalArea(m_chBtmWth[i], m_chDepth[i], m_chSideSlope[i]);
    // Hydraulic radius
    float radius = cross_area / wet_perimeter;
    // The storage time constant calculated for the reach segment with bankfull flows.
    float k_bankfull = StorageTimeConstant(m_chMan[i], m_chSlope[i], m_chLen[i], radius); // Hour

    // The storage time constant calculated for the reach segment with one-tenth of the bankfull flows.
    float wet_perimeter2 = ChannelWettingPerimeter(m_chBtmWth[i], 0.1f * m_chDepth[i], m_chSideSlope[i]);
    float cross_area2 = ChannelCrossSectionalArea(m_chBtmWth[i], 0.1f * m_chDepth[i], m_chSideSlope[i]);
    float radius2 = cross_area2 / wet_perimeter2;
    float k_bankfull2 = StorageTimeConstant(m_chMan[i], m_chSlope[i], m_chLen[i], radius2); // Hour

    float xkm = k_bankfull * m_mskCoef1 + k_bankfull2 * m_mskCoef2;
    // Eq. 7:1.4.9 in SWAT Theory 2009.
    // Check Muskingum numerical stability
    float detmax = 2.f * xkm * (1.f - m_mskX);
    float detmin = 2.f * xkm * m_mskX;
    // Discretize time interval to meet the stability criterion
    float det = 24.f; // hours, time step
    int nn = 0;       // number of subdaily computation points for stable routing
    if (det > detmax) {
        if (det / 2.f <= detmax) {
            det = 12.f;
            nn = 2;
        } else if (det / 4.f <= detmax) {
            det = 6.f;
            nn = 4;
        } else {
            det = 1;
            nn = 24;
        }
    } else {
        det = 24;
        nn = 1;
    }
    // get coefficients of Muskingum
    float temp = detmax + det;
    float c1 = (det - detmin) / temp;
    float c2 = (det + detmin) / temp;
    float c3 = (detmax - det) / temp;
    // make sure any coefficient is positive. Not sure whether this is needed.
    //if (c1 < 0) {
    //    c2 += c1;
    //    c1 = 0.f;
    //}
    //if (c3 < 0) {
    //    c2 += c1;
    //    c3 = 0.f;
    //}

#ifdef PRINT_DEBUG
    cout << " chStorage before routing " << m_chStorage[i] << endl;
#endif
    m_rteWtrOut[i] = qIn * m_dt;   // m^3
    float wtrin = qIn * m_dt / nn; // Inflow during a sub time interval, m^3
    float vol = 0.f;               // volume of water in reach, m^3
    float volrt = 0.f;             // flow rate, m^3/s
    float max_rate = 0.f;          // maximum flow capacity of the channel at bank full (m^3/s)
    float sdti = 0.f;              // average flow on day in reach, m^3/s, i.e., m_qRchOut
    float vc = 0.f;                // average flow velocity in channel, m/s
    float rchp = 0.f;              // wet perimeter, m
    float rcharea = 0.f;           // cross-sectional area, m^2
    float rchradius = 0.f;         // hydraulic radius
    float rtwtr = 0.f;             // water leaving reach on day, m^3, i.e., m_rteWtrOut
    float rttlc = 0.f;             // transmission losses from reach on day, m^3
    float qinday = 0.f;            // m^3
    float qoutday = 0.f;           // m^3
    //ljj++
    // Calculate River to hband
    float rchareaa = (m_chSto[i]+qIn * m_dt)/m_chLen[i];
    m_rch_ht[i] = sqrt(rchareaa/m_chSideSlope[i] + pow(m_chBtmWth[i]/(2*m_chSideSlope[i]),2))-m_chBtmWth[i]/(2*m_chSideSlope[i]); //m     
    float he = gw_height[i] - m_GWMAX*0.001f; //head of Ele GW
    float dh = m_rch_ht[i]-he;
    if (m_ispermafrost[i]==1) dh = -1.f*gw_height[i];   //河道深切峡谷
    if (abs(dh) <= 1e-6) dh = 0.f;
    float K = m_Kbank[i];
    float Q = 0.f;   //flux  River to Element
    //Q = 1.f*dh* K * m_chArea[i];//m3

    float qsep = 0.f;
    float rto = m_gw_sh[i]/gw_height[i];
    //the height from topsoil to gw bottom
    float gw_max = 0.f; 
    //for permafrost, the groundwater is very shallow
    if (m_ispermafrost[i]==1) gw_max = 2.5f*rto;
    if (dh > 0.f){
        Q = m_Kg * pow(abs(dh), m_Base_ex) * rto; // m3
        qsep = Min(Q, m_chSto[i]+qIn * m_dt);   //m3
        qgUp -= qsep/m_dt*(qgUp/(qIn+UTIL_ZERO));
        qsUp -= qsep/m_dt*(qsUp/(qIn+UTIL_ZERO));
        qiUp -= qsep/m_dt*(qiUp/(qIn+UTIL_ZERO));
        qgSub -= qsep/m_dt*(qgSub/(qIn+UTIL_ZERO));
        m_olQ2Rch[i] -= qsep/m_dt*(m_olQ2Rch[i]/(qIn+UTIL_ZERO));
        qiSub -= qsep/m_dt*(qiSub/(qIn+UTIL_ZERO));
        qIn  -= qsep/m_dt;
        if(qIn<=0){
            m_chSto[i] +=qIn/m_dt;
            qIn = 0.f;
            qgUp = 0.f;
            qsUp = 0.f;
            qiUp = 0.f;
            qgSub = 0.f;
            m_olQ2Rch[i] = 0.f;
            qiSub = 0.f;
            m_chSto[i] = Max(m_chSto[i],0.f);
        }
    }else{  
        Q = -1.f *m_Kg* pow(abs(dh), m_Base_ex) * rto; // m3
        qsep = Max(Q, -1.f*m_gw_sh[i]);   //m3
        if(m_gw_sh[i]<=UTIL_ZERO) qsep = 0.f;
       //qsep = Q;
        qgsSub -= qsep/m_dt;
        qIn  -= qsep/m_dt;
    }
    float gw_sh_o = m_gw_sh[i];
    if (nullptr != m_gw_sh) {
        m_gw_sh[i] += qsep; // m3
        m_gw_sh[i] = Max(m_gw_sh[i],0.f);
    }
    //if higher than gwmax, water could flow out immediately
    // if (m_ispermafrost[i]==1) {
    //     if(m_gw_sh[i]>=gw_max){
    //         float ul_water = m_gw_sh[i] - gw_max;
    //         qgsSub += ul_water/m_dt;
    //         qIn  += ul_water/m_dt;
    //         m_gw_sh[i] = gw_max;
    //     }
    // }
    m_seepage[i] = m_gw_sh[i] - gw_sh_o;
    m_rteWtrOut[i] = qIn * m_dt;
    wtrin = qIn * m_dt / nn;
    //if (i==11) cout<<dh<<"  "<<m_seepage[i]<<"  "<< gw_height[i]<<endl;

    // Iterate for the day
    for (int ii = 0; ii < nn; ii++) {
        // Calculate volume of water in reach
        vol = m_chSto[i] + wtrin; // m^3
        // Find average flowrate in a sub time interval, m^3/s
        volrt = vol / (86400.f / nn);
        // Find maximum flow capacity of the channel at bank full, m^3/s
        max_rate = manningQ(cross_area, radius, m_chMan[i], m_chSlope[i]);
        sdti = 0.f;
        m_chWtrDepth[i] = 0.f;
        // If average flowrate is greater than the channel capacity at bank full,
        //   then simulate flood plain flow, else simulate the regular channel flow.
        if (volrt > max_rate) {
            m_chWtrDepth[i] = m_chDepth[i];
            sdti = max_rate;
            //a = b * d + chsslope * d * d
            rcharea = m_chBtmWth[i] * m_chDepth[i] + m_chSideSlope[i] * m_chDepth[i] * m_chDepth[i];
            rchp = m_chDepth[i];
            float p = m_chBtmWth[i] + 2. * m_chDepth[i] * sqrt(1. + m_chSideSlope[i] * m_chSideSlope[i]);
            // Find the cross-sectional area and depth for volrt by iteration method at 1cm interval depth.
            // Find the depth until the discharge rate is equal to volrt
            while (sdti < volrt) {
                m_chWtrDepth[i] += 0.01f; // Increase 1cm at each interation
                float addarea = rcharea + ((m_chWth[i] * 5) + 4 * m_chWtrDepth[i]) * m_chWtrDepth[i];
                float addp = p + (m_chWth[i] * 4) + 2. * m_chWtrDepth[i] * sqrt(1. + 4 * 4);
                radius = addarea / addp;
                rcharea = addarea;
	            rchp = addp;
                // rcharea = ChannelCrossSectionalArea(m_chBtmWth[i], m_chDepth[i], m_chWtrDepth[i],
                //                                     m_chSideSlope[i], m_chWth[i], 4.f);
                // rchp = ChannelWettingPerimeter(m_chBtmWth[i], m_chDepth[i], m_chWtrDepth[i],
                //                                m_chSideSlope[i], m_chWth[i], 4.f);
                // radius = rcharea / rchp;
                sdti = manningQ(rcharea, radius, m_chMan[i], m_chSlope[i]);
            }
            sdti = volrt;
        } else {
            // Find the cross-sectional area and depth for volrt by iteration method at 1cm interval depth
            // Find the depth until the discharge rate is equal to volrt.
            while (sdti < volrt) {
                m_chWtrDepth[i] += 0.01f;
                rcharea = (m_chBtmWth[i] + m_chSideSlope[i] * m_chWtrDepth[i]) * m_chWtrDepth[i];
                rchp = m_chBtmWth[i] + 2. * m_chWtrDepth[i] * sqrt(1. + m_chSideSlope[i] * m_chSideSlope[i]);
                // rcharea = ChannelCrossSectionalArea(m_chBtmWth[i], m_chWtrDepth[i], m_chSideSlope[i]);
                // rchp = ChannelWettingPerimeter(m_chBtmWth[i], m_chWtrDepth[i], m_chSideSlope[i]);
                rchradius = rcharea / rchp;
                sdti = manningQ(rcharea, rchradius, m_chMan[i], m_chSlope[i]);
            }
            sdti = volrt;
        }
        // Calculate top width of channel at water level, topw in SWAT
        if (m_chWtrDepth[i] <= m_chDepth[i]) {
            m_chWtrWth[i] = m_chBtmWth[i] + 2.f * m_chWtrDepth[i] * m_chSideSlope[i];
        } else {
            m_chWtrWth[i] = 5.f * m_chWth[i] + 2.f * (m_chWtrDepth[i] - m_chDepth[i]) * 4.f;
        }
        if (sdti > 0.f) {
            // Calculate velocity and travel time
            vc = sdti / rcharea;                       // vel_chan(:) in SWAT
            float rttime = m_chLen[i] / (3600.f * vc); // reach travel time, hr
            // Compute water leaving reach on day
            rtwtr = c1 * wtrin + c2 * m_flowIn[i] + c3 * m_flowOut[i];
            if (rtwtr < 0.f) rtwtr = 0.f;
            rtwtr = Min(rtwtr, wtrin + m_chSto[i]);
            // Calculate amount of water in channel at end of day
            m_chSto[i] += wtrin - rtwtr;
            // Add if statement to keep m_chStorage from becoming negative
            if (m_chSto[i] < 0.f) m_chSto[i] = 0.f;

            // Transmission and evaporation losses are proportionally taken from the channel storage
            //   and from volume flowing out
            if (rtwtr > 0.f) {
                // Total time in hours to clear the water
                rttlc = det * m_Kchb[i] * 0.001f * m_chLen[i] * rchp; // m^3
                float rttlc2 = rttlc * m_chSto[i] / (rtwtr + m_chSto[i]);
                float rttlc1 = 0.f;
                if (m_chSto[i] <= rttlc2) {
                    rttlc2 = Min(rttlc2, m_chSto[i]);
                }
                m_chSto[i] -= rttlc2;
                rttlc1 = rttlc - rttlc2;
                if (rtwtr <= rttlc1) {
                    rttlc1 = Min(rttlc1, rtwtr);
                }
                rtwtr -= rttlc1;
                rttlc = rttlc1 + rttlc2; // Total water loss by transmission
            }
            // Calculate evaporation
            float rtevp = 0.f;
            float rtevp1 = 0.f;
            float rtevp2 = 0.f;
            if (rtwtr > 0.f) {
                /// In SWAT source code, line 306 of rtmusk.f, I think aaa should be divided by nn! By lj.
                float aaa = m_Epch * m_petSubbsn[i] * 0.001f / nn; // m
                if (m_chWtrDepth[i] <= m_chDepth[i]) {
                    rtevp = aaa * m_chLen[i] * m_chWtrWth[i]; // m^3
                } else {
                    if (aaa <= m_chWtrDepth[i] - m_chDepth[i]) {
                        rtevp = aaa * m_chLen[i] * m_chWtrWth[i];
                    } else {
                        rtevp = aaa;
                        m_chWtrWth[i] = m_chBtmWth[i] + 2.f * m_chDepth[i] * m_chSideSlope[i];
                        rtevp *= m_chLen[i] * m_chWtrWth[i]; // m^3
                    }
                }
                rtevp2 = rtevp * m_chSto[i] / (rtwtr + m_chSto[i]);
                if (m_chSto[i] <= rtevp2) {
                    rtevp2 = Min(rtevp2, m_chSto[i]);
                }
                m_chSto[i] -= rtevp2;
                rtevp1 = rtevp - rtevp2;
                if (rtwtr <= rtevp1) {
                    rtevp1 = Min(rtevp1, rtwtr);
                }
                rtwtr -= rtevp1;
                rtevp = rtevp1 + rtevp2; // Total water loss by evaporation
            }
            // Define flow parameters for current iteration
            m_flowIn[i] = wtrin;
            m_flowOut[i] = rtwtr;
            // Define flow parameters for current day
            qinday += wtrin;
            qoutday += rtwtr;
            // Total outflow for the day
            rtwtr = qoutday;
        } else {
            rtwtr = 0.f;
            sdti = 0.f;
            m_chSto[i] = 0.f;
            m_flowIn[i] = 0.f;
            m_flowOut[i] = 0.f;
        }
    } /* Iterate for the day */
    if (rtwtr < 0.f) rtwtr = 0.f;
    if (m_chSto[i] < 0.f) m_chSto[i] = 0.f;
    if (m_chSto[i] < 10.f) {
        rtwtr += m_chSto[i];
        m_chSto[i] = 0.f;
    }
    m_qRchOut[i] = sdti;
    m_rteWtrOut[i] = rtwtr;
    m_chCrossArea[i] = rcharea;
    m_qRchOut[i] = rtwtr / m_dt;
    if(m_qRchOut[i]<=UTIL_ZERO) m_rteWtrOut[i] = 0.f;

    float qInSum = m_olQ2Rch[i] + qiSub + qgSub + qsUp + qiUp + qgUp +qgsSub + qgsUp;
    if (qInSum < UTIL_ZERO) {
        // In case of divided by zero.
        m_qsRchOut[i] = 0.f;
        m_qiRchOut[i] = 0.f;
        m_qgRchOut[i] = 0.f;
        m_qgsRchOut[i] = 0.f;
        m_qRchOut[i] = 0.f;
    } else {
        // In my opinion, these lines should use `qIn` instead of `qInSum`. By lj.
        m_qsRchOut[i] = m_qRchOut[i] * (m_olQ2Rch[i] + qsUp) / (qIn+UTIL_ZERO);
        m_qiRchOut[i] = m_qRchOut[i] * (qiSub + qiUp) / (qIn+UTIL_ZERO);
        m_qgRchOut[i] = m_qRchOut[i] * (qgSub + qgUp) / (qIn+UTIL_ZERO);
        m_qgsRchOut[i] = m_qRchOut[i] * (qgsSub + qgsUp) / (qIn+UTIL_ZERO);
    }

    // Add transmission losses to bank storage/deep aquifer (i.e., groundwater in current version)
    if (rttlc > 0.f) {
        float trnsrch = 0.5f;
        if (rchp > 0.f) {
            trnsrch = m_chBtmWth[i] / rchp; // Use bottom width / wetting perimeter to estimate.
        }
        m_bankSto[i] += rttlc * (1.f - trnsrch); // m^3
        if (nullptr != m_gwSto) {
            m_gwSto[i] += rttlc * trnsrch / m_chArea[i] * 1000.f; // mm
        }
    }

    // todo, compute revap from bank storage. In SWAT, revap coefficient is equal to gw_revap.

#ifdef PRINT_DEBUG
    cout << " chStorage after routing " << m_chStorage[i] << endl;
    cout << " surfq: " << m_qsCh[i] << ", ifluq: " << m_qiCh[i] << ", groudq: " << m_qgCh[i] << endl;
#endif
    return true;
}
//ljj++
bool MUSK_CH::LakeBudget(const int i) {
    m_chWtrDepth[i] = 0.f;
    //! 1. add all the inflow water
    float qIn = 0.f; /// Water entering reach on current day from both current subbasin and upstreams
    // 1.1. water from this subbasin
    qIn += m_olQ2Rch[i]; /// surface flow
    float qiSub = 0.f;   /// interflow flow
    if (nullptr != m_ifluQ2Rch && m_ifluQ2Rch[i] >= 0.f) {
        qiSub = m_ifluQ2Rch[i];
        qIn += qiSub;
    }
    float qgSub = 0.f; /// groundwater flow
    if (nullptr != m_gndQ2Rch && m_gndQ2Rch[i] >= 0.f) {
        qgSub = m_gndQ2Rch[i];
        qIn += qgSub;
    }
    float qgsSub = 0.f;
    float dh = gw_height[i];   //lake 
    float rto = m_gw_sh[i]/gw_height[i];
    if (abs(dh) <= 1e-6) dh = 0.f;
    float Q = m_Kg * pow(abs(dh), m_Base_ex) * rto; // m3
    Q = Max(Q,1.e-6f);
    qgsSub += Q/m_dt;
    qIn += qgsSub;
    m_gw_sh[i] -= Q; // m3

    // 1.2. water from upstream reaches
    float qsUp = 0.f;
    float qiUp = 0.f;
    float qgUp = 0.f;
    float qgsUp = 0.f;
    for (auto upRchID = m_reachUpStream.at(i).begin(); upRchID != m_reachUpStream.at(i).end(); ++upRchID) {
        if (m_qsRchOut[*upRchID] != m_qsRchOut[*upRchID]) {
            cout << "DayOfYear: " << m_dayOfYear << ", rchID: " << i << ", upRchID: " << *upRchID <<
                    ", surface part illegal!" << endl;
            return false;
        }
        if (m_qiRchOut[*upRchID] != m_qiRchOut[*upRchID]) {
            cout << "DayOfYear: " << m_dayOfYear << ", rchID: " << i << ", upRchID: " << *upRchID <<
                    ", subsurface part illegal!" << endl;
            return false;
        }
        if (m_qgRchOut[*upRchID] != m_qgRchOut[*upRchID]) {
            cout << "DayOfYear: " << m_dayOfYear << ", rchID: " << i << ", upRchID: " << *upRchID <<
                    ", groundwater part illegal!" << endl;
            return false;
        }
        if (m_qgRchOut[*upRchID] != m_qgRchOut[*upRchID]) {
            cout << "DayOfYear: " << m_dayOfYear << ", rchID: " << i << ", upRchID: " << *upRchID <<
                    ", groundwater part illegal!" << endl;
            return false;
        }
        if (m_qsRchOut[*upRchID] > 0.f) qsUp += m_qsRchOut[*upRchID];
        if (m_qiRchOut[*upRchID] > 0.f) qiUp += m_qiRchOut[*upRchID];
        if (m_qgRchOut[*upRchID] > 0.f) qgUp += m_qgRchOut[*upRchID];
        if (m_qgsRchOut[*upRchID] > 0.f) qgsUp += m_qgsRchOut[*upRchID];
    }
    qIn += qsUp + qiUp + qgUp + qgsUp;   

    float pre_Sto = m_chSto[i];
    //add precipitation
    m_chSto[i] += m_prec[i]* m_dt;  
    float wtrin = qIn * m_dt;  //m3  
    float rtwtr = 0.f;    //flow out

    //!!! important: for lake unit, chSto == lakeSto
    //m_chSto[i] += wtrin;
    
    //! 2. minus all the outflow water
    // linear reservior
    rtwtr = 0.f;  //m3/s
    float h0 = 0.f;
    h0 = m_minvol;
	float thwl = m_lakedpini[i] * h0; 
    float max_outflow = Max(0.f, (m_lakedp[i] - thwl) * m_lakearea[i]);
	// lake evaporation
	float rtevp = 0.f; 
	rtevp = m_evlake * m_petSubbsn[i]/m_petFactor * 0.001f * m_lakearea[i]; //m3
	rtevp = Min(rtevp, m_chSto[i]);
    m_chSto[i] -= rtevp;
	// lake groundwater
	float LakeOutGw = m_chSto[i] * m_lakeseep;
    LakeOutGw = Min(LakeOutGw, m_chSto[i]);
    if (nullptr != m_gwSto) {
        m_gwSto[i] += LakeOutGw / m_chArea[i] * 1000.f; // updated groundwater storage
        m_chSto[i] -= LakeOutGw;
    }
	
	// add qIn
	m_chSto[i] += qIn * m_dt;

    float SI = m_chSto[i] / m_dt + (qIn + m_qin1[i])/2 - m_qout1[i]/2;

    // Lake parameter A (suggested  value equal to outflow width in [m])
    // float lakefactor = m_lakearea[i] / m_dt / sqrt(alpha[i]);
    float lakefactor = m_lakearea[i] / m_dt / sqrt(m_chWth[i]*m_lakealpha[i]);
    if (m_lakedp[i] > thwl && m_chSto[i]>0.f){
        rtwtr = pow(sqrt((lakefactor*lakefactor) + 2*SI)-lakefactor,2);
        if ((m_qout1[i]+rtwtr)*0.5 * m_dt > max_outflow){
            rtwtr = 2*(max_outflow / m_dt) -m_qout1[i];
            rtwtr = Max(rtwtr,0.f);
        }
        m_chSto[i] -= (m_qout1[i]+rtwtr)*0.5 * m_dt;
    }
    else{
        rtwtr = 0.f;
    }

	// update lake water level
    if(m_chSto[i]>=m_lakevol[i]) {
        rtwtr +=  (m_chSto[i]- m_lakevol[i])/m_dt;
        //m_qRchOut[i]+= (m_chSto[i]- m_lakevol[i])/m_dt;  
        m_chSto[i]  = m_lakevol[i];
    }
    m_qRchOut[i] = rtwtr*0.5;
    m_lakedp[i] = m_chSto[i] / m_lakearea[i];
    m_rteWtrOut[i] = m_qRchOut[i] * m_dt;   // m^3

    m_qin1[i] = qIn;
    m_qout1[i] = rtwtr*0.5;
    m_chWtrDepth[i] = m_lakedp[i];

    float qInSum = m_olQ2Rch[i] + qiSub + qgSub + qsUp + qiUp + qgUp;
    if (qInSum < UTIL_ZERO) {
        // In case of divided by zero.
        // m_qsRchOut[i] = 0.f;
        // m_qiRchOut[i] = 0.f;
        // m_qgRchOut[i] = 0.f;
        // m_qgsRchOut[i] = 0.f;
        // m_qRchOut[i] = 0.f;
    } else {
        // In my opinion, these lines should use `qIn` instead of `qInSum`. By lj.
        // m_qsRchOut[i] = m_qRchOut[i] * (m_olQ2Rch[i] + qsUp) / qIn;
        // m_qiRchOut[i] = m_qRchOut[i] * (qiSub + qiUp) / qIn;
        // m_qgRchOut[i] = m_qRchOut[i] * (qgSub + qgUp) / qIn;
        m_qsRchOut[i] = m_qRchOut[i];
        m_qiRchOut[i] = 0.f;
        m_qgRchOut[i] = 0.f;
        m_qgsRchOut[i] = 0.f;
    }
    return true;
}

bool MUSK_CH::ResBudget(const int i) {
    //! 1. add all the inflow water
    float qIn = 0.f; /// Water entering reach on current day from both current subbasin and upstreams
    // 1.1. water from this subbasin
    qIn += m_olQ2Rch[i]; /// surface flow
    float qiSub = 0.f;   /// interflow flow
    if (nullptr != m_ifluQ2Rch && m_ifluQ2Rch[i] >= 0.f) {
        qiSub = m_ifluQ2Rch[i];
        qIn += qiSub;
    }
    float qgSub = 0.f; /// groundwater flow
    if (nullptr != m_gndQ2Rch && m_gndQ2Rch[i] >= 0.f) {
        qgSub = m_gndQ2Rch[i];
        qIn += qgSub;
    }
    float qgsSub = 0.f;
    float dh = gw_height[i];   //lake 
    float rto = m_gw_sh[i]/gw_height[i];
    if (abs(dh) <= 1e-6) dh = 0.f;
    float Q = m_Kg * pow(abs(dh), m_Base_ex) * rto; // m3
    Q = Max(Q,1.e-6f);
    qgsSub += Q/m_dt;
    qIn += qgsSub;
    m_gw_sh[i] -= Q; // m3
    // 1.2. water from upstream reaches
    float qsUp = 0.f;
    float qiUp = 0.f;
    float qgUp = 0.f;
    float qgsUp = 0.f;
    for (auto upRchID = m_reachUpStream.at(i).begin(); upRchID != m_reachUpStream.at(i).end(); ++upRchID) {
        if (m_qsRchOut[*upRchID] != m_qsRchOut[*upRchID]) {
            cout << "DayOfYear: " << m_dayOfYear << ", rchID: " << i << ", upRchID: " << *upRchID <<
                    ", surface part illegal!" << endl;
            return false;
        }
        if (m_qiRchOut[*upRchID] != m_qiRchOut[*upRchID]) {
            cout << "DayOfYear: " << m_dayOfYear << ", rchID: " << i << ", upRchID: " << *upRchID <<
                    ", subsurface part illegal!" << endl;
            return false;
        }
        if (m_qgRchOut[*upRchID] != m_qgRchOut[*upRchID]) {
            cout << "DayOfYear: " << m_dayOfYear << ", rchID: " << i << ", upRchID: " << *upRchID <<
                    ", groundwater part illegal!" << endl;
            return false;
        }
        if (m_qsRchOut[*upRchID] > 0.f) qsUp += m_qsRchOut[*upRchID];
        if (m_qiRchOut[*upRchID] > 0.f) qiUp += m_qiRchOut[*upRchID];
        if (m_qgRchOut[*upRchID] > 0.f) qgUp += m_qgRchOut[*upRchID];
        if (m_qgsRchOut[*upRchID] > 0.f) qgsUp += m_qgsRchOut[*upRchID];
    }
    qIn += qsUp + qiUp + qgUp + qgsUp;     

    float pre_Sto = m_chSto[i];

    //add precipitation
    m_chSto[i] += m_prec[i]*m_dt;  

    float wtrin = qIn * m_dt;  //m3  
    float rtwtr = 0.f;    //flow out

	// evaporation
	float rtevp = 0.f; 
	rtevp = m_evlake * m_petSubbsn[i]/m_petFactor * 0.001f * m_lakearea[i]; //m3
	rtevp = Min(rtevp, m_chSto[i]);
    m_chSto[i] -= rtevp;

	// add qIn
	m_chSto[i] += qIn * m_dt;

    //initial 
    float TotalReservoirStorageM3CC = m_lakevol[i];
    float ConservativeStorageLimitCC = m_ResLc[i]; //minimum storage, 0.1 as defult
    float NormalStorageLimitCC = m_ResLn[i]; //normal storage, 0.3 as defult
    float FloodStorageLimitCC = m_ResLf[i]; //maximum storage, 0.97 as defult
    float adjust_Normal_FloodCC = m_ResAdjust[i]; //adjust parameter, 0.01~0.99, 0.7 as defult
    float Normal_FloodStorageLimitCC = NormalStorageLimitCC + adjust_Normal_FloodCC * (FloodStorageLimitCC - NormalStorageLimitCC);

    float InvDtSecDay = m_dt;
    float MinReservoirOutflowCC = 0.05* m_natural_flow[i];  //minimum outflow 
    float NormalReservoirOutflowCC = 0.3* m_natural_flow[i]; //normal outflow 
    float NonDamagingReservoirOutflowCC = 0.97* m_natural_flow[i]; //non-damaging outflow 
    float ReservoirRnormqMultCC=1.f;  //!!todo
    NormalReservoirOutflowCC = NormalReservoirOutflowCC * ReservoirRnormqMultCC;
    if(NormalReservoirOutflowCC > MinReservoirOutflowCC){
        NormalReservoirOutflowCC = NormalReservoirOutflowCC;
    }else{
        NormalReservoirOutflowCC = MinReservoirOutflowCC+0.01;
    }
    if(NormalReservoirOutflowCC < NonDamagingReservoirOutflowCC){
        NormalReservoirOutflowCC = NormalReservoirOutflowCC;
    }else{
        NormalReservoirOutflowCC = NonDamagingReservoirOutflowCC-0.01;
    }

    //New reservoir fill (fraction)
    float ReservoirFillCC = m_chSto[i] / TotalReservoirStorageM3CC;
    
    //below 2Lc
    float ReservoirOutflow1 = Min(MinReservoirOutflowCC, m_chSto[i] / InvDtSecDay);

    //2Lc<F<=Ln
    float DeltaO = NormalReservoirOutflowCC - MinReservoirOutflowCC;
    float DeltaLN = NormalStorageLimitCC - 2 * ConservativeStorageLimitCC;
    float ReservoirOutflow2 = MinReservoirOutflowCC + DeltaO * (ReservoirFillCC - 2 * ConservativeStorageLimitCC) / DeltaLN;
    
    //Ln<F<Lf
    float DeltaNFL = FloodStorageLimitCC - Normal_FloodStorageLimitCC;
    float ReservoirOutflow3a = NormalReservoirOutflowCC;
    float ReservoirOutflow3b = NormalReservoirOutflowCC + ((ReservoirFillCC - Normal_FloodStorageLimitCC) / DeltaNFL) 
                                * (NonDamagingReservoirOutflowCC - NormalReservoirOutflowCC);

    //F>Lf
    float temp = Min(NonDamagingReservoirOutflowCC, Max(qIn * 1.2, NormalReservoirOutflowCC));
    float ReservoirOutflow4 = Min(Max((ReservoirFillCC - FloodStorageLimitCC-0.01) *
                                TotalReservoirStorageM3CC / InvDtSecDay,NonDamagingReservoirOutflowCC), temp);
    
    //Reservoir outflow [m3/s] 
    rtwtr = 0.f;
    rtwtr = ReservoirOutflow1;
    if(ReservoirFillCC > 2 * ConservativeStorageLimitCC) rtwtr = ReservoirOutflow2;
    if(ReservoirFillCC > NormalStorageLimitCC) rtwtr = ReservoirOutflow3a;
    if(ReservoirFillCC > Normal_FloodStorageLimitCC) rtwtr = ReservoirOutflow3b;
    if(ReservoirFillCC > FloodStorageLimitCC) rtwtr = ReservoirOutflow4;


    temp = Min(rtwtr,Max(qIn, NormalReservoirOutflowCC));

    if((rtwtr > 1.2 * qIn) & (rtwtr > NormalReservoirOutflowCC) & (ReservoirFillCC < FloodStorageLimitCC)){
        rtwtr = temp;
    }
    m_qRchOut[i] = rtwtr;
    m_rteWtrOut[i] = m_qRchOut[i] * m_dt;   // m^3
    m_rteWtrOut[i] = Min(m_rteWtrOut[i], m_chSto[i]);
    m_rteWtrOut[i] = Max(m_rteWtrOut[i], m_chSto[i] - TotalReservoirStorageM3CC);
    m_qRchOut[i] = m_rteWtrOut[i] / m_dt;   //m3 s-1
    m_chSto[i] = m_chSto[i] - m_rteWtrOut[i];
    ReservoirFillCC = m_chSto[i] /TotalReservoirStorageM3CC;
    if(ReservoirFillCC<=0.f) ReservoirFillCC = 0.f;
    m_chWtrDepth[i] = ReservoirFillCC;

    float qInSum = m_olQ2Rch[i] + qiSub + qgSub + qsUp + qiUp + qgUp;
    if (qInSum < UTIL_ZERO) {
        // In case of divided by zero.
        // m_qsRchOut[i] = 0.f;
        // m_qiRchOut[i] = 0.f;
        // m_qgRchOut[i] = 0.f;
        // m_qRchOut[i] = 0.f;
    } else {
        // In my opinion, these lines should use `qIn` instead of `qInSum`. By lj.
        // m_qsRchOut[i] = m_qRchOut[i] * (m_olQ2Rch[i] + qsUp) / qIn;
        // m_qiRchOut[i] = m_qRchOut[i] * (qiSub + qiUp) / qIn;
        // m_qgRchOut[i] = m_qRchOut[i] * (qgSub + qgUp) / qIn;
        m_qsRchOut[i] = m_qRchOut[i];
        //if(m_temp[i] <0) m_qsRchOut[i] +=m_prec[i] +qsUp+ qiUp + qgUp+ m_olQ2Rch[i]+m_ifluQ2Rch[i]+m_gndQ2Rch[i];
        m_qiRchOut[i] = 0.f;
        m_qgRchOut[i] = 0.f;
    }

    return true;
}

void MUSK_CH::SetSubbasins(clsSubbasins* subbsns) {
    if (m_subbasinsInfo == nullptr) {
        m_subbasinsInfo = subbsns;
        // m_nSubbasins = m_subbasinsInfo->GetSubbasinNumber(); // Set in SetValue()! lj
        m_subbasinIDs = m_subbasinsInfo->GetSubbasinIDs();
    }
}