#include "DOCGroundwater.h"

#include "text.h"

DOCGroundwater::DOCGroundwater():
	//input
	m_nSubbsns(-1), m_recessionCoefficient(-1.f),  m_nCells(-1),curBasinArea(nullptr),
	m_subbsnID(nullptr),m_soilPerco(nullptr),  gw_delaye(nullptr),m_hlife_docgw(NODATA_VALUE),
	m_subbasinsInfo(nullptr),m_delay(NODATA_VALUE),m_dp_co(NODATA_VALUE),m_area(nullptr),
	m_gw_shallow(nullptr),
	
	//output
	m_gwS_DOCSto(nullptr), m_gwS_DOCconc(nullptr), m_gwD_DOCSto(nullptr), m_Deepgrndwtr_DOC(nullptr),m_recharge1(nullptr),
	m_DeepDOCtoCH(nullptr)
{
}

DOCGroundwater::~DOCGroundwater() {
	if (curBasinArea != nullptr) Release1DArray(curBasinArea);
	if (gw_delaye != nullptr) Release1DArray(gw_delaye);
	if (m_recharge1 != nullptr) Release1DArray(m_recharge1);
	if (m_gwS_DOCSto != nullptr) Release1DArray(m_gwS_DOCSto);
	if (m_gwD_DOCSto != nullptr) Release1DArray(m_gwD_DOCSto);
	if (m_Deepgrndwtr_DOC != nullptr) Release1DArray(m_Deepgrndwtr_DOC);
	if (m_DeepDOCtoCH != nullptr) Release1DArray(m_DeepDOCtoCH);
	if (m_gwS_DOCconc != nullptr) Release1DArray(m_gwS_DOCconc);

}

void DOCGroundwater::SetSubbasins(clsSubbasins* subbasins) {
	if (nullptr == m_subbasinsInfo) {
		m_subbasinsInfo = subbasins;
		m_subbasinIDs = m_subbasinsInfo->GetSubbasinIDs();
	}
}

bool DOCGroundwater::CheckInputData() {
	CHECK_POSITIVE(MID_CarbonGW, m_nSubbsns);
	CHECK_POSITIVE(MID_CarbonGW, m_nCells);
    return true;
}

void DOCGroundwater::SetValue(const char* key, const float value) {
    string sk(key);
	if (StringMatch(sk, VAR_SUBBSNID_NUM)) m_nSubbsns = CVT_INT(value);
	if (StringMatch(sk, Tag_CellSize)) m_nCells = CVT_INT(value);
	if (StringMatch(sk, VAR_DELAY)) m_delay = value;
	if (StringMatch(sk, VAR_DF_COEF)) m_dp_co = value;
	if (StringMatch(sk, VAR_HLDOCGW)) m_hlife_docgw = value;
		
	// if (StringMatch(sk, VAR_KG)) m_recessionCoefficient = value;
	
	if (StringMatch(sk, VAR_KDOC)) m_kdoc = value;
}

void DOCGroundwater::Set1DData(const char* key, const int n, float* data) {
	string sk(key);  
	if (StringMatch(sk, VAR_SUBBSN)) m_subbsnID = data;
	if (StringMatch(sk, VAR_AHRU)) m_area = data;
	if (StringMatch(sk, VAR_PERC_LOWEST_DOC)) m_soilPerco = data;
	if (StringMatch(sk, VAR_GW_SH)) {
        CheckInputSize(MID_MUSK_CH, key, n - 1, m_nSubbsns);
        m_gw_shallow = data;
    }

}

void DOCGroundwater::Set2DData(const char* key, const int nrows, const int ncols, float** data) {
	string sk(key);

}

void DOCGroundwater::InitialOutputs() {
	CHECK_POSITIVE(MID_CarbonGW, m_nSubbsns);
	CHECK_POSITIVE(MID_CarbonGW, m_nCells);
	if (curBasinArea == nullptr) {
		Initialize1DArray(m_nSubbsns + 1, curBasinArea, 0.f);
		for (auto it = m_subbasinIDs.begin(); it != m_subbasinIDs.end(); ++it) {
		int subID = *it;
		Subbasin* curSub = m_subbasinsInfo->GetSubbasinByID(subID);
		// get percolation from the bottom soil layer at the subbasin scale
		int curCellsNum = curSub->GetCellCount();
		int* curCells = curSub->GetCells();
			for (int i = 0; i < curCellsNum; i++) {
				curBasinArea[subID] += m_area[curCells[i]];
			}
		}
	}
	if (gw_delaye == nullptr) Initialize1DArray(m_nSubbsns + 1, gw_delaye, exp(-1./m_delay));
	if (m_recharge1 == nullptr) Initialize1DArray(m_nSubbsns + 1, m_recharge1, 0.f);
	if (m_gwD_DOCSto == nullptr)Initialize1DArray(m_nSubbsns + 1, m_gwD_DOCSto, 0.f);
	if (m_Deepgrndwtr_DOC == nullptr)Initialize1DArray(m_nSubbsns + 1, m_Deepgrndwtr_DOC, 0.f);
	if (m_DeepDOCtoCH == nullptr)Initialize1DArray(m_nSubbsns + 1, m_DeepDOCtoCH, 0.f);
	if (m_gwS_DOCSto == nullptr)Initialize1DArray(m_nSubbsns + 1, m_gwS_DOCSto, 0.f);
	if (m_gwS_DOCconc == nullptr)Initialize1DArray(m_nSubbsns + 1, m_gwS_DOCconc, 0.f);

}

int DOCGroundwater::Execute() {
    CheckInputData();
    InitialOutputs();
	for (auto it = m_subbasinIDs.begin(); it != m_subbasinIDs.end(); ++it) {
		int subID = *it;
		m_DeepDOCtoCH[subID] = 0.f;
		m_gwS_DOCconc[subID] = 0.f;

		Subbasin* curSub = m_subbasinsInfo->GetSubbasinByID(subID);
		// get percolation from the bottom soil layer at the subbasin scale
		int curCellsNum = curSub->GetCellCount();
		int* curCells = curSub->GetCells();
		float perco = 0.f;
		float rchrg1 = 0.f;
		if (curSub->GetPerco() > 0.f) rchrg1 = m_recharge1[subID] ;
		for (int i = 0; i < curCellsNum; i++) {
            int index = curCells[i];
            float tmp_perc = m_soilPerco[index];
            if (tmp_perc > 0) {
                perco += tmp_perc * (m_area[index] / curBasinArea[subID])*(1-gw_delaye[subID])+rchrg1* (m_area[index] / curBasinArea[subID])*gw_delaye[subID];// + rchrg1*gw_delaye[subID];
			} else {
                m_soilPerco[index] = 0.f;
            }
        }
		m_recharge1[subID] = perco;
		
		float ratio2gw = 1.f;
        perco *= ratio2gw;
        float percoDeep = perco * m_dp_co; ///< deep percolation

        m_gwS_DOCSto[subID] += (perco - percoDeep);   //kg/ha
		m_gwS_DOCSto[subID] = m_gwS_DOCSto[subID] * exp(-0.693f /m_hlife_docgw); 
		m_gwS_DOCconc[subID] = 100* m_gwS_DOCSto[subID] / (m_gw_shallow[subID]/curBasinArea[subID]*1000.f);  //mg/L

		m_gwD_DOCSto[subID] += percoDeep;
		//if(subID == 1) cout<<m_dp_co<<"  "<<percoDeep<<"  "<<m_gwD_DOCSto[subID]<<"  "<<endl;
		float m_deepWaterDepth = curSub->GetGw();
		float m_RG = curSub->GetRg();
		float xx = (m_deepWaterDepth + m_RG);
		if (xx > 0.f) {
			xx = m_gwD_DOCSto[subID] / (m_deepWaterDepth + m_RG);
		}
		else {
			xx = 0.f;
		}

		if (xx < 1.e-6f) xx = 0.f;
		m_Deepgrndwtr_DOC[subID] = 0.f;
		m_Deepgrndwtr_DOC[subID] = xx * m_RG;
		//subtract DOC transport losses from the shallow aquifer
		m_gwD_DOCSto[subID] = m_gwD_DOCSto[subID] - m_Deepgrndwtr_DOC[subID];
		m_gwD_DOCSto[subID] = Max(0.f, m_gwD_DOCSto[subID]);

		//compute DOC reaction losses in the groundwater
		m_gwD_DOCSto[subID] = m_gwD_DOCSto[subID] * exp(-0.693f /m_hlife_docgw); 
		m_gwD_DOCSto[subID] = Max(0.f, m_gwD_DOCSto[subID]);

		m_gwS_DOCSto[subID] = Max(0.f, m_gwS_DOCSto[subID]);
	}

#pragma omp parallel
    {
        float* tmp_dDOCtoCH = new(nothrow) float[m_nSubbsns + 1];
        for (int i = 0; i <= m_nSubbsns; i++) {
            tmp_dDOCtoCH[i] = 0.f;

        }
#pragma omp for
        for (int i = 1; i < m_nSubbsns; i++) {
            m_DeepDOCtoCH[i] += m_Deepgrndwtr_DOC[i] * curBasinArea[i] * 0.0001f;  //kg
        }

		delete[] tmp_dDOCtoCH;
        tmp_dDOCtoCH = nullptr;

    } /* END of #pragma omp parallel */
    // sum all the subbasins and put the sum value in the zero-index of the array
    for (int i = 1; i < m_nSubbsns + 1; i++) {
        m_DeepDOCtoCH[0] += m_DeepDOCtoCH[i];   //units: kg

    }
		return 0;
}

void DOCGroundwater::Get1DData(const char* key, int* n, float** data) {
    InitialOutputs();
    string sk(key);
	if   (StringMatch(sk, VAR_GWD_RDOCtoCH)) {
		*data = m_DeepDOCtoCH;
		*n = m_nSubbsns + 1;
	}
	else if (StringMatch(sk, VAR_GWS_RDOCsto)) {
		*data = m_gwS_DOCSto;
		*n = m_nSubbsns + 1;
	}
	else if (StringMatch(sk, VAR_GWS_RDOCconc)) {
		*data = m_gwS_DOCconc;
		*n = m_nSubbsns + 1;
	}
   }


void DOCGroundwater::Get2DData(const char* key, int* nrows, int* ncols, float*** data) {
	InitialOutputs();
	string sk(key);
}