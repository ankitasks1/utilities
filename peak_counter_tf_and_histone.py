import os,sys
path = os.getcwd()

#print(path)


list1_data = ('''
CRAMP1_NTC_IgG_NTC_peaks.narrowPeak
CRAMP1_shSUZ12_IgG_shSUZ12_peaks.narrowPeak
FLAG-H14_FLAG_control_peaks.narrowPeak
H14_NTC_IgG_NTC_peaks.narrowPeak
H14_shCRAMP1_IgG_shCRAMP1_peaks.narrowPeak
H14_shSUZ12_IgG_shSUZ12_peaks.narrowPeak
H1_NTC_IgG_NTC_peaks.narrowPeak
H1_shCRAMP1_IgG_shCRAMP1_peaks.narrowPeak
H1_shSUZ12_IgG_shSUZ12_peaks.narrowPeak
H3K27me3_NTC_IgG_NTC_peaks.narrowPeak
H3K27me3_shCRAMP1_IgG_shCRAMP1_peaks.narrowPeak
H3K27me3_shSUZ12_IgG_shSUZ12_peaks.narrowPeak
H3K9me3_NTC_IgG_NTC_peaks.narrowPeak
H3K9me3_shCRAMP1_IgG_shCRAMP1_peaks.narrowPeak
H3K9me3_shSUZ12_IgG_shSUZ12_peaks.narrowPeak
MTF2_NTC_IgG_NTC_peaks.narrowPeak
MTF2_shCRAMP1_IgG_shCRAMP1_peaks.narrowPeak
MTF2_shSUZ12_IgG_shSUZ12_peaks.narrowPeak
S12_shCRAMP1_IgG_shCRAMP1_peaks.narrowPeak
SUZ12_NTC_IgG_NTC_peaks.narrowPeak
V5-CRAMP1_V5_control_peaks.narrowPeak
ADNP_ENCFF739AJO_peaks_id.bed
AFF1_ENCFF195YGC_peaks_id.bed
AFF1_ENCFF674XTY_peaks_id.bed
AGO1_ENCFF100VYA_peaks_id.bed
AGO1_ENCFF794IRP_peaks_id.bed
AGO1_ENCFF836VRV_peaks_id.bed
ARHGAP35_ENCFF952WKN_peaks_id.bed
ARID1B_ENCFF879NTL_peaks_id.bed
ARID2_ENCFF913WRW_peaks_id.bed
ARID3A_ENCFF891OQP_peaks_id.bed
ARID3B_ENCFF270TSN_peaks_id.bed
ARID4B_ENCFF086FAZ_peaks_id.bed
ARID4B_ENCFF528IDR_peaks_id.bed
ARNT_ENCFF075OIT_peaks_id.bed
ARNT_ENCFF553BVI_peaks_id.bed
ARNT_ENCFF728ITJ_peaks_id.bed
ASH1L_ENCFF553IUR_peaks_id.bed
ATF1_ENCFF206BGR_peaks_id.bed
ATF1_ENCFF627RSK_peaks_id.bed
ATF1_ENCFF793HSJ_peaks_id.bed
ATF1_ENCFF924GPK_peaks_id.bed
ATF2_ENCFF121HYT_peaks_id.bed
ATF2_ENCFF192ASP_peaks_id.bed
ATF3_ENCFF107DBQ_peaks_id.bed
ATF3_ENCFF556SMM_peaks_id.bed
ATF3_ENCFF718VHT_peaks_id.bed
ATF3_ENCFF895RNS_peaks_id.bed
ATF4_ENCFF086TAD_peaks_id.bed
ATF4_ENCFF250MUC_peaks_id.bed
ATF6_ENCFF863ZFH_peaks_id.bed
ATF7_ENCFF951BFN_peaks_id.bed
BACH1_ENCFF026BDA_peaks_id.bed
BACH1_ENCFF678TXX_peaks_id.bed
BCL6_ENCFF941EDY_peaks_id.bed
BCLAF1_ENCFF070ZTX_peaks_id.bed
BCLAF1_ENCFF725JXZ_peaks_id.bed
BCLAF1_ENCFF941SRW_peaks_id.bed
BCOR_ENCFF124FTW_peaks_id.bed
BDP1_ENCFF284FTY_peaks_id.bed
BHLHE40_ENCFF154IVU_peaks_id.bed
BMI1_ENCFF600KEF_peaks_id.bed
BRCA1_ENCFF258GDS_peaks_id.bed
BRCA1_ENCFF620FIH_peaks_id.bed
BRD4_ENCFF130JVF_peaks_id.bed
BRD9_ENCFF775WSL_peaks_id.bed
BRF1_ENCFF577LSK_peaks_id.bed
BRF2_ENCFF944LWS_peaks_id.bed
C11orf30_ENCFF038WQY_peaks_id.bed
CAMTA2_ENCFF123IHX_peaks_id.bed
CBFA2T2_ENCFF543GET_peaks_id.bed
CBFA2T3_ENCFF082DOH_peaks_id.bed
CBFB_ENCFF802NHC_peaks_id.bed
CBX1_ENCFF470CSE_peaks_id.bed
CBX2_ENCFF258XBJ_peaks_id.bed
CBX3_ENCFF068OEJ_peaks_id.bed
CBX3_ENCFF386ZWO_peaks_id.bed
CBX5_ENCFF066OJB_peaks_id.bed
CBX8_ENCFF522HZT_peaks_id.bed
CC2D1A_ENCFF156GHD_peaks_id.bed
CCAR2_ENCFF583GKE_peaks_id.bed
CCAR2_ENCFF657LFS_peaks_id.bed
CCAR2_ENCFF704PGT_peaks_id.bed
CCNT2_ENCFF642ZYO_peaks_id.bed
CDC5L_ENCFF408RSJ_peaks_id.bed
CEBPB_ENCFF022KBK_peaks_id.bed
CEBPB_ENCFF712ZNR_peaks_id.bed
CEBPB_ENCFF882ARK_peaks_id.bed
CEBPG_ENCFF456USL_peaks_id.bed
CEBPG_ENCFF558AJI_peaks_id.bed
CEBPG_ENCFF998LIR_peaks_id.bed
CEBPZ_ENCFF588GNU_peaks_id.bed
CGGBP1_ENCFF771XZZ_peaks_id.bed
CHAMP1_ENCFF515GUE_peaks_id.bed
CHAMP1_ENCFF541VYN_peaks_id.bed
CHCHD3_ENCFF645IVF_peaks_id.bed
CHD1_ENCFF408NUX_peaks_id.bed
CHD2_ENCFF947AEO_peaks_id.bed
CHD4_ENCFF985QBS_peaks_id.bed
CHD7_ENCFF722UJW_peaks_id.bed
CLOCK_ENCFF648UJW_peaks_id.bed
COPS2_ENCFF942JZR_peaks_id.bed
CREB1_ENCFF193LLN_peaks_id.bed
CREB1_ENCFF970QKS_peaks_id.bed
CREB3_ENCFF950IYB_peaks_id.bed
CREB3L1_ENCFF118GMS_peaks_id.bed
CREB5_ENCFF875JMR_peaks_id.bed
CREBBP_ENCFF532VPN_peaks_id.bed
CREM_ENCFF324ELP_peaks_id.bed
CSDE1_ENCFF079ERC_peaks_id.bed
CSDE1_ENCFF746OAD_peaks_id.bed
CTBP1_ENCFF782PIZ_peaks_id.bed
CTCF_ENCFF221SKA_peaks_id.bed
CTCF_ENCFF582SNT_peaks_id.bed
CTCF_ENCFF660GHM_peaks_id.bed
CTCF_ENCFF736NYC_peaks_id.bed
CTCF_ENCFF769AUF_peaks_id.bed
CTCFL_ENCFF630YVJ_peaks_id.bed
CUX1_ENCFF136KLM_peaks_id.bed
CUX1_ENCFF838GFC_peaks_id.bed
CXXC5_ENCFF798LLI_peaks_id.bed
DACH1_ENCFF729DNM_peaks_id.bed
DDIT3_ENCFF869RFC_peaks_id.bed
DDX20_ENCFF746TUQ_peaks_id.bed
DDX20_ENCFF748EAO_peaks_id.bed
DEAF1_ENCFF030AXK_peaks_id.bed
DEAF1_ENCFF813TRY_peaks_id.bed
DIDO1_ENCFF928BHE_peaks_id.bed
DLX4_ENCFF832DBU_peaks_id.bed
DMBX1_ENCFF359ZQQ_peaks_id.bed
DMTF1_ENCFF505INE_peaks_id.bed
DNMT1_ENCFF412WOK_peaks_id.bed
DPF2_ENCFF542HJS_peaks_id.bed
DPF2_ENCFF716PXH_peaks_id.bed
E2F1_ENCFF053HPQ_peaks_id.bed
E2F1_ENCFF104VFF_peaks_id.bed
E2F3_ENCFF710UAZ_peaks_id.bed
E2F4_ENCFF064CWA_peaks_id.bed
E2F4_ENCFF221TRR_peaks_id.bed
E2F5_ENCFF082QHP_peaks_id.bed
E2F5_ENCFF173QUY_peaks_id.bed
E2F6_ENCFF085GHX_peaks_id.bed
E2F6_ENCFF831NDQ_peaks_id.bed
E2F7_ENCFF938ZPZ_peaks_id.bed
E2F8_ENCFF925GLW_peaks_id.bed
E4F1_ENCFF439KTZ_peaks_id.bed
E4F1_ENCFF683WRK_peaks_id.bed
EGR1_ENCFF136DDO_peaks_id.bed
EGR1_ENCFF566PRZ_peaks_id.bed
EGR1_ENCFF640AKG_peaks_id.bed
EHMT2_ENCFF189OHQ_peaks_id.bed
ELF1_ENCFF133TSU_peaks_id.bed
ELF1_ENCFF350OMH_peaks_id.bed
ELF1_ENCFF688TNZ_peaks_id.bed
ELF1_ENCFF763MXW_peaks_id.bed
ELF2_ENCFF695PDY_peaks_id.bed
ELF4_ENCFF518EGY_peaks_id.bed
ELF4_ENCFF616FCV_peaks_id.bed
ELF4_ENCFF832HNJ_peaks_id.bed
ELK1_ENCFF586MZX_peaks_id.bed
ELK1_ENCFF715WGN_peaks_id.bed
ELK3_ENCFF198NGY_peaks_id.bed
EP300_ENCFF702XPO_peaks_id.bed
EP300_ENCFF821MKR_peaks_id.bed
EP400_ENCFF167CQF_peaks_id.bed
ERF_ENCFF330EGV_peaks_id.bed
ERF_ENCFF462ZIG_peaks_id.bed
ESRRA_ENCFF588QRH_peaks_id.bed
ESRRB_ENCFF106ECV_peaks_id.bed
ETS1_ENCFF886BDQ_peaks_id.bed
ETS2_ENCFF772QLT_peaks_id.bed
ETV1_ENCFF065RZP_peaks_id.bed
ETV5_ENCFF604LXR_peaks_id.bed
ETV6_ENCFF405QTW_peaks_id.bed
ETV6_ENCFF584QFY_peaks_id.bed
ETV6_ENCFF965PER_peaks_id.bed
EWSR1_ENCFF924FYI_peaks_id.bed
EZH2_ENCFF804RVA_peaks_id.bed
FIP1L1_ENCFF084DTV_peaks_id.bed
FIP1L1_ENCFF285TMA_peaks_id.bed
FIP1L1_ENCFF883WYX_peaks_id.bed
FOS_ENCFF258PLH_peaks_id.bed
FOSL1_ENCFF004HXL_peaks_id.bed
FOSL1_ENCFF256NXT_peaks_id.bed
FOXA1_ENCFF497OQD_peaks_id.bed
FOXA3_ENCFF197SXI_peaks_id.bed
FOXA3_ENCFF504WWL_peaks_id.bed
FOXJ2_ENCFF957CVJ_peaks_id.bed
FOXJ3_ENCFF567CPM_peaks_id.bed
FOXK1_ENCFF516ZWP_peaks_id.bed
FOXK2_ENCFF702AQS_peaks_id.bed
FOXK2_ENCFF768QMM_peaks_id.bed
FOXM1_ENCFF711ZED_peaks_id.bed
FOXM1_ENCFF820MJD_peaks_id.bed
FOXO4_ENCFF750AGR_peaks_id.bed
FOXP1_ENCFF491EEI_peaks_id.bed
FOXP4_ENCFF974IPP_peaks_id.bed
FUS_ENCFF142CPK_peaks_id.bed
FUS_ENCFF581MLB_peaks_id.bed
FUS_ENCFF688ARM_peaks_id.bed
GABPA_ENCFF489EME_peaks_id.bed
GABPA_ENCFF889XXJ_peaks_id.bed
GABPB1_ENCFF009RFC_peaks_id.bed
GABPB1_ENCFF764HIS_peaks_id.bed
GABPB2_ENCFF103YZQ_peaks_id.bed
GATA1_ENCFF509ZLE_peaks_id.bed
GATA1_ENCFF657CTC_peaks_id.bed
GATA2_ENCFF165ZEP_peaks_id.bed
GATA2_ENCFF242YZU_peaks_id.bed
GATA2_ENCFF497ISV_peaks_id.bed
GATA2_ENCFF772OKO_peaks_id.bed
GATAD2A_ENCFF295YRT_peaks_id.bed
GATAD2B_ENCFF388KGW_peaks_id.bed
GMEB1_ENCFF365ETH_peaks_id.bed
GMEB1_ENCFF863WFD_peaks_id.bed
GTF2A2_ENCFF920CRL_peaks_id.bed
GTF2B_ENCFF987FJR_peaks_id.bed
GTF2E2_ENCFF564KIU_peaks_id.bed
GTF2F1_ENCFF075BRS_peaks_id.bed
GTF2F1_ENCFF225ZPU_peaks_id.bed
GTF2F1_ENCFF246VJH_peaks_id.bed
GTF2F1_ENCFF837ZJJ_peaks_id.bed
GTF2F1_ENCFF843UHP_peaks_id.bed
GTF2I_ENCFF866OZW_peaks_id.bed
GTF3C2_ENCFF496PIQ_peaks_id.bed
HBP1_ENCFF681YLP_peaks_id.bed
HCFC1_ENCFF139FMU_peaks_id.bed
HDAC1_ENCFF245OJC_peaks_id.bed
HDAC1_ENCFF432KJA_peaks_id.bed
HDAC1_ENCFF507QLB_peaks_id.bed
HDAC1_ENCFF669MJX_peaks_id.bed
HDAC2_ENCFF150FJT_peaks_id.bed
HDAC2_ENCFF182MPT_peaks_id.bed
HDAC2_ENCFF630HGY_peaks_id.bed
HDAC2_ENCFF652ZZF_peaks_id.bed
HDAC3_ENCFF775PIY_peaks_id.bed
HDAC6_ENCFF736FGN_peaks_id.bed
HDAC8_ENCFF566PEY_peaks_id.bed
HDAC8_ENCFF611FTI_peaks_id.bed
HDGF_ENCFF621FPC_peaks_id.bed
HDGF_ENCFF633XFQ_peaks_id.bed
HES1_ENCFF671ZGP_peaks_id.bed
HEY1_ENCFF180KAV_peaks_id.bed
HINFP_ENCFF826QOG_peaks_id.bed
HIVEP1_ENCFF449JMB_peaks_id.bed
HLTF_ENCFF614IBI_peaks_id.bed
HMBOX1_ENCFF558DSF_peaks_id.bed
HMBOX1_ENCFF672ZQW_peaks_id.bed
HMBOX1_ENCFF718DFX_peaks_id.bed
HMG20A_ENCFF118ATD_peaks_id.bed
HMG20A_ENCFF866ARI_peaks_id.bed
HMG20B_ENCFF099OSI_peaks_id.bed
HMGN3_ENCFF615AKB_peaks_id.bed
HNRNPH1_ENCFF292JRY_peaks_id.bed
HNRNPH1_ENCFF797TVO_peaks_id.bed
HNRNPH1_ENCFF844QFF_peaks_id.bed
HNRNPK_ENCFF505RNR_peaks_id.bed
HNRNPL_ENCFF010STZ_peaks_id.bed
HNRNPL_ENCFF854WAP_peaks_id.bed
HNRNPL_ENCFF984ESZ_peaks_id.bed
HNRNPLL_ENCFF122QSN_peaks_id.bed
HNRNPLL_ENCFF662WPN_peaks_id.bed
HNRNPLL_ENCFF833ZNA_peaks_id.bed
HNRNPUL1_ENCFF700VSW_peaks_id.bed
HNRNPUL1_ENCFF768TJI_peaks_id.bed
HNRNPUL1_ENCFF991ZSC_peaks_id.bed
HOMEZ_ENCFF601LMD_peaks_id.bed
HOXB6_ENCFF285MHB_peaks_id.bed
HSF4_ENCFF023JGT_peaks_id.bed
ID3_ENCFF375SIS_peaks_id.bed
ID3_ENCFF794SCQ_peaks_id.bed
IFI16_ENCFF518SJY_peaks_id.bed
IKZF1_ENCFF637SIR_peaks_id.bed
IKZF1_ENCFF711ILQ_peaks_id.bed
ILF3_ENCFF382CUT_peaks_id.bed
ILK_ENCFF143MEF_peaks_id.bed
IRF1_ENCFF093BBG_peaks_id.bed
IRF1_ENCFF557FUM_peaks_id.bed
IRF2_ENCFF430YXJ_peaks_id.bed
IRF2_ENCFF559ODJ_peaks_id.bed
IRF9_ENCFF724CHN_peaks_id.bed
JUNB_ENCFF328OCA_peaks_id.bed
JUNB_ENCFF362UGH_peaks_id.bed
JUNB_ENCFF932KRN_peaks_id.bed
JUND_ENCFF273KIA_peaks_id.bed
JUND_ENCFF306SZL_peaks_id.bed
JUN_ENCFF190CGV_peaks_id.bed
JUN_ENCFF589QXC_peaks_id.bed
JUN_ENCFF865UPM_peaks_id.bed
KAT2B_ENCFF349VSP_peaks_id.bed
KAT7_ENCFF555IKI_peaks_id.bed
KAT8_ENCFF113LBV_peaks_id.bed
KDM1A_ENCFF054XCG_peaks_id.bed
KDM1A_ENCFF346UGW_peaks_id.bed
KDM1A_ENCFF829FZW_peaks_id.bed
KDM2B_ENCFF702THZ_peaks_id.bed
KDM4B_ENCFF304PGC_peaks_id.bed
KDM4B_ENCFF928UFW_peaks_id.bed
KDM5B_ENCFF807FNB_peaks_id.bed
KHSRP_ENCFF525XXS_peaks_id.bed
KLF10_ENCFF142ZTD_peaks_id.bed
KLF13_ENCFF453MMH_peaks_id.bed
KLF16_ENCFF488OTN_peaks_id.bed
KLF1_ENCFF674KVR_peaks_id.bed
KLF6_ENCFF563ZJM_peaks_id.bed
L3MBTL2_ENCFF063NIH_peaks_id.bed
LARP7_ENCFF304PQQ_peaks_id.bed
LCOR_ENCFF266MUT_peaks_id.bed
LEF1_ENCFF591TLJ_peaks_id.bed
LEF1_ENCFF689HWD_peaks_id.bed
MAFF_ENCFF119AHD_peaks_id.bed
MAFG_ENCFF493DXA_peaks_id.bed
MAFK_ENCFF439TJM_peaks_id.bed
MAX_ENCFF171FVU_peaks_id.bed
MAX_ENCFF266YHW_peaks_id.bed
MAX_ENCFF493DCP_peaks_id.bed
MAX_ENCFF822FKQ_peaks_id.bed
MAZ_ENCFF389FLV_peaks_id.bed
MAZ_ENCFF704UVS_peaks_id.bed
MAZ_ENCFF837NNR_peaks_id.bed
MBD1_ENCFF313VEF_peaks_id.bed
MBD2_ENCFF788YHU_peaks_id.bed
MCM2_ENCFF147IAH_peaks_id.bed
MCM2_ENCFF239WEU_peaks_id.bed
MCM3_ENCFF583XAW_peaks_id.bed
MCM5_ENCFF282KHW_peaks_id.bed
MCM5_ENCFF344LTH_peaks_id.bed
MCM7_ENCFF239OHL_peaks_id.bed
MCM7_ENCFF519LSA_peaks_id.bed
MCM7_ENCFF674ZQI_peaks_id.bed
MECOM_ENCFF188PLS_peaks_id.bed
MEF2A_ENCFF031NTF_peaks_id.bed
MEF2D_ENCFF242ULW_peaks_id.bed
MEIS2_ENCFF861ZJL_peaks_id.bed
MGA_ENCFF524ZER_peaks_id.bed
MIER1_ENCFF562XKK_peaks_id.bed
MITF_ENCFF076GSK_peaks_id.bed
MITF_ENCFF605XYC_peaks_id.bed
MLLT1_ENCFF355OYM_peaks_id.bed
MLLT1_ENCFF556MUV_peaks_id.bed
MLX_ENCFF128YWS_peaks_id.bed
MLX_ENCFF171ZNN_peaks_id.bed
MNT_ENCFF440AQF_peaks_id.bed
MNT_ENCFF461QQT_peaks_id.bed
MNT_ENCFF986SXD_peaks_id.bed
MTA1_ENCFF058ZHN_peaks_id.bed
MTA2_ENCFF080IBW_peaks_id.bed
MTA2_ENCFF478LQC_peaks_id.bed
MTA3_ENCFF280AIK_peaks_id.bed
MTF1_ENCFF478MAJ_peaks_id.bed
MXD1_ENCFF559MTF_peaks_id.bed
MXI1_ENCFF068IGH_peaks_id.bed
MYBL2_ENCFF227LKX_peaks_id.bed
MYBL2_ENCFF637PFI_peaks_id.bed
MYC_ENCFF114VAI_peaks_id.bed
MYC_ENCFF566CTX_peaks_id.bed
MYC_ENCFF608CXN_peaks_id.bed
MYC_ENCFF675PHY_peaks_id.bed
MYNN_ENCFF516ZEQ_peaks_id.bed
MYNN_ENCFF603EAX_peaks_id.bed
MYNN_ENCFF749GMT_peaks_id.bed
NBN_ENCFF670POK_peaks_id.bed
NCOA1_ENCFF514WJW_peaks_id.bed
NCOA1_ENCFF660DEP_peaks_id.bed
NCOA1_ENCFF904YZI_peaks_id.bed
NCOA2_ENCFF057SID_peaks_id.bed
NCOA2_ENCFF872RQG_peaks_id.bed
NCOA4_ENCFF067BCD_peaks_id.bed
NCOA6_ENCFF996OEL_peaks_id.bed
NCOR1_ENCFF321NXM_peaks_id.bed
NCOR1_ENCFF816AEF_peaks_id.bed
NCOR1_ENCFF836PGX_peaks_id.bed
NCOR1_ENCFF898CGY_peaks_id.bed
NELFE_ENCFF265KKE_peaks_id.bed
NEUROD1_ENCFF625QHR_peaks_id.bed
NFATC3_ENCFF778UHP_peaks_id.bed
NFATC3_ENCFF972OBR_peaks_id.bed
NFE2_ENCFF023IFO_peaks_id.bed
NFE2_ENCFF474PRJ_peaks_id.bed
NFE2_ENCFF714MVQ_peaks_id.bed
NFE2L1_ENCFF366XDW_peaks_id.bed
NFIC_ENCFF169MAE_peaks_id.bed
NFIX_ENCFF726LLI_peaks_id.bed
NFRKB_ENCFF057NNR_peaks_id.bed
NFRKB_ENCFF321LMQ_peaks_id.bed
NFXL1_ENCFF201WFM_peaks_id.bed
NFXL1_ENCFF341CNM_peaks_id.bed
NFYA_ENCFF786CRW_peaks_id.bed
NFYA_ENCFF908HSL_peaks_id.bed
NFYB_ENCFF718ZFY_peaks_id.bed
NKRF_ENCFF068MGV_peaks_id.bed
NONO_ENCFF134NMW_peaks_id.bed
NONO_ENCFF211TTD_peaks_id.bed
NONO_ENCFF318GYM_peaks_id.bed
NONO_ENCFF823CQK_peaks_id.bed
NR0B1_ENCFF931LUC_peaks_id.bed
NR1H2_ENCFF798KXN_peaks_id.bed
NR2C1_ENCFF199GFU_peaks_id.bed
NR2C1_ENCFF469ZBB_peaks_id.bed
NR2C2_ENCFF263VIC_peaks_id.bed
NR2C2_ENCFF495XYS_peaks_id.bed
NR2C2_ENCFF946CTR_peaks_id.bed
NR2F1_ENCFF746GDG_peaks_id.bed
NR2F2_ENCFF847ZHF_peaks_id.bed
NR2F6_ENCFF529SRS_peaks_id.bed
NR2F6_ENCFF647RGI_peaks_id.bed
NR2F6_ENCFF835FPJ_peaks_id.bed
NR3C1_ENCFF201BGD_peaks_id.bed
NR3C1_ENCFF960SMD_peaks_id.bed
NR3C1_ENCFF961ADR_peaks_id.bed
NR4A1_ENCFF137TCL_peaks_id.bed
NR4A1_ENCFF772UYI_peaks_id.bed
NR4A1_ENCFF852ZIK_peaks_id.bed
NRF1_ENCFF259YUE_peaks_id.bed
NRF1_ENCFF410RJD_peaks_id.bed
NRF1_ENCFF762CYI_peaks_id.bed
NRF1_ENCFF777PKJ_peaks_id.bed
NUFIP1_ENCFF669RLC_peaks_id.bed
OTX1_ENCFF245JPK_peaks_id.bed
PATZ1_ENCFF743JJF_peaks_id.bed
PBX2_ENCFF345BUG_peaks_id.bed
PBX2_ENCFF841RMH_peaks_id.bed
PCBP1_ENCFF467RYH_peaks_id.bed
PCBP1_ENCFF941RVL_peaks_id.bed
PCBP1_ENCFF956MGE_peaks_id.bed
PCBP2_ENCFF732IRT_peaks_id.bed
PCBP2_ENCFF881XQF_peaks_id.bed
PCBP2_ENCFF941XZW_peaks_id.bed
PHB2_ENCFF785QJM_peaks_id.bed
PHB_ENCFF424OWM_peaks_id.bed
PHF20_ENCFF053PDX_peaks_id.bed
PHF21A_ENCFF215JWS_peaks_id.bed
PHF8_ENCFF981ISM_peaks_id.bed
PHTF2_ENCFF130WQT_peaks_id.bed
PKNOX1_ENCFF099RDJ_peaks_id.bed
PML_ENCFF051FNO_peaks_id.bed
POLR2A_ENCFF107SJD_peaks_id.bed
POLR2A_ENCFF355MNE_peaks_id.bed
POLR2A_ENCFF410CPY_peaks_id.bed
POLR2A_ENCFF634JRD_peaks_id.bed
POLR2A_ENCFF842JME_peaks_id.bed
POLR2A_ENCFF921FKB_peaks_id.bed
POLR2AphosphoS2_ENCFF950YZE_peaks_id.bed
POLR2AphosphoS2_ENCFF951KHS_peaks_id.bed
POLR2AphosphoS5_ENCFF060MMW_peaks_id.bed
POLR2B_ENCFF821QOS_peaks_id.bed
POLR2G_ENCFF273EYJ_peaks_id.bed
POLR2G_ENCFF283CUY_peaks_id.bed
POLR2G_ENCFF413LLO_peaks_id.bed
POLR2H_ENCFF084KHS_peaks_id.bed
POLR3A_ENCFF627DJP_peaks_id.bed
POLR3G_ENCFF307GHJ_peaks_id.bed
POU5F1_ENCFF724DNU_peaks_id.bed
PPARD_ENCFF893YKY_peaks_id.bed
PRDM10_ENCFF405BKB_peaks_id.bed
PRMT5_ENCFF018TNP_peaks_id.bed
PRPF4_ENCFF417RQZ_peaks_id.bed
PRPF4_ENCFF785ACI_peaks_id.bed
PRPF4_ENCFF886UMM_peaks_id.bed
PTBP1_ENCFF694BIA_peaks_id.bed
PTBP1_ENCFF835NOD_peaks_id.bed
PTBP1_ENCFF917HXV_peaks_id.bed
PTRF_ENCFF491UBF_peaks_id.bed
PTTG1_ENCFF918BZG_peaks_id.bed
PURB_ENCFF500BWO_peaks_id.bed
PYGO2_ENCFF232TVE_peaks_id.bed
PYGO2_ENCFF335VOJ_peaks_id.bed
RAD21_ENCFF057JFH_peaks_id.bed
RAD21_ENCFF258VXX_peaks_id.bed
RAD21_ENCFF330SHG_peaks_id.bed
RAD21_ENCFF930WPG_peaks_id.bed
RAD51_ENCFF536WBT_peaks_id.bed
RB1_ENCFF479SOJ_peaks_id.bed
RBBP5_ENCFF942VGF_peaks_id.bed
RBFOX2_ENCFF120IDE_peaks_id.bed
RBFOX2_ENCFF232ASB_peaks_id.bed
RBFOX2_ENCFF538ACI_peaks_id.bed
RBM14_ENCFF092DLN_peaks_id.bed
RBM14_ENCFF465UMU_peaks_id.bed
RBM14_ENCFF682SIY_peaks_id.bed
RBM15_ENCFF091TCH_peaks_id.bed
RBM15_ENCFF563WDZ_peaks_id.bed
RBM15_ENCFF971VJZ_peaks_id.bed
RBM17_ENCFF056OIG_peaks_id.bed
RBM17_ENCFF861YKK_peaks_id.bed
RBM17_ENCFF937QYU_peaks_id.bed
RBM22_ENCFF420IBN_peaks_id.bed
RBM22_ENCFF522JUV_peaks_id.bed
RBM22_ENCFF987DPX_peaks_id.bed
RBM25_ENCFF102XVH_peaks_id.bed
RBM25_ENCFF680WBN_peaks_id.bed
RBM25_ENCFF849VEO_peaks_id.bed
RBM34_ENCFF071GJH_peaks_id.bed
RBM34_ENCFF670ILH_peaks_id.bed
RBM34_ENCFF782GWS_peaks_id.bed
RBM39_ENCFF209JJD_peaks_id.bed
RBM39_ENCFF503DIK_peaks_id.bed
RBM39_ENCFF983LFS_peaks_id.bed
RBPJ_ENCFF691JZW_peaks_id.bed
RCOR1_ENCFF224YBQ_peaks_id.bed
RCOR1_ENCFF627PLM_peaks_id.bed
RELA_ENCFF936WES_peaks_id.bed
REST_ENCFF118ECK_peaks_id.bed
REST_ENCFF539MIO_peaks_id.bed
REST_ENCFF707MDI_peaks_id.bed
REST_ENCFF761YYL_peaks_id.bed
RFX1_ENCFF768VRG_peaks_id.bed
RFX1_ENCFF820GPR_peaks_id.bed
RFX5_ENCFF661EUX_peaks_id.bed
RFX7_ENCFF774NZZ_peaks_id.bed
RHOXF2B_ENCFF727JXN_peaks_id.bed
RLF_ENCFF752VQB_peaks_id.bed
RNF219_ENCFF417FFL_peaks_id.bed
RNF2_ENCFF060JVZ_peaks_id.bed
RNF2_ENCFF271MJE_peaks_id.bed
RNF2_ENCFF408XPG_peaks_id.bed
RNF2_ENCFF409SDU_peaks_id.bed
RNF2_ENCFF923GDR_peaks_id.bed
RREB1_ENCFF057RJK_peaks_id.bed
RUNX1_ENCFF003LPE_peaks_id.bed
RUNX1_ENCFF374EFU_peaks_id.bed
SAFB2_ENCFF087DKT_peaks_id.bed
SAFB2_ENCFF176MDV_peaks_id.bed
SAFB2_ENCFF624NUZ_peaks_id.bed
SAFB_ENCFF411YVY_peaks_id.bed
SAFB_ENCFF478MPX_peaks_id.bed
SAFB_ENCFF537RNU_peaks_id.bed
SAP30_ENCFF383IEP_peaks_id.bed
SETDB1_ENCFF267QRI_peaks_id.bed
SETDB1_ENCFF319ZYX_peaks_id.bed
SETDB1_ENCFF752EMN_peaks_id.bed
SFPQ_ENCFF652ZEN_peaks_id.bed
SFPQ_ENCFF757ODD_peaks_id.bed
SIN3A_ENCFF469KAH_peaks_id.bed
SIN3A_ENCFF884MNF_peaks_id.bed
SIN3B_ENCFF257NOX_peaks_id.bed
SIRT6_ENCFF217ACS_peaks_id.bed
SIRT6_ENCFF821XJU_peaks_id.bed
SIX5_ENCFF742ZHT_peaks_id.bed
SIX5_ENCFF754ZWP_peaks_id.bed
SKIL_ENCFF067KCK_peaks_id.bed
SLC30A9_ENCFF441QXF_peaks_id.bed
SMAD1_ENCFF120PGJ_peaks_id.bed
SMAD2_ENCFF730QQD_peaks_id.bed
SMAD3_ENCFF335ZTU_peaks_id.bed
SMAD3_ENCFF718DSJ_peaks_id.bed
SMAD4_ENCFF232ZHD_peaks_id.bed
SMAD4_ENCFF934NEX_peaks_id.bed
SMAD5_ENCFF971YTJ_peaks_id.bed
SMARCA4_ENCFF267OGF_peaks_id.bed
SMARCA4_ENCFF486QCV_peaks_id.bed
SMARCA4_ENCFF793FMG_peaks_id.bed
SMARCA5_ENCFF722RWS_peaks_id.bed
SMARCB1_ENCFF217LXF_peaks_id.bed
SMARCC2_ENCFF908MWB_peaks_id.bed
SMARCE1_ENCFF148YMC_peaks_id.bed
SMC3_ENCFF289LLT_peaks_id.bed
SNAPC5_ENCFF656THP_peaks_id.bed
SNIP1_ENCFF020XNM_peaks_id.bed
SNRNP70_ENCFF033KXY_peaks_id.bed
SNRNP70_ENCFF206MJS_peaks_id.bed
SNRNP70_ENCFF306TAD_peaks_id.bed
SOX6_ENCFF125GZU_peaks_id.bed
SP1_ENCFF171NEU_peaks_id.bed
SP1_ENCFF412ERW_peaks_id.bed
SP1_ENCFF553GPK_peaks_id.bed
SP2_ENCFF741FZG_peaks_id.bed
SPI1_ENCFF888CKG_peaks_id.bed
SREBF1_ENCFF029RBI_peaks_id.bed
SREBF2_ENCFF710YFH_peaks_id.bed
SRF_ENCFF087EVW_peaks_id.bed
SRF_ENCFF101SEZ_peaks_id.bed
SRF_ENCFF716OTE_peaks_id.bed
SRSF1_ENCFF792SXS_peaks_id.bed
SRSF1_ENCFF798URE_peaks_id.bed
SRSF1_ENCFF862DEO_peaks_id.bed
SRSF3_ENCFF048BKZ_peaks_id.bed
SRSF3_ENCFF722XRW_peaks_id.bed
SRSF3_ENCFF926XGK_peaks_id.bed
SRSF7_ENCFF059WVE_peaks_id.bed
SRSF7_ENCFF217MBF_peaks_id.bed
SRSF7_ENCFF550VUN_peaks_id.bed
SRSF9_ENCFF137IBM_peaks_id.bed
SRSF9_ENCFF217HAW_peaks_id.bed
SRSF9_ENCFF746IEZ_peaks_id.bed
STAG1_ENCFF921BXP_peaks_id.bed
STAT5A_ENCFF187BVL_peaks_id.bed
STAT5B_ENCFF122FTW_peaks_id.bed
SUPT5H_ENCFF738XMN_peaks_id.bed
SUZ12_ENCFF647TBC_peaks_id.bed
SUZ12_ENCFF889QYR_peaks_id.bed
TAF15_ENCFF547PES_peaks_id.bed
TAF15_ENCFF617CAZ_peaks_id.bed
TAF15_ENCFF710LLF_peaks_id.bed
TAF1_ENCFF784BTI_peaks_id.bed
TAF7_ENCFF199CQK_peaks_id.bed
TAF7_ENCFF792CKI_peaks_id.bed
TAF9B_ENCFF432HRL_peaks_id.bed
TAL1_ENCFF101DBG_peaks_id.bed
TAL1_ENCFF852ZRK_peaks_id.bed
TARDBP_ENCFF448YOS_peaks_id.bed
TARDBP_ENCFF564QOL_peaks_id.bed
TARDBP_ENCFF880BPX_peaks_id.bed
TARDBP_ENCFF905VXX_peaks_id.bed
TARDBP_ENCFF986BYT_peaks_id.bed
TBL1XR1_ENCFF034SWM_peaks_id.bed
TBL1XR1_ENCFF967URJ_peaks_id.bed
TBP_ENCFF370YGS_peaks_id.bed
TBPL1_ENCFF927NEU_peaks_id.bed
TBX18_ENCFF405EGZ_peaks_id.bed
TCF12_ENCFF504PFQ_peaks_id.bed
TCF12_ENCFF808QUD_peaks_id.bed
TCF15_ENCFF415MKA_peaks_id.bed
TCF3_ENCFF773RNU_peaks_id.bed
TCF7_ENCFF736XUU_peaks_id.bed
TCF7L2_ENCFF732ZNB_peaks_id.bed
TCFL5_ENCFF410SWS_peaks_id.bed
TEAD1_ENCFF214SNH_peaks_id.bed
TEAD1_ENCFF692HUD_peaks_id.bed
TEAD2_ENCFF209WPT_peaks_id.bed
TEAD4_ENCFF303CVC_peaks_id.bed
TEAD4_ENCFF501XJP_peaks_id.bed
TFAM_ENCFF269PAM_peaks_id.bed
TFAP4_ENCFF438IYI_peaks_id.bed
TFCP2_ENCFF392MAR_peaks_id.bed
TFDP1_ENCFF044IDV_peaks_id.bed
TFDP1_ENCFF362JHQ_peaks_id.bed
TFE3_ENCFF592NJN_peaks_id.bed
TGIF2_ENCFF238XIA_peaks_id.bed
THAP12_ENCFF206LMB_peaks_id.bed
THAP1_ENCFF823RYG_peaks_id.bed
THAP7_ENCFF372YFU_peaks_id.bed
THRA_ENCFF992HUS_peaks_id.bed
THRAP3_ENCFF067FJF_peaks_id.bed
THRB_ENCFF906PUB_peaks_id.bed
TOE1_ENCFF563WUP_peaks_id.bed
TOE1_ENCFF917HUH_peaks_id.bed
TRIM24_ENCFF527MLG_peaks_id.bed
TRIM24_ENCFF597YMB_peaks_id.bed
TRIM25_ENCFF110HVI_peaks_id.bed
TRIM25_ENCFF448LKL_peaks_id.bed
TRIM25_ENCFF514DDG_peaks_id.bed
TRIM25_ENCFF574LAO_peaks_id.bed
TRIM28_ENCFF172CRL_peaks_id.bed
TRIM28_ENCFF634OOU_peaks_id.bed
TRIM28_ENCFF770JQM_peaks_id.bed
TRIP13_ENCFF955QCD_peaks_id.bed
TSC22D4_ENCFF705SFR_peaks_id.bed
TSHZ1_ENCFF144ZLB_peaks_id.bed
U2AF1_ENCFF482DRO_peaks_id.bed
U2AF1_ENCFF730URE_peaks_id.bed
U2AF1_ENCFF752DQV_peaks_id.bed
U2AF2_ENCFF134HBP_peaks_id.bed
U2AF2_ENCFF164CTH_peaks_id.bed
U2AF2_ENCFF398EQF_peaks_id.bed
UBTF_ENCFF326CIT_peaks_id.bed
UBTF_ENCFF568ZPW_peaks_id.bed
USF1_ENCFF255QDL_peaks_id.bed
USF1_ENCFF310CCS_peaks_id.bed
USF2_ENCFF289ZIR_peaks_id.bed
USF2_ENCFF640ZIN_peaks_id.bed
USF2_ENCFF744HVD_peaks_id.bed
VEZF1_ENCFF433ERT_peaks_id.bed
WHSC1_ENCFF862UUR_peaks_id.bed
XRCC3_ENCFF137OGC_peaks_id.bed
XRCC4_ENCFF540QPM_peaks_id.bed
XRCC5_ENCFF708IYN_peaks_id.bed
XRCC5_ENCFF929TWP_peaks_id.bed
XRCC5_ENCFF998KKJ_peaks_id.bed
YBX1_ENCFF286IPW_peaks_id.bed
YBX3_ENCFF877ZKU_peaks_id.bed
YY1_ENCFF074OAM_peaks_id.bed
YY1_ENCFF398UQZ_peaks_id.bed
YY1_ENCFF589PZO_peaks_id.bed
ZBED1_ENCFF392INU_peaks_id.bed
ZBTB11_ENCFF227HVU_peaks_id.bed
ZBTB11_ENCFF539GWU_peaks_id.bed
ZBTB11_ENCFF554GSE_peaks_id.bed
ZBTB11_ENCFF991VAW_peaks_id.bed
ZBTB12_ENCFF830MTX_peaks_id.bed
ZBTB17_ENCFF292HSS_peaks_id.bed
ZBTB26_ENCFF435SNP_peaks_id.bed
ZBTB2_ENCFF209DVO_peaks_id.bed
ZBTB33_ENCFF146GZZ_peaks_id.bed
ZBTB33_ENCFF917RIN_peaks_id.bed
ZBTB34_ENCFF683VXW_peaks_id.bed
ZBTB40_ENCFF132DPY_peaks_id.bed
ZBTB40_ENCFF624PTB_peaks_id.bed
ZBTB40_ENCFF803CYY_peaks_id.bed
ZBTB49_ENCFF308VAI_peaks_id.bed
ZBTB5_ENCFF079EGS_peaks_id.bed
ZBTB5_ENCFF460GXA_peaks_id.bed
ZBTB7A_ENCFF346AOR_peaks_id.bed
ZBTB8A_ENCFF328SSL_peaks_id.bed
ZBTB8A_ENCFF423CWQ_peaks_id.bed
ZBTB8A_ENCFF589EVD_peaks_id.bed
ZBTB9_ENCFF588CXX_peaks_id.bed
ZC3H11A_ENCFF528SIV_peaks_id.bed
ZC3H4_ENCFF270TMW_peaks_id.bed
ZC3H8_ENCFF045AOZ_peaks_id.bed
ZC3H8_ENCFF230DZT_peaks_id.bed
ZEB2_ENCFF242AOL_peaks_id.bed
ZEB2_ENCFF475CIS_peaks_id.bed
ZFP1_ENCFF795AWO_peaks_id.bed
ZFP30_ENCFF141HPL_peaks_id.bed
ZFP36_ENCFF839GAS_peaks_id.bed
ZFP91_ENCFF718IGO_peaks_id.bed
ZFP91_ENCFF835PDO_peaks_id.bed
ZFPM2_ENCFF549NPY_peaks_id.bed
ZFX_ENCFF840NZE_peaks_id.bed
ZFX_ENCFF910JTR_peaks_id.bed
ZHX1_ENCFF175YPM_peaks_id.bed
ZKSCAN1_ENCFF560UGR_peaks_id.bed
ZKSCAN3_ENCFF697DRN_peaks_id.bed
ZKSCAN8_ENCFF908SNB_peaks_id.bed
ZMIZ1_ENCFF316IYH_peaks_id.bed
ZMIZ1_ENCFF881DAT_peaks_id.bed
ZMYM3_ENCFF117XRE_peaks_id.bed
ZNF121_ENCFF650DWZ_peaks_id.bed
ZNF124_ENCFF881MMD_peaks_id.bed
ZNF12_ENCFF021JCJ_peaks_id.bed
ZNF133_ENCFF810OHB_peaks_id.bed
ZNF134_ENCFF169AWH_peaks_id.bed
ZNF140_ENCFF207XCM_peaks_id.bed
ZNF143_ENCFF540TWL_peaks_id.bed
ZNF143_ENCFF596HNP_peaks_id.bed
ZNF146_ENCFF296DCP_peaks_id.bed
ZNF146_ENCFF742XBE_peaks_id.bed
ZNF146_ENCFF932FVX_peaks_id.bed
ZNF148_ENCFF283UWH_peaks_id.bed
ZNF165_ENCFF169QYL_peaks_id.bed
ZNF174_ENCFF822FXR_peaks_id.bed
ZNF175_ENCFF192MEM_peaks_id.bed
ZNF184_ENCFF058VGZ_peaks_id.bed
ZNF184_ENCFF266FUJ_peaks_id.bed
ZNF184_ENCFF882XTP_peaks_id.bed
ZNF197_ENCFF987UBO_peaks_id.bed
ZNF212_ENCFF435MHH_peaks_id.bed
ZNF215_ENCFF922RHM_peaks_id.bed
ZNF217_ENCFF515LWL_peaks_id.bed
ZNF224_ENCFF905RSN_peaks_id.bed
ZNF232_ENCFF672PJF_peaks_id.bed
ZNF239_ENCFF569TJP_peaks_id.bed
ZNF23_ENCFF063VRD_peaks_id.bed
ZNF24_ENCFF225GCF_peaks_id.bed
ZNF24_ENCFF268YOM_peaks_id.bed
ZNF24_ENCFF649OPX_peaks_id.bed
ZNF24_ENCFF683VTL_peaks_id.bed
ZNF253_ENCFF019NGO_peaks_id.bed
ZNF257_ENCFF738IQL_peaks_id.bed
ZNF263_ENCFF264BDD_peaks_id.bed
ZNF263_ENCFF295XBK_peaks_id.bed
ZNF274_ENCFF045CFL_peaks_id.bed
ZNF274_ENCFF742DPO_peaks_id.bed
ZNF277_ENCFF382QUI_peaks_id.bed
ZNF280A_ENCFF932ZCO_peaks_id.bed
ZNF280B_ENCFF253UAF_peaks_id.bed
ZNF281_ENCFF648ORA_peaks_id.bed
ZNF282_ENCFF110GDP_peaks_id.bed
ZNF311_ENCFF080KWE_peaks_id.bed
ZNF316_ENCFF451AEQ_peaks_id.bed
ZNF316_ENCFF840SLB_peaks_id.bed
ZNF318_ENCFF545WAL_peaks_id.bed
ZNF318_ENCFF990SFN_peaks_id.bed
ZNF319_ENCFF873VFI_peaks_id.bed
ZNF324_ENCFF265MQC_peaks_id.bed
ZNF347_ENCFF843NBV_peaks_id.bed
ZNF354B_ENCFF649UCX_peaks_id.bed
ZNF354B_ENCFF712NHB_peaks_id.bed
ZNF354C_ENCFF604RBX_peaks_id.bed
ZNF384_ENCFF864XZP_peaks_id.bed
ZNF395_ENCFF430IPX_peaks_id.bed
ZNF397_ENCFF193BGB_peaks_id.bed
ZNF398_ENCFF380KJB_peaks_id.bed
ZNF3_ENCFF180XUM_peaks_id.bed
ZNF3_ENCFF510PVQ_peaks_id.bed
ZNF3_ENCFF570YIC_peaks_id.bed
ZNF3_ENCFF798AFV_peaks_id.bed
ZNF407_ENCFF015FXW_peaks_id.bed
ZNF407_ENCFF255HEF_peaks_id.bed
ZNF408_ENCFF207MLR_peaks_id.bed
ZNF41_ENCFF629EAD_peaks_id.bed
ZNF431_ENCFF442RYH_peaks_id.bed
ZNF436_ENCFF318ZUN_peaks_id.bed
ZNF444_ENCFF295XCB_peaks_id.bed
ZNF445_ENCFF416TFM_peaks_id.bed
ZNF449_ENCFF443BGE_peaks_id.bed
ZNF507_ENCFF072JDK_peaks_id.bed
ZNF507_ENCFF878SVX_peaks_id.bed
ZNF511_ENCFF653ABW_peaks_id.bed
ZNF512_ENCFF184IOY_peaks_id.bed
ZNF512_ENCFF817EJA_peaks_id.bed
ZNF518B_ENCFF096AYV_peaks_id.bed
ZNF551_ENCFF913YMX_peaks_id.bed
ZNF561_ENCFF701RWL_peaks_id.bed
ZNF57_ENCFF945IST_peaks_id.bed
ZNF583_ENCFF053ABJ_peaks_id.bed
ZNF584_ENCFF005MBI_peaks_id.bed
ZNF586_ENCFF887ZDT_peaks_id.bed
ZNF589_ENCFF345IHK_peaks_id.bed
ZNF592_ENCFF847QJI_peaks_id.bed
ZNF639_ENCFF217EJG_peaks_id.bed
ZNF639_ENCFF643AUL_peaks_id.bed
ZNF639_ENCFF645OCW_peaks_id.bed
ZNF639_ENCFF669VJB_peaks_id.bed
ZNF644_ENCFF773XPT_peaks_id.bed
ZNF655_ENCFF156EZP_peaks_id.bed
ZNF668_ENCFF885HJI_peaks_id.bed
ZNF695_ENCFF493LZB_peaks_id.bed
ZNF696_ENCFF899BJR_peaks_id.bed
ZNF699_ENCFF854SNX_peaks_id.bed
ZNF700_ENCFF208OME_peaks_id.bed
ZNF707_ENCFF959ZKZ_peaks_id.bed
ZNF717_ENCFF159END_peaks_id.bed
ZNF740_ENCFF004YCK_peaks_id.bed
ZNF740_ENCFF219NIA_peaks_id.bed
ZNF740_ENCFF959SMC_peaks_id.bed
ZNF75A_ENCFF307UZC_peaks_id.bed
ZNF764_ENCFF049MDN_peaks_id.bed
ZNF766_ENCFF466PKS_peaks_id.bed
ZNF76_ENCFF305XQM_peaks_id.bed
ZNF778_ENCFF165JQS_peaks_id.bed
ZNF77_ENCFF809JNC_peaks_id.bed
ZNF780A_ENCFF794OSY_peaks_id.bed
ZNF785_ENCFF352NGK_peaks_id.bed
ZNF79_ENCFF670RLH_peaks_id.bed
ZNF7_ENCFF018HWM_peaks_id.bed
ZNF830_ENCFF150ZBY_peaks_id.bed
ZNF830_ENCFF672NBD_peaks_id.bed
ZNF830_ENCFF896IUI_peaks_id.bed
ZNF83_ENCFF451LLC_peaks_id.bed
ZNF84_ENCFF014HYS_peaks_id.bed
ZSCAN29_ENCFF151WYQ_peaks_id.bed
ZSCAN29_ENCFF407STM_peaks_id.bed
ZSCAN32_ENCFF537HHU_peaks_id.bed
ZZZ3_ENCFF797VEK_peaks_id.bed
''').strip().split('\n')

#Add only the the files corrected for id and scaffolds below i.e. _chr.txt from public and new cutnrun data. For chip-seq it is ok not to add _chr.txt 
#some files which needs to be intersect with itself will be generated above and there _chr.txt form is added below
list2_data = ('''
CRAMP1_NTC_IgG_NTC_peaks.narrowPeak_chr.txt
CRAMP1_shSUZ12_IgG_shSUZ12_peaks.narrowPeak_chr.txt
V5-CRAMP1_V5_control_peaks.narrowPeak_chr.txt
SETDB1_ENCFF319ZYX_peaks_id.bed_chr.txt
SUZ12_NTC_IgG_NTC_peaks.narrowPeak_chr.txt
SUZ12_ENCFF647TBC_peaks_id.bed_chr.txt
H3K9me3_NTC_IgG_NTC_peaks.narrowPeak_chr.txt
H3K27me3_NTC_IgG_NTC_peaks.narrowPeak_chr.txt
H3K27ac_peaks.broadPeak.chr.txt
H3K27me3_k562_intersected_peaks.broadPeak.chr.txt
H3K27me3_noInput_peaks.broadPeak.chr.txt
H3K36me3_peaks.broadPeak.chr.txt
H3K4me1_peaks.broadPeak.chr.txt
H3K4me2_peaks.broadPeak.chr.txt
H3K4me3_peaks.broadPeak.chr.txt
H3K79me2_peaks.broadPeak.chr.txt
H3K9ac_peaks.broadPeak.chr.txt
H3K9me3_peaks.broadPeak.chr.txt
H4K16ac_peaks.broadPeak.chr.txt
H2AFZ_ENCFF213OTI_peaks_id.bed
H3K27ac_ENCFF544LXB_peaks_id.bed
H3K27me3_ENCFF323WOT_peaks_id.bed
H3K27me3_ENCFF801AHF_peaks_id.bed
H3K36me3_ENCFF193ERO_peaks_id.bed
H3K36me3_ENCFF561OUZ_peaks_id.bed
H3K4me1_ENCFF135ZLM_peaks_id.bed
H3K4me1_ENCFF540NGG_peaks_id.bed
H3K4me2_ENCFF749KLQ_peaks_id.bed
H3K4me3_ENCFF122CSI_peaks_id.bed
H3K4me3_ENCFF689QIJ_peaks_id.bed
H3K4me3_ENCFF706WUF_peaks_id.bed
H3K4me3_ENCFF885FQN_peaks_id.bed
H3K79me2_ENCFF209OQD_peaks_id.bed
H3K9ac_ENCFF148UQI_peaks_id.bed
H3K9ac_ENCFF891CHI_peaks_id.bed
H3K9me1_ENCFF462AVD_peaks_id.bed
H3K9me3_ENCFF963GZJ_peaks_id.bed
H4K20me1_ENCFF909RKY_peaks_id.bed
''').strip().split('\n')

for target in list1_data:

    if os.path.getsize(target + '_chr.txt') > 0:
        mypeakfile = open( target + '_chr.txt', 'r')
        mypeakfile = mypeakfile.read().strip().split('\n')
        peaklist = []
        for peakid in mypeakfile:
            peakid = peakid.split('\t')
            peaklist.append(peakid[3])
        uniq_peaklist = set(peaklist)
        #print(uniq_peaklist)
        #print(len(peaklist))
        #print(len(uniq_peaklist))

        for histone in list2_data:
            print(target + '_' + histone)
            os.system('bedtools intersect -a ' + target + '_chr.txt -b ' + histone + ' > ' + target + '_' + histone + '_intersect.txt')
            if os.path.getsize(target + '_' + histone + '_intersect.txt') > 0:
                myintersectsfile = open(target + '_' + histone + '_intersect.txt', 'r')
                myintersectsfile = myintersectsfile.read().strip().split('\n')
                myintersectlist = []
                for myintersect in myintersectsfile:
                    myintersect = myintersect.split('\t')
                    myintersectlist.append(myintersect[3])
                uniq_myintersectlist = set(myintersectlist)
                #print(uniq_myintersectlist)
                #print(len(myintersectlist))
                #print(len(uniq_myintersectlist))
                mycount = len(uniq_myintersectlist) / len(uniq_peaklist)
                #print(mycount)
                if os.path.exists('peak_counter_tf_and_histone_out.txt'):
                    os.remove('peak_counter_tf_and_histone_out.txt')
                    with open("peak_counter_tf_and_histone_out.txt", "a") as myintersectedfiles:
                        myintersectedfiles.write(''.join(target + '\t' + histone + '\t'  + str(len(uniq_myintersectlist)) + '\t' + str(len(uniq_peaklist)) + '\t' + str(mycount) + '\n'))
                        myintersectedfiles.close()
                    print(target + '\t' + histone + '\t'  + str(len(uniq_myintersectlist)) + '\t' + str(len(uniq_peaklist)) + '\t' + str(mycount))
print('Analysis finished')


