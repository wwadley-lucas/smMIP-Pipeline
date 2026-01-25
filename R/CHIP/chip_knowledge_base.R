# CHIP Knowledge Base
# Gene-level clinical reference data for CHIP mutations
# Based on published literature and clinical guidelines

# ============================================================================
# CHIP GENE INFORMATION DATABASE
# ============================================================================

#' Get comprehensive CHIP gene information
#'
#' Returns a list containing clinical information for all genes in the
#' CHIP/myeloid panel, including function, clinical significance, and references.
#'
#' @return Named list of gene information
get_chip_gene_info <- function() {
  list(
    # ========================================================================
    # DNA METHYLATION / EPIGENETIC REGULATORS
    # ========================================================================
    DNMT3A = list(
      full_name = "DNA Methyltransferase 3 Alpha",
      gene_function = "Encodes a DNA methyltransferase that establishes de novo DNA methylation patterns. Critical for hematopoietic stem cell self-renewal and differentiation.",
      pathway = "Epigenetic regulation / DNA methylation",
      chip_frequency = "Most common (~40% of CHIP cases)",
      risk_category = "Intermediate",
      clinical_significance = "Associated with increased cardiovascular disease risk (2-4 fold), modest AML progression risk (~0.5-1% per year). Single mutations generally lower risk than multiple mutations or co-occurring mutations.",
      hotspot_mutations = c("R882H", "R882C", "R882S"),
      hotspot_note = "R882 codon mutations are most common (~60% of DNMT3A mutations), cause dominant-negative effects",
      cardiovascular_risk = "Elevated risk of coronary heart disease, early-onset myocardial infarction, ischemic stroke",
      progression_to_malignancy = "~0.5-1% annual risk of progression to hematologic malignancy",
      monitoring_recommendation = "Annual CBC with differential; consider echocardiography for cardiovascular assessment",
      references = "Jaiswal et al., NEJM 2014; Genovese et al., NEJM 2014; Buscarlet et al., Blood 2017"
    ),

    TET2 = list(
      full_name = "Tet Methylcytosine Dioxygenase 2",
      gene_function = "Encodes enzyme that converts 5-methylcytosine to 5-hydroxymethylcytosine. Involved in DNA demethylation and regulation of hematopoietic stem cell function.",
      pathway = "Epigenetic regulation / DNA demethylation",
      chip_frequency = "Second most common (~20-25% of CHIP cases)",
      risk_category = "Intermediate",
      clinical_significance = "Associated with increased cardiovascular disease risk and inflammation. Loss-of-function mutations promote clonal expansion. May enhance response to vitamin C therapy.",
      hotspot_mutations = c("Frameshift mutations throughout gene"),
      hotspot_note = "No single hotspot; loss-of-function mutations distributed across gene",
      cardiovascular_risk = "Elevated risk of atherosclerosis, heart failure; linked to IL-1beta-mediated inflammation",
      progression_to_malignancy = "~0.5-1% annual risk; higher when combined with other mutations",
      monitoring_recommendation = "Annual CBC; cardiovascular risk assessment; consider inflammatory markers (hsCRP)",
      references = "Fuster et al., Science 2017; Jaiswal et al., NEJM 2017; Cimmino et al., Cell 2017"
    ),

    ASXL1 = list(
      full_name = "ASXL Transcriptional Regulator 1",
      gene_function = "Chromatin-binding protein that regulates gene expression through interaction with polycomb repressive complex 2 (PRC2). Involved in histone modification.",
      pathway = "Epigenetic regulation / Chromatin modification",
      chip_frequency = "Common (~10-15% of CHIP cases)",
      risk_category = "High",
      clinical_significance = "Frameshift mutations (especially in exon 12) are pathogenic. Higher risk of progression to myeloid malignancy compared to DNMT3A/TET2. Associated with adverse prognosis in MDS/AML.",
      hotspot_mutations = c("G646Wfs*12", "E635Rfs*15", "Y591*"),
      hotspot_note = "Most mutations are truncating mutations in exon 12",
      cardiovascular_risk = "Elevated cardiovascular risk similar to other CHIP genes",
      progression_to_malignancy = "Higher risk (~1-2% annual); particularly concerning when combined with other mutations",
      monitoring_recommendation = "More frequent monitoring (every 6 months); prompt evaluation if cytopenias develop",
      references = "Abdel-Wahab et al., Nature Genetics 2012; Inoue et al., Cancer Cell 2013"
    ),

    IDH1 = list(
      full_name = "Isocitrate Dehydrogenase 1",
      gene_function = "Cytoplasmic enzyme that converts isocitrate to alpha-ketoglutarate. Mutations produce oncometabolite 2-hydroxyglutarate (2-HG) that inhibits TET2.",
      pathway = "Epigenetic regulation / Metabolism",
      chip_frequency = "Uncommon in CHIP (~1-2%)",
      risk_category = "High",
      clinical_significance = "Gain-of-function mutations produce 2-HG oncometabolite. More commonly seen in AML than CHIP. FDA-approved targeted therapy available (ivosidenib).",
      hotspot_mutations = c("R132H", "R132C", "R132S"),
      hotspot_note = "All pathogenic mutations occur at R132 codon",
      cardiovascular_risk = "Less well-characterized than DNMT3A/TET2",
      progression_to_malignancy = "Higher progression risk; therapeutic target in AML",
      monitoring_recommendation = "Close monitoring; consider early referral to hematology/oncology",
      references = "DiNardo et al., Blood 2015; Papaemmanuil et al., NEJM 2016"
    ),

    IDH2 = list(
      full_name = "Isocitrate Dehydrogenase 2",
      gene_function = "Mitochondrial enzyme that converts isocitrate to alpha-ketoglutarate. Mutations produce oncometabolite 2-hydroxyglutarate.",
      pathway = "Epigenetic regulation / Metabolism",
      chip_frequency = "Uncommon in CHIP (~1-2%)",
      risk_category = "High",
      clinical_significance = "Similar to IDH1. Gain-of-function mutations. FDA-approved targeted therapy available (enasidenib). R140 and R172 are distinct mutation sites with different clinical features.",
      hotspot_mutations = c("R140Q", "R172K", "R172S"),
      hotspot_note = "Mutations at R140 or R172 codons only",
      cardiovascular_risk = "Less well-characterized",
      progression_to_malignancy = "Higher progression risk; therapeutic target in AML",
      monitoring_recommendation = "Close monitoring; consider early referral to hematology/oncology",
      references = "Stein et al., Blood 2017; DiNardo et al., Blood 2015"
    ),

    EZH2 = list(
      full_name = "Enhancer of Zeste Homolog 2",
      gene_function = "Catalytic subunit of PRC2 complex. Methylates histone H3 at lysine 27 (H3K27me3) to repress gene transcription.",
      pathway = "Epigenetic regulation / Chromatin modification",
      chip_frequency = "Uncommon (~2-3%)",
      risk_category = "High",
      clinical_significance = "Loss-of-function mutations associated with myeloid malignancies. Can be gain-of-function in lymphoid malignancies. EZH2 inhibitors available for certain contexts.",
      hotspot_mutations = c("Various throughout gene"),
      hotspot_note = "Both loss-of-function (myeloid) and gain-of-function (lymphoid) mutations occur",
      cardiovascular_risk = "Not well-characterized in CHIP context",
      progression_to_malignancy = "Associated with MDS and AML development",
      monitoring_recommendation = "Regular hematologic monitoring; CBC every 6 months",
      references = "Ernst et al., Nature Genetics 2010; Nikoloski et al., Nature Genetics 2010"
    ),

    # ========================================================================
    # SPLICEOSOME GENES
    # ========================================================================
    SF3B1 = list(
      full_name = "Splicing Factor 3b Subunit 1",
      gene_function = "Core component of U2 small nuclear ribonucleoprotein (snRNP). Essential for pre-mRNA splicing at branch point recognition.",
      pathway = "RNA splicing / Spliceosome",
      chip_frequency = "Common in CHIP (~5-10%)",
      risk_category = "Intermediate-High",
      clinical_significance = "Mutations cause altered splice site selection. Strongly associated with ring sideroblasts in MDS. May confer better prognosis in MDS compared to other spliceosome mutations.",
      hotspot_mutations = c("K700E", "K666N", "H662Q"),
      hotspot_note = "Most mutations cluster in HEAT repeats; K700E most common",
      cardiovascular_risk = "Elevated cardiovascular risk",
      progression_to_malignancy = "Associated with MDS (especially MDS-RS); lower transformation risk to AML compared to other spliceosome genes",
      monitoring_recommendation = "Monitor for anemia; iron studies if ring sideroblasts suspected",
      references = "Papaemmanuil et al., NEJM 2011; Yoshida et al., Nature 2011"
    ),

    SRSF2 = list(
      full_name = "Serine and Arginine Rich Splicing Factor 2",
      gene_function = "Member of SR protein family. Regulates constitutive and alternative splicing by binding to exonic splicing enhancers.",
      pathway = "RNA splicing / Spliceosome",
      chip_frequency = "Common (~5%)",
      risk_category = "High",
      clinical_significance = "P95 hotspot mutations alter RNA binding specificity. Associated with worse prognosis in MDS. Often co-occurs with TET2 mutations.",
      hotspot_mutations = c("P95H", "P95L", "P95R"),
      hotspot_note = "Nearly all mutations at P95 codon",
      cardiovascular_risk = "Elevated cardiovascular risk",
      progression_to_malignancy = "Higher risk of progression; poor prognosis in MDS/CMML",
      monitoring_recommendation = "Close monitoring every 6 months; low threshold for bone marrow evaluation",
      references = "Yoshida et al., Nature 2011; Kim et al., Leukemia 2015"
    ),

    U2AF1 = list(
      full_name = "U2 Small Nuclear RNA Auxiliary Factor 1",
      gene_function = "Splicing factor that recognizes the 3' splice site AG dinucleotide. Required for spliceosome assembly.",
      pathway = "RNA splicing / Spliceosome",
      chip_frequency = "Less common (~2-3%)",
      risk_category = "High",
      clinical_significance = "Mutations at S34 and Q157 alter splice site recognition. Associated with MDS and secondary AML. Poor prognosis.",
      hotspot_mutations = c("S34F", "S34Y", "Q157P", "Q157R"),
      hotspot_note = "Two mutational hotspots: S34 and Q157",
      cardiovascular_risk = "Elevated risk",
      progression_to_malignancy = "High risk of progression to AML",
      monitoring_recommendation = "Close monitoring; consider early bone marrow evaluation",
      references = "Graubert et al., Nature Genetics 2012; Yoshida et al., Nature 2011"
    ),

    ZRSR2 = list(
      full_name = "Zinc Finger CCCH-Type, RNA Binding Motif and Serine/Arginine Rich 2",
      gene_function = "Component of U12-type spliceosome. Involved in recognition of 3' splice site in minor introns.",
      pathway = "RNA splicing / Minor spliceosome",
      chip_frequency = "Uncommon (~1-2%)",
      risk_category = "Intermediate-High",
      clinical_significance = "X-linked gene; mutations more common in males. Associated with MDS. Affects minor intron splicing.",
      hotspot_mutations = c("Truncating mutations throughout"),
      hotspot_note = "Loss-of-function mutations; no single hotspot",
      cardiovascular_risk = "Less well-characterized",
      progression_to_malignancy = "Associated with MDS development",
      monitoring_recommendation = "Regular CBC monitoring",
      references = "Madan et al., Nature Communications 2015"
    ),

    # ========================================================================
    # SIGNALING PATHWAY GENES
    # ========================================================================
    JAK2 = list(
      full_name = "Janus Kinase 2",
      gene_function = "Non-receptor tyrosine kinase essential for cytokine signaling. Mediates signaling from erythropoietin, thrombopoietin, and other receptors.",
      pathway = "JAK-STAT signaling",
      chip_frequency = "Common (~5-10%)",
      risk_category = "High",
      clinical_significance = "V617F mutation is pathognomonic for myeloproliferative neoplasms (MPNs). Even at low VAF, indicates risk of MPN development. JAK inhibitors (ruxolitinib) available.",
      hotspot_mutations = c("V617F"),
      hotspot_note = "V617F accounts for >95% of JAK2 mutations in CHIP/MPN",
      cardiovascular_risk = "Significantly elevated thrombotic risk; venous and arterial thrombosis",
      progression_to_malignancy = "High risk of MPN (PV, ET, PMF); may already indicate early MPN",
      monitoring_recommendation = "Urgent hematology referral; CBC with smear; consider bone marrow biopsy if VAF >2%",
      references = "James et al., Nature 2005; Levine et al., Cancer Cell 2005"
    ),

    MPL = list(
      full_name = "MPL Proto-Oncogene, Thrombopoietin Receptor",
      gene_function = "Thrombopoietin receptor that regulates megakaryocyte and platelet production.",
      pathway = "JAK-STAT signaling",
      chip_frequency = "Uncommon (~1%)",
      risk_category = "High",
      clinical_significance = "W515 mutations cause constitutive activation. Associated with essential thrombocythemia and primary myelofibrosis.",
      hotspot_mutations = c("W515L", "W515K", "W515R", "S505N"),
      hotspot_note = "W515 codon mutations most common",
      cardiovascular_risk = "Elevated thrombotic risk",
      progression_to_malignancy = "Indicates MPN or pre-MPN state",
      monitoring_recommendation = "Hematology referral; platelet count monitoring; bone marrow evaluation",
      references = "Pikman et al., PLoS Medicine 2006"
    ),

    CALR = list(
      full_name = "Calreticulin",
      gene_function = "Endoplasmic reticulum chaperone protein. Mutant forms activate MPL/JAK-STAT signaling.",
      pathway = "JAK-STAT signaling (via mutant-mediated MPL activation)",
      chip_frequency = "Uncommon (~1-2%)",
      risk_category = "High",
      clinical_significance = "Frameshift mutations in exon 9 create novel C-terminus that activates MPL. Type 1 (52bp del) and Type 2 (5bp ins) have different prognostic implications.",
      hotspot_mutations = c("L367fs (Type 1)", "K385fs (Type 2)"),
      hotspot_note = "All pathogenic mutations are frameshifts in exon 9",
      cardiovascular_risk = "Lower thrombotic risk compared to JAK2",
      progression_to_malignancy = "Indicates ET or PMF; Type 1 generally better prognosis",
      monitoring_recommendation = "Hematology referral; bone marrow evaluation",
      references = "Klampfl et al., NEJM 2013; Nangalia et al., NEJM 2013"
    ),

    NRAS = list(
      full_name = "NRAS Proto-Oncogene, GTPase",
      gene_function = "Small GTPase that activates RAS/MAPK signaling pathway. Regulates cell proliferation and survival.",
      pathway = "RAS/MAPK signaling",
      chip_frequency = "Less common (~1-2%)",
      risk_category = "High",
      clinical_significance = "Activating mutations at G12, G13, Q61. Associated with AML and CMML. May respond to MEK inhibitors.",
      hotspot_mutations = c("G12D", "G12V", "G13D", "Q61R", "Q61K"),
      hotspot_note = "Mutations only at codons 12, 13, or 61 are pathogenic",
      cardiovascular_risk = "Not well-characterized in CHIP context",
      progression_to_malignancy = "Higher risk; associated with AML",
      monitoring_recommendation = "Close hematologic monitoring",
      references = "Bacher et al., Blood 2010"
    ),

    KRAS = list(
      full_name = "KRAS Proto-Oncogene, GTPase",
      gene_function = "Small GTPase in RAS/MAPK pathway. Similar function to NRAS.",
      pathway = "RAS/MAPK signaling",
      chip_frequency = "Uncommon (~1%)",
      risk_category = "High",
      clinical_significance = "Similar to NRAS. Activating mutations at G12, G13, Q61.",
      hotspot_mutations = c("G12D", "G12V", "G12C", "G13D", "Q61H"),
      hotspot_note = "Same hotspots as NRAS",
      cardiovascular_risk = "Not well-characterized",
      progression_to_malignancy = "Associated with CMML and AML",
      monitoring_recommendation = "Close hematologic monitoring",
      references = "Bacher et al., Blood 2010"
    ),

    KIT = list(
      full_name = "KIT Proto-Oncogene, Receptor Tyrosine Kinase",
      gene_function = "Receptor tyrosine kinase for stem cell factor (SCF). Critical for hematopoietic stem cell maintenance.",
      pathway = "Receptor tyrosine kinase signaling",
      chip_frequency = "Rare in CHIP",
      risk_category = "High",
      clinical_significance = "D816V mutation associated with systemic mastocytosis and AML (especially core-binding factor AML). Imatinib-resistant; avapritinib approved.",
      hotspot_mutations = c("D816V", "D816Y", "N822K"),
      hotspot_note = "D816V most common in myeloid malignancies",
      cardiovascular_risk = "Not well-characterized",
      progression_to_malignancy = "Associated with mastocytosis and AML",
      monitoring_recommendation = "Hematology referral if detected",
      references = "Nagata et al., Nature 1995"
    ),

    FLT3 = list(
      full_name = "Fms Related Receptor Tyrosine Kinase 3",
      gene_function = "Receptor tyrosine kinase expressed on hematopoietic progenitors. Regulates proliferation and survival.",
      pathway = "Receptor tyrosine kinase signaling",
      chip_frequency = "Rare in CHIP (common in AML)",
      risk_category = "High",
      clinical_significance = "ITD (internal tandem duplication) and TKD (D835) mutations. More commonly seen in de novo AML. FLT3 inhibitors (midostaurin, gilteritinib) available.",
      hotspot_mutations = c("ITD", "D835Y", "D835V"),
      hotspot_note = "ITD in juxtamembrane domain; TKD at D835",
      cardiovascular_risk = "Not well-characterized",
      progression_to_malignancy = "High risk if detected; may indicate occult AML",
      monitoring_recommendation = "Urgent hematology referral; bone marrow evaluation",
      references = "Nakao et al., Leukemia 1996; Yamamoto et al., Blood 2001"
    ),

    CBL = list(
      full_name = "Cbl Proto-Oncogene",
      gene_function = "E3 ubiquitin ligase that negatively regulates receptor tyrosine kinase signaling by promoting receptor degradation.",
      pathway = "Ubiquitin-mediated signaling regulation",
      chip_frequency = "Uncommon (~1%)",
      risk_category = "Intermediate-High",
      clinical_significance = "Loss-of-function mutations lead to prolonged RTK signaling. Associated with CMML and juvenile myelomonocytic leukemia (JMML).",
      hotspot_mutations = c("Mutations in linker and RING finger domains"),
      hotspot_note = "Cluster in exons 8-9",
      cardiovascular_risk = "Not well-characterized",
      progression_to_malignancy = "Associated with CMML",
      monitoring_recommendation = "Monitor monocyte counts; hematology referral if monocytosis",
      references = "Sanada et al., Nature 2009"
    ),

    PTPN11 = list(
      full_name = "Protein Tyrosine Phosphatase Non-Receptor Type 11",
      gene_function = "Encodes SHP-2 phosphatase. Positive regulator of RAS signaling.",
      pathway = "RAS/MAPK signaling",
      chip_frequency = "Uncommon (~1%)",
      risk_category = "Intermediate-High",
      clinical_significance = "Gain-of-function mutations. Germline mutations cause Noonan syndrome. Somatic mutations associated with JMML and AML.",
      hotspot_mutations = c("E76K", "G60V", "A72V", "D61Y"),
      hotspot_note = "Mutations cluster in N-SH2 and PTP domains",
      cardiovascular_risk = "Associated with congenital heart disease (germline)",
      progression_to_malignancy = "Associated with JMML and AML",
      monitoring_recommendation = "Hematology referral",
      references = "Tartaglia et al., Nature Genetics 2001"
    ),

    GNAS = list(
      full_name = "GNAS Complex Locus",
      gene_function = "Encodes Gs-alpha subunit of heterotrimeric G proteins. Regulates cAMP signaling.",
      pathway = "G-protein coupled receptor signaling",
      chip_frequency = "Rare",
      risk_category = "Intermediate",
      clinical_significance = "Activating mutations associated with McCune-Albright syndrome (germline) and various tumors (somatic).",
      hotspot_mutations = c("R201C", "R201H", "Q227L"),
      hotspot_note = "Mutations impair GTPase activity",
      cardiovascular_risk = "Not well-characterized in CHIP",
      progression_to_malignancy = "Uncommon in myeloid malignancies",
      monitoring_recommendation = "Standard CHIP monitoring",
      references = "Weinstein et al., NEJM 1991"
    ),

    # ========================================================================
    # TRANSCRIPTION FACTORS
    # ========================================================================
    RUNX1 = list(
      full_name = "RUNX Family Transcription Factor 1",
      gene_function = "Master regulator of hematopoiesis. Required for definitive hematopoiesis and megakaryocyte/platelet development.",
      pathway = "Transcriptional regulation",
      chip_frequency = "Less common (~2-3%)",
      risk_category = "High",
      clinical_significance = "Germline mutations cause familial platelet disorder with predisposition to AML (FPD-AML). Somatic mutations associated with AML and MDS.",
      hotspot_mutations = c("Mutations throughout Runt domain and TAD"),
      hotspot_note = "Concentrated in Runt homology domain",
      cardiovascular_risk = "Not well-characterized",
      progression_to_malignancy = "High risk; associated with therapy-related AML and MDS",
      monitoring_recommendation = "Close monitoring; family history assessment for germline consideration",
      references = "Song et al., Nature Genetics 1999; Osato et al., Blood 1999"
    ),

    CEBPA = list(
      full_name = "CCAAT Enhancer Binding Protein Alpha",
      gene_function = "Transcription factor essential for granulocyte differentiation.",
      pathway = "Transcriptional regulation",
      chip_frequency = "Rare in CHIP",
      risk_category = "Intermediate-High",
      clinical_significance = "Biallelic mutations define a distinct AML subtype with favorable prognosis. N-terminal and C-terminal mutations have different effects.",
      hotspot_mutations = c("N-terminal frameshifts", "C-terminal in-frame insertions"),
      hotspot_note = "Biallelic mutations (N-term + C-term) most significant",
      cardiovascular_risk = "Not well-characterized",
      progression_to_malignancy = "Associated with AML; biallelic favorable prognosis",
      monitoring_recommendation = "Hematology referral if detected",
      references = "Pabst et al., Nature Genetics 2001"
    ),

    WT1 = list(
      full_name = "WT1 Transcription Factor",
      gene_function = "Zinc finger transcription factor. Regulates cell proliferation and differentiation.",
      pathway = "Transcriptional regulation",
      chip_frequency = "Uncommon",
      risk_category = "High",
      clinical_significance = "Germline mutations cause Wilms tumor predisposition. Somatic mutations associated with AML, often with poor prognosis.",
      hotspot_mutations = c("Exon 7 and exon 9 mutations"),
      hotspot_note = "Cluster in zinc finger domains",
      cardiovascular_risk = "Not well-characterized",
      progression_to_malignancy = "Associated with AML; adverse prognosis",
      monitoring_recommendation = "Hematology referral",
      references = "Paschka et al., JCO 2008"
    ),

    GATA1 = list(
      full_name = "GATA Binding Protein 1",
      gene_function = "Transcription factor essential for erythroid and megakaryocyte development.",
      pathway = "Transcriptional regulation",
      chip_frequency = "Rare",
      risk_category = "High",
      clinical_significance = "Germline mutations cause Diamond-Blackfan anemia and other cytopenias. Somatic mutations in Down syndrome-associated acute megakaryoblastic leukemia.",
      hotspot_mutations = c("Exon 2 mutations creating GATA1s"),
      hotspot_note = "N-terminal truncating mutations",
      cardiovascular_risk = "Not well-characterized",
      progression_to_malignancy = "Context-dependent; associated with DS-AMKL",
      monitoring_recommendation = "Specialist referral",
      references = "Wechsler et al., Nature Genetics 2002"
    ),

    GATA2 = list(
      full_name = "GATA Binding Protein 2",
      gene_function = "Transcription factor critical for hematopoietic stem cell maintenance and lymphoid development.",
      pathway = "Transcriptional regulation",
      chip_frequency = "Rare",
      risk_category = "High",
      clinical_significance = "Germline mutations cause GATA2 deficiency syndrome with immunodeficiency and MDS/AML predisposition. Important to distinguish germline from somatic.",
      hotspot_mutations = c("Zinc finger domain mutations"),
      hotspot_note = "ZF2 domain hotspot",
      cardiovascular_risk = "Not well-characterized",
      progression_to_malignancy = "High risk in germline carriers; MDS/AML predisposition",
      monitoring_recommendation = "Genetic counseling; family screening consideration",
      references = "Hahn et al., Nature Genetics 2011; Spinner et al., Blood 2014"
    ),

    NPM1 = list(
      full_name = "Nucleophosmin 1",
      gene_function = "Nucleolar phosphoprotein involved in ribosome biogenesis, centrosome duplication, and genomic stability.",
      pathway = "Nucleolar function / Ribosome biogenesis",
      chip_frequency = "Rare in CHIP (common in AML)",
      risk_category = "High",
      clinical_significance = "Exon 12 mutations cause cytoplasmic localization. Defines a distinct AML subtype with favorable prognosis (when not co-mutated with FLT3-ITD).",
      hotspot_mutations = c("W288fs (Type A)", "W288fs variants"),
      hotspot_note = "All mutations are 4bp insertions in exon 12",
      cardiovascular_risk = "Not well-characterized",
      progression_to_malignancy = "If detected, likely indicates AML or pre-AML state",
      monitoring_recommendation = "Urgent hematology referral; bone marrow evaluation",
      references = "Falini et al., NEJM 2005"
    ),

    # ========================================================================
    # TUMOR SUPPRESSORS
    # ========================================================================
    TP53 = list(
      full_name = "Tumor Protein P53",
      gene_function = "Master tumor suppressor. Regulates cell cycle arrest, DNA repair, and apoptosis in response to cellular stress.",
      pathway = "DNA damage response / Cell cycle regulation",
      chip_frequency = "Less common (~3-5%)",
      risk_category = "Very High",
      clinical_significance = "Highest risk CHIP gene. Associated with therapy-related MDS/AML, complex karyotype, and very poor prognosis. Often indicates prior cytotoxic therapy exposure.",
      hotspot_mutations = c("R175H", "R248Q", "R273H", "R282W", "Y220C"),
      hotspot_note = "Hotspots in DNA-binding domain",
      cardiovascular_risk = "Elevated; germline mutations associated with Li-Fraumeni syndrome",
      progression_to_malignancy = "Highest risk (~10-15% at 10 years); therapy-related myeloid neoplasms",
      monitoring_recommendation = "Urgent hematology referral; close monitoring; avoid unnecessary cytotoxic exposures",
      references = "Wong et al., Nature Medicine 2015; Coombs et al., Cell Stem Cell 2017"
    ),

    PPM1D = list(
      full_name = "Protein Phosphatase, Mg2+/Mn2+ Dependent 1D",
      gene_function = "Serine/threonine phosphatase that negatively regulates p53 and DNA damage response.",
      pathway = "DNA damage response",
      chip_frequency = "Uncommon (~2%)",
      risk_category = "High",
      clinical_significance = "Truncating mutations in exon 6 cause gain-of-function by removing regulatory domain. Strongly associated with prior cytotoxic therapy. Increases with age.",
      hotspot_mutations = c("Truncating mutations in exon 6"),
      hotspot_note = "All pathogenic mutations are truncating in exon 6",
      cardiovascular_risk = "Not well-characterized",
      progression_to_malignancy = "Associated with therapy-related MDS/AML",
      monitoring_recommendation = "Review prior therapy exposure; hematology referral",
      references = "Hsu et al., Cell 2018; Kahn et al., Nature Medicine 2018"
    ),

    PHF6 = list(
      full_name = "PHD Finger Protein 6",
      gene_function = "Chromatin-associated protein involved in transcriptional regulation. Tumor suppressor.",
      pathway = "Chromatin regulation",
      chip_frequency = "Rare",
      risk_category = "High",
      clinical_significance = "X-linked gene. More common in T-ALL than myeloid malignancies. Germline mutations cause Borjeson-Forssman-Lehmann syndrome.",
      hotspot_mutations = c("Truncating mutations throughout"),
      hotspot_note = "Loss-of-function mutations",
      cardiovascular_risk = "Not well-characterized",
      progression_to_malignancy = "Associated with T-ALL and AML",
      monitoring_recommendation = "Specialist referral",
      references = "Van Vlierberghe et al., Nature Genetics 2010"
    ),

    # ========================================================================
    # COHESIN COMPLEX
    # ========================================================================
    STAG2 = list(
      full_name = "Stromal Antigen 2",
      gene_function = "Core component of cohesin complex. Regulates sister chromatid cohesion and gene expression.",
      pathway = "Cohesin complex / Chromosome segregation",
      chip_frequency = "Uncommon (~2%)",
      risk_category = "Intermediate-High",
      clinical_significance = "X-linked gene; loss-of-function mutations. Associated with MDS and AML. May affect differentiation rather than proliferation.",
      hotspot_mutations = c("Truncating mutations throughout"),
      hotspot_note = "Loss-of-function mutations distributed across gene",
      cardiovascular_risk = "Not well-characterized",
      progression_to_malignancy = "Associated with MDS and AML",
      monitoring_recommendation = "Regular CBC monitoring",
      references = "Kon et al., Nature Genetics 2013"
    ),

    SMC1A = list(
      full_name = "Structural Maintenance of Chromosomes 1A",
      gene_function = "Core component of cohesin complex.",
      pathway = "Cohesin complex / Chromosome segregation",
      chip_frequency = "Rare",
      risk_category = "Intermediate",
      clinical_significance = "Less commonly mutated than STAG2 in myeloid malignancies. X-linked; germline mutations cause Cornelia de Lange syndrome.",
      hotspot_mutations = c("Various"),
      hotspot_note = "Missense mutations",
      cardiovascular_risk = "Not well-characterized",
      progression_to_malignancy = "Associated with MDS",
      monitoring_recommendation = "Standard CHIP monitoring",
      references = "Kon et al., Nature Genetics 2013"
    ),

    SMC3 = list(
      full_name = "Structural Maintenance of Chromosomes 3",
      gene_function = "Core component of cohesin complex.",
      pathway = "Cohesin complex / Chromosome segregation",
      chip_frequency = "Rare",
      risk_category = "Intermediate",
      clinical_significance = "Less commonly mutated than STAG2. Germline mutations cause Cornelia de Lange syndrome.",
      hotspot_mutations = c("Various"),
      hotspot_note = "Missense mutations",
      cardiovascular_risk = "Not well-characterized",
      progression_to_malignancy = "Associated with MDS",
      monitoring_recommendation = "Standard CHIP monitoring",
      references = "Kon et al., Nature Genetics 2013"
    ),

    RAD21 = list(
      full_name = "RAD21 Cohesin Complex Component",
      gene_function = "Kleisin subunit of cohesin complex. Required for sister chromatid cohesion.",
      pathway = "Cohesin complex / Chromosome segregation",
      chip_frequency = "Rare",
      risk_category = "Intermediate",
      clinical_significance = "Associated with MDS and AML. Less commonly mutated than STAG2.",
      hotspot_mutations = c("Various"),
      hotspot_note = "Truncating and missense mutations",
      cardiovascular_risk = "Not well-characterized",
      progression_to_malignancy = "Associated with MDS",
      monitoring_recommendation = "Standard CHIP monitoring",
      references = "Thota et al., Blood 2014"
    ),

    # ========================================================================
    # OTHER GENES
    # ========================================================================
    BCOR = list(
      full_name = "BCL6 Corepressor",
      gene_function = "Transcriptional corepressor that interacts with BCL6. Involved in PRC1-mediated gene silencing.",
      pathway = "Transcriptional repression",
      chip_frequency = "Uncommon (~2%)",
      risk_category = "Intermediate-High",
      clinical_significance = "Loss-of-function mutations. Associated with MDS and AML. X-linked; germline mutations cause oculofaciocardiodental syndrome.",
      hotspot_mutations = c("Truncating mutations throughout"),
      hotspot_note = "Loss-of-function",
      cardiovascular_risk = "Not well-characterized",
      progression_to_malignancy = "Associated with MDS and AML",
      monitoring_recommendation = "Regular CBC monitoring",
      references = "Grossmann et al., Blood 2011"
    ),

    SETBP1 = list(
      full_name = "SET Binding Protein 1",
      gene_function = "Binds SET protein; mutations cause protein stabilization and increased proliferation.",
      pathway = "Cell proliferation / PP2A regulation",
      chip_frequency = "Uncommon (~1%)",
      risk_category = "High",
      clinical_significance = "Hotspot mutations in SKI homologous region. Associated with atypical CML and CMML. Poor prognosis.",
      hotspot_mutations = c("D868N", "G870S", "I871T"),
      hotspot_note = "Cluster in SKI homologous domain (codons 858-871)",
      cardiovascular_risk = "Not well-characterized",
      progression_to_malignancy = "High risk; associated with aCML and CMML",
      monitoring_recommendation = "Hematology referral",
      references = "Piazza et al., Nature Genetics 2013"
    ),

    UTP23 = list(
      full_name = "UTP23 Small Subunit Processome Component",
      gene_function = "Component of small subunit processome involved in 18S rRNA maturation.",
      pathway = "Ribosome biogenesis",
      chip_frequency = "Rare",
      risk_category = "Low-Intermediate",
      clinical_significance = "Less well-characterized in CHIP context. May be a passenger mutation in some cases.",
      hotspot_mutations = c("Not well-defined"),
      hotspot_note = "Unknown",
      cardiovascular_risk = "Not characterized",
      progression_to_malignancy = "Unknown",
      monitoring_recommendation = "Standard CHIP monitoring",
      references = "Limited data"
    )
  )
}

#' Get clinical interpretation for a specific gene
#'
#' @param gene Gene symbol (e.g., "DNMT3A")
#' @return List containing gene information or NULL if not found
get_gene_info <- function(gene) {
  chip_info <- get_chip_gene_info()
  if (gene %in% names(chip_info)) {
    return(chip_info[[gene]])
  }
  return(NULL)
}

#' Get risk category for a gene
#'
#' @param gene Gene symbol
#' @return Character string: "Very High", "High", "Intermediate-High", "Intermediate", "Low-Intermediate", or "Unknown"
get_gene_risk <- function(gene) {
  info <- get_gene_info(gene)
  if (!is.null(info)) {
    return(info$risk_category)
  }
  return("Unknown")
}

#' Check if a protein change is a known hotspot mutation
#'
#' @param gene Gene symbol
#' @param protein_change Protein change (e.g., "V617F", "R882H")
#' @return Logical TRUE if hotspot
is_hotspot_mutation <- function(gene, protein_change) {
  info <- get_gene_info(gene)
  if (is.null(info) || is.null(info$hotspot_mutations)) {
    return(FALSE)
  }
  # Clean protein change
  clean_change <- gsub("^p\\.", "", protein_change)
  return(clean_change %in% info$hotspot_mutations)
}

#' Generate summary interpretation for a mutation
#'
#' @param gene Gene symbol
#' @param protein_change Protein change
#' @param vaf Variant allele frequency
#' @return Character string with clinical interpretation
generate_mutation_summary <- function(gene, protein_change, vaf) {
  info <- get_gene_info(gene)

  if (is.null(info)) {
    return(paste0("Mutation in ", gene, " detected at ", round(vaf * 100, 2),
                  "% VAF. Clinical significance unknown - gene not in standard CHIP panel."))
  }

  is_hotspot <- is_hotspot_mutation(gene, protein_change)
  risk <- info$risk_category

  summary_parts <- c()

  # Opening statement
  if (is_hotspot) {
    summary_parts <- c(summary_parts,
                       paste0("**Known pathogenic hotspot mutation** in ", info$full_name,
                              " (", gene, ") detected at ", round(vaf * 100, 2), "% VAF."))
  } else {
    summary_parts <- c(summary_parts,
                       paste0("Mutation in ", info$full_name, " (", gene,
                              ") detected at ", round(vaf * 100, 2), "% VAF."))
  }

  # Risk statement
  summary_parts <- c(summary_parts,
                     paste0("Risk category: **", risk, "**."))

  # Clinical significance
  summary_parts <- c(summary_parts, info$clinical_significance)

  # Monitoring recommendation
  summary_parts <- c(summary_parts,
                     paste0("Recommendation: ", info$monitoring_recommendation))

  return(paste(summary_parts, collapse = " "))
}

#' Get overall risk assessment based on mutation profile
#'
#' @param genes Vector of mutated gene symbols
#' @return List with overall_risk and explanation
assess_overall_risk <- function(genes) {
  risk_levels <- sapply(genes, get_gene_risk)

  # Priority: Very High > High > Intermediate-High > Intermediate > Low-Intermediate > Unknown
  if ("Very High" %in% risk_levels) {
    return(list(
      overall_risk = "Very High",
      explanation = "TP53 or equivalent high-risk mutation detected. Elevated risk of progression to myeloid malignancy. Close monitoring and hematology referral recommended."
    ))
  } else if ("High" %in% risk_levels) {
    return(list(
      overall_risk = "High",
      explanation = "High-risk gene mutation(s) detected. Increased risk of progression to myeloid malignancy and cardiovascular disease. Consider hematology referral and more frequent monitoring."
    ))
  } else if ("Intermediate-High" %in% risk_levels) {
    return(list(
      overall_risk = "Intermediate-High",
      explanation = "Intermediate-to-high risk mutation profile. Moderate elevation in risk of progression. Regular monitoring recommended."
    ))
  } else if ("Intermediate" %in% risk_levels) {
    return(list(
      overall_risk = "Intermediate",
      explanation = "Intermediate risk mutation profile (DNMT3A/TET2-like). Associated with modest increase in cardiovascular and malignancy risk. Annual monitoring generally sufficient."
    ))
  } else if (length(genes) > 0) {
    return(list(
      overall_risk = "Low-Intermediate",
      explanation = "Low-to-intermediate risk mutation profile. Standard CHIP monitoring recommended."
    ))
  } else {
    return(list(
      overall_risk = "None",
      explanation = "No CHIP mutations detected above reporting threshold."
    ))
  }
}
