import warnings
from pathlib import Path
from typing import List, Optional

import pandas as pd

from data_loading.utils import load_dataframe_from_file


def merge_transcriptom_data_to_raw_hospital(transcriptome_dataset: pd.DataFrame,
                                            raw_hospital_dataset: pd.DataFrame,
                                            filter_transcriptome_dataset_by_col: Optional[str] = "Transcriptom",
                                            transcriptome_dataset_patient_id_col_name: Optional[
                                                str] = "PID") -> pd.DataFrame:
    dataset = transcriptome_dataset
    if filter_transcriptome_dataset_by_col is not None:
        dataset = transcriptome_dataset[transcriptome_dataset[filter_transcriptome_dataset_by_col].fillna(False)]
    if "Unnamed: 0" in dataset.columns:
        dataset = dataset.drop(columns=["Unnamed: 0"])
    if transcriptome_dataset_patient_id_col_name is not None:
        dataset = dataset.set_index("PID")

    # add post treatment columns
    post_treatment_cols = [col for col in raw_hospital_dataset.columns if ".2" in col]
    # add fish_columns
    fish_cols = [col for col in raw_hospital_dataset.columns if
                 "t(" in col or "del(" in col or col in ['1q21+', 'IGH rearrangement',
                                                         'Cytogenetics Risk (0=standard risk, 1=single hit, 2=2+ hits)']]
    post_treatment_data = raw_hospital_dataset[raw_hospital_dataset["Time"] != "Post"][
        ["Code"] + post_treatment_cols + fish_cols].set_index("Code")
    dataset = dataset.merge(post_treatment_data, how="left", left_index=True, right_index=True, validate="one_to_one")

    # add data to TAL_3 patients
    TAL3_patients = [code for code in dataset.index if "P" == code[0]]
    dataset.loc[TAL3_patients, "Stage"] = 3
    dataset.loc[TAL3_patients, "Lenalidomide"] = 2

    return dataset


def generate_refracrotines_dataset(dataset: pd.DataFrame, treatment: str, non_ref_policy: str, feats: List[str]):
    pre_biopsy_ref_mask = dataset[treatment] == 2
    post_biopsy_ref_mask = dataset[f"{treatment}.2" if treatment != "DARA" else "Daratumumab.2"] == 2
    ref_mask = pre_biopsy_ref_mask | post_biopsy_ref_mask

    NDMM_STAGE_LIST = [1, 2]
    if non_ref_policy == "NDMM":
        newly_diagnosed = dataset["Stage"].apply(lambda x: x in NDMM_STAGE_LIST)
        non_ref_mask = newly_diagnosed
    elif non_ref_policy == "NDMM-POST_TREATMENT_REF":
        newly_diagnosed = dataset["Stage"].apply(lambda x: x in NDMM_STAGE_LIST)
        non_ref_mask = (newly_diagnosed) | (dataset[f"{treatment}.2"] == 4)
    elif non_ref_policy == "NON_EXPOSED":
        non_ref_mask = dataset[treatment].isna()
    else:
        raise NotImplementedError

    X = pd.concat([dataset[ref_mask][feats], dataset[non_ref_mask][feats]], axis=0)
    y = pd.concat([ref_mask[ref_mask].astype(int), non_ref_mask[non_ref_mask].astype(int) - 1])

    return X, y


DEFAULT_RESPONSE_CODING = {
    "response": "R",
    "no_response": "NR",
    "no_data": None
}
ALL_TREATMENTS = ["Bortezomib", "Ixazomib", "Carfilzomib", "Lenalidomide", "Thalidomide", "Pomalidomide",
                  "Cyclophosphamide", "Chemotherapy", "Venetoclax", "Dexamethasone", "Prednisone",
                  "Daratumumab", "Elotuzumab", "Belantamab", "Talquetamab", "Teclistamab", "Cevostamab",
                  "Selinexor", "Auto-SCT", "CART"]  


def load_and_process_clinical_data(clinical_data_path: Path, code_lower_case: bool, get_treatment_history: bool,
                                   get_hospital_stage: bool, get_post_treatment: bool, get_combination_exposure: bool,
                                   get_pfs_data: bool, get_fish_data: bool = False,
                                   additional_cols: Optional[List[str]] = None,
                                   treatment_names: Optional[List[str]] = None):
    clinical_data = load_dataframe_from_file(clinical_data_path)
    if code_lower_case:
        clinical_data['Code'] = clinical_data['Code'].str.lower()

    requested_cols = []
    requested_cols += ["Code", "Biopsy sequence No."]
    if treatment_names is None:
        treatment_names = ALL_TREATMENTS

    if get_hospital_stage:
        hospital_stage = 'Plasma cell dyscrasia at Bx time(0=NDMM, 1=RRMM, 2=SMM 3=MGUS,4=NDAL, 5=RRAL, 6=NDSPC, 7=MGRS, 8=None)'
        generated_hospital_stage = 'Disease Stage Hospital'
        if hospital_stage in clinical_data.columns:
            hospital_stage_map = {0: 'NDMM',
                                  1: "RRMM",
                                  2: "SMM",
                                  3: "MGUS",
                                  4: "AL", 5: "AL",
                                  6: "MGUS", 7: "MGUS",
                                  8: None}
        else:  # new version of clinical data
            hospital_stage = 'Plasma cell dyscrasia at Bx time(0=NDMM, 1=RRMM, 2=SMM 3=MGUS,4=NDAL, 5=RRAL, 6=NDSPC, 7=MGRS, 8=None, 10=Naïve_NDMM)'
            hospital_stage_map = {0: 'non_naive_NDMM',
                                  1: "RRMM",
                                  2: "SMM",
                                  3: "MGUS",
                                  4: "AL", 5: "AL",
                                  6: "MGUS", 7: "MGUS",
                                  8: 'Healthy',
                                  10: 'NDMM'}
        clinical_data[generated_hospital_stage] = clinical_data[hospital_stage].map(hospital_stage_map)

        requested_cols += [generated_hospital_stage]
    if get_treatment_history:
        requested_cols += treatment_names
    if get_post_treatment:
        requested_cols += [f"{treatment}.2" for treatment in treatment_names]
    if get_combination_exposure:
        requested_cols += [
            "Triple Ref.", "Triple Exp.",
            "Penta Ref.", "Penta Exp.",
            "Belantamab", "Belantamab Ref"
        ]
    if get_pfs_data:
        requested_cols += [' PFS (Months)', 'PFS Event', 'Biopsy date', 'Progression date', 'Last FU Date']

    if get_fish_data:
        fish_cols = [col for col in clinical_data.columns
                     if "t(" in col or "del(" in col or col in ['1q21+', 'IGH rearrangement',
                                                                'Cytogenetics Risk (0=standard risk, 1=single hit, 2=2+ hits)']]
        requested_cols += fish_cols
    if additional_cols is not None:
        for col in additional_cols:
            if col not in clinical_data.columns:
                raise ValueError(f"{col} is requested but not in clinical data")
            else:
                requested_cols.append(col)
    requested_clinical_data = clinical_data[requested_cols]
    return requested_clinical_data


def add_response_columns_to_drug_combination(dataset: pd.DataFrame, combination: str,
                                             response_policy: str, no_response_policy: str,
                                             coding: Optional[dict] = None,
                                             raise_errors=True) -> pd.DataFrame:
    df = dataset.copy()
    if coding is None:
        coding = DEFAULT_RESPONSE_CODING

    no_response_mask = pd.Series([False] * len(df), index=df.index)
    if 'refractory' in no_response_policy:
        refractory_mask = df[f"{combination} Ref."] == 1
        no_response_mask = no_response_mask | refractory_mask

    response_mask = pd.Series([False] * len(df), index=df.index)
    if 'exposed_non_refractory' in response_policy:
        exposed_non_refractory_mask = (df[f"{combination} Ref."] != 1) & (df[f"{combination} Exp."] == 1)
        response_mask = response_mask | exposed_non_refractory_mask
    if 'NDMM_SMM' in response_policy:
        # NDMM_mask = df['Disease Stage (MGUS=0, SMM=1, NDMM=2, RRMM=3, None=4)'].apply(lambda x: x in [1, 2])
        NDMM_mask = df['Disease'].apply(lambda x: x in ['NDMM', 'SMM'])
        response_mask = response_mask | NDMM_mask

    if raise_errors:
        assert sum(no_response_mask) > 0, f"no patients follow no response policy for combination: {combination}"
        assert sum(response_mask) > 0, f"no patients follow response policy for combination: {combination}"
        assert sum(response_mask & no_response_mask) == 0, \
            f"some patients follow both response and no response policies for combination: {combination}"

    response_series = pd.Series([coding["no_data"]] * len(df), index=df.index)
    response_series[response_mask] = coding["response"]
    response_series[no_response_mask] = coding["no_response"]

    df[f"{combination}_response"] = response_series

    return df


def add_response_columns_to_specific_treatment(dataset: pd.DataFrame, treatment: str,
                                               response_policy: str, no_response_policy: str,
                                               coding: Optional[dict] = None,
                                               raise_errors=True) -> pd.DataFrame:
    df = dataset.copy()
    if coding is None:
        coding = DEFAULT_RESPONSE_CODING

    # no_response_policy like 'pre_exposed|post_refractory'
    no_response_mask = pd.Series([False] * len(df), index=df.index)
    if 'pre_exposed' in no_response_policy:
        pre_exposed_mask = df[treatment] == 2 | df[treatment] == 1
        no_response_mask = no_response_mask | pre_exposed_mask
    elif 'pre_refractory' in no_response_policy:
        pre_refractory_mask = df[treatment] == 2
        no_response_mask = no_response_mask | pre_refractory_mask

    if 'last_line_exposed' in no_response_policy:
        last_line_exposed_mask = df[f"{treatment}.1"] == 2 | df[f"{treatment}.1"] == 1
        no_response_mask = no_response_mask | last_line_exposed_mask
    elif 'last_line_refractory' in no_response_policy:
        last_line_refractory_mask = df[f"{treatment}.1"] == 2
        no_response_mask = no_response_mask | last_line_refractory_mask

    if 'post_refractory' in no_response_policy:
        post_refractory_mask = df[f"{treatment}.2"] == 2
        no_response_mask = no_response_mask | post_refractory_mask

    # response_policy like 'NDMM_SMM|post_sensitive|not_exposed'
    response_mask = pd.Series([False] * len(df), index=df.index)
    if 'NDMM_SMM' in response_policy:
        # NDMM_mask = df['Disease Stage (MGUS=0, SMM=1, NDMM=2, RRMM=3, None=4)'].apply(lambda x: x in [1, 2])
        NDMM_mask = df['Disease'].apply(lambda x: x in ['NDMM', 'SMM'])
        response_mask = response_mask | NDMM_mask
    if 'post_sensitive' in response_policy:
        post_sensitive = dataset[f"{treatment}.2"] == 4
        treated_alone_mask = (~ dataset[dataset.columns[dataset.columns.str.contains("\.2")]].isna()).sum(axis=1) == 1
        response_mask = response_mask | (post_sensitive & treated_alone_mask)
    if 'not_exposed' in response_policy:
        not_exposed_mask = df[treatment].isna()
        response_mask = response_mask | not_exposed_mask

    if raise_errors:
        assert sum(no_response_mask) > 0, f"no patients follow no response policy for treatment: {treatment}"
        assert sum(response_mask) > 0, f"no patients follow response policy for treatment: {treatment}"

    response_series = pd.Series([coding["no_data"]] * len(df), index=df.index)
    response_series[response_mask] = coding["response"]
    response_series[no_response_mask] = coding["no_response"]

    df[f"{treatment}_response"] = response_series
    if sum(response_mask & no_response_mask) != 0:
        cols_to_print = ['Hospital.Code', 'Biopsy.Sequence', 'CD45', 'PC', 'Disease', 'Project', 'Cohort', 'Method',
                         f'{treatment}', f'{treatment}.2', f"{treatment}_response",
                         'Plasma cell dyscrasia at Bx time(0=NDMM, 1=RRMM, 2=SMM 3=MGUS,4=NDAL, 5=RRAL, 6=NDSPC, 7=MGRS, 8=None)']
        cols_to_print = list(set(cols_to_print).intersection(df.columns))
        warn = f"some patients follow both response and no response policies for treatment: {treatment}\n" \
               f"will consider them as non responders\n\n" \
               f"{df[response_mask & no_response_mask][cols_to_print]} "
        warnings.warn(warn)

    return df


def add_Kydar_response(df_with_clinical_data: pd.DataFrame, number_of_months: int = 4, coding: Optional[dict] = None):
    if coding is None:
        coding = DEFAULT_RESPONSE_CODING
    kydar_mask = df_with_clinical_data["Project"] == "Kydar"
    enough_pfs_mask = df_with_clinical_data[" PFS (Months)"] >= number_of_months

    kydar_response = 'Kydar_response'
    df_with_clinical_data[kydar_response] = coding['no_data']
    df_with_clinical_data[kydar_response][kydar_mask & enough_pfs_mask] = coding['response']
    df_with_clinical_data[kydar_response][kydar_mask & (~enough_pfs_mask)] = coding['no_response']
    return df_with_clinical_data


def add_CART_response(patient_df: pd.DataFrame, full_clinical_df: pd.DataFrame, pfs_policy="9M PFS",
                      coding: Optional[dict] = None):
    if coding is None:
        coding = DEFAULT_RESPONSE_CODING

    if pfs_policy not in ("9M PFS", "6M PFS", '3M PFS'):
        raise ValueError()
    df = patient_df.copy()

    CART_response_col = 'CART_response'
    full_clinical_df["Hospital.Code"] = full_clinical_df["Pt No."].apply(
        lambda x: f"CART_P{x}" if len(str(x)) == 2 else f"CART_P0{x}").str.lower()
    full_clinical_df[CART_response_col] = full_clinical_df[pfs_policy]
    full_clinical_df[CART_response_col] = full_clinical_df[CART_response_col].map(
        {1: coding['response'], 0: coding['no_response']})

    if CART_response_col in df.columns:
        df = df.drop(columns=CART_response_col)
    df = df.merge(full_clinical_df[["Hospital.Code", CART_response_col]], how='left', on="Hospital.Code")
    return df


def add_general_response(patient_df: pd.DataFrame, add_best_response: bool = True, add_pfs: bool = True,
                         pfs_thresh_months=9, coding: Optional[dict] = None):
    patient_df = patient_df.copy()
    _DAYS_A_MONTH = 30.4

    if coding is None:
        coding = DEFAULT_RESPONSE_CODING
    if not (add_best_response or add_pfs):
        raise ValueError("must generate at least one of the respnses: best_response and/or pfs")
    if add_best_response:
        general_response_col = 'general_response'
        post_treatment_cols = [f'{treatment}.2' for treatment in ALL_TREATMENTS]
        patient_df[general_response_col] = (patient_df[post_treatment_cols] == 4).any(axis=1).map(
            {True: 'response', False: 'no_response'})
        patients_with_post_data = patient_df[post_treatment_cols].any(axis=1)
        patient_df[general_response_col][~patients_with_post_data] = 'no_data'

        patient_df[general_response_col] = patient_df[general_response_col].map(coding)
    if add_pfs:
        general_pfs_response_col = 'general_pfs_response'
        progression_date_map = {"25/02//2019": "25/02/2019", '13/082018': '13/08/2018'}  # manual typos correction
        patient_df['Progression date'] = patient_df['Progression date'].apply(
            lambda x: progression_date_map[x] if x in progression_date_map else x)

        time_from_biopsy_to_FU = pd.to_datetime(patient_df['Last FU Date']) - pd.to_datetime(patient_df['Biopsy date'])
        time_from_biopsy_to_FU_months = time_from_biopsy_to_FU.apply(lambda x: x.days / _DAYS_A_MONTH)

        pfs_higher_then_thresh = patient_df[' PFS (Months)'] > pfs_thresh_months
        patient_df[general_pfs_response_col] = pfs_higher_then_thresh.map(
            {True: 'response', False: 'no_response'})

        not_reliable_pfs = (~pfs_higher_then_thresh) & (time_from_biopsy_to_FU_months < pfs_thresh_months)
        patient_df[general_pfs_response_col][not_reliable_pfs] = 'no_data'
        patient_df[general_pfs_response_col][patient_df[' PFS (Months)'].isna()] = 'no_data'

        patient_df[general_pfs_response_col] = patient_df[general_pfs_response_col].map(coding)

    return patient_df

def add_general_response_cens(patient_df: pd.DataFrame, pfs_thresh_months=9, coding: Optional[dict] = None):
    """
    Determine response based on PFS data and follow-up time.
    
    Args:
        patient_df: DataFrame containing:
            - 'PFS Event': 1 for progression/death, 0 for censoring
            - 'PFS (Months)': Time to progression/censoring
            - 'Last FU Date': Last follow-up date
            - 'Biopsy date': Initial biopsy date
    """
    if coding is None:
        coding = DEFAULT_RESPONSE_CODING
        
    general_pfs_response_col = f'general_pfs_response_cens_{str(pfs_thresh_months)}M'
    
    time_from_biopsy_to_FU = pd.to_datetime(patient_df['Last FU Date']) - pd.to_datetime(patient_df['Biopsy date'])
    time_from_biopsy_to_FU_months = time_from_biopsy_to_FU.apply(lambda x: x.days / 30.4)
    
    patient_df[general_pfs_response_col] = 'no_data'
    
    # get masks
    has_pfs_data = ~patient_df[' PFS (Months)'].isna()
    is_progression = patient_df['PFS Event'] == 1
    is_censored = patient_df['PFS Event'] == 0
    pfs_higher_than_thresh = patient_df[' PFS (Months)'] > pfs_thresh_months
    adequate_followup = time_from_biopsy_to_FU_months >= pfs_thresh_months
    
    # clear pesponses
    # - progressed after threshold
    # - censored after threshold with adequate follow-up
    response_mask = (
        has_pfs_data & 
        pfs_higher_than_thresh & 
        (is_progression | (is_censored & adequate_followup))
    )
    patient_df.loc[response_mask, general_pfs_response_col] = 'response'
    
    # clear non-responses
    # - who definitely progressed before threshold
    no_response_mask = has_pfs_data & is_progression & ~pfs_higher_than_thresh
    patient_df.loc[no_response_mask, general_pfs_response_col] = 'no_response'
    
    # no data
    # - missing PFS data
    # - censored before threshold
    # - censored after threshold but without adequate follow-up
    no_data_mask = (
        ~has_pfs_data | 
        (is_censored & ~pfs_higher_than_thresh)
    )
    patient_df.loc[no_data_mask, general_pfs_response_col] = 'no_data'
    
    patient_df[general_pfs_response_col] = patient_df[general_pfs_response_col].map(coding)
    print(patient_df[general_pfs_response_col].value_counts(dropna = False))
    return patient_df

def generate_gen_info(row, columns):
    info = []
    for col in columns:
        if row[col] == 'Yes':
            info.append(f'{col}: Yes')
        elif row[col] == 'No':
            info.append(f'{col}: No')
    return ', '.join(info)
