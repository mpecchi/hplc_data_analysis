# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 08:57:55 2022

@author: mp933
"""
import numpy as np

import pandas as pd
import pubchempy as pcp
import operator as oper
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib as plib

# =============================================================================
# # graphical paramters    
# =============================================================================
sns.set_style("whitegrid")
clrs = sns.color_palette('deep', 10)   # list with colors
# list with linestyles for plotting
lnstls = [(0, ()),  # solid
          (0, (1, 1)),  # 'densely dotted'
          (0, (5, 1)),  # 'densely dashed'
          (0, (3, 1, 1, 1)),  # 'densely dashdotted'
          (0, (3, 1, 1, 1, 1, 1)),  # 'densely dashdotdotted'
          (0, (5, 5)),  # 'dashed'
          (0, (3, 5, 1, 5)),  # 'dashdotted'
          (0, (1, 5)),  # dotted
          (0, (3, 5, 1, 5, 1, 5)),  # 'dashdotdotted'
          (0, (1, 10)),  # 'loosely dotted'
          (0, (5, 10)),  # 'loosely dashed'
          (0, (3, 10, 1, 10)),  # 'loosely dashdotted'
          (0, (3, 10, 1, 10, 1, 10))]  # 'loosely dashdotdotted'
# list with letters for plotting
lttrs = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm']
# list with markers for plotting
mrkrs = ["o", "v", "X", "s", "p", "^", "P", "<", ">", "*", "d"]
# list with hatches for plotting
htchs = (None, '//', '...', '--', 'O', '\\\\', 'oo', '\\\\\\', '/////', '.....')


def PathsCreate(subfolder=''):
    ''' This function creates 2 folder paths (independently of the os). \
        For each path checks the folder existence, if False it creates it.
    in_path : pathlib object
        path to the Input folder.
    out_path : pathlib object
        path to the Output folder.''' 
    try:  # start from script folder
        script_path = plib.Path(__file__).parents[0]  # script folder path
    except NameError:  # if a cell is run alone and the __file__ is not available
        script_path = plib.Path.cwd()  # the current working directory
    # create all necessary paths to subfolders starting from the script one
    in_path = plib.Path(script_path, "Input", subfolder)
    out_path = plib.Path(script_path, 'Output', subfolder)
    # check existence of each folders and create them if missing
    plib.Path(in_path).mkdir(parents=True, exist_ok=True)
    plib.Path(out_path).mkdir(parents=True, exist_ok=True)
    return in_path, out_path  # returns the two object-paths




def PrintLog(out_path, text):
    # used to both print info on screen and save them on a report file in the 
    # provided path, prints are appended on the old file
    report_path = plib.Path(out_path, 'LogFile.txt')
    print(text)
    with open(report_path, "a") as f:
        print(text, file=f)  


def CheckSubstCAS(cas, CASsubstDB=None):
    # the CASsubsDB contains misudentified CAS to substitute and known contaminants
    # that are substitute with 0-0-0 so they get ignored
    if not isinstance(cas, str):  # if cas is not a str, so its nan
        cas = '0-0-0'
        return cas
    if cas == '0-0-0':  # if CAS is unidentified nothing we can do
        pass
    else:  # check that cas are not mispelled
        # CAS = xxxxxxx-yy-z,  xxxxx (2 to 7 digits, y 2 digidts, z one)
        parts = cas.split('-')
        if len(parts[1]) == 1:  # a zero is missing on the left digit
            parts[1] = '0' + parts[1]
        cas =  '-'.join(parts)
        if not CASsubstDB is None:
            wrongCASs = CASsubstDB['WrongCAS'].tolist()
            if cas in wrongCASs:
                cas = CASsubstDB.loc[CASsubstDB['WrongCAS'] == 
                                          cas, 'RightCAS'].iloc[0]
    return cas


def CAStoDB(cas2search, name2search, comp_db, ElemenFunGroups,
            in_path=None, out_path=None, filename='CompDB'):
   
    """ this function does...
    """
    elements = ElemenFunGroups['elements']['elements'].tolist()
    el_masses = ElemenFunGroups['elements']['mass_elements'].tolist()
    functional_groups = ElemenFunGroups['functional_groups'
                                        ]['functional_groups'].tolist()
    
    if comp_db is None: # if No DB is provided, start from scratch
        columns = ['CAS', 'iupac_name', 'molecular_formula', 'canonical_smiles', 
                   'molecular_weight'] \
            + elements + ['mp_'+el for el in elements] + \
                ['fg_'+fr for fr in functional_groups] + ['GC_name', 'str-CAS']
        comp_db = pd.DataFrame(columns=columns).set_index('CAS')
    if CheckSubstCAS(cas2search) in comp_db.index.to_list():
        return comp_db  # since the cas is already there
    else:
        cas = CheckSubstCAS(cas2search)
        try:  # try to retrieve compound info using cas (uses pubchempy)
            comp = pcp.get_compounds(cas, 'name')[0]    
            CasFound = True            
        except IndexError:  # if cas does not work, uses the names
            try: 
                comp = pcp.get_compounds(name2search, 'name')[0]
                CasFound = True
            except IndexError:  # if name does not work, uses the names
                CasFound = False  # CAS AND NAME DID NOT WORK
        # once the CAS search ended (positively or negatively)
        df = pd.DataFrame()  #create a df to store the info  
        df.loc[0, 'CAS'] = cas  # add to the df
        if CasFound:  # data will populate the df, otherwise only nan  
            if comp.iupac_name:  # if it is available
                df.loc[0, 'iupac_name'] = comp.iupac_name
            else:
                df.loc[0, 'iupac_name'] = name2search
            
            df.loc[0, 'molecular_weight'] = float(comp.molecular_weight)
            df.loc[0, 'molecular_formula'] = comp.molecular_formula
            df.loc[0, 'canonical_smiles'] = comp.canonical_smiles        
            for e, el in enumerate(elements):  # count all atoms presence            
                df.loc[0, el] = comp.to_dict()['elements'].count(el)
                if df.loc[0, el] > 0:
                    df.loc[0, 'mp_' + el] = (df.loc[0, el]*el_masses[e]/
                                             df.loc[0, 'molecular_weight'])
                else:
                    df.loc[0, 'mp_' + el] = 0 
          
            # FUNCTIONAL GROUPS IDENTIFIER
            if functional_groups is not None:                           
                fun_gr_values = ['['+str(fr)+']' for fr in functional_groups]
                fun_gr_labels = ['fg_' + fr for fr in functional_groups]
    
                fragmentation_scheme = dict(zip(fun_gr_labels, fun_gr_values))
    
                frg = fragmenter(fragmentation_scheme, 
                                  fragmentation_scheme_order=fun_gr_labels, 
                                  algorithm='simple')
                # df = pd.DataFrame(index=smiles, columns=frg_sch_ord)
                smi = comp.canonical_smiles
                fragmentation, success, fragmentation_matches = frg.fragment(smi)
                for key in fragmentation_scheme.keys():
                    try:
                        df.loc[0, key] = fragmentation[key]
                    except KeyError:
                        df.loc[0, key] = 0
                                        
        else:  # save a list with only the cas and name of non-found cas
            report_path = plib.Path(out_path, 'UnfoundCAS.txt')
            with open(report_path, "a") as f:
                print(cas + '\t' + name2search, file=f)  
            PrintLog(out_path, cas + ' not found and added to the UnfoundCAS')
                    
        df.loc[0, 'str-CAS'] = "'" + cas  # to avoid excel conversion to date
        df.loc[0, 'GC_name'] = name2search  # the name from GC file
        df.set_index('CAS', inplace=True)  # set the cas number as the index 
        # create and save the update comp_db and save the old one for safety
        UpdatedDB = pd.concat([comp_db, df], axis=0)     
        db_path = plib.Path(PathsCreate()[0], filename + '.xlsx')
        UpdatedDB.to_excel(db_path)   
        PrintLog(out_path, cas + ' added to CompDB')
        return UpdatedDB


def ClassificToDB(cas, classif_db, ClassifInstr, CompDB,
                  in_path=None, out_path=None, filename='ClassifDB'):
        
    clsf_keys = ClassifInstr.keys()
    
    if classif_db is None: 
        classif_db = pd.DataFrame(columns=['CAS', 'iupac_name'] + 
                                  [k for k in clsf_keys] + 
                                  ['str-CAS']).set_index('CAS')
    if cas in classif_db.index.to_list():
        return classif_db  # since the cas is already there
    df = pd.DataFrame()
    df.loc[0, 'CAS'] = cas  # add to the df
    
    df.loc[0, 'str-CAS'] = "'" + cas  # to avoid excel conversion to d
    df.loc[0, 'iupac_name'] = CompDB.loc[cas, 'iupac_name']
    for key in clsf_keys:
        clsfc = ClassifInstr[key] 
        df.loc[0, key] = 'n.a.'
        if key.split('-')[1] == 'count':
            for cl in list(clsfc):  # classifications
                if df.loc[0, key] != 'n.a.':
                    break  # if a class was found
                if CompDB.loc[cas, clsfc[cl].dropna().tolist()].sum() > 0:  
                    df.loc[0, key] = cl
                    # prefixes given for each classification
        elif key.split('-')[1] == 'str':
            if CompDB.loc[cas, 'iupac_name']:  # must be there to classify
                for cl in list(clsfc):  # classifications
                    if df.loc[0, key] != 'n.a.':
                        break  # if a class was found
                    # prefixes given for each classification
                    for prefix in clsfc[cl].dropna():
                        try:
                            if prefix in CompDB.loc[cas, 'iupac_name']:  # then set the class
                                df.loc[0, key] = cl
                        except TypeError:
                            df.loc[0, key] = 'n.a.'
    df.set_index('CAS', inplace=True)  # set the cas number as the index
    # create and save the update DB and save the old one for safety
    UpdatedDB = pd.concat([classif_db, df], axis=0)   
    
    UpdatedDB.to_excel(plib.Path(plib.Path(in_path, filename + '.xlsx')))   
    PrintLog(out_path, cas + ' classified')
    return UpdatedDB
 

def GCMS_file_analysis(in_path, out_path, file, ElemenFunGroups, CompDB, 
                       ClassifInstr, ClassifDB, CASsubstDB, RunsInfo, CalDB,
                       RRFsDB, cal_method, 
                       cols_cal_area=['Area 1', 'Area 2', 'Area 3', 'Area 4', 'Area 5'],
                       cols_cal_ppms=['PPM 1', 'PPM 2', 'PPM 3', 'PPM 4', 'PPM 5']
                       ):
    # =========================================================================
    # general calculations for the same sample (valid for all CAS in it)
    # =========================================================================
    PrintLog(out_path, '\n'+ file + " sample started.")     
    # import GC output for sample, it contains unidentified compounds
    sample_000 = pd.read_csv(plib.Path(in_path, file + '.txt'),
                             delimiter='\t', index_col=0, skiprows=8)
    # rename columns to avoid confusion
    sample_000.rename(columns={'CAS #': 'CAS', 'Name':'GC_name'}, inplace=True)
    # correct all CAS that are missing the zero in xxxxx-0y-z and add it
    # and change wrong ones with right ones
    sample_000['CAS'] = sample_000['CAS'].apply(CheckSubstCAS,
                                                CASsubstDB=CASsubstDB)
    # drop all CAS == '0-0-0'
    sample = sample_000[sample_000['CAS'] != '0-0-0'].copy() 
    # areas attributed to the same CAS are summed 
    sample_sum_areas = sample.groupby('CAS')['Area'].sum()
    # the first ret time is kept for each duplicated CAS  
    sample.drop_duplicates(subset='CAS', keep='first', inplace=True)
    # if Excel convert CAS entry into date format, str-CAS is unaffected
    sample['str-CAS'] = "'" + sample['CAS']  # safe from Excel data autoformat
    sample.set_index('CAS', inplace=True)  # set the cas as the index
    sample['Area'] = sample_sum_areas  # used summed areas as areas 
    RunInfo = RunsInfo.loc[file, :]  # keep only the file row
    if cal_method == 'ExtStd':  # if external calibration is available
        conc_sample_mg_L = (RunInfo['mass_sample_mg']/
                            RunInfo['vol_Solv_mL']*1000)
    elif cal_method == 'IntStd':       
        # cas of internal standards 1 and 2
        cas_IS1, cas_IS2 = RunInfo['CAS_IS1'], RunInfo['CAS_IS2']
        # ensures the IS1 and IS2 information is in the CompDB and ClassifDB
        for c in [cas_IS1, cas_IS2]:
            CompDB = CAStoDB(c, None, CompDB, ElemenFunGroups,
                             in_path, out_path)
            ClassifDB = ClassificToDB(c, ClassifDB, ClassifInstr, CompDB, 
                                      in_path, out_path)
        # areas of internal standards 1 and 2 
        area_IS1 = sample.loc[cas_IS1, 'Area']
        area_IS2 = sample.loc[cas_IS2, 'Area']
        # ideal concentration of IS1 assuming no loss
        conc_IS1_id = (RunInfo['conc_IS1_ug_mL']*RunInfo['V_IS1_initial_mL']/
                    RunInfo['V_solvent_mL'])
        # real conc of IS1 from its RRF with IS2
        conc_IS1 = (RRFsDB.loc[cas_IS1, 'RRF_'+ cas_IS2] *
                    RunInfo['conc_IS2_ug_mL'] * area_IS1/area_IS2)  #ug/mL
        # conc correction factor
        CCF = conc_IS1_id*RunInfo['aliquot_factor']/conc_IS1
        if CCF < 0.95:  # it shoudl be > 1
            print('WARNING: CCF = ', round(CCF, 3))
    # =========================================================================
    # specific calculations for the each CAS in the sample
    # =========================================================================
    for cas in sample.index: # iterate over all peaks detected by the GC sample
        # check for the presence of the CAS, add it to the CompDb if not present
        CompDB = CAStoDB(cas, sample.loc[cas, 'GC_name'], CompDB, 
                         ElemenFunGroups, in_path, out_path)
        CompInfo = CompDB.loc[cas, :]
        MW_cas = CompInfo['molecular_weight']
        if cal_method == 'ExtStd':
            # check if there is a calibration for the specific compound
            if cas in CalDB.index.tolist():  
                cas_clbr = cas  # cas used for the calibration
                CalibrComp = 'Self'  # type of calibration used for this cas
            else:  # semi-calibration with the nearest MW
                # get the nearest MW from the overall list of compounds
                cas_clbr = CalDB.index[np.abs(CalDB['MW'] - MW_cas).argmin()]
                # note the type of calibration and the compound used for it
                CalibrComp = cas_clbr + ' (MW=' + \
                    str(round(CalDB.loc[cas_clbr, 'MW'], 0)) + ')'                            
            # areas and ppms for the calibration are taken from df_clbr
            cal_areas = CalDB.loc[cas_clbr, cols_cal_area].to_numpy(dtype=float)
            cal_ppms = CalDB.loc[cas_clbr, cols_cal_ppms].to_numpy(dtype=float)
            # linear fit of calibration curve (exclude nan), get ppm from area 
            fit = np.polyfit(cal_areas[~np.isnan(cal_areas)], 
                             cal_ppms[~np.isnan(cal_ppms)], 1)
            # concentration of CAS at the injection solution (GC vial)
            conc_inj_mg_L = np.poly1d(fit)(sample.loc[cas, 'Area'])      
            if conc_inj_mg_L < 0:  # conc cannot be < 0 
                conc_inj_mg_L = 0
            # fraction of cas over sample (using conc in GC vial) [mg/mg]
            f_sample = conc_inj_mg_L/(conc_sample_mg_L)  # [mg/mg]
            # yield refers to the original feedstock to produce the sample
            yield_m = f_sample*RunInfo['yield_sample_f']  # [mg/mg]
        elif cal_method == 'IntStd':  # Internal standards calibration
            if cas in RRFsDB.index.tolist():  # if RRF is available for cas
                cas_rrf = cas  # cas used for the calibration
                CalibrComp = 'Self'  # type of calibration used for this cas
            else:  # calibration using RRF from the nearest MW
                cas_rrf = RRFsDB.index[np.abs(RRFsDB['MW']-MW_cas).argmin()]
                # note the type of calibration and the compound used for it
                CalibrComp = cas_rrf + ' (MW=' + \
                    str(round(RRFsDB.loc[cas_rrf, 'MW'], 0)) + ')'
            # uncorrected concentration of cas in the injection vial rel to IS2
            unc_conc_inj_mg_L = (RRFsDB.loc[cas_rrf, 'RRF_' + cas_IS2]*
                                 RunInfo['conc_IS2_ug_mL']
                                 *sample.loc[cas,'Area']/area_IS2)
            # corrected concentration using the CCF [mg/L]
            conc_inj_mg_L = unc_conc_inj_mg_L*CCF/RunInfo['aliquot_factor']  
            # fraction of cas in the analyzed sample
            f_sample = (conc_inj_mg_L*RunInfo['V_solvent_mL']
                        /RunInfo['mass_sample_mg']/1000)  # [ug/ug]
            # here refers to mass of ceramic
            yield_m = (conc_inj_mg_L*RunInfo['V_solvent_mL']
                          /RunInfo['m_ceramic_g']/1000)   # [ug/g]            
        elif cal_method == 'None':  # areas are used
            CalibrComp = 'None'  # no calibration is applied
            # no concentration is available, conc_inj_mg_L is the areas
            conc_inj_mg_L = sample.loc[cas, 'Area%']/100  # fraction  
            # since its %area, dilution does not matter
            f_sample, yield_m = np.nan, np.nan  # no way of knowing
        # store information in the sample df for each identified cas (compound)
        sample.loc[cas, 'iupac_name'] =  CompInfo['iupac_name']       
        sample.loc[cas, 'molecular_weight'] = MW_cas    
        # make sure compound is classified in the ClassifDB
        ClassifDB = ClassificToDB(cas, ClassifDB, ClassifInstr, CompDB, 
                                  in_path, out_path)
        classifications = list(ClassifInstr.keys())  # list of classifications
        for cl in classifications:  # add classification tags in the sample df
            sample.loc[cas, cl] = ClassifDB.loc[cas, cl]
        sample.loc[cas, 'CalibrComp'] = CalibrComp
        sample.loc[cas, 'conc_inj_mg_L'] = conc_inj_mg_L  # mg/L
        sample.loc[cas, 'f_sample'] = f_sample # [-]
        sample.loc[cas, 'yield_m'] = yield_m # [g/g]  
    sample.to_excel(plib.Path(out_path, 'SingleSamples', file + '.xlsx'))
    PrintLog(out_path, 'f_sample detected: ' + 
             str(round(sample['f_sample'].sum(), 2)))
    PrintLog(out_path, file + " sample completed.") 
    return sample
    
    
def GCMS_folder_analysis(folder, f_sample_filt_thresh=0.01):
    shared_path = PathsCreate()[0]  # path to the common input folder (for all proj)
    ElemenFunGroups = pd.read_excel(plib.Path(shared_path, 'ElemenFunGroups.xlsx'),
                                    engine='openpyxl', sheet_name=None)
    try: # import the CompDB from the common input folder
        CompDB = pd.read_excel(plib.Path(shared_path, 'CompDB.xlsx'),
                               index_col='CAS')
    except FileNotFoundError:  # if the file is not there create an empty one
        CompDB = None  # so the import function can create a new one
    in_path, out_path = PathsCreate(folder)  # set path to project folders
    try:  # load DB with CAS to be substituted for the specific project
        CASsubstDB = pd.read_excel(plib.Path(in_path, 'CASsubstDB.xlsx'),
                                   engine='openpyxl') 
    except FileNotFoundError:  # use the shared one in the input folder
        CASsubstDB = pd.read_excel(plib.Path(shared_path, 'CASsubstDB.xlsx'),
                                   engine='openpyxl') 
    try:  # import the CalDB
        CalDB = pd.read_excel(plib.Path(in_path, 'CalDB.xlsx'), engine='openpyxl',
                              index_col='CAS') 
        cols_cal_area = ['Area 1', 'Area 2', 'Area 3', 'Area 4', 'Area 5']
        cols_cal_ppms = ['PPM 1', 'PPM 2', 'PPM 3', 'PPM 4', 'PPM 5']
        CalDB[cols_cal_area + cols_cal_ppms] = \
            CalDB[cols_cal_area + cols_cal_ppms].apply(pd.to_numeric, 
                                                       errors='coerce')
        cal_method = 'ExtStd'  # set the method
        RRFsDB = None
    except FileNotFoundError:  # if the file is not there create an empty 
        CalDB = None 
        try:  # check if the relative response factor DB is available
            RRFsDB = pd.read_excel(plib.Path(in_path, 'RRFsDB.xlsx'), 
                                   engine='openpyxl') 
            RRFsDB.set_index('CAS', inplace=True)  # set the CAS as index
            cal_method = 'IntStd'
        except FileNotFoundError:
            cal_method = 'None'
            RRFsDB = None
    
    # Classif instructions used to create and update ClassifDB 
    ClassifInstr = pd.read_excel(plib.Path(in_path, 'ClassifInstr.xlsx'), 
                                 engine='openpyxl', sheet_name=None)
    classifications = list(ClassifInstr.keys())  # list of classifications
    try: # if a ClassificationDB is available load it and retrieve classif names
        ClassifDB = pd.read_excel(plib.Path(in_path, 'ClassifDB.xlsx'),
                                  index_col='CAS')    
    except FileNotFoundError:  # if not available, load instructions to create it
        ClassifDB = None
    
    # load info about yield and dilutions from the specific folder in the input
    RunsInfo = pd.read_excel(plib.Path(in_path, folder + '.xlsx'), 
                             engine='openpyxl', index_col='Filename')
    filenames = RunsInfo.index.tolist()  # extract list of filenames
    samples = []  # create empty list and df to store results during looping 
    # create folder where single reports are stored
    plib.Path(out_path, 'SingleSamples').mkdir(parents=True, exist_ok=True)
    # =============================================================================
    # # loop for each sample
    # =============================================================================
    for f, file in enumerate(filenames):
        sample = GCMS_file_analysis(in_path, out_path, file, ElemenFunGroups, 
                                    CompDB, ClassifInstr, ClassifDB, CASsubstDB, 
                                    RunsInfo, CalDB, RRFsDB, cal_method)         
        samples.append(sample)
    # =============================================================================
    # # FINAL REPORTS WITH ALL SAMPLES TOGETHER
    # =============================================================================
    # Create Report (original, sorted, filtered) for all samples.
    # for each of these report, also create subreport with only 'conc_inj', 'f_sample'
    # 'yield'. save them all in the AllSample folder
    # create the AllSamples folder
    try:
        files = RunsInfo.loc[:, 'plt_label'].tolist()
    except KeyError:
        files = filenames
    # create the AllSamples folder
    plib.Path(out_path, 'AllSamples').mkdir(parents=True, exist_ok=True)
    # define the columns that will be put in the reports
    cols_init = ['iupac_name'] + classifications # left columns
    cols_conc_inj = ['conc_inj_mg_L_' + file for file in files]
    cols_f_sample = ['f_sample_' + file for file in files] 
    cols_yield_m = ['yield_m_' + file for file in files] 
    cols_conc = cols_conc_inj + cols_f_sample + cols_yield_m
    cols_final = ['str-CAS', 'CalibrComp', 'molecular_weight', 'Ret.Time'] 
    cols_all = cols_init + cols_conc + cols_final
    # create the report df, add common columns and renamed sample columns
    Rep = pd.DataFrame(columns=cols_all, index=pd.concat(samples, axis=1).index)
    for file, sample in zip(files, samples):
        Rep.update(sample.rename(columns={'conc_inj_mg_L': 'conc_inj_mg_L_'+ file,
                                          'f_sample': 'f_sample_'+ file,
                                          'yield_m': 'yield_m_'+ file}))
    # Original report    
    Rep.loc[:, cols_conc] = Rep.loc[:, cols_conc].fillna(0)  # put nan to zero
    Rep = Rep[Rep[cols_conc].sum(axis=1) > 0]  # remove all 0 rows
    # SortedReport sort report based on sum of f_sample (descending)
    RepSort = Rep.sort_index(key=Rep[cols_f_sample].sum(1).get, 
                             ascending=False, inplace=False)
    # FilterReport filter based on the minimum threshold of f_sample (default=0.01)
    RepFilt = RepSort.copy()  # copy the filtered df
    RepFilt.loc[(RepFilt[cols_f_sample] < f_sample_filt_thresh).all(axis=1), 
                   cols_conc] = 0  #  set all values in f_sample_cols to 0
    RepFilt = RepFilt[RepFilt[cols_conc].sum(axis=1) > 0] # remove all 0 rows
    # save all full Reps (original, sorted, filtered)
    Rep.to_excel(plib.Path(out_path, 'AllSamples', 'Rep.xlsx'))
    RepSort.to_excel(plib.Path(out_path, 'AllSamples', 'RepSort.xlsx'))
    RepFilt.to_excel(plib.Path(out_path, 'AllSamples', 'RepFilt.xlsx'))
    # produce Reps with only one of the concentration values for all types of rep
    for rep, lab in zip([Rep, RepSort, RepFilt], ['Rep_', 'RepSort_', 'RepFilt_']):
        for param, tag in zip(['conc_inj_mg_L', 'f_sample', 'yield_m'], 
                              ['inj', 'sample', 'yield']):
            cols = cols_init + [cc for cc in cols_conc if param in cc] + cols_final
            rep.loc[:, cols].to_excel(plib.Path(out_path, 'AllSamples', lab + tag 
                                                + '.xlsx'))
    # create the Aggreg folder for storing aggregated df
    plib.Path(out_path, 'Aggreg').mkdir(parents=True, exist_ok=True)
    # Compute aggregated results based on classes
    Aggs = []  # to store df for different classes
    for c, clsfc in enumerate(classifications):
        # for cols_conc, groupby class keys and sum to get aggregated values
        temp = Rep[[clsfc] + cols_conc
                   ].set_index(clsfc).groupby(clsfc)[cols_conc].sum()
        Agg = temp.sort_index(key=temp[cols_f_sample].sum(1).get, 
                                ascending=False, inplace=False)
        Aggs.append(Agg)  # stored Agg in Aggs (one Agg for each classific)
        Agg.to_excel(plib.Path(out_path, 'Aggreg', 'Agg_' + clsfc + '.xlsx'))
        # save aggregated df with only one parameter at a time
        for param, tag in zip(['conc_inj_mg_L', 'f_sample', 'yield_m'], 
                              ['inj', 'sample', 'yield']):
            cols = [cc for cc in cols_conc if param in cc]
            Agg.loc[:, cols].to_excel(plib.Path(out_path, 'Aggreg', 'Agg_' 
                                                + clsfc + '_' + tag + '.xlsx'))
    return files, Rep, RepSort, RepFilt, Aggs

    
def PlotAggreg(Aggs, folder, filenames, classification_idx, 
               param='f_sample', 
               transpose=False, rot=0, yLim=None, 
               threshold=0.01, leg_ncol=1, fig_name=None):
    """
    This function plots aggregate data based on a given parameter.
    
    Parameters
    ----------
    aggs : list of pandas DataFrames
        A list of dataframes containing the aggregate data.
    folder : str
        The name of the folder where the plots will be saved.
    filenames : list of str
        A list of strings containing the names to use for the plots.
    classification_idx : int
        The index of the dataframe in `aggs` to use for the plot.
    param : str, optional
        The parameter to use for the plot. Default is 'f_sample'.
    transpose : bool, optional
        If True, the data will be transposed before plotting. Default is False.
    rot : int, optional
        The rotation angle for the x-axis labels. Default is 0.
    ylim : tuple of float, optional
        The y-axis limits for the plot. Default is None.
    threshold : float, optional
        A threshold value used to filter the data. Default is 0.01.
    leg_ncol : int, optional
        The number of columns for the legend. Default is 1.
    fig_name : str, optional
        The name to use for the figure. If not provided, a default name will be used.
    
    Returns
    -------
    None
    """
    
    
    in_path, out_path = PathsCreate(folder)  # set path to project folders
    # list of tags for labeling plots based on parameters
    param2tag = dict(zip(['conc_inj_mg_L', 'f_sample', 'yield_m'], 
                         ['inj', 'sample', 'yield']))
    tag = param2tag[param]  # create tag
    
    cols = [cc for cc in list(Aggs[0]) if param in cc]  # select cols based on param
    classifications = [a.index.name for a in Aggs]
    classification = classifications[classification_idx]
    Agg = Aggs[classification_idx]
    AggPlot = Agg.loc[:, cols].copy()
    AggPlot.rename(index={'carboxylic acid': 'carb. acid'}, inplace=True)
    AggPlot.columns = filenames  # use nicer plt_labels if available
    if fig_name is None:         
        fig_name = 'Plt_' + classification + '_' + tag
    if not transpose:  # if True, samples names go into the legend
        AggPlot = AggPlot.transpose()
        
        fig_name += '_T'
    else:
        AggPlot.loc['total', :] = AggPlot.sum(axis=0)
    # filter values below a certain throshold     
    AggPlot = AggPlot.loc[:, AggPlot.max(axis=0) >= threshold]    
    # create figure
    fig, ax, axt, fig_par = FigCreate(rows=1, cols=1, plot_type=0,
                                      paper_col=1.5)
    # plot bars
    AggPlot.plot(ax=ax[0], kind='bar', rot=rot, width=.9, edgecolor='k', 
                 legend=False)
    if not transpose:  # samples are on the x-axis so plot total as scatter
        ax[0].scatter(AggPlot.index, AggPlot.sum(axis=1), color='k', 
                      linestyle='None', edgecolor='k', s=50, label='total')
    bars = ax[0].patches  # add patches to the bars
    hatches = [p for p in htchs for i in range(len(AggPlot))]
    for bar, hatch in zip(bars, hatches):
        bar.set_hatch(hatch)
    ax[0].legend(bbox_to_anchor=(1, 1.0), ncol=leg_ncol)
    ax[0].set(xlabel=None)
    FigSave(fig_name, out_path, fig, ax, axt, fig_par,
            yLab='Yield [-]', yLim=yLim, legend=False, tight_layout=True)

class fragmenter:
    # tested with Python 3.8.8 and RDKit version 2021.09.4
    
    from rdkit import Chem
    import marshal as marshal
    from rdkit.Chem import rdmolops

    # does a substructure match and then checks whether the match 
    # is adjacent to previous matches
    @classmethod
    def get_substruct_matches(cls, mol_searched_for, mol_searched_in, atomIdxs_to_which_new_matches_have_to_be_adjacent):
        
        valid_matches = []
        
        if mol_searched_in.GetNumAtoms() >= mol_searched_for.GetNumAtoms():
            matches = mol_searched_in.GetSubstructMatches(mol_searched_for)
            
            if matches:
                for match in matches:
                        add_this_match = True
                        if len(atomIdxs_to_which_new_matches_have_to_be_adjacent) > 0:
                            add_this_match = False
                            
                            for i in match:
                                for neighbor in mol_searched_in.GetAtomWithIdx(i).GetNeighbors():
                                    if neighbor.GetIdx() in atomIdxs_to_which_new_matches_have_to_be_adjacent:
                                        add_this_match = True
                                        break
                                
                        if add_this_match:
                            valid_matches.append(match)
                
        return valid_matches
    
    # count heavier isotopes of hydrogen correctly
    @classmethod
    def get_heavy_atom_count(cls, mol):
        heavy_atom_count = 0
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 1:
                heavy_atom_count += 1
        
        return heavy_atom_count
    
    def __init__(self, fragmentation_scheme = {}, fragmentation_scheme_order = None, match_hydrogens = False, algorithm = '', n_atoms_cuttoff = -1, function_to_choose_fragmentation = False, n_max_fragmentations_to_find = -1):

        if not type(fragmentation_scheme) is dict:
            raise TypeError('fragmentation_scheme must be a dctionary with integers as keys and either strings or list of strings as values.')
            
        if len(fragmentation_scheme) == 0:
            raise ValueError('fragmentation_scheme must be provided.')
        
        if not algorithm in ['simple', 'complete', 'combined']:
            raise ValueError('Algorithm must be either simple ,complete or combined.')
            
        if algorithm == 'simple':
            if n_max_fragmentations_to_find != -1:
                raise ValueError('Setting n_max_fragmentations_to_find only makes sense with complete or combined algorithm.')
        
        self.algorithm = algorithm
        
        if algorithm in ['combined', 'complete']:
            if n_atoms_cuttoff == -1:
                raise ValueError('n_atoms_cuttoff needs to be specified for complete or combined algorithms.')
                
            if function_to_choose_fragmentation == False:
                raise ValueError('function_to_choose_fragmentation needs to be specified for complete or combined algorithms.')
                
            if not callable(function_to_choose_fragmentation):
                raise TypeError('function_to_choose_fragmentation needs to be a function.')
            else:
                if type(function_to_choose_fragmentation([{}, {}])) != dict:
                    raise TypeError('function_to_choose_fragmentation needs to take a list of fragmentations and choose one of it')
                
            if n_max_fragmentations_to_find != -1:
                if n_max_fragmentations_to_find < 1:
                    raise ValueError('n_max_fragmentations_to_find has to be 1 or higher.')

        if fragmentation_scheme_order is None:
            fragmentation_scheme_order = []

        if algorithm in ['simple', 'combined']:
            assert len(fragmentation_scheme) == len(fragmentation_scheme_order)
        else:
            fragmentation_scheme_order = [key for key in fragmentation_scheme.keys()]
            
        self.n_max_fragmentations_to_find = n_max_fragmentations_to_find
        
        self.n_atoms_cuttoff = n_atoms_cuttoff

        self.match_hydrogens = match_hydrogens
        
        self.fragmentation_scheme = fragmentation_scheme
        
        self.function_to_choose_fragmentation = function_to_choose_fragmentation
        
        # create a lookup dictionaries to faster finding a group number
        self._fragmentation_scheme_group_number_lookup = {}
        self._fragmentation_scheme_pattern_lookup = {}
        self.fragmentation_scheme_order = fragmentation_scheme_order

        for group_number, list_SMARTS in fragmentation_scheme.items():
            
            if type(list_SMARTS) is not list:
                list_SMARTS = [list_SMARTS]
                
            for SMARTS in list_SMARTS:
                if SMARTS != '':
                    self._fragmentation_scheme_group_number_lookup[SMARTS] = group_number
                    
                    mol_SMARTS = fragmenter.Chem.MolFromSmarts(SMARTS)
                    self._fragmentation_scheme_pattern_lookup[SMARTS] = mol_SMARTS

    def fragment(self, SMILES_or_molecule):
        
        if type(SMILES_or_molecule) is str:
            mol_SMILES = fragmenter.Chem.MolFromSmiles(SMILES_or_molecule)
            mol_SMILES = fragmenter.Chem.AddHs(mol_SMILES) if self.match_hydrogens else mol_SMILES
            is_valid_SMILES = mol_SMILES is not None
            
            if not is_valid_SMILES:
                raise ValueError('Following SMILES is not valid: ' + SMILES_or_molecule)

        else:
            mol_SMILES = SMILES_or_molecule

        # iterate over all separated molecules
        success = []
        fragmentation = {}
        fragmentation_matches = {}
        for mol in fragmenter.rdmolops.GetMolFrags(mol_SMILES, asMols = True):

            this_mol_fragmentation, this_mol_success = self.__get_fragmentation(mol)
    
            for SMARTS, matches in this_mol_fragmentation.items():
                group_number = self._fragmentation_scheme_group_number_lookup[SMARTS]
                
                if not group_number in fragmentation:
                    fragmentation[group_number] = 0
                    fragmentation_matches[group_number] = []
                
                fragmentation[group_number] += len(matches)
                fragmentation_matches[group_number].extend(matches)   
                
            success.append(this_mol_success)

        return fragmentation, all(success), fragmentation_matches
    
    def fragment_complete(self, SMILES_or_molecule):
        
        if type(SMILES_or_molecule) is str:
            mol_SMILES = fragmenter.Chem.MolFromSmiles(SMILES_or_molecule)
            mol_SMILES = fragmenter.Chem.AddHs(mol_SMILES) if self.match_hydrogens else mol_SMILES
            is_valid_SMILES = mol_SMILES is not None
            
            if not is_valid_SMILES:
                raise ValueError('Following SMILES is not valid: ' + SMILES_or_molecule)

        else:
            mol_SMILES = SMILES_or_molecule

        if len(fragmenter.rdmolops.GetMolFrags(mol_SMILES)) != 1:
            raise ValueError('fragment_complete does not accept multifragment molecules.')

        temp_fragmentations, success = self.__complete_fragmentation(mol_SMILES)

        fragmentations = []
        fragmentations_matches = []
        for temp_fragmentation in temp_fragmentations:
            fragmentation = {}
            fragmentation_matches = {}
            for SMARTS, matches in temp_fragmentation.items():
                group_number = self._fragmentation_scheme_group_number_lookup[SMARTS]

                fragmentation[group_number] = len(matches)
                fragmentation_matches[group_number] = matches

            fragmentations.append(fragmentation)
            fragmentations_matches.append(fragmentation_matches)

        return fragmentations, success, fragmentations_matches


    def __get_fragmentation(self, mol_SMILES):
        
        success = False
        fragmentation = {}
        if self.algorithm in ['simple', 'combined']:
            fragmentation, success = self.__simple_fragmentation(mol_SMILES)
        
        if success:
            return fragmentation, success
        
        if self.algorithm in ['combined', 'complete']:
            fragmentations, success = self.__complete_fragmentation(mol_SMILES)
            
            if success:
                fragmentation = self.function_to_choose_fragmentation(fragmentations)
        
        return fragmentation, success
    
    def __simple_fragmentation(self, mol_SMILES):

        if self.match_hydrogens:
            target_atom_count = len(mol_SMILES.GetAtoms())
        else:
            target_atom_count = fragmenter.get_heavy_atom_count(mol_SMILES)
    
        success = False
        fragmentation = {}
        
        fragmentation, atomIdxs_included_in_fragmentation = self.__search_non_overlapping_solution(mol_SMILES, {}, set(), set())
        success = len(atomIdxs_included_in_fragmentation) == target_atom_count
        
        # if not successful, clean up molecule and search again
        level = 1
        while not success:
            fragmentation_so_far , atomIdxs_included_in_fragmentation_so_far = fragmenter.__clean_molecule_surrounding_unmatched_atoms(mol_SMILES, fragmentation, atomIdxs_included_in_fragmentation, level)
            level += 1
            
            if len(atomIdxs_included_in_fragmentation_so_far) == 0:
                break
            
            fragmentation_so_far, atomIdxs_included_in_fragmentation_so_far = self.__search_non_overlapping_solution(mol_SMILES, fragmentation_so_far, atomIdxs_included_in_fragmentation_so_far, atomIdxs_included_in_fragmentation_so_far)
            
            success = len(atomIdxs_included_in_fragmentation_so_far) == target_atom_count
            
            if success:
                fragmentation = fragmentation_so_far
            
        return fragmentation, success
    
    def __search_non_overlapping_solution(self, mol_searched_in, fragmentation, atomIdxs_included_in_fragmentation, atomIdxs_to_which_new_matches_have_to_be_adjacent):
        
        n_atomIdxs_included_in_fragmentation = len(atomIdxs_included_in_fragmentation) - 1
        
        while n_atomIdxs_included_in_fragmentation != len(atomIdxs_included_in_fragmentation):
            n_atomIdxs_included_in_fragmentation = len(atomIdxs_included_in_fragmentation)
            
             
            for group_number in self.fragmentation_scheme_order:
                list_SMARTS = self.fragmentation_scheme[group_number]
                if type(list_SMARTS) is not list:
                    list_SMARTS = [list_SMARTS]
                
                for SMARTS in list_SMARTS:
                    if SMARTS != "":  
                        fragmentation, atomIdxs_included_in_fragmentation = self.__get_next_non_overlapping_match(mol_searched_in, SMARTS, fragmentation, atomIdxs_included_in_fragmentation, atomIdxs_to_which_new_matches_have_to_be_adjacent)

        return fragmentation, atomIdxs_included_in_fragmentation

    def __get_next_non_overlapping_match(self, mol_searched_in, SMARTS, fragmentation, atomIdxs_included_in_fragmentation, atomIdxs_to_which_new_matches_have_to_be_adjacent):
        
        mol_searched_for = self._fragmentation_scheme_pattern_lookup[SMARTS]
        
        if atomIdxs_to_which_new_matches_have_to_be_adjacent:
            matches = fragmenter.get_substruct_matches(mol_searched_for, mol_searched_in, atomIdxs_to_which_new_matches_have_to_be_adjacent)
        else:
            matches = fragmenter.get_substruct_matches(mol_searched_for, mol_searched_in, set())
        
        if matches:
            for match in matches:
                all_atoms_of_new_match_are_unassigned = atomIdxs_included_in_fragmentation.isdisjoint(match)
        
                if all_atoms_of_new_match_are_unassigned:
                    if not SMARTS in fragmentation:
                        fragmentation[SMARTS] = []
                        
                    fragmentation[SMARTS].append(match)
                    atomIdxs_included_in_fragmentation.update(match)   
                
        return fragmentation, atomIdxs_included_in_fragmentation

    @classmethod
    def __clean_molecule_surrounding_unmatched_atoms(cls, mol_searched_in, fragmentation, atomIdxs_included_in_fragmentation, level):
    
        for i in range(0, level):
            
            atoms_missing = set(range(0, fragmenter.get_heavy_atom_count(mol_searched_in))).difference(atomIdxs_included_in_fragmentation)
                        
            new_fragmentation = fragmenter.marshal.loads(fragmenter.marshal.dumps(fragmentation))
            
            for atomIdx in atoms_missing:
                for neighbor in mol_searched_in.GetAtomWithIdx(atomIdx).GetNeighbors():
                    for smart, atoms_found in fragmentation.items():
                        for atoms in atoms_found:
                            if neighbor.GetIdx() in atoms:
                                if smart in new_fragmentation:
                                    if new_fragmentation[smart].count(atoms) > 0:
                                        new_fragmentation[smart].remove(atoms)
                                
                        if smart in new_fragmentation:
                            if len(new_fragmentation[smart]) == 0:
                                new_fragmentation.pop(smart)
                                
                                
            new_atomIdxs_included_in_fragmentation = set()
            for i in new_fragmentation.values():
                for j in i:
                    new_atomIdxs_included_in_fragmentation.update(j)
                    
            atomIdxs_included_in_fragmentation = new_atomIdxs_included_in_fragmentation
            fragmentation = new_fragmentation
            
        return fragmentation, atomIdxs_included_in_fragmentation
       
    def __complete_fragmentation(self, mol_SMILES):
    
        heavy_atom_count = fragmenter.get_heavy_atom_count(mol_SMILES)
        
        if heavy_atom_count > self.n_atoms_cuttoff:
            return {}, False
        
        completed_fragmentations = []
        groups_leading_to_incomplete_fragmentations = []        
        completed_fragmentations, groups_leading_to_incomplete_fragmentations, incomplete_fragmentation_found = self.__get_next_non_overlapping_adjacent_match_recursively(mol_SMILES, heavy_atom_count, completed_fragmentations, groups_leading_to_incomplete_fragmentations, {}, set(), set(), self.n_max_fragmentations_to_find)
        success = len(completed_fragmentations) > 0
        
        return completed_fragmentations, success
        
    def __get_next_non_overlapping_adjacent_match_recursively(self, mol_searched_in, heavy_atom_count, completed_fragmentations, groups_leading_to_incomplete_fragmentations, fragmentation_so_far, atomIdxs_included_in_fragmentation_so_far, atomIdxs_to_which_new_matches_have_to_be_adjacent, n_max_fragmentations_to_find = -1):
      
        n_completed_fragmentations = len(completed_fragmentations)
        incomplete_fragmentation_found = False
        complete_fragmentation_found = False
        
        if len(completed_fragmentations) == n_max_fragmentations_to_find:
            return completed_fragmentations, groups_leading_to_incomplete_fragmentations, incomplete_fragmentation_found
                    
                
        for group_number in self.fragmentation_scheme_order:
            list_SMARTS = self.fragmentation_scheme[group_number]
            
            if complete_fragmentation_found:
                break
            
            if type(list_SMARTS) is not list:
                list_SMARTS = [list_SMARTS]
            
            for SMARTS in list_SMARTS:
                if complete_fragmentation_found:
                    break
                
                if SMARTS != "":  
                    matches = fragmenter.get_substruct_matches(self._fragmentation_scheme_pattern_lookup[SMARTS], mol_searched_in, atomIdxs_included_in_fragmentation_so_far)
                    
                    for match in matches:
                        
                        # only allow non-overlapping matches
                        all_atoms_are_unassigned = atomIdxs_included_in_fragmentation_so_far.isdisjoint(match)
                        if not all_atoms_are_unassigned:
                            continue
                        
                        # only allow matches that do not contain groups leading to incomplete matches
                        for groups_leading_to_incomplete_fragmentation in groups_leading_to_incomplete_fragmentations:
                            if fragmenter.__is_fragmentation_subset_of_other_fragmentation(groups_leading_to_incomplete_fragmentation, fragmentation_so_far):
                                return completed_fragmentations, groups_leading_to_incomplete_fragmentations, incomplete_fragmentation_found
                        
                        # only allow matches that will lead to new fragmentations
                        use_this_match = True
                        n_found_groups = len(fragmentation_so_far)
                        
                        for completed_fragmentation in completed_fragmentations:
                            
                            if not SMARTS in completed_fragmentation:
                                continue
                            
                            if n_found_groups == 0:
                                use_this_match = not fragmenter.__is_match_contained_in_fragmentation(match, SMARTS, completed_fragmentation)
                            else:
                                if fragmenter.__is_fragmentation_subset_of_other_fragmentation(fragmentation_so_far, completed_fragmentation):
                                    use_this_match = not fragmenter.__is_match_contained_in_fragmentation(match, SMARTS, completed_fragmentation)
                            
                            if not use_this_match:
                                break
                                
                        if not use_this_match:
                            continue
                        
                        # make a deepcopy here, otherwise the variables are modified down the road
                        # marshal is used here because it works faster than copy.deepcopy
                        this_SMARTS_fragmentation_so_far = fragmenter.marshal.loads(fragmenter.marshal.dumps(fragmentation_so_far))
                        this_SMARTS_atomIdxs_included_in_fragmentation_so_far = atomIdxs_included_in_fragmentation_so_far.copy()
                        
                        if not SMARTS in this_SMARTS_fragmentation_so_far:
                            this_SMARTS_fragmentation_so_far[SMARTS] = []
                            
                        this_SMARTS_fragmentation_so_far[SMARTS].append(match)
                        this_SMARTS_atomIdxs_included_in_fragmentation_so_far.update(match)
                        
                        # only allow matches that do not contain groups leading to incomplete matches
                        for groups_leading_to_incomplete_match in groups_leading_to_incomplete_fragmentations:
                            if fragmenter.__is_fragmentation_subset_of_other_fragmentation(groups_leading_to_incomplete_match, this_SMARTS_fragmentation_so_far):
                                use_this_match = False
                                break
                            
                        if not use_this_match:
                            continue
                        
                        # if the complete molecule has not been fragmented, continue to do so
                        if len(this_SMARTS_atomIdxs_included_in_fragmentation_so_far) < heavy_atom_count:
                            completed_fragmentations, groups_leading_to_incomplete_fragmentations, incomplete_fragmentation_found = self.__get_next_non_overlapping_adjacent_match_recursively(mol_searched_in, heavy_atom_count, completed_fragmentations, groups_leading_to_incomplete_fragmentations, this_SMARTS_fragmentation_so_far, this_SMARTS_atomIdxs_included_in_fragmentation_so_far, this_SMARTS_atomIdxs_included_in_fragmentation_so_far, n_max_fragmentations_to_find)
                            break
                        
                        # if the complete molecule has been fragmented, save and return
                        if len(this_SMARTS_atomIdxs_included_in_fragmentation_so_far) == heavy_atom_count:                 
                            completed_fragmentations.append(this_SMARTS_fragmentation_so_far)
                            complete_fragmentation_found = True
                            break
                        
        # if until here no new fragmentation was found check whether an incomplete fragmentation was found
        if n_completed_fragmentations == len(completed_fragmentations):      
            
            if not incomplete_fragmentation_found:
                
                incomplete_matched_groups = {}
                
                if len(atomIdxs_included_in_fragmentation_so_far) > 0:
                    unassignes_atom_idx = set(range(0, heavy_atom_count)).difference(atomIdxs_included_in_fragmentation_so_far)
                    for atom_idx in unassignes_atom_idx:
                        neighbor_atoms_idx = [i.GetIdx() for i in mol_searched_in.GetAtomWithIdx(atom_idx).GetNeighbors()]
                        
                        for neighbor_atom_idx in neighbor_atoms_idx:
                            for found_smarts, found_matches in fragmentation_so_far.items():
                                for found_match in found_matches:
                                    if neighbor_atom_idx in found_match:
                                        if not found_smarts in incomplete_matched_groups:
                                            incomplete_matched_groups[found_smarts] = []
                                            
                                        if found_match not in incomplete_matched_groups[found_smarts]:
                                            incomplete_matched_groups[found_smarts].append(found_match)
                                    
                    is_subset_of_groups_already_found = False
                    indexes_to_remove = []
                    
                    for idx, groups_leading_to_incomplete_match in enumerate(groups_leading_to_incomplete_fragmentations):
                        is_subset_of_groups_already_found = fragmenter.__is_fragmentation_subset_of_other_fragmentation(incomplete_matched_groups, groups_leading_to_incomplete_match)
                        if is_subset_of_groups_already_found:
                            indexes_to_remove.append(idx)
                        
                    for index in sorted(indexes_to_remove, reverse=True):
                        del groups_leading_to_incomplete_fragmentations[index]
                        
                    groups_leading_to_incomplete_fragmentations.append(incomplete_matched_groups)
                    groups_leading_to_incomplete_fragmentations = sorted(groups_leading_to_incomplete_fragmentations, key = len)
                    
                    incomplete_fragmentation_found =  True
    
        return completed_fragmentations, groups_leading_to_incomplete_fragmentations, incomplete_fragmentation_found
    
    @classmethod
    def __is_fragmentation_subset_of_other_fragmentation(cls, fragmentation, other_fragmentation):
        n_found_groups = len(fragmentation)
        n_found_other_groups = len(other_fragmentation)
        
        if n_found_groups == 0:
            return False
            
        if n_found_other_groups < n_found_groups:
            return False
        
        n_found_SMARTS_that_are_subset = 0
        for found_SMARTS, _ in fragmentation.items():
            if found_SMARTS in other_fragmentation:
                found_matches_set = set(frozenset(i) for i in fragmentation[found_SMARTS])
                found_other_matches_set =  set(frozenset(i) for i in other_fragmentation[found_SMARTS])
                
                if found_matches_set.issubset(found_other_matches_set):
                    n_found_SMARTS_that_are_subset += 1
            else:
                return False
            
        return n_found_SMARTS_that_are_subset == n_found_groups
    
    @classmethod
    def __is_match_contained_in_fragmentation(cls, match, SMARTS, fragmentation):
        if not SMARTS in fragmentation:
            return False
            
        found_matches_set = set(frozenset(i) for i in fragmentation[SMARTS])
        match_set = set(match)
        
        return match_set in found_matches_set
  
    
# =============================================================================
# # GRAPHICAL FUNCTIONS
# =============================================================================

sns.set_style("whitegrid")
clrs = sns.color_palette('deep', 10)   # list with colors
# list with linestyles for plotting
lnstls = [(0, ()),  # solid
          (0, (1, 1)),  # 'densely dotted'
          (0, (5, 1)),  # 'densely dashed'
          (0, (3, 1, 1, 1)),  # 'densely dashdotted'
          (0, (3, 1, 1, 1, 1, 1)),  # 'densely dashdotdotted'
          (0, (5, 5)),  # 'dashed'
          (0, (3, 5, 1, 5)),  # 'dashdotted'
          (0, (1, 5)),  # dotted
          (0, (3, 5, 1, 5, 1, 5)),  # 'dashdotdotted'
          (0, (1, 10)),  # 'loosely dotted'
          (0, (5, 10)),  # 'loosely dashed'
          (0, (3, 10, 1, 10)),  # 'loosely dashdotted'
          (0, (3, 10, 1, 10, 1, 10))]  # 'loosely dashdotdotted'
# list with letters for plotting
lttrs = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm']
# list with markers for plotting
mrkrs = ["o", "v", "X", "s", "p", "^", "P", "<", ">", "*", "d"]
# list with hatches for plotting
htchs = (None, '//', '...', '--', 'O', '\\\\', 'oo', '\\\\\\', '/////', '.....')      

def FigCreate(rows=1, cols=1, plot_type=0, paper_col=1,
              gridspec_hgt_fr=None,
              gridspec_wdt_fr=None, hgt_mltp=1, font='DejaVu Sans'):
    """
    This function creates all the necessary objects to produce plots with
    replicable characteristics.

    Parameters
    ----------
    rows : int, optional
        number of plot rows in the grid. The default is 1.
    cols : int, optional
        number of plot colums in the grid. The default is 1.
    plot_type : int, optional
        one of the different plots available. The default is 0.
        Plot types and their labels:
        0. Std: standard plot (single or grid rows x cols)
        1. Twin-x: secondary axis plot (single or grid rows x cols)
        2. Cmap: colormap plot (single or grid rows x cols)
        3. Ir: Infra-red curves plot (no brokenXaxis)
        4. Ir: Infra-red curves plot (WITH brokenXaxis)
    paper_col : int, optional
        single or double column size for the plot, meaning the actual space
        it will fit in a paper. The default is 1.
    gridspec_wdt_ratio : list of float, optional
        for multiple cols, list of relative width of each subplot.
        The default is None.
    no_par : bool, optional
        if True, no size setting parameters are passed. The default is None.
    font: str
        if the str'Times' is given, it sets Times New Roman as the default
        font for the plot, otherwise the default Dejavu Sans is maintained.
        Default is 'Dejavu Sans'
    hgt_mltp: float
        multiplies the fig height. default is to 1. best using values between
        .65 and 2. may not work with multiplot and paper_col=1 or out of the
        specified range
    Returns
    -------
    fig : object
        the figure object to be passed to FigSave.
    lst_ax : list of axis
        list of axis (it is a list even with 1 axis) on which to plot.
    lst_axt : list of axis
        list of secondary axis (it is a list even with 1 axis)
    fig_par : list of float
        ist of parameters to reserve space around the plot canvas.

    """

    if font == 'Times':  # set Times New Roman as the plot font fot text
        sns.set("paper", font_scale=1.25)
        sns.set_palette('deep', 10)
        sns.set_style("ticks", {'font.family': 'Times New Roman'})
    else:  # leave Dejavu Sans (default) as the plot font fot text
        sns.set("paper", font_scale=1.25)
        sns.set_palette('deep', 10)
        sns.set_style("ticks")
    # single or double column in paperthat the figure will occupy
    if cols > 3:  # numer of columns (thus of plots in the figure)
        raise ValueError('\n FigCreate: cols>2 not supported')

    # width of the figure in inches, it's fixed to keep the same text size
    # is 6, 9, 12 for 1, 1.5, and 3 paper_col (columns in paper)
    fig_wdt = 6*paper_col  # width of the plot in inches
    fig_hgt = 4*paper_col*rows/cols*hgt_mltp  # heigth of the figure in inches
    px = 0.06*(6/fig_wdt)*cols  # set px so that (A) fits the square
    py = px*fig_wdt/fig_hgt/cols*rows/hgt_mltp  # set py so that (A) fits
    # if more rows are added, it increases, but if cols areadded it decreases
    # to maintain the plot ratio
    # set plot margins
    sp_lab_wdt = 0.156/paper_col  # hor. space for labels
    sp_nar_wdt = 0.02294/paper_col  # space narrow no labels (horiz)
    sp_lab_hgt = 0.147/paper_col/rows*cols/hgt_mltp  # space for labels (vert)
    sp_nar_hgt = 0.02/paper_col/rows*cols/hgt_mltp  # space narrow no labels
    # (vert)
    # =========================================================================
    # # 0. Std: standard plot (single or grid rows x cols)
    # =========================================================================
    if plot_type == 0:
        fig, ax = plt.subplots(rows, cols, figsize=(fig_wdt, fig_hgt))
        if rows*cols == 1:  # only 1 plot
            lst_ax = [ax]  # create ax list for uniform iterations over 1 obj.
        elif rows*cols > 1:  # more than one plot
            lst_ax = [axs for axs in ax.flatten()]  # create list of axis
        lst_axt = None  # no secondary axis in this plot_type
        # horizontal space between plot in percentage
        sp_btp_wdt = (0.26*paper_col**2 - 1.09*paper_col + 1.35)
        # vertical space between plot in percentage !!! needs DEBUG
        sp_btp_hgt = .2/paper_col*cols/hgt_mltp
        # left, bottom, right, top, widthspace, heightspace
        fig_par = [sp_lab_wdt, sp_lab_hgt, 1-sp_nar_wdt, 1-sp_nar_hgt,
                   sp_btp_wdt, sp_btp_hgt, px, py]
    # =========================================================================
    # # 1. Twin-x: secondary axis plot (single or grid rows x cols)
    # =========================================================================
    elif plot_type == 1:
        fig, ax = plt.subplots(rows, cols, figsize=(fig_wdt, fig_hgt))
        if rows*cols == 1:  # only 1 plot
            lst_ax = [ax]  # create ax list for uniform iterations over 1 obj.
            lst_axt = [ax.twinx()]  # create a list with secondary axis object
        elif rows*cols > 1:  # more than one plot
            lst_ax = [axs for axs in ax.flatten()]  # create list of axis
            # create list of secondary twin axis
            lst_axt = [axs.twinx() for axs in ax.flatten()]
        # horizontal space between plot in percentage !!! needs DEBUG
        sp_btp_wdt = 1.36*paper_col**2 - 5.28*paper_col + 5.57
        # vertical space between plot in percentage !!! needs DEBUG
        sp_btp_hgt = .2/paper_col*cols/hgt_mltp
        # left, bottom, right(DIFFERENT FROM STD), top, widthspace, heightspace
        fig_par = [sp_lab_wdt, sp_lab_hgt, 1-sp_lab_wdt, 1-sp_nar_hgt,
                   sp_btp_wdt, sp_btp_hgt, px, py]

    return fig, lst_ax, lst_axt, fig_par


def FigSave(fig_name, out_path, fig, lst_ax, lst_axt, fig_par,
            xLab=None, yLab=None, ytLab=None,
            xLim=None, yLim=None, ytLim=None,
            xTicks=None, yTicks=None, ytTicks=None,
            xTickLabels=None, yTickLabels=None, ytTickLabels=None,
            legend=None, ncol_leg=1, 
            annotate_lttrs=False, annotate_lttrs_loc='down',
            pdf=False, svg=False, transparency=False,
            subfolder=None, tight_layout=False, title=False,
            ):
    '''
    FIXES:
        1. px, py moved to FIgCreate

    This function takes the obects created in FigCreate and allows to modify
    their appeareance and saving the results.

    Parameters
    ----------
    fig_name : str
        name of figure. It is the name of the png od pfd file to be saved
    out_path : pathlib.Path object. path to the output folder.
    fig : figure object. created in FigSave.
    lst_ax : list of axis. Created in FigCreate
    lst_axt : list of twin (secondary) axis. Created in FigCreate
    fig_par : list of figure parameters for space settings
        left, bottom, right, top, widthspace, heightspace, px, py.
        Created in FigCreate
    tight_layout : bool
        If True, ignore fig_par[0:6] and fit the figure to the tightest layout
        possible. Avoids to lose part of figure, but loses control of margins
    xLab : str.list, optional
        label of the x axis. The default is None.
        can be given as
        0. xLab=None: no axis gets an xlabel
        1. xLab='label': only one str, all axis get the same xlabel
        2. xLab=['label1', None, Label2, ...]: the list must have the size of
            lst_ax and contain labels and or None values. Each axis is
            assigned its label, where None is given, no label is set.
    yLab : str, optional
        label of the y axis. The default is None. Same options as xLab
    ytLab : str, optional
        label of the secondary y-axis. The default is None.
        Same options as xLab
    xLim : list of two values, list of lists, optional
        limits of x axis. The default is None.
        can be given as
        0. xLim=None: no axis gets a xlim
        1. xLab=[a,b]: all axis get the same xlim
        2. xLab=[[a,b], None, [c,d], ...]: the list must have the size of
            lst_ax and contain [a,b] and or None values. Each axis is
            assigned its limit, where None is given, no llimit is set.
    yLim : list of two values, optional
        limits of y axis. The default is None. Same options as xLim
    ytLim : list of two values, optional
        limits of secondary y axis. The default is None.
        Same options as xLim
    xTicks : list of int or float, optional
        list of tiks value to be shown on the axis. The default is None.
    yTicks : list of int or float, optional
        list of tiks value to be shown on the axis. The default is None.
    ytTicks : TYPE, optional
        list of tiks value to be shown on the axis. The default is None.
    legend : str, optional
        contains info on the legend location. To avoid printing the legend
        (also in case it is empty) set it to None.
        The default is 'best'.
    ncol_leg : int, optional
        number of columns in the legend. The default is 1.
    annotate_lttrs : bool, optional
        if True, each plot is assigned a letter between () in the lower left
        corner. The default is False. If a string is given, the string is used
        as the letter in the plot even for single plots.
    annotate_lttrs_loc: str.
        default is 'down', if 'up' is given, the letters are placed on the left
        top corner.
    pdf : bool, optional
        if True, the figure is saved also in pdf in the output folder.
        The default is False, so only a png file with 300dpi is saved
    transparency : bool, optional
        if True, background of PNG figure is transparent, defautls is False.
    subfolder : str, optional
        name of the subfolder inside the output folder where the output will
        be saved. If the folder does not exists, it is created.
        The default is None.
    '''

    fig_adj_par = fig_par[0:6]
    if not any(fig_par[0:6]):  # True if all element in fig_par[0:6] are False
        tight_layout = True
    px = fig_par[6]
    py = fig_par[7]
    n_ax = len(lst_ax)  # number of ax objects
    # for xLab, yLab, ytLab creates a list with same length as n_ax.
    # only one value is given all axis are given the same label
    # if a list is given, each axis is given a different value, where False
    # is specified, no value is given to that particular axis
    vrbls = [xLab, yLab, ytLab, legend]  # collect variables for iteration
    lst_xLab, lst_yLab, lst_ytLab, lst_legend \
        = [], [], [], []  # create lists for iteration
    lst_vrbls = [lst_xLab, lst_yLab, lst_ytLab, lst_legend]  # collect lists
    for vrbl, lst_vrbl in zip(vrbls, lst_vrbls):
        if vrbl is None:  # label is not given for any axis
            lst_vrbl[:] = [None]*n_ax
        else:  # label is given
            if np.size(vrbl) == 1:  # only one value is given
                if type(vrbl) == str:  # create a list before replicating it
                    lst_vrbl[:] = [vrbl]*n_ax  # each axis gets same label
                elif type(vrbl) == list:  # replicate the list
                    lst_vrbl[:] = vrbl*n_ax  # each axis gets same label
            elif np.size(vrbl) == n_ax:  # each axis has been assigned its lab
                lst_vrbl[:] = vrbl  # copy the label inside the list
            else:
                print(vrbl)
                print('Labels/legend size does not match axes number')
    # for xLim, yLim, ytLim creates a list with same length as n_ax.
    # If one list like [a,b] is given, all axis have the same limits, if a list
    # of the same length of the axis is given, each axis has its lim. Where
    # None is given, no lim is set on that axis
    vrbls = [xLim, yLim, ytLim, xTicks, yTicks, ytTicks, xTickLabels,
             yTickLabels, ytTickLabels]  # collect variables for iteration
    lst_xLim, lst_yLim, lst_ytLim, lst_xTicks, lst_yTicks, lst_ytTicks, \
        lst_xTickLabels, lst_yTickLabels, lst_ytTickLabels = \
            [], [], [], [], [], [], [], [], [] # create lists for iteration
    lst_vrbls = [lst_xLim, lst_yLim, lst_ytLim, lst_xTicks, lst_yTicks,
                 lst_ytTicks, lst_xTickLabels, lst_yTickLabels,
                 lst_ytTickLabels]  # collect lists
    for vrbl, lst_vrbl in zip(vrbls, lst_vrbls):
        if vrbl is None:  # limit is not given for any axis
            lst_vrbl[:] = [None]*n_ax
        else:
            # if only list and None are in vrbl, it is [[], None, [], ..]
            # each axis has been assigned its limits
            if any([isinstance(v, (int, float, np.int32, str))
                    for v in vrbl]):
                temporary = []  # necessary to allow append on [:]
                for i in range(n_ax):
                    temporary.append(vrbl)  # give it to all axis
                lst_vrbl[:] = temporary
            else:  # xLim=[[a,b], None, ...] = [list, bool] # no float
                lst_vrbl[:] = vrbl  # a lim for each axis is already given
    # loops over each axs in the ax array and set the different properties
    for i, axs in enumerate(lst_ax):
        # for each property, if the variable is not false, it is set
        if lst_xLab[i] is not None:
            axs.set_xlabel(lst_xLab[i])
        if lst_yLab[i] is not None:
            axs.set_ylabel(lst_yLab[i])
        if lst_xLim[i] is not None:
            axs.set_xlim([lst_xLim[i][0]*(1 + px) - px*lst_xLim[i][1],
                          lst_xLim[i][1]*(1 + px) - px*lst_xLim[i][0]])
        if lst_yLim[i] is not None:
            axs.set_ylim([lst_yLim[i][0]*(1 + py) - py*lst_yLim[i][1],
                          lst_yLim[i][1]*(1 + py) - py*lst_yLim[i][0]])
        if lst_xTicks[i] is not None:
            axs.set_xticks(lst_xTicks[i])
        if lst_yTicks[i] is not None:
            axs.set_yticks(lst_yTicks[i])
        if lst_xTickLabels[i] is not None:
            axs.set_xticklabels(lst_xTickLabels[i])
        if lst_yTickLabels[i] is not None:
            axs.set_yticklabels(lst_yTickLabels[i])
        if annotate_lttrs is not False:
            if annotate_lttrs_loc == 'down':
                y_lttrs = (py/px*.02)
            elif annotate_lttrs_loc == 'up':
                y_lttrs = 1 - py
            if n_ax == 1:  # if only one plot is given, do not put the letters
                axs.annotate('(' + annotate_lttrs + ')',
                              xycoords='axes fraction',
                              xy=(0, 0), rotation=0, size='large',
                              xytext=(0, y_lttrs), weight='bold')
            elif n_ax > 1:  # if only one plot is given, do not put the letters
                axs.annotate('(' + lttrs[i] + ')', xycoords='axes fraction',
                              xy=(0, 0), rotation=0, size='large',
                              xytext=(0, y_lttrs), weight='bold')
    # if secondary (twin) axis are given, set thier properties
    if lst_axt is not None:
        for i, axst in enumerate(lst_axt):
            axst.grid(False)  # grid is always false on secondaty axis
            # for each property, if the variable is not false, it is set
            if lst_ytLab[i] is not None:
                axst.set_ylabel(lst_ytLab[i])
            if lst_ytLim[i] is not None:
                axst.set_ylim([lst_ytLim[i][0]*(1 + py) - py*lst_ytLim[i][1],
                              lst_ytLim[i][1]*(1 + py) - py*lst_ytLim[i][0]])
            if lst_ytTicks[i] is not None:
                axst.set_yticks(lst_ytTicks[i])
            if lst_ytTickLabels[i] is not None:
                axst.set_yticklabels(lst_ytTickLabels[i])
    # create a legend merging the entries for each couple of ax and axt
    if any(lst_legend):
        if lst_axt is None:  # with no axt, only axs in ax needs a legend
            for i, axs in enumerate(lst_ax):
                axs.legend(loc=lst_legend[i], ncol=ncol_leg)
        else:  # merge the legend for each couple of ax and axt
            i = 0
            for axs, axst in zip(lst_ax, lst_axt):
                hnd_ax, lab_ax = axs.get_legend_handles_labels()
                hnd_axt, lab_axt = axst.get_legend_handles_labels()
                axs.legend(hnd_ax + hnd_axt, lab_ax + lab_axt, loc=lst_legend[i],
                           ncol=ncol_leg)
                i += 1
    try:
        fig.align_labels()  # align labels of subplots, needed only for multi plot
    except AttributeError:
        print('align_labels not performed')
    # if a subfolder is specified, create the subfolder inside the output
    # folder if not already there and save the figure in it
    if subfolder is not None:
        out_path = plib.Path(out_path, subfolder)  # update out_path
        plib.Path(out_path).mkdir(parents=True, exist_ok=True)  # check if
        # folder is there, if not create it
    # set figure margins and save the figure in the output folder
    if tight_layout is False:  # if margins are given sets margins and save
        fig.subplots_adjust(*fig_adj_par[0:6])  # set margins
        plt.savefig(plib.Path(out_path, fig_name + '.png'), dpi=300,
                    transparent=transparency)
        if pdf is not False:  # save also as pdf
            plt.savefig(plib.Path(out_path, fig_name + '.pdf'))
        if svg is not False:  # save also as pdf
            plt.savefig(plib.Path(out_path, fig_name + '.svg'))
    else:  # margins are not given, use a tight layout option and save
        plt.savefig(plib.Path(out_path, fig_name + '.png'),
                    bbox_inches="tight", dpi=300, transparent=transparency)
        if pdf is not False:  # save also as pdf
            plt.savefig(plib.Path(out_path, fig_name + '.pdf'),
                        bbox_inches="tight")
        if svg is not False:  # save also as pdf
            plt.savefig(plib.Path(out_path, fig_name + '.svg'),
                        bbox_inches="tight")
    # add the title after saving, so it's only visible in the console
    if title is True:
        lst_ax[0].annotate(fig_name, xycoords='axes fraction', size='small',
                            xy=(0, 0), xytext=(0.05, .95), clip_on=True)
    return