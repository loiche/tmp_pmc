import re
import spacy
import collections
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def make_bar_plot(target_col_name, figsize_=(20, 3)):
    plt.figure(figsize=figsize_)
    plt.bar(REFOLDdb_df.loc[:, target_col_name].value_counts().index, 
            REFOLDdb_df.loc[:, target_col_name].value_counts().values,
           )
    plt.xticks(rotation=-90)
    plt.title(target_col_name)
    plt.show()


REFOLDdb_df = pd.read_csv("../../data/REFOLDdb/REFOLDdb_summary.csv", index_col=0)
REFOLDdb_df.dtypes

REFOLDdb_index_category_dict = {"Protein": ['Protein Name', 'Abbreviated Name', 'SCOP Family', 'Structure Notes', 'Organism', 'UniProt Accession', 'SCOP Unique ID', 'Structure Solved','Class', 'Molecularity'],
                                "Construct": ['Full Length', 'Domain', 'Chimera', 'Variants', 'Chain Length', 'Molecular Weight', 'Pi', 'Disulphides', 'Full Sequence', 'Notes'],
                               "Expression": ['Report', 'Project Aim', 'Fusion','Protein Expression and Production', 'Expression Host','Expression Strain', 'Expression Temp', 'Expression Time','Expression Vector', 'Expression Protocol', 'Method of Induction','Cell Density at Induction', 'Cell Disruption Method', 'Lytic Agent','Pre-Refolding Purification', 'Solubility'],
                               "Refolding":['Refolding Method', 'Wash Buffer', 'Solubilization Buffer', 'Refolding Buffer', 'Tag Cleaved', 'Refolding pH', 'Refolding Temperature', 'Protein Concentration', 'Refolding Time', 'Redox Agent', 'Redox Agent Concentration', 'Refolding Protocol', 'Refolding Assay', 'Refolding Chaperones', 'Refolding Additives', 'Additives Concentration', 'Refolding Yield', 'Purity']}

# construct
#- 数値っぽいもの
#    - Chain Length
#    - Molecular Weight

# Expression
#- 数値っぽいもの
#    - "Expression Temp": ok
#    - "Expression Time": Expression Time min, Expression Time maxに分解
#    - "Cell Density at Induction": 直す必要なし？
    
# "Expression Temp"
REFOLDdb_df.loc[:, "Expression Temp"]
REFOLDdb_df.loc[:, "Expression Temp"].describe()

### "Expression Time"
#5-6daysが気になる
REFOLDdb_df.loc[:, "Expression Time"].describe()
make_bar_plot("Expression Time")

def fix_expression_time(target_value):
    flag_bar = False
    flag_min = False
    
    if type(target_value)==str:
        #特殊系
        if "16-h" == target_value:
            return 16, 16
        if "not stated" == target_value or "unknown" == target_value or "Unknown" == target_value or "Not stated" == target_value or "?" == target_value or "ON" == target_value:
            return False, False
        if "overnight" == target_value or "4h/overnight" == target_value:
            return 12, 12
        if "3g" == target_value:
            return 3, 3
        if "up to 16h" == target_value:
            return 16, 16
        if "5-6 days" == target_value:
            return 24 * 5, 24 * 6
        
        # 空白の削除
        target_value = target_value.replace(" ", "")
        # hoursの削除
        target_value = target_value.replace("hours", "")
        # hrsの削除
        target_value = target_value.replace("hrs", "")
        # hrの削除
        target_value = target_value.replace("hr", "")
        # hの削除
        target_value = target_value.replace("h", "")
        # H
        target_value = target_value.replace("H", "")
        
        #削除
        if "-" in target_value:
            flag_bar = True
        if "min" in target_value:
            flag_min = True
        
        if flag_bar == 1 and flag_min == 0:
        # -への対応
            target_value_min, target_value_max = target_value.split("-")
            
        # minへの対応
        if flag_bar == 0 and flag_min == 1:
            target_value = target_value.replace("min", "")
            try:
                target_value = float(target_value) / 60
                target_value_min, target_value_max = target_value, target_value
            except:
                raise ValueError("minが入っている文字列に数値以外の文字列も入っています。")
                
        if flag_bar == 0 and flag_min == 0:
            target_value_min, target_value_max = target_value, target_value
                
    else:
        target_value_min, target_value_max = target_value, target_value
    
    return target_value_min, target_value_max

# expression timeのテスト
test_expression_time_list = REFOLDdb_df["Expression Time"].value_counts().keys()
print('{0:18} | {1:18} | {2:18}'.format("original", "target_value_min", "target_value_max"))
for _, target_value in enumerate(test_expression_time_list):
    target_value_min, target_value_max = fix_expression_time(target_value)
    try:
        float(target_value_min), float(target_value_max)
    except:
        print('{0:18} | {1:18} | {2:18}'.format(target_value, target_value_min, target_value_max), _)

# 記録
expression_time_min_list = []
expression_time_max_list = []
for key, value in REFOLDdb_df.loc[:, "Expression Time"].items():
    target_value_min, target_value_max = fix_expression_time(value)
    if target_value_min and target_value_max:
        expression_time_min_list.append(float(target_value_min))
        expression_time_max_list.append(float(target_value_max))
    else:
        expression_time_min_list.append(target_value_min)
        expression_time_max_list.append(target_value_max)

expression_time_min = "Expression Time min"
expression_time_max = "Expression Time max"
REFOLDdb_df[expression_time_min] = expression_time_min_list
REFOLDdb_df[expression_time_max] = expression_time_max_list

plt.hist(REFOLDdb_df.loc[:, expression_time_min], bins=100);
plt.hist(REFOLDdb_df.loc[:, expression_time_max], bins=100);

make_bar_plot("Cell Density at Induction")

def fix_cell_density_at_induction(target_value):
    if type(target_value)==str and len(target_value.split("=")) == 2:
        target_value = target_value.replace("\n", "").replace("\t", "")
        OD_set, OD_value = target_value.split("=")
        
    else:
        pass
    
    return OD_set, OD_value

# Refolding
#- 数値っぽいもの
#    - Refolding Temperature
#    - Protein Concentration
#    - Refolding Time
#    - Refolding Yield
#    - Purity 

REFOLDdb_df.loc[:, REFOLDdb_index_category_dict["Refolding"]]

plt.hist(REFOLDdb_df.loc[:, "Refolding Temperature"], bins=10);

REFOLDdb_df.loc[:, "Refolding Temperature"].describe()

#外れ値の確認
###正しかった
REFOLDdb_df[REFOLDdb_df.loc[:, "Refolding Temperature"] > 50].loc[227, "Refolding Protocol"]

#正しい値へ変更
REFOLDdb_df[REFOLDdb_df.loc[:, "Refolding Temperature"] > 50].loc[100, "Refolding Protocol"]
REFOLDdb_df.loc[100, "Refolding Temperature"] = 4

#削除
REFOLDdb_df[REFOLDdb_df.loc[:, "Refolding Temperature"] > 50].loc[851, "Refolding Protocol"]
REFOLDdb_df = REFOLDdb_df[REFOLDdb_df["Protein Name"] != "prion-like yeast protein"]

REFOLDdb_df.loc[:, "Refolding Temperature"].describe()
plt.hist(REFOLDdb_df.loc[:, "Refolding Temperature"], bins=60);

#Protein Concentration
make_bar_plot("Protein Concentration")
def extract_value(target_value):
    target_value = target_value.replace("&#8722;", "")
    target_value = target_value.replace("&#956;", "")
    extracted_value_list = re.findall("\d+(?:\.\d+)?", target_value)
    if type(extracted_value_list) == list and len(extracted_value_list) == 2:
        max_protein_concentration_value , min_protein_concentration_value = extracted_value_list[1], extracted_value_list[0]
    elif len(extracted_value_list) == 1:
        max_protein_concentration_value , min_protein_concentration_value = extracted_value_list[0], extracted_value_list[0]
    else:
        max_protein_concentration_value , min_protein_concentration_value = np.nan, np.nan
        
    return max_protein_concentration_value , min_protein_concentration_value

max_protein_concentration_value_list = []
min_protein_concentration_value_list = []
for _ in REFOLDdb_df.loc[:, "Protein Concentration"].values:
    if type(_) == str:
        max_protein_concentration_value , min_protein_concentration_value = extract_value(_)
        max_protein_concentration_value_list.append(float(max_protein_concentration_value))
        min_protein_concentration_value_list.append(float(min_protein_concentration_value))
    else:
        max_protein_concentration_value_list.append(_)
        min_protein_concentration_value_list.append(_)

REFOLDdb_df.loc[:, "Protein_Concentration_max"] = max_protein_concentration_value_list
REFOLDdb_df.loc[:, "Protein_Concentration_min"] = min_protein_concentration_value_list

def extract_unit(target_value):
    target_value = target_value.replace("&#", "")
    target_value = target_value.replace(";", "")
    target_value = target_value.replace(" ", "")
    target_value = target_value.replace(".", "")
    target_value = target_value.replace("-", "")
    unit_list = re.findall("\D+", target_value)
    target_unit = "_".join(unit_list)
    
    unit_convert_dict = {"mg/ml": "mg/ml",
                         "microM": "mM",
                         "microg/ml":"mg/ml",
                         "micrograms/ml": "mg/ml",
                         "micromolar": "mM",
                         "mgml": "mg/ml",
                         "upto_mg/ml": "mg/ml",
                         "lessthan_mg/ml": "mg/ml",
                         "<_mg/ml": "mg/ml",
                         "mg/mlorless": "mg/ml",
                         "<_micrograms/ml": "mg/ml",
                         "/mg/ml": "mg/ml",
                         'micrograms/mL': "mg/ml",
                         'nd': np.nan,
                         "h": np.nan,
                         "uMfinal": "uM",
                         '>_mg/L': "ug/ml",
                         'g/L': "mg/ml",
                         "micrograms/m": "mg/ml",
                         "microgram/ml": "mg/ml",
                         "dilute": np.nan,
                         'mg/mL': "mg/ml",
                         '/_mgperml': "mg/ml",
                         ',_µg/ml': "mg/ml",
                         "mgper_mlrefoldingbuffer": "mg/ml",
                         "overnight": np.nan,
                         'microg/l': "mg/ml",
                         'µg/mL': "ug/ml",
                         "<_ug/ml": "ug/ml",
                         'mg·mL': "mg/ml",
                         'mg/L': "ug/ml",
                         'min': np.nan,
                         '–_M': "M",
                         'mgto_mg/ml': "mg/ml",
                         "": np.nan,
                         "µM": "uM",
                         "notstated": np.nan,
                         'mg': "mg",
                         'µg': "ug",
                         'A_/ml': np.nan,
                         '%': '%',
                         "ug/ml": "ug/ml",
                         "uM": "uM",
                         "mM": "mM"
                        }
    
    if target_unit in unit_convert_dict.keys():
        target_unit = unit_convert_dict[target_unit]
    else:
        raise ValueError()
        
    return target_unit

protein_concentration_unit_list = []
for _ in REFOLDdb_df.loc[:, "Protein Concentration"].values:
    if type(_) == str:
        protein_concentration_unit = extract_unit(_)
        protein_concentration_unit_list.append(protein_concentration_unit)
    else:
        protein_concentration_unit_list.append(_)

collections.Counter(protein_concentration_unit_list)
REFOLDdb_df.loc[:, "Protein_Concentration_unit"] = protein_concentration_unit_list
make_bar_plot("Protein_Concentration_min")
make_bar_plot("Protein_Concentration_max")
make_bar_plot("Protein_Concentration_unit")
REFOLDdb_df["Protein_Concentration_unit"].value_counts()

# Refolding Time
make_bar_plot("Refolding Time")

def extract_refolding_time(target_value):
    #target_value = target_value.replace("&#8722;", "")
    #target_value = target_value.replace("&#956;", "")
    extracted_value_list = re.findall("\d+(?:\.\d+)?", target_value)
    if type(extracted_value_list) == list and len(extracted_value_list) == 2:
        max_refolding_time_value , min_refolding_time_value = extracted_value_list[1], extracted_value_list[0]
    elif len(extracted_value_list) == 1:
        max_refolding_time_value , min_refolding_time_value = extracted_value_list[0], extracted_value_list[0]
    else:
        max_refolding_time_value , min_refolding_time_value = np.nan, np.nan
        
    return max_refolding_time_value, min_refolding_time_value

max_refolding_time_value_list = []
min_refolding_time_value_list = []
for _ in REFOLDdb_df.loc[:, "Refolding Time"].values:
    if type(_) == str:
        max_refolding_time_value , min_refolding_time_value = extract_value(_)
        max_refolding_time_value_list.append(float(max_refolding_time_value))
        min_refolding_time_value_list.append(float(min_refolding_time_value))
        #print(max_refolding_time_value, type(max_refolding_time_value), min_refolding_time_value, type(min_refolding_time_value))
    else:
        max_refolding_time_value_list.append(_)
        min_refolding_time_value_list.append(_)

REFOLDdb_df.loc[:, "refolding_time_max"] = max_refolding_time_value_list
REFOLDdb_df.loc[:, "refolding_time_min"] = min_refolding_time_value_list

def extract_refolding_time_unit(target_value):
    target_value = target_value.replace("&#", "")
    target_value = target_value.replace(";", "")
    target_value = target_value.replace(" ", "")
    target_value = target_value.replace(".", "")
    target_value = target_value.replace("-", "")
    unit_list = re.findall("\D+", target_value)
    target_unit = "_".join(unit_list)
    
    unit_convert_dict = {"weeks": "week", 
                         "h": "h",
                         "": np.nan,
                         "min": "min",
                         "mins": "min",
                         "atleast_h": "h",
                         "overnight": "overnight",
                         ">_h": "h",
                         "week": "week",
                         "min_h": "h",
                         "days": "day",
                         "hr": "h",
                         ">_hr": "h",
                         "hrs": "h",
                         "~_hours": "h",
                         "hours": "h",
                         "hs": "h",
                         "d": "day",
                         "hour": "h",
                         "overnig": "overnight",
                         "h+": "h",
                         "o/night": "overnight",
                         "Overnight": "overnight",
                         "x_h": "h",
                         "notstated": np.nan,
                         "day": "day",
                         "m": "min",
                         "heachstep": "h",
                         "onehour": "onehour"
                        }
    
    if target_unit in unit_convert_dict.keys():
        target_unit = unit_convert_dict[target_unit]
    else:
        print(target_unit)
        raise ValueError()
        
    return target_unit
refolding_time_unit_list = []
for _ in REFOLDdb_df.loc[:, "Refolding Time"].values:
    if type(_) == str:
        refolding_time_unit = extract_refolding_time_unit(_)
        refolding_time_unit_list.append(refolding_time_unit)
    else:
        refolding_time_unit_list.append(_)

collections.Counter(refolding_time_unit_list)
REFOLDdb_df.loc[:, "Refolding_Time_unit"] = refolding_time_unit_list

make_bar_plot("Refolding_Time_unit")
plt.hist(REFOLDdb_df["refolding_time_max"].values)
plt.hist(REFOLDdb_df["refolding_time_min"].values)

make_bar_plot("Refolding Yield")

def extract_refolding_yield(target_value):
    extracted_value_list = re.findall("\d+(?:\.\d+)?", target_value)
    if type(extracted_value_list) == list and len(extracted_value_list) == 2:
        max_refolding_yield_value , min_refolding_yield_value = extracted_value_list[1], extracted_value_list[0]
    elif len(extracted_value_list) == 1:
        max_refolding_yield_value , min_refolding_yield_value = extracted_value_list[0], extracted_value_list[0]
    else:
        max_refolding_yield_value , min_refolding_yield_value = np.nan, np.nan
        
    return max_refolding_yield_value, min_refolding_yield_value

max_refolding_yield_value_list = []
min_refolding_yield_value_list = []
for _ in REFOLDdb_df.loc[:, "Refolding Yield"].values:
    if type(_) == str:
        max_refolding_yield_value , min_refolding_yield_value = extract_refolding_yield(_)
        max_refolding_yield_value_list.append(float(max_refolding_yield_value))
        min_refolding_yield_value_list.append(float(min_refolding_yield_value))
    else:
        max_refolding_yield_value_list.append(_)
        min_refolding_yield_value_list.append(_)

REFOLDdb_df.loc[:, "refolding_yield_max"] = max_refolding_yield_value_list
REFOLDdb_df.loc[:, "refolding_yield_min"] = min_refolding_yield_value_list

def extract_refolding_yield_unit(target_value):
    target_value = target_value.replace("&#", "")
    target_value = target_value.replace(";", "")
    target_value = target_value.replace(" ", "")
    target_value = target_value.replace(".", "")
    target_value = target_value.replace("-", "")
    unit_list = re.findall("\D+", target_value)
    target_unit = "_".join(unit_list)
    
    unit_convert_dict = {'mgfrom_mlculture': 'mgfrom_mlculture',
                         "%": "%",
                         "mg/_Lculture": "mg/l",
                         "mg/_Lmedia": "mg/l",
                         "mg/Lcultureapprox_%(w/w)ofinclusionbody": '9mg/L culture; approx. 10% (w/w) of inclusion body',
                         "%_mg/Lculture": "%_mg/Lculture",
                         "mg/Lcultur": "mg/l",
                         "mg/Lculture": "mg/l",
                         ">_mg/Lmedia": "mg/l",
                         'mg/Lculture(fermentation),_mg/L(shakeflask)': '25mg/L culture (fermentation), 2mg/L (shake-flask)',
                         "milligramquantities": np.nan,
                         "mg/L": "mg/l",
                         "mgfrom_mgoffrozencells": '0.4mg from 0.5 mg of frozen cells',
                         "mg/Lculture_%": '37mg/L culture - 36%',
                         "<_%": "%",
                         ">_%": "%",
                         "mgperlitrecellculture": '6 mg per litre cell culture',
                         "microgramspermilligram": '23 micrograms per milligram',
                         "mg": "mg",
                         "Notstated": np.nan,
                         "g/L": "g/l",
                         "mgperL": "mg/l",
                         "mg/Lrefold": "mg/l",
                         "mg/LinLBbroth": "mg/l",
                         "mg/ml": "g/l",
                         "nearly_%": "%",
                         "%recovery": "%",
                         "mg/_mgIBs": '35mg/100mg IBs',
                         ">_mg/L": "mg/l",
                         "mg/_L": '0.5mg/4L',
                         "": "%",
                         "mg/gpackedcells": '0.5mg/g packed cells',
                         "mg/mLcellculture": "mg/ml", 
                         "%(ofmass)": "%",
                         "~_mg/litreculture": "mg",
                         "approx_%": "%",
                         "~_%": "%",
                         ">_mggp_fusionprotein/gEcoliwetweight": '> 20 mg gp41 fusion protein/g E. coli wet weight',
                         "ug": "ug",
                         "mg/litreofculture": "mg/l",
                         "mg/Lcells": "mg/l",
                         "microg/_mlculture": "mg/ml",
                         "U/_mlculture": '0.1U/100ml culture',
                         "microg/ml": "mg/ml",
                         "mgprotein/L": "mg/l",
                         "upto_%": "%",
                         "nm": "nm",
                         "mg/_gcellss": '3mg/10g cellss',
                         "mg/l": "mg/l",
                         "mg/gIB": "mg/g",
                         "%,_mgperL": '20%, 20mg per L',
                         "mg/gcell": "mg/g",
                         "g/Lculture": "mg/ml",
                         "mgperlit": "mg/l",
                         "mg/litre": "mg/l",
                         "mg(_%)": '15 mg (30%)',
                         "U/mg": "mg/g",
                         "mgprot/_mg": '80mg prot/110mg',
                         "mgperliter": "mg/l", 
                         "%(±_%)": '55% (±6%)',
                         "to_mg/l": '0.75 to 1 mg/l',
                         ",_%": '2,0%',
                         "ug/l": "ug/l",
                         "%_%": "%",
                         "to_mg/L": 'mg/l',
                         "to_mg/liter": "mg/l",
                         "–_mg/L": "mg/l",
                         "mg/lculture": "mg/l",
                         "mg/_g": "mg/g",
                         "mg/Lcul": "mg/l",
                         "~_mg": "mg",
                         "–_mgper": ' 5–10 mg per 1',
                         "mgper_ml": '50 mg per 500ml',
                         "noactivity": np.nan
                        }
    
    if target_unit in unit_convert_dict.keys():
        target_unit = unit_convert_dict[target_unit]
    else:
        print(target_unit)
        raise ValueError()
        
    return target_unit

refolding_yield_unit_list = []
for _ in REFOLDdb_df.loc[:, "Refolding Yield"].values:
    if type(_) == str:
        refolding_yield_unit = extract_refolding_yield_unit(_)
        refolding_yield_unit_list.append(refolding_yield_unit)
    else:
        refolding_yield_unit_list.append(_)

collections.Counter(refolding_yield_unit_list)
REFOLDdb_df.loc[:, "refolding_yield_unit"] = refolding_yield_unit_list

make_bar_plot("refolding_yield_unit")
plt.hist(REFOLDdb_df["refolding_yield_max"], bins=100)
plt.hist(REFOLDdb_df["refolding_yield_min"], bins=100)

# purity
make_bar_plot("Purity")
def extract_purity(target_value):
    extracted_value_list = re.findall("\d+(?:\.\d+)?", target_value)
    if type(extracted_value_list) == list and len(extracted_value_list) == 2:
        max_purity_value , min_purity_value = extracted_value_list[1], extracted_value_list[0]
    elif len(extracted_value_list) == 1:
        max_purity_value , min_purity_value = extracted_value_list[0], extracted_value_list[0]
    else:
        max_purity_value , min_purity_value = np.nan, np.nan
        
    return max_purity_value , min_purity_value

max_purity_value_list = []
min_purity_value_list = []
for _ in REFOLDdb_df.loc[:, "Purity"].values:
    if type(_) == str:
        max_purity_value , min_purity_value = extract_purity(_)
        max_purity_value_list.append(float(max_purity_value))
        min_purity_value_list.append(float(min_purity_value))
    else:
        max_purity_value_list.append(_)
        min_purity_value_list.append(_)
        
REFOLDdb_df.loc[:, "purity_max"] = max_purity_value_list
REFOLDdb_df.loc[:, "purity_min"] = min_purity_value_list

def extract_purity_unit(target_value):
    target_value = target_value.replace("&#", "")
    target_value = target_value.replace(";", "")
    target_value = target_value.replace(" ", "")
    target_value = target_value.replace(".", "")
    target_value = target_value.replace("-", "")
    unit_list = re.findall("\D+", target_value)
    target_unit = "_".join(unit_list)
    
    unit_convert_dict = {"seegelattached": np.nan,
                         ">_%": "%",
                         "Notstated": np.nan,
                         "%prochymosin,asjudgedbySDS?PAGE": "%",
                         "Seegel": np.nan,
                         "pure": "pure",
                         "%": "%",
                         'afewc': np.nan,
                         "morethan_%": "%",
                         "singlebandonSDSgel": 'single band on SDS gel',
                         "manyvisiblebandsonSDS": np.nan,
                         "to_veryfaintimpuritybandsonSDSPAGE": np.nan,
                         "singlebandonSDSPAGE": 'single band on SDS PAGE',
                         "somecontaminatingbandsonSDS": np.nan,
                         "fewveryfaintbandsonSDSPAGE": np.nan,
                         "fewfeintbandsonSDSPAGE": np.nan,
                         "manystrongbandsonSDSPAGE": np.nan,
                         "manycontaminatingbandsonSDSPAGE": np.nan,
                         "homogeneous": np.nan,
                         "onecontaminatingbandonSDSPAGE": np.nan,
                         "fewcontaminatingbandsonSDSPAGE": np.nan,
                         "%pure": "%",
                         "nd": np.nan,
                         "NotStated": np.nan,
                         "": "%",
                         "~_%": "%",
                         ">_%(CoomassiestainedSDSPAGE)": "%",
                         ">_%aftercolumnchromatography": "%",
                         "%+": "%",
                         "mixture": np.nan,
                         "NMRgrade": np.nan,
                         "Homogen": np.nan,
                         "High": np.nan,
                         ">": "%",
                         "±": ' 88 ± 3',
                         "mg/": "mg/",
                         "mg/l": "mg/l",
                         "mg/L": "mg/L",
                         "mg": "mg"
                        }
    
    if target_unit in unit_convert_dict.keys():
        target_unit = unit_convert_dict[target_unit]
    else:
        print(target_unit)
        raise ValueError()
        
    return target_unit

purity_unit_list = []
for _ in REFOLDdb_df.loc[:, "Purity"].values:
    if type(_) == str:
        purity_unit = extract_purity_unit(_)
        purity_unit_list.append(purity_unit)
    else:
        purity_unit_list.append(_)
REFOLDdb_df.loc[:, "purity_unit"] = purity_unit_list

collections.Counter(purity_unit_list)

make_bar_plot("purity_unit")

plt.hist(REFOLDdb_df["purity_max"].values)
plt.hist(REFOLDdb_df["purity_min"].values)

numeric_cols_list = ["Chain Length", "Molecular Weight", "Pi", "Expression Temp", "Refolding pH", "Refolding Temperature"]
REFOLDdb_df.loc[:, numeric_cols_list].describe()