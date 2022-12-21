import re
import spacy
import collections
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def make_bar_plot(target_df, target_col_name, figsize_=(20, 3)):
    plt.figure(figsize=figsize_)
    plt.bar(target_df.loc[:, target_col_name].value_counts().index, 
            target_df.loc[:, target_col_name].value_counts().values,
           )
    plt.xticks(rotation=-90)
    plt.title(target_col_name)
    plt.show()

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

def extract_protein_concentration_value(target_value):
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

def extract_protein_concentration_unit(target_value):
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

def annotate_protein_concentration_df(target_df):
    """
    target_df: REFOLD_db
    """
    protein_concentration_dict = {"Protein_Concentration_max_mg/ml": [], 
                                  "Protein_Concentration_min_mg/ml": [], 
                                  "Protein_Concentration_max_mM": [], 
                                  "Protein_Concentration_min_mM": [], 
                                  }
    for key, value in target_df.loc[: ,["Protein_Concentration_max", "Protein_Concentration_min", "Protein_Concentration_unit"]].iterrows():
        if value["Protein_Concentration_unit"] == "mg/ml":
            protein_concentration_dict["Protein_Concentration_max_mg/ml"].append(value["Protein_Concentration_max"])
            protein_concentration_dict["Protein_Concentration_min_mg/ml"].append(value["Protein_Concentration_min"])
            protein_concentration_dict["Protein_Concentration_max_mM"].append(np.nan)
            protein_concentration_dict["Protein_Concentration_min_mM"].append(np.nan)
        elif value["Protein_Concentration_unit"] == "mM":
            protein_concentration_dict["Protein_Concentration_max_mg/ml"].append(np.nan)
            protein_concentration_dict["Protein_Concentration_min_mg/ml"].append(np.nan)
            protein_concentration_dict["Protein_Concentration_max_mM"].append(value["Protein_Concentration_max"])
            protein_concentration_dict["Protein_Concentration_min_mM"].append(value["Protein_Concentration_min"])
        elif value["Protein_Concentration_unit"] == "uM":
            protein_concentration_dict["Protein_Concentration_max_mg/ml"].append(np.nan)
            protein_concentration_dict["Protein_Concentration_min_mg/ml"].append(np.nan)
            protein_concentration_dict["Protein_Concentration_max_mM"].append(value["Protein_Concentration_max"] / 1000)
            protein_concentration_dict["Protein_Concentration_min_mM"].append(value["Protein_Concentration_min"] / 1000)
        elif value["Protein_Concentration_unit"] == "ug/ml":
            protein_concentration_dict["Protein_Concentration_max_mg/ml"].append(value["Protein_Concentration_max"] / 1000)
            protein_concentration_dict["Protein_Concentration_min_mg/ml"].append(value["Protein_Concentration_min"] / 1000)
            protein_concentration_dict["Protein_Concentration_max_mM"].append(np.nan)
            protein_concentration_dict["Protein_Concentration_min_mM"].append(np.nan)
        elif value["Protein_Concentration_unit"] == "M":
            protein_concentration_dict["Protein_Concentration_max_mg/ml"].append(np.nan)
            protein_concentration_dict["Protein_Concentration_min_mg/ml"].append(np.nan)
            protein_concentration_dict["Protein_Concentration_max_mM"].append(value["Protein_Concentration_max"] * 1000)
            protein_concentration_dict["Protein_Concentration_min_mM"].append(value["Protein_Concentration_min"] * 1000)
        else:
            protein_concentration_dict["Protein_Concentration_max_mg/ml"].append(np.nan)
            protein_concentration_dict["Protein_Concentration_min_mg/ml"].append(np.nan)
            protein_concentration_dict["Protein_Concentration_max_mM"].append(np.nan)
            protein_concentration_dict["Protein_Concentration_min_mM"].append(np.nan)
            if value["Protein_Concentration_unit"] in target_df["Protein_Concentration_unit"].value_counts().keys():
                print(value)
                
    protein_concentration_df = pd.DataFrame(protein_concentration_dict, index=target_df.index)
    target_df = pd.concat([target_df, protein_concentration_df], axis=1)
    return target_df

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

def annotate_refolding_time_df(target_df):
    """
    target_df: REFOLD_db
    """
    
    refolding_time_dict = {"refolding_time_max_h": [], 
                           "refolding_time_min_h": [], 
                          }

    for key, value in target_df.loc[: ,["refolding_time_max", "refolding_time_min", "Refolding_Time_unit"]].iterrows():
        if value["Refolding_Time_unit"] == "h":
            refolding_time_dict["refolding_time_max_h"].append(value["refolding_time_max"])
            refolding_time_dict["refolding_time_min_h"].append(value["refolding_time_min"])
        elif value["Refolding_Time_unit"] == "week":
            refolding_time_dict["refolding_time_max_h"].append(value["refolding_time_max"] * 168)
            refolding_time_dict["refolding_time_min_h"].append(value["refolding_time_min"] * 168)
        elif value["Refolding_Time_unit"] == "min":
            refolding_time_dict["refolding_time_max_h"].append(value["refolding_time_max"] / 60)
            refolding_time_dict["refolding_time_min_h"].append(value["refolding_time_min"]/ 60)
        elif value["Refolding_Time_unit"] == "overnight":
            refolding_time_dict["refolding_time_max_h"].append(12)
            refolding_time_dict["refolding_time_min_h"].append(12)
        elif value["Refolding_Time_unit"] == "day":
            refolding_time_dict["refolding_time_max_h"].append(value["refolding_time_max"] * 24)
            refolding_time_dict["refolding_time_min_h"].append(value["refolding_time_min"] * 24)
        elif value["Refolding_Time_unit"] == "onehour":
            refolding_time_dict["refolding_time_max_h"].append(1)
            refolding_time_dict["refolding_time_min_h"].append(1)
        else:
            refolding_time_dict["refolding_time_max_h"].append(np.nan)
            refolding_time_dict["refolding_time_min_h"].append(np.nan)
            if value["Refolding_Time_unit"] in target_df["Refolding_Time_unit"].value_counts().keys():
                print(value)
                
    refolding_time_df = pd.DataFrame(refolding_time_dict, index=target_df.index)
    target_df = pd.concat([target_df, refolding_time_df], axis=1)
    return target_df

def extract_refolding_yield(target_value):
    extracted_value_list = re.findall("\d+(?:\.\d+)?", target_value)
    if type(extracted_value_list) == list and len(extracted_value_list) == 2:
        max_refolding_yield_value , min_refolding_yield_value = extracted_value_list[1], extracted_value_list[0]
    elif len(extracted_value_list) == 1:
        max_refolding_yield_value , min_refolding_yield_value = extracted_value_list[0], extracted_value_list[0]
    else:
        max_refolding_yield_value , min_refolding_yield_value = np.nan, np.nan
        
    return max_refolding_yield_value, min_refolding_yield_value

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

def annotate_refolding_yield_df(target_df):
    """
    target_df: REFOLD_db
    """
    
    
    refolding_yield_dict = {"refolding_yield_max_%": [], 
                       "refolding_yield_min_%": [], 
                        "refolding_yield_max_mg/l": [], 
                       "refolding_yield_min_mg/l": [], 
                        "refolding_yield_max_mg": [], 
                       "refolding_yield_min_mg": [], 
                      }

    for key, value in target_df.loc[: ,["refolding_yield_max", "refolding_yield_min", "refolding_yield_unit"]].iterrows():
        if value["refolding_yield_unit"] == "%":
            refolding_yield_dict["refolding_yield_max_%"].append(value["refolding_yield_max"])
            refolding_yield_dict["refolding_yield_min_%"].append(value["refolding_yield_min"])
            refolding_yield_dict["refolding_yield_max_mg/l"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg/l"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
        elif value["refolding_yield_unit"] == "mg/l":
            refolding_yield_dict["refolding_yield_max_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg/l"].append(value["refolding_yield_max"])
            refolding_yield_dict["refolding_yield_min_mg/l"].append(value["refolding_yield_min"])
            refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
        elif value["refolding_yield_unit"] == "mg":
            refolding_yield_dict["refolding_yield_max_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg/l"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg/l"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg"].append(value["refolding_yield_max"])
            refolding_yield_dict["refolding_yield_min_mg"].append(value["refolding_yield_min"])
        elif value["refolding_yield_unit"] == "g/l":
            refolding_yield_dict["refolding_yield_max_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg/l"].append(value["refolding_yield_max"] * 1000)
            refolding_yield_dict["refolding_yield_min_mg/l"].append(value["refolding_yield_min"] * 1000)
            refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
        elif value["refolding_yield_unit"] == "mg/g":
            refolding_yield_dict["refolding_yield_max_%"].append(value["refolding_yield_max"] / 1000)
            refolding_yield_dict["refolding_yield_min_%"].append(value["refolding_yield_min"] / 1000)
            refolding_yield_dict["refolding_yield_max_mg/l"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg/l"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
        elif value["refolding_yield_unit"] == "mg/ml":
            refolding_yield_dict["refolding_yield_max_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg/l"].append(value["refolding_yield_max"] * 1000)
            refolding_yield_dict["refolding_yield_min_mg/l"].append(value["refolding_yield_min"] * 1000)
            refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
        elif value["refolding_yield_unit"] == "mgfrom_mlculture":
            refolding_yield_dict["refolding_yield_max_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg/l"].append(value["refolding_yield_max"])
            refolding_yield_dict["refolding_yield_min_mg/l"].append(value["refolding_yield_min"])
            refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
        elif value["refolding_yield_unit"] == '9mg/L culture; approx. 10% (w/w) of inclusion body':
            refolding_yield_dict["refolding_yield_max_%"].append(10)
            refolding_yield_dict["refolding_yield_min_%"].append(10)
            refolding_yield_dict["refolding_yield_max_mg/l"].append(9)
            refolding_yield_dict["refolding_yield_min_mg/l"].append(9)
            refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
        elif value["refolding_yield_unit"] == "%_mg/Lculture":
            if target_df.loc[key, "Refolding Yield"] == "5.7% - 4-5mg/L culture":
                refolding_yield_dict["refolding_yield_max_%"].append(5.7)
                refolding_yield_dict["refolding_yield_min_%"].append(5.7)
                refolding_yield_dict["refolding_yield_max_mg/l"].append(5)
                refolding_yield_dict["refolding_yield_min_mg/l"].append(4)
                refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
                refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
            elif target_df.loc[key, "Refolding Yield"] == "50%- 6mg/L culture":
                refolding_yield_dict["refolding_yield_max_%"].append(50)
                refolding_yield_dict["refolding_yield_min_%"].append(50)
                refolding_yield_dict["refolding_yield_max_mg/l"].append(6)
                refolding_yield_dict["refolding_yield_min_mg/l"].append(6)
                refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
                refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
            elif target_df.loc[key, "Refolding Yield"] == "80% - 360mg/L culture":
                refolding_yield_dict["refolding_yield_max_%"].append(80)
                refolding_yield_dict["refolding_yield_min_%"].append(80)
                refolding_yield_dict["refolding_yield_max_mg/l"].append(360)
                refolding_yield_dict["refolding_yield_min_mg/l"].append(360)
                refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
                refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
            elif target_df.loc[key, "Refolding Yield"] == "4.4%- 0.15mg/L culture":
                refolding_yield_dict["refolding_yield_max_%"].append(4.4)
                refolding_yield_dict["refolding_yield_min_%"].append(4.4)
                refolding_yield_dict["refolding_yield_max_mg/l"].append(0.15)
                refolding_yield_dict["refolding_yield_min_mg/l"].append(0.15)
                refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
                refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
            else:
                print(str(target_df.loc[key, "Refolding Yield"]))
                raise ValueError(f"%_mg/Lcultureのうち予期せぬ値があります。")
        elif value["refolding_yield_unit"] == "25mg/L culture (fermentation), 2mg/L (shake-flask)":
            refolding_yield_dict["refolding_yield_max_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg/l"].append(2)
            refolding_yield_dict["refolding_yield_min_mg/l"].append(2)
            refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
        elif value["refolding_yield_unit"] == '0.4mg from 0.5 mg of frozen cells':
            refolding_yield_dict["refolding_yield_max_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg/l"].append(0.5)
            refolding_yield_dict["refolding_yield_min_mg/l"].append(0.4)
            refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
        elif value["refolding_yield_unit"] == '37mg/L culture - 36%':
            refolding_yield_dict["refolding_yield_max_%"].append(36)
            refolding_yield_dict["refolding_yield_min_%"].append(36)
            refolding_yield_dict["refolding_yield_max_mg/l"].append(37)
            refolding_yield_dict["refolding_yield_min_mg/l"].append(37)
            refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
        elif value["refolding_yield_unit"] == '6 mg per litre cell culture':
            refolding_yield_dict["refolding_yield_max_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg/l"].append(6)
            refolding_yield_dict["refolding_yield_min_mg/l"].append(6)
            refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
        elif value["refolding_yield_unit"] == '23 micrograms per milligram':
            refolding_yield_dict["refolding_yield_max_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg/l"].append(23 * 1000)
            refolding_yield_dict["refolding_yield_min_mg/l"].append(23 * 1000)
            refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
        elif value["refolding_yield_unit"] == "35mg/100mg IBs":
            refolding_yield_dict["refolding_yield_max_%"].append(0.35)
            refolding_yield_dict["refolding_yield_min_%"].append(0.35)
            refolding_yield_dict["refolding_yield_max_mg/l"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg/l"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
        elif value["refolding_yield_unit"] == "0.5mg/4L":
            refolding_yield_dict["refolding_yield_max_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg/l"].append(0.5 / 4)
            refolding_yield_dict["refolding_yield_min_mg/l"].append(0.5 / 4)
            refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
        elif value["refolding_yield_unit"] == '0.5mg/g packed cells':
            refolding_yield_dict["refolding_yield_max_%"].append(0.5 / 1000)
            refolding_yield_dict["refolding_yield_min_%"].append(0.5 / 1000)
            refolding_yield_dict["refolding_yield_max_mg/l"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg/l"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
        elif value["refolding_yield_unit"] == '> 20 mg gp41 fusion protein/g E. coli wet weight':
            refolding_yield_dict["refolding_yield_max_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg/l"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg/l"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg"].append(20)
            refolding_yield_dict["refolding_yield_min_mg"].append(20)
        elif value["refolding_yield_unit"] == "ug":
            refolding_yield_dict["refolding_yield_max_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg/l"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg/l"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg"].append(value["refolding_yield_max"] / 1000)
            refolding_yield_dict["refolding_yield_min_mg"].append(value["refolding_yield_min"] / 1000)
        elif value["refolding_yield_unit"] == "0.1U/100ml culture":
            refolding_yield_dict["refolding_yield_max_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg/l"].append(value["refolding_yield_max"] / 1000 * 10)
            refolding_yield_dict["refolding_yield_min_mg/l"].append(value["refolding_yield_min"] / 1000 * 10)
            refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
        elif value["refolding_yield_unit"] == "nm":
            refolding_yield_dict["refolding_yield_max_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg/l"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg/l"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
        elif value["refolding_yield_unit"] == "3mg/10g cellss":
            refolding_yield_dict["refolding_yield_max_%"].append(value["refolding_yield_max"] / 1000)
            refolding_yield_dict["refolding_yield_min_%"].append(value["refolding_yield_min"] / 1000)
            refolding_yield_dict["refolding_yield_max_mg/l"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg/l"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
        elif value["refolding_yield_unit"] == "20%, 20mg per L":
            refolding_yield_dict["refolding_yield_max_%"].append(20)
            refolding_yield_dict["refolding_yield_min_%"].append(20)
            refolding_yield_dict["refolding_yield_max_mg/l"].append(20)
            refolding_yield_dict["refolding_yield_min_mg/l"].append(20)
            refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
        elif value["refolding_yield_unit"] == "15 mg (30%)":
            refolding_yield_dict["refolding_yield_max_%"].append(30)
            refolding_yield_dict["refolding_yield_min_%"].append(30)
            refolding_yield_dict["refolding_yield_max_mg/l"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg/l"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg"].append(15)
            refolding_yield_dict["refolding_yield_min_mg"].append(15)
        elif value["refolding_yield_unit"] == "80mg prot/110mg":
            refolding_yield_dict["refolding_yield_max_%"].append(80 / 110)
            refolding_yield_dict["refolding_yield_min_%"].append(80 / 110)
            refolding_yield_dict["refolding_yield_max_mg/l"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg/l"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
        elif value["refolding_yield_unit"] == '55% (±6%)':
            refolding_yield_dict["refolding_yield_max_%"].append(61)
            refolding_yield_dict["refolding_yield_min_%"].append(49)
            refolding_yield_dict["refolding_yield_max_mg/l"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg/l"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
        elif value["refolding_yield_unit"] == "0.75 to 1 mg/l":
            refolding_yield_dict["refolding_yield_max_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg/l"].append(1)
            refolding_yield_dict["refolding_yield_min_mg/l"].append(0.75)
            refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
        elif value["refolding_yield_unit"] == "2,0%":
            refolding_yield_dict["refolding_yield_max_%"].append(2)
            refolding_yield_dict["refolding_yield_min_%"].append(2)
            refolding_yield_dict["refolding_yield_max_mg/l"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg/l"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
        elif value["refolding_yield_unit"] == "ug/l":
            refolding_yield_dict["refolding_yield_max_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg/l"].append(value["refolding_yield_max"] / 1000)
            refolding_yield_dict["refolding_yield_min_mg/l"].append(value["refolding_yield_min"] / 1000)
            refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
        elif value["refolding_yield_unit"] == " 5–10 mg per 1":
            refolding_yield_dict["refolding_yield_max_%"].append(0.10)
            refolding_yield_dict["refolding_yield_min_%"].append(0.05)
            refolding_yield_dict["refolding_yield_max_mg/l"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg/l"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
        elif value["refolding_yield_unit"] == "50 mg per 500ml":
            refolding_yield_dict["refolding_yield_max_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg/l"].append(100 / 1000)
            refolding_yield_dict["refolding_yield_min_mg/l"].append(100 / 1000)
            refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
        elif value["refolding_yield_unit"] is np.nan:
            refolding_yield_dict["refolding_yield_max_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_%"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg/l"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg/l"].append(np.nan)
            refolding_yield_dict["refolding_yield_max_mg"].append(np.nan)
            refolding_yield_dict["refolding_yield_min_mg"].append(np.nan)
        else:
            raise ValueError()

    refolding_yield_df = pd.DataFrame(refolding_yield_dict, index=target_df.index)
    target_df = pd.concat([target_df, refolding_yield_df], axis=1)
    return target_df

def extract_purity(target_value):
    extracted_value_list = re.findall("\d+(?:\.\d+)?", target_value)
    if type(extracted_value_list) == list and len(extracted_value_list) == 2:
        max_purity_value , min_purity_value = extracted_value_list[1], extracted_value_list[0]
    elif len(extracted_value_list) == 1:
        max_purity_value , min_purity_value = extracted_value_list[0], extracted_value_list[0]
    else:
        max_purity_value , min_purity_value = np.nan, np.nan
        
    return max_purity_value , min_purity_value

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

def annotate_purity(target_df):
    """
    target_df: target_df
    """
    purity_dict = {"purity_max_%": [], 
                   "purity_min_%": [], 
                  }

    for key, value in target_df.loc[: ,["purity_max", "purity_min", "purity_unit"]].iterrows():
        if value["purity_unit"] == "%":
            purity_dict["purity_max_%"].append(value["purity_max"])
            purity_dict["purity_min_%"].append(value["purity_min"])
        elif value["purity_unit"] == 'pure':
            purity_dict["purity_max_%"].append(100)
            purity_dict["purity_min_%"].append(100)
        elif value["purity_unit"] == 'single band on SDS gel':
            purity_dict["purity_max_%"].append(100)
            purity_dict["purity_min_%"].append(100)
        elif value["purity_unit"] == 'single band on SDS PAGE':
            purity_dict["purity_max_%"].append(100)
            purity_dict["purity_min_%"].append(100)
        elif value["purity_unit"] == ' 88 ± 3':
            purity_dict["purity_max_%"].append(91)
            purity_dict["purity_min_%"].append(85)
        elif value["purity_unit"] == 'mg/' or value["purity_unit"] == 'mg/l' or value["purity_unit"] == 'mg/L' or value["purity_unit"] == 'mg':
            purity_dict["purity_max_%"].append(np.nan)
            purity_dict["purity_min_%"].append(np.nan)
        elif value["purity_unit"] is np.nan:
            purity_dict["purity_max_%"].append(np.nan)
            purity_dict["purity_min_%"].append(np.nan)
        else:
            raise ValueError()

    purity_df = pd.DataFrame(purity_dict, index=target_df.index)
    target_df = pd.concat([target_df, purity_df], axis=1)
    return target_df

def fix_refolding_temperature(target_df):
    #修正
    target_df[target_df.loc[:, "Refolding Temperature"] > 50].loc[100, "Refolding Protocol"]
    target_df.loc[100, "Refolding Temperature"] = 4
    #削除
    target_df = target_df[target_df["Protein Name"] != "prion-like yeast protein"]
    return target_df


