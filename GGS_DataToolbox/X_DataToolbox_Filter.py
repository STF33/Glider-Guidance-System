# =========================
# IMPORTS
# =========================

import pandas as pd

# =========================

### FUNCTION:
def x_data_filter(df, conditions, operator='AND'):

    '''
    Filter a DataFrame based on specified conditions and logical operator.
    
    Args:
    - df (pd.DataFrame): The DataFrame to filter.
    - conditions (dict): A dictionary where keys are column names and values are tuples containing an operator and a value.
    - operator (str): Logical operator to combine conditions ('AND' or 'OR').
    
    Returns:
    - df_filtered (pd.DataFrame): The filtered DataFrame.
    '''

    if operator not in ['AND', 'OR']:
        raise ValueError("Operator must be 'AND' or 'OR'")
    if operator == 'AND':
        mask = pd.Series([True] * len(df))
        for column, (op, value) in conditions.items():
            if op == '==':
                mask &= (df[column] == value)
            elif op == '!=':
                mask &= (df[column] != value)
            elif op == '>':
                mask &= (df[column] > value)
            elif op == '<':
                mask &= (df[column] < value)
            elif op == '>=':
                mask &= (df[column] >= value)
            elif op == '<=':
                mask &= (df[column] <= value)
            else:
                raise ValueError(f"Unsupported operator: {op}")
    elif operator == 'OR':
        mask = pd.Series([False] * len(df))
        for column, (op, value) in conditions.items():
            if op == '==':
                mask |= (df[column] == value)
            elif op == '!=':
                mask |= (df[column] != value)
            elif op == '>':
                mask |= (df[column] > value)
            elif op == '<':
                mask |= (df[column] < value)
            elif op == '>=':
                mask |= (df[column] >= value)
            elif op == '<=':
                mask |= (df[column] <= value)
            else:
                raise ValueError(f"Unsupported operator: {op}")
    df_filtered = df[mask]

    return df_filtered

# conditions = {
#     'm_de_oil_vol': ('<', 0),
#     'c_air_pump': ('==', 0)
# }
