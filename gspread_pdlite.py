# Light wrapper for using google sheets
# mostly using gspread and pandas
# lighter weight, less checks and validation than gspread-pandas
# To see how gspread_pandas interacts with gspread, see:
# https://github.com/burnash/gspread/blob/master/gspread/spreadsheet.py
#
# Main thing to watch for: need to sanitize data for gsheet API before saving.
# Here I provide some convenience sanitizers, but in general 
# one should be careful to do that themselves beforehand.
#
import gspread
import pandas as pd
from io import StringIO
import collections

# ===== Datatype utilities
def to_num(x_str, else_str=False):
    """cast to number

    Args:
        x_str (str): value to try to cast. preferably string.
        else_str (bool): whether to force non numbers into strings.

    Returns:
        (): number if castable, otherwise original string

    Notes: 
        don't try to fill in empty strings here
        maybe I should try something else...

        e.g. what if there are numpy or pandas ints/floats in there?
        hard to check...
    """
    number_types = (int,float,pd.np.int64,pd.pd.int32,pd.np.float64,pd.np.float32)
    if isinstance(x_str,number_types):
        return x_str
    elif isinstance(x_str,str):
        if "." in x_str: #is float
            try:
                x = float(x_str)
                return x
            except:
                if else_str:
                    return str(x_str)
                else:
                    return x_str
        else: #is int
            try:
                x = int(x_str)
                return x
            except:
                if else_str:
                    return str(x_str)
                else:
                    return x_str
    else:
        return str(x_str)

def sanitize(df,index=False):
    """Use pandas' to_csv() export to sanitize a dataframe for google API.

    Args:
        df (DataFrame): 

    Returns:
        DataFrame: sanitized data frame.

    Notes:
        Essentially turns non-number values into strings that GoogleAPI can accept.
        Alternative is to use applymap with to_num().
    """
    # remember to set na_filter = False. For some reason 
    tmp_bytes = df.to_csv(index=index).encode("utf-8") 
    # remember to set na_filter = False. For some reason 
    tmp_df = pd.read_csv( StringIO(str(tmp_bytes,"utf-8")), index_col=False, na_filter=False)
    return tmp_df


# ===== Utilities for reading public sheets, not using gspread
def urlfy(sheet_id,sheet_name=None,gid=None):
    url_base = f"https://docs.google.com/spreadsheets/d/{sheet_id}/"

    if sheet_name is not None:
        # sheet_name takes precedence over gid
        url = url_base + f"gviz/tq?tqx=out:csv&sheet={sheet_name}"
    elif gid is not None:
        url = url_base + f"export?format=csv&gid={gid}"
    else:
        url = url_base + f"export?format=csv"

    return url

def url_to_df(url):
    df = pd.read_csv(url, dtype=str).fillna("")
    return df

# usage is two lines:
#   url = urlfy(sheet_id, sheet_name)
#   df = url_to_df(url)
# or one line:
#   df = url_to_df( urlfy(sheet_id, sheet_name) )

# ===== Using gspread
# Note, NEED to have some kind of oauth certification!
# can use own certification object, or use built-in

# Note, gspread_pandas creates a custom Spread() object 
#   that itself requires a Client() object

def open_sheet(name, cert, by_title=False):
    if isinstance(cert,collections.Mapping):
        gc = gspread.service_account_from_dict(cert)
    else: #certification is a certification object
        # see https://docs.gspread.org/en/latest/oauth2.html#service-account
        gc = gspread.authorize(cert)

    if by_title:
        sheets = gc.open(name)
    elif "google.com" in name:
        sheets = gc.open_by_url(sheets)
    else:
        sheets = gc.open_by_key(name)

    return sheets

# subsequently, have to retrieve worksheets using:
#  
# `sheet.worksheet(name)`
# `sheet.get_worksheet(idx)`
# `sheet.get_worksheet_by_id(gid)`
#
# later can think about a wrapper object, like gspread_pandas

def worksheet_to_df(worksheet):
    """_summary_

    Args:
        worksheet (gspread worksheet): 

    Notes:
        ASSUME, first row is columns, rest are values
    """
    vals = worksheet.get_all_values()
    df = pd.DataFrame(vals[1:],columns=vals[0])

    #cast values to number if possible
    # alternative is to use worksheet.get_all_records(), 
    # but I think that might be an expensive 
    # and redundant list of dictionaries
    # consider using my wrapper using pands.to_csv to sanitize
    df = df.applymap(to_num, na_action="ignore")

    return df

def worksheet_append( worksheet, data ):
    """append df to a worksheet

    Args:
        worksheet (gspread worksheet): 
        df (several): 
        
    Notes:
        Series should use series.to_frame().transpose() first
        Note, this is really simple
        or, just do worksheet.append_row(Series.)

        if data is not DataFrame or Series, should santizie manually beforehand!
    """
    print("HERE")
    print(data)
    if isinstance(data,pd.DataFrame):
        # consider using my wrapper using pands.to_csv to sanitize
        #data = data.applymap( lambda x: to_num(x,else_str=True) )
        data = sanitize(data)
        worksheet.append_rows( data.values.tolist() )
    elif isinstance(data,pd.Series):
        # consider using my wrapper using pands.to_csv to sanitize
        #data = data.map( lambda x: to_num(x,else_str=True))
        data = sanitize(data.to_frame().transpose()) #casts to dataframe
        worksheet.append_rows( data.values.tolist() ) #dataframe format
    else:
        try: #treat as one row first
            worksheet.append_row( data )
        except: #should be list of lists
            worksheet.append_rows( data )

def worksheet_update( worksheet, df ):
    """Update values to worksheet.

    Args:
        worksheet (_type_): _description_
        df (DataFrame): Must be data frame, or else can just use simple update

    Notes:
        Default update() behavior doesn't clear the sheet!

        Dataframe values should be numbers or strings!
        Or else Google API will complain.
        no validation done right now. 

        One approach is to use an applymap() to sanitize the dataframe

        Alternative would be to export with csv first, 
        i.e. let pandas do the string exporting

        Shouldn't use series, otherwise just use gspread defaults.
    """
    #worksheet.update([df.columns.values.tolist()] + df.values.tolist())
    tmp_df = sanitize(df)
    tmp_data = [tmp_df.columns.values.tolist()] + tmp_df.values.tolist()
    
    worksheet.clear()
    worksheet.update(tmp_data)

