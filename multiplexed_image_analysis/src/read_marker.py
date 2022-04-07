from pathlib import Path
import datatable as dt

def read_marker(path, marker):
    # Generate a list with all the files on the folder
    files = Path(path).rglob("*.csv")

    # iread returns an iterator of all files, rbind combines them together
    # pd.concat() works similar to dr.rbind()
    all_csvs = dt.rbind(dt.iread(tuple(files), columns={marker:marker, ...:None})).to_pandas()

    return all_csvs