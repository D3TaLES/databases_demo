import pandas as pd


class ParseExcel:
    def __init__(self, file_path):
        self.file_path = file_path
        self.parse_file()

    def parse_file(self):
        try:
            df = pd.read_excel(self.file_path, header=None, names=['col1', 'col2'])
        except:
            df = pd.read_csv(self.file_path, header=None, names=['col1', 'col2'])

        self.data_df = df.iloc[4:, :].astype(float, errors='ignore').rename(columns={'col1': 'wavelength', 'col2': 'absorbance'})
        self.string_data = df.iloc[:3, :]

    @property
    def integration_time(self):
        query = self.string_data[self.string_data["col1"].str.contains('Integration Time')]['col2'].values
        return query[0] if query else None

    @property
    def date_recorded(self):
        query = self.string_data[self.string_data["col1"].str.contains('Timestamp')]['col2'].values
        return query[0] if query else ''

    @property
    def absorbance_data(self):
        return self.data_df.to_dict('list')

