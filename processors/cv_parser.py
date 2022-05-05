import numpy as np
import dateutil.parser


class ParseChiCV:
    def __init__(self, file_path):
        self.file_path = file_path

        self.parse_file()

    def parse_file(self):
        with open(self.file_path, "r") as f:
            line = f.readline()
            self.date_recorded = dateutil.parser.parse(line.strip()).isoformat()

            segment = 1
            line = f.readline()
            while not line.startswith("Potential"):
                if line.startswith("File"):
                    self.file_name = self.extract_value_unit(line, value_break=":")
                elif line.startswith("Header"):
                    self.header = self.extract_value_unit(line, value_break=":")
                elif line.startswith("Note"):
                    self.note = self.extract_value_unit(line, value_break=":")
                elif line.startswith("Init E"):
                    self.init_e = self.extract_value_unit(line, value_type='float')
                elif line.startswith("High E"):
                    self.high_e = self.extract_value_unit(line, value_type='float')
                elif line.startswith("Low E"):
                    self.low_e = self.extract_value_unit(line, value_type='float')
                elif line.startswith("Init P/N"):
                    self.init_p_n = self.extract_value_unit(line)
                elif line.startswith("Scan Rate"):
                    self.scan_rate = self.extract_value_unit(line, value_type='float')
                elif line.startswith("Segment ="):
                    self.segment = self.extract_value_unit(line, value_type='int')
                elif line.startswith("Sample Interval"):
                    self.sample_interval = self.extract_value_unit(line, value_type='float')
                elif line.startswith("Quiet Time"):
                    self.quiet_time = self.extract_value_unit(line, value_type='float')
                elif line.startswith("Sensitivity"):
                    self.sensitivity = self.extract_value_unit(line, value_type='float')
                elif line.startswith("Comp R"):
                    self.comp_R = self.extract_value_unit(line, value_type='float')

                line = f.readline()

            line = f.readline()
            while not line.strip().split():
                line = f.readline()
            scan_data = []
            self.data_points_per_scan = int(abs(self.high_e["value"] - self.low_e["value"]) / self.sample_interval["value"])
            while line.strip().split():
                scan = []
                data_point = 0
                while data_point < self.data_points_per_scan:
                    if line.strip().split():
                        scan.append([float(x) for x in line.strip().split(",")])
                    data_point += 1
                    line = f.readline()
                arr_data = np.array(scan).tolist()
                scan_data.append(arr_data)
                line = f.readline()
        self.num_scans = len(scan_data)
        self.scan_data = scan_data

    @staticmethod
    def extract_value_unit(line, value_break="=", value_type='str'):
        if value_type == 'float':
            value = float("".join(line.split(value_break)[1:]).strip())
        elif value_type == 'int':
            value = int("".join(line.split(value_break)[1:]).strip())
        else:
            value = ("".join(line.split(value_break)[1:]).strip())
        return value
