import subprocess
from common.validation_report import ValidationReport


class Validator:

    def validate(self, file_path):
        report = ValidationReport()

        process = subprocess.Popen(["fastq_info", "-r", "-s", file_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        err_lines = stderr.decode().split('\n')
        for line in err_lines:
            if "ERROR" in line:
                report.log_error(line.rstrip())

        report.state = "INVALID" if report.errors else "VALID"

        return report
