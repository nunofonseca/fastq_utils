import subprocess
from common.validation_report import ValidationReport


class Validator:

    def validate(self, file_path):
        report = ValidationReport()

        process = subprocess.Popen(["fastq_info", file_path], stdout=subprocess.PIPE)
        process.wait()
        for line in process.stderr:
            report.log_error(line)

        report.state = "INVALID" if report.errors else "VALID"

        return report
