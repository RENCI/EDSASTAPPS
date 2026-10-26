import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


REPO_ROOT=Path(__file__).resolve().parents[1]
RUN_HARVESTER=REPO_ROOT / 'run_harvester'
sys.path.insert(0, str(RUN_HARVESTER))

from station_validation import archived_log_name, validate_station_sets


class StationValidationTests(unittest.TestCase):
    def setUp(self):
        self.index=pd.date_range('2026-10-21', periods=2, freq='h')
        self.metadata=pd.DataFrame(index=['CONTROL', 'CTIN7'])

    def test_recovered_ctin7_passes_validation(self):
        data=pd.DataFrame(
            {'CONTROL': [1.0, 2.0], 'CTIN7': [3.0, 4.0]},
            index=self.index,
        )
        report=validate_station_sets(
            data,
            self.metadata,
            requested={'CONTROL', 'CTIN7'},
            policy='fail',
        )
        self.assertFalse(report['has_mismatch'])

    def test_missing_ctin7_fails_before_csv_output(self):
        data=pd.DataFrame({'CONTROL': [1.0, 2.0]}, index=self.index)
        with self.assertRaisesRegex(RuntimeError, 'CTIN7'):
            validate_station_sets(data, self.metadata, policy='fail')

    def test_warn_policy_reports_missing_and_all_fill_stations(self):
        data=pd.DataFrame(
            {'CONTROL': [1.0, 2.0], 'CTIN7': [-99999, -99999]},
            index=self.index,
        )
        report=validate_station_sets(data, self.metadata, policy='warn')
        self.assertTrue(report['has_mismatch'])
        self.assertEqual(['CTIN7'], report['stations_without_valid_data'])

    def test_log_name_contains_job_identity(self):
        name=archived_log_name(
            '4406-2026102106-gfsforecast',
            'CONTRAILS_COASTAL',
            'FORECAST',
            'pod/1',
        )
        self.assertEqual(
            'run_harvester_adcirc_4406-2026102106-gfsforecast__CONTRAILS_COASTAL__FORECAST__pod-1.log',
            name,
        )

    def test_launcher_stops_before_nowcast_when_forecast_fails(self):
        script=RUN_HARVESTER / 'execute_multiprovider_run_harvest_adcirc.sh'
        with tempfile.TemporaryDirectory() as temp_dir:
            temp=Path(temp_dir)
            bin_dir=temp / 'bin'
            bin_dir.mkdir()
            calls=temp / 'calls.txt'
            git_stub=bin_dir / 'git'
            python_stub=bin_dir / 'python'
            git_stub.write_text('#!/bin/sh\nexit 0\n')
            python_stub.write_text('#!/bin/sh\nprintf "%s\\n" "$*" >> "$CALLS_FILE"\nexit 1\n')
            git_stub.chmod(0o755)
            python_stub.chmod(0o755)
            env=os.environ.copy()
            env.update(
                {
                    'PATH': f'{bin_dir}:{env["PATH"]}',
                    'PYTHONPATH': 'fixture',
                    'LOG_PATH': str(temp / 'logs'),
                    'CALLS_FILE': str(calls),
                }
            )
            result=subprocess.run(
                [str(script), 'forecast-url', str(temp / 'output'), 'model-id', 'CONTRAILS_COASTAL'],
                cwd=temp,
                env=env,
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertNotEqual(0, result.returncode)
            call_lines=calls.read_text().splitlines()
            self.assertEqual(1, len(call_lines))
            self.assertIn('--station_mismatch_policy fail', call_lines[0])
            self.assertNotIn('--ensemble nowcast', call_lines[0])


if __name__ == '__main__':
    unittest.main()
